#include "calculatedialog.h"
#include "ui_calculatedialog.h"
#include "resultplot.h"
#include "Calc/BaseComputeField.hpp"
#include "Algo/Algoritm.h"
#include "processthread.h"

QString EConductionToString(EConduction::TYPE value)
{
    if(value == EConduction::Metal)
        return "metal";
    if(value == EConduction::Dielectric)
        return "dielectric";
    return "";
}

QString complexToString(TComplex val)
{
    return QString("(%1; %2)").arg(QString::number(val.real())).arg(QString::number(val.imag()));
}

CalculateDialog::CalculateDialog(ModelView *modelView, QWidget *parent) :
    QWidget(parent), mModelView(modelView),
    ui(new Ui::CalculateDialog)
{
    ui->setupUi(this);

    InitDialog();

}

CalculateDialog::~CalculateDialog()
{
    delete ui;
}


void CalculateDialog::InitDialog()
{
    ui->grpbxInputParam->setTitle(tr("Input parametrs"));


    ui->lblWavenumber->setTextFormat(Qt::RichText);



    mMaterial[EConductionToString(EConduction::Dielectric)] = EConduction::Dielectric;
    mMaterial[EConductionToString(EConduction::Metal)] = EConduction::Metal;


    QList<QString> mapKeys = mMaterial.keys();
    for(int i=0; i<mapKeys.size(); ++i)
        ui->cmbbxMaterial->addItem(mapKeys[i]);

    SetMaterial(EConduction::Metal);

    connect(ui->btnCalculate,SIGNAL(clicked()),this,SLOT(OnBtnCalculate()));
    connect(ui->btnClear,SIGNAL(clicked()),this,SLOT(OnBtnClearLog()));



    //ui->verticalLayout_2->addWidget(new ResultPlot(0,0,1,1,0,0));

}

double CalculateDialog::Wavelength()
{
    return ui->dspnbxWavelength->value();
}

double CalculateDialog::AzimuthAngle()
{
    return ui->dspnbxAzimuthAngle->value();
}

double CalculateDialog::ZenithAngle()
{
    return ui->dspnbxZenithAngle->value();
}

TComplex CalculateDialog::M1()
{
    return TComplex(ui->dspnbxM1Real->value(),ui->dspnbxM1Img->value());
}

TComplex CalculateDialog::E1()
{
    return TComplex(ui->dspnbxE1Real->value(),ui->dspnbxE1Img->value());
}

EConduction::TYPE CalculateDialog::Material()
{
    return mMaterial[ui->cmbbxMaterial->currentText()];
}


void CalculateDialog::SetWavelength(double value)
{
    ui->dspnbxWavelength->setValue(value);
}


void CalculateDialog::SetM1(const TComplex &value)
{
    ui->dspnbxM1Real->setValue(value.real());
    ui->dspnbxM1Img->setValue(value.imag());
}

void CalculateDialog::SetE1(const TComplex &value)
{
    ui->dspnbxE1Real->setValue(value.real());
    ui->dspnbxE1Img->setValue(value.imag());
}

void CalculateDialog::SetMaterial(const EConduction::TYPE &value)
{
    ui->cmbbxMaterial->setCurrentIndex(ui->cmbbxMaterial->findText(EConductionToString(value)));
}



void CalculateDialog::OnBtnCalculate()
{
    // Process calculate


    QCursor curs = cursor();
    setCursor(Qt::WaitCursor);

    QVector<TriangleShared> triangles =  mModelView->GetModelScene()->GetModel()->GetTriangles();
    BaseComputeField<double> compute;

    TComplex result(0);

    ProcessThread pthread;
    Algoritm m_algo;
    m_algo.SetModel(mModelView->GetModelScene()->GetModel());

    TVector polVec;

    polVec.Set(ZenithAngle()*M_PI/180, AzimuthAngle()*M_PI/180, 1000);

    m_algo.SetSightVector(polVec);
    mModelView->GetModelScene()->setDrawNormal(true);
    mModelView->GetModelScene()->setNormal(polVec.Theta(), polVec.Phi());

    pthread.setAlgs(&m_algo);

    QEventLoop loop;

    connect(&pthread,SIGNAL(finished()),&loop,SLOT(quit()));
    pthread.start();

    loop.exec();


    for(double j=ui->dspnbxPhinabFrom->value(); j<=ui->dspnbxPhinabTo->value(); j+=ui->dspnbxPhinabStep->value())
    {

        for(int i=0; i< triangles.size(); ++i)
        {
            const TriangleShared &tr = triangles[i];

            if(!tr.IsVisible())
                continue;

            TTriangle sTr(tr.p1(), tr.p2(), tr.p3());
            result += compute.CalculateTriangleField(sTr, AzimuthAngle(),ZenithAngle(),j,Wavelength(),M1(),E1(),Material());

        }

        QString msg;


        msg = QString("Tetha %1, Phi %2, PhiNab %3, Wavelength %4, Material %5, E1 %6, M1 %7")
                .arg(QString::number(AzimuthAngle()))
                .arg(QString::number(ZenithAngle()))
                .arg(QString::number(j))
                .arg(QString::number(Wavelength()))
                .arg(EConductionToString(Material()))
                .arg(complexToString(E1()))
                .arg(complexToString(M1()));

        Print("===================");
        Print(msg);
        Print("Result: "+complexToString(result));
        Print("");

    }


    setCursor(curs);

}




void CalculateDialog::Print(QString text)
{
    ui->resultBrowser->append(text);
}




void CalculateDialog::OnBtnClearLog()
{
    ui->resultBrowser->clear();
}
