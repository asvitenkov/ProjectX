#include "calculatedialog.h"
#include "ui_calculatedialog.h"
#include "resultplot.h"


QString EConductionToString(EConduction::TYPE value)
{
    if(value == EConduction::Metal)
        return "metal";
    if(value == EConduction::Dielectric)
        return "dielectric";
    return "";
}

CalculateDialog::CalculateDialog(QWidget *parent) :
    QWidget(parent),
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



    mMaterial[EConductionToString(EConduction::Metal)] = EConduction::Metal;
    mMaterial[EConductionToString(EConduction::Dielectric)] = EConduction::Dielectric;

    QList<QString> mapKeys = mMaterial.keys();
    for(int i=0; i<mapKeys.size(); ++i)
        ui->cmbbxMaterial->addItem(mapKeys[i]);


    connect(ui->btnCalculate,SIGNAL(clicked()),this,SLOT(OnBtnCalculate()));



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


double CalculateDialog::PhiNab()
{
    return ui->dspnbxPhinab->value();
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

}




void CalculateDialog::Print(QString text)
{
    ui->resultBrowser->append(text);
}
