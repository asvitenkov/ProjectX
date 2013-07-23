#include "calculatedialog.h"
#include "ui_calculatedialog.h"


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


    QString waveNumberFormatStr = QString("<html>%1 (m<sup>-1</sup>)</html>").arg(tr("Wavenumber"));
    ui->lblWavenumber->setTextFormat(Qt::RichText);
    ui->lblWavenumber->setText(waveNumberFormatStr);



    mMaterial[EConductionToString(EConduction::Metal)] = EConduction::Metal;
    mMaterial[EConductionToString(EConduction::Dielectric)] = EConduction::Dielectric;

    QList<QString> mapKeys = mMaterial.keys();
    for(int i=0; i<mapKeys.size(); ++i)
        ui->cmbbxMaterial->addItem(mapKeys[i]);


    connect(ui->btnCalculate,SIGNAL(clicked()),this,SLOT(OnBtnCalculate()));

}

double CalculateDialog::Wavenumber()
{
    return ui->dspnbxWavenumber->value();
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


void CalculateDialog::SetWavenumber(double value)
{
    ui->dspnbxWavenumber->setValue(value);
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
