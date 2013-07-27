#ifndef CALCULATEDIALOG_H
#define CALCULATEDIALOG_H

#include <QWidget>
#include "Calc/IComputeField.hpp"
#include "Math/MathDefines.h"
#include <QMap>

namespace Ui {
class CalculateDialog;
}

class CalculateDialog : public QWidget
{
    Q_OBJECT
    
public:
    explicit CalculateDialog(QWidget *parent = 0);
    ~CalculateDialog();
    
    double Wavelength();
    double ZenithAngle();
    double AzimuthAngle();
    double PhiNab();
    TComplex E1();
    TComplex M1();
    EConduction::TYPE Material();

    void SetWavelength(double value);
    void SetE1(const TComplex &value);
    void SetM1(const TComplex &value);
    void SetMaterial(const EConduction::TYPE &value);




private:
    void Print(QString text);
    void InitDialog();

    Ui::CalculateDialog *ui;
    QMap<QString,EConduction::TYPE> mMaterial;

protected slots:
    void OnBtnCalculate();
};

#endif // CALCULATEDIALOG_H
