#ifndef BASECOMPUTEFIELD_H
#define BASECOMPUTEFIELD_H


#include "Calc/IComputeField.hpp"
#include <math.h>
#include <QDebug>


template <typename T>
class BaseComputeField : public IComputeField<T>
{
    typedef IComputeField<T> BaseClass;
    typedef typename BaseClass::TriangleType TriangleType;
    typedef typename BaseClass::ComplexType ComplexType;
    typedef typename TriangleType::VecType VecType;

public:
    BaseComputeField(){}

    ComplexType CalculateTriangleField(const TriangleType &tr, T phi, T theta, T phiNab, T lambda, ComplexType m1, ComplexType e1, EConduction::TYPE type) const
    {
        ComplexType field(0,0);

        // переводим градусы к радианам
        theta = theta * M_PI / 180.0;
        phi = phi * M_PI / 180.0;
        phiNab = phiNab * M_PI / 180.0;

        // определение направлений облучения и приёма и поляризации
        T CP, SP, CF, SF;
        CP = cos(theta);
        SP = sin(theta);
        CF = cos(phi);
        SF = sin(phi);

        T R0[3], RR[3];

        R0[0] =  CP;
        R0[1] = -SP * SF;
        R0[2] = -SP * CF;

        RR[0] = -CP;
        RR[1] =  SP * SF;
        RR[2] =  SP * CF;

        T CFG, SFG;
        CFG = cos(phi + phiNab);
        SFG = sin(phi + phiNab);

        T R1[3];
        R1[0] = -CP;
        R1[1] =  SP * SFG;
        R1[2] =  SP * CFG ;

        T P0[3], P1[3];
        P0[0] =   0;
        P0[1] =  CF;
        P0[2] = -SF;

        P1[0] =    0;
        P1[1] = -CFG;
        P1[2] = -SFG;


        T PK0[3], PK1[3];
        PK0[0] = SP;
        PK0[1] = CP * SF;
        PK0[2] = CP * CF;

        PK1[0] = -SP;
        PK1[1] = -CP * SFG;
        PK1[2] = -CP * CFG;

        // k - волновое число
        T k = 2 * M_PI / lambda;

        T Z1[3], Z2[3], Z3[3];

        Z1[0] = tr.p1().x;
        Z1[1] = tr.p1().y;
        Z1[2] = tr.p1().z;

        Z2[0] = tr.p2().x;
        Z2[1] = tr.p2().y;
        Z2[2] = tr.p2().z;

        Z3[0] = tr.p3().x;
        Z3[1] = tr.p3().y;
        Z3[2] = tr.p3().z;

        T Z[3][3] = {
                      Z1[0], Z1[1], Z1[2],
                      Z2[0], Z2[1], Z2[2],
                      Z3[0], Z3[1], Z3[2]
                    };

        T DR0[3], DR1[3], DR2[3];

        for(int i=0; i<3; i++)
        {
            DR0[i] = R0[i] + R1[i];
            DR1[i] = Z2[i] - Z3[i];
            DR2[i] = Z1[i] - Z3[i];
        }

        T dd = tr.Square() / 2.0;

        T a1 = 0, a2 = 0;

        for(int i=0; i<3; i++)
        {
            a1+=DR0[i]*DR1[i];
            a2+=DR0[i]*DR2[i];
        }

        a1 = -a1 / k;
        a2 = -a2 / k;

        //qDebug() << a1;
        //qDebug() << a2;

        return field;
    }
};

#endif // BASECOMPUTEFIELD_H
