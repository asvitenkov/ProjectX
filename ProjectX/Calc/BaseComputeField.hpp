#ifndef BASECOMPUTEFIELD_H
#define BASECOMPUTEFIELD_H


#include "Calc/IComputeField.hpp"
#include <math.h>
#include <QDebug>
#include <limits.h>


#define ImOne ComplexType(1,0)
#define ImI ComplexType(0,1)
#define ImZero ComplexType(0)
#define ImTwo ComplexType(2)
#define VecToCVec(v) vector_cast<T,ComplexType>(v)

#define PrintComplex(n) qDebug() << #n << QString("%1 + %2i").arg(QString::number(n.real(),'e',3)).arg(QString::number(n.imag(),'e',3))
#define PrintComplexVector(n) \
    qDebug() << "===============================";  \
    qDebug() << "    "#n; \
    PrintComplex(n.x); \
    PrintComplex(n.y); \
    PrintComplex(n.z); \
    qDebug() << "===============================\n"

double round(double number)
{
    return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
}

template <typename T>
class BaseComputeField : public IComputeField<T>
{
    typedef IComputeField<T> BaseClass;
    typedef typename BaseClass::TriangleType TriangleType;
    typedef typename BaseClass::ComplexType ComplexType;
    typedef typename TriangleType::VecType VecType;
    typedef Vec3<ComplexType> CVecType;
    typedef Point3<T> PointType;


//    CVecType VecToCVec(const VecType& vec) const
//    {
//        CVecType cvec(ComplexType(vec.x,0),ComplexType(vec.y,0),ComplexType(vec.z,0) );

//        return cvec;
//    }

public:
    BaseComputeField()
    {

    }



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

        VecType R0(CP, -SP * SF, -SP * CF);
        VecType RR(-CP, SP * SF, SP * CF);

        T CFG, SFG;
        CFG = cos(phi + phiNab);
        SFG = sin(phi + phiNab);

        VecType R1(-CP, SP*SFG, SP*CFG);
        VecType P0(  0,     CF,    -SF);
        VecType P1(  0,   -CFG,   -SFG);


        VecType PKO(SP, CP*SF, CP*CF);
        VecType PK1(-SP, -CP*SFG, -CP*CFG);


        // k - волновое число
        T k = 2 * M_PI / lambda;

        PointType p1 = tr.p1();
        PointType p2 = tr.p2();
        PointType p3 = tr.p3();

        VecType Z1(p1.x, p1.y, p1.z);
        VecType Z2(p2.x, p2.y, p2.z);
        VecType Z3(p3.x, p3.y, p3.z);

        T Z[3][3] = {
                      Z1[0], Z1[1], Z1[2],
                      Z2[0], Z2[1], Z2[2],
                      Z3[0], Z3[1], Z3[2]
                    };

        VecType DR0, DR1, DR2;


        DR0 = R0 + R1;
        DR1 = Z2 - Z3;
        DR2 = Z1 - Z3;


        T dd = tr.Square() / 2.0;

        T a1 = 0, a2 = 0;

        a1 = (DR0*DR1).ManhattanLength();
        a2 = (DR0*DR2).ManhattanLength();

        a1 = -a1 / k;
        a2 = -a2 / k;

        VecType vNormal = tr.Normal();

        T tes = 3 * (RR*vNormal).ManhattanLength();
        VecType n(vNormal);
        if( tes >0 )
            n = n*(-1);

        tes = (P0*n).ManhattanLength();
        VecType P0t = P0 - n*tes;

        tes = (RR * n).ManhattanLength();

        T tes2= 0;

        if( qAbs(tes) > BaseClass::mZero )
            tes2 = tes;

        ComplexType c0 = sqrt(ImOne-(ImOne - tes2*tes2)/e1*m1);
        ComplexType c2 = c0;


        ComplexType c1;
        if(type == EConduction::Metal)
        {
            c1 =  ComplexType(0,0);
        }
        else if(type == EConduction::Dielectric)
        {
            // ?????
        }

        ComplexType c = sqrt(m1/e1) * c0 * (sin(c1)*cos(c1));
        c0 = 2*tes2;

        c1 = (ImI * c0 - ImOne)/(ImI * c0 + ImOne);

        VecType R0p = n.CrossProduct(RR);
        VecType R0t = RR - n*tes2;
        CVecType P1t;


        if(type == EConduction::Metal )
        {
            ComplexType tmp =  (e1/m1 * (ImI*c + c2*c2/tes2 ));
            T length = (R0p*P0).ManhattanLength();
            ComplexType cTmp(length);
            P1t = VecToCVec(P0t)*c1 + ImZero + VecToCVec(R0p) * (cTmp/tmp) ;
        }
        else if( type == EConduction::Dielectric)
        {
            ComplexType lengthR0t((R0t*P0).ManhattanLength(),0);
            ComplexType tmp = lengthR0t / (ImI * c * tes2);
            ComplexType lengthR0p((R0p*P0).ManhattanLength(),0);
            ComplexType tmp2 = e1/m1*(ImI*c+(c2*c2/tes2));
            P1t = VecToCVec(P0t)*c1 + VecToCVec(R0t)*tmp*((ImI*ImTwo*c)/(ImI*c0+ImOne)) + VecToCVec(R0p)*lengthR0p/tmp2;
        }

        ComplexType tes5 = (VecToCVec(RR)*P1t).ManhattanLength();

        CVecType Pr1 = P1t + VecToCVec(n) * tes5/ComplexType(tes);

        VecType RR1 = RR - n*tes*2;

        VecType EtTmp = n.CrossProduct(P0);

        CVecType cn = VecToCVec(n)+ImI;

        CVecType Et = VecToCVec(EtTmp) + (cn.CrossProduct(Pr1));

        tes = (RR*n).ManhattanLength();

        VecType HtTmp = P0*tes - RR*( (P0*n).ManhattanLength() );

        tes = (RR1 * n).ManhattanLength();

        CVecType cdr3 = Pr1 * tes;

        tes5 = (Pr1*VecToCVec(n)).ManhattanLength();

        cdr3 = cdr3 - VecToCVec(RR1)*tes5;

        tes = (RR1*Z3).ManhattanLength();


        CVecType Ht = VecToCVec(HtTmp) + cdr3;

        ComplexType a = (Et * VecToCVec(P1.CrossProduct(R1))).ManhattanLength();

        tes = (P1 * R1).ManhattanLength();

        VecType dr3 = R1 * tes;

        dr3 = P1 - dr3;

        a += (Ht * VecToCVec(dr3)).ManhattanLength();

        CVecType f(a,a,a);


        // Вычисление Bi
        VecType Bi0(1.0/6, 1.0/6, 1.0/2);


        // Bi1
        CVecType Bi1;
        {
            ComplexType v11 = -1/(a1*a1);
            ComplexType v12 = (ImOne/(ComplexType(a1)-ImI))/a1;
            ComplexType v13 = round(1/(a1*a1));
            ComplexType v14 = (ImOne / (a1 - ImI)) * (ImOne / (a1 - ImI));
            ComplexType v15 = (std::exp(ImI*a1)-ImOne)/(ImI*a1);
            ComplexType Bi11 = -0.5 * (  v11 - v12 + ( v13 + v14 ) * v15 );
            ComplexType Bi12(1.0/6);
            ComplexType Bi13 = ( (ImOne - ImI*a1)*(std::exp(ImI*a1)) - ImOne)/(a1*a1);

            Bi1 = CVecType(Bi11, Bi12, Bi13);
        }

        // Bi2
        CVecType Bi2;
        {

            ComplexType Bi21, Bi22, Bi23;

            ComplexType ia2 = ImI*a2;
            // Bi21
            {
                ComplexType numNum = std::exp(ia2)-ImOne;
                ComplexType numDenum = ia2;

                ComplexType num = numNum/numDenum - ImOne - ia2/ImTwo;
                ComplexType denum = a2*a2;
                Bi21 = - num/denum;
            }
            // Bi22
            {
                ComplexType numV1 = ImOne - (std::exp(ia2)-ImOne)/(ia2);
                PrintComplex(numV1);
                ComplexType numV2Num = ia2 * std::exp(ia2) - std::exp(ia2) + ImOne;
                PrintComplex(numV2Num);
                ComplexType numV2Denum = ia2;
                PrintComplex(numV2Denum);
                ComplexType numV2 = (numV2Num)/numV2Denum;
                PrintComplex(numV2);
                ComplexType num =numV1 - numV2;
                PrintComplex(num);
                ComplexType denum = a2*a2;
                PrintComplex(denum);

                Bi22 =  - num/denum;
            }
            // Bi23
            {
                ComplexType  num = (std::exp(ia2)-ImOne)/ia2 - ImOne;
                ComplexType denum = ia2;

                Bi23 = num/denum;
            }
            Bi2 = CVecType(Bi21,Bi22,Bi23);
        }

        // Bi3
        CVecType Bi3;
        {
            ComplexType Bi31, Bi32, Bi33;

            ComplexType ia1 = ImI * a1;
            // Bi31
            {
                ComplexType numV1 = ImOne - (std::exp(ia1)-ImOne)/(ia1);
                PrintComplex(numV1);
                ComplexType numV2Num = ia1 * std::exp(ia1) - std::exp(ia1) + ImOne;
                PrintComplex(numV2Num);
                ComplexType numV2Denum = ia1;
                PrintComplex(numV2Denum);
                ComplexType numV2 = (numV2Num)/numV2Denum;
                PrintComplex(numV2);
                ComplexType num =numV1 + numV2;
                PrintComplex(num);
                ComplexType denum = a1*a1;
                PrintComplex(denum);

                Bi31 =  - num/denum;
            }
            // Bi32
            {
                ComplexType numNum = std::exp(ia1)-ImOne;
                ComplexType numDenum = ia1;

                ComplexType num = numNum/numDenum - ImOne - ia1/ImTwo;
                ComplexType denum = a1*a1;
                Bi32 = - num/denum;
            }
            // Bi33
            {
                ComplexType  num = (std::exp(ia1)-ImOne)/(ia1 - ImOne);
                ComplexType denum = ia1;

                Bi33 = num/denum;
            }

            Bi3 = CVecType(Bi31,Bi32,Bi33);
        }

        // Bi4
        CVecType Bi4;
        {
            ComplexType Bi41,Bi42,Bi43;
            ComplexType ia2 = ImI * a2;
            ComplexType ia1 = ImI * a1;
            ComplexType eia2 = std::exp(ia2);
            ComplexType eia1 = std::exp(ia1);
            // Bi41
            {
                ComplexType num1 = (eia2 - ImOne) / ia2;
                PrintComplex(num1);
                ComplexType num2 = (eia1 - ImOne) / ia1;
                PrintComplex(num2);
                ComplexType num3 = (a2-a1);
                PrintComplex(num3);
                ComplexType num4 =  (ia1 * eia1 - eia1 + ImOne)/(ImI * a1*a1);
                PrintComplex(num4);
                ComplexType num = num1 - num2 - num3 * num4;
                PrintComplex(num);
                ComplexType denum = (a2-a1)*(a2-a1);

                Bi41 = -num/denum;
            }
            // Bi42
            {
                ComplexType num1 = (eia1 - ImOne) / ia1;
                PrintComplex(num1);
                ComplexType num2 = (eia2 - ImOne) / ia2;
                PrintComplex(num2);
                ComplexType num3 = (a1-a2);
                PrintComplex(num3);
                ComplexType num4 =  (ia2 * eia2 - eia2 + ImOne)/(ImI * a2*a2);
                PrintComplex(num4);
                ComplexType num = num1 - num2 - num3 * num4;
                PrintComplex(num);
                ComplexType denum = (a1-a2)*(a1-a2);

                Bi42 = -num/denum;
            }
            // Bi43
            {
                ComplexType num1 = (eia2-ImOne)/ia2;
                ComplexType num2 = (eia1-ImOne)/ia1;
                ComplexType num = num1 - num2;
                ComplexType denum = ImI*(a2-a1);

                Bi43 = num/denum;
            }
            PrintComplex(Bi41);
            PrintComplex(Bi42);
            PrintComplex(Bi43);
        }

//        PrintComplexVector(Bi1);
//        PrintComplexVector(Bi2);
//        PrintComplexVector(Bi3);

        return field;
    }
};

#endif // BASECOMPUTEFIELD_H
