//#ifndef COMPUTEFIELDDLL_H
//#define COMPUTEFIELDDLL_H

//#include "IComputeField.hpp"
//#include "Math/Triangle.hpp"
//#include "libfortran/libfortran.h"

//template<typename T>
//class ComputeFieldDll : public IComputeField<T>
//{

//    typedef IComputeField<T> BaseClass;
//    typedef typename BaseClass::TriangleType TriangleType;
//    typedef typename BaseClass::ComplexType ComplexType;
//    typedef typename TriangleType::VecType VecType;

//public:
//    ComputeFieldDll()
//    {

//    }

//    ComplexType CalculateTriangleFieldCommon( const TriangleType &tr, T k0, ComplexType& e1, ComplexType &m1, EConduction::TYPE mType, T p0[3], T p1[3], T r0[3], T r1[3], T rr[3]) const
//    {

//        T z[3][3];
//        T n0[3][3];
//        T _e1[2] = { e1.real(), e1.imag() };
//        T _m1[2] = { m1.real(), m1.imag() };
//        double ep[2] = { 0.0, 0.0 };
//        T dd = 0.0;
//        T tol = (mType == EConduction::Metal)? 0: -1;
//        // P1
//        z[0][0]=tr.p1().x;
//        z[0][1]=tr.p1().y;
//        z[0][2]=tr.p1().z;

//        //P2
//        z[1][0]=tr.p2().x;
//        z[1][1]=tr.p2().y;
//        z[1][2]=tr.p2().z;

//        //P3
//        z[2][0]=tr.p3().x;
//        z[2][1]=tr.p3().y;
//        z[2][2]=tr.p3().z;


//        // TODO  сделать нормаль к вершинам треугольника
//        VecType nVec = tr.Normal();

//        for(int i=0; i<3; ++i)
//        {
////            n0[i][0] = nVec.x;
////            n0[i][1] = nVec.y;
////            n0[i][2] = nVec.z;

//            n0[i][0] = 0;
//            n0[i][1] = 0;
//            n0[i][2] = 1;
//        }


//        AARIST_DLL_COMMON(z, k0, n0, dd, _e1, _m1, tol, p0, p1, r0, r1, rr, ep);

//        ComplexType calcField = ComplexType(ep[0],ep[1]);

//        return calcField;
//    }

//    ComplexType CalculateTriangleField(const TriangleType &tr, T k0, ComplexType &e1, ComplexType &m1, EConduction::TYPE mType, T p0[3], T p1[3]) const
//    {

//        T r0[3] = {0.0, 0.0, 0.0};
//        T r1[3] = {0.0, 0.0, 0.0};
//        T rr[3] = {0.0, 0.0, 0.0};

//        return CalculateTriangleFieldCommon(tr, k0,e1,m1,mType,p0,p1,r0,r1,rr);
//    }

//    ComplexType CalculateEdgeField() const
//    {
//        return ComplexType(0,0);
//    }
//};

//#endif // COMPUTEFIELDDLL_H
