#ifndef ICOMPUTEFIELD_H
#define ICOMPUTEFIELD_H

#include "Math/Vector3.hpp"
#include "Math/Point3.hpp"
#include "Math/Triangle.hpp"
#include "Math/Complex.hpp"


namespace EConduction
{
    enum TYPE
    {
         Metal
        ,Dielectric
    };
}

template<typename T>
class IComputeField
{


public:    
    typedef IComputeField<T> ThisType;
    typedef Point3<T>	PointType;
    typedef Vec3<T>		VecType;
    typedef Triangle<T> TriangleType;
    typedef std::complex<T> ComplexType;

    /*
            Функция расчета поля на элементарном треугольнике
            !tr - расчитываемый треугольник
            !k0 - волновое число
            !e1,m1 - относительные проницаемости комплексные числа
            !tol - признак (0-металл, -1 - диэлектрик)
            !P0,P1 - поляризации передатчика и приемника
            return расчитанное поле
        */
    virtual ComplexType CalculateTriangleField(const IComputeField<T>::TriangleType &tr, T k0, ComplexType &e1, ComplexType &m1, EConduction::TYPE mType, T p0[3], T p1[3]) const = 0;

    /*
        Функция расчета поля на элементарном треугольнике (вынесены общие переменные)
        ! z - координаты точек
        !k0 - волновое число
        !n0 - нормали к вершинам треугольника
        !ep - расчитанное поле
        !dd - площадь треугольника
        !e1,m1 - относительные проницаемости
        !tol - признак (0-металл, -1 - диэлектрик)
        !P0,P1 - поляризации передатчика и приемника
        !R0[3],R1[3],RR[3] - общие переменные (назначение не знаем!!!!)
    */
    virtual ComplexType CalculateTriangleFieldCommon(const IComputeField<T>::TriangleType &tr, T k0, ComplexType &e1, ComplexType &m1, EConduction::TYPE mType, T p0[3], T p1[3], T r0[3], T r1[3], T rr[3]) const = 0;

    //virtual ComplexType CalculateEdgeField() const = 0;

};

#endif // ICOMPUTEFIELD_H
