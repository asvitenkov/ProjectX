#ifndef ICOMPUTEFIELD_H
#define ICOMPUTEFIELD_H

#include "Math/Vector3.hpp"
#include "Math/Point3.hpp"
#include "Math/Triangle.hpp"
#include "Math/Complex.hpp"
#include "Math/MathDefines.h"


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

    IComputeField()
    {
        mZero = std::numeric_limits<T>::epsilon();
    }

    /*
            Функция расчета поля на элементарном треугольнике
            !tr - расчитываемый треугольник
            !phi, theta - азимутальный и зенитный углы в сферических координатах в градусах
            !pniNab - пока не известно что
            !lambda -  длина волны
            !m1, e1 абсолютная комплексная проницаемость
            !type - тип поверхности треугольника
            return расчитанное поле
        */
    virtual ComplexType CalculateTriangleField(const IComputeField<T>::TriangleType &tr, T phi, T theta, T phiNab, T lambda, ComplexType m1, ComplexType e1, EConduction::TYPE type) const = 0;

protected:
    T mZero;
};

#endif // ICOMPUTEFIELD_H
