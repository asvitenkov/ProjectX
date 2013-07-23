#ifndef MATHDEFINES_H
#define MATHDEFINES_H

#include "Complex.hpp"
#include "Point3.hpp"
#include "Triangle.hpp"
#include "Vector3.hpp"
#include "LineShared.h"

typedef double _T_;

typedef std::complex<_T_> TComplex;
typedef Point3<_T_> TPoint3;
typedef Point3<_T_> TPoint2;
typedef Triangle<_T_> TTriangle;
typedef Vec3<_T_> TVector;
typedef Vec3<TComplex> TCVector;

namespace MathHelper
{
    inline TPoint3 Redirect(const TPoint3* p1, const TPoint3* p2)
    {
        _T_ scalar = p1->x * p2->x + p1->y * p2->y + p1->z * p2->z;

        return (*p1) - (*p2)*(scalar);
    }
}

#endif // MATHDEFINES_H
