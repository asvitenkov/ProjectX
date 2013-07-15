#ifndef __CALC_H__
#define __CALC_H__ 

#include "Math/Triangle.hpp"
#include "Math/Point3.hpp"
#include "Math/Vector3.hpp"
#include "Math/Complex.hpp"

// Электрическая проводимость
namespace ElCondType
{
	enum TYPE
	{
		Metall = 0,
		Dielectric = -1
	};
}

template<typename T>
class Calc
{
public:
	Calc();
	~Calc();

	typedef std::complex<T>	ComplexType;
	typedef Vec3<T>		VecType;
	typedef Point3<T>	PointType;
	typedef Triangle<T> TriangleType;
	typedef Vec3<ComplexType> CVecType;

	void FieldOfTriangle(	
		const TriangleType& tr,		//!z  - координаты точек 
		const VecType& p0,			//!P0,P1 - поляризации передатчика и приемника
		const VecType& p1,			//!P0,P1 - поляризации передатчика и приемника
		const T& k0,				//!k0 - волновое число
		const ComplexType& e1,		//e1,m1 - относительные проницаемости 
		const ComplexType& m1,		//e1,m1 - относительные проницаемости
		const ElCondType::TYPE tol, //!tol - признак (0-металл, -1 - диэлектрик)
		CVecType& ep				// !ep - расчитанное поле
);
private:
	void Init();
};

#include "Calc.inl"

#endif //__CALC_H__