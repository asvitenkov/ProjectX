#ifndef __CALC_H__
#define __CALC_H__ 

#include "Math/Triangle.hpp"
#include "Math/Point3.hpp"
#include "Math/Vector3.hpp"
#include "Math/Complex.hpp"

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

	void FieldOfTriangle(const TriangleType& tr);
private:
	void Init();
};

#include "Calc.inl"

#endif //__CALC_H__