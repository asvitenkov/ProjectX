#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

#include "Vector3.hpp"
#include "Point3.hpp"

#include <iostream>
#include <math.h>
#include <limits>

template<typename T>
class Triangle
{
public:
	typedef Triangle<T> ThisType;
	typedef Point3<T>	PointType;
	typedef Vec3<T>		VecType;

	Triangle()
	{
		Assign(PointType(), PointType(), PointType());
	};

	Triangle( const ThisType& v )
	{
		Assign(v);
	};

	Triangle( PointType fx, PointType fy, PointType fz )
	{
		Assign(fx, fy, fz);
	};

	Triangle( const PointType vecArray[3] )
	{
		Assign(vecArray[0], vecArray[1], vecArray[2]);
	};

	~Triangle()
	{
		Assign(PointType(), PointType(), PointType());
	};

	friend std::ostream& operator << (std::ostream& output, const ThisType& v)
	{
		output << "Triangle(" << p1 << ", " << p2 << ", " << p3 << ")";
		return output;
	};

	ThisType operator = ( const ThisType& v )
	{
		Assign(v);
		return (*this);
	};

	ThisType operator = ( const PointType* pf )
	{
		Assign(pf[0], pf[1], pf[2]);
		return (*this);
	};

	ThisType operator - ( const ThisType& v ) const
	{
		return ThisType(p1 - v.p1, p2 - v.p2, p3 - v.p3);
	};

	ThisType operator + ( const ThisType& v ) const
	{
		return ThisType(p1 + v.p1, p2 + v.p2, p3 + v.p3);
	};

	ThisType operator * ( const T& f ) const
	{
		return ThisType(p1 * f, p2 * f, p3 * f);
	};

	ThisType operator / ( const T& f ) const
	{
		return ThisType(p1 / f, p2 / f, p3 / f);
	};
	
	ThisType operator += ( const ThisType& v )
	{
		Assign(p1 + v.p1, p2 + v.p2, p3 + v.p3);
		return (*this);
	}

	ThisType operator -= ( const ThisType& v )
	{
		Assign(p1 - v.p1, p2 - v.p2, p3 - v.p3);
		return (*this);
	};

	ThisType operator /= ( const T& f )
	{
		Assign(p1 / f, p2 / f, p3 / f);
		return (*this);
	};

	ThisType operator *= ( const T& f )
	{
		Assign(p1 * f, p2 * f, p3 * f);
		return (*this);
	};

	const PointType& operator [] ( const unsigned int index ) const
	{
		if(index == 0) return p1;
		if(index == 1) return p2;
		if(index == 2) return p3;

		assert("Out of bounds.");

		return p1;
	};

	PointType operator [] ( const unsigned int index )
	{
		if(index == 0) return p1;
		if(index == 1) return p2;
		if(index == 2) return p3;

		assert("Out of bounds.");

		return p1;
	};


	bool operator == ( const ThisType& v ) const
	{
		bool bRes 
			=  ((p1 == v.p1) || (p1 == v.p2) || (p1 == v.p3))
			&& ((p2 == v.p1) || (p2 == v.p2) || (p2 == v.p3))
			&& ((p3 == v.p1) || (p3 == v.p2) || (p3 == v.p3));

		return bRes; 
	};

	VecType Normal() const
	{
		const T nx = p1.y * p2.z - p1.z * p2.y;
		const T ny = p1.z * p2.x - p1.x * p2.z;
		const T nz = p1.x * p2.y - p1.y * p2.x;
		
		VecType norm(nx, ny, nz);
		norm.Normalize();

		return norm;
	};

	T Square() const
	{
		VecType v1(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
		VecType v2(p2.x - p3.x, p2.y - p3.y, p2.z - p3.z);

		VecType v3 = v1.CrossProduct(v2);

		const T square = v3.Length()/2;

		return square;
	};

	PointType p1, p2, p3;
private:
	inline void Assign(PointType p31, PointType p32, PointType p33)
	{
		p1 = p31;
		p2 = p32;
		p3 = p33;
	};

	inline void Assign(const ThisType& v)
	{
		Assign(v.p1, v.p2, v.p3);
	};
};

typedef Triangle<float> TriangleF;
typedef Triangle<double> TriangleD;

#endif
