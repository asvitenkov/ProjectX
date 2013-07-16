#ifndef __VEC3_H__
#define __VEC3_H__

#include "Point3.hpp"
#include "Complex.hpp"

#include <iostream>
#include <math.h>
#include <limits>

template<typename T>
class Vec3 : public Point3<T>
{
public:
	typedef Vec3<T> ThisType;
	typedef Point3<T> Point3Type;
	typedef std::complex<T> ComplexType;

	Vec3()
	{
		Assign(T(0), T(0), T(0));
	};

	Vec3( const ThisType& v )
	{
		Assign(v.x, v.y, v.z);
	};

	Vec3( const Point3Type& p )
	{
		Assign(p);
	};

	Vec3( T fx, T fy, T fz )
	{
		Assign(fx, fy, fz);
	};

	Vec3( const T xyzArray[3] )
	{
		Assign(xyzArray[0], xyzArray[1], xyzArray[2]);
	};

	~Vec3()
	{
		Assign(T(0), T(0), T(0));
	};

	friend std::ostream& operator << (std::ostream& output, const ThisType& v)
	{
		output << "Vector(" << v.x << ", " << v.y << ", " << v.z << ")";
		return output;
	};

	ThisType operator = ( const ThisType& v )
	{
		Assign(v);
		return (*this);
	};

	ThisType operator = ( const T& f)
	{
		Assign(f, f, f);
		return (*this);
	};

	ThisType operator = ( const T* pf )
	{
		Assign(pf[0], pf[1], pf[2]);
		return (*this);
	};

	ThisType operator - ( const ThisType& v ) const
	{
		return ThisType(x - v.x, y - v.y, z - v.z);
	};

	ThisType operator + ( const ThisType& v ) const
	{
		return ThisType(x + v.x, y + v.y, z + v.z);
	};

	ThisType operator * ( const ThisType& v ) const
	{
		return ThisType(x * v.x, y * v.y, z * v.z);
	};

	ThisType operator / ( const ThisType& v ) const
	{
		return ThisType(x / v.x, y / v.y, z / v.z);
	};

	ThisType operator - ( const T& f ) const
	{
		return ThisType(x - v.x, y - v.y, z - v.z);
	};

	ThisType operator + ( const T& f ) const
	{
		return ThisType(x + v.x, y + v.y, z + v.z);
	};

	ThisType operator * ( const T& f ) const
	{
		return ThisType(x * f, y * f, z * f);
	};

	ThisType operator / ( const T& f ) const
	{
		return ThisType(x / f, y / f, z / f);
	};
	
	ThisType operator += ( const ThisType& v )
	{
		Assign(x + v.x, y + v.y, z + v.z);
		return (*this);
	}

	ThisType operator -= ( const ThisType& v )
	{
		Assign(x - v.x, y - v.y, z - v.z);
		return (*this);
	};

	ThisType operator /= ( const T& f )
	{
		Assign(x / f, y / f, z / f);
		return (*this);
	};

	ThisType operator *= ( const T& f )
	{
		Assign(x * f, y * f, z * f);
		return (*this);
	};

	bool operator == ( const ThisType& v ) const
	{
		bool bRes 
			=  abs(x - v.x) < std::numeric_limits<T>::epsilon()
			&& abs(y - v.y) < std::numeric_limits<T>::epsilon()
			&& abs(z - v.z) < std::numeric_limits<T>::epsilon();

		return bRes; 
	};

	ThisType	Revert()
	{
		return ThisType(-x,-y,-z);
	}

	ThisType	Normalize()
	{
		const T d = Length();

		if( d > std::numeric_limits<T>::epsilon() )
		{
			(*this) /= d;
		}
		else
		{
			(*this) = T(0);
		}

		return (*this);
	}

	ThisType	Normal() const
	{
		ThisType ret(*this);
		ret.Normalize();
		return ret;
	};

	T	LengthSquared() const
	{
		const T lengthSquared = x*x + y*y+ z*z;
		return lengthSquared;
	};

	T	Length() const
	{
		const T length = sqrt(LengthSquared());
		return length;
	};

	T	DotProduct( const ThisType& v ) const
	{
		const T ret = x * v.x + y * v.y + z * v.z;
		return ret;
	};

	ThisType CrossProduct( const ThisType& v ) const
	{
		const T fx = y * v.z - z * v.y;
		const T fy = z * v.x - x * v.z;
		const T fz = x * v.y - y * v.x;

		return ThisType(fx, fy, fz);
	}
};

template<typename From, typename To>
Vec3<To> vector_cast(const Vec3<From>& v)
{
	const Vec3<To> toVector(To(v.x), To(v.y), To(v.z));

	return toVector;
}

typedef Vec3<float> Vec3F;
typedef Vec3<double> Vec3D;

#endif
