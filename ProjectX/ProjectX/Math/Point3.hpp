#ifndef __POINT3_H__
#define __POINT3_H__

#include <iostream>
#include <math.h>
#include <limits>

template<typename T>
class Point3
{
public:
	typedef Point3<T> ThisType;

	Point3()
	{
		Assign(T(0), T(0), T(0));
	};

	Point3( const ThisType& v )
	{
		Assign(v);
	};

	Point3( T fx, T fy, T fz )
	{
		Assign(fx, fy, fz);
	};

	Point3( const T xyzArray[3] )
	{
		Assign(xyzArray[0], xyzArray[1], xyzArray[2]);
	};

	~Point3()
	{
		Assign(T(0), T(0), T(0));
	};

	friend std::ostream& operator << (std::ostream& output, const ThisType& v)
	{
		output << "Point(" << v.x << ", " << v.y << ", " << v.z << ")";
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
			=  (x - v.x) < std::numeric_limits<T>::epsilon()
			&& (y - v.y) < std::numeric_limits<T>::epsilon()
			&& (z - v.z) < std::numeric_limits<T>::epsilon();

		return bRes; 
	};
		
	T x, y, z;
protected:
	inline void Assign(T fx, T fy, T fz)
	{
		x = fx;
		y = fy;
		z = fz;
	};

	inline void Assign(const ThisType& v)
	{
		Assign(v.x, v.y, v.z);
	};
};

typedef Point3<float> Point3F;
typedef Point3<double> Point3D;

#endif
