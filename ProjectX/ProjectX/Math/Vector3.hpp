#ifndef __VEC3_H__
#define __VEC3_H__

#include <iostream>
#include <math.h>
#include <limits>

template<typename T>
class Vec3
{
public:
	typedef Vec3<T> ThisType;

	Vec3()
	{
		Assign(T(0), T(0), T(0));
	};

	Vec3( const ThisType& v )
	{
		Assign(v);
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
		output << "(" << v.x << ", " << v.y << ", " << v.z << ")";
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
		return ThisType(x * v.x, y * v.y, z * v.z);
	};

	ThisType operator / ( const T& f ) const
	{
		return ThisType(x / v.x, y / v.y, z / v.z);
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
		Assign(x * v.x, y * v.y, z * v.z);
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

	T	Length() const
	{
		return T(sqrt( x*x + y*y+ z*z ));
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
		
	T x, y, z;
private:
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

typedef Vec3<float> Vec3F;
typedef Vec3<double> Vec3D;

#endif
