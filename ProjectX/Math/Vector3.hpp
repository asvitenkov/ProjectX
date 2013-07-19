#ifndef __VEC3_H__
#define __VEC3_H__

#include "Point3.hpp"

#include <iostream>
#include <math.h>
#include <limits>

template<typename T>
class Vec3 : public Point3<T>
{
public:
	typedef Vec3<T> ThisType;
	typedef Point3<T> Point3Type;

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
        return ThisType(this->x - v.x, this->y - v.y, this->z - v.z);
	};

	ThisType operator + ( const ThisType& v ) const
	{
        return ThisType(this->x + v.x, this->y + v.y, this->z + v.z);
	};

	ThisType operator * ( const ThisType& v ) const
	{
        return ThisType(this->x * v.x, this->y * v.y, this->z * v.z);
	};

	ThisType operator / ( const ThisType& v ) const
	{
        return ThisType(this->x / v.x, this->y / v.y, this->z / v.z);
	};

	ThisType operator - ( const T& f ) const
	{
        return ThisType(this->x - f, this->y - f, this->z - f);
	};

	ThisType operator + ( const T& f ) const
	{
        return ThisType(this->x + f, this->y + f, this->z + f);
	};

	ThisType operator * ( const T& f ) const
	{
        return ThisType(this->x * f, this->y * f, this->z * f);
	};

	ThisType operator / ( const T& f ) const
	{
        return ThisType(this->x / f, this->y / f, this->z / f);
	};
	
	ThisType operator += ( const ThisType& v )
	{
        Assign(this->x + v.x, this->y + v.y, this->z + v.z);
		return (*this);
	}

	ThisType operator -= ( const ThisType& v )
	{
        Assign(this->x - v.x, this->y - v.y, this->z - v.z);
		return (*this);
	};

	ThisType operator /= ( const T& f )
	{
        Assign(this->x / f, this->y / f, this->z / f);
		return (*this);
	};

	ThisType operator *= ( const T& f )
	{
        Assign(this->x * f, this->y * f, this->z * f);
		return (*this);
	};

	bool operator == ( const ThisType& v ) const
	{
		bool bRes 
            =  abs(this->x - v.x) < std::numeric_limits<T>::epsilon()
            && abs(this->y - v.y) < std::numeric_limits<T>::epsilon()
            && abs(this->z - v.z) < std::numeric_limits<T>::epsilon();

		return bRes; 
	};

	ThisType	Revert()
	{
        return ThisType(-this->x,-this->y,-this->z);
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
        const T lengthSquared = this->x*this->x + this->y*this->y+ this->z*this->z;
		return lengthSquared;
	};

	T	Length() const
	{
		const T length = sqrt(LengthSquared());
		return length;
	};

	T	DotProduct( const ThisType& v ) const
	{
        const T ret = this->x * v.x + this->y * v.y + this->z * v.z;
		return ret;
	};

	ThisType CrossProduct( const ThisType& v ) const
	{
        const T fx = this->y * v.z - this->z * v.y;
        const T fy = this->z * v.x - this->x * v.z;
        const T fz = this->x * v.y - this->y * v.x;

		return ThisType(fx, fy, fz);
	}

	T Theta() const
	{
        return atan2(sqrt(this->x*this->x+this->y*this->y), this->z);
	}
	
	T Phi() const
	{
        return atan2(this->y, this->x);
	}

	void Set(const T theta, const T phi, const T length)
	{
		const T xx = length*cos(phi)*sin(theta);
		const T yy = length*sin(phi)*sin(theta);
		const T zz = cos(theta);

		Assign(xx, yy, zz);
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
