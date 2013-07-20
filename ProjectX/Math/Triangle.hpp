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
        Assign(PointType::Zero, PointType::Zero, PointType::Zero);
	};

	Triangle( const ThisType& v )
	{
        Assign(v.p1(), v.p2(), v.p3());
	};

    Triangle( const PointType& fx, const PointType& fy, const PointType& fz )
	{
		Assign(fx, fy, fz);
	};

	Triangle( const PointType vecArray[3] )
	{
		Assign(vecArray[0], vecArray[1], vecArray[2]);
	};

	~Triangle()
	{

	};

    void Set( const PointType& fx, const PointType& fy, const PointType& fz )
    {
        Assign(fx, fy, fz);
    };

    friend std::ostream& operator << (std::ostream& output, const ThisType& v)
	{
        output << "Triangle(" << v.p1() << ", " << v.p2() << ", " << v.p3() << ")";
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
                =  ((p1() == v.p1()) || (p1() == v.p2()) || (p1() == v.p3()))
            && ((p2() == v.p1()) || (p2() == v.p2()) || (p2() == v.p3()))
            && ((p3() == v.p1()) || (p3() == v.p2()) || (p3() == v.p3()));

        return bRes;
	};

    bool HasSamePoint(const ThisType& v) const
    {
        bool bRes
                =  ((p1() == v.p1()) || (p1() == v.p2()) || (p1() == v.p3()))
            || ((p2() == v.p1()) || (p2() == v.p2()) || (p2() == v.p3()))
            || ((p3() == v.p1()) || (p3() == v.p2()) || (p3() == v.p3()));

        return true;
    }

	VecType Normal() const
	{
		return m_normal;
	};

	T Square() const
	{
		return m_square;
	};

#define DECL(p) private: \
    PointType m_##p; \
 	public: \
    const PointType& p () const { return m_##p;} const \
    PointType& r##p() const { return m_##p;}  const \
    void Set##p(const PointType& p_t ) \
    { if(p_t==m_##p) return; m_##p = p_t; Assign(m_p1, m_p2, m_p3); };

    DECL(p1)
    DECL(p2)
    DECL(p3)
private:
	inline void Assign(const PointType& p31, const PointType& p32, const PointType& p33)
	{
		m_p1 = p31;
		m_p2 = p32;
		m_p3 = p33;

		Recalc();
	};

	inline void Assign(const ThisType& v)
	{
        Assign(v.p1(), v.p2(), v.p3());
	};

	void Recalc()
	{
		CalcNormal();
		CalcSquare();
		CalcCenter();
	}

	void CalcNormal()
	{
		m_normal.x = m_p1.y * m_p2.z - m_p1.z * m_p2.y;
		m_normal.y = m_p1.z * m_p2.x - m_p1.x * m_p2.z;
		m_normal.z = m_p1.x * m_p2.y - m_p1.y * m_p2.x;

		m_normal.Normalize();
	}

	void CalcCenter()
	{
		m_center.x = (m_p1.x + m_p2.x + m_p3.x)/3;
		m_center.y = (m_p1.y + m_p2.y + m_p3.y)/3;
		m_center.z = (m_p1.z + m_p2.z + m_p3.z)/3;
	}

	void CalcSquare()
	{
		const VecType v1(m_p1.x - m_p2.x, m_p1.y - m_p2.y, m_p1.z - m_p2.z);
		const VecType v2(m_p2.x - m_p3.x, m_p2.y - m_p3.y, m_p2.z - m_p3.z);

		const VecType v3 = v1.CrossProduct(v2);

		m_square = v3.Length()/2;
	}

	PointType	m_center;
	VecType		m_normal;
	T			m_square;
};

typedef Triangle<float> TriangleF;
typedef Triangle<double> TriangleD;

#endif
