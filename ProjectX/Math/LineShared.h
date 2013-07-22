#ifndef __LineShared_H__
#define __LineShared_H__

#include <iostream>

#include <assert.h>

class LineShared
{
public:
	typedef LineShared ThisType;
	static ThisType Zero;

	LineShared()
	{
		Assign(int(0), int(0));
	};

	LineShared( const ThisType& v )
	{
		Assign(v);
	};

	LineShared( int a, int b)
	{
		Assign(a, b);
	};

	LineShared( const int xyzArray[2] )
	{
		Assign(xyzArray[0], xyzArray[1]);
	};

	~LineShared()
	{
	};

	friend std::ostream& operator << (std::ostream& output, const ThisType& v)
	{
        output << "LineShared(" << v.a << "->" << v.b << ")";
		return output;
	};

	ThisType operator = ( const ThisType& v )
	{
		Assign(v);
		return (*this);
	};

	ThisType operator = ( const int& f)
	{
		Assign(f, f);
		return (*this);
	};

	ThisType operator = ( const int* pf )
	{
		Assign(pf[0], pf[1]);
		return (*this);
	};

    int operator [] ( const unsigned int index )
	{
		if(index == 0) return a;
		if(index == 1) return b;

		assert("Out of bounds.");

        return a;
	};

    int operator [] ( const unsigned int index ) const
	{
		if(index == 0) return a;
		if(index == 1) return b;

		assert("Out of bounds.");

        return a;
	};

	bool operator == ( const ThisType& v ) const
	{
		bool bRes 
			=  ((a == v.a) && (b == v.b)) ||
				((a == v.a) && (b == v.a));

		return bRes; 
	};

	int a, b;
protected:
	inline void Assign(int fx, int fy)
	{
		a = fx;
		b = fy;
	};

	inline void Assign(const ThisType& v)
	{
        Assign(v.a, v.b);
	};
};

#endif
