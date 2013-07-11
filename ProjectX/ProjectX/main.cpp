// main.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "Math/Vector3.hpp"
#include "Math/Complex.hpp"

#include <iostream>
#include <conio.h>

#include "Radiolocation/Calc.hpp"

int _tmain(int argc, _TCHAR* argv[])
{
	Vec3F v(3, 4, 1);

	std::cout << v;
	std::cout << v.Normal();
	
	v.Normalize();
	std::cout << v << " " << v.Revert();

	Point3F p1(-1, 0, 2);
	Point3F p2( 1,-2, 5);
	Point3F p3( 3, 0, 4);

	TriangleF tr(p1,p2,p3);

	std::cout << tr.Square();

	Calc<float> cur;
	cur.FieldOfTriangle(tr);

	_getch();
	return 0;
}

