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
	typedef float _Ty;

	Point3<_Ty> point1(-1, 0, 2);
	Point3<_Ty> point2( 1,-2, 5);
	Point3<_Ty> point3( 3, 0, 4);

	Triangle<_Ty> tr(point1, point2, point3);	//!z  - координаты точек 
	Vec3<_Ty> p0;				//!P0,P1 - поляризации передатчика и приемника
	Vec3<_Ty> p1;				//!P0,P1 - поляризации передатчика и приемника
	_Ty k0;						//!k0 - волновое число
	std::complex<_Ty> e1;		//e1,m1 - относительные проницаемости 
	std::complex<_Ty> m1;		//e1,m1 - относительные проницаемости
	ElCondType::TYPE tol(ElCondType::Metall); //!tol - признак (0-металл, -1 - диэлектрик)
	Vec3<std::complex<_Ty>> ep;	// !ep - расчитанное поле
	
	Calc<_Ty> cur;
	cur.FieldOfTriangle(tr, p0, p1, k0, e1, m1, tol, ep);

	std::cout << ep;

	_getch();
	return 0;
}

