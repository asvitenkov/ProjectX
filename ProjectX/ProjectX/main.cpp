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
	Point3F point1(-1, 0, 2);
	Point3F point2( 1,-2, 5);
	Point3F point3( 3, 0, 4);

	typedef float _Ty;
	Triangle<_Ty> tr(point1, point2, point3);	//!z  - ���������� ����� 
	Vec3<_Ty> p0;				//!P0,P1 - ����������� ����������� � ���������
	Vec3<_Ty> p1;				//!P0,P1 - ����������� ����������� � ���������
	_Ty k0;						//!k0 - �������� �����
	std::complex<_Ty> e1;		//e1,m1 - ������������� ������������� 
	std::complex<_Ty> m1;		//e1,m1 - ������������� �������������
	ElCondType::TYPE tol(ElCondType::Metall); //!tol - ������� (0-������, -1 - ����������)
	Vec3<std::complex<_Ty>> ep;	// !ep - ����������� ����
	
	Calc<float> cur;
	cur.FieldOfTriangle(tr, p0, p1, k0, e1, m1, tol, ep);

	std::cout << ep;

	_getch();
	return 0;
}

