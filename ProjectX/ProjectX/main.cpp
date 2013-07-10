// main.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "Math/Vector3.hpp"
#include "Math/Complex.hpp"

#include <iostream>
#include <conio.h>

int _tmain(int argc, _TCHAR* argv[])
{
	Vec3F v(3, 4, 1);

	std::cout << v;
	std::cout << v.Normal();
	
	v.Normalize();
	std::cout << v << " " << v.Revert();

	_getch();
	return 0;
}

