// rcs_test.cpp : main project file.

#include "stdafx.h"

#include <windows.h>
#include <iostream>

using namespace System;

extern "C"
{
	/*
		‘ункци€ расчета пол€ на элементарном треугольнике
		! z - координаты точек 
		!k0 - волновое число
		!n0 - нормали к вершинам треугольника
		!ep - расчитанное поле
		!dd - площадь треугольника
		!e1,m1 - относительные проницаемости
		!tol - признак (0-металл, -1 - диэлектрик)
		!P0,P1 - пол€ризации передатчика и приемника
	*/
	int __declspec(dllimport) __stdcall AARIST_DLL(float z[3][3],float & k0, float n0[3][3],
		float & dd, float e1[2], float m1[2], float & tol, 
		float P0[3], float P1[3], double ep[2]);
	/*
		‘ункци€ расчета пол€ на элементарном треугольнике (вынесены общие переменные)
		! z - координаты точек 
		!k0 - волновое число
		!n0 - нормали к вершинам треугольника
		!ep - расчитанное поле
		!dd - площадь треугольника
		!e1,m1 - относительные проницаемости
		!tol - признак (0-металл, -1 - диэлектрик)
		!P0,P1 - пол€ризации передатчика и приемника
		!R0[3],R1[3],RR[3] - общие переменные (назначение не знаем!!!!)
	*/
	int __declspec(dllimport) __stdcall AARIST_DLL_COMMON(float z[3][3],float & k0, float n0[3][3],
		float & dd, float e1[2], float m1[2], float & tol, 
		float P0[3], float P1[3], float R0[3], float R1[3], float RR[3], double ep[2]);

	/*
		‘ункци€ расчета пол€, рассе€нного кромкой
		!P0, P1 - пол€ризации передатчика и приемника
		!lym - длина волны
		!ev - поле, рассе€нное кромкой
		!epr1 - Ёѕ– кромки
		!ei0,mi0 - относительные проницаемости материала покрыти€ кромки
		!tol0 - радиус покрыти€ на кромках (здесь равен 0)
		!tet00 - угол поворота кромки (здесь pi/2)
		!tau0 - параметр раствора кромки (здесь 3/2)
		!acr0 - радиус цилиндра
		!dleen0 - длина цилиндра
	*/
	int __declspec(dllimport) __stdcall POINTBLE_DLL(float P0[3],float P1[3], float & lym,
		float ev[2], float ei0[2], float mi0[2], float & tol0, float & tet00, float & tau0,
		float & acr0, float & dleen0, float & epr1);
	/*
		‘ункци€ расчета пол€, рассе€нного кромкой (вынесены общие переменные)
		!P0, P1 - пол€ризации передатчика и приемника
		!lym - длина волны
		!ev - поле, рассе€нное кромкой
		!epr1 - Ёѕ– кромки
		!ei0,mi0 - относительные проницаемости материала покрыти€ кромки
		!tol0 - радиус покрыти€ на кромках (здесь равен 0)
		!tet00 - угол поворота кромки (здесь pi/2)
		!tau0 - параметр раствора кромки (здесь 3/2)
		!acr0 - радиус цилиндра
		!dleen0 - длина цилиндра
		!Rr0[3], Rr1[3], index, pola, polb - общие переменные (назначение не знаем!!!!)
	*/
	int __declspec(dllimport) __stdcall POINTBLE_DLL_COMMON(float P0[3],float P1[3], float & lym,
		float ev[2], float ei0[2], float mi0[2], float & tol0, float & tet00, float & tau0,
		float & acr0, float & dleen0, float Rr0[3], float Rr1[3], int & index,
		float & pola, float & polb, float & epr1);

}

int main(array<System::String ^> ^args)
{
	float z[3][3], k0 = 0.0, n0[3][3], dd = 0.0, e1[2] = {0.0};
	float m1[2] = {0.0}, tol = 0.0, P0[3] = {0.0}, P1[3] = {0.0};
	float R0[3] = {0.0}, R1[3] = {0.0}, RR[3] = {0.0};

	float Rr0[3] = {0.0}, Rr1[3] = {0.0}, pola = 0.0, polb = 0.0;
	int index = 0;
	double ep[2] = {0.0};
	for (int i = 0; i<3; i++)
		for (int j = 0; j<3; j++)
		{
			z[i][j] = 1.0;
			n0[i][j] = 1.0;
		}
	float lym = 0.0, epr1 = 0.0, tol0 = 0.0; 
	float tet00 = 0.0, tau0 = 0.0, acr0 = 0.0, dleen0 = 0.0;
	float ev[2] = {0.0}, ei0[2] = {0.0}, mi0[2] = {0.0};


	AARIST_DLL(z, k0, n0, dd, e1, m1, tol, P0, P1, ep);

	Console::WriteLine("AARIST_DLL: \n"+"ep = "+Convert::ToString(ep[0]));

	AARIST_DLL_COMMON(z, k0, n0, dd, e1, m1, tol, P0, P1, R0, R1, RR, ep);
	Console::WriteLine("AARIST_DLL_COMMON: \n"+"ep = "+Convert::ToString(ep[0]));
	Console::WriteLine("R0 = "+Convert::ToString(R0[0]));
	Console::WriteLine("R1 = "+Convert::ToString(R1[0]));
	Console::WriteLine("RR = "+Convert::ToString(RR[0]));

	POINTBLE_DLL(P0, P1, lym, ev, ei0, mi0, tol0, tet00, tau0, acr0, dleen0, epr1);
	Console::WriteLine("POINTBLE_DLL: \n"+"epr1 = "+Convert::ToString(epr1));

	POINTBLE_DLL_COMMON(P0, P1, lym, ev, ei0, mi0, tol0, tet00, tau0, acr0, dleen0, 
		Rr0, Rr1, index, pola, polb, epr1);
	Console::WriteLine("POINTBLE_DLL_COMMON: \n"+"epr1 = "+Convert::ToString(epr1));
	Console::WriteLine("Rr0 = "+Convert::ToString(Rr0[0]));
	Console::WriteLine("Rr1 = "+Convert::ToString(Rr1[0]));
	Console::WriteLine("index = "+Convert::ToString(index));
	Console::WriteLine("Pola = "+Convert::ToString(pola));
	Console::WriteLine("Polb = "+Convert::ToString(polb));
	
	system("pause");
    return 0;
}
