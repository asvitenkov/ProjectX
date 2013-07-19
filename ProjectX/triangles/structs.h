#ifndef _STRUCTS_H
#define _STRUCTS_H

#include "Math/MathDefines.h"

class triangle_t
{
public:
    int A, B, C;
    double distance, radius;
    TPoint3 nVector,cVector,bVector,rVector;
    TPoint3 aSystem,bSystem,cSystem,mSystem;
    TPoint3 center_3D, redirectedCenter_3D;
    bool dead;
};

struct polarVector
{
	double r; //long
	double fi;//between oX and proect
	double te;//between oZ and QVector
};

struct line
{
    int A;
    int B;
};

#endif
