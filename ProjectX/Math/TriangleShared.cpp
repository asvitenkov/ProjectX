#include "TriangleShared.h"

QVector<TPoint3> TriangleShared::Points;
QVector<TPoint3> TriangleShared::RedirectPoints;

TriangleShared::TriangleShared()
{

}

TriangleShared::TriangleShared(const TriangleShared& tr)
    : A(tr.A), B(tr.B), C(tr.C)
    , distance(tr.distance), radius(tr.radius)
    , nVector(tr.nVector), cVector(tr.cVector), bVector(tr.bVector)
    , aSystem(tr.aSystem), bSystem(tr.bSystem)
    , cSystem(tr.cSystem), mSystem(tr.mSystem)
    , center_3D(tr.center_3D), redirectedCenter_3D(tr.redirectedCenter_3D)
    , dead(tr.dead)
{

}

void TriangleShared::calcCenter()
{
    center_3D.x =(Points[A].x + Points[B].x + Points[C].x)/3;
    center_3D.y =(Points[A].y + Points[B].y + Points[C].y)/3;
    center_3D.z =(Points[A].z + Points[B].z + Points[C].z)/3;
}

bool TriangleShared::crossedTriangles(const TriangleShared& t2)
{
    if(crossedLines(Points[A], Points[B], Points[t2.A], Points[t2.B])) return true;
    if(crossedLines(Points[A], Points[C], Points[t2.A], Points[t2.B])) return true;
    if(crossedLines(Points[B], Points[C], Points[t2.A], Points[t2.B])) return true;

    if(crossedLines(Points[A], Points[B], Points[t2.A], Points[t2.C])) return true;
    if(crossedLines(Points[A], Points[C], Points[t2.A], Points[t2.C])) return true;
    if(crossedLines(Points[B], Points[C], Points[t2.A], Points[t2.C])) return true;

    if(crossedLines(Points[A], Points[B], Points[t2.C], Points[t2.B])) return true;
    if(crossedLines(Points[A], Points[C], Points[t2.C], Points[t2.B])) return true;
    if(crossedLines(Points[B], Points[C], Points[t2.C], Points[t2.B])) return true;
    return false;
}

bool TriangleShared::crossedLines(TPoint2 p11, TPoint2 p12, TPoint2 p21, TPoint2 p22)
{
    double k1 = (p11.y - p12.y)/(p11.x - p12.x);
    double k2 = (p21.y - p22.y)/(p21.x - p22.x);
    double x = (p21.y - p11.y + k1*p11.x - k2*p21.x)/(k1 - k2);
    if(((x <= qMax(p11.x,p12.x))&&(x >= qMin(p12.x,p11.x)))&&((x <= qMax(p22.x,p21.x))&&(x >= qMin(p21.x,p22.x))))
        return true;
    return false;
}


void TriangleShared::findVectors()
{
    float d = sqrt((TriangleShared::Points[A].x - center_3D.x)*(TriangleShared::Points[A].x - center_3D.x)+
                   (TriangleShared::Points[A].y - center_3D.y)*(TriangleShared::Points[A].y - center_3D.y)+
                   (TriangleShared::Points[A].z - center_3D.z)*(TriangleShared::Points[A].z - center_3D.z));

    nVector = countNormal();
    cVector.x = (TriangleShared::Points[A].x - center_3D.x)/d;
    cVector.y = (TriangleShared::Points[A].y - -center_3D.y)/d;
    cVector.z = (TriangleShared::Points[A].z - center_3D.z)/d;
    bVector = cVector.CrossProduct(nVector);
}

void TriangleShared::findDistance(const TVector& v)
{
    distance = (v.x*center_3D.x+v.y*center_3D.y+v.z*center_3D.z);
}

TVector TriangleShared::countNormal() const
{
    TVector a, b;
    a.x = TriangleShared::Points[B].x - TriangleShared::Points[A].x;
    a.y = TriangleShared::Points[B].y - TriangleShared::Points[A].y;
    a.z = TriangleShared::Points[B].z - TriangleShared::Points[A].z;
    b.x = TriangleShared::Points[C].x - TriangleShared::Points[A].x;
    b.y = TriangleShared::Points[C].y - TriangleShared::Points[A].y;
    b.z = TriangleShared::Points[C].z - TriangleShared::Points[A].z;
    return a.CrossProduct(b);
}

bool TriangleShared::crossedRadius(const TriangleShared& t2) const
{
    double d;
    d = (redirectedCenter_3D.x - t2.redirectedCenter_3D.x)*(redirectedCenter_3D.x - t2.redirectedCenter_3D.x)+
            (redirectedCenter_3D.y - t2.redirectedCenter_3D.y)*(redirectedCenter_3D.y - t2.redirectedCenter_3D.y)+
            (redirectedCenter_3D.z - t2.redirectedCenter_3D.z)*(redirectedCenter_3D.z - t2.redirectedCenter_3D.z);
    if(d < ( radius + t2.radius ))
        return true;
    return false;
}

bool TriangleShared::sameTriangle(const TriangleShared& t2)
{
    if(A==t2.A) return true;
    if(A==t2.B) return true;
    if(A==t2.C) return true;
    if(B==t2.A) return true;
    if(B==t2.B) return true;
    if(B==t2.C) return true;
    if(C==t2.A) return true;
    if(C==t2.B) return true;
    if(C==t2.C) return true;
    return false;
}

double TriangleShared::findRadius()
{
    double Radius = 0;
    double tempA=0;
    double tempB=0;
    double tempC=0;
    tempA=( RedirectPoints[A].x - redirectedCenter_3D.x)*( RedirectPoints[ A].x- redirectedCenter_3D.x)+
            ( RedirectPoints[A].y- redirectedCenter_3D.y)*( RedirectPoints[ A].y- redirectedCenter_3D.y)+
            ( RedirectPoints[A].z- redirectedCenter_3D.z)*( RedirectPoints[ A].z- redirectedCenter_3D.z);
    tempB=( RedirectPoints[B].x- redirectedCenter_3D.x)*( RedirectPoints[ B].x- redirectedCenter_3D.x)+
            ( RedirectPoints[B].y- redirectedCenter_3D.y)*( RedirectPoints[ B].y- redirectedCenter_3D.y)+
            ( RedirectPoints[B].z- redirectedCenter_3D.z)*( RedirectPoints[ B].z- redirectedCenter_3D.z);
    tempC=( RedirectPoints[C].x- redirectedCenter_3D.x)*( RedirectPoints[ C].x- redirectedCenter_3D.x)+
            ( RedirectPoints[C].y- redirectedCenter_3D.y)*( RedirectPoints[ C].y- redirectedCenter_3D.y)+
            ( RedirectPoints[C].z- redirectedCenter_3D.z)*( RedirectPoints[ C].z- redirectedCenter_3D.z);
    if(tempA>=tempB){
        if(tempA>=tempC)
            Radius = sqrt(tempA);
        else
            Radius = sqrt(tempC);
    }
    else
    {
        if(tempB>=tempC)
            Radius = sqrt(tempB);
        else
            Radius = sqrt(tempC);
    }

    radius = Radius;
}
