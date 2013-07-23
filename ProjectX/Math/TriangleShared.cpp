#include "TriangleShared.h"

QVector<TPoint2> TriangleShared::Points2D;
QVector<TPoint3> TriangleShared::Points3D;
QVector<TPoint3> TriangleShared::RedirectPoints;

TriangleShared::TriangleShared()
    : dead(false)
{

}

TriangleShared::TriangleShared(const TriangleShared& tr)
    : A(tr.A), B(tr.B), C(tr.C)
    , distance(tr.distance), radius(tr.radius)
    , m_Center3D(tr.m_Center3D), m_RedirectedCenter3D(tr.m_RedirectedCenter3D)
    , dead(tr.dead)
{

}

void TriangleShared::Set(int _A, int _B, int _C)
{
    A = _A;
    B = _B;
    C = _C;
}

void TriangleShared::CaclFileds(const TVector& v)
{
    m_Center3D.x =(p1().x + p2().x + p3().x)/3;
    m_Center3D.y =(p1().y + p2().y + p3().y)/3;
    m_Center3D.z =(p1().z + p2().z + p3().z)/3;

    m_RedirectedCenter3D = MathHelper::Redirect(&m_Center3D, &v);
}

bool TriangleShared::crossedTriangles(const TriangleShared& t2)
{
    if(crossedLines(Points2D[A], Points2D[B], Points2D[t2.A], Points2D[t2.B])) return true;
    if(crossedLines(Points2D[A], Points2D[C], Points2D[t2.A], Points2D[t2.B])) return true;
    if(crossedLines(Points2D[B], Points2D[C], Points2D[t2.A], Points2D[t2.B])) return true;

    if(crossedLines(Points2D[A], Points2D[B], Points2D[t2.A], Points2D[t2.C])) return true;
    if(crossedLines(Points2D[A], Points2D[C], Points2D[t2.A], Points2D[t2.C])) return true;
    if(crossedLines(Points2D[B], Points2D[C], Points2D[t2.A], Points2D[t2.C])) return true;

    if(crossedLines(Points2D[A], Points2D[B], Points2D[t2.C], Points2D[t2.B])) return true;
    if(crossedLines(Points2D[A], Points2D[C], Points2D[t2.C], Points2D[t2.B])) return true;
    if(crossedLines(Points2D[B], Points2D[C], Points2D[t2.C], Points2D[t2.B])) return true;
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

void TriangleShared::findDistance(const TVector& v)
{
    distance = (v.x*m_Center3D.x+v.y*m_Center3D.y+v.z*m_Center3D.z);
}

TVector TriangleShared::countNormal() const
{
    TVector a, b;
    a.x = TriangleShared::Points2D[B].x - TriangleShared::Points2D[A].x;
    a.y = TriangleShared::Points2D[B].y - TriangleShared::Points2D[A].y;
    a.z = TriangleShared::Points2D[B].z - TriangleShared::Points2D[A].z;
    b.x = TriangleShared::Points2D[C].x - TriangleShared::Points2D[A].x;
    b.y = TriangleShared::Points2D[C].y - TriangleShared::Points2D[A].y;
    b.z = TriangleShared::Points2D[C].z - TriangleShared::Points2D[A].z;
    return a.CrossProduct(b);
}

bool TriangleShared::crossedRadius(const TriangleShared& t2) const
{
    double d;
    d = (m_RedirectedCenter3D.x - t2.m_RedirectedCenter3D.x)*(m_RedirectedCenter3D.x - t2.m_RedirectedCenter3D.x)+
            (m_RedirectedCenter3D.y - t2.m_RedirectedCenter3D.y)*(m_RedirectedCenter3D.y - t2.m_RedirectedCenter3D.y)+
            (m_RedirectedCenter3D.z - t2.m_RedirectedCenter3D.z)*(m_RedirectedCenter3D.z - t2.m_RedirectedCenter3D.z);

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
    tempA=( RedirectPoints[A].x - m_RedirectedCenter3D.x)*( RedirectPoints[ A].x- m_RedirectedCenter3D.x)+
            ( RedirectPoints[A].y- m_RedirectedCenter3D.y)*( RedirectPoints[ A].y- m_RedirectedCenter3D.y)+
            ( RedirectPoints[A].z- m_RedirectedCenter3D.z)*( RedirectPoints[ A].z- m_RedirectedCenter3D.z);
    tempB=( RedirectPoints[B].x- m_RedirectedCenter3D.x)*( RedirectPoints[ B].x- m_RedirectedCenter3D.x)+
            ( RedirectPoints[B].y- m_RedirectedCenter3D.y)*( RedirectPoints[ B].y- m_RedirectedCenter3D.y)+
            ( RedirectPoints[B].z- m_RedirectedCenter3D.z)*( RedirectPoints[ B].z- m_RedirectedCenter3D.z);
    tempC=( RedirectPoints[C].x- m_RedirectedCenter3D.x)*( RedirectPoints[ C].x- m_RedirectedCenter3D.x)+
            ( RedirectPoints[C].y- m_RedirectedCenter3D.y)*( RedirectPoints[ C].y- m_RedirectedCenter3D.y)+
            ( RedirectPoints[C].z- m_RedirectedCenter3D.z)*( RedirectPoints[ C].z- m_RedirectedCenter3D.z);
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

    return Radius;
}
