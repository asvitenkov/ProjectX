#include "blockanalizator.h"

BlockAnalizator::BlockAnalizator(QVector<triangle_t> &tr, QVector<TPoint2> *pt, QObject *parent) :
    QThread(parent), triangles(tr), points(pt)
{
}

bool BlockAnalizator::crossedTriangles(triangle_t t1, triangle_t t2)//просматривает каждую сторону с каждой см. ниже
{
    if(crossedLines(points->at(t1.A), points->at(t1.B), points->at(t2.A), points->at(t2.B))) return true;
    if(crossedLines(points->at(t1.A), points->at(t1.C), points->at(t2.A), points->at(t2.B))) return true;
    if(crossedLines(points->at(t1.B), points->at(t1.C), points->at(t2.A), points->at(t2.B))) return true;

    if(crossedLines(points->at(t1.A), points->at(t1.B), points->at(t2.A), points->at(t2.C))) return true;
    if(crossedLines(points->at(t1.A), points->at(t1.C), points->at(t2.A), points->at(t2.C))) return true;
    if(crossedLines(points->at(t1.B), points->at(t1.C), points->at(t2.A), points->at(t2.C))) return true;

    if(crossedLines(points->at(t1.A), points->at(t1.B), points->at(t2.C), points->at(t2.B))) return true;
    if(crossedLines(points->at(t1.A), points->at(t1.C), points->at(t2.C), points->at(t2.B))) return true;
    if(crossedLines(points->at(t1.B), points->at(t1.C), points->at(t2.C), points->at(t2.B))) return true;
    return false;
}

bool BlockAnalizator::crossedLines(TPoint2 p11, TPoint2 p12, TPoint2 p21, TPoint2 p22)////пересекаются ли линии треугольников?
{
    double k1 = (p11.y - p12.y)/(p11.x - p12.x);
    double k2 = (p21.y - p22.y)/(p21.x - p22.x);
    double x = (p21.y - p11.y + k1*p11.x - k2*p21.x)/(k1 - k2);
    if(((x <= qMax(p11.x,p12.x))&&(x >= qMin(p12.x,p11.x)))&&((x <= qMax(p22.x,p21.x))&&(x >= qMin(p21.x,p22.x))))
        return true;
    return false;
}

bool BlockAnalizator::crossedRadius(triangle_t t1, triangle_t t2)//пересекаются ли радиусы треугольников?
{
    double d;
    d = (t1.redirectedCenter_3D.x - t2.redirectedCenter_3D.x)*(t1.redirectedCenter_3D.x - t2.redirectedCenter_3D.x)+
            (t1.redirectedCenter_3D.y - t2.redirectedCenter_3D.y)*(t1.redirectedCenter_3D.y - t2.redirectedCenter_3D.y)+
            (t1.redirectedCenter_3D.z - t2.redirectedCenter_3D.z)*(t1.redirectedCenter_3D.z - t2.redirectedCenter_3D.z);
    if(d < ( t1.radius + t2.radius ))
        return true;
    return false;
}

bool BlockAnalizator::sameTriangle(triangle_t t1, triangle_t t2)//является ли соседним треугольник?
{
    if(t1.A==t2.A) return true;
    if(t1.A==t2.B) return true;
    if(t1.A==t2.C) return true;
    if(t1.B==t2.A) return true;
    if(t1.B==t2.B) return true;
    if(t1.B==t2.C) return true;
    if(t1.C==t2.A) return true;
    if(t1.C==t2.B) return true;
    if(t1.C==t2.C) return true;
    return false;
}

void BlockAnalizator::run()
{
    int a = 0;
    int n = triangles.size() - 1;

    for(int i = 0; i < n; ++i)
    {

        //emit processStateChanged((100*i)/(n-1));

        if(!triangles[i].dead)
            for(int j = i+1; j < triangles.size(); ++j)
            {
                if((!triangles[j].dead) && (!sameTriangle(triangles[i],triangles[j])))
                {
                    if(crossedRadius(triangles[i],triangles[j]))
                        if(crossedTriangles(triangles[i],triangles[j]))
                        {
                            if(triangles[i].distance < triangles[j].distance)
                                triangles[i].dead = true;
                            else
                                triangles[j].dead = true;
                        }
                }
            }
    }
}
