#ifndef TRIANGLESHARED_H
#define TRIANGLESHARED_H

#include "Math/MathDefines.h"

#include <QVector>

class TriangleShared
{
public:
    TriangleShared();
    TriangleShared(const TriangleShared&  );

    static QVector<TPoint3> Points;
    static QVector<TPoint3> RedirectPoints;

    void calcCenter();
    bool crossedTriangles(const TriangleShared& t2);
    TVector countNormal() const;
    void findDistance(const TVector& v);
    bool crossedRadius(const TriangleShared& t2) const;
    bool sameTriangle(const TriangleShared& t2);
    double findRadius();

    int A, B, C;
    double distance, radius;
    TPoint3 center_3D, redirectedCenter_3D;
    bool dead;
private:
    bool crossedLines(TPoint2 p11, TPoint2 p12, TPoint2 p21, TPoint2 p22);
};

#endif // TRIANGLESHARED_H
