#ifndef TRIANGLESHARED_H
#define TRIANGLESHARED_H

#include "Math/MathDefines.h"

#include <QVector>

class TriangleShared
{
public:
    TriangleShared();
    TriangleShared(const TriangleShared&  );

    static QVector<TPoint2> Points2D;
    static QVector<TPoint3> Points3D;
    static QVector<TPoint3> RedirectPoints;

    void Set(int A, int B, int C);

    void CaclFileds(const TVector& v);

    void calcCenter();
    bool crossedTriangles(const TriangleShared& t2);
    void findVectors();
    TVector countNormal() const;
    void findDistance(const TVector& v);
    bool crossedRadius(const TriangleShared& t2) const;
    bool sameTriangle(const TriangleShared& t2);
    double findRadius();

    const TPoint3& p1() const { return Points3D[A]; }
    const TPoint3& p2() const { return Points3D[B]; }
    const TPoint3& p3() const { return Points3D[C]; }

    const TPoint3& Center() const { return m_Center3D; }
    const TPoint3& RedirectCenter() const { return m_RedirectedCenter3D; }

    double distance, radius;
    bool dead;
private:
    bool crossedLines(TPoint2 p11, TPoint2 p12, TPoint2 p21, TPoint2 p22);
    int A, B, C;
    TPoint3 m_Center3D;
    TPoint3 m_RedirectedCenter3D;
};

#endif // TRIANGLESHARED_H
