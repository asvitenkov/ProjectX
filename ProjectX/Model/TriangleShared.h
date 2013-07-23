#ifndef TRIANGLESHARED_H
#define TRIANGLESHARED_H

#include "Math/MathDefines.h"
#include <QVector>

class Model;

class TriangleShared
{
public:
    TriangleShared(Model* parent = 0);
    TriangleShared(const TriangleShared&);

    static QVector<TPoint2> Points2D;
    static QVector<TPoint3> Points3D;

    void Set(int A, int B, int C);

    const TPoint3& p1() const;
    const TPoint3& p3() const;
    const TPoint3& p2() const;

    const int& A() const;
    const int& B() const;
    const int& C() const;

    void SetVisible(bool visible) {m_Visible = visible;}
    bool IsVisible() const {return m_Visible;}


    TVector Normal() const
    {
        double xx = p1().y * p2().z - p1().z * p2().y;
        double yy = p1().z * p2().x - p1().x * p2().z;
        double zz = p1().x * p2().y - p1().y * p2().x;

        return TVector(xx, yy, zz);
    }

    double Square() const
    {
        const TVector v1(p1().x - p2().x, p1().y - p2().y, p1().z - p2().z);
        const TVector v2(p2().x - p3().x, p2().y - p3().y, p2().z - p3().z);

        const TVector v3 = v1.CrossProduct(v2);

        return v3.Length()/2;
    };

    TPoint3 Center() const
    {
        return (p1()+p2()+p3())/3;
    }

private:
    Model* m_parent;

    int m_A, m_B, m_C;

    bool m_Visible;

    TPoint3	m_center;
    TVector	m_normal;
    double	m_square;
};

#endif // TRIANGLESHARED_H
