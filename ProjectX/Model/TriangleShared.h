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

    inline const int& A() const {return m_A;}
    inline const int& B() const {return m_B;}
    inline const int& C() const {return m_C;}

    void SetVisible(bool visible) {m_Visible = visible;}
    bool IsVisible() const {return m_Visible;}

    double Square() const
    {
        return m_square;
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
