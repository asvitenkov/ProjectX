#include "TriangleShared.h"
#include "Model.h"

QVector<TPoint2> m_point2D;
QVector<TPoint3> TriangleShared::Points3D;

TriangleShared::TriangleShared(Model* parent)
    : m_parent(parent), m_Visible(true), m_square(-1)
{

}

TriangleShared::TriangleShared(const TriangleShared& tr)
    : m_parent(tr.m_parent)
    , m_A(tr.m_A), m_B(tr.m_B), m_C(tr.m_C)
    ,  m_Visible(tr.m_Visible)
{

}

void TriangleShared::Set(int _A, int _B, int _C)
{
    m_A = _A;
    m_B = _B;
    m_C = _C;

    const TVector v1(p1().x - p2().x, p1().y - p2().y, p1().z - p2().z);
    const TVector v2(p2().x - p3().x, p2().y - p3().y, p2().z - p3().z);

    const TVector v3 = v1.CrossProduct(v2);

    m_square = abs(v3.Length()/2);
}

const TPoint3& TriangleShared::p1() const
{
    return m_parent->m_points[m_A];
}

const TPoint3& TriangleShared::p3() const
{
    return m_parent->m_points[m_C];
}

const TPoint3& TriangleShared::p2() const
{
    return m_parent->m_points[m_B];
}
