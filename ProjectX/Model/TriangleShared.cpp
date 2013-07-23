#include "TriangleShared.h"
#include "Model.h"

QVector<TPoint2> m_point2D;
QVector<TPoint3> TriangleShared::Points3D;

TriangleShared::TriangleShared(Model* parent)
    : m_parent(parent), m_Visible(true)
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

const int& TriangleShared::A() const
{
    return m_A;
}

const int& TriangleShared::B() const
{
    return m_B;
}

const int& TriangleShared::C() const
{
    return m_C;
}
