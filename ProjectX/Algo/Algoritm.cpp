#define _USE_MATH_DEFINES

#include "Algo/Algoritm.h"

#include <iostream>
#include <fstream>
#include <ostream>
#include <cstdlib>
#include <QVector>
#include <cmath>
#include <ctime>

#include <QDebug>

Model* Algoritm::manipulate()
{
    m_redirectPoints.clear();
    m_point2D.clear();
    checkedTrianglesCount = 0x0;
    time_t start_timestamp = time(NULL);

    //process
    {
        emit processStateChanged(0);

        emit statusChanged("Preparing triangles:");

        Prepare();
        Calc(0, m_model->GetTriangles().size());

        emit statusChanged("done\n");

        cleanPhase2();
    }

    for(int i=0; i < m_model->GetTriangles().size(); ++i)
    {
        m_model->GetTriangles()[i].SetVisible(!m_descr[i].IsDead);
    }

    time_t end_timestamp = time(NULL);
    std::cout << start_timestamp << " " << end_timestamp << std::endl;

    double workTime = difftime(end_timestamp, start_timestamp);

    emit statusChanged("Work time: " +QString::number(workTime) +" s or " +QString::number(workTime/60) +" min");
    emit statusChanged("Count of projection triangles: " +QString::number(m_model->GetTriangles().size()));

    emit processStateChanged(100);

    return m_model;
}

void Algoritm::SetModel(Model* model)
{
    m_model = model;
}

void Algoritm::SetSightVector(const TVector &v)
{
    m_vec = v;
}

void Algoritm::cleanPhase2()
{
    //clean - phase 2
    emit statusChanged("Clean phase 2: ");

    killNonShown(0, m_model->GetTriangles().size());

    emit processStateChanged(98);
    emit statusChanged("done\n");
}

void Algoritm::killNonShown(int start, int end)
{
    const int n = m_model->GetTriangles().size();

    if(end > n)
    {
        end = n;
    }

    for(int i = start; i < end; ++i)
    {
        if(m_descr[i].IsDead == false)
        {
            for(int j = i+1; j < m_model->GetTriangles().size(); ++j)
            {
                if(!HaveSamePoint(m_model->GetTriangles()[j], m_model->GetTriangles()[i]))
                {
                    if(IsCrossedTriangles(m_model->GetTriangles()[j], m_model->GetTriangles()[i]))
                    {
                        if(m_descr[i].Distance < m_descr[j].Distance)
                            m_descr[i].IsDead = true;
                        else
                            m_descr[j].IsDead= true;
                    }
                }
            }
        }

        checkedTrianglesCount++;
    }
}

void Algoritm::Prepare()
{
    const size_t N = m_model->GetTriangles().size();

    m_descr.resize(N);
    m_redirectPoints.resize(N);
    m_point2D.resize(N);
}

void Algoritm::Calc(int start, int end)
{
    const int N = m_model->GetTriangles().size();
    if(end > N)
    {
        end = N;
    }

    FindRedirectPoint(start, end);

    for(int i = start; i < N; i++)
    {
        m_descr[i].IsDead = false;
        FindCenter(m_model->GetTriangles()[i], m_descr[i]);
        FindRadius(m_model->GetTriangles()[i], m_descr[i]);
        FindDistance(m_model->GetTriangles()[i], m_descr[i]);
        FindNormal(m_model->GetTriangles()[i], m_descr[i]);
    }
}

void Algoritm::FindRadius(const TriangleShared& tr, TriangleDescr& descr)
{
    double Radius = 0;
    double tempA=0;
    double tempB=0;
    double tempC=0;

    TPoint3& p1 = m_redirectPoints[tr.A()];
    TPoint3& p2 = m_redirectPoints[tr.B()];
    TPoint3& p3 = m_redirectPoints[tr.C()];

    TPoint3& cntr = descr.RedirectedCenter3D;

    tempA=( p1.x - cntr.x)*( p1.x- cntr.x)+
            ( p1.y- cntr.y)*( p1.y- cntr.y)+
            ( p1.z- cntr.z)*( p1.z- cntr.z);

    tempB=( p2.x- cntr.x)*( p2.x- cntr.x)+
            ( p2.y- cntr.y)*( p2.y- cntr.y)+
            ( p2.z- cntr.z)*( p2.z- cntr.z);

    tempC=(p3.x- cntr.x)*( p3.x- cntr.x)+
            (p3.y- cntr.y)*( p3.y- cntr.y)+
            (p3.z- cntr.z)*( p3.z- cntr.z);

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

    descr.Radius = Radius;
}


void Algoritm::FindCenter(const TriangleShared& tr, TriangleDescr& descr)
{
    descr.Center3D = (tr.p1() + tr.p2() + tr.p3())/3;
    descr.RedirectedCenter3D = (MathHelper::Redirect(&descr.Center3D, &m_vec));
}

void Algoritm::FindRedirectPoint(int start, int end)
{
    const int N = m_model->GetPoints().size();

    if(end > N)
    {
        end = N;
    }

    for(int i = start; i < end; ++i)
    {
        m_redirectPoints[i] = (MathHelper::Redirect(&m_model->GetPoints()[i], &m_vec));
    }

    const double tcos = cos(m_vec.Theta());
    const double tsin = sin(m_vec.Theta());
    const double pcos = cos(m_vec.Phi());
    const double psin = sin(m_vec.Phi());

    for(int i = start; i < end; ++i)
    {
        m_point2D[i].x = -psin*m_redirectPoints[i].x + pcos*m_redirectPoints[i].y;
        m_point2D[i].y = -tcos*pcos*m_redirectPoints[i].x + -tcos*psin*m_redirectPoints[i].y + tsin*m_redirectPoints[i].z;
    }
}

bool Algoritm::IsCrossedTriangles(const TriangleShared& t1, const TriangleShared& t2) const
{
    const TPoint3& p11 = m_point2D[t1.A()];
    const TPoint3& p12 = m_point2D[t1.B()];
    const TPoint3& p13 = m_point2D[t1.C()];

    const TPoint3& p21 = m_point2D[t2.A()];
    const TPoint3& p22 = m_point2D[t2.B()];
    const TPoint3& p23 = m_point2D[t2.C()];

    if(IsCrossedLines(p11, p12, p21, p22)) return true;
    if(IsCrossedLines(p11, p13, p21, p22)) return true;
    if(IsCrossedLines(p12, p13, p21, p22)) return true;

    if(IsCrossedLines(p11, p12, p21, p23)) return true;
    if(IsCrossedLines(p11, p13, p21, p23)) return true;
    if(IsCrossedLines(p12, p13, p21, p23)) return true;

    if(IsCrossedLines(p11, p12, p23, p22)) return true;
    if(IsCrossedLines(p11, p13, p23, p22)) return true;
    if(IsCrossedLines(p12, p13, p23, p22)) return true;

    return false;
}

bool Algoritm::IsCrossedLines(const TPoint2& p11, const TPoint2& p12, const TPoint2& p21, const TPoint2& p22) const
{
    double k1 = (p11.y - p12.y)/(p11.x - p12.x);
    double k2 = (p21.y - p22.y)/(p21.x - p22.x);

    double x = (p21.y - p11.y + k1*p11.x - k2*p21.x)/(k1 - k2);

    if(((x <= qMax(p11.x,p12.x))&&(x >= qMin(p12.x,p11.x)))&&((x <= qMax(p22.x,p21.x))&&(x >= qMin(p21.x,p22.x))))
        return true;

    return false;
}

void Algoritm::FindDistance(const TriangleShared& tr, TriangleDescr& descr) const
{
    descr.Distance = (m_vec.x * descr.Center3D.x
                + m_vec.y * descr.Center3D.y
                + m_vec.z * descr.Center3D.z);
}

void Algoritm::FindNormal(const TriangleShared& tr, TriangleDescr& descr) const
{
    TVector a, b;

    a = m_point2D[tr.B()] - m_point2D[tr.A()];
    b = m_point2D[tr.C()] - m_point2D[tr.A()];

    descr.Normal = a.CrossProduct(b);
}

bool Algoritm::IsCrossedRadius(const TriangleDescr& t1, const TriangleDescr& t2) const
{
    double d;
    const TPoint3& p1 = t1.RedirectedCenter3D;
    const TPoint3& p2 = t2.RedirectedCenter3D;

    // TODO: Can refactor
    d = (p1.x - p2.x)*(p1.x - p2.x)+
            (p1.y - p2.y)*(p1.y - p2.y)+
            (p1.z - p2.z)*(p1.z - p2.z);

    if(d < ( d + t2.Radius ))
        return true;

    return false;
}

bool Algoritm::HaveSamePoint(const TriangleShared& t1, const TriangleShared& t2) const
{
    if(t1.A() == t2.A()) return true;
    if(t1.A() == t2.B()) return true;
    if(t1.A() == t2.C()) return true;
    if(t1.B() == t2.A()) return true;
    if(t1.B() == t2.B()) return true;
    if(t1.B() == t2.C()) return true;
    if(t1.C() == t2.A()) return true;
    if(t1.C() == t2.B()) return true;
    if(t1.C() == t2.C()) return true;
    return false;
}
