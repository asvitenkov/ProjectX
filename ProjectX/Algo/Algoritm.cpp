#define _USE_MATH_DEFINES

#include "Algo/Algoritm.h"

#include <iostream>
#include <fstream>
#include <ostream>
#include <cstdlib>
#include <QVector>
#include <QTimer>
#include <cmath>

#include "TimeDiff.h"

#include <QDebug>

Algoritm::Algoritm(QObject* obj)
    : QObject(obj)
    , m_Epsilon(0.01)
{
}

Algoritm::~Algoritm()
{
}

Model* Algoritm::Calculate()
{
    TimeDiff time;
    time.Start();

    m_finalProcessedTriangles = m_model->GetTriangles().size()*2;
    emit processStateChanged(0);
    //process
    {
        emit statusChanged("Sight vector: Theta-" + QString::number(m_vec.Theta())+", Phi-" + QString::number(m_vec.Phi()));
        emit statusChanged("Preparing triangles..." + QString::number(m_model->GetTriangles().size()));

        Prepare();
        Calc(0, m_model->GetTriangles().size());

        emit statusChanged("Done... Time: " + QString::number(time.Diff()));

        emit statusChanged("Clean. Phase 1...");

        CleanPhase1();

        emit statusChanged("Done... Time: " + QString::number(time.Diff()));

        //CleanPhase2();

        emit statusChanged("Done... Time: " + QString::number(time.Diff()));
    }
    emit processStateChanged(98);

    for(int i=0; i < m_model->GetTriangles().size(); ++i)
    {
        m_model->GetTriangles()[i].SetVisible(!m_descr[i].IsDead);
    }

    emit statusChanged("Time: " + QString::number(time.Diff()));
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

void Algoritm::CleanPhase1()
{
    m_processedTriangles  = 0;
    QThreadPool* pool = new QThreadPool();

    const int step = m_model->GetTriangles().size()/THREAD_COUNT;
    for(int i = 0; i<THREAD_COUNT; i++)
    {
        pool->start(new MultiTask(this, &Algoritm::CleanByNormal, i*step, step + i*step));
    }
    pool->waitForDone();
    delete pool;
    //CleanByNormal(0, m_model->GetTriangles().size());
}

void Algoritm::CleanPhase2()
{
    m_processedTriangles  = 0;
    QThreadPool* pool = new QThreadPool();

    const int step = m_model->GetTriangles().size()/THREAD_COUNT;
    for(int i = 0; i<THREAD_COUNT; i++)
    {
        pool->start(new MultiTask(this, &Algoritm::KillNonShown, i*step, step + i*step));
    }
    pool->waitForDone();
    delete pool;
    //KillNonShown(0, m_model->GetTriangles().size());
}

void Algoritm::CleanByNormal(int start, int end)
{
    const int n = m_model->GetTriangles().size();

    if(end > n)
    {
        end = n;
    }

    for(int i = start; i < end; ++i)
    {
        TriangleShared& t = m_model->GetTriangles()[i];
        TVector v1(t.p1()-t.p2());
        TVector v2(t.p3()-t.p2());

        TVector face(v1.CrossProduct(v2));
        double distance = -(face.x*v1.x + face.y*v1.y + face.z*v1.z);
        TPoint3 center =t.Center();

        double m = -sgn(face.x*center.x + face.y*center.y + face.z*center.z + distance);

        face *= m;
        distance *= m;

        const double value = (face.DotProduct(m_vec)+distance);
        if(!(value < -m_Epsilon))
        {
            m_descr[i].IsDead = true;
        }
        m_processedTriangles++;
    }
}

void Algoritm::KillNonShown(int start, int end)
{
    const int N = m_model->GetTriangles().size();

    if(end > N)
    {
        end = N;
    }

    for(int i = start; i < end; ++i)
    {
        TriangleDescr& descr_i = m_descr[i];
        const TriangleShared& tr_i = m_model->GetTriangles()[i];

        if(!descr_i.IsDead)
        {
            for(int j = i+1; j < N; ++j)
            {
                TriangleDescr& descr_j = m_descr[j];
                const TriangleShared& tr_j = m_model->GetTriangles()[j];

                if(!HaveSamePoint(tr_j, tr_i))
                {
                    if(IsCrossedTriangles(tr_j, tr_i))
                    {
                        if(descr_i.Distance < descr_j.Distance)
                        {
                            descr_i.IsDead = true;
                        }
                        else
                        {
                            descr_j.IsDead= true;
                        }
                    }
                }
            }
        }
        m_processedTriangles++;
    }
}

void Algoritm::Prepare()
{
    const size_t N = m_model->GetTriangles().size();

    m_redirectPoints.clear();
    m_point2D.clear();
    m_descr.clear();

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
        const TriangleShared& tr = m_model->GetTriangles()[i];
        TriangleDescr& descr = m_descr[i];

        descr.IsDead = false;
        FindCenter(tr, descr);
        FindRadius(tr, descr);
        FindDistance(tr, descr);
        FindNormal(tr, descr);
    }
}

void Algoritm::FindRadius(const TriangleShared& tr, TriangleDescr& descr)
{
    TPoint2& p1 = m_point2D[tr.A()];
    TPoint2& p2 = m_point2D[tr.B()];
    TPoint2& p3 = m_point2D[tr.C()];

    TPoint3& cntr = descr.Center2D;

    const double tempA = ( p1.x- cntr.x)*( p1.x- cntr.x)+
            ( p1.y- cntr.y)*( p1.y- cntr.y);

    const double tempB = ( p2.x- cntr.x)*( p2.x- cntr.x)+
            ( p2.y- cntr.y)*( p2.y- cntr.y);

    const double tempC = (p3.x- cntr.x)*( p3.x- cntr.x)+
            (p3.y- cntr.y)*( p3.y- cntr.y);

    descr.SquaredRadius = std::max(tempA, std::max(tempB, tempC));
}

void Algoritm::FindCenter(const TriangleShared& tr, TriangleDescr& descr)
{
    descr.Center3D = tr.Center();
    descr.RedirectedCenter3D = (MathHelper::Redirect(&descr.Center3D, &m_vec));

    const double tcos = cos(m_vec.Theta());
    const double tsin = sin(m_vec.Theta());
    const double pcos = cos(m_vec.Phi());
    const double psin = sin(m_vec.Phi());

    descr.Center2D.x = -psin*descr.RedirectedCenter3D.x + pcos*descr.RedirectedCenter3D.y;
    descr.Center2D.y = -tcos*pcos*descr.RedirectedCenter3D.x + -tcos*psin*descr.RedirectedCenter3D.y + tsin*descr.RedirectedCenter3D.z;
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
    if(IsTrianglePointInTriangle(t1,t2) || IsTrianglePointInTriangle(t2,t1)) return true;

    const TPoint2& p11 = m_point2D[t1.A()];
    const TPoint2& p12 = m_point2D[t1.B()];
    const TPoint2& p13 = m_point2D[t1.C()];

    const TPoint2& p21 = m_point2D[t2.A()];
    const TPoint2& p22 = m_point2D[t2.B()];
    const TPoint2& p23 = m_point2D[t2.C()];

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

bool Algoritm::IsCrossedLines(const TPoint2& p1, const TPoint2& p2, const TPoint2& p3, const TPoint2& p4) const
{
    const double a1 = -(p2.y - p1.y);
    const double b1 = +(p2.x - p1.x);
    const double d1 = -(a1*p1.x + b1*p1.y);

    const double a2 = -(p4.y - p3.y);
    const double b2 = +(p4.x - p3.x);
    const double d2 = -(a2*p3.x + b2*p3.y);

    const double seg1_line2_start = a2*p1.x + b2*p1.y + d2;
    const double seg1_line2_end = a2*p2.x + b2*p2.y + d2;

    const double seg2_line1_start = a1*p3.x + b1*p3.y + d1;
    const double seg2_line1_end = a1*p4.x + b1*p4.y + d1;

    return !(seg1_line2_start * seg1_line2_end >= 0 || seg2_line1_start * seg2_line1_end >= 0);
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
    const TPoint2& p1 = t1.Center2D;
    const TPoint2& p2 = t2.Center2D;

    const double distSqCenter =
            ( p1.x - p2.x)*( p1.x - p2.x)+
            ( p1.y - p2.y)*( p1.y - p2.y);

    return (t1.SquaredRadius + t2.SquaredRadius) < distSqCenter;
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

bool Algoritm::IsTrianglePointInTriangle(const TriangleShared& t1, const TriangleShared& t2) const
{
    const TPoint2& p11 = m_point2D[t1.A()];
    const TPoint2& p12 = m_point2D[t1.B()];
    const TPoint2& p13 = m_point2D[t1.C()];

    const TPoint2* p[3] = {&m_point2D[t2.A()],
                           &m_point2D[t2.B()],
                           &m_point2D[t2.C()]};

    const double xA = p12.x - p11.x;
    const double yA = p12.y - p11.y;

    const double xB = p13.x - p11.x;
    const double yB = p13.y - p11.y;

    bool crossed = false;
    double xC, yC, m, l;
    int i = 0;
    while(!(crossed || i>2))
    {
        xC = p[i]->x - p11.x;
        yC = p[i]->y - p11.y;

        m = (xC*yA - xA*yC)/(xB*yA - xA*yB);

        if(m > 0 && m <= 1.0)
        {
            l = (xC - m*xB) / xA;
            crossed = ( (l >= 0) && ((m + l) <= 1));
        }
        ++i;
    }

    return crossed;
}

void Algoritm::NotifyAboutState()
{
    emit processStateChanged(m_processedTriangles/m_finalProcessedTriangles*100);
}
