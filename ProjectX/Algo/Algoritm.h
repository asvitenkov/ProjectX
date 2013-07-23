#ifndef ALGS_H
#define ALGS_H

#include "exception.h"

#include "Model/TriangleShared.h"
#include "Model/Model.h"

#include <QVector>

#include <QString>
#include <QDateTime>
#include <QObject>

#define START_PROCCESS 15
#define THREAD_COUNT 16
#define NET_SIZE 10

struct TriangleDescr
{
    bool IsDead;

    double Distance;
    double Radius;

    TPoint3 RedirectedCenter3D;
    TPoint3 Center3D;
    TVector Normal;
};

class Algoritm: public QObject{
    Q_OBJECT
public:
    void SetModel(Model* model);
    void SetSightVector(const TVector& v);

    Model* manipulate();

    QVector<TPoint3> findEdgePoints() throw (NoLinesException);
private:
    Model* m_model;

    QVector<TPoint3>        m_redirectPoints;
    QVector<TriangleDescr>  m_descr;
    QVector<TPoint3>        m_point2D;

    TVector m_vec;
    std::string status;
    unsigned int checkedTrianglesCount;
    int m_amout;

    void Prepare();
    void Calc(int start, int end);

    /* main action groups */
    void cleanPhase1();
    void cleanPhase2();
    void cleanByNormal(int start, int end);
    void killNonShown(int start, int end);

    void FindRadius(const TriangleShared& tr, TriangleDescr& descr);
    void FindDistance(const TriangleShared& tr, TriangleDescr& descr) const;
    void FindNormal(const TriangleShared& tr, TriangleDescr& descr) const;
    void FindCenter(const TriangleShared& tr, TriangleDescr& descr);
    void FindRedirectPoint(int start, int end);

    bool HaveSamePoint(const TriangleShared& t1, const TriangleShared& t2) const;
    bool IsCrossedTriangles(const TriangleShared& t1, const TriangleShared& t2) const;
    bool IsCrossedLines(const TPoint2& p11, const TPoint2& p12, const TPoint2& p21, const TPoint2& p22) const;
    bool IsCrossedRadius(const TriangleDescr& t1, const TriangleDescr& t2) const;

    void prepareTriangles();
    void redirectRadiuses();
signals:
    void statusChanged(QString);
    void processStateChanged(int);
};

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
#endif // ALGS_H
