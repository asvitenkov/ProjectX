#ifndef ALGS_H
#define ALGS_H

#include "structs.h"
#include "exception.h"

#include <QVector>

#include <QString>
#include <QDateTime>
#include <QObject>

#define START_PROCCESS 15
#define THREAD_COUNT 16
#define NET_SIZE 10


class Algoritm: public QObject{
   /* QT lib use */
    Q_OBJECT
    /* Threads which use private alg functions */
    friend class ReditrectThread;
    friend class CleanThread;

public:
    QVector<triangle_t> manipulate();

    void setPolVector(polarVector& v);
    QVector<TPoint3> findEdgePoints() throw (NoLinesException);

    /* Work with files */
    QVector<triangle_t> readMailFile(const char *path);
    void writeFile(const char* path, QVector<triangle_t>& triangles);

    QVector<triangle_t> srcTriangles;
    QVector<TPoint3> points_3D;
    QVector<TPoint3> redirectedPoints_3D;
    QVector<TPoint2> points_2D;
    QVector<triangle_t> outTriangles;

private:
    polarVector polVec;
    std::string status;
    unsigned int checkedTrianglesCount;
    double sx, sy;

    /* main action groups */
    void cleanPhase2();

    TPoint3 countNormal(triangle_t tr);
    TPoint3 countVectorProduct(TPoint3 a, TPoint3 b);
    double countScalarProduct(TPoint3 a, TPoint3 b);
    double countLength(TPoint3 vec);
    double countAngle(TPoint3 a, TPoint3 b);
    /* */
    void cleanByNormal(const QVector<triangle_t>& triangles);
    triangle_t findMiddle(triangle_t &t);
    triangle_t findDistance(triangle_t &t, polarVector v);
    triangle_t findLocals(triangle_t &t);
    TPoint3 translateLocal(TPoint3 p, TPoint3 cVector,TPoint3 nVector,TPoint3 bVector);
    void translate();
    double findRadius(triangle_t tr);

    /* For multithread */
    void killNonShown(int b, int e);


    /* TODO delete this methods */
    bool crossedTriangles(triangle_t& t1, triangle_t& t2);
    bool crossedLines(TPoint2 p11, TPoint2 p12, TPoint2 p21, TPoint2 p22);
    bool crossedRadius(triangle_t& t1, triangle_t& t2);
    bool sameTriangle(triangle_t& t1, triangle_t& t2);
    TPoint3 redirectPoint(TPoint3 p, polarVector v);

    bool redirectPoints();
    /* For multithread */
    bool redirectPoints(int b, int e);

    void findRadiuses();
    void prepareTriangles();
    float raschet(triangle_t t,float R, int tol,float e1,float m1,float k0,float P0,float P1);
    void redirectRadiuses();

    double findMaxRadius();
    void calculateWindowSize();

    TComplex reflectionCoefficient(triangle_t t,float R, int tol,float e1,float m1,float k0,float P0,float P1);
    TPoint3 countOrt(TPoint3 p);
    triangle_t findVectors(triangle_t &t);
    TPoint3 countSum(TPoint3 p1, TPoint3 p2);
    TPoint3 countMultiple(TPoint3 p, float m);
    TPoint3 countEVector(TPoint3 p);
    TPoint3 countProectionVector(TPoint3 p);
    TPoint3 countProectionComplexVector(TPoint3 p);
    TPoint3 countDecartVector(polarVector v);

signals:
    void statusChanged(QString);
    void processStateChanged(int);
};

#include <QThread>
class ReditrectThread : public QThread
{
    Q_OBJECT
public:
    explicit ReditrectThread(QObject *parent = 0);
    void run();
signals:

public slots:

};
class CleanThread: public QThread
{
     Q_OBJECT
public:
    CleanThread();
    void run();
    void setJob(Algoritm* alg, int b, int e);
    static void sleep(unsigned long time)
    {
        QThread::sleep(time);
    }

private:
    int b, e;
    Algoritm *alg;
};
#endif // ALGS_H
