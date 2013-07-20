#ifndef ALGS_H
#define ALGS_H

#include "Math/TriangleShared.h"
#include "structs.h"
#include "exception.h"

#include <QVector>

#include <QThread>
#include <QString>
#include <QDateTime>
#include <QObject>

#define START_PROCCESS 15
#define THREAD_COUNT 1
#define NET_SIZE 10


class Algoritm: public QObject{
   /* QT lib use */
    Q_OBJECT
    /* Threads which use private alg functions */
    friend class ReditrectThread;
    friend class CleanThread;

public:
    QVector<TriangleShared> manipulate();

    void SetViewVector(const TVector& v);
    QVector<TPoint3> findEdgePoints() throw (NoLinesException);

    /* Work with files */
    QVector<TriangleShared> readMailFile(const char *path);
    void writeFile(const char* path, QVector<TriangleShared>&  iangles);

    QVector<TriangleShared> srcTriangles;
    QVector<TPoint2> points_2D;
    QVector<TriangleShared> outTriangles;

private:
    TVector m_viewVector;
    std::string status;
    unsigned int checkedTrianglesCount;
    double sx, sy;

    /* main action groups */
    void cleanPhase2();

    /* */
    void cleanByNormal(const QVector<TriangleShared>&  iangles);
    TriangleShared findLocals(TriangleShared &t);
    TPoint3 translateLocal(TPoint3 p, TPoint3 cVector,TPoint3 nVector,TPoint3 bVector);
    void  translate();
    double findRadius(TriangleShared  );

    /* For multithread */
    void killNonShown(int b, int e);

    /* TODO delete this methods */
    TPoint3 redirectPoint(TPoint3 p, TVector v);

    bool redirectPoints();

    void findRadiuses();
    void prepareTriangles();
    void redirectRadiuses();

    double findMaxRadius();
    void calculateWindowSize();

    TPoint3 countEVector(TPoint3 p);
    TPoint3 countProectionVector(TPoint3 p);
    TPoint3 countProectionComplexVector(TPoint3 p);
    TPoint3 countDecartVector(polarVector v);

signals:
    void statusChanged(QString);
    void processStateChanged(int);
};

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
