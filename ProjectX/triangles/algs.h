#ifndef ALGS_H
#define ALGS_H

#include "exception.h"

#include "Math/TriangleShared.h"
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
    QVector<TriangleShared> manipulate();

    void setPolVector(TVector& v);
    QVector<TPoint3> findEdgePoints() throw (NoLinesException);

    /* Work with files */
    QVector<TriangleShared> readMailFile(const char *path);
    void writeFile(const char* path, QVector<TriangleShared>& triangles);

    QVector<TriangleShared> srcTriangles;
    QVector<TriangleShared> outTriangles;

private:
    TVector polVec;
    std::string status;
    unsigned int checkedTrianglesCount;

    /* main action groups */
    void cleanPhase2();

    TPoint3 countVectorProduct(TPoint3 a, TPoint3 b);
    double countScalarProduct(TPoint3 a, TPoint3 b);
    double countLength(TPoint3 vec);
    double countAngle(TPoint3 a, TPoint3 b);
    /* */
    void cleanByNormal(const QVector<TriangleShared>& triangles);
    void translate();

    /* For multithread */
    void killNonShown(int b, int e);


    /* TODO delete this methods */
    TPoint3 redirectPoint(TPoint3 p, TVector v);

    bool redirectPoints();
    /* For multithread */
    bool redirectPoints(int b, int e);

    void findRadiuses();
    void prepareTriangles();
    void redirectRadiuses();
signals:
    void statusChanged(QString);
    void processStateChanged(int);
};
#endif // ALGS_H
