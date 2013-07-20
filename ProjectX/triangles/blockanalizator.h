#ifndef BLOCKANALIZATOR_H
#define BLOCKANALIZATOR_H

#include <QThread>
#include <QVector>
#include "Math/TriangleShared.h"

class BlockAnalizator : public QThread
{
    Q_OBJECT
public:
    explicit BlockAnalizator(QVector<TriangleShared> & , QVector<TPoint2> *pt, QObject *parent = 0);
    void run();

private:
    QVector<TriangleShared>  triangles;
    QVector<TPoint2>* points;

    bool sameTriangle(TriangleShared t1, TriangleShared t2);
    bool crossedTriangles(TriangleShared t1, TriangleShared t2);
    bool crossedLines(TPoint2 p11, TPoint2 p12, TPoint2 p21, TPoint2 p22);
    bool crossedRadius(TriangleShared t1, TriangleShared t2);
signals:
    
public slots:
    
};

#endif // BLOCKANALIZATOR_H
