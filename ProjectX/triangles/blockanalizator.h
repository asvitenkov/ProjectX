#ifndef BLOCKANALIZATOR_H
#define BLOCKANALIZATOR_H

#include <QThread>
#include <QVector>
#include <triangles/structs.h>

class BlockAnalizator : public QThread
{
    Q_OBJECT
public:
    explicit BlockAnalizator(QVector<triangle_t> &tr, QVector<TPoint2> *pt, QObject *parent = 0);
    void run();

private:
    QVector<triangle_t> triangles;
    QVector<TPoint2>* points;

    bool sameTriangle(triangle_t t1, triangle_t t2);
    bool crossedTriangles(triangle_t t1, triangle_t t2);
    bool crossedLines(TPoint2 p11, TPoint2 p12, TPoint2 p21, TPoint2 p22);
    bool crossedRadius(triangle_t t1, triangle_t t2);
signals:
    
public slots:
    
};

#endif // BLOCKANALIZATOR_H
