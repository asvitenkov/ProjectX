#ifndef EDGESELECTOR_H
#define EDGESELECTOR_H

#include <QVector>
#include "triangles/structs.h"
#include "triangles/exception.h"

#define ACCURACY 0.000001
#define DABS(x) (x>0)?(x):(-x)

/* Written by Denis Vashchuk */
/* Class for searching border points */

class EdgeSelector
{
public:
    EdgeSelector();

    void readFile(const char* path);

    /* method written by Denis Novogrodski and Denis Vashchuk*/
    QVector<unsigned int> findEdgePoints(const QVector<TPoint3> &srcPoints) throw (NoLinesException);
    /* Вектор линий */
    QVector<line> lines;
    /* Вектор точек*/
    QVector<TPoint3> points;
};

#endif // EDGESELECTOR_H
