#ifndef EDGESELECTOR_H
#define EDGESELECTOR_H

#include "Math/MathDefines.h"

#include <QVector>
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

    QVector<unsigned int> findEdgePoints(const QVector<TPoint3> &srcPoints) throw (NoLinesException);
    QVector<LineShared> lines;
    QVector<TPoint3> points;
};

#endif // EDGESELECTOR_H
