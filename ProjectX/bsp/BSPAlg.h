#ifndef BSPALG_H
#define BSPALG_H

#include "Math/TriangleShared.h"

#include <QVector>
#include <QObject>

class BSPAlg : public QObject
{
    Q_OBJECT
public:
    explicit BSPAlg(QVector<TriangleShared>*  iangles, QVector<TPoint3>* points, QObject *parent = 0);

private:
    QVector<TriangleShared> * iangles;
    QVector<TPoint3> *points;

    void split();
    void makeTree();
signals:
    
public slots:
    
};

#endif // BSPALG_H
