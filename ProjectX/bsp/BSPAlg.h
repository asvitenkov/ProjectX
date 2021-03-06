#ifndef BSPALG_H
#define BSPALG_H

#include "Model/TriangleShared.h"

#include <QVector>
#include <QObject>

class BSPAlg : public QObject
{
    Q_OBJECT
public:
    explicit BSPAlg(QVector<TriangleShared>* triangles, QVector<TPoint3>* points, QObject *parent = 0);

private:
    QVector<TriangleShared> *triangles;
    QVector<TPoint3> *points;

    void split();
    void makeTree();
signals:
    
public slots:
    
};

#endif // BSPALG_H
