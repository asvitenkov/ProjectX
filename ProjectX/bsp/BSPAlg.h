#ifndef BSPALG_H
#define BSPALG_H

#include <triangles/structs.h>
#include <QVector>
#include <QObject>

class BSPAlg : public QObject
{
    Q_OBJECT
public:
    explicit BSPAlg(QVector<triangle_t>* triangles, QVector<TPoint3>* points, QObject *parent = 0);

private:
    QVector<triangle_t> *triangles;
    QVector<TPoint3> *points;

    void split();
    void makeTree();
signals:
    
public slots:
    
};

#endif // BSPALG_H
