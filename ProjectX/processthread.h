#ifndef PROCESSTHREAD_H
#define PROCESSTHREAD_H

#include <QThread>
#include "triangles/algs.h"

class ProcessThread : public QThread
{
    Q_OBJECT
public:
    explicit ProcessThread(QObject *parent = 0);
    virtual void run();
    void setAlgs(Algoritm* alg);
private:
    Algoritm *alg;
signals:
    
public slots:
    
};

#endif // PROCESSTHREAD_H
