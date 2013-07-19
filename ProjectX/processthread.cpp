#include "processthread.h"

ProcessThread::ProcessThread(QObject *parent) :
    QThread(parent)
{
}

void ProcessThread::run()
{
    alg->manipulate();
}

void ProcessThread::setAlgs(Algoritm *alg)
{
    this->alg = alg;
}
