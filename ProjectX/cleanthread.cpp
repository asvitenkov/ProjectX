#include "triangles/algs.h"

CleanThread::CleanThread()
{
}

void CleanThread::run()
{
    alg->killNonShown(b,e);
}

void CleanThread::setJob(Algoritm *alg, int b, int e)
{
    this->alg = alg;
    this->b = b;
    this->e = e;
}
