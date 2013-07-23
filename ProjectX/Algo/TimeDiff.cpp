#include "TimeDiff.h"

TimeDiff::TimeDiff()
{

}

TimeDiff::~TimeDiff()
{

}

void    TimeDiff::Start()
{
    m_start = QTime::currentTime();
}

int   TimeDiff::Diff()
{
    return QTime::currentTime().msec() - m_start.msec();
}

void    TimeDiff::End()
{

}
