#include "TimeDiff.h"

TimeDiff::TimeDiff()
{

}

TimeDiff::~TimeDiff()
{

}

void TimeDiff::Start()
{
    m_start = QTime::currentTime();
}

unsigned int TimeDiff::Diff()
{
    return  m_start.msecsTo(QTime::currentTime());
}

void TimeDiff::End()
{

}
