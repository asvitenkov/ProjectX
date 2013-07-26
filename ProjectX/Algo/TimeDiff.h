#ifndef TIMEDIFF_H
#define TIMEDIFF_H

#include <QTime>

class TimeDiff
{
public:
    TimeDiff();
    ~TimeDiff();

    void    Start();
    unsigned int Diff();
    void    End();
private:
    QTime m_start;
};

#endif // TIMEDIFF_H
