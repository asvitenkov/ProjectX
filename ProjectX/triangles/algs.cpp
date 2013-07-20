#define _USE_MATH_DEFINES
#include "triangles/algs.h"

#include <iostream>
#include <fstream>
#include <ostream>
#include <cstdlib>
#include <QVector>
#include <cmath>
#include <ctime>

#include <QDebug>


QVector<TriangleShared> Algoritm::manipulate()
{
    outTriangles.clear();
    TriangleShared::RedirectPoints.clear();
    points_2D.clear();
    checkedTrianglesCount = 0x0;
    time_t start_timestamp = time(NULL);

    outTriangles = srcTriangles;

    emit processStateChanged(0);
    /*
        emit statusChanged("Clean phase 1. ");
        //clean - phase 1;
        cleanByNormal(srcTriangles);
        emit statusChanged("done\n");
    */

    //find 3D redirect points in  iangles
    emit statusChanged("Preparing  iangles:");
    prepareTriangles();//v
    emit statusChanged("done\n");

    /* redirect all points from TriangleShared::Points to RedirectPoints vectors  AND RADIUSES*/
    redirectPoints();//v
    redirectRadiuses();
    emit statusChanged("Reditect points: RedirectPoints size is " + QString::number(TriangleShared::RedirectPoints.size()) +"\n");

    //find radiuses and distances from real  iangle to its redirect
    emit statusChanged("Finding radiuses: ");
    findRadiuses();
    emit statusChanged("done\n");

    translate();
    emit statusChanged("Translate: Amount equal to " + QString::number(points_2D.size()));
    emit statusChanged("done\n");

    //finding window size
    this->calculateWindowSize();
    qDebug() << "Window size " << sx << " " << sy;

    cleanPhase2();

    //PROFIT
    QVector<TriangleShared> tempTriangles;
    // making new QVector with result  iangles
    for(size_t i=0; i<outTriangles.size(); ++i)
    {
        if(!outTriangles[i].dead)
            tempTriangles.push_back(outTriangles[i]);
    }
    outTriangles = tempTriangles;

    time_t end_timestamp = time(NULL);
    std::cout << start_timestamp << " " << end_timestamp << std::endl;

    double workTime = difftime(end_timestamp, start_timestamp);

    emit statusChanged("Work time: " +QString::number(workTime) +" s or " +QString::number(workTime/60) +" min");
    emit statusChanged("Count of projection  iangles: " +QString::number(outTriangles.size()));
    emit processStateChanged(100);

    return tempTriangles;
}

void Algoritm::SetViewVector(const TVector& v)
{
    m_viewVector = v;
}

void Algoritm::cleanByNormal(const QVector<TriangleShared>&  iangles)
{
    TPoint3 ptVec(m_viewVector);

    //cleaning backout  iangles by normal
    for(size_t i = 0; i <  iangles.size(); i++)
    {
        if( iangles[i].countNormal().Angle(m_viewVector) > 90)
        {
            outTriangles.push_back( iangles[i]);
        }
    }
}

QVector<TriangleShared> Algoritm::readMailFile(const char* path)
{
    TriangleShared::Points.clear();
    srcTriangles.clear();

    TPoint3 tmp_pt;
    TriangleShared tmp_ ;
    int points_num,  aingles_num;

    std::ifstream ifs;
    ifs.open(path);
    ifs >> points_num >>  aingles_num;

    for(int i = 0; i < points_num; i++)
    {
        ifs >> tmp_pt.x >> tmp_pt.y >> tmp_pt.z;
        TriangleShared::Points.push_back(tmp_pt);
    }
    for(int i = 0; i <  aingles_num; i++)
    {
        ifs >> tmp_ .A >> tmp_ .B >> tmp_ .C;
        tmp_ .A--;
        tmp_ .B--;
        tmp_ .C--;
        if(tmp_ .A == tmp_ .B || tmp_ .A == tmp_ .C || tmp_ .B == tmp_ .C) continue;
        srcTriangles.push_back(tmp_ );
    }

    return srcTriangles;
}

void Algoritm::writeFile(const char* path, QVector<TriangleShared>&  iangles)
{
    std::ofstream f;
    f.open(path);

    if(!f.is_open())
        return;

    f <<   iangles.size() << std::endl;

    for(size_t i=0; i<  iangles.size(); i++)
    {
        TriangleShared t =  iangles.at(i);

        for(int j=0; j<3; j++)
        {
            //f << t.points[j].x << " " << t.points[j].y << " " << t.points[j].z << " ";
        }

        f << std::endl;
    }

    f.close();
}

void Algoritm::cleanPhase2()
{
    //clean - phase 2
    emit statusChanged("Clean phase 2: ");

    /*Multithread now! killNonShown();*/
//    CleanThread cleanThreads[THREAD_COUNT];
//    emit statusChanged( "using " + QS ing::number(THREAD_COUNT) +" threads: ");

//    int gap = outTriangles.size() / THREAD_COUNT;

//    for(int i = 0; i<THREAD_COUNT; i++)
//    {
//        cleanThreads[i].setJob(this, i*gap, i*gap+gap);
//        cleanThreads[i].start();
//    }

//    /* Waiting all threads */
//    bool jobFinished = false;
//    while(!jobFinished)
//    {
//        jobFinished =  ue;
//        for(int i=0; i<THREAD_COUNT; i++)
//        {
//            if( cleanThreads[i].isRunning() )
//            {
//                jobFinished = false;
//                break;
//            }
//        }

//        CleanThread::sleep(1);
//        emit processStateChanged(checkedTrianglesCount*98/srcTriangles.size());
//    }

    killNonShown(0, outTriangles.size());
    emit processStateChanged(98);
    emit statusChanged("done\n");
}


void Algoritm::killNonShown(int b, int e)//// ÑƒÐ±Ð¸Ð²Ð°ÐµÑ‚ Ð½ÐµÐ²Ð¸Ð´Ð¸Ð¼Ñ‹Ðµ Ñ‚Ñ€ÐµÑƒÐ³Ð¾Ð»ÑŒÐ½Ð¸ÐºÐ¸
{
    const int n = outTriangles.size()-1;

    for(int i = b; i < e+1; ++i)
    {
        if(i > n)
        {
            break;
        }
        if(outTriangles[i].dead == false)
            for(int j = i+1; j < outTriangles.size(); ++j)
            {
                if((!outTriangles[j].dead) && (!outTriangles[i].sameTriangle(outTriangles[j])))
                {
                    if(outTriangles[i].crossedRadius(outTriangles[j]))
                        if(outTriangles[i].crossedTriangles(outTriangles[j]))
                        {
                            if(outTriangles[i].distance < outTriangles[j].distance)
                                outTriangles[i].dead =  true;
                            else
                                outTriangles[j].dead =  true;
                        }
                }
            }

        checkedTrianglesCount++;
    }
}

//Ð´Ð¾Ð±Ð°Ð²Ð¸Ñ‚ÑŒ Ð² Ñ‚Ñ€ÐµÑƒÐ³Ð¾Ð»ÑŒÐ½Ð¸ÐºÐ¸ Ð¿Ð»Ð¾Ñ‰Ð°Ð´ÑŒ
//Ð½Ð¾Ñ€Ð¼Ð°Ð»Ð¸ Ðº Ð²ÐµÑ€ÑˆÐ¸Ð½Ð°Ð¼ Ñ‚Ñ€ÐµÑƒÐ³Ð¾Ð»ÑŒÐ½Ð¸ÐºÐ°

void Algoritm::translate()          //// Ñ„Ð¾Ñ€Ð¼Ð¸Ñ€ÑƒÐµÑ‚ Ð¼Ð°ÑÑÐ¸Ð² Ñ‚Ð¾Ñ‡ÐµÐº Ð² Ð¿Ñ€Ð¾ÑÑ‚Ñ€Ð°Ð½ÑÑ‚Ð²Ðµ 2D
{
    points_2D.clear();
    TPoint2 p;

    double a = cos(m_viewVector.Theta());
    double b = sin(m_viewVector.Theta());
    double c = cos(m_viewVector.Phi());
    double d = sin(m_viewVector.Phi());

    for(int i = 0; i < TriangleShared::RedirectPoints.size(); ++i)
    {
        p.x = -d*TriangleShared::RedirectPoints[i].x + c*TriangleShared::RedirectPoints[i].y;
        p.y = -a*c*TriangleShared::RedirectPoints[i].x + -a*d*TriangleShared::RedirectPoints[i].y + b*TriangleShared::RedirectPoints[i].z;
        points_2D.push_back(p);
    }
}

TPoint3 Algoritm::redirectPoint(TPoint3 p, TVector v)
{
    double d;
    TPoint3 tempPoint;

    d = v.x*p.x + v.y*p.y + v.z*p.z;

    tempPoint.x = p.x - d*v.x;
    tempPoint.y = p.y - d*v.y;
    tempPoint.z = p.z - d*v.z;

    return tempPoint;
}

bool Algoritm::redirectPoints()         //Ð¿Ñ€Ð¾ÐµÑ� Ð¸Ñ€ÑƒÐµÑ‚ Ñ‚Ð¾Ñ‡ÐºÐ¸ Ð² Ð¿Ð»Ð¾ÑÐºÐ¾ÑÑ‚ÑŒ Ñ‡ÐµÑ€ÐµÐ· (0,0,0) Ñ‚Ð¾Ñ‡ÐºÑƒ
{
    for(int i = 0; i < TriangleShared::Points.size(); ++i)
    {
        TriangleShared::RedirectPoints.push_back(redirectPoint(TriangleShared::Points[i], m_viewVector));
        TriangleShared::RedirectPoints.push_back(redirectPoint(TriangleShared::Points[i], m_viewVector));
    }
    return  true;
}

void Algoritm::findRadiuses()
{
    for(int i = 0; i< outTriangles.size(); ++i)
    {
        outTriangles[i].findRadius();
    }
}

void Algoritm::prepareTriangles()
{
    for(int i = 0; i < outTriangles.size(); i++)
    {
        outTriangles[i].dead = false;
        outTriangles[i].calcCenter();
        outTriangles[i].findDistance(m_viewVector);
    }
}

void Algoritm::redirectRadiuses()////Ð½Ð°Ñ…Ð¾Ð´Ð¸Ñ‚ Ñ� ÐµÐ½Ñ‚Ñ€ Ñ‚Ñ€ÐµÑƒÐ³Ð¾Ð»ÑŒÐ½Ð¸ÐºÐ° Ð´Ð»Ñ Ð½Ð°Ñ…Ð¾Ð¶Ð´ÐµÐ½Ð¸Ñ Ñ€Ð°Ð´Ð¸ÑƒÑÐ°
{
    for(int i = 0; i<outTriangles.size(); ++i)
    {
        outTriangles[i].redirectedCenter_3D = redirectPoint(outTriangles[i].center_3D,m_viewVector);
    }
}

double Algoritm::findMaxRadius()
{
    double max = 0;
    foreach(TriangleShared  iangle, outTriangles)
    {
        max = std::max(max,  iangle.radius);
    }

    return max;
}

void Algoritm::calculateWindowSize()
{

    double maxX, maxY, minX, minY;

    maxX=minX=points_2D[0].x;
    maxY=minY=points_2D[0].y;

    foreach(TPoint2 point, points_2D)
    {
        maxX = std::max(point.x, maxX);
        maxY = std::max(point.y, maxY);
        minX = std::min(point.x, minX);
        minY = std::min(point.y, minY);
    }

    double dx = maxX - minX;
    double dy = maxY - minY;
    double maxRadius = findMaxRadius();

    qDebug() << "Max radius is " << maxRadius << "Block size " << dx << "x" << dy;

           sx = (dx/maxRadius*NET_SIZE);

           sx = dx/sx;

           sy = (dy/maxRadius*NET_SIZE);
           sy = dy/sy;
}

TPoint3 Algoritm::countEVector(TPoint3 p)//e QVector (Ð¸Ð· ÑƒÑÐ»Ð¾Ð²Ð¸Ñ Ð¿ÐµÑ€Ð¿ÐµÐ½Ð´Ð¸ÐºÑƒÐ»ÑÑ€Ð½Ð¾ÑÑ‚Ð¸, Ð³Ð´Ðµ y=1 z=1)
{
    TPoint3 tp;
    tp.x = -(p.y + p.z)/p.x;
    tp.y = 1;
    tp.z = 1;
    return tp;
}

TPoint3 Algoritm::countProectionVector(TPoint3 p)//q Ð¿ÐµÑ€Ð¿ÐµÐ½Ð´Ð¸ÐºÑƒÐ»ÑÑ€
{
    TPoint3 tp;
    tp.x = p.x;
    tp.y = 0;
    tp.z = p.z;
    return tp;
}

TPoint3 Algoritm::countProectionComplexVector(TPoint3 p)//q*
{
    TPoint3 tp;
    tp.x = p.z;
    tp.y = 0;
    tp.z = -p.x;
    return tp;
}

TPoint3 Algoritm::countDecartVector(polarVector v)
{
    TPoint3 p;
    p.x=cos(v.fi)*sin(v.te);
    p.y=sin(v.fi)*sin(v.te);
    p.z=cos(v.te);
    return p;
}


