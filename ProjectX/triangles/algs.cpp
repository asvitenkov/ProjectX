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


QVector<triangle_t> Algoritm::manipulate()
{
    outTriangles.clear();
    redirectedPoints_3D.clear();
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

    //find 3D redirect points in triangles
    emit statusChanged("Preparing triangles:");
    prepareTriangles();//v
    emit statusChanged("done\n");

    /* redirect all points from points_3D to redirectedPoints_3D vectors  AND RADIUSES*/
    redirectPoints();//v
    redirectRadiuses();
    emit statusChanged("Reditect points: redirectedPoints_3D size is " + QString::number(redirectedPoints_3D.size()) +"\n");

    //find radiuses and distances from real triangle to its redirect
    emit statusChanged("Finding radiuses: ");
    findRadiuses();
    emit statusChanged("done\n");

    //find 2D redirect points ( pushing in points_2D )
    translate();
    emit statusChanged("Translate: points_2D size is " + QString::number(points_2D.size()));
    emit statusChanged("done\n");

    //finding window size
    this->calculateWindowSize();
    qDebug() << "Window size " << sx << " " << sy;

    cleanPhase2();

    //PROFIT
    QVector<triangle_t> tempTriangles;
    // making new QVector with result triangles
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
    emit statusChanged("Count of projection triangles: " +QString::number(outTriangles.size()));
    emit processStateChanged(100);

    return tempTriangles;
}

void Algoritm::setPolVector(polarVector &v)
{
    polVec = v;
}

void Algoritm::cleanByNormal(const QVector<triangle_t>& triangles)
{
    TPoint3 ptVec;
    ptVec.x = cos(polVec.fi)*sin(polVec.te);
    ptVec.y = sin(polVec.fi)*sin(polVec.te);
    ptVec.z = cos(polVec.te);

    //cleaning backout triangles by normal
    for(size_t i = 0; i < triangles.size(); i++)
    {
        if(countAngle(countNormal(triangles[i]), ptVec) > 90)
        {
            outTriangles.push_back(triangles[i]);
        }
    }
}

TPoint3 Algoritm::countNormal(triangle_t tr)
{
    TPoint3 a, b;
    a.x = points_3D[tr.B].x - points_3D[tr.A].x;
    a.y = points_3D[tr.B].y - points_3D[tr.A].y;
    a.z = points_3D[tr.B].z - points_3D[tr.A].z;
    b.x = points_3D[tr.C].x - points_3D[tr.A].x;
    b.y = points_3D[tr.C].y - points_3D[tr.A].y;
    b.z = points_3D[tr.C].z - points_3D[tr.A].z;
    return countVectorProduct(a, b);
}

TPoint3 Algoritm::countVectorProduct(TPoint3 a, TPoint3 b)
{
    TPoint3 c;
    c.x = a.y*b.z - a.z*b.y;
    c.y = a.z*b.x - a.x*b.z;
    c.z = a.x*b.y - a.y*b.x;
    return c;
}

double Algoritm::countScalarProduct(TPoint3 a, TPoint3 b)
{
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

double Algoritm::countLength(TPoint3 vec)
{
    return sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
}

double Algoritm::countAngle(TPoint3 a, TPoint3 b)
{
    return acos( countScalarProduct(a, b)/(countLength(a)*countLength(b)) )*180/M_PI;
}

QVector<triangle_t> Algoritm::readMailFile(const char* path)
{
    srcTriangles.clear();

    TPoint3 tmp_pt;
    triangle_t tmp_tr;
    int points_num, traingles_num;

    std::ifstream ifs;
    ifs.open(path);
    ifs >> points_num >> traingles_num;

    for(int i = 0; i < points_num; i++)
    {
        ifs >> tmp_pt.x >> tmp_pt.y >> tmp_pt.z;
        points_3D.push_back(tmp_pt);
    }
    for(int i = 0; i < traingles_num; i++)
    {
        ifs >> tmp_tr.A >> tmp_tr.B >> tmp_tr.C;
        tmp_tr.A--;
        tmp_tr.B--;
        tmp_tr.C--;
        if(tmp_tr.A == tmp_tr.B || tmp_tr.A == tmp_tr.C || tmp_tr.B == tmp_tr.C) continue;
        srcTriangles.push_back(tmp_tr);
    }

    return srcTriangles;
}

void Algoritm::writeFile(const char* path, QVector<triangle_t>& triangles)
{
    std::ofstream f;
    f.open(path);

    if(!f.is_open())
        return;

    f <<  triangles.size() << std::endl;

    for(size_t i=0; i< triangles.size(); i++)
    {
        struct triangle_t t = triangles.at(i);

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

//    /*Multithread now! killNonShown();*/
//    CleanThread cleanThreads[THREAD_COUNT];
//    emit statusChanged( "using " + QString::number(THREAD_COUNT) +" threads: ");

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
//        jobFinished = true;
//        for(int i=0; i<THREAD_COUNT; i++)
//        {
//            if( cleanThreads[i].isRunning() )
//            {
//                jobFinished = false;
//                break;
//            }
//        }

//        CleanThread::sleep(0.3);
//        emit processStateChanged(checkedTrianglesCount*98/srcTriangles.size());
//    }
    killNonShown(0, outTriangles.size());

    emit processStateChanged(98);
    emit statusChanged("done\n");
}


void Algoritm::killNonShown(int b, int e)
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
                if((!outTriangles[j].dead) && (!sameTriangle(outTriangles[i],outTriangles[j])))
                {
                    if(crossedRadius(outTriangles[i],outTriangles[j]))
                        if(crossedTriangles(outTriangles[i],outTriangles[j]))
                        {
                            if(outTriangles[i].distance < outTriangles[j].distance)
                                outTriangles[i].dead = true;
                            else
                                outTriangles[j].dead = true;
                        }
                }
            }

        checkedTrianglesCount++;
    }
}

bool Algoritm::crossedTriangles(triangle_t& t1, triangle_t& t2)
{
    if(crossedLines(points_2D[t1.A], points_2D[t1.B], points_2D[t2.A], points_2D[t2.B])) return true;
    if(crossedLines(points_2D[t1.A], points_2D[t1.C], points_2D[t2.A], points_2D[t2.B])) return true;
    if(crossedLines(points_2D[t1.B], points_2D[t1.C], points_2D[t2.A], points_2D[t2.B])) return true;

    if(crossedLines(points_2D[t1.A], points_2D[t1.B], points_2D[t2.A], points_2D[t2.C])) return true;
    if(crossedLines(points_2D[t1.A], points_2D[t1.C], points_2D[t2.A], points_2D[t2.C])) return true;
    if(crossedLines(points_2D[t1.B], points_2D[t1.C], points_2D[t2.A], points_2D[t2.C])) return true;

    if(crossedLines(points_2D[t1.A], points_2D[t1.B], points_2D[t2.C], points_2D[t2.B])) return true;
    if(crossedLines(points_2D[t1.A], points_2D[t1.C], points_2D[t2.C], points_2D[t2.B])) return true;
    if(crossedLines(points_2D[t1.B], points_2D[t1.C], points_2D[t2.C], points_2D[t2.B])) return true;
    return false;
}

bool Algoritm::crossedLines(TPoint2 p11, TPoint2 p12, TPoint2 p21, TPoint2 p22)////Ð¿ÐµÑ€ÐµÑÐµÐºÐ°ÑŽÑ‚ÑÑ Ð»Ð¸ Ð»Ð¸Ð½Ð¸Ð¸ Ñ‚Ñ€ÐµÑƒÐ³Ð¾Ð»ÑŒÐ½Ð¸ÐºÐ¾Ð²?
{
    double k1 = (p11.y - p12.y)/(p11.x - p12.x);
    double k2 = (p21.y - p22.y)/(p21.x - p22.x);
    double x = (p21.y - p11.y + k1*p11.x - k2*p21.x)/(k1 - k2);
    if(((x <= qMax(p11.x,p12.x))&&(x >= qMin(p12.x,p11.x)))&&((x <= qMax(p22.x,p21.x))&&(x >= qMin(p21.x,p22.x))))
        return true;
    return false;
}

bool Algoritm::crossedRadius(triangle_t& t1, triangle_t& t2)//Ð¿ÐµÑ€ÐµÑÐµÐºÐ°ÑŽÑ‚ÑÑ Ð»Ð¸ Ñ€Ð°Ð´Ð¸ÑƒÑÑ‹ Ñ‚Ñ€ÐµÑƒÐ³Ð¾Ð»ÑŒÐ½Ð¸ÐºÐ¾Ð²?
{
    double d;
    d = (t1.redirectedCenter_3D.x - t2.redirectedCenter_3D.x)*(t1.redirectedCenter_3D.x - t2.redirectedCenter_3D.x)+
            (t1.redirectedCenter_3D.y - t2.redirectedCenter_3D.y)*(t1.redirectedCenter_3D.y - t2.redirectedCenter_3D.y)+
            (t1.redirectedCenter_3D.z - t2.redirectedCenter_3D.z)*(t1.redirectedCenter_3D.z - t2.redirectedCenter_3D.z);
    if(d < ( t1.radius + t2.radius ))
        return true;
    return false;
}

bool Algoritm::sameTriangle(triangle_t& t1, triangle_t& t2)//ÑÐ²Ð»ÑÐµÑ‚ÑÑ Ð»Ð¸ ÑÐ¾ÑÐµÐ´Ð½Ð¸Ð¼ Ñ‚Ñ€ÐµÑƒÐ³Ð¾Ð»ÑŒÐ½Ð¸Ðº?
{
    if(t1.A==t2.A) return true;
    if(t1.A==t2.B) return true;
    if(t1.A==t2.C) return true;
    if(t1.B==t2.A) return true;
    if(t1.B==t2.B) return true;
    if(t1.B==t2.C) return true;
    if(t1.C==t2.A) return true;
    if(t1.C==t2.B) return true;
    if(t1.C==t2.C) return true;
    return false;
}
//Ð´Ð¾Ð±Ð°Ð²Ð¸Ñ‚ÑŒ Ð² Ñ‚Ñ€ÐµÑƒÐ³Ð¾Ð»ÑŒÐ½Ð¸ÐºÐ¸ Ð¿Ð»Ð¾Ñ‰Ð°Ð´ÑŒ
//Ð½Ð¾Ñ€Ð¼Ð°Ð»Ð¸ Ðº Ð²ÐµÑ€ÑˆÐ¸Ð½Ð°Ð¼ Ñ‚Ñ€ÐµÑƒÐ³Ð¾Ð»ÑŒÐ½Ð¸ÐºÐ°

void Algoritm::translate()          //// Ñ„Ð¾Ñ€Ð¼Ð¸Ñ€ÑƒÐµÑ‚ Ð¼Ð°ÑÑÐ¸Ð² Ñ‚Ð¾Ñ‡ÐµÐº Ð² Ð¿Ñ€Ð¾ÑÑ‚Ñ€Ð°Ð½ÑÑ‚Ð²Ðµ 2D
{
    points_2D.clear();
    TPoint2 p;
    double a = cos(polVec.te);
    double b = sin(polVec.te);
    double c = cos(polVec.fi);
    double d = sin(polVec.fi);
    for(int i = 0; i < redirectedPoints_3D.size(); ++i)
    {
        p.x = -d*redirectedPoints_3D[i].x + c*redirectedPoints_3D[i].y;
        p.y = -a*c*redirectedPoints_3D[i].x + -a*d*redirectedPoints_3D[i].y + b*redirectedPoints_3D[i].z;
        points_2D.push_back(p);
    }
}

double Algoritm::findRadius(triangle_t tr) //Ð½Ð°Ñ…Ð¾Ð´Ð¸Ñ‚ Ð¿ÑÐµÐ²Ð´Ð¾-Ñ€Ð°Ð´Ð¸ÑƒÑ Ð´Ð»Ñ Ð¿Ð¾Ð»Ð¸Ð³Ð¾Ð½Ð°
{
    double tempA=0;
    double tempB=0;
    double tempC=0;
    tempA=(redirectedPoints_3D[tr.A].x-tr.redirectedCenter_3D.x)*(redirectedPoints_3D[tr.A].x-tr.redirectedCenter_3D.x)+
            (redirectedPoints_3D[tr.A].y-tr.redirectedCenter_3D.y)*(redirectedPoints_3D[tr.A].y-tr.redirectedCenter_3D.y)+
            (redirectedPoints_3D[tr.A].z-tr.redirectedCenter_3D.z)*(redirectedPoints_3D[tr.A].z-tr.redirectedCenter_3D.z);
    tempB=(redirectedPoints_3D[tr.B].x-tr.redirectedCenter_3D.x)*(redirectedPoints_3D[tr.B].x-tr.redirectedCenter_3D.x)+
            (redirectedPoints_3D[tr.B].y-tr.redirectedCenter_3D.y)*(redirectedPoints_3D[tr.B].y-tr.redirectedCenter_3D.y)+
            (redirectedPoints_3D[tr.B].z-tr.redirectedCenter_3D.z)*(redirectedPoints_3D[tr.B].z-tr.redirectedCenter_3D.z);
    tempC=(redirectedPoints_3D[tr.C].x-tr.redirectedCenter_3D.x)*(redirectedPoints_3D[tr.C].x-tr.redirectedCenter_3D.x)+
            (redirectedPoints_3D[tr.C].y-tr.redirectedCenter_3D.y)*(redirectedPoints_3D[tr.C].y-tr.redirectedCenter_3D.y)+
            (redirectedPoints_3D[tr.C].z-tr.redirectedCenter_3D.z)*(redirectedPoints_3D[tr.C].z-tr.redirectedCenter_3D.z);
    if(tempA>=tempB){
        if(tempA>=tempC)
            return sqrt(tempA);
        else
            return sqrt(tempC);
    }
    else
    {
        if(tempB>=tempC)
            return sqrt(tempB);
        else
            return sqrt(tempC);
    }
}

TPoint3 Algoritm::redirectPoint(TPoint3 p, polarVector v) // Ð¿Ñ€Ð¾ÐµÑ� Ð¸Ñ€ÑƒÐµÑ‚ Ñ‚Ð¾Ñ‡ÐºÑƒ ÑÐ¼. Ð½Ð¸Ð¶Ðµ
{
    double d,A,B,C;
    TPoint3 tempPoint;
    A=cos(v.fi)*sin(v.te);
    B=sin(v.fi)*sin(v.te);
    C=cos(v.te);
    d=A*p.x+B*p.y+C*p.z;
    tempPoint.x = p.x - d*A;
    tempPoint.y = p.y - d*B;
    tempPoint.z = p.z - d*C;
    return tempPoint;
}

bool Algoritm::redirectPoints()         //Ð¿Ñ€Ð¾ÐµÑ� Ð¸Ñ€ÑƒÐµÑ‚ Ñ‚Ð¾Ñ‡ÐºÐ¸ Ð² Ð¿Ð»Ð¾ÑÐºÐ¾ÑÑ‚ÑŒ Ñ‡ÐµÑ€ÐµÐ· (0,0,0) Ñ‚Ð¾Ñ‡ÐºÑƒ
{
    for(int i = 0; i < points_3D.size(); ++i)
    {
        redirectedPoints_3D.push_back(redirectPoint(points_3D[i], polVec));
    }
    return true;
}

bool Algoritm::redirectPoints(int b, int e)         //Ð¿Ñ€Ð¾ÐµÑ� Ð¸Ñ€ÑƒÐµÑ‚ Ñ‚Ð¾Ñ‡ÐºÐ¸ Ð² Ð¿Ð»Ð¾ÑÐºÐ¾ÑÑ‚ÑŒ Ñ‡ÐµÑ€ÐµÐ· (0,0,0) Ñ‚Ð¾Ñ‡ÐºÑƒ
{
    for(int i = b; i < e+1; ++i)
    {
        redirectedPoints_3D.push_back(redirectPoint(points_3D[i], polVec));
    }
    return true;
}

void Algoritm::findRadiuses() // Ð½Ð°Ñ…Ð¾Ð´Ð¸Ñ‚ Ñ€Ð°Ð´Ð¸ÑƒÑÑ‹ Ð´Ð»Ñ Ñ‚Ñ€ÐµÑƒÐ³Ð¾Ð»ÑŒÐ½Ð¸ÐºÐ¾Ð²
{
    for(int i = 0; i< outTriangles.size(); ++i)
    {
        outTriangles[i].radius = findRadius(outTriangles[i]);
    }
}

float Algoritm::raschet(triangle_t t,	//Ñ‚Ñ€ÐµÑƒÐ³Ð¾Ð»ÑŒÐ½Ð¸Ðº                                   .... Ñ€Ð°ÑÑ‡Ñ‘Ñ‚
                        float R,    //Ñ€Ð°ÑÑÑ‚Ð¾ÑÐ½Ð¸Ðµ Ð´Ð¾ ÐºÐ¾Ð¾Ñ€Ð´Ð¸Ð½Ð°Ñ‚Ñ‹ 0.0 Ð¾Ñ‚ Ð°Ð½Ñ‚ÐµÐ½Ð½Ñ‹
                        int tol,	//Ð¿Ñ€Ð¸Ð·Ð½Ð°Ðº (0-Ð¼ÐµÑ‚Ð°Ð»Ð», -1 - Ð´Ð¸ÑÐ»ÐµÐºÑ‚Ñ€Ð¸Ðº)
                        float e1,	//ÑÐ»ÐµÐºÑ‚Ñ€Ð¸Ñ‡ÐµÑÐºÐ°Ñ Ð¿Ñ€Ð¾Ð½Ð¸Ñ� Ð°ÐµÐ¼Ð¾ÑÑ‚ÑŒ ÑÑ€ÐµÐ´Ñ‹
                        float m1,	//Ð¼Ð°Ð³Ð½Ð¸Ñ‚Ð½Ð°Ñ Ð¿Ñ€Ð¾Ð½Ð¸Ñ� Ð°ÐµÐ¼Ð¾ÑÑ‚ÑŒ ÑÑ€ÐµÐ´Ñ‹
                        float k0,	//Ð²Ð¾Ð»Ð½Ð¾Ð²Ð¾Ðµ Ñ‡Ð¸ÑÐ»Ð¾
                        float P0,	//Ð¿Ð¾Ð»ÑÑ€Ð¸Ð·Ð°Ñ� Ð¸Ð¸ Ð¿ÐµÑ€ÐµÐ´Ð°Ñ‚Ñ‡Ð¸ÐºÐ°
                        float P1)	//Ð¿Ð¾Ð»ÑÑ€Ð¸Ð·Ð°Ñ� Ð¸Ð¸ Ð¿Ñ€Ð¸ÐµÐ¼Ð½Ð¸ÐºÐ°
{
    float T;
    //COMPLEX j
    int j = -1;
    TPoint3 q,q1,a1,a2,a3,a21,a32,a13,a21_,a32_,a13_;
    q.x = -2*t.rVector.x;
    q.y = 0;
    q.z = -2*t.rVector.z;
    q1.x = -2*t.rVector.z;
    q1.y = 0;
    q1.z = 2*t.rVector.x;

    a1.x = t.aSystem.x - t.mSystem.x;
    a1.y = t.aSystem.y - t.mSystem.y;
    a1.z = t.aSystem.z - t.mSystem.z;
    a2.x = t.bSystem.x - t.mSystem.x;
    a2.y = t.bSystem.y - t.mSystem.y;
    a2.z = t.bSystem.z - t.mSystem.z;
    a3.x = t.cSystem.x - t.mSystem.x;
    a3.y = t.cSystem.y - t.mSystem.y;
    a3.z = t.cSystem.z - t.mSystem.z;

    a21.x = a2.x - a1.x;
    a21.y = a2.y - a1.y;
    a21.z = a2.z - a1.z;
    a32.x = a3.x - a2.x;
    a32.y = a3.y - a2.y;
    a32.z = a3.z - a2.z;
    a13.x = a1.x - a3.x;
    a13.y = a1.y - a3.y;
    a13.z = a1.z - a3.z;

    a21_.x = a2.x + a1.x;
    a21_.y = a2.y + a1.y;
    a21_.z = a2.z + a1.z;
    a32_.x = a3.x + a2.x;
    a32_.y = a3.y + a2.y;
    a32_.z = a3.z + a2.z;
    a13_.x = a1.x + a3.x;
    a13_.y = a1.y + a3.y;
    a13_.z = a1.z + a3.z;

    float Ei = 0;
    float D1, D2, D3, D;
    D1 = countScalarProduct(q1,a21)*((sin(0.5*k0*countScalarProduct(q1,a21)))/(0,5*k0*countScalarProduct(q1,a21)))
            *exp(-j*0.5*k0*countScalarProduct(q1,a21_));
    D2 = countScalarProduct(q1,a32)*((sin(0.5*k0*countScalarProduct(q1,a32)))/(0,5*k0*countScalarProduct(q1,a32)))
            *exp(-j*0.5*k0*countScalarProduct(q1,a32_));
    D3 = countScalarProduct(q1,a13)*((sin(0.5*k0*countScalarProduct(q1,a13)))/(0,5*k0*countScalarProduct(q1,a13)))
            *exp(-j*0.5*k0*countScalarProduct(q1,a13_));
    D = D1 + D2 + D3;
    Ei = exp(-j*2*k0*(R+t.distance))*(T)*D/
            (4*3.14*(R+t.distance)*(R+t.distance)*(t.rVector.x*t.rVector.x+t.rVector.z*t.rVector.z));
    return Ei;
}

void Algoritm::prepareTriangles()////Ð¿Ð¾Ð´Ð³Ð¾Ñ‚Ð°Ð²Ð»Ð¸Ð²Ð°ÐµÑ‚ Ñ‚Ñ€ÐµÑƒÐ³Ð¾Ð»ÑŒÐ½Ð¸ÐºÐ¸
{
    for(int i = 0; i < outTriangles.size(); i++)
    {
        outTriangles[i].dead = false;
        findMiddle(outTriangles[i]);
        findDistance(outTriangles[i], polVec);
        findVectors(outTriangles[i]);
        findLocals(outTriangles[i]);
    }
}

triangle_t Algoritm::findDistance(triangle_t& tr, polarVector v)///Ð½Ð°Ñ…Ð¾Ð´Ð¸Ñ‚ Ð´Ð¸ÑÑ‚Ð°Ð½Ñ� Ð¸ÑŽ Ñ‚Ñ€ÐµÑƒÐ³Ð¾Ð»ÑŒÐ½Ð¸ÐºÐ° (Ð´ÐµÐ»ÑŒÑ‚Ð°)!!!!Ð½ÑƒÐ¶Ð½Ð¾ +R Ð¾Ñ‚ Ð Ð›Ð¡!!!
{
    double A,B,C;
    A=cos(v.fi)*sin(v.te);
    B=sin(v.fi)*sin(v.te);
    C=cos(v.te);
    tr.distance=(A*tr.center_3D.x+B*tr.center_3D.y+C*tr.center_3D.z);
    return tr;
}

triangle_t Algoritm::findLocals(triangle_t &t)// Ð½Ð°Ñ…Ð¾Ð´Ð¸Ñ‚ Ñ‚Ð¾Ñ‡ÐºÐ¸ Ð¾Ñ‚Ð½Ð¾ÑÐ¸Ñ‚ÐµÐ»ÑŒÐ½Ð¾ Ð±Ð°Ð·Ð¸ÑÐ° Ð¿Ð¾Ð»Ð¸Ð³Ð¾Ð½Ð°
{
    t.aSystem = translateLocal(points_3D[t.A],t.cVector,t.nVector,t.bVector);
    t.bSystem = translateLocal(points_3D[t.B],t.cVector,t.nVector,t.bVector);
    t.cSystem = translateLocal(points_3D[t.C],t.cVector,t.nVector,t.bVector);
    t.mSystem.x = (t.aSystem.x + t.bSystem.x + t.cSystem.x)/3;
    t.mSystem.y = (t.aSystem.y + t.bSystem.y + t.cSystem.y)/3;
    t.mSystem.z = (t.aSystem.z + t.bSystem.z + t.cSystem.z)/3;
    TPoint3 pVec;
    pVec.x = cos(polVec.fi)*sin(polVec.te);
    pVec.y = sin(polVec.fi)*sin(polVec.te);
    pVec.z = cos(polVec.te);
    t.rVector = translateLocal(pVec,t.cVector,t.nVector,t.bVector);
    return t;
}

TPoint3 Algoritm::translateLocal(TPoint3 p, TPoint3 cVector, TPoint3 nVector, TPoint3 bVector)///Ð¿ÐµÑ€ÐµÐ²Ð¾Ð´Ð¸Ñ‚ Ñ‚Ð¾Ñ‡ÐºÑƒ p Ð² Ð±Ð°Ð·Ð¸Ñ Ñ‚Ñ€ÐµÑƒÐ³Ð¾Ð»ÑŒÐ½Ð¸ÐºÐ°
{
    TPoint3 tempP;
    tempP.x = cVector.x*p.x + cVector.y*p.y + cVector.z*p.z;
    tempP.y = nVector.x*p.x + nVector.y*p.y + nVector.z*p.z;
    tempP.z = bVector.x*p.x + bVector.y*p.y + bVector.z*p.z;
    return tempP;
}

void Algoritm::redirectRadiuses()////Ð½Ð°Ñ…Ð¾Ð´Ð¸Ñ‚ Ñ� ÐµÐ½Ñ‚Ñ€ Ñ‚Ñ€ÐµÑƒÐ³Ð¾Ð»ÑŒÐ½Ð¸ÐºÐ° Ð´Ð»Ñ Ð½Ð°Ñ…Ð¾Ð¶Ð´ÐµÐ½Ð¸Ñ Ñ€Ð°Ð´Ð¸ÑƒÑÐ°
{
    for(int i = 0; i<outTriangles.size(); ++i)
    {
        outTriangles[i].redirectedCenter_3D = redirectPoint(outTriangles[i].center_3D,polVec);
    }
}

double Algoritm::findMaxRadius()
{
    double max = 0;
    foreach(triangle_t triangle, outTriangles)
    {
        max = std::max(max, triangle.radius);
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

triangle_t Algoritm::findMiddle(triangle_t &t)//Ð½Ð°Ñ…Ð¾Ð´Ð¸Ñ‚ Ñ� ÐµÐ½Ñ‚Ñ€ Ð² Ñ‚Ñ€ÐµÑƒÐ³Ð¾Ð»ÑŒÐ½Ð¸ÐºÐµ
{
    t.center_3D.x=(points_3D[t.A].x+points_3D[t.B].x+points_3D[t.C].x)/3;
    t.center_3D.y=(points_3D[t.A].y+points_3D[t.B].y+points_3D[t.C].y)/3;
    t.center_3D.z=(points_3D[t.A].z+points_3D[t.B].z+points_3D[t.C].z)/3;
    return t;
}


TComplex Algoritm::reflectionCoefficient(triangle_t t,	//Ñ‚Ñ€ÐµÑƒÐ³Ð¾Ð»ÑŒÐ½Ð¸Ðº                                   .... Ñ€Ð°ÑÑ‡Ñ‘Ñ‚
                                   float R,    //Ñ€Ð°ÑÑÑ‚Ð¾ÑÐ½Ð¸Ðµ Ð´Ð¾ ÐºÐ¾Ð¾Ñ€Ð´Ð¸Ð½Ð°Ñ‚Ñ‹ 0.0 Ð¾Ñ‚ Ð°Ð½Ñ‚ÐµÐ½Ð½Ñ‹
                                   int   tol,	//Ð¿Ñ€Ð¸Ð·Ð½Ð°Ðº (0-Ð¼ÐµÑ‚Ð°Ð»Ð», -1 - Ð´Ð¸ÑÐ»ÐµÐºÑ‚Ñ€Ð¸Ðº)
                                   float e1,	//ÑÐ»ÐµÐºÑ‚Ñ€Ð¸Ñ‡ÐµÑÐºÐ°Ñ Ð¿Ñ€Ð¾Ð½Ð¸Ñ� Ð°ÐµÐ¼Ð¾ÑÑ‚ÑŒ ÑÑ€ÐµÐ´Ñ‹
                                   float m1,	//Ð¼Ð°Ð³Ð½Ð¸Ñ‚Ð½Ð°Ñ Ð¿Ñ€Ð¾Ð½Ð¸Ñ� Ð°ÐµÐ¼Ð¾ÑÑ‚ÑŒ ÑÑ€ÐµÐ´Ñ‹
                                   float k0,	//Ð²Ð¾Ð»Ð½Ð¾Ð²Ð¾Ðµ Ñ‡Ð¸ÑÐ»Ð¾
                                   float P0,	//Ð¿Ð¾Ð»ÑÑ€Ð¸Ð·Ð°Ñ� Ð¸Ð¸ Ð¿ÐµÑ€ÐµÐ´Ð°Ñ‚Ñ‡Ð¸ÐºÐ°
                                   float P1)	//Ð¿Ð¾Ð»ÑÑ€Ð¸Ð·Ð°Ñ� Ð¸Ð¸ Ð¿Ñ€Ð¸ÐµÐ¼Ð½Ð¸ÐºÐ°
{
    int Fv = 1;
    int Fh = -1;
    TPoint3 Tv, Th;
    TPoint3 y0,z0,p,rs;
    TPoint3 q, qComplex;
    float T;

    rs = countMultiple(t.rVector,-2);
    z0 = countOrt(countVectorProduct(t.rVector,t.nVector));
    y0 = countVectorProduct(z0,t.rVector);
    p = countVectorProduct(t.nVector,z0);
    q = countProectionVector(rs);
    qComplex = countProectionComplexVector(rs);

    Tv = countSum(countMultiple(p,(1 + Fv)*countScalarProduct(countVectorProduct(t.rVector,countEVector(t.rVector)),z0)),
    countMultiple(countVectorProduct(rs,z0),(-1)*(1 - Fv)*(countScalarProduct(countEVector(t.rVector),y0)*countScalarProduct(y0,p))));

    Th =countSum(
    countMultiple(z0,(Fh - 1)*countScalarProduct(countVectorProduct(t.rVector,countEVector(t.rVector)),y0)*countScalarProduct(y0,p)),
    countMultiple(countVectorProduct(p,rs),(-1)*(Fh + 1)*countScalarProduct(countEVector(t.rVector),z0)));

    T = countScalarProduct(countSum(Tv, Th), countOrt(rs));

    TPoint3 a21, a32, a13, a21_, a32_, a13_;
    a21 = countSum(t.bSystem,t.aSystem);
    a21_ = countSum(t.bSystem,countMultiple(t.aSystem,-1));
    a32 = countSum(t.cSystem,t.bSystem);
    a32_ = countSum(t.cSystem,countMultiple(t.bSystem,-1));
    a13 = countSum(t.aSystem,t.cSystem);
    a13_ = countSum(t.aSystem,countMultiple(t.cSystem,-1));

    float d1R, d2R, d3R, d1arg, d2arg, d3arg;
    d1R = -countScalarProduct(qComplex, a21_)*(sin(k0/2*countScalarProduct(q,a21_))/(k0/2*countScalarProduct(q,a21_)));
    d2R = -countScalarProduct(qComplex, a32_)*(sin(k0/2*countScalarProduct(q,a32_))/(k0/2*countScalarProduct(q,a32_)));
    d3R = -countScalarProduct(qComplex, a13_)*(sin(k0/2*countScalarProduct(q,a13_))/(k0/2*countScalarProduct(q,a13_)));
    d1arg = (-1)*k0/2*countScalarProduct(q,a21);
    d2arg = (-1)*k0/2*countScalarProduct(q,a32);
    d3arg = (-1)*k0/2*countScalarProduct(q,a13);

    TComplex D1(d1R * cos(d1arg), d1R * sin(d1arg));

    TComplex D2(d2R * cos(d2arg), d2R * sin(d2arg));

    TComplex D3(d3R * cos(d3arg), d3R * sin(d3arg));
    TComplex D(0, 0);

    D += D1;
    D += D2;
    D += D3;

    float argEi = (-1)*k0*(R + t.distance);
    TComplex Ei(cos(argEi), sin(argEi));

    Ei += 1/(4*M_PI)/(countLength(q)*countLength(q))*T;
    Ei += D;
    return Ei;
}

TPoint3 Algoritm::countOrt(TPoint3 p)
{
    TPoint3 tempPoint;
    float length = countLength(p);
    tempPoint.x = p.x / length;
    tempPoint.y = p.y / length;
    tempPoint.z = p.z / length;
    return tempPoint;
}
TPoint3 Algoritm::countSum(TPoint3 p1, TPoint3 p2)
{
    TPoint3 tp;
    tp.x = p1.x + p2.x;
    tp.y = p1.y + p2.y;
    tp.z = p1.z + p2.z;
    return tp;
}
TPoint3 Algoritm::countMultiple(TPoint3 p, float m)
{
    TPoint3 tp;
    tp.x = p.x * m;
    tp.y = p.y * m;
    tp.z = p.z * m;
    return tp;
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

triangle_t Algoritm::findVectors(triangle_t &t)//Ð½Ð°Ñ…Ð¾Ð´Ð¸Ñ‚ Ð²ÐµÐºÑ‚Ð¾Ñ€Ð° Ð´Ð»Ñ Ð½Ð°Ñ…Ð¾Ð¶Ð´ÐµÐ½Ð¸Ñ Ð±Ð°Ð·Ð¸ÑÐ° Ñ‚Ñ€ÐµÑƒÐ³Ð¾Ð»ÑŒÐ½Ð¸ÐºÐ°
{
    float d = sqrt((points_3D[t.A].x - t.center_3D.x)*(points_3D[t.A].x - t.center_3D.x)+
                   (points_3D[t.A].y - t.center_3D.y)*(points_3D[t.A].y - t.center_3D.y)+
                   (points_3D[t.A].z - t.center_3D.z)*(points_3D[t.A].z - t.center_3D.z));
    t.nVector = countNormal(t);
    t.cVector.x = (points_3D[t.A].x - t.center_3D.x)/d;
    t.cVector.y = (points_3D[t.A].y - -t.center_3D.y)/d;
    t.cVector.z = (points_3D[t.A].z - t.center_3D.z)/d;
    t.bVector = countVectorProduct(t.cVector,t.nVector);
    return t;
}

