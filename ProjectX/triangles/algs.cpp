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
    TriangleShared::Points2D.clear();
    checkedTrianglesCount = 0x0;
    time_t start_timestamp = time(NULL);

    outTriangles = srcTriangles;

    emit processStateChanged(0);

    //find 3D redirect points in triangles
    emit statusChanged("Preparing triangles:");
    prepareTriangles();//v
    emit statusChanged("done\n");

    /* redirect all points from TriangleShared::Points3D to TriangleShared::RedirectPoints vectors  AND RADIUSES*/
    redirectPoints();//v
    emit statusChanged("Reditect points: TriangleShared::RedirectPoints size is " + QString::number(TriangleShared::RedirectPoints.size()) +"\n");

    //find radiuses and distances from real triangle to its redirect
    emit statusChanged("Finding radiuses: ");
    findRadiuses();
    emit statusChanged("done\n");

    //find 2D redirect points ( pushing in TriangleShared::Points )
    translate();
    emit statusChanged("Translate: TriangleShared::Points size is " + QString::number(TriangleShared::Points2D.size()));
    emit statusChanged("done\n");

    cleanPhase2();

    //PROFIT
    QVector<TriangleShared> tempTriangles;
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

void Algoritm::setPolVector(TVector &v)
{
    polVec = v;
}

void Algoritm::cleanByNormal(const QVector<TriangleShared>& triangles)
{
    //cleaning backout triangles by normal
    for(size_t i = 0; i < triangles.size(); i++)
    {
        if(triangles[i].countNormal().Angle(polVec) > M_PI_2)
        {
            outTriangles.push_back(triangles[i]);
        }
    }
}

QVector<TriangleShared> Algoritm::readMailFile(const char* path)
{
    srcTriangles.clear();

    TPoint3 tmp_pt;
    TriangleShared tmp_tr;
    int A, B, C;
    int points_num, traingles_num;

    std::ifstream ifs;
    ifs.open(path);
    ifs >> points_num >> traingles_num;

    for(int i = 0; i < points_num; i++)
    {
        ifs >> tmp_pt.x >> tmp_pt.y >> tmp_pt.z;
        TriangleShared::Points3D.push_back(tmp_pt);
    }
    for(int i = 0; i < traingles_num; i++)
    {
        ifs >> A >> B >> C;

        if(A == B || A == C || B == C)
            continue;

        tmp_tr.Set(A-1, B-1, C-1);
        srcTriangles.push_back(tmp_tr);
    }

    return srcTriangles;
}

void Algoritm::writeFile(const char* path, QVector<TriangleShared>& triangles)
{
    std::ofstream f;
    f.open(path);

    if(!f.is_open())
        return;

    f <<  triangles.size() << std::endl;

    f.close();
}

void Algoritm::cleanPhase2()
{
    //clean - phase 2
    emit statusChanged("Clean phase 2: ");

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
                if((!outTriangles[j].dead) && (!outTriangles[i].sameTriangle(outTriangles[j])))
                {
                    if(outTriangles[i].crossedRadius(outTriangles[j]))
                        if(outTriangles[i].crossedTriangles(outTriangles[j]))
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

void Algoritm::translate()
{
    TriangleShared::Points2D.clear();
    TPoint2 p;
    double a = cos(polVec.Theta());
    double b = sin(polVec.Theta());
    double c = cos(polVec.Phi());
    double d = sin(polVec.Phi());
    for(int i = 0; i < TriangleShared::RedirectPoints.size(); ++i)
    {
        p.x = -d*TriangleShared::RedirectPoints[i].x + c*TriangleShared::RedirectPoints[i].y;
        p.y = -a*c*TriangleShared::RedirectPoints[i].x + -a*d*TriangleShared::RedirectPoints[i].y + b*TriangleShared::RedirectPoints[i].z;
        TriangleShared::Points2D.push_back(p);
    }
}

bool Algoritm::redirectPoints()
{
    for(int i = 0; i < TriangleShared::Points3D.size(); ++i)
    {
        TriangleShared::RedirectPoints.push_back(MathHelper::Redirect(&TriangleShared::Points3D[i], &polVec));
    }
    return true;
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
        outTriangles[i].CaclFileds(polVec);
        outTriangles[i].findDistance(polVec);
    }
}
