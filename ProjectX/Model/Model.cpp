#include "Model.h"

#include <fstream>

Model::Model(QObject* widget)
    : QObject(widget)
{
}

void Model::Create(const QString& filepath)
{
    m_triangles.clear();
    m_points.clear();

    TPoint3 pointTmp;
    TriangleShared triangleTmp(this);

    int A, B, C;
    int pointAmount, trainglesAmount;

    std::ifstream ifs;
    ifs.open(filepath.toStdString().c_str());
    ifs >> pointAmount >> trainglesAmount;

    for(int i = 0; i < pointAmount; ++i)
    {
        ifs >> pointTmp.x >> pointTmp.y >> pointTmp.z;
        AddPoint(pointTmp);
    }

    for(int i = 0; i < trainglesAmount; ++i)
    {
        ifs >> A >> B >> C;

        if(A == B || A == C || B == C)
            continue;

        triangleTmp.Set(A-1, B-1, C-1);
        m_triangles.push_back(triangleTmp);
    }
}

void Model::AddPoint(const TPoint3& p)
{
    m_points.push_back(p);
}

QVector<TriangleShared>& Model::GetTriangles()
{
    return m_triangles;
}

const QVector<TPoint3>& Model::GetPoints() const
{
    return m_points;
}

void Model::SetTriangles(const QVector<TriangleShared>& newTr)
{
    m_triangles = newTr;
}

void Model::SetPoints(const QVector<TPoint3>& newP)
{
    m_points = newP;
}
