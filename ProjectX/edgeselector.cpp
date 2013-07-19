#include "edgeselector.h"

#include <fstream>
#include <iostream>

#include <QRegExp>
#include <QDebug>

EdgeSelector::EdgeSelector()
{
}

QVector<unsigned int> EdgeSelector::findEdgePoints(const QVector<TPoint3> &srcPoints) throw (NoLinesException)
{
    QVector<unsigned int> outPoints;

    if(lines.size() == 0 )
        throw NoLinesException("Edge lines QVector is empty.");

    for(int i=0; i<srcPoints.size(); i++)
    {
        for(int j = 0; j < lines.size(); j++)
        {

            double x = points[i].x;
            double x1 = points[lines[j].A].x;
            double x2 = points[lines[j].B].x;

            double y = points[i].y;
            double y1 = points[lines[j].A].y;
            double y2 = points[lines[j].B].y;

            double z = points[i].z;
            double z1 = points[lines[j].A].z;
            double z2 = points[lines[j].B].z;

            double f1 = (x-x1)*(y2-y1)*(z2-z1);
            double f2 = (y-y1)*(x2-x1)*(z2-z1);
            double f3 = (z-z1)*(x2-x1)*(y2-y1);

            qDebug() << f1 << ' ' << f2 << ' ' << f3;

            if(DABS(f1-f2) < ACCURACY && DABS(f1-f3) < ACCURACY && DABS(f2-f3) < ACCURACY && (f1 || f2 || f3 ) )
            {
                outPoints.push_back(i);
                break;
            }
        }
    }

    return outPoints;
}

void EdgeSelector::readFile(const char *path)
{
    TPoint3 tmp_pt;

    std::fstream file;
    file.open(path);

    std::cout << "Opening file: " << path << std::endl;

    if(!file.is_open())
        throw Exception("Can't open extended file.");

    std::string line;

    QRegExp pointExp("^Point\\((\\d+)\\)\\s=\\s\\{(-?[0-9\\.]+),\\s(-?[0-9\\.]+),\\s(-?[0-9\\.]+),\\s(-?[0-9\\.]+)\\}");
    QRegExp lineExp("^Line\\((\\d+)\\)\\s=\\s\\{([0-9]+),\\s([0-9]+)\\}");

    while( !file.eof() )
    {
        std::getline(file, line);

        if(pointExp.indexIn( line.data()) != -1 )
        {
            int num = pointExp.cap(1).toInt();
            double x = pointExp.cap(2).toDouble();
            double y = pointExp.cap(3).toDouble();
            double z = pointExp.cap(4).toDouble();
            double a = pointExp.cap(5).toDouble();

            TPoint3 currentPoint(x, y, z);

            points.push_back(currentPoint);
        }
        else if( lineExp.indexIn(line.data()) != -1)
        {
            int num = lineExp.cap(1).toInt();
            double a = lineExp.cap(2).toInt()-1;
            double b = lineExp.cap(3).toInt()-1;

            struct line currentLine;
            currentLine.A = a;
            currentLine.B = b;
            lines.push_back(currentLine);

            std::cout << "Parsed line " << num << " point " << a << " " << b << std::endl;
        }
    }
}
