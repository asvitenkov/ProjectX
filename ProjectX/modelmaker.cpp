#include "modelmaker.h"
#include "math.h"
#include <QDebug>
#include <QMouseEvent>

static void qNormalizeAngle(int &angle)
{
    while (angle < 0)
        angle += 360 * 16;
    while (angle > 360 * 16)
        angle -= 360 * 16;

}

ModelMaker::ModelMaker(QWidget *parent) :
    QGLWidget(parent), xOffset(0), yOffset(0)
{
    xRot = 0;
    yRot = 0;
    zRot = 0;

    drawNormal = false;

    scaling = 1;
}

void ModelMaker::setTriangles(const QVector<TriangleShared>& t)
{
    UsingTriangles = t;
    if(UsingTriangles.isEmpty()){
        updateGL();
        return;
    }

    TPoint3 p;
    for(int i=0; i<UsingTriangles.size(); i++)
    {
        p = getPointFun(UsingTriangles.at(i).A);
        this->max = qMax(p.x, this->max);
        this->max = qMax(p.y, this->max);
        this->max = qMax(p.z, this->max);

        p = getPointFun(UsingTriangles.at(i).B);
        this->max = qMax(p.x, this->max);
        this->max = qMax(p.y, this->max);
        this->max = qMax(p.z, this->max);

        this->max = qMax(getPointFun(UsingTriangles.at(i).C).x, this->max);
        this->max = qMax(getPointFun(UsingTriangles.at(i).C).y, this->max);
        this->max = qMax(getPointFun(UsingTriangles.at(i).C).z, this->max);
    }
    updateGL();
}

void ModelMaker::setPoints(QVector<TPoint3> *p)
{
    points = p;
}

void ModelMaker::setBorderPointsIndexes(QVector<unsigned int> &indexes)
{
    borderPointsIndexes = indexes;
    updateGL();
}

void ModelMaker::setMax(double _max)
{
    max = _max;
}

int ModelMaker::getXRot()
{
    return xRot;
}

int ModelMaker::getYRot()
{
    return yRot;
}

int ModelMaker::getZRot()
{
    return zRot;
}

void ModelMaker::setNormal(double _tetta, double _fi)
{
    tetta = _tetta;
    fi = _fi;
    updateGL();
}

void ModelMaker::setDrawNormal(const bool enable)
{
    drawNormal = enable;
}

void ModelMaker::drawAxis()
{

    //draw axis
    //z-axis
    qglColor(Qt::blue);
    glLineWidth(3.0);
    glBegin(GL_LINES);
    glVertex3d(0,0, 1.0);
    glVertex3d(0,0,0);
    glEnd();
    qglColor(Qt::red);
    //y-axis
    glBegin(GL_LINES);
    glVertex3d(0, 1.0 ,0);
    glVertex3d(0,0,0);
    glEnd();
    //x-axis
    qglColor(Qt::green);
    glBegin(GL_LINES);
    glVertex3d(1.0,0,0);
    glVertex3d(0,0,0);
    glEnd();
}

void ModelMaker::drawBorderPoints()
{
    if(borderPointsIndexes.size() <  0)
        return;

    for(int i=0; i<borderPointsIndexes.size(); i++)
    {
        unsigned int index = borderPointsIndexes.at(i);
        TPoint3 point = points->at(index);

        glPointSize(5.0f);
        qglColor(Qt::red);

        glBegin(GL_POINTS);
        glVertex3f(point.x/max, point.y/max, point.z/max);
        glEnd();
    }
}

TPoint3 ModelMaker::getPointFun(int index)
{
    return points->at(index);
}

void ModelMaker::initializeGL()
{
}

void ModelMaker::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glTranslatef(xOffset, yOffset , -10.0f);
    glRotatef(xRot / 16.0, 1.0, 0.0, 0.0);
    glRotatef(yRot / 16.0, 0.0, 1.0, 0.0);
    glRotatef(zRot / 16.0, 0.0, 0.0, 1.0);

    if(drawNormal) {
        qglColor(Qt::magenta);
        glBegin(GL_LINES);
        glVertex3d(sin(tetta)*cos(fi), sin(tetta)*sin(fi), cos(tetta));
        glVertex3d(0,0,0);
        glEnd();
    }

    qglColor(Qt::cyan);

    glLineWidth(1.0);
    glScalef(scaling, scaling, scaling);

    for(int i=0; i<UsingTriangles.size(); i++)
    {
        TriangleShared t = UsingTriangles.at(i);

        qglColor((t.dead)?(Qt::cyan):(Qt::blue));

        glBegin(GL_LINE_STRIP);
        glVertex3d((GLdouble) (getPointFun(t.A).x/max),
                   (GLdouble) (getPointFun(t.A).y / max), (GLdouble) (getPointFun(t.A).z / max));
        glVertex3d((GLdouble) (getPointFun(t.B).x/max),
                   (GLdouble) (getPointFun(t.B).y / max), (GLdouble) (getPointFun(t.B).z / max));
        glVertex3d((GLdouble) (getPointFun(t.C).x/max),
                   (GLdouble) (getPointFun(t.C).y / max), (GLdouble) (getPointFun(t.C).z / max));
        glVertex3d((GLdouble) (getPointFun(t.A).x/max),
                   (GLdouble) (getPointFun(t.A).y / max), (GLdouble) (getPointFun(t.A).z / max));
        glEnd();
    }
    drawBorderPoints();
    drawAxis();
}

void ModelMaker::resizeGL(int width, int height)
{
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float aspect = (float)this->width() / (float)this->height();
#ifdef QT_OPENGL_ES_1
    glOrthof(-aspect, aspect, -1.0, 1.0, 4.0, 15.0);
#else
    glOrtho(-aspect, aspect, -1.0, 1.0, 4.0, 15.0);
#endif
    glMatrixMode(GL_MODELVIEW);

}

void ModelMaker::mousePressEvent(QMouseEvent *event)
{
    lastPos = event->pos();
}

void ModelMaker::mouseMoveEvent(QMouseEvent *event)
{
    int dx = event->x() - lastPos.x();
    int dy = event->y() - lastPos.y();

    if (event->buttons() & Qt::LeftButton) {
        setXRotation(xRot + 8 * dy);
        setYRotation(yRot + 8 * dx);
    } else if (event->buttons() & Qt::RightButton) {
        setXRotation(xRot + 8 * dy);
        setZRotation(zRot + 8 * dx);
    }
    lastPos = event->pos();
}

void ModelMaker::setXRotation(int angle)
{
    qNormalizeAngle(angle);

    if (angle != xRot) {
        xRot = angle;
        emit xRotationChanged(angle);
        updateGL();
    }
}

void ModelMaker::setYRotation(int angle)
{
    qNormalizeAngle(angle);

    if (angle != yRot) {
        yRot = angle;
        emit yRotationChanged(angle);
        updateGL();
    }
}

void ModelMaker::setZRotation(int angle)
{
    qNormalizeAngle(angle);

    if (angle != zRot) {
        zRot = angle;
        emit zRotationChanged(angle);
        updateGL();
    }
}

void ModelMaker::setScaling(int _scaling)
{
    scaling = ( _scaling ) / 100.0;
    updateGL();
}

void ModelMaker::setXOffset(int x)
{
    xOffset = (double) x / OFFSET_RATIO;
    updateGL();
}

void ModelMaker::setYOffset(int y)
{
    yOffset = (float) y / OFFSET_RATIO;
    updateGL();
}
