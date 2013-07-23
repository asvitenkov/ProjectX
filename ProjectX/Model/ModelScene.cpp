#include "ModelScene.h"
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

ModelScene::ModelScene(QWidget *parent) :
    QGLWidget(parent), xOffset(0), yOffset(0)
{
    xRot = 0;
    yRot = 0;
    zRot = 0;

    drawNormal = false;

    scaling = 1;

    m_model = new Model();
}

ModelScene::~ModelScene()
{
    if(m_model)
    {
        delete m_model;
    }
}

void ModelScene::SetModel(Model* model)
{
    if(m_model)
    {
        delete m_model;
    }

    m_model = model;
    OnChange();
}


void ModelScene::Update()
{
    updateGL();
}

void ModelScene::FindMax()
{
    for(int i=0; i<m_model->GetTriangles().size(); i++)
    {
        const TPoint3& p1 = m_model->GetTriangles()[i].p1();
        const TPoint3& p2 = m_model->GetTriangles()[i].p2();
        const TPoint3& p3 = m_model->GetTriangles()[i].p3();

        this->max = qMax(p1.x, this->max);
        this->max = qMax(p1.y, this->max);
        this->max = qMax(p1.z, this->max);

        this->max = qMax(p2.x, this->max);
        this->max = qMax(p2.y, this->max);
        this->max = qMax(p2.z, this->max);

        this->max = qMax(p3.x, this->max);
        this->max = qMax(p3.y, this->max);
        this->max = qMax(p3.z, this->max);
    }
}

void ModelScene::OnChange()
{
    FindMax();
    updateGL();
}

void ModelScene::setBorderPointsIndexes(QVector<unsigned int> &indexes)
{
    borderPointsIndexes = indexes;
    updateGL();
}

int ModelScene::getXRot()
{
    return xRot;
}

int ModelScene::getYRot()
{
    return yRot;
}

int ModelScene::getZRot()
{
    return zRot;
}

void ModelScene::setNormal(double _tetta, double _fi)
{
    tetta = _tetta;
    fi = _fi;
    updateGL();
}

void ModelScene::setDrawNormal(const bool enable)
{
    drawNormal = enable;
}

void ModelScene::drawAxis()
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

void ModelScene::drawBorderPoints()
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

TPoint3 ModelScene::getPointFun(int index)
{
    return points->at(index);
}

void ModelScene::initializeGL()
{
}

void ModelScene::paintGL()
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

    for(int i=0; i<m_model->m_triangles.size(); i++)
    {
        TriangleShared& t = m_model->m_triangles[i];

        if(!t.IsVisible())
            continue;

        glBegin(GL_LINE_STRIP);
        glVertex3d((GLdouble) (t.p1().x/max),
                   (GLdouble) (t.p1().y / max), (GLdouble) (t.p1().z / max));
        glVertex3d((GLdouble) (t.p2().x/max),
                   (GLdouble) (t.p2().y / max), (GLdouble) (t.p2().z / max));
        glVertex3d((GLdouble) (t.p3().x/max),
                   (GLdouble) (t.p3().y / max), (GLdouble) (t.p3().z / max));
        glVertex3d((GLdouble) (t.p1().x/max),
                   (GLdouble) (t.p1().y / max), (GLdouble) (t.p1().z / max));
        glEnd();
    }
    drawBorderPoints();
    drawAxis();
}

void ModelScene::resizeGL(int width, int height)
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

void ModelScene::mousePressEvent(QMouseEvent *event)
{
    lastPos = event->pos();
}

void ModelScene::mouseMoveEvent(QMouseEvent *event)
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

void ModelScene::setXRotation(int angle)
{
    qNormalizeAngle(angle);

    if (angle != xRot) {
        xRot = angle;
        emit xRotationChanged(angle);
        updateGL();
    }
}

void ModelScene::setYRotation(int angle)
{
    qNormalizeAngle(angle);

    if (angle != yRot) {
        yRot = angle;
        emit yRotationChanged(angle);
        updateGL();
    }
}

void ModelScene::setZRotation(int angle)
{
    qNormalizeAngle(angle);

    if (angle != zRot) {
        zRot = angle;
        emit zRotationChanged(angle);
        updateGL();
    }
}

void ModelScene::setScaling(int _scaling)
{
    scaling = ( _scaling ) / 100.0;
    updateGL();
}

void ModelScene::setXOffset(int x)
{
    xOffset = (double) x / OFFSET_RATIO;
    updateGL();
}

void ModelScene::setYOffset(int y)
{
    yOffset = (float) y / OFFSET_RATIO;
    updateGL();
}
