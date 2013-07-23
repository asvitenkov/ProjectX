#ifndef MODELMAKER_H
#define MODELMAKER_H

#include <QGLWidget>
#include <QVector>

#include "Model/TriangleShared.h"
#include "Model/Model.h"

#define OFFSET_RATIO 100

class ModelScene : public QGLWidget
{
    Q_OBJECT
    friend class Model;
public:
    explicit ModelScene(QWidget *parent = 0);
    ~ModelScene();

    void Update();
    void SetModel(Model* model);
    void setBorderPointsIndexes(QVector<unsigned int> &indexes);

    int getXRot();
    int getYRot();
    int getZRot();
    void setNormal(double tetta, double fi);
    void setDrawNormal(const bool enable);

private:
    void drawAxis();

    void drawBorderPoints();


    int xRot;
    int yRot;
    int zRot;


    float xOffset;
    float yOffset;

    double tetta;
    double fi;

    bool drawNormal;

    double scaling;
    double max;

    Model* m_model;

    void OnChange();
    void FindMax();

    QVector<TPoint3> *points;
    QVector<unsigned int> borderPointsIndexes;

    QPoint lastPos;

protected:
    /* Get 3d point */
    TPoint3 getPointFun(int index);
    /* Init OpenGL Enviroment */
    void initializeGL();
    /*
     *The paintGL() function is used to paint the contents of the scene onto the widget.
     *For widgets that only need to be decorated with
     *pure OpenGL content, we reimplement QGLWidget::paintGL() instead of reimplementing QWidget::paintEvent():
     */
    void paintGL();
    /*
     *The resizeGL() function is used to ensure that the OpenGL implementation renders the scene onto a viewport
     *that matches the size of the widget, using the correct transformation from 3D coordinates to 2D
     *viewport coordinates.
     */
    void resizeGL(int width, int height);
    /*
     *The mousePressEvent() function simply records the position of the mouse when a button is initially pressed:
     **/
    void mousePressEvent(QMouseEvent *event);
    /*
     *The mouseMoveEvent() function uses the previous location of the mouse
     *cursor to determine how much the object in the scene should be rotated, and in which direction
    */
    void mouseMoveEvent(QMouseEvent *event);
public slots:
    void setXRotation(int angle);
    void setYRotation(int angle);
    void setZRotation(int angle);
    void setScaling(int scaling);
    void setXOffset(int x);
    void setYOffset(int y);
signals:
    void xRotationChanged(int angle);
    void yRotationChanged(int angle);
    void zRotationChanged(int angle);
};

#endif // MODELMAKER_H
