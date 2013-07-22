#ifndef MODELMAKER_H
#define MODELMAKER_H

#include <QGLWidget>
#include <QVector>
#include "triangles/structs.h"

#define OFFSET_RATIO 100

class ModelMaker : public QGLWidget
{
    Q_OBJECT
public:
    explicit ModelMaker(QWidget *parent = 0);

    /* Установка вектора треугольников для прорисовки */
    void setTriangles(const QVector<triangle_t>& t);
    /* Установка точек\вершин треугольника */
    void setPoints(QVector<TPoint3> *p);
    /* Установка точек лежащих на гранях для прорисовки */
    void setBorderPointsIndexes(QVector<unsigned int> &indexes);
    /* Установка максимального значения координат при прорисовке */
    void setMax(double _max);

    //получение угла поворота относительно осей
    int getXRot();
    int getYRot();
    int getZRot();
    //установка координат нормали, с которой происходит просмотр
    void setNormal(double tetta, double fi);
    //установка разрешения прорисовки нормали
    void setDrawNormal(const bool enable);

private:
    /* Прорисовка  координатных осей */
    void drawAxis();
    /* Прорисовка точек лежащих на гранях */
    void drawBorderPoints();

    /* Коэффициенты поворота по осям */
    int xRot;
    int yRot;
    int zRot;

    /* Коэффициент сдвига центра по осям Ox и Oy */
    float xOffset;
    float yOffset;
    /* Координаты нормали зрения */
    double tetta;
    double fi;
    /* Флаг разрешения прорисовки нормали */
    bool drawNormal;

    /* Коэффициент маштаба */
    double scaling;
    double max;


    /* Вектора треугольников, вершин, вершин принадлежащих граням */
    QVector<triangle_t> UsingTriangles;
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
