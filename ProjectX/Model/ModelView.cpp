#include "ModelView.h"
#include "ModelScene.h"
#include <QSlider>
#include <QHBoxLayout>

ModelView::ModelView(QWidget* widget)
    : QWidget(widget)
{
    QSlider *xSlider = new QSlider(Qt::Vertical);
    QSlider *ySlider = new QSlider(Qt::Vertical);
    QSlider *zSlider = new QSlider(Qt::Vertical);
    QSlider *sSlider = new QSlider(Qt::Vertical);
    QSlider *axSlider = new QSlider(Qt::Vertical);
    QSlider *aySlider = new QSlider(Qt::Vertical);

    m_scene = new ModelScene();

    xSlider->setStyleSheet("QSlider { background-color: #33aa33; }");
    ySlider->setStyleSheet("QSlider { background-color: #aa3333; }");
    zSlider->setStyleSheet("QSlider { background-color: #3333aa; }");
    sSlider->setMaximum(+300.0);
    sSlider->setMinimum(-001.0);
    sSlider->setValue(100);

    xSlider->setMaximum(5740);
    ySlider->setMaximum(5740);
    zSlider->setMaximum(5740);

    axSlider->setMaximum(100);
    axSlider->setMinimum(-100);
    axSlider->setValue(0);

    aySlider->setMinimum(-100);
    aySlider->setMaximum(100);
    aySlider->setValue(0);

    this->setLayout(new QHBoxLayout());
    this->layout()->addWidget(axSlider);
    this->layout()->addWidget(aySlider);

    this->layout()->addWidget(m_scene);
    this->layout()->addWidget(xSlider);
    this->layout()->addWidget(ySlider);
    this->layout()->addWidget(zSlider);
    this->layout()->addWidget(sSlider);

    connect(xSlider, SIGNAL(valueChanged(int)), m_scene, SLOT(setXRotation(int)));
    connect(ySlider, SIGNAL(valueChanged(int)), m_scene, SLOT(setYRotation(int)));
    connect(zSlider, SIGNAL(valueChanged(int)), m_scene, SLOT(setZRotation(int)));
    connect(sSlider, SIGNAL(valueChanged(int)), m_scene, SLOT(setScaling(int)));

    connect(m_scene, SIGNAL(xRotationChanged(int)), xSlider, SLOT(setValue(int)));
    connect(m_scene, SIGNAL(yRotationChanged(int)), ySlider, SLOT(setValue(int)));
    connect(m_scene, SIGNAL(zRotationChanged(int)), zSlider, SLOT(setValue(int)));

    connect(axSlider, SIGNAL(valueChanged(int)), m_scene, SLOT(setXOffset(int)));
    connect(aySlider, SIGNAL(valueChanged(int)), m_scene, SLOT(setYOffset(int)));
}

ModelScene* ModelView::GetModelScene()
{
    return m_scene;
}
