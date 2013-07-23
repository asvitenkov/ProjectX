#ifndef MODELVIEWER_H
#define MODELVIEWER_H

#include <QWidget>

class ModelScene;

class ModelView: public QWidget
{
    Q_OBJECT
public:
    ModelView(QWidget* parent = 0x0);
    ModelScene* GetModelScene();
private:
    ModelScene* m_scene;
};

#endif // MODELVIEWER_H
