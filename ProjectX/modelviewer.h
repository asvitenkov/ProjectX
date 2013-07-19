#ifndef MODELVIEWER_H
#define MODELVIEWER_H

#include <QWidget>
#include "modelmaker.h"

class ModelViewer: public QWidget
{
    Q_OBJECT
public:
    ModelViewer(QWidget* parent = 0x0);
    ModelMaker& getModelMaker();
private:
    ModelMaker modelMaker;
};

#endif // MODELVIEWER_H
