#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include "processthread.h"
#include "Algo/Algoritm.h"

#include "Model/ModelView.h"
#include "Model/Model.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    ModelView m_view, outputModelMaker;

private:
    QVector<TriangleShared> triangles;
    Ui::MainWindow *ui;
    QVector<TriangleShared> data;
    Algoritm *m_algo;
    ProcessThread *pthread;

    Model* m_model;

public slots:
    void openHandler();
    void processHandler();
    void procesBorderLinesFile(QString filePath);
    void saveHandler();
    void rotationAngleChanged();
    void stateChanged(QString state);
    void processStatusChanged(int);
};

#endif // MAINWINDOW_H
