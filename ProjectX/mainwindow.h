#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "modelviewer.h"
#include "processthread.h"
#include "triangles/algs.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    ModelViewer inputModelMaker;
private:
    QVector<TriangleShared> triangles;
    Ui::MainWindow *ui;
    QVector<TriangleShared> data;
    Algoritm *alg;
    ProcessThread *pthread;

public slots:
    /* Обработка открытых файлов */
    void openHandler();
    /* Процес обработки проекции */
    void processHandler();
    /* Поиск вершин, принадлежащих граням */
    void procesBorderLinesFile(QString filePath);
    /* Сохранение файла */
    void saveHandler();
    /* Угол поворота изменен */
    void rotationAngleChanged();
    /*Изменен статус работы*/
    void stateChanged(QString state);
    /*Изменился процент выполнения*/
    void processStatusChanged(int);
};

#endif // MAINWINDOW_H
