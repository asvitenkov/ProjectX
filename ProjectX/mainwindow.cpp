#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "triangles/algs.h"
#include "edgeselector.h"
#include "math.h"


#include <QString>
#include <QHBoxLayout>
#include <QFileDialog>
#include <QDebug>
#include <QMessageBox>
#include <QDate>

#define win32 true


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    alg(0x0),
    pthread(0x0)
{
    ui->setupUi(this);
    ui->originalViewTab->setLayout(new QHBoxLayout());
    ui->originalViewTab->layout()->addWidget(&inputModelMaker);
#ifndef win32
    ui->originalViewTab->layout()->addWidget(&outputModelMaker);
#endif

    connect(ui->actionOpen, SIGNAL(triggered()),
            this, SLOT(openHandler()));
    connect(ui->showButton, SIGNAL(clicked()), this, SLOT(processHandler()));
    connect(ui->showButton, SIGNAL(clicked()), this, SLOT(processHandler()));
    connect(ui->actionSave, SIGNAL(triggered()), this, SLOT(saveHandler()));

    connect( &inputModelMaker.getModelMaker(), SIGNAL(xRotationChanged(int)), this, SLOT(rotationAngleChanged()) );
    connect( &inputModelMaker.getModelMaker(), SIGNAL(zRotationChanged(int)), this, SLOT(rotationAngleChanged()) );
    connect( &inputModelMaker.getModelMaker(), SIGNAL(yRotationChanged(int)), this, SLOT(rotationAngleChanged()) );

    pthread = new ProcessThread();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::openHandler()
{
    QString filename = QFileDialog::getOpenFileName(this,
                                                    "Select triangles file",
                                                     QDir::homePath(),
                                                    "Triangles file (*.mail);; Extended information (*)");
    if(filename.isEmpty())
        return;

    if(filename.endsWith(".mail"))
    {
        if(alg != 0x0) {
            disconnect(alg, SIGNAL(processStateChanged(int)));
            disconnect(alg, SIGNAL(statusChanged(QString)));

            inputModelMaker.getModelMaker().setPoints(0x0);
            inputModelMaker.getModelMaker().setTriangles(QVector<triangle_t>());
            delete alg;
        }

        alg =  new Algoritm();

        connect(alg, SIGNAL(processStateChanged(int)), this, SLOT(processStatusChanged(int)), Qt::QueuedConnection);
        connect(alg, SIGNAL(statusChanged(QString)), this, SLOT(stateChanged(QString)), Qt::QueuedConnection);

        alg->readMailFile(filename.toStdString().data());
        inputModelMaker.getModelMaker().setPoints(&alg->points_3D);
        inputModelMaker.getModelMaker().setTriangles(alg->srcTriangles);
    }
    else
    {
        if(alg != 0x0)
            procesBorderLinesFile(filename);
        else
        {
            QMessageBox::warning(this, tr("Warning"), "Open please mail file for first");
            return;
        }
    }
}

void MainWindow::processHandler()
{
    if(alg->srcTriangles.size() == 0 )
        QMessageBox::warning(this, tr("Warning"), tr("Empty triangles array"));

    this->setCursor(Qt::WaitCursor);

    QDate startDate = QDate::currentDate();

    polarVector polVec;

    polVec.r = 1;
    polVec.te = ui->aSpinBox->value();
    polVec.fi = ui->fiSphinBox->value();

    alg->setPolVector(polVec);

    pthread->setAlgs(alg);
    pthread->start();

    inputModelMaker.getModelMaker().setDrawNormal(true);
    inputModelMaker.getModelMaker().setNormal(ui->aSpinBox->value(), ui->fiSphinBox->value());
}

void MainWindow::procesBorderLinesFile(QString filename)
{
    EdgeSelector selector;
    try {
        selector.readFile(filename.toStdString().data());
        QVector<unsigned int> edgePoints = selector.findEdgePoints(alg->points_3D);
        qDebug() << "Selected " << edgePoints.size() << " points";
        QMessageBox::information(this, tr("Information"), tr("Border selecting done. Selected: ")
                                 + QString::number(edgePoints.size()) +" points.");

        inputModelMaker.getModelMaker().setBorderPointsIndexes(edgePoints);
    } catch (Exception e)
    {
        QMessageBox::warning(this, "Exception", e.getMessage().data());
    }
}

void MainWindow::saveHandler()
{
    QString filePath = QFileDialog::getSaveFileName(this, "Saving triangles file", QDir::homePath(), "Triangles file (*.mail)");
    if(filePath.isEmpty())
        return;
    this->setCursor(Qt::WaitCursor);
    alg->writeFile(filePath.toStdString().data(), data);
    this->setCursor(Qt::ArrowCursor);
}

void MainWindow::rotationAngleChanged()
{
    int xRot = inputModelMaker.getModelMaker().getXRot() / 16;
    int yRot = inputModelMaker.getModelMaker().getYRot() / 16;
    int zRot = inputModelMaker.getModelMaker().getZRot() / 16;
}

void MainWindow::stateChanged(QString state)
{
    ui->logTextEdit->append(state);
}

void MainWindow::processStatusChanged(int x)
{
    ui->progressBar->setValue(x);

    if(x > 10)
    {
#ifdef win32
        inputModelMaker.getModelMaker().setPoints(&alg->points_3D);
        inputModelMaker.getModelMaker().setTriangles(alg->outTriangles);
#else
        outputModelMaker.getModelMaker().setPoints(&alg->points_3D);
        outputModelMaker.getModelMaker().setTriangles(alg->outTriangles);
#endif
        this->setCursor(Qt::ArrowCursor);
    }
}

