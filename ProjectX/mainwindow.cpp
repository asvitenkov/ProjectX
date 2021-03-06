#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "Algo/Algoritm.h"
#include "edgeselector.h"
#include "math.h"


#include <QString>
#include <QHBoxLayout>
#include <QFileDialog>
#include <QDebug>
#include <QMessageBox>
#include <QDate>
#include "calculatedialog.h"

#define win32 true


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
    , m_algo(NULL)
    , pthread(NULL)
{
    ui->setupUi(this);
    ui->originalViewTab->setLayout(new QHBoxLayout());
    ui->originalViewTab->layout()->addWidget(&m_view);

    CalculateDialog *dialog = new CalculateDialog(&m_view);

    ui->projectionViewWidget->addTab(dialog, tr("Calculate"));

    connect(ui->actionOpen, SIGNAL(triggered()),
            this, SLOT(openHandler()));
    connect(ui->showButton, SIGNAL(clicked()), this, SLOT(processHandler()));
    connect(ui->showButton, SIGNAL(clicked()), this, SLOT(processHandler()));
    connect(ui->actionSave, SIGNAL(triggered()), this, SLOT(saveHandler()));
    connect(ui->actionReset, SIGNAL(triggered()), this, SLOT(resetHandler()));

    connect( m_view.GetModelScene(), SIGNAL(xRotationChanged(int)), this, SLOT(rotationAngleChanged()) );
    connect( m_view.GetModelScene(), SIGNAL(zRotationChanged(int)), this, SLOT(rotationAngleChanged()) );
    connect( m_view.GetModelScene(), SIGNAL(yRotationChanged(int)), this, SLOT(rotationAngleChanged()) );

    ui->aSpinBox->setMinimum(0);
    ui->aSpinBox->setMaximum(360);
    ui->fiSphinBox->setMinimum(0);
    ui->fiSphinBox->setMaximum(360.0);

    pthread = new ProcessThread();

    m_algo = new Algoritm(pthread);
    m_algo->SetModel(m_view.GetModelScene()->GetModel());

    connect(m_algo, SIGNAL(processStateChanged(int)), this, SLOT(processStatusChanged(int)), Qt::QueuedConnection);
    connect(m_algo, SIGNAL(statusChanged(QString)), this, SLOT(stateChanged(QString)), Qt::QueuedConnection);

    m_timer = new QTimer(this);
    m_timer->setInterval(1000);
    m_timer->setSingleShot(false);
    connect(m_timer, SIGNAL(timeout()), m_algo, SLOT(NotifyAboutState()));

    m_timer->start();
}

MainWindow::~MainWindow()
{
    delete pthread;
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
        Model* model = m_view.GetModelScene()->GetModel();
        model->Create(filename);
        m_view.GetModelScene()->ReChange();
    }
    else
    {
        procesBorderLinesFile(filename);
    }
}

void MainWindow::processHandler()
{
    setCursor(Qt::WaitCursor);

    TVector polVec;

    polVec.Set(ui->aSpinBox->value()*M_PI/180, ui->fiSphinBox->value()*M_PI/180, 1);

    m_algo->SetSightVector(polVec);
    m_view.GetModelScene()->setDrawNormal(true);
    m_view.GetModelScene()->setNormal(polVec.Theta(), polVec.Phi());

    pthread->setAlgs(m_algo);
    pthread->start();
}

void MainWindow::procesBorderLinesFile(QString filename)
{
    EdgeSelector selector;
    try {
        selector.readFile(filename.toStdString().data());
        QVector<unsigned int> edgePoints = selector.findEdgePoints(TriangleShared::Points3D);
        qDebug() << "Selected " << edgePoints.size() << " points";
        QMessageBox::information(this, tr("Information"), tr("Border selecting done. Selected: ")
                                 + QString::number(edgePoints.size()) +" points.");

        m_view.GetModelScene()->setBorderPointsIndexes(edgePoints);
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
    //m_algo->writeFile(filePath.toStdString().data(), data);
    this->setCursor(Qt::ArrowCursor);
}

void MainWindow::rotationAngleChanged()
{
    int xRot = m_view.GetModelScene()->getXRot() / 16;
    int yRot = m_view.GetModelScene()->getYRot() / 16;
    int zRot = m_view.GetModelScene()->getZRot() / 16;
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
        this->setCursor(Qt::ArrowCursor);
        m_view.GetModelScene()->Update();
    }
}

void MainWindow::resetHandler()
{
    Model* model = m_view.GetModelScene()->GetModel();
    for(int i = 0; i<model->GetTriangles().size(); i++)
    {
        model->GetTriangles()[i].SetVisible(true);
    }
}
