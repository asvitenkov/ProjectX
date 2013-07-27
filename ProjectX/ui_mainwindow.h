/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created: Sat 27. Jul 15:31:40 2013
**      by: Qt User Interface Compiler version 4.8.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QGroupBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QProgressBar>
#include <QtGui/QPushButton>
#include <QtGui/QStatusBar>
#include <QtGui/QTabWidget>
#include <QtGui/QTextEdit>
#include <QtGui/QToolBar>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *actionOpen;
    QAction *actionSave;
    QAction *actionReset;
    QWidget *centralWidget;
    QVBoxLayout *verticalLayout;
    QTabWidget *projectionViewWidget;
    QWidget *originalViewTab;
    QWidget *doubleViewTab;
    QVBoxLayout *verticalLayout_2;
    QTextEdit *logTextEdit;
    QGroupBox *parametersGroupBOx;
    QHBoxLayout *horizontalLayout_2;
    QProgressBar *progressBar;
    QLabel *label_2;
    QDoubleSpinBox *aSpinBox;
    QLabel *label_3;
    QDoubleSpinBox *fiSphinBox;
    QPushButton *showButton;
    QMenuBar *menuBar;
    QMenu *menuFile;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;
    QToolBar *toolBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(910, 652);
        MainWindow->setWindowOpacity(1);
        actionOpen = new QAction(MainWindow);
        actionOpen->setObjectName(QString::fromUtf8("actionOpen"));
        actionSave = new QAction(MainWindow);
        actionSave->setObjectName(QString::fromUtf8("actionSave"));
        actionReset = new QAction(MainWindow);
        actionReset->setObjectName(QString::fromUtf8("actionReset"));
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        verticalLayout = new QVBoxLayout(centralWidget);
        verticalLayout->setSpacing(6);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        projectionViewWidget = new QTabWidget(centralWidget);
        projectionViewWidget->setObjectName(QString::fromUtf8("projectionViewWidget"));
        originalViewTab = new QWidget();
        originalViewTab->setObjectName(QString::fromUtf8("originalViewTab"));
        projectionViewWidget->addTab(originalViewTab, QString());
        doubleViewTab = new QWidget();
        doubleViewTab->setObjectName(QString::fromUtf8("doubleViewTab"));
        verticalLayout_2 = new QVBoxLayout(doubleViewTab);
        verticalLayout_2->setSpacing(6);
        verticalLayout_2->setContentsMargins(11, 11, 11, 11);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        logTextEdit = new QTextEdit(doubleViewTab);
        logTextEdit->setObjectName(QString::fromUtf8("logTextEdit"));
        logTextEdit->setReadOnly(true);

        verticalLayout_2->addWidget(logTextEdit);

        projectionViewWidget->addTab(doubleViewTab, QString());

        verticalLayout->addWidget(projectionViewWidget);

        parametersGroupBOx = new QGroupBox(centralWidget);
        parametersGroupBOx->setObjectName(QString::fromUtf8("parametersGroupBOx"));
        parametersGroupBOx->setSizeIncrement(QSize(0, 1));
        horizontalLayout_2 = new QHBoxLayout(parametersGroupBOx);
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        progressBar = new QProgressBar(parametersGroupBOx);
        progressBar->setObjectName(QString::fromUtf8("progressBar"));
        progressBar->setValue(0);

        horizontalLayout_2->addWidget(progressBar);

        label_2 = new QLabel(parametersGroupBOx);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        horizontalLayout_2->addWidget(label_2);

        aSpinBox = new QDoubleSpinBox(parametersGroupBOx);
        aSpinBox->setObjectName(QString::fromUtf8("aSpinBox"));
        aSpinBox->setMaximum(360);

        horizontalLayout_2->addWidget(aSpinBox);

        label_3 = new QLabel(parametersGroupBOx);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        horizontalLayout_2->addWidget(label_3);

        fiSphinBox = new QDoubleSpinBox(parametersGroupBOx);
        fiSphinBox->setObjectName(QString::fromUtf8("fiSphinBox"));
        fiSphinBox->setMaximum(360);

        horizontalLayout_2->addWidget(fiSphinBox);

        showButton = new QPushButton(parametersGroupBOx);
        showButton->setObjectName(QString::fromUtf8("showButton"));
        showButton->setFocusPolicy(Qt::ClickFocus);
        showButton->setAutoDefault(true);
        showButton->setDefault(true);
        showButton->setFlat(false);

        horizontalLayout_2->addWidget(showButton);


        verticalLayout->addWidget(parametersGroupBOx);

        MainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 910, 21));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        MainWindow->setMenuBar(menuBar);
        mainToolBar = new QToolBar(MainWindow);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        MainWindow->setStatusBar(statusBar);
        toolBar = new QToolBar(MainWindow);
        toolBar->setObjectName(QString::fromUtf8("toolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, toolBar);

        menuBar->addAction(menuFile->menuAction());
        menuFile->addAction(actionOpen);
        menuFile->addAction(actionSave);
        menuFile->addAction(actionReset);

        retranslateUi(MainWindow);

        projectionViewWidget->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "ProjectX", 0, QApplication::UnicodeUTF8));
        actionOpen->setText(QApplication::translate("MainWindow", "Open", 0, QApplication::UnicodeUTF8));
        actionOpen->setShortcut(QApplication::translate("MainWindow", "Ctrl+O", 0, QApplication::UnicodeUTF8));
        actionSave->setText(QApplication::translate("MainWindow", "Save", 0, QApplication::UnicodeUTF8));
        actionSave->setShortcut(QApplication::translate("MainWindow", "Ctrl+S", 0, QApplication::UnicodeUTF8));
        actionReset->setText(QApplication::translate("MainWindow", "Reset", 0, QApplication::UnicodeUTF8));
        projectionViewWidget->setTabText(projectionViewWidget->indexOf(originalViewTab), QApplication::translate("MainWindow", "Model viewss", 0, QApplication::UnicodeUTF8));
        projectionViewWidget->setTabText(projectionViewWidget->indexOf(doubleViewTab), QApplication::translate("MainWindow", "Log", 0, QApplication::UnicodeUTF8));
        parametersGroupBOx->setTitle(QApplication::translate("MainWindow", "Parameters", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("MainWindow", "<p>&#920;</p>\n"
"\n"
"", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("MainWindow", "<p>&#946;</p>", 0, QApplication::UnicodeUTF8));
        showButton->setText(QApplication::translate("MainWindow", "Show", 0, QApplication::UnicodeUTF8));
        menuFile->setTitle(QApplication::translate("MainWindow", "File", 0, QApplication::UnicodeUTF8));
        toolBar->setWindowTitle(QApplication::translate("MainWindow", "toolBar", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
