#-------------------------------------------------
#
# Project created by QtCreator 2013-02-21T09:11:56
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = ProjectX
TEMPLATE = app

QMAKE_CFLAGS_DEBUG += -Itriangles

SOURCES += main.cpp\
        mainwindow.cpp \
        Algo/exception.cpp \
        edgeselector.cpp \
        processthread.cpp \
        bsp/BSPAlg.cpp \
        Model/TriangleShared.cpp \
        Model/Model.cpp \
        Model/ModelView.cpp \
        Model/ModelScene.cpp \
        Algo/Algoritm.cpp \
        calculatedialog.cpp

HEADERS  += mainwindow.h \
            Algo/exception.h \
            edgeselector.h \
            processthread.h \
            bsp/BSPAlg.h \
            Math/Vector3.hpp \
            Math/Triangle.hpp \
            Math/Point3.hpp \
            Math/Complex.hpp \
            Math/MathDefines.h \
            Math/LineShared.h \
            Model/TriangleShared.h \
            Model/Model.h \
            Model/ModelView.h \
            Model/ModelScene.h \
            Algo/Algoritm.h \
            Calc/ComputeFieldDll.hpp \
            Calc/IComputeField.hpp \
            calculatedialog.h


FORMS    += mainwindow.ui \
    calculatedialog.ui

LIBS += -L$$PWD/libfortran/ -lRCS_functions_cyl
INCLUDEPATH += $$PWD/libfortran
DEPENDPATH += $$PWD/libfortran
PRE_TARGETDEPS += $$PWD/libfortran/RCS_functions_cyl.lib
