#-------------------------------------------------
#
# Project created by QtCreator 2013-02-21T09:11:56
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Hawk
TEMPLATE = app

QMAKE_CFLAGS_DEBUG += -Itriangles

SOURCES += main.cpp\
        mainwindow.cpp \
        modelmaker.cpp \
        triangles/algs.cpp \
        modelviewer.cpp \
        triangles/exception.cpp \
        edgeselector.cpp \
        processthread.cpp \
        bsp/BSPAlg.cpp \
        Math/TriangleShared.cpp \
        calculatedialog.cpp

HEADERS  += mainwindow.h \
            modelmaker.h \
            triangles/algs.h \
            modelviewer.h \
            triangles/exception.h \
            edgeselector.h \
            processthread.h \
            bsp/BSPAlg.h \
            Math/Vector3.hpp \
            Math/Triangle.hpp \
            Math/Point3.hpp \
            Math/Complex.hpp \
            Math/MathDefines.h \
            Math/LineShared.h \            
            Calc/ComputeFieldDll.hpp \
            Calc/IComputeField.hpp \
            calculatedialog.h


FORMS    += mainwindow.ui \
    calculatedialog.ui

LIBS += -L$$PWD/libfortran/ -lRCS_functions_cyl
INCLUDEPATH += $$PWD/libfortran
DEPENDPATH += $$PWD/libfortran
PRE_TARGETDEPS += $$PWD/libfortran/RCS_functions_cyl.lib
