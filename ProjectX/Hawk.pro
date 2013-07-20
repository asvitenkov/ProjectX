#-------------------------------------------------
#
# Project created by QtCreator 2013-02-21T09:11:56
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Hawk
TEMPLATE = app

QMAKE_CFLAGS_DEBUG += -I iangles

SOURCES += main.cpp\
        mainwindow.cpp \
        modelmaker.cpp \
        triangles/algs.cpp \
        modelviewer.cpp \
        triangles/exception.cpp \
        edgeselector.cpp \
        processthread.cpp \
        reditrectthread.cpp \
        cleanthread.cpp \
        bsp/BSPAlg.cpp \
        triangles/blockanalizator.cpp \
        Math/TriangleShared.cpp

HEADERS  += mainwindow.h \
            modelmaker.h \
            triangles/s ucts.h \
            triangles/algs.h \
            modelviewer.h \
            triangles/exception.h \
            edgeselector.h \
            processthread.h \
            bsp/BSPAlg.h \
            triangles/blockanalizator.h \
            Math/Vector3.hpp \
            Math/Triangle.hpp \
            Math/Point3.hpp \
            Math/Complex.hpp \
            Math/MathDefines.h \
            Math/TriangleShared.h

FORMS    += mainwindow.ui
