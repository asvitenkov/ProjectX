#include "resultplot.h"

#include <QGridLayout>

ResultPlot::ResultPlot(double startRVal, double startCVal, double rowCount, double colCount, double rowStep, double colStep, QWidget *parent):
    QWidget(parent), mStartRVal(startRVal), mStartCVal(startCVal), mColCount(colCount), mRowCount(rowCount), mRowStep(rowStep), mColStep(colStep), mTableWidget(NULL)
{
    //mTableWidget.setParent(this, Qt::Widget);
    init();
}


void ResultPlot::init()
{

    if(mTableWidget != NULL)
    {
        delete mTableWidget;
        mTableWidget = NULL;
    }

    if(mColCount <= 0 || mRowCount <= 0)
        return;

    mTableWidget = new QTableWidget(mColCount, mRowCount);

    double tmpVal;
    for(int i=0; i< mRowCount; i++)
    {
        tmpVal = i * mRowStep + mStartRVal;
        mTableWidget->setVerticalHeaderItem(i, new QTableWidgetItem(QString::number(tmpVal)));
    }

    for(int i=0; i< mColCount; i++)
    {
        tmpVal = i * mColStep + mStartCVal;
        mTableWidget->setHorizontalHeaderItem(i, new QTableWidgetItem(QString::number(tmpVal)));
    }

    //mTableWidget->setHorizontalHeaderItem();

    QGridLayout *pLayout = new QGridLayout(this);
    pLayout->addWidget(mTableWidget);
    this->setLayout(pLayout);
}
bool  ResultPlot::addItem(double rVal, double cVal, double value)
{

}


bool ResultPlot::addItem(int row, int col, double value)
{

}

int ResultPlot::valueToCol(double value)
{

}

int ResultPlot::valueToRow(double value)
{

}
