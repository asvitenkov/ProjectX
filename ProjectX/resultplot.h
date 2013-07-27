#ifndef RESULTPLOT_H
#define RESULTPLOT_H

#include <QWidget>
#include <QTableWidget>

class ResultPlot : public QWidget
{
    Q_OBJECT
public:
    explicit ResultPlot(double startRVal, double startCVal, double rowCount, double colCount, double rowStep, double colStep,  QWidget *parent = 0);
    bool addItem(double rVal, double cVal, double value);
    bool addItem(int row, int col, double value);

signals:
    
public slots:
    
private:
    void init();
    int valueToRow(double value);
    int valueToCol(double value);

    double mStartRVal, mStartCVal;
    double mRowStep, mColStep;
    double mMaxVal, mMinVal;

    int mRowCount;
    int mColCount;

    QTableWidget *mTableWidget;

};

#endif // RESULTPLOT_H
