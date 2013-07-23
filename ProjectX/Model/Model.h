#ifndef MODEL_H
#define MODEL_H

#include <QWidget>
#include <QVector>
#include <QString>
#include "TriangleShared.h"
#include "ModelScene.h"
#include "ModelView.h"

class Model : public QObject
{
    Q_OBJECT

    friend class TriangleShared;
    friend class ModelScene;
public:
    explicit Model(QObject* widget = 0);

    void Create(const QString& filepath);
    QVector<TriangleShared>& GetTriangles();
    const QVector<TPoint3>& GetPoints() const;

    QVector<TriangleShared>& GetVisibleTringle();

    void SetTriangles(const QVector<TriangleShared>& newTr);
    void SetPoints(const QVector<TPoint3>& newP);
protected:
    void AddPoint(const TPoint3& p);

    QVector<TriangleShared> m_triangles;
    QVector<TPoint3>        m_points;
};

#endif // MODEL_H
