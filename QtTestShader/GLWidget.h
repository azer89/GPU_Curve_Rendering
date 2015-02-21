
#ifndef __GLWIDGET_H__
#define __GLWIDGET_H__

#include "stdafx.h"
#include "CurveRenderer.h"
#include "MyPoint.h"

//#include <QGLBuffer>


class GLWidget : public QGLWidget
{
	Q_OBJECT

public:
	GLWidget(QGLFormat format, QWidget *parent = 0);
	~GLWidget();

	void DrawCircle(CVSystem::MyPoint pt, float radius, double r, double g, double b, double a);
	void DrawLine(float x1, float y1, float x2, float y2);

	int GetCurveType()
	{
		return cRenderer->curveTypeDebug;
	}

protected:
	bool event( QEvent * event );
	void initializeGL();	
	void paintGL();
	void resizeGL(int w, int h);

	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event);
	
private:

	CurveRenderer* cRenderer;
	std::vector<CVSystem::MyPoint> points;

	QMatrix4x4 pMatrix;
	QGLShaderProgram shaderProgram;

	int pressIdx;

private:
};

#endif