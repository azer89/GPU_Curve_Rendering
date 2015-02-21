#ifndef QTTESTSHADER_H
#define QTTESTSHADER_H

#include <QtWidgets/QMainWindow>
#include "ui_QtTestShader.h"
#include "GLWidget.h"

class QtTestShader : public QMainWindow
{
	Q_OBJECT

public:
	QtTestShader(QWidget *parent = 0);
	~QtTestShader();

private:
	Ui::QtTestShaderClass ui;

	QGridLayout *myGridLayout;
	GLWidget *myGLWidget;

private slots:
	void UpdateTimer();
};

#endif // QTTESTSHADER_H
