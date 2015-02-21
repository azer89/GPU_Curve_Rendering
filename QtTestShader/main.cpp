#include "stdafx.h"
#include "QtTestShader.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	QtTestShader w;
	w.show();
	return a.exec();
}
