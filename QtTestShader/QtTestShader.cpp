#include "stdafx.h"
#include "QtTestShader.h"

QtTestShader::QtTestShader(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	
	QGLFormat glFormat;
	glFormat.setVersion(2, 0);
	glFormat.setProfile(QGLFormat::CompatibilityProfile);
	//glFormat.setProfile( QGLFormat::CoreProfile );
	glFormat.setSampleBuffers( true );
	
	myGridLayout = new QGridLayout(ui.centralWidget);
	myGLWidget = new GLWidget(glFormat, ui.centralWidget);

	myGridLayout->setSpacing(6);
	myGridLayout->setContentsMargins(11, 11, 11, 11);
	myGridLayout->setObjectName(QStringLiteral("gridLayout"));	
	myGLWidget->setObjectName(QStringLiteral("widget"));

	myGridLayout->addWidget(myGLWidget, 0, 0, 1, 1);

	// timer
	QTimer *timer = new QTimer(this);
	connect(timer, SIGNAL(timeout()), this, SLOT(UpdateTimer()));
	timer->start(100);
}

QtTestShader::~QtTestShader()
{
}

void QtTestShader::UpdateTimer()
{
	//std::cout << "update\n";
	int cType = myGLWidget->GetCurveType();

	/*
	CURVE_TYPE_UNKNOWN = 0,
	CURVE_TYPE_SERPENTINE,	
	CURVE_TYPE_LOOP,	
	CURVE_TYPE_CUSP_WITH_INFLECTION_AT_INFINITY,
	CURVE_TYPE_CUPS_WITH_CUPS_AT_INFINITY,
	CURVE_TYPE_QUADRATIC,
	CURVE_TYPE_LINE,
	*/

	QString str = "CURVE_TYPE_UNKNOWN";

	if(cType == 1)		str = "CURVE_TYPE_SERPENTINE";
	else if(cType == 2) str = "CURVE_TYPE_LOOP";
	else if(cType == 3) str = "CURVE_TYPE_CUSP";
	else if(cType == 4) str = "CURVE_TYPE_QUADRATIC";
	else if(cType == 5) str = "CURVE_TYPE_LINE";

	ui.statusBar->showMessage(str, 100);
}
