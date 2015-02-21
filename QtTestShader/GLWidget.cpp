
#include "stdafx.h"
#include "GLWidget.h"
#include "MyPoint.h"

GLWidget::GLWidget(QGLFormat format, QWidget *parent) :
	QGLWidget(format, parent),
	cRenderer(0),
	pressIdx(-1)
{
	cRenderer = new CurveRenderer();

	/*
	// Case loop artifact #1
	points.push_back(CVSystem::MyPoint(108, 530));
	points.push_back(CVSystem::MyPoint(117, 74));
	points.push_back(CVSystem::MyPoint(639, 227));
	points.push_back(CVSystem::MyPoint(896, 368));
	*/
	
	// Case loop artifact #2
	points.push_back(CVSystem::MyPoint(154, 670));
	points.push_back(CVSystem::MyPoint(368, 511));
	points.push_back(CVSystem::MyPoint(722, 354));
	points.push_back(CVSystem::MyPoint(932, 551));
	
}

GLWidget::~GLWidget()
{
	if(cRenderer) delete cRenderer;
}

void GLWidget::initializeGL() 
{	
	glClearColor(0.75, 0.75, 0.75, 1.0); 

	//shaderProgram.addShaderFromSourceFile(QGLShader::Vertex,	"VertexShader.vsh");
	shaderProgram.addShaderFromSourceFile(QGLShader::Fragment,	"Cubic.fsh");
	shaderProgram.link();

	cRenderer->shaderProgram = &shaderProgram;
}

void GLWidget::resizeGL( int width, int height )
{
}

bool GLWidget::event( QEvent * event ) { return QGLWidget::event(event); }

/* void GLWidget::resizeGL(int width, int height)  { QGLWidget::resizeGL(width, height); }*/

void GLWidget::paintGL() 
{
	glViewport(0, 0, this->width(),  this->height());	

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluOrtho2D(0, this->width(),  this->height(),  0);	// flip the y axis

	glMatrixMode(GL_MODELVIEW);	
	glLoadIdentity();
	glScalef(1, 1, 1);

	glClear(GL_COLOR_BUFFER_BIT);

	/*GLfloat m1[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, m1);
	QMatrix4x4 mv_mat( m1[0],m1[1],m1[2],m1[3],
					   m1[4],m1[5],m1[6],m1[7],
					   m1[8],m1[9],m1[10],m1[11],
					   m1[12],m1[13],m1[14],m1[15]);

	GLfloat m2[16];
	glGetFloatv(GL_PROJECTION_MATRIX, m2);
	QMatrix4x4 p_mat( m2[0], m2[1], m2[2], m2[3],
					  m2[4], m2[5], m2[6], m2[7],
					  m2[8], m2[9], m2[10],m2[11],
					  m2[12],m2[13],m2[14],m2[15]);*/
	
	shaderProgram.bind();

	//shaderProgram.setUniformValue("mat", mv_mat * p_mat);
	//shaderProgram.setUniformValue("color", QColor(Qt::red

	/*
	// Demo for quadratic curve
	glBegin(GL_TRIANGLES);

	glTexCoord2f(0.0f, 0.0f);
	glVertex3f(15,  15, 0);

	glTexCoord2f(0.5f, 0.0f);
	glVertex3f(25,  90, 0);

	glTexCoord2f(1.0f, 1.0f);
	glVertex3f(90,  25, 0);

	glEnd();*/

	cRenderer->ComputeCubic( points[0].x, points[0].y, 
						  points[1].x, points[1].y, 
						  points[2].x, points[2].y, 
						  points[3].x, points[3].y);

	shaderProgram.release();

	if(cRenderer->drawAdditionalTri)
		cRenderer->DrawPlainTriangle(cRenderer->atri[0].x, cRenderer->atri[0].y, cRenderer->atri[1].x, cRenderer->atri[1].y, cRenderer->atri[2].x, cRenderer->atri[2].y);

	/*
	DrawLine(cRenderer->t1[0].x, cRenderer->t1[0].y, cRenderer->t2[0].x, cRenderer->t2[0].y);
	DrawLine(cRenderer->t2[0].x, cRenderer->t2[0].y, cRenderer->t3[0].x, cRenderer->t3[0].y);
	DrawLine(cRenderer->t3[0].x, cRenderer->t3[0].y, cRenderer->t1[0].x, cRenderer->t1[0].y);

	DrawLine(cRenderer->t1[1].x, cRenderer->t1[1].y, cRenderer->t2[1].x, cRenderer->t2[1].y);
	DrawLine(cRenderer->t2[1].x, cRenderer->t2[1].y, cRenderer->t3[1].x, cRenderer->t3[1].y);
	DrawLine(cRenderer->t3[1].x, cRenderer->t3[1].y, cRenderer->t1[1].x, cRenderer->t1[1].y);

	DrawLine(cRenderer->t1[2].x, cRenderer->t1[2].y, cRenderer->t2[2].x, cRenderer->t2[2].y);
	DrawLine(cRenderer->t2[2].x, cRenderer->t2[2].y, cRenderer->t3[2].x, cRenderer->t3[2].y);
	DrawLine(cRenderer->t3[2].x, cRenderer->t3[2].y, cRenderer->t1[2].x, cRenderer->t1[2].y);*/

	DrawLine(points[0].x, points[0].y, points[1].x, points[1].y);
	DrawLine(points[2].x, points[2].y, points[3].x, points[3].y);

	DrawCircle(points[0], 15, 0, 1.0,   0.0, 1.0);
	DrawCircle(points[1], 15, 0, 1.0, 0,   1.0);
	DrawCircle(points[2], 15, 0, 1.0, 0,   1.0);
	DrawCircle(points[3], 15, 0, 1.0,   0.0, 1.0);

	
}

void GLWidget::DrawCircle(CVSystem::MyPoint pt, float radius, double r, double g, double b, double a)
{
	glDepthMask(GL_FALSE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	float delta_theta = 0.01;

	glColor4f(r, g, b, a);
	glBegin( GL_POLYGON );

	for( float angle = 0; angle < 2* M_PI; angle += delta_theta )
		glVertex3f( pt.x + radius * cos(angle), pt.y + radius * sin(angle), 0 );

	glEnd();
	glDisable(GL_BLEND);
	glDepthMask(GL_TRUE);
}

void GLWidget::DrawLine(float x1, float y1, float x2, float y2)
{
	glDepthMask(GL_FALSE);	
	glLineWidth(3.0); 
	glColor4f(1, 0, 0, 1.0);

	glBegin(GL_LINES);

	glVertex3f(x1, y1, 0);	
	glVertex3f(x2, y2, 0);

	glEnd();
	glDepthMask(GL_TRUE);  
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
	CVSystem::MyPoint cP(event->x(), event->y());
	for(int a = 0; a < points.size(); a++)
	{
		float d = points[a].Distance(cP);
		if(d < 10.0f)
		{
			pressIdx = a;
			break;
		}
	}

	if(pressIdx >= 0)
	{		
		points[pressIdx].x = event->x();
		points[pressIdx].y = event->y();
		this->updateGL();
	}
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
	if(pressIdx >= 0)
	{
		points[pressIdx].x = event->x();
		points[pressIdx].y =  event->y();
		this->updateGL();
	}
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
	if(pressIdx >= 0)
	{
		points[pressIdx].x = event->x();
		points[pressIdx].y = event->y();
		this->updateGL();
	}
	pressIdx = -1;
}
