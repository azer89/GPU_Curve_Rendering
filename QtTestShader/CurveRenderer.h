
#ifndef __VECTORRENDERER_H__
#define __VECTORRENDERER_H__

#include "stdafx.h"

#include "MyPoint.h"
#include "MyVector3.h"
#include "MyVertex.h"

#include <iostream>

#include <Eigen/Dense>

typedef enum _CurveType
{
	CURVE_TYPE_UNKNOWN		= 0,
	CURVE_TYPE_SERPENTINE	= 1,	
	CURVE_TYPE_LOOP			= 2,	
	CURVE_TYPE_CUSP			= 3,
	CURVE_TYPE_QUADRATIC    = 4,
	CURVE_TYPE_LINE			= 5,

} CurveType;

class CurveRenderer
{
public:
	CurveRenderer() {
		drawAdditionalTri = false;
		atri = std::vector<CVSystem::MyPoint>(3);
	}
	~CurveRenderer(){}

	QGLShaderProgram* shaderProgram;
	int curveTypeDebug;
		
	void ComputeCubic (float x0, float y0,
					float x1, float y1,
					float x2, float y2,
					float x3, float y3,
					int recursiveType = -1);

	CurveType DetermineType (CVSystem::MyPoint    v0, CVSystem::MyPoint    v1,
							 CVSystem::MyPoint    v2, CVSystem::MyPoint    v3,
							 double& d0, double& d1,
							 double& d2, double& d3);

	void Triangulation (double  x0, double  y0,
						   double  x1, double  y1,
						   double  x2, double  y2,
						   double  x3, double  y3,
						   std::vector<double> klm);

	void DrawTriangle (MyVertex v0, MyVertex v1, MyVertex v2);

	void DrawPlainTriangle (double x0, double y0,
							double x1, double y1,
							double x2, double y2);

	void DrawPlainTriangle (double x0, double y0,
							double x1, double y1,
							double x2, double y2,
							bool isInside);

	bool ApproxEqual(CVSystem::MyPoint& v0, CVSystem::MyPoint& v1);
	bool PointInTriangle(CVSystem::MyPoint& point, CVSystem::MyPoint& a, CVSystem::MyPoint& b, CVSystem::MyPoint& c);
	int Orientation(CVSystem::MyPoint& p1, CVSystem::MyPoint& p2, CVSystem::MyPoint& p3);
	bool LinesIntersect(CVSystem::MyPoint& p1, CVSystem::MyPoint& q1, CVSystem::MyPoint& p2, CVSystem::MyPoint& q2);

public:
	bool drawAdditionalTri;
	std::vector<CVSystem::MyPoint> atri;
};

#endif