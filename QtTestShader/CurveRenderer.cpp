
#include "stdafx.h"
#include "CurveRenderer.h"

void CurveRenderer::ComputeCubic(float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3, int recursiveType)
{	
	double d0 = 0; double d1 = 0; double d2 = 0; double d3 = 0;
	CurveType curve_type = DetermineType (CVSystem::MyPoint(x0, y0),
										  CVSystem::MyPoint(x1, y1),
										  CVSystem::MyPoint(x2, y2),
										  CVSystem::MyPoint(x3, y3),
										  d0, d1, d2, d3);

	// debug
	curveTypeDebug = (int)curve_type;

	double OneThird = 1.0f / 3.0f;
	double TwoThirds = 2.0f / 3.0f;

	float t1;
	float ls;
	float lt;
	float ms;
	float mt;

	float ltMinusLs;
	float mtMinusMs;
	float lsMinusLt;
	float ql;
	float qm;

	// bool flip
	bool flip = false;

	// artifact on loop
	int errorLoop = -1;
	double splitParam = 0;

	std::vector<double> klm(12);
		

	switch (curve_type)
	{
	case CURVE_TYPE_UNKNOWN:
		break;

	case CURVE_TYPE_SERPENTINE:	

		t1 = sqrtf(9.0f * d2 * d2 - 12 * d1 * d3);
		ls = 3.0f * d2 - t1;
		lt = 6.0f * d1;
		ms = 3.0f * d2 + t1;
		mt = lt;
		ltMinusLs = lt - ls;
		mtMinusMs = mt - ms;

		klm[0] =  ls * ms;
		klm[1] =  	ls * ls * ls;
		klm[2] =  	ms * ms * ms;

		klm[3] =  OneThird * (3.0f * ls * ms - ls * mt - lt * ms);
		klm[4]  = 	ls * ls * (ls - lt);
		klm[5]  = 	ms * ms * (ms - mt);

		klm[6] =  OneThird * (lt * (mt - 2.0f * ms) + ls * (3.0f * ms - 2.0f * mt));
		klm[7] =  	ltMinusLs * ltMinusLs * ls;
		klm[8] =  	mtMinusMs * mtMinusMs * ms;

		klm[9] =  ltMinusLs * mtMinusMs;
		klm[10] =  	-(ltMinusLs * ltMinusLs * ltMinusLs);
		klm[11]  = 	-(mtMinusMs * mtMinusMs * mtMinusMs);

		if (d1 < 0.0f)
			flip = true;

		break;

	case CURVE_TYPE_LOOP:

		t1 = sqrtf(4.0f * d1 * d3 - 3.0f * d2 * d2);
	    ls = d2 - t1;
		lt = 2.0f * d1;
		ms = d2 + t1;
		mt = lt;

		// Figure out whether there is a rendering artifact requiring
		// the curve to be subdivided by the caller.
		 ql = ls / lt;
		 qm = ms / mt;
		if (0.0f < ql && ql < 1.0f) 
		{
			errorLoop = 1;
			splitParam = ql;
			//std::cout << "error loop 1\n";
		}

		if (0.0f < qm && qm < 1.0f) 
		{
			errorLoop = 2;
			splitParam = qm;

			//std::cout << "error loop 2\n";
		}
				 
		ltMinusLs = lt - ls;
		mtMinusMs = mt - ms;
		klm[0] =  ls * ms;
		klm[1] =  	ls * ls * ms;
		klm[2] =  	ls * ms * ms;

		klm[3] =  OneThird * (-ls * mt - lt * ms + 3.0f * ls * ms);
		klm[4] =  	-OneThird * ls * (ls * (mt - 3.0f * ms) + 2.0f * lt * ms);
		klm[5] =  	-OneThird * ms * (ls * (2.0f * mt - 3.0f * ms) + lt * ms);

		klm[6] =  OneThird * (lt * (mt - 2.0f * ms) + ls * (3.0f * ms - 2.0f * mt));
		klm[7] =  	OneThird * (lt - ls) * (ls * (2.0f * mt - 3.0f * ms) + lt * ms);
		klm[8] =  	OneThird * (mt - ms) * (ls * (mt - 3.0f * ms) + 2.0f * lt * ms);

		klm[9] =  ltMinusLs * mtMinusMs;
		klm[10] =  	-(ltMinusLs * ltMinusLs) * mtMinusMs;
		klm[11] =  	-ltMinusLs * mtMinusMs * mtMinusMs;

		if(recursiveType == -1)
			flip =  ((d1 > 0.0f && klm[0] < 0.0f) || (d1 < 0.0f && klm[0] > 0.0f));

		break;

	case CURVE_TYPE_CUSP:
		ls = d3;
		lt = 3.0f * d2;
		lsMinusLt = ls - lt;
		klm[0] =  ls;
		klm[1] =  ls * ls * ls;
		klm[2] =  1.0f;

		klm[3] =  ls - OneThird * lt;
		klm[4] =  	ls * ls * lsMinusLt;
		klm[5] =  	1.0f;

		klm[6] =  ls - TwoThirds * lt;
		klm[7] =  	lsMinusLt * lsMinusLt * ls;
		klm[8] =  	1.0f;

		klm[9] =  lsMinusLt;
		klm[10] =  	lsMinusLt * lsMinusLt * lsMinusLt;
		klm[11] =  	1.0f;

		break;

	case CURVE_TYPE_QUADRATIC:
		klm[0] =  0;
		klm[1] =  0;
		klm[2] =  0;

		klm[3] =  OneThird;
		klm[4] =  0;
		klm[5] =  OneThird;

		klm[6] =  TwoThirds;
		klm[7] =  OneThird;
		klm[8] =  TwoThirds;

		klm[9] =  1;
		klm[10] = 1;
		klm[11] = 1;

		if (d3 < 0)
			flip = true;

		break;

	case CURVE_TYPE_LINE:
		break;
	}

	if(errorLoop != -1 && recursiveType == -1)
	{
		double x01 = (x1 - x0) * splitParam + x0;		double x12 = (x2 - x1) * splitParam + x1;		double x23 = (x3 - x2) * splitParam + x2;
		double y01 = (y1 - y0) * splitParam + y0;		double y12 = (y2 - y1) * splitParam + y1;		double y23 = (y3 - y2) * splitParam + y2;		

		double x012 = (x12 - x01) * splitParam + x01;	double x123 = (x23 - x12) * splitParam + x12;
		double y012 = (y12 - y01) * splitParam + y01;	double y123 = (y23 - y12) * splitParam + y12;

		double x0123 = (x123 - x012) * splitParam + x012;
		double y0123 = (y123 - y012) * splitParam + y012;

		

		drawAdditionalTri = true;
		atri[0] = CVSystem::MyPoint(x0, y0);
		atri[1] = CVSystem::MyPoint(x0123, y0123);
		atri[2] = CVSystem::MyPoint(x3, y3);
		//DrawPlainTriangle(atri[0].x, atri[0].y, atri[1].x, atri[1].y, atri[2].x, atri[2].y);


		if(errorLoop == 1)	// flip second
		{
			ComputeCubic(x0, y0, x01, y01, x012, y012, x0123, y0123, 0);
			ComputeCubic(x0123,  y0123, x123, y123, x23, y23, x3, y3, 1);
		}
		else if(errorLoop == 2) // flip first
		{
			ComputeCubic(x0, y0, x01, y01, x012, y012, x0123, y0123, 1);
			ComputeCubic(x0123,  y0123, x123, y123, x23, y23, x3, y3, 0);
		}

		/*if(flip) 
		{
			DrawPlainTriangle(x0, y0,  x0123, y0123, x3, y3, false);
		}
		else
		{
			DrawPlainTriangle(x0, y0,  x0123, y0123, x3, y3, true);
		}*/

		return;
	}
	else if(errorLoop == -1 && recursiveType == -1)
	{
		drawAdditionalTri = false;
	}

	if(recursiveType == 1)
		flip = !flip;

	if(flip)
	{
		klm[0] = -klm[0]; klm[1] = -klm[1];
		klm[3] = -klm[3]; klm[4] = -klm[4];
		klm[6] = -klm[6]; klm[7] = -klm[7];
		klm[9] = -klm[9]; klm[10] = -klm[10];
	}

	Triangulation (x0, y0, x1, y1, x2, y2, x3, y3, klm);

	//if(recursiveType != -1)
	//{
	//	DrawPlainTriangle(atri[0].x, atri[0].y, atri[1].x, atri[1].y, atri[2].x, atri[2].y);
	//}

	//if(drawAdditionalTri)
	//{
	//	std::cout << "error loop\n";
	//	DrawPlainTriangle(atri[0].x, atri[0].y, atri[1].x, atri[1].y, atri[2].x, atri[2].y, false);
	//}


}

void CurveRenderer::Triangulation (double  x0, double  y0,
									  double  x1, double  y1,
									  double  x2, double  y2,
									  double  x3, double  y3,
									  std::vector<double> klm)
{
	std::vector<MyVertex> vertices;
	vertices.push_back(MyVertex(x0, y0, klm[0], klm[1], klm[2]));
	vertices.push_back(MyVertex(x1, y1, klm[3], klm[4], klm[5]));
	vertices.push_back(MyVertex(x2, y2, klm[6], klm[7], klm[8]));
	vertices.push_back(MyVertex(x3, y3, klm[9], klm[10], klm[11]));

	// First test for degenerate cases.
	for (int i = 0; i < 4; ++i) 
	{
		for (int j = i + 1; j < 4; ++j) 
		{
			if (ApproxEqual(vertices[i].xyCoor, vertices[j].xyCoor)) 
			{
				// Two of the vertices are coincident, so we can eliminate at
				// least one triangle. We might be able to eliminate the other
				// as well, but this seems sufficient to avoid degenerate
				// triangulations.
				int indices[3] = { 0 };
				int index = 0;
				for (int k = 0; k < 4; ++k)
					if (k != j)
						indices[index++] = k;

				DrawTriangle(vertices[indices[0]],
					vertices[indices[1]],
					vertices[indices[2]]);

				return;
			}
		}
	}

	// See whether any of the points are fully contained in the
	// triangle defined by the other three.
	for (int i = 0; i < 4; ++i) 
	{
		int indices[3] = { 0 };
		int index = 0;
		for (int j = 0; j < 4; ++j)
			if (i != j)
				indices[index++] = j;

		if (PointInTriangle(vertices[i].xyCoor, vertices[indices[0]].xyCoor, vertices[indices[1]].xyCoor, vertices[indices[2]].xyCoor)) 
		{
			// Produce three triangles surrounding this interior vertex.
			for (int j = 0; j < 3; ++j)
				DrawTriangle(vertices[indices[j % 3]], vertices[indices[(j + 1) % 3]], vertices[i]);

			// Mark the interior vertex so we ignore it if trying to trace
			// the interior edge ?
			//m_vertices[i].setInterior(true);
			return;
		}
	}

	// There are only a few permutations of the vertices, ignoring
	// rotations, which are irrelevant:
	//
	//  0--3  0--2  0--3  0--1  0--2  0--1
	//  |  |  |  |  |  |  |  |  |  |  |  |
	//  |  |  |  |  |  |  |  |  |  |  |  |
	//  1--2  1--3  2--1  2--3  3--1  3--2
	//
	// Note that three of these are reflections of each other.
	// Therefore there are only three possible triangulations:
	//
	//  0--3  0--2  0--3
	//  |\ |  |\ |  |\ |
	//  | \|  | \|  | \|
	//  1--2  1--3  2--1
	//
	// From which we can choose by seeing which of the potential
	// diagonals intersect. Note that we choose the shortest diagonal
	// to split the quad.
	if (LinesIntersect(vertices[0].xyCoor, vertices[2].xyCoor, vertices[1].xyCoor, vertices[3].xyCoor)) 
	{
		if ((vertices[2].xyCoor - vertices[0].xyCoor).DiagonalLengthSquared() < (vertices[3].xyCoor - vertices[1].xyCoor).DiagonalLengthSquared()) 
		{
			DrawTriangle(vertices[0], vertices[1], vertices[2]);
			DrawTriangle(vertices[0], vertices[2], vertices[3]);
		}
		else 
		{
			DrawTriangle(vertices[0], vertices[1], vertices[3]);
			DrawTriangle(vertices[1], vertices[2], vertices[3]);
		}
	} 
	else if (LinesIntersect(vertices[0].xyCoor, vertices[3].xyCoor, vertices[1].xyCoor, vertices[2].xyCoor)) 
	{
		if ((vertices[3].xyCoor - vertices[0].xyCoor).DiagonalLengthSquared() < (vertices[2].xyCoor - vertices[1].xyCoor).DiagonalLengthSquared()) 
		{
			DrawTriangle(vertices[0], vertices[1], vertices[3]);
			DrawTriangle(vertices[0], vertices[3], vertices[2]);
		} 
		else 
		{
			DrawTriangle(vertices[0], vertices[1], vertices[2]);
			DrawTriangle(vertices[2], vertices[1], vertices[3]);
		}
	} 
	else 
	{
		// Lines (0->1), (2->3) intersect -- or should, modulo numerical
		// precision issues
		if ((vertices[1].xyCoor - vertices[0].xyCoor).DiagonalLengthSquared() < (vertices[3].xyCoor - vertices[2].xyCoor).DiagonalLengthSquared())
		{
			DrawTriangle(vertices[0],vertices[2], vertices[1]);
			DrawTriangle(vertices[0], vertices[1], vertices[3]);
		} 
		else 
		{
			DrawTriangle(vertices[0],vertices[2], vertices[3]);
			DrawTriangle(vertices[3], vertices[2],vertices[1]);
		}
	}
}

void CurveRenderer::DrawPlainTriangle (double x0, double y0,
						double x1, double y1,
						double x2, double y2)
{
	//glDepthMask(GL_FALSE);
	//glEnable(GL_BLEND);
	//glBlendFunc(GL_DST_COLOR, GL_ONE_MINUS_SRC_ALPHA);

	glColor4f(0.0, 0.0, 0.0, 1.0);
	glBegin (GL_TRIANGLES);
	glVertex3f (x0, y0, 0.0f);
	glVertex3f (x1, y1, 0.0f);
	glVertex3f (x2, y2, 0.0f);
	glEnd ();

	//glDisable(GL_BLEND);
	//glDepthMask(GL_TRUE);
}

void CurveRenderer::DrawPlainTriangle (double x0, double y0, double x1, double y1, double x2, double y2, bool isInside)
{
	if(isInside)
	{
		glColor4f(0.0, 0.0, 0.0, 1.0);
	}
	else
	{
		glColor4f(1.0, 1.0, 0.0, 1.0);
	}

	glBegin (GL_TRIANGLES);
	glVertex3f (x0, y0, 0.0f);
	glVertex3f (x1, y1, 0.0f);
	glVertex3f (x2, y2, 0.0f);
	glEnd ();
}

void CurveRenderer::DrawTriangle (MyVertex v0,
								  MyVertex v1,
								  MyVertex v2)
{
	//std::cout << "draw\n";

	shaderProgram->setUniformValue("insideColor",  QColor(Qt::black));
	shaderProgram->setUniformValue("outsideColor", QColor(Qt::white));

	glDepthMask(GL_FALSE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_DST_COLOR, GL_ONE_MINUS_SRC_ALPHA);

	

	glBegin (GL_TRIANGLES);
	glTexCoord3f (v0.k, v0.l, v0.m);
	glVertex3f (v0.xyCoor.x, v0.xyCoor.y, 0.0f);

	glTexCoord3f (v1.k, v1.l, v1.m);
	glVertex3f (v1.xyCoor.x, v1.xyCoor.y, 0.0f);

	glTexCoord3f (v2.k, v2.l, v2.m);
	glVertex3f (v2.xyCoor.x, v2.xyCoor.y, 0.0f);
	glEnd ();


	glDisable(GL_BLEND);
	glDepthMask(GL_TRUE);
}


CurveType CurveRenderer::DetermineType (CVSystem::MyPoint    v0,
										CVSystem::MyPoint    v1,
										CVSystem::MyPoint    v2,
										CVSystem::MyPoint    v3,
										double& d0,
										double& d1,
										double& d2,
										double& d3)
{
	d0 = 0;

	Eigen:: Vector3d b0(v0.x, v0.y, 1.0f);
	Eigen:: Vector3d b1(v1.x, v1.y, 1.0f);
	Eigen:: Vector3d b2(v2.x, v2.y, 1.0f);
	Eigen:: Vector3d b3(v3.x, v3.y, 1.0f);

	double a1 = b0.dot(b3.cross(b2));
	double a2 = b1.dot(b0.cross(b3));
	double a3 = b2.dot(b1.cross(b0));

	d1 = a1 - 2 * a2 + 3 * a3;
	d2 = -a2 + 3 * a3;
	d3 = 3 * a3;

	// normalize?
	/*double length = sqrt(d0 * d0 + d1 * d1 + d2* d2 + d3 * d3);
	d0 /= length;
	d1 /= length;
	d2 /= length;
	d3 /= length;*/

	double D = 3.0 * d2 * d2 - 4.0 * d1 * d3;
	double disc = d1 * d1 * D;

	/*if (disc > -M_EPS && disc < M_EPS)
		disc = 0;

	if (d0 > -M_EPS && d0 < M_EPS)
		d0 = 0;

	if (d1 > -M_EPS && d1 < M_EPS)
		d1 = 0;

	if (d2 > -M_EPS && d2 < M_EPS)
		d2 = 0;

	if (d3 > -M_EPS && d3 < M_EPS)
		d3 = 0;
	*/
	//curve_type = CURVE_TYPE_UNKNOWN;

	if (!disc) 
	{
		if (!d1 && !d2) 
		{
			if (!d3)
			{
				return CURVE_TYPE_LINE;
			}

			return CURVE_TYPE_QUADRATIC;
		}

		if (d1 != 0)
		{
			return CURVE_TYPE_CUSP;
		}

		if (D < 0)
			return CURVE_TYPE_LOOP;

		return CURVE_TYPE_SERPENTINE;
	}

	if (disc > 0) return CURVE_TYPE_SERPENTINE;

	// discriminant < 0
	return CURVE_TYPE_LOOP;
}

bool CurveRenderer::ApproxEqual(CVSystem::MyPoint& v0, CVSystem::MyPoint& v1)
{
	return (v0 - v1).DiagonalLengthSquared() < M_EPS * M_EPS;
}

bool CurveRenderer::PointInTriangle(CVSystem::MyPoint& point,
					 CVSystem::MyPoint& a,
					 CVSystem::MyPoint& b,
					 CVSystem::MyPoint& c)
{
	// Algorithm from http://www.blackpawn.com/texts/pointinpoly/default.html
	float x0 = c.x - a.x;
	float y0 = c.y - a.y;
	float x1 = b.x - a.x;
	float y1 = b.y - a.y;
	float x2 = point.x - a.x;
	float y2 = point.y - a.y;

	float dot00 = x0 * x0 + y0 * y0;
	float dot01 = x0 * x1 + y0 * y1;
	float dot02 = x0 * x2 + y0 * y2;
	float dot11 = x1 * x1 + y1 * y1;
	float dot12 = x1 * x2 + y1 * y2;
	float denominator = dot00 * dot11 - dot01 * dot01;
	if (!denominator)
		// Triangle is zero-area. Treat query point as not being inside.
			return false;
	// Compute
	float inverseDenominator = 1.0f / denominator;
	float u = (dot11 * dot02 - dot01 * dot12) * inverseDenominator;
	float v = (dot00 * dot12 - dot01 * dot02) * inverseDenominator;

	return (u > 0.0f) && (v > 0.0f) && (u + v < 1.0f);
}

// Utility functions local to this file.
int CurveRenderer::Orientation(CVSystem::MyPoint& p1,
				CVSystem::MyPoint& p2,
				CVSystem::MyPoint& p3)
{
	float crossProduct = (p2.y - p1.y) * (p3.x - p2.x) - (p3.y - p2.y) * (p2.x - p1.x);
	return (crossProduct < 0.0f) ? -1 : ((crossProduct > 0.0f) ? 1 : 0);
}

bool CurveRenderer::LinesIntersect(CVSystem::MyPoint& p1,
				   CVSystem::MyPoint& q1,
				   CVSystem::MyPoint& p2,
				   CVSystem::MyPoint& q2)
{
	return (Orientation(p1, q1, p2) != Orientation(p1, q1, q2)
		&& Orientation(p2, q2, p1) != Orientation(p2, q2, q1));
}

/*

void CurveRenderer::DrawCubic(float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3, bool flip, bool isRecursiveCall)
{
//std::vector<double> coord1(12);

double d0 = 0; double d1 = 0; double d2 = 0; double d3 = 0;
CurveType curve_type = DetermineType (CVSystem::MyPoint(x0, y0),
CVSystem::MyPoint(x1, y1),
CVSystem::MyPoint(x2, y2),
CVSystem::MyPoint(x3, y3),
d0, d1, d2, d3);

// debug
curveTypeDebug = (int)curve_type;

// serpentine
float tl;
float sl;	
float tm;
float sm;
//float ls, lt, ms, mt;

// loop
float td;
float sd;
float te;
float se;
double td_sd;
double te_se;

Eigen::MatrixXd invM3(4,4);
invM3 << 1, 0,         0,         0,
1, 1.0 / 3.0, 0,         0,
1, 2.0 / 3.0, 1.0 / 3.0, 0,
1, 1,         1,         1;

Eigen::MatrixXd F(4,4);

double splitLoop = -1.0;

switch (curve_type)
{
case CURVE_TYPE_UNKNOWN:
break;

case CURVE_TYPE_SERPENTINE:	

tl = d2 + ((1.0 / sqrt(3.0)) * sqrt(3.0 * d2 * d2 - 4.0 * d1 * d3));
sl = 2.0 * d1;				
tm = d2 - ((1.0 / sqrt(3.0)) * sqrt(3.0 * d2 * d2 - 4.0 * d1 * d3));
sm = 2.0 * d1;

F << tl * tm,					tl * tl * tl,			tm * tm * tm,			1,
-(sm * tl) -(sl * tm),		-(3.0 * sl * tl * tl),	-(3.0 * sm * tm * tm),	0,
sl * sm,					3.0 * sl * sl * tl,		3.0 * sm * sm *	 tm,	0,
0,							-(sl * sl * sl),		-(sm * sm * sm),		0;	




break;

case CURVE_TYPE_LOOP:
td = d2 + sqrt(4.0 * d1 * d3 - 3.0 * d2 *d2);
sd = 2.0 * d1;
te = d2 - sqrt(4.0 * d1 * d3 - 3.0 * d2 * d2);
se = 2.0 * d1;

// Get artifact, should subdivide the curve by getting t
td_sd = td / sd;
te_se = te / se;
if(td_sd > 0.0 && td_sd < 1.0) splitLoop = td / sd;
else if(te_se > 0.0 && te_se < 1.0) splitLoop = te / se;

F << td * te,					td * td * te,								td * te * te,							1,		
(-se * td) - (se * te),	(-se * td * td) - (2.0 * sd * te * td),		(-sd * te * te) - (2.0 * se * td * te),	0,
sd * se,					te * sd * sd + 2.0 * se * td* sd,			td * se * se + 2 * sd * te * se,		0,
0,							-sd * sd * se,								-sd * se * se,							0;

break;

case CURVE_TYPE_CUSP_WITH_INFLECTION_AT_INFINITY:
break;

case CURVE_TYPE_CUPS_WITH_CUPS_AT_INFINITY:
break;

case CURVE_TYPE_QUADRATIC:
break;

case CURVE_TYPE_LINE:
break;

}

if(splitLoop > 0.0 && splitLoop < 1.0 && !isRecursiveCall)
{
// SPLIT
double x01 = (x1 - x0) * splitLoop + x0;		double x12 = (x2 - x1) * splitLoop + x1;		double x23 = (x3 - x2) * splitLoop + x2;
double y01 = (y1 - y0) * splitLoop + y0;		double y12 = (y2 - y1) * splitLoop + y1;		double y23 = (y3 - y2) * splitLoop + y2;		

double x012 = (x12 - x01) * splitLoop + x01;	double x123 = (x23 - x12) * splitLoop + x12;
double y012 = (y12 - y01) * splitLoop + y01;	double y123 = (y23 - y12) * splitLoop + y12;

double x0123 = (x123 - x012) * splitLoop + x012;
double y0123 = (y123 - y012) * splitLoop + y012;

// Curve A
DrawCubic(x0, y0, x01, y01, x012, y012, x0123, y0123, false, true);

// Curve B
DrawCubic(x0123,  y0123, x123, y123, x23, y23, x3, y3, true, true);
}
else
{

// Draw as usual...
Eigen::MatrixXd rMat = invM3 * F;

bool serpFlip = false;
if(curve_type == CurveType::CURVE_TYPE_SERPENTINE)
{
std::cout << d0 << " " <<  d1 << " " << d2 << " " << d3 << "\n";
}

std::vector<double> cds(12);
cds[0] = rMat(0,0);		cds[1] =  rMat(0,1);	cds[2] =  rMat(0,2);
cds[3] = rMat(1,0);		cds[4] =  rMat(1,1);	cds[5] =  rMat(1,2);
cds[6] = rMat(2,0);		cds[7] =  rMat(2,1);	cds[8] =  rMat(2,2);
cds[9] = rMat(3,0);		cds[10] = rMat(3,1);	cds[11] = rMat(3,2);
DrawDelaunayMesh (x0, y0, x1, y1, x2, y2, x3, y3, cds, flip, serpFlip);
}
}

*/


	/*
	void CurveRenderer::DrawDelaunayMesh (double  x0, double  y0,
	double  x1, double  y1,
	double  x2, double  y2,
	double  x3, double  y3,
	std::vector<double> st,
	bool flip, bool serpFlip)
	{
	int val = 0;
	if (!IsClockwise(x0, y0, x1, y1, x2, y2))
	val = 1;
	if(flip) val = 1 - val;

	if(val == 0)
	{
	st[0] = -st[0];		st[1] = -st[1];
	st[3] = -st[3];		st[4] = -st[4];
	st[6] = -st[6];		st[7] = -st[7];
	st[9] = -st[9];		st[10] = -st[10];
	}

	//shaderProgram->setUniformValue("flag", val);

	// 1
	if (!IsInsideCircle(x0, y0, x1, y1, x2, y2, x3, y3))
	DrawTriangle (x0, y0, x1, y1, x2, y2,
	st[0], st[1], st[2],
	st[3], st[4], st[5],
	st[6], st[7], st[8]);

	// 2
	if (!IsInsideCircle (x0, y0, x2, y2, x3, y3, x1, y1))
	DrawTriangle (x0, y0, x2, y2, x3, y3,
	st[0], st[1], st[2],
	st[6], st[7], st[8],
	st[9], st[10], st[11]);


	// 3
	if (!IsInsideCircle (x1, y1, x2, y2, x3, y3, x0, y0))
	DrawTriangle (x1, y1, x2, y2, x3, y3,
	st[3], st[4], st[5],
	st[6], st[7], st[8],
	st[9], st[10], st[11]);

	// 4
	if (!IsInsideCircle (x0, y0, x1, y1, x3, y3, x2, y2))
	DrawTriangle (x0, y0, x1, y1, x3, y3,
	st[0], st[1], st[2],
	st[3], st[4], st[5],
	st[9], st[10], st[11]);
	}
	*/