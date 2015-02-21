

#ifndef __VERTEX_H__
#define __VERTEX_H__

#include <iostream>
#include "MyPoint.h"


	typedef struct MyVertex
	{
	public:
		CVSystem::MyPoint xyCoor;

		double k;
		double l;
		double m;

		MyVertex()
		{
			this->xyCoor.x = -1;
			this->xyCoor.y = -1;

			this->k = -1;
			this->l = -1;
			this->m = -1;
		}

		MyVertex(float x, float y, double k, double l, double m)
		{
			this->xyCoor.x = x;
			this->xyCoor.y = y;

			this->k = k;
			this->l = l;
			this->m = m;
		}
	};

#endif