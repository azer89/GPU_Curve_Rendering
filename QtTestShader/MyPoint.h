
#ifndef __MYPOINT_H__
#define __MYPOINT_H__

#include <iostream>

namespace CVSystem
{
	typedef struct MyPoint
	{
	public:
		float x;
		float y;

		MyPoint()
		{
			this->x = -1;
			this->y = -1;
		}

		MyPoint(float x, float y)
		{
			this->x = x;
			this->y = y;
		}

		bool Invalid()
		{
			if(((int)x) == -1 && ((int)y) == -1) 
				return true;
			return false;
		}

		MyPoint Norm()
		{
			float vlength = std::sqrt( x * x + y * y );
			return MyPoint(this->x / vlength, this->y / vlength);
		}

		float Distance(MyPoint other)
		{
			float xDist = x - other.x;
			float yDist = y - other.y;
			return sqrt(xDist * xDist + yDist * yDist);
		}

		MyPoint MyPoint::operator+ (const MyPoint& other) { return MyPoint(x + other.x, y + other.y); }

		MyPoint MyPoint::operator- (const MyPoint& other) { return MyPoint(x - other.x, y - other.y); }

		bool MyPoint::operator== (const MyPoint& other) { return ((int)x == (int)other.x && (int)y == (int)other.y); }

		MyPoint MyPoint::operator+= (const MyPoint& other)
		{
			x += other.x;
			y += other.y;
			return *this;
		}

		MyPoint MyPoint::operator-= (const MyPoint& other)
		{
			x -= other.x;
			y -= other.y;
			return *this;
		}

		MyPoint MyPoint::operator/ (const float& val) { return MyPoint(x / val, y / val); }

		MyPoint MyPoint::operator* (const float& val) { return MyPoint(x * val, y * val); }

		MyPoint MyPoint::operator*= (const float& val)
		{
			x *= val;
			y *= val;
			return *this;
		}

		MyPoint MyPoint::operator/= (const float& val)
		{
			x /= val;
			y /= val;
			return *this;
		}

		double DiagonalLengthSquared()
		{
			return x * x + y * y;
		}
	};
}
#endif