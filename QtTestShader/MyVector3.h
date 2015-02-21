
#ifndef __MYVECTOR3_H__
#define __MYVECTOR3_H__

#include <iostream>

namespace CVSystem
{
	typedef struct MyVector3
	{
	public:
		double x;
		double y;
		double z;

		MyVector3()
		{
			this->x = -1;
			this->y = -1;
			this->z = -1;
		}

		MyVector3(double x, double y, double z)
		{
			this->x = x;
			this->y = y;
			this->z = z;
		}

		bool Invalid()
		{
			if(((int)x) == -1 && ((int)y) == -1 && ((int)z) == -1) 
				return true;
			return false;
		}

		MyVector3 Norm()
		{
			float vlength = std::sqrt(x * x + y * y + z * z);
			return MyVector3(this->x / vlength, this->y / vlength, this->z / vlength);
		}

		float Distance(MyVector3 other)
		{
			float xDist = x - other.x;
			float yDist = y - other.y;
			float zDist = z - other.z;
			return sqrt(xDist * xDist + yDist * yDist + zDist * zDist);
		}

		MyVector3 MyVector3::operator+ (const MyVector3& other) { return MyVector3(x + other.x, y + other.y, z + other.z); }

		MyVector3 MyVector3::operator- (const MyVector3& other) { return MyVector3(x - other.x, y - other.y, z - other.z); }

		//bool MyVector3::operator== (const MyVector3& other) { return ((int)x == (int)other.x && (int)y == (int)other.y && (int)z == (int)other.z); }

		MyVector3 MyVector3::operator+= (const MyVector3& other)
		{
			x += other.x;
			y += other.y;
			z += other.z;
			return *this;
		}

		MyVector3 MyVector3::operator-= (const MyVector3& other)
		{
			x -= other.x;
			y -= other.y;
			z -= other.z;
			return *this;
		}

		MyVector3 MyVector3::operator/ (const double& val) { return MyVector3(x / val, y / val, z / val); }

		MyVector3 MyVector3::operator* (const double& val) { return MyVector3(x * val, y * val, z / val); }

		/*MyVector3 MyVector3::operator*= (const double& val)
		{
			x *= val;
			y *= val;
			z *= val;
			return *this;
		}*/

		/*MyVector3 MyVector3::operator/= (const double& val)
		{
			x /= val;
			y /= val;
			z /= val;
			return *this;
		}*/
	};
}
#endif