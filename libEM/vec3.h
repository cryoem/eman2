/**
 * $Id$
 */

#ifndef eman__vec3_h__
#define eman__vec3_h__ 1

#include <vector>
using std::vector;

#include <math.h>

#define SQR(X) ((X)*(X))

namespace EMAN
{
	/** Vec3i defines a 3-element integer vector and various vector
     * operations. 
     *
     * Some usage examples:
	 *
	 *    Vec3i v1(1,2,3);
	 *    Vec3i v2(2,4,6);
	 *    Vec3i v3 = v1 + v2;
	 *    int dot = v1 * v2;
     */
	
	class Vec3i
	{
	public:
		Vec3i() {
			vec[0] = 0;
			vec[1] = 0;
			vec[2] = 0;
		}

		Vec3i(int x, int y, int z = 0)
		{
			vec[0] = x;
			vec[1] = y;
			vec[2] = z;
		}

		Vec3i(const vector < int > &v)
		{
			vec[0] = v[0];
			vec[1] = v[1];
			vec[2] = v[2];
		}
		
		
		virtual ~ Vec3i() {
		}

		float normalize()
		{
			float len = length();
			if (len != 0) {
				vec[0] = (int) (vec[0] / len);
				vec[1] = (int) (vec[1] / len);
				vec[2] = (int) (vec[2] / len);
			}
			else {
				set_value(0, 0, 0);
			}
			return len;
		}

		float length() const
		{
			float t = (float)(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
			return (float)sqrt(t);
		}

		int dot(const Vec3i & v) const
		{
			return (vec[0] * v[0] + vec[1] * v[1] + vec[2] * v[2]);
		}


		Vec3i cross(const Vec3i  & v) const
		{
			return Vec3i((vec[1] * v[2] - vec[2] * v[1]),
						 (vec[2] * v[0] - vec[0] * v[2]),
						 (vec[0] * v[1] - vec[1] * v[0]));
		}

		vector < int > as_list()const
		{
			vector < int > v(3);
			v[0] = vec[0];
			v[1] = vec[1];
			v[2] = vec[2];
			return v;
		}
		
		void set_value(const vector < int > &v)
		{
			vec[0] =  v[0];
			vec[1] =  v[1];
			vec[2] =  v[2];
		}
		
		void set_value(int x, int y, int z)
		{
			vec[0] =  x;
			vec[1] =  y;
			vec[2] =  z;
		}

		int operator[] (int i) const
		{
			return vec[i];
		}
		
		int & operator[] (int i)
		{
			return vec[i];
		}
		
		Vec3i &operator +=(const Vec3i &v) {
			vec[0] += v[0];
			vec[1] += v[1];
			vec[2] += v[2];
			return *this;
		}
		
		Vec3i  &operator +=(int d) {
			vec[0] +=   d;
			vec[1] +=   d;
			vec[2] +=   d;
			return *this;
		}

		Vec3i  &operator -=(const Vec3i &v) {
			vec[0] -=  v[0];
			vec[1] -=  v[1];
			vec[2] -=  v[2];
			return *this;
		}
		
		Vec3i  &operator -=(int d) {
			vec[0] -=  d;
			vec[1] -=  d;
			vec[2] -=  d;
			return *this;
		}
		
		Vec3i  &operator *=(int d) {
			vec[0] *=  d;
			vec[1] *=  d;
			vec[2] *=  d;
			return *this;
		}
		
		
		Vec3i &operator /=(int d) {
			if (d != 0) {
				vec[0] /=  d;
				vec[1] /=  d;
				vec[2] /=  d;
			}
			return *this;
		}

	  private:
		int vec[3];
	};

	
	inline Vec3i operator +(const Vec3i &v1, const Vec3i &v2)
	{
		return Vec3i(v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]);
	}
	
	inline Vec3i operator +(const Vec3i &v, int n)
	{
		Vec3i v1 = v;
		v1 += n;
		return v1;
	}
	
	inline Vec3i operator +(int n, const Vec3i &v)
	{
		Vec3i v1 = v;
		v1 += n;
		return v1;
	}

	
	inline Vec3i operator -(const Vec3i &v1, const Vec3i &v2)
	{
		return Vec3i(v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]);
	}	

	inline Vec3i operator -(const Vec3i &v, int n)
	{
		Vec3i v1 = v;
		v1 -= n;
		return v1;
	}
	
	inline Vec3i operator -(int n, const Vec3i &v)
	{
		Vec3i v1 = v;
		v1 -= n;
		return v1;
	}

	inline int operator *(const Vec3i &v1, const Vec3i &v2)
	{
		return v1.dot(v2);
	}
	
	inline Vec3i operator *(int d, const Vec3i & v)
	{
		Vec3i v1 = v;
		v1 *= d;
		return v1;
	}
	
	inline Vec3i operator *(const Vec3i & v, int d) {
		Vec3i v1 = v;
		v1 *= d;
		return v1;
	}
	
	inline Vec3i operator /(const Vec3i & v, int d) {
		Vec3i v1 = v;
		if (d != 0) {
			v1 /= d;
		}
		return v1;
	}

	
	inline bool operator ==(const Vec3i &v1, const Vec3i &v2) {
		if (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2]) {
			return true;
		}
		return false;
	}

	
	inline bool operator !=(const Vec3i &v1, const Vec3i &v2) {
		if (v1[0] != v2[0] || v1[1] != v2[1] || v1[2] != v2[2]) {
			return true;
		}
		return false;
	}

	/** Vec3f defines a 3-element float vector and various vector
     * operations. 
     *
     * Some usage examples:
	 *
	 *    Vec3f v1(1.1,2.2,3);
	 *    Vec3f v2(2,4.12,6.2);
	 *    Vec3f v3 = v1 + v2;
	 *    float dot = v1 * v2;
     */
	
	class Vec3f
	{
	public:
		Vec3f() {
			vec[0] = 0;
			vec[1] = 0;
			vec[2] = 0;
		}

		
		Vec3f(float x, float y, float z = 0)
		{
			vec[0] = x;
			vec[1] = y;
			vec[2] = z;
		}

		Vec3f(const vector < float > &v)
		{
			vec[0] = v[0];
			vec[1] = v[1];
			vec[2] = v[2];
		}
		
		Vec3f(const Vec3i &v)
		{
			vec[0] = (float)v[0];
			vec[1] = (float)v[1];
			vec[2] = (float)v[2];
		}

		virtual ~ Vec3f() {
		}

		float normalize()
		{
			float len = length();
			if (len != 0) {
				(*this) *= (1.0f / len);
			}
			else {
				set_value(0, 0, 0);
			}
			return len;
		}

		float length() const
		{
			float t = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
			return (float)sqrt(t);
		}

		float dot(const Vec3f & v) const
		{
			return (vec[0] * v[0] + vec[1] * v[1] + vec[2] * v[2]);
		}
#if 0
		float dot(const Vec3i & v) const
		{
			return (vec[0] * v[0] + vec[1] * v[1] + vec[2] * v[2]);
		}
		
		Vec3f cross(const Vec3i & v) const
		{
			return Vec3f((vec[1] * v[2] - vec[2] * v[1]),
						(vec[2] * v[0] - vec[0] * v[2]),
						(vec[0] * v[1] - vec[1] * v[0]));
		}
		
#endif
		Vec3f cross(const Vec3f & v) const
		{
			return Vec3f((vec[1] * v[2] - vec[2] * v[1]),
						(vec[2] * v[0] - vec[0] * v[2]),
						(vec[0] * v[1] - vec[1] * v[0]));
		}
		
		vector < float > as_list()const
		{
			vector < float > v(3);
			v[0] = vec[0];
			v[1] = vec[1];
			v[2] = vec[2];
			return v;
		}
		
		void set_value(const vector < float > &v)
		{
			vec[0] =  v[0];
			vec[1] =  v[1];
			vec[2] =  v[2];
		}
		
		void set_value(float x, float y, float z)
		{
			vec[0] =  x;
			vec[1] =  y;
			vec[2] =  z;
		}

		inline float at(int i) { return vec[i]; }
		
		float operator[] (int i) const
		{
			return vec[i];
		}
		
		float & operator[] (int i)
		{
			return vec[i];
		}
		
		Vec3f &operator +=(const Vec3f &v) {
			vec[0] += v[0];
			vec[1] += v[1];
			vec[2] += v[2];
			return *this;
		}
		
		Vec3f &operator +=(const Vec3i &v) {
			vec[0] += v[0];
			vec[1] += v[1];
			vec[2] += v[2];
			return *this;
		}
		
		Vec3f &operator +=(float d) {
			vec[0] +=   d;
			vec[1] +=   d;
			vec[2] +=   d;
			return *this;
		}

		Vec3f &operator -=(const Vec3f &v) {
			vec[0] -=  v[0];
			vec[1] -=  v[1];
			vec[2] -=  v[2];
			return *this;
		}
		
		Vec3f &operator -=(const Vec3i &v) {
			vec[0] -=  v[0];
			vec[1] -=  v[1];
			vec[2] -=  v[2];
			return *this;
		}
		
		Vec3f &operator -=(float d) {
			vec[0] -=  d;
			vec[1] -=  d;
			vec[2] -=  d;
			return *this;
		}
		
		Vec3f &operator *=(float d) {
			vec[0] *=  d;
			vec[1] *=  d;
			vec[2] *=  d;
			return *this;
		}
		
		
		Vec3f &operator /=(float d) {
			if (d != 0) {
				vec[0] /=  d;
				vec[1] /=  d;
				vec[2] /=  d;
			}
			return *this;
		}

	  private:
		float vec[3];
	};

	
	inline Vec3f operator +(const Vec3f &v1, const Vec3f &v2)
	{
		return Vec3f(v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]);
	}

	inline Vec3f operator -(const Vec3f &v1, const Vec3f &v2)
	{
		return Vec3f(v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]);
	}	
	
	inline float operator *(const Vec3f &v1, const Vec3f &v2)
	{
		return v1.dot(v2);
	}
		
	inline Vec3f operator +(const Vec3i &v1, const Vec3f &v2)
	{
		return Vec3f(v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]);
	}

	inline Vec3f operator +(const Vec3f &v1, const Vec3i &v2)
	{
		return Vec3f(v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]);
	}

	inline Vec3f operator -(const Vec3i &v1, const Vec3f &v2)
	{
		return Vec3f(v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]);
	}	
	
	inline Vec3f operator -(const Vec3f &v1, const Vec3i &v2)
	{
		return Vec3f(v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]);
	}
	
	inline Vec3f operator -(const Vec3f &v)
	{
		return Vec3f(-v[0],-v[1],-v[2]);
	}
	
	inline float operator *(const Vec3f &v1, const Vec3i &v2)
	{
		return v1.dot(v2);
	}
		
	inline float operator *(const Vec3i &v1, const Vec3f &v2)
	{
		return v2.dot(v1);
	}
	
	inline Vec3f operator +(float d, const Vec3f & v)
	{
		Vec3f v1 = v;
		v1 += d;
		return v1;
	}
	
	inline Vec3f operator +(const Vec3f & v, float d) {
		Vec3f v1 = v;
		v1 += d;
		return v1;
	}
	
	inline Vec3f operator -(float d, const Vec3f & v)
	{
		Vec3f v1 = v;
		v1 -= d;
		return v1;
	}
	
	inline Vec3f operator -(const Vec3f & v, float d) {
		Vec3f v1 = v;
		v1 -= d;
		return v1;
	}
	
	inline Vec3f operator *(float d, const Vec3f & v)
	{
		Vec3f v1 = v;
		v1 *= d;
		return v1;
	}
	
	inline Vec3f operator *(const Vec3f & v, float d) {
		Vec3f v1 = v;
		v1 *= d;
		return v1;
	}
	
	inline Vec3f operator /(const Vec3f & v, float d) {
		Vec3f v1 = v;
		if (d != 0) {
			v1 /= d;
		}
		return v1;
	}

	
	inline bool operator ==(const Vec3f &v1, const Vec3f &v2) {
		if (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2]) {
			return true;
		}
		return false;
	}

	
	inline bool operator !=(const Vec3f &v1, const Vec3f &v2) {
		if (v1[0] != v2[0] || v1[1] != v2[1] || v1[2] != v2[2]) {
			return true;
		}
		return false;
	}
	
}
#endif
