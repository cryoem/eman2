/**
 * $Id$
 */

#ifndef eman__vec3_h__
#define eman__vec3_h__ 1

#include <vector>
#include <math.h>
using std::vector;

namespace EMAN
{
	/** Template Vec3 defines a 3-element vector and various vector
     * operations. The vector may store any numeric data type
     * including int, float, double, etc.
     *
     * Some usage examples:
	 *
	 *    Vec3<float> v1(1,2,3);
	 *    Vec3<float> v2(2,4,6);
	 *    Vec3<float> v3 = v1 + v2;
	 *    double dot = v1 * v2;
     */
	template < class T > class Vec3
	{
	  public:
		Vec3() {
			vec[0] = 0;
			vec[1] = 0;
			vec[2] = 0;
		}

		template <class T1, class T2, class T3>
		Vec3(T1 x, T2 y, T3 z)
		{
			vec[0] = (T)x;
			vec[1] = (T)y;
			vec[2] = (T)z;
		}
		
		template <class T2>
		Vec3(const vector < T2 > &v)
		{
			vec[0] = (T) v[0];
			vec[1] = (T) v[1];
			vec[2] = (T) v[2];
		}
		
		template <class T2>
		Vec3(const Vec3 < T2 > &v)
		{
			vec[0] = (T) v[0];
			vec[1] = (T) v[1];
			vec[2] = (T) v[2];
		}
		
		template <class T2>
		Vec3<T> & operator=(const Vec3 < T2 > &v)
		{
			if (this != &v) {
				vec[0] = (T) v[0];
				vec[1] = (T) v[1];
				vec[2] = (T) v[2];
			}
			return *this;
		}


		virtual ~ Vec3() {
		}

		double normalize()
		{
			double len = length();
			if (len != 0) {
				(*this) *= (1.0 / len);
			}
			else {
				set_value(0, 0, 0);
			}
			return len;
		}

		double length() const
		{
			T t = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
			return sqrt((double) t);
		}

		template <class T2>
		double dot(const Vec3 < T2 > &v) const
		{
			return (vec[0] * v[0] + vec[1] * v[1] + vec[2] * v[2]);
		}

		template <class T2>
		Vec3 < T > cross(const Vec3 < T2 > &v) const
		{
			return Vec3((T) (vec[1] * v[2] - vec[2] * v[1]),
						(T) (vec[2] * v[0] - vec[0] * v[2]),
						(T) (vec[0] * v[1] - vec[1] * v[0]));
		}

		Vec3 < T > &negate()
		{
			vec[0] = -vec[0];
			vec[1] = -vec[1];
			vec[2] = -vec[2];
			return (*this);
		}

		vector < T > get_as_list()const
		{
			vector < T > v(3);
			v[0] = vec[0];
			v[1] = vec[1];
			v[2] = vec[2];
			return v;
		}
		
		template <class T2>
		void set_value(const vector < T2 > &v)
		{
			vec[0] = (T) v[0];
			vec[1] = (T) v[1];
			vec[2] = (T) v[2];
		}
		
		template <class T1, class T2, class T3>
		void set_value(T1 x, T2 y, T3 z)
		{
			vec[0] = (T) x;
			vec[1] = (T) y;
			vec[2] = (T) z;
		}

		T operator[] (int i) const
		{
			return vec[i];
		}
		T & operator[] (int i)
		{
			return vec[i];
		}
		
		template <class T2>
		Vec3 < T > &operator +=(const Vec3 < T2 > &v) {
			vec[0] +=  (T) v[0];
			vec[1] +=  (T) v[1];
			vec[2] +=  (T) v[2];
			return *this;
		}
		
		template <class T2>
		Vec3 < T > &operator +=(T2 d) {
			vec[0] +=  (T) d;
			vec[1] +=  (T) d;
			vec[2] +=  (T) d;
			return *this;
		}

		template <class T2>
		Vec3 < T > &operator -=(const Vec3 < T2 > &v) {
			vec[0] -= (T) v[0];
			vec[1] -= (T) v[1];
			vec[2] -= (T) v[2];
			return *this;
		}
		
		template <class T2>
		Vec3 < T > &operator -=(T2 d) {
			vec[0] -= (T) d;
			vec[1] -= (T) d;
			vec[2] -= (T) d;
			return *this;
		}
		
		
		template <class T2>
		Vec3 < T > &operator *=(T2 d) {
			vec[0] *= (T) d;
			vec[1] *= (T) d;
			vec[2] *= (T) d;
			return *this;
		}
		
		template <class T2>
		Vec3 < T > &operator /=(T2 d) {
			if (d != 0) {
				vec[0] /= (T) d;
				vec[1] /= (T) d;
				vec[2] /= (T) d;
			}
			return *this;
		}

	  private:
		T vec[3];
	};


	template < class T1, class T2 >
	Vec3 < float >operator +(const Vec3 < T1 > &v1, const Vec3 < T2 > &v2)
	{
		return Vec3 < float >(v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]);
	}

	template < class T1, class T2 >
	Vec3 < float >operator -(const Vec3 < T1 > &v1, const Vec3 < T2 > &v2)
	{
		return Vec3 < float >(v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]);
	}
	
	template < class T1, class T2 >
	double operator *(const Vec3 < T1 > &v1, const Vec3 < T2 > &v2)
	{
		return v1.dot(v2);
	}

	template < class T1, class T2 >
	Vec3 < T1 > operator *(T2 d, const Vec3 < T1 >& v)
	{
		Vec3 < T1 > v1 = v;
		v1 *= d;
		return v1;
	}

	template < class T1, class T2 >
	Vec3 < T1 > operator /(T2 d, const Vec3 < T1 >& v) {
		Vec3 < T1 > v1 = v;
		if (d != 0) {
			v1 /= d;
		}
		return v1;
	}

	
	template < class T1, class T2 >
	Vec3 < T1 > operator *(const Vec3 < T1 >& v, T2 d) {
		Vec3 < T1 > v1 = v;
		v1 *= d;
		return v1;
	}

	template < class T1, class T2 >
	Vec3 < T1 > operator /(const Vec3 < T1 >& v, T2 d) {
		Vec3 < T1 > v1 = v;
		if (d != 0) {
			v1 /= d;
		}
		return v1;
	}

	template < class T1, class T2 >
	bool operator ==(const Vec3 < T1 > &v1, const Vec3 < T2 > &v2) {
		if (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2]) {
			return true;
		}
		return false;
	}

	template < class T1, class T2 >
	bool operator !=(const Vec3 < T1 > &v1, const Vec3 < T2 > &v2) {
		if (v1[0] != v2[0] || v1[1] != v2[1] || v1[2] != v2[2]) {
			return true;
		}
		return false;
	}
}

#endif
