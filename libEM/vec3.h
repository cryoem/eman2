/**
 * $Id$
 */

/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
 * Copyright (c) 2000-2006 Baylor College of Medicine
 * 
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 * 
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * 
 * */

#ifndef eman__vec3_h__
#define eman__vec3_h__ 1

#include <vector>
using std::vector;

#include <cmath>

namespace EMAN
{
	/** @ingroup tested3c */
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

		/** contruct a Vec3i object with all elements equal to 0.
		 */
		Vec3i() {
			vec[0] = 0;
			vec[1] = 0;
			vec[2] = 0;
		}
		
		/** contruct a Vec3i object given (x,y) or (x,y,z) values.
		 * @param x  Value of the first item.
		 * @param y  Value of the second item.
		 * @param z  Value of the third item. If not specified,
		 * default to 0.
		 */
		Vec3i(int x, int y, int z = 0)
		{
			vec[0] = x;
			vec[1] = y;
			vec[2] = z;
		}

		/** Construct a Vec3i object given a std::vector object. The
		 * std::vector object should have at least 3 items.
		 * @param v The std::vector object. It should have at least 3 items.
		 */
		Vec3i(const vector < int > &v)
		{
			vec[0] = v[0];
			vec[1] = v[1];
			vec[2] = v[2];
		}

		Vec3i(const Vec3i &v)
		{
			vec[0] = v[0];
			vec[1] = v[1];
			vec[2] = v[2];
		}
		
		
		virtual ~ Vec3i() {
		}

		/** Normalize the vector and return its length before the
		 * normalization.
		 * @return The length of the Vec before normalization.
		 */
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
		
		/** Calculate its length.
		 * @return The vector's length.
		 */
		float length() const
		{
			float t = (float)(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
			return (float)sqrt(t);
		}
		
		/** Calculate its squared length.
		 * no sqrt called 
		 * @return The vector's length squared.
		 */
		float squared_length() const
		{
			float t = static_cast<float>( vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2] );
			return t;
		}
		
		/** Calculate the dot product of 'this' vector with a second
		 * vector.
		 * @param v  The second vector to do the dot product.
		 * @return The dot product.
		 */
		int dot(const Vec3i & v) const
		{
			return (vec[0] * v[0] + vec[1] * v[1] + vec[2] * v[2]);
		}
		
		/** Calculate the cross product of 'this' vector with a second
		 * vector.
		 * @param v  The second vector to do the cross product.
		 * @return The cross product.
		 */
		Vec3i cross(const Vec3i  & v) const
		{
			return Vec3i((vec[1] * v[2] - vec[2] * v[1]),
						 (vec[2] * v[0] - vec[0] * v[2]),
						 (vec[0] * v[1] - vec[1] * v[0]));
		}

		/** Return the values of this vector as a std::vector.
		 * @return The std::vector version of this vector.
		 */
		vector < int > as_list()const
		{
			vector < int > v(3);
			v[0] = vec[0];
			v[1] = vec[1];
			v[2] = vec[2];
			return v;
		}
		
		/** Set new values using a std::vector object.
		 * @param v A std::vector object used to set 'this' vector's value.
		 *  It should have at least 3 items.
		 */
		void set_value(const vector < int > &v)
		{
			vec[0] =  v[0];
			vec[1] =  v[1];
			vec[2] =  v[2];
		}

		/** Set values at a particular index
		 * @param index  The index to be set
		 * @param value  The value to be set
		 */
		void set_value_at(int index, int value)
		{
			vec[index] =  value;
		}

		/** Set new values to this vector object.
		 * @param x Value of the first item.
		 * @param y Value of the second item.
		 * @param z Value of the third item.
		 */
		void set_value(int x, int y, int z)
		{
			vec[0] =  x;
			vec[1] =  y;
			vec[2] =  z;
		}

		/** Get the ith item of the vector. Used in the right side of
		 * the assignment.
		 *
		 * @param i The index of the item to get. Its validality is
		 * not checked.
		 * @return The ith item of the vector.
		 */
		int operator[] (int i) const
		{
			return vec[i];
		}

		/** Get the ith item of the vector. Used in the left side of
		 * the assignment.
		 *
		 * @param i The index of the item to get. Its validality is
		 * not checked.
		 * @return The ith item of the vector.
		 */
		int & operator[] (int i)
		{
			return vec[i];
		}

		/** 'this' += v; Add the 2 vectors by adding item by item.
		 * @param v The vector used to be added to 'this' vector.
		 * @return The new 'this' as a result of add.
		 */
		Vec3i &operator +=(const Vec3i &v) {
			vec[0] += v[0];
			vec[1] += v[1];
			vec[2] += v[2];
			return *this;
		}

		/** 'this' += d. Add d to each item of this vector.
		 * @param d The number used to be added to this vector.
		 * @return The new 'this' as a result of add.
		 */
		Vec3i  &operator +=(int d) {
			vec[0] +=   d;
			vec[1] +=   d;
			vec[2] +=   d;
			return *this;
		}
		
		/** 'this' -= v; Minus the 2 vectors item by item.
		 * @param v The vector used to be substracted from 'this' vector.
		 * @return The new 'this' as a result of substraction.
		 */
		Vec3i  &operator -=(const Vec3i &v) {
			vec[0] -=  v[0];
			vec[1] -=  v[1];
			vec[2] -=  v[2];
			return *this;
		}
		
		/** 'this' -= d; Minus a number from each item of 'this' vector.
		 * @param d The number used to be substracted from 'this' vector.
		 * @return The new 'this' as a result of substraction.
		 */
		Vec3i  &operator -=(int d) {
			vec[0] -=  d;
			vec[1] -=  d;
			vec[2] -=  d;
			return *this;
		}

		/** 'this' *= d; Multiply a number on each item of 'this' vector.
		 * @param d The number to multiply.
		 * @return The new 'this' as a result of multiplication.
		 */
		Vec3i  &operator *=(int d) {
			vec[0] *=  d;
			vec[1] *=  d;
			vec[2] *=  d;
			return *this;
		}
		
		/** 'this' /= d; Divide a number on each item of 'this' vector.
		 * @param d The number to divide.
		 * @return The new 'this' as a result of division.
		 */
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

	/** @ingroup tested3c */
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
		/** contruct a Vec3f object with all elements equal to 0.
		 */
		Vec3f() {
			vec[0] = 0;
			vec[1] = 0;
			vec[2] = 0;
		}

			
		/** contruct a Vec3f object given (x,y) or (x,y,z) values.
		 * @param x  Value of the first item.
		 * @param y  Value of the second item.
		 * @param z  Value of the third item. If not specified,
		 * default to 0.
		 */
		Vec3f(float x, float y, float z = 0)
		{
			vec[0] = x;
			vec[1] = y;
			vec[2] = z;
		}

		/** Construct a Vec3f object given a std::vector object. The
		 * std::vector object should have at least 3 items.
		 * @param v The std::vector object. It should have at least 3 items.
		 */
		Vec3f(const vector < float > &v)
		{
			vec[0] = v[0];
			vec[1] = v[1];
			vec[2] = v[2];
		}

		/** convert a Vec3i vector to a Vec3f vector.
		 */
		Vec3f(const Vec3i &v)
		{
			vec[0] = (float)v[0];
			vec[1] = (float)v[1];
			vec[2] = (float)v[2];
		}

		/** copy constructor
		 */
		Vec3f(const Vec3f &v)
		{
			vec[0] = v[0];
			vec[1] = v[1];
			vec[2] = v[2];
		}

		virtual ~ Vec3f() {
		}
		
		/** Normalize the vector and return its length before the
		 * normalization.
		 * @return The length of the Vec before normalization.
		 */
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

		/** Calculate its length.
		 * @return The vector's length.
		 */
		float length() const
		{
			float t = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
			return (float)sqrt(t);
		}
		
		/** Calculate its squared length.
		 * The same as calling length except sqrt is not called.
		 * @return The vector's length squared.
		 */
		float squared_length() const
		{
			float t = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
			return t;
		}

		/** Calculate the dot product of 'this' vector with a second
		 * vector.
		 * @param v  The second vector to do the dot product.
		 * @return The dot product.
		 */
		float dot(const Vec3f & v) const
		{
			return (vec[0] * v[0] + vec[1] * v[1] + vec[2] * v[2]);
		}

		/** Calculate the cross product of 'this' vector with a second
		 * vector.
		 * @param v  The second vector to do the cross product.
		 * @return The cross product.
		 */
		Vec3f cross(const Vec3f & v) const
		{
			return Vec3f((vec[1] * v[2] - vec[2] * v[1]),
						 (vec[2] * v[0] - vec[0] * v[2]),
						 (vec[0] * v[1] - vec[1] * v[0]));
		}
		
		/** Return the values of this vector as a std::vector.
		 * @return The std::vector version of this vector.
		 */
		vector < float > as_list()const
		{
			vector < float > v(3);
			v[0] = vec[0];
			v[1] = vec[1];
			v[2] = vec[2];
			return v;
		}
		
		/** Set new values using a std::vector object.
		 * @param v A std::vector object used to set 'this' vector's value.
		 *  It should have at least 3 items.
		 */
		void set_value(const vector < float > &v)
		{
			vec[0] =  v[0];
			vec[1] =  v[1];
			vec[2] =  v[2];
		}
		/** Set values at a particular index
		 * @param index  The index to be set
		 * @param value  The value to be set
		 */
		void set_value_at(int index, float value)
		{
			vec[index] =  value;
		}

		
		/** Set new values to this vector object.
		 * @param x Value of the first item.
		 * @param y Value of the second item.
		 * @param z Value of the third item.
		 */
		void set_value(float x, float y, float z)
		{
			vec[0] =  x;
			vec[1] =  y;
			vec[2] =  z;
		}
		
		/** Get the ith item of the vector.
		 */
		inline float at(int i) { return vec[i]; }
		
		/** Get the ith item of the vector. Used in the right side of
		 * the assignment.
		 *
		 * @param i The index of the item to get. Its validality is
		 * not checked.
		 * @return The ith item of the vector.
		 */
		float operator[] (int i) const
		{
			return vec[i];
		}
		
		/** Get the ith item of the vector. Used in the left side of
		 * the assignment.
		 *
		 * @param i The index of the item to get. Its validality is
		 * not checked.
		 * @return The ith item of the vector.
		 */
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
