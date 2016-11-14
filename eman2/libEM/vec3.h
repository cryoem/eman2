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

#include <iostream>
using std::cout;
using std::endl;

namespace EMAN
{

	/** The Vec4 object is a templated object, intended to instantiated
	 * with basic types such as int, float, double etc. You may try
	 * to use other more generic types such as classes but you may get
	 * bad results.
	 * Note that the normalize and length operations are precise only to 32 bits
	 * Note there are convenient typedef so one needn't bother about using template
	 * terminology
	 *
	 * WARNING This class is not as complete as Vec4f, in that I have not yet implmented object operators
	 *
	 * For now, only functionality is for normalization
	 * typedef Vec4<float> Vec4f;
	 * typedef Vec4<int> Vec4i;
	 * typedef Vec4<double> Vec4d; // Not recommended for use unless precision is addressed in this class
	 * @author John Flanagan
	 * @date July 2011
	 */
	template<typename Type>
	class Vec4
	{
		public:
		/** One can always cast to the type of a Vec4 by accessing Vec4<Type>::type
		*/
		typedef Type type;
		
		Vec4() 
		{
			vec[0] = static_cast<Type>(0);
			vec[1] = static_cast<Type>(0);
			vec[2] = static_cast<Type>(0);
			vec[3] = static_cast<Type>(0);
		}
		
		template<typename Type2, typename Type3, typename Type4, typename Type5>
		Vec4(const Type2& a, const Type3& b, const Type4& c, const Type5& d)
		{
			vec[0] = static_cast<Type>(a);
			vec[1] = static_cast<Type>(b);
			vec[2] = static_cast<Type>(c);
			vec[3] = static_cast<Type>(d);
		}
		
		/** Construct a Vec3 object given a std::vector object. The
		 * std::vector object should have at least 3 items.
		 * @param v The std::vector object. It should have at least 3 items.
		 */
		template<typename Type2>
		Vec4(const vector < Type2 > &v)
		{
			vec[0] = static_cast<Type>(v[0]);
			vec[1] = static_cast<Type>(v[1]);
			vec[2] = static_cast<Type>(v[2]);
			vec[3] = static_cast<Type>(v[3]);
		}

		/** Copy constructor copies vector elements
		*/
		template<typename Type2>
		Vec4(const Vec4<Type2> &v)
		{
			vec[0] = v[0];
			vec[1] = v[1];
			vec[2] = v[2];
			vec[2] = v[3];
		}
		
		/** Destructor
		*/
		~Vec4() {}
		
		/** Calculate its length.
		 * @return The vector's length
		 * Warning - float precision
		 */
		float length() const
		{
			float t = (float)(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2] + vec[3] * vec[3]);
			return (float)sqrt(t);
		}
		
		/** Normalize the vector and return its length before the
		 * normalization.
		 * @return The length of the Vec before normalization.
		 */
		float normalize()
		{
			// Warning - float precision
			float len = length();
			if (len != 0) {
				vec[0] = static_cast<Type> (vec[0] / len);
				vec[1] = static_cast<Type> (vec[1] / len);
				vec[2] = static_cast<Type> (vec[2] / len);
				vec[3] = static_cast<Type> (vec[3] / len);
			}
			else {
				set_value(0, 0, 0, 0);	
			}
			return len;
		}
		
		/** Set new values using a std::vector object.
		 * @param v A std::vector object used to set 'this' vector's value.
		 *  It should have at least 3 items.
		 */
		template<typename Type2>
		void set_value(const vector < Type2 > &v)
		{
			vec[0] =  static_cast<Type>(v[0]);
			vec[1] =  static_cast<Type>(v[1]);
			vec[2] =  static_cast<Type>(v[2]);
			vec[3] =  static_cast<Type>(v[3]);
		}

		/** Set values at a particular index
		 * @param index  The index to be set
		 * @param value  The value to be set
		 */
		template<typename Type2>
		void set_value_at(int index, const Type2& value)
		{
			vec[index] = static_cast<Type>(value);
		}

		/** Set new values to this vector object.
		 * @param x Value of the first item.
		 * @param y Value of the second item.
		 * @param z Value of the third item.
		 */
		void set_value(const Type& a, const Type& b, const Type& c, const Type& d)
		{
			vec[0] =  a;
			vec[1] =  b;
			vec[2] =  c;
			vec[3] =  d;
		}

		/** Get the ith item of the vector. Used in the left side of
		* the assignment.
		*
		* @param i The index of the item to get. Its validality is
		* not checked.
		* @return The ith item of the vector.
		*/
		inline Type at(int i)
		{
			return vec[i];
		}
		
		/** For python __len__
		 *
		 * @return the number of elements in this container. it's always 3.
		 * */
		int number_of_element()
		{
			return 4;
		}

		/** Add this function to make it iterable in Python,
		 * so we can call list() or tuple() to convert Vec3f in python to a list or tuple.
		 *
		 * @return the iterator (here is the pointer) of the first element
		 * */
		Type * begin()
		{
			return &vec[0];
		}

		/** Add this function to make it iterable in Python,
		 * so we can call list() or tuple() to convert Vec3f in python to a list or tuple.
		 *
		 * @return the iterator (here is the pointer) of the one beyond the last element.
		 * */
		Type * end()
		{
			return &vec[4];
		}
		
		/** Get the ith item of the vector. Used in the right side of
		 * the assignment.
		 *
		 * @param i The index of the item to get. Its validality is
		 * not checked.
		 * @return The ith item of the vector.
		 */
		inline Type operator[] (int i) const
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
		inline Type& operator[] (int i)
		{
			return vec[i];
		}
		
		private:
			Type vec[4];
	};
	
	typedef Vec4<float> Vec4f;
	typedef Vec4<int> Vec4i;
	typedef Vec4<double> Vec4d;
	
	/** The Vec3 object is a templated object, intended to instantiated
	 * with basic types such as int, float, double etc. You may try
	 * to use other more generic types such as classes but you may get
	 * bad results.
	 * Note that the normalize and length operations are precise only to 32 bits
	 * Note there are convenient typedef so one needn't bother about using template
	 * terminology
	 * typedef Vec3<float> Vec3f;
	 * typedef Vec3<int> Vec3i;
	 * typedef Vec3<double> Vec3d; // Not recommended for use unless precision is addressed in this class
	 * @author David Woolford (based on the work of who ever wrote the original Vec3f and Vec3i classes - extracted into a template)
	 * @date August 2008
	 * @ingroup tested3c REALLY?
	 */
	template<typename Type>
	class Vec3
	{
		public:
		/** One can always cast to the type of a Vec3 by accessing Vec3<Type>::type
		*/
		typedef Type type;

		/** contruct a Vec3 object with all elements equal to 0.
		 */
		Vec3() /*: vec[0](0),vec[1](0),vec[2](0)*/ {

			vec[0] = static_cast<Type>(0);
			vec[1] = static_cast<Type>(0);
			vec[2] = static_cast<Type>(0);
		}

		/** contruct a Vec3 object given (x,y) or (x,y,z) values.
		 * @param x  Value of the first item.
		 * @param y  Value of the second item.
		 * @param z  Value of the third item. If not specified,
		 * default to 0.
		 */
		template<typename Type2, typename Type3, typename Type4>
		Vec3(const Type2& x, const Type3& y, const Type4& z = 0)
		{
			vec[0] = static_cast<Type>(x);
			vec[1] = static_cast<Type>(y);
			vec[2] = static_cast<Type>(z);
		}

		/** Construct a Vec3 object given a std::vector object. The
		 * std::vector object should have at least 3 items.
		 * @param v The std::vector object. It should have at least 3 items.
		 */
		template<typename Type2>
		Vec3(const vector < Type2 > &v)
		{
			vec[0] = static_cast<Type>(v[0]);
			vec[1] = static_cast<Type>(v[1]);
			vec[2] = static_cast<Type>(v[2]);
		}

		/** Copy constructor copies vector elements
		*/
		template<typename Type2>
		Vec3(const Vec3<Type2> &v)
		{
			vec[0] = v[0];
			vec[1] = v[1];
			vec[2] = v[2];
		}

		/** Destructor
		*/
		~Vec3() {}


		/** Normalize the vector and return its length before the
		 * normalization.
		 * @return The length of the Vec before normalization.
		 */
		float normalize()
		{
			// Warning - float precision
			float len = length();
			if (len != 0) {
				vec[0] = static_cast<Type> (vec[0] / len);
				vec[1] = static_cast<Type> (vec[1] / len);
				vec[2] = static_cast<Type> (vec[2] / len);
			}
			else {
				set_value(0, 0, 0);
			}
			return len;
		}


		/** Calculate its length.
		 * @return The vector's length
		 * Warning - float precision
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
		Type squared_length() const
		{
			return  vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2] ;
		}

		/** Calculate the dot product of 'this' vector with a second
		 * vector.
		 * @param v  The second vector to do the dot product.
		 * @return The dot product.
		 */
		template<typename Type2>
		Type dot(const Vec3<Type2> & v) const
		{
			return static_cast<Type>((vec[0] * v[0] + vec[1] * v[1] + vec[2] * v[2]));
		}

		/** Calculate the cross product of 'this' vector with a second
		 * vector.
		 * @param v  The second vector to do the cross product.
		 * @return The cross product.
		 */
		template<typename Type2>
		Vec3<Type> cross(const Vec3<Type2>  & v) const
		{
			return Vec3<Type>((vec[1] * v[2] - vec[2] * v[1]),
						  (vec[2] * v[0] - vec[0] * v[2]),
						   (vec[0] * v[1] - vec[1] * v[0]));
		}

		/** Return the values of this vector as a std::vector.
		 * @return The std::vector version of this vector.
		 */
		vector<Type> as_list() const
		{
			vector < Type > v(3);
			v[0] = vec[0];
			v[1] = vec[1];
			v[2] = vec[2];
			return v;
		}

		/** Set new values using a std::vector object.
		 * @param v A std::vector object used to set 'this' vector's value.
		 *  It should have at least 3 items.
		 */
		template<typename Type2>
		void set_value(const vector < Type2 > &v)
		{
			vec[0] =  static_cast<Type>(v[0]);
			vec[1] =  static_cast<Type>(v[1]);
			vec[2] =  static_cast<Type>(v[2]);
		}

		/** Set values at a particular index
		 * @param index  The index to be set
		 * @param value  The value to be set
		 */
		template<typename Type2>
		void set_value_at(int index, const Type2& value)
		{
			vec[index] = static_cast<Type>(value);
		}

		/** Set new values to this vector object.
		 * @param x Value of the first item.
		 * @param y Value of the second item.
		 * @param z Value of the third item.
		 */
		void set_value(const Type& x, const Type& y,  const Type& z)
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
		inline Type operator[] (int i) const
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
		inline Type& operator[] (int i)
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
		inline Type at(int i)
		{
			return vec[i];
		}

		/** For python __len__
		 *
		 * @return the number of elements in this container. it's always 3.
		 * */
		int number_of_element()
		{
			return 3;
		}

		/** Add this function to make it iterable in Python,
		 * so we can call list() or tuple() to convert Vec3f in python to a list or tuple.
		 *
		 * @return the iterator (here is the pointer) of the first element
		 * */
		Type * begin()
		{
			return &vec[0];
		}

		/** Add this function to make it iterable in Python,
		 * so we can call list() or tuple() to convert Vec3f in python to a list or tuple.
		 *
		 * @return the iterator (here is the pointer) of the one beyond the last element.
		 * */
		Type * end()
		{
			return &vec[3];
		}

		/** 'this' += v; Add the 2 vectors by adding item by item.
		 * @param v The vector used to be added to 'this' vector.
		 * @return The new 'this' as a result of add.
		 */
		template<typename Type2>
		Vec3<Type>& operator +=(const Vec3<Type2> &v) {
			vec[0] = static_cast<Type>(vec[0]+v[0]);
			vec[1] = static_cast<Type>(vec[1]+v[1]);
			vec[2] = static_cast<Type>(vec[2]+v[2]);
			return *this;
		}

		/** 'this' += d. Add d to each item of this vector.
		 * @param d The number used to be added to this vector.
		 * @return The new 'this' as a result of add.
		 */
		template<typename Type2>
		Vec3<Type>  &operator +=(const Type2& d) {
			vec[0] = static_cast<Type>(vec[0]+d);
			vec[1] = static_cast<Type>(vec[1]+d);
			vec[2] = static_cast<Type>(vec[2]+d);
			return *this;
		}

		/** 'this' -= v; Minus the 2 vectors item by item.
		 * @param v The vector used to be substracted from 'this' vector.
		 * @return The new 'this' as a result of substraction.
		 */
		template<typename Type2>
		Vec3<Type>  &operator -=(const Vec3<Type2> &v) {
			vec[0] = static_cast<Type>(vec[0]-v[0]);
			vec[1] = static_cast<Type>(vec[1]-v[1]);
			vec[2] = static_cast<Type>(vec[2]-v[2]);
			return *this;
		}

		/** 'this' -= d; Minus a number from each item of 'this' vector.
		 * @param d The number used to be substracted from 'this' vector.
		 * @return The new 'this' as a result of substraction.
		 */
		template<typename Type2>
		Vec3<Type>& operator -=(const Type2& d) {
			vec[0] = static_cast<Type>(vec[0]-d);
			vec[1] = static_cast<Type>(vec[1]-d);
			vec[2] = static_cast<Type>(vec[2]-d);
			return *this;
		}

		/** 'this' *= d; Multiply a number on each item of 'this' vector.
		 * @param d The number to multiply.
		 * @return The new 'this' as a result of multiplication.
		 */
		template<typename Type2>
		Vec3<Type>  &operator *=(const Type2& d) {
			vec[0] = static_cast<Type>(vec[0]*d);
			vec[1] = static_cast<Type>(vec[1]*d);
			vec[2] = static_cast<Type>(vec[2]*d);
			return *this;
		}

		/** 'this' /= d; Divide a number on each item of 'this' vector.
		 * @param d The number to divide.
		 * @return The new 'this' as a result of division.
		 */
		template<typename Type2>
		Vec3<Type>  &operator /=(const Type2& d) {
			vec[0] = static_cast<Type>(vec[0]/d);
			vec[1] = static_cast<Type>(vec[1]/d);
			vec[2] = static_cast<Type>(vec[2]/d);
			return *this;
		}

		template<typename Type2>
		operator vector<Type2>() const {
			vector<Type2> v(vec,vec+3);
			return v;
		}


		private:
			Type vec[3];
	};

	template<typename Type,typename Type2>
	inline Vec3<Type> operator +(const Vec3<Type> &v1, const Vec3<Type2> &v2)
	{

		return Vec3<Type>(static_cast<Type>(v1[0] + v2[0]), static_cast<Type>(v1[1] + v2[1]),static_cast<Type>(v1[2] + v2[2]));;
	}

	template<typename Type,typename Type2>
	inline Vec3<Type> operator +(const Vec3<Type> &v, const Type2& n)
	{
		Vec3<Type> v1(v);
		v1 += n;
		return v1;
	}

// 	template<typename Type,typename Type2>
// 	inline Vec3<Type> operator +(const Type2& n, const Vec3<Type> &v)
// 	{
// 		Vec3<Type> v1(v);
// 		v1 += n;
// 		return v1;
// 	}

	template<typename Type,typename Type2>
	inline Vec3<Type> operator -(const Vec3<Type> &v1, const Vec3<Type2> &v2)
	{
		return Vec3<Type>(static_cast<Type>(v1[0] - v2[0]),
						  static_cast<Type>(v1[1] - v2[1]),
						  static_cast<Type>(v1[2] - v2[2]));
	}

	template<typename Type,typename Type2>
	inline Vec3<Type> operator -(const Vec3<Type> &v, const Type2& n)
	{
		Vec3<Type> v1(v);
		v1 -= n;
		return v1;
	}
	template<typename Type>
	inline Vec3<Type> operator -(const Vec3<Type> &v)
	{
		return Vec3<Type>(-v[0],-v[1],-v[2]);
	}

// 	template<typename Type,typename Type2>
// 	inline Vec3<Type> operator -(const Type2& n, const Vec3<Type> &v)
// 	{
// 		Vec3<Type> v1(v);
// 		v1 -= n;
// 		return v1;
// 	}

	template<typename Type,typename Type2>
	inline Type operator *(const Vec3<Type> &v1, const Vec3<Type2> &v2)
	{
		return v1.dot(v2);
	}

	template<typename Type,typename Type2>
	inline Vec3<Type2> operator *(const Type& d, const Vec3<Type2> & v)
	{
		// Preserve the vector type
		Vec3<Type2> v1(v);
		v1 *= d;
		return v1;
	}

	template<typename Type,typename Type2>
	inline Vec3<Type> operator *(const Vec3<Type> & v,const Type2& d) {
		// Preserve the vector type
		Vec3<Type> v1(v);
		v1 *= d;
		return v1;
	}

	template<typename Type,typename Type2>
	inline Vec3<Type2> operator /(const Type& d, const Vec3<Type2> & v)
	{
		// Preserve the vector type
		Vec3<Type2> v1(v);
		v1 /= d;
		return v1;
	}

	template<typename Type,typename Type2>
	inline Vec3<Type> operator /(const Vec3<Type> & v,const Type2& d) {
		// Preserve the vector type
		Vec3<Type> v1(v);
		v1 /= d;
		return v1;
	}

	template<typename Type,typename Type2>
	inline bool operator ==(const Vec3<Type> &v1, const Vec3<Type2> &v2) {
		if (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2]) {
			return true;
		}
		return false;
	}

	template<typename Type,typename Type2>
	inline bool operator !=(const Vec3<Type> &v1, const Vec3<Type2> &v2) {
		if (v1[0] != v2[0] || v1[1] != v2[1] || v1[2] != v2[2]) {
			return true;
		}
		return false;
	}

	typedef Vec3<float> Vec3f;
	typedef Vec3<int> Vec3i;
	typedef Vec3<double> Vec3d;


	/** The Vec2 is precisely the same as Vec3 except it works exclusively in 2D
	 * Note there are convenient typedef so one needn't bother about using template
	 * terminology
	 * typedef Vec2<float> Vec2f;
	 * typedef Vec2<int> Vec2i;
	 * typedef Vec2double> Vec2d; // Not recommended for use unless precision is addressed in this class
	 * @author David Woolford (based on the work of who ever wrote the original Vec3f and Vec3i classes - extracted into a template)
	 * @date August 2008
	 * @ingroup tested3c REALLY?
	 */
	template<typename Type>
	class Vec2
	{
		public:
		/** One can always cast to the type of a Vec2 by accessing Vec2<Type>::type
		 */
		typedef Type type;

		/** contruct a Vec2 object with all elements equal to 0.
		 */
		Vec2() /*: vec[0](0),vec[1](0),vec[2](0)*/ {
			vec[0] = static_cast<Type>(0);
			vec[1] = static_cast<Type>(0);
		}

		/** contruct a Vec2 object given (x,y) or (x,y,z) values.
		* @param x  Value of the first item.
		* @param y  Value of the second item.
		* default to 0.
		 */
		template<typename Type2, typename Type3>
		Vec2(const Type2& x, const Type3& y)
		{
			vec[0] = static_cast<Type>(x);
			vec[1] = static_cast<Type>(y);
		}

		/** Construct a Vec2 object given a std::vector object. The
		* std::vector object should have at least 3 items.
		* @param v The std::vector object. It should have at least 3 items.
		 */
		template<typename Type2>
		Vec2(const vector < Type2 > &v)
		{
			vec[0] = static_cast<Type>(v[0]);
			vec[1] = static_cast<Type>(v[1]);
		}

		/** Copy constructor copies vector elements
		 */
		template<typename Type2>
		Vec2(const Vec2<Type2> &v)
		{
			vec[0] = static_cast<Type>(v[0]);
			vec[1] = static_cast<Type>(v[1]);
		}

		/** Destructor
		 */
		~Vec2() {}


		/** Normalize the vector and return its length before the
		 * normalization.
		 * @return The length of the Vec before normalization.
		 */
		float normalize()
		{
		// Warning - float precision
			float len = length();
			if (len != 0) {
				vec[0] = static_cast<Type> (vec[0] / len);
				vec[1] = static_cast<Type> (vec[1] / len);
			}
			else {
				set_value(0, 0);
			}
			return len;
		}


		/** Calculate its length.
		* @return The vector's length
		* Warning - float precision
		 */
		float length() const
		{
			float t = (float)(vec[0] * vec[0] + vec[1] * vec[1]);
			return (float)sqrt(t);
		}

		/** Calculate its squared length.
		 * no sqrt called
		 * @return The vector's length squared.
		 */
		Type squared_length() const
		{
			return  vec[0] * vec[0] + vec[1] * vec[1] ;
		}

		/** Calculate the dot product of 'this' vector with a second
		* vector.
		* @param v  The second vector to do the dot product.
		* @return The dot product.
		 */
		template<typename Type2>
		Type dot(const Vec2<Type2> & v) const
		{
			return static_cast<Type>((vec[0] * v[0] + vec[1] * v[1]));
		}

		/** Return the values of this vector as a std::vector.
		 * @return The std::vector version of this vector.
		 */
		vector<Type> as_list() const
		{
			vector < Type > v(2);
			v[0] = vec[0];
			v[1] = vec[1];
			return v;
		}

		/** Set new values using a std::vector object.
		 * @param v A std::vector object used to set 'this' vector's value.
		 *  It should have at least 3 items.
		 */
		template<typename Type2>
		void set_value(const vector < Type2 > &v)
		{
			vec[0] =  static_cast<Type>(v[0]);
			vec[1] =  static_cast<Type>(v[1]);;
		}

		/** Set values at a particular index
		 * @param index  The index to be set
		 * @param value  The value to be set
		 */
		template<typename Type2>
		void set_value_at(int index, const Type2& value)
		{
			vec[index] = static_cast<Type>(value);
		}

		/** Set new values to this vector object.
		 * @param x Value of the first item.
		 * @param y Value of the second item.
		 */
		void set_value(const Type& x, const Type& y)
		{
			vec[0] =  x;
			vec[1] =  y;
		}

		/** Get the ith item of the vector. Used in the right side of
		 * the assignment.
		 *
		 * @param i The index of the item to get. Its validality is
		 * not checked.
		 * @return The ith item of the vector.
		 */
		inline Type operator[] (int i) const { return vec[i]; }

		/** Get the ith item of the vector. Used in the left side of
		 * the assignment.
		 *
		 * @param i The index of the item to get. Its validality is
		  * not checked.
		* @return The ith item of the vector.
		 */
		inline Type& operator[] (int i) { return vec[i]; }

		/** Get the ith item of the vector. Used in the left side of
		 * the assignment.
		 *
		 * @param i The index of the item to get. Its validality is
		 * not checked.
		 * @return The ith item of the vector.
		 */
		inline Type at(int i) { return vec[i]; }


		/** For python __len__
		 *
		 * @return the number of elements in this container. it's always 2.
		 * */
		int number_of_element()
		{
			return 2;
		}

		/** Add this function to make it iterable in Python,
		 * so we can call list() or tuple() to convert Vec3f in python to a list or tuple.
		 *
		 * @return the iterator (here is the pointer) of the first element
		 * */
		Type * begin()
		{
			return &vec[0];
		}

		/** Add this function to make it iterable in Python,
		 * so we can call list() or tuple() to convert Vec3f in python to a list or tuple.
		 *
		 * @return the iterator (here is the pointer) of the one beyond the last element.
		 * */
		Type * end()
		{
			return &vec[2];
		}


		/** 'this' += v; Add the 2 vectors by adding item by item.
		 * @param v The vector used to be added to 'this' vector.
		 * @return The new 'this' as a result of add.
		 */
		template<typename Type2>
		Vec2<Type>& operator +=(const Vec2<Type2> &v) {
			vec[0] = static_cast<Type>(vec[0]+v[0]);
			vec[1] = static_cast<Type>(vec[1]+v[1]);
			return *this;
		}

		/** 'this' += d. Add d to each item of this vector.
		 * @param d The number used to be added to this vector.
		 * @return The new 'this' as a result of add.
		 */
		template<typename Type2>
		Vec2<Type>& operator +=(const Type2& d) {
			vec[0] = static_cast<Type>(vec[0]+d);
			vec[1] = static_cast<Type>(vec[1]+d);
			return *this;
		}

		/** 'this' -= v; Minus the 2 vectors item by item.
		 * @param v The vector used to be substracted from 'this' vector.
		 * @return The new 'this' as a result of substraction.
		 */
		template<typename Type2>
		Vec2<Type>& operator -=(const Vec2<Type2> &v) {
			vec[0] = static_cast<Type>(vec[0]-v[0]);
			vec[1] = static_cast<Type>(vec[1]-v[1]);
			return *this;
		}

		/** 'this' -= d; Minus a number from each item of 'this' vector.
		 * @param d The number used to be substracted from 'this' vector.
		 * @return The new 'this' as a result of substraction.
		 */
		template<typename Type2>
		Vec2<Type>& operator -=(const Type2& d) {
			vec[0] = static_cast<Type>(vec[0]-d);
			vec[1] = static_cast<Type>(vec[1]-d);
			return *this;
		}

		/** 'this' *= d; Multiply a number on each item of 'this' vector.
		 * @param d The number to multiply.
		 * @return The new 'this' as a result of multiplication.
		 */
		template<typename Type2>
		Vec2<Type>& operator *=(const Type2& d) {
			vec[0] = static_cast<Type>(vec[0]*d);
			vec[1] = static_cast<Type>(vec[1]*d);
			return *this;
		}

		/** 'this' /= d; Divide a number on each item of 'this' vector.
		 * @param d The number to divide.
		 * @return The new 'this' as a result of division.
		 */
		template<typename Type2>
		Vec2<Type>& operator /=(const Type2& d) {
			vec[0] = static_cast<Type>(vec[0]/d);
			vec[1] = static_cast<Type>(vec[1]/d);
			return *this;
		}


		private:
			Type vec[2];
	};

	template<typename Type,typename Type2>
	inline Vec2<Type> operator +(const Vec2<Type> &v1, const Vec2<Type2> &v2)
	{
		return Vec2<Type>(static_cast<Type>(v1[0] + v2[0]), static_cast<Type>(v1[1] + v2[1]));;
	}

	template<typename Type,typename Type2>
	inline Vec2<Type> operator +(const Vec2<Type> &v, const Type2& n)
	{
		Vec2<Type> v1(v);
		v1 += n;
		return v1;
	}

	template<typename Type,typename Type2>
	inline Vec2<Type> operator -(const Vec2<Type> &v1, const Vec2<Type2> &v2)
	{
		return Vec2<Type>(static_cast<Type>(v1[0] - v2[0]),	static_cast<Type>(v1[1] - v2[1]));
	}

	template<typename Type,typename Type2>
	inline Vec2<Type> operator -(const Vec2<Type> &v, const Type2& n)
	{
		Vec2<Type> v1(v);
		v1 -= n;
		return v1;
	}

	template<typename Type>
	inline Vec2<Type> operator -(const Vec2<Type> &v)
	{
		return Vec2<Type>(-v[0],-v[1]);
	}


	template<typename Type,typename Type2>
	inline Type operator *(const Vec2<Type> &v1, const Vec2<Type2> &v2)
	{
		return v1.dot(v2);
	}

	template<typename Type,typename Type2>
	inline Vec2<Type2> operator *(const Type& d, const Vec2<Type2> & v)
	{
		// Preserve the vector type
		Vec2<Type2> v1(v);
		v1 *= d;
		return v1;
	}

	template<typename Type,typename Type2>
	inline Vec2<Type> operator *(const Vec2<Type> & v,const Type2& d) {
	// Preserve the vector type
		Vec2<Type> v1(v);
		v1 *= d;
		return v1;
	}

	template<typename Type,typename Type2>
	inline Vec2<Type2> operator /(const Type& d, const Vec2<Type2> & v)
	{
	// Preserve the vector type
		Vec2<Type2> v1(v);
		v1 /= d;
		return v1;
	}

	template<typename Type,typename Type2>
	inline Vec2<Type> operator /(const Vec2<Type> & v,const Type2& d) {
		// Preserve the vector type
		Vec2<Type> v1(v);
		v1 /= d;
		return v1;
	}

	template<typename Type,typename Type2>
	inline bool operator ==(const Vec2<Type> &v1, const Vec2<Type2> &v2) {
		if (v1[0] == v2[0] && v1[1] == v2[1] ) {
			return true;
		}
		return false;
	}

	template<typename Type,typename Type2>
	inline bool operator !=(const Vec2<Type> &v1, const Vec2<Type2> &v2) {
		if (v1[0] != v2[0] || v1[1] != v2[1] ) {
			return true;
		}
		return false;
	}

	typedef Vec2<float> Vec2f;
	typedef Vec2<int> Vec2i;
	typedef Vec2<double> Vec2d;
}
#endif
