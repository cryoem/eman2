/**
 * $Id$
 */
#ifndef eman__transform_h__
#define eman__transform_h__ 1

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <string>
#include <vector>
#include <math.h>
#include "geometry.h"

using std::vector;
using std::string;

namespace EMAN
{
	/** Template Vec3 defines a 3-element vector and various vector
     * operations. The vector may store any numeric data type
     * including int, float, double, etc.
     *
     * For example, to define a 3-element float vector: Vec3<float> v;
     */
	template < class T > class Vec3
	{
	  public:
		Vec3() {
			vec[0] = 0;
			vec[1] = 0;
			vec[2] = 0;
		}

		Vec3(T x, T y, T z)
		{
			vec[0] = (T)x;
			vec[1] = (T)y;
			vec[2] = (T)z;
		}

		Vec3(const vector < T > &v)
		{
			vec[0] = v[0];
			vec[1] = v[1];
			vec[2] = v[2];
		}

		virtual ~ Vec3() {
		}

		T normalize()
		{
			T len = length();
			if (len != 0) {
				(*this) *= (1 / len);
			}
			else {
				set_value(0, 0, 0);
			}
			return len;
		}

		T length() const
		{
			T t = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
			return (T) sqrt((float) t);
		}

		T dot(const Vec3 < T > &v) const
		{
			return (vec[0] * v[0] + vec[1] * v[1] + vec[2] * v[2]);
		}

		Vec3 < T > cross(const Vec3 < T > &v) const
		{
			return Vec3(vec[1] * v.vec[2] - vec[2] * v.vec[1],
						vec[2] * v.vec[0] - vec[0] * v.vec[2],
						vec[0] * v.vec[1] - vec[1] * v.vec[0]);
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

		void set_value(const vector < T > &v)
		{
			vec[0] = v[0];
			vec[1] = v[1];
			vec[2] = v[2];
		}

		void set_value(T x, T y, T z)
		{
			vec[0] = x;
			vec[1] = y;
			vec[2] = z;
		}

		T operator[] (int i) const
		{
			return vec[i];
		}
		T & operator[] (int i)
		{
			return vec[i];
		}

		Vec3 < T > &operator +=(const Vec3 < T > &v) {
			vec[0] += v[0];
			vec[1] += v[1];
			vec[2] += v[2];
			return *this;
		}

		Vec3 < T > &operator -=(const Vec3 < T > &v) {
			vec[0] -= v[0];
			vec[1] -= v[1];
			vec[2] -= v[2];
			return *this;
		}

		Vec3 < T > &operator *=(T d) {
			vec[0] *= d;
			vec[1] *= d;
			vec[2] *= d;
			return *this;
		}

		Vec3 < T > &operator /=(T d) {
			if (d != 0) {
				vec[0] /= d;
				vec[1] /= d;
				vec[2] /= d;
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


	template < class T > Vec3 < T > operator *(T d, const Vec3 < T > v)
	{
		Vec3 < T > v1 = v;
		v1 *= d;
		return v1;
	}

	template < class T > Vec3 < T > operator /(T d, const Vec3 < T > v) {
		Vec3 < T > v1 = v;
		if (d != 0) {
			Vec3 < T > v1 = v;
			v1 /= d;
		}
		return v1;
	}


	template < class T > Vec3 < T > operator *(const Vec3 < T > v, T d) {
		Vec3 < T > v1 = v;
		v1 *= d;
		return v1;
	}

	template < class T > Vec3 < T > operator /(const Vec3 < T > v, T d) {
		Vec3 < T > v1 = v;
		if (d != 0) {
			Vec3 < T > v1 = v;
			v1 /= d;
		}
		return v1;
	}

	template < class T > bool operator ==(const Vec3 < T > &v1, const Vec3 < T > &v2) {
		if (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2]) {
			return true;
		}
		return false;
	}

	template < class T > bool operator !=(const Vec3 < T > &v1, const Vec3 < T > &v2) {
		if (v1[0] != v2[0] || v1[1] != v2[1] || v1[2] != v2[2]) {
			return true;
		}
		return false;
	}

	/** Matrix3f defines a 3x3 floating-point matrix and various
     * matrix operation. It uses gsl matrix as its internal
     * representation and implementation.
     */
	class Matrix3f
	{
	  public:
		Matrix3f();
		Matrix3f(float m0, float m1, float m2, float m3, float m4, float m5, float m6, float m7,
				 float m8);
		Matrix3f(const vector < float >&m);
		Matrix3f(const Matrix3f & m);
		virtual ~ Matrix3f();

		Matrix3f & operator=(const Matrix3f & m);

		void make_identity();

		/** this = this * m
	 */
		Matrix3f & mult_right(const Matrix3f & m);

		/** this = m * this
		 */
		Matrix3f & mult_left(const Matrix3f & m);

		void set_value(const vector < float >&m);
		vector < float >get_as_list() const;

		/** inverse this matrix and return it.
		 */
		Matrix3f & inverse();
		Matrix3f & transpose();

		/** return the inverse of this matrix without changing this matrix.
		 */
		Matrix3f create_inverse() const;

		Vec3 < float >get_vector(int i) const;

		double *operator[] (int i);
		const double *operator[] (int i) const;

		Matrix3f & operator+=(float f);
		Matrix3f & operator-=(float f);
		Matrix3f & operator*=(float f);
		Matrix3f & operator/=(float f);

		Matrix3f & operator+=(const Matrix3f & m);
		Matrix3f & operator-=(const Matrix3f & m);
		Matrix3f & operator*=(const Matrix3f & m);
		Matrix3f & operator/=(const Matrix3f & m);

	  private:
		Matrix3f(gsl_matrix * m);
		gsl_matrix *get_gsl_matrix() const;

		gsl_matrix *matrix;
	};

	inline Matrix3f::Matrix3f()
	{
		matrix = gsl_matrix_alloc(3, 3);
		gsl_matrix_set_identity(matrix);
	}

	inline Matrix3f & Matrix3f::operator+=(float f)
	{
		gsl_matrix_add_constant(matrix, f);
		return *this;
	}

	inline Matrix3f & Matrix3f::operator-=(float f)
	{
		gsl_matrix_add_constant(matrix, -f);
		return *this;
	}

	inline Matrix3f & Matrix3f::operator*=(float scale)
	{
		gsl_matrix_scale(matrix, scale);
		return *this;
	}

	inline Matrix3f & Matrix3f::operator/=(float scale)
	{
		if (scale != 0) {
			gsl_matrix_scale(matrix, scale);
		}
		return *this;
	}


	inline Matrix3f & Matrix3f::operator+=(const Matrix3f & m)
	{
		gsl_matrix_add(matrix, m.get_gsl_matrix());
		return *this;
	}

	inline Matrix3f & Matrix3f::operator-=(const Matrix3f & m)
	{
		gsl_matrix_sub(matrix, m.get_gsl_matrix());
		return *this;
	}

	inline Matrix3f & Matrix3f::operator*=(const Matrix3f & m)
	{
		int n=3;
		gsl_matrix* temp = gsl_matrix_alloc(n, n);
		gsl_matrix_memcpy(temp, matrix);
		gsl_matrix_set_zero(matrix);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, 
					   temp, m.get_gsl_matrix(), 0.0, matrix);
		gsl_matrix_free(temp);
		temp = 0;
		return *this;
	}

	inline Matrix3f & Matrix3f::operator/=(const Matrix3f & m)
	{
		Matrix3f m_inverse = m.create_inverse();
		(*this) *= m_inverse;
		return *this;
	}

	inline double *Matrix3f::operator[] (int i)
	{
		return gsl_matrix_ptr(matrix, i, 0);
	}

	inline const double *Matrix3f::operator[] (int i) const
	{
		return gsl_matrix_const_ptr(matrix, i, 0);
	}

	Matrix3f operator+(float f, const Matrix3f & m2);
	Matrix3f operator-(float f, const Matrix3f & m2);
	Matrix3f operator*(float f, const Matrix3f & m2);
	Matrix3f operator/(float f, const Matrix3f & m2);

	Matrix3f operator+(const Matrix3f & m1, float f);
	Matrix3f operator-(const Matrix3f & m1, float f);
	Matrix3f operator*(const Matrix3f & m1, float f);
	Matrix3f operator/(const Matrix3f & m1, float f);

	Matrix3f operator+(const Matrix3f & m1, const Matrix3f & m2);
	Matrix3f operator-(const Matrix3f & m1, const Matrix3f & m2);
	Matrix3f operator*(const Matrix3f & m1, const Matrix3f & m2);
	Matrix3f operator/(const Matrix3f & m1, const Matrix3f & m2);

	bool operator==(const Matrix3f & m1, const Matrix3f & m2);
	bool operator!=(const Matrix3f & m1, const Matrix3f & m2);

	Vec3 < float >operator*(const Vec3 < int >&v, const Matrix3f & m2);
	Vec3 < float >operator*(const Matrix3f & m1, const Vec3 < int >&v);

	Vec3 < float >operator*(const Vec3 < float >&v, const Matrix3f & m2);
	Vec3 < float >operator*(const Matrix3f & m1, const Vec3 < float >&v);



	/** Matrix4f defines a 4x4 floating-point matrix and various
     * matrix operation. It uses gsl matrix as its internal
     * representation and implementation.
     */
	class Matrix4f
	{
	  public:
		Matrix4f();
		Matrix4f(const vector < float >&m);
		Matrix4f(const Matrix3f & m3);
		Matrix4f(const Matrix4f & m);
		virtual ~ Matrix4f();

		Matrix4f & operator=(const Matrix4f & m);

		/** this = this * m
		 */
		Matrix4f & mult_right(const Matrix4f & m);

		/** this = m * this
		 */
		Matrix4f & mult_left(const Matrix4f & m);

		void make_identity();

		void set_value(const vector < float >&m);
		vector < float >get_as_list() const;
		Matrix3f get_matrix3() const;

		/** inverse this matrix and return it.
		 */
		Matrix4f & inverse();
		Matrix4f & transpose();

		/** return the inverse of this matrix without changing this matrix.
		 */
		Matrix4f create_inverse() const;

		Matrix4f & operator+=(float f);
		Matrix4f & operator-=(float f);
		Matrix4f & operator*=(float f);
		Matrix4f & operator/=(float f);

		Matrix4f & operator+=(const Matrix4f & m);
		Matrix4f & operator-=(const Matrix4f & m);
		Matrix4f & operator*=(const Matrix4f & m);
		Matrix4f & operator/=(const Matrix4f & m);

		double *operator[] (int i);
		const double *operator[] (int i) const;

	  private:
		Matrix4f(gsl_matrix * m);
		gsl_matrix *get_gsl_matrix() const;
		gsl_matrix *matrix;
	};


	Matrix4f operator+(float f, const Matrix4f & m2);
	Matrix4f operator-(float f, const Matrix4f & m2);
	Matrix4f operator*(float f, const Matrix4f & m2);
	Matrix4f operator/(float f, const Matrix4f & m2);

	Matrix4f operator+(const Matrix4f & m1, float f);
	Matrix4f operator-(const Matrix4f & m1, float f);
	Matrix4f operator*(const Matrix4f & m1, float f);
	Matrix4f operator/(const Matrix4f & m1, float f);

	Matrix4f operator+(const Matrix4f & m1, const Matrix4f & m2);
	Matrix4f operator-(const Matrix4f & m1, const Matrix4f & m2);
	Matrix4f operator*(const Matrix4f & m1, const Matrix4f & m2);
	Matrix4f operator/(const Matrix4f & m1, const Matrix4f & m2);

	bool operator==(const Matrix4f & m1, const Matrix4f & m2);
	bool operator!=(const Matrix4f & m1, const Matrix4f & m2);


	/** Quaternion is used in Rotation and Transformation to replace Euler angles.
     *
     * Quaternions extend the concept of rotation in three dimensions to
     * rotation in four dimensions. This avoids the problem of "gimbal-lock"
     * and allows for the implementation of smooth and continuous rotation.
     *
     * Euler angles have the disadvantage of being
     * susceptible to "Gimbal lock" where attempts to rotate an
     * object fail due to the order in which the rotations are performed.
     * 
     * Quaternions are a solution to this problem. Instead of rotating an
     * object through a series of successive rotations, a quaternion allows
     * the programmer to rotate an object through a single arbitary rotation
     * axis.
     * 
     * Because the rotation axis is specifed as a unit direction vector,
     * it may be calculated through vector mathematics or from spherical
     * coordinates ie (longitude/latitude).
     * 
     * Quaternions offer another advantage in that they be interpolated.
     * This allows for smooth and predictable rotation effects.
     */
	class Quaternion
	{
	  public:
		Quaternion();
		Quaternion(float e0, float e1, float e2, float e3);
		Quaternion(float radians, const Vec3 < float >&axis);
		Quaternion(const Vec3 < float >&axis, float radians);
		Quaternion(const Matrix3f & m);
		Quaternion(const Matrix4f & m);
		 ~Quaternion()
		{
		}

		float norm() const
		{
			return (e0 * e0 + e1 * e1 + e2 * e2 + e3 * e3);
		}

		Quaternion conj() const
		{
			return (Quaternion(e0, -e1, -e2, -e3));
		}

		float abs() const
		{
			return sqrt(norm());
		}

		void normalize();
		Quaternion & inverse();
		Quaternion create_inverse() const;

		Vec3 < float >rotate(const Vec3 < float >&v) const;

		float to_angle() const;
		Vec3 < float >to_axis() const;

		Matrix3f to_matrix3() const;
		Matrix4f to_matrix4() const;

		float real() const;
		Vec3 < float >unreal() const;

		vector < float >get_as_list() const;

		Quaternion & operator+=(const Quaternion & q);
		Quaternion & operator-=(const Quaternion & q);
		Quaternion & operator*=(const Quaternion & q);
		Quaternion & operator*=(float s);
		Quaternion & operator/=(const Quaternion & q);
		Quaternion & operator/=(float s);

		static Quaternion interpolate(const Quaternion & from, const Quaternion & to,
									  float percent);
	  private:
		float e0;
		float e1;
		float e2;
		float e3;
	};

	Quaternion operator+(const Quaternion & q1, const Quaternion & q2);
	Quaternion operator-(const Quaternion & q1, const Quaternion & q2);

	Quaternion operator*(const Quaternion & q1, const Quaternion & q2);
	Quaternion operator*(const Quaternion & q, float s);
	Quaternion operator*(float s, const Quaternion & q);
	Quaternion operator/(const Quaternion & q1, const Quaternion & q2);

	bool operator==(const Quaternion & q1, const Quaternion & q2);
	bool operator!=(const Quaternion & q1, const Quaternion & q2);

	/** Rotation class defines various conventions of Euler angles.
     * It can also convert any convention to any other convention.
     *
     * Any rotation may be described using three angles called Euler Angles.
     * There are several conventions for Euler angles, depending on
     * the axes about which the rotations are carried out.    
     *   - x-convention: (z, x, z')
     *   - y-convention
     *   - xyz-convention
     *
     * Currently the following conventions are supported: EMAN,
     * IMAGIC, SPIN,  QUATERNION, MATRIX, SGIROT, SPIDER, MRC.
     *
     * EMAN is a (z,x,z') -> (az, alt, phi), usually specified as
     * (alt, az, phi).
     */
	class Rotation
	{
	  public:
		static const float ERR_LIMIT;
		enum Type
		{ EMAN, IMAGIC, SPIN, QUATERNION, MATRIX, SGIROT, SPIDER, MRC, UNKNOWN };

	  public:
		Rotation();
		Rotation(float a1, float a2, float a3, Type type);
		Rotation(float e1, float e2, float e3, float e0, Type type);
		Rotation(const Quaternion & q);
		Rotation(const Matrix3f & m);

		 ~Rotation();

		Rotation & inverse();
		Rotation create_inverse() const;

		float diff(const Rotation & r);
		Rotation & rotate_from_left(const Rotation & left_rotate);

		void set_sym(string symname);
		int get_max_nsym() const;
		Rotation get_sym(int sym_index);

		void set_angle(float a1, float a2, float a3, Type type);
		void set_angle(float e1, float e2, float e3, float e0, Type type);
		void set_angle(const Matrix3f & m);
		void set_angle(const Quaternion & q);

		bool is_valid() const;
		void rectify();

		float eman_alt() const;
		float eman_az() const;
		float eman_phi() const;

		float mrc_theta() const;
		float mrc_phi() const;
		float mrc_omega() const;

		float imagic_alpha() const;
		float imagic_beta() const;
		float imagic_gamma() const;

		float spider_phi() const;
		float spider_theta() const;
		float spider_gamma() const;

		vector < float >get_spin_axis() const;
		void get_spin_axis(float *q, float *n1, float *n2, float *n3) const;

		vector < float >get_sgi() const;
		void get_sgi(float *q, float *n1, float *n2, float *n3) const;

		Quaternion get_quaternion() const;
		Matrix3f get_matrix3() const;
		Matrix4f get_matrix4() const;

		Rotation & operator*=(const Rotation & e);
		Rotation & operator/=(const Rotation & e);


		static Rotation interpolate(const Rotation & from, const Rotation & to, float percent);

	  private:
		enum SymType
		{ CSYM, DSYM, ICOS_SYM, OCT_SYM, ISYM, UNKNOWN_SYM };
		SymType get_sym_type(string symname);

	  private:
		Matrix3f matrix;
		Type type;
		string symname;
		int nsym;
		int cur_sym;
	};

	Rotation operator*(const Rotation & e1, const Rotation & e2);
	Rotation operator/(const Rotation & e1, const Rotation & e2);

	bool operator==(const Rotation & e1, const Rotation & e2);
	bool operator!=(const Rotation & e1, const Rotation & e2);


	/** Transform defines a transformation, which can be rotation,
     * translation, scale, and their combinations.
	 *
	 * Internally a transformation is stored in a 4x4 matrix. With the
	 * left-top 3x3 submatrix as the rotation part.
	 *
	 * Center is stored in (m[3][0], m[3][1], m[3][2])
	 * post translation is stored in (m[0][3], m[1][3], m[2][3])
	 * pre translation is stored out side of this matrix separately.
	 *
	 *   | R R R T |
	 *   | R R R T |
	 *   | R R R T |
	 *   | C C C 1 |
	 *
     */
	class Transform
	{
	  public:
		enum TransformType
		{
			ROTATION = 1 << 1,
			TRANSLATION = 1 << 2,
			SCALE = 1 << 3,
			UNIFORM_SCALE = 1 << 4,
			IDENTITY = 1 << 5,
			TRANSFORM = 1 << 6
		};

	  public:
		Transform()
		{
			matrix.make_identity();
		}

		Transform(const Matrix4f & m)
		{
			matrix = m;
		}

		Transform(const Vec3 < float >&translation)
		{
			for (int i = 0; i < 3; i++) {
				matrix[i][3] = translation[i];
			}
		}

		Transform(const Rotation & rotation)
		{
			matrix = rotation.get_matrix4();
		}

		Transform(const Rotation & rotation, const Vec3 < float >&post_translation)
		{
			Matrix3f rotate_m = rotation.get_matrix3();

			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					matrix[i][j] = rotate_m[i][j];
				}
			}

			for (int i = 0; i < 3; i++) {
				matrix[i][3] = post_translation[i];
			}
		}

		virtual ~ Transform() {
		}

		Transform & set_rotate_instance(const Rotation & r)
		{
			matrix = r.get_matrix4();
			return (*this);
		}

		Transform & set_rotate_instance(const Matrix3f & m) {
			matrix = Matrix4f(m);
			return (*this);
		}

		Transform & set_translate_instance(const Vec3 < float >&s) {
			matrix.make_identity();
			return set_post_translate(s);
		}

		Transform & set_scale_instance(const Vec3 < float >&s) {
			matrix.make_identity();
			for (int i = 0; i < 3; i++) {
				matrix[i][i] = s[i];
			}
			return (*this);
		}

		Transform & set_transform_instance(const Vec3 < float >&translation,
										   const Rotation & roration,
										   const Vec3 < float >&scale_factor,
										   const Rotation & scale_orientation,
										   const Vec3 < float >&center);

		Transform & set_transform_instance(const Vec3 < float >&translation,
										   const Rotation & rotation,
										   const Vec3 < float >&scale_factor)
		{
			return set_transform_instance(translation, rotation, scale_factor,
										  Rotation(1, 0, 0, 0, Rotation::QUATERNION),
										  Vec3 < float >(0, 0, 0));
		}

		Transform & set_center(const Vec3 < float >&c)
		{
			pre_trans = (float) -1.0 * c;

			for (int i = 0; i < 3; i++) {
				matrix[3][i] += c[i];
			}
			return (*this);
		}


		Transform & set_matrix(const Matrix4f & m) {
			matrix = m;
			return (*this);
		}

		Transform & set_post_translate(const Vec3 < float >&s) {
			for (int i = 0; i < 3; i++) {
				matrix[i][3] = s[i];
			}
			return (*this);
		}

		Transform & inverse() {
			matrix.inverse();
			return (*this);
		}

		Transform create_inverse()
		{
			return Transform(matrix.create_inverse());
		}

		Transform & transpose() {
			matrix.transpose();
			return (*this);
		}

	/** [this] = [this] x [t]
	 */
		Transform & post_concatenate(const Transform & t) {
			(*this) *= t;
			return (*this);
		}

	/** [this] = [t] x [this]
	 */
		Transform & pre_concatenate(const Transform & t) {
			Transform t1 = t;
			t1 *= (*this);
			(*this) = t1;
			return (*this);
		}

	/** Concatenates this transform with a translation transformation.
	*/
		Transform & translate(const Vec3 < float >&v) {
			if (v != Vec3 < float >(0, 0, 0)) {
				Matrix4f m;
				for (int i = 0; i < 3; i++) {
					m[i][3] = v[i];
				}
				matrix *= m;
			}
			return (*this);
		}


	/** Concatenates this transform with a rotation transformation.
	 */
		Transform & rotate(const Rotation & r) {
			if (r != Rotation(1, 0, 0, 0, Rotation::QUATERNION)) {
				Matrix4f m = r.get_matrix4();
				(*this) *= m;
			}
			return (*this);
		}

		Transform & rotate_center(const Rotation & r, const Vec3 < float >center) {
			translate((float) -1.0 * center);
			rotate(r);
			translate(center);
			return (*this);
		}

		Transform & rotate_scale(const Rotation & rotation,
								 const Vec3 < float >&scale_factor, const Vec3 < float >&center) {
			scale(scale_factor);
			rotate_center(rotation, center);
			return (*this);
		}

		Transform & pre_translate_rotate(const Vec3 < float >t, const Rotation & r) {
			translate(t);
			rotate(r);
			return (*this);
		}

		Transform & post_translate_rotate(const Rotation & r, const Vec3 < float >t) {
			rotate(r);
			translate(t);
			return (*this);
		}

	/** Concatenates this transform with a scaling transformation.
	 */
		Transform & scale(const Vec3 < float >&scale_factor) {
			if (scale_factor != Vec3 < float >(1, 1, 1)) {
				Matrix4f m;
				m.make_identity();
				for (int i = 0; i < 3; i++) {
					m[i][i] = scale_factor[i];
				}
				matrix *= m;
			}
			return (*this);
		}


		Vec3 < float >transform(const Vec3 < float >&v);
		Vec3 < float >inverse_transform(const Vec3 < float >&v);

		Rotation get_rotation() const
		{
			return Rotation(matrix.get_matrix3());
		}

		float get_scale(int i) const;		
		Vec3 < float >get_scale() const;
		
		Vec3 < float >get_center() const
		{
			return pre_trans;
		}

		Matrix4f get_matrix() const
		{
			return matrix;
		}

		Vec3 < float >get_pre_translate() const
		{
			return pre_trans;
		}

		Vec3 < float >get_post_translate() const
		{
			return Vec3 < float >(matrix[0][3], matrix[1][3], matrix[2][3]);
		}

		int get_type() const;

		Transform & operator+=(const Transform & t)
		{
			matrix += t.get_matrix();
			return (*this);
		}

		Transform & operator-=(const Transform & t) {
			matrix -= t.get_matrix();
			return (*this);
		}

		Transform & operator*=(const Transform & t) {
			matrix *= t.get_matrix();
			return (*this);
		}

		Transform & operator*=(float scalar) {
			matrix *= scalar;
			return (*this);
		}

		Transform & operator/=(const Transform & t) {
			matrix /= t.get_matrix();
			return (*this);
		}

		Transform & operator/=(float scalar) {
			matrix /= scalar;
			return (*this);
		}

	  private:
		Matrix4f matrix;
		Vec3 < float >pre_trans;
	};

	Transform operator+(const Transform & t1, const Transform & t2);
	Transform operator-(const Transform & t1, const Transform & t2);

	Transform operator*(const Transform & t1, const Transform & t2);
	Transform operator*(const Transform & t, float s);
	Transform operator*(float s, const Transform & t);

	Transform operator/(const Transform & t1, const Transform & t2);
	Transform operator/(const Transform & t, float s);
	Transform operator/(float s, const Transform & t);

}


#endif
