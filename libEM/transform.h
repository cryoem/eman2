#ifndef eman__transform_h__
#define eman__transform_h__ 1

#include <gsl/gsl_linalg.h>
#include <string>
#include <vector>
#include <math.h>
#include "geometry.h"

using std::vector;
using std::string;

namespace EMAN
{
    class Vec3f
    {
    public:
	Vec3f()
	{
	    vec[0] = 0;
	    vec[1] = 0;
	    vec[2] = 0;
	}

	Vec3f(float x, float y, float z)
	{
	    vec[0] = x;
	    vec[1] = y;
	    vec[2] = z;
	}

	Vec3f(const vector<float> &v)
	{
	    vec[0] = v[0];
	    vec[1] = v[1];
	    vec[2] = v[2];
	}
	
	virtual ~Vec3f()
	{
	}

	float normalize();
	float length() const;
	
	float dot(const Vec3f & v) const
	{
	    return (vec[0] * v[0] + vec[1] * v[1] + vec[2] * v[2]);
	}

	Vec3f cross(const Vec3f & v) const
	{
	    return Vec3f(vec[1] * v.vec[2] - vec[2] * v.vec[1],
			 vec[2] * v.vec[0] - vec[0] * v.vec[2],
			 vec[0] * v.vec[1] - vec[1] * v.vec[0]);
	}

	Vec3f & negate()
	{
	    vec[0] = -vec[0];
	    vec[1] = -vec[1];
	    vec[2] = -vec[2];
	    return (*this);
	}

	vector<float> get_value() const
	{
	    vector<float> v(3);
	    v[0] = vec[0];
	    v[1] = vec[1];
	    v[2] = vec[2];
	    return v;
	}

	void set_value(const vector<float> &v)
	{
	    vec[0] = v[0];
	    vec[1] = v[1];
	    vec[2] = v[2];
	}
	void set_value(float x, float y, float z)
	{
	    vec[0] = x;
	    vec[1] = y;
	    vec[2] = z;
	}

	float operator[] (int i) const
	{
	    return vec[i];
	}
	float &operator[] (int i)
	{
	    return vec[i];
	}

	Vec3f & operator +=(const Vec3f & v)
	{
	    vec[0] += v[0];
	    vec[1] += v[1];
	    vec[2] += v[2];
	    return *this;
	}
	
	Vec3f & operator -=(const Vec3f & v)
	{
	    vec[0] -= v[0];
	    vec[1] -= v[1];
	    vec[2] -= v[2];
	    return *this;
	}

	Vec3f & operator *=(float d)
	{
	    vec[0] *= d;
	    vec[1] *= d;
	    vec[2] *= d;
	    return *this;
	}
	
	Vec3f & operator /=(float d)
	{
	    if (d != 0) {
		vec[0] /= d;
		vec[1] /= d;
		vec[2] /= d;
	    }
	    return *this;
	}

	friend Vec3f operator +(const Vec3f & v1, const Vec3f & v2);
	friend Vec3f operator -(const Vec3f & v1, const Vec3f & v2);

	friend Vec3f operator *(float d, const Vec3f v);
	friend Vec3f operator /(float d, const Vec3f v);

	friend Vec3f operator *(const Vec3f v, float d);
	friend Vec3f operator /(const Vec3f v, float d);

	friend bool operator ==(const Vec3f & v1, const Vec3f & v2);
	friend bool operator !=(const Vec3f & v1, const Vec3f & v2);

    private:
	float vec[3];
    };


    class Matrix3f
    {
    public:
	Matrix3f();
	Matrix3f(float m0, float m1, float m2, float m3, float m4, float m5, float m6, float m7,
		 float m8);
	Matrix3f(const vector<float> &m);
	virtual ~Matrix3f();

	void make_identity();

	Matrix3f & mult_right(const Matrix3f & m);
	Matrix3f & mult_left(const Matrix3f & m);

	void set_value(const vector<float> &m);
	vector<float> get_value() const;

	Matrix3f & inverse();
	Matrix3f & transpose();
	Matrix3f create_inverse() const;

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

	friend Matrix3f operator+(float f, const Matrix3f & m2);
	friend Matrix3f operator-(float f, const Matrix3f & m2);
	friend Matrix3f operator*(float f, const Matrix3f & m2);
	friend Matrix3f operator/(float f, const Matrix3f & m2);

	friend Matrix3f operator+(const Matrix3f & m1, float f);
	friend Matrix3f operator-(const Matrix3f & m1, float f);
	friend Matrix3f operator*(const Matrix3f & m1, float f);
	friend Matrix3f operator/(const Matrix3f & m1, float f);

	friend Matrix3f operator+(const Matrix3f & m1, const Matrix3f & m2);
	friend Matrix3f operator-(const Matrix3f & m1, const Matrix3f & m2);
	friend Matrix3f operator*(const Matrix3f & m1, const Matrix3f & m2);
	friend Matrix3f operator/(const Matrix3f & m1, const Matrix3f & m2);

	friend bool operator==(const Matrix3f & m1, const Matrix3f & m2);
	friend bool operator!=(const Matrix3f & m1, const Matrix3f & m2);

	friend Vec3f operator*(const Vec3f & v, const Matrix3f & m2);
	friend Vec3f operator*(const Matrix3f & m1, const Vec3f & v);

    private:
	Matrix3f(gsl_matrix * m);
	gsl_matrix *get_gsl_matrix() const;

	gsl_matrix *matrix;
    };


    class Matrix4f
    {
    public:
	Matrix4f();
	Matrix4f(const vector<float> &m);
	Matrix4f(const Matrix3f & m3);
	virtual ~Matrix4f();

	Matrix4f & mult_right(const Matrix4f & m);
	Matrix4f & mult_left(const Matrix4f & m);

	void make_identity();

	void set_value(const vector<float> &m);
	vector<float> get_value() const;
	Matrix3f get_matrix3() const;

	Matrix4f & inverse();
	Matrix4f & transpose();
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

	friend Matrix4f operator+(float f, const Matrix4f & m2);
	friend Matrix4f operator-(float f, const Matrix4f & m2);
	friend Matrix4f operator*(float f, const Matrix4f & m2);
	friend Matrix4f operator/(float f, const Matrix4f & m2);

	friend Matrix4f operator+(const Matrix4f & m1, float f);
	friend Matrix4f operator-(const Matrix4f & m1, float f);
	friend Matrix4f operator*(const Matrix4f & m1, float f);
	friend Matrix4f operator/(const Matrix4f & m1, float f);

	friend Matrix4f operator+(const Matrix4f & m1, const Matrix4f & m2);
	friend Matrix4f operator-(const Matrix4f & m1, const Matrix4f & m2);
	friend Matrix4f operator*(const Matrix4f & m1, const Matrix4f & m2);
	friend Matrix4f operator/(const Matrix4f & m1, const Matrix4f & m2);

	friend bool operator==(const Matrix4f & m1, const Matrix4f & m2);
	friend bool operator!=(const Matrix4f & m1, const Matrix4f & m2);

    private:
	Matrix4f(gsl_matrix * m);
	gsl_matrix *get_gsl_matrix() const;
	gsl_matrix *matrix;
    };

    /*
      Euler angles have the disadvantage of being
      susceptible to "Gimbal lock" where attempts to rotate an
      object fail due to the order in which the rotations are performed.

      Quaternions are a solution to this problem. Instead of rotating an
      object through a series of successive rotations, a quaternion allows
      the programmer to rotate an object through a single arbitary rotation
      axis.

      Because the rotation axis is specifed as a unit direction vector,
      it may be calculated through vector mathematics or from spherical
      coordinates ie (longitude/latitude).

      Quaternions offer another advantage in that they be interpolated.
      This allows for smooth and predictable rotation effects.
    */

    class Quaternion
    {
    public:
	Quaternion();
	Quaternion(float e0, float e1, float e2, float e3);
	Quaternion(float radians, const Vec3f & axis);
	Quaternion(const Vec3f & axis, float radians);
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

	Vec3f rotate(const Vec3f & v) const;

	void interpolate(const Quaternion & q1, float alpha);
	void interpolate(const Quaternion & q1, const Quaternion & q2, float alpha);

	float to_angle() const;
	Vec3f to_axis() const;

	Matrix3f to_matrix3() const;
	Matrix4f to_matrix4() const;

	float real() const;
	Vec3f unreal() const;

	vector<float> get_value() const;

	Quaternion & operator+=(const Quaternion & q);
	Quaternion & operator-=(const Quaternion & q);
	Quaternion & operator*=(const Quaternion & q);
	Quaternion & operator*=(float s);
	Quaternion & operator/=(const Quaternion & q);
	Quaternion & operator/=(float s);

	friend Quaternion operator+(const Quaternion & q1, const Quaternion & q2);
	friend Quaternion operator-(const Quaternion & q1, const Quaternion & q2);

	friend Quaternion operator*(const Quaternion & q1, const Quaternion & q2);
	friend Quaternion operator*(const Quaternion & q, float s);
	friend Quaternion operator*(float s, const Quaternion & q);
	friend Quaternion operator/(const Quaternion & q1, const Quaternion & q2);

	friend bool operator==(const Quaternion & q1, const Quaternion & q2);
	friend bool operator!=(const Quaternion & q1, const Quaternion & q2);

    private:
	float e0;
	float e1;
	float e2;
	float e3;
    };


    class Rotation
    {
    public:
	static const float ERR_LIMIT;
	enum Type { EMAN, IMAGIC, SPIN, QUATERNION, MATRIX, SGIROT, SPIDER, MRC, UNKNOWN };

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
	float spider_psi() const;

	vector<float> get_spin_axis() const;
	void get_spin_axis(float *q, float *n1, float *n2, float *n3) const;

	vector<float> get_sgi() const;
	void get_sgi(float *q, float *n1, float *n2, float *n3) const;

	Quaternion get_quaternion() const;
	Matrix3f get_matrix3() const;
	Matrix4f get_matrix4() const;

	Rotation & operator*=(const Rotation & e);
	Rotation & operator/=(const Rotation & e);

	friend Rotation operator*(const Rotation & e1, const Rotation & e2);
	friend Rotation operator/(const Rotation & e1, const Rotation & e2);

	friend bool operator==(const Rotation & e1, const Rotation & e2);
	friend bool operator!=(const Rotation & e1, const Rotation & e2);

    private:
	enum SymType { CSYM, DSYM, ICOS_SYM, OCT_SYM, ISYM, UNKNOWN_SYM };
	SymType get_sym_type(string symname);

    private:
	Matrix3f matrix;
	Type type;
	string symname;
	int nsym;
	int cur_sym;


    };

    // composite transform: = (-Translate) * Scale * Rotate * Translate
    // Translate = center
    class Transform
    {
    public:
	enum TransformType {
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
	}

	Transform(const Matrix4f & m)
	{
	    matrix = m;
	}

	Transform(const Rotation & rotation)
	{
	    matrix = rotation.get_matrix4();
	}

	Transform(const Rotation & rotation, const Vec3f & post_translation)
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

	virtual ~Transform() {
	}

	Transform & set_rotate_instance(const Rotation & r)
	{
	    matrix = r.get_matrix4();
	    return (*this);
	}

	Transform & set_rotate_instance(const Matrix3f & m)
	{
	    matrix = Matrix4f(m);
	    return (*this);
	}

	Transform & set_translate_instance(const Vec3f & s)
	{
	    matrix.make_identity();
	    return set_post_translate(s);
	}

	Transform & set_scale_instance(const Vec3f & s)
	{
	    matrix.make_identity();
	    for (int i = 0; i < 3; i++) {
		matrix[i][i] = s[i];
	    }
	    return (*this);
	}

	Transform & set_transform_instance(const Vec3f & translation, const Rotation & roration,
					   const Vec3f & scale_factor,
					   const Rotation & scale_orientation,
					   const Vec3f & center);

	Transform & set_transform_instance(const Vec3f & translation,
					   const Rotation & rotation, const Vec3f & scale_factor)
	{
	    return set_transform_instance(translation, rotation, scale_factor,
					  Rotation(1, 0, 0, 0, Rotation::QUATERNION),
					  Vec3f(0, 0, 0));
	}

	Transform & set_center(const Vec3f & c)
	{
	    pre_trans = -1.0 * c;

	    for (int i = 0; i < 3; i++) {
		matrix[i][3] += c[i];
	    }
	    return (*this);
	}


	Transform & set_matrix(const Matrix4f & m)
	{
	    matrix = m;
	    return (*this);
	}

	Transform & set_post_translate(const Vec3f & s)
	{
	    for (int i = 0; i < 3; i++) {
		matrix[3][i] = s[i];
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

	// [this] = [this] x [t]
	Transform & post_concatenate(const Transform & t)
	{
	    (*this) = (*this) * t;
	    return (*this);
	}

	// [this] = [t] x [this]
	Transform & pre_concatenate(const Transform & t)
	{
	    (*this) = t * (*this);
	    return (*this);
	}

	// Concatenates this transform with a translation transformation.
	Transform & translate(const Vec3f & v)
	{
	    if (v != Vec3f(0, 0, 0)) {
		Matrix4f m;
		for (int i = 0; i < 3; i++) {
		    m[3][i] = v[i];
		}
		matrix *= m;
	    }
	    return (*this);
	}


	// Concatenates this transform with a rotation transformation.
	Transform & rotate(const Rotation & r)
	{
	    if (r != Rotation(1, 0, 0, 0, Rotation::QUATERNION))
		{
		    Matrix4f m = r.get_matrix4();
		    (*this) = (*this) * m;
		}
	    return (*this);
	}

	Transform & rotate_center(const Rotation & r, const Vec3f center)
	{
	    translate(-1.0 * center);
	    rotate(r);
	    translate(center);
	    return (*this);
	}

	Transform & rotate_scale(const Rotation & rotation,
				 const Vec3f & scale_factor, const Vec3f & center)
	{
	    scale(scale_factor);
	    rotate_center(rotation, center);
	    return (*this);
	}

	Transform & pre_translate_rotate(const Vec3f t, const Rotation & r)
	{
	    translate(t);
	    rotate(r);
	    return (*this);
	}

	Transform & post_translate_rotate(const Rotation & r, const Vec3f t)
	{
	    rotate(r);
	    translate(t);
	    return (*this);
	}

	// Concatenates this transform with a scaling transformation.
	Transform & scale(const Vec3f & scale_factor)
	{
	    if (scale_factor != Vec3f(1, 1, 1)) {
		Matrix4f m;
		m.make_identity();
		for (int i = 0; i < 3; i++) {
		    m[i][i] = scale_factor[i];
		}
		matrix *= m;
	    }
	    return (*this);
	}


	Vec3f transform(const Vec3f & v);
	Vec3f inverse_transform(const Vec3f & v);

	Rotation get_rotation(Rotation & r) const
	{
	    return Rotation(matrix.get_matrix3());
	}

	Vec3f get_scale() const
	{
	    return Vec3f(matrix[0][0], matrix[1][1], matrix[2][2]);
	}
	
	Vec3f get_center() const
	{
	    return pre_trans;
	}

	Matrix4f get_matrix() const
	{
	    return matrix;
	}

	Vec3f get_pre_translate() const
	{
	    return pre_trans;
	}
	
	Vec3f get_post_translate() const
	{
	    return Vec3f(matrix[3][0], matrix[3][1], matrix[3][2]);
	}

	int get_type() const;

	Transform & operator+=(const Transform & t)
	{
	    matrix += t.get_matrix();
	    return (*this);
	}

	Transform & operator-=(const Transform & t)
	{
	    matrix -= t.get_matrix();
	    return (*this);
	}
	
	Transform & operator*=(const Transform & t)
	{
	    matrix *= t.get_matrix();
	    return (*this);
	}
	
	Transform & operator*=(float scalar)
	{
	    matrix *= scalar;
	    return (*this);
	}
	
	Transform & operator/=(const Transform & t)
	{
	    matrix /= t.get_matrix();
	    return (*this);
	}
	
	Transform & operator/=(float scalar)
	{
	    matrix /= scalar;
	    return (*this);
	}

	friend Transform operator+(const Transform & t1, const Transform & t2);
	friend Transform operator-(const Transform & t1, const Transform & t2);

	friend Transform operator*(const Transform & t1, const Transform & t2);
	friend Transform operator*(const Transform & t, float s);
	friend Transform operator*(float s, const Transform & t);

	friend Transform operator/(const Transform & t1, const Transform & t2);
	friend Transform operator/(const Transform & t, float s);
	friend Transform operator/(float s, const Transform & t);

	static Transform interpolate(const Transform & t1, const Transform & t2, float percent);

    private:
	Matrix4f matrix;
	Vec3f pre_trans;
    };




}


#endif
