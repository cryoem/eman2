/**
 * $Id$
 */
#include "quaternion.h"


using namespace EMAN;

Quaternion::Quaternion()
:	e0(0), e1(0), e2(0), e3(0)
{
}

Quaternion::Quaternion(float radians, const Vec3 < float >&axis)
{
	Vec3 < float >q = axis;
	//normalize();

	q *= sin(radians / 2.0f);
	e0 = cos(radians / 2.0f);

	e1 = q[0];
	e2 = q[1];
	e3 = q[2];
}

Quaternion::Quaternion(const Vec3 < float >&axis, float radians)
{
	Vec3 < float >q = axis;
	//normalize();

	q *= sin(radians / 2.0f);
	e0 = cos(radians / 2.0f);

	e1 = q[0];
	e2 = q[1];
	e3 = q[2];
}

Quaternion::Quaternion(float ee0, float ee1, float ee2, float ee3)
:e0(ee0), e1(ee1), e2(ee2), e3(ee3)
{
	//normalize();
}

Quaternion::Quaternion(const Matrix3f & m)
{
	int i = 0;

	if (m[0][0] > m[1][1]) {
		if (m[0][0] > m[2][2]) {
			i = 0;
		}
		else {
			i = 2;
		}
	}
	else {
		if (m[1][1] > m[2][2]) {
			i = 1;
		}
		else {
			i = 2;
		}
	}

	if (m[0][0] + m[1][1] + m[2][2] > m[i][i]) {
		e0 = sqrt(m[0][0] + m[1][1] + m[2][2] + 1) / 2.0f;
		e1 = (m[1][2] - m[2][1]) / (4 * e0);
		e2 = (m[2][0] - m[0][2]) / (4 * e0);
		e3 = (m[0][1] - m[1][0]) / (4 * e0);
	}
	else {
		float quat[3];
		int j = (i + 1) % 3;
		int k = (i + 2) % 3;

		quat[i] = sqrt(m[i][i] - m[j][j] - m[k][k] + 1) / 2.0f;
		quat[j] = (m[i][j] + m[j][i]) / (4 * quat[i]);
		quat[k] = (m[i][k] + m[k][i]) / (4 * quat[i]);

		e0 = (m[j][k] - m[k][j]) / (4 * quat[i]);
		e1 = quat[0];
		e2 = quat[1];
		e3 = quat[2];
	}

	//normalize();
}


Quaternion::Quaternion(const Matrix4f & m)
{
	int i = 0;

	if (m[0][0] > m[1][1]) {
		if (m[0][0] > m[2][2]) {
			i = 0;
		}
		else {
			i = 2;
		}
	}
	else {
		if (m[1][1] > m[2][2]) {
			i = 1;
		}
		else {
			i = 2;
		}
	}

	if (m[0][0] + m[1][1] + m[2][2] > m[i][i]) {
		e0 = sqrt(m[0][0] + m[1][1] + m[2][2] + m[3][3]) / 2.0f;
		e1 = (m[1][2] - m[2][1]) / (4 * e0);
		e2 = (m[2][0] - m[0][2]) / (4 * e0);
		e3 = (m[0][1] - m[1][0]) / (4 * e0);
	}
	else {
		float quat[3];
		int j = (i + 1) % 3;
		int k = (i + 2) % 3;

		quat[i] = sqrt(m[i][i] - m[j][j] - m[k][k] + m[3][3]) / 2.0f;
		quat[j] = (m[i][j] + m[j][i]) / (4 * quat[i]);
		quat[k] = (m[i][k] + m[k][i]) / (4 * quat[i]);

		e0 = (m[j][k] - m[k][j]) / (4 * quat[i]);
		e1 = quat[0];
		e2 = quat[1];
		e3 = quat[2];
	}

	//normalize();
}


void Quaternion::normalize()
{
	float dist = 1.0f / sqrt(norm());
	e0 *= dist;
	e1 *= dist;
	e2 *= dist;
	e3 *= dist;
}

Quaternion & Quaternion::inverse()
{
	float f = 1.0f / norm();
	e0 *= f;
	e1 *= -f;
	e2 *= -f;
	e3 *= -f;
	return (*this);
}

Quaternion Quaternion::create_inverse() const
{
	Quaternion q = *this;
	return q.inverse();
}


Vec3 < float >Quaternion::rotate(const Vec3 < float >&v) const
{
	Vec3 < float >i(e1, e2, e3);
	Vec3 < float >v1 = i.cross(v) * (2 * e0);
	Vec3 < float >v2 = v1.cross(i) * (float) 2;
	Vec3 < float >rotated = v + v1 - v2;

	return rotated;
}


float Quaternion::to_angle() const
{
	Vec3 < float >q(e1, e2, e3);
	float len = q.length();
	float radians = 0;

	if (len > 0.00001f) {
		radians = 2.0f * acos(e0);
	}
	else {
		radians = 0;
	}
	return radians;
}

Vec3 < float >Quaternion::to_axis() const
{
	Vec3 < float >q(e1, e2, e3);
	float len = q.length();
	Vec3 < float >axis;

	if (len > 0.00001f) {
		axis = q * ((float) (1.0f / len));
	}
	else {
		axis.set_value(0.0f, 0.0f, 1.0f);
	}
	return axis;
}


Matrix3f Quaternion::to_matrix3() const
{
	vector < float >m(9);

	m[0] = e0 * e0 + e1 * e1 - e2 * e2 - e3 * e3;
	m[1] = 2.0f * (e1 * e2 + e0 * e3);
	m[2] = 2.0f * (e1 * e3 - e0 * e2);

	m[3] = 2.0f * (e1 * e2 - e0 * e3);
	m[4] = e0 * e0 + e1 * e1 + e2 * e2 - e3 * e3;
	m[5] = 2.0f * (e2 * e3 + e0 * e1);

	m[6] = 2.0f * (e1 * e3 + e0 * e2);
	m[7] = 2.0f * (e2 * e3 - e0 * e1);
	m[8] = e0 * e0 - e1 * e1 - e2 * e2 + e3 * e3;

	return Matrix3f(m);
}

Matrix4f Quaternion::to_matrix4() const
{
	vector < float >m(16);

	m[0] = e0 * e0 + e1 * e1 - e2 * e2 - e3 * e3;
	m[1] = 2.0f * (e1 * e2 + e0 * e3);
	m[2] = 2.0f * (e1 * e3 - e0 * e2);
	m[3] = 0;

	m[4] = 2.0f * (e1 * e2 - e0 * e3);
	m[5] = e0 * e0 + e1 * e1 + e2 * e2 - e3 * e3;
	m[6] = 2.0f * (e2 * e3 + e0 * e1);
	m[7] = 0;

	m[8] = 2.0f * (e1 * e3 + e0 * e2);
	m[9] = 2.0f * (e2 * e3 - e0 * e1);
	m[10] = e0 * e0 - e1 * e1 - e2 * e2 + e3 * e3;
	m[11] = 0;

	m[12] = 0;
	m[13] = 0;
	m[14] = 0;
	m[15] = e0 * e0 + e1 * e1 + e2 * e2 + e3 * e3;

	return Matrix4f(m);
}



float Quaternion::real() const
{
	return e0;
}

Vec3 < float >Quaternion::unreal() const
{
	return Vec3 < float >(e1, e2, e3);
}

vector < float >Quaternion::get_as_list() const
{
	vector < float >v(4);
	v[0] = e0;
	v[1] = e1;
	v[2] = e2;
	v[3] = e3;

	return v;
}

Quaternion & Quaternion::operator+=(const Quaternion & q)
{
	e0 += q.e0;
	e1 += q.e1;
	e2 += q.e2;
	e3 += q.e3;
	return *this;
}

Quaternion & Quaternion::operator-=(const Quaternion & q)
{
	e0 -= q.e0;
	e1 -= q.e1;
	e2 -= q.e2;
	e3 -= q.e3;
	return *this;
}

Quaternion & Quaternion::operator*=(const Quaternion & q)
{
	float a = e0 * q.e0 - e1 * q.e1 - e2 * q.e2 - e3 * q.e3;
	float b = e0 * q.e1 + e1 * q.e0 + e2 * q.e3 - e3 * q.e2;
	float c = e0 * q.e2 - e1 * q.e3 + e2 * q.e0 + e3 * q.e1;
	float d = e0 * q.e3 + e1 * q.e2 - e2 * q.e1 + e3 * q.e0;

	e0 = a;
	e1 = b;
	e2 = c;
	e3 = d;


	// normalize();

	return (*this);
}

Quaternion & Quaternion::operator*=(float s)
{
	e0 *= s;
	e1 *= s;
	e2 *= s;
	e3 *= s;
	return (*this);
}

Quaternion & Quaternion::operator/=(const Quaternion & q)
{
	float qn = q.norm();

	float a = (+e0 * q.e0 + e1 * q.e1 + e2 * q.e2 + e3 * q.e3) / qn;
	float b = (-e0 * q.e1 + e1 * q.e0 - e2 * q.e3 + e3 * q.e2) / qn;
	float c = (-e0 * q.e2 + e1 * q.e3 + e2 * q.e0 - e3 * q.e1) / qn;
	float d = (-e0 * q.e3 - e1 * q.e2 + e2 * q.e1 + e3 * q.e0) / qn;

	e0 = a;
	e1 = b;
	e2 = c;
	e3 = d;

	return (*this);
}

Quaternion & Quaternion::operator/=(float s)
{
	if (s != 0) {
		e0 /= s;
		e1 /= s;
		e2 /= s;
		e3 /= s;
	}

	return (*this);
}


Quaternion EMAN::operator+(const Quaternion & q1, const Quaternion & q2)
{
	Quaternion q = q1;
	q += q2;
	return q;
}

Quaternion EMAN::operator-(const Quaternion & q1, const Quaternion & q2)
{
	Quaternion q = q1;
	q -= q2;
	return q;
}


Quaternion EMAN::operator*(const Quaternion & q1, const Quaternion & q2)
{
	Quaternion q = q1;
	q *= q2;
	return q;
}

Quaternion EMAN::operator*(const Quaternion & q, float s)
{
	Quaternion q1 = q;
	q1 *= s;
	return q1;
}

Quaternion EMAN::operator*(float s, const Quaternion & q)
{
	Quaternion q1 = q;
	q1 *= s;
	return q1;
}

Quaternion EMAN::operator/(const Quaternion & q1, const Quaternion & q2)
{
	Quaternion q = q1;
	q /= q2;
	return q;
}


bool EMAN::operator==(const Quaternion & q1, const Quaternion & q2)
{
	bool result = true;
	vector < float >v1 = q1.get_as_list();
	vector < float >v2 = q2.get_as_list();

	for (size_t i = 0; i < v1.size(); i++) {
		if (v1[i] != v2[i]) {
			result = false;
			break;
		}
	}

	return result;
}

bool EMAN::operator!=(const Quaternion & q1, const Quaternion & q2)
{
	return (!(q1 == q2));
}


Quaternion Quaternion::interpolate(const Quaternion & from, const Quaternion & to, float t)
{
	const double epsilon = 0.00001f;
	double cosom = from.e1 * to.e1 + from.e2 * to.e2 + from.e3 * to.e3 + from.e0 * to.e0;

	Quaternion q;
	if (cosom < 0.0f) {
		cosom = -cosom;
		q = q - to;
	}
	else {
		q = to;
	}

	double scale0 = 1 - t;
	double scale1 = t;

	if ((1 - cosom) > epsilon) {
		double omega = acos(cosom);
		double sinom = sin(omega);
		scale0 = sin((1 - t) * omega) / sinom;
		scale1 = sin(t * omega) / sinom;
	}

	return (scale0 * from + scale1 * q);
}
