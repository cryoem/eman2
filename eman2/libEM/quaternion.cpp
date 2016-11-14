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

#include "quaternion.h"

using namespace EMAN;

Quaternion::Quaternion()
:	e0(0), e1(0), e2(0), e3(0)
{
}

Quaternion::Quaternion(float radians, const Vec3f &axis)
{
	Vec3f q = axis;
	//normalize();

	q *= sin(radians / 2.0f);
	e0 = cos(radians / 2.0f);

	e1 = q[0];
	e2 = q[1];
	e3 = q[2];
}

Quaternion::Quaternion(const Vec3f &axis, float radians)
{
	Vec3f q = axis;
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

Quaternion::Quaternion(const vector<float> & m)
{
	int i = 0;

	if (m[0] > m[4]) {
		if (m[0] > m[8]) {
			i = 0;
		}
		else {
			i = 2;
		}
	}
	else {
		if (m[4] > m[8]) {
			i = 1;
		}
		else {
			i = 2;
		}
	}

	if (m[0] + m[4] + m[8] > m[i*3+i]) {
		e0 = (float) (sqrt(m[0] + m[4] + m[8] + 1) / 2.0);
		e1 = (float) ((m[5] - m[7]) / (4 * e0));
		e2 = (float) ((m[6] - m[2]) / (4 * e0));
		e3 = (float) ((m[1] - m[3]) / (4 * e0));
	}
	else {
		float quat[3];
		int j = (i + 1) % 3;
		int k = (i + 2) % 3;

		quat[i] = (float) (sqrt(m[i*3+i] - m[j*3+j] - m[k*3+k] + 1) / 2.0);
		quat[j] = (float) ((m[i*3+j] + m[j*3+i]) / (4 * quat[i]));
		quat[k] = (float) ((m[i*3+k] + m[k*3+i]) / (4 * quat[i]));

		e0 = (float) ((m[j*3+k] - m[k*3+j]) / (4 * quat[i]));
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


Vec3f Quaternion::rotate(const Vec3f &v) const
{
	Vec3f i(e1, e2, e3);
	Vec3f v1 = i.cross(v) * (2 * e0);
	Vec3f v2 = v1.cross(i) * (float) 2;
	Vec3f rotated = v + v1 - v2;

	return rotated;
}


float Quaternion::to_angle() const
{
	Vec3f q(e1, e2, e3);
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

Vec3f Quaternion::to_axis() const
{
	Vec3f q(e1, e2, e3);
	float len = q.length();
	Vec3f axis;

	if (len > 0.00001f) {
		axis = q * ((float) (1.0f / len));
	}
	else {
		axis.set_value(0.0f, 0.0f, 1.0f);
	}
	return axis;
}


vector<float> Quaternion::to_matrix3() const
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

	return m;
}



float Quaternion::real() const
{
	return e0;
}

Vec3f Quaternion::unreal() const
{
	return Vec3f(e1, e2, e3);
}

vector < float >Quaternion::as_list() const
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
	const float err_limit = 0.00001f;
	
	vector < float >v1 = q1.as_list();
	vector < float >v2 = q2.as_list();

	for (size_t i = 0; i < v1.size(); i++) {
		if (fabs(v1[i] - v2[i]) > err_limit) {
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


Quaternion Quaternion::interpolate(const Quaternion & from,
								   const Quaternion & to, float t)
{
	const double epsilon = 0.00001;
	double cosom = from.e1 * to.e1 + from.e2 * to.e2 + from.e3 * to.e3 + from.e0 * to.e0;

	Quaternion q;
	if (cosom < 0) {
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

	float scale0f = (float) scale0;
	float scale1f = (float) scale1;
	
	return (scale0f * from + scale1f * q);
}
