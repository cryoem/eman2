/**
 * $Id$
 */
#include "transform.h"

using namespace EMAN;

Transform & Transform::set_transform_instance(const Vec3f &translation,
											  const Rotation & rotation,
											  const Vec3f &scale_factor,
											  const Rotation & scale_orientation,
											  const Vec3f &center)
{
	Rotation so = scale_orientation;

	matrix.make_identity();
	translate((-1.0f) * center);

	if (so != Rotation(1, 0, 0, 0, Rotation::QUATERNION)) {
		Rotation inverse_so = so.create_inverse();
		rotate(inverse_so);
	}
	scale(scale_factor);
	rotate(so);
	rotate(rotation);
	translate(center);
	translate(translation);

	return (*this);
}



Vec3f Transform::transform(const Vec3f &x)
{
	Vec3f x1 = x - pre_trans;
	Matrix3f r = matrix.get_matrix3();
	x1 = r * x1;
	x1 += get_post_translate();

	return x1;
}


Vec3f Transform::inverse_transform(const Vec3f &v)
{
	Transform t = create_inverse();
	return t.transform(v);
}

float Transform::get_scale(int i) const
{
	float s = (float)(matrix[i][0] *  matrix[i][0] + matrix[i][1]*matrix[i][1]
		+ matrix[i][2] * matrix[i][2]);
	return s;
}

Vec3f  Transform::get_scale() const
{
	float s1 = (float)(matrix[0][0] *  matrix[0][0] + matrix[0][1]*matrix[0][1]
		+ matrix[0][2] * matrix[0][2]);
	float s2 = (float)(matrix[1][0] *  matrix[1][0] + matrix[1][1]*matrix[1][1]
		+ matrix[1][2] * matrix[1][2]);

	float s3 = (float)(matrix[2][0] *  matrix[2][0] + matrix[2][1]*matrix[2][1]
		+ matrix[2][2] * matrix[2][2]);
	
	return Vec3f(s1, s2, s3);
}


int Transform::get_type() const
{
	return TRANSFORM;
}



Transform EMAN::operator+(const Transform & t1, const Transform & t2)
{
	Transform t = t1;
	t += t2;
	return t;
}

Transform EMAN::operator-(const Transform & t1, const Transform & t2)
{
	Transform t = t1;
	t -= t2;
	return t;
}


Transform EMAN::operator*(const Transform & t1, const Transform & t2)
{
	Transform t = t1;
	t *= t2;
	return t;
}

Transform EMAN::operator*(const Transform & t, float s)
{
	Transform new_t = t;
	new_t *= s;
	return new_t;
}

Transform EMAN::operator*(float s, const Transform & t)
{
	Transform new_t = t;
	new_t *= s;
	return new_t;
}


Transform EMAN::operator/(const Transform & t1, const Transform & t2)
{
	Transform t = t1;
	t /= t2;
	return t;
}

Transform EMAN::operator/(const Transform & t, float s)
{
	Transform new_t = t;
	new_t /= s;
	return new_t;
}

Transform EMAN::operator/(float s, const Transform & t)
{
	Transform new_t = t;
	new_t /= s;
	return new_t;
}

Rotation Rotation::interpolate(const Rotation & from, const Rotation & to, float percent)
{
	Quaternion q1 = from.get_quaternion();
	Quaternion q2 = to.get_quaternion();
	Quaternion q = Quaternion::interpolate(q1, q2, percent);
	return Rotation(q);
}

