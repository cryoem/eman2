/**
 * $Id$
 */

#ifndef eman__quaternion_h__
#define eman__quaternion_h__ 1

#include "vec3.h"
#include "matrix.h"
#include <math.h>

namespace EMAN
{
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
		Quaternion(float radians, const Vec3f &axis);
		Quaternion(const Vec3f &axis, float radians);
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

		Vec3f rotate(const Vec3f &v) const;

		float to_angle() const;
		Vec3f to_axis() const;

		Matrix3f to_matrix3() const;
		Matrix4f to_matrix4() const;

		float real() const;
		Vec3f unreal() const;

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
}

#endif
