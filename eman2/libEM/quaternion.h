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

#ifndef eman__quaternion_h__
#define eman__quaternion_h__ 1

#include "vec3.h"

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
		explicit Quaternion(const vector<float> & matrix3);
		
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

		vector<float> to_matrix3() const;

		float real() const;
		Vec3f unreal() const;

		vector < float > as_list() const;

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

#endif	//eman__quaternion_h__
