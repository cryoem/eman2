/**
 * $Id$
 */
#ifndef eman__transform_h__
#define eman__transform_h__ 1

#include "rotation.h"
#include "vec3.h"
#include "matrix.h"

namespace EMAN
{
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

		Transform(const Vec3f &translation)
		{
			for (int i = 0; i < 3; i++) {
				matrix[i][3] = translation[i];
			}
		}

		Transform(const Rotation & rotation)
		{
			matrix = rotation.get_matrix4();
		}

		Transform(const Rotation & rotation, const Vec3f &post_translation)
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

		Transform & set_translate_instance(const Vec3f &s) {
			matrix.make_identity();
			return set_post_translate(s);
		}

		Transform & set_scale_instance(const Vec3f &s) {
			matrix.make_identity();
			for (int i = 0; i < 3; i++) {
				matrix[i][i] = s[i];
			}
			return (*this);
		}

		Transform & set_transform_instance(const Vec3f &translation,
										   const Rotation & roration,
										   const Vec3f &scale_factor,
										   const Rotation & scale_orientation,
										   const Vec3f &center);

		Transform & set_transform_instance(const Vec3f &translation,
										   const Rotation & rotation,
										   const Vec3f &scale_factor)
		{
			return set_transform_instance(translation, rotation, scale_factor,
										  Rotation(1, 0, 0, 0, Rotation::QUATERNION),
										  Vec3f (0, 0, 0));
		}

		Transform & set_center(const Vec3f &c)
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

		Transform & set_post_translate(const Vec3f &s) {
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
		Transform & translate(const Vec3f &v) {
			if (v != Vec3f (0, 0, 0)) {
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

		Transform & rotate_center(const Rotation & r, const Vec3f center) {
			translate((float) -1.0 * center);
			rotate(r);
			translate(center);
			return (*this);
		}

		Transform & rotate_scale(const Rotation & rotation,
								 const Vec3f &scale_factor, const Vec3f &center) {
			scale(scale_factor);
			rotate_center(rotation, center);
			return (*this);
		}

		Transform & pre_translate_rotate(const Vec3f t, const Rotation & r) {
			translate(t);
			rotate(r);
			return (*this);
		}

		Transform & post_translate_rotate(const Rotation & r, const Vec3f t) {
			rotate(r);
			translate(t);
			return (*this);
		}

	/** Concatenates this transform with a scaling transformation.
	 */
		Transform & scale(const Vec3f &scale_factor) {
			if (scale_factor != Vec3f (1, 1, 1)) {
				Matrix4f m;
				m.make_identity();
				for (int i = 0; i < 3; i++) {
					m[i][i] = scale_factor[i];
				}
				matrix *= m;
			}
			return (*this);
		}


		Vec3f transform(const Vec3f &v);
		Vec3f inverse_transform(const Vec3f &v);

		Rotation get_rotation() const
		{
			return Rotation(matrix.get_matrix3());
		}

		float get_scale(int i) const;		
		Vec3f get_scale() const;
		
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
			return Vec3f (matrix[0][3], matrix[1][3], matrix[2][3]);
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
		Vec3f pre_trans;
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
