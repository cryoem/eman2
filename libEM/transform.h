/**
 * $Id$
 */
#ifndef eman__transform_h__
#define eman__transform_h__ 1

#include "vec3.h"
#include <vector>
#include <map>
#include <string>

using std::string;
using std::map;
using std::vector;

namespace EMAN
{
	/** Transform defines a transformation, which can be rotation,
     * translation, scale, and their combinations.
	 *
	 * Internally a transformation is stored in a 4x3 matrix.
	 *    a b c
	 *    e f g
	 *    i j k
	 *    m n o
	 *
	 * The left-top 3x3 submatrix
	 *	  a b c
	 *    e f g
	 *    i j k
	 * provides rotation, scaling and skewing.
	 *
	 * post translation is stored in (m,n,o).
	 *
	 * a separate vector containing the pretranslation, with an implicit
	 * column [0 0 0 1] at the end when 4x4 multiplies are required.
	 *
	 * The 'center of rotation' is NOT implemented as a separate vector, 
	 * but as a combination of pre and post translations.
	 *
	 * a matrix  M is called orthogonal M * transpose(M) = 1. All
	 * Orthogonal Matrices have determinants of 1
	 *
	 *
     */
	class Transform
	{
	public:
		static const float ERR_LIMIT;
		enum EulerType
		{
			UNKNOWN,
			EMAN,
			IMAGIC,
			SPIN,
			QUATERNION,
			MATRIX,
			SGIROT,
			SPIDER,
			MRC
		};
	public:
		Transform();

		
		Transform(const vector<float> & rotation, EulerType euler_type);
		
		Transform(const vector<float> & rotation, EulerType euler_type, 
				  const Vec3f& posttrans);
		
		Transform(const vector<float> & rotation, EulerType euler_type, 
				  const Vec3f & pretrans, const Vec3f& posttrans);
		
		Transform(map<string, float>& rotation, EulerType euler_type);

		Transform(map<string, float>& rotation, EulerType euler_type, 
				  const Vec3f& posttrans);
		
		Transform(map<string, float>& rotation, EulerType euler_type, 
				  const Vec3f & pretrans, const Vec3f& posttrans);
		
		virtual ~ Transform();

		void to_identity();

		// reorthogonalize the matrix
		void orthogonalize();	
		Transform inverse();
	
		void set_pretrans(const Vec3f & pretrans);
		void set_posttrans(const Vec3f & posttrans);
		void set_center(const Vec3f & center);
		void set_rotation(map<string, float>& rotation, EulerType euler_type);
		void set_rotation(const vector<float>& rotation, EulerType euler_type);
		void set_scale(float scale);

		Vec3f get_pretrans() const;
		Vec3f get_posttrans() const;
		Vec3f get_center() const;
		map<string, float> get_rotation(EulerType euler_type) const;
		
		float get_scale() const;
		
		// returns orthogonality coefficient 0-1 range;
		float is_orthogonal() const;
		
		float * operator[] (int i);
		const float * operator[] (int i) const;

		// This should return an iterator which will produce
		// all of the symmetrically equivalent Transform objects
		// for the given symmetry including the current Transform
		//iterator get_symmetry(string symname);
		
	  private:
		void init();
		float eman_alt() const;
		float eman_az() const;
		float eman_phi() const;
		vector<float> matrix2quaternion() const;
		vector<float> matrix2sgi() const;
		void quaternion2matrix(float e0, float e1, float e2, float e3);
		
		float **matrix;
		Vec3f pre_trans;
	};
#if 0
	Transform operator+(const Transform & t1, const Transform & t2);
	Transform operator-(const Transform & t1, const Transform & t2);

	Transform operator*(const Transform & t1, const Transform & t2);
	Transform operator*(const Transform & t, float s);
	Transform operator*(float s, const Transform & t);

	Transform operator/(const Transform & t1, const Transform & t2);
	Transform operator/(const Transform & t, float s);
	Transform operator/(float s, const Transform & t);
#endif
}



#endif
