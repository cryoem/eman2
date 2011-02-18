/**
 * $Id$
 */

/*
 * Author: Steven Ludtke (sludtke@bcm.edu)
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


#ifndef eman__transform_h__
#define eman__transform_h__ 1

#ifdef _WIN32
	#pragma warning(disable:4819)
#endif	//_WIN32

#include "vec3.h"
#include "emobject.h"
#include <cstdio>

// #include <vector>
// using std::vector;
//
// #include <map>
// using std::map;
//
// #include <string>
// using std::string;



namespace EMAN
{
	/** A Transform object is a somewhat specialized object designed specifically for EMAN2/Sparx storage of
	 * alignment parameters and euler orientations.
	 * It's designed to store four transformations in a specific order, namely
	 * Transform = MTSR
	 * Where M is a mirroring operation (about the x-axis) or the identity
	 * T is a Translation matrix
	 * S is a uniform scaling matrix
	 * R is a rotation matrix
	 * This means you can call set_scale, set_trans, set_rotation in any order but still have the operations arranged
	 * internally in the order of MTSR. This is somewhat restrictive, for example in the context of how OpenGL handles
	 * transformations, but in practice is nicely suited to the situations that arise in EMAN2 - namely, alignment and
	 * projection orientation characterization.
	 *
	 * Note that you can fool the Transform object into storing any matrix by using the constructors that take array arguments.
	 * This can useful, for example, for shearing your image.
	 *
	 * See http://blake.bcm.tmc.edu/emanwiki/Eman2TransformInPython for using it from Python and detailed discussion
	 * See test_transform.py for examples of the way it is unit tested
	 * See http://blake.bcm.tmc.edu/emanwiki/EMAN2/Tutorials/RotateTranslate for examples showing how to transform EMDatas with it.
	 * @author David Woolford (based on Philip Baldwin's original Transform3D code)
	 * @date September 2008
	 * @ingroup tested3c
	 */
	class Transform {
// 		friend Transform EMAN::operator*(const Transform & M2, const Transform & M1);
		public:
			static const float ERR_LIMIT;

			/** Default constructor
			* Internal matrix is the identity
			*/
			Transform();

			/** Copy constructor
			* @param rhs the object to be copied
			*/
			Transform( const Transform& rhs );

			/** Assignment operator
			 * @param that that which this will become
			 */
			Transform& operator=(const Transform& that);

			/** Equality comparision operator
			 * @param rhs the Transform object compared to
			 * @return a boolean indicate if the pass in transform is equal this object
			 * */
			bool operator==(const Transform& rhs) const;

			/** Unequality comparision operator
			 * @param rhs the Transform object compared to
			 * @return a boolean indicate if the pass in transform is not equal this object
			 * */
			bool operator!=(const Transform& rhs) const;

			/** Construction using a dictionary
			 * @param d the dictionary containing the parameters
			 */
			Transform(const Dict& d);

			/** Construction using an array of floats
			 * @param array the array of values that will become the internal matrix. row order (3 rows of 4)
			 */
			Transform(const float array[12]);

			/** Construction using a vector of size 12
			 * @param array the array of values that will become the internal matrix. row order (3 rows of 4)
			 */
			Transform(const vector<float> array);


			~Transform() { }

			/** Set a rotation using a specific Euler type and the dictionary interface
			* Works for all Euler types
			* @param rotation a dictionary containing all key-entry pair required of the associated Euler type
			*/
			void set_rotation( const Dict &rotation );

			/** Determine the rotation that would transform a vector pointing in the Z direction
			 * so that it points in the direction of the argument vector
			 * Automatically normalizes the vector
			 * @param v the direction you want to solve for
			 */
			void set_rotation(const Vec3f & v);

			/** Get a rotation in any Euler format
			 * @param euler_type the requested Euler type
			 * @return a dictionary containing the key-entry pairs describing the rotations in terms of the requested Euler type
			 */
			Dict get_rotation(const string& euler_type = "eman") const;


			/** Get the rotation part of the tranformation matrix as a Transform object
			 * @return a Transform describing the rotation only
			 */
			Transform get_rotation_transform() const;

			/** How do I get the transform that will yield the horizontally flipped projection?
			 * @return the transform that will yield the horizontally flipped projection
			 */
			Transform get_hflip_transform() const;

			/** How do I get the transform that will yield the vertically flipped projection?
			 * @return the transform that will yield the vertically flipped projection
			 */
			Transform get_vflip_transform() const;

			/** Set the parameters of the entire transform.
			 * keys acted upon are "type" - if this exists then the correct euler angles need to be included -
			 * also "tx","ty","tz", "scale", and "mirror"
			 * @param d the dictionary containing the parameters
			 */
			void set_params(const Dict& d);

			/** Set the parameters of the entire transform as though they there in the inverse format

			 * in other words, calling set_params_inverse(get_params_inverse()) should essentially leave
			 * the object unchanged.
			 *  @param d the dictionary containing the inverse parameters
			 */
			void set_params_inverse(const Dict& d);

			/** Get the parameters of the entire transform, using a specific euler convention
			 * @param euler_type the euler type of the retrieved rotation
			 * @return a dictionary containing the parameters
			 */
			Dict get_params(const string& euler_type) const;

			/** Get the parameters of the inverse of the transform as though it were in RSMT order not MTSR
			 * @param euler_type the euler type of the retrieved rotation
			 * @return a dictionary containing the parameters
			 */
			Dict get_params_inverse(const string& euler_type) const;

			//=============== set and get post trans =============
			/** Set the post translation component
			 * @param x the x translation
			 * @param y the y translation
			 * @param z the z translation
			 */
			void set_trans(const float& x, const float& y, const float& z=0);

			/** Set the post translation component using a Vec3f
			 * @param v the 3D translation vector
			 */
			inline void set_trans(const Vec3f& v) { set_trans(v[0],v[1],v[2]); }

			/** Set the post translation component using a Vec2f
			 * @param v the 2D translation vector
			 */
			inline void set_trans(const Vec2f& v) { set_trans(v[0],v[1]); }

			/** Get the post trans as a vec3f
			 * @return the  translation vector
			 */
			Vec3f get_trans() const;

			/** Get the degenerant 2D post trans as a vec2f
			 * @return the 2D translation vector
			 */
			Vec2f get_trans_2d() const;

			//================= get pre trans is supported =========

			/** Get the translation vector as though this object was MSRT_ not MTSR, where T_ is what you want
			 * Note M means post x mirror, T means translation, S means scale, and R means rotaiton
			 * @return the pre translation vector
			 */
			Vec3f get_pre_trans() const;

			/** 2D version of getting the translation vector as though this object was MSRT_ not MTSR, where T_ is what you want
			 * Note M means post x mirror, T means translation, S means scale, and R means rotation
			 * @return the pre translation vector
			 */
			Vec2f get_pre_trans_2d() const;

			/** Set the translational component of the matrix as though it was MSRT_ not MTSR, where
			 * T_ is the pre translation. Internally the correct form of MTSR is computed.
			 * @param v the vector (Vec3f or Vec2f) that is the pre trans
			 */
			template<typename type>
			void set_pre_trans(const type& v);

			//=============== set and get scale =============
			/** Set the scale
			 * @param scale the amount to scale by
			 */
			void set_scale(const float& scale);

			/** Get the scale that was applied
			 * @return the scale factor
			 */
			float get_scale() const;

			//=============== set and get post x mirror =============
			/** Query whether x_mirroring is occuring
			 * @return whether x_mirroring is occuring
			 */
			bool get_mirror() const;

			/** Set whether or not x_mirroring is occuring
			 * @param x_mirror whether x_mirroring should be applied
			 */
			void set_mirror(const bool x_mirror);

			//=============== other stuff ============================
			/** Get scale and x_mirror with 1 function call.
			 * More efficient than calling get_scale and get_x_mirror separately
			 * @param scale a reference to the value that will be assigned the scale value
			 * @param x_mirror a reference to the value that will be assigned the x_mirror value
			 */
			void get_scale_and_mirror(float& scale, bool& x_mirror) const;

			/** Force the internal matrix to become the identity
			*/
			void to_identity();

			/** Returns whethers or this matrix is the identity
			 */
			bool is_identity() const;

			/** Reorthogonalize the rotation part of the matrix in place.
			 * Does this by performing the SVD decomposition of the rotation matrix R
			 * such that R = USV^T - since the eigenvalues of a rotation matrix are all 1 we
			 * enforce that S should be the identity and produce a corrected matrix R' = UV^T
			 *
			 */
			void orthogonalize();

			/** Get the determinant of the matrix
			 * @return the determinant
			 */
			float get_determinant() const;

			/** Print the contents of the internal matrix verbatim to standard out
			 */
			void printme() const {
				printf("%8.6f %8.6f %8.6f %8.6f\n",matrix[0][0],matrix[0][1],matrix[0][2],matrix[0][3]);
				printf("%8.6f %8.6f %8.6f %8.6f\n",matrix[1][0],matrix[1][1],matrix[1][2],matrix[1][3]);
				printf("%8.6f %8.6f %8.6f %8.6f\n",matrix[2][0],matrix[2][1],matrix[2][2],matrix[2][3]);
				printf("%8.6f %8.6f %8.6f %8.6f\n",0.0,0.0,0.0,1.0);

			}

			//=============== get set matrix ============================
			/** Set the transformation matrix using a vector. Must be of length 12.
			 * @param v the transformation matrix stored as a vector - 3 rows of 4.
			 */
			void set_matrix(const vector<float>& v);

			/** Get the transformation matrix using a vector.
			 * @return a vector - 3 rows of 4 - that stores the values of the transformation matrix
			 */
			vector<float> get_matrix() const;

			/** Get the inverse of this transformation matrix
			 * @return the inverse of this transformation matrix
			 */
			void invert();

			/** Get the inverse of this transformation matrix
			 * @return the inverse of this transformation matrix
			 */
			Transform inverse() const;

			/** Get the transpose of this transformation matrix
		 	 */
			void transpose_inplace();

			/** Get the transpose of this transformation matrix
			 * @return the transpose of this transformation matrix
			 */
			Transform transpose() const;

			/** Get the value stored in the internal transformation matrix at at coordinate (r,c)
			 */
			inline float at(int r,int c) const { return matrix[r][c]; }

			/** Set the value stored in the internal transformation matrix at at coordinate (r,c) to value
			 */
			inline void set(int r, int c, float value) { matrix[r][c] = value; }

			/** Operator[] convenience
			 * so Transform3D[2][2] etc terminology can be used
			 */
			inline float * operator[] (int i) { return matrix[i]; }

			/** Operator[] convenience
			 * so Transform3D[2][2] etc terminology can be used
			 */
			inline const float * operator[] (int i) const { return matrix[i]; }

			/** Transform 2D coordinates using the internal transformation matrix
			* @param x the x coordinate of the transformed point
			* @param y the y coordinate of the transformed point
			* @return the transformed vector
			 */
			inline Vec2f transform(const float& x, const float& y) const {
//				assert_valid_2d();
				Vec2f ret;
				ret[0] = matrix[0][0]*x + matrix[0][1]*y + matrix[0][3];
				ret[1] = matrix[1][0]*x + matrix[1][1]*y + matrix[1][3];
				return ret;
			}

			/** Transform a 2D vector using the internal transformation matrix
			* @param v a two dimensional vector to be transformed
			* @return the transformed vector
			 */
			template<typename Type>
			inline Vec2f transform(const Vec2<Type>& v) const {
				return transform(v[0],v[1]);
			}

			/** Transform 3D coordinates using the internal transformation matrix
			 * @param x the x coordinate of the transformed point
			 * @param y the y coordinate of the transformed point
			 * @param z the z coordinate of the transformed point
			 * @return the transformed vector
			 */
			inline Vec3f transform(const float& x, const float& y, const float& z) const {
// 				assert_consistent_type(THREED);
				Vec3f ret;
				ret[0] = matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z + matrix[0][3];
				ret[1] = matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z + matrix[1][3];
				ret[2] = matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z + matrix[2][3];
				return ret;
			}

			/** Transform a 3D vector using the internal transformation matrix
			 * @param v a three dimensional vector to be transformed
			 * @return the transformed vector
			 */
			template<typename Type>
			inline Vec3f transform(const Vec3<Type>& v) const {
// 				assert_consistent_type(THREED); // Transform does the assertion
				return transform(v[0],v[1],v[2]);
			}


			/** Get a matrix row as a Vec3f
			 * required for back compatibility with Tranform3D - see PawelProjector
			* @param i the row number (starting at 0)
			* @return the ith row
			*/
			inline Vec3f get_matrix3_row(int i) const {
				return Vec3f(matrix[i][0], matrix[i][1], matrix[i][2]);
			}

			/** get the number of symmetries associated with the given symmetry name
			 */
			static int get_nsym(const string & sym);

			/** Apply the symmetry deduced from the function arguments to this Transform and
			 * return the result
			 *
			 */
			Transform get_sym(const string & sym, int n) const;

			/**
			*/
			void copy_matrix_into_array(float* const) const;

			/** Negates the Transform - a useful way, for example, for getting an orientation on the opposite side
			* of the sphere
			* @return a transform that has been negated
			*/
			Transform negate() const;

			/** Get the transform that moves any tetrahedron generated by eman2 so that it matches the 2-2-2 (MRC, FREALIGN) convention
			 * @return a Transforms, alt=54.73561, phi=45
			 */
			static Transform tet_3_to_2();

			/** Get the transform that moves any icosahedron generated by eman2 so that it matches the 2-2-2 (MRC, FREALIGN) convention
			 * @return a Transforms, alt=58.282523, az=270
			*/
			static Transform icos_5_to_2();

		private:
			float matrix[3][4];

			void assert_valid_2d() const;

			/// This map is used to validate keys in the argument maps - e.g. if the type is 2d and the angle is not "alpha" then we should throw
			static vector<string> permissable_2d_not_rot;
			static vector<string> permissable_3d_not_rot;
			static map<string,vector<string> >  permissable_rot_keys;

			/** Called internally to initialize permissable_2d_not_rot, permissable_3d_not_rot, and permissable_rot_keys static members
			*/
			void init_permissable_keys();

			/** Test to ensure the parametes in the given dictionary are valid
			* Throws if an error is detected
			* Generic - works in every circumstance (set_params, set_rotation, set_params_inv)
			* Uses static members permissable_2d_not_rot, permissable_3d_not_rot, and permissable_rot_keys as basis of decision
			* @param d the dictionary that was the function argument of the set_params, set_rotation or the  set_params_inv function
			* @exception InvalidParameterException if the dictionary is invalid in anyway
			*/
			void detect_problem_keys(const Dict& d);

	};
	/// Matrix times Matrix, a pure mathematical operation
	Transform operator*(const Transform & M2, const Transform & M1);

	/// Matrix times Vector, a pure mathematical operation
	template<typename Type>
	Vec3f operator*( const Transform& M, const Vec3<Type> & v)
	{
		return M.transform(v);
	}

	/// Matrix times Vector, a pure mathematical operation
	template<typename Type>
	Vec2f operator*( const Transform& M, const Vec2<Type> & v)
	{
		return M.transform(v);
	}

	/// Vector times a matrix. Highly specialized. Useful when
	/// the upper 3x3 only contains rotations and you want to quickly
	/// multiply by the rotation matrix inverse (transpose)
	template<typename Type>
	Vec3f operator*(const Vec3<Type> & v, const Transform & M)
	{
		float x = v[0] * M[0][0] + v[1] * M[1][0] + v[2] * M[2][0] ;
		float y = v[0] * M[0][1] + v[1] * M[1][1] + v[2] * M[2][1];
		float z = v[0] * M[0][2] + v[1] * M[1][2] + v[2] * M[2][2];
		return Vec3f(x, y, z);
	}

	template<typename type>
	void Transform::set_pre_trans(const type& v) {

		Transform tmp;
		Dict rot = get_rotation("eman");
		tmp.set_rotation(rot);

		float scale = get_scale();
		if (scale != 1.0 ) tmp.set_scale(scale);

		Transform trans;
		trans.set_trans(v);

		trans = tmp*trans;

		Transform tmp2;
		tmp2.set_rotation(rot);
		tmp2.invert(); // invert
		if (scale != 1.0 ) tmp2.set_scale(1.0f/scale);


		trans = trans*tmp2;

		set_trans(trans.get_trans());
	}

	/** Transform3D
	 * These are  a collection of transformation tools: rotation, translation,
	 * and construction of symmetric objects
	 * @author Philip Baldwin <Philip.Baldwin@uth.tmc.edu> and Steven Ludtke
	 * @date $Date: 2005/04/04 17:41pm
	 * @see Phil's article
	 * Transform defines a transformation, which can be rotation,
	 * translation, scale, and their combinations.
	 *
	 * Internally a transformation is stored in a 4x4 matrix.
	 *         a b c d
	 *         e f g h           R        v
	 *  M=     j k m n    =      vpre     1    , where R is 3by3, v is 3by1
	 *         p q r 1
	 *  The standard computer graphics convention is identical to ours after setting vpre to
	 *  zero and can be found in many references including Frackowiak et al; Human Brain Function
	 *
	 *
	 * The left-top 3x3 submatrix
	 *
	 *        a b c
	 *   R =  e f g
	 *        j k m
	 *
	 * provides rotation, scaling and skewing (not yet implimented).
	 *
	 * The cumulative translation is stored in (d, h, n).
	 * We put the post-translation into (p, q, r), since it is convenient
	 * to carry along at times. When matrices are multiplied or printed, these
	 * are hidden to the user. They can only be found by applying the post_translation
	 * method, and these elements are non-zero. Otherwise the post-translation method returns
	 * the cumulative translationmlb
	 *
	 * If rotations need to be found
	 * around alternate origins, then brief calculations need to be performed
	 * Pre and Post Translations should be kept as separate vectors
	 *
	 * a matrix  R is called orthogonal if
	 *           R * transpose(R) = 1.
	 * All Real Orthogonal Matrices have eigenvalues with unit modulus and determinant
	 * therefore equal to  "\pm 1"
	 * @ingroup tested3c
     */
//	class Transform3D
//	{
//	public:
//		static const float ERR_LIMIT;
//		enum EulerType
//		{
//			UNKNOWN,
//			EMAN,
//			IMAGIC,
//			SPIN,
//			QUATERNION,
//			SGIROT,
//			SPIDER,
//			MRC,
//			XYZ,
//			MATRIX
//		};
//
//		/** Default constructor
//		 * Internal matrix is the identity
//		 */
//		Transform3D();
//
//		/** Copy constructor
//		 * @param rhs the object to be copied
//		 */
//	    Transform3D( const Transform3D& rhs );
//
//		 /** Construct a Transform3D object describing a rotation, assuming the EMAN Euler type
//		 * @param az EMAN - az
//		 * @param alt EMAN - alt
//		 * @param phi EMAN - phi
//		 */
//		Transform3D(const float& az,const float& alt,const float& phi); // EMAN by default
//
//		 /** Construct a Transform3D object describing a rotation (assuming the EMAN Euler type) and
//		 * a post translation
//		 * @param az EMAN - az
//		 * @param alt EMAN - alt
//		 * @param phi EMAN - phi
//		 * @param posttrans the post translation vector
//		 */
//		Transform3D(const float& az, const float& alt, const float& phi, const Vec3f& posttrans);
//
//		 /** Construct a Transform3D object describing a pre trans, a rotation
//		 * assuming the EMAN Euler type) and a post translation
//		 * @param pretrans the pre translation vector
//		 * @param az EMAN - az
//		 * @param alt EMAN - alt
//		 * @param phi EMAN - phi
//		 * @param posttrans the post translation vector
//		  */
//		Transform3D(const Vec3f & pretrans, const float& az,const float& alt,const float& phi, const Vec3f& posttrans);
//
//		/** Construct a Transform3D object describing a rotation, using a specific Euler type.
//		 * works for EMAN, SPIDER, MRC, IMAGIC and  XYZ (Euler types required 3 parameters)
//		 * @param euler_type the Euler type either  EMAN, SPIDER, MRC, IMAGIC and  XYZ
//		 * @param a1 EMAN - az, SPIDER - phi, MRC - phi, IMAGIC - alpha, XYZ - xtilt
//		 * @param a2 EMAN - alt, SPIDER - theta, MRC - theta, IMAGIC - beta, XYZ - ytilt
//		 * @param a3 EMAN - phi, SPIDER - psi, MRC - omega, IMAGIC - gamma, XYZ - ztilt
//		 */
//		Transform3D(EulerType euler_type, const float& a1, const float& a2, const float& a3) ;
//
//		/** Construct a Transform3D object describing a rotation, using a specific Euler type.
//		 * Works for Euler types that require 4 parameters
//		 * @param euler_type the Euler type either QUATERNION, SPIN or SGIROT
//		 * @param a1 QUATERNION - e0, SPN and SGIROT - Omega
//		 * @param a2 QUATERNION - e1, SPN and SGIROT - n1
//		 * @param a3 QUATERNION - e2, SPN and SGIROT - n2
//		 * @param a4 QUATERNION - e3, SPN and SGIROT - n3
//		*/
//		Transform3D(EulerType euler_type, const float& a1, const float& a2, const float& a3, const float& a4) ; // o
//
//		/** Construct a Transform3D object consisting only of a rotation, using a specific Euler type.
//		 * Works for all Euler types
//		 * @param euler_type any Euler type
//		 * @param rotation a dictionary containing all key-entry pair required of the associated Euler type
//		 */
//		Transform3D(EulerType euler_type, const Dict& rotation);
//
//		/** Construct a Transform3D object consisting only of a rotation by initializing the internal
//		 * rotation matrix component wise
//		 * @param mij the value to be placed in the internal transformation matrix at coordinate (i-1,j-1)
//		 */
//		Transform3D(const float& m11, const float& m12, const float& m13,
//					const float& m21, const float& m22, const float& m23,
//	 				const float& m31, const float& m32, const float& m33);
//
//		/** Destructor
//		 */
//		~Transform3D();
//
//		/** FIXME insert comments
//		 */
//		void apply_scale(const float& scale);
//
//		/** FIXME insert comments
//		 */
//		void set_scale(const float& scale);
//
//		/** Reorthogonalize the matrix
//		 */
//		void orthogonalize();
//
//		/**  create the transpose in place
//		 */
//		void transpose();
//
//		void set_rotation(const float& az, const float& alt,const float& phi);
//
//		/** Sets the rotation as defined by the EulerType
//		 * works for EMAN, SPIDER, MRC, IMAGIC and  XYZ
//		 * @param euler_type the Euler type either  EMAN, SPIDER, MRC, IMAGIC and  XYZ
//		 * @param a1 EMAN - az, SPIDER - phi, MRC - phi, IMAGIC - alpha, XYZ - xtilt
//		 * @param a2 EMAN - alt, SPIDER - theta, MRC - theta, IMAGIC - beta, XYZ - ytilt
//		 * @param a3 EMAN - phi, SPIDER - psi, MRC - omega, IMAGIC - gamma, XYZ - ztilt
//		*/
//		void set_rotation(EulerType euler_type,const float& a1,const float& a2,const float& a3);
//
//		/** Set quaternion-based rotations
//		 * Works for QUATERNION, SPIN and SGIROT
//		 * @param euler_type the Euler type either QUATERNION, SPIN or SGIROT
//		 * @param a1 QUATERNION - e0, SPN and SGIROT - Omega
//		 * @param a2 QUATERNION - e1, SPN and SGIROT - n1
//		 * @param a3 QUATERNION - e2, SPN and SGIROT - n2
//		 * @param a4 QUATERNION - e3, SPN and SGIROT - n3
//		 */
//		void set_rotation(EulerType euler_type,const float& a1,const float& a2,const float& a3, const float& a4);
//
//		/** set the internal rotation matrix component wise
//		 * @param mij the value to be placed in the internal transformation matrix at coordinate (i-1,j-1)
//		 */
//		void set_rotation(const float& m11, const float& m12, const float& m13,
//						  const float& m21, const float& m22, const float& m23,
//						  const float& m31, const float& m32, const float& m33);
//
//// 		void set_params(const Dict& params, const string& parameter_convention, const EulerType euler_type = EMAN);
////
//// 		Dict get_params(const string& parameter_convention, const EulerType euler_type = EMAN) const;
//
//		/** Set a rotation using a specific Euler type and the dictionary interface
//		 * Works for all Euler types
//		 * @param euler_type any Euler type
//		 * @param rotation a dictionary containing all key-entry pair required of the associated Euler type
//		 */
//		void set_rotation(EulerType euler_type, const Dict &rotation );
//
//		/** Get a rotation in any Euler format
//		 * @param euler_type the requested Euler type
//		 * @return a dictionary containing the key-entry pairs describing the rotations in terms of the requested Euler type
//		 */
//		Dict get_rotation(EulerType euler_type=EMAN) const;
//
//		/** returns a rotation that maps a pair of unit vectors, a,b to a second  pair A,B
//		 * @param eahat, ebhat, eAhat, eBhat are all unit vectors
//		 * @return  a transform3D rotation
//		 */
//		void set_rotation(const Vec3f & eahat, const Vec3f & ebhat,
//		                                    const Vec3f & eAhat, const Vec3f & eBhat);
//		/** returns the magnitude of the rotation
//		*/
//		float get_mag() const;
//
//		/** returns the spin-axis (or finger) of the rotation
//		*/
//		Vec3f get_finger() const;
//
//		/** Gets one of two pre translation vectors
//		 * @param flag if 0 returns the pre translation vector, if 1 all translation is treated as pre
//		 * @return the translation vector
//		 */
//		Vec3f get_pretrans( int flag=0) const; // flag=1 => all trans is pre
//
//		/** Gets one of two post translation vectors.
//		 * when the flag is 1 then the contents of the Transform3D matrix right column are returned
//		 * @param flag if 0 returns the post translation vector, if 1 all translation is treated as post
//		 * @return the translation vector
//		 */
//		Vec3f get_posttrans(int flag=0) const; // flag=1 => all trans is post
//
//		/** Get the total translation as a post translation. Calls get_postrans(1)
//		 * @return the translation vector
//		 */
//		Vec3f get_total_posttrans() const;
//
//		/** Get the total translation as a pre translation. Calls get_pretrans(1)
//		 * @return the translation vector
//		 */
//		Vec3f get_total_pretrans() const;
//
//		/** This doesn't do anything, it returns an empty vector. Why is it being used?
//		 */
// 		Vec3f get_center() const; // This doesn't do anything
//
//		/** Get a matrix column as a Vec3f
//		 * @param i the column number (starting at 0)
//		 * @return the ith column
//		 */
//		Vec3f get_matrix3_col(int i) const;
//
//		/** Get a matrix row as a Vec3f
//		 * @param i the row number (starting at 0)
//		 * @return the ith row
//		 */
//		Vec3f get_matrix3_row(int i) const;
//
//		/** Perform a full transform a Vec3f using the internal transformation matrix
//		 * @param v3f the vector to be transformed
//		 * @return the transformed vector
//		 */
//		Vec3f transform(const Vec3f & v3f) const;
//
//		/** Rotate a Vec3f using the internal rotation matrix
//		 * @param v3f the vector to be rotated
//		 * @return the rotated vector
//		 */
//		Vec3f rotate(const Vec3f & v3f) const;
//
//		/** FIXME insert comments
//		 */
//		Transform3D inverseUsingAngs() const;
//
//		/** FIXME insert comments
//		 */
//		Transform3D inverse() const;
//
//
//		/** Print the Transform3D matrix
//		 */
//		void printme() const {
//			for (int i=0; i<3; i++) {
//				printf("%6.15f\t%6.15f\t%6.15f\t%6.1f\n",
//					   matrix[i][0],matrix[i][1],matrix[i][2],matrix[i][3]);
//			}
//			printf("%6.3f\t%6.3f\t%6.3f\t%6.3f\n",0.0,0.0,0.0,1.0);
//			printf("\n");
//		}
//
//		/** Get the value stored in the internal transformation matrix at at coordinate (r,c)
//		 */
//		inline float at(int r,int c) const { return matrix[r][c]; }
//
//		/** Set the value stored in the internal transformation matrix at at coordinate (r,c) to value
//		 */
//		void set(int r, int c, float value) { matrix[r][c] = value; }
//
//		/** Operator[] convenience
//		 * so Transform3D[2][2] etc terminology can be used
//		 */
//		inline float * operator[] (int i) { return matrix[i]; }
//
//		/** Operator[] convenience
//		 * so Transform3D[2][2] etc terminology can be used
//		 */
//		inline const float * operator[] (int i) const { return matrix[i]; }
//
//		static int get_nsym(const string & sym);
//		Transform3D get_sym(const string & sym, int n) const;
//
//		/** Set functions FIXME insert more comments from here down
//		 */
//		void set_center(const Vec3f & center);
//		void set_pretrans(const Vec3f & pretrans);
//		void set_pretrans(const float& dx, const float& dy, const float& dz);
//		void set_pretrans(const float& dx, const float& dy);
//		void set_pretrans(const Vec2f& pretrans);
//		void set_posttrans(const Vec3f & posttrans);
//		void set_posttrans(const float& dx, const float& dy, const float& dz);
//		void set_posttrans(const float& dx, const float& dy);
//		void set_posttrans(const Vec2f& posttrans);
//
//		void set_post_x_mirror(const bool b) { post_x_mirror = b; }
//		bool get_post_x_mirror() const { return post_x_mirror; }
//
//		float get_scale() const;
//
//		void to_identity();
//        bool is_identity();
//
//		/** Convert a list of euler angles to a vector of Transform3D objects.
//		 *
//		 *	@param[in] eulertype The type of Euler angles that is being passed in.
//		 *	@param[in] angles A flat vector of angles.
//		 *
//		 *	@return Vector of pointers to Transform3D objects.
//		 */
//		static vector<Transform3D*>
//			angles2tfvec(EulerType eulertype, const vector<float> angles);
//
//		/** This added by d.woolford, will eventually be removed by author
//		*/
//		void dsaw_zero_hack(){
//			for (int j=0; j<4; ++j) {
//				for (int i=0; i<4; i++) {
//				if ( fabs(matrix[j][i]) < 0.000001 )
//					matrix[j][i] = 0.0;
//				}
//			}
//
//		}
//
//	protected:
//		enum SymType
//		{      CSYM,
//			DSYM,
//			TET_SYM,
//			ICOS_SYM,
//			OCT_SYM,
//			ISYM,
//			UNKNOWN_SYM
//		};
//
//		void init();
//
//		static SymType get_sym_type(const string & symname);
//
//		float matrix[4][4];
//
//		bool post_x_mirror;
//
//		Transform3D::EulerType s;
//	}; // ends Class
//
//	Transform3D operator*(const Transform3D & M1, const Transform3D & M2);
// 	Vec3f operator*(const Vec3f & v    , const Transform3D & M);
// 	Vec3f operator*(const Transform3D & M, const Vec3f & v    );

//	template<typename Type>
//	Vec3f operator*(const Vec3<Type> & v, const Transform3D & M)   // YYY
//	{
////               This is the right multiplication of a row vector, v by a transform3D matrix M
//		float x = v[0] * M[0][0] + v[1] * M[1][0] + v[2] * M[2][0] ;
//		float y = v[0] * M[0][1] + v[1] * M[1][1] + v[2] * M[2][1];
//		float z = v[0] * M[0][2] + v[1] * M[1][2] + v[2] * M[2][2];
//		return Vec3f(x, y, z);
//	}
//
//	template<typename Type>
//	Vec3f operator*( const Transform3D & M, const Vec3<Type> & v)      // YYY
//	{
////      This is the  left multiplication of a vector, v by a matrix M
//		float x = M[0][0] * v[0] + M[0][1] * v[1] + M[0][2] * v[2] + M[0][3];
//		float y = M[1][0] * v[0] + M[1][1] * v[1] + M[1][2] * v[2] + M[1][3];
//		float z = M[2][0] * v[0] + M[2][1] * v[1] + M[2][2] * v[2] + M[2][3];
//		return Vec3f(x, y, z);
//	}
//
//
//	template<typename Type>
//	Vec2f operator*( const Transform3D & M, const Vec2<Type> & v)      // YYY
//	{
////      This is the  left multiplication of a vector, v by a matrix M
//		float x = M[0][0] * v[0] + M[0][1] * v[1] + M[0][3] ;
//		float y = M[1][0] * v[0] + M[1][1] * v[1] + M[1][3];
//		return Vec2f(x, y);
//	}

}  // ends NameSpace EMAN



#endif


/* vim: set ts=4 noet: */
