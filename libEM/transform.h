/**
 * $Id$
 */
#ifndef eman__transform_h__
#define eman__transform_h__ 1

#include "vec3.h"
#include "emobject.h"
//#include <vector>
//#include <map>
//#include <string>

using std::string;
using std::map;
using std::vector;

namespace EMAN
{
	/** @file transform.h
	 *   These are  a collection of transformation tools: rotation, translation,
	 *            and construction of symmetric objects
         *  @author Philip Baldwin and Steve Ludtke
	 *    <Philip.Baldwin@uth.tmc.edu>
	 *	Transform defines a transformation, which can be rotation,
         *         translation, scale, and their combinations.
	 *
	 *  @date $Date: 2005/04/04 17:41pm
	 *
	 *  @see Phil's article
	 *
	 * Internally a transformation is stored in a 4x4 matrix.
	 *         a b c d
	 *         e f g h           R     v
	 *  M=     j k m n    =      0     1    , where R is 3by3, v is 3by1
	 *         0 0 0 1
	 *  This is a standard computer graphics convention and can be found in many
	 *    references including Frackowiak et al; Human Brain Function
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
	 * The (post) translation is stored in (d, h, n).
	 *
	 *
	 * If rotations need to be found
	 *    around alternate origins, then brief calculations need to be performed
	 * Pre and Post Translations should be kept as separate vectors
	 *
	 * a matrix  R is called orthogonal if
	 *           R * transpose(R) = 1.
	 * All Real Orthogonal Matrices have eigenvalues with unit modulus and determinant
	 *  therefore equal to  "\pm 1"
	 *
     */
     /** @ingroup tested3c */
	class Transform3D
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
			SGIROT,
			SPIDER,
			MRC,
			XYZ,
			MATRIX
		};
		
	     // C1
		Transform3D();

	     // C2   
		Transform3D(float az,float alt, float phi); // EMAN by default


             //  C3  Usual Constructor: Post Trans, after appying Rot
		Transform3D(const Vec3f& posttrans, float az,float alt, float phi);
                

 	     // C4
	     	Transform3D(EulerType euler_type, float a1, float a2, float a3) ; // only EMAN: az alt phi
									      // SPIDER     phi theta psi
		
	     // C5   First apply pretrans: Then rotation
		Transform3D(EulerType euler_type, const Dict& rotation);
		

	     // C6   First apply pretrans: Then rotation: Then posttrans
		Transform3D(const Vec3f & pretrans, const Vec3f& posttrans, float az, float alt, float phi);

              
		Transform3D(float m11, float m12, float m13,
					float m21, float m22, float m23,
					float m31, float m32, float m33);

		virtual ~ Transform3D();  // COmega   Deconstructor

		void set_posttrans(const Vec3f & posttrans);
		void apply_scale(float scale);
		void set_scale(float scale);
		void orthogonalize();	// reorthogonalize the matrix
		void transpose();	// create the transpose in place

		void set_rotation(float az, float alt,float phi);
		void set_rotation(EulerType euler_type, float a1, float a2, float a3); // just SPIDER and EMAN
		
		void set_rotation(float m11, float m12, float m13,
                                  float m21, float m22, float m23,
			          float m31, float m32, float m33);

		void set_rotation(EulerType euler_type, const Dict &rotation );
		

		/** returns a rotation that maps a pair of unit vectors, a,b to a second  pair A,B
		 * @param eahat, ebhat, eAhat, eBhat are all unit vectors
		 * @return a transform3D rotation
		 */
		void set_rotation(const Vec3f & eahat, const Vec3f & ebhat,
		                                    const Vec3f & eAhat, const Vec3f & eBhat); 


		/** returns the magnitude of the rotation
		*/
		float get_mag() const;
		/** returns the spin-axis (or finger) of the rotation
		*/
		Vec3f get_finger() const;
		Vec3f get_posttrans(Vec3f &pretrans) const;
		Vec3f get_posttrans() const;
 		Vec3f get_center() const;
		Vec3f get_matrix3_col(int i) const;
		Vec3f get_matrix3_row(int i) const;
		
		Transform3D inverse();
					
		Dict get_rotation(EulerType euler_type=EMAN) const;

		void printme() const {
			for (int i=0; i<4; i++) {
				printf("%6.3f\t%6.3f\t%6.3f\t%6.3f\n",
					   matrix[i][0],matrix[i][1],matrix[i][2],matrix[i][3]);
			}
			printf("\n");
		}

		inline float at(int r,int c) { return matrix[r][c]; }
		float * operator[] (int i);
		const float * operator[] (int i) const;

            //
		static int get_nsym(const string & sym);
            Transform3D get_sym(const string & sym, int n) const;
		void set_center(const Vec3f & center);
		void set_pretrans(const Vec3f & pretrans);
            //
		float get_scale() const;   

		void to_identity();
        bool is_identity();

		/** Convert a list of euler angles to a vector of Transform3D objects.
		 *
		 *	@param[in] eulertype The type of Euler angles that is being passed in.
		 *	@param[in] angles A flat vector of angles.
		 *
		 *	@return Vector of pointers to Transform3D objects.
		 */
		static vector<Transform3D*>
			angles2tfvec(EulerType eulertype, const vector<float> angles);


	private:
		enum SymType
		{
			CSYM,
			DSYM,
			TET_SYM,
			ICOS_SYM,
			OCT_SYM,
			ISYM,
			UNKNOWN_SYM
		};

		void init();

		static SymType get_sym_type(const string & symname);

		float matrix[4][4];

		static map<string, int> symmetry_map;
	}; // ends Class

	Transform3D operator*(const Transform3D & M1, const Transform3D & M2);
	Vec3f operator*(const Vec3f & v    , const Transform3D & M);
	Vec3f operator*(const Transform3D & M, const Vec3f & v    );



}  // ends NameSpace EMAN


#endif


/* vim: set ts=4 noet: */
