/**
 * $Id$
 */
#ifndef eman__transform_h__
#define eman__transform_h__ 1

#include "vec3.h"
#include "emobject.h"
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
			SGIROT,
			SPIDER,
			MRC
		};
	public:
		Transform();
		
		Transform(EulerType euler_type, float a1,float a2,float a3, float a4=0);
		
		Transform(const Vec3f& posttrans, EulerType euler_type,
				  float a1,float a2, float a3, float a4=0);
		
		Transform(const Vec3f & pretrans, const Vec3f& posttrans, EulerType euler_type, 
				  float a1,float a2, float a3, float a4=0);
		
		Transform(EulerType euler_type, Dict& rotation);

		Transform(const Vec3f& posttrans, EulerType euler_type, Dict& rotation);
		
		Transform(const Vec3f & pretrans, const Vec3f& posttrans,   
				  EulerType euler_type, Dict& rotation);
		
		virtual ~ Transform();

		void to_identity();
		bool is_identity();
		
		void orthogonalize();	// reorthogonalize the matrix
		Transform inverse();
	
		void set_pretrans(const Vec3f & pretrans);
		void set_posttrans(const Vec3f & posttrans);
		void set_center(const Vec3f & center);
		void set_rotation(EulerType euler_type, float a0, float a1,float a2, float a3=0);
		void set_rotation(EulerType euler_type, Dict &rotation );
		void set_scale(float scale);

		Vec3f get_pretrans() const;
		Vec3f get_posttrans() const;
		Vec3f get_center() const;
		Dict get_rotation(EulerType euler_type) const;
		Vec3f get_matrix3_col(int i) const;
		Vec3f get_matrix3_row(int i) const;
		float get_scale() const;

		// returns orthogonality coefficient 0-1 range;
		float orthogonality() const;
		
		static int get_nsym(const string & sym);
		Transform get_sym(const string & sym, int n);
		
		// A few utilities used for matrix operations
		inline void row_mult(int n,float f)
		{
			matrix[n][0]*=f;
			matrix[n][1]*=f;
			matrix[n][2]*=f;
		}
		
		inline void row_add(int n,float a,float b,float c)
		{
			matrix[n][0]+=a;
			matrix[n][1]+=b;
			matrix[n][2]+=c;
		}

		void row_exch(int n1,int n2) { 
			float t; 
			t=matrix[n1][0]; matrix[n1][0]=matrix[n2][0]; matrix[n2][0]=t; 
			t=matrix[n1][1]; matrix[n1][1]=matrix[n2][1]; matrix[n2][1]=t; 
			t=matrix[n1][2]; matrix[n1][2]=matrix[n2][2]; matrix[n2][2]=t; 
		}
		void print() {
			for (int i=0; i<3; i++) {
				printf("%6.3f\t%6.3f\t%6.3f\n",
					   matrix[i][0],matrix[i][1],matrix[i][2]);
			}
			printf("\n");
		}
				
		inline float at(int r,int c) { return matrix[r][c]; }
		float * operator[] (int i);
		const float * operator[] (int i) const;

		Transform operator*=(const Transform& t);
	
		
	  private:
		enum SymType
		{
			CSYM,
			DSYM,
			ICOS_SYM,
			OCT_SYM,
			ISYM,
			UNKNOWN_SYM
		};

		void init();

		float eman_alt() const;
		float eman_az() const;
		float eman_phi() const;
		
		vector<float> matrix2quaternion() const;
		vector<float> matrix2sgi() const;
		void quaternion2matrix(float e0, float e1, float e2, float e3);
		static SymType get_sym_type(const string & symname);
		
		float matrix[4][3];
		Vec3f pre_trans;

		static map<string, int> symmetry_map;
	};

	Vec3f operator*(const Vec3f & v, const Transform & t);
	Transform operator*(const Transform & t1, const Transform & t2);
	
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
