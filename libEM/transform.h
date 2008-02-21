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
	 *         e f g h           R        v
	 *  M=     j k m n    =      vpre     1    , where R is 3by3, v is 3by1
	 *         p q r 1
	 *  The standard computer graphics convention is identical to ours after setting vpre to 
	 *     zero and can be found in many
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
	 * The cumulative translation is stored in (d, h, n).
	 * We put the post-translation into (p, q, r), since it is convenient
	 *   to carry along at times. When matrices are multiplied or printed, these
	 *   are hidden to the user. They can only be found by applying the post_translation
	 *    method, and these elements are non-zero. Otherwise the post-translation method returns
	 *         the cumulative translationmlb
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

             // copy constructor
	    Transform3D( const Transform3D& rhs );

	     // C2   
		Transform3D(float az,float alt, float phi); // EMAN by default


             //  C3  Usual Constructor: Post Trans, after appying Rot
		Transform3D(float az, float alt, float phi, const Vec3f& posttrans);
 
 	     // C4
		Transform3D(EulerType euler_type, float a1, float a2, float a3) ; // only EMAN: az alt phi
								                            // SPIDER     phi theta psi
		
	     // C5   First apply pretrans: Then rotation
		Transform3D(EulerType euler_type, const Dict& rotation);
		

	     // C6   First apply pretrans: Then rotation: Then posttrans
		Transform3D(const Vec3f & pretrans, float az, float alt, float phi,  const Vec3f& posttrans);

		Transform3D(float m11, float m12, float m13,
					float m21, float m22, float m23,
					float m31, float m32, float m33);

		~ Transform3D();  // COmega   Deconstructor


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
		 * @return  a transform3D rotation
		 */
		void set_rotation(const Vec3f & eahat, const Vec3f & ebhat,
		                                    const Vec3f & eAhat, const Vec3f & eBhat); 


		/** returns the magnitude of the rotation
		*/
		float get_mag() const;
		/** returns the spin-axis (or finger) of the rotation
		*/
		Vec3f get_finger() const;
		Vec3f get_pretrans( int flag=0) const; // flag=1 => all trans is pre
		Vec3f get_posttrans(int flag=0) const; // flag=1 => all trans is post
 		Vec3f get_center() const;
		Vec3f get_matrix3_col(int i) const;
		Vec3f get_matrix3_row(int i) const;
        Vec3f transform(Vec3f & v3f) const ; // This applies the full tranform to the vec
        Vec3f rotate(Vec3f & v3f) const ;  // This just applies the rotation to the vec
		
		Transform3D inverseUsingAngs() const;
		Transform3D inverse() const;
					
		Dict get_rotation(EulerType euler_type=EMAN) const;

		void printme() const {
			for (int i=0; i<3; i++) {
				printf("%6.15f\t%6.15f\t%6.15f\t%6.1f\n",
					   matrix[i][0],matrix[i][1],matrix[i][2],matrix[i][3]);
			}
			printf("%6.3f\t%6.3f\t%6.3f\t%6.3f\n",0.0,0.0,0.0,1.0);
			printf("\n");
		}

		inline float at(int r,int c) const { return matrix[r][c]; }
		void set(int r, int c, float value) { matrix[r][c] = value; }
		float * operator[] (int i);
		const float * operator[] (int i) const;
		
		static int get_nsym(const string & sym);
		Transform3D get_sym(const string & sym, int n) const;
		void set_center(const Vec3f & center);
		void set_pretrans(const Vec3f & pretrans); // flag=1 means count all translation as pre
		void set_pretrans(float dx, float dy, float dz);
		void set_pretrans(float dx, float dy);
		void set_posttrans(const Vec3f & posttrans);// flag=1 means count all translation as post
		void set_posttrans(float dx, float dy, float dz);
		void set_posttrans(float dx, float dy);

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

		void dsaw_zero_hack(){
			for (int j=0; j<4; ++j) {
				for (int i=0; i<4; i++) {
				if ( fabs(matrix[j][i]) < 0.000001 )
					matrix[j][i] = 0.0;
				}
			}
			
		}

	private:
		enum SymType
		{      CSYM,
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

	class FactoryBase
	{
	public:
		FactoryBase() {}
		virtual ~FactoryBase() {};
		
		// Factory dependent functionality
		virtual string get_name() const = 0;

		virtual string get_desc() const = 0;
		
		Dict get_params() const	{ return params; }

		void set_params(const Dict & new_params)
		{
			// This commented out by d.woolford. 
			//params.clear();
			insert_params(new_params);
		}
		
		virtual TypeDict get_param_types() const = 0;
		
		void insert_params(const Dict & new_params)
		{
			// this is really inserting OR individually replacing...
			// the old data will be kept if it is not written over
			TypeDict permissable_params = get_param_types();
			for ( Dict::const_iterator it = new_params.begin(); it != new_params.end(); ++it )
			{
			
				if ( !permissable_params.find_type(it->first) )
				{
					throw InvalidParameterException(it->first);
				}
				params[it->first] = it->second;
			}	
		}
		// end factory dependent functionality
		
		
		protected:
			mutable Dict params;
	};
	
	/** Symmetry3D 
	* an abstract (meaning it has pure virtual functions) base class for Symmetry3D objects
	* Objects of this type must provide delimiters for the asymmetric unit (get_delimiters), and
	* must also provide all of the rotational symmetric operations (get_sym(int n)).
	* get_delimiter returns a dictionary with "alt_max" and "az_max" keys, which correspond to the
	* encompassing azimuth and altitude angles of the asymmetric unit. These can be interpreted in a
	* relatively straight forward fashion when dealing with C and D symmetries to demarcate the asymmetric
	* unit, however when dealing with Platonic symmetries the story is not so basic.... to be continued
	*
	* see http://blake.bcm.edu/emanwiki/EMAN2/Symmetry for figures and description of what we're doing
	* here, for all the symmetries
	* @author David Woolford
	* @date Feb 2008
	*/
	class Symmetry3D : public FactoryBase
	{
		public:
		Symmetry3D() {};
		virtual  ~Symmetry3D() {};
		
		virtual TypeDict get_param_types() const
		{
			// NOTE - child classes should call
			// TypeDict d = Symmetry3D::get_params_types();
			// to initialize their type dict before inserting their own parameters
			TypeDict d;
			d.put("inc_mirror", EMObject::BOOL, "Include mirror portion of asymmetric unit. Default faulse.");
			return d;
		}
		// end factory dependent functionality
		
		// Symmetry virtual behavior
		virtual Dict get_delimiters() const = 0;
		
		virtual Transform3D get_sym(int n) const = 0;
		
		// This is a hack, because this functionality is only relevant to platonic symmetries. But it could
		// grow into functionality for the other symmetries.
		virtual bool is_in_asym_unit(const float& altitude, const float& azimuth) const { return true; }
		// This is a hack, because this functionality is only relevant to platonic symmetries. But it could
		// grow into functionality for the other symmetries.
		virtual float get_az_alignment_offset() const { return 0.0; }
		
		virtual float get_h(const float& prop,const float& altitude) const = 0;
		
		virtual float get_h_base(const float& prop,const float& altitude, const int maxcsym) const;
		
		virtual int get_nsym() const = 0;
		
		virtual bool is_platonic() const { return false; }
		
		bool is_h_sym() const;
		
		vector<Transform3D> gen_orientations(const string& generatorname="standard", const Dict& parms=Dict());

	};
	
	class CSym : public Symmetry3D
	{
		public:
		CSym() {};
		virtual  ~CSym() {};
		
		static Symmetry3D *NEW()
		{
			return new CSym();
		}
		
		virtual string get_name() const { return NAME; }

		virtual string get_desc() const { return "C symmetry support"; }
		
		virtual TypeDict get_param_types() const
		{
			TypeDict d = Symmetry3D::get_param_types();
			d.put("nsym", EMObject::INT, "The symmetry number");
			return d;
		}
		
		virtual Dict get_delimiters() const;
		
		virtual Transform3D get_sym(int n) const;
		
		virtual float get_h(const float& prop,const float& altitude) const;
		
		virtual int get_nsym() const { return params["nsym"]; };
		
		static const string NAME;
	};
	
	class DSym : public Symmetry3D
	{
		public:
		DSym() {};
		virtual  ~DSym() {};
		
		static Symmetry3D *NEW()
		{
			return new DSym();
		}
		
		virtual string get_name() const { return NAME; }

		virtual string get_desc() const { return "D symmetry support"; }
		
		virtual TypeDict get_param_types() const
		{
			TypeDict d = Symmetry3D::get_param_types();
			d.put("nsym", EMObject::INT, "The symmetry number");
			return d;
		}
		
		virtual Dict get_delimiters() const;
		
		virtual Transform3D get_sym(int n) const;
		
		virtual float get_h(const float& prop,const float& altitude) const;
		
		
		virtual int get_nsym() const { return 2*(int)params["nsym"]; };
		
		static const string NAME;
	};
	
	class HSym : public Symmetry3D
	{
		public:
			HSym() {};
			virtual  ~HSym() {};
		
			static Symmetry3D *NEW()
			{
				return new HSym();
			}
		
			virtual string get_name() const { return NAME; }

			virtual string get_desc() const { return "H symmetry support"; }
		
			virtual TypeDict get_param_types() const
			{
				TypeDict d = Symmetry3D::get_param_types();
				d.put("nsym", EMObject::INT, "The symmetry number");
				d.put("equator_range", EMObject::FLOAT, "The amount altitude angles are allowed to vary above and below the equator. Default is 5");
				d.put("dz", EMObject::FLOAT, "The translational distance (along z) between succesive identical subunits in angstrom (default a/pix is 1)");
				d.put("daz", EMObject::FLOAT, "The rotational angle (about z) between successive identical subunits in degrees");
				d.put("apix", EMObject::FLOAT, "Angstrom per pixel, default is one.");
				return d;
			}
		
			virtual Dict get_delimiters() const;
		
			virtual Transform3D get_sym(int n) const;
		
			virtual float get_h(const float& prop,const float& altitude) const;
		
		
			virtual int get_nsym() const { return (int)params["nsym"]; };
		
			static const string NAME;
	};
	
	/**
	* see www.math.utah.edu/~alfeld/math/polyhedra/polyhedra.html for pictures of platonic solids
	* By default we will put largest symmetry along z-axis.
	*
	*
	*
	*

    * The ICOS symmetry group has the face along z-axis

	* IMPORTANT NOTE is anything that derives from this class must call init() in its constructor
	*/
	class PlatonicSym : public Symmetry3D
	{
		public:
		PlatonicSym() {};
		virtual  ~PlatonicSym() {};
		
		virtual string get_name() const = 0;

		virtual string get_desc() const = 0;
		
		virtual Dict get_delimiters() const;
		
		virtual Transform3D get_sym(int n) const = 0;
		
		virtual float get_h(const float& prop,const float& altitude) const;
		
		// A virtual function particular to PlatonicSym
		virtual int get_max_csym() const = 0;
		
		// This function works for icosahedral and octahedral, but not tetrahedral
		virtual bool is_in_asym_unit(const float& altitude, const float& azimuth) const;
		
		virtual bool is_platonic() const { return true; }
		protected:
		Dict platonic_params;

		// Called to initialize platonic_params, should be called in the constructor of all
		// Platonic solids that derive from this
		void init();
		
		// Returns the lower bound of the asymmetric unit, as dependent on azimuth, and on alpha -
		// alpha is alt_max for icos and oct, but may be alt_max/2.0 for tet depending on mirror
		// symmetry etc
		float platonic_alt_lower_bound(const float& azimuth, const float& alpha) const;
			
	};
	
	class TetrahedralSym : public PlatonicSym
	{
		public:
		TetrahedralSym()  {init();}
		virtual  ~TetrahedralSym() {}
		
		static Symmetry3D *NEW()
		{
			return new TetrahedralSym();
		}
		
		virtual string get_name() const { return NAME; };

		virtual string get_desc() const { return "Tetrahedral symmetry support"; }

		virtual int get_max_csym() const { return 3; }
		
		virtual Transform3D get_sym(int n) const;
		
		virtual bool is_in_asym_unit(const float& altitude, const float& azimuth) const;
		
		virtual int get_nsym() const { return 12; };
		
		virtual float get_az_alignment_offset() const { return  180.0; }
		
		static const string NAME;
		
		
	};
	
	class OctahedralSym : public PlatonicSym
	{
		public:
		OctahedralSym()  {init();}
		virtual  ~OctahedralSym() {}
		
		static Symmetry3D *NEW()
		{
			return new OctahedralSym();
		}
		
		virtual string get_name() const { return NAME; };

		virtual string get_desc() const { return "Octahedral symmetry support"; }

		virtual int get_max_csym() const { return 4; }
		
		virtual Transform3D get_sym(int n) const;
		
		virtual int get_nsym() const { return 24; };
		
		static const string NAME;
	};
	
	/**
	* Doctor Phil has this to say about icosahedral symmetry:
	* " Each Platonic Solid has 2E symmetry elements.
	* An icosahedron has   m=5, n=3, F=20 E=30=nF/2, V=12=nF/m,since vertices shared by 5 triangles
	* It is composed of 20 triangles. E=3*20/2
	* An dodecahedron has m=3, n=5   F=12 E=30  V=20
	* It is composed of 12 pentagons. E=5*12/2;   V= 5*12/3, since vertices shared by 3 pentagons"
	*
	*
	*/
	class IcosahedralSym : public PlatonicSym
	{
		public:
		IcosahedralSym() {init(); }
		virtual  ~IcosahedralSym() { }
		
		static Symmetry3D *NEW()
		{
			return new IcosahedralSym();
		}
			
		virtual string get_name() const { return NAME; };

		virtual string get_desc() const { return "Icosahedral symmetry support"; }

		virtual int get_max_csym() const { return 5; }// 5 is the greatest symmetry
		
		virtual Transform3D get_sym(int n) const;
		
		virtual int get_nsym() const { return 60; };
		
		virtual float get_az_alignment_offset() const { return  234.0; /*270.0-36.0;*/ }
		
		static const string NAME;
	};
	
	template <> Factory < Symmetry3D >::Factory();
	void dump_symmetries();
	map<string, vector<string> > dump_symmetries_list();
	
	class SymOrientationGenerator : public FactoryBase
	{
	public:
		SymOrientationGenerator() {};
		virtual ~SymOrientationGenerator() {};
		
		virtual vector<Transform3D> gen_orientations(const Symmetry3D* const)  = 0;
		
	};
	
	class AsymmUnitCoverer : public SymOrientationGenerator
	{
	public:
		AsymmUnitCoverer() {};
		virtual  ~AsymmUnitCoverer() {};
		
		static SymOrientationGenerator *NEW()
		{
			return new AsymmUnitCoverer();
		}
		
		virtual string get_name() const { return "asymmunit"; }

		virtual string get_desc() const { return "Generate projections distributed quasi-uniformaly over the asymmetric unit"; }
		
		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("prop", EMObject::FLOAT, "The angular separation of projections");
			d.put("perturb", EMObject::BOOL, "Whether or not to perturb the generated orientations, default is false");
			d.put("smear", EMObject::BOOL, "Whether or not smear the projections in plane");
			d.put("phitoo", EMObject::FLOAT, "Specifying a non zero value for this argument will cause phi rotations to be included");
			return d;
		}
		
		virtual vector<Transform3D> gen_orientations(const Symmetry3D* const);
		
	};
	
	template <> Factory < SymOrientationGenerator >::Factory();
}  // ends NameSpace EMAN



#endif


/* vim: set ts=4 noet: */
