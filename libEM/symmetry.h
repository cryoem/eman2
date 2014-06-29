
/*
 * Author: David Woolford, 09/23/2008 (woolford@bcm.edu)
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
#ifndef eman__symmetry_h__
#define eman__symmetry_h__ 1


#include "emobject.h"
#include "vec3.h"
#include "transform.h"

namespace EMAN {

	/** Symmetry3D - A base class for 3D Symmetry objects.
	* Objects of this type must provide delimiters for the asymmetric unit (get_delimiters), and
	* must also provide all of the rotational symmetric operations (get_sym(const int n)). They must also
	* provide the total number of unique symmetric operations with get_nsym (except in helical symmetry).
	* get_delimiter returns a dictionary with "alt_max" and "az_max" keys, which correspond to the
	* encompassing azimuth and altitude angles of the asymmetric unit. These can be interpreted in a
	* relatively straight forward fashion when dealing with C and D symmetries to demarcate the asymmetric
	* unit, however when dealing with Platonic symmetries the asymmetric unit is not so trivial.
	* see http://blake.bcm.edu/emanwiki/EMAN2/Symmetry for figures and description of what we're doing
	* here, for all the symmetries, and look in the comments of the PlatonicSym classes themselves.
	* It inherits from a factory base, making it amenable to incorporation in EMAN2 style factories
	* @author David Woolford with Philip Baldwin and Steven Ludtke
	* @date Feb 2008
	 */
	class Symmetry3D : public FactoryBase
	{
	public:
		typedef vector<vector<Vec3f> >::const_iterator cit;
		typedef vector<vector<Vec3f> >::iterator ncit;
		Symmetry3D();
		virtual  ~Symmetry3D();

		/** Every Symmetry3D object must return a dictionary containing the delimiters
		 * that define its asymmetric unit (this is not strictly true in the case of the PlatonicSym class)
		 * @param inc_mirror whether or not the mirror part of the asymmetric unit should be included in the consideration
		 * @return a dictionary containing atleast "alt_max" and "az_max"
		 */
		virtual Dict get_delimiters(const bool inc_mirror=false) const = 0;

		/** Every Symmetry3D object must provide access to the full set of its symmetry operators
		 * via this function
		 * @param n the symmetry operator number
		 * @return a Transform object describing the symmetry operation
		 */
		virtual Transform get_sym(const int n) const = 0;

		/** The total number of unique symmetry operations that will be return by this object when
		 * a calling program access Symmetry3D::get_sym. 
		 */
		virtual int get_nsym() const = 0;

		/**This functionality is only relevant to platonic symmetries. But it could
		 * grow into functionality for the other symmetries.
		 */
		virtual float get_az_alignment_offset() const { return 0.0; }



		/** A function that is used to determine if this is a platonic symmetry object
		 * This function is only virtually overidden by the PlatonicSym symmetry, which returns true, not false
		 * @return false - indicating that this is not a platonic symmetry object
		 */
		virtual bool is_platonic_sym() const { return false; }

		/** A function that is used to determine if this is a Helical symmetry object
		 * This function is only virtually overidden by the HSym symmetry, which returns true, not false
		 * @return false - indicating that this is not a helical symmetry object
		 */
		virtual bool is_h_sym() const { return false; }

		/** A function that is used to determine if this is a c symmetry object
		 * This function is only virtually overidden by the CSym object, which returns true
		 * @return false - indicating that this is not a helical symmetry object
		 */
		virtual bool is_c_sym() const { return false; }

		/** A function that is used to determine if this is a d symmetry object
		 * This function is only virtually overidden by the DSym object, which returns true
		 * @return false - indicating that this is not a helical symmetry object
		 */
		virtual bool is_d_sym() const { return false; }

		/** A function that is used to determine if this is the tetrahedral symmetry object
		 * This function is only virtually overidden by the TetSym object, which returns true
		 * @return false - indicating that this is not a tetrahedral symmetry object
		 */
		virtual bool is_tet_sym() const { return false; }



		/** The Symmetry3D object must return the maximum degree of symmetry it exhibits about any one axis.
		 * This function is only called in the AsymmUnitOrientationGenerator.
		 */
		virtual int get_max_csym() const = 0;

		/** The Symmetry3D object must be capable of returning an ordered list of points on the unit
		 * sphere that define its asymmetric unit (with mirror considerations). The list can
		 * be any length, and must be constructed carefully. If the list consists of points A B and C,
		 * then arcs on the unit sphere connecting A to B, then B to C, then C to A must define the
		 * asymmetric unit (with or without its mirror portion). i.e. it is a cyclic list, on any
		 * length
		 * @param inc_mirror whether or not to include the mirror portion of the asymmetric unit
		 * @return a vector or points which define a cyclic set of great arcs on the unit sphere. Each point may be connected to the point that proceeds it, and the last point may be connected to the first point, and this demarcates the asymmetric unit.
		 */
		virtual vector<Vec3f> get_asym_unit_points(bool inc_mirror) const = 0;
		/** Ask the Symmetry3D object to generate a set of orientations in its asymmetric unit
		 * using an OrientationGenerator constructed from the given parameters (using a Factory).
		 * This is reminiscent of the strategy design pattern
		 * @param generatorname the string name of the OrientationGenerator, as accessed for the OrientationGenerator factory
		 * @param parms the parameters handed to OrientationGenerator::set_params after initial construction
		 * @return a set of orientations in the unit sphere
		 */
		vector<Transform> gen_orientations(const string& generatorname="eman", const Dict& parms=Dict());

		/** A function to be used when generating orientations over portion of the unit sphere
		 * defined by parameters returned by get_delimiters. In platonic symmetry altitude and azimuth
		 * alone are not enough to correctly demarcate the asymmetric unit. See the get_delimiters comments.
		 * @param altitude the EMAN style altitude of the 3D orientation in degrees
		 * @param azimuth the EMAN style azimuth of the 3D orientation in degrees
		 * @param inc_mirror whether or not to include orientations if they are in the mirror portion of the asymmetric unit
		 * @return true or false, depending on whether or not the orientation is within the asymmetric unit
		 */
		virtual bool is_in_asym_unit(const float& altitude, const float& azimuth, const bool inc_mirror) const = 0;

		/** A function that will reduce an orientation, as characterized by Euler anges, into a specific asymmetric unit.
		 * Default behavior is to map the given orientation into the default asymmetric unit of the symmetry (n=0). This
		 * is a concrete implementation that works for all symmetries, relying on a concrete instance of the get_asym_unit_triangles
		 * function
		 * @param t3d a Transform characterizing an orientation
		 * @param n the number of the asymmetric unit you wish to map the given orientation into. There is a strong relationship between n and to Symmetry3D::get_sym
		 * @return the orientation the specified asymmetric unit (by default this is the default asymmetric unit of the symmetry)
		 * @ingroup tested3c
		 */
		virtual Transform reduce(const Transform& t3d, int n=0) const;


		/** A function that will determine in which asymmetric unit a given orientation resides
		 * The asymmetric unit 'number' will depend entirely on the order in which different symmetry operations are return by the
		 * Symmetry3D::get_sym function
		 * @param t3d a Transform characterizing an orientation
		 * @return the asymmetric unit number the the orientation is in
		 */
		virtual int in_which_asym_unit(const Transform& t3d) const;

		/** A function that will determine in which asymmetric unit a given vector resides
		 * The asymmetric unit 'number' will depend entirely on the order in which different symmetry operations are return by the
		 * Symmetry3D::get_sym function
		 * The vector is a point
		 * @param v a Vec3f characterizing a point
		 * @return the asymmetric unit number the the orientation is in
		 */
		virtual int point_in_which_asym_unit(const Vec3f& v) const;

		/** Get triangles that precisely occlude the projection area of the default asymmetric unit. This will be used
		 * for collision detection in Symmetry3D::reduce
		 * @param inc_mirror whether to include the mirror portion of the asymmetric unit
		 */
		virtual vector<vector<Vec3f> > get_asym_unit_triangles(bool inc_mirror) const = 0;

		/** Gets a vector of Transform objects that define the set of asymmetric units that touch the default
		 * asymmetric unit. The 'default asymmetric unit' is defined by the results of Symmetry3d::get_asym_unit_points
		 * and is sensitive to whether or not you want to include the mirror part of the asymmetric unit.
		 * This function is useful when used in conjunction with Symmetry3D::reduce, and particularly when finding
		 * the angular deviation of particles through different stages of iterative Single Particle Reconstruction
		 * This function could be expanded to work for an asymmetric unit number supplied by the user.
		 * @param inc_mirror whether or not to include the mirror portion of the asymmetric unit
		 * @return a vector of Transform objects that map the default asymmetric unit to the neighboring asymmetric unit
		 */
		virtual vector<Transform> get_touching_au_transforms(bool inc_mirror = true) const;

		virtual vector<Transform> get_syms() const;
		static vector<Transform> get_symmetries(const string& symmetry);
	protected:
		/// The asymmetric unit planes are cached to provide a great speed up
		/// the point_in_which_asym_unit function, which is called by reduce and by in_which_asym_unit
		mutable float** cached_au_planes;

		/// Have to remember the cache size
		mutable int cache_size;
		/// This is stores the number of triangles returned by get_asym_unit_triangles(true)
		mutable int num_triangles;
		/// This cache is of size cache_size
		mutable vector< vector<Vec3f> > au_sym_triangles;
		/** Establish the asymmetric unit planes cache
		*/
		void cache_au_planes() const;

		/** Clear the asymmetric unit planes cache
		*/
		void delete_au_planes();
	private:
		/** Disallow copy construction */
		Symmetry3D(const Symmetry3D&);
		/** Disallow assignment */
		Symmetry3D& operator=(const Symmetry3D&);
};

	/** An encapsulation of cyclic 3D symmetry
 * @author David Woolford (based on previous work by Phil Baldwin and Steve Ludtke)
 * @date Feb 2008
	 */
	class CSym : public Symmetry3D
{
	public:
		CSym() {};
		virtual  ~CSym() {};

		/** Factory support function NEW
		 * @return a newly instantiated class of this type
		 */
		static Symmetry3D *NEW()
		{
			return new CSym();
		}

		/** Return CSym::NAME
		 * @return the unique name of this class
		 */
		virtual string get_name() const { return NAME; }

		/** Get a description
		 * @return a clear desciption of this class
		 */
		virtual string get_desc() const { return "C symmetry support"; }

		/** Get a dictionary containing the permissable parameters of this class
		 * @return a dictionary containing the permissable parameters of this class
		 */
		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("nsym", EMObject::INT, "The symmetry number");
			return d;
		}

		/** Get the altitude and phi angle of the c symmetry, which depends on nysm.
		 * The "alt_max" value in the return dicts is 180 or 90 degrees, depending inc_mirror
		 * The "az_max" is 360/nsym degrees.
		 * @param inc_mirror whether or not to include the part of the asymmetric unit which contains the mirror projections of the other half
		 * @return a dictionary containing the keys "alt_max" and "az_max"
		 * @exception InvalidValueException if nsym is less than or equal to zero
		 */
		virtual Dict get_delimiters(const bool inc_mirror=false) const;

		/** Provides access to the complete set of rotational symmetry operations associated with this symmetry.
		 * Rotational symmetry operations for C symmetry are always about the z-axis (in the EMAN convention), and
		 * therefore the only non zero return angle is azimuth. Specifically, it is n*360/nsym degrees.
		 * @param n the rotational symmetry operation number. If n is greater than nsym we take n modulo nsym
		 * @return a transform containing the correct rotational symmetric operation.
		 * @exception InvalidValueException if nsym is less than or equal to zero
		 */
		virtual Transform get_sym(const int n) const;

		/** Gets the total number of unique roational symmetry operations associated with this symmetry
		 * For C symmetry, this is simply nsym
		 * @return the degree of of cyclic symmetry (nsym)
		 */
		virtual int get_nsym() const { return params["nsym"]; };


		/** Gets the maximum symmetry of this object. This is used by OrientationGenerators, and is
		 * probably not something a general user would utilize.
		 * @return the degree of of cyclic symmetry (nsym) - this is the maximum symmetry
		 */
		virtual int get_max_csym() const { return params["nsym"]; }

		/// The name of this class - used to access it from factories etc. Should be "c"
		static const string NAME;

		/** to demarcate the asymmetric unit. The last should may be connected to the first.
		 * @param inc_mirror whether or not to include the mirror portion of the asymmetric unit
		 * @return a cyclic set of points which can be connected using great arcs on the unit sphere
		 */
		virtual vector<Vec3f> get_asym_unit_points(bool inc_mirror = false) const;

		/** A function to be used when generating orientations over portion of the unit sphere
		 * defined by parameters returned by get_delimiters. In platonic symmetry altitude and azimuth
		 * alone are not enough to correctly demarcate the asymmetric unit. See the get_delimiters comments.
		 * @param altitude the EMAN style altitude of the 3D orientation in degrees
		 * @param azimuth the EMAN style azimuth of the 3D orientation in degrees
		 * @param inc_mirror whether or not to include orientations if they are in the mirror portion of the asymmetric unit
		 * @return true or false, depending on whether or not the orientation is within the asymmetric unit
		 */
		virtual bool is_in_asym_unit(const float& altitude, const float& azimuth, const bool inc_mirror) const;

		/** Returns true - this is indeed a c symmetry object
		 * @return true - indicating that this is a c symmetry object
		 */
		virtual bool is_c_sym() const { return  true; }


		/** Get triangles that precisely occlude the projection area of the default asymmetric unit. This is used
		 * for collision detection in Symmetry3D::reduce
		 * @param inc_mirror whether to include the mirror portion of the asymmetric unit
		 */
		virtual vector<vector<Vec3f> > get_asym_unit_triangles(bool inc_mirror) const;
	private:
		/** Disallow copy construction */
		CSym(const CSym&);
		/** Disallow assignment */
		CSym& operator=(const CSym&);

};

	/** An encapsulation of dihedral 3D symmetry
 * @author David Woolford (based on previous work by Phil Baldwin and Steve Ludtke)
 * @date Feb 2008
	 */
	class DSym : public Symmetry3D
{
	public:
		DSym() {};
		virtual  ~DSym() {};

		/** Factory support function NEW
		 * @return a newly instantiated class of this type
		 */
		static Symmetry3D *NEW()
		{
			return new DSym();
		}

		/** Return DSym::NAME
		 * @return the unique name of this class
		 */
		virtual string get_name() const { return NAME; }

		/** Get a description
		 * @return a clear desciption of this class
		 */
		virtual string get_desc() const { return "D symmetry support"; }

		/** Get a dictionary containing the permissable parameters of this class
		 * @return a dictionary containing the permissable parameters of this class
		 */
		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("nsym", EMObject::INT, "The symmetry number");
			return d;
		}

		/** Get the altitude and phi angle of the d symmetry, which depends on nysm.
		 * The "alt_max" is always 90 degrees
		 * The "az_max" is 360/nsym degrees of 180/nsym, depending the inc_mirror argument
		 * @param inc_mirror whether or not to include the part of the asymmetric unit which contains the mirror projections of the other half
		 * @return a dictionary containing the keys "alt_max" and "az_max"
		 * @exception InvalidValueException if nsym is less than or equal to zero
		 */
		virtual Dict get_delimiters(const bool inc_mirror=false) const;

		/** Provides access to the complete set of rotational symmetry operations associated with this symmetry.
		 * The first half symmetry operations returned by this function are all about the z axis (i.e. only azimuth
		 * is non zero. The second half of the symmetry operations are replicas of the first half, except that they
		 * have an additional 180 degree rotation about x (in EMAN terms, the altitude angle is 180).
		 * @param n the rotational symmetry operation number. If n is greater than nsym we take n modulo nsym
		 * @return a transform containing the correct rotational symmetric operation.
		 * @exception InvalidValueException if nsym is less than or equal to zero
		 */
		virtual Transform get_sym(const int n) const;

		/** Gets the total number of unique roational symmetry operations associated with this symmetry
		 * For D symmetry, this is simply 2*nsym
		 * @return two times nsym
		 */
		virtual int get_nsym() const { return 2*(int)params["nsym"]; };


		/** Gets the maximum symmetry of this object. This is used by OrientationGenerators, and is
		 * probably not something a general user would utilize.
		 * @return nsym - this is the maximum symmetry about a given any axis for D symmetry
		 */
		virtual int get_max_csym() const { return params["nsym"]; }

		/** @param inc_mirror whether or not to include the mirror portion of the asymmetric unit
		 * @return a cyclic set of points which can be connected using great arcs on the unit sphere
		 * to demarcate the asymmetric unit. The last should may be connected to the first.
		 */
		virtual vector<Vec3f> get_asym_unit_points(bool inc_mirror = false) const;

		/** A function to be used when generating orientations over portion of the unit sphere
		 * defined by parameters returned by get_delimiters. In platonic symmetry altitude and azimuth
		 * alone are not enough to correctly demarcate the asymmetric unit. See the get_delimiters comments.
		 * @param altitude the EMAN style altitude of the 3D orientation in degrees
		 * @param azimuth the EMAN style azimuth of the 3D orientation in degrees
		 * @param inc_mirror whether or not to include orientations if they are in the mirror portion of the asymmetric unit
		 * @return true or false, depending on whether or not the orientation is within the asymmetric unit
		 */
		virtual bool is_in_asym_unit(const float& altitude, const float& azimuth, const bool inc_mirror) const;

		/// The name of this class - used to access it from factories etc. Should be "d"
		static const string NAME;

		/** Returns true - this is indeed a c symmetry object
		 * @return true - indicating that this is a c symmetry object
		 */
		virtual bool is_d_sym() const { return  true; }

		/** Get triangles that precisely occlude the projection area of the default asymmetric unit. This is used
		 * for collision detection in Symmetry3D::reduce
		 * @param inc_mirror whether to include the mirror portion of the asymmetric unit
		 */
		virtual vector<vector<Vec3f> > get_asym_unit_triangles(bool inc_mirror) const;
	private:
		/** Disallow copy construction */
		DSym(const DSym&);
		/** Disallow assignment */
		DSym& operator=(const DSym&);
};

	/** An encapsulation of helical 3D symmetry
 * @author David Woolford (based on previous work by Phil Baldwin and Steve Ludtke)
 * @date Feb 2008
	 */
	class HSym : public Symmetry3D
{
	public:
		HSym() {};
		virtual  ~HSym() {};

			/** Factory support function NEW
		 * @return a newly instantiated class of this type
			 */
		static Symmetry3D *NEW()
		{
			return new HSym();
		}

			/** Return HSym::NAME
		 * @return the unique name of this class
			 */
		virtual string get_name() const { return NAME; }

			/** Get a description
		 * @return a clear desciption of this class
			 */
		virtual string get_desc() const { return "Helical symmetry, with support for N-start, pitch and limited tilt range. Specify as H<nsym>:<nstart>:<daz>:<tz in pix>[:<maxtilt>]"; }

			/** Get a dictionary containing the permissable parameters of this class
		 * Of all the symmetries, helical has the most options. This is because
		 * different approaches have to taken in regard to defining an asymmetric unit
		 * and to returning the set of rotational and translational symmetry operations
		 * @return a dictionary containing the permissable parameters of this class
			 */
		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("nsym", EMObject::INT, "The number of asymmetric units to generate. This could be infinite for helical symmetry. Normally a multiple of nstart.");
			d.put("nstart", EMObject::INT, "The Cn symmetry of a single Z-slice of the helix.");
			d.put("tz", EMObject::FLOAT, "The translational distance (along z) between successive identical subunits in angstroms (default A/pix is 1)");
			d.put("daz", EMObject::FLOAT, "The rotational angle (about z) between successive identical subunits in degrees");
			d.put("apix", EMObject::FLOAT, "Angstroms per pixel, default is 1.0, used only for tz");
			d.put("maxtilt", EMObject::FLOAT, "When generating projections, normally only 'side views' are created (3-D Z along Y in 2-D). This is the maximum out of plane tilt in degrees.");
			return d;
		}

			/** Get the altitude and phi angle of the d symmetry, which depends on nysm.
		 * The "alt_max" is always 90
		 * The "alt_min" 90-maxtilt
		 * The "az_max" is always 360/nsym degrees
		 * Helical symmetry argument is the only symmetry not to act on the inc_mirror argument. The current policy
		 * is the orientation generator using this symmetry should make its own accomodation for the inclusion of
		 * mirror orientations if the symmetry is helical (hence the presence of the is_h_sym function in
		 * the Symmetry3D class).
		 * @param inc_mirror this parameter is not specifically acted upon in this class
		 * @return a dictionary containing the keys "alt_max" and "az_max" and "alt_min"
		 * @exception InvalidValueException if nsym is less than or equal to zero
			 */
		virtual Dict get_delimiters(const bool inc_mirror=false) const;


			/** Provides access to the complete set of rotational and translational symmetry operations
		 * associated with helical symmetry. This symmetry operations are generated in a straightforward
		 * way from the parameters of this class, specifically the return Transform object has an
		 * azimuth of n times the "d_az" (as specified in the parameters of this class), and has a post
		 * translation of "tz" in the z direction.
		 * @param n the helical symmetry operation number.
		 * @return a transform containing the correct rotational and translational symmetry operation.
			 */
		virtual Transform get_sym(const int n) const;

			/** For symmetries in general this function is supposed to return the number
		 * of unique symmetric operations that can be applied for the given Symmetry3D object.
		 * For helical symmetries this is provided by the user as a parameter when setting up the helical symmetry. Generally a multiple of nstart.
		 * @return the number of symmetric rotations that can be applied without going beyond 360 degrees
		 * @exception InvalidValueException if d_az (as stored internally in parms) is less than or equal to zero
			 */
		virtual int get_nsym() const { return (int)params["nsym"]; }; 
		/*virtual int get_nsym() const {
			float daz = params.set_default("daz",0.0f);
			if ( daz <= 0 ) throw InvalidValueException(daz,"Error, you must specify a positive non zero d_az");
			return static_cast<int>(360.0/daz);
		};*/

			/** Gets the maximum cylcic symmetry exhibited by this object. This is used by OrientationGenerators, and is
		 * probably not something a general user would utilize.
		 * @return nsym - this is the symmetry of the helix
			 */
		virtual int get_max_csym() const { return (int)params["nstart"]; }	// may not be 

			/// The name of this class - used to access it from factories etc. Should be "h"
		static const string NAME;

			/** Determines whether or not this Symmetry3D is the helical type - returns true
		 * @return true - indicating that this is a helical symmetry object
			 */
		virtual bool is_h_sym() const { return true; }

			/** A function to be used when generating orientations over portion of the unit sphere
		 * defined by parameters returned by get_delimiters. In platonic symmetry altitude and azimuth
		 * alone are not enough to correctly demarcate the asymmetric unit. See the get_delimiters comments.
		 * @param altitude the EMAN style altitude of the 3D orientation in degrees
		 * @param azimuth the EMAN style azimuth of the 3D orientation in degrees
		 * @param inc_mirror whether or not to include orientations if they are in the mirror portion of the asymmetric unit
		 * @return true or false, depending on whether or not the orientation is within the asymmetric unit
			 */
		virtual bool is_in_asym_unit(const float& altitude, const float& azimuth, const bool inc_mirror) const;

			/** @param inc_mirror whether or not to include the mirror portion of the asymmetric unit
		 * @return a cyclic set of points which can be connected using great arcs on the unit sphere
		 * to demarcate the asymmetric unit. The last should may be connected to the first.
			 */
		virtual vector<Vec3f> get_asym_unit_points(bool inc_mirror = false) const;

			/** Get triangles that precisely occlude the projection area of the default asymmetric unit. This is used
		 * for collision detection in Symmetry3D::reduce
		 * @param inc_mirror whether to include the mirror portion of the asymmetric unit
			 */
		virtual vector<vector<Vec3f> > get_asym_unit_triangles(bool inc_mirror) const;
	private:
		/** Disallow copy construction */
		HSym(const HSym&);
		/** Disallow assignment */
		HSym& operator=(const HSym&);
};

	/** A base (or parent) class for the Platonic symmetries. It cannot be instantieted on its own.
	 * Doctor Phil says:
 * "see www.math.utah.edu/~alfeld/math/polyhedra/polyhedra.html for pictures of platonic solids"
 * Also, see http://blake.bcm.edu/emanwiki/EMAN2/Symmetry for a good pictorial description of what's going on here
 * This class has a fundamental role to play in terms of the Platonic symmetries that derive from it.
 * It is based heavily on the manuscript  Baldwin and Penczek, 2007. The Transform Class in SPARX and
 * EMAN2. JSB 157(250-261), where the important angles of the asymmetric units in Platonic solids are
 * described. The MOST IMPORTANT THING TO NOTE is anything that derives from this class must call init()
 * in its constructor. However, because it is unlikey that any class will inherit from this one seeing
 * as the set of Platonic symmetries is finite.
 * @author David Woolford (based on previous work by Phil Baldwin and Steve Ludtke)
 * @date Feb 2008
	 */
	class PlatonicSym : public Symmetry3D
{
	public:
		PlatonicSym() {};
		virtual  ~PlatonicSym() {};

		/** Get a dictionary containing the permissable parameters of this class
		 * Platonic symmetries actually have no parameters.
		 * @return a dictionary containing the permissable parameters of this class  ( which is none)
		 */
		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			return d;
		}

		/** Returns the range of altitude and azimuth angles which encompass
		 * the asymmetric unit of the Platonic symmetry (and more). As a general
		 * rule you may generate your orientations evenly over the range altitude
		 * range as accessed by "alt_max" key in the return dictionary, and over
		 * the azimuth range as accessed by the "az_max", but your must call
		 * the function is_in_asym_unit as you do it, to accomodate for orientations
		 * in the range that are actually beyond the asymmetric unit. See
		 * http://blake.bcm.edu/emanwiki/EMAN2/Symmetry for pictures and descriptions.
		 * If the inc_mirror is true, the return "az_max" key is twice as large as if not,
		 * but only if the platonic symmetry is Icos or Oct. If the symmetry is Tet, the
		 * mirror considerations are taken into account in is_in_asym_unit. This is a bit
		 * of a design flaw, but it works.
		 * @param inc_mirror whether or not to consider the mirror portion of the asymmetric unit (only changes the return values if the symmetry is Icos or Oct)
		 * @return a dictionary containing the "az_max" and "alt_max" keys which define angles, in degrees
		 */
		virtual Dict get_delimiters(const bool inc_mirror=false) const;

		/** A function to be used when generating orientations over portion of the unit sphere
		 * defined by parameters returned by get_delimiters. altitude and azimuth alone are not
		 * enough to correctly demarcate the asymmetric unit. See the get_delimiters comments.
		 * @param altitude the EMAN style altitude of the 3D orientation in degrees
		 * @param azimuth the EMAN style azimuth of the 3D orientation in degrees
		 * @param inc_mirror whether or not to include orientations if they are in the mirror portion of the asymmetric unit
		 * @return true or false, depending on whether or not the orientation is within the asymmetric unit
		 */
		virtual bool is_in_asym_unit(const float& altitude, const float& azimuth, const bool inc_mirror) const;

		/** Determines whether or not this Symmetry3D is the platonic type - returns true
		 * @return true - indicating that this is a platonic symmetry object
		 */
		virtual bool is_platonic_sym() const { return true; }

	protected:
		/// A dictionary that stores important angles, in radians
		Dict platonic_params;

		/** Init - Called to initialize platonic_params, should be called in the constructor of all
		 * Platonic solids that derive from this. This function generates the important angles
		 * of the platonic symmetries which is dependent only on the function get_max_csym ( which
		 * must be defined in all classes that inherit from this class)
		 */
		void init();

		/** Returns the lower bound of the asymmetric unit, as dependent on azimuth, and on alpha -
		 * alpha is alt_max for icos and oct, but may be alt_max/2.0 for tet depending on mirror
		 * symmetry etc
		 * @param azimuth an EMAN style 3D azimuth angle, in radians
		 * @param alpha an EMAN style altitude angle that helps to define arcs on the unit sphere. See Baldwin and Penczek, 2007. The Transform Class in SPARX and EMAN2. JSB 157(250-261) where the angle alpha is described
		 * @return the altitude corresponding to the lower bound for the given azimuth, in radians
		 */
		float platonic_alt_lower_bound(const float& azimuth, const float& alpha) const;

		/** @param inc_mirror whether or not to include the mirror portion of the asymmetric unit
		 * @return a cyclic set of points which can be connected using great arcs on the unit sphere
		 * to demarcate the asymmetric unit. The last should may be connected to the first.
		 */
		virtual vector<Vec3f> get_asym_unit_points(bool inc_mirror = false) const;

		/** Get triangles that precisely occlude the projection area of the default asymmetric unit. This is used
		 * for collision detection in Symmetry3D::reduce
		 * @param inc_mirror whether to include the mirror portion of the asymmetric unit
		 */
		virtual vector<vector<Vec3f> > get_asym_unit_triangles(bool inc_mirror) const;
	private:
		/** Disallow copy construction */
		PlatonicSym(const PlatonicSym&);
		/** Disallow assignment */
		PlatonicSym& operator=(const PlatonicSym&);
};

	/** An encapsulation of tetrahedral symmetry
	 * Doctor Phil has this to say about tetrahedral symmetry:
 * " Each Platonic Solid has 2E symmetry elements.
 *	 The tetrahedron has n=m=3; F=4, E=6=nF/2, V=4=nF/m.
 *   It is composed of four triangles."
 * @author David Woolford (based on previous work by Phil Baldwin and Steve Ludtke)
 * @date Feb 2008
	 */
	class TetrahedralSym : public PlatonicSym
{
	public:
		/** Constructor calls PlatonicSym::init
		 */
		TetrahedralSym()  {init();}
		virtual  ~TetrahedralSym() {}

		/** Factory support function NEW
		 * @return a newly instantiated class of this type
		 */
		static Symmetry3D *NEW()
		{
			return new TetrahedralSym();
		}

		/** Return TetrahedralSym::NAME
		 * @return the unique name of this class
		 */
		virtual string get_name() const { return NAME; }


		/** Get a description
		 * @return a clear desciption of this class
		 */
		virtual string get_desc() const { return "Tetrahedral symmetry support"; }

		/** Gets the maximum symmetry of this object. This is used by OrientationGenerators, and is
		 * probably not something a general user would utilize.
		 * @return for tetrahedral symmetry, this number is 3
		 */
		virtual int get_max_csym() const { return 3; }

		/** This function provides access to the unique rotational symmetries of a tetrahedron.
		 * In this implementation, the tetrahedral symmetry group has a face along the z-axis. In all, there are
		 * 12 (accessed by get_nysm) unique rotational symmetric operations for the tetrahedron.
		 * In the terminology defined Append A (titled Symmetry Elements) in the manuscript  Baldwin and Penczek, 2007.
		  * The Transform Class in SPARX and EMAN2. JSB 157(250-261), Doctor Phil has this to say:
		 * "B^3=A^3=1;  BABA=1; implies   A^2=BAB, ABA=B^2 , AB^2A = B^2AB^2 and
		 *  12 words with at most a single A
		 *  1 B BB  A  BA AB BBA BAB ABB BBAB BABB BBABB
		 *  at most one A is necessary"
		 * @param n the symmetric operation number
		 * @return a transform containing the correct rotational symmetry operation.
		 */
		virtual Transform get_sym(const int n) const;

		/** In tetrahedral symmetry special consideration must be taken when generating orientations
		 * in the asymmetric unit. This function is a specialization of the functionality in
		 * PlatonicSym::is_in_asym_unit
		 * @param altitude the EMAN style altitude of the 3D orientation in degrees
		 * @param azimuth the EMAN style azimuth of the 3D orientation in degrees
		 * @param inc_mirror whether or not to include orientations if they are in the mirror portion of the asymmetric unit
		 * @return true or false, depending on whether or not the orientation is within the asymmetric unit
		 */
		virtual bool is_in_asym_unit(const float& altitude, const float& azimuth, const bool inc_mirror) const;

		/** Gets the total number of unique roational symmetry operations associated with this symmetry
		 * For tetrahedral symmetry symmetry, this is 12
		 * @return 12
		 */
		virtual int get_nsym() const { return 12; };

		/** Get the azimuth alignment offset required to ensure that orientations align correctly
		 * with symmetric axes of the tetrahedron. This offset is directly related to the way
		 * the symmetric operations are generated by get_sym. All orientations generated as a
		 * result of using the delimiters supplied by this class should by offset by this azimuth
		 * to ensure proper alignment with tetrahedral objects in EMAN2
		 */
		virtual float get_az_alignment_offset() const;

		/// The name of this class - used to access it from factories etc. Should be "tet"
		static const string NAME;

		/** @param inc_mirror whether or not to include the mirror portion of the asymmetric unit
		 * @return a cyclic set of points which can be connected using great arcs on the unit sphere
		 * to demarcate the asymmetric unit. The last should may be connected to the first.
		 */
		virtual vector<Vec3f> get_asym_unit_points(bool inc_mirror = false) const;

		/** A function that is used to determine if this is the tetrahedral symmetry object
		 * @return true - indicating that this is not a tetrahedral symmetry object
		 */
		virtual bool is_tet_sym() const { return true; }

	private:
		/** Disallow copy construction */
		TetrahedralSym(const TetrahedralSym&);
		/** Disallow assignment */
		TetrahedralSym& operator=(const TetrahedralSym&);


};

	/** An encapsulation of octahedral symmetry
	* Doctor Phil has this to say about the octahedral symmetry:
 * "Each Platonic Solid has 2E symmetry elements.
 * "A cube has   m=3, n=4, F=6 E=12=nF/2, V=8=nF/m,since vertices shared by 3 squares;
 *  It is composed of 6 squares.
 *  An octahedron has   m=4, n=3, F=8 E=12=nF/2, V=6=nF/m,since vertices shared by 4 triangles"
 * @author David Woolford (based on previous work by Phil Baldwin and Steve Ludtke)
 * @date Feb 2008
	 */

	class OctahedralSym : public PlatonicSym
{
	public:
		/** Constructor calls PlatonicSym::init
		 */
		OctahedralSym()  {init();}
		virtual  ~OctahedralSym() {}

		/** Factory support function NEW
		 * @return a newly instantiated class of this type
		 */
		static Symmetry3D *NEW()
		{
			return new OctahedralSym();
		}

		/** Return OctahedralSym::NAME
		 * @return the unique name of this class
		 */
		virtual string get_name() const { return NAME; };

		/** Get a description
		 * @return a clear desciption of this class
		 */
		virtual string get_desc() const { return "Octahedral symmetry support"; }

		/** Gets the maximum symmetry of this object. This is used by OrientationGenerators, and is
		 * probably not something a general user would utilize.
		 * @return for octahedral symmetry, this number is 4
		 */
		virtual int get_max_csym() const { return 4; }

		/** This function provides access to the unique rotational symmetries of an octahedron.
		 * We have placed the octahedral symmetry group with a face along the z-axis. In all, there are
		 * 24 (accessed by get_nysm) unique rotational symmetric operations for the octahedron.
		 * In the terminology defined Append A (titled Symmetry Elements) in the manuscript  Baldwin and Penczek, 2007.
		 * The Transform Class in SPARX and EMAN2. JSB 157(250-261), Doctor Phil has this to say:
		 * "B^4=A^3=1;  BABA=1; implies   AA=BAB, ABA=B^3 , AB^2A = BBBABBB and
		 *  20 words with at most a single A
		 *  1 B BB BBB A  BA AB BBA BAB ABB BBBA BBAB BABB ABBB BBBAB BBABB BABBB
		 *  BBBABB BBABBB BBBABBB
		 *  also ABBBA is distinct yields 4 more words
		 *  ABBBA BABBBA BBABBBA BBBABBBA
		 *  for a total of 24 words
		 *  Note A BBB A BBB A  reduces to BBABB
		 *  and B A BBB A is the same as A BBB A BBB etc."
		 * @param n the symmetric operation number.
		 * @return a transform containing the correct rotational symmetry operation.
		 */
		virtual Transform get_sym(const int n) const;

		/** Gets the total number of unique roational symmetry operations associated with this symmetry
		 * For octahedral symmetry this is 24
		 * @return 24
		 */
		virtual int get_nsym() const { return 24; };

		/// The name of this class - used to access it from factories etc. Should be "oct"
		static const string NAME;
	private:
		/** Disallow copy construction */
		OctahedralSym(const OctahedralSym&);
		/** Disallow assignment */
		OctahedralSym& operator=(const OctahedralSym&);
};

	/** An encapsulation of icosahedral symmetry
	* Doctor Phil has this to say about icosahedral symmetry:
 * "Each Platonic Solid has 2E symmetry elements.
 * An icosahedron has   m=5, n=3, F=20 E=30=nF/2, V=12=nF/m,since vertices shared by 5 triangles
 * It is composed of 20 triangles. E=3*20/2
 * A  dodecahedron has m=3, n=5   F=12 E=30  V=20
 * It is composed of 12 pentagons. E=5*12/2;   V= 5*12/3, since vertices shared by 3 pentagons"
 * @author David Woolford (based on previous work by Phil Baldwin and Steve Ludtke)
 * @date Feb 2008
	 */
	class IcosahedralSym : public PlatonicSym
{
	public:
		/** Constructor calls PlatonicSym::init
		 */
		IcosahedralSym() {init(); }
		virtual  ~IcosahedralSym() { }

		/** Factory support function NEW
		 * @return a newly instantiated class of this type
		 */
		static Symmetry3D *NEW()
		{
			return new IcosahedralSym();
		}

		/** Return IcosahedralSym::NAME
		 * @return the unique name of this class
		 */
		virtual string get_name() const { return NAME; };

		/** Get a description
		 * @return a clear desciption of this class
		 */
		virtual string get_desc() const { return "Icosahedral symmetry support"; }

		/** Gets the maximum symmetry of this object. This is used by OrientationGenerators, and is
		 * probably not something a general user would utilize.
		 * @return for icosahedral symmetry, this number is 5
		 */
		virtual int get_max_csym() const { return 5; }// 5 is the greatest symmetry

		/** This function provides access to the unique rotational symmetries of an icosahedron.
		 * We have placed the icosahedral symmetry group with a face along the z-axis. In all, there are
		 * 60 (accessed by get_nysm) unique rotational symmetric operations for the icosahedron.
		 * @param n the symmetric operation number.
		 * @return a transform containing the correct rotational symmetry operation.
		 */
		virtual Transform get_sym(const int n) const;

		/** Gets the total number of unique rotational symmetry operations associated with this symmetry
		 * For icosahedral symmetry, this is 60
		 * @return 60
		 */
		virtual int get_nsym() const { return 60; };

		/** Get the azimuth alignment offset required to ensure that orientations align correctly
		 * with symmetric axes of the icosahedron. This offset is directly related to the way
		 * the symmetric operations are generated by get_sym. All orientations generated as a
		 * result of using the delimiters supplied by this class should by offset by this azimuth
		 * to ensure proper alignment with tetrahedral objects in EMAN2
		 */
		virtual float get_az_alignment_offset() const;

		/// The name of this class - used to access it from factories etc. Should be "icos"
		static const string NAME;
	private:
		/** Disallow copy construction */
		IcosahedralSym(const IcosahedralSym&);
		/** Disallow assignment */
		IcosahedralSym& operator=(const IcosahedralSym&);
};

/** An encapsulation of icosahedral symmetry 222 */
class Icosahedral2Sym : public PlatonicSym
{
	public:
		/** Constructor calls PlatonicSym::init
		 */
		Icosahedral2Sym() { }
		virtual  ~Icosahedral2Sym() { init(); }

		/** Factory support function NEW
		 * @return a newly instantiated class of this type
		 */
		static Symmetry3D *NEW()
		{
			return new Icosahedral2Sym();
		}

		/** Return IcosahedralSym::NAME
		 * @return the unique name of this class
		 */
		virtual string get_name() const { return NAME; };

		/** Get a description
		 * @return a clear desciption of this class
		 */
		virtual string get_desc() const { return "Icosahedral 2 symmetry support"; }

		/** Gets the maximum symmetry of this object. This is used by OrientationGenerators, and is
		 * probably not something a general user would utilize.
		 * @return for icosahedral symmetry, this number is 2
		 */
		virtual int get_max_csym() const { return 2; }

		/** This function provides access to the unique rotational symmetries of an icosahedron.
		 * We have placed the icosahedral symmetry group with a face along the z-axis. In all, there are
		 * 60 (accessed by get_nysm) unique rotational symmetric operations for the icosahedron.
		 * @param n the symmetric operation number.
		 * @return a transform containing the correct rotational symmetry operation.
		 */
		virtual Transform get_sym(const int n) const;

		/** Gets the total number of unique roational symmetry operations associated with this symmetry
		 * For icosahedral symmetry, this is 60
		 * @return 60
		 */
		virtual int get_nsym() const { return 60; };

		/** Get the azimuth alignment offset required to ensure that orientations align correctly
		 * with symmetric axes of the icosahedron. This offset is directly related to the way
		 * the symmetric operations are generated by get_sym. All orientations generated as a
		 * result of using the delimiters supplied by this class should by offset by this azimuth
		 * to ensure proper alignment with tetrahedral objects in EMAN2
		 */
		virtual float get_az_alignment_offset() const;

		/// The name of this class - used to access it from factories etc. Should be "icos2"
		static const string NAME;
	private:
		/** Disallow copy construction */
		Icosahedral2Sym(const Icosahedral2Sym&);
		/** Disallow assignment */
		Icosahedral2Sym& operator=(const Icosahedral2Sym&);
};

	/// A template specialization for the Symmetry3D factory. Adds all of the symmetries
	template <> Factory < Symmetry3D >::Factory();
	/// A template specialization for get - so people can call get with strings like "c1","d4" etc - this avoids have to use Dicts to specify the nsym
	template <> Symmetry3D* Factory < Symmetry3D >::get(const string & instancename);
	/// dump symmetries, useful for obtaining symmetry information
	void dump_symmetries();
	/// dump_symmetries_list, useful for obtaining symmetry information
	map<string, vector<string> > dump_symmetries_list();

	/** An orientation generator is a kind of class that will generate orientations for a given symmetry
	 * If one needs to generate orientations in the unit sphere, one simply uses the C1 symmetry.
	 * It inherits from a factory base, making it amenable to incorporation in EMAN2 style factories.
	 * Objects that inherit from this class must write a gen_orientations function, in addition to
	 * fulfilling the responsibilities of the FactoryBase class
	 * @author David Woolford
	 * @date Feb 2008
	 */
	class OrientationGenerator : public FactoryBase
{
	public:
		OrientationGenerator() {};
		virtual ~OrientationGenerator() {};

		/** generate orientations given some symmetry type
		 * @param sym the symmetry which defines the interesting asymmetric unit
		 * @return a vector of Transform objects containing the generated set of orientations
		 */
		virtual vector<Transform> gen_orientations(const Symmetry3D* const sym) const  = 0;

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("phitoo", EMObject::FLOAT,  "Specifying a non zero value for this argument will cause phi rotations to be included. The value specified is the angular spacing of the phi rotations in degrees. The default for this value is 0, causing no extra phi rotations to be included.");
			d.put("random_phi", EMObject::BOOL,  "Causes the orientations to have a random phi. This occurs before the phitoo parameter is considered.");
			return d;
		}

		/** This functions adds one or more Transform objects to the vector v, depending
		 * on the parameters stored in the dictionary (which the inheriting class may
		 * include by initializing the typedict in get_param_types by calling
		 *
		 * @code
		 * TypeDict d = OrientationGenerator::get_param_types();
		 * @endcode
		 *
		 * to initialize. If phitoo is no zero, this cause extra orientations to be included in phi (in steps of phitoo).
		 * If random_phi is true, the phi of the Transform object is randomized.
		 * This function is for internal convenience of child classes.
		 * @param v the vector to add Transform objects to
		 * @param az the azimuth to be used as a basis for generated Transform objects (in degrees)
		 * @param alt the altitude to be used as a basis for generated Transform objects (in degrees)
		 * @return and indication of success (true or false). False is only ever return if phitoo is less than 0.
		 */
		bool add_orientation(vector<Transform>& v, const float& az, const float& alt) const;



		/** This function gets the optimal value of the delta (or angular spacing) of the orientations
		 * based on a desired total number of orientations (n). It does this using a bifurcation strategy,
		 * calling get_orientations_tally (which must be supplied by the child class)
		 * using the next best guess etc. The solution may not exist (simply  because the
		 * orientation generation strategy does not contain it), so a best guess may be returned.
		 *
		 * The inheriting class must supply the get_orientations_tally function, which returns
		 * the number of orientations generated (an int), for a given delta.
		 *
		 * @param sym the symmetry which defines the interesting asymmetric unit
		 * @param n the desired number of orientations
		 * @return the optimal value of delta to ensure as near to the desired number of orientations is generated
		 */
		float get_optimal_delta(const Symmetry3D* const sym, const int& n) const;

		/** This function returns how many orientations will be generated for a given delta (angular spacing)
		 * It should general do this by simulating the function gen_orientations
		 * @param sym the symmetry which defines the interesting asymmetric unit
		 * @param delta the desired angular spacing of the orientations
		 * @return the number of orientations that will be generated using these parameters
		 */
		virtual int get_orientations_tally(const Symmetry3D* const sym, const float& delta) const = 0;

	protected:
		void get_az_max(const Symmetry3D* const sym, const float& altmax, const bool inc_mirror, const float& alt_iterator,const float& h,bool& d_odd_mirror_flag, float& azmax_adjusted) const;

	private:
		/** Disallow copy construction */
		OrientationGenerator(const OrientationGenerator&);
		/** Disallow assignment */
		OrientationGenerator& operator=(const OrientationGenerator&);
};



 /** EmanOrientationGenerator generates orientations quasi-evenly distributed in the asymmetric unit.
 * Historically, it is an adaptation of the method first used in EMAN1 and developed by Steve
 * Ludtke. In EMAN2 it is more or less the same thing, but with more precise treatmeant of the
 * platonic symmetries. In terms of approach, the altitude angles in the asymmetric unit are traversed
 * constantly in steps of "prop" (a parameter of this class). However, the azimuth steps vary according
 * to altitude, and this helps to achieve a more even distribution of orientations.
 * @author David Woolford (based on previous work by Phil Baldwin and Steve Ludtke)
 * @date Feb 2008
  */
class EmanOrientationGenerator : public OrientationGenerator
{
	public:
		EmanOrientationGenerator() {};
		virtual  ~EmanOrientationGenerator() {};

		/** Factory support function NEW
		 * @return a newly instantiated class of this type
		 */
		static OrientationGenerator *NEW()
		{
			return new EmanOrientationGenerator();
		}

		/** Return 	"eman"
		 * @return the unique name of this class
		 */
		virtual string get_name() const { return NAME; }

		/** Get a description
		 * @return a clear desciption of this class
		 */
		virtual string get_desc() const { return "Generate orientations distributed quasi-uniformaly over the asymmetric unit using an altitude-proportional strategy"; }

		/** Get a dictionary containing the permissable parameters of this class
		 * @return a dictionary containing the permissable parameters of this class
		 * parameters are explained in the dictionary itself
		 */
		virtual TypeDict get_param_types() const
		{
			TypeDict d = OrientationGenerator::get_param_types();
			d.put("delta", EMObject::FLOAT, "The angular separation of orientations in degrees. This option is mutually exclusively of the n argument.");
			d.put("perturb", EMObject::BOOL, "Whether or not to perturb the generated orientations in a small local area, default is true.");
			d.put("n", EMObject::INT, "The number of orientations to generate. This option is mutually exclusively of the delta argument.Will attempt to get as close to the number specified as possible.");
			d.put("inc_mirror", EMObject::BOOL, "Indicates whether or not to include the mirror portion of the asymmetric unit. Default is false.");
			d.put("alt_min", EMObject::FLOAT, "Minimum altitude value to include (alt=0 is Z axis). Default=0");
			d.put("alt_max", EMObject::FLOAT, "Maximum altitude value to include (alt=90 is X-Y plane). Default=no limit");
			d.put("breaksym", EMObject::BOOL, "If specified, still generates orientations filling the unit (hemi)sphere, but does it by filling one asymmetric unit, then generating all symmetric equivalents.");
			return d;
		}

		/** generate orientations given some symmetry type
		 * @param sym the symmetry which defines the interesting asymmetric unit
		 * @return a vector of Transform objects containing the set of evenly distributed orientations
		 */
		virtual vector<Transform> gen_orientations(const Symmetry3D* const sym) const;

		/// The name of this class - used to access it from factories etc. Should be "icos"
		static const string NAME;
	private:
		/** Disallow copy construction */
		EmanOrientationGenerator(const EmanOrientationGenerator&);
		/** Disallow assignment */
		EmanOrientationGenerator& operator=(const EmanOrientationGenerator&);

		/** This function returns how many orientations will be generated for a given delta (angular spacing)
		 * It does this by simulated gen_orientations.
		 * @param sym the symmetry which defines the interesting asymmetric unit
		 * @param delta the desired angular spacing of the orientations
		 * @return the number of orientations that will be generated using these parameters
		 */
		virtual int get_orientations_tally(const Symmetry3D* const sym, const float& delta) const;


		/** Gets the optimum azimuth delta (angular step) for a given altitude, delta and
		 * maximum symmetry. This function is important for the generation of evenly distributed
		 * orientations
		 * @param delta - the angular spacing of the altitude angles, this is usually the "delta" parameter
		 * @param altitude the altitude along which the azimuth is going to be varied
		 * @param maxcsym the maximum csym of the Symmetry3D object - this is usually Symmetry3D::get_max_csym
		 * @return the optimal azimuth angular spacing
		 */
		float get_az_delta(const float& delta,const float& altitude, const int maxcsym) const;

};

/** Random Orientation Generator - carefully generates uniformly random orientations in any asymmetric unit.
 *  For points distributed in the unit sphere, just use the CSym type with nysm = 1.
 * (i.e. c1 symmetry)
 * @author David Woolford
 * @date March 2008
 */
class RandomOrientationGenerator : public OrientationGenerator
{
	public:
		RandomOrientationGenerator() {}
		virtual ~RandomOrientationGenerator() {}

		/** Factory support function NEW
		 * @return a newly instantiated class of this type
		 */
		static OrientationGenerator *NEW()
		{
			return new RandomOrientationGenerator();
		}

		/** Return 	"random"
		 * @return the unique name of this class
		 */
		virtual string get_name() const { return NAME; }

		/** Get a description
		 * @return a clear desciption of this class
		 */
		virtual string get_desc() const { return "Generate random orientations within an asymmetric unit"; }

		/** Get a dictionary containing the permissable parameters of this class
		 * @return a dictionary containing the permissable parameters of this class
		 * parameters are explained in the dictionary itself
		 */
		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("n", EMObject::INT, "The number of orientations to generate.");
			d.put("inc_mirror", EMObject::BOOL, "Indicates whether or not to include the mirror portion of the asymmetric unit. Default is false.");
			d.put("phitoo", EMObject::BOOL, "Makes phi random as well");
			return d;
		}

		/** Generate random orientations in the asymmetric unit of the symmetry
		 * @param sym the symmetry which defines the interesting asymmetric unit
		 * @return a vector of Transform objects containing the set of evenly distributed orientations
		 */
		virtual vector<Transform> gen_orientations(const Symmetry3D* const sym) const;

		/// The name of this class - used to access it from factories etc.
		static const string NAME;

		virtual int get_orientations_tally(const Symmetry3D* const sym, const float& delta) const { (void)sym; (void)delta; return 0; }
	private:
		/** Disallow copy construction */
		RandomOrientationGenerator(const RandomOrientationGenerator&);
		/** Disallow assignment */
		RandomOrientationGenerator& operator=(const RandomOrientationGenerator&);
};

/**Sparx even orientation generator - see util_sparx.cpp - Util::even_angles(...)
 * This orientation generator is based on work presented in Penczek et al., 1994 P.A. Penczek, R.A.
 * Grassucci and J. Frank, The ribosome at improved resolution: new techniques for merging and
 * orientation refinement in 3D cryo-electron microscopy of biological particles, Ultramicroscopy 53 (1994).
 *
 * This is a proportional approach very similar to the eman approach - the differences between these
 * two approaches is mostly visible near altitude=0
 *
 * @author David Woolford (ported directly from Sparx utilities.py, which is written by Pawel Penczek)
 * @date March 2008
 */
	class EvenOrientationGenerator : public OrientationGenerator
{
	public:
		EvenOrientationGenerator() {}
		virtual ~EvenOrientationGenerator() {}

		/** Factory support function NEW
		 * @return a newly instantiated class of this type
		 */
		static OrientationGenerator *NEW()
		{
			return new EvenOrientationGenerator();
		}


		/** Return 	"even"
		 * @return the unique name of this class
		 */
		virtual string get_name() const { return NAME; }

		/** Get a description
		 * @return a clear desciption of this class
		 */
		virtual string get_desc() const { return "Generate quasi-evenly distributed orientations within an asymmetric unit using Penczek's (94) approach"; }

		/** Get a dictionary containing the permissable parameters of this class
		 * @return a dictionary containing the permissable parameters of this class
		 * parameters are explained in the dictionary itself
		 */
		virtual TypeDict get_param_types() const
		{
			TypeDict d = OrientationGenerator::get_param_types();
			d.put("delta", EMObject::FLOAT, "The angular separation of orientations in degrees. This option is mutually exclusively of the n argument.");
			d.put("inc_mirror", EMObject::BOOL, "Indicates whether or not to include the mirror portion of the asymmetric unit. Default is false.");
			d.put("n", EMObject::INT, "The number of orientations to generate. This option is mutually exclusively of the delta argument.Will attempt to get as close to the number specified as possible.");
			return d;
		}

		/** Generate even distributed orientations in the asymmetric unit of the symmetry
		 * @param sym the symmetry which defines the interesting asymmetric unit
		 * @return a vector of Transform objects containing the set of evenly distributed orientations
		 */
		virtual vector<Transform> gen_orientations(const Symmetry3D* const sym) const;

		/// The name of this class - used to access it from factories etc. Should be "icos"
		static const string NAME;
	private:
		/** Disallow copy construction */
		EvenOrientationGenerator(const EvenOrientationGenerator&);
		/** Disallow assignment */
		EvenOrientationGenerator& operator=(const EvenOrientationGenerator&);
		/** This function returns how many orientations will be generated for a given delta (angular spacing)
		 * It does this by simulated gen_orientations.
		 * @param sym the symmetry which defines the interesting asymmetric unit
		 * @param delta the desired angular spacing of the orientations
		 * @return the number of orientations that will be generated using these parameters
		 */
		virtual int get_orientations_tally(const Symmetry3D* const sym, const float& delta) const;

};

/** Saff orientation generator - based on the work of Saff and Kuijlaars, 1997 E.B. Saff and A.B.J. Kuijlaars, Distributing many points on a sphere,
 * Mathematical Intelligencer 19 (1997), pp. 5–11. This is a spiral based approach
 * @author David Woolford (ported directly from Sparx utilities.py, which is written by Pawel Penczek)
 * @date March 2008
 */
class SaffOrientationGenerator : public OrientationGenerator
{
	public:
		SaffOrientationGenerator() {}
		virtual ~SaffOrientationGenerator() {}

		/** Factory support function NEW
		 * @return a newly instantiated class of this type
		 */
		static OrientationGenerator *NEW()
		{
			return new SaffOrientationGenerator();
		}

		/** Return 	"saff"
		 * @return the unique name of this class
		 */
		virtual string get_name() const { return NAME; }

		/** Get a description
		 * @return a clear desciption of this class
		 */
		virtual string get_desc() const { return "Generate quasi-evenly distributed orientations within an asymmetric unit using a spiraling method attributed to Saff"; }

		/** Get a dictionary containing the permissable parameters of this class
		 * @return a dictionary containing the permissable parameters of this class
		 * parameters are explained in the dictionary itself
		 */
		virtual TypeDict get_param_types() const
		{
			TypeDict d = OrientationGenerator::get_param_types();
			d.put("n", EMObject::INT, "The number of orientations to generate. This option is mutually exclusively of the delta argument.Will attempt to get as close to the number specified as possible.");
			d.put("inc_mirror", EMObject::BOOL, "Indicates whether or not to include the mirror portion of the asymmetric unit. Default is false.");
			d.put("delta", EMObject::FLOAT, "The angular separation of orientations in degrees. This option is mutually exclusively of the n argument.");
			return d;
		}

		/** Generate Saff orientations in the asymmetric unit of the symmetry
		 * @param sym the symmetry which defines the interesting asymmetric unit
		 * @return a vector of Transform objects containing the set of evenly distributed orientations
		 */
		virtual vector<Transform> gen_orientations(const Symmetry3D* const sym) const;

		/// The name of this class - used to access it from factories etc. Should be "icos"
		static const string NAME;
	private:
		/** Disallow copy construction */
		SaffOrientationGenerator(const SaffOrientationGenerator&);
		/** Disallow assignment */
		SaffOrientationGenerator& operator=(const SaffOrientationGenerator&);
		/** This function returns how many orientations will be generated for a given delta (angular spacing)
		 * It does this by simulated gen_orientations.
		 * @param sym the symmetry which defines the interesting asymmetric unit
		 * @param delta the desired angular spacing of the orientations
		 * @return the number of orientations that will be generated using these parameters
		 */
		virtual int get_orientations_tally(const Symmetry3D* const sym, const float& delta) const;

		// This was a function that paid special considerations to the overall algorithm in the
		// case of the Platonic symmetries, which have non trivial asymmetric units. But unfortunately
		// it was bug-prone, and the approach in place already seemed good enough
		//vector<Transform> gen_platonic_orientations(const Symmetry3D* const sym, const float& delta) const;
};


/** Optimum orientation generator. Optimally distributes points on a unit sphere, then slices out
 * a correctly sized asymmetric unit, depending on the symmetry type. The approach relies on an initial
 * distribution of points on the unit sphere, which may be generated using any of the other orientation
 * generators. By default, the Saff orientation generator is used.
 *
 * @author David Woolford
 * @date March 2008
 */
class OptimumOrientationGenerator : public OrientationGenerator
{
	public:
		OptimumOrientationGenerator() {}
		virtual ~OptimumOrientationGenerator() {}

		/** Factory support function NEW
		 * @return a newly instantiated class of this type
		 */
		static OrientationGenerator *NEW()
		{
			return new OptimumOrientationGenerator();
		}

		/** Return 	"opt"
		 * @return the unique name of this class
		 */
		virtual string get_name() const { return NAME; }

		/** Get a description
		 * @return a clear desciption of this class
		*/
		virtual string get_desc() const { return "Generate optimally distributed orientations within an asymmetric using a basic optimization technique"; }

		/** Get a dictionary containing the permissable parameters of this class
		 * @return a dictionary containing the permissable parameters of this class
		 * parameters are explained in the dictionary itself
		 */
		virtual TypeDict get_param_types() const
		{
			TypeDict d = OrientationGenerator::get_param_types();
			d.put("n", EMObject::INT, "The number of orientations to generate. This option is mutually exclusively of the delta argument.Will attempt to get as close to the number specified as possible.");
			d.put("inc_mirror", EMObject::BOOL, "Indicates whether or not to include the mirror portion of the asymmetric unit. Default is false.");
			d.put("delta", EMObject::FLOAT, "The angular separation of orientations in degrees. This option is mutually exclusively of the n argument.");
			d.put("use", EMObject::STRING, "The orientation generation technique used to generate the initial distribution on the unit sphere.");
			return d;
		}

		/** Generate Saff orientations in the asymmetric unit of the symmetry
		 * @param sym the symmetry which defines the interesting asymmetric unit
		 * @return a vector of Transform objects containing the set of evenly distributed orientations
		 */
		virtual vector<Transform> gen_orientations(const Symmetry3D* const sym) const;

		/// The name of this class - used to access it from factories etc. Should be "icos"
		static const string NAME;
	private:
		/** Disallow copy construction */
		OptimumOrientationGenerator(const OptimumOrientationGenerator&);
		/** Disallow assignment */
		OptimumOrientationGenerator& operator=(const OptimumOrientationGenerator&);
		/** This function returns how many orientations will be generated for a given delta (angular spacing)
		 * It does this by simulated gen_orientations.
		 * @param sym the symmetry which defines the interesting asymmetric unit
		 * @param delta the desired angular spacing of the orientations
		 * @return the number of orientations that will be generated using these parameters
		 */
		virtual int get_orientations_tally(const Symmetry3D* const sym, const float& delta) const;


		/// Optimize the distances in separating points on the unit sphere, as described by the
		/// the rotations in Transform objects.
		vector<Vec3f> optimize_distances(const vector<Transform>& v) const;
};

	/// Template specialization for the OrientationGenerator class
	template <> Factory < OrientationGenerator >::Factory();
	/// Dumps useful information about the OrientationGenerator factory
	void dump_orientgens();
	/// Can be used to get useful information about the OrientationGenerator factory
	map<string, vector<string> > dump_orientgens_list();
} // namespace EMAN

#endif // eman__symmetry_h__

