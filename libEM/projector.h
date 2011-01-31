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

#ifndef eman__projector_h__
#define eman__projector_h__ 1

#include "transform.h"

using std::string;

namespace EMAN
{
	class EMData;

	/** Projector class defines a method to generate 2D projections
	 * from a 3D model. Projector class is the base class for all projectors.
	 * Each specific projector has a unique name and should be called
     * through the name.
     *
	 * All Projector classes in EMAN are managed by a Factory
	 * pattern. So each Projector class must define:
	 *   - a unique name to idenfity itself in the factory.
	 *   - a static method to register itself in the factory.
	 *
     * Typical usage of Projectors:
     *
     *  - How to get all the Projector types
     @code
     *    vector<string> all_projectors = Factory<Projector>::get_list();
     @endcode
	 *
     *  - How to use a Projector
     @code
     *    EMData* img = ...;
     *    Projector* proj = Factory<Projector>::get("fft");
     *    EMData* result = proj->project3d(img);
     @endcode
	 *
     *  - How to define a new Projector
     *    A new projector type "XYZProjector" should implement at
     *    least the following 3 functions:
     @code
     *        EMData *project3d(EMData * em) const;
     *        string get_name() const { return "xyz"; }
     *        static Projector* NEW() { return new XYZProjector(); }
	 @endcode
	*/

	class Projector
	{
	  public:
		virtual ~ Projector()
		{
		}

		/** Project an 3D image into a 2D image.
		 * @return A 2D image from the projection.
		 */
		virtual EMData *project3d(EMData * image) const = 0;

		/** Back-project a 2D image into a 3D image.
		 * @return A 3D image from the backprojection.
                 */
		virtual EMData *backproject3d(EMData * image) const = 0;

		/** Get the projector's name. Each projector is indentified by
		 * unique name.
		 * @return The projector's name.
		 */
		virtual string get_name() const = 0;

		virtual string get_desc() const = 0;

		/** Get the projector parameters in a key/value dictionary.
		 * return A key/value pair dictionary containing the
		 *         parameters.
		 */
		virtual Dict get_params() const
		{
			return params;
		}
		/** Set the projector parameters using a key/value dictionary */
		void set_params(const Dict & new_params)
		{
			params = new_params;
		}
		/** Get processor parameter information in a dictionary. Each
		 * parameter has one record in the dictionary. Each record
		 * contains its name, data-type, and description.
		 *
		 * @return A dictionary containing the parameter info.
		 */
		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			return d;
		}

	  protected:
		Dict params;
	};

	/** Gaussian FFT 3D projection.
     * use integer 'mode' to determine the gaussian width and the way
     * to interpolate a point in a 3d complex image.
     * valid mode range: [1,7].
     * the gauss widths are as follows with mode from 1 to 7:
     *
     * mode 1:  0;
     * mode 2:  4.0 / (M_PI * M_PI);
     * mode 3:  6.4 / (M_PI * M_PI);
     * mode 4:  8.8 / (M_PI * M_PI);
     * mode 5:  0;
     * mode 6:  10.4 / (M_PI * M_PI);
     * mode 7:  10.4 / (M_PI * M_PI);
     */

	class GaussFFTProjector:public Projector
	{
	  public:
		GaussFFTProjector():alt(0), az(0), phi(0)
		{
		}

		EMData *project3d(EMData * image) const;
                // no implementation yet
		EMData *backproject3d(EMData * image) const;


		void set_params(const Dict & new_params)
		{
			Projector::set_params(new_params);
			alt = params["alt"];
			az = params["az"];
			phi = params["phi"];
		}

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Projections using a Gaussian kernel in Fourier space. Produces artifacts, not recommended.";
		}

		static Projector *NEW()
		{
			return new GaussFFTProjector();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("transform", EMObject::TRANSFORM);
			d.put("mode", EMObject::INT);
			return d;
		}
		
		static const string NAME;

	  private:
		float alt, az, phi;
		bool interp_ft_3d(int mode, EMData * image, float x, float y,
						  float z, float *data, float gauss_width) const;
	};

	/** Fourier gridding projection routine.
	 *
	 *  @see P. A. Penczek, R. Renka, and H. Schomberg,
	 *       J. Opt. Soc. Am. A _21_, 499-509 (2004)
     */
	class FourierGriddingProjector:public Projector
	{
	  public:
		EMData * project3d(EMData * image) const;
                // no implementation yet
		EMData * backproject3d(EMData * image) const;


		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Fourier-space projection using gridding.";
		}

		static Projector *NEW()
		{
			return new FourierGriddingProjector();
		}
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("transform", EMObject::TRANSFORM);
			d.put("kb_alpha", EMObject::FLOAT);
			d.put("kb_K", EMObject::FLOAT);
			d.put("angletype", EMObject::STRING);
			d.put("anglelist", EMObject::FLOATARRAY);
			d.put("theta", EMObject::FLOAT);
			d.put("psi", EMObject::FLOAT);
			d.put("npad", EMObject::INT);
			return d;
		}
		
		static const string NAME;
	};


	/** Pawel Penczek's optimized projection routine.
     */
	class PawelProjector:public Projector
	{
	  public:
		EMData * project3d(EMData * image) const;
		EMData * backproject3d(EMData * image) const;

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Pawel Penczek's optimized real-space projection generation. Minor interpolation artifacts.";
		}

		static Projector *NEW()
		{
			return new PawelProjector();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("transform", EMObject::TRANSFORM);
			d.put("origin_x", EMObject::INT);
			d.put("origin_y", EMObject::INT);
			d.put("origin_z", EMObject::INT);
			d.put("radius", EMObject::INT);
			d.put("anglelist", EMObject::FLOATARRAY);
			d.put("angletype", EMObject::STRING);
			d.put("theta", EMObject::FLOAT);
			d.put("psi", EMObject::FLOAT);
			return d;
		}
		
		static const string NAME;

	  private:
		// Same structure as the IPCUBE structure in Spider
		struct IPCube
		{
			int start;
			int end;
			Vec3i loc;
		};
		// Process the number of valid x-lines (rows)
		// within the radius
		void prepcubes(int nx, int ny, int nz, int ri, Vec3i origin,
				       int& nn, IPCube* ipcube=NULL) const;
	};

	/** Fast real-space 3D projection.
	 * @param Transform object used for projection
     */
	class StandardProjector:public Projector
	{
	  public:
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("transform", EMObject::TRANSFORM, "Transform object used for projection");
			return d;
		}

		EMData * project3d(EMData * image) const;
                // no implementation yet
		EMData * backproject3d(EMData * image) const;

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Simple real-space projection. Most accurate.";
		}

		static Projector *NEW()
		{
			return new StandardProjector();
		}
		
		static const string NAME;
	};

	/** Fast real space projection using Bi-Linear interpolation. (C. Yang)
        */
	class ChaoProjector:public Projector
	{
	  public:
		EMData * project3d(EMData * vol) const;
		EMData * backproject3d(EMData * imagestack) const;

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Fast real space projection generation with Bi-Linear interpolation.";
		}

		static Projector *NEW()
		{
			return new ChaoProjector();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("transform", EMObject::TRANSFORM);
			d.put("origin_x",  EMObject::INT);
			d.put("origin_y",  EMObject::INT);
			d.put("origin_z",  EMObject::INT);
			d.put("anglelist", EMObject::FLOATARRAY);
			d.put("radius",    EMObject::FLOAT);
			return d;
		}
		
		static const string NAME;

	  private:
        	int getnnz(Vec3i volsize, int ri, Vec3i origin, int *nray, int *nnz) const;
        	int cb2sph(float *cube, Vec3i volsize, int ri, Vec3i origin, int nnz0, int *ptrs,
        			int *cord  , float *sphere) const;
        	int sph2cb(float *sphere, Vec3i volsize, int nray, int ri, int nnz0,
        			int   *ptrs  , int *cord, float *cube) const;
        	int fwdpj3(Vec3i volsize, int nray, int nnz  , float *dm,
        	           Vec3i origin, int ri  , int *ptrs,
        	           int *cord, float *x, float *y) const;
        	int bckpj3(Vec3i volsize, int nray, int nnz, float *dm,
        	           Vec3i origin, int ri, int *ptrs, int *cord,
        	           float *x, float *y) const;
        	int ifix(float a) const;
        	void setdm(vector<float> anglelist, string const angletype, float *dm) const;
	};

	template <> Factory < Projector >::Factory();

	void dump_projectors();
	map<string, vector<string> > dump_projectors_list();
}

#endif	//eman__projector_h__

/* vim: set ts=4 noet nospell: */
