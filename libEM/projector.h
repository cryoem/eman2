/**
 * $Id$
 */
#ifndef eman__projector_h__
#define eman__projector_h__ 1

#include "emobject.h"
#include "transform.h"
#include <string>
#include <math.h>

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

		void set_params(const Dict & new_params)
		{
			Projector::set_params(new_params);
			alt = params["alt"];
			az = params["az"];
			phi = params["phi"];
		}

		string get_name() const
		{
			return "gauss_fft";
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
			d.put("alt", EMObject::FLOAT);
			d.put("az", EMObject::FLOAT);
			d.put("phi", EMObject::FLOAT);
			d.put("mode", EMObject::INT);
			return d;
		}

	  private:
		float alt, az, phi;
		void interp_ft_3d(int mode, EMData * image, float x, float y,
						  float z, float *data, float gauss_width) const;
	};


	/** Pawel Penczek's optimized projection routine.
     */
	class PawelProjector:public Projector
	{
	  public:
		EMData * project3d(EMData * image) const;

		string get_name() const
		{
			return "pawel";
		}

		string get_desc() const
		{
			return "Pawel Penczek's optimized real-space projection generation. Minor interpolation artifacts.";
		}

		static Projector *NEW()
		{
			return new PawelProjector();
		}

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

	/** fast real-space isosurface 3D proejection.
     */
	class SimpleIsoSurfaceProjector:public Projector
	{
	  public:
		EMData * project3d(EMData * image) const;

		string get_name() const
		{
			return "simple_iso_surface";
		}

		string get_desc() const
		{
			return "Simple isosurface rendering, not projections.";
		}

		static Projector *NEW()
		{
			return new SimpleIsoSurfaceProjector();
		}
	};

	/** Fast real-space 3D projection.
     */
	class StandardProjector:public Projector
	{
	  public:
		EMData * project3d(EMData * image) const;

		string get_name() const
		{
			return "standard";
		}

		string get_desc() const
		{
			return "Simple real-space projection. Most accurate.";
		}

		static Projector *NEW()
		{
			return new StandardProjector();
		}
	};

	/** Real-space 3D projection with trilinear interpolation.
     * It accumulates the results directly in the 2D image
     * instead of storing in another rotated copy then accumulating
     * (ie: less memory requirement). It does not modify the self
     * data either.
     */
	class StandardFastProjector:public Projector
	{
	  public:
		EMData * project3d(EMData * image) const;

		string get_name() const
		{
			return "standardfast";
		}

		string get_desc() const
		{
			return "Standard with some optimizations. Might be minor accuracy difference.";
		}

		static Projector *NEW()
		{
			return new StandardFastProjector();
		}
	};

	template <> Factory < Projector >::Factory();

	void dump_projectors();
}

#endif
