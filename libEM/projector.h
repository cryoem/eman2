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
	 *   a) a unique name to idenfity itself in the factory.
	 *   b) a static method to register itself in the factory.
	 *
     * Typical usage of Projectors:
     *
     * 1. How to get all the Projector types
     *
     *    vector<string> all_projectors = Factory<Projector>::get_list();
     *
     * 2. How to use a Projector
     *
     *    EMData* img = ...;
     *    Projector* proj = Factory<Projector>::get("FFT");
     *    EMData* result = proj->project3d(img);
     *
     * 3. How to define a new Projector
     *
     *    A new projector type "XYZProjector" should implement at
     *    least the following 3 functions:
     *
     *        EMData *project3d(EMData * em) const;
     *        string get_name() const { return "XYZ"; }
     *        static Projector* NEW() { return new XYZProjector(); }
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
		/** Get filter parameter information in a dictionary. Each
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
			return "GaussFFT";
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
     * Subject to some small artifacts due to interpolation scheme used.
     */
	class PawelProjector:public Projector
	{
	  public:
		EMData * project3d(EMData * image) const;

		string get_name() const
		{
			return "Pawel";
		}

		static Projector *NEW()
		{
			return new PawelProjector();
		}

	  private:

		struct Pointers
		{
			Vec3 < int >location;
			int start;
			int end;
		};

	};

	/** fast real-space isosurface 3D proejection.
     */
	class SimpleIsoSurfaceProjector:public Projector
	{
	  public:
		EMData * project3d(EMData * image) const;

		string get_name() const
		{
			return "SimpleIsoSurface";
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
			return "Standard";
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
	class StandardBigProjector:public Projector
	{
	  public:
		EMData * project3d(EMData * image) const;

		string get_name() const
		{
			return "StandardBig";
		}

		static Projector *NEW()
		{
			return new StandardBigProjector();
		}
	};

	template <> Factory < Projector >::Factory();

	void dump_projectors();
}

#endif
