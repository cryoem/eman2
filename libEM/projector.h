/**
 * $Id$
 */
#ifndef eman__projector_h__
#define eman__projector_h__ 1


#include "emobject.h"
#include <map>
#include <string>
#include <math.h>

using std::map;
using std::string;


namespace EMAN
{
    class EMData;

    /** Projector class is the base class for all 3D projectors.
     * Each specific projector has a unique name and should be called
     * through the name.
     *
     * Typical usage of Projectors:
     *
     * 1. How to get all the Projector types
     *
     *    vector<string> all_projectors = Factory<Projector>.instance()->get_list();
     *
     * 2. How to use a Projector
     *
     *    EMData* img = ...;
     *    Projector* proj = Factory<Projector>.instance()->get("FFT");
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
	virtual ~Projector()
	{
	}

	virtual EMData *project3d(EMData * em) const = 0;
	
	virtual string get_name() const = 0;
	
	virtual Dict get_params() const
	{
	    return params;
	}
	
	void set_params(const Dict & new_params)
	{
	    params = new_params;
	}

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
     * 0;
     * 4.0 / (M_PI * M_PI);
     * 6.4 / (M_PI * M_PI);
     * 8.8 / (M_PI * M_PI);
     * 0;
     * 10.4 / (M_PI * M_PI);
     * 10.4 / (M_PI * M_PI);
     */

    class GaussFFTProjector : public Projector
    {
    public:
	GaussFFTProjector() : alt(0), az(0), phi(0)
	{
	}

	EMData *project3d(EMData * em) const;
	
	void set_params(const Dict & new_params)
	{
	    Projector::set_params(new_params);
	    alt = params["alt"].get_float();
	    az = params["az"].get_float();
	    phi = params["phi"].get_float();
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
    class PawelProjector : public Projector
    {
    public:
	EMData *project3d(EMData * em) const;
	
	string get_name() const
	{
	    return "Pawel";
	}

	static Projector *NEW()
	{
	    return new PawelProjector();
	}	
    };

    /** fast real-space isosurface 3D proejection.
     */
    class SimpleIsoSurfaceProjector : public Projector
    {
    public:
	EMData *project3d(EMData * em) const;
	
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
    class StandardProjector : public Projector
    {
    public:
	EMData *project3d(EMData * em) const;
	
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
    class StandardBigProjector : public Projector
    {
    public:
	EMData *project3d(EMData * em) const;
	
	string get_name() const
	{
	    return "StandardBig";
	}
	
	static Projector *NEW()
	{
	    return new StandardBigProjector();
	}
    };

    template<> Factory<Projector>::Factory();
}

#endif
