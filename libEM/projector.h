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

    // mode = 1,2,3,4,5,6,7
    class FFTProjector : public Projector
    {
    public:
	FFTProjector() : alt(0), az(0), phi(0)
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

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("alt", EMObject::FLOAT);
	    d.put("az", EMObject::FLOAT);
	    d.put("phi", EMObject::FLOAT);
	    return d;
	}
	
    protected:
	float alt, az, phi;
	
	virtual float get_gaussian_width() const = 0;
	virtual void interp_ft_3d(EMData * image, float x, float y, float z, float *data) const = 0;
    };

    class RoundFFTProjector : public FFTProjector
    {
    public:
	string get_name() const
	{
	    return "RoundFFTProjector";
	}
	
	static Projector *NEW()
	{
	    return new RoundFFTProjector();
	}
	
    protected:
	void interp_ft_3d(EMData * image, float x, float y, float z, float *data) const;
	
	float get_gaussian_width() const
	{
	    return 0;
	}
    };

    class Gaussian2FFTProjectorr : public FFTProjector
    {
    public:
	string get_name() const
	{
	    return "Gaussian2FFTProjectorr";
	}
	
	static Projector *NEW()
	{
	    return new Gaussian2FFTProjectorr();
	}
	
    protected:
	void interp_ft_3d(EMData * image, float x, float y, float z, float *data) const;
	
	float get_gaussian_width() const
	{
	    return 4.0 / (M_PI * M_PI);
	}
    };


    class Gaussian3FFTProjectorr : public FFTProjector
    {
    public:
	string get_name() const
	{
	    return "Gaussian3FFTProjectorr";
	}
	
	static Projector *NEW()
	{
	    return new Gaussian3FFTProjectorr();
	}
	
    protected:
	void interp_ft_3d(EMData * image, float x, float y, float z, float *data) const;
	
	float get_gaussian_width() const
	{
	    return 6.4 / (M_PI * M_PI);
	}
    };


    class Gaussian4FFTProjectorr : public FFTProjector
    {
    public:
	string get_name() const
	{
	    return "Gaussian4FFTProjectorr";
	}
	
	static Projector *NEW()
	{
	    return new Gaussian4FFTProjectorr();
	}
	
    protected:
	float get_gaussian_width() const
	{
	    return 8.8 / (M_PI * M_PI);
	}
	
	void interp_ft_3d(EMData * image, float x, float y, float z, float *data) const;
    };


    class Gaussian5FFTProjectorr : public FFTProjector
    {
    public:
	string get_name() const
	{
	    return "Gaussian5FFTProjectorr";
	}
	
	static Projector *NEW()
	{
	    return new Gaussian5FFTProjectorr();
	}
	
    protected:
	float get_gaussian_width() const
	{
	    return 0;
	}
	
	void interp_ft_3d(EMData * image, float x, float y, float z, float *data) const;
    };


    class Gaussian6FFTProjectorr : public FFTProjector
    {
    public:
	string get_name() const
	{
	    return "Gaussian6FFTProjectorr";
	}
	
	static Projector *NEW()
	{
	    return new Gaussian6FFTProjectorr();
	}
    protected:
	float get_gaussian_width() const
	{
	    return 10.4 / (M_PI * M_PI);
	}
	
	void interp_ft_3d(EMData * image, float x, float y, float z, float *data) const;
    };


    class Gaussian7FFTProjectorr : public FFTProjector
    {
    public:
	string get_name() const
	{
	    return "Gaussian7FFTProjectorr";
	}
	
	static Projector *NEW()
	{
	    return new Gaussian7FFTProjectorr();
	}
	
    protected:
	float get_gaussian_width() const
	{
	    return 10.4 / (M_PI * M_PI);
	}
	
	void interp_ft_3d(EMData * image, float x, float y, float z, float *data) const;
    };


    // mode = -1
    class OptimizedProjector : public Projector
    {
    public:
	EMData *project3d(EMData * em) const;
	
	string get_name() const
	{
	    return "OptimizedProjector";
	}

	static Projector *NEW()
	{
	    return new OptimizedProjector();
	}	
    };

    // mode = -2
    class FastProjector : public Projector
    {
    public:	
	EMData *project3d(EMData * em) const;
	
	string get_name() const
	{
	    return "FastProjector";
	}
	
	static Projector *NEW()
	{
	    return new FastProjector();
	}	
    };

    // mode = -3
    class FastSurfaceProjector : public Projector
    {
    public:
	EMData *project3d(EMData * em) const;
	
	string get_name() const
	{
	    return "FastSurfaceProjector";
	}
	
	static Projector *NEW()
	{
	    return new FastSurfaceProjector();
	}	
    };

    // mode = -4
    class SlowAccurateProjector : public Projector
    {
    public:
	EMData *project3d(EMData * em) const;
	
	string get_name() const
	{
	    return "SlowAccurateProjector";
	}
	
	static Projector *NEW()
	{
	    return new SlowAccurateProjector();
	}	
    };

    // mode = -5
    class SlowAccurateYProjector : public Projector
    {
    public:
	EMData *project3d(EMData * em) const;
	
	string get_name() const
	{
	    return "SlowAccurateYProjector";
	}
	
	static Projector *NEW()
	{
	    return new SlowAccurateYProjector();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    return d;
	}
    };

    // mode = -6
    class SlowAccurate2DProjector : public Projector
    {
    public:
	EMData *project3d(EMData * em) const;
	
	string get_name() const
	{
	    return "SlowAccurate2DProjector";
	}
	
	static Projector *NEW()
	{
	    return new SlowAccurate2DProjector();
	}
    };

}

#endif
