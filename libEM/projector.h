#ifndef eman__projector_h__
#define eman__projector_h__ 1


#include "emobject.h"
#include <map>
#include <string>
#include <math.h>

using std::map;
using std::string;


namespace EMAN {

    class EMData;

    class Projector {
    public:
	virtual ~Projector() {}
	
	virtual Dict get_params() const { return params; }
	void set_params(const Dict& new_params) { params = new_params; }
	
	virtual EMData* project3d(EMData* em) const = 0;
	virtual TypeDict get_param_types() const = 0;
	virtual string get_name() const = 0;
	
    protected:
	Dict params;
    };

    // mode = 1,2,3,4,5,6,7
    class FFTProjector : public Projector {
    public:
	FFTProjector() : alt(0), az(0), phi(0) {}

	void set_params(const Dict& new_params)
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
	
	EMData* project3d(EMData* em) const;
    protected:
	virtual float get_gaussian_width() const = 0;
	virtual void interp_ft_3d(EMData* image, float x, float y, float z, float* data) const = 0;
	
	float alt, az, phi;
    };

    class RoundFFTProjector : public FFTProjector {
    public:
	string get_name() const { return "RoundFFTProjector"; }
	static Projector* NEW() { return new RoundFFTProjector(); }
    protected:
	float get_gaussian_width() const { return 0; }
	void interp_ft_3d(EMData* image, float x, float y, float z, float* data) const;
    };

    class Gaussian2FFTProjectorr : public FFTProjector {
    public:
	string get_name() const { return "Gaussian2FFTProjectorr"; }
	static Projector* NEW() { return new Gaussian2FFTProjectorr(); }
    protected:
	float get_gaussian_width() const { return 4.0/(M_PI*M_PI); }
	void interp_ft_3d(EMData* image, float x, float y, float z, float* data) const;
    };


    class Gaussian3FFTProjectorr : public FFTProjector {
    public:
	string get_name() const { return "Gaussian3FFTProjectorr"; }
	static Projector* NEW() { return new Gaussian3FFTProjectorr(); }
    protected:
	float get_gaussian_width() const { return 6.4/(M_PI*M_PI); }
	void interp_ft_3d(EMData* image, float x, float y, float z, float* data) const;
    };


    class Gaussian4FFTProjectorr : public FFTProjector {
    public:
	string get_name() const { return "Gaussian4FFTProjectorr"; }
	static Projector* NEW() { return new Gaussian4FFTProjectorr(); }
    protected:
	float get_gaussian_width() const { return 8.8/(M_PI*M_PI); }
	void interp_ft_3d(EMData* image, float x, float y, float z, float* data) const;
    };

    
    class Gaussian5FFTProjectorr : public FFTProjector {
    public:
	string get_name() const { return "Gaussian5FFTProjectorr"; }
	static Projector* NEW() { return new Gaussian5FFTProjectorr(); }
    protected:
	float get_gaussian_width() const { return 0; }
	void interp_ft_3d(EMData* image, float x, float y, float z, float* data) const;
    };


    class Gaussian6FFTProjectorr : public FFTProjector {
    public:
	string get_name() const { return "Gaussian6FFTProjectorr"; }
	static Projector* NEW() { return new Gaussian6FFTProjectorr(); }
    protected:
	float get_gaussian_width() const { return 10.4/(M_PI*M_PI); }
	void interp_ft_3d(EMData* image, float x, float y, float z, float* data) const;
    };


    class Gaussian7FFTProjectorr : public FFTProjector {
    public:
	string get_name() const { return "Gaussian7FFTProjectorr"; }
	static Projector* NEW() { return new Gaussian7FFTProjectorr(); }
    protected:
	float get_gaussian_width() const { return 10.4/(M_PI*M_PI); }
	void interp_ft_3d(EMData* image, float x, float y, float z, float* data) const;
    };

    
    // mode = -1
    class OptimizedProjector :  public Projector {
    public:
	string get_name() const { return "OptimizedProjector"; }
	static Projector* NEW() { return new OptimizedProjector(); }
	TypeDict get_param_types() const
	{
	    TypeDict d;
	    return d;
	}
	
	EMData* project3d(EMData* em) const;
    };

    // mode = -2
    class FastProjector :  public Projector {
    public:
	string get_name() const { return "FastProjector"; }
	static Projector* NEW() { return new FastProjector(); }
	TypeDict get_param_types() const
	{
	    TypeDict d;
	    return d;
	}
	
	EMData* project3d(EMData* em) const;
    };

    // mode = -3
    class FastSurfaceProjector :  public Projector {
    public:
	string get_name() const { return "FastSurfaceProjector"; }
	static Projector* NEW() { return new FastSurfaceProjector(); }
	TypeDict get_param_types() const
	{
	    TypeDict d;
	    return d;
	}
	
	EMData* project3d(EMData* em) const;
    };
    
    // mode = -4
    class SlowAccurateProjector :  public Projector {
    public:
	string get_name() const { return "SlowAccurateProjector"; }
	static Projector* NEW() { return new SlowAccurateProjector(); }
	TypeDict get_param_types() const
	{
	    TypeDict d;
	    return d;
	}
	
	EMData* project3d(EMData* em) const;
    };

 // mode = -5
    class SlowAccurateYProjector :  public Projector {
    public:
	string get_name() const { return "SlowAccurateYProjector"; }
	static Projector* NEW() { return new SlowAccurateYProjector(); }
	TypeDict get_param_types() const
	{
	    TypeDict d;
	    return d;
	}
	
	EMData* project3d(EMData* em) const;
    };
    
// mode = -6
    class SlowAccurate2DProjector :  public Projector {
    public:
	string get_name() const { return "SlowAccurate2DProjector"; }
	static Projector* NEW() { return new SlowAccurate2DProjector(); }
	TypeDict get_param_types() const
	{
	    TypeDict d;
	    return d;
	}
	
	EMData* project3d(EMData* em) const;
    };



    ///////////////////////// 
    
    typedef Projector* (*ProjectorType)();
    
    class ProjectorFactory {
    public:
	static ProjectorFactory* instance();

	void add(ProjectorType projector);
	Projector* get(string projector_name);
	Projector* get(string projector_name, const Dict& params);
	vector<string> get_list();
	
    private:
	ProjectorFactory();
	ProjectorFactory(const ProjectorFactory& f);
	~ProjectorFactory();

	void force_add(ProjectorType projector);
	
	static ProjectorFactory* my_instance;
	map<string, ProjectorType> my_dict;
    };

}
    
#endif
