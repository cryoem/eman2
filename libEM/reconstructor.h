#ifndef eman_reconstructor_h__
#define eman_reconstructor_h__ 1

#include <vector>
#include <string>
#include <math.h>
#include "emobject.h"

using std::vector;
using std::string;

namespace EMAN {

    // these magic numbers and macros need to be removed later
    
#define I5G	(10.4/(M_PI*M_PI))	// used for 5x5x5 interpolation
#define I4G	(8.8/(M_PI*M_PI))	// used for 4 interpolation
#define I3G	(6.4/(M_PI*M_PI))	// used for 3 and 5x5x5 interpolation
#define I2G	(4.0/(M_PI*M_PI))	// 2 interpolation
#define CTFOS	5			// curve oversampling,
					// generates ny*CTFOS/2 points
    
    class EMData;
    class Rotation;

    class Reconstructor {
    public:
	virtual ~Reconstructor() {}
	virtual int setup() = 0;
	virtual int insert_slice(EMData* slice, const Rotation& euler) = 0;
	virtual EMData* finish() = 0;
	
	virtual string get_name() const = 0;
	
	virtual Dict get_params() const { return params; }
	virtual void set_params(const Dict& new_params) { params = new_params; }
	virtual TypeDict get_param_types() const  = 0;
	
    protected:
	mutable Dict params;
    };

    class FftReconstructor : public Reconstructor {
    public:
	FftReconstructor();
	~FftReconstructor();
	
	int setup();
	int insert_slice(EMData* slice, const Rotation& euler);
	EMData* finish();

	string get_name() const { return "FftReconstructor"; }
	static Reconstructor* NEW() { return new FftReconstructor(); }

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("size", EMObject::INT);
	    d.put("mode", EMObject::INT);
	    d.put("weight", EMObject::FLOAT);
	    d.put("dlog", EMObject::INT);
	    return d;
	}
	
    private:
	EMData* image;
	int nx;
	int ny;
	int nz;
    };


    class WfFftReconstructor : public Reconstructor {
    public:
	WfFftReconstructor();
	~WfFftReconstructor();
	
	int setup();
	int insert_slice(EMData* slice, const Rotation& euler);
	EMData* finish();

	string get_name() const { return "WfFftReconstructor"; }
	static Reconstructor* NEW() { return new WfFftReconstructor(); }

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("size", EMObject::INT);
	    d.put("mode", EMObject::INT);
	    d.put("padratio", EMObject::FLOAT);
	    d.put("snr", EMObject::FLOATARRAY);
	    return d;
	}
	
    private:
	EMData* image;
	int nx;
	int ny;
	int nz;
    };

    class BpRealReconstructor : public Reconstructor {
    public:
	BpRealReconstructor();
	~BpRealReconstructor();
	
	int setup();
	int insert_slice(EMData* slice, const Rotation& euler);
	EMData* finish();

	string get_name() const { return "BpRealReconstructor"; }
	static Reconstructor* NEW() { return new BpRealReconstructor(); }

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("size", EMObject::INT);
	    d.put("weight", EMObject::FLOAT);
	    return d;
	}
	
    private:
	EMData* image;
	int nx;
	int ny;
	int nz;
    };

    
    
    typedef Reconstructor* (ReconstructorType)();

    class ReconstructorFactory {
    public:
	static ReconstructorFactory* instance();
	void add(ReconstructorType reconstructor);
	Reconstructor* get(string reconstructor_name);
	vector<string> get_list();
    };
    
    
};

#endif
