#ifndef eman_averager_h__
#define eman_averager_h__ 1

#include <vector>
using std::vector;

#include "emobject.h"

namespace EMAN {

    class EMData;
    class XYData;
    
    class Averager {
    public:
	virtual ~Averager() {}
	virtual EMData* average(const vector<EMData*>& image_list) const = 0;
	virtual void set_params(const Dict& new_params) { params = new_params; }
	virtual TypeDict get_param_types() const = 0;
	virtual string get_name() const = 0;
	
    protected:
	mutable Dict params;
    };


    class ImageAverager : public Averager {
    public:
	EMData* average(const vector<EMData*>& image_list) const;
	string get_name() const { return "Image"; }
	static Averager* NEW() { return new ImageAverager(); }

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("sigma", EMObject::EMDATA);
	    d.put("ignore0", EMObject::INT);
	    return d;
	}
    };

    class IterationAverager : public Averager {
    public:
	EMData* average(const vector<EMData*>& image_list) const;
	string get_name() const { return "Iteration"; }
	static Averager* NEW() { return new IterationAverager(); }

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    return d;
	}
    };


    class CtfAverager : public Averager {
    public:
	EMData* average(const vector<EMData*>& image_list) const;
    protected:
	CtfAverager() : sf(0), curves(0) {}
	XYData* sf;
	EMData* curves;
	vector<float> snr;
	string outfile;
    };

    
    class WeightingAverager : public CtfAverager {
    public:
	
	string get_name() const { return "Weighting"; }
	static Averager* NEW() { return new WeightingAverager(); }

	void set_params(const Dict& new_params)
	{
	    params = new_params;
	    curves = params["curves"].get_emdata();
	    sf = params["sf"].get_xydata();
	}
	
	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("curves", EMObject::EMDATA);
	    d.put("sf", EMObject::XYDATA);
	    return d;
	}
    };

    

    class CtfCAverager : public CtfAverager {
    public:
	string get_name() const { return "CtfC"; }
	static Averager* NEW() { return new CtfCAverager(); }

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    return d;
	}
    };
    


    class CtfCWAverager : public CtfAverager {
    public:
	string get_name() const { return "CtfCW"; }
	static Averager* NEW() { return new CtfCWAverager(); }
	
	void set_params(const Dict& new_params)
	{
	    params = new_params;
	    snr = params["snr"].get_farray();
	}
	
	TypeDict get_param_types() const
	{
	    TypeDict d;
	    return d;
	}
    };

    
    class CtfCWautoAverager : public CtfAverager {
    public:
	string get_name() const { return "CtfCWauto"; }
	static Averager* NEW() { return new CtfCWautoAverager(); }

	void set_params(const Dict& new_params)
	{
	    params = new_params;
	    outfile = params["outfile"].get_string();
	}
	
	TypeDict get_param_types() const
	{
	    TypeDict d;
	    return d;
	}
    };
    
}


#endif
