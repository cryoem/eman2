#ifndef eman_averager_h__
#define eman_averager_h__ 1

namespace EMAN {

    class EMData;
    
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
	static Cmp* NEW() { return new Averager(); }

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    return d;
	}
    };

    class IterationAverager : public Averager {
    public:
	EMData* average(const vector<EMData*>& image_list) const;
	string get_name() const { return "Iteration"; }
	static Cmp* NEW() { return new Averager(); }

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
	EMData* sf;
	EMData* curves;
	vector<float> snr;
	string outfile;
    };

    
    class WeightingAverager : public CtfAverager {
    public:
	
	string get_name() const { return "Weighting"; }
	static Cmp* NEW() { return new WeightingAverager(); }

	void set_params(const Dict& new_params)
	{
	    params = new_params;
	    curves = params["curves"].get_EMData();
	    sf = params["sf"].get_XYData();
	}
	
	TypeDict get_param_types() const
	{
	    TypeDict d;
	    return d;
	}
    };

    

    class CtfCAverager : public CtfAverager {
    public:
	string get_name() const { return "CtfC"; }
	static Cmp* NEW() { return new Averager(); }

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    return d;
	}
    };
    


    class CtfCWAverager : public CtfAverager {
    public:
	string get_name() const { return "CtfCW"; }
	static Cmp* NEW() { return new CtfCWAverager(); }
	
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
	static Cmp* NEW() { return new CtfCWautoAverager(); }

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
