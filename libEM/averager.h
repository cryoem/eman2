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


    class IterationAverager : public Averager {
    public:
	EMData* average(const vector<EMData*>& image_list) const;
	string get_name() const { return "Averager"; }
	static Cmp* NEW() { return new Averager(); }

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    return d;
	}
    };

    
    class WeightingAverager : public Averager {
    public:
	EMData* average(const vector<EMData*>& image_list) const;
	string get_name() const { return "Averager"; }
	static Cmp* NEW() { return new Averager(); }

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    return d;
	}
    };

    

    class CtfCAverager : public Averager {
    public:
	EMData* average(const vector<EMData*>& image_list) const;
	string get_name() const { return "Averager"; }
	static Cmp* NEW() { return new Averager(); }

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    return d;
	}
    };

    

    class CtfCWautoAverager : public Averager {
    public:
	EMData* average(const vector<EMData*>& image_list) const;
	string get_name() const { return "Averager"; }
	static Cmp* NEW() { return new Averager(); }

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    return d;
	}
    };

    

    class CtfCWAverager : public Averager {
    public:
	EMData* average(const vector<EMData*>& image_list) const;
	string get_name() const { return "Averager"; }
	static Cmp* NEW() { return new Averager(); }

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    return d;
	}
    };

    
    
}


#endif
