#ifndef eman__projector_h__
#define eman__projector_h__ 1


#include "emobject.h"
#include <map>
#include <string>


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

    // mode = 1
    class FFTProjector : public Projector {
    public:
	string get_name() const { return "FFTProjector"; }
	static Projector* NEW() { return new FFTProjector(); }
	TypeDict get_param_types() const
	{
	    TypeDict d;
	    return d;
	}
	
	EMData* project3d(EMData* em) const;
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
