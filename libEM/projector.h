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
	Projector();
	virtual ~Projector();

	virtual EMData* project3d(EMData* em) const = 0;
	Dict get_params() const;
	void set_params(const Dict& new_params);
	
	virtual string get_name() const = 0;
	
    private:
	Dict params;
    };


#define DEFINE_PROJECTOR(T) \
class T: public Projector { \
public: \
    EMData* project3d(EMData* em) const; \
    string get_name() const; \
    static Projector* NEW(); \
}; 

    DEFINE_PROJECTOR(FFTProjector); // 0
    DEFINE_PROJECTOR(OptimizedProjector);  // -1
    DEFINE_PROJECTOR(FastRealProjector); // -2
    DEFINE_PROJECTOR(SlowAccurateProjector); // -4
    DEFINE_PROJECTOR(FastSurfaceProjector); // -3

    typedef Projector* (*ProjectorType)();
    class ProjectorFactory {
    public:
	static ProjectorFactory* instance();

	void add(ProjectorType projector);
	Projector* get(string projector_name);
	Projector* get(string projector_name, const Dict& params);
    private:
	ProjectorFactory();
	ProjectorFactory(const ProjectorFactory& f);

	static ProjectorFactory* f_instance;
	static map<string, ProjectorType> projector_dict;
    };

}
    
#endif
