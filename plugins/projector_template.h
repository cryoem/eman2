#ifndef eman_projector_template_h__
#define eman_projector_template_h__ 1

#include "projector.h"


namespace EMAN {
    
    /** XYZProjector is an projector template for defining new
     * projectors. Please add your own code at the proper place.
     *
     * 1) Replace all 'XYZ' with your new projector name.
     * 2) Define the projector parameter names and types in get_param_types().
     * 3) Implement the projector in XYZProjector::project3d().
     */
    class XYZProjector : public Projector
    {
    public:
	EMData *project3d(EMData * em) const;
	
	string get_name() const
	{
	    return "XYZ";
	}

	static Projector *NEW()
	{
	    return new XYZProjector();
	}
	/** Add your projector parameter names and types in
	 * get_param_types(). For available parameter types, please
	 * refer class EMObject.
	 * 
	 * As an example, XYZProjector has 2 parameters:
	 *    float param1;
	 *    int param2;
	 */
	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("param1", EMObject::FLOAT);
	    d.put("param2", EMObject::INT);
	    return d;
	}
    };
    
    /** Add your new projector to ProjectorFactoryExt().
     */
    class ProjectorFactoryExt
    {
    public:
	ProjectorFactoryExt() {
	    Factory<Projector> *projectors = Factory<Projector>::instance();
	    projectors->add(&XYZProjector::NEW);
	}
    };


    static ProjectorFactoryExt pf_ext;
}

#endif
