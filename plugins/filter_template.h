#ifndef eman_filter_template_h__
#define eman_filter_template_h__ 1

#include "filter.h"

namespace EMAN {
    
    /** XYZFilter is a filter template for defining new
     * filters. Please add your own code at the proper place.
     *
     * 1) Replace all 'XYZ' with your new filter name.
     * 2) Define the filter parameter names and types in get_param_types().
     * 3) Implement the filter in XYZFilter::process().
     */
    class XYZFilter : public Filter {
    public:
	void process(EMData *image);

	string get_name() const
	{
	    return "XYZ";
	}

	static Filter *NEW()
	{
	    return new XYZFilter();
	}

	/** Add your filter parameter names and types in
	 * get_param_types(). For available parameter types, please
	 * refer class EMObject.
	 * 
	 * As an example, XYZFilter has 2 parameters:
	 *    int value1;
	 *    float value2;
	 */
	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("value1", EMObject::INT);
	    d.put("value2", EMObject::FLOAT);
	    return d;
	}
    };


    /** Add your new filter to FilterFactoryExt().
     */
    class FilterFactoryExt {
    public:
	FilterFactoryExt() {
	    Factory<Filter>::add(&XYZFilter::NEW);
	}
    };

    static FilterFactoryExt filter_factory_ext;

}


#endif
