#ifndef eman_filter_template_h__
#define eman_filter_template_h__ 1

#include "filter.h"

namespace EMAN
{

	/** XYZProcessor is a processor template for defining new
     * processors. Please add your own code at the proper place.
     *
     * 1) Replace all 'XYZ' with your new processor name.
     * 2) Define the processor parameter names and types in get_param_types().
     * 3) Implement the processor in XYZProcessor::process().
     */
	class XYZProcessor:public Processor
	{
	public:
		void process(EMData * image);

		string get_name() const
		{
			return "XYZ";
		}

		static Processor *NEW()
		{
			return new XYZProcessor();
		}

		string get_desc() const
		{
			return "add your documentation here.";
		}
		
		/** Add your processor parameter names and types in
		 * get_param_types(). For available parameter types, please
		 * refer class EMObject.
		 * 
		 * As an example, XYZProcessor has 2 parameters:
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


	/** Add your new processor to FilterFactoryExt().
     */
	class FilterFactoryExt
	{
	public:
		FilterFactoryExt()
		{
			//Factory < Processor >::add(&XYZProcessor::NEW);
		}
	};

	static FilterFactoryExt filter_factory_ext;

}


#endif
