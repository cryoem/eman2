#ifndef eman_averager_template_h__
#define eman_averager_template_h__ 1

#include "averager.h"

namespace EMAN
{

	/** XYZAverager is an averager  template for defining new
	 * averagers. Please add your own code at the proper place.
	 *
	 * 1) Replace all 'XYZ' with your new averager name.
	 * 2) Define the averager parameter names and types in get_param_types().
	 * 3) Implement the averager in XYZAverager::align().
	 */
	class XYZAverager:public Averager
	{
	public:
		XYZAverager();
		void add_image( EMData * image);
		EMData * finish();

		string get_name() const
		{
			return "XYZ";
		}

		static Averager *NEW()
		{
			return new XYZAverager();
		}

		/** Add your averager parameter names and types in
		 * get_param_types(). For available parameter types, please
		 * refer class EMObject.
		 * 
		 * As an example, XYZAverager has 3 parameters:
		 *    EMData *param1;
		 *    int param2;
		 *    float param3;
		 */
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("param1", EMObject::EMDATA);
			d.put("param2", EMObject::INT);
			d.put("param3", EMObject::FLOAT);
			return d;
		}
	private:
		EMData *result;
	};

	/** Add your new averager to AveragerFactoryExt().
     */
	class AveragerFactoryExt
	{
	public:
		AveragerFactoryExt()
		{
			Factory < Averager >::add(&XYZAverager::NEW);
		}
	};

	static AveragerFactoryExt averager_factory_ext;
}

#endif
