#ifndef eman_cmp_template_h__
#define eman_cmp_template_h__ 1

#include "cmp.h"

namespace EMAN
{

	/** XYZCmp is a cmp template for defining new
     * cmps. Please add your own code at the proper place.
     *
     * 1) Replace all 'XYZ' with your new cmp name.
     * 2) Define the cmp parameter names and types in get_param_types().
     * 3) Implement the cmp in XYZCmp::cmp().
     */
	class XYZCmp:public Cmp
	{
	  public:
		float cmp(EMData * image, EMData * with) const;

		string get_name() const
		{
			return "xYZ";
		}

		string get_desc() const
		{
			return "XYZ Description";
		}

		static Cmp *NEW()
		{
			return new XYZCmp();
		}

		/** Add your cmp parameter names and types in
		 * get_param_types(). For available parameter types, please
		 * refer class EMObject.
		 * 
		 * As an example, XYZCmp has 3 parameters:
		 *    EMData *with;
		 *    int param1;
		 *    float param2;
		 */
		TypeDict get_param_types() const
		{
			TypeDict d;
			  d.put("param1", EMObject::INT);
			  d.put("param2", EMObject::FLOAT);
			  return d;
		}
	};

	/** Add your new cmp to CmpFactoryExt().
     */
	class CmpFactoryExt
	{
	  public:
		CmpFactoryExt()
		{
			Factory < Cmp >::add(&XYZCmp::NEW);
		}
	};

	static CmpFactoryExt cf_ext;

}

#endif
