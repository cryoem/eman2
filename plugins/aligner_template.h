#ifndef eman_aligner_template_h__
#define eman_aligner_template_h__ 1

#include "aligner.h"

namespace EMAN
{

	/** XYZAligner is an aligner  template for defining new
     * aligners. Please add your own code at the proper place.
     *
     * 1) Replace all 'XYZ' with your new aligner name.
     * 2) Define the aligner parameter names and types in get_param_types().
     * 3) Implement the aligner in XYZAligner::align().
     */
	class XYZAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img, const string & cmp_name) const;

		EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img);
		}
		
		string get_name() const
		{
			return "XYZ";
		}
		
		string get_desc() const
		{
			return "XYZ description";
		}

		static Aligner *NEW()
		{
			return new XYZAligner();
		}

		/** Add your aligner parameter names and types in
		 * get_param_types(). For available parameter types, please
		 * refer class EMObject.
		 * 
		 * As an example, XYZAligner has 3 parameters:
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


	/** Add your new aligner to AlignerFactoryExt().
     */
	class AlignerFactoryExt
	{
	  public:
		AlignerFactoryExt()
		{
			Factory < Aligner >::add(&XYZAligner::NEW);
		}
	};

	static AlignerFactoryExt aligner_factory_ext;
}

#endif
