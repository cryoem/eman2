/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
 * Copyright (c) 2000-2006 Baylor College of Medicine
 * 
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 * 
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * 
 * */

#ifndef eman_processor_template_h__
#define eman_processor_template_h__ 1

#include "processor.h"
#include "emdata.h"

namespace EMAN
{

	/** XYZProcessor is a processor template for defining new
     * processors. Please add your own code at the proper place.
     *
     * 1) Replace all 'XYZ' with your new processor name.
     * 2) Define the processor parameter names and types in get_param_types().
     * 3) Implement the processor in XYZProcessor::process_inplace().
     */
	class XYZProcessor:public Processor
	{
	public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return NAME;
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
		
		static const string NAME;
	};

	
// 	class SubstituteZeroPixelsProcessor:public Processor
// 	{
// 		public:
// 			void process_inplace(EMData * image);
// 
// 			string get_name() const
// 			{
// 				return "substitute.zeropixels";
// 			}
// 
// 			static Processor *NEW()
// 			{
// 				return new SubstituteZeroPixelsProcessor();
// 			}
// 
// 			string get_desc() const
// 			{
// 				return "Replaces zero pixels in the image with corresponding (coordinate-wise) pixel values in the EMData object specified in the parameters.";
// 			}
// 
// 			TypeDict get_param_types() const
// 			{
// 				TypeDict d;
// 				d.put("image", EMObject::EMDATA, "EMData object that contains useful image data");
// 				return d;
// 			}
// 	};

	/** Add your new processor to FilterFactoryExt().
     */
//	class FilterFactoryExt
//	{
//	public:
//		FilterFactoryExt()
//		{
//			Factory < Processor >::add<XYZProcessor>();
//// 			Factory < Processor >::add(&SubstituteZeroPixelsProcessor::NEW);
//		}
//	};
//
//	static FilterFactoryExt filter_factory_ext;

}


#endif	//eman_processor_template_h__
