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
			return NAME;
		}

		string get_desc() const
		{
			return "XYZ description";
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
		
		static const string NAME;
	private:
		EMData *result;
	};

	/** Add your new averager to AveragerFactoryExt().
     */
//	class AveragerFactoryExt
//	{
//	public:
//		AveragerFactoryExt()
//		{
//			Factory < Averager >::add<XYZAverager>();
//		}
//	};
//
//	static AveragerFactoryExt averager_factory_ext;
}

#endif
