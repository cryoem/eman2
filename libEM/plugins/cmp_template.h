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

#ifndef eman_cmp_template_h__
#define eman_cmp_template_h__ 1

#include "emdata.h"
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
			return NAME;
		}

		string get_desc() const
		{
			return "dot product using the center of the image";
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
			  d.put("radius", EMObject::INT);
			  d.put("param2", EMObject::FLOAT);
			  return d;
		}
		
		static const string NAME;
	};

	/** Add your new cmp to CmpFactoryExt().
     */
//	class CmpFactoryExt
//	{
//	  public:
//		CmpFactoryExt()
//		{
//			Factory < Cmp >::add<XYZCmp>();
//		}
//	};
//
//	static CmpFactoryExt cf_ext;

}

#endif
