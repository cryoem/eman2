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

#ifndef eman_projector_template_h__
#define eman_projector_template_h__ 1

#include "projector.h"


namespace EMAN
{

	/** XYZProjector is an projector template for defining new
     * projectors. Please add your own code at the proper place.
     *
     * 1) Replace all 'XYZ' with your new projector name.
     * 2) Define the projector parameter names and types in get_param_types().
     * 3) Implement the projector in XYZProjector::project3d().
     */
	class XYZProjector:public Projector
	{
	  public:
		EMData * project3d(EMData * em) const;
// chao added
		EMData * backproject3d(EMData * em) const;

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "XYZ Description";
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
		
		static const string NAME;
	};

	/** Add your new projector to ProjectorFactoryExt().
     */
//	class ProjectorFactoryExt
//	{
//	  public:
//		ProjectorFactoryExt()
//		{
//			Factory < Projector >::add<XYZProjector>();
//		}
//	};
//
//
//	static ProjectorFactoryExt pf_ext;
}

#endif
