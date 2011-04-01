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

#ifndef eman_reconstructor_template_h__
#define eman_reconstructor_template_h__ 1

#include "reconstructor.h"

namespace EMAN
{

	/** XYZReconstructor is a reconstructor template for defining new
     * reconstructors. Please add your own code at the proper place.
     *
     * 1) Replace all 'XYZ' with your new reconstructor name.
     * 2) Define the reconstructor parameter names and types in get_param_types().
     * 3) Implement the reconstructor in setup(), insert_slice(), and finish();
     */
	class XYZReconstructor:public Reconstructor
	{
	  public:
		XYZReconstructor();
		~XYZReconstructor();

	/** initialize the reconstructor
	 */
		void setup();

	/** insert each image slice to the reconstructor. You may call
	 * this function multiple times.
	 */
		int insert_slice(const EMData * const slice, const Transform & euler,const float weight);

	/** finish reconstruction and return the complete model.
	 */
		EMData *finish(bool doift);

		string get_name() const
		{
			return NAME;
		}
		string get_desc() const
		{
			return "XYZ Description";
		}

		static Reconstructor *NEW()
		{
			return new XYZReconstructor();
		}

	/** Add your reconstructor parameter names and types in
	 * get_param_types(). For available parameter types, please
	 * refer class EMObject.
	 * 
	 * As an example, XYZReconstructor has 3 parameters:
	 *    int size;
	 *    float patratio;
	 *    vector<float> snr;
	 */
		TypeDict get_param_types() const
		{
			TypeDict d;
			  d.put("size", EMObject::INT);
			  d.put("padratio", EMObject::FLOAT);
			  d.put("snr", EMObject::FLOATARRAY);
			  return d;
		}

		static const string NAME;
	  private:
		  EMData * image;
		int nx;
		int ny;
		int nz;
	};

	/** Add your new reconstructor to ReconstructorFactoryExt().
     */
//	class ReconstructorFactoryExt
//	{
//	  public:
//		ReconstructorFactoryExt()
//		{
//			Factory < Reconstructor >::add<XYZReconstructor>();
//		}
//	};
//
//	static ReconstructorFactoryExt rf_ext;
}


#endif
