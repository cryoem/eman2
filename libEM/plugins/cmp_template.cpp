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

#include "cmp_template.h"
#include "transform.h"


using namespace EMAN;

const string XYZCmp::NAME = "middle";

float XYZCmp::cmp(EMData * image, EMData * with) const
{
	validate_input_args(image, with);

#if 0
	int param1 = params["param1"];
	float param2 = params["param2"];
	// do something
#endif
	int nx=image->get_xsize();
	int ny=image->get_ysize();
	
	double ret=0;
	for (int i=nx/4; i<nx*3/4; i++) {
		for (int j=ny/4; j<ny*3/4; j++) {
			ret+=image->get_value_at(i,j)*with->get_value_at(i,j);
		}
	}

	return (float) ret;
}
