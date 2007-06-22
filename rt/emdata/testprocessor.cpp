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

#include "emdata.h"
using namespace EMAN;

/*
 *@author Wen Jiang
 *@date 2006-7-18
 */
 
int main()
{
	EMData *d;
	int size=512;
	d = new EMData();
	d->set_size(size,size);
	
	// create a random plane
	srand(time(0));	// get_frand() generates the same sequence of numbers everytime
	float xc=Util::get_frand(0,size), yc=Util::get_frand(0,size), zc=Util::get_frand(0,size);
	float norm[3] = {Util::get_frand(0,1),Util::get_frand(0,1),Util::get_frand(0,1)};
	for(int i=0; i<3; i++) norm[i]/=Util::hypot3(norm[0], norm[1], norm[2]); 
	
	for(int j=0; j<size; j++){
		for(int i=0; i<size; i++){
			d->set_value_at(i,j,0,-((i-xc)*norm[0]+(j-yc)*norm[1])/norm[2]+zc);
		}
	}
	float mean = d->get_attr("mean");
	
	float epsilon = 1e-5; // percent

	EMData* dmask=d->copy();
	dmask->process_inplace("mask.sharp", Dict("inner_radius", size/3-5, "outer_radius", size/3+5, "value", 0));
	
	d->process_inplace("filter.gradientPlaneRemover", Dict("mask", dmask));
	
	int err = 0;
	float max_residual = fabs(d->get_attr("maximum"));
	if ( max_residual > epsilon * mean) {
		err = 1;
		printf("FAILED: max residual=%g percent > error limit=%g\n", max_residual / mean, epsilon);
	}
	else {
		err = 0;
		printf("SUCCESS: max residual=%g percent < error limit=%g\n", max_residual / mean, epsilon);
	}

	return err;
}
	

	
		
		
