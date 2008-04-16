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


#include "boxingtools.h"

using namespace EMAN;


vector<float> BoxingTools::get_min_delta_profile(const EMData* const image, int x, int y, int radius)
{
	float peakval = image->get_value_at(x,y);
	
	vector<float> profile(radius,0); // this sets the vectors size to radius, and the values to 0
	int radius_squared = radius*radius;
	for(int k = -radius; k <= radius; ++k) {
		for(int j = -radius; j <= radius; ++j) {
			// Translate coordinates
			int xx = x+j;
			int yy = y+k;
			
			// Protect against accessing pixel out of bounds
			if ( xx >= image->get_xsize() || xx < 0 ) continue;
			if ( yy >= image->get_ysize() || yy < 0 ) continue;
			
			// We don't need to pay attention to the origin
			if ( xx == x and yy == y) continue;
			
			// The idx is the radius, rounded down. This creates a certain type of pattern that
			// can only really be explained visually...
			
			// Protect against vector accesses beyond the boundary
			int square_length = k*k + j*j;
			if (square_length > radius_squared ) continue;
			int idx = (int) sqrtf(k*k + j*j);
			// decrement the idx, because the origin information is redundant
			idx -= 1;
			
			// Finally, get the drop relative to the origin
			float val = peakval - image->get_value_at(xx,yy);
			
			// Store it if the drop is smaller than the current value (or if there is no value)
			if ( profile[idx] > val || profile[idx] == 0 ) profile[idx] = val;
			
		}
	}
	
	return profile;
}

bool BoxingTools::is_local_maximum(const EMData* const image, int x, int y, int radius,EMData* efficiency_map)
{
	float peakval = image->get_value_at(x,y);
	int radius_squared = radius*radius;
	for(int k = -radius; k <= radius; ++k) {
		for(int j = -radius; j <= radius; ++j) {
			// Translate coordinates
			int xx = x+j;
			int yy = y+k;
			
			
			
// 			// Protect against accessing pixel out of bounds
// 			if ( xx >= image->get_xsize() || xx < 0 ) continue;
// 			if ( yy >= image->get_ysize() || yy < 0 ) continue;
			
			// We don't need to pay attention to the origin
			if ( xx == x and yy == y) continue;
			
			if ((k*k+j*j)>radius_squared) continue;
			
			if ( image->get_value_at(xx,yy) > peakval)  return false;
		}
	}
	
	for(int k = -radius; k <= radius; ++k) {
		for(int j = -radius; j <= radius; ++j) {
			// Translate coordinates
			int xx = x+j;
			int yy = y+k;
			
			if ((k*k+j*j)>radius_squared) continue;
			// Protect against accessing pixel out of bounds
// 			if ( xx >= image->get_xsize() || xx < 0 ) continue;
// 			if ( yy >= image->get_ysize() || yy < 0 ) continue;
// 			
			efficiency_map->set_value_at(xx,yy,0);
		}
	}
	
	return true;
	
}

