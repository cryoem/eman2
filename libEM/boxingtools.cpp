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
			
			// Protect against accessing pixels out of bounds
			if ( xx >= image->get_xsize() || xx < 0 ) continue;
			if ( yy >= image->get_ysize() || yy < 0 ) continue;
			
			// We don't need to pay attention to the origin
			if ( xx == x and yy == y) continue;
			
			// Protect against vector accesses beyond the boundary
			int square_length = k*k + j*j;
			if (square_length > radius_squared ) continue;
			
			// The idx is the radius, rounded down. This creates a certain type of pattern that
			// can only really be explained visually...
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
			if ( xx >= image->get_xsize() || xx < 0 ) continue;
			if ( yy >= image->get_ysize() || yy < 0 ) continue;
			
			// We don't need to pay attention to the origin
			if ( xx == x and yy == y) continue;
			
			if ((k*k+j*j)>radius_squared) continue;
			
			if ( image->get_value_at(xx,yy) > peakval)  return false;
		}
	}
	
	set_radial_zero(efficiency_map,x,y,2*radius);
	
	return true;
	
}

vector<IntPoint> BoxingTools::auto_correlation_pick(const EMData* const image, float threshold, int radius, const vector<float>& profile, EMData* const efficiency)
{
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	
	vector<IntPoint> solution;

	int r = radius+1;
	
	for(int j = r; j < ny-r;++j) {
		for(int k = r; k < nx-r;++k) {
			
			if (efficiency->get_value_at(k,j) == 0) continue;
			
			if (image->get_value_at(k,j) < threshold) continue;
			
			vector<float> p(r,0);
						
			if (hi_brid(image,k,j,r,efficiency,p)) {
				if (p[radius] >= profile[radius]) {
					solution.push_back(IntPoint(k,j));
				}
			}
					
			// remove if smart is used
// 			set_radial_zero(efficiency,k,j,2*radius);
		}
	}
	
	return solution;
}


bool BoxingTools::hi_brid(const EMData* const image, int x, int y, int radius,EMData* const efficiency_map, vector<float>& profile)
{
	float peakval = image->get_value_at(x,y);
	
	int radius_squared = radius*radius;
	for(int k = -radius; k <= radius; ++k) {
		for(int j = -radius; j <= radius; ++j) {
			// Translate coordinates
			int xx = x+j;
			int yy = y+k;
			
			// Protect against accessing pixels out of bounds
			if ( xx >= image->get_xsize() || xx < 0 ) continue;
			if ( yy >= image->get_ysize() || yy < 0 ) continue;
			
			// We don't need to pay attention to the origin
			if ( xx == x and yy == y) continue;
			
			// Protect against vector accesses beyond the boundary
			int square_length = k*k + j*j;
			if (square_length > radius_squared ) continue;
			
			// It's not a local maximum!
			if ( image->get_value_at(xx,yy) > peakval)  return false;
			
			// The idx is the radius, rounded down. This creates a certain type of pattern that
			// can only really be explained visually...
			int idx = (int) sqrtf(k*k + j*j);
			// decrement the idx, because the origin information is redundant
			idx -= 1;
			
			// Finally, get the drop relative to the origin
			float val = peakval - image->get_value_at(xx,yy);
			
			// Store it if the drop is smaller than the current value (or if there is no value)
			if ( profile[idx] > val || profile[idx] == 0 ) profile[idx] = val;
			
		}
	}
	
	set_radial_zero(efficiency_map,x,y,2*radius);
	
	return true;
}


void BoxingTools::set_radial_zero(EMData* const efficiency, int x, int y, int radius)
{
	int radius_squared = radius*radius;
	for(int k = -radius; k <= radius; ++k) {
		for(int j = -radius; j <= radius; ++j) {
			// Translate coordinates
			int xx = x+j;
			int yy = y+k;
			
			if ((k*k+j*j)>radius_squared) continue;
			// Protect against accessing pixel out of bounds
			if ( xx >= efficiency->get_xsize() || xx < 0 ) continue;
			if ( yy >= efficiency->get_ysize() || yy < 0 ) continue;
			
			efficiency->set_value_at(xx,yy,0);
		}
	}
}

IntPoint BoxingTools::find_radial_max(const EMData* const map, int x, int y, int radius)
{
	float currentmax = map->get_value_at(x,y);
	
	IntPoint soln(x,y);
	
	int radius_squared = radius*radius;
	for(int k = -radius; k <= radius; ++k) {
		for(int j = -radius; j <= radius; ++j) {
			// Translate coordinates
			int xx = x+j;
			int yy = y+k;
			
			// Protect against accessing pixels out of bounds
			if ( xx >= map->get_xsize() || xx < 0 ) continue;
			if ( yy >= map->get_ysize() || yy < 0 ) continue;
			
			// Protect against vector accesses beyond the boundary
			int square_length = k*k + j*j;
			if (square_length > radius_squared ) continue;
			
			float val = map->get_value_at(xx,yy);
			
			if (val > currentmax) {
				currentmax = val;
				soln[0] = xx;
				soln[1] = yy;
			}
		}
	}
	
	return soln;
}
