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
 
#ifndef eman_boxingtools_h__
#define eman_boxingtools_h__ 1

#include "emdata.h"
#include "geometry.h"
#include <vector>
using std::vector;

namespace EMAN
{
	/** BoxingTools is class for encapsulating common 
	 * boxing operations that may become expensive if they are implemented in 
	 * python.
	 */
	class BoxingTools
	{
	public:
		BoxingTools() {}
		~BoxingTools() {}
		
		/** Gets a pixel minimum delta radial profile about some pixel focal point.
		 * Useful for automated correlation-based boxing.
		 *
		 * @param image the image containing the interesting pixels values (typically a correlation image)
		 * @param x the x coordinate of the pixel focal point
		 * @param y the y coordinate of the pixel focal point
		 * @param radius the constraining radius of the minimum delta profile
		 * @author David Woolford
		 * @date April 2008
		*/
		static vector<float> get_min_delta_profile(const EMData* const image, int x, int y, int radius);
		
		/** Determines if a given pixel is the maximum in local radial neighborhood
		 * Useful for automated correlation-based boxing.
		 *
		 * @param image the image containing the interesting pixels values (typically a correlation image)
		 * @param x the x coordinate of the candidate pixel maximum 
		 * @param y the y coordinate of the candidate pixel maximum 
		 * @param radius the constraining radius of the local neighborhood
		 * @param efficiency an EMData object with the same dimensions as the main image. If a local maximum is found, the pixels in the radial neighborhood that was queried are set to 0. This can be exploited by the calling function to minimize queries. 
		 * @author David Woolford
		 * @date April 2008
		 */
		static bool is_local_maximum(const EMData* const image, int x, int y, int radius, EMData* const efficiency_map);
		
		
		static vector<IntPoint> auto_correlation_pick(const EMData* const image, float threshold, int radius, const vector<float>& profile, EMData* const efficiency);
		
		static bool hi_brid(const EMData* const image, int x, int y, int radius,EMData* const efficiency_map, vector<float>& profile);
		
		static void set_radial_zero(EMData* const efficiency, int x, int y, int radius);
		
		static IntPoint find_radial_max(const EMData* const map, int x, int y, int radius);
	private:
		
	};

} // namespace EMAN

#endif // eman_boxingtools_h__
