/**
 * $Id$
 */

/*
 * Author: David Woolford, 04/15/2008 (woolford@bcm.edu)
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

#include <vector>
using std::vector;

#include <map>
using std::map;

#include <gsl/gsl_matrix.h>

namespace EMAN
{
	/** BoxingTools is class for encapsulating common
	 * boxing operations that may become expensive if
	 * they are implemented in python.
	 * @author David Woolford
	 * @date April 2008
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
		 * @param exclusion_map an EMData object with the same dimensions as the main image. If a local maximum is found, the pixels in the radial neighborhood that was queried are set to 0. This can be exploited by the calling function to minimize queries.
		 * @author David Woolford
		 * @date April 2008
		 */
		static bool is_local_maximum(const EMData* const image, int x, int y, int radius, EMData* const exclusion_map);


		static vector<IntPoint> auto_correlation_pick(const EMData* const image, float threshold, int radius, const vector<float>& profile, EMData* const exclusion,const int cradius, int mode=1);

		static bool hi_brid(const EMData* const image, int x, int y, int radius,EMData* const exclusion_map, vector<float>& profile);

		static void set_radial_non_zero(EMData* const exclusion, int x, int y, int radius);

		static IntPoint find_radial_max(const EMData* const map, int x, int y, int radius);

		static map<unsigned int, unsigned int> classify(const vector<vector<float> >& data, const unsigned int& classes = 4);

		static Vec3f get_color( const unsigned int index );

		static void set_region( EMData* const image, const EMData* const mask, const int x, const int y, const float& val );

		enum CmpMode {
			SWARM_DIFFERENCE,
			SWARM_RATIO,
   			SWARM_AVERAGE_RATIO
		};

		static void set_mode( const CmpMode m ) { mode = m; }
	private:
		// A vector to store the "seed" starting colors, which exist at the corners of the cube.
		// Then the vector is grown as more colours are asked for.
		static vector<Vec3f> colors;
		static CmpMode mode;
	};

	class BoxSVDClassifier
	{
		public:
			BoxSVDClassifier(const vector<vector<float> >& data, const unsigned int& classes = 4);

			~BoxSVDClassifier();

			map< unsigned int, unsigned int> go();

	// Because the results of the classification always give a random class number, I thought it
	// be easier for the user if the classes with the greatest numbers always had the same colour.
	//  This is achieved using this function
			static map< unsigned int, unsigned int> colorMappingByClassSize( const map< unsigned int, unsigned int>& grouping );
		private:
			const vector<vector<float> >& mData;

			unsigned int mColumns;
			unsigned int mRows;

			unsigned int mClasses;


			map< unsigned int, unsigned int> randomSeedCluster(const gsl_matrix* const svd_coords, unsigned int matrix_dims);
			map< unsigned int, unsigned int> getIterativeCluster(const gsl_matrix* const svd_coords, const map< unsigned int, unsigned int>& current_grouping);

			bool setDims( const vector<vector<float> >& data );

			vector<vector<float> > getDistances( const gsl_matrix* const svd_coords, const gsl_matrix* const ref_coords);

			map< unsigned int, unsigned int> getMapping(const vector<vector<float> >& distances);
	};

} // namespace EMAN

#endif // eman_boxingtools_h__
