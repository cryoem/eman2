/**
* $Id$
*/

/*
 * Author: Wen Jiang, 08/11/2004 (jiang12@purdue.edu)
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

#ifndef eman_pointarray_h_
#define eman_pointarray_h_

#include "emdata.h"
#include "transform.h"

#if defined NFFT || NFFT2
extern "C"
{
	#include "nfft.h"
	#include "utils.h"
}
#endif

#include <sys/stat.h>

namespace EMAN
{
	/** PointArray defines a double array of points with values in a 3D space. */
	class PointArray
	{
		public:
		enum Density2PointsArrayAlgorithm
		{
			PEAKS_SUB, PEAKS_DIV, KMEANS
		};

		public:
		PointArray();
		explicit PointArray(int nn);
		~PointArray();
		void zero();
		PointArray *copy();
		PointArray & operator=(PointArray & pa);
		size_t get_number_points() const;
		void set_number_points(size_t nn);
		bool read_from_pdb(const char *file);
		void save_to_pdb(const char *file);
		FloatPoint get_center();
		void center_to_zero();	// shift center to zero
		Region get_bounding_box();
		double *get_points_array();
		Vec3f get_vector_at(int i);	// get the ith vector in the list
		double get_value_at(int i);  // get the amplitude of the ith vector
		void set_vector_at(int i,Vec3f vec,double value);
		void set_vector_at(int i,vector<double>);
		void set_points_array(double *p);
		/** Returns all x,y,z triplets packed into a vector<float>
		*
		* @return All points packed into a vector<float>
		*/
		vector<float> get_points();

		/** Calculates a (symmetrized) distance matrix for the current PointArray
		*
		* @param sortrows if set, will sort the values in each row. The return will no longer be a true similarity matrix.
		* @return An EMData object containing the similarity matrix
		*/
		EMData *distmx(int sortrows=0);	// computes a distance matrix

		/** Will try to establish a 1-1 correspondence between points in
		*   two different PointArray objects (this and to). Returns
		*   a vector<int> where the index is addresses the points in 'this'
		*   and the value addresses points in 'to'. A value of -1 means there
		*   was no match for that point.
		*
		* @return A vector<int> with the mapping of points from 'this' to 'to'. e.g. - ret[2]=5 means point 2 in 'this' matches point 5 in 'to'
		*/
		vector<int> match_points(PointArray *to,float max_miss=-1.0);


		/** Aligns one PointArray to another in 2 dimensions
		*
		* @param to Another PointArray to align to
		* @param max_dist
		* @return a Transform3D to map 'this' to 'to'
		*/
		Transform3D *align_2d(PointArray *to,float max_dist);	// computes the optimal alignment between two (non-identical) sets of points

		/** Translationally aligns one PointArray to another in 2 dimensions
		*
		* @param to Another PointArray to align to
		* @param flags
		* @param dxhint
		* @param dyhint
		* @return A Pixel containing dx,dy and a quality factor (smaller better)
		*/
		vector<float> align_trans_2d(PointArray *to,int flags=0,float dxhint=0,float dyhint=0 );		// translational alignment of point sets in 2-D

		void mask(double rmax, double rmin = 0.0);
		void mask_asymmetric_unit(const string & sym);
		void transform(Transform3D transform);
		
		/** Does Transform*v as opposed to v*Transform (as in the transform function)
		 * @param transform an EMAN2 Transform object 
		*/
		void right_transform(const Transform& transform);
		
		void set_from(vector<float>);
		void set_from(PointArray * source, const string & sym = "", Transform3D *transform=0);
		void set_from(double *source, int num, const string & sym = "", Transform3D *transform=0);
		void set_from_density_map(EMData * map, int num, float thresh, float apix,
								  Density2PointsArrayAlgorithm mode = PEAKS_DIV);
		void sort_by_axis(int axis = 1);	// x,y,z axes = 0, 1, 2
		EMData *pdb2mrc_by_nfft(int map_size, float apix, float res);	// return real space 3-D map
		EMData *pdb2mrc_by_summation(int map_size, float apix, float res);	// return real space 3-D map
		EMData *projection_by_nfft(int image_size, float apix, float res = 0);	// return 2-D Fourier Transform
		EMData *projection_by_summation(int image_size, float apix, float res);	// return 2-D real space image
		void replace_by_summation(EMData *image, int i, Vec3f vec, float amp, float apix, float res); // changes a single Gaussian from the projection

		/** Optimizes a pointarray based on a set of projection images (EMData objects)
		 * This is effectively a 3D reconstruction algorithm.
		 *
		 * @author Steve Ludtke  11/27/2004
		 * @param proj A vector of EMData objects containing projections with orientations
		 * @param pixres Size of each Gaussian in pixels
		 */
		void opt_from_proj(const vector<EMData*> & proj,float pixres);

		private:
		double *points;
		size_t n;
	};
}

#endif
