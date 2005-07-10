/**
* $Id$
*/
#ifndef eman_pointarray_h_
#define eman_pointarray_h_

#include "ctf.h"
#include "emdata.h"

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
		PointArray(unsigned int nn);
		~PointArray();
		void zero();
		PointArray *copy();
		PointArray & operator=(PointArray & pa);
		unsigned int get_number_points();
		void set_number_points(unsigned int nn);
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
		void mask(double rmax, double rmin = 0.0);
		void mask_asymmetric_unit(const string & sym);
		void set_from(PointArray * source, const string & sym = "", Transform3D *transform=0);
		void set_from(double *source, unsigned int num, const string & sym = "", Transform3D *transform=0);
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
		 * @param proj A vector of EMData objects conatining projections with orientations
		 * @param pixres Size of each Gaussian in pixels
		 */
		void opt_from_proj(const vector<EMData*> & proj,float pixres);

		private:
		double *points;
		unsigned int n;
	};
}

#endif
