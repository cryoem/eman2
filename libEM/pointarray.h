/**
* $Id$
*/
#ifndef eman_pointarray_h_
#define eman_pointarray_h_

#include "transform.h"
#include "geometry.h"
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
		void mask_asymmetric_unit(string sym);
		void set_from(PointArray * source, string sym = "", Transform * transform = 0);
		void set_from(double *source, unsigned int num, string sym = "", Transform * transform = 0);
		void set_from_density_map(EMData * map, int num, float thresh, float apix,
								  Density2PointsArrayAlgorithm mode = PEAKS_DIV);
		void sort_by_axis(int axis = 1);	// x,y,z axes = 0, 1, 2
		EMData *pdb2mrc_by_nfft(int map_size, float apix, float res);	// return real space 3-D map
		EMData *pdb2mrc_by_summation(int map_size, float apix, float res);	// return real space 3-D map
		EMData *projection_by_nfft(int image_size, float apix, float res = 0);	// return 2-D Fourier Transform
		EMData *projection_by_summation(int image_size, float apix, float res);	// return 2-D real space image
		
		private:
		double *points;
		unsigned int n;
	};
}

#endif
