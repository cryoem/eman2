/**
 * $Id$
 */
#ifndef eman_pointarray_h_
#define eman_pointarray_h_

#include "transform.h"
#include "geometry.h"
#include "emdata.h"
extern "C"
{
#include "nfft.h"
#include "utils.h"
}
#include <sys/stat.h>

namespace EMAN
{
	/** PointArray defines a double array of points with values in a 3D space. */
	class PointArray
	{
	  public:
		enum Optimizer
		{
			ConjugateGradientFletcherReeves, ConjugateGradientPolakRibiere,
			ConjugateGradientBFGS, SimplexNelderMead
		};

		enum OptimizedParameters
		{
			Map = 1 << 0, Orientation = 1 << 1, Center = 1 << 2, Defocus = 1 << 3,
			Astigmatism = 1 << 4, BFactor = 1 << 5, Drift = 1 << 6, Scale = 1 << 7,
			Distortion = 1 << 8, BeamTilt = 1 << 9, DepthOfView = 1 << 10
		};

	  public:
		  PointArray();
		  PointArray(unsigned int nn);
		 ~PointArray();
		unsigned int get_number_points();
		void set_number_points(unsigned int nn);
		bool read_from_pdb(const char *file);
		FloatPoint get_center();
		void center_to_zero();	// shift center to zero
		Region get_bounding_box();
		double *get_points_array();
		void set_points_array(double *p);
		void set_from(PointArray * source, Transform * transform, string sym = "");
		void set_from(double *source, unsigned int num, Transform * transform, string sym = "");
		void sort_by_axis(int axis = 1);	// x,y,z axes = 0, 1, 2
		EMData *pdb2mrc_by_nfft(int map_size, float apix, float res);	// return real space 3-D map
		EMData *pdb2mrc_by_summation(int map_size, float apix, float res);	// return real space 3-D map
		EMData *projection_by_nfft(int image_size, float apix, float res = 0);	// return 2-D Fourier Transform
		EMData *projection_by_summation(int image_size, float apix, float res);	// return 2-D real space image
		bool refine(vector < EMData * >images, string sym = "", OptimizedParameters optparam =
					Map, Optimizer optimizer = ConjugateGradientBFGS);
	  private:
		double *points;
		unsigned int n;
	};
}

#endif
