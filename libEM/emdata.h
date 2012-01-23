/**
 * $Id$
 */

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

#ifndef eman__emdata_h__
#define eman__emdata_h__ 1

#ifdef _WIN32
	#pragma warning(disable:4819)
#endif	//_WIN32

#include <cfloat>
#include <complex>
#include <fstream>

#include "sparx/fundamentals.h"
#include "emutil.h"
#include "util.h"
#include "sparx/emarray.h"
#include "geometry.h"
#include "transform.h"
#ifdef EMAN2_USING_CUDA
#include <cuda_runtime_api.h>
#include "cuda/cuda_util.h"
#endif // EMAN2_USING_CUDA
using std::string;
using std::vector;
using std::map;
//using std::complex;	//comment this out for conflict with ACML
using std::ostream;

#include <utility>
using std::pair;

namespace EMAN
{
	class ImageIO;
	class Ctf;
	class XYData;
	class Transform;
	class GLUtil;

	typedef boost::multi_array_ref<float, 2> MArray2D;
	typedef boost::multi_array_ref<float, 3> MArray3D;
	typedef boost::multi_array_ref<std::complex<float>, 2> MCArray2D;
	typedef boost::multi_array_ref<std::complex<float>, 3> MCArray3D;
	typedef boost::multi_array<int, 2> MIArray2D;
	typedef boost::multi_array<int, 3> MIArray3D;

	/** @ingroup tested3c */
	/** EMData stores an image's data and defines core image processing routines.
     * The image is 1D, 2D or 3D, in real space or fourier space (complex image).
	 *
	 * Data are ordered with x increasing fastest, then y, then z.
	*/
	class EMData
	{
		friend class GLUtil;

		/** For all image I/O */
		#include "emdata_io.h"

		/** For anything read/set image's information */
		#include "emdata_metadata.h"

		/** For modular class functions, process, align, etc. */
		#include "emdata_modular.h"

		/** For fft, wavelet, insert data */
		#include "emdata_transform.h"

		/** For get/set values, basic math operations, operators */
		#include "emdata_core.h"

		/** This is the header of EMData stay in sparx directory */
		#include "sparx/emdata_sparx.h"
#ifdef EMAN2_USING_CUDA
		/** This is for CUDA related functionality */
		#include "emdata_cuda.h"
#endif // EMAN2_USING_CUDA

	public:
		enum FFTPLACE { FFT_OUT_OF_PLACE, FFT_IN_PLACE };
		enum WINDOWPLACE { WINDOW_OUT_OF_PLACE, WINDOW_IN_PLACE };

		/** Construct an empty EMData instance. It has no image data. */
		EMData();
		~ EMData();

		/** Construct from an image file.
		 * @param filename the image file name
		 * @param image_index the image index for stack image file, default 0 */
		explicit EMData(const string& filename, int image_index=0);

		/**# makes an image of the specified size, either real or complex.
		 * For complex image, the user would specify the real-space dimensions.
		 * @ingroup CUDA_ENABLED
		 * @param nx size for x dimension
		 * @param ny size for y dimension
		 * @param nz size for z dimension, default 1
		 * @param is_real boolean to specify real(true) or complex(false) image, default real */
		EMData(int nx, int ny, int nz=1, bool is_real=true);

		/** Construction from a data pointer, dimensions must be supplied.
		 * Takes possession of the pointer.
		 * data pointer must be allocated using malloc!
		 * @param data a pointer to the pixel data which is stored in memory. Takes possession
		 * @param nx the number of pixels in the x direction
		 * @param ny the number of pixels in the y direction
		 * @param nz the number of pixels in the z direction
		 * @param attr_dict attribute dictionary for this image
		 */
		EMData(float* data, const int nx, const int ny, const int nz, const Dict& attr_dict = Dict());
		
		/** Construction from a data pointer for usage in cuda, dimensions must be supplied.
		 * Takes possession of the pointer.
		 * data pointer must be allocated using malloc!
		 * @ingroup CUDA_ENABLED
		 * @param data a pointer to the pixel data which is stored in memory. Takes possession
		 * @param cudadata a pointer to the pixel data which is stored in cudamemory. Takes possession
		 * @param nx the number of pixels in the x direction
		 * @param ny the number of pixels in the y direction
		 * @param nz the number of pixels in the z direction
		 */
		EMData(float* data, float* cudadata, const int nx, const int ny, const int nz, const Dict& attr_dict = Dict());
		//I do not wish to use a default dict if none is provided. Copying Dicts arround is really expensive in CUDA.

		/** Construct from an EMData (copy constructor).
		 * Performs a deep copy
		 * @ingroup CUDA_ENABLED
		 * @param that the EMData to copy
		*/
		EMData(const EMData& that);

		/** EMData assignment operator
		 * Performs a deep copy
		 * @ingroup CUDA_ENABLED
		 * @param that the EMData to copy
		*/
		EMData& operator=(const EMData& that);


		/** Get an inclusive clip. Pads to fill if larger than this image.
		 * .
		 * @param area The clip area, can be 2D/3D
		 * @param fill the value to assign new pixels outside the area of the original image
		 * @exception ImageDimensionException if any of the dimensions of the argument region are negative
		 * @return The clip image.
		 */
		EMData *get_clip(const Region & area, const float fill = 0) const;

		/** Clip the image inplace - clipping region must be smaller than the current region
		 * internally memory is reallocated
		 * @exception ImageDimensionException if any of the dimensions of the argument region are negative
		 * @param area The clip area, can be 2D/3D.
		 * @param fill_value the value that new region
		 */
		void clip_inplace(const Region & area,const float& fill_value=0);

		/** Get the top half of this 3D image.
		 * @exception ImageDimensionException If this image is not 3D.
		 * @return The top half of this image.
		 */
		EMData *get_top_half() const;


		/** This will extract an arbitrarily oriented and sized region from the
		 *  image.
		 *
		 *  @param xform The transformation of the region.
		 *  @param size Size of the clip.
		 *  @param scale Scaling put on the returned image.
		 *  @return The clip image.
		 */
		EMData *get_rotated_clip(const Transform & xform, const IntSize &size, float scale=1.0);

		/** Window the center of an image.
		 *  Often an image is padded with zeros for fourier interpolation.  In
		 *  that case the desired lxlxl volume (or lxl area) lies in the center
		 *  of a larger volume (or area).  This routine creates a new object
		 *  that contains only the desired window.  (This routine is a thin
		 *  wrapper around get_clip.)
		 *
		 * @param l Length of the window.
		 * @return An image object that has been windowed.
		 */
		EMData* window_center(int l);


		/** Set up for fftslice operations.
		 * When interpolating in fourier space there is a little
		 * problem when we get close to x=0, since f(-x,-y,-z) = f(x,y,z)* .
		 * So this makes a supplementary array that allows for up to +-2
		 * point interpolation all the way to the origin in x.
		 *
		 * 3D only; complex image only
		 *
		 * @param redo If true,  recalculate the supplementary array.
		 * @exception ImageFormatException If the image is not a
		 * complex image.
		 * @exception ImageDimensionException If the image is not 3D.
		 * @return The supplementary array.
		 */
		float *setup4slice(bool redo = true);


		/** scale the image by a factor.
		 * @param scale_factor scale factor.
		 */
		void scale(float scale_factor);


		/** Translate this image.
		 * @param dx Translation distance in x direction.
		 * @param dy Translation distance in y direction.
		 * @param dz Translation distance in z direction.
		 */
		void translate(float dx, float dy, float dz);


		/** Translate this image.
		 * @param translation The translation distance vector.
		 */
		void translate(const Vec3f &translation);


		/** Translate this image. integer only translation
		 *  could be done faster, without interpolation.
		 * @param dx Translation distance in x direction.
		 * @param dy Translation distance in y direction.
		 * @param dz Translation distance in z direction.
		 */
		void translate(int dx, int dy, int dz);


		/** Translate this image. integer only translation
		 *  could be done faster, without interpolation.
		 * @param translation The translation distance vector.
		 */
		void translate(const Vec3i &translation);


		/** Rotate this image.
		 * DEPRECATED USE EMData::transform()
		 * @param t Transformation rotation.
		 */
		void rotate(const Transform & t);

		float max_3D_pixel_error(const Transform &t1, const Transform &t2, float r);

		/** Rotate this image.
		 * DEPRECATED USE EMData::Transform
		 * @param az  Rotation euler angle az  in EMAN convention.
		 * @param alt Rotation euler angle alt in EMAN convention.
		 * @param phi Rotation euler angle phi in EMAN convention.
		 */
		void rotate(float az, float alt, float phi);


		/** Rotate then translate the image.
		 * DEPRECATED USE EMData::Transform
		 * @param t The rotation and translation transformation to be done.
		 */
//		void rotate_translate(const Transform3D & t);

		/**  Transform the image
		 * @param t the transform object that describes the transformation to be applied to the image.
		 */
		inline void transform(const Transform& t) {
			ENTERFUNC;
			process_inplace("xform",Dict("transform",(Transform*)(&t)));
			//update(); no need, process_inplace did it
			EXITFUNC;
		}

		/** Apply a transformation to the image.
		 * DEPRECATED USE EMData::Transform
		 * @param t transform object that describes the transformation to be applied to the image.
		 */
		inline void rotate_translate(const Transform & t) {
			cout << "Deprecation warning. Please consider using EMData::transform() instead " << endl;
			transform(t); }

		/** Rotate then translate the image.
		 * DEPRECATED USE EMData::Transform
		 * @param az  Rotation euler angle az  in EMAN convention.
		 * @param alt Rotation euler angle alt in EMAN convention.
		 * @param phi Rotation euler angle phi in EMAN convention.
		 * @param dx Translation distance in x direction.
		 * @param dy Translation distance in y direction.
		 * @param dz Translation distance in z direction.
		 */
		void rotate_translate(float az, float alt, float phi, float dx, float dy, float dz);


		/** Rotate then translate the image.
		 * DEPRECATED USE EMData::Transform
		 * @param az  Rotation euler angle az  in EMAN convention.
		 * @param alt Rotation euler angle alt in EMAN convention.
		 * @param phi Rotation euler angle phi in EMAN convention.
		 * @param dx Translation distance in x direction.
		 * @param dy Translation distance in y direction.
		 * @param dz Translation distance in z direction.
		 * @param pdx Pretranslation distance in x direction.
		 * @param pdy Pretranslation distance in y direction.
		 * @param pdz Pretranslation distance in z direction.
		 */
		void rotate_translate(float az, float alt, float phi, float dx, float dy,
							  float dz, float pdx, float pdy, float pdz);


		/** This performs a translation of each line along x with wraparound.
		 *  This is equivalent to a rotation when performed on 'unwrapped' maps.
		 *  @param dx Translation distance align x direction.
		 *  @exception ImageDimensionException If the image is 3D.
		 */
		void rotate_x(int dx);


		/** Fast rotation by 180 degrees. Square 2D image only.
		 *  @exception ImageFormatException If the image is not square.
		 *  @exception ImageDimensionException If the image is not 2D.
		 */
		inline void rotate_180() {
			ENTERFUNC;
			process_inplace("math.rotate.180",Dict());
			EXITFUNC;
		}


		/** dot product of 2 images. Then 'this' image is rotated/translated.
		 * It is much faster than Rotate/Translate then dot product.
		 * 2D images only.
		 *
		 * @param with The image used to do the dot product.
		 * @param dx Translation distance in x direction.
		 * @param dy Translation distance in y direction.
		 * @param da Rotation euler angle in degrees
		 * @param mirror
		 * @exception ImageFormatException If the 2 images are not the same size.
		 * @exception ImageDimensionException If the image is 3D.
		 * @return
		 */
		double dot_rotate_translate(EMData * with, float dx, float dy, float da,const bool mirror=false);


		/** This does a normalized dot product of a little image with a big image
		 * using real-space methods. The result is the same size as 'this',
		 * but a border 1/2 the size of 'little_img' will be zero.
		 * This routine is only efficient when 'little_img' is fairly small.
		 *
		 * @param little_img A small image.
		 * @param do_sigma Calculate sigma or not.
		 * @exception ImageDimensionException If the image is not 1D/2D.
		 * @return normalized dot product image.
		 */
		EMData *little_big_dot(EMData * little_img, bool do_sigma = false);


		/** Radon Transform: an algorithm that transforms an original
		 * image into a series of equiangular projections. When
		 * applied to a 2D object, the output of the Radon transform is a
		 * series of 1D lines.
		 *
		 * Do radon transformation on this image. This image must be
		 * 2D square.

		 * @exception ImageFormatException If the image is not square.
		 * @exception ImageDimensionException If the image is not 2D.
		 * @return Radon transform image in square.
		 */
		EMData *do_radon();


		/** Calculate Cross-Correlation Function (CCF).
		 * Calculate the correlation of two 1-, 2-, or 3-dimensional
		 * images.  Note: this method internally just calls the
		 * correlation function from fundamentals.h.
		 *
		 * @param[in] with The image used to calculate the CCF. If 'with' is
		 * NULL, the autocorrelation function is computed instead.
		 * @param[in] fpflag Specify how periodicity (or normalization) should
		 * be handled.  See fundamentals.h  The default is "CIRCULANT".  for
		 * specific flags.
		 * @param center whether or not to center the image (bring bottom left corner to center)
		 * @return Real-space image.
		 * @exception ImageDimensionException if nx > 1 and nx < 2*radius + 1
		 * @ingroup CUDA_ENABLED
		 */
		EMData *calc_ccf(EMData * with = 0, fp_flag fpflag = CIRCULANT, bool center=false);

		/** Zero the pixels in the bottom left corner of the image
		 *  If radius is greater than 1, than circulant zeroing occurs
		 *  assuming that the center of operation starts in the bottom left
		 *  corner and proceed outwards to the NE and backwards in a circulant
		 *  fashion towards the SW.
		 *  Intended to zero the area corresponding to the middle of the image,
		 *  as generated by calc_ccf
		 * @param radius the radius to the zeroing operation
		 * @exception ImageDimensionException if nx > 1 and nx < 2*radius + 1
		 * @exception ImageDimensionException if ny > 1 and ny < 2*radius + 1
		 * @exception ImageDimensionException if nz > 1 and nz < 2*radius + 1
		 */
		void zero_corner_circulant(const int radius = 0);

		/** Calculate Cross-Correlation Function (CCF) in the x-direction
		 * and adds them up, result in 1D.
		 * WARNING: this routine will modify the 'this' and 'with' to contain
		 * 1D fft's without setting some flags. This is an optimization
		 * for rotational alignment.
		 * @ingroup CUDA_ENABLED
		 * @param with The image used to calculate CCF.
		 * @param y0 Starting position in x-direction.
		 * @param y1 Ending position in x-direction. '-1' means the
		 *        end of the row.
		 * @param nosum If true, returns an image y1-y0+1 pixels high.
		 * @see #calc_ccf()
		 * @exception NullPointerException If input image 'with' is NULL.
		 * @exception ImageFormatException If 'with' and 'this' are
		 * not same size.
		 * @exception ImageDimensionException If 'this' image is 3D.
		 * @return The result image containing the CCF.
		 */
		EMData *calc_ccfx( EMData * const with, int y0 = 0, int y1 = -1, bool nosum = false, bool flip = false);


		/** Makes a 'rotational footprint', which is an 'unwound'
		 * autocorrelation function. generally the image should be
		 * edge-normalized and masked before using this.
		 * @ingroup CUDA_ENABLED
		 * @param unwrap RFP undergoes polar->cartesian x-form
		 * @exception ImageFormatException If image size is not even.
		 * @return The rotaional footprint image.
		 */
		EMData *make_rotational_footprint(bool unwrap = true);
		EMData *make_rotational_footprint_e1(bool unwrap = true);
		EMData *make_rotational_footprint_cmc(bool unwrap = true);

		/** Makes a 'footprint' for the current image. This is image containing
		 * a rotational & translational invariant of the parent image. The size of the
		 * resulting image depends on the selected type.
		 *
		 * type 0- The original, default footprint derived from the rotational footprint
		 * types 1-6 - bispectrum-based
		 * types 1,3,5 - returns Fouier-like images
		 * types 2,4,6 - returns real-space-like images
		 * type 1,2 - simple r1,r2, 2-D footprints
		 * type 3,4 - r1,r2,anle 3D footprints
		 * type 5,6 - same as 1,2 but with the cube root of the final products used
		 *
		 * @param type Select one of several possible algorithms for producing the invariants
		 * @exception ImageFormatException If image size is not even.
		 * @return The footprint image.
		 */
		EMData *make_footprint(int type=0);


		/** Calculates mutual correlation function (MCF) between 2 images.
		 * If 'with' is NULL, this does mirror ACF.
		 *
		 * @param with The image used to calculate MCF.
		 * @param tocorner Set whether to translate the result image
		 *        to the corner.
		 * @param filter The filter image used in calculating MCF.
		 * @exception ImageFormatException If 'with' is not NULL and
		 * it doesn't have the same size to 'this' image.
		 * @exception NullPointerException If FFT returns NULL image.
		 * @return Mutual correlation function image.
		 */
		EMData *calc_mutual_correlation(EMData * with, bool tocorner = false, EMData * filter = 0);


		/** Maps to polar coordinates from Cartesian coordinates. Optionaly radially weighted.
		 * When used with RFP, this provides 1 pixel accuracy at 75% radius.
		 * 2D only.
		 * @ingroup CUDA_ENABLED
		 * @param r1 inner ring (all rings less than r1 are discarded) Default = 4
		 * @param r2 outer ring (all rings > r2 are discarded) Default = ny/2 - 2 - floor(hyp(dx,dy))
		 * @param xs Number of angular bins. Default = (2 if do360) PI * ny/4 - xs % 8
		 * @param dx origin offset in x
		 * @param dy origin offest in y
		 * @param do360  If true, do 0-360 degree mapping. Otherwise, do 0-180 degree mapping.
		 * @param weight_radial if true (default) reights the pixel value by its radius
		 * @exception ImageDimensionException If 'this' image is not 2D.
		 * @exception UnexpectedBehaviorException if the dimension of this image and the function arguments are incompatibale - i.e. the return image is less than 0 in some dimension.
		 * @return The image in Cartesian coordinates.
		 */
		EMData *unwrap(int r1 = -1, int r2 = -1, int xs = -1, int dx = 0,
							   int dy = 0, bool do360 = false, bool weight_radial=true) const;

		EMData * unwrap_largerR(int r1,int r2,int xs, float rmax_f);

		EMData *oneDfftPolar(int size, float rmax, float MAXR);

		
		/** multiplies by a radial function in fourier space.
		 *
		 * @param x0  starting point x coordinate.
		 * @param dx  step of x.
		 * @param array radial function data array.
		 * @param interp Do the interpolation or not.
		 */
		void apply_radial_func(float x0, float dx, vector < float >array, bool interp = true);


		/** calculates radial distribution. works for real and imaginary images.
		 * Returns mean radial amplitude, or intensity if inten is set. Note that the complex
		 * origin is at (0,0), with periodic boundaries. Note that the inten option is NOT
		 * equivalent to returning amplitude and squaring the result.
		 *
		 * @param n number of points.
		 * @param x0 starting point x coordinate.
		 * @param dx step of x.
		 * @param inten returns intensity (amp^2) rather than amplitude if set
		 * @return The radial distribution in an array.
		 */
		vector < float >calc_radial_dist(int n, float x0, float dx,bool inten);


		/** calculates radial distribution subdivided by angle. works for real and imaginary images.
		 * 2-D only. The first returns a single vector of n*nwedge points, with radius varying first.
		 * That is, the first n points represent the radial profile in the first wedge.
		 * @param n number of points.
		 * @param x0 starting x coordinate.
		 * @param dx step of x.
		 * @param nwedge int number of wedges to divide the circle into
		 * @param inten returns intensity (amp^2) rather than amplitude if set
		 * @exception ImageDimensionException If 'this' image is not 2D.
		 * @return nwedge radial distributions packed into a single vector<float>
		 */
		vector < float >calc_radial_dist(int n, float x0, float dx, int nwedge, bool inten);


		/** Replace the image its complex conjugate.
		 * @exception ImageFormatException Image must be complex (and RI)
		 */
		void cconj();


		/** Adds 'obj' to 'this' incoherently. 'obj' and 'this' should
		 * be same size. Both images should be complex images.
		 *
		 * @param obj The image added to 'this' image.
		 * @exception ImageFormatException If the 2 images are not
		 * same size; or if the 2 images are not complex images.
		 */
		void add_incoherent(EMData * obj);


		/** Calculates the histogram of 'this' image. The result is
		 * stored in float array 'hist'. If hist_min = hist_max, use
		 * image data min as hist_min; use image data max as hist_max.
		 *
		 * @param hist_size Histogram array's size.
		 * @param hist_min Minimum histogram value.
		 * @param hist_max Maximum histogram value.
		 * @param brt
		 * @param cont
		 * @return histogram array of this image.
		 */
		vector <float> calc_hist(int hist_size = 128, float hist_min = 0, float hist_max = 0, const float& brt = 0.0f, const float& cont = 1.0f);


		/** Caculates the azimuthal distributions.
		 * works for real or complex images, 2D only.
		 *
		 * @param n  Number of elements.
		 * @param a0 Starting angle.
		 * @param da Angle step.
		 * @param rmin Minimum radius.
		 * @param rmax  Maximum radius.
		 * @exception ImageDimensionException If image is 3D.
		 * @return Float array to store the data.
		 */
		vector<float> calc_az_dist(int n, float a0, float da, float rmin,
								   float rmax);

#if 0
		void calc_rcf(EMData * with, vector < float >&sum_array);
#endif
		/** Calculates the distance between 2 vectors. 'this' image is
		 * 1D, which contains a vector; 'second_img' may be nD. One of
		 * its row is used as the second vector. 'second_img' and
		 * 'this' must have the same x size.
		 *
		 * @param second_img The image used to caculate the distance.
		 * @param y_index Specifies which row in 'second_img' is used to do
		 * the caculation.
		 * @exception ImageDimensionException If 'this' image is not 1D.
		 * @exception ImageFormatException If the 2 images don't have
		 * same xsize.
		 * @return The distance between 2 vectors.
		 */
		float calc_dist(EMData * second_img, int y_index = 0) const;

		/** Calculates the cross correlation with local normalization
		 * between 2 images. This is a faster version of local correlation
		 * that make use of Fourier convolution and correlation.
		 * With is the template - the thing that you wish to find in the this image.
		 * It should not be unecessarily padded. The function uses the size of with
		 * to determine the extent of the local neighborhood used in the local
		 * normalization (for technical details, see calc_fast_sigma_image).
		 * Note that this function circularly masks the template at its radius
		 * so the calling function need not do this beforehand.
		 * Works in 1,2 and 3D.
		 *
		 * @param with The image used to calculate cross correlation (the template)
		 * @return the local normalized cross correlation image - the phase origin is at the corner of the image
		 * @author David Woolford
		 * @date April 2008
		 */
		EMData *calc_flcf(EMData * with);

		/** Calculates the local standard deviation (sigma) image using the given
		 * mask image. The mask image is typically much smaller than this image,
		 * and consists of ones, or is a small circle consisting of ones. The extent
		 * of the non zero neighborhood explicitly defines the range over which
		 * the local standard deviation is determined.
		 * Fourier convolution is used to do the math, ala Roseman (2003, Ultramicroscopy)
		 * However, Roseman was just working on methods Van Heel had presented earlier.
		 * The normalize flag causes the mask image to be processed so that it has a
		 * unit sum.
		 * Works in 1,2 and 3D
		 *
		 * @param mask the image that will be used to define the neighborhood for determine the local standard deviation
		 * @return the sigma image, the phase origin is at the corner (not the center)
		 * @exception ImageDimensionException if the dimensions of with do not match those of this
		 * @exception ImageDimensionException if any of the dimensions sizes of with exceed of this image's.
		 * @author David Woolford
		 * @date April 2008
		*/
		EMData *calc_fast_sigma_image( EMData* mask);

		/** Convolutes 2 data sets. The 2 images must be of the same size.
		 * @param with One data set. 'this' image is the other data set.
		 * @exception NullPointerException If FFT resturns NULL image.
		 * @return The result image.
		 */
		EMData *convolute(EMData * with);

#if 0
		void create_ctf_map(CtfMapType type, XYData * sf = 0);
#endif


		/** @ingroup tested2 */
		/** Finds common lines between 2 complex images.
		 *
		 * This function does not assume any symmetry, just blindly
		 * compute the "sinogram" and the user has to take care how
		 * to interpret the returned "sinogram". it only considers
		 * inplane rotation and assumes prefect centering and identical
		 * scale.
		 *
		 * @param image1 The first complex image.
		 * @param image2 The second complex image.
		 * @param mode Either 0 or 1 or 2. mode 0 is a summed
		 *   dot-product, larger value means better match; mode 1 is
		 *   weighted phase residual, lower value means better match.
		 * @param steps 1/2 of the resolution of the map.
		 * @param horizontal In horizontal way or not.
		 * @exception NullPointerException If 'image1' or 'image2' is NULL.
		 * @exception OutofRangeException If 'mode' is invalid.
		 * @exception ImageFormatException If 'image1' 'image2' are
		 * not same size.
		 */
		void common_lines(EMData * image1, EMData * image2, int mode = 0,
						  int steps = 180, bool horizontal = false);

		/** @ingroup tested2 */
		/** Finds common lines between 2 real images.
		 *
		 * @param image1 The first image.
		 * @param image2 The second image.
		 * @param steps 1/2 of the resolution of the map.
		 * @param horizontal In horizontal way or not.
		 * @exception NullPointerException If 'image1' or 'image2' is NULL.
		 * @exception ImageFormatException If 'image1' 'image2' are not same size.
		 */
		void common_lines_real(EMData * image1, EMData * image2,
							   int steps = 180, bool horizontal = false);

		/** cut a 2D slice out of a real 3D map. Put slice into 'this' image.
		 *
		 * @param map The real 3D map.
		 * @param tr orientation of the slice as encapsulated in a Transform object.
		 * @param interpolate Do interpolation or not.
		 * @exception NullPointerException If map is NULL.
		 * @exception ImageDimensionException If this image is not 2D.
		 * @exception ImageDimensionException If map image is not 3D.
		 * @exception ImageFormatException If this image is complex
		 * @exception ImageFormatException If map is complex
		 * @author David Woolford (adapted from an original version by Steve Ludtke)
		 * @date Feb 2008
		 */
		void cut_slice(const EMData * const map, const Transform& tr, bool interpolate = true);

		/** Opposite of the cut_slice(). It will take a slice and insert
		 * the data into a real 3D map. It does not interpolate, it uses
		 * the nearest neighbor.
		 *
		 * @param map  The real 3D map.
		 * @param tr Orientation of the slice.
		 * @exception NullPointerException If map is NULL.
		 * @exception ImageDimensionException If this image is not 2D.
		 * @exception ImageDimensionException If map image is not 3D.
		 * @exception ImageFormatException If this image is complex
		 * @exception ImageFormatException If map is complex
		 * @author David Woolford (adapted from an original version by Steve Ludtke)
		 * @date Feb 2008
		 */
		void uncut_slice(EMData * const map, const Transform& tr) const;

		/** Extract a box from EMData in an abritrary orrientation. Used for extracting helix boxes from tomograms
		
		EMData extract_box(const Transform& cs, const int ix, const int fx, const int iy, const int yf, const int zi, const int zf);
		* @param cs Transform describing the coordinate system of the box realative to the standard coord system
		* @param r Region describe the volume to extract, in the local coordinate system
		* @author John Flanagan
		* @date Aug 2011
		**/
		EMData *extract_box(const Transform& cs, const Region& r);
		
		/** function for MarchingCubes, for 3D image display
		 * @return the resolution
		 * */
		int getResolution() const {
			int resolution = 0;
			int num = 1;
			while(num < get_xsize()) {
				resolution++;
				num = 1 << resolution;
			}

			return resolution;
		}
		
		/** Printing EMData params for debugging purpose.
		 * */
		void debug_print_parms()
		{
			std::cout << "Printing EMData params" << std::endl;
			for ( Dict::const_iterator it = attr_dict.begin(); it != attr_dict.end(); ++it )
			{
				std::cout << (it->first) << " " << (it->second).to_str() << std::endl;
			}
			std::cout << "Done printing EMData params" << std::endl;
		}

		/**Set the x,y,z origin of the image
		 * @param origin_x  the x origin
		 * @param origin_y  the y origin
		 * @param origin_z  the z origin
		 * */
		void set_xyz_origin(float origin_x, float origin_y, float origin_z);
		
		/**Find the mean and variance of voxels in the missing wedge
		 * @param wedgeangle the angle of the missing wedge
		 * @param influnce the region of influnce in fourier space. This is a fudge factor between 0 and 0.5
		 * @param wedgedirection the direction of the wedge, so far only a wedge along Z is supported (set wedgedirection to 0)
		 * */
		EMData* compute_missingwedge_stats(float wedgeangle, float influnce = 0.05, int wedgedirection = 0);
		
		static int totalalloc;
	private:
		/** This EMDataFlags is deprecated. For anything which is currently handled by setting a
		 * bit in 'flags', instead, set or read an appropriately named attribute
		 * in the attributes dictionary. While there is a small overhead in the
		 * string lookup, none of these things should be called in the innermost
		 * loop anywhere, so it should be fine. --Grant
		 * */
		enum EMDataFlags {
//			EMDATA_COMPLEX = 1 << 1,
//			EMDATA_RI = 1 << 2,	       // real/imaginary or amp/phase
			EMDATA_BUSY = 1 << 3,	   // someone is modifying data, NO LONGER USED
			EMDATA_HASCTFF = 1 << 4,   // has CTF info in the image file
			EMDATA_NEEDUPD = 1 << 5,   // needs a real update
//			EMDATA_COMPLEXX = 1 << 6,  // 1D fft's in X
			EMDATA_FLIP = 1 << 7,	   // is the image flipped
			EMDATA_PAD = 1 << 8,       // is the image fft padded
			EMDATA_FFTODD = 1 << 9,	   // is the (real-space) nx odd
			EMDATA_SHUFFLE = 1 << 10,  // fft been shuffled? (so O is centered) PRB
			EMDATA_FH = 1 << 11,        // is the complex image a FH image
			EMDATA_CPU_NEEDS_UPDATE = 1 << 12, // CUDA related: is the CPU version of the image out out data
			EMDATA_GPU_NEEDS_UPDATE = 1 << 13, // CUDA related: is the GPU version of the image out out data
			EMDATA_GPU_RO_NEEDS_UPDATE = 1 << 14 // // CUDA related: is the GPU RO version of the image out out data
		};

		void update_stat() const;
		void save_byteorder_to_dict(ImageIO * imageio);

	private:
		/** to store all image header info */
		mutable Dict attr_dict;
		/** image real data */
		mutable float *rdata;
		/** supplementary data array */
		float *supp;

		/** CTF data
		 * All CTF data become attribute ctf(vector<float>) in attr_dict  --Grant Tang*/
		//Ctf *ctf;

		/** flags */
		mutable int flags;
		// Incremented every time the image changes
		int changecount;
		/** image size */
		int nx, ny, nz, nxy;
		size_t nxyz;
		/** array index offsets */
		int xoff, yoff, zoff;

		/** translation from the original location */
		Vec3f all_translation;
//		Vec3f all_rotation;    /** rotation (az, alt, phi) from the original locaton*/

		string path;
		int pathnum;

		/** This is a cached rotational footprint, can save much time */
		mutable EMData* rot_fp;

		// Clip inplace variables is a local class used from convenience in EMData::clip_inplace
		// Added by d.woolford
		class ClipInplaceVariables
		{
			public:
				ClipInplaceVariables(const int p_nx, const int p_ny, const int p_nz, const int n_nx, const int n_ny, const int n_nz,const int xtrans, const int ytrans, const int ztrans) :
					prv_nx(p_nx), prv_ny(p_ny), prv_nz(p_nz), new_nx(n_nx), new_ny(n_ny), new_nz(n_nz), xshift(xtrans), yshift(ytrans), zshift(ztrans),
				 x_iter(prv_nx), y_iter(prv_ny), z_iter(prv_nz), new_z_top(0), new_z_bottom(0),  new_y_back(0), new_y_front(0),new_x_left(0), new_x_right(0),
				prv_z_top(0), prv_z_bottom(0), prv_y_back(0), prv_y_front(0), prv_x_left(0), prv_x_right(0)
			{
				if ( xtrans > 0 ) x_iter -= xtrans;
				if ( x_iter < 0 ) x_iter = 0;
				if ( ytrans > 0 ) y_iter -= ytrans;
				if ( y_iter < 0 ) y_iter = 0;
				if ( ztrans > 0 ) z_iter -= ztrans;
				if ( z_iter < 0 ) z_iter = 0;

				// Get the depth in the new volume where slices are inserted
				// if this value is zero it means that the last z-slice in the new
				// volume contains image data
				if ( (new_nz + ztrans) > prv_nz ) new_z_top = new_nz + ztrans - prv_nz;
				if ( (new_ny + ytrans) > prv_ny ) new_y_back = new_ny + ytrans - prv_ny;
				if ( (new_nx + xtrans) > prv_nx ) new_x_right = new_nx + xtrans - prv_nx;

				if ( (new_nz + ztrans) < prv_nz )
				{
					prv_z_top = prv_nz - new_nz - ztrans;
					z_iter -= prv_z_top;
				}
				if ( (new_ny + ytrans) < prv_ny )
				{
					prv_y_back = prv_ny - new_ny - ytrans;
					y_iter -= prv_y_back;
				}
				if ( (new_nx + xtrans) < prv_nx )
				{
					prv_x_right = prv_nx - new_nx - xtrans;
					x_iter -= prv_x_right;
				}

				if ( xtrans > 0 ) prv_x_left = xtrans;
				if ( ytrans > 0 ) prv_y_front = ytrans;
				if ( ztrans > 0 ) prv_z_bottom = ztrans;

				if ( xtrans < 0 ) new_x_left = -xtrans;
				if ( ytrans < 0 ) new_y_front = -ytrans;
				if ( ztrans < 0 ) new_z_bottom = -ztrans;

			}
			~ClipInplaceVariables() {}

			int prv_nx, prv_ny, prv_nz, new_nx, new_ny, new_nz;
			int xshift, yshift, zshift;
			int x_iter, y_iter, z_iter;
			int new_z_top, new_z_bottom, new_y_back, new_y_front, new_x_left, new_x_right;
			int prv_z_top, prv_z_bottom,  prv_y_back, prv_y_front, prv_x_left, prv_x_right;
		};



		//		/**  Do the Fourier Harmonic Transform  PRB
		//		 * Takes a real image, returns the FH
		//		 * Sets the EMDATA_FH switch to indicate that it is an FH image
		//		 * @exception ImageFormatException If the image is not a square real odd image.
		//		 * @return the FH image.
		//		 */
		// 		EMData* do_FH();

		//		/**   Do the Inverse Fourier Harmonic Transform   PRB
		//		 * Takes an FH image, returns a square  complex image with odd sides
		//		 * @exception ImageFormatException If the image is not the FH of something
		//		 * @return a square complex image with odd sides
		//		 */
		// 		EMData* do_FH2F();



		//		/** Caclulates normalization and phase residual for a slice in
		//		 * an already existing volume. phase residual is calculated
		//		 * using only the inner 1/2 of the fourier sphere. Both the
		//		 * slice image and this image must be in complex image format.
		//		 *
		//		 * @param slice An slice image to be normalized.
		//		 * @param orient Orientation of the slice.
		//		 * @exception ImageFormatException If the images are not complex.
		//		 * @exception ImageDimensionException If the image is 3D.
		//		 * @return A float number pair (result, phase-residual).
		//		 */
		//		FloatPoint normalize_slice(EMData * slice, const Transform3D & orient);

		//		/** Caclulates normalization and phase residual for a slice in
		//		 * an already existing volume. phase residual is calculated
		//		 * using only the inner 1/2 of the fourier sphere. Both the
		//		 * slice image and this image must be in complex image format.
		//		 *
		//		 * @param slice An slice image to be normalized.
		//		 * @param alt Orientation euler angle alt (in EMAN convention).
		//		 * @param az  Orientation euler angle az  (in EMAN convention).
		//		 * @param phi Orientation euler angle phi (in EMAN convention).
		//		 * @exception ImageFormatException If the images are not complex.
		//		 * @exception ImageDimensionException If the image is 3D.
		//		 * @return A float number pair (result, phase-residual).
		//		 */
		//		FloatPoint normalize_slice(EMData * slice, float az, float alt, float phi);

		//		/** Get the normalization and phase residual values
		//		 * Used for normalizaton and error measurement when 2D slices are inserted into a 3D volume of Fourier pixels
		//		 * Originally added for use by the FourierReconstructor object
		//		 * @return the normalization const (pair.first) and the phase residual (pair.second)
		//		 * @param slice -the slice to be inserted into the 3D volume
		//		 * @param euler - the euler angle orientation of the slice
		//		 * @exception ImageDimensionException If this image is not 3D.ImageFormatException
		//		 * @exception ImageFormatException If this image is not complex
		//		 * @exception ImageFormatException If the slice not complex
		//		 */
		//		pair<float, float> get_normalization_and_phaseres( const EMData* const slice, const Transform3D& euler );

		//		/** cut a 2D slice out of a this 3D image and return it
		//		 * An alternative to cut_slice
		//		 * @param tr orientation of the slice as encapsulated in a Transform object
		//		 * @exception ImageDimensionException If this image is not 3D.
		//		 * @exception ImageFormatException If this image is complex
		//		 * @author David Woolford (adapted from an original version by Steve Ludtke)
		//		 * @date Feb 2009
		//		 */
		// 		EMData* get_cut_slice(const Transform& tr);

	};


	EMData * operator+(const EMData & em, float n);
	EMData * operator-(const EMData & em, float n);
	EMData * operator*(const EMData & em, float n);
	EMData * operator/(const EMData & em, float n);

	EMData * operator+(float n, const EMData & em);
	EMData * operator-(float n, const EMData & em);
	EMData * operator*(float n, const EMData & em);
	EMData * operator/(float n, const EMData & em);

	EMData * rsub(const EMData & em, float n);
	EMData * rdiv(const EMData & em, float n);

	EMData * operator+(const EMData & a, const EMData & b);
	EMData * operator-(const EMData & a, const EMData & b);
	EMData * operator*(const EMData & a, const EMData & b);
	EMData * operator/(const EMData & a, const EMData & b);


/*   Next  is Modified by PRB      Transform3D::EMAN,
	inline Transform3D EMData::get_transform() const
	{
		return Transform3D((float)attr_dict["euler_alt"],
				   (float)attr_dict["euler_az"],
				   (float)attr_dict["euler_phi"]);
	}
*/


}


#endif

/* vim: set ts=4 noet nospell: */
