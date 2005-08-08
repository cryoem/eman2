/**
 * $Id$
 */
#ifndef eman__emdata_h__
#define eman__emdata_h__ 1

#include <float.h>
//#include <boost/multi_array.hpp>
#include <complex>
#include <fstream>

#include "fundamentals.h"
#include "emutil.h"
#include "util.h"
//#include "vec3.h"
#include "geometry.h"
#include "transform.h"

using std::string;
using std::vector;
using std::map;
using std::complex;
using std::ostream;

namespace EMAN
{
	class ImageIO;
	class Ctf;
	class XYData;
	class Transform3D;
	
	typedef boost::multi_array_ref<float, 2> MArray2D;
	typedef boost::multi_array_ref<float, 3> MArray3D;
	typedef boost::multi_array_ref<complex<float>, 2> MCArray2D;
	typedef boost::multi_array_ref<complex<float>, 3> MCArray3D;
	typedef boost::multi_array<int, 2> MIArray2D;
	typedef boost::multi_array<int, 3> MIArray3D;
	
	/** EMData stores an image's data and defines core image processing routines.
     * The image is 1D, 2D or 3D, in real space or fourier space (complex image).
	 *
	 * Data are ordered with x increasing fastest, then y, then z.
	*/
	class EMData
	{
	static int totalalloc;
	public:
		
		enum FFTPLACE { FFT_OUT_OF_PLACE, FFT_IN_PLACE };
		enum WINDOWPLACE { WINDOW_OUT_OF_PLACE, WINDOW_IN_PLACE };
		
		/** Construct an empty EMData instance. It has no image data. */
		EMData();
		virtual ~ EMData();

		/** read an image file and stores its information to this
		 * EMData object.
		 *
		 * If a region is given, then only read a
		 * region of the image file. The region will be this
		 * EMData object. The given region must be inside the given
		 * image file. Otherwise, an error will be created.
		 *
		 * @param filename The image file name.
		 * @param img_index The nth image you want to read.
		 * @param header_only To read only the header or both header and data.
		 * @param region To read only a region of the image.
		 * @param is_3d  Whether to treat the image as a single 3D or a
		 *   set of 2Ds. This is a hint for certain image formats which
		 *   has no difference between 3D image and set of 2Ds.
		 * @exception ImageFormatException
		 * @exception ImageReadException
		 */
		void read_image(const string & filename, int img_index = 0,
						bool header_only = false,
						const Region * region = 0, bool is_3d = false);

		/** write the header and data out to an image.
		 *
		 * If the img_index = -1, append the image to the given image file.
		 *
		 * If the given image file already exists, this image
		 * format only stores 1 image, and no region is given, then
		 * truncate the image file  to  zero length before writing
		 * data out. For header writing only, no truncation happens.
		 *
		 * If a region is given, then write a region only.
		 *
		 * @param filename The image file name.
		 * @param img_index The nth image to write as.
		 * @param imgtype Write to the given image format type. if not
		 *        specified, use the 'filename' extension to decide.
		 * @param header_only To write only the header or both header and data.
		 * @param region Define the region to write to.
		 * @param filestoragetype The image data type used in the output file.
		 * @param use_host_endian To write in the host computer byte order.
		 * 
		 * @exception ImageFormatException
		 * @exception ImageWriteException
		 */
		void write_image(const string & filename,
						 int img_index = 0,
						 EMUtil::ImageType imgtype = EMUtil::IMAGE_UNKNOWN,
						 bool header_only = false,
						 const Region * region = 0,
						 EMUtil::EMDataType filestoragetype = EMUtil::EM_FLOAT,
						 bool use_host_endian = true);

		/** append to an image file; If the file doesn't exist, create one.
		 *
		 * @param filename The image file name.
		 * @param imgtype Write to the given image format type. if not
		 *        specified, use the 'filename' extension to decide.
		 * @param header_only To write only the header or both header and data.
		 */
		void append_image(const string & filename,
						  EMUtil::ImageType imgtype = EMUtil::IMAGE_UNKNOWN,
						  bool header_only = false);

		/** Append data to a LST image file.
		 * @param filename The LST image file name.
		 * @param reffile Reference file name.
		 * @param refn The reference file number.
		 * @param comment The comment to the added reference file.
		 * @see lstio.h
		 */
		void write_lst(const string & filename, 
					   const string & reffile="", int refn=-1,
					   const string & comment="");
		   
		
		/** Print the image data to a file stream (standard out by default).
		 * @param out Output stream; cout by default.
		 * @param str Message string to be printed.
		 */
		void print_image(const string str = string(""), 
				ostream& out = std::cout);

		/** Apply a processor with its parameters on this image.
		 * @param processorname Processor Name.
		 * @param params Processor parameters in a keyed dictionary.
		 * @exception NotExistingObjectError If the processor doesn't exist.
		 */
		void process(const string & processorname, const Dict & params = Dict());

		/** Compare this image with another image.
		 * @param cmpname Comparison algorithm name.
		 * @param with The image you want to compare to.
		 * @param params Comparison parameters in a keyed dictionary.
		 * @exception NotExistingObjectError If the comparison algorithm doesn't exist.
		 * @return comparison score. The bigger, the better.
		 */
		float cmp(const string & cmpname, EMData * with, const Dict & params);

		/** Align this image with another image and return the result image.
		 *
		 * @param aligner_name Alignment algorithm name.
		 * @param to_img The image 'this' image aligns to.
		 * @param params  Alignment algorithm parameters in a keyed dictionary.
		 * @param comp_name Comparison algorithm used in alignment.
		 * @param cmp_params Parameter dictionary for comparison algorithm.
		 * @exception NotExistingObjectError If the alignment algorithm doesn't exist.
		 * @return The result image.
		 */
		EMData *align(const string & aligner_name, EMData * to_img,
					  const Dict & params = Dict(), const string & comp_name = "", 
					  const Dict& cmp_params = Dict());

		/** Calculate the projection of this image and return the result.
		 * @param projector_name Projection algorithm name.
		 * @param params Projection Algorithm parameters.
		 * @exception NotExistingObjectError If the projection algorithm doesn't exist.
		 * @return The result image.
		 */
		EMData *project(const string & projector_name, const Dict & params = Dict());

		/** Make a copy of this image including both data and header.
		 * @return A copy of this image including both data and header.
		 */
		EMData *copy() const;
		
		/** Make an image with a copy of the current image's header.
		 * @return An image with a copy of the current image's header.
		 */
		EMData *copy_head() const;

		/** Get an inclusive clip. Pads 0 if larger than this image.
		 * area can be 2D/3D.
		 * @param area The clip area.
		 * @return The clip image.
		 */
		EMData *get_clip(const Region & area);

		/** Insert a clip into this image.
		 * @param block An image block.
		 * @param origin The origin location to insert the clip.
		 * @exception ImageFormatException If clip is outside of the
		 * destination image (i.e., this image).
		 */
		void insert_clip(EMData * block, const IntPoint & origin);

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
		EMData *get_rotated_clip(const Transform3D & xform, const IntSize &size, float scale=1.0);
				
		/** Add a scaled image into another image at a specified location.
		 *  This is used, for example, to accumulate gaussians in
		 *  programs like pdb2mrc.py. The center of 'block' will be positioned at
		 *  'center' with scale factor 'scale'. Densities will be interpolated in
		 *  'block' and multiplied by 'mult'.
		 *
		 * @param block The image to inserted.
		 * @param center The center of the inserted block in 'this'.
		 * @param scale  Scale factor.
		 * @param mult_factor Number used to multiply the block's densities.
		 * @exception ImageDimensionException If 'this' image is not 2D/3D.
		 */
		void insert_scaled_sum(EMData *block, const FloatPoint & center,
							   float scale=1.0, float mult_factor=1.0);


		/** During image reconstruction the image may have been
		 *  padded with zeros for fourier interpolation.  In that
		 *  case the desired lxlxl image lies in the center of a larger
		 *  volume.  This routine creates a new
		 *  object that contains only the desired lxlxl volume.
		 *  (This routine is a thin wrapper around get_clip.)
		 *
		 * @return An image object that has been windowed.
		 */
		EMData* window_padded(int l);


		
		/** Multiply a real image by (-1)**(ix+iy+iz) to center
		 *  the fft version.
		 *
		 */
		void center_origin();

		/** Multiply a Fourier image by (-1)**(ix+iy+iz) to center it.
		 *
		 */
		void center_origin_fft();

		/** return an image object that has been padded with zeros to
		 *  an integral (npad) factor (e.g., npad=4 means the new 
		 *  image will be 2x larger in each direction).  The default
		 *  is to pad 4x.
		 * The current image is not changed.
		 *
		 * @return An image object that has been padded npad-times.
		 */
		EMData *zeropad_ntimes(int npad=4);

		/** return an image object that has been fft-padded/unpadded.
		 * The current image is not changed.
		 *
		 * @return An image object that has been fft-padded/unpadded.
		 */
		EMData *pad_fft(int npad = 1);

		/** Remove padding, leaving a single corner of the image.
		 *  The current image is changed in place.
		 *
		 *  The assumption is that after an in-place inverse fft
		 *  the real-space image contains too much information 
		 *  because it may have been zero-padded some integer factor
		 *  of times and it has also been extended slightly along x
		 *  for the fft.  Here we keep only the data corresponding 
		 *  to ix=0,...,nxold-1, iy=0,...,nyold-1, iz=0,...,nzold-1,
		 *  where nxold, nyold, nzold are the sizes of the original
		 *  image.
		 *
		 * @return Pointer to the depadded image.
		 */
		void postift_depad_corner_inplace();



		/** returns the fourier harmonic transform (FH) image of the current
		 * image (in real space). The current image is not changed. The result is in
		 * real/imaginary format. The FH switch is set on.
		 *
		 */
		EMData *real2FH(float OverSamplekB);

		/** returns the fourier version of the image 
		 * from the FH version. The current image is not changed. The result is in
		 * real/imaginary format. The FH switch is set off.
		 *
		 */
		EMData *FH2F(int Size, float OverSamplekB);

		/** return the fast fourier transform (FFT) image of the current
		 * image. the current image is not changed. The result is in
		 * real/imaginary format.
		 *
		 */
		EMData *do_fft();

		/** Do FFT inplace. And return the FFT image.
		 * @return The FFT of the current image in real/imaginary format.
		 */
		EMData* do_fft_inplace();

		/** return the inverse fourier transform (IFT) image of the current
		 * image. the current image is not changed.
		 *
		 * @exception ImageFormatException If the image is not a complex image.
		 * @return The current image's inverse fourier transform image.
		 */
		EMData *do_ift();

		/* Do IFT inplace. And return the IFT image.
		 * @return The IFT image.
		 */
		EMData* do_ift_inplace();
		
		/**  Do the Fourier Harmonic Transform  PRB
		 * Takes a real image, returns the FH
		 * Sets the EMDATA_FH switch to indicate that it is an FH image
		 * @exception ImageFormatException If the image is not a square real odd image.
		 * @return the FH image.
		 */
//		EMData* do_FH();
		
		/**   Do the Inverse Fourier Harmonic Transform   PRB
		 * Takes an FH image, returns a square  complex image with odd sides
		 * @exception ImageFormatException If the image is not the FH of something
		 * @return a square complex image with odd sides
		 */
//		EMData* do_FH2F();

		/** return the amplitudes of the FFT including the left half
		 *
		 * @exception ImageFormatException If the image is not a complex image.
		 * @return The current FFT image's amplitude image.
		 */
		EMData *get_fft_amplitude();

		/** return the amplitudes of the 2D FFT including the left half
		 *     PRB
		 * @exception ImageFormatException If the image is not a complex image.
		 * @return The current FFT image's amplitude image.
		 */
		EMData *get_fft_amplitude2D();

		/** return the phases of the FFT including the left half
		 *
		 * @exception ImageFormatException If the image is not a complex image.
		 * @return The current FFT image's phase image.
		 */
		EMData *get_fft_phase();

		/** Caclulates normalization and phase residual for a slice in
		 * an already existing volume. phase residual is calculated
		 * using only the inner 1/2 of the fourier sphere. Both the
		 * slice image and this image must be in complex image format.
		 *
		 * @param slice An slice image to be normalized.
		 * @param orient Orientation of the slice.
		 * @exception ImageFormatException If the images are not complex.
		 * @exception ImageDimensionException If the image is 3D.
		 * @return A float number pair (result, phase-residual).
		 */
//		FloatPoint normalize_slice(EMData * slice, const Transform3D & orient);

		/** Caclulates normalization and phase residual for a slice in
		 * an already existing volume. phase residual is calculated
		 * using only the inner 1/2 of the fourier sphere. Both the
		 * slice image and this image must be in complex image format.
		 *
		 * @param slice An slice image to be normalized.
		 * @param alt Orientation euler angle alt (in EMAN convention).
		 * @param az  Orientation euler angle az  (in EMAN convention).
		 * @param phi Orientation euler angle phi (in EMAN convention).
		 * @exception ImageFormatException If the images are not complex.
		 * @exception ImageDimensionException If the image is 3D.
		 * @return A float number pair (result, phase-residual).
		 */
//		FloatPoint normalize_slice(EMData * slice, float az, float alt, float phi);

		/** Render the image into an 8-bit image. 2D image only.
		 *
		 * @param x	origin of the area to render
		 * @param y
		 * @param xsize	size of the area to render in output pixels
		 * @param ysize
		 * @param bpl	bytes per line, if asrgb remember *3
		 * @param scale	scale factor for rendering
		 * @param min_gray	minimum gray value to render (0-255)
		 * @param max_gray	maximum gray value to render (0-255)
		 * @param min_render	float image density corresponding to min_gray
		 * @param max_render	float image density corresponding to max_gray
		 * @param asrgb	duplicate each output pixel 3x for RGB rendering
		 * @exception ImageDimensionException If the image is not 2D.
		 */
		std::string render_amp8(int x, int y, int xsize, int ysize,
						 int bpl, float scale, int min_gray, int max_gray,
						 float min_render, float max_render,int asrgb);
		
		/** Render the image into a 24-bit image. 2D image only.
		 * @param x
		 * @param y
		 * @param xsize
		 * @param ysize
		 * @param bpl
		 * @param scale
		 * @param min_gray
		 * @param max_gray
		 * @param min_render
		 * @param max_render
		 * @param ref
		 * @param cmap
		 * @exception ImageDimensionException If the image is not 2D.
		 */		 
		void render_amp24(int x, int y, int xsize, int ysize,
						  int bpl, float scale, int min_gray, int max_gray,
						  float min_render, float max_render,
						  void *ref, void cmap(void *, int coord, unsigned char *tri));

		/** convert the complex image from real/imaginary to amplitude/phase */
		void ri2ap();

		/** convert the complex image from amplitude/phase to real/imaginary */
		void ap2ri();

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
		 * @param t Transformation rotation.
		 */
		void rotate(const Transform3D & t);

		/** Rotate this image.
		 * @param az  Rotation euler angle az  in EMAN convention.
		 * @param alt Rotation euler angle alt in EMAN convention.		 
		 * @param phi Rotation euler angle phi in EMAN convention.
		 */
		void rotate(float az, float alt, float phi);

		/** Rotate then translate the image.
		 * @param t The rotation and translation transformation to be done.
		 */
		void rotate_translate(const Transform3D & t);

		/** Rotate then translate the image.
		 * @param alt Rotation euler angle alt in EMAN convention.
		 * @param az  Rotation euler angle az  in EMAN convention.
		 * @param phi Rotation euler angle phi in EMAN convention.
		 * @param dx Translation distance in x direction.
		 * @param dy Translation distance in y direction.
		 * @param dz Translation distance in z direction.
		 */
		void rotate_translate(float az, float alt, float phi, float dx, float dy, float dz);
		
		/** Rotate then translate the image.
		 * @param alt Rotation euler angle alt in EMAN convention.
		 * @param az  Rotation euler angle az  in EMAN convention.
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
		void rotate_180();

		/** dot product of 2 images. Then 'this' image is rotated/translated.
		 * It is much faster than Rotate/Translate then dot product. 
		 * 2D images only.
		 *
		 * @param with The image used to do the dot product.
		 * @param dx Translation distance in x direction.
		 * @param dy Translation distance in y direction.
		 * @param da Rotation euler angle.
		 * @exception ImageFormatException If the 2 images are not the
		 * same size.
		 * @exception ImageDimensionException If the image is 3D.
		 * @return
		 */
		double dot_rotate_translate(EMData * with, float dx, float dy, float da);

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
		 *
		 * Calculate the correlation of two 1-, 2-, or 3-dimensional
		 * images.  Note: this method internally just calls the
		 * correlation function from fundamentals.h.
		 *
		 * @param[in] with The image used to calculate the CCF. If 'with' is
		 * NULL, the autocorrelation function is computed instead.
		 * @param[in] fpflag Specify how periodicity (or normalization) should
		 * be handled.  See fundamentals.h  The default is "CIRCULANT".  for
		 * specific flags.
		 * @return Real-space image.
		 */
		EMData *calc_ccf(EMData * with, fp_flag fpflag = CIRCULANT);

		/** Calculate Cross-Correlation Function (CCF) in the x-direction 
		 * and adds them up, result in 1D.
		 * WARNING: this routine will modify the 'this' and 'with' to contain
		 * 1D fft's without setting some flags. This is an optimization
		 * for rotational alignment.
		 *
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
		EMData *calc_ccfx(EMData * with, int y0 = 0, int y1 = -1, bool nosum = false);

		/** Makes a 'rotational footprint', which is an 'unwound'
		 * autocorrelation function. generally the image should be
		 * edge-normalized and masked before using this.
		 *
		 * @param unwrap To cache the rfp or not. false means not cached.
		 * @param premasked Is the image pre-masked?
		 * @exception ImageFormatException If image size is not even.
		 * @return The rotaional footprint image.
		 */
		EMData *make_rotational_footprint(bool premasked = false, bool unwrap = true);

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

		/** maps polar coordinates to Cartesian coordinates. radially weighted.
		 * When used with RFP, this provides 1 pixel accuracy at 75% radius.
		 * 2D only.
		 *
		 * @param r1
		 * @param r2
		 * @param xs
		 * @param dx
		 * @param dy
		 * @param do360  If true, do 0-360 degree mapping. Otherwise,
		 * do 0-180 degree mapping.
		 * @exception ImageDimensionException If 'this' image is not 2D.
		 * @return The image in Cartesian coordinates.
		 */		 
		EMData *unwrap(int r1 = -1, int r2 = -1, int xs = -1, int dx = 0,
					   int dy = 0, bool do360 = false);

		/** Reduces the size of the image by a factor 
		 * using the average value of the pixels in a block.
		 * @param shrink_factor Image shrink factor.
		 * @exception InvalidValueException If shrink factor is invalid.
		 */
		void mean_shrink(float shrink_factor);
		
		/* Reduces the size of the image by a factor using a local median processor.
		 *
		 * @param shrink_factor Image shrink factor.
		 * @exception InvalidValueException If shrink factor is invalid.
		 */
		void median_shrink(int shrink_factor);

		/** multiplies by a radial function in fourier space.
		 *
		 * @param x0  starting point x coordinate.
		 * @param dx  step of x.
		 * @param array radial function data array.
		 * @param interp Do the interpolation or not.
		 */
		void apply_radial_func(float x0, float dx, vector < float >array, bool interp = true);

		/** calculates radial distribution. works for real and imaginary images. 
		 * 
		 * @param n number of points.
		 * @param x0 starting point x coordinate.
		 * @param dx step of x.
		 * @return The radial distribution in an array.
		 */					
		vector < float >calc_radial_dist(int n, float x0, float dx);

		/** calculates radial distribution. works for real and imaginary images. 
		 * @param n number of points.
		 * @param x0 starting x coordinate.
		 * @param dx step of x.
		 * @param acen The direction.
		 * @param arange The angular range around the direction in radians.
		 * @exception ImageDimensionException If 'this' image is not 2D.
		 * @return The radial distribution in an array.
		 */
		vector < float >calc_radial_dist(int n, float x0, float dx, float acen, float arange);

		/** Create a (1-D) rotationally averaged image.
		 * 
		 * @return 1-D rotationally-averaged image
		 */					
		EMData* rotavg();

		/** add a number to each pixel value of the image. Image may be real or complex.
		 * @param f The number added to 'this' image.
		 * @param keepzero If set will not modify pixels that are exactly zero
		 */
		void add(float f,int keepzero=0);
		
		/** add a same-size image to this image pixel by pixel.
		 *
		 * @param image The image added to 'this' image.
		 * @exception ImageFormatException If the 2 images are not same size.
		 */
		void add(const EMData & image);
		
		/** subtract a number to each pixel value of the image.
		 * @param f The number subtracted from 'this' image.		 
		 */
		void sub(float f);
		
		/** subtract a same-size image from this image pixel by pixel.
		 * @param image The image subtracted  from 'this' image.
		 * @exception ImageFormatException If the 2 images are not same size.
		 */
		void sub(const EMData & image);

		/** multiply a number to each pixel value of the image.
		 * @param n The number multiplied to 'this' image.
		 */
		void mult(int n)
		{
			mult((float)n);
		}

		/** multiply a number to each pixel value of the image.
		 * @param f The number multiplied to 'this' image.
		 */
		void mult(float f);
		
		/** multiply each pixel of this image with each pixel of some
		 * other same-size image.
		 *
		 * @param image The image multiplied to 'this' image.
		 * @exception ImageFormatException If the 2 images are not same size.
		 */
		void mult(const EMData & image);

		/** make each pixel value divided by a number.
		 * @param f The number 'this' image divided by.
		 */
		void div(float f);

		/** make each pixel value divided by pixel value of another
		 * same-size image.
		 * @param image The image 'this' image divided by.
		 * @exception ImageFormatException If the 2 images are not same size.
		 */
		void div(const EMData & image);

		/** Get the image pixel density data in a 1D float array.
		 * @return The image pixel density data.
		 */
		float *get_data() const;
		
		/** Done with data manipulation. It marks EMData as changed.
		 *
		 * This function is used together with 'get_data()'.
		 * A typical case is 1) call
		 * get_data(); 2) work on the data. 3) if data is changed, then
		 * call done_data. If not changed, no need to call done_data.
		 */
		void done_data();

		/** Mark EMData as changed, statistics, etc will be updated at need.*/
		void update();

		/** Make all the pixel value = 0. */
		void to_zero();
		
		/** Make all the pixel value = 1. */
		void to_one();


		/** Adds 'obj' to 'this' incoherently. 'obj' and 'this' should
		 * be same size. Both images should be complex images.
		 *
		 * @param obj The image added to 'this' image.
		 * @exception ImageFormatException If the 2 images are not
		 * same size; or if the 2 images are not complex images.
		 */
		void add_incoherent(EMData * obj);

		/** Caculates the fourier ring/shell correlation coefficients
		 * as an array with ysize/2 elements (corners not calculated).
		 * The input image 'with' must have the same size to 'this' image.
		 *
		 * @param with The image used to caculate the fourier shell
		 * correlation together with 'this' image.
		 * @exception ImageFormatException If the 2 images are not
		 * same size.
		 * @return The fourier shell correlation coefficients array.
		 */
		vector < float >calc_fourier_shell_correlation(EMData * with);

		/** Calculates the histogram of 'this' image. The result is
		 * stored in float array 'hist'. If 'add' is true, the new
		 * histogram data will be added to existing 'hist' array.
		 * If hist_min = hist_max, use image data min as hist_min; use
		 * image data max as hist_max.
		 *
		 * @param hist Float array storing histogram data.
		 * @param hist_min Minimum histogram value.
		 * @param hist_max Maximum histogram value.
		 * @param add If true, the new histogram data will be added to
		 * existing 'hist' array.
		 */
		void calc_hist(vector < float >&hist, float hist_min = 0, float hist_max = 0,
					   bool add = false);

		/** Caculates the azimuthal distributions.
		 * works for real or complex images, 2D only.
		 *
		 * @param n  Number of elements.
		 * @param a0 Starting angle.
		 * @param da Angle step.
		 * @param data Float array to store the data.
		 * @param rmin Minimum radius.
		 * @param rmax  Maximum radius.
		 * @exception ImageDimensionException If image is 3D.
		 */
		void calc_az_dist(int n, float a0, float da, float *data,
						  float rmin, float rmax);
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
		 * between 2 images. This is a fater version of local correlation.
		 *
		 * This is the fast local correlation program based upon local
		 * real space correlation. The traditional real-space technies
		 * are know to be senisitve for finding small objects in a large
		 * field due to local optimization of numerical
		 * scaling. However, they are slow to compute. This technique is
		 * based upon Fourier transform and claimed to be two times
		 * faster than the explicit real-space formula.
		 * 
		 * The technique is published in Ultramicroscopy by Alan M. Roseman, 
		 * MRC Cambridge.
		 *
		 * Imlemented by htet khant, 02-2003
		 *
		 * @param with The image used to calculate cross correlation.
		 * @param radius
		 * @param maskfilter
		 * @return the cross correlation image.
		 */
		EMData *calc_flcf(EMData * with, int radius = 50,
						  const string & maskfilter = "mask.sharp");

		/** Convolutes 2 data sets. The 2 images must be of the same size.
		 * @param with One data set. 'this' image is the other data set.
		 * @exception NullPointerException If FFT resturns NULL image.
		 * @return The result image.
		 */
		EMData *convolute(EMData * with);

		/** check whether the image physical file has the CTF info or not.
		 * @return True if it has the CTF information. Otherwise, false.
		*/
		bool has_ctff() const;

#if 0
		void create_ctf_map(CtfMapType type, XYData * sf = 0);
#endif

		/** Dot product 2 images. The 2 images must be of same size.
		 * If 'evenonly' is true, only calculates pixels with even
		 * positions assuming all pixels are in a single array. If
		 * 'evenonly' is false, calculates all pixels. Shortcut for
		 * cmp("dot")
		 *
		 * @param with The image to do dot product with.
		 * @return The dot product result.
		 */
		float dot(EMData * with);

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
		 * @param steps: 1/2 of the resolution of the map.
		 * @param horizontal In horizontal way or not.
		 * @exception NullPointerException If 'image1' or 'image2' is NULL.
		 * @exception OutofRangeException If 'mode' is invalid.
		 * @exception ImageFormatException If 'image1' 'image2' are
		 * not same size.
		 */
		void common_lines(EMData * image1, EMData * image2, int mode = 0,
						  int steps = 180, bool horizontal = false);

		/** Finds common lines between 2 real images.
		 *
		 * @param image1 The first image.
		 * @param image2 The second image.
		 * @param steps: 1/2 of the resolution of the map.
		 * @param horizontal In horizontal way or not.
		 * @exception NullPointerException If 'image1' or 'image2' is NULL.
		 * @exception ImageFormatException If 'image1' 'image2' are
		 * not same size.
		 */
		void common_lines_real(EMData * image1, EMData * image2,
							   int steps = 180, bool horizontal = false);

		/** cut a 2D slice out of a real 3D map. Put slice into 'this' image.
		 *
		 * @param map The real 3D map.
		 * @param dz 
		 * @param orientation Orientation of the slice.
		 * @param interpolate Do interpolation or not.
		 * @param dx
		 * @param dy
		 */
		void cut_slice(const EMData * map, float dz, Transform3D * orientation = 0,
					   bool interpolate = true, float dx = 0, float dy = 0);

		/** Opposite of the cut_slice(). It will take a slice and insert
		 * the data into a real 3D map. It does not interpolate, it uses
		 * the nearest neighbor.
		 *
		 * @param map  The real 3D map.
		 * @param dz
		 * @param orientation Orientation of the slice.
		 * @param dx
		 * @param dy
		 */		 
		void uncut_slice(EMData * map, float dz, Transform3D * orientation = 0,
						 float dx = 0, float dy = 0);

		/** Calculates the density value at the peak of the
		 * image histogram, sort of like the mode of the density.
		 * @return The density value at the peak of the image histogram.
		*/
		float calc_center_density();

		/** Calculates sigma above and below the mean and returns the
		 * difference between them.
		 * @return The difference between sigma above and below the mean.
		 */
		float calc_sigma_diff();

		/** Calculates the coordinates of the minimum-value pixel.
		 * @return The coordinates of the minimum-value pixel.
		 */
		IntPoint calc_min_location() const;
		
		/** Calculates the coordinates of the maximum-value pixel.
		 * @return The coordinates of the maximum-value pixel.
		 */
		IntPoint calc_max_location() const;

		/** Calculates the index of minimum-value pixel when assuming
		 * all pixels are in a 1D array.
		 * @return Index of the minimum-value pixel.
		 */
		int calc_min_index() const;

		/** Calculates the index of maximum-value pixel when assuming
		 * all pixels are in a 1D array.
		 * @return Index of the maximum-value pixel.
		 */
		int calc_max_index() const;

		/** Calculate and return a sorted list of pixels whose values
		 * are above a specified threshold. The pixels are sorted 
		 * from high to low.
		 *
		 * @param threshold The specified pixel value. Returned pixels
		 *        should have higher values than it.
		 * @return A sorted list of pixels with their values, and
		 *     locations. Their values are higher than threshold.
		 */
		vector<Pixel> calc_highest_locations(float threshold);

		/** Calculates the mean pixel values around the (1 pixel) edge
		 * of the image.
		 *
		 * @return The mean pixel values around the (1 pixel) edge.
		*/
		float get_edge_mean() const;

		/** Calculates the circular edge mean by applying a circular
		 * mask on 'this' image.
		 * @return The circular edge mean.
		 */
		float get_circle_mean();

		/** Get ctf parameter of this image.
		 * @return The ctf parameter.
		 */
		Ctf *get_ctf() const;

		/** Set the CTF parameter of this image.
		 * @param ctf The CTF parameter object.
		 */
		void set_ctf(Ctf * ctf);

		/** Get 'this' image's translation vector from the original
		 * location.
		 * @return 'this' image's translation vector from the original
		 * location.
		 */
		Vec3f get_translation() const;

		/** Set 'this' images' translation vector from the original
		 * location.
		 * @param new_translation The new translation vector.
		 */
		void set_translation(const Vec3f &new_translation);

		/** Set 'this' images' translation vector from the original
		 * location.
		 * @param dx The translation distance in x direction.
		 * @param dy The translation distance in y direction.
		 * @param dz The translation distance in z direction.
		 */
		void set_translation(float dx, float dy, float dz);

		/** Get the 3D orientation of 'this' image.
		 * @return The 3D orientation of 'this' image.
		 */
		Transform3D get_transform() const; 

		/** Define the 3D orientation of this particle, also
		 * used to indicate relative rotations for reconstructions
		 *
		 * @param alt 'alt' Euler angle in EMAN convention.
		 * @param az  'az' Euler angle in EMAN convention.
		 * @param phi 'phi' Euler angle in EMAN convention.
		 */
		void set_rotation(float az, float alt, float phi);

		/** Resize 'this' image.
		 *
		 * @param nx  x size of this image.
		 * @param ny  y size of this image.
		 * @param nz  z size of this image.
		 */
		void set_size(int nx, int ny=1, int nz=1);

		/** Resize 'this' complex image.
		 *
		 * @param nx  x size of this image.
		 * @param ny  y size of this image.
		 * @param nz  z size of this image.
		 */
		void set_complex_size(int nx, int ny=1, int nz=1) {
			set_size(nx*2, ny, nz); 
		}

		/** Set the path
		 * @param path The new path.
		 */
		void set_path(const string & path);

		/** Set the number of paths.
		 * @param n The number of paths.
		 */
		void set_pathnum(int n);

		/** Get image raw pixel data in a 2D multi-array format. The
		 * array shares the memory space with the image data.
		 * Notice: the subscription order is d[y][x]
		 *
		 * It should be used on 2D image only.
		 *
		 * @return 2D multi-array format of the raw data.
		 */
		MArray2D get_2dview() const;
		
		/** Get image raw pixel data in a 3D multi-array format. The
		 * array shares the memory space with the image data.
		 * Notice: the subscription order is d[z][y][x]
		 *
		 * It should be used on 3D image only.
		 *
		 * @return 3D multi-array format of the raw data.
		 */
		MArray3D get_3dview() const;
		
		/** Get complex image raw pixel data in a 2D multi-array format. 
		 * The array shares the memory space with the image data.
		 *
		 * It should be used on 2D image only.
		 *
		 * @return 2D multi-array format of the raw data.
		 */
		MCArray2D get_2dcview() const;
		
		/** Get complex image raw pixel data in a 3D multi-array format. 
		 * The array shares the memory space with the image data.
		 *
		 * It should be used on 3D image only.
		 *
		 * @return 3D multi-array format of the raw data.
		 */
		MCArray3D get_3dcview() const;
		
		/** Get pointer to a complex image raw pixel data in a 3D multi-array format. 
		 * The array shares the memory space with the image data.
		 *
		 * It should be used on 3D image only.
		 *
		 * @return Pointer to a 3D multi-array format of the raw data.
		 */
		MCArray3D* get_3dcviewptr() const;
		
		/** Get image raw pixel data in a 2D multi-array format. The
		 * data coordinates is translated by (x0,y0) such that
		 * array[y0][x0] points to the pixel at the origin location.
		 * the data coordiates translated by (x0,y0). The
		 * array shares the memory space with the image data.
		 *
		 * It should be used on 2D image only.
		 *
		 * @param x0 X-axis translation amount.
		 * @param y0 Y-axis translation amount.
		 * @return 2D multi-array format of the raw data.
		 */
		MArray2D get_2dview(int x0, int y0) const;

		/** Get image raw pixel data in a 3D multi-array format. The
		 * data coordinates is translated by (x0,y0,z0) such that
		 * array[z0][y0][x0] points to the pixel at the origin location.
		 * the data coordiates translated by (x0,y0,z0). The
		 * array shares the memory space with the image data.
		 *
		 * It should be used on 3D image only.
		 *
		 * @param x0 X-axis translation amount.
		 * @param y0 Y-axis translation amount.
		 * @param z0 Z-axis translation amount.
		 * @return 3D multi-array format of the raw data.
		 */
		MArray3D get_3dview(int x0, int y0, int z0) const;

		/** Get complex image raw pixel data in a 2D multi-array format. The
		 * data coordinates is translated by (x0,y0) such that
		 * array[y0][x0] points to the pixel at the origin location.
		 * the data coordiates translated by (x0,y0). The
		 * array shares the memory space with the image data.
		 *
		 * It should be used on 2D image only.
		 *
		 * @param x0 X-axis translation amount.
		 * @param y0 Y-axis translation amount.
		 * @return 2D multi-array format of the raw data.
		 */
		MCArray2D get_2dcview(int x0, int y0) const;

		/** Get complex image raw pixel data in a 3D multi-array format. The
		 * data coordinates is translated by (x0,y0,z0) such that
		 * array[z0][y0][x0] points to the pixel at the origin location.
		 * the data coordiates translated by (x0,y0,z0). The
		 * array shares the memory space with the image data.
		 *
		 * It should be used on 3D image only.
		 *
		 * @param x0 X-axis translation amount.
		 * @param y0 Y-axis translation amount.
		 * @param z0 Z-axis translation amount.
		 * @return 3D multi-array format of the raw data.
		 */
		MCArray3D get_3dcview(int x0, int y0, int z0) const;
		
		/** Get one row of a 1D/2D image.
		 *
		 * @param row_index Index of the row.
		 * @exception ImageDimensionException If this image is 3D.
		 * @return A 1D image with the row data.
		 */
		EMData *get_row(int row_index) const;

		/** Set one row of a 1D/2D image.
		 *
		 * @param data The row image data.
		 * @param row_index Index of the row.
		 * @exception ImageDimensionException If this image is 3D.
		 */
		void set_row(const EMData * data, int row_index);

		/** Get one column of a 2D images.
		 *
		 * @param col_index Index of the column.
		 * @exception ImageDimensionException If this image is not 2D.
		 * @return A 1D image with the column data.
		 */
		EMData *get_col(int col_index) const;

		/** Set one column of a 2D image.
		 *
		 * @param data The column image data.
		 * @param col_index Index of the column.
		 * @exception ImageDimensionException If this image is not 2D.
		 */
		void set_col(const EMData * data, int col_index);

		/** The generic way to get any image header information
		 * given a header attribute name.
		 *
		 * @param attr_name The header attribute name.
		 * @return The attribute value.
		 */
		EMObject get_attr(const string & attr_name);
		
		/** Set a header attribute's value.
		 *
		 * @param attr_name The header attribute name.
		 * @param val The attribute value.
		 */
		void set_attr(const string & attr_name, EMObject val);

		/** Get the image attribute dictionary containing all the
		 * image attribute names and attribute values.
		 *
		 * @return The image attribute dictionary containing all
		 * attribute names and values.
		 */
		Dict get_attr_dict();

		/** Set the attribute dictionary to a new dictioanary.
		 *
		 * @param new_dict The new attribute dictionary.
		 */
		void set_attr_dict(const Dict & new_dict);

		/** Get the image x-dimensional size.
		 * @return Image x-dimensional size.
		 */
		int get_xsize() const;
		
		/** Get the image y-dimensional size.
		 * @return Image y-dimensional size.
		 */
		int get_ysize() const;
		
		/** Get the image z-dimensional size.
		 * @return Image z-dimensional size.
		 */
		int get_zsize() const;

		/** Get image dimension.
		 * @return image dimension.
		 */
		int get_ndim() const;

		/** Get the pixel density value at coordinates (x,y,z).
		 * The validity of x, y, and z is not checked.
		 *
		 * @param x The x cooridinate.
		 * @param y The y cooridinate.
		 * @param z The z cooridinate.
		 * @return The pixel density value at coordinates (x,y,z).
		 */
		float get_value_at(int x, int y, int z) const;

		
		/** Get the pixel density value at coordinates (x,y). 2D only.
		 * The validity of x, y is not checked.
		 *
		 * @param x The x cooridinate.
		 * @param y The y cooridinate.
		 * @return The pixel density value at coordinates (x,y).
		 */
		float get_value_at(int x, int y) const;

		/** Get the pixel density value given an index 'i' assuming
		 * the pixles are stored in a 1D array. The validity of i
		 * is not checked.
		 * 
		 * @param i  1D data array index.
		 * @return The pixel density value
		 */
		float get_value_at(size_t i) const;

		/** A safer, slower way to get the pixel density value at
		 * coordinates (x,y,z). The validity of x, y, and z is checked.
		 * If the coordinates are out of range, return 0;
		 *
		 * @param x The x cooridinate.
		 * @param y The y cooridinate.
		 * @param z The z cooridinate.
		 * @return The pixel density value at coordinates (x,y,z).
		 */
		float sget_value_at(int x, int y, int z) const;

		/** A safer, slower way to get the pixel density value at
		 * coordinates (x,y). 2D only. The validity of x, y is checked.
		 * If the coordinates are out of range, return 0;
		 *
		 * @param x The x cooridinate.
		 * @param y The y cooridinate.
		 * @return The pixel density value at coordinates (x,y).
		 */
		float sget_value_at(int x, int y) const;

		/** A safer, slower way to get the pixel density value
		 * given an index 'i' assuming
		 * the pixles are stored in a 1D array. The validity of i
		 * is checked. If i is out of range, return 0;
		 * 
		 * @param i  1D data array index.
		 * @return The pixel density value
		 */
		float sget_value_at(size_t i) const;

		/** Get pixel density value at interpolation of (x,y).
		 * The validity of x, y is checked.2D image only. 
		 *
		 * @param x The x cooridinate.
		 * @param y The y cooridinate.
		 * @return The pixel density value at coordinates (x,y).
		 */
		float sget_value_at_interp(float x, float y) const;
		
		/** Get the pixel density value at interpolation of (x,y,z).
		 * The validity of x, y, and z is checked.
		 *
		 * @param x The x cooridinate.
		 * @param y The y cooridinate.
		 * @param z The z cooridinate.
		 * @return The pixel density value at coordinates (x,y,z).
		 */
		float sget_value_at_interp(float x, float y, float z) const;

		/** Set the pixel density value at coordinates (x,y,z).
		 * The validity of x, y, and z is not checked.
		 * This implementation does bounds checking.
		 * 
		 * @param x The x cooridinate.
		 * @param y The y cooridinate.
		 * @param z The z cooridinate.
		 * @param v The pixel density value at coordinates (x,y,z).
		 */
		void set_value_at(int x, int y, int z, float v);
		
		/** Set the pixel density value at coordinates (x,y,z).
		 * The validity of x, y, and z is not checked.
		 * This implementation has no bounds checking.
		 * 
		 * @param x The x cooridinate.
		 * @param y The y cooridinate.
		 * @param z The z cooridinate.
		 * @param v The pixel density value at coordinates (x,y,z).
		 */
		void set_value_at_fast(int x, int y, int z, float v);
		
		/** Set the pixel density value at coordinates (x,y).
		 * 2D image only. The validity of x, y, and z is not checked.
		 *
		 * @param x The x cooridinate.
		 * @param y The y cooridinate.
		 * @param v The pixel density value at coordinates (x,y).
		 */
		void set_value_at(int x, int y, float v);
		
		/** Set the pixel density value at coordinates (x,y).
		 * 2D image only. The validity of x, y, and z is not checked.
		 *
		 * @param x The x cooridinate.
		 * @param y The y cooridinate.
		 * @param v The pixel density value at coordinates (x,y).
		 */
		void set_value_at_fast(int x, int y, float v);

		/** Has this image been shuffled?
		 * @return Whether this image has been shuffled to put origin in the center.
		 */
		bool is_shuffle() const;
		
		/** Is this a FH image?
		 * @return Whether this is a FH image or not.
		 */
		bool is_FH() const;
		
		/** Is this a complex image?
		 * @return Whether this is a complex image or not.
		 */
		bool is_complex() const;

		/** Is this a real image?
		 * @return Whether this is image is real (not complex) or not.
		 */
		bool is_real() const;

		/** Mark this image as a shuffled image.
		 * @param is_shuffle If true, a shuffled image. If false, not 
		 *          a shuffled image.
		 */
		void set_shuffle(bool is_shuffle);

		/** Mark this complex image as a FH image.
		 * @param is_FH If true, a FH image. If false, 
		 *        not a FH image.
		 */
		void set_FH(bool is_FH);

		/** Mark this image as a complex image.
		 * @param is_complex If true, a complex image. If false, a real
		 * image.
		 */
		void set_complex(bool is_complex);

		/** Is this image a 1D FFT image in X direction?
		 * @return Whether this image is a 1D FFT image in X
		 * direction.
		 */
		bool is_complex_x() const;

		/** Marks this image a 1D FFT image in X direction.
		 * @param is_complex_x If true, a 1D FFT image in X direction;
		 * If false, not such an image.
		 */
		void set_complex_x(bool is_complex_x);

		/** Is this image flipped?
		 * @return Whether this image is flipped or not.
		 */
		bool is_flipped() const;

		/** Mark this image as flipped.
		 * @param is_flipped If true, mark this image as flipped;
		 * If false, mark this image as not flipped.
		 */
		void set_flipped(bool is_flipped);

		/** Is this image a real/imaginary format complex image?
		 * @return Whether this image is real/imaginary format complex
		 * image.
		 */
		bool is_ri() const;

		/** Mark this image as a real/imaginary format complex image.
		 * @param is_ri If true, mark as real/imaginary format; If
		 * false, mark as amp/phase format.
		 */
		void set_ri(bool is_ri);

		/** Is this image already extended along x for ffts?
		 * @return Whether this image is extended along x for ffts.
		 */
		bool is_fftpadded() const;

		/** Mark this image as already extended along x for ffts.
		 * @param is_padded If true, mark as padded along x; If
		 * false, mark as not padded along x.
		 */
		void set_fftpad(bool is_padded);

		/** Does this image correspond to a (real-space) odd nx?
		 * @return Whether this image has a (real-space) odd nx.
		 */
		bool is_fftodd() const;

		/** Mark this image as having (real-space) odd nx.
		 * @param is_fftodd If true, mark as nx odd; If
		 * false, mark as nx not odd.
		 */
		void set_fftodd(bool is_fftodd);

		/** Set the number of complex elements along x.
		 * @param nxc is the number of complex elements along x.
		 */
		void set_nxc(int nxc);

		EMData & operator+=(float n);
		EMData & operator-=(float n);
		EMData & operator*=(float n);
		EMData & operator/=(float n);

		EMData & operator+=(const EMData & em);
		EMData & operator-=(const EMData & em);
		EMData & operator*=(const EMData & em);
		EMData & operator/=(const EMData & em);
		
		/** return a image to the power of n
		 * @param n
		 */
		EMData & power(int n);		
		
		/** return real and imaginary part of a complex image as a real image format
		 */
		EMData & real();
		EMData & imag();
		
		/** create a complex image from a real image, this complex image is in real/imaginary format
		 * @param img give an artificial imaginary part
		 */
		EMData & real2complex(float img = 0.0f);
		
		/** Read a set of images from file specified by 'filename'.
		 * Which images are read is set by 'img_indices'.
		 * @param filename The image file name.
		 * @param img_indices Which images are read. If it is empty,
		 *     all images are read. If it is not empty, only those 
		 *     in this array are read.
		 * @param header_only If true, only read image header. If
		 *     false, read both data and header.
		 * @return The set of images read from filename.
		 */
		static vector < EMData * >read_images(const string & filename,
											  vector < int >img_indices = vector < int >(),
											  bool header_only = false);

		/** Read a set of images from file specified by 'filename'. If
		 * the given 'ext' is not empty, replace 'filename's extension it.
		 * Images with index from img_index_start to img_index_end are read.
		 *
		 * @param filename The image file name.
		 * @param img_index_start Starting image index.
		 * @param img_index_end Ending image index.
		 * @param header_only If true, only read image header. If
		 *     false, read both data and header.
		 * @param ext The new image filename extension.
		 * @return The set of images read from filename.
		 */
		static vector < EMData * >read_images_ext(const string & filename,
												  int img_index_start,
												  int img_index_end,
												  bool header_only = false,
												  const string & ext = "");


		/** Helper function for method nn.
		 *
		 * @param j y fourier index (frequency)
		 * @param n number of real elements.
		 * @param n2 Number of complex elements.
		 * @param x  Complex matrix of [0:n2][1:n][1:n]
		 * @param nr Normalization matrix [0:n2][1:n][1:n]
		 * @param bi Fourier transform matrix [0:n2][1:n]
		 * @param tf Transform3D reference
		 * @return The set of images read from filename.
		 */
		static void onelinenn(int j, int n, int n2, MCArray3D& x,
				              MIArray3D& nr, MCArray2D& bi, 
							  const Transform3D& tf);

		/** Nearest Neighbor interpolation.
		 *  Modifies the current object.
		 *
		 * @param norm Normalization data.
		 * @param myfft FFT data.
		 * @param tf Transform3D reference
		 */
		void nn(MIArray3D& norm, EMData* myfft, const Transform3D& tf);

		/** Symmetrize plane 0
		 *  Modifies the current object.
		 *
		 * @param norm Normalization data.
		 */
		void symplane0(MIArray3D& norm);

		/** Symmetrize volume in real space.
		 *  
		 *  @param[in] symmetry Point group of the target volume.
		 *  
		 *  @return New symmetrized volume object.
		 */
		EMData* symvol(string symmetry);

		/** Generate Rotated-Circulantly-Translated image 
		 *
		 *  @param[in] ang Rotation angle in degrees.
		 *  @param[in] delx Translation along x
		 *  @param[in] dely Translation along y
		 *  
		 *  @return New rotated/translated/scaled image
		 */
		EMData* 
		rot_trans2D(float ang, float delx=0.f, float dely=0.f);
		/** Generate Rotated-Scaled-Circulantly-Translated image 
		 *  (or image slice).
		 *
		 *  If the image is a volume, then only the specified slice
		 *  is rotated/translated/scaled.
		 *  
		 *  @param[in] ang Rotation angle in degrees.
		 *  @param[in] scale Scaling factor
		 *  @param[in] delx Amount to translate rotation origin along x
		 *  @param[in] dely Amount to translate rotation origin along y
		 *  
		 *  @return New rotated/translated/scaled image
		 */
		EMData* 
		rot_scale_trans2D(float ang, float scale, float delx, 
			 			float dely, int zslice = 1);
		/** Value of 2-D analytic masking (or 2-D convolution) at off-grid point.
		 *  
		 *  The only requirement for the window function object is that
		 *  it overload operator()(const float) and return a float.
		 *
		 *  @param[in] x x-value of the desired (potentially off-grid) point
		 *  @param[in] y y-value of the desired (potentially off-grid) point
		 *  @param[in] win Window (mask/kernel) function object.
		 *  @param[in] size Size of real-space kernel/mask.
		 *
		 *  @return Value of masked/convolved image at (x,y)
		 */
		//template<class Win>
		//float getconvpt2d(float x, float y, Win win, int size = 7);
		float getconvpt2d_kbi0(float x, float y, 
				Util::KaiserBessel::kbi0_win win, int size = 7);
		/** 2-D rotation using gridding convolution.
		 *  
		 *  The only requirement for the window function object is that
		 *  it overload operator()(const float) and return a float.
		 *
		 *  This routine does _not_ deconvolve out the window function
		 *  after rotation.
		 *
		 *  @param[in] x x-value of the desired (potentially off-grid) point
		 *  @param[in] y y-value of the desired (potentially off-grid) point
		 *  @param[in] win Window (mask/kernel) function object.
		 *  @param[in] size Size of real-space kernel/mask.
		 *
		 *  @return Rotated/convolved EMData image.
		 */
		//template<class Win>
		//EMData* rotconvtrunc2d(float ang, Win win, int size = 7);
		//EMData* rotconvtrunc2d_kbi0(float ang, float alpha, int size) {
		//	Util::KaiserBessel kb(alpha, size-1);
		//	return rotconvtrunc2d(ang, kb.get_kbi0_win(), size);
		//}
		EMData* rotconvtrunc2d_kbi0(float ang, float alpha, int size);


	private:
		enum EMDataFlags {
			EMDATA_COMPLEX = 1 << 1,
			EMDATA_RI = 1 << 2,	       // real/imaginary or amp/phase
			EMDATA_BUSY = 1 << 3,	   // someone is modifying data
			EMDATA_HASCTFF = 1 << 4,   // has CTF info in the image file
			EMDATA_NEEDUPD = 1 << 5,   // needs a real update			
			EMDATA_COMPLEXX = 1 << 6,  // 1D fft's in X
			EMDATA_FLIP = 1 << 7,	   // is the image flipped
			EMDATA_PAD = 1 << 8,       // is the image fft padded 
			EMDATA_FFTODD = 1 << 9,	   // is the (real-space) nx odd
			EMDATA_SHUFFLE = 1 << 10,  // fft been shuffled? (so O is centered) PRB
			EMDATA_FH = 1 << 11        // is the complex image a FH image
		};

		void update_stat();
		void set_xyz_origin(float origin_x, float origin_y, float origin_z);
		void scale_pixel(float scale_factor) const;
		void save_byteorder_to_dict(ImageIO * imageio);
		
	private:
		/** to store all image header info */
		mutable Dict attr_dict; 
		/** image real data */
		float *rdata;	 
		/** supplementary data array */       
		float *supp;    
		/** CTF data */        
		Ctf *ctf;		   
		/** rotational foot print */     
		EMData *rfp;          
		/** flags */  
		int flags;       
		/** image size */       
		int nx, ny, nz;	        

		/** translation from the original location */
		Vec3f all_translation; 
//		Vec3f all_rotation;    /** rotation (az, alt, phi) from the original locaton*/

		string path;
		int pathnum;
	};


	EMData * operator+(const EMData & em, float n);
	EMData * operator-(const EMData & em, float n);
	EMData * operator*(const EMData & em, float n);
	EMData * operator/(const EMData & em, float n);

	EMData * operator+(float n, const EMData & em);
	EMData * operator-(float n, const EMData & em);
	EMData * operator*(float n, const EMData & em);
	EMData * operator/(float n, const EMData & em);

	EMData * operator+(const EMData & a, const EMData & b);
	EMData * operator-(const EMData & a, const EMData & b);
	EMData * operator*(const EMData & a, const EMData & b);
	EMData * operator/(const EMData & a, const EMData & b);

	inline int EMData::get_xsize() const
	{
		return nx;
	}

	inline int EMData::get_ysize() const
	{
		return ny;
	}

	inline int EMData::get_zsize() const
	{
		return nz;
	}

	inline int EMData::get_ndim() const
	{
		if (nz <= 1) {
			if (ny <= 1) {
				return 1;
			}
			else {
				return 2;
			}
		}

		return 3;
	}


	inline float EMData::get_value_at(int x, int y, int z) const
	{
		return rdata[x + y * nx + z * nx * ny];
	}


	inline float EMData::get_value_at(int x, int y) const
	{
		return rdata[x + y * nx];
	}

	inline float EMData::get_value_at(size_t i) const
	{
		return rdata[i];
	}


	inline void EMData::set_value_at(int x, int y, int z, float v)
	{
		if( x>=nx && x<0 )
		{
			throw OutofRangeException(0, nx-1, x, "x dimension index");
		}
		else if( y>=ny && y<0 )
		{
			throw OutofRangeException(0, ny-1, y, "y dimension index");
		}
		else if( z>=nz && z<0 )
		{
			throw OutofRangeException(0, nz-1, z, "z dimension index");
		}
		else
		{
			rdata[x + y * nx + z * nx * ny] = v;
			flags |= EMDATA_NEEDUPD;
		}
	}

	inline void EMData::set_value_at_fast(int x, int y, int z, float v)
	{
		rdata[x + y * nx + z * nx * ny] = v;
		flags |= EMDATA_NEEDUPD;
	}
	
	inline void EMData::set_value_at(int x, int y, float v)
	{
		if( x>=nx && x<0 )
		{
			throw OutofRangeException(0, nx-1, x, "x dimension index");
		}
		else if( y>=ny && y<0 )
		{
			throw OutofRangeException(0, ny-1, y, "y dimension index");
		}
		else
		{
			rdata[x + y * nx] = v;
			flags |= EMDATA_NEEDUPD;
		}
	}

	inline void EMData::set_value_at_fast(int x, int y, float v)
	{
		rdata[x + y * nx] = v;
		flags |= EMDATA_NEEDUPD;
	}

	inline void EMData::update()
	{
		flags |= EMDATA_NEEDUPD;
	}


	inline bool EMData::is_shuffle() const
	{  //     PRB
		if (flags & EMDATA_SHUFFLE) {
			return true;
		}
		else {
			return false;
		}
	}

	inline bool EMData::is_FH() const
	{  //     PRB
		if (flags & EMDATA_FH) {
			return true;
		}
		else {
			return false;
		}
	}


	inline bool EMData::is_complex() const
	{
		if (flags & EMDATA_COMPLEX) {
			return true;
		}
		else {
			return false;
		}
	}

	inline bool EMData::is_real() const
	{
		return !is_complex();
	}

	inline bool EMData::is_complex_x() const
	{
		if (flags & EMDATA_COMPLEXX) {
			return true;
		}
		else {
			return false;
		}
	}
	inline bool EMData::is_ri() const
	{
		if (flags & EMDATA_RI) {
			return true;	
		}
		else {
			return false;
		}
	}
	inline bool EMData::is_fftpadded() const
	{
		if (flags & EMDATA_PAD) {
			return true;	
		}
		else {
			return false;
		}
	}
	inline bool EMData::is_fftodd() const
	{
		if (flags & EMDATA_FFTODD) {
			return true;	
		}
		else {
			return false;
		}
	}


	inline void EMData::set_nxc(int nxc)
	{
		attr_dict["nxc"] = nxc;
	}

	inline void EMData::set_shuffle(bool is_shuffle)
	{ // PRB
		if (is_shuffle) {
//			printf("entered correct part of set_shuffle \n");
			flags |=  EMDATA_SHUFFLE;
		}
		else {
			flags &= ~EMDATA_SHUFFLE;
		}
	}

	inline void EMData::set_FH(bool is_FH)
	{ // PRB
		if (is_FH) {
			flags |=  EMDATA_FH;
		}
		else {
			flags &= ~EMDATA_FH;
		}
	}

	inline void EMData::set_complex(bool is_complex)
	{
		if (is_complex) {
			flags |= EMDATA_COMPLEX;
		}
		else {
			flags &= ~EMDATA_COMPLEX;
		}
	}

	inline void EMData::set_complex_x(bool is_complex_x)
	{
		if (is_complex_x) {
			flags |= EMDATA_COMPLEXX;
		}
		else {
			flags &= ~EMDATA_COMPLEXX;
		}
	}


	inline void EMData::set_ri(bool is_ri)
	{
		if (is_ri) {
			flags |= EMDATA_RI;
		}
		else {
			flags &= ~EMDATA_RI;
		}
	}

	inline void EMData::set_fftpad(bool is_fftpadded)
	{
		if (is_fftpadded) {
			flags |= EMDATA_PAD;
		}
		else {
			flags &= ~EMDATA_PAD;
		}
	}

	inline void EMData::set_fftodd(bool is_fftodd)
	{
		if (is_fftodd) {
			flags |= EMDATA_FFTODD;
		}
		else {
			flags &= ~EMDATA_FFTODD;
		}
	}

	inline bool EMData::is_flipped() const
	{
		if (flags & EMDATA_FLIP) {
			return true;
		}
		else {
			return false;
		}
	}

	inline void EMData::set_flipped(bool is_flipped)
	{
		if (is_flipped) {
			flags |= EMDATA_FLIP;
		}
		else {
			flags &= ~EMDATA_FLIP;
		}
	}

	inline void EMData::set_path(const string & new_path)
	{
		path = new_path;
	}

	inline void EMData::set_pathnum(int n)
	{
		pathnum = n;
	}

	inline bool EMData::has_ctff() const
	{
		if ((flags & EMDATA_HASCTFF)) {
			return true;
		}
		
		return false;
	}

	inline Ctf *EMData::get_ctf() const
	{
		return ctf;
	}

	inline Vec3f EMData::get_translation() const
	{
		return all_translation;
	}

	inline void EMData::set_translation(const Vec3f &t)
	{
		all_translation = t;
	}

	inline void EMData::set_translation(float dx, float dy, float dz)
	{
		all_translation = Vec3f(dx, dy, dz);
	}

	inline Transform3D EMData::get_transform() const
	{
		return Transform3D( (float)attr_dict["euler_alt"],
				    (float)attr_dict["euler_az"],
				    (float)attr_dict["euler_phi"]);
	}
/*   Next  is Modified by PRB      Transform3D::EMAN,
	inline Transform3D EMData::get_transform() const 
	{
		return Transform3D((float)attr_dict["euler_alt"],
				   (float)attr_dict["euler_az"],
				   (float)attr_dict["euler_phi"]);
	}
*/
	inline void EMData::set_rotation(float az, float alt, float phi)
	{
        attr_dict["orientation_convention"] = "EMAN";
		attr_dict["euler_alt"]=alt;
		attr_dict["euler_az"]=az;
		attr_dict["euler_phi"]=phi;
    }
	
	inline void EMData::scale_pixel(float scale) const
	{
		attr_dict["apix_x"] = ((float) attr_dict["apix_x"]) * scale;
		attr_dict["apix_y"] = ((float) attr_dict["apix_y"]) * scale;
		attr_dict["apix_z"] = ((float) attr_dict["apix_z"]) * scale;
	}

	inline void EMData::set_attr(const string & key, EMObject val)
	{
		attr_dict[key] = val;
	}
	
	
}

			
#endif

/* vim: set ts=4 noet: */
