/**
 * $Id$
 */
#ifndef eman__emdata_h__
#define eman__emdata_h__ 1

#include <string>
#include <vector>
#include <map>
#include <math.h>
#include <float.h>

#include "emobject.h"
#include "emutil.h"
#include "util.h"
#include "transform.h"

using std::string;
using std::vector;
using std::map;


namespace EMAN
{
	class ImageIO;
	class Ctf;
	class XYData;

	/** EMData stores an image's data and defines core image processing routines.
     * The image is 2D or 3D, in real space or fourier space (complex image).
     */
	class EMData
	{
	public:
		EMData();
		virtual ~ EMData();

		/** read an image file and stores its information.
		 *
		 * @param filename The image file name.
		 * @param img_index The nth image you want to read.
		 * @param header_only To read only the header or both header and data.
		 * @param region To read only a region of the image.
		 * @param is_3d  Whether to treat the image as a single 3D or a
		 *   set of 2Ds. This is a hint for certain image formats which
		 *   has no difference between 3D image and set of 2Ds.
		 */
		void read_image(string filename, int img_index = 0, bool header_only = false,
						const Region * region = 0, bool is_3d = false);

		/** write the header and data out to an image.
		 *
		 * @param filename The image file name.
		 * @param img_index The nth image to write as.
		 * @param imgtype Write to the given image format type. if not
		 *        specified, use the 'filename' extension to decide.
		 * @param header_only To write only the header or both header and data.
		 * @param use_host_endian To write in the host computer byte order.
		 */
		void write_image(string filename, int img_index = 0,
						 EMUtil::ImageType imgtype = EMUtil::IMAGE_UNKNOWN, bool header_only = false,
						 bool use_host_endian = true);

		/** append to an image file; If the file doesn't exist, create one.
		 *
		 * @param filename The image file name.
		 * @param imgtype Write to the given image format type. if not
		 *        specified, use the 'filename' extension to decide.
		 * @param header_only To write only the header or both header and data.
		 */
		void append_image(string filename, EMUtil::ImageType imgtype = EMUtil::IMAGE_UNKNOWN,
						  bool header_only = false);

		/** Apply a filter with its parameters on this image.
		 * @param filtername Filter Name.
		 * @param params Filter parameters in a keyed dictionary.
		 * @exception NotExistingObjectError If the filter doesn't exist.
		 */
		void filter(string filtername, const Dict & params = Dict());

		/** Compare this image with another image.
		 * @param cmpname Comparison algorithm name.
		 * @param params Comparison parameters in a keyed dictionary.
		 * @exception NotExistingObjectError If the comparison algorithm doesn't exist.
		 * @return comparison score. The bigger, the better.
		 */
		float cmp(string cmpname, const Dict & params = Dict());

		/** Align this image with another image and return the result image.
		 *
		 * @param aligner_name Alignment algorithm name.
		 * @param params  Alignment algorithm parameters in a keyed dictionary.
		 * @param comp_name Comparison algorithm used in alighment.
		 * @exception NotExistingObjectError If the alignment algorithm doesn't exist.
		 * @return The result image.
		 */
		EMData *align(string aligner_name, const Dict & params = Dict(), string comp_name = "");

		/** Calculate the projection of this image and return the result.
		 * @param projector_name Projection algorithm name.
		 * @param params Projection Algorithm parameters.
		 * @exception NotExistingObjectError If the projection algorithm doesn't exist.
		 * @return The result image.
		 */
		EMData *project(string projector_name, const Dict & params = Dict());

		/** Return a copy of this image including both data and header*/
		EMData *copy(bool withparent = true);
		
		/** Return a image with only a copy of the header */
		EMData *copy_head();

		/** Inclusive clip. Pads 0 if larger than this image. */
		EMData *get_clip(const Region & area);
		void insert_clip(EMData * block, const Point < int >&originn);
		EMData *get_top_half() const;
		
		/** This will exctract an arbitrarily oriented and sized region from the
		 *  image. The orientation is defined by rotating the target image into
		 *  the source image ('this'). 'center' defines the location of the
		 *  center of the returned image in 'this' */ 
		EMData *get_rotated_clip(Point <float>&center, Rotation &orient, Size &size, float scale=1.0);
				
		/** Add a scaled image into another image at a specified location
		 *  This is used, for example, to accumulate gaussians in
		 *  programs like pdb2mrc.py. The center of 'block' will be positioned at
		 *  'center' with scale factor 'scale. Densities will be interpolated in
		 *  'block' and multiplied by 'mult' */
		void insert_scaled_sum(EMData *block, const Point <float>&center, float scale=1.0, float mult=1.0);

		/** return the fast fourier transform image of the current
		 * image. the current image is not changed.
		 * the result is in real/imaginary format.
		 */
		EMData *do_fft();

		/** return the inverse fourier transform image of the current
		 * image. the current image is not changed. 
		 */
		EMData *do_ift();

		Point < float >normalize_slice(EMData * slice, float alt, float az, float phi);

		void render_amp8(unsigned char *data, int x, int y, int xsize, int ysize,
						 int bpl, float scale, int min_gray, int max_gray,
						 float min_render, float max_render);
		void render_amp8_wrapper(long data, int x, int y, int xsize, int ysize,
								 int bpl, float scale, int min_gray, int max_gray,
								 float min_render, float max_render);
		void render_amp24(unsigned char *data, int x, int y, int xsize, int ysize,
						  int bpl, float scale, int min_gray, int max_gray,
						  float min_render, float max_render,
						  void *ref, void cmap(void *, int coord, unsigned char *tri));

		/** convert the complex image from real/imaginary to amplitude/phase */
		void ri2ap();
		/** convert the complex image from amplitude/phase to real/imaginary */
		void ap2ri();

		float *setup4slice(bool redo = false);

		void scale(float s);

		void translate(float dx, float dy, float dz);
		void translate(const Vec3 < float >&translation);

		void rotate(float alt, float az, float phi);
		void rotate(const Rotation & r);

		void rotate_translate(float alt, float az, float phi, float dx, float dy, float dz);
		void rotate_translate(const Rotation & rotation, const Vec3 < float >&translation);
		void rotate_translate(const Transform & xform);

		void rotate_x(int dx);
		void rotate_180();

		double dot_rotate_translate(EMData * data, float dx, float dy, float da);

		EMData *little_big_dot(EMData * little_img, bool do_sigma = false);

		/** Radon Transform: an algorithm that transforms an original
		 * image into a series of equiangular projections. When
		 * applied to a 2D object, the output of the Radon transform is a
		 * series of 1D lines.
		 *
		 * Do radon transformation on this image. This image must be
		 * 2D square.
		 * @return Radon transform image in square.
		 */
		EMData *do_radon();

		/** Calculate Cross-Correlation Function (CCF).
		 *
		 * CCF is a 2D function that is obtained by forming the scalar
		 * cross-product of two images (i.e., the sum of products of
		 * equivalent pixels) as a function of a 2D shift vector. The
		 * CCF is often used to achieve alignment between two images,
		 * since it displays a high value (a peak) at the place where
		 * a motif contained in both images come into register.
		 *
		 * @param with The image used to calculate CCF. If 'with' is
		 *   NULL, does mirror ACF.
		 * @param tocorner Set whether to translate the result image
		 *   to the corner.
		 * @param filter The filter image used in calculating CCF.
		 * @return The result image containing the CCF.
		 */
		EMData *calc_ccf(EMData * with, bool tocorner = false, EMData * filter = 0);

		/** Calculate Cross-Correlation Function (CCF) in the
		 * x-direction and adds them up, result in 1D.
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
		 */
		EMData *calc_ccfx(EMData * with, int y0 = 0, int y1 = -1, bool nosum = false);

		/** Makes a 'rotational footprint', which is an 'unwound'
		 * autocorrelation function. generally the image should be
		 * edgenormalized and masked before using this.
		 *
		 * @param unwrap To cache the rfp or not. false means not cached.
		 * @param premasked
		 */
		EMData *make_rotational_footprint(bool premasked = false, bool unwrap = true);

		/** Calculate mutual correlation function (MCF) between 2 images.
		 *
		 * @param with The image used to calculate MCF.
		 * @param tocorner Set whether to translate the result image
		 *        to the corner.
		 * @param filter The filter image used in calculating MCF.
		 */
		EMData *calc_mutual_correlation(EMData * with, bool tocorner = false, EMData * filter = 0);

		/** maps polar coordinates to Cartesian coordinates */
		EMData *unwrap(int r1 = -1, int r2 = -1, int xs = -1, int dx = 0,
					   int dy = 0, bool do360 = false);

		/** Reduces the size of the image by a factor of 'shrink_factor'
		 * using the average value of the pixels in a block.
		 */
		void mean_shrink(int shrink_factor);
		
		/* Reduces the size of the image by a factor of 'shrink_factor'
		 * using a local median filter.
		 */
		void median_shrink(int shrink_factor);

		/** multiplies by a radial function in fourier space */
		void apply_radial_func(float x0, float dx, vector < float >array, bool interp = true);

		/** calculates radial distribution. works for real and imaginary images. 
		 * 
		 * @param n number of points.
		 * @param x0 starting x coordinate.
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
		 * @return The radial distribution in an array.
		 */
		vector < float >calc_radial_dist(int n, float x0, float dx, float acen, float arange);

		void add(float f);
		void add(const EMData & em);

		void sub(float f);
		void sub(const EMData & em);

		void mult(float f);
		void mult(const EMData & em);

		void div(float f);
		void div(const EMData & em);

		float *get_data() const;
		void done_data();

		void update();
		void to_zero();
		void to_one();

		void dump_data(string filename);

		void add_incoherent(EMData * obj);

		vector < float >calc_fourier_shell_correlation(EMData * with);
		void calc_hist(vector < float >&hist, float hist_min = 0, float hist_max = 0,
					   bool add = false);
		void calc_az_dist(int n, float a0, float da, float *d, float rmin, float rmax);
#if 0
		void calc_rcf(EMData * with, vector < float >&sum_array);
#endif
		float calc_dist(EMData * second_img, int y_index = 0) const;
		EMData *calc_flcf(EMData * with, int radius = 50, string maskfilter = "SharpMask");

		EMData *convolute(EMData * with);

		bool has_ctff() const;

#if 0
		void clear_ctf();
		void create_ctf_map(CtfMapType type, XYData * sf = 0);
#endif

		float dot(EMData * with, bool evenonly = false);

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
		 */
		void common_lines(EMData * image1, EMData * image2, int mode = 0,
						  int steps = 180, bool horizontal = false);

		/** Finds common lines between 2 real images.
		 */
		void common_lines_real(EMData * image1, EMData * image2,
							   int steps = 180, bool horizontal = false);

		/** cut a slice out of a real 3D map */
		void cut_slice(EMData * map, float dz, Rotation * orientation = 0,
					   bool interpolate = true, float dx = 0, float dy = 0);

		/** Opposite of the cut_slice(). It will take a slice and insert
		 * the data into a real 3D map. It does not interpolate, it uses
		 * the nearest neighbor.
		 */
		void uncut_slice(EMData * map, float dz, Rotation * orientation = 0,
						 float dx = 0, float dy = 0);


		float calc_density_center();
		float calc_sigma_diff();
		
		Point < int >calc_min_location() const;
		Point < int >calc_max_location() const;

		int calc_min_index() const;
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
		
		
		float get_edge_mean() const;
		float get_circle_mean();

		void setup_insert_slice(int size);

		Ctf *get_ctf() const;
		void set_ctf(Ctf * ctf);

		Vec3 < float >get_translation() const;
		void set_translation(const Vec3 < float >&new_translation);
		void set_translation(float dx, float dy, float dz);

		Rotation get_rotation() const;
		void set_rotation(float alt, float az, float phi);

		void set_size(int nx, int ny, int nz);
		void set_path(const string & path);
		void set_pathnum(int n);

		EMData *get_row(int row_index) const;
		void set_row(const EMData * data, int row_index);

		EMData *get_col(int col_index) const;
		void set_col(const EMData * data, int n);


		EMObject get_attr(string key);
		void set_attr(string key, EMObject val);		
		Dict get_attr_dict();

		int get_xsize() const;
		int get_ysize() const;
		int get_zsize() const;
		int get_ndim() const;

		EMData *get_parent() const;
		void set_parent(EMData * new_parent);

		string get_name() const;
		void set_name(const string & name);

		float get_value_at(int x, int y, int z) const;
		float get_value_at(int x, int y) const;
		float get_value_at(size_t i) const;

		float sget_value_at(int x, int y, int z) const;
		float sget_value_at(int x, int y) const;

		float get_value_at_interp(float x, float y) const;
		float get_value_at_interp(float x, float y, float z) const;

		void set_value_at(int x, int y, int z, float v);
		void set_value_at(int x, int y, float v);

		bool is_complex() const;
		void set_complex(bool is_complex);

		bool is_complex_x() const;
		void set_complex_x(bool is_complex_x);

		bool is_flipped() const;
		void set_flipped(bool is_flipped);

		bool is_ri() const;
		void set_ri(bool is_ri);

		EMData & operator+=(float n);
		EMData & operator-=(float n);
		EMData & operator*=(float n);
		EMData & operator/=(float n);

		EMData & operator+=(const EMData & em);
		EMData & operator-=(const EMData & em);
		EMData & operator*=(const EMData & em);
		EMData & operator/=(const EMData & em);

		static vector < EMData * >read_images(string filename,
											  vector < int >img_indices = vector < int >(),
											  bool header_only = false);

		static vector < EMData * >read_images_ext(string filename, int img_index_start,
												  int img_index_end, bool header_only =
												  false, string ext = "");


	private:
		enum EMDataFlags {
			EMDATA_COMPLEX = 1 << 0,
			EMDATA_RI = 1 << 1,	       // real/imaginary or amp/phase
			EMDATA_BUSY = 1 << 2,	   // someone is modifying data
			EMDATA_SHARED = 1 << 3,	   // Stored in shared memory
			EMDATA_SWAPPED = 1 << 4,   // Data is swapped = may be offloaded if memory is tight,
			EMDATA_HASCTF = 1 << 6,	   // has CTF info
			EMDATA_NEEDUPD = 1 << 7,   // needs a realupdate= ,
			EMDATA_NEEDHIST = 1 << 8,  // histogram needs update
			EMDATA_NEWRFP = 1 << 9,	   // needs new rotational footprint
			EMDATA_NODATA = 1 << 10,   // no actual data
			EMDATA_COMPLEXX = 1 << 11, // 1D fft's in X
			EMDATA_FLIP = 1 << 12,	   // a flag only
			EMDATA_CHANGED = (EMDATA_NEEDUPD + EMDATA_NEEDHIST + EMDATA_NEWRFP)
		};

		void update_stat();
		void set_xyz_origin(float origin_x, float origin_y, float origin_z);
		void scale_pixel(float scale_factor) const;

	private:
		/** to store all image header information */
		mutable Dict attr_dict;
		float *rdata;	  /** image real data */
		float *supp;
		Ctf *ctf;		  /** CTF data */
		EMData *parent;

		EMData *rfp;
		int flags;

		int nx, ny, nz;	  /** image size */

		Vec3 < float >all_translation;
		/** translation from the original location */
		Vec3 < float >all_rotation;
		/** rotation (alt, az, phi) from the original locaton*/

		string name;
		string path;
		int pathnum;
	};


	EMData operator+(const EMData & em, float n);
	EMData operator-(const EMData & em, float n);
	EMData operator*(const EMData & em, float n);
	EMData operator/(const EMData & em, float n);

	EMData operator+(float n, const EMData & em);
	EMData operator-(float n, const EMData & em);
	EMData operator*(float n, const EMData & em);
	EMData operator/(float n, const EMData & em);

	EMData operator+(const EMData & a, const EMData & b);
	EMData operator-(const EMData & a, const EMData & b);
	EMData operator*(const EMData & a, const EMData & b);
	EMData operator/(const EMData & a, const EMData & b);

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
	
	inline float EMData::sget_value_at(int x, int y, int z) const
	{
		if (x < 0 || y < 0 || z < 0 || x >= nx || y >= ny || z >= nz) {
			return 0;
		}
		return rdata[x + y * nx + z * nx * ny];
	}

	inline float EMData::sget_value_at(int x, int y) const
	{
		if (x < 0 || y < 0 || x >= nx || y >= ny) {
			return 0;
		}
		return rdata[x + y * nx];
	}


	inline float EMData::get_value_at_interp(float xx, float yy) const
	{
		int x = static_cast < int >(floor(xx));
		int y = static_cast < int >(floor(yy));

		float p1 = sget_value_at(x, y);
		float p2 = sget_value_at(x + 1, y);
		float p3 = sget_value_at(x + 1, y + 1);
		float p4 = sget_value_at(x, y + 1);

		return Util::bilinear_interpolate(p1, p2, p3, p4, xx - x, yy - y);
	}

	inline float EMData::get_value_at_interp(float xx, float yy, float zz) const
	{
		int x = (int) floor(xx);
		int y = (int) floor(yy);
		int z = (int) floor(zz);
		float p1 = sget_value_at(x, y, z);
		float p2 = sget_value_at(x + 1, y, z);
		float p3 = sget_value_at(x, y + 1, z);
		float p4 = sget_value_at(x + 1, y + 1, z);

		float p5 = sget_value_at(x, y, z + 1);
		float p6 = sget_value_at(x + 1, y, z + 1);
		float p7 = sget_value_at(x, y + 1, z + 1);
		float p8 = sget_value_at(x + 1, y + 1, z + 1);

		return Util::trilinear_interpolate(p1, p2, p3, p4, p5, p6, p7, p8, xx - x, yy - y, zz - z);
	}

	inline void EMData::set_value_at(int x, int y, int z, float v)
	{
		rdata[x + y * nx + z * nx * ny] = v;
		flags |= EMDATA_NEEDUPD;
	}


	inline void EMData::set_value_at(int x, int y, float v)
	{
		rdata[x + y * nx] = v;
		flags |= EMDATA_NEEDUPD;
	}

	inline void EMData::update()
	{
		flags |= EMDATA_CHANGED;
	}


	inline bool EMData::is_complex() const
	{
		return (flags & EMDATA_COMPLEX);
	}

	inline bool EMData::is_complex_x() const
	{
		return (flags & EMDATA_COMPLEXX);
	}
	inline bool EMData::is_ri() const
	{
		return (flags & EMDATA_RI);
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

	inline bool EMData::is_flipped() const
	{
		return (flags & EMDATA_FLIP);
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

	inline void EMData::set_name(const string & new_name)
	{
		name = new_name;
	}

	inline string EMData::get_name() const
	{
		return name;
	}

	inline bool EMData::has_ctff() const
	{
		if (name[0] == '!' && name[1] == '$') {
			return true;
		}
		else {
			return false;
		}
	}

	inline Ctf *EMData::get_ctf() const
	{
		return ctf;
	}

	inline EMData *EMData::get_parent() const
	{
		return parent;
	}
	inline void EMData::set_parent(EMData * new_parent)
	{
		parent = new_parent;
	}


	inline Vec3 < float >EMData::get_translation() const
	{
		return all_translation;
	}

	inline void EMData::set_translation(const Vec3 < float >&t)
	{
		all_translation = t;
	}

	inline void EMData::set_translation(float dx, float dy, float dz)
	{
		all_translation = Vec3 < float >(dx, dy, dz);
	}

	inline Rotation EMData::get_rotation() const
	{
		return Rotation(all_rotation[0], all_rotation[1], all_rotation[2], Rotation::EMAN);
	}

	inline void EMData::set_rotation(float alt, float az, float phi)
	{
		all_rotation = Vec3 < float >(alt, az, phi);
	}

	inline void EMData::render_amp8_wrapper(long data, int x, int y, int xsize, int ysize,
											int bpl, float scale, int min_gray, int max_gray,
											float min_render, float max_render)
	{

		render_amp8((unsigned char *) data, x, y, xsize, ysize, bpl,
					scale, min_gray, max_gray, min_render, max_render);
	}
	
	inline void EMData::scale_pixel(float scale) const
	{
		attr_dict["spacing_row"] = EMObject((float) attr_dict["spacing_row"] * scale);
		attr_dict["spacing_col"] = EMObject((float) attr_dict["spacing_col"] * scale);
		attr_dict["spacing_sec"] = EMObject((float) attr_dict["spacing_sec"] * scale);
	}

	inline EMObject EMData::get_attr(string key)
	{
		update_stat();
		return attr_dict[key];
	}
	
	inline void EMData::set_attr(string key, EMObject val)
	{
		attr_dict[key] = val;
	}
	
	
}

#endif
