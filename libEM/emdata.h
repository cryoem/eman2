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
	static bool HEADER_ONLY;
	static bool HEADER_AND_DATA;
	static bool IS_3D;
	static bool NOT_3D;
	static bool DATA_READ_ONLY;
	static bool DATA_READ_WRITE;


    public:
	EMData();
	virtual ~EMData();

	int read_image(string filename, int img_index = 0, bool header_only = false,
		       const Region * r = 0, bool is_3d = NOT_3D);

	int write_image(string filename, int img_index = 0,
			EMUtil::ImageType imgtype = EMUtil::IMAGE_UNKNOWN, bool header_only = false);

	int append_image(string filename, EMUtil::ImageType imgtype = EMUtil::IMAGE_UNKNOWN,
			 bool header_only = false);

	int filter(string filtername, const Dict & params = Dict());
	float cmp(string cmpname, const Dict & params);
	EMData *align(string aligner_name, const Dict & params, string comp_name = "");
	EMData *project(string projector_name, const Dict & params);
	
	EMData *copy(bool withfft = false, bool withparent = true);
	EMData *copy_head();

	EMData *get_clip(const Region & area);
	void insert_clip(EMData * block, const Point<int> & originn);

	EMData *do_fft();
	EMData *do_ift();

	Point<float> normalize_slice(EMData * slice, float alt, float az, float phi);

	int render_amp8(unsigned char *data, int x, int y, int xsize, int ysize,
			int bpl, float scale, int min_gray, int max_gray,
			float min_render, float max_render);
	int render_amp8_wrapper(int data, int x, int y, int xsize, int ysize,
				int bpl, float scale, int min_gray, int max_gray,
				float min_render, float max_render);
	int render_amp24(unsigned char *data, int x, int y, int xsize, int ysize,
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
	void translate(const Vec3<float> & translation);
	
	void rotate(float alt, float az, float phi);
	void rotate(const Rotation & r);

	void rotate_translate(float alt, float az, float phi, float dx, float dy, float dz);
	void rotate_translate(const Rotation& rotation, const Vec3<float>& translation);
	void rotate_translate(const Transform & xform);
	
	void rotate_x(int dx);
	int rotate_180();
	
	double dot_rotate_translate(EMData * data, float dx, float dy, float da);

	EMData *little_big_dot(EMData * little_img, bool do_sigma = false);
	EMData *do_radon();
	
	/** calculate Cross-correlation function (CCF).
	 *
	 * CCF is a 2D function that is obtained by forming the scalar
	 * cross-product of two images (i.e., the sum of products of
	 * equivalent pixels) as a function of a 2D shift vector. The
	 * CCF is often used to achieve alignment between two images,
	 * since it displays a high value (a peak) at the place where
	 * a motif contained in both images come into register.
	 */
	EMData *calc_ccf(EMData * with, bool tocorner = false, EMData * filter = 0);
	EMData *make_rotational_footprint(bool premasked = false, bool unwrap = true);
	EMData *calc_ccfx(EMData * with, int y0 = 0, int y1 = -1, bool nosum = false);
	EMData *calc_mutual_correlation(EMData * with, bool tocorner = false,
					EMData * filter = 0);
	EMData *unwrap(int r1 = -1, int r2 = -1, int xs = -1, int dx = 0,
		       int dy = 0, bool do360 = false);

	int mean_shrink(int shrink_factor);
	int median_shrink(int shrink_factor);

	void apply_radial_func(float x0, float dx, vector<float> array, bool interp = true);
	vector<float> calc_radial_dist(int n, float x0, float dx);
	vector<float> calc_radial_dist(int n, float x0, float dx, float acen, float amwid);
	
	int add(float f);
	int add(const EMData & em);

	int sub(float f);
	int sub(const EMData & em);

	int mult(float f);
	int mult(const EMData & em);

	int div(float f);
	int div(const EMData & em);

	float *get_data() const;
	void done_data();

	void update();
	void to_zero();
	void to_one();

	void dump_data(string filename);

	int add_incoherent(EMData * obj);

	vector<float> calc_fourier_shell_correlation(EMData * with);
	void calc_hist(vector<float> & hist, float hist_min = 0, float hist_max = 0, bool add =
		       false);
	int calc_az_dist(int n, float a0, float da, float *d, float rmin, float rmax);
#if 0
	void calc_rcf(EMData * with, vector<float> & sum_array);
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

	int common_lines(EMData * image1, EMData * image2, int mode = 0,
			 int steps = 180, bool horizontal = false);

	int common_lines_real(EMData * image1, EMData * image2,
			      int steps = 180, bool horizontal = false);

	int cut_slice(EMData * map, float dz, Rotation * orientation = 0,
		      bool interpolate = true, float dx = 0, float dy = 0);

	int uncut_slice(EMData * map, float dz, Rotation * orientation = 0,
			float dx = 0, float dy = 0);


	float get_edge_mean() const;
	float get_circle_mean();

	void setup_insert_slice(int size);

	Ctf *get_ctf() const;
	void set_ctf(Ctf * ctf);

	Vec3<float> get_translation() const;
	void set_translation(const Vec3<float> & new_translation);
	void set_translation(float dx, float dy, float dz);
	
	Rotation get_rotation() const;
	void set_rotation(float alt, float az, float phi);
	
	void set_size(int nx, int ny, int nz);
	void set_path(const string & path);
	void set_pathnum(int n);

	EMData *get_row(int row_index) const;
	void set_row(const EMData * d, int row_index);

	EMData *get_col(int col_index) const;
	void set_col(const EMData * d, int n);

	float get_align_score() const;
	void set_align_score(float score);

	float get_density_center();
	float get_sigma_diff();

	Dict get_attr_dict();
	void set_attr_dict(string key, EMObject val);

	int get_average_nimg() const;
	void set_average_nimg(int n);
	
	float get_max();
	float get_min();
	float get_mean();
	float get_sigma();
	float get_skewness();
	float get_kurtosis();
	
	Point<int> get_min_location() const;
	Point<int> get_max_location() const;

	int get_min_index() const;
	int get_max_index() const;

	int get_xsize() const;
	int get_ysize() const;
	int get_zsize() const;
	int get_ndim() const;

	float get_xorigin() const;
	float get_yorigin() const;
	float get_zorigin() const;

	float get_xpixel() const;
	float get_ypixel() const;
	float get_zpixel() const;
	
	EMData *get_parent() const;
	void set_parent(EMData * new_parent);

	string get_name() const;
	void set_name(const string & name);

	float get_value_at(int x, int y, int z) const;
	float get_value_at(int x, int y) const;

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

	static vector<EMData *> read_images_by_index(string filename, vector<int>img_indices,
						     bool header_only = false);

	static vector<EMData *> read_images_by_ext(string filename, int img_index_start,
						   int img_index_end, bool header_only =
						   false, string ext = "");


    private:
	enum EMDataFlags {
	    EMDATA_COMPLEX = 1 << 0,
	    EMDATA_RI = 1 << 1,	        // real/imaginary or amp/phase
	    EMDATA_BUSY = 1 << 2,	// someone is modifying data
	    EMDATA_SHARED = 1 << 3,	// Stored in shared memory
	    EMDATA_SWAPPED = 1 << 4,	// Data is swapped = may be offloaded if memory is tight,
	    EMDATA_HASCTF = 1 << 6,	// has CTF info
	    EMDATA_NEEDUPD = 1 << 7,	// needs a realupdate= ,
	    EMDATA_NEEDHIST = 1 << 8,	// histogram needs update
	    EMDATA_NEWRFP = 1 << 9,	// needs new rotational footprint
	    EMDATA_NODATA = 1 << 10,	// no actual data
	    EMDATA_COMPLEXX = 1 << 11,	// 1D fft's in X
	    EMDATA_FLIP = 1 << 12,	// a flag only
	    EMDATA_CHANGED = (EMDATA_NEEDUPD + EMDATA_NEEDHIST + EMDATA_NEWRFP)
	};

	int update_stat();
	void set_xyz_origin(float origin_x, float origin_y, float origin_z);
	void scale_pixel(float scale_factor) const;

    private:
	mutable Dict attr_dict;
	float *rdata;
	float *supp;
	Ctf *ctf;
	EMData *parent;
	
	EMData *rfp;
	int flags;

	int nx;
	int ny;
	int nz;

	int average_nimg;
	
	Vec3<float> all_translation;  // from the original location
	Vec3<float> all_rotation;     // (alt, az, phi) from the original locaton
	
	float align_score;
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

    inline float EMData::get_xorigin() const
    {
	return attr_dict["origin_row"].get_float();
    }

    
    inline float EMData::get_yorigin() const
    {
	return attr_dict["origin_col"].get_float();
    }

    inline float EMData::get_zorigin() const
    {
	return attr_dict["origin_sec"].get_float();
    }

    
    inline float EMData::get_xpixel() const
    {
	return attr_dict["spacing_row"].get_float();
    }

    
    inline float EMData::get_ypixel() const
    {
	return attr_dict["spacing_col"].get_float();
    }

    
    inline float EMData::get_zpixel() const
    {
	return attr_dict["spacing_sec"].get_float();
    }

    
    inline float EMData::get_value_at(int x, int y, int z) const
    {
	return rdata[x + y * nx + z * nx * ny];
    }


    inline float EMData::get_value_at(int x, int y) const
    {
	return rdata[x + y * nx];
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
	int x = static_cast<int>(floor(xx));
	int y = static_cast<int>(floor(yy));

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

    inline int EMData::get_average_nimg() const
    {
	return average_nimg;
    }

    inline void EMData::set_average_nimg(int n)
    {
	average_nimg = n;
    }
    
    inline float EMData::get_max()
    {
	update_stat();
	return attr_dict["maximum"].get_float();
    }

    inline float EMData::get_min()
    {
	update_stat();
	return attr_dict["minimum"].get_float();
    }
    
    inline float EMData::get_mean()
    {
	update_stat();
	return attr_dict["mean"].get_float();
    }

    inline float EMData::get_sigma()
    {
	update_stat();
	return attr_dict["sigma"].get_float();
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

    inline void EMData::set_align_score(float score)
    {
	align_score = score;
    }
    inline float EMData::get_align_score() const
    {
	return align_score;
    }

    inline Vec3<float> EMData::get_translation() const
    {
	return all_translation;
    }

    inline void EMData::set_translation(const Vec3<float> & t)
    {
	all_translation = t;
    }
    
    inline void EMData::set_translation(float dx, float dy, float dz)
    {
	all_translation = Vec3<float>(dx, dy, dz);
    }
    
    inline Rotation EMData::get_rotation() const
    {
	return Rotation(all_rotation[0], all_rotation[1], all_rotation[2],
			Rotation::EMAN);
    }
    
    inline void EMData::set_rotation(float alt, float az, float phi) 
    {
	all_rotation = Vec3<float>(alt, az, phi);
    }

    
    inline int EMData::render_amp8_wrapper(int data, int x, int y, int xsize, int ysize,
					   int bpl, float scale, int min_gray, int max_gray,
					   float min_render, float max_render)
    {
#if 0
	return render_amp8((unsigned char*)data, x, y, xsize, ysize, bpl,
			   scale, min_gray, max_gray, min_render, max_render);
#endif
	return 0;
    }
    
    inline void EMData::scale_pixel(float scale) const
    {	
	attr_dict["spacing_row"] = EMObject(attr_dict["spacing_row"].get_float() * scale);
	attr_dict["spacing_col"] = EMObject(attr_dict["spacing_col"].get_float() * scale);
	attr_dict["spacing_sec"] = EMObject(attr_dict["spacing_sec"].get_float() * scale);
    }

}

#endif
