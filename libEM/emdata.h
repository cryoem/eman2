#ifndef eman__emdata_h__
#define eman__emdata_h__ 1

#include <string>
#include <vector>
#include <map>
#include <math.h>

#include "emobject.h"
#include "emutil.h"
#include "util.h"
#include "transform.h"

using std::string;
using std::vector;
using std::map;

namespace EMAN {
    class ImageIO;
    class SimpleCtf;
    class XYData;
    
    class EMData {
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
		       Region* r = 0, bool is_3d = NOT_3D);
	
	int write_image(string filename, int img_index = 0,
			EMUtil::ImageType imgtype = EMUtil::IMAGE_UNKNOWN, bool header_only = false);

	int filter(string filtername, const Dict& params);

	float cmp(string cmpname, const Dict& params);

	EMData* align(string aligner_name, const Dict& params, string comp_name = "");
	
	EMData* copy(bool withfft = false, bool withparent = true) const;
	EMData* copy_head();
	
	EMData* get_clip(const Region& area);
	void insert_clip(EMData* block, const Point<int>& originn);

	void normalize();

	void calc_hist(vector<float>& hist, float hist_min = 0, float hist_max = 0, bool add = false);

	int render_amp8(unsigned char* data, int x, int y, int xsize, int ysize,
			int bpl, float scale, int min_gray, int max_gray,
			float min_render, float max_render);
    
	int render_amp24(unsigned char* data, int x, int y, int xsize, int ysize,
			 int bpl, float scale, int min_gray, int max_gray,
			 float min_render, float max_render,
			 void *ref, void cmap(void*, int coord, unsigned char *tri));
    
	int render_pha24(unsigned char* data, int x, int y, int xsize, int ysize, 
			 int bpl, float scale, float min_render, float max_render);
	
	int calc_az_dist(int n, float a0, float da, float *d, float rmin, float rmax);
	
	bool is_complex() const;
	void set_complex(bool is_complex);

	bool is_complex_x() const;
	void set_complex_x(bool is_complex_x);	
	
	void set_ri(bool is_ri);
	void ri2ap();
	void ap2ri();

	EMData* get_parent() const;
	void set_parent(EMData* new_parent) { parent = new_parent; }

	void setup4slice(bool redo = true);
	void to_corner();

	float* get_supp() const { return 0; }
	
	EMData* do_fft();
	EMData* do_ift();
	void gimme_fft();
	
	void rotate_x(int dx);
	void rotate_translate(float scale = 1.0, float dxc = 0,
			      float dyc = 0, float dzc = 0, int r = 0);
	void fast_rotate_translate(float scale = 1.0);
	double dot_rotate_translate(EMData* data, float dx, float dy, float da);

	int fast_translate(bool inplace = true);
	int rotate_180();
	
	EMData* do_radon();
	EMData* vertical_acf(int maxdy);
	
	void vertical_flip();
	void horizontal_flip();
	
	void set_flipped(bool flipped)
	{
	    if (flipped) {
		flags |= EMDATA_FLIP;
	    }
	    else {
		flags&=~EMDATA_FLIP;
	    }
	}
	
	float dot(EMData* with, bool evenonly = false) { return 0; }
	
	EMData* calc_ccf(EMData* with, bool tocorner = false, EMData* filter = 0);
	EMData* make_rotational_footprint(bool premasked = false, bool unwrap = true);
	EMData* calc_ccfx(EMData* with, int y0 = 0, int y1 = -1, bool nosum = false);
#if 0
	void calc_rcf(EMData *with, vector<float>& sum_array);
#endif
	EMData* unwrap(int r1 = -1, int r2 = -1, int xs = -1, int dx = 0, int dy = 0,  bool do360 = false);

	vector<float> calc_fourier_shell_correlation(EMData *with);
	
	int mean_shrink(int shrink_factor);
	int median_shrink(int shrink_factor);
	
	void apply_radial_func(int, float, vector<float> array);
	    
	void set_talign_params(float dx, float dy);
	void set_talign_params(float dx, float dy, float dz);

	void set_ralign_params(float alt, float az, float phi);
	void set_ralign_params(const Rotation& r);

	void set_align_score(float score) { align_score = score; }
	float get_align_score() const { return align_score; }
	
	int add(float f);
	int add(const EMData& em);
	int sub(const EMData& em);
	
	int mult(float f);
	int mult(const EMData& em);

	int div(float f);
	int div(const EMData& em);
	
	float* get_data() const;
	void done_data();

	void update();
	void to_zero();
	void to_one();
	
	SimpleCtf* get_ctf();
	void set_ctf(const SimpleCtf& ctf);

	void set_size(int nx, int ny, int nz);

	void set_pixel_size(float pixel_size);
	float get_pixel_size() const;
	
	Point<int> get_min_location() const;
	Point<int> get_max_location() const;

	int get_min_index() const;
	int get_max_index() const;
	
	Dict get_attr_dict();
	
	float get_value_at(int x, int y, int z) const;
	float get_value_at(int x, int y) const;
	
	float sget_value_at(int x, int y, int z) const;
	float sget_value_at(int x, int y) const;

	float get_value_at_interp(float x, float y) const;
	float get_value_at_interp(float x, float y, float z) const;
	
	void set_value_at(int x, int y, int z, float v);
	void set_value_at(int x, int y, float v);

	int get_x() const;
	int get_y() const;
	int get_z() const;

	Vec3f get_translation() const;
	void set_translation(const Vec3f& t);

	Rotation get_rotation() const { return rotation; }
	Vec3f get_trans_align() const { return trans_align; }
	
	void EMData::set_name(string name);
	void EMData::set_path(string path);

	void dump_data(string filename);

	static vector<EMData*> read_images_by_index(string filename, vector<int> img_indices,
						    bool header_only=false);
	static vector<EMData*> read_images_by_ext(string filename, int img_index_start, int img_index_end,
						  bool header_only = false, string ext="");

	EMData& operator+=(float n);
        EMData& operator-=(float n);
        EMData& operator*=(float n);
        EMData& operator/=(float n);

        EMData& operator+=(const EMData& em);
        EMData& operator-=(const EMData& em);
        EMData& operator*=(const EMData& em);
        EMData& operator/=(const EMData& em);

        friend EMData operator+(const EMData& em, float n);
        friend EMData operator-(const EMData& em, float n);
        friend EMData operator*(const EMData& em, float n);
        friend EMData operator/(const EMData& em, float n);

        friend EMData operator+(float n, const EMData& em);
        friend EMData operator-(float n, const EMData& em);
        friend EMData operator*(float n, const EMData& em);
        friend EMData operator/(float n, const EMData& em);

        friend EMData operator+(const EMData& a, const EMData& b);
        friend EMData operator-(const EMData& a, const EMData& b);
        friend EMData operator*(const EMData& a, const EMData& b);
        friend EMData operator/(const EMData& a, const EMData& b);

	
    private:
	enum EMDataFlags {
	    EMDATA_COMPLEX      = 1<<0, 
	    EMDATA_RI	        = 1<<1,		// real/imaginary or amp/phase
	    EMDATA_BUSY	        = 1<<2,		// someone is modifying data
	    EMDATA_SHARED	= 1<<3,		// Stored in shared memory
	    EMDATA_SWAPPED	= 1<<4,		// Data is swapped = may be offloaded if memory is tight,
	    EMDATA_NEWFFT	= 1<<5,		// Data has changed, redo fft
	    EMDATA_HASCTF	= 1<<6,		// has CTF info
	    EMDATA_NEEDUPD	= 1<<7,		// needs a realupdate= ,
	    EMDATA_NEEDHIST	= 1<<8,		// histogram needs update
	    EMDATA_NEWRFP	= 1<<9,		// needs new rotational footprint
	    EMDATA_NODATA	= 1<<10,	// no actual data
	    EMDATA_COMPLEXX	= 1<<11,       	// 1D fft's in X
	    EMDATA_FLIP         = 1<<12,
	    EMDATA_CHANGED      = (EMDATA_NEWFFT+EMDATA_NEEDUPD+EMDATA_NEEDHIST+EMDATA_NEWRFP)
	};

	int update_stat();
	void set_xyz_origin(float origin_x, float origin_y, float origin_z);
	
    private:
	mutable map<string, EMObject> attr_dict;
	float* rdata;
	float* supp;
	SimpleCtf* ctf;
	EMData* parent;
	EMData* fft;
	int flags;
	float pixel_size;
	
	int nx;
	int ny;
	int nz;
	
	Vec3f translation;
	Rotation rotation;
	Vec3f trans_align;
	float align_score;
    };


    inline int EMData::get_x() const
    {
	return nx;
    }
    
    inline int EMData::get_y() const
    {	
	return ny;
    }    

    inline int EMData::get_z() const
    {
	return nz;
    }

    inline Vec3f EMData::get_translation() const { return translation; }
    
    inline void EMData::set_translation(const Vec3f& t) { translation = t; }
    
    
    inline float EMData::get_value_at(int x, int y, int z) const
    {
	return rdata[x+y*nx+z*nx*ny]; 
    }


    inline float EMData::get_value_at(int x, int y) const
    {
	return rdata[x+y*nx];
    }

	
    inline float EMData::sget_value_at(int x, int y, int z) const
    {
	if (x < 0 || y < 0 || z < 0 || x >= nx || y >= ny || z >= nz) {
	    return 0;
	}
	return rdata[x+y*nx+z*nx*ny];
    }

    inline float EMData::sget_value_at(int x, int y) const
    {
	if (x < 0 || y < 0 || x >= nx || y >= ny) {
	    return 0;
	}
	return rdata[x+y*nx];
    }


    inline float EMData::get_value_at_interp(float xx, float yy) const
    {
	int x = static_cast<int>(floor(xx));
	int y = static_cast<int>(floor(yy));

	float p1 = sget_value_at(x, y);
	float p2 = sget_value_at(x+1, y);
	float p3 = sget_value_at(x+1, y+1);
	float p4 = sget_value_at(x, y+1);
	
	return Util::bilinear_interpolate(p1, p2, p3, p4, xx-x, yy-y);
    }

    inline float EMData::get_value_at_interp(float xx, float yy, float zz) const
    {
	int x = (int)floor(xx);
	int y = (int)floor(yy);
	int z = (int)floor(zz);
	float p1 = sget_value_at(x, y, z);
	float p2 = sget_value_at(x+1, y, z);
	float p3 = sget_value_at(x, y+1, z);
	float p4 = sget_value_at(x+1, y+1, z);
	
	float p5 = sget_value_at(x, y, z+1);
	float p6 = sget_value_at(x+1, y, z+1);
	float p7 = sget_value_at(x, y+1, z+1);
	float p8 = sget_value_at(x+1, y+1, z+1);
	
	return Util::trilinear_interpolate(p1, p2, p3, p4, p5, p6, p7, p8, xx-x, yy-y, zz-z);
    }
    
    inline void EMData::set_value_at(int x, int y, int z, float v)
    {	
	rdata[x+y*nx+z*nx*ny] = v;
	flags |= EMDATA_NEEDUPD;
    }


    inline void EMData::set_value_at(int x, int y, float v)
    {	
	rdata[x+y*nx] = v;
	flags |= EMDATA_NEEDUPD;
    }

    inline void EMData::update()
    {
	flags |= EMDATA_CHANGED;
    }

    inline void EMData::set_pixel_size(float new_pixel_size) { pixel_size = new_pixel_size; }
    inline float EMData::get_pixel_size() const { return pixel_size; }


    inline bool EMData::is_complex() const { return (flags & EMDATA_COMPLEX); }
    inline bool EMData::is_complex_x() const { return (flags & EMDATA_COMPLEXX); }
    
    inline void EMData::set_ri(bool is_ri) 
    {
	if (is_ri) {
	    flags |= EMDATA_RI;
	}
	else {
	    flags &= ~EMDATA_RI;
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

    
}

#endif
