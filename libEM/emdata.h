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
	
	EMData* copy(bool withfft = false, bool withparent = true) const;
	
	EMData* get_clip(const Region& area);
	void insert_clip(EMData* block, const Point<int>& originn);

	int read_image(string filename, int img_index = 0, bool header_only = false,
		       Region* r = 0, bool is_3d = NOT_3D);
	int write_image(string filename, int img_index = 0,
			EMUtil::ImageType imgtype = EMUtil::IMAGE_UNKNOWN, bool header_only = false);

	int filter(string filtername, const Dict& params);
	
	void normalize();

	bool is_complex() const;
	void ri2ap();
	void ap2ri();

	EMData* get_parent() const;
	
	EMData* do_fft();
	EMData* do_ift();
	void gimme_fft();

	EMData* calc_ccf(EMData* with, bool tocorner = false, EMData* filter = 0);
	EMData* make_rotational_footprint(bool premasked = false, bool unwrap = true);
	EMData* calc_ccfx(EMData* with, int y0 = 0, int y1 = -1, bool nosum = false);
	EMData* unwrap(int r1 = -1, int r2 = -1, int xs = -1, int dx = 0, int dy = 0,  bool do360 = false);

	float* calc_fourier_shell_correlation(EMData *with);
	
	void apply_radial_func(int, float, vector<float> array);
	    
	void set_talign_params(float dx, float dy);
	void set_talign_params(float dx, float dy, float dz);

	void set_ralign_params(float alt, float az, float phi);
	void set_ralign_params(const Rotation& r);

	int add(float f);
	int add(const EMData& em);
	int sub(const EMData& em);
	
	int mult(float f);
	int mult(const EMData& em);

	int div(float f);
	int div(const EMData& em);
	
	float* get_data() const;
	void done_data();
	
	SimpleCtf* get_ctf();
	void set_ctf(const SimpleCtf& ctf);

	void set_size(int nx, int ny, int nz);
	
	Dict get_attr_dict();
	
	float get_value_at(int x, int y, int z) const;
	float get_value_at(int x, int y) const;
	
	float sget_value_at(int x, int y, int z) const;
	float sget_value_at(int x, int y) const;

	float get_value_at_interp(float x, float y) const;
	
	void set_value_at(int x, int y, int z, float v);
	void set_value_at(int x, int y, float v);

	int get_x() const;
	int get_y() const;
	int get_z() const;

	Vec3f get_translation() const;
	void set_translation(const Vec3f& t);

	

	
	void dump_data(string filename);

	static vector<EMData*> read_images_by_index(string filename, vector<int> img_indices, bool header_only=false);
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
	    EMDATA_COMPLEXX	= 1<<11		// 1D fft's in X
	};

	int update_stat();

    private:
	mutable map<string, EMObject> attr_dict;
	float* rdata;
	SimpleCtf* ctf;
	int flags;
	int rocount;
	int nx;
	int ny;
	int nz;
	Vec3f translation;
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

	float v = Util::bilinear_interpolate(sget_value_at(x, y),
					   sget_value_at(x+1, y),
					   sget_value_at(x+1, y+1),
					   sget_value_at(x, y+1),
					   xx-x,
					   yy-y);
	return v;
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
}

#endif
