#ifndef __emdata_h__
#define __emdata_h__

#include <string>
#include <vector>
#include <map>
#include <math.h>

#include "emobject.h"
#include "emutil.h"
#include "util.h"

using std::string;
using std::vector;
using std::map;

namespace EMAN {
    class ImageIO;
    class Region;
    class SimpleCtf;
    
    class EMData {
    public:
	static bool HEADER_ONLY;
	static bool HEADER_AND_DATA;
	static bool IS_3D;
	static bool NOT_3D;
	
    public:
	EMData();
	virtual ~EMData();
	
	int read_image(string filename, int img_index = 0, bool header_only = false,
		       Region* r = 0, bool is_3d = NOT_3D);
	int write_image(string filename, int img_index = 0,
			EMUtil::ImageType imgtype = EMUtil::IMAGE_UNKNOWN, bool header_only = false);

	map<string, EMObject> get_attr_dict() const;
	float* get_data() const;

	SimpleCtf* get_ctf();
	void set_ctf(const SimpleCtf& ctf);
	
	float get_value_at(int x, int y, int z);
	float get_value_at(int x, int y);
	
	float sget_value_at(int x, int y, int z);
	float sget_value_at(int x, int y);

	float get_value_at_interp(float x, float y);
	
	void set_value_at(int x, int y, int z, float v);
	void set_value_at(int x, int y, float v);

	void dump_data(string filename);

	static vector<EMData*> read_images_by_index(string filename, vector<int> img_indices, bool header_only = false);
	static vector<EMData*> read_images_by_ext(string filename, int img_index_start, int img_index_end,
						  bool header_only = false, string ext = "");
#if 0
	EMData& operator+=(float n);
        EMData& operator-=(float n);
        EMData& operator*=(float n);
        EMData& operator/=(float n);

        EMData& operator+=(const EMData& em);
        EMData& operator-=(const EMData& em);
        EMData& operator*=(const EMData& em);
        EMData& operator/=(const EMData& em);

        friend EMData& operator+(const EMData& em, float n);
        friend EMData& operator-(const EMData& em, float n);
        friend EMData& operator*(const EMData& em, float n);
        friend EMData& operator/(const EMData& em, float n);

        friend EMData& operator+(float n, const EMData& em);
        friend EMData& operator-(float n, const EMData& em);
        friend EMData& operator*(float n, const EMData& em);
        friend EMData& operator/(float n, const EMData& em);

        friend EMData& operator+(const EMData& a, const EMData& b);
        friend EMData& operator-(const EMData& a, const EMData& b);
        friend EMData& operator*(const EMData& a, const EMData& b);
        friend EMData& operator/(const EMData& a, const EMData& b);
#endif
	
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
	void set_size(int nx, int ny, int nz);
	
    private:
	map<string, EMObject> attr_dict;
	float* rdata;
	SimpleCtf* ctf;
	int flags;
    };

    
    inline float EMData::get_value_at(int x, int y, int z)
    {
	int nx = attr_dict["nx"].get_int();
	int ny = attr_dict["ny"].get_int();
	return rdata[x+y*nx+z*nx*ny]; 
    }


    inline float EMData::get_value_at(int x, int y)
    {
	int nx = attr_dict["nx"].get_int();
	return rdata[x+y*nx];
    }

	
    inline float EMData::sget_value_at(int x, int y, int z)
    {
	int nx = attr_dict["nx"].get_int();
	int ny = attr_dict["ny"].get_int();
	int nz = attr_dict["nz"].get_int();
    
	if (x < 0 || y < 0 || z < 0 || x >= nx || y >= ny || z >= nz) {
	    return 0;
	}
	return rdata[x+y*nx+z*nx*ny];
    }

    inline float EMData::sget_value_at(int x, int y)
    {
	int nx = attr_dict["nx"].get_int();
	int ny = attr_dict["ny"].get_int();
    
	if (x < 0 || y < 0 || x >= nx || y >= ny) {
	    return 0;
	}
	return rdata[x+y*nx];
    }


    inline float EMData::get_value_at_interp(float xx, float yy)
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
	int nx = attr_dict["nx"].get_int();
	int ny = attr_dict["ny"].get_int();
	
	rdata[x+y*nx+z*nx*ny] = v;
	flags |= EMDATA_NEEDUPD;
    }


    inline void EMData::set_value_at(int x, int y, float v)
    {
	int nx = attr_dict["nx"].get_int();
	
	rdata[x+y*nx] = v;
	flags |= EMDATA_NEEDUPD;
    }
    
}

#endif
