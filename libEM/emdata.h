#ifndef __emdata_h__
#define __emdata_h__

#include <string>
#include <vector>
#include <map>
#include "emobject.h"
#include "emutil.h"

using std::string;
using std::vector;
using std::map;

namespace EMAN {
    class ImageIO;
    class Region;
    class SimpleCtf;
    
    class EMData {
    public:
	EMData();
	virtual ~EMData();

	int read_image(string filename, int img_index = 0, bool nodata = false,
		       Region* r = 0, bool is_3d = false);
	int write_image(string filename, int img_index = 0,
			EMUtil::ImageType imgtype = EMUtil::IMAGE_UNKNOWN, bool nodata = false);

	map<string, EMObject> get_attr_dict() const;
	float* get_data() const;

	void dump_data(string filename);

	static vector<EMData*> read_images(string filename, int img_indices[], int nimg, bool nodata = false);
	static vector<EMData*> read_images(string filename, int img_index_start,
					   int img_index_end, bool nodata = false);
	static vector<EMData*> read_images_ext(string filename, string ext, int img_index_start,
					       int img_index_end, bool nodata = false);
	
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

    
}

#endif
