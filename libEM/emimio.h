#ifndef __emimio_h__
#define __emimio_h__


#include "imageio.h"
#include <stdio.h>

namespace EMAN {
    
    class EmimIO : public ImageIO {
    public:
	EmimIO(string filename, IOMode rw_mode = READ_ONLY);
	~EmimIO();
	
	DEFINE_IMAGEIO_FUNC;
	
    private:
	static const char* MAGIC;
	enum { EMIM_COMPLEX = 1, PARM_NONE = 0, PARM_NORM = 1, NUM_INT_IN_FILE_HEADER = 16 };
	
	struct EmimFileHeader {
	    char magic[4];		// 'EMIM'
	    int order;			// This defines the byte order, value is 2 on correct machine
	    int count;			// # images in the file (24 bits max)
	    int nx;
	    int ny;
	    int nz;		        // image size (all images in 1 file same size)
	    int flag;			// flags are the same as well
	    float pixel;		// pixel/voxel size in A
	    int misc[7];		// misc usage
	    int headlen;		// length of individual image headers (currently 166 x 4)
	    char name[64];		// optional description of entire set of images
	};

	struct EmimImageHeader {
	    char name[80];		// data set name
	    unsigned int time;	        // time at which image was acquired
	    int mgnum;			// micrograph number or CCD frame
	    int clipnum;		// id of clip within micrograph/frame
	    int id;			// another id number
	    float dx;
	    float dy;
	    float dz;		       // translation
	    float alt,az,phi;	       // orientation
	    short nparm;	       // number of parameter sets
	    int ptype[4];	       // type of each parameter set
	    float param[4][32];	       // Parameterized CTF/SNR, etc...
	};

    private:
	string filename;
	IOMode rw_mode;
	EmimFileHeader efh;
	FILE* emim_file;
	bool is_big_endian;
	bool initialized;
    };
}


#endif
