/**
 * $Id$
 */
#ifndef eman__pifio_h__
#define eman__pifio_h__ 1


#include "imageio.h"
#include <stdio.h>

namespace EMAN
{
    class PifIO : public ImageIO
    {
    public:
	PifIO(string filename, IOMode rw_mode = READ_ONLY);
	~PifIO();

	DEFINE_IMAGEIO_FUNC;
	static bool is_valid(const void *first_block);
	
    private:
	enum {
	    PIF_MAGIC_NUM = 8
	};

	enum PifDataMode {
	    PIF_CHAR = 0,
	    PIF_SHORT = 1,
	    PIF_FLOAT_INT = 2,
	    PIF_SHORT_COMPLEX = 3,
	    PIF_FLOAT_INT_COMPLEX = 4,
	    PIF_BOXED_DATA = 6,
	    PIF_SHORT_FLOAT = 7,
	    PIF_SHORT_FLOAT_COMPLEX = 8,
	    PIF_FLOAT = 9,
	    PIF_FLOAT_COMPLEX = 10,
	    PIF_MAP_FLOAT_SHORT = 20,
	    PIF_MAP_FLOAT_INT = 21,
	    PIF_INVALID
	};

	// note there are no floats in these files. Floats are stored as ints with
	// a scaling factor.
	struct PifFileHeader
	{
	    int magic[2];	// magic number; identify PIF file
	    char scalefactor[16];	// to convert float ints -> floats
	    int nimg;		// # images in file
	    int endian;		// endianness, 0 -> vax,intel (little), 1 -> SGI, PowerPC (big)
	    char program[32];	// program which generated image
	    int htype;		// 1 - all images same number of pixels and depth, 0 - otherwise
	    int nx;		// number of columns
	    int ny;		// number of rows
	    int nz;		// number of sections
	    int mode;		// image data type                              
	    int pad[107];
	};

	struct PifColorMap
	{			// color map for depthcued images
	    short r[256];
	    short g[256];
	    short b[256];
	};

	struct PifImageHeader
	{
	    int nx;
	    int ny;
	    int nz;		// size of this image
	    int mode;		// image data type
	    int bkg;		// background value
	    int radius;		// boxed image radius
	    int xstart;
	    int ystart;
	    int zstart;		// starting number of each axis
	    int mx;
	    int my;
	    int mz;		// intervals along x,y,z
	    int xlen;
	    int ylen;
	    int zlen;		// cell dimensions (floatints)
	    int alpha;
	    int beta;
	    int gamma;		// angles (floatints)
	    int mapc;
	    int mapr;
	    int maps;		// axes->sections (1,2,3=x,y,z)
	    int min;
	    int max;
	    int mean;
	    int sigma;		// statistics (floatints)
	    int ispg;		// spacegroup
	    int nsymbt;		// bytes for symmetry ops
	    int xorigin;
	    int yorigin;	// origin (floatint)
	    char title[80];
	    char time[32];
	    char imagenum[16];	// unique micrograph number
	    char scannum[8];	// scan number of micrograph
	    int aoverb;
	    int mapabang;
	    int pad[63];	// later use
	};
    private:
	int get_mode_size(PifDataMode mode);
	bool is_float(int mode);
	void fseek_to(int image_index);
	int to_em_datatype(int pif_datatype);
	int to_pif_datatype(int em_datatype);

    private:
	string filename;
	IOMode rw_mode;
	PifFileHeader pfh;
	FILE *pif_file;
	int mode_size;
	bool is_big_endian;
	bool initialized;
	bool is_new_file;
	float real_scale_factor;
    };

}


#endif
