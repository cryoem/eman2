#ifndef eman__tiffio_h__
#define eman__tiffio_h__ 1

#ifdef EM_TIFF

#include "imageio.h"

typedef struct tiff TIFF;

namespace EMAN
{

    class TiffIO : public ImageIO
    {
    public:
	TiffIO(string filename, IOMode rw_mode = READ_ONLY);
	~TiffIO();

	DEFINE_IMAGEIO_FUNC;
	static bool is_valid(const void *first_block);
	
    private:
	enum {
	    TIFF_LITTLE_ENDIAN = 0x49,
	    TIFF_BIG_ENDIAN = 0x4d
	};

	string filename;
	IOMode rw_mode;
	TIFF *tiff_file;
	unsigned short bitspersample;
	bool is_big_endian;
	bool initialized;
    };
}

#endif

#endif
