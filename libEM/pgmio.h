#ifndef eman__pgmio_h__
#define eman__pgmio_h__ 1

#include "imageio.h"
#include <stdio.h>

namespace EMAN
{
    class PgmIO : public ImageIO
    {
    public:
	PgmIO(string filename, IOMode rw_mode = READ_ONLY);
	~PgmIO();

	DEFINE_IMAGEIO_FUNC;
	static bool is_valid(const void *first_block);
	
    private:
	static const char *MAGIC_ASCII;
	static const char *MAGIC_BINARY;

	enum FileType {
	    PGM_ASCII,
	    PGM_BINARY,
	    PGM_UNKNOWN_FILE
	};
	
	enum DataType {
	    PGM_UCHAR,
	    PGM_USHORT,
	    PGM_UNKNOWN_TYPE
	};
	
    private:
	string filename;
	IOMode rw_mode;
	FILE *pgm_file;

	bool is_big_endian;
	bool initialized;

	int nx;
	int ny;
	int maxval;
	int minval;

	DataType datatype;
	long file_offset;
    };
}

#endif
