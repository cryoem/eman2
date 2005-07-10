/**
 * $Id$
 */
#ifndef eman__pgmio_h__
#define eman__pgmio_h__ 1

#include "imageio.h"

namespace EMAN
{
	/** A PGM file = header + data. Header is always in ASCII
	 * format. Data can be in either ASCII or BINARY format. Only
	 * Binary format is supported in EMAN so far.
	 * 
	 * A PGM file contains 1 2D image.
	 *
	 */
	class PgmIO:public ImageIO
	{
	  public:
		PgmIO(const string & filename, IOMode rw_mode = READ_ONLY);
		~PgmIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);

	  private:
		static const char *MAGIC_ASCII;
		static const char *MAGIC_BINARY;

		enum FileType
		{
			PGM_ASCII,
			PGM_BINARY,
			PGM_UNKNOWN_FILE
		};

		enum DataType
		{
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
		off_t file_offset;
	};
}

#endif	//eman__pgmio_h__
