/**
 * $Id$
 */
#ifndef eman__pngio_h__
#define eman__pngio_h__ 1

#ifdef EM_PNG

#include <png.h>
#include "imageio.h"
#include <stdio.h>


namespace EMAN
{
	/** PngIO reads/writes a 2D PNG image. Currently 8-bit and 16-bit
	 * PNG read/write are supported.
	 */
	class PngIO:public ImageIO
	{
	  public:
		PngIO(const string & filename, IOMode rw_mode = READ_ONLY);
		~PngIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);

	  private:
		enum
		{
			PNG_BYTES_TO_CHECK = 8
		};

		enum BitDepthType
		{
			PNG_CHAR_DEPTH,
			PNG_SHORT_DEPTH,
			PNG_INVALID_DEPTH
		};

	  private:
		  string filename;
		IOMode rw_mode;
		FILE *png_file;

		bool initialized;

		png_structp png_ptr;
		png_infop info_ptr;
		png_infop end_info;

		png_uint_32 nx;
		png_uint_32 ny;
		BitDepthType depth_type;
		int number_passes;
	};

}

#endif

#endif
