/**
 * $Id$
 */
#ifndef eman__amiraio_h__
#define eman__amiraio_h__ 1


#include "imageio.h"
#include <stdio.h>

namespace EMAN
{

	/** Amira file = ASCII header + binary data.
	 * Its first line has some magic
	 * name to label it as an Amira image. The first few lines of the
	 * file is the ASCII header. Followed the header is the data in
	 * binary format. The data has nx x ny x nz pixels.
	 *
	 * An Amira file has only 1 2D or 3D image.
	 */
	class AmiraIO:public ImageIO
	{
	  public:
		AmiraIO(const string & filename, IOMode rw_mode = READ_ONLY);
		~AmiraIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);

	  private:
		string filename;
		IOMode rw_mode;
		FILE *amira_file;

		bool is_big_endian;
		bool initialized;
		int nx;
		int ny;
		int nz;

		static const char *MAGIC;
	};

}


#endif
