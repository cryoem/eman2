/**
 * $Id$
 */
#ifndef eman__salio_h__
#define eman__salio_h__ 1


#include "imageio.h"
#include <stdio.h>

namespace EMAN
{
	// to read images from Perkin Elmer PDS Microdensitometer

	class SalIO:public ImageIO
	{
	  public:
		SalIO(string filename, IOMode rw_mode = READ_ONLY);
		~SalIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);

	  private:
		static const char *HDR_EXT;
		static const char *IMG_EXT;
		static const char *MAGIC;

		enum ScanAxis
		{
			X_SCAN_AXIS,
			Y_SCAN_AXIS
		};

		enum ScanMode
		{
			NON_RASTER_SCAN,
			RASTER_SCAN
		};

	  private:
		  string filename;
		IOMode rw_mode;
		FILE *sal_file;
		bool initialized;
		int nx;
		int ny;
		int record_length;
		ScanMode scan_mode;
		float pixel;
	};

}


#endif
