/**
 * $Id$
 */
#ifndef eman__salio_h__
#define eman__salio_h__ 1

#include "imageio.h"

namespace EMAN
{
	/** A SAL image is an image from Perkin Elmer PDS Microdensitometer.
	 * A SAL image consists of 2 files: 1 header file "X.hdr" and a
	 * data file "X.img". Header file is in ASCII format. Data file is
	 * in binary format.
	 *
	 * Each pair of hdr/img SAL files contains 1 2D image.
	*/

	class SalIO:public ImageIO
	{
	  public:
		SalIO(const string & filename, IOMode rw_mode = READ_ONLY);
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


#endif	//eman__salio_h__
