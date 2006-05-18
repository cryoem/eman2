/**
 * $Id$
 */
#ifndef eman__jpegio_h__
#define eman__jpegio_h__ 1

#ifdef EM_JPEG

extern "C" {
#include <jpeglib.h>
}
#include "imageio.h"


namespace EMAN
{
	/** JpegIO reads/writes a 2D JPEG image. Currently write-only */
	class JpegIO:public ImageIO
	{
	  public:
		JpegIO(const string & filename, IOMode rw_mode = WRITE_ONLY);
		~JpegIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);

	  private:
		string filename;
		IOMode rw_mode;
		FILE *jpeg_file;

		bool initialized;
		float rendermin;
		float rendermax;

		struct jpeg_compress_struct cinfo;
		struct jpeg_error_mgr jerr;

	};

}

#endif	//EM_JPEG

#endif	//eman__jpegio_h__
