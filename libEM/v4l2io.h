/**
 * $Id$
 */
#ifndef eman__v4l2_h__
#define eman__v4l2_h__ 1

#include "imageio.h"
#include <asm/types.h>
#include <linux/videodev2.h>

namespace EMAN
{
	/** Read-only. Acquires images from the V4L2 interface in real-time
		(video4linux). ie - this will read from a framegrabber, etc.
	 */
	class V4L2IO:public ImageIO
	{
	  public:
		V4L2IO(const string & filename, IOMode rw_mode = READ_ONLY);
		~V4L2IO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);

	  private:
		char *filename;
		int v4l_file;

		bool initialized;

		int nx;
		int ny;

		struct v4l2_capability cap;
		struct v4l2_cropcap cropcap;
		struct v4l2_crop crop;
		struct v4l2_format fmt;
	};
}

#endif	//eman__pgmio_h__
