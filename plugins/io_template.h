#ifndef eman__xyzio_h__
#define eman__xyzio_h__1

#include "imageio.h"
#include <stdio.h>

namespace EMAN
{
	/** XYZIO is a sample Image IO class. It defines all required API
     * that you may need to implement.
     */

	class XYZIO:public ImageIO
	{
	  public:
		XYZIO(string filename, IOMode rw_mode = READ_ONLY);
		~XYZIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);

	  private:
		  string filename;
		IOMode rw_mode;
		FILE *xyz_file;

		bool is_big_endian;
		bool initialized;
	};

}


#endif
