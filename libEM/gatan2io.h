/**
 * $Id$
 */
#ifndef eman__gatan2io_h__
#define eman__gatan2io_h__ 1

#include "imageio.h"

namespace EMAN
{

	/** Gatan2 Image file = header + data.
	 * header is defined in Gatan2Header. data is nx by ny.
	 * A Gatan2 file contains 1 2D image.
	 */
	   
	class Gatan2IO:public ImageIO
	{
	  public:
		Gatan2IO(const string & filename, IOMode rw_mode = READ_ONLY);
		~Gatan2IO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);
	  private:
		enum DataType
		{
			GATAN2_SHORT = 1,
			GATAN2_FLOAT = 2,
			GATAN2_COMPLEX = 3,
			GATAN2_PACKED_COMPLEX = 5,
			GATAN2_CHAR = 6,
			GATAN2_INT = 7,
			GATAN2_INVALID
		};

		struct Gatan2Header
		{
			short version;
			short un1;
			short un2;
			short nx;
			short ny;
			short len;
			short type;
		};

		int to_em_datatype(int gatan_type);

	  private:
		string filename;
		IOMode rw_mode;
		FILE *gatan2_file;
		Gatan2Header gatanh;

		bool is_big_endian;
		bool initialized;
	};

}


#endif	//eman__gatan2io_h__
