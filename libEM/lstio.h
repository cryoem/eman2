/**
 * $Id$
 */
#ifndef eman__lstio_h__
#define eman__lstio_h__ 1

#include "imageio.h"

namespace EMAN
{
	/** A LST file is an ASCII file that contains a list of image
	 * file names. Each line of a LST file has the following format:
	 * reference_image_index  reference-image-filename comments
	 */
		
		
	class LstIO:public ImageIO
	{
	  public:
		LstIO(const string & filename, IOMode rw_mode = READ_ONLY);
		~LstIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);

		bool is_single_image_format() const
		{
			return false;
		}
		int get_nimg();
	  private:
		string filename;
		IOMode rw_mode;
		FILE *lst_file;

		bool is_big_endian;
		bool initialized;
		int nimg;

		ImageIO *imageio;
		string ref_filename;

		int last_lst_index;
		int last_ref_index;

		int calc_ref_image_index(int image_index);
		static const char *MAGIC;
	};

}


#endif	//eman__lstio_h__
