#ifndef __lstio_h__
#define __lstio_h__


#include "imageio.h"
#include <stdio.h>

namespace EMAN {

    class LstIO: public ImageIO {
    public:
	LstIO(string filename, IOMode rw_mode = READ_ONLY);
	~LstIO();
	
	DEFINE_IMAGEIO_FUNC;
	
    private:
	int calc_ref_image_index(int image_index);
	static const char* MAGIC;
    private:
	string filename;
	IOMode rw_mode;
	FILE* lst_file;
	
	bool is_big_endian;
	bool initialized;
	int nimg;

	ImageIO* imageio;
	string ref_filename;
	
	int last_lst_index;
	int last_ref_index;
    };

}


#endif
