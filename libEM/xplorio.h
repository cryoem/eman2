#ifndef __xplorio_h__
#define __xplorio_h__


#include "imageio.h"
#include <stdio.h>

namespace EMAN {

    class XplorIO: public ImageIO {
    public:
	XplorIO(string filename, IOMode rw_mode = READ_ONLY);
	~XplorIO();
	
	DEFINE_IMAGEIO_FUNC;

    private:
	string filename;
	IOMode rw_mode;
	FILE* xplor_file;
	
	bool is_big_endian;
	bool initialized;

	int nx;
	int ny;
	int nz;
    };
    
}


#endif
