#ifndef eman__amiraio_h__
#define eman__amiraio_h__ 1


#include "imageio.h"
#include <stdio.h>

namespace EMAN {

    class AmiraIO: public ImageIO {
    public:
	AmiraIO(string filename, IOMode rw_mode = READ_ONLY);
	~AmiraIO();
	
	DEFINE_IMAGEIO_FUNC;

    private:
	static const char* MAGIC;
	
    private:
	string filename;
	IOMode rw_mode;
	FILE* amira_file;
	
	bool is_big_endian;
	bool initialized;
	int nx;
	int ny;
	int nz;
    };
    
}


#endif
