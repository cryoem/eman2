/**
 * $Id$
 */
#ifndef eman__icosio_h__
#define eman__icosio_h__ 1


#include "imageio.h"
#include <stdio.h>

namespace EMAN
{

    /* file format:
       1. header
       2. data: ny*nz of rows. Each row is: 
       int nx*sizeof(float), nx*float, int nx*sizeof(float)
    */


    class IcosIO : public ImageIO
    {
    public:
	IcosIO(string filename, IOMode rw_mode = READ_ONLY);
	~IcosIO();

	DEFINE_IMAGEIO_FUNC;
	static bool is_valid(const void *first_block);

    private:
	enum { STAMP = 72, STAMP1 = 72, STAMP2 = 20, STAMP3 = 20 };

	struct IcosHeader
	{
	    int stamp;		/* = 72 */
	    char title[72];	/* title of the map */
	    int stamp1;		/* = 72 */
	    int stamp2;		/* = 20 */
	    int nx;		/* number of rows */
	    int ny;		/* number of columnss */
	    int nz;		/* number of sections */
	    float min;
	    float max;
	    int stamp3;		/* = 20 */
	};

    private:
	string filename;
	IOMode rw_mode;
	IcosHeader icosh;
	FILE *icos_file;
	bool is_big_endian;
	bool initialized;
	bool is_new_file;
    };


}


#endif
