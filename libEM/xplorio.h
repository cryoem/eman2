/**
 * $Id$
 */
#ifndef eman__xplorio_h__
#define eman__xplorio_h__ 1


#include "imageio.h"
#include <stdio.h>

namespace EMAN
{
	/** XPLOR image format is in ASCII:
	 *
	 * 1. header:  (note: each integer takes 8 Bytes space, each float
	 *               is 12.5E format.)
	 *    line 1: empty
	 *    line 2: int: number of lines for title (ntitle)
	 *    next ntitle lines: string: titles
	 *    line ntitle+3: 9 int: nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax
	 *    line ntitle+4: 6 float: cell size (x, y, z), cell angles (alpha, beta, gamma)
	 *    line ntitle+5: string: ZYX (this is the section mode. always ZYX)
	 *    
	 *
	 * 2. data
	 *
	 *    for zindex = 1 to nz:
	 *      zindex
	 *      nx*ny floats. each line has 6 floats. each float is in 12.5E format.
	 *
	 * A XPLOR file contains 1 2D or 3D image.
	 */ 
	
	class XplorIO:public ImageIO
	{
	  public:
		XplorIO(string filename, IOMode rw_mode = READ_ONLY);
		~XplorIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);

	  private:
		string filename;
		IOMode rw_mode;
		FILE *xplor_file;

		bool is_big_endian;
		bool initialized;

		int nx;
		int ny;
		int nz;

		float apix_x;
		float apix_y;
		float apix_z;

		float cell_alpha;
		float cell_beta;
		float cell_gama;
		
		static string SECTION_MODE;
		static int NFLOAT_PER_LINE;
		static int INTEGER_SIZE;
		static int FLOAT_SIZE;
	};

}


#endif
