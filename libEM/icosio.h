/**
 * $Id$
 */
#ifndef eman__icosio_h__
#define eman__icosio_h__ 1


#include "imageio.h"
#include <stdio.h>

namespace EMAN
{

	/** ICOS file = header + data.
	 *  1. header, defined in IcosHeader.
	 *  2. data: ny*nz of rows. Each row = n1 + row-data + n2
	 *     where n1: an integer. n1 = nx*sizeof(float).
	 *           n2: an integer. n2 = n1.
	 *     row-data: nx numbers of float points.
	 *
	 * An Icos file stores 1 2D or 3D image.
	 * 
	 */
	class IcosIO:public ImageIO
	{
	  public:
		IcosIO(const string & filename, IOMode rw_mode = READ_ONLY);
		~IcosIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);

	  private:
		enum
		{ STAMP = 72, STAMP1 = 72, STAMP2 = 20, STAMP3 = 20 };

		struct IcosHeader
		{
			int stamp;			/* = 72 */
			char title[72];		/* title of the map */
			int stamp1;			/* = 72 */
			int stamp2;			/* = 20 */
			int nx;				/* number of rows */
			int ny;				/* number of columnss */
			int nz;				/* number of sections */
			float min;
			float max;
			int stamp3;			/* = 20 */
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
