/**
 * $Id$
 */

/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
 * Copyright (c) 2000-2006 Baylor College of Medicine
 * 
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 * 
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * 
 * */

#ifndef eman__pngio_h__
#define eman__pngio_h__ 1

#ifdef EM_PNG

#include <png.h>
#include "imageio.h"


namespace EMAN
{
	/** PngIO reads/writes a 2D PNG image. Currently 8-bit and 16-bit
	 * PNG read/write are supported.
	 */
	class PngIO:public ImageIO
	{
	  public:
		explicit PngIO(const string & filename, IOMode rw_mode = READ_ONLY);
		~PngIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);

	  private:
		enum
		{
			PNG_BYTES_TO_CHECK = 8
		};

		enum BitDepthType
		{
			PNG_CHAR_DEPTH,
			PNG_SHORT_DEPTH,
			PNG_INVALID_DEPTH
		};

	  private:
		string filename;
		IOMode rw_mode;
		FILE *png_file;

		bool initialized;

		png_structp png_ptr;
		png_infop info_ptr;
		png_infop end_info;

		png_uint_32 nx;
		png_uint_32 ny;
		BitDepthType depth_type;
		int number_passes;
		
		float rendermin;
		float rendermax;
	};

}

#endif	//EM_PNG

#endif	//eman__pngio_h__
