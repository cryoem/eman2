/**
 * $Id$
 */

/*
 * Author: Grant Tang, 08/06/2010 (gtang@bcm.edu)
 * Copyright (c) 2000-2010 Baylor College of Medicine
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

#ifndef __df3io_h__
#define __df3io_h__

#include "imageio.h"
#include <cstdio>

/** df3 file format (http://www.povray.org/documentation/view/3.6.1/374)
 * Header:
 * The df3 format consists of a 6 byte header of three 16-bit integers with
 * high order byte first. These three values give the x,y,z size of the data
 * in pixels (or more appropriately called voxels ).
 * Data:
 * The header is followed by x*y*z unsigned integer bytes of data with a
 * resolution of 8, 16 or 32 bit. The data are written with high order byte
 * first (big-endian). The resolution of the data is determined by the size
 * of the df3-file. That is, if the file is twice (minus header, of course)
 * as long as an 8 bit file then it is assumed to contain 16 bit ints and
 * if it is four times as long 32 bit ints.
 * */

namespace EMAN
{
	class Df3IO:public ImageIO
	{
	public:
		explicit Df3IO(const string & filename, IOMode rw_mode = READ_ONLY);
		~Df3IO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block, off_t file_size = 0);
	private:
		unsigned short nx, ny, nz;
		string filename;
		IOMode rw_mode;
		FILE *df3file;
		bool initialized;
		bool is_new_file;
	};
}


#endif	//__df3io_h__
