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
		explicit Gatan2IO(const string & filename, IOMode rw_mode = READ_ONLY);
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
