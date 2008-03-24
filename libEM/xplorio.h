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

#ifndef eman__xplorio_h__
#define eman__xplorio_h__ 1

#include "imageio.h"

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
	 * A XPLOR file contains one  2D or 3D image.
	 */ 
	
	class XplorIO:public ImageIO
	{
	  public:
		explicit XplorIO(const string & filename, IOMode rw_mode = READ_ONLY);
		~XplorIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);

	  private:		
		string filename;
		IOMode rw_mode;
		FILE *xplor_file;

		bool is_big_endian;
		bool initialized;
		bool is_new_file;
		int nlines_in_header;
		
		int nx;
		int ny;
		int nz;

		float apix_x;
		float apix_y;
		float apix_z;

		float cell_alpha;
		float cell_beta;
		float cell_gama;
		
		static const string SECTION_MODE;
		static const int NFLOAT_PER_LINE;
		static const int INTEGER_SIZE;
		static const int FLOAT_SIZE;
		static const char * OUTFORMAT;
	};

}


#endif	//eman__xplorio_h__
