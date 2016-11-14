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

#ifndef eman__amiraio_h__
#define eman__amiraio_h__ 1

#include "imageio.h"

namespace EMAN
{

	/** Amira file = ASCII header + binary data.
	 * Its first line has some magic
	 * name to label it as an Amira image. The first few lines of the
	 * file is the ASCII header. Followed the header is the data in
	 * binary format. The data has nx x ny x nz pixels.
	 *
	 * An Amira file has only 1 2D or 3D image.
	 */
	class AmiraIO:public ImageIO
	{
	  public:
		AmiraIO(const string & filename, IOMode rw_mode = READ_ONLY);
		~AmiraIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);

	  private:
		string filename;
		IOMode rw_mode;
		FILE *amira_file;

		bool is_big_endian;
		bool initialized;
		EMUtil::EMDataType dt;
		int nx;
		int ny;
		int nz;
		float pixel;		// A/pixel (or voxel)
		float xorigin, yorigin, zorigin;	// the origin coordinates in Angstrom for the first pixel (0,0,0).

		static const char *MAGIC;
	};

}


#endif	//eman__amiraio_h__
