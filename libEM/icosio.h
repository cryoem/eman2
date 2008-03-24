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

#ifndef eman__icosio_h__
#define eman__icosio_h__ 1

#include "imageio.h"

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
		explicit IcosIO(const string & filename, IOMode rw_mode = READ_ONLY);
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


#endif	//eman__icosio_h__
