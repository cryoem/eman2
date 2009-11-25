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

#ifndef eman__xyzio_h__
#define eman__xyzio_h__1

#include "imageio.h"

namespace EMAN
{
	/** XYZIO is a sample Image IO class. It defines all required API
     * that you may need to implement.
     */

	class XYZIO:public ImageIO
	{
	  public:
		XYZIO(const string & filename, IOMode rw_mode = READ_ONLY);
		~XYZIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);
		int get_nimg();

	  private:
		string filename;
		IOMode rw_mode;
		FILE *xyz_file;

		bool is_big_endian;
		bool initialized;
	};

}


#endif
