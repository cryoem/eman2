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
 
#ifndef eman__emio_h__
#define eman__emio_h__ 1

#include "imageio.h"

namespace EMAN
{
	/** EmIO defines I/O operations on EM image format.
     * EM image = header + data with (data = nx * ny * nz).
	 *
	 * An EM image file stores 1 single 2D or 3D image.
     */
	class EmIO:public ImageIO
	{
	  public:
		explicit EmIO(const string & filename, IOMode rw_mode = READ_ONLY);
		~EmIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block, off_t file_size = 0);
		static size_t get_mode_size(char data_type);
		static int get_machine_type();
		static int to_em_datatype(char t);

	  private:
		struct EMHeader
		{
			char machine;		// 0=OS-8, 1=VAX; 2=Convex; 3=SGI; 5=Mac; 6=PC
			char is_new_ver;	// OS-9 only
			char not_used1;
			char data_type;		// 1=byte, 2=short; 4=int; 5=float; 8=complex; 9=double
			int nx;
			int ny;
			int nz;
			char comment[80];
			int parameters[40];
			char username[20];
			char date[8];
			char userdata[228];
		};

		enum DataType
		{
			EM_EM_CHAR = 1,
			EM_EM_SHORT = 2,
			EM_EM_INT = 4,
			EM_EM_FLOAT = 5,
			EM_EM_COMPLEX = 8,
			EM_EM_DOUBLE = 9,
			EM_EM_UNKNOWN
		};

		enum MachineType
		{
			EM_OS8 = 0,
			EM_VAX = 1,
			EM_CONVEX = 2,
			EM_SGI = 3,
			EM_MAC = 5,
			EM_PC = 6,
			EM_UNKNOWN_MACHINE
		};

	  private:
		string filename;
		IOMode rw_mode;
		FILE *em_file;
		EMHeader emh;

		size_t mode_size;
		DataType mode;
		bool is_big_endian;
		bool initialized;
		bool is_new_file;
	};

}


#endif	//eman__emio_h__
