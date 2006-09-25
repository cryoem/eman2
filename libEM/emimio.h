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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * 
 * */
 
#ifndef eman__emimio_h__
#define eman__emimio_h__ 1

#include "imageio.h"

namespace EMAN
{
	/** EMIM image format = 1 EmimFileHeader + n x (EmimImageHeader + data).
	 * EmimFileHeader defines the number of images 'n' in the
	 * file. All following images have the same dimensions:nx x ny x nz.
	 */
	
	class EmimIO:public ImageIO
	{
	  public:
		EmimIO(const string & filename, IOMode rw_mode = READ_ONLY);
		~EmimIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);

		
		bool is_single_image_format() const
		{
			return false;
		}
		int get_nimg();
	  private:
		static const char *MAGIC;

		enum
		{
			EMIM_COMPLEX = 1,
			PARM_NONE = 0,
			PARM_NORM = 1,
			NUM_INT_IN_FILE_HEADER = 16
		};

		struct EmimFileHeader
		{
			char magic[4];		// 'EMIM'
			int order;			// This defines the byte order, value is 2 on correct machine
			int count;			// # images in the file (24 bits max)
			int nx;
			int ny;
			int nz;				// image size (all images in 1 file have same size)
			int flag;			// flags are the same as well
			float pixel;		// pixel/voxel size in A
			int misc[7];		// misc usage
			int headlen;		// length of individual image headers (currently 166 x 4)
			char name[64];		// optional description of entire set of images
		};

		struct EmimImageHeader
		{
			char name[80];		// data set name
			unsigned int time;	// time at which image was acquired
			int mgnum;			// micrograph number or CCD frame
			int clipnum;		// id of clip within micrograph/frame
			int id;				// another id number
			float dx;
			float dy;
			float dz;			// translation
			float alt, az, phi;	// orientation
			short nparm;		// number of parameter sets
			int ptype[4];		// type of each parameter set
			float param[4][32];	// Parameterized CTF/SNR, etc...
		};

	  private:
		string filename;
		IOMode rw_mode;
		EmimFileHeader efh;
		FILE *emim_file;
		bool is_big_endian;
		bool initialized;
	};
}


#endif	//eman__emimio_h__
