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

#ifndef eman__pifio_h__
#define eman__pifio_h__ 1

#include "imageio.h"

namespace EMAN
{
	/** PIF(Portable Image Format for EM Data) is an image format from
	 * Purdue University.
	 *
	 * A PIF file = file header + (image header + image data) + (image header + image data) ...
	 *
	 * A PIF file has a overall file header followed by n images. Each
	 * image has a header and data block.
	 *
	 * EMAN only supports homogeneous PIF file, which means all images
	 * n a PIF should have the same dimensions. We also assume the
	 * (nx,ny,nz) in PifFileHeader are equal to (nx,ny,nz) in each
	 * Image header.
	 */
	class PifIO:public ImageIO
	{
	  public:
		explicit PifIO(const string & filename, IOMode rw_mode = READ_ONLY);
		~PifIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);

		bool is_single_image_format() const
		{
			return false;
		}
		int get_nimg();
		
	  private:
		enum
		{
			PIF_MAGIC_NUM = 8
		};

		enum PifDataMode
		{
			PIF_CHAR = 0,
			PIF_SHORT = 1,
			PIF_FLOAT_INT = 2,
			PIF_SHORT_COMPLEX = 3,
			PIF_FLOAT_INT_COMPLEX = 4,
			PIF_BOXED_DATA = 6,			// byte, read as 0
			PIF_SHORT_FLOAT = 7,
			PIF_SHORT_FLOAT_COMPLEX = 8,
			PIF_FLOAT = 9,
			PIF_FLOAT_COMPLEX = 10,
			PIF_MAP_FLOAT_SHORT = 20,
			PIF_MAP_FLOAT_INT = 21,
			PIF_MAP_FLOAT_INT_2 = 40,			// 4 byte floatint, read as 2
			PIF_BOXED_FLOAT_INT = 46,		// 4 byte floatint, read as 2
			PIF_INVALID
		};

		// note there are no floats in these files. Floats are stored as ints with
		// a scaling factor.
		struct PifFileHeader
		{
			int magic[2];		// magic number; identify PIF file
			char scalefactor[16];	// to convert float ints -> floats
			int nimg;			// number of images in file
			int endian;			// endianness, 0 -> vax,intel (little), 1 -> SGI, PowerPC (big)
			char program[32];	// program which generated image
			int htype;			// 1 - all images same number of pixels and depth, 0 - otherwise
			int nx;				// number of columns
			int ny;				// number of rows
			int nz;				// number of sections
			int mode;			// image data type

			int even;
			int mrcX;			// MRC X dim
			int mrcY;			// MRC Y dim
			int mrcZ;			// MRC Z dim
			char scale_fac[16];
			char ctf_a0[13];
			char ctf_a1[13];
			char ctf_a2[13];
			char ctf_a3[13];
			char ctf_a4[13];
			char ctf_a5[13];
			char pad[318];
		};

		struct PifColorMap
		{						// color map for depthcued images
			short r[256];
			short g[256];
			short b[256];
		};

		struct PifImageHeader
		{
			int nx;
			int ny;
			int nz;				// size of this image
			int mode;			// image data type
			int bkg;			// background value
			int radius;			// boxed image radius
			int xstart;
			int ystart;
			int zstart;			// starting number of each axis
			int mx;
			int my;
			int mz;				// intervals along x,y,z
			int xlen;
			int ylen;
			int zlen;			// cell dimensions (floatints)
			int alpha;
			int beta;
			int gamma;			// angles (floatints)
			int mapc;
			int mapr;
			int maps;			// axes->sections (1,2,3=x,y,z)
			int min;
			int max;
			int mean;
			int sigma;			// statistics (floatints)
			int ispg;			// spacegroup
			int nsymbt;			// bytes for symmetry ops
			int xorigin;
			int yorigin;		// origin (floatint)
			char title[80];
			char time[32];
			char imagenum[16];	// unique micrograph number
			char scannum[8];	// scan number of micrograph
			int aoverb;
			int mapabang;
			int pad[63];		// later use
		};
	  
		int get_mode_size(PifDataMode mode);
		bool is_float_int(int mode);
		void fseek_to(int image_index);
		int to_em_datatype(int pif_datatype);
		int to_pif_datatype(int em_datatype);

	    string filename;
		IOMode rw_mode;
		PifFileHeader pfh;
		FILE *pif_file;
		int mode_size;
		bool is_big_endian;
		bool initialized;
		bool is_new_file;
		float real_scale_factor;
	};

}


#endif	//eman__pifio_h__
