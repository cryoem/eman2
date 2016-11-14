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

#ifndef eman__imagicio_h__
#define eman__imagicio_h__ 1

#include "imageio.h"

namespace EMAN
{
	/** IMAGIC-5 Header File Format
	 *
	 * An IMAGIC-5 file has 2 files:
	 * a) a header file with extension ".hed". It contains information
	 *    for every image.
	 * b) an image file with extension ".img". It contains raw data.
	 *
	 * The header file contains one (fixed-size) record per image
	 * stored. Every header record consists of 256 REAL/float
	 * for every image.
	 * 
	 * The image file contains only the raw data. Depending on the
	 * internal IMAGIC-5 format used, which can be REAL, INTG, PACK
	 * or COMP, the data is stored as REAL/float, INTEGER/int,
	 * INTEGER*1/byte or 2x REAL/float, respectively. The first pixel
	 * stored is the upper left one. The data is stored line
	 * by line, section by section, volume by volume.
	 * 
	 * 3D imagic uses the same format to 2D. it is a bunch of 2D slices.
	 * use the 'hint' IS_3D to treat "2D slices" as 3D volume.
	 * 
	 * imagic doesn't store multiple 3D images in one file
	 * (header/data pair).
	 * 
	 */

	class ImagicIO:public ImageIO
	{
	  public:
		static const char *HED_EXT;
		static const char *IMG_EXT;
		
		explicit ImagicIO(string filename, IOMode rw_mode = READ_ONLY);
		~ImagicIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);
		
		bool is_single_image_format() const
		{
			return false;
		}
		
		int get_nimg();
		
	  private:
		static const char *REAL_TYPE_MAGIC;
		static const char *CTF_MAGIC;


		enum DataType
		{
			IMAGIC_UCHAR,
			IMAGIC_USHORT,
			IMAGIC_FLOAT,
			IMAGIC_FLOAT_COMPLEX,
			IMAGIC_FFT_FLOAT_COMPLEX,
			IMAGIC_UNKNOWN_TYPE
		};

		enum
		{
			NUM_4BYTES_PRE_IXOLD = 14,
			NUM_4BYTES_AFTER_IXOLD = 14,
			NUM_4BYTES_AFTER_SPACE = 207
		};

		struct ImagicHeader
		{
			int imgnum;			// image number, [1,n]
			int count;			// total number of images - 1 (only first image), [0,n-1]
			int error;			// Error code for this image
			int headrec;		// # of header records/image (always 1)
			int mday;			// image creation time
			int month;
			int year;
			int hour;
			int minute;
			int sec;
			int reals;			// image size in reals
			int pixels;			// image size in pixels
			int ny;				// # of lines / image
			int nx;				// # of pixels / line
			char type[4];		// PACK, INTG, REAL, COMP, RECO
			int ixold;			// Top left X-coord. in image before windowing 
			int iyold;			// Top left Y-coord. in image before windowing 
			float avdens;		// average density
			float sigma;		// deviation of density
			float varia;		// variance of density
			float oldav;		// old average density
			float max;			// max density
			float min;			// min density
			int complex;		// not used
			float cellx;		// not used
			float celly;		// not used
			float cellz;		// not used
			float cella1;		// not used
			float cella2;		// not used
			char label[80];		// image id string
			int space[8];
			float mrc1[4];
			int mrc2;
			int space2[7];
			int lbuf;			// effective buffer len = nx
			int inn;			// lines in buffer = 1
			int iblp;			// buffer lines/image = ny
			int ifb;			// 1st line in buf = 0
			int lbr;			// last buf line read = -1
			int lbw;			// last buf line written = 0
			int lastlr;			// last line called for read = -1
			int lastlw;			// last line called for write = 1
			int ncflag;			// decode to complex = 0
			int num;			// file number = 40 (?)
			int nhalf;			// leff/2
			int ibsd;			// record size for r/w (words) = nx*2
			int ihfl;			// file # = 8
			int lcbr;			// lin count read buf = -1
			int lcbw;			// lin count wr buf = 1
			int imstr;			// calc stat on rd = -1
			int imstw;			// calc stat on wr = -1
			int istart;			// begin line in buf = 1
			int iend;			// end line in buf = nx
			int leff;			// eff line len = nx
			int linbuf;			// line len (16 bit) nx *2
			int ntotbuf;		// total buf in pgm = -1
			int space3[5];
			int icstart;		// complex line start = 1
			int icend;			// complex line end = nx/2
			int rdonly;			// read only = 0
			int misc[157];		// Remainder of header (EMAN1 specific settings not supported)
		};

		size_t get_datatype_size(DataType t);
		int to_em_datatype(DataType t);
		void make_header_host_endian(ImagicHeader & hed);
		void swap_header(ImagicHeader & hed);
		DataType get_datatype_from_name(const char *name);

		/** the Ctf object is a EMAN1Ctf object. */
		Ctf * read_ctf(const ImagicHeader& hed) const;
		void write_ctf(const Ctf * const ctf, int image_index = 0);
		
	  private:
		string filename;
		string hed_filename;
		string img_filename;

		IOMode rw_mode;
		FILE *hed_file;
		FILE *img_file;

		ImagicHeader imagich;
		bool is_big_endian;
		bool initialized;
		bool is_new_hed;
		bool is_new_img;

		DataType datatype;
		int nz;
	};

}


#endif	//eman__imagicio_h__
