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

#ifndef eman__mrcio_h__
#define eman__mrcio_h__ 1

#include "imageio.h"

namespace EMAN
{
	/** MRC file = header + data (nx x ny x nz).
	 * A MRC image file stores 1D, 2D or 3D image. The image's
	 * dimensions and pixel type are defined in the header.
	 */
	
	class MrcIO:public ImageIO
	{
	public:
		explicit MrcIO(const string & filename, IOMode rw_mode = READ_ONLY);
		~MrcIO();

		DEFINE_IMAGEIO_FUNC;

		int read_ctf(Ctf & ctf, int image_index = 0);
		void write_ctf(const Ctf & ctf, int image_index = 0);

		static bool is_valid(const void *first_block, off_t file_size = 0);
		static int get_mode_size(int mm);
		static int to_em_datatype(int mrcmode);
		static int to_mrcmode(int em_datatype, int is_complex);

	private:
		enum MrcMode {
			MRC_UCHAR = 0,
			MRC_SHORT,
			MRC_FLOAT,
			MRC_SHORT_COMPLEX,
			MRC_FLOAT_COMPLEX,
			MRC_USHORT = 6,		//non-standard
			MRC_UCHAR3 = 16,	//unsigned char * 3, for rgb data, non-standard
			MRC_UNKNOWN
		};

		enum {
			MRC_NUM_LABELS = 10,
			MRC_LABEL_SIZE = 80,
			NUM_4BYTES_PRE_MAP = 52,
			NUM_4BYTES_AFTER_MAP = 3
		};

		/* updated to MRC Image2000 format which is compatible with CCP4 format */
		struct MrcHeader
		{
			int nx;				/* number of columns */
			int ny;				/* number of rows */
			int nz;				/* number of sections */

			int mode;			/* See modes above. */

			int nxstart;		/* No. of first column in map, default 0. */
			int nystart;		/* No. of first row in map, default 0. */
			int nzstart;		/* No. of first section in map,default 0. */

			int mx;				/* Number of intervals along X. */
			int my;				/* Number of intervals along Y. */
			int mz;				/* Number of intervals along Z. */

			/* Cell: treat a whole 2D image as a cell */
			float xlen;			/* Cell dimensions (Angstroms). */
			float ylen;			/* Cell dimensions (Angstroms). */
			float zlen;			/* Cell dimensions (Angstroms). */

			float alpha;		/* Cell angles (Degrees). */
			float beta;			/* Cell angles (Degrees). */
			float gamma;		/* Cell angles (Degrees). */

			/* axis X => 1, Y => 2, Z => 3 */
			int mapc;			/* Which axis corresponds to Columns.  */
			int mapr;			/* Which axis corresponds to Rows.     */
			int maps;			/* Which axis corresponds to Sections. */

			float amin;			/* Minimum density value. */
			float amax;			/* Maximum density value. */
			float amean;		/* Mean density value.    */

			int ispg;			/* Space group number (0 for images). */

			int nsymbt;			/* Number of chars used for storing symmetry operators. */

			int user[25];

			float xorigin;		/* X origin. */
			float yorigin;		/* Y origin. */
			float zorigin;		/* Z origin. */

			char map[4];		/* constant string "MAP "  */
			int machinestamp;	/* machine stamp in CCP4 convention:
								   big endian=0x11110000 little endian=0x44440000 */
						/* There is an ambiguity in the specification, using 0x11111111 & 4 instead */

			float rms;			/* rms deviation of map from mean density */

			int nlabels;		/* Number of labels being used. */
			char labels[MRC_NUM_LABELS][MRC_LABEL_SIZE];
		};

		static const char *CTF_MAGIC;
		static const char *SHORT_CTF_MAGIC;


	private:
		string filename;
		IOMode rw_mode;
		FILE *mrcfile;
		MrcHeader mrch;
		int mode_size;

		int is_ri;
		bool is_big_endian;
		bool is_new_file;
		bool initialized;
		
		/** generate the machine stamp used in MRC image format. */
		static int generate_machine_stamp();
		void swap_header(MrcHeader& mrch);
		
		/** This is a utility function used to calculate new min/max/mean/sigma
		 * when write MRC file as 16 bit or 8 bit. It will write new set of 
		 * min/max/mean/sigma to mrch.
		 * this function needs get the output data storage type from mrch.mode.*/
		void update_stat(void* data);
	};
}

#endif	//eman__mrcio_h__
