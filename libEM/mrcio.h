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

		int get_nimg();

	private:
		enum MrcMode {
			MRC_UCHAR = 0,
			MRC_SHORT = 1,
			MRC_FLOAT = 2,
			MRC_SHORT_COMPLEX=3,
			MRC_FLOAT_COMPLEX=4,
			MRC_USHORT = 6,
			MRC_UCHAR3 = 16,	// unsigned char * 3, for rgb data, non-standard
			MRC_CHAR = 17,		// non-standard - signed char
			MRC_UHEX = 101,	// 2 4-bit values packed into each 8-bit byte
			MRC_UNKNOWN = 18
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
			int nx;				/* 0 - number of columns */
			int ny;				/* 1 - number of rows */
			int nz;				/* 2 - number of sections */

			int mode;			/* 3 - See modes above. */

			int nxstart;		/* 4 - No. of first column in map, default 0. */
			int nystart;		/* 5 - No. of first row in map, default 0. */
			int nzstart;		/* 6 - No. of first section in map,default 0. */

			int mx;				/* 7 - Number of intervals along X. */
			int my;				/* 8 - Number of intervals along Y. */
			int mz;				/* 9 - Number of intervals along Z. */

			/* Cell: treat a whole 2D image as a cell */
			float xlen;			/* 10 - Cell dimensions (Angstroms). */
			float ylen;			/* 11 - Cell dimensions (Angstroms). */
			float zlen;			/* 12 - Cell dimensions (Angstroms). */

			float alpha;		/* 13 - Cell angles (Degrees). */
			float beta;			/* 14 - Cell angles (Degrees). */
			float gamma;		/* 15 - Cell angles (Degrees). */

			/* axis X => 1, Y => 2, Z => 3 */
			int mapc;			/* 16 - Which axis corresponds to Columns.  */
			int mapr;			/* 17 - Which axis corresponds to Rows.     */
			int maps;			/* 18 - Which axis corresponds to Sections. */

			float amin;			/* 19 - Minimum density value. */
			float amax;			/* 20 - Maximum density value. */
			float amean;		/* 21 - Mean density value.    */

			int ispg;			/* 22 - Space group number (0 for images). */

			int nsymbt;			/* 23 - Number of chars used for storing symmetry operators. */

			int user1[15];			// 24 - 38
			int imod_flags;			/* 39 - bit flags used by IMOD - >= 16 for 8 bit packed */
			int user2[9];

			float xorigin;		/* 40 - X origin. */
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

		/** Extended MRC format for tomography
		 * As used by Fei; original definition of extended header by Dave Agard and Bram Koster
		 * Contact Dustin Morado <Dustin.Morado@uth.tmc.edu> for details.
		 *
		 * The extended MRC format consists of three parts:
		 * 1. the header (same as in the original MRC format definition)
		 * 2. the extended header
		 * 3. the data
		 * The MRC file is in little-endian (PC) format (IMOD can handle this). If you need to swap
		 * to big-endian (UNIX or Mac), all data must be swapped according to their data format.
		 *
		 * All image data will be 2-byte integer. */
		struct FeiMrcHeader
		{
			int nx;				/* The number of pixels in the x direction of the image */
			int ny;				/* The number of pixels in the y direction of the image */
			int nz;				/* The number of pixels in the z direction of the image
			 	 	 	 	 	 Effectively in the tomography tilt series this means the
			 	 	 	 	 	 number of images in the tilt series. */

			int mode;			/* Defines the data type. Should always be 1 (2-byte integer)
			 	 	 	 	 	 in the case of tomography*/

			int nxstart;		/* set to 0: not used; lower bound of columns */
			int nystart;		/* set to 0: not used; lower bound of rows */
			int nzstart;		/* set to 0: not used; lower bound of sections */

			int mx;				/* set to nx: not used; grid size x */
			int my;				/* set to ny: not used; grid size y */
			int mz;				/* set to nz: not used; grid size z */

			float xlen;			/* set to mx: not used; cell size in x Angstroms (pixel spacing=3Dxlen/mx) */
			float ylen;			/* set to mx: not used; cell size in y Angstroms (pixel spacing=3Dxlen/my) */
			float zlen;			/* set to mx: not used; cell size in z Angstroms (pixel spacing=3Dxlen/mz) */

			float alpha;		/* set to 90: not used; cell angles in degrees */
			float beta;			/* set to 90: not used; cell angles in degrees */
			float gamma;		/* set to 90: not used; cell angles in degrees */

			/* axis X => 1, Y => 2, Z => 3 */
			int mapc;			/* set to 1: not used; mapping columns, rows, sections on axis (x=3D1, y=3D2, z=3D3) */
			int mapr;			/* set to 2: not used; mapping columns, rows, sections on axis (x=3D1, y=3D2, z=3D3)     */
			int maps;			/* set to 3: not used; mapping columns, rows, sections on axis (x=3D1, y=3D2, z=3D3) */

			float amin;			/* minimum pixel value of all images in file */
			float amax;			/* maximum pixel value of all images in file */
			float amean;		/* mean pixel value of all images in file */

			short ispg;			/* set to 0: not used; space group number (0 for images) */

			short nsymbt;			/* set to 0: not used; number of bytes used for storing symmetry operators */

			int	next;			/* This value gives the offset (in bytes) from the end
								of the file header to the first dataset (image).
								Thus you will find the first image at 1024 + next bytes. */
			short dvid;			/* set to 0: not used; creator id */
			char extra[30];		/* set to 0: not used, extra 30 bytes data */
			short numintegers;	/* set to 0: not used */
			short numfloats;	/* set to 32; we always expect a extended header of 32 floats */

			short sub;
			short zfac;
			float min2;
			float max2;
			float min3;
			float max3;
			float min4;
			float max4;
			short idtype;
			short lens;
			short nd1;
			short nd2;
			short vd1;
			short vd2;
			float tiltangles[9];	/* set to 0; not used; used to rotate model to match rotated image */

			float zorg;			/* set to 0: not used; origin of image */
			float xorg;
			float yorg;

			int nlabl;			/* number of labels */
			char labl[MRC_NUM_LABELS][MRC_LABEL_SIZE]; 	/* Arrays of characters that can be used for description.
			 	 	 	 	 	 	 	 	 	 	 	 	 Label0 is used for copyright information, always start with "Fei"*/
		};

		/** The extended header used by Fei MRC image. It contains the information about a maximum of 1024 images.
		 * Each section is 128 bytes long. The extended header is thus 1024*128 bytes (always the same length,
		 * reagrdless of how many images are present.)
		 * Not always 1024*128 bytes, but at least 128*nz. the length of extended header is defined in the regular header field "next" */
		struct FeiMrcExtHeader
		{
			float a_tilt;		/* Alpha tilt, in degrees */
			float b_tilt;		/* beta tilt, in degrees */

			float x_stage;		/* Stage x position. Normally in SI units (meters), but some older files may be
			 	 	 	 	 	 in micrometers. Check by looking at values for x,y,z. If one of these exceeds 1,
			 	 	 	 	 	 it will be micrometers. */
			float y_stage;		/* Stage y position. For testing of units see x_stage */
			float z_stage;		/* Stage z position. For testing of units see x_stage */

			float x_shift;		/* Image shift x. For testing of units see x_stage */
			float y_shift;		/* Image shift y. For testing of units see x_stage */

			float defocus;		/* Defocus as read from microscope. For testing of units see x_stage */
			float exp_time;		/* Exposure time in seconds */
			float mean_int;		/* Mean value of image */

			float tilt_axis;	/* The orientation of the tilt axis in the image in degrees.
			 	 	 	 	 	 Vertical to the top is 0=B0, the direction of positive rotation is anti-clockwise */
			float pixel_size;	/* The pixel size of the images in SI units (meters) */
			float magnification;	/*The magnification used in SI units (volts) */
			float ht;			/* Value of the high tension in SI units (volts) */
			float binning;		/* The binning of the CCD or STEM acquisition */
			float appliedDefocus;	/* The intended application defocus in SI units (meters),
			 	 	 	 	 	 as defined for example in the tomography parameters view. */

			float remainder[16];	/* not used */
		};

		static const char *CTF_MAGIC;
		static const char *SHORT_CTF_MAGIC;

	private:
		string filename;
		IOMode rw_mode;
		FILE *mrcfile;
		int mode_size;

		union {
			MrcHeader mrch;
			FeiMrcHeader feimrch;
		};

		/* the extended MRC format for tomography, used by FEI */
		bool isFEI;

		/* for MRCS (MRC stack) format */
		bool is_stack;
		int stack_size;

		/* for MRC 8 bit packed format (2 4-bit values in each byte) */
		bool is_8_bit_packed;
		bool use_given_dimensions;

		int is_ri;
		bool is_big_endian;
		bool is_new_file;
		bool initialized;
		bool is_transpose;
		float rendermin;
		float rendermax;
		
		/** generate the machine stamp used in MRC image format. */
		static int generate_machine_stamp();
		void swap_header(MrcHeader& mrch);
		
		/** This is a utility function used to calculate new min/max/mean/sigma
		 * when write MRC file as 16 bit or 8 bit. It will write new set of 
		 * min/max/mean/sigma to mrch.
		 * this function needs get the output data storage type from mrch.mode.*/
		void update_stats(void* data, size_t size);

		/** This is a utility routine to tell whether to byte swap MRC header. */
		void check_swap(int * data, char * filnam, bool show_errors,
							 bool & do_swap, bool & have_err);

		int read_mrc_header(Dict & dict, int image_index = 0, const Region * area = 0, bool is_3d = false);
		int read_fei_header(Dict & dict, int image_index = 0, const Region * area = 0, bool is_3d = false);

		//utility funciton to tranpose x and y dimension in case the source mrc image is mapc=2,mapr=1
		int transpose(float *data, int nx, int ny, int nz) const;
	};
}

#endif	//eman__mrcio_h__
