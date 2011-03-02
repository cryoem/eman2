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

#ifndef eman__imagicio2_h__
#define eman__imagicio2_h__ 1

#include "imageio.h"

namespace EMAN
{
	/** IMAGIC-5 Header File Format
	 * 
	 * Renewed 4D version: http://www.imagescience.de/formats/formats.htm
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
	 * This renewed version support storing multiple 3D images in one file
	 * (header/data pair).
	 * 
	 * EMAN2 will read both old and new Imagic header,
	 * but write to new format from - Grant Tang.
	 */

	class ImagicIO2:public ImageIO
	{
	  public:
		static const char *HED_EXT;
		static const char *IMG_EXT;
		
		explicit ImagicIO2(string filename, IOMode rw_mode = READ_ONLY);
		~ImagicIO2();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);
		
		 /* If returns -1, then we have an old-style Imagic file
		 *
		 * @return 0 for new Imagic4D file, -1 for old-style Imagic file*/
		int init_test();

		/**If this file format is only for single iamge
		 *
		 * @return false this file format support stack*/
		bool is_single_image_format() const
		{
			return false;
		}
		
		/** Get number of images in this file
		 *
		 * @return number of images*/
		int get_nimg();
		
	  private:
		static const char *REAL_TYPE_MAGIC;
		static const char *CTF_MAGIC;


		enum DataType
		{
			IMAGIC_CHAR,
			IMAGIC_SHORT,
			IMAGIC_FLOAT,
			IMAGIC_FLOAT_COMPLEX,
			IMAGIC_FFT_FLOAT_COMPLEX,
			IMAGIC_UNKNOWN_TYPE
		};

		/*This enum is used to skip those char type when swap bytes between big endian and little endian*/
		enum
		{
			NUM_4BYTES_PRE_IYLP = 14,
			NUM_4BYTES_AFTER_IXOLD = 14,
			NUM_4BYTES_AFTER_NAME = 150	//before HISTORY
		};
		
		enum RealType
		{
			VAX_VMS = 16777216,			//for VAX/VMS
			LINUX_WINDOWS = 33686018,	//for OSF, ULTRIX, LINUX, MS WINDOWS
			SGI_IBM = 67372036			//for SiliconGraphics, SUN, HP, IBM
		};

		/** IMAGIC-4D file format
		 * http://www.imagescience.de/formats/formats.htm
		 * */
		struct Imagic4D {
			int imgnum;	//1 image location number (1,2,3,...)
			int count;	//2 total # images in file (1st record only), total number of images - 1 (0,1,2,...)
			int error;	//3 error code for this image during IMAGIC run
			int headrec;//4 number of header records per image (=1 always)
			int month;	//5 creation month
			int mday;	//6	creation day
			int year;	//7 creation year
			int hour;	//8 creation hour
			int minute; //9	creation minute
			int sec;	//10 creation second
			int rsize;	//11 image size in bytes
			int izold;	//12 top left Z co-ordinate before THREED-CUT
			int ny;		//13 number of lines per image (for 1D data IXLP1=1)
			int nx;		//14  number of pixels per line (IYLP)
			char type[4];	//15 only 'REAL', or 'INTG' are implemented for EMAN
							// 4 characters determining the image type
							// REAL : REAL/float
							// INTG : INTEGER*2/short
							// PACK : PACK/byte
							// COMP : 2 REAL / 2 float
							// RECO : complex format with 0 in imaginary part)
			int ixold;	//16  top left X co-ordinate before CUT-IMAGE (boxing)
			int iyold;	//17  top left Y co-ordinate before CUT-IMAGE (boxing)
			float avdens;	//18 average density in image
			float sigma;	//19 standard deviation of density
			float user1;	//20 at user's own disposal
			float user2;	//21 at user's own disposal
			float densmax;	//22 highest density in image
			float densmin;	//23 minimal density in image
			int complex;	//24 label indicating that data is always complex
			float defocus1;	//25 defocus value 1
			float defocus2;	//26 defocus value 2
			float defangle;	//27 defocus angle
			float sinostart;	//28 start angle if image is a sinogram
			float sinoend;	//29 end angle if image is a sinogram
			char label[80];	//30-49 coded NAME/TITLE of the image (80 characters)
			float ccc3d;	//50 3D simularity criteria
			int ref3d;		//51 3D membership
			int mident;		//52 micrograph identification number
			int ezshift;	//53 equivalent shift in Z direction 
			int ealpha;		//54 equivalent Euler angle alpha
			int ebeta;		//55 equivalent Euler angle beta 
			int egamma;		//56 equivalent Euler angle gamma
			int unused1;	//57 currently not used
			int unused2;		//58 currently not used
			int nalisum;	//59 number of image summed
			int pgroup;		//60 point-group symmetry in international notation (622, for example) 
			int izlp;		//61 number of 2D planes in 3D data (for 1D/2D: IZLP=1)
			int i4lp;		//62 number of objects in file:
							// 1D (IXLP=1): number of 1D spectra
							// 2D (IZLP=1): number of 2D images
							// 3D (IZLP>1): number of 3D volumes
			int i5lp;		//63
			int i6lp;		//64
			float alpha;	//65 Euler angle alpha (3D and Angular Reconst.)
			float beta;		//66 Euler angle beta (3D and Angular Reconst.)
			float gamma;	//67 Euler angle gamma (3D and Angular Reconst.)
			int imavers;	//68 IMAGIC version, which created the file (yyyymmdd)
			int realtype;	//69  floating point type, machine stamp
							// 16777216 for VAX/VMS
							// 33686018 for OSF,ULTRIX, LINUX, MS Windows
							// 67372036 for SiliconGraphics, SUN, HP, IBM
			char buffer[120];		//70-99	 Variables that control the buffering
									// during read/write in IMAGIC-5 programs.
									// PLEASE DO NOT TOUCH !
			float angle;	//100 last rotation angle
			float voltage;	//101 acceleration voltage (kv)
			int spaberr;	//102 sperical aberration (mm)
			int pcoher;		//103 partial coherence
			float ccc;		//104 cross correlation peak hight
			float errar;	//105 error in angular reconstitution,  if -1.0: the file is a special file (FABOSA)
			float err3d;	//106 error in 3D reconstruction
			int ref;		//107 (multi-) reference number
			float classno;	//108 class number in MSA classification
			float locold;	//109 location number before CUT-IMAGE (boxing), or before copy in ANG-RECONST and EX-COPY
			float repqual;	//110 representation quality, used in MSA-RUN and MSA (eigen) filtering
			float zshift;	//111 last shift in Z direction
			float xshift;	//112 last shift in X direction
			float yshift;	//113 last shift in Y direction
			float numcls;	//114 number members in the class specified in CLASSNO, if this image represents a class average (class-sum image)
			float ovqual;	//115 overall quality of the class in CLASSNO
			float eangle;	//116 equivalent angle
			float exshift;	//117 equivalent shift in X direction
			float eyshift;	//118 equivalent shift in Y direction
			float cmtotvar;	//119 total variance in data matrix relative to center of mass (MSA calculations)
			float informat;	//120 Gauss norm / real*FT Space information of the data set
			int numeigen;	//121 number of eigen values in MSA
			int niactive;	//122 number of active images in MSA calculations
			float resolx;	//123 Angstrom per pixel/voxel in X direction,  if DAT1(105) = -1.0 (FABOSA): mm per pixel
			float resoly;	//124 Angstrom per pixel/voxel in Y direction
			float resolz;	//125 Angstrom per pixel/voxel in Z direction
			float alpha2;	//126 Euler angle alpha (from projection matching), Special FABOSA variables if DAT1(105) = -1.0
			float beta2;	//127 Euler angle beta (from projection matching), Special FABOSA variables if DAT1(105) = -1.0
			float gamma2;	//128 Euler angle gamma (from projection matching), Special FABOSA variables if DAT1(105) = -1.0
			float nmetric;	//129 Metric used in MSA calculations
			float actmsa;	//130 a flag indicating whether the "image" is active or not. Used during MSA calculations
			float coosmsa[69];	//131-199  co-ordinates of "image" along factorial axis
								// number 1 through 69 (maximum possible).
								// Used during MSA calculations.
//			float eigval;					//150 eigval, eigenvalues if the "images" represent eigenimages (eigenvalue #1 into loc#1 etc.)
			char history[228];	//220-256 coded history of image (228 characters)
		};

		size_t get_datatype_size(DataType t) const;
		int to_em_datatype(DataType t) const;
		void make_header_host_endian(Imagic4D & hed) const;
		void swap_header(Imagic4D & hed) const;
		DataType get_datatype_from_name(const char *name) const;

		/** the Ctf object is a EMAN1Ctf object. */
		Ctf * read_ctf(const Imagic4D& hed) const;
		void write_ctf(const Ctf * const ctf, int image_index = 0);
		
		int generate_machine_stamp() const;
	  private:
		string filename;
		string hed_filename;
		string img_filename;

		IOMode rw_mode;
		FILE *hed_file;
		FILE *img_file;

		Imagic4D imagich;
		bool is_big_endian;
		bool initialized;
		bool is_new_hed;
		bool is_new_img;

		DataType datatype;
		int nz;
	};

}


#endif	//eman__imagicio2_h__
