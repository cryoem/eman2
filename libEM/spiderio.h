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

#ifndef eman__spiderio_h__
#define eman__spiderio_h__ 1

#include "imageio.h"

namespace EMAN
{
	/** SPIDER: (System for Processing Image Data from Electron microscopy 
	 * and Related fields) is an image processing system for electron
	 * microscopy. SpiderIO reads/writes images used in spider format.
	 * (reference: http://www.wadsworth.org/spider_doc/spider/docs/index.html)
     *
	 * A single image = header + data.
	 *
	 * header-length = ( ceiling(1024/(nx*4)) * (nx*4) )
	 * where nx is the number of pixels per row. 
	 * The data is (nx * ny * nslice) of floating numbers,
	 * where ny is number of rows per slice. nslice is number of slices.
	 *
	 * There are 2 spider image formats:
	 *
	 *  - single-spider.
	 *       It is simple 2D/3D non-stack file. image = header + data.
	 *
	 *  - image stack.
	 *        There is one overall stack header followed by a stack 
	 *        of images in which each image has its own image header.
	 *        EMAN currently only supports homogeneous image stack,
	 *        which means all images have the same sizes. 
     *        
     *        if there is only 1 image in the file, it can be overwritten
     *        by a different-size image.
     *
	 * Record: In spider terminology, each row is called a record.
	 *
	 * Note: To read the overall image header in a stacked spider
	 * file, use image_index = -1.
     *
	 */
	class SpiderIO:public ImageIO
	{
	  public:
		explicit SpiderIO(const string & filename, IOMode rw_mode = READ_ONLY);
		~SpiderIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);
		
		bool is_single_image_format() const { return false;	}
		
		/** get the number of images in this stacked SPIDER image
		 * @return the number of images
		 * */
		int get_nimg();

	protected:
		struct SpiderHeader
		{
			float nslice;		// number of slices in volume; 1 for a 2D image.
			float nrow;         // nrow, number of rows per slice 
			float irec;			// total number of records in the file (unused)
			float nhistrec;     // obsolete, unused
			
			/** file type: 
			 *   1	: 2D image; 
			 *   3	: 3D volume;
			 * -11	: 2D Fourier, mixed radix odd. 
			 * -12	: 2D Fourier, mixed radix even.
			 * -21	: 3D Fourier, mixed radix odd.
			 * -22	: 3D Fourier, mixed radix even.
			 */
			float type;	        //iform, file type
			
			/** max/min flag. Is set at 0 when the file is created,
			 * and at 1 when the maximum, minimum, average, and
			 * standard deviation have been computed and stored into this
			 * header.
			 */
			float mmvalid;		// imami, max/min flag.
			float max;          // fmax, max value
			float min;          // fmin, min value
			float mean;         // av, average value
			float sigma;		// sig, std dev, -1=unknown
			float ihist;        // obsolete, no longer used
			float nsam;         // nsam, number of pixels per row
			float headrec;		// labrec, number of records in header
			float angvalid;		// iangle, flag, 1 if tilt angles have been computed
			float phi;          // tilt angle
			float theta;        // tilt angle
			float gamma;        // tilt angle  (also called psi).
			float dx;           // xoff, x translation
			float dy;           // yoff, y translation
			float dz;			// zoff, z translation
			float scale;        // scale factor
			float headlen;		// labbyt, header length in bytes
			float reclen;       // lenbyt, record length in bytes

			/** istack = 0 for simple 2D or 3D (non-stack) files.
			 * In an "image stack" there is one overall stack header followed by 
			 * a stack of images in which each image has its own image header. 
			 * (An image stack differs from a simple 3D image in that each 
			 * stacked image has its own header.) A value >0 in this position 
			 * in the overall stack header indicates a stack of images. 
			 * A value of <0 inthis position in the overall stack header 
			 * indicates an indexed stack of images and gives the maximum image
			 * number allowed in the index.
			 * for stacked image, istack=2 in overall header, istack =-1
			 * in following individual images.
			*/
			float istack;
			float inuse;        // not used
			/** maxim is only used in the overall header for a stacked
			 * image file. It is the number of the highest image
			 * currently used in the stack. The number is updated, if necessary,
			 * when an image is added or deleted from the stack.
			 */
			float maxim;
			/** imgnum is only used in a stacked image header. It is the number 
			 *	of the current image or zero if the image is unused.
			 */
			float imgnum;
			
			/**This position is only used in the overall header of indexed stacks.
			 * There, this position is the highest index currently in use.*/
			float lastindx;    
			
			float u6;		   // unused
			float u7;          // unused
			
			/** flag that additional angles are present in header. 
			 * 1 = one additional rotation is present, 
			 * 2 = additional rotation that preceeds the rotation that was stored in
			 * words 15..20.
			 */
			float Kangle;	
			float phi1;
			float theta1;
			float psi1;
			float phi2;
			float theta2;
			float psi2;
			char  u8[48];		//unused
			float xf[27];		// reserved for Jose Maria's transforms
			float u9[135];      // unused
			char date[11];      // creation date e.g. 27-MAY-1999 
			char time[8];       // creation time e.g. 09:43:19 
			char title[160];    
		};

		enum SpiderType
		{
			IMAGE_2D 			=   1,
			IMAGE_3D 			=   3,
			IMAGE_2D_FFT_ODD 	= -11,
			IMAGE_2D_FFT_EVEN 	= -12,
			IMAGE_3D_FFT_ODD 	= -21,
			IMAGE_3D_FFT_EVEN 	= -22
		};

		enum
		{
			SINGLE_IMAGE_HEADER		=   0,
			OVERALL_STACK_HEADER	= 	2,
			NUM_FLOATS_IN_HEADER 	= 211
		};
		
		/** write a SPIDER header to spider_file
		 * @param dict the dictionary contain header information
		 * @param area the region we want to write
		 * @param image_index the image index inthe stack, it's 0-indexed
		 * @param offset the offset in the spider_file
		 * @param hp the SpiderHeader pointer all header info will write to, then the content of hp will write to spider file
		 * @param ISTACK the istack field for SPIDER's header, 0 for dingle image, 2 for stacked image
		 * @param MAXIM maxim field for header, only used for overall header
		 * @param IMGNUM imgnum for header, only used for stacked image header
		 * @param use_host_endian use host machine's endian, set to false will swap the byte order
		 * @exception ImageWriteException
		 * @return 0 for sucess
		 * */
		//ISTACK, MAXIM, IMGNUM only useful for overall stack header
		int write_single_header(const Dict & dict, const Region* area, int image_index, size_t offset,
								SpiderHeader *& hp,	int ISTACK, int MAXIM=1, int IMGNUM=1, bool use_host_endian=true);	
		
		/** write a single image data
		 * @param data the data block to be written
		 * @param area the data region we want to write
		 * @param hp the SpiderHeader pointer all header info will write to
		 * @param offset the offset in the spider_file
		 * @param img_index the image index inthe stack, it's 0-indexed
		 * @param max_nimg max image number in a stack
		 * @param use_host_endian use host machine's endian, set to false will swap the byte order
		 * @exception ImageWriteException
		 * @return 0 for success
		 * */
		int write_single_data(float *data, const Region * area, SpiderHeader *& hp,
								size_t offset, int img_index, int max_nimg, bool use_host_endian=true);
		
		/**check the data block to see if it represents valid stacked SPIDER image file header
		 * @param first_block the pointer to first block of the file
		 * @return boolean result of the check, true for valid stacked SPIDER image 
		 * */					  
		virtual bool is_valid_spider(const void *first_block);

	  protected:
		bool need_swap() const;
		void swap_data(float *data, size_t nitems);
		void swap_header(SpiderHeader * header);

	  protected:
		string filename;
		IOMode rw_mode;

		FILE *spider_file;
		SpiderHeader *first_h; // overall image header
		SpiderHeader *cur_h;   // the current reading/writing image header
		
		bool is_big_endian;
		bool initialized;
		bool is_new_file;
	};
}

#endif	//eman__spiderio_h__
