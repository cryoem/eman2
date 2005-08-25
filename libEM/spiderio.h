/**
 * $Id$
 */
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
		SpiderIO(const string & filename, IOMode rw_mode = READ_ONLY);
		~SpiderIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);
		
		bool is_single_image_format() const
		{
			return false;
		}

		int get_nimg();

	protected:
		struct SpiderHeader
		{
			float nslice;		// number of slices in volume; 1 for a 2D image.
			float ny;           // number of rows per slice
			float u1;			// unused
			float u2;           // unused
			
			/** file type: 1: 2d real; 3: 3D real;
			 * -11: 2d Fourier, mixed  radix odd. 
			 * -12: 2D Fourier, mixed radix even.
			 * -21: 3D Fourier, mixed radix odd.
			 * -22: 3D Fourier, mixed radix even.
			 */
			float type;	        // file type
			
			/** max/min flag. Is set at 0 when the file is created,
			 * and at 1 when the maximum, minimum, average, and
			 * standard deviation have been computed and stored into this
			 * header.
			 */
			float mmvalid;		// max/min flag.
			float max;          // max value
			float min;          // min value
			float mean;         // average value
			float sigma;		// std dev, -1=unknown
			float u3;           // unused
			float nx;           // number of pixels per row
			float headrec;		// number of records in header
			float angvalid;		// 1 if tilt angles have been computed
			float phi;          // tilt angle
			float theta;        // tilt angle
			float gamma;        // tilt angle  (also called psi).
			float dx;           // x translation
			float dy;           // y translation
			float dz;			// z translation
			float scale;        // scale factor
			float headlen;		// header length in bytes
			float reclen;       // record length in bytes

			/** istack = 0 for simple 2D or 3D (non-stack) files.
			 * for stacked image, istack=2 in overall header, istack =-1
			 * in following individual images.
			*/
			float istack;
			float inuse;        // not used
			/** maxim is only used in the overall header for a stacked
			 * image file. It is the number of the highest image
			 * currently used in the stack.
			 */
			float maxim;
			/** imgnum is only used in a stacked image header. It is the number 
			 *	of the current image or zero if the image is unused.
			 */
			float imgnum;
			
			float u6;          // unused
			float u7;          // unused
			/** flag that additional angles are present in header. 1 =
			 *	one additional rotation is present, 2 = additional
			 * rotation that preceeds the rotation that was stored in
			 * words 15..20.
			 */
			float ang2valid;	// additional angles
			float phi1;
			float theta1;
			float psi1;
			float phi2;
			float theta2;
			float psi2;
			float xf[14];		// reserved for Jose Maria's transforms
			float u8[161];      // unused
			char date[12];      // creation date e.g. 27-MAY-1999 
			char time[8];       // creation time e.g. 09:43:19 
			char name[160];     // title
		};

		enum SpiderType
		{
			IMAGE_2D = 1,
			IMAGE_3D = 3,
			IMAGE_2D_FFT_ODD = -11,
			IMAGE_2D_FFT_EVEN = -12,
			IMAGE_3D_FFT_ODD = -21,
			IMAGE_3D_FFT_EVEN = -21
		};

		enum
		{
			SINGLE_IMAGE_HEADER = 0,
			STACK_IMAGE_HEADER = -1,
			STACK_OVERALL_HEADER = 2,
			NUM_FLOATS_IN_HEADER = 211
		};

		int write_single_header(const Dict & dict, const Region* area,
								EMUtil::EMDataType filestoragetype, bool use_host_endian);
		int write_single_data(float *data, const Region * area,
							  EMUtil::EMDataType filestoragetype, bool use_host_endian);
		virtual bool is_valid_spider(const void *first_block);

	  private:
		bool need_swap() const;
		void swap_data(float *data, size_t nitems);
		void swap_header(SpiderHeader * header);

	  private:
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
