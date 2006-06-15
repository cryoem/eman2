/**
 * $Id$
 */
#ifndef eman__single_spider_io_h__
#define eman__single_spider_io_h__ 1

#include "spiderio.h"


namespace EMAN
{
	/** Single Spider Image I/O class. For Spider and Single Spider
	 * image format, please refer spiderio.h.
	 * @see spiderio.h
	 */	
	class SingleSpiderIO:public SpiderIO
	{
	  public:
		/** SingleSpiderIO constructor.
		 * @param filename The filename of a single spider image file.
		 * @param rw_mode Read/Write file open mode.
		 */
		SingleSpiderIO(const string & filename, IOMode rw_mode = READ_ONLY);
		~SingleSpiderIO();
		
		/** Write a single SPIDER format header to an image.
		 *
		 * @param dict A keyed-dictionary storing the header information.
		 * @param image_index The index of the image to write.
		 * @param area The region to write data to.
		 * @param filestoragetype The image data type used in the output file.
		 * @param use_host_endian Whether to use the host machine
		 *        endian to write out or not. If false, use the
		 *        endian opposite to the host machine's endian.
		 * @return 0 if OK; 1 if error.
		 */
		int write_header(const Dict & dict, int image_index = 0, const Region* area = 0,
						 EMUtil::EMDataType filestoragetype = EMUtil::EM_FLOAT,
						 bool use_host_endian = true);
		
		/** Write data to an image.
		 *
		 * @param data An array storing the data.
		 * @param image_index The index of the image to write.
		 * @param area The region to write data to.
		 * @param filestoragetype The image data type used in the output file.
		 * @param use_host_endian Whether to use the host machine
		 *        endian to write out or not. If false, use the
		 *        endian opposite to the host machine's endian.
		 * @return 0 if OK; 1 if error.
		 */
		int write_data(float *data, int image_index = 0, const Region* area = 0,
					   EMUtil::EMDataType filestoragetype = EMUtil::EM_FLOAT,
					   bool use_host_endian = true);

		static bool is_valid(const void *first_block);
		
		bool is_single_image_format() const
		{
			return true;
		}
		
	  protected:
		  bool is_valid_spider(const void *first_block);
	};

}

#endif
