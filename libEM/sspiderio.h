/**
 * $Id$
 */
#ifndef eman__single_spider_io_h__
#define eman__single_spider_io_h__ 1

#include "spiderio.h"


namespace EMAN
{
	class SingleSpiderIO:public SpiderIO
	{
	  public:
		SingleSpiderIO(const string & filename, IOMode rw_mode = READ_ONLY);
		~SingleSpiderIO();

		int write_header(const Dict & dict, int image_index = 0, const Region* area = 0,
						 EMUtil::EMDataType filestoragetype = EMUtil::EM_FLOAT,
						 bool use_host_endian = true);
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
