#ifndef eman__single_spider_io_h__
#define eman__single_spider_io_h__ 1

#include "spiderio.h"


namespace EMAN
{
    class SingleSpiderIO : public SpiderIO
    {
    public:
	SingleSpiderIO(string filename, IOMode rw_mode = READ_ONLY);
	~SingleSpiderIO();

	int write_header(map<string, EMObject> & dict, int img_index = 0);
	int write_data(float *p_data, int img_index = 0);

	static bool is_valid(const void *first_block);

    protected:
	bool is_valid_spider(const void *first_block);
    };

}

#endif
