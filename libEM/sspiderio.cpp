/**
 * $Id$
 */
#include "sspiderio.h"
#include "log.h"
#include "geometry.h"
#include "portable_fileio.h"
#include <assert.h>

using namespace EMAN;

SingleSpiderIO::SingleSpiderIO(string file, IOMode rw)
    : SpiderIO(file, rw)
{
}

SingleSpiderIO::~SingleSpiderIO()
{
}


int SingleSpiderIO::write_header(const Dict & dict, int image_index)
{
    Log::logger()->log("SingleSpiderIO::write_header()");
    assert(image_index == 0);
    return write_single_header(dict);
}

int SingleSpiderIO::write_data(float *data, int image_index)
{
    Log::logger()->log("SingleSpiderIO::write_data()");
    assert(image_index == 0);
    return write_single_data(data);
}

bool SingleSpiderIO::is_valid(const void *first_block)
{
    Log::logger()->log("SingleSpiderIO::is_valid()");
    if (!first_block) {
	return false;
    }

    const float *data = static_cast<const float *>(first_block);
    int nslice = static_cast<int>(data[0]);
    int type = static_cast<int>(data[4]);
    int ny = static_cast<int>(data[1]);
    int istack = static_cast<int>(data[23]);

    bool data_big_endian = ByteOrder::is_data_big_endian(&nslice);

    if (data_big_endian != ByteOrder::is_machine_big_endian()) {
	ByteOrder::swap_bytes(&nslice);
	ByteOrder::swap_bytes(&type);
	ByteOrder::swap_bytes(&ny);
	ByteOrder::swap_bytes(&istack);
    }

    if ((int (nslice)) !=nslice) {
	return false;
    }

    const int max_dim = 1 << 20;
    int itype = static_cast<int>(type);

    if ((itype == IMAGE_2D || itype == IMAGE_3D) && (istack == SINGLE_IMAGE_HEADER) &&
	(nslice > 0 && nslice < max_dim) && (ny > 0 && ny < max_dim)) {
	return true;
    }

    return false;
}


bool SingleSpiderIO::is_valid_spider(const void *first_block)
{
    return SingleSpiderIO::is_valid(first_block);
}
