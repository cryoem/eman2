/**
 * $Id$
 */
#include "sspiderio.h"
#include "log.h"
#include "geometry.h"
#include "portable_fileio.h"


using namespace EMAN;

SingleSpiderIO::SingleSpiderIO(const string & file, IOMode rw)
:	SpiderIO(file, rw)
{
}

SingleSpiderIO::~SingleSpiderIO()
{
}


int SingleSpiderIO::write_header(const Dict & dict, int , const Region* area,
								 bool use_host_endian)
{
	return write_single_header(dict, area, use_host_endian);
}

int SingleSpiderIO::write_data(float *data, int , const Region* area,
							   bool use_host_endian)
{
	return write_single_data(data, area, use_host_endian);
}

bool SingleSpiderIO::is_valid(const void *first_block)
{
	ENTERFUNC;
	bool result = false;
	
	if (!first_block) {
		result = false;
	}

	const float *data = static_cast < const float *>(first_block);
	int nslice = static_cast < int >(data[0]);
	int type = static_cast < int >(data[4]);
	int ny = static_cast < int >(data[1]);
	int istack = static_cast < int >(data[23]);

	bool data_big_endian = ByteOrder::is_data_big_endian(&nslice);

	if (data_big_endian != ByteOrder::is_host_big_endian()) {
		ByteOrder::swap_bytes(&nslice);
		ByteOrder::swap_bytes(&type);
		ByteOrder::swap_bytes(&ny);
		ByteOrder::swap_bytes(&istack);
	}

	if ((int (nslice)) !=nslice) {
		result = false;
	}

	const int max_dim = 1 << 20;
	int itype = static_cast < int >(type);

	if ((itype == IMAGE_2D || itype == IMAGE_3D) && (istack == SINGLE_IMAGE_HEADER) &&
		(nslice > 0 && nslice < max_dim) && (ny > 0 && ny < max_dim)) {
		result = true;
	}

	EXITFUNC;
	return result;
}


bool SingleSpiderIO::is_valid_spider(const void *first_block)
{
	return SingleSpiderIO::is_valid(first_block);
}
