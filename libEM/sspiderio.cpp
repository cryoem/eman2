/**
 * $Id$
 */
#include "sspiderio.h"
#include "geometry.h"
#include "portable_fileio.h"

#include <iostream>

using namespace EMAN;

SingleSpiderIO::SingleSpiderIO(const string & file, IOMode rw)
:	SpiderIO(file, rw)
{
}


SingleSpiderIO::~SingleSpiderIO()
{
}


int SingleSpiderIO::write_header(const Dict & dict, int , const Region* area,
								 EMUtil::EMDataType filestoragetype, bool use_host_endian)
{
	size_t offset = 0;
	int image_index = 0;
	return write_single_header(dict, area, image_index, offset, first_h, SINGLE_IMAGE_HEADER);
}


int SingleSpiderIO::write_data(float *data, int , const Region* area,
							   EMUtil::EMDataType filestoragetype, bool use_host_endian)
{
	size_t offset = (int) first_h->headlen;
	return write_single_data(data, area, first_h, offset, 0, 1);
}


bool SingleSpiderIO::is_valid(const void *first_block)
{
	ENTERFUNC;
	bool result = false;
	
	if (first_block) {
		const float *data = static_cast < const float *>(first_block);
		float nslice = data[0];
		float nrow = data[1];
		float iform = data[4];
		float nsam = data[11];
		float labrec = data[12];	//NO. of records in file header
		float labbyt = data[21];	//total NO. of bytes in header
		float lenbyt = data[22];	//record length in bytes
		float istack = data[23];
		
		bool big_endian = ByteOrder::is_float_big_endian(nslice);
		if (big_endian != ByteOrder::is_host_big_endian()) {
			ByteOrder::swap_bytes(&nslice);
			ByteOrder::swap_bytes(&nrow);
			ByteOrder::swap_bytes(&iform);
			ByteOrder::swap_bytes(&nsam);
			ByteOrder::swap_bytes(&labrec);
			ByteOrder::swap_bytes(&labbyt);
			ByteOrder::swap_bytes(&lenbyt);
			ByteOrder::swap_bytes(&istack);
		}
		
		if( int(nslice) != nslice || int(nrow) != nrow
				|| int(iform) != iform || int(nsam) != nsam 
				|| int(labrec) != labrec || int(labbyt) != labbyt
				|| int(lenbyt) != lenbyt ) {
			result =  false;
		}
		else {
			int itype = static_cast < int >(iform);
			if( int(istack) != SINGLE_IMAGE_HEADER ) {
				result = false; //istack>0 for overall header, istack<0 for indexed stack of image
			}
			else if( itype == IMAGE_2D_FFT_ODD || itype == IMAGE_2D_FFT_EVEN 
					|| itype == IMAGE_3D_FFT_ODD || itype == IMAGE_3D_FFT_EVEN ) {
				result = false;	//Complex SPIDER image not supported in EMAN2
			}
			else if (itype == IMAGE_2D || itype == IMAGE_3D) {
				result = true;
			}
		}
		
		int ilabrec = static_cast<int>(labrec);
		int ilabbyt = static_cast<int>(labbyt);
		int ilenbyt = static_cast<int>(lenbyt);
		if( ilabbyt != ilabrec * ilenbyt ) {
			result = false;
		}
	}

	EXITFUNC;
	return result;
}


bool SingleSpiderIO::is_valid_spider(const void *first_block)
{
	return SingleSpiderIO::is_valid(first_block);
}
