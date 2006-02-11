/**
 * $Id$
 */
#include "spiderio.h"
#include "geometry.h"
#include "portable_fileio.h"
#include "Assert.h"

using namespace EMAN;

SpiderIO::SpiderIO(const string & spider_filename, IOMode rw)
:	filename(spider_filename), rw_mode(rw)
{
	spider_file = 0;
	first_h = 0;
	cur_h = 0;

	is_big_endian = ByteOrder::is_host_big_endian();
	initialized = false;
	is_new_file = false;
}

SpiderIO::~SpiderIO()
{
	if (spider_file) {
		fclose(spider_file);
		spider_file = 0;
	}

	if (first_h) {
		free(first_h);
		first_h = 0;
	}
	if (cur_h) {
		free(cur_h);
		cur_h = 0;
	}

}

void SpiderIO::init()
{
	if (initialized) {
		return;
	}
	ENTERFUNC;
	initialized = true;
	spider_file = sfopen(filename, rw_mode, &is_new_file);

	if (!is_new_file) {
		first_h = static_cast < SpiderHeader * >(calloc(1, sizeof(SpiderHeader)));

		if (fread(first_h, sizeof(SpiderHeader), 1, spider_file) != 1) {
			throw ImageReadException(filename, "SPIDER header");
		}

		if (!is_valid_spider(first_h)) {
			throw ImageReadException(filename, "invalid SPIDER");
		}

		float nslice = first_h->nslice;

		is_big_endian = ByteOrder::is_float_big_endian(nslice);
		become_host_endian((float *) first_h, NUM_FLOATS_IN_HEADER);
		
		if (first_h->istack == SINGLE_IMAGE_HEADER && rw_mode == WRITE_ONLY) {
			fclose(spider_file);
			spider_file = 0;
			
			spider_file = fopen(filename.c_str(), "wb");
		}
	}

	EXITFUNC;
}

bool SpiderIO::is_valid_spider(const void *first_block)
{
	return SpiderIO::is_valid(first_block);
}

bool SpiderIO::is_valid(const void *first_block)
{
	ENTERFUNC;
	bool result = false;
	
	if (!first_block) {
		 result = false;
	}
	else {
		const float *data = static_cast < const float *>(first_block);
		float nslice = data[0];
		float type = data[4];
		float ny = data[1];
		float istack = data[23];

		bool bigendian = ByteOrder::is_float_big_endian(nslice);
				
		if (bigendian != ByteOrder::is_host_big_endian()) {
			ByteOrder::swap_bytes(&nslice);
			ByteOrder::swap_bytes(&type);
			ByteOrder::swap_bytes(&ny);
			ByteOrder::swap_bytes(&istack);
		}

		if ((int (nslice)) !=nslice) {
			 result = false;
		}
		else {
			const int max_dim = 1 << 20;
			int itype = static_cast < int >(type);
			
			if ((itype == IMAGE_2D || itype == IMAGE_3D) &&
				(istack == STACK_OVERALL_HEADER) &&
				(nslice > 0 && nslice < max_dim) && (ny > 0 && ny < max_dim)) {
				 result = true;
			}
		}
	}

	EXITFUNC;
	return result;
}

/** If image_index = -1, read the overall spider header.
 */
int SpiderIO::read_header(Dict & dict, int image_index, const Region * area, bool)
{
	ENTERFUNC;

    bool is_read_overall_header = false;
    if (image_index == -1) {
        is_read_overall_header = true;
        image_index = 0;
    }
    
	check_read_access(image_index);
	
	if (!first_h) {
		throw ImageReadException(filename, "empty spider header");
	}
	
	check_region(area, FloatSize(first_h->nx, first_h->ny,first_h->nslice), 
				 is_new_file);
	
	float size = first_h->nx * first_h->ny * first_h->nslice;
	int overall_headlen = 0;

    if (!is_read_overall_header) {
        if (first_h->istack == STACK_OVERALL_HEADER) {
            overall_headlen = (int) first_h->headlen;
        }
        else {
            if(image_index != 0) {
                char desc[1024];
                sprintf(desc, "image index must be 0. Your image index = %d.", image_index);
                throw ImageReadException(filename, desc);
            }
        }
    }
    
	size_t single_image_size = (size_t) (first_h->headlen + size * sizeof(float));
	size_t offset = overall_headlen + single_image_size * image_index;

	SpiderHeader *cur_image_hed;
	if (offset == 0) {
		cur_image_hed = first_h;
	}
	else {
		cur_image_hed = static_cast < SpiderHeader * >(calloc(1, sizeof(SpiderHeader)));
		portable_fseek(spider_file, offset, SEEK_SET);

		if (fread(cur_image_hed, sizeof(SpiderHeader), 1, spider_file) != 1) {
			char desc[1024];
			sprintf(desc, "read spider header with image_index = %d failed", image_index);
			throw ImageReadException(filename, desc);
		}

		become_host_endian((float *) cur_image_hed, NUM_FLOATS_IN_HEADER);

		if (cur_image_hed->nx != first_h->nx || cur_image_hed->ny != first_h->ny
			|| cur_image_hed->nslice != first_h->nslice) {
			char desc[1024];
			sprintf(desc, "%dth image size %dx%dx%d != overall size %dx%dx%d",
					image_index, (int)cur_image_hed->nx, (int)cur_image_hed->ny,
					(int)cur_image_hed->nslice, 
					(int)first_h->nx, (int)first_h->ny, (int)first_h->nslice);
			throw ImageReadException(filename, desc);
		}
	}


	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, (int) cur_image_hed->nx, &xlen, (int) cur_image_hed->ny,
							&ylen, (int) cur_image_hed->nslice, &zlen);

	dict["nx"] = xlen;
	dict["ny"] = ylen;
	dict["nz"] = zlen;

	dict["datatype"] = EMUtil::EM_FLOAT;

	dict["minimum"] = cur_image_hed->min;
	dict["maximum"] = cur_image_hed->max;
	dict["mean"] = cur_image_hed->mean;
	dict["sigma"] = cur_image_hed->sigma;

	dict["SPIDER.nslice"] = (int) cur_image_hed->nslice;
	dict["SPIDER.type"] = (int) cur_image_hed->type;
	dict["SPIDER.headrec"] = (int) cur_image_hed->headrec;
	dict["SPIDER.angvalid"] = cur_image_hed->angvalid;

	dict["SPIDER.headlen"] = (int) cur_image_hed->headlen;
	dict["SPIDER.dx"] = cur_image_hed->dx;
	dict["SPIDER.dy"] = cur_image_hed->dy;
	dict["SPIDER.dz"] = cur_image_hed->dz;

	dict["SPIDER.reclen"] = (int) cur_image_hed->reclen;
	dict["SPIDER.istack"] = (int) cur_image_hed->istack;
	dict["SPIDER.date"] = cur_image_hed->date;
	dict["SPIDER.name"] = cur_image_hed->name;
	dict["SPIDER.maxim"] = (int)cur_image_hed->maxim;
	dict["SPIDER.imgnum"] = (int)cur_image_hed->imgnum;
	// complex check is in is_complex_mode, but even/odd check needs to be here
	int type = dict["SPIDER.type"];
	if ( type == IMAGE_2D_FFT_ODD
			|| type == IMAGE_3D_FFT_ODD ) 
		dict["is_fftodd"] = 1;
	else dict["is_fftodd"] = 0;
	// complex spider files are always ri, so we just set it unconditionally
	dict["is_complex_ri"] = 1;

	if (offset != 0) {
		if( cur_image_hed )
		{
			free(cur_image_hed);
			cur_image_hed = 0;
		}
	}
	EXITFUNC;
	return 0;
}


int SpiderIO::write_header(const Dict & dict, int image_index, const Region* area,
						   EMUtil::EMDataType, bool)
{
	ENTERFUNC;
	
	check_write_access(rw_mode, image_index);
	if (area) {
		if (!cur_h) {
			if (image_index == 0) {
				cur_h = static_cast < SpiderHeader * >(calloc(1, sizeof(SpiderHeader)));
				memcpy(cur_h, first_h, sizeof(SpiderHeader));
			}
			else {
				cur_h = static_cast < SpiderHeader * >(calloc(1, sizeof(SpiderHeader)));
				float size = first_h->nx * first_h->ny * first_h->nslice;
				size_t single_image_size = (int) (first_h->headlen + size * sizeof(float));
				size_t offset = (int) first_h->headlen + single_image_size * image_index;
				portable_fseek(spider_file, offset, SEEK_SET);

				if (fread(cur_h, sizeof(SpiderHeader), 1, spider_file) != 1) {
					char desc[1024];
					sprintf(desc, "read spider header with image_index = %d failed",
							image_index);
					throw ImageReadException(filename, desc);
				}
			}
		}
		
		check_region(area, FloatSize(cur_h->nx, cur_h->ny,cur_h->nslice), 
					 is_new_file);
		EXITFUNC;
		return 0;
	}
	
	int nx1 = dict["nx"];
	int ny1 = dict["ny"];
	int nz1 = dict["nz"];

	size_t header_size = sizeof(SpiderHeader);
	size_t record_size = nx1 * sizeof(float);
	size_t num_records = header_size / record_size;

	if ((header_size % record_size) != 0) {
		num_records++;
	}

	size_t header_length = num_records * record_size;
	size_t irec = nx1 * nz1 + num_records;

	if (!is_new_file) {
		if (first_h->istack != STACK_OVERALL_HEADER) {
			throw ImageWriteException(filename,
									  "cannot write to single-image Spider files");
		}

        // if there is only 1 image in the file and image_index = 0,
        // file overwriting is permitted for different-size image.
        if ( ((int) first_h->maxim <= 1) && (image_index == 0)) {
            first_h->nx = (float)nx1;
            first_h->ny = (float)ny1;
            first_h->nslice = (float)nz1;
        }
        else {        
            if (nx1 != first_h->nx || ny1 != first_h->ny || nz1 != first_h->nslice) {
                char desc[256];
                sprintf(desc, "new size %dx%dx%d != existing size %dx%dx%d",
                        nx1, ny1, nz1, (int)first_h->nx, 
                        (int)first_h->ny, (int)first_h->nslice);
                throw ImageWriteException(filename, desc);
            }
        }
        
		int old_maxim = static_cast < int >(first_h->maxim);

		if (image_index >= first_h->maxim) {
			first_h->maxim = (float) (image_index + 1);
		}

		first_h->max = dict["maximum"];
		first_h->min = dict["minimum"];
		first_h->mean = dict["mean"];
		first_h->sigma = dict["sigma"];

		if (need_swap()) {
			ByteOrder::swap_bytes((float *) first_h, NUM_FLOATS_IN_HEADER);
		}

		rewind(spider_file);

		if (fwrite(first_h, header_size, 1, spider_file) != 1) {
			throw ImageWriteException(filename, "cannot write to Spider file");
		}

		if (need_swap()) {
			ByteOrder::swap_bytes((float *) first_h, NUM_FLOATS_IN_HEADER);
		}

		first_h->maxim = (float) old_maxim;
		portable_fseek(spider_file, 0, SEEK_END);
	}
	else {
		if (!first_h) {
			first_h = (SpiderHeader *) calloc(1, header_length);
			first_h->mmvalid = 0;
			first_h->sigma = -1.0;
			first_h->angvalid = 0;
			first_h->scale = 1.0;
			first_h->nslice = (float) nz1;
			first_h->nx = (float) nx1;
			first_h->ny = (float) ny1;

			if (1 == int(dict["is_complex"])) {
				if (1 == int(dict["is_fftodd"])) {
					if (1 == nz1) {
						first_h->type = IMAGE_2D_FFT_ODD;
					} else {
						first_h->type = IMAGE_3D_FFT_ODD;
					}
				} else {
					if (1 == nz1) {
						first_h->type = IMAGE_2D_FFT_EVEN;
					} else {
						first_h->type = IMAGE_3D_FFT_EVEN;
					}
				}
			} else {
				if (nz1 == 1) {
					first_h->type = IMAGE_2D;
				}
				else {
					first_h->type = IMAGE_3D;
				}
			}

			first_h->reclen = (float)record_size;
			first_h->headrec = (float)num_records;
			first_h->headlen = (float)header_length;
			first_h->inuse = -1;
			first_h->istack = STACK_OVERALL_HEADER;
			first_h->maxim = 1;
			first_h->u1 = (float)irec;

			first_h->max = dict["maximum"];
			first_h->min = dict["minimum"];
			first_h->mean = dict["mean"];
			first_h->sigma = dict["sigma"];
		}
		else {
			first_h->maxim++;
		}

		rewind(spider_file);

		if (fwrite(first_h, header_length, 1, spider_file) != 1) {
			throw ImageWriteException(filename,
									  "cannot write a new header to Spider file");
		}
	}

	if (!cur_h) {
		cur_h = (SpiderHeader *) calloc(1, header_size);
	}

	cur_h->scale = 1.0;
	cur_h->nslice = (float)nz1;
	cur_h->nx = (float)nx1;
	cur_h->ny = (float)ny1;

	if (1 == int(dict["is_complex"])) {
		if (1 == int(dict["is_fftodd"])) {
			if (1 == nz1) {
				cur_h->type = IMAGE_2D_FFT_ODD;
			} else {
				cur_h->type = IMAGE_3D_FFT_ODD;
			}
		} else {
			if (1 == nz1) {
				cur_h->type = IMAGE_2D_FFT_EVEN;
			} else {
				cur_h->type = IMAGE_3D_FFT_EVEN;
			}
		}
	} else {
		if (nz1 == 1) {
			cur_h->type = IMAGE_2D;
		}
		else {
			cur_h->type = IMAGE_3D;
		}
	}

	cur_h->mmvalid = 1;

	cur_h->max = dict["maximum"];
	cur_h->min = dict["minimum"];
	cur_h->mean = dict["mean"];
	cur_h->sigma = dict["sigma"];

	cur_h->reclen = static_cast < float >(nx1 * sizeof(float));
	cur_h->inuse = 1.0;
	cur_h->istack = STACK_IMAGE_HEADER;
	cur_h->imgnum = (float)(image_index + 1);
	cur_h->headrec = (float)num_records;
	cur_h->headlen = (float)header_length;
	cur_h->angvalid = 1.0;
	cur_h->u1 = (float)irec;

	size_t data_size = (size_t)(cur_h->nx) * (size_t)(cur_h->ny) * (size_t)(cur_h->nslice);
	off_t cur_offset = (off_t) (cur_h->headlen +
								(data_size*sizeof(float) + cur_h->headlen) * image_index);
	portable_fseek(spider_file, cur_offset, SEEK_SET);
	fwrite(cur_h, header_size, 1, spider_file);
	
	int pad_size =  (int)cur_h->headlen - header_size;
	if (pad_size > 0) {
		char *pad = static_cast < char *>(calloc(pad_size, 1));
		fwrite(pad, pad_size, 1, spider_file);
        if( pad )
        {
        	free(pad);
        	pad = 0;
        }
	}
	
	EXITFUNC;
	return 0;

}

int SpiderIO::write_single_header(const Dict & dict, const Region *area,
								  EMUtil::EMDataType, bool)
{
	ENTERFUNC;

	check_write_access(rw_mode, 0);

	if (area) {
		check_region(area, FloatSize(first_h->nx, first_h->ny, first_h->nslice),
					 is_new_file);
		EXITFUNC;
		return 0;
	}
	
	int nx = dict["nx"];
	int ny = dict["ny"];
	int nz = dict["nz"];

	size_t header_size = sizeof(SpiderHeader);
	size_t record_size = nx * sizeof(float);
	size_t num_records = header_size / record_size;

	if ((header_size % record_size) != 0) {
		num_records++;
	}

	size_t header_length = num_records * record_size;

	if (!first_h) {
		first_h = static_cast < SpiderHeader * >(calloc(1, header_size));
	}

	first_h->mmvalid = 0;
	first_h->sigma = -1.0;
	first_h->angvalid = 0;
	first_h->scale = 1.0;
	first_h->istack = SINGLE_IMAGE_HEADER;
	first_h->nslice = (float)nz;
	first_h->nx = (float)nx;
	first_h->ny = (float)ny;

	first_h->max = dict["maximum"];
	first_h->min = dict["minimum"];
	first_h->mean = dict["mean"];
	first_h->sigma = dict["sigma"];


	if (1 == int(dict["is_complex"])) {
		if (1 == int(dict["is_fftodd"])) {
			if (1 == nz) {
				first_h->type = IMAGE_2D_FFT_ODD;
			} else {
				first_h->type = IMAGE_3D_FFT_ODD;
			}
		} else {
			if (1 == nz) {
				first_h->type = IMAGE_2D_FFT_EVEN;
			} else {
				first_h->type = IMAGE_3D_FFT_EVEN;
			}
		}
	} else {
		if (nz == 1) {
			first_h->type = IMAGE_2D;
		}
		else {
			first_h->type = IMAGE_3D;
		}
	}

	first_h->reclen = (float)record_size;
	first_h->headrec = (float)num_records;
	first_h->headlen = (float)header_length;
	first_h->inuse = 1;

	if (fwrite(first_h, header_size, 1, spider_file) != 1) {
		throw ImageWriteException(filename, "write single spider header failed");
	}

	size_t pad_size = header_length - header_size;
	char *pad = static_cast < char *>(calloc(pad_size, 1));
	fwrite(pad, pad_size, 1, spider_file);
	if( pad )
	{
		free(pad);
		pad = 0;
	}
	
	EXITFUNC;
	return 0;
}


int SpiderIO::read_data(float *data, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	check_read_access(image_index, data);

	check_region(area, FloatSize((int) first_h->nx, (int) first_h->ny,
								 (int) first_h->nslice), is_new_file);

	int overall_headlen = 0;
	if (first_h->istack == STACK_OVERALL_HEADER) {
		overall_headlen = (int) first_h->headlen;
	}
	else {
		if(image_index != 0) {
			char desc[256];
			sprintf(desc, "image index must be 0. Your image index = %d.", image_index);
			throw ImageReadException(filename, desc);
		}
	}

	size_t size = static_cast < size_t > (first_h->nx * first_h->ny * first_h->nslice);
	size_t single_image_size = static_cast < size_t > (first_h->headlen +
													   size * sizeof(float));
	off_t offset = overall_headlen + single_image_size * image_index;
	portable_fseek(spider_file, offset, SEEK_SET);
        
	portable_fseek(spider_file, (int) first_h->headlen, SEEK_CUR);
   
#if 1
	EMUtil::process_region_io(data, spider_file, READ_ONLY, 0, sizeof(float), 
							  (int) first_h->nx, (int) first_h->ny, 
							  (int) first_h->nslice, area);
#endif
#if 0
	unsigned int nz = static_cast < unsigned int >(first_h->nslice);
	int sec_size = static_cast < int >(first_h->nx * first_h->ny * sizeof(float));

	if (fread(data, sec_size, nz, spider_file) != nz) {
		LOGERR("Incomplete SPIDER data read");
		return 1;
	}
#endif

	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, (int) first_h->nx, &xlen, (int) first_h->ny, &ylen,
							(int) first_h->nslice, &zlen);

	int data_size = xlen * ylen * zlen;
	become_host_endian(data, data_size);
	EXITFUNC;
	return 0;
}

int SpiderIO::write_data(float *data, int image_index, const Region* area,
						 EMUtil::EMDataType, bool)
{
	ENTERFUNC;

	check_write_access(rw_mode, image_index, 0, data);
	if (area) {
		check_region(area, FloatSize(cur_h->nx, cur_h->ny, cur_h->nslice),
					 is_new_file);
	}
	
	if (!cur_h) {
		throw ImageWriteException(filename, "please write header before write data");
	}

	if (cur_h->istack == SINGLE_IMAGE_HEADER) {
		throw ImageWriteException(filename,
								  "cannot mix single spider and stack spider.");
	}

	int hl = static_cast < int >(cur_h->headlen);
	int head_size = sizeof(SpiderHeader);
	int pad_size = hl - head_size;

	char *pad = static_cast < char *>(calloc(pad_size, 1));
	size_t data_size = static_cast <size_t>(cur_h->nx * cur_h->ny * cur_h->nslice);

	int sec_size = static_cast < int >(cur_h->nx * cur_h->ny * sizeof(float));
	int nz = static_cast < int >(cur_h->nslice);

	swap_header(cur_h);
	swap_data(data, data_size);
	
	while (image_index > first_h->maxim) {
		first_h->maxim++;
		first_h->inuse = 1.0;
		first_h->imgnum = first_h->maxim;

		portable_fseek(spider_file, 0, SEEK_END);
		fwrite(cur_h, head_size, 1, spider_file);
		fwrite(pad, pad_size, 1, spider_file);
		fwrite(data, sec_size, nz, spider_file);
	}

	off_t cur_offset = (off_t) (hl + (data_size * sizeof(float) + hl) * image_index + hl);
	portable_fseek(spider_file, cur_offset, SEEK_SET);

#if 1
    // with region I/O support
	EMUtil::process_region_io(data, spider_file, WRITE_ONLY, image_index,
							  sizeof(float), (int)cur_h->nx, (int)cur_h->ny,
							  (int)cur_h->nslice, area);
#endif
#if 0
    // no region I/O support. This part should be removed if the
    // region I/O version works.
	if (fwrite(data, sec_size, nz, spider_file) != (unsigned int) nz) {
		throw ImageWriteException(filename, "spider writing failed");
	}
#endif
	swap_header(cur_h);
	swap_data(data, data_size);

	if( pad )
	{
		free(pad);
		pad = 0;
	}
	EXITFUNC;
	return 0;
}


void SpiderIO::flush()
{
	fflush(spider_file);
}

int SpiderIO::write_single_data(float *data, const Region * area, EMUtil::EMDataType, bool)
{
	ENTERFUNC;
	
	check_write_access(rw_mode, 0, 1, data);

	if (area) {
		check_region(area, FloatSize(first_h->nx, first_h->ny, first_h->nslice),
					 is_new_file);
	}
	
	if (!first_h) {
		throw ImageWriteException(filename, "NULL image header");
	}
	
	if (first_h->istack != SINGLE_IMAGE_HEADER) {
		throw ImageWriteException(filename, "no mix of single spider with stack spider.");
	}

	portable_fseek(spider_file, (int) first_h->headlen, SEEK_SET);

	// remove the stuff in #if 0 ... #endif if the following works.
	EMUtil::process_region_io(data, spider_file, WRITE_ONLY,0, sizeof(float), 
							  (int) first_h->nx, (int) first_h->ny, 
							  (int) first_h->nslice, area);
	
#if 0
	size_t sec_size = (size_t)first_h->nx * (size_t)first_h->ny * sizeof(float));
	size_t nz = static_cast < size_t >(first_h->nslice);

	if (fwrite(data, sec_size, nz, spider_file) != nz) {
		throw ImageWriteException(filename, "single spider image writing failed");
	}
#endif

	EXITFUNC;
	return 0;
}

void SpiderIO::swap_data(float *data, size_t size)
{
	if (data && need_swap()) {
		ByteOrder::swap_bytes(data, size);
	}
}

void SpiderIO::swap_header(SpiderHeader * header)
{
	if (header && need_swap()) {
		ByteOrder::swap_bytes((float *) header, NUM_FLOATS_IN_HEADER);
	}
}

bool SpiderIO::is_complex_mode()
{
	int type = static_cast<int>(first_h->type);
	if (type == IMAGE_2D_FFT_ODD
		|| type == IMAGE_2D_FFT_EVEN
		|| type == IMAGE_3D_FFT_ODD
		|| type == IMAGE_3D_FFT_EVEN
	   ) return true;
	return false;
}

bool SpiderIO::is_image_big_endian()
{
	init();
	return is_big_endian;
}

int SpiderIO::get_nimg()
{
	init();
	if (!first_h) {
		Assert(is_new_file);
		return 0;
	}
	
	if (first_h->istack == STACK_OVERALL_HEADER) {
		return static_cast < int >(first_h->maxim);
	}

	return 1;
}

bool SpiderIO::need_swap() const
{
	if (!is_new_file && (is_big_endian != ByteOrder::is_host_big_endian())) {
		return true;
	}
	return false;
}
/* vim: set ts=4 noet: */
