/**
 * $Id$
 */
#include "spiderio.h"
#include "log.h"
#include "geometry.h"
#include "portable_fileio.h"
#include "emutil.h"
#include <assert.h>

using namespace EMAN;


SpiderIO::SpiderIO(string spider_filename, IOMode rw)
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

int SpiderIO::init()
{
	ENTERFUNC;

	static int err = 0;
	if (initialized) {
		return err;
	}

	initialized = true;

	spider_file = sfopen(filename, rw_mode, &is_new_file);
	if (!spider_file) {
		err = 1;
		return err;
	}

	if (!is_new_file) {
		first_h = static_cast < SpiderHeader * >(calloc(sizeof(SpiderHeader), 1));

		if (fread(first_h, sizeof(SpiderHeader), 1, spider_file) != 1) {
			LOGERR("read header failed on file: '%s'", filename.c_str());
			err = 1;
			return err;
		}

		if (!is_valid_spider(first_h)) {
			LOGERR("'%s' is not a valid Spider file", filename.c_str());
			err = 1;
			return err;
		}
		int nslice = static_cast < int >(first_h->nslice);
		is_big_endian = ByteOrder::is_data_big_endian(&nslice);
		become_host_endian((float *) first_h, NUM_FLOATS_IN_HEADER);

		if (first_h->istack == SINGLE_IMAGE_HEADER && rw_mode != READ_ONLY) {
			fclose(spider_file);
			spider_file = 0;
			free(first_h);
			first_h = 0;
			spider_file = fopen(filename.c_str(), "wb");
		}
	}


	return 0;
}

bool SpiderIO::is_valid_spider(const void *first_block)
{
	return SpiderIO::is_valid(first_block);
}

bool SpiderIO::is_valid(const void *first_block)
{
	ENTERFUNC;

	if (!first_block) {
		return false;
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
		return false;
	}

	const int max_dim = 1 << 20;
	int itype = static_cast < int >(type);

	if ((itype == IMAGE_2D || itype == IMAGE_3D) && (istack == STACK_OVERALL_HEADER) &&
		(nslice > 0 && nslice < max_dim) && (ny > 0 && ny < max_dim)) {
		return true;
	}

	return false;
}

int SpiderIO::read_header(Dict & dict, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	if (check_read_access(image_index) != 0) {
		return 1;
	}

	if (!first_h) {
		LOGDEBUG("empty header. nothing to read");
		return 1;
	}

	if (check_region(area, IntSize((int) first_h->nx, (int) first_h->ny,
								   (int) first_h->nslice)) != 0) {
		return 1;
	}

	float size = first_h->nx * first_h->ny * first_h->nslice;
	int overall_headlen = 0;

	if (first_h->istack == STACK_OVERALL_HEADER) {
		overall_headlen = (int) first_h->headlen;
	}
	else {
		assert(image_index == 0);
	}

	int single_image_size = (int) (first_h->headlen + size * sizeof(float));
	int offset = overall_headlen + single_image_size * image_index;

	SpiderHeader *cur_image_hed;
	if (offset == 0) {
		cur_image_hed = first_h;
	}
	else {
		cur_image_hed = static_cast < SpiderHeader * >(calloc(sizeof(SpiderHeader), 1));
		portable_fseek(spider_file, offset, SEEK_SET);

		if (fread(cur_image_hed, sizeof(SpiderHeader), 1, spider_file) != 1) {
			LOGERR("read spider header with image_index = %d failed", image_index);
			return 1;
		}

		become_host_endian((float *) &cur_image_hed, NUM_FLOATS_IN_HEADER);

		if (cur_image_hed->nx != first_h->nx || cur_image_hed->ny != first_h->ny
			|| cur_image_hed->nslice != first_h->nslice) {
			LOGERR("%dth image size %dx%dx%d != overall size %dx%dx%d",
								 image_index, cur_image_hed->nx, cur_image_hed->ny,
								 cur_image_hed->nslice, first_h->nx, first_h->ny, first_h->nslice);
			return 1;
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

	dict["SPIDER.nslice"] = cur_image_hed->nslice;
	dict["SPIDER.type"] = cur_image_hed->type;
	dict["SPIDER.headrec"] = cur_image_hed->headrec;
	dict["SPIDER.angvalid"] = cur_image_hed->angvalid;

	dict["SPIDER.headlen"] = cur_image_hed->headlen;
	dict["SPIDER.dx"] = cur_image_hed->dx;
	dict["SPIDER.dy"] = cur_image_hed->dy;
	dict["SPIDER.dz"] = cur_image_hed->dz;

	dict["SPIDER.reclen"] = cur_image_hed->reclen;
	dict["SPIDER.istack"] = cur_image_hed->istack;
	dict["SPIDER.date"] = cur_image_hed->date;
	dict["SPIDER.name"] = cur_image_hed->name;
	dict["SPIDER.maxim"] = cur_image_hed->maxim;
	dict["SPIDER.imgnum"] = cur_image_hed->imgnum;

	if (offset != 0) {
		free(cur_image_hed);
		cur_image_hed = 0;
	}

	return 0;
}


int SpiderIO::write_header(const Dict & dict, int image_index, const Region* area, bool)
{
	ENTERFUNC;
	
	if (check_write_access(rw_mode, image_index) != 0) {
		return 1;
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
			LOGERR("cannot write to single-image Spider files");
			return 1;
		}

		if (nx1 != first_h->nx || ny1 != first_h->ny || nz1 != first_h->nslice) {
			LOGERR("new size %dx%dx%d != existing %s size %dx%dx%d",
								 filename.c_str(), nx1, ny1, nz1,
								 first_h->nx, first_h->ny, first_h->nslice);
			return 1;
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
			LOGERR("cannot write to Spider file '%s'", filename.c_str());
			return 1;
		}

		if (need_swap()) {
			ByteOrder::swap_bytes((float *) first_h, NUM_FLOATS_IN_HEADER);
		}

		first_h->maxim = (float) old_maxim;
		portable_fseek(spider_file, 0, SEEK_END);
	}
	else {
		if (!first_h) {
			first_h = (SpiderHeader *) calloc(header_length, 1);
			first_h->mmvalid = 0;
			first_h->sigma = -1.0;
			first_h->angvalid = 0;
			first_h->scale = 1.0;
			first_h->nslice = (float) nz1;
			first_h->nx = (float) nx1;
			first_h->ny = (float) ny1;

			if (nz1 == 1) {
				first_h->type = IMAGE_2D;
			}
			else {
				first_h->type = IMAGE_3D;
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
			LOGERR("cannot write a new header to Spider file '%s'", filename.c_str());
			return 1;
		}
	}

	if (!cur_h) {
		cur_h = (SpiderHeader *) calloc(header_size, 1);
	}

	cur_h->scale = 1.0;
	cur_h->nslice = (float)nz1;
	cur_h->nx = (float)nx1;
	cur_h->ny = (float)ny1;

	if (nz1 == 1) {
		cur_h->type = IMAGE_2D;
	}
	else {
		cur_h->type = IMAGE_3D;
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

	return 0;

}

int SpiderIO::write_single_header(const Dict & dict)
{
	ENTERFUNC;

	if (check_write_access(rw_mode, 0) != 0) {
		return 1;
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
		first_h = static_cast < SpiderHeader * >(calloc(header_size, 1));
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


	if (nz == 1) {
		first_h->type = IMAGE_2D;
	}
	else {
		first_h->type = IMAGE_3D;
	}

	first_h->reclen = (float)record_size;
	first_h->headrec = (float)num_records;
	first_h->headlen = (float)header_length;
	first_h->inuse = 1;

	if (fwrite(first_h, header_size, 1, spider_file) != 1) {
		LOGERR("write single spider header failed");
		return 1;
	}

	size_t pad_size = header_length - header_size;
	char *pad = static_cast < char *>(calloc(pad_size, 1));
	fwrite(pad, pad_size, 1, spider_file);
	free(pad);

	return 0;
}


int SpiderIO::read_data(float *data, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	if (check_read_access(image_index, true, data) != 0) {
		return 1;
	}

	if (check_region(area, IntSize((int) first_h->nx, (int) first_h->ny,
								   (int) first_h->nslice)) != 0) {
		return 1;
	}

	int overall_headlen = 0;
	if (first_h->istack == STACK_OVERALL_HEADER) {
		overall_headlen = (int) first_h->headlen;
	}
	else {
		assert(image_index == 0);
	}

	size_t size = static_cast < size_t > (first_h->nx * first_h->ny * first_h->nslice);
	size_t single_image_size = static_cast < size_t > (first_h->headlen + size * sizeof(float));
	off_t offset = overall_headlen + single_image_size * image_index;
	portable_fseek(spider_file, offset, SEEK_SET);

	size_t pad_size = static_cast < size_t > (first_h->headlen - sizeof(SpiderHeader));
	portable_fseek(spider_file, pad_size, SEEK_CUR);

	int err = EMUtil::get_region_data((unsigned char *) data, spider_file, 0, sizeof(float),
									  (int) first_h->nx, (int) first_h->ny, (int) first_h->nslice,
									  area);
	if (err) {
		return 1;
	}

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

	return 0;
}

int SpiderIO::write_data(float *data, int image_index, const Region* area, bool)
{
	ENTERFUNC;

	if (check_write_access(rw_mode, image_index, true, data) != 0) {
		return 1;
	}

	if (!cur_h) {
		LOGERR("please write header before write data");
		return 1;
	}

	if (cur_h->istack == SINGLE_IMAGE_HEADER) {
		LOGERR("cannot mix your single spider and stack spider.");
		return 1;
	}

	int hl = static_cast < int >(cur_h->headlen);
	int head_size = sizeof(SpiderHeader);
	int pad_size = hl - head_size;

	assert(pad_size >= 0);

	char *pad = static_cast < char *>(calloc(pad_size, 1));
	int img_size = static_cast < int >(cur_h->nx * cur_h->ny * cur_h->nslice);

	int sec_size = static_cast < int >(cur_h->nx * cur_h->ny * sizeof(float));
	int nz = static_cast < int >(cur_h->nslice);

	swap_header(cur_h);
	swap_data(data, img_size);

	while (image_index > first_h->maxim) {
		first_h->maxim++;
		first_h->inuse = 1.0;
		first_h->imgnum = first_h->maxim;

		portable_fseek(spider_file, 0, SEEK_END);
		fwrite(cur_h, head_size, 1, spider_file);
		fwrite(pad, pad_size, 1, spider_file);
		fwrite(data, sec_size, nz, spider_file);
	}

	off_t cur_offset = (off_t) (hl + (img_size * 1.0 * sizeof(float) + hl) * image_index);
	portable_fseek(spider_file, cur_offset, SEEK_SET);

	fwrite(cur_h, head_size, 1, spider_file);
	fwrite(pad, pad_size, 1, spider_file);

	int err = 0;
	if (fwrite(data, sec_size, nz, spider_file) != (unsigned int) nz) {
		LOGERR("cannot write to spider file '%s'", filename.c_str());
		err = 1;
	}

	swap_header(cur_h);
	swap_data(data, img_size);

	free(pad);
	return err;
}

int SpiderIO::write_single_data(float *data)
{
	ENTERFUNC;
	
	if (check_write_access(rw_mode, 0, true, data) != 0) {
		return 1;
	}

	assert(first_h != 0);

	if (first_h->istack != SINGLE_IMAGE_HEADER) {
		LOGERR("cannot mix your single spider and stack spider.");
		return 1;
	}

	portable_fseek(spider_file, (int) first_h->headlen, SEEK_SET);

	int sec_size = static_cast < int >(first_h->nx * first_h->ny * sizeof(float));
	unsigned int nz = static_cast < unsigned int >(first_h->nslice);

	if (fwrite(data, sec_size, nz, spider_file) != nz) {
		LOGERR("write data to single spider image failed");
		return 1;
	}
	return 0;
}

void SpiderIO::swap_data(float *data, int size)
{
	if (need_swap()) {
		ByteOrder::swap_bytes(data, size);
	}
}

void SpiderIO::swap_header(SpiderHeader * header)
{
	if (need_swap()) {
		ByteOrder::swap_bytes((float *) header, NUM_FLOATS_IN_HEADER);
	}
}

bool SpiderIO::is_complex_mode()
{
	return false;
}

bool SpiderIO::is_image_big_endian()
{
	return is_big_endian;
}

int SpiderIO::get_nimg()
{
	if (init() != 0) {
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
