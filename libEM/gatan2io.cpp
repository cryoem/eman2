/**
 * $Id$
 */
#include "gatan2io.h"
#include "log.h"
#include "emutil.h"
#include "geometry.h"
#include "portable_fileio.h"
#include <assert.h>

using namespace EMAN;

Gatan2IO::Gatan2IO(string file, IOMode rw)
:	filename(file), rw_mode(rw), gatan2_file(0), initialized(false)
{
	is_big_endian = ByteOrder::is_host_big_endian();
}

Gatan2IO::~Gatan2IO()
{
	if (gatan2_file) {
		fclose(gatan2_file);
		gatan2_file = 0;
	}
}

int Gatan2IO::init()
{
	ENTERFUNC;
	
	static int err = 0;
	if (initialized) {
		return err;
	}
	
	initialized = true;

	bool is_new_file = false;
	gatan2_file = sfopen(filename, rw_mode, &is_new_file);

	if (!gatan2_file) {
		err = 1;
		return err;
	}

	if (!is_new_file) {
		if (fread(&gatanh, sizeof(Gatan2Header), 1, gatan2_file) != 1) {
			LOGERR("cannot read header from Gatan2 file '%s'", filename.c_str());
			err = 1;
			return err;
		}

		if (!is_valid(&gatanh)) {
			LOGERR("'%s' is not a valid Gatan2 file", filename.c_str());
			err = 1;
			return err;
		}

		is_big_endian = ByteOrder::is_data_big_endian(&gatanh.len);
		become_host_endian((short *) &gatanh, sizeof(Gatan2Header) / sizeof(short));
	}

	return 0;
}

bool Gatan2IO::is_valid(const void *first_block)
{
	ENTERFUNC;

	if (!first_block) {
		return false;
	}

	const short *data = static_cast < const short *>(first_block);
	short len = data[5];
	short type = data[6];

	bool data_big_endian = ByteOrder::is_data_big_endian(&len);

	if (data_big_endian != ByteOrder::is_host_big_endian()) {
		ByteOrder::swap_bytes(&len);
		ByteOrder::swap_bytes(&type);
	}

	int double_size = sizeof(double);
	if (len > 0 && len <= double_size && type > 0 && type <= GATAN2_INVALID) {
		return true;
	}
	return false;
}

int Gatan2IO::read_header(Dict & dict, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	if (check_read_access(image_index) != 0) {
		return 1;
	}

	if (is_complex_mode()) {
		LOGERR("Cannot read complex Gatan2 files");
		return 1;
	}
	if (check_region(area, IntSize(gatanh.nx, gatanh.ny)) != 0) {
		return 1;
	}
	int xlen = 0, ylen = 0;
	EMUtil::get_region_dims(area, gatanh.nx, &xlen, gatanh.ny, &ylen);

	dict["nx"] = xlen;
	dict["ny"] = ylen;
	dict["nz"] = 1;
	dict["datatype"] = to_em_datatype(gatanh.type);
	EXITFUNC;
	return 0;
}

int Gatan2IO::write_header(const Dict &, int, const Region* area, bool)
{
	ENTERFUNC;
	LOGWARN("Gatan2 write is not supported.");
	EXITFUNC;
	return 1;
}

int Gatan2IO::read_data(float *data, int image_index, const Region * area, bool)
{
	ENTERFUNC;
	if (check_read_access(image_index, true, data) != 0) {
		return 1;
	}

	if (is_complex_mode()) {
		LOGERR("Cannot read complex Gatan2 files");
		return 1;
	}
	if (check_region(area, IntSize(gatanh.nx, gatanh.ny)) != 0) {
		return 1;
	}

	portable_fseek(gatan2_file, sizeof(Gatan2Header), SEEK_SET);

#if 0
	if (fread(data, gatanh.nx * gatanh.len, gatanh.ny, gatan2_file) != (unsigned int) gatanh.ny) {
		LOGDEBUG("Data read incomplete in Gatan file '%s'", filename.c_str());
		return 1;
	}
#endif

	int size = gatanh.nx * gatanh.ny;
	short *sdata = (short *) data;
	unsigned char *cdata = (unsigned char *) data;
	int *ldata = (int *) data;

	int err = EMUtil::get_region_data(cdata, gatan2_file, image_index, gatanh.len,
									  gatanh.nx, gatanh.ny, 1, area);
	if (err) {
		return 1;
	}

	int i = 0;

	switch (gatanh.type) {
	case GATAN2_SHORT:
		become_host_endian((short *) data, size);
		for (i = size - 1; i >= 0; i--) {
			data[i] = static_cast < float >(sdata[i]);
		}
		break;
	case GATAN2_FLOAT:
		become_host_endian(data, size);
		break;
	case GATAN2_CHAR:
		for (i = size - 1; i >= 0; i--) {
			data[i] = static_cast < float >(cdata[i]);
		}
		break;
	case GATAN2_INT:
		become_host_endian((int *) data, size);
		for (i = size - 1; i >= 0; i--) {
			data[i] = static_cast < float >(ldata[i]);
		}
		break;
	default:
		LOGERR("don't know how to handle this type");
		return 1;
	}
	EXITFUNC;
	return 0;
}

int Gatan2IO::write_data(float *, int, const Region* area, bool)
{
	ENTERFUNC;
	LOGWARN("Gatan2 write is not supported.");
	EXITFUNC;
	return 1;
}

void Gatan2IO::flush()
{
}

bool Gatan2IO::is_complex_mode()
{
	if (gatanh.type == GATAN2_COMPLEX || gatanh.type == GATAN2_PACKED_COMPLEX) {
		return true;
	}
	return false;
}

bool Gatan2IO::is_image_big_endian()
{
	return is_big_endian;
}

int Gatan2IO::get_nimg()
{
	if (init() != 0) {
		return 0;
	}

	return 1;
}

int Gatan2IO::to_em_datatype(int gatan_type)
{
	switch (gatan_type) {
	case GATAN2_SHORT:
		return EMUtil::EM_SHORT;

	case GATAN2_FLOAT:
		return EMUtil::EM_FLOAT;

	case GATAN2_CHAR:
		return EMUtil::EM_CHAR;

	case GATAN2_INT:
		return EMUtil::EM_INT;
	}

	return EMUtil::EM_UNKNOWN;
}
