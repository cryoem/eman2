/**
 * $Id$
 */
#include "gatan2io.h"
#include "log.h"
#include "emutil.h"
#include "geometry.h"
#include "portable_fileio.h"


using namespace EMAN;

Gatan2IO::Gatan2IO(const string & file, IOMode rw)
:	filename(file), rw_mode(rw), gatan2_file(0), initialized(false)
{
	is_big_endian = ByteOrder::is_host_big_endian();
	memset(&gatanh, 0, sizeof(Gatan2Header));
}

Gatan2IO::~Gatan2IO()
{
	if (gatan2_file) {
		fclose(gatan2_file);
		gatan2_file = 0;
	}
}

void Gatan2IO::init()
{
	ENTERFUNC;
	
	if (initialized) {
		return;
	}
	
	initialized = true;

	bool is_new_file = false;
	gatan2_file = sfopen(filename, rw_mode, &is_new_file);

	if (!is_new_file) {
		if (fread(&gatanh, sizeof(Gatan2Header), 1, gatan2_file) != 1) {
			throw ImageReadException(filename, "Gatan2 Header");
		}

		if (!is_valid(&gatanh)) {
			throw ImageReadException(filename, "invalid Gatan2 file");
		}

		is_big_endian = ByteOrder::is_data_big_endian(&gatanh.len);
		become_host_endian((short *) &gatanh, sizeof(Gatan2Header) / sizeof(short));
	}
	EXITFUNC;
}

bool Gatan2IO::is_valid(const void *first_block)
{
	ENTERFUNC;
	bool result = false;
	
	if (!first_block) {
		result = false;
	}
	else {
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
			result = true;
		}
	}
	EXITFUNC;
	return result;
}

int Gatan2IO::read_header(Dict & dict, int image_index, const Region * area, bool)
{
	ENTERFUNC;
	check_read_access(image_index);
	
	if (is_complex_mode()) {
		throw ImageReadException(filename, "Cannot read complex Gatan2 files");
	}
	else {
		check_region(area, IntSize(gatanh.nx, gatanh.ny));

		int xlen = 0, ylen = 0;
		EMUtil::get_region_dims(area, gatanh.nx, &xlen, gatanh.ny, &ylen);
				
		dict["nx"] = xlen;
		dict["ny"] = ylen;
		dict["nz"] = 1;
		dict["datatype"] = to_em_datatype(gatanh.type);	
	}
	
	EXITFUNC;
	return 0;
}

int Gatan2IO::write_header(const Dict &, int, const Region* , bool)
{
	ENTERFUNC;
	LOGWARN("Gatan2 write is not supported.");
	EXITFUNC;
	return 1;
}

int Gatan2IO::read_data(float *data, int image_index, const Region * area, bool )
{
	ENTERFUNC;
	check_read_access(image_index, data);

	if (is_complex_mode()) {
		throw ImageReadException(filename, "Cannot read complex Gatan2 files");
	}
	
	check_region(area, IntSize(gatanh.nx, gatanh.ny));
	
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

	EMUtil::process_region_io(cdata, gatan2_file, READ_ONLY, image_index, gatanh.len,
							  gatanh.nx, gatanh.ny, 1, area);

	switch (gatanh.type) {
	case GATAN2_SHORT:
		become_host_endian((short *) data, size);
		for (int i = size - 1; i >= 0; i--) {
			data[i] = static_cast < float >(sdata[i]);
		}
		break;
	case GATAN2_FLOAT:
		become_host_endian(data, size);
		break;
	case GATAN2_CHAR:
		for (int i = size - 1; i >= 0; i--) {
			data[i] = static_cast < float >(cdata[i]);
		}
		break;
	case GATAN2_INT:
		become_host_endian((int *) data, size);
		for (int i = size - 1; i >= 0; i--) {
			data[i] = static_cast < float >(ldata[i]);
		}
		break;
	default:
		throw ImageReadException(filename, "unsupported Gatan2 data type");
	}
	EXITFUNC;
	return 0;
}

int Gatan2IO::write_data(float *, int, const Region*, bool)
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
