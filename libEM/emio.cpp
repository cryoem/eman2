/**
 * $Id$
 */
#include "emio.h"
#include "log.h"
#include "portable_fileio.h"
#include "emutil.h"
#include "geometry.h"
#include <assert.h>

using namespace EMAN;

EmIO::EmIO(string file, IOMode rw)
:	filename(file), rw_mode(rw), em_file(0), initialized(false)
{
	mode_size = 0;
	mode = EM_EM_UNKNOWN;
	is_big_endian = ByteOrder::is_host_big_endian();
}

EmIO::~EmIO()
{
	if (em_file) {
		fclose(em_file);
		em_file = 0;
	}
}

int EmIO::init()
{
	ENTERFUNC;
	
	static int err = 0;
	if (initialized) {
		return err;
	}

	
	initialized = true;

	bool is_new_file = false;
	em_file = sfopen(filename, rw_mode, &is_new_file);

	if (!em_file) {
		err = 1;
		return err;
	}

	if (!is_new_file) {
		if (fread(&emh, sizeof(EMHeader), 1, em_file) != 1) {
			LOGERR("cannot read EM image file '%s'", filename.c_str());
			err = 1;
			return err;
		}
		if (!is_valid(&emh)) {
			LOGERR("'%s' is not a valid EM image", filename.c_str());
			err = 1;
			return err;
		}

		is_big_endian = ByteOrder::is_data_big_endian(&emh.nz);
		become_host_endian(&emh.nx);
		become_host_endian(&emh.ny);
		become_host_endian(&emh.nz);

		mode = (DataType) emh.data_type;

		if (mode == EM_EM_DOUBLE) {
			LOGERR("DOUBLE data type is not supported for EM image");
			err = 1;
			return err;
		}

		mode_size = get_mode_size(emh.data_type);
		if (is_complex_mode()) {
			emh.nx *= 2;
		}
	}

	return 0;
}

bool EmIO::is_valid(const void *first_block, off_t file_size)
{
	ENTERFUNC;

	if (!first_block) {
		return false;
	}

	const char *data = static_cast < const char *>(first_block);
	char machine = data[0];
	char is_new_ver = data[1];
	char data_type = data[3];

	const int *data1 = static_cast < const int *>(first_block);
	int nx = data1[1];
	int ny = data1[2];
	int nz = data1[3];

	bool data_big_endian = ByteOrder::is_data_big_endian(&nz);
	if (data_big_endian != ByteOrder::is_host_big_endian()) {
		ByteOrder::swap_bytes(&nx);
		ByteOrder::swap_bytes(&ny);
		ByteOrder::swap_bytes(&nz);
	}

	const int max_dim = 1 << 20;

	if (((int) machine >= EM_OS8 && machine <= EM_PC) &&
		(is_new_ver == 0 || is_new_ver == 1) &&
		(data_type >= EM_EM_CHAR && data_type <= EM_EM_DOUBLE) &&
		(nx > 1 && nx < max_dim) && (ny > 0 && ny < max_dim) && (nz > 0 && nz < max_dim)) {
		if (file_size > 0) {
			off_t file_size1 = nx * ny * nz * get_mode_size(data_type) + sizeof(EMHeader);
			if (file_size == file_size1) {
				return true;
			}
		}
		else {
			return true;
		}
	}

	return false;
}

int EmIO::read_header(Dict & dict, int image_index, const Region * area, bool is_3d)
{
	ENTERFUNC;

	if (check_read_access(image_index, area) != 0) {
		return 1;
	}
	if (check_region(area, IntSize(emh.nx, emh.ny, emh.nz)) != 0) {
		return 1;
	}

	int nimg = 1;
	if (is_3d) {
		nimg = emh.nz;
	}

	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, emh.nx, &xlen, emh.ny, &ylen, nimg, &zlen);

	dict["nx"] = xlen;
	dict["ny"] = ylen;
	dict["nz"] = 1;
	dict["datatype"] = to_em_datatype(emh.data_type);
	EXITFUNC;
	return 0;
}

int EmIO::write_header(const Dict & dict, int image_index, const Region* area, bool)
{
	ENTERFUNC;
	int err = 0;
	
	if (check_write_access(rw_mode, image_index) != 0) {
		err = 1;
	}
	else {
		emh.machine = static_cast < char >(get_machine_type());
		emh.nx = dict["nx"];
		emh.ny = dict["ny"];
		emh.nz = dict["nz"];
		emh.data_type = EM_EM_FLOAT;

		if (fwrite(&emh, sizeof(EMHeader), 1, em_file) != 1) {
			LOGERR("cannot write header to file '%s'", filename.c_str());
			err = 1;
		}
	}
	EXITFUNC;
	return err;
}

int EmIO::read_data(float *data, int image_index, const Region * area, bool is_3d)
{
	ENTERFUNC;

	if (check_read_access(image_index, true, data) != 0) {
		return 1;
	}

	if (check_region(area, IntSize(emh.nx, emh.ny, emh.nz)) != 0) {
		return 1;
	}

	int nimg = 1;
	if (is_3d) {
		image_index = 0;
		nimg = emh.nz;
	}

	size_t header_sz = sizeof(EMHeader);
	off_t file_offset = header_sz + emh.nx * emh.ny * mode_size * image_index;
	portable_fseek(em_file, file_offset, SEEK_SET);

	unsigned char *cdata = (unsigned char *) data;
	int err = EMUtil::get_region_data(cdata, em_file, image_index, mode_size,
									  emh.nx, emh.ny, nimg, area);
	if (err) {
		return 1;
	}

	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, emh.nx, &xlen, emh.ny, &ylen, nimg, &zlen);

	int total_sz = xlen * ylen * zlen;

	if (mode_size == sizeof(short)) {
		become_host_endian((short *) cdata, total_sz);
	}
	else if (mode_size == sizeof(int)) {
		become_host_endian((int *) cdata, total_sz);
	}
	else if (mode_size == sizeof(double)) {
		LOGERR("double type image is not supported");
		return 1;
	}

	for (int k = total_sz - 1; k >= 0; k--) {
		float curr_data = 0;

		if (mode == EM_EM_CHAR) {
			curr_data = static_cast < float >(cdata[k]);
		}
		else if (mode == EM_EM_SHORT) {
			curr_data = static_cast < float >(((short *) cdata)[k]);
		}
		else if (mode == EM_EM_INT) {
			curr_data = static_cast < float >(((int *) cdata)[k]);
		}
		else if (mode == EM_EM_FLOAT || mode == EM_EM_COMPLEX) {
			curr_data = ((float *) cdata)[k];
		}
		else if (mode_size == sizeof(double)) {
			LOGERR("double type image is not supported");
			return 1;
		}

		data[k] = curr_data;
	}

	EXITFUNC;
	return 0;
}

int EmIO::write_data(float *data, int image_index, const Region* area, bool)
{
	ENTERFUNC;

	if (check_write_access(rw_mode, image_index, true, data) != 0) {
		return 1;
	}

	int sec_size = emh.nx * emh.ny;
	int row_size = sizeof(float) * emh.nx;

	for (int i = 0; i < emh.nz; i++) {
		int k = i * sec_size;
		for (int j = 0; j < emh.ny; j++) {
			fwrite(&data[k + j * emh.nx], row_size, 1, em_file);
		}
	}
	EXITFUNC;
	return 0;
}

void EmIO::flush()
{
	fflush(em_file);
}

bool EmIO::is_complex_mode()
{
	if (emh.data_type == EM_EM_COMPLEX) {
		return true;
	}
	return false;
}

bool EmIO::is_image_big_endian()
{
	return is_big_endian;
}

int EmIO::get_nimg()
{
	if (init() != 0) {
		return 0;
	}

	return emh.nz;
}

int EmIO::to_em_datatype(char t)
{
	DataType type = static_cast < DataType > (t);
	switch (type) {
	case EM_EM_CHAR:
		return EMUtil::EM_CHAR;
	case EM_EM_SHORT:
		return EMUtil::EM_SHORT;
	case EM_EM_INT:
		return EMUtil::EM_INT;
	case EM_EM_FLOAT:
		return EMUtil::EM_FLOAT;
	case EM_EM_DOUBLE:
		return EMUtil::EM_DOUBLE;
	case EM_EM_COMPLEX:
		return EMUtil::EM_FLOAT_COMPLEX;
	default:
		break;
	}
	return EMUtil::EM_UNKNOWN;
}

int EmIO::get_machine_type()
{
	int m = EM_UNKNOWN_MACHINE;
#ifdef __sgi
	m = EM_SGI;
#elif defined __linux__
	m = EM_PC;
#elif defined __CYGWIN__
	m = EM_PC;
#elif defined WIN32
	m = EM_PC;
#else
	m = EM_UNKNOWN_MACHINE;
#endif
	return m;
}

size_t EmIO::get_mode_size(char data_type)
{
	int mode = (int) data_type;
	switch (mode) {
	case EM_EM_CHAR:
		return sizeof(char);
	case EM_EM_SHORT:
		return sizeof(short);
	case EM_EM_INT:
	case EM_EM_FLOAT:
	case EM_EM_COMPLEX:
		return sizeof(int);
	case EM_EM_DOUBLE:
		return sizeof(double);
	}
	return 0;
}
