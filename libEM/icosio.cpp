/**
 * $Id$
 */
#include "icosio.h"
#include "log.h"
#include "emutil.h"
#include "portable_fileio.h"
#include "geometry.h"

using namespace EMAN;


IcosIO::IcosIO(string file, IOMode rw)
    :  filename(file), rw_mode(rw), icos_file(0), initialized(false)
{
    is_big_endian = ByteOrder::is_machine_big_endian();
    is_new_file = false;
}

IcosIO::~IcosIO()
{
    if (icos_file) {
	fclose(icos_file);
	icos_file = 0;
    }
}

int IcosIO::init()
{
    static int err = 0;
    if (initialized) {
	return err;
    }
    Log::logger()->log("IcosIO::init()");
    initialized = true;

    icos_file = sfopen(filename, rw_mode, &is_new_file);
    if (!icos_file) {
	err = 1;
	return err;
    }

    if (!is_new_file) {
	if (fread(&icosh, sizeof(IcosHeader), 1, icos_file) != 1) {
	    Log::logger()->error("cannot read ICOS file '%s'", filename.c_str());
	    err = 1;
	    return err;
	}

	if (!is_valid(&icosh)) {
	    Log::logger()->error("invalid ICOS file");
	    err = 1;
	    return err;
	}

	become_platform_endian((int *) &icosh, sizeof(IcosHeader) / sizeof(int));
	is_big_endian = ByteOrder::is_data_big_endian(&icosh.stamp);
    }

    return 0;
}

bool IcosIO::is_valid(const void *first_block)
{
    Log::logger()->log("IcosIO::is_valid()");
    if (!first_block) {
	return false;
    }

    const int *data = static_cast<const int *>(first_block);
    int stamp = data[0];
    int stamp1 = data[19];
    int stamp2 = data[20];
    int stamp3 = data[26];

    bool data_big_endian = ByteOrder::is_data_big_endian(&stamp);

    if (data_big_endian != ByteOrder::is_machine_big_endian()) {
	ByteOrder::swap_bytes(&stamp);
	ByteOrder::swap_bytes(&stamp1);
	ByteOrder::swap_bytes(&stamp2);
	ByteOrder::swap_bytes(&stamp3);
    }

    if (stamp == STAMP && stamp1 == STAMP1 && stamp2 == STAMP2 && stamp3 == STAMP3) {
	return true;
    }
    return false;
}

int IcosIO::read_header(Dict & dict, int image_index,const Region * area, bool is_3d)
{
    Log::logger()->log("IcosIO::read_header() from file '%s'", filename.c_str());

    if (check_read_access(image_index) != 0) {
	return 1;
    }

    int nimg = 1;
    if (is_3d) {
	nimg = icosh.nz;
    }

    if (check_region(area, Size(icosh.nx, icosh.ny, nimg)) != 0) {
	return 1;
    }

    int xlen = 0, ylen = 0, zlen = 0;
    EMUtil::get_region_dims(area, icosh.nx, &xlen, icosh.ny, &ylen, nimg, &zlen);

    dict["nx"] = xlen;
    dict["ny"] = ylen;
    dict["nz"] = zlen;

    dict["datatype"] = EMUtil::EM_FLOAT;
    dict["minimum"] = icosh.min;
    dict["maximum"] = icosh.max;

    return 0;
}

int IcosIO::write_header(const Dict & dict, int image_index)
{
    Log::logger()->log("IcosIO::write_header() to file '%s'", filename.c_str());
    if (check_write_access(rw_mode, image_index) != 0) {
	return 1;
    }

    if (image_index > 0) {
	Log::logger()->error("image_index = %d. appending to ICOS file is not supported.",
			     image_index);
	return 1;
    }

    icosh.stamp = STAMP;
    icosh.stamp1 = STAMP1;
    icosh.stamp2 = STAMP2;
    icosh.stamp3 = STAMP3;

    icosh.nx = dict["nx"].get_int();
    icosh.ny = dict["ny"].get_int();
    icosh.nz = dict["nz"].get_int();

    icosh.min = dict["minimum"].get_float();
    icosh.max = dict["maximum"].get_float();

    if (fwrite(&icosh, sizeof(IcosHeader), 1, icos_file) != 1) {
	Log::logger()->error("cannot write header to file '%s'", filename.c_str());
	return 1;
    }

    return 0;
}

int IcosIO::read_data(float *data, int image_index, const Region * area, bool is_3d)
{
    Log::logger()->log("IcosIO::read_data() from file '%s'", filename.c_str());

    if (check_read_access(image_index, true, data) != 0) {
	return 1;
    }

    if (is_3d && image_index != 0) {
	Log::logger()->warn("read whole 3D image. reset image index from %d to 0", image_index);
	image_index = 0;
    }

    int nimg = 1;
    if (is_3d) {
	nimg = icosh.nz;
    }

    if (check_region(area, Size(icosh.nx, icosh.ny, nimg)) != 0) {
	return 1;
    }

    portable_fseek(icos_file, sizeof(IcosHeader), SEEK_SET);

    int err = EMUtil::get_region_data((unsigned char *) data, icos_file, image_index,
				      sizeof(float), icosh.nx, icosh.ny, nimg, area,
				      false, sizeof(int), sizeof(int));
    if (err) {
	return 1;
    }

#if 0
    int row_size = icosh.nx * sizeof(float) + 2 * sizeof(int);
    off_t img_size = row_size * icosh.ny;
    off_t offset = sizeof(IcosHeader) + img_size * image_index;
    int int_size = sizeof(int);
    int nrows = icosh.ny * icosh.nz;
    if (!is_3d) {
	nrows = icosh.ny;
    }
    int data_row_size = icosh.nx * sizeof(float);

    for (int j = 0; j < nrows; j++) {
	portable_fseek(icos_file, int_size, SEEK_CUR);
	if (fread(&data[j * icosh.nx], data_row_size, 1, icos_file) != 1) {
	    Log::logger()->error("incomplete data read %d/%d blocks on file '%s'", j, nrows,
				 filename.c_str());
	}
	portable_fseek(icos_file, int_size, SEEK_CUR);
    }
#endif


    int xlen = 0, ylen = 0, zlen = 0;
    EMUtil::get_region_dims(area, icosh.nx, &xlen, icosh.ny, &ylen, nimg, &zlen);
    become_platform_endian(data, xlen * ylen * zlen);

    return 0;
}

int IcosIO::write_data(float *data, int image_index)
{
    Log::logger()->log("IcosIO::write_data() to file '%s'", filename.c_str());
    if (check_write_access(rw_mode, image_index, true, data) != 0) {
	return 1;
    }

    if (icosh.nx <= 0) {
	Log::logger()->error("Please write a valid ICOS header before you write the data");
	return 1;
    }

    if (image_index > 0) {
	Log::logger()->error("image_index = %d. appending to ICOS file is not supported.",
			     image_index);
	return 1;
    }

    int float_size = sizeof(float);
    int nx = icosh.nx;
    float *buf = new float[nx + 2];
    buf[0] = float_size * nx;
    buf[nx + 1] = buf[0];
    int nrows = icosh.ny * icosh.nz;

    int row_size = (nx + 2) * float_size;
    int err = 0;

    for (int j = 0; j < nrows; j++) {
	memcpy(&buf[1], &data[nx * j], nx * float_size);
	if (fwrite(buf, row_size, 1, icos_file) != 1) {
	    Log::logger()->error("writing ICOS data out failed");
	    err = 1;
	    break;
	}
    }

    delete [] buf;
    buf = 0;

    return err;
}

bool IcosIO::is_complex_mode()
{
    return false;
}

bool IcosIO::is_image_big_endian()
{
    return is_big_endian;
}

int IcosIO::get_nimg()
{
    if (init() != 0) {
	return 0;
    }

    return icosh.nz;
}
