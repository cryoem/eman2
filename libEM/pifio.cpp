/**
 * $Id$
 */
#include "pifio.h"
#include "log.h"
#include "emutil.h"
#include "portable_fileio.h"
#include "geometry.h"
#include <assert.h>

using namespace EMAN;

PifIO::PifIO(string pif_filename, IOMode rw)
    : filename(pif_filename), rw_mode(rw)
{
    pif_file = 0;
    mode_size = 0;
    is_big_endian = ByteOrder::is_machine_big_endian();
    initialized = false;
    real_scale_factor = 1;
    is_new_file = false;
    memset(&pfh, 0, sizeof(PifFileHeader));
}

PifIO::~PifIO()
{
    if (pif_file) {
	fclose(pif_file);
	pif_file = 0;
    }
}


int PifIO::get_mode_size(PifDataMode mode)
{
    int size = 0;

    switch (mode) {
    case PIF_CHAR:
    case PIF_BOXED_DATA:
	size = sizeof(char);
	break;
    case PIF_SHORT:
    case PIF_SHORT_FLOAT:
    case PIF_SHORT_COMPLEX:
    case PIF_SHORT_FLOAT_COMPLEX:
    case PIF_MAP_FLOAT_SHORT:
	size = sizeof(short);
	break;
    case PIF_FLOAT:
    case PIF_FLOAT_INT:
    case PIF_FLOAT_COMPLEX:
    case PIF_FLOAT_INT_COMPLEX:
    case PIF_MAP_FLOAT_INT:
	size = sizeof(int);
	break;
    default:
	break;
    }
    return size;
}

bool PifIO::is_float(int m)
{
    PifDataMode mode = static_cast<PifDataMode>(m);
    switch (mode) {
    case PIF_SHORT_FLOAT:
    case PIF_SHORT_FLOAT_COMPLEX:
    case PIF_FLOAT:
    case PIF_FLOAT_INT:
    case PIF_FLOAT_COMPLEX:
    case PIF_FLOAT_INT_COMPLEX:
    case PIF_MAP_FLOAT_SHORT:
    case PIF_MAP_FLOAT_INT:
	return true;
    default:
	break;
    }
    return false;
}

int PifIO::init()
{
    static int err = 0;
    if (initialized) {
	return err;
    }
    Log::logger()->log("PifIO::init()");
    initialized = true;

    pif_file = sfopen(filename, rw_mode, &is_new_file);

    if (!pif_file) {
	err = 1;
	return err;
    }

    if (!is_new_file) {
	if (fread(&pfh, sizeof(PifFileHeader), 1, pif_file) != 1) {
	    Log::logger()->error("read header failed on PIF file: '%s'", filename.c_str());
	    err = 1;
	    return err;
	}

	if (!is_valid(&pfh)) {
	    Log::logger()->error("'%s' is not a valid PIF file", filename.c_str());
	    err = 1;
	    return err;
	}

	is_big_endian = ByteOrder::is_data_big_endian(&pfh.nz);
	become_platform_endian(&pfh.htype);

	if (pfh.htype != 1) {
	    Log::logger()->error("only handle PIF that all projects have the same dimensions");
	    return 1;
	}

	become_platform_endian(&pfh.mode);
	become_platform_endian(&pfh.nx);
	become_platform_endian(&pfh.ny);
	become_platform_endian(&pfh.nz);
	become_platform_endian(&pfh.nimg);

	if (is_float(pfh.mode)) {
	    real_scale_factor = atof(pfh.scalefactor);
	}

	mode_size = get_mode_size(static_cast<PifDataMode>(pfh.mode));

	if (is_complex_mode()) {
	    pfh.nx *= 2;
	}
    }

    return 0;
}

bool PifIO::is_valid(const void *first_block)
{
    Log::logger()->log("PifIO::is_valid()");
    if (!first_block) {
	return false;
    }

    const int *data = static_cast<const int *>(first_block);
    int m1 = data[0];
    int m2 = data[1];
    int endian = data[7];
    bool data_big_endian = static_cast<bool>(endian);

    if (data_big_endian != ByteOrder::is_machine_big_endian()) {
	ByteOrder::swap_bytes(&m1);
	ByteOrder::swap_bytes(&m2);
    }

    if (m1 == PIF_MAGIC_NUM && m2 == PIF_MAGIC_NUM) {
	return true;
    }

    return false;
}

void PifIO::fseek_to(int image_index)
{
    int pih_sz = sizeof(PifImageHeader);
    int image_size = 0;

    if (pfh.nimg == 1) {
	image_size = pfh.nx * pfh.ny * pfh.nz;
    }
    else {
	assert(pfh.nz == pfh.nimg);
	image_size = pfh.nx * pfh.ny;
    }

    size_t file_offset = sizeof(PifFileHeader) + (pih_sz + image_size * mode_size) * image_index;

    portable_fseek(pif_file, file_offset, SEEK_SET);
}


int PifIO::read_header(Dict & dict, int image_index, const Region * area, bool)
{
    Log::logger()->log("PifIO::read_header() from file '%s'", filename.c_str());
    if (check_read_access(image_index) != 0) {
	return 1;
    }

    int pih_sz = sizeof(PifImageHeader);
    fseek_to(image_index);

    PifImageHeader pih;

    if (fread(&pih, pih_sz, 1, pif_file) != 1) {
	Log::logger()->error("read pif file '%s' failed", filename.c_str());
	return 1;
    }
    if (check_region(area, Size(pih.nx, pih.ny, pih.nz)) != 0) {
	return 1;
    }
    int xlen = 0, ylen = 0, zlen = 0;
    EMUtil::get_region_dims(area, pih.nx, &xlen, pih.ny, &ylen, pih.nz, &zlen);

    dict["nx"] = EMObject(xlen);
    dict["ny"] = EMObject(ylen);
    dict["nz"] = EMObject(zlen);

    dict["datatype"] = EMObject(to_em_datatype(pih.mode));

    dict["spacing_row"] = EMObject(static_cast<float>(pih.xlen));
    dict["spacing_col"] = EMObject(static_cast<float>(pih.ylen));
    dict["spacing_sec"] = EMObject(static_cast<float>(pih.zlen));

    dict["minimum"] = EMObject(static_cast<float>(pih.min));
    dict["maximum"] = EMObject(static_cast<float>(pih.max));
    dict["mean"] = EMObject(static_cast<float>(pih.mean));
    dict["sigma"] = EMObject(static_cast<float>(pih.sigma));

    dict["origin_row"] = EMObject(static_cast<float>(pih.xorigin));
    dict["origin_col"] = EMObject(static_cast<float>(pih.yorigin));

    return 0;
}

int PifIO::write_header(const Dict & dict, int image_index)
{
    Log::logger()->log("PifIO::write_header() to file '%s'", filename.c_str());

    if (check_write_access(rw_mode, image_index) != 0) {
	return 1;
    }

    time_t t0 = time(0);
    struct tm *t = localtime(&t0);

    if (!is_new_file) {
	if (is_big_endian != ByteOrder::is_machine_big_endian()) {
	    Log::logger()->
		error("cannot write to existing file '%s' which is in opposite byte order",
		      filename.c_str());
	    return 1;
	}

	int nx1 = dict["nx"].get_int();
	int ny1 = dict["ny"].get_int();
	int nz1 = dict["nz"].get_int();

	int mode1 = to_pif_datatype(dict["datatype"].get_int());

	if (nx1 != pfh.nx || ny1 != pfh.ny || nz1 != pfh.nz) {
	    Log::logger()->error("cannot write to different size file %s (%dx%dx%d != %dx%dx%d)",
				 filename.c_str(), nx1, ny1, nz1, pfh.nx, pfh.ny, pfh.nz);
	    return 1;
	}

	if (mode1 != pfh.mode) {
	    Log::logger()->error("cannot write to different data type file %s",
				 filename.c_str());
	    return 1;
	}

	fseek_to(image_index);
    }
    else {
	pfh.magic[0] = PIF_MAGIC_NUM;
	pfh.magic[1] = PIF_MAGIC_NUM;
	sprintf(pfh.scalefactor, "1.0");

	pfh.mode = PIF_FLOAT_INT;
	sprintf(pfh.program, "EMAN %d/%02d/%02d", t->tm_mon, t->tm_mday, t->tm_year);

	pfh.htype = 1;
	pfh.nx = dict["nx"].get_int();
	pfh.ny = dict["ny"].get_int();
	pfh.nz = dict["nz"].get_int();
	pfh.nimg += 1;
	pfh.endian = (int) ByteOrder::is_machine_big_endian();
	fwrite(&pfh, sizeof(PifFileHeader), 1, pif_file);
    }

    PifImageHeader pih;
    memset(&pih, 0, sizeof(PifImageHeader));
    pih.nx = dict["nx"].get_int();
    pih.ny = dict["ny"].get_int();
    pih.nz = dict["nz"].get_int();

    pih.mode = PIF_FLOAT;
    pih.xlen = dict["spacing_row"].get_int();
    pih.ylen = dict["spacing_col"].get_int();
    pih.zlen = dict["spacing_sec"].get_int();
    pih.alpha = 90;
    pih.beta = 90;
    pih.gamma = 90;
    pih.mapc = 1;
    pih.mapr = 2;
    pih.maps = 3;
    pih.min = dict["minimum"].get_int();
    pih.max = dict["maximum"].get_int();
    pih.mean = dict["mean"].get_int();
    pih.sigma = dict["sigma"].get_int();

    pih.xorigin = dict["origin_row"].get_int();
    pih.yorigin = dict["origin_col"].get_int();

    sprintf(pih.time, "%d/%02d/%02d %d:%02d",
	    t->tm_mon, t->tm_mday, t->tm_year, t->tm_hour, t->tm_min);
    fwrite(&pih, sizeof(PifImageHeader), 1, pif_file);

    return 0;
}

int PifIO::read_data(float *data, int image_index, const Region *, bool)
{
    Log::logger()->log("PifIO::read_data() from file '%s'", filename.c_str());

    if (check_read_access(image_index, true, data) != 0) {
	return 1;
    }

    int pih_sz = sizeof(PifImageHeader);
    fseek_to(image_index);
    portable_fseek(pif_file, pih_sz, SEEK_CUR);

    int buf_size = pfh.nx * mode_size;
    unsigned char *buf = new unsigned char[buf_size];

    int num_layers = pfh.nz;
    if (pfh.nz == pfh.nimg) {
	num_layers = 1;
    }

    for (int l = 0; l < num_layers; l++) {
	int offset1 = l * pfh.nx * pfh.ny;
	for (int j = 0; j < pfh.ny; j++) {
	    if (fread(buf, mode_size, pfh.nx, pif_file) != (unsigned int) pfh.nx) {
		Log::logger()->error("read PIF image file '%s' failed", filename.c_str());
		delete [] buf;
		buf = 0;
		return 1;
	    }

	    if (mode_size == sizeof(short)) {
		become_platform_endian((short *) buf, pfh.nx);
	    }
	    else if (mode_size == sizeof(int)) {
		become_platform_endian((int *) buf, pfh.nx);
	    }

	    int offset2 = offset1 + j * pfh.nx;

	    for (int k = 0; k < pfh.nx; k++) {
		float curr_data = 0;

		if (mode_size == sizeof(char)) {
		    curr_data = (float) buf[k];
		}
		else if (mode_size == sizeof(short)) {
		    curr_data = (float) ((short *) buf)[k];
		}
		else if (mode_size == sizeof(int)) {
		    curr_data = (float) ((int *) buf)[k];
		}
		data[offset2 + k] = curr_data * real_scale_factor;
	    }
	}
    }

    delete [] buf;
    buf = 0;
    return 0;
}


int PifIO::write_data(float *data, int image_index)
{

    Log::logger()->log("PifIO::write_data() to file '%s'", filename.c_str());
    if (check_write_access(rw_mode, image_index, true, data) != 0) {
	return 1;
    }

    fseek_to(image_index);

    int nx = pfh.nx;
    int ny = pfh.ny;
    int nz = pfh.nz;

    int *buf = new int[nx];
    for (int i = 0; i < nz; i++) {
	for (int j = 0; j < ny; j++) {
	    for (int k = 0; k < pfh.nx; k++) {
		buf[k] = (int) data[i * nx * ny + j * nx + k];
	    }
	    fwrite(buf, sizeof(int) * nx, 1, pif_file);
	}
    }

    delete [] buf;
    buf = 0;
    return 0;
}


bool PifIO::is_complex_mode()
{
    if (pfh.mode == PIF_SHORT_COMPLEX ||
	pfh.mode == PIF_FLOAT_INT_COMPLEX ||
	pfh.mode == PIF_FLOAT_COMPLEX || pfh.mode == PIF_SHORT_FLOAT_COMPLEX) {
	return true;
    }
    return false;
}

bool PifIO::is_image_big_endian()
{
    return is_big_endian;
}

int PifIO::get_nimg()
{
    if (init() != 0) {
	return 0;
    }

    return pfh.nimg;
}


int PifIO::to_em_datatype(int p)
{
    PifDataMode mode = static_cast<PifDataMode>(p);
    EMUtil::EMDataType e = EMUtil::EM_UNKNOWN;

    switch (mode) {
    case PIF_CHAR:
    case PIF_BOXED_DATA:
	e = EMUtil::EM_CHAR;
	break;

    case PIF_SHORT:
    case PIF_SHORT_FLOAT:
    case PIF_MAP_FLOAT_SHORT:
	e = EMUtil::EM_SHORT;
	break;

    case PIF_SHORT_COMPLEX:
    case PIF_SHORT_FLOAT_COMPLEX:
	e = EMUtil::EM_SHORT_COMPLEX;
	break;

    case PIF_FLOAT:
    case PIF_FLOAT_INT:
    case PIF_MAP_FLOAT_INT:
	e = EMUtil::EM_FLOAT;
	break;
    case PIF_FLOAT_COMPLEX:
    case PIF_FLOAT_INT_COMPLEX:
	e = EMUtil::EM_FLOAT_COMPLEX;
	break;
    case PIF_INVALID:
	e = EMUtil::EM_UNKNOWN;
	break;
    }
    return e;
}

int PifIO::to_pif_datatype(int e)
{
    PifDataMode m = PIF_INVALID;

    switch (e) {
    case EMUtil::EM_CHAR:
	m = PIF_BOXED_DATA;
	break;
    case EMUtil::EM_SHORT:
	m = PIF_SHORT;
	break;
    case EMUtil::EM_SHORT_COMPLEX:
	m = PIF_SHORT_COMPLEX;
	break;
    case EMUtil::EM_FLOAT:
	m = PIF_FLOAT_INT;
	break;
    case EMUtil::EM_FLOAT_COMPLEX:
	m = PIF_FLOAT_COMPLEX;
	break;
    default:
	Log::logger()->error("unknown PIF mode: %d", e);
    }

    return m;
}
