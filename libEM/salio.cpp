/**
 * $Id$
 */
#include "salio.h"
#include "log.h"
#include "emutil.h"
#include "util.h"
#include "geometry.h"

#ifndef WIN32
#include <sys/param.h>
#endif
#include <limits.h>
#include <assert.h>

using namespace EMAN;

const char *SalIO::HDR_EXT = ".hdr";
const char *SalIO::IMG_EXT = ".img";
const char *SalIO::MAGIC = " IDENTIFICATION";


SalIO::SalIO(string file, IOMode rw)
    :  filename(file), rw_mode(rw), sal_file(0), initialized(false)
{
    nx = 0;
    ny = 0;
    record_length = 512;
    scan_mode = NON_RASTER_SCAN;
    pixel = 4.6667;
}

SalIO::~SalIO()
{
    if (sal_file) {
	fclose(sal_file);
	sal_file = 0;
    }
}

int SalIO::init()
{
    static int err = 0;
    if (initialized) {
	return err;
    }
    Log::logger()->log("SalIO::init()");
    initialized = true;

    string hdr_filename = Util::get_filename_by_ext(filename, HDR_EXT);
    string img_filename = Util::get_filename_by_ext(filename, IMG_EXT);

    bool is_new_file = false;
    sal_file = sfopen(hdr_filename, rw_mode, &is_new_file);
    if (!sal_file) {
	err = 1;
	return err;
    }


    char scan_type[MAXPATHLEN];

    ScanAxis axis = X_SCAN_AXIS;

    if (!is_new_file) {
	char buf[MAXPATHLEN];
	if (fgets(buf, MAXPATHLEN, sal_file)) {
	    if (!is_valid(buf)) {
		Log::logger()->error("'%s' is not a valid SAL file", filename.c_str());
		err = 1;
		return err;
	    }
	}

	while (fgets(buf, MAXPATHLEN, sal_file)) {
	    const char *buf1 = buf + 1;

	    if (Util::sstrncmp(buf1, "NXP")) {
		sscanf(strchr(buf, '=') + 1, " %d", &nx);
	    }
	    else if (Util::sstrncmp(buf1, "NYP")) {
		sscanf(strchr(buf, '=') + 1, " %d", &ny);
	    }
	    else if (Util::sstrncmp(buf1, "AXSCAN")) {
		char *t = strrchr(buf, '\'');
		if (t[-1] == 'Y') {
		    axis = Y_SCAN_AXIS;
		}
	    }
	    else if (Util::sstrncmp(buf1, "FILE REC LEN")) {
		sscanf(strchr(buf, '=') + 1, " %d", &record_length);
	    }
	    else if (Util::sstrncmp(buf1, "SCAN TYPE")) {
		sscanf(strchr(buf, '\'') + 1, " %s", scan_type);
		if (scan_type[0] == 'R') {
		    scan_mode = RASTER_SCAN;
		}
	    }
	    else if (Util::sstrncmp(buf1, "DELTAX")) {
		sscanf(strchr(buf, '=') + 1, " %f", &pixel);
		pixel /= 3.0;
	    }
	}


	if (axis == Y_SCAN_AXIS) {
	    int t = nx;
	    nx = ny;
	    ny = t;
	}
    }
    fclose(sal_file);
    sal_file = sfopen(img_filename, rw_mode);
    if (!sal_file) {
	err = 1;
	return err;
    }

    return 0;
}

bool SalIO::is_valid(const void *first_block)
{
    Log::logger()->log("SalIO::is_valid()");
    if (!first_block) {
	return false;
    }

    return Util::check_file_by_magic(first_block, MAGIC);
}

int SalIO::read_header(Dict & dict, int image_index, const Region * area,
		       bool is_3d)
{
    Log::logger()->log("SalIO::read_header() from file '%s'", filename.c_str());

    if (check_read_access(image_index) != 0) {
	return 1;
    }

    assert(!is_3d);

    if (check_region(area, Size(nx, ny)) != 0) {
	return 1;
    }
    int xlen = 0, ylen = 0;
    EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen);

    dict["nx"] = xlen;
    dict["ny"] = ylen;
    dict["nz"] = 1;
    dict["datatype"] = EMUtil::EM_SHORT;
    dict["pixel"] = pixel;

    return 0;
}

int SalIO::write_header(const Dict &, int)
{
    Log::logger()->log("SalIO::write_header() to file '%s'", filename.c_str());
    Log::logger()->warn("SAL write is not supported.");
    return 1;
}

int SalIO::read_data(float *data, int image_index, const Region * area, bool )
{
    Log::logger()->log("SalIO::read_data() from file '%s'", filename.c_str());

    if (check_read_access(image_index, true, data) != 0) {
	return 1;
    }
    if (check_region(area, Size(nx, ny)) != 0) {
	return 1;
    }

    assert(scan_mode == NON_RASTER_SCAN);

    rewind(sal_file);

    size_t mode_size = sizeof(short);
    unsigned char *cdata = (unsigned char *) data;
    short *sdata = (short *) data;
    size_t row_size = nx * mode_size;
    size_t block_size = (((row_size - 1) / record_length) + 1) * record_length;
    size_t post_row = block_size - row_size;

    int err = EMUtil::get_region_data(cdata, sal_file, image_index, mode_size,
				      nx, ny, 1, area, false, 0, post_row);
    if (err) {
	return 1;
    }

#if 0
    int row_size = nx * mode_size;
    int block_size = (((row_size - 1) / record_length) + 1) * record_length;

    for (int j = 0; j < ny; j++) {
	if (fread(&cdata[j * row_size], block_size, 1, sal_file) != 1) {
	    Log::logger()->error("Incomplete SAL data read %d/%d blocks", j, ny);
	    return 1;
	}
    }
#endif

    int xlen = 0, ylen = 0;
    EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen);

    if (scan_mode == NON_RASTER_SCAN) {
	become_platform_endian(sdata, xlen * ylen);

	for (int i = 0; i < ylen; i += 2) {
	    for (int j = 0; j < xlen / 2; j++) {
		short sw = sdata[i * xlen + j];
		sdata[i * xlen + j] = sdata[i * xlen + xlen - j - 1];
		sdata[i * xlen + xlen - j - 1] = sw;
	    }
	}
    }

    for (int i = xlen * ylen - 1; i >= 0; i--) {
	data[i] = static_cast<float>((cdata[i * 2 + 1] * UCHAR_MAX) + cdata[i * 2]);
    }
    return 0;
}

int SalIO::write_data(float *, int)
{
    Log::logger()->log("SalIO::write_data() to file '%s'", filename.c_str());
    Log::logger()->warn("SAL write is not supported.");
    return 1;
}


bool SalIO::is_complex_mode()
{
    return false;
}

bool SalIO::is_image_big_endian()
{
    return false;
}

int SalIO::get_nimg()
{
    if (init() != 0) {
	return 0;
    }

    return 1;
}
