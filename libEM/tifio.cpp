/**
 * $Id$
 */
#ifdef EM_TIFF

#include "tifio.h"
#include "log.h"
#include "emutil.h"
#include "util.h"
#include "geometry.h"
#include <tiffio.h>
#include <limits.h>
#include <assert.h>

using namespace EMAN;

TiffIO::TiffIO(string tiff_filename, IOMode rw)
:	filename(tiff_filename), rw_mode(rw), tiff_file(0), bitspersample(0), initialized(false)
{
	is_big_endian = ByteOrder::is_host_big_endian();
}


TiffIO::~TiffIO()
{
	if (tiff_file) {
		TIFFClose(tiff_file);
	}
}

int TiffIO::init()
{
	static int err = 0;
	if (initialized) {
		return err;
	}
	//initialized = true;

	LOGDEBUG("TiffIO::init()");

	if (rw_mode != READ_ONLY) {
		LOGERR("wrong rw mode. Only reading is supported for TIFF.");
		err = 1;
		return err;
	}

	FILE *tmp_in = fopen(filename.c_str(), "rb");
	if (!tmp_in) {
		LOGERR("cannot open TIFF image '%s'", filename.c_str());
		err = 1;
		return err;
	}
	char buf[64];
	if (fread(buf, sizeof(buf), 1, tmp_in) != 1) {
		LOGERR("cannot read TIFF image '%s'", filename.c_str());
		err = 1;
	}

	if (!err && !is_valid(&buf)) {
		LOGERR("'%s' is an invalid TIFF image", filename.c_str());
		err = 1;
	}

	fclose(tmp_in);
	tmp_in = 0;

	if (err) {
		return err;
	}

	TIFFSetWarningHandler(0);

	tiff_file = TIFFOpen(filename.c_str(), "r");
	if (!tiff_file) {
		LOGERR("cannot open TIFF image '%s'", filename.c_str());
		err = 1;
		return err;
	}

	if (buf[0] == TIFF_BIG_ENDIAN) {
		is_big_endian = true;
	}
	else {
		is_big_endian = false;
	}

	TIFFGetField(tiff_file, TIFFTAG_BITSPERSAMPLE, &bitspersample);

	if (bitspersample != CHAR_BIT && bitspersample != (CHAR_BIT * sizeof(short))) {
		LOGERR("invalid %d bits. only %d-bit and %d-bit TIFF are supported",
							 bitspersample, CHAR_BIT, (CHAR_BIT * sizeof(short)));
		err = 1;
	}

	return err;
}

bool TiffIO::is_valid(const void *first_block)
{
	LOGDEBUG("TiffIO::is_valid()");

	if (!first_block) {
		return false;
	}

	const char *data = static_cast < const char *>(first_block);

	if ((data[0] == data[1]) && (data[0] == TIFF_LITTLE_ENDIAN || data[1] == TIFF_BIG_ENDIAN)) {
		return true;
	}
	return false;
}

int TiffIO::read_header(Dict & dict, int img_index, const Region * area, bool)
{
	LOGDEBUG("TiffIO::read_header() on file '%s'", filename.c_str());
	if (check_read_access(img_index) != 0) {
		return 1;
	}

	int nx = 0;
	int ny = 0;
	TIFFGetField(tiff_file, TIFFTAG_IMAGEWIDTH, &nx);
	TIFFGetField(tiff_file, TIFFTAG_IMAGELENGTH, &ny);

	if (check_region(area, Size(nx, ny)) != 0) {
		return 1;
	}

	float min = 0;
	float max = 0;
	TIFFDataType data_type = TIFF_NOTYPE;
	float resolution_x = 0;
	float resolution_y = 0;

	TIFFGetField(tiff_file, TIFFTAG_MINSAMPLEVALUE, &min);
	TIFFGetField(tiff_file, TIFFTAG_MAXSAMPLEVALUE, &max);

	TIFFGetField(tiff_file, TIFFTAG_SAMPLEFORMAT, &data_type);
	TIFFGetField(tiff_file, TIFFTAG_XRESOLUTION, &resolution_x);
	TIFFGetField(tiff_file, TIFFTAG_YRESOLUTION, &resolution_y);

	int xlen = 0, ylen = 0;
	EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen);

	dict["nx"] = xlen;
	dict["ny"] = ylen;
	dict["nz"] = 1;

	dict["minimum"] = min;
	dict["maximum"] = max;

	if (bitspersample == CHAR_BIT) {
		dict["datatype"] = EMUtil::EM_UCHAR;
	}
	else if (bitspersample == sizeof(unsigned short) * CHAR_BIT) {
		dict["datatype"] = EMUtil::EM_USHORT;
	}

	dict["bitspersample"] = bitspersample;
	dict["resolution_x"] = resolution_x;
	dict["resolution_y"] = resolution_y;

	return 0;
}

int TiffIO::read_data(float *rdata, int img_index, const Region * area, bool is_3d)
{
	LOGDEBUG("TiffIO::read_data() on file '%s'", filename.c_str());

	if (check_read_access(img_index, true, rdata) != 0) {
		return 1;
	}
	assert(!is_3d);

	int nx = 0;
	int ny = 0;
	TIFFGetField(tiff_file, TIFFTAG_IMAGEWIDTH, &nx);
	TIFFGetField(tiff_file, TIFFTAG_IMAGELENGTH, &ny);

	if (check_region(area, Size(nx, ny)) != 0) {
		return 1;
	}

	int xlen = 0, ylen = 0, x0 = 0, y0 = 0;
	EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen);
	EMUtil::get_region_origins(area, &x0, &y0);

	int err = 0;
	int strip_size = TIFFStripSize(tiff_file);
	uint32 num_strips = TIFFNumberOfStrips(tiff_file);

	unsigned char *cdata = static_cast < unsigned char *>(_TIFFmalloc(strip_size));

	int k = 0;
	int num_read = 0;
	int mode_size = bitspersample / CHAR_BIT;
	int total_rows = 0;

	for (uint32 i = 0; i < num_strips; i++) {
		if ((num_read = TIFFReadEncodedStrip(tiff_file, i, cdata, strip_size)) == -1) {
			LOGERR("reading stripped TiFF image '%s' failed", filename.c_str());
			err = 1;
			break;
		}

		int nitems = num_read / mode_size;
		assert(nitems % nx == 0);
		int nrows = nitems / nx;
		total_rows += nrows;

		int y_start = 0;
		int y_end = nrows;

		if (area) {
			if (total_rows >= y0 && total_rows < y0 + nrows) {
				y_start = nrows - (total_rows - y0);
			}
			else if (total_rows >= (y0 + ylen) && total_rows < (y0 + ylen + nrows)) {
				y_end = y0 + ylen - total_rows + nrows;
			}
			else if (total_rows >= (y0 + ylen + nrows)) {
				break;
			}
		}

		for (int l = y_start; l < y_end; l++) {
			for (int j = x0; j < x0 + xlen; j++) {
				if (bitspersample == CHAR_BIT) {
					rdata[k] = static_cast < float >(cdata[l * nx + j]);
				}
				else if (bitspersample == sizeof(unsigned short) * CHAR_BIT) {
					rdata[k] = static_cast < float >(((unsigned short *) cdata)[l * nx + j]);
				}
				k++;
			}
		}
	}

	Util::flip_image(rdata, xlen, ylen);

	_TIFFfree(cdata);

	return err;
}


int TiffIO::write_header(const Dict &, int, bool)
{
	LOGDEBUG("TiffIO::write_header()");
	LOGERR("TIFF write is not supported");
	return 0;
}

int TiffIO::write_data(float *, int, bool)
{
	LOGDEBUG("TiffIO::write_data()");
	LOGERR("TIFF write is not supported");
	return 0;
}



bool TiffIO::is_complex_mode()
{
	return false;
}

bool TiffIO::is_image_big_endian()
{
	init();
	return is_big_endian;
}

int TiffIO::get_nimg()
{
	if (init() != 0) {
		return 0;
	}

	return 1;
}

#endif
