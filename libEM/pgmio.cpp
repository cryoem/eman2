/**
 * $Id$
 */
#include "pgmio.h"
#include "log.h"
#include "geometry.h"
#include "util.h"
#include "emutil.h"
#include "portable_fileio.h"
#include <limits.h>
#include <ctype.h>

using namespace EMAN;

const char *PgmIO::MAGIC_BINARY = "P5";
const char *PgmIO::MAGIC_ASCII = "P2";

PgmIO::PgmIO(string file, IOMode rw)
    :  filename(file), pgm_file(0), initialized(false)
{
    rw_mode = rw;
    is_big_endian = true;

    nx = 0;
    ny = 0;
    maxval = 0;
    minval = 0;
    datatype = PGM_UNKNOWN_TYPE;
    file_offset = 0;
}

PgmIO::~PgmIO()
{
    if (pgm_file) {
	fclose(pgm_file);
	pgm_file = 0;
    }
}

static int read_int_and_space(FILE * in)
{
    char buf[32];
    int c = 0;

    int i = 0;
    while (!isspace(c = getc(in))) {
	buf[i] = static_cast<char>(c);
	i++;
    }

    return atoi(buf);
}

int PgmIO::init()
{
    static int err = 0;
    if (initialized) {
	return err;
    }
    Log::logger()->log("PgmIO::init()");
    initialized = true;

    bool is_new_file = false;
    pgm_file = sfopen(filename, rw_mode, &is_new_file, true);

    if (!is_new_file) {
	const int bufsz = 1024;
	char buf[bufsz];

	buf[0] = static_cast<char>(getc(pgm_file));
	buf[1] = static_cast<char>(getc(pgm_file));
	buf[2] = '\0';
	getc(pgm_file);

	if (!is_valid(&buf)) {
	    Log::logger()->error("not a valid PGM file");
	    err = 1;
	    return 1;
	}

	char c = '\0';

	while ((c = static_cast<char>(getc(pgm_file))) == '#') {
	    fgets(buf, bufsz, pgm_file);
	}
	ungetc(c, pgm_file);

	nx = read_int_and_space(pgm_file);
	ny = read_int_and_space(pgm_file);
	maxval = read_int_and_space(pgm_file);

	if (nx <= 0 || ny <= 0) {
	    Log::logger()->error("invalid file size: %dx%d", nx, ny);
	    err = 1;
	    return 1;
	}
	if (maxval > USHRT_MAX) {
	    Log::logger()->error("not a valid PGM file: max gray value '%d' cannot > $d",
				 maxval, USHRT_MAX);
	    err = 1;
	    return 1;
	}
	else if (maxval > UCHAR_MAX) {
	    datatype = PGM_USHORT;
	}
	else if (maxval > 0) {
	    datatype = PGM_UCHAR;
	}
	else {
	    Log::logger()->error("not a valid PGM file. max gray value '%d' cannot <= 0", maxval);
	    err = 1;
	    return 1;
	}
	file_offset = portable_ftell(pgm_file);
    }

    return 0;

}

bool PgmIO::is_valid(const void *first_block)
{
    Log::logger()->log("PgmIO::is_valid()");
    return Util::check_file_by_magic(first_block, MAGIC_BINARY);
}

int PgmIO::read_header(Dict & dict, int image_index, const Region * area, bool)
{
    Log::logger()->log("PgmIO::read_header() from file '%s'", filename.c_str());

    if (check_read_access(image_index) != 0) {
	return 1;
    }

    if (check_region(area, Size(nx, ny)) != 0) {
	return 1;
    }

    int xlen = 0, ylen = 0;
    EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen);

    dict["nx"] = xlen;
    dict["ny"] = ylen;
    dict["nz"] = 1;

    if (datatype == PGM_UCHAR) {
	dict["datatype"] = EMUtil::EM_UCHAR;
    }
    else {
	dict["datatype"] = EMUtil::EM_USHORT;
    }
    
    dict["max_gray"] = maxval;
    dict["min_gray"] = minval;

    return 0;
}

int PgmIO::write_header(const Dict & dict, int image_index)
{
    Log::logger()->log("PgmIO::write_header() to file '%s'", filename.c_str());
    if (check_write_access(rw_mode, image_index) != 0) {
	return 1;
    }

    nx = dict["nx"];
    ny = dict["ny"];

    minval = dict["min_gray"];
    maxval = dict["max_gray"];

    int nz = dict["nz"];
    if (nz != 1) {
	Log::logger()->error("Cannot write 3d image as PGM. Your image nz = %d", nz);
	return 1;
    }

    fprintf(pgm_file, "%s\n%d %d\n%d\n", MAGIC_BINARY, nx, ny, maxval);

    return 0;
}

int PgmIO::read_data(float *data, int image_index, const Region * area, bool)
{
    Log::logger()->log("PgmIO::read_data() from file '%s'", filename.c_str());

    if (check_read_access(image_index, true, data) != 0) {
	return 1;
    }

    if (check_region(area, Size(nx, ny)) != 0) {
	return 1;
    }

    portable_fseek(pgm_file, file_offset, SEEK_SET);

    unsigned char *cdata = (unsigned char *) (data);
    unsigned short *sdata = (unsigned short *) (data);

    size_t mode_size = 0;
    if (datatype == PGM_UCHAR) {
	mode_size = sizeof(unsigned char);
    }
    else if (datatype == PGM_USHORT) {
	mode_size = sizeof(unsigned short);
    }

    int err = EMUtil::get_region_data(cdata, pgm_file, image_index, mode_size,
				      nx, ny, 1, area, true);
    if (err) {
	return 1;
    }

#if 0
    int ushort_size = sizeof(unsigned short);

    for (int j = 0; j < ny; j++) {
	int n = 0;
	if (datatype == PGM_UCHAR) {
	    n = fread(&cdata[(ny - j - 1) * nx], nx, 1, pgm_file);
	}
	else {
	    n = fread(&sdata[(ny - j - 1) * nx], nx * ushort_size, 1, pgm_file);
	}

	if (n != 1) {
	    Log::logger()->error("Incomplete data read in PGM file '%s'", filename.c_str());
	    return 1;
	}
    }
#endif

    int xlen = 0, ylen = 0;
    EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen);

    if (datatype == PGM_USHORT) {
	become_platform_endian(sdata, xlen * ylen);
    }

    for (int k = xlen * ylen - 1; k >= 0; k--) {
	if (datatype == PGM_UCHAR) {
	    data[k] = static_cast<float>(cdata[k]);
	}
	else {
	    data[k] = static_cast<float>(sdata[k]);
	}
    }

    return 0;
}

int PgmIO::write_data(float *data, int image_index)
{
    Log::logger()->log("PgmIO::write_data() to file '%s'", filename.c_str());
    if (check_write_access(rw_mode, image_index, true, data) != 0) {
	return 1;
    }
    portable_fseek(pgm_file, file_offset, SEEK_SET);

    Log::logger()->error("not working yet. need to normalize data before write");

    //fwrite(data, nx, ny, pgm_file);

    return 0;
}


bool PgmIO::is_complex_mode()
{
    return false;
}

bool PgmIO::is_image_big_endian()
{
    return is_big_endian;
}

int PgmIO::get_nimg()
{
    if (init() != 0) {
	return 0;
    }

    return 1;
}
