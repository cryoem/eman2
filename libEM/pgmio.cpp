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
:	filename(file), pgm_file(0), initialized(false)
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
		buf[i] = static_cast < char >(c);
		i++;
	}

	return atoi(buf);
}

int PgmIO::init()
{
	ENTERFUNC;
	
	static int err = 0;
	if (initialized) {
		return err;
	}

	initialized = true;

	bool is_new_file = false;
	pgm_file = sfopen(filename, rw_mode, &is_new_file, true);

	if (!is_new_file) {
		const int bufsz = 1024;
		char buf[bufsz];

		buf[0] = static_cast < char >(getc(pgm_file));
		buf[1] = static_cast < char >(getc(pgm_file));
		buf[2] = '\0';
		getc(pgm_file);

		if (!is_valid(&buf)) {
			LOGERR("not a valid PGM file");
			err = 1;
			return 1;
		}

		char c = '\0';

		while ((c = static_cast < char >(getc(pgm_file))) == '#') {
			fgets(buf, bufsz, pgm_file);
		}
		ungetc(c, pgm_file);

		nx = read_int_and_space(pgm_file);
		ny = read_int_and_space(pgm_file);
		maxval = read_int_and_space(pgm_file);

		if (nx <= 0 || ny <= 0) {
			LOGERR("invalid file size: %dx%d", nx, ny);
			err = 1;
			return 1;
		}
		if (maxval > USHRT_MAX) {
			LOGERR("not a valid PGM file: max gray value '%d' cannot > $d",
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
			LOGERR("not a valid PGM file. max gray value '%d' cannot <= 0", maxval);
			err = 1;
			return 1;
		}
		file_offset = portable_ftell(pgm_file);
	}
	EXITFUNC;
	return 0;

}

bool PgmIO::is_valid(const void *first_block)
{
	ENTERFUNC;
	bool result = Util::check_file_by_magic(first_block, MAGIC_BINARY);
	EXITFUNC;
	return result;
}

int PgmIO::read_header(Dict & dict, int image_index, const Region * area, bool)
{
	ENTERFUNC;
	int err = 0;
	
	if (check_read_access(image_index) != 0) {
	    err = 1;
	}
	else {
		if (check_region(area, IntSize(nx, ny)) != 0) {
			err = 1;
		}
		else {
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
		}
	}
	EXITFUNC;
	return err;
}

int PgmIO::write_header(const Dict & dict, int image_index, const Region*, bool)
{
	ENTERFUNC;
	int err = 0;
	
	if (check_write_access(rw_mode, image_index) != 0) {
		err = 1;
	}
	else {
		int nz = dict["nz"];
		if (nz != 1) {
			LOGERR("Cannot write 3D image as PGM. Your image nz = %d", nz);
			err = 1;
		}
		else {
			rewind(pgm_file);
			
			nx = dict["nx"];
			ny = dict["ny"];
			minval = dict["min_gray"];
			maxval = dict["max_gray"];
			fprintf(pgm_file, "%s\n%d %d\n%d\n", MAGIC_BINARY, nx, ny, maxval);
		}
	}
	EXITFUNC;
	return err;
}

int PgmIO::read_data(float *data, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	if (check_read_access(image_index, data) != 0) {
		return 1;
	}

	if (check_region(area, IntSize(nx, ny)) != 0) {
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

	EMUtil::process_region_io(cdata, pgm_file, READ_ONLY, image_index, 
							  mode_size, nx, ny, 1, area, true);

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
			LOGERR("Incomplete data read in PGM file '%s'", filename.c_str());
			return 1;
		}
	}
#endif

	int xlen = 0, ylen = 0;
	EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen);

	if (datatype == PGM_USHORT) {
		become_host_endian(sdata, xlen * ylen);
	}

	for (int k = xlen * ylen - 1; k >= 0; k--) {
		if (datatype == PGM_UCHAR) {
			data[k] = static_cast < float >(cdata[k]);
		}
		else {
			data[k] = static_cast < float >(sdata[k]);
		}
	}
	EXITFUNC;
	return 0;
}

int PgmIO::write_data(float *data, int image_index, const Region* area, bool)
{
	ENTERFUNC;

	if (check_write_access(rw_mode, image_index, 1, data) != 0) {
		return 1;
	}
	portable_fseek(pgm_file, file_offset, SEEK_SET);

	LOGERR("not working yet. need to normalize data before write");

	//fwrite(data, nx, ny, pgm_file);
	EXITFUNC;
	return 1;
}

void PgmIO::flush()
{
	fflush(pgm_file);
}


bool PgmIO::is_complex_mode()
{
	return false;
}

bool PgmIO::is_image_big_endian()
{
	return is_big_endian;
}

