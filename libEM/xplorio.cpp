/**
 * $Id$
 */
#include "xplorio.h"
#include "log.h"
#ifdef WIN32
#include <time.h>
#endif


using namespace EMAN;

XplorIO::XplorIO(string file, IOMode rw)
:	filename(file), rw_mode(rw), xplor_file(0), initialized(false)
{
	is_big_endian = ByteOrder::is_host_big_endian();
	nx = 0;
	ny = 0;
	nz = 0;
}

XplorIO::~XplorIO()
{
	if (xplor_file) {
		fclose(xplor_file);
		xplor_file = 0;
	}
}

int XplorIO::init()
{
	static int err = 0;
	if (initialized) {
		return err;
	}
	Log::logger()->log("XplorIO::init()");
	initialized = true;

	bool is_new_file = false;
	xplor_file = sfopen(filename, rw_mode, &is_new_file);

	if (!xplor_file) {
		err = 1;
		return err;
	}


	return 0;
}

bool XplorIO::is_valid(const void *first_block)
{
	Log::logger()->log("XplorIO::is_valid()");
	if (!first_block) {
		return false;
	}
	return false;
}

int XplorIO::read_header(Dict &, int, const Region *, bool)
{
	Log::logger()->log("XplorIO::read_header() from file '%s'", filename.c_str());
	Log::logger()->warn("XPLOR read is not supported.");
	return 1;
}

int XplorIO::write_header(const Dict & dict, int image_index, bool)
{
	Log::logger()->log("XplorIO::write_header() to file '%s'", filename.c_str());
	if (check_write_access(rw_mode, image_index) != 0) {
		return 1;
	}

	nx = dict["nx"];
	ny = dict["ny"];
	nz = dict["nz"];
	float pixel = dict["pixel"];

	time_t t0 = time(0);
	struct tm *t = localtime(&t0);

	fprintf(xplor_file, "\n%8d\n\"%s\" written by EMAN at %s", 1, filename.c_str(), asctime(t));

	int z0 = -nz / 2;
	int z1 = (nz - 1) / 2;

	if (2 * nz - 1 == nx && 2 * nz - 1 == ny) {
		z0 = 0;
		z1 = nz - 1;
	}

	fprintf(xplor_file, "%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
			nx, -nx / 2, nx % 2 ? nx / 2 : nx / 2 - 1, ny, -ny / 2,
			ny % 2 ? ny / 2 : ny / 2 - 1, nz, z0, z1);
	fprintf(xplor_file, "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\nZYX\n", nx * pixel, ny * pixel,
			nz * pixel, 90.0, 90.0, 90.0);
	return 0;
}

int XplorIO::read_data(float *, int, const Region *, bool)
{
	Log::logger()->log("XplorIO::read_data() from file '%s'", filename.c_str());
	Log::logger()->warn("XPLOR read is not supported.");
	return 1;
}

int XplorIO::write_data(float *data, int image_index, bool)
{
	Log::logger()->log("XplorIO::write_data() to file '%s'", filename.c_str());
	if (check_write_access(rw_mode, image_index, true, data) != 0) {
		return 1;
	}

	int nsecs = nx * ny;
	int step = 6;

	for (int k = 0; k < nz; k++) {
		fprintf(xplor_file, "%8d\n", k);

		for (int i = 0; i < nsecs - step; i += step) {
			for (int j = 0; j < step; j++) {
				fprintf(xplor_file, "%12.5E", data[k * nsecs + i + j]);
			}
			fprintf(xplor_file, "\n");
		}

		for (int l = (nsecs - 1) / step * step; l < nsecs; l++) {
			fprintf(xplor_file, "%12.5E", data[k * nsecs + l]);
		}

		fprintf(xplor_file, "\n");
	}

	fprintf(xplor_file, "%8d\n", -9999);

	return 0;
}

bool XplorIO::is_complex_mode()
{
	return false;
}

bool XplorIO::is_image_big_endian()
{
	return is_big_endian;
}

int XplorIO::get_nimg()
{
	if (init() != 0) {
		return 0;
	}

	return 1;
}
