#include "io_template.h"
#include "log.h"
#include "portable_fileio.h"

using namespace EMAN;

XYZIO::XYZIO(string file, IOMode rw)
:	filename(file), rw_mode(rw), xyz_file(0), initialized(false)
{
	is_big_endian = ByteOrder::is_host_big_endian();
}

XYZIO::~XYZIO()
{
	if (xyz_file) {
		fclose(xyz_file);
		xyz_file = 0;
	}
}

int XYZIO::init()
{
	static int err = 0;
	if (initialized) {
		return err;
	}
	LOGDEBUG("XYZIO::init()");
	initialized = true;

	bool is_new_file = false;
	xyz_file = sfopen(filename, rw_mode, &is_new_file);

	if (!xyz_file) {
		err = 1;
		return err;
	}

	if (!is_new_file) {

	}

	return 0;
}

bool XYZIO::is_valid(const void *first_block)
{
	LOGDEBUG("XYZIO::is_valid()");
	if (!first_block) {
		return false;
	}
	return false;
}

int XYZIO::read_header(Dict & dict, int image_index, const Region * area, bool is_3d)
{
	LOGDEBUG("XYZIO::read_header() from file '%s'", filename.c_str());

	if (check_read_access(image_index) != 0) {
		return 1;
	}

	return 0;
}

int XYZIO::write_header(const Dict & dict, int image_index, const Region *region, bool)
{
	LOGDEBUG("XYZIO::write_header() to file '%s'", filename.c_str());
	if (check_write_access(rw_mode, image_index) != 0) {
		return 1;
	}

	return 0;
}

int XYZIO::read_data(float *data, int image_index, const Region * area, bool is_3d)
{
	LOGDEBUG("XYZIO::read_data() from file '%s'", filename.c_str());

	if (check_read_access(image_index, true, data) != 0) {
		return 1;
	}

	return 0;
}

int XYZIO::write_data(float *data, int image_index, const Region *region, bool)
{
	LOGDEBUG("XYZIO::write_data() to file '%s'", filename.c_str());
	if (check_write_access(rw_mode, image_index, true, data) != 0) {
		return 1;
	}
	return 0;
}


bool XYZIO::is_complex_mode()
{
	return false;
}

bool XYZIO::is_image_big_endian()
{
	return is_big_endian;
}

int XYZIO::get_nimg()
{
	if (init() != 0) {
		return 0;
	}

	return 1;
}
