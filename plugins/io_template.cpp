#include "io_template.h"
#include "portable_fileio.h"

using namespace EMAN;

XYZIO::XYZIO(const string & file, IOMode rw)
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

void XYZIO::init()
{
	if (initialized) {
		return ;
	}

	ENTERFUNC;
	
	initialized = true;
	bool is_new_file = false;
	xyz_file = sfopen(filename, rw_mode, &is_new_file);

	if (!is_new_file) {

	}

	EXITFUNC;
}

bool XYZIO::is_valid(const void *first_block)
{
	ENTERFUNC;
	bool result = false;
	if (!first_block) {
		result = false;
	}

	// check image format validality here
	
	EXITFUNC;
	return result;
}

int XYZIO::read_header(Dict & , int image_index, const Region * , bool )
{
	ENTERFUNC;
	check_read_access(image_index);

	// read header info here
	
	EXITFUNC;
	return 0;
}

int XYZIO::write_header(const Dict & , int image_index, const Region *,
						EMUtil::EMDataType , bool)
{
	ENTERFUNC;
	check_write_access(rw_mode, image_index);
	// write header info here
	EXITFUNC;
	return 0;
}

int XYZIO::read_data(float *data, int image_index, const Region * , bool )
{
	ENTERFUNC;
	check_read_access(image_index, data);

	// read image data here

	EXITFUNC;
	return 0;
}

int XYZIO::write_data(float *data, int image_index, const Region *,
					  EMUtil::EMDataType , bool)
{
	ENTERFUNC;
	check_write_access(rw_mode, image_index, 0, data);

	// write image data here
	
	EXITFUNC;
	return 0;
}

void XYZIO::flush()
{
	if (xyz_file) {
		fflush(xyz_file);
	}
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
	init();

	return 1;
}
