/**
 * $Id$
 */
#include "emimio.h"
#include "log.h"
#include "emutil.h"
#include "portable_fileio.h"

using namespace EMAN;

const char *EmimIO::MAGIC = "EMIM";

EmimIO::EmimIO(string file, IOMode rw)
:	filename(file), rw_mode(rw), emim_file(0), initialized(false)
{
	is_big_endian = ByteOrder::is_host_big_endian();
	memset(&efh, 0, sizeof(EmimFileHeader));
}

EmimIO::~EmimIO()
{
	if (emim_file) {
		fclose(emim_file);
		emim_file = 0;
	}
}

void EmimIO::init()
{
	ENTERFUNC;
	if (initialized) {
		return;
	}
	
	initialized = true;
	bool is_new_file = false;
	emim_file = sfopen(filename, rw_mode, &is_new_file);

	if (!is_new_file) {
		if (fread(&efh, sizeof(EmimFileHeader), 1, emim_file) != 1) {
			throw ImageReadException(filename, "EMIM header");
		}

		if (!is_valid(&efh)) {
			throw ImageReadException(filename, "invalid EMIM file");
		}

		become_host_endian((int *) &efh, NUM_INT_IN_FILE_HEADER);
		is_big_endian = ByteOrder::is_data_big_endian(&efh.count);
	}

	EXITFUNC;
}

bool EmimIO::is_valid(const void *first_block)
{
	ENTERFUNC;
	
	if (!first_block) {
		return false;
	}

	const char *data = static_cast < const char *>(first_block);
	const int *idata = static_cast < const int *>(first_block);
	int count = idata[2];

	if (strncmp(data, MAGIC, sizeof(MAGIC)) == 0) {
		bool data_big_endian = ByteOrder::is_data_big_endian(&count);

		if (data_big_endian != ByteOrder::is_host_big_endian()) {
			ByteOrder::swap_bytes(&count);
		}

		if (count >= 0 && count <= 1 << 20) {
			return true;
		}
	}
	EXITFUNC;
	return false;
}

int EmimIO::read_header(Dict & dict, int image_index, const Region * area, bool)
{
	ENTERFUNC;
	
	check_read_access(image_index);

	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, efh.nx, &xlen, efh.ny, &ylen, efh.nz, &zlen);

	dict["nx"] = xlen;
	dict["ny"] = ylen;
	dict["nz"] = zlen;

	dict["datatype"] = EMUtil::EM_FLOAT;
	dict["pixel"] = efh.pixel;

	off_t imgsize = efh.nx * efh.ny * efh.nz * sizeof(float) + sizeof(EmimImageHeader);
	off_t offset = sizeof(EmimFileHeader) + imgsize * image_index;

	portable_fseek(emim_file, offset, SEEK_SET);

	EmimImageHeader eih;
	fread(&eih, sizeof(EmimImageHeader), 1, emim_file);

	int n = eih.mgnum;
	become_host_endian(&n);

	char mgnum[32];
	sprintf(mgnum, "%d", n);

	dict["micrograph_id"] = mgnum;
	EXITFUNC;
	return 0;

}

int EmimIO::write_header(const Dict &, int, const Region* , bool)
{
	ENTERFUNC;
	LOGWARN("EMIM write header is not supported.");
	EXITFUNC;
	return 1;
}

int EmimIO::read_data(float *data, int image_index, const Region * area, bool)
{
	ENTERFUNC;
	int err = 0;
	check_read_access(image_index, data);
	
	off_t imgsize = efh.nx * efh.ny * efh.nz * sizeof(float) + sizeof(EmimImageHeader);
	off_t offset = sizeof(EmimFileHeader) + imgsize * (image_index+1);
	portable_fseek(emim_file, offset, SEEK_SET);
		
	unsigned char *cdata = (unsigned char *) data;
	EMUtil::process_region_io(cdata, emim_file, READ_ONLY, 0, sizeof(float),
							  efh.nx, efh.ny, efh.nz, area);
		
	become_host_endian(data, efh.nx * efh.ny * efh.nz);
		
	
	EXITFUNC;
	return err;
}

int EmimIO::write_data(float *, int, const Region* , bool)
{
	ENTERFUNC;
	LOGWARN("EMIM write data is not supported.");
	EXITFUNC;
	return 1;
}

void EmimIO::flush()
{
}

bool EmimIO::is_complex_mode()
{
	if (efh.flag & EMIM_COMPLEX) {
		return true;
	}
	return false;
}

bool EmimIO::is_image_big_endian()
{
	return is_big_endian;
}

int EmimIO::get_nimg()
{
	init();
	return efh.count;
}
