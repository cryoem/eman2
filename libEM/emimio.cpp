/**
 * $Id$
 */

/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
 * Copyright (c) 2000-2006 Baylor College of Medicine
 *
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * */

#include <cstring>
#include "emimio.h"
#include "portable_fileio.h"

using namespace EMAN;

const char *EmimIO::MAGIC = "EMIM";

EmimIO::EmimIO(const string & file, IOMode rw)
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

	//single image format, index can only be zero
	if(image_index == -1) {
		image_index = 0;
	}
	image_index = 0;
	check_read_access(image_index);

	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, efh.nx, &xlen, efh.ny, &ylen, efh.nz, &zlen);

	dict["nx"] = xlen;
	dict["ny"] = ylen;
	dict["nz"] = zlen;

	dict["datatype"] = EMUtil::EM_FLOAT;
	dict["pixel"] = efh.pixel;

	off_t imgsize = (off_t)efh.nx * (off_t)efh.ny * (off_t)efh.nz * (off_t)sizeof(float) + (off_t)sizeof(EmimImageHeader);
	off_t offset = (off_t)sizeof(EmimFileHeader) + imgsize * (off_t)image_index;

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

int EmimIO::write_header(const Dict &, int, const Region* , EMUtil::EMDataType, bool)
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
	//single image format, index can only be zero
	image_index = 0;
	check_read_access(image_index, data);

	off_t imgsize = (off_t)efh.nx * (off_t)efh.ny * (off_t)efh.nz * (off_t)sizeof(float) + (off_t)sizeof(EmimImageHeader);
	off_t offset = (off_t)sizeof(EmimFileHeader) + imgsize * (off_t)(image_index+1);
	portable_fseek(emim_file, offset, SEEK_SET);

	unsigned char *cdata = (unsigned char *) data;
	EMUtil::process_region_io(cdata, emim_file, READ_ONLY, 0, sizeof(float),
							  efh.nx, efh.ny, efh.nz, area);

	become_host_endian(data, efh.nx * efh.ny * efh.nz);


	EXITFUNC;
	return err;
}

int EmimIO::write_data(float *, int, const Region* , EMUtil::EMDataType, bool)
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
	init();
	if (efh.flag & EMIM_COMPLEX) {
		return true;
	}
	return false;
}

bool EmimIO::is_image_big_endian()
{
	init();
	return is_big_endian;
}

int EmimIO::get_nimg()
{
	init();
	return efh.count;
}
