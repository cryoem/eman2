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
#include "icosio.h"
#include "portable_fileio.h"
#include "geometry.h"

using namespace EMAN;


IcosIO::IcosIO(const string & file, IOMode rw)
:	filename(file), rw_mode(rw), icos_file(0), initialized(false)
{
	is_big_endian = ByteOrder::is_host_big_endian();
	is_new_file = false;
	memset(&icosh, 0, sizeof(IcosHeader));
}

IcosIO::~IcosIO()
{
	if (icos_file) {
		fclose(icos_file);
		icos_file = 0;
	}
}

void IcosIO::init()
{
	ENTERFUNC;
	if (initialized) {
		return ;
	}

	initialized = true;
	icos_file = sfopen(filename, rw_mode, &is_new_file);

	if (!is_new_file) {
		if (fread(&icosh, sizeof(IcosHeader), 1, icos_file) != 1) {
			throw ImageReadException(filename, "ICOS header");
		}

		if (!is_valid(&icosh)) {
			throw ImageReadException(filename, "invalid ICOS file");
		}
		is_big_endian = ByteOrder::is_data_big_endian(&icosh.stamp);
		become_host_endian((int *) &icosh, sizeof(IcosHeader) / sizeof(int));
	}
	EXITFUNC;
}

bool IcosIO::is_valid(const void *first_block)
{
	ENTERFUNC;
	bool result = false;
	if (!first_block) {
		result = false;
	}
	else {
		const int *data = static_cast < const int *>(first_block);
		int stamp = data[0];
		int stamp1 = data[19];
		int stamp2 = data[20];
		int stamp3 = data[26];

		bool data_big_endian = ByteOrder::is_data_big_endian(&stamp);

		if (data_big_endian != ByteOrder::is_host_big_endian()) {
			ByteOrder::swap_bytes(&stamp);
			ByteOrder::swap_bytes(&stamp1);
			ByteOrder::swap_bytes(&stamp2);
			ByteOrder::swap_bytes(&stamp3);
		}

		if (stamp == STAMP && stamp1 == STAMP1 && stamp2 == STAMP2 && stamp3 == STAMP3) {
			result = true;
		}
	}
	EXITFUNC;
	return result;
}

int IcosIO::read_header(Dict & dict, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	//single image format, index can only be zero
	image_index = 0;
	check_read_access(image_index);
	check_region(area, IntSize(icosh.nx, icosh.ny, icosh.nz),false,false);

	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, icosh.nx, &xlen, icosh.ny, &ylen, icosh.nz, &zlen);

	dict["nx"] = xlen;
	dict["ny"] = ylen;
	dict["nz"] = zlen;
	dict["datatype"] = EMUtil::EM_FLOAT;
	dict["minimum"] = icosh.min;
	dict["maximum"] = icosh.max;

	EXITFUNC;
	return 0;
}

int IcosIO::write_header(const Dict & dict, int image_index, const Region* ,
						 EMUtil::EMDataType, bool)
{
	ENTERFUNC;
	int err = 0;
	//single image format, index can only be zero
	image_index = 0;
	check_write_access(rw_mode, image_index, 1);

	icosh.stamp = STAMP;
	icosh.stamp1 = STAMP1;
	icosh.stamp2 = STAMP2;
	icosh.stamp3 = STAMP3;

	icosh.nx = dict["nx"];
	icosh.ny = dict["ny"];
	icosh.nz = dict["nz"];

	icosh.min = dict["minimum"];
	icosh.max = dict["maximum"];

	rewind(icos_file);

	if (fwrite(&icosh, sizeof(IcosHeader), 1, icos_file) != 1) {
		throw ImageWriteException(filename, "ICOS header write");
	}

	EXITFUNC;
	return err;
}

int IcosIO::read_data(float *data, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	//single image format, index can only be zero
	image_index = 0;
	check_read_access(image_index, data);
	check_region(area, IntSize(icosh.nx, icosh.ny, icosh.nz),false,false);

	portable_fseek(icos_file, sizeof(IcosHeader), SEEK_SET);

	EMUtil::process_region_io((unsigned char *) data, icos_file,
							  READ_ONLY, image_index,
							  sizeof(float), icosh.nx, icosh.ny, icosh.nz,
							  area, false, EMUtil::IMAGE_ICOS, sizeof(int), sizeof(int));

	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, icosh.nx, &xlen, icosh.ny, &ylen, icosh.nz, &zlen);
	become_host_endian(data, xlen * ylen * zlen);


	EXITFUNC;
	return 0;
}

int IcosIO::write_data(float *data, int image_index, const Region* area,
					   EMUtil::EMDataType, bool)
{
	ENTERFUNC;
	//single image format, index can only be zero
	image_index = 0;
	check_write_access(rw_mode, image_index, 1, data);
	portable_fseek(icos_file, sizeof(IcosHeader), SEEK_SET);

	EMUtil::process_region_io(data, icos_file, rw_mode, image_index,
							  sizeof(float), icosh.nx, icosh.ny, icosh.nz, area,
							  false, EMUtil::IMAGE_ICOS, sizeof(int), sizeof(int));

#if 0
	int float_size = sizeof(float);
	int nx = icosh.nx;
	float *buf = new float[nx + 2];
	buf[0] = (float)float_size * nx;
	buf[nx + 1] = buf[0];
	int nrows = icosh.ny * icosh.nz;

	int row_size = (nx + 2) * float_size;


	for (int j = 0; j < nrows; j++) {
		memcpy(&buf[1], &data[nx * j], nx * float_size);
		if (fwrite(buf, row_size, 1, icos_file) != 1) {
			throw ImageWriteException(filename, "ICOS data");
		}
	}

	if( buf )
	{
		delete[]buf;
		buf = 0;
	}
#endif
	EXITFUNC;
	return 0;
}

void IcosIO::flush()
{
	fflush(icos_file);
}

bool IcosIO::is_complex_mode()
{
	return false;
}

bool IcosIO::is_image_big_endian()
{
	init();
	return is_big_endian;
}

