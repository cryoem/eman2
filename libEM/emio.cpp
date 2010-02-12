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
#include "emio.h"
#include "portable_fileio.h"
#include "geometry.h"

using namespace EMAN;

EmIO::EmIO(const string & file, IOMode rw)
:	filename(file), rw_mode(rw), em_file(0), initialized(false)
{
	mode_size = 0;
	mode = EM_EM_UNKNOWN;
	is_big_endian = ByteOrder::is_host_big_endian();
	is_new_file = false;
	memset(&emh, 0, sizeof(EMHeader));
}

EmIO::~EmIO()
{
	if (em_file) {
		fclose(em_file);
		em_file = 0;
	}
}

void EmIO::init()
{
	ENTERFUNC;

	if (initialized) {
		return;
	}


	initialized = true;
	em_file = sfopen(filename, rw_mode, &is_new_file);

	if (!is_new_file) {
		if (fread(&emh, sizeof(EMHeader), 1, em_file) != 1) {
			throw ImageReadException(filename, "EM header");
		}
		if (!is_valid(&emh)) {
			throw ImageReadException(filename, "invalid EM image");
		}

		is_big_endian = ByteOrder::is_data_big_endian(&emh.nz);
		become_host_endian(&emh.nx);
		become_host_endian(&emh.ny);
		become_host_endian(&emh.nz);

		mode = (DataType) emh.data_type;

		if (mode == EM_EM_DOUBLE) {
			throw ImageReadException(filename, "DOUBLE data type not supported for EM image");
		}

		mode_size = get_mode_size(emh.data_type);
		if (is_complex_mode()) {
			emh.nx *= 2;
		}
	}
	EXITFUNC;
}

bool EmIO::is_valid(const void *first_block, off_t file_size)
{
	ENTERFUNC;

	if (!first_block) {
		return false;
	}

	const char *data = static_cast < const char *>(first_block);
	char machine = data[0];
	char is_new_ver = data[1];
	char data_type = data[3];

	const int *data1 = static_cast < const int *>(first_block);
	int nx = data1[1];
	int ny = data1[2];
	int nz = data1[3];

	bool data_big_endian = ByteOrder::is_data_big_endian(&nz);
	if (data_big_endian != ByteOrder::is_host_big_endian()) {
		ByteOrder::swap_bytes(&nx);
		ByteOrder::swap_bytes(&ny);
		ByteOrder::swap_bytes(&nz);
	}

	const int max_dim = 1 << 20;

	if (((int) machine >= EM_OS8 && machine <= EM_PC) &&
		(is_new_ver == 0 || is_new_ver == 1) &&
		(data_type >= EM_EM_CHAR && data_type <= EM_EM_DOUBLE) &&
		(nx > 1 && nx < max_dim) && (ny > 0 && ny < max_dim) && (nz > 0 && nz < max_dim)) {
		if (file_size > 0) {
			off_t file_size1 = (off_t)nx * (off_t)ny * (off_t)nz * (off_t)get_mode_size(data_type) + (off_t)sizeof(EMHeader);
			if (file_size == file_size1) {
				return true;
			}
		}
		else {
			return true;
		}
	}

	return false;
}

int EmIO::read_header(Dict & dict, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	//single image format, index can only be zero
	if(image_index == -1) {
		image_index = 0;
	}

	if(image_index != 0) {
		throw ImageReadException(filename, "no stack allowed for MRC image. For take 2D slice out of 3D image, read the 3D image first, then use get_clip().");
	}

	init();
	check_region(area, IntSize(emh.nx, emh.ny, emh.nz),false,false);

	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, emh.nx, &xlen, emh.ny, &ylen, emh.nz, &zlen);

	dict["nx"] = xlen;
	dict["ny"] = ylen;
	dict["nz"] = zlen;
	dict["datatype"] = to_em_datatype(emh.data_type);
	EXITFUNC;
	return 0;
}

int EmIO::write_header(const Dict & dict, int image_index, const Region* area,
					   EMUtil::EMDataType, bool)
{
	ENTERFUNC;
	//single image format, index can only be zero
	if(image_index != 0) {
		throw ImageWriteException(filename, "MRC file does not support stack.");
	}
	check_write_access(rw_mode, image_index, 1);
	if (area) {
		check_region(area, FloatSize(emh.nx, emh.ny, emh.nz), is_new_file);
		EXITFUNC;
		return 0;
	}

	emh.machine = static_cast < char >(get_machine_type());
	emh.nx = dict["nx"];
	emh.ny = dict["ny"];
	emh.nz = dict["nz"];
	emh.data_type = EM_EM_FLOAT;

	rewind(em_file);
	if (fwrite(&emh, sizeof(EMHeader), 1, em_file) != 1) {
		throw ImageWriteException(filename, "EM Header");
	}

	EXITFUNC;
	return 0;
}

int EmIO::read_data(float *data, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	//single image format, index can only be zero
	image_index = 0;
	check_read_access(image_index, data);
	check_region(area, IntSize(emh.nx, emh.ny, emh.nz),false,false);

	portable_fseek(em_file, sizeof(EMHeader), SEEK_SET);

	unsigned char *cdata = (unsigned char *) data;
	EMUtil::process_region_io(cdata, em_file, READ_ONLY, image_index, mode_size,
							  emh.nx, emh.ny, emh.nz, area);

	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, emh.nx, &xlen, emh.ny, &ylen, emh.nz, &zlen);

	int total_sz = xlen * ylen * zlen;

	if (mode_size == sizeof(short)) {
		become_host_endian((short *) cdata, total_sz);
	}
	else if (mode_size == sizeof(int)) {
		become_host_endian((int *) cdata, total_sz);
	}
	else if (mode_size == sizeof(double)) {
		throw ImageReadException(filename, "double type image is not supported");
	}

	for (int k = total_sz - 1; k >= 0; k--) {
		float curr_data = 0;

		if (mode == EM_EM_CHAR) {
			curr_data = static_cast < float >(cdata[k]);
		}
		else if (mode == EM_EM_SHORT) {
			curr_data = static_cast < float >(((short *) cdata)[k]);
		}
		else if (mode == EM_EM_INT) {
			curr_data = static_cast < float >(((int *) cdata)[k]);
		}
		else if (mode == EM_EM_FLOAT || mode == EM_EM_COMPLEX) {
			curr_data = ((float *) cdata)[k];
		}
		else if (mode_size == sizeof(double)) {
			throw ImageReadException(filename, "double type image is not supported");
		}

		data[k] = curr_data;
	}

	EXITFUNC;
	return 0;
}

int EmIO::write_data(float *data, int image_index, const Region* area, EMUtil::EMDataType, bool)
{
	ENTERFUNC;
	//single image format, index can only be zero
	image_index = 0;
	check_write_access(rw_mode, image_index, 1, data);
	portable_fseek(em_file, sizeof(EMHeader), SEEK_SET);

	EMUtil::process_region_io(data, em_file, rw_mode,
							  image_index, sizeof(float),
							  emh.nx, emh.ny, emh.nz, area);
#if 0
	int sec_size = emh.nx * emh.ny;
	int row_size = sizeof(float) * emh.nx;

	for (int i = 0; i < emh.nz; i++) {
		int k = i * sec_size;
		for (int j = 0; j < emh.ny; j++) {
			fwrite(&data[k + j * emh.nx], row_size, 1, em_file);
		}
	}
#endif

	EXITFUNC;
	return 0;
}

void EmIO::flush()
{
	fflush(em_file);
}

bool EmIO::is_complex_mode()
{
	init();
	if (emh.data_type == EM_EM_COMPLEX) {
		return true;
	}
	return false;
}

bool EmIO::is_image_big_endian()
{
	init();
	return is_big_endian;
}


int EmIO::to_em_datatype(char t)
{
	DataType type = static_cast < DataType > (t);
	switch (type) {
	case EM_EM_CHAR:
		return EMUtil::EM_CHAR;
	case EM_EM_SHORT:
		return EMUtil::EM_SHORT;
	case EM_EM_INT:
		return EMUtil::EM_INT;
	case EM_EM_FLOAT:
		return EMUtil::EM_FLOAT;
	case EM_EM_DOUBLE:
		return EMUtil::EM_DOUBLE;
	case EM_EM_COMPLEX:
		return EMUtil::EM_FLOAT_COMPLEX;
	default:
		break;
	}
	return EMUtil::EM_UNKNOWN;
}

int EmIO::get_machine_type()
{
	int m = EM_UNKNOWN_MACHINE;
#ifdef __sgi
	m = EM_SGI;
#elif defined __linux__
	m = EM_PC;
#elif defined __CYGWIN__
	m = EM_PC;
#elif defined WIN32
	m = EM_PC;
#elif defined MACOS
	m = EM_MAC;
#elif defined macintosh
	m = EM_MAC;
#elif defined __darwin__
	m = EM_MAC;
#elif defined __APPLE__
	m = EM_MAC;
#else
	m = EM_UNKNOWN_MACHINE;
#endif
	return m;
}

size_t EmIO::get_mode_size(char data_type)
{
	int mode = (int) data_type;
	switch (mode) {
	case EM_EM_CHAR:
		return sizeof(char);
	case EM_EM_SHORT:
		return sizeof(short);
	case EM_EM_INT:
	case EM_EM_FLOAT:
	case EM_EM_COMPLEX:
		return sizeof(int);
	case EM_EM_DOUBLE:
		return sizeof(double);
	}
	return 0;
}
