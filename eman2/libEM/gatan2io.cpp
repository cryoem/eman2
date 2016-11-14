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

#include "gatan2io.h"
#include "geometry.h"
#include "portable_fileio.h"
#include <cstring>

using namespace EMAN;

Gatan2IO::Gatan2IO(const string & file, IOMode rw)
:	filename(file), rw_mode(rw), gatan2_file(0), initialized(false)
{
	is_big_endian = ByteOrder::is_host_big_endian();
	memset(&gatanh, 0, sizeof(Gatan2Header));
}

Gatan2IO::~Gatan2IO()
{
	if (gatan2_file) {
		fclose(gatan2_file);
		gatan2_file = 0;
	}
}

void Gatan2IO::init()
{
	ENTERFUNC;
	
	if (initialized) {
		return;
	}
	
	initialized = true;

	bool is_new_file = false;
	gatan2_file = sfopen(filename, rw_mode, &is_new_file);

	if (!is_new_file) {
		if (fread(&gatanh, sizeof(Gatan2Header), 1, gatan2_file) != 1) {
			throw ImageReadException(filename, "Gatan2 Header");
		}

		if (!is_valid(&gatanh)) {
			throw ImageReadException(filename, "invalid Gatan2 file");
		}

		is_big_endian = ByteOrder::is_data_big_endian(&gatanh.len);
		become_host_endian((short *) &gatanh, sizeof(Gatan2Header) / sizeof(short));
	}
	EXITFUNC;
}

bool Gatan2IO::is_valid(const void *first_block)
{
	ENTERFUNC;
	bool result = false;
	
	if (!first_block) {
		result = false;
	}
	else {
		const short *data = static_cast < const short *>(first_block);
		short len = data[5];
		short type = data[6];

		bool data_big_endian = ByteOrder::is_data_big_endian(&len);

		if (data_big_endian != ByteOrder::is_host_big_endian()) {
			ByteOrder::swap_bytes(&len);
			ByteOrder::swap_bytes(&type);
		}

		int double_size = sizeof(double);
		if (len > 0 && len <= double_size && type > 0 && type <= GATAN2_INVALID) {
			result = true;
		}
	}
	EXITFUNC;
	return result;
}

int Gatan2IO::read_header(Dict & dict, int image_index, const Region * area, bool)
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
	
	if (is_complex_mode()) {
		throw ImageReadException(filename, "Cannot read complex Gatan2 files");
	}
	else {
		check_region(area, IntSize(gatanh.nx, gatanh.ny));

		int xlen = 0, ylen = 0;
		EMUtil::get_region_dims(area, gatanh.nx, &xlen, gatanh.ny, &ylen);
				
		dict["nx"] = xlen;
		dict["ny"] = ylen;
		dict["nz"] = 1;
		dict["datatype"] = to_em_datatype(gatanh.type);	
	}
	
	EXITFUNC;
	return 0;
}

int Gatan2IO::write_header(const Dict &, int, const Region* , EMUtil::EMDataType, bool)
{
	ENTERFUNC;
	LOGWARN("Gatan2 write is not supported.");
	EXITFUNC;
	return 1;
}

int Gatan2IO::read_data(float *data, int image_index, const Region * area, bool )
{
	ENTERFUNC;
	//single image format, index can only be zero
	image_index = 0;
	check_read_access(image_index, data);

	if (is_complex_mode()) {
		throw ImageReadException(filename, "Cannot read complex Gatan2 files");
	}
	
	check_region(area, IntSize(gatanh.nx, gatanh.ny));
	
	portable_fseek(gatan2_file, sizeof(Gatan2Header), SEEK_SET);

#if 0
	if (fread(data, gatanh.nx * gatanh.len, gatanh.ny, gatan2_file) != (unsigned int) gatanh.ny) {
		LOGDEBUG("Data read incomplete in Gatan file '%s'", filename.c_str());
		return 1;
	}
#endif

	int size = gatanh.nx * gatanh.ny;
	short *sdata = (short *) data;
	unsigned char *cdata = (unsigned char *) data;
	int *ldata = (int *) data;

	EMUtil::process_region_io(cdata, gatan2_file, READ_ONLY, image_index, gatanh.len,
							  gatanh.nx, gatanh.ny, 1, area);

	switch (gatanh.type) {
	case GATAN2_SHORT:
		become_host_endian((short *) data, size);
		for (int i = size - 1; i >= 0; i--) {
			data[i] = static_cast < float >(sdata[i]);
		}
		break;
	case GATAN2_FLOAT:
		become_host_endian(data, size);
		break;
	case GATAN2_CHAR:
		for (int i = size - 1; i >= 0; i--) {
			data[i] = static_cast < float >(cdata[i]);
		}
		break;
	case GATAN2_INT:
		become_host_endian((int *) data, size);
		for (int i = size - 1; i >= 0; i--) {
			data[i] = static_cast < float >(ldata[i]);
		}
		break;
	default:
		throw ImageReadException(filename, "unsupported Gatan2 data type");
	}
	EXITFUNC;
	return 0;
}

int Gatan2IO::write_data(float *, int, const Region*, EMUtil::EMDataType, bool)
{
	ENTERFUNC;
	LOGWARN("Gatan2 write is not supported.");
	EXITFUNC;
	return 1;
}

void Gatan2IO::flush()
{
}

bool Gatan2IO::is_complex_mode()
{
	init();
	if (gatanh.type == GATAN2_COMPLEX || gatanh.type == GATAN2_PACKED_COMPLEX) {
		return true;
	}
	return false;
}

bool Gatan2IO::is_image_big_endian()
{
	init();
	return is_big_endian;
}



int Gatan2IO::to_em_datatype(int gatan_type)
{
	switch (gatan_type) {
	case GATAN2_SHORT:
		return EMUtil::EM_SHORT;

	case GATAN2_FLOAT:
		return EMUtil::EM_FLOAT;

	case GATAN2_CHAR:
		return EMUtil::EM_CHAR;

	case GATAN2_INT:
		return EMUtil::EM_INT;
	}

	return EMUtil::EM_UNKNOWN;
}
