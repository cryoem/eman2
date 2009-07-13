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

#include "pifio.h"
#include "portable_fileio.h"
#include "geometry.h"
#include "imageio.h"
#include <cstring>

#ifdef WIN32
#include <time.h>
#endif

using namespace EMAN;

PifIO::PifIO(const string & pif_filename, IOMode rw)
:	filename(pif_filename), rw_mode(rw)
{
	pif_file = 0;
	mode_size = 0;
	is_big_endian = ByteOrder::is_host_big_endian();
	initialized = false;
	real_scale_factor = 1;
	is_new_file = false;
	memset(&pfh, 0, sizeof(PifFileHeader));
}

PifIO::~PifIO()
{
	if (pif_file) {
		fclose(pif_file);
		pif_file = 0;
	}
}


int PifIO::get_mode_size(PifDataMode mode)
{
	int size = 0;

	switch (mode) {
	case PIF_CHAR:
	case PIF_BOXED_DATA:
		size = sizeof(char);
		break;
	case PIF_SHORT:
	case PIF_SHORT_FLOAT:
	case PIF_SHORT_COMPLEX:
	case PIF_SHORT_FLOAT_COMPLEX:
	case PIF_MAP_FLOAT_SHORT:
		size = sizeof(short);
		break;
	case PIF_FLOAT:
	case PIF_FLOAT_COMPLEX:
		size = sizeof(float);
		break;
	case PIF_FLOAT_INT:
	case PIF_FLOAT_INT_COMPLEX:
	case PIF_MAP_FLOAT_INT:
		size = sizeof(int);
		break;
	default:
		break;
	}
	return size;
}

bool PifIO::is_float_int(int m)
{
	PifDataMode mode = static_cast < PifDataMode > (m);
	switch (mode) {
	case PIF_SHORT_FLOAT:
	case PIF_SHORT_FLOAT_COMPLEX:
	//case PIF_FLOAT:
	case PIF_FLOAT_INT:
	//case PIF_FLOAT_COMPLEX:
	case PIF_FLOAT_INT_COMPLEX:
	case PIF_MAP_FLOAT_SHORT:
	case PIF_MAP_FLOAT_INT:
		return true;
	default:
		break;
	}
	return false;
}

void PifIO::init()
{
	ENTERFUNC;
	if (initialized) {
		return;
	}

	initialized = true;
	pif_file = sfopen(filename, rw_mode, &is_new_file);

	if (!is_new_file) {
		if (fread(&pfh, sizeof(PifFileHeader), 1, pif_file) != 1) {
			throw ImageReadException(filename, "PIF file header");
		}

		if (!is_valid(&pfh)) {
			throw ImageReadException(filename, "invalid PIF file");
		}

		is_big_endian = ByteOrder::is_data_big_endian(&pfh.nz);
		become_host_endian(&pfh.htype);

		if (pfh.htype != 1) {
			string desc = "only support PIF with all projects having the same dimensions";
			throw ImageReadException(filename, desc);
		}

		become_host_endian(&pfh.mode);
		become_host_endian(&pfh.nx);
		become_host_endian(&pfh.ny);
		become_host_endian(&pfh.nz);
		become_host_endian(&pfh.nimg);

		if (is_float_int(pfh.mode)) {
			real_scale_factor = (float) atof(pfh.scalefactor);
		}

		mode_size = get_mode_size(static_cast < PifDataMode > (pfh.mode));

		if (is_complex_mode()) {
			pfh.nx *= 2;
		}
	}
	EXITFUNC;
}

bool PifIO::is_valid(const void *first_block)
{
	ENTERFUNC;
	bool result = false;

	if (first_block) {
		const int *data = static_cast < const int *>(first_block);
		int m1 = data[0];
		int m2 = data[1];
		int endian = data[7];
		bool data_big_endian = false;
		if (endian) {
			data_big_endian = true;
		}

		if (data_big_endian != ByteOrder::is_host_big_endian()) {
			ByteOrder::swap_bytes(&m1);
			ByteOrder::swap_bytes(&m2);
		}

		if (m1 == PIF_MAGIC_NUM && m2 == PIF_MAGIC_NUM) {
			 result = true;
		}
	}

	EXITFUNC;
	return result;
}

void PifIO::fseek_to(int image_index)
{
	int pih_sz = sizeof(PifImageHeader);
	int image_size = 0;

#if 0
	// this works for some images that PURDUE people gave to me.
	// But those images don't follow the PIF specification. So
	// I believe they are in wrong format.
	if (pfh.nimg == 1) {
		image_size = pfh.nx * pfh.ny * pfh.nz;
	}
	else {
		image_size = pfh.nx * pfh.ny;
	}
#endif
	image_size = pfh.nx * pfh.ny * pfh.nz;

	size_t file_offset = sizeof(PifFileHeader) +
		(pih_sz + image_size * mode_size) * image_index;

	portable_fseek(pif_file, file_offset, SEEK_SET);
}


int PifIO::read_header(Dict & dict, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	check_read_access(image_index);
	fseek_to(image_index);

	int pih_sz = sizeof(PifImageHeader);
	PifImageHeader pih;

	if (fread(&pih, pih_sz, 1, pif_file) != 1) {
		throw ImageReadException(filename, "PIF Image header");
	}
	else {
		check_region(area, FloatSize(pih.nx, pih.ny, pih.nz), is_new_file);
		int xlen = 0, ylen = 0, zlen = 0;
		EMUtil::get_region_dims(area, pih.nx, &xlen, pih.ny, &ylen, pih.nz, &zlen);

		dict["nx"] = xlen;
		dict["ny"] = ylen;
		dict["nz"] = zlen;

		dict["datatype"] = to_em_datatype(pih.mode);

		dict["apix_x"] = static_cast < float >(pih.xlen);
		dict["apix_y"] = static_cast < float >(pih.ylen);
		dict["apix_z"] = static_cast < float >(pih.zlen);

		dict["minimum"] = static_cast < float >(pih.min);
		dict["maximum"] = static_cast < float >(pih.max);
		dict["mean"] = static_cast < float >(pih.mean);
		dict["sigma"] = static_cast < float >(pih.sigma);

		dict["origin_x"] = static_cast < float >(pih.xorigin);
		dict["origin_y"] = static_cast < float >(pih.yorigin);
	}

	EXITFUNC;

	return 0;
}

int PifIO::write_header(const Dict & dict, int image_index, const Region* area,
						EMUtil::EMDataType, bool)
{
	ENTERFUNC;

	check_write_access(rw_mode, image_index);

	if (area) {
		check_region(area, FloatSize(pfh.nx, pfh.ny, pfh.nz), is_new_file);
		EXITFUNC;
		return 0;
	}

	time_t t0 = time(0);
	struct tm *t = localtime(&t0);

	if (!is_new_file) {
		if (is_big_endian != ByteOrder::is_host_big_endian()) {
			throw ImageWriteException(filename, "writing to opposite byteorder file");
		}

		int nx1 = dict["nx"];
		int ny1 = dict["ny"];
		int nz1 = dict["nz"];

		int mode1 = to_pif_datatype(dict["datatype"]);

		if (nx1 != pfh.nx || ny1 != pfh.ny || nz1 != pfh.nz) {
			LOGERR("cannot write to different size file %s (%dx%dx%d != %dx%dx%d)",
				   filename.c_str(), nx1, ny1, nz1, pfh.nx, pfh.ny, pfh.nz);
			throw ImageWriteException(filename, "different image size");
		}

		if (mode1 != pfh.mode) {
			throw ImageWriteException(filename, "different data type");
		}

		fseek_to(image_index);
	}
	else {
		pfh.magic[0] = PIF_MAGIC_NUM;
		pfh.magic[1] = PIF_MAGIC_NUM;
		sprintf(pfh.scalefactor, "1.0");

		pfh.mode = PIF_FLOAT_INT;
		sprintf(pfh.program, "EMAN %d/%02d/%02d", t->tm_mon, t->tm_mday, t->tm_year);

		pfh.htype = 1;
		pfh.nx = dict["nx"];
		pfh.ny = dict["ny"];
		pfh.nz = dict["nz"];
		pfh.nimg += 1;
		pfh.endian = (int) ByteOrder::is_host_big_endian();
		fwrite(&pfh, sizeof(PifFileHeader), 1, pif_file);
	}

	PifImageHeader pih;
	memset(&pih, 0, sizeof(PifImageHeader));
	pih.nx = dict["nx"];
	pih.ny = dict["ny"];
	pih.nz = dict["nz"];

	pih.mode = PIF_FLOAT;
	pih.xlen = dict["apix_x"];
	pih.ylen = dict["apix_y"];
	pih.zlen = dict["apix_z"];
	pih.alpha = 90;
	pih.beta = 90;
	pih.gamma = 90;
	pih.mapc = 1;
	pih.mapr = 2;
	pih.maps = 3;
	pih.min = dict["minimum"];
	pih.max = dict["maximum"];
	pih.mean = dict["mean"];
	pih.sigma = dict["sigma"];

	pih.xorigin = dict["origin_x"];
	pih.yorigin = dict["origin_y"];

	sprintf(pih.time, "%d/%02d/%02d %d:%02d",
			t->tm_mon, t->tm_mday, t->tm_year, t->tm_hour, t->tm_min);
	fwrite(&pih, sizeof(PifImageHeader), 1, pif_file);

	EXITFUNC;
	return 0;
}

int PifIO::read_data(float *data, int image_index, const Region *area, bool)
{
	ENTERFUNC;

	check_read_access(image_index, data);
	fseek_to(image_index);

	int pih_sz = sizeof(PifImageHeader);
	PifImageHeader pih;

	if (fread(&pih, pih_sz, 1, pif_file) != 1) {
		throw ImageReadException(filename, "PIF Image header");
	}

	if (area) {
		check_region(area, FloatSize(pih.nx, pih.ny, pih.nz), is_new_file);
	}

	PifDataMode data_mode = static_cast < PifDataMode > (pih.mode);
	int num_layers = pih.nz;
#if 0
	if (pfh.nz == pfh.nimg) {
		num_layers = 1;
	}
#endif
	// new way to read PIF data. The new way includes region reading.
	// If it is tested to be OK, remove the code in #if 0 ... #endif
	unsigned char * cdata = (unsigned char*)data;
	short *sdata = (short*) data;

	EMUtil::process_region_io(cdata, pif_file, READ_ONLY,
							  0, mode_size, pih.nx, pih.ny,
							  num_layers, area);

	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, pih.nx, &xlen, pih.ny, &ylen, pih.nz, &zlen);
	size_t size = xlen * ylen * zlen;

	if(data_mode == PIF_FLOAT || data_mode == PIF_FLOAT_COMPLEX)
	{
		become_host_endian< float >(data, size);
	}
	else {
		if (mode_size == sizeof(short)) {
			become_host_endian((short *) sdata, size);
		}
		else if (mode_size == sizeof(int)) {
			become_host_endian((int *) data, size);
		}

		if (mode_size == sizeof(char)) {
			for (size_t i = 0; i < size; i++) {
				size_t j = size - 1 - i;
				data[j] = (float)(cdata[j]) * real_scale_factor;
			}
		}
		else if (mode_size == sizeof(short)) {
			for (size_t i = 0; i < size; i++) {
				size_t j = size - 1 - i;
				data[j] = (float)(sdata[j]) * real_scale_factor;
			}
		}
		else if (mode_size == sizeof(int)) {
			for (size_t i = 0; i < size; i++) {
				size_t j = size - 1 - i;
				data[j] = (float) ((int *)data)[j] * real_scale_factor;
			}
		}
	}

	// end of new way for region reading

#if 0
	int buf_size = pih.nx * mode_size;
	unsigned char *buf = new unsigned char[buf_size];

	for (int l = 0; l < num_layers; l++) {
		int offset1 = l * pfh.nx * pfh.ny;
		for (int j = 0; j < pfh.ny; j++) {
			if (fread(buf, mode_size, pfh.nx, pif_file) != (unsigned int) pfh.nx) {
				if( buf )
				{
					delete[]buf;
					buf = 0;
				}
				throw ImageReadException(filename, "incomplete PIF read");
			}

			if (mode_size == sizeof(short)) {
				become_host_endian((short *) buf, pfh.nx);
			}
			else if (mode_size == sizeof(int)) {
				become_host_endian((int *) buf, pfh.nx);
			}

			int offset2 = offset1 + j * pfh.nx;

			for (int k = 0; k < pfh.nx; k++) {
				float curr_data = 0;

				if (mode_size == sizeof(char)) {
					curr_data = (float) buf[k];
				}
				else if (mode_size == sizeof(short)) {
					curr_data = (float) ((short *) buf)[k];
				}
				else if (mode_size == sizeof(int)) {
					curr_data = (float) ((int *) buf)[k];
				}
				data[offset2 + k] = curr_data * real_scale_factor;
			}
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


int PifIO::write_data(float *data, int image_index, const Region* area,
					  EMUtil::EMDataType, bool)
{
	ENTERFUNC;

	check_write_access(rw_mode, image_index, 0, data);
	fseek_to(image_index);

	int nx = pfh.nx;
	int ny = pfh.ny;
	int nz = pfh.nz;

	check_region(area, FloatSize(nx, ny, nz), is_new_file);

	// to write a region; if it works, please remove the code
	// in #if 0 ... #endif
	EMUtil::process_region_io(data, pif_file, WRITE_ONLY, 0,
							  mode_size, nx, ny, nz, area);

#if 0
	size_t idx;
	int *buf = new int[nx];
	for (int i = 0; i < nz; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < pfh.nx; k++) {
				idx = i * nx * ny + j * nx + k;
				buf[k] = (int) data[idx];
			}
			fwrite(buf, sizeof(int) * nx, 1, pif_file);
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

void PifIO::flush()
{
	fflush(pif_file);
}

bool PifIO::is_complex_mode()
{
	init();
	if (pfh.mode == PIF_SHORT_COMPLEX ||
		pfh.mode == PIF_FLOAT_INT_COMPLEX ||
		pfh.mode == PIF_FLOAT_COMPLEX || pfh.mode == PIF_SHORT_FLOAT_COMPLEX) {
		return true;
	}
	return false;
}

bool PifIO::is_image_big_endian()
{
	init();
	return is_big_endian;
}

int PifIO::get_nimg()
{
	init();
	return pfh.nimg;
}


int PifIO::to_em_datatype(int p)
{
	PifDataMode mode = static_cast < PifDataMode > (p);
	EMUtil::EMDataType e = EMUtil::EM_UNKNOWN;

	switch (mode) {
	case PIF_CHAR:
	case PIF_BOXED_DATA:
		e = EMUtil::EM_CHAR;
		break;

	case PIF_SHORT:
	case PIF_SHORT_FLOAT:
	case PIF_MAP_FLOAT_SHORT:
		e = EMUtil::EM_SHORT;
		break;

	case PIF_SHORT_COMPLEX:
	case PIF_SHORT_FLOAT_COMPLEX:
		e = EMUtil::EM_SHORT_COMPLEX;
		break;

	case PIF_FLOAT:
	case PIF_FLOAT_INT:
	case PIF_MAP_FLOAT_INT:
		e = EMUtil::EM_FLOAT;
		break;
	case PIF_FLOAT_COMPLEX:
	case PIF_FLOAT_INT_COMPLEX:
		e = EMUtil::EM_FLOAT_COMPLEX;
		break;
	case PIF_INVALID:
		e = EMUtil::EM_UNKNOWN;
		break;
	}
	return e;
}

int PifIO::to_pif_datatype(int e)
{
	PifDataMode m = PIF_INVALID;

	switch (e) {
	case EMUtil::EM_CHAR:
		m = PIF_BOXED_DATA;
		break;
	case EMUtil::EM_SHORT:
		m = PIF_SHORT;
		break;
	case EMUtil::EM_SHORT_COMPLEX:
		m = PIF_SHORT_COMPLEX;
		break;
	case EMUtil::EM_FLOAT:
		m = PIF_FLOAT_INT;
		break;
	case EMUtil::EM_FLOAT_COMPLEX:
		m = PIF_FLOAT_COMPLEX;
		break;
	default:
		LOGERR("unknown PIF mode: %d", e);
	}

	return m;
}
