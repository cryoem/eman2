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

#ifdef EM_PNG

#include <climits>
#include "pngio.h"
#include "geometry.h"
#include "util.h"


using namespace EMAN;

PngIO::PngIO(const string & file, IOMode rw)
:	filename(file), rw_mode(rw), png_file(0), initialized(false),
	png_ptr(0), info_ptr(0), end_info(0), nx(0), ny(0),
	depth_type(PNG_INVALID_DEPTH), number_passes(0), rendermin(0), rendermax(0)
{}

PngIO::~PngIO()
{
	if (png_file) {
		fclose(png_file);
		png_file = 0;
	}

	png_ptr = 0;
	info_ptr = 0;
	end_info = 0;
}

void PngIO::init()
{
	ENTERFUNC;
	if (initialized) {
		return;
	}

	initialized = true;

	bool is_new_file = false;
	png_file = sfopen(filename, rw_mode, &is_new_file, true);

	if (!is_new_file) {
		png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);
	}
	else {
		png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);
	}

	if (!png_ptr) {
		throw ImageReadException(filename, "cannot initialize libpng data structure");
	}

	info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr) {
		throw ImageReadException(filename, "cannot create png info data structure");
	}

	end_info = png_create_info_struct(png_ptr);
	if (!end_info) {
		throw ImageReadException(filename, "cannot create png end info structure");
	}

	if (setjmp(png_ptr->jmpbuf)) {
		throw ImageReadException(filename, "an error occurs within png");
	}

	png_init_io(png_ptr, png_file);

	if (!is_new_file) {
		unsigned char header[PNG_BYTES_TO_CHECK];
		fread(header, sizeof(unsigned char), PNG_BYTES_TO_CHECK, png_file);
		if (!is_valid(header)) {
			throw ImageReadException(filename, "invalid PNG format");
		}

		png_set_sig_bytes(png_ptr, PNG_BYTES_TO_CHECK);

		png_read_info(png_ptr, info_ptr);

		nx = png_get_image_width(png_ptr, info_ptr);
		ny = png_get_image_height(png_ptr, info_ptr);
		int bit_depth = png_get_bit_depth(png_ptr, info_ptr);
		int color_type = png_get_color_type(png_ptr, info_ptr);

		if (nx == 0 || ny == 0) {
			throw ImageReadException(filename, "PNG file size = 0");
		}

		if (bit_depth == CHAR_BIT) {
			depth_type = PNG_CHAR_DEPTH;
		}
		else if (bit_depth == CHAR_BIT * sizeof(short)) {
			depth_type = PNG_SHORT_DEPTH;
		}
		else {
			depth_type = PNG_INVALID_DEPTH;
			char desc[256];
			sprintf(desc, "not support png with depth = %d bit", bit_depth);
			throw ImageReadException(filename, desc);
		}

		png_set_packing(png_ptr);

		if ((color_type == PNG_COLOR_TYPE_GRAY) && (bit_depth < CHAR_BIT)) {
			png_set_expand(png_ptr);
		}

		number_passes = png_set_interlace_handling(png_ptr);

		if (bit_depth > CHAR_BIT) {
			png_set_swap(png_ptr);
		}

		png_read_update_info(png_ptr, info_ptr);
	}
	EXITFUNC;
}

bool PngIO::is_valid(const void *first_block)
{
	ENTERFUNC;
	bool result = false;

	if (!first_block) {
		result = false;
	}
	else {
		if (png_sig_cmp((png_byte *) first_block, (png_size_t) 0, PNG_BYTES_TO_CHECK) == 0) {
			result = true;
		}
	}
	EXITFUNC;
	return result;
}

int PngIO::read_header(Dict & dict, int image_index, const Region * area, bool)
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
	
	int nx1 = static_cast < int >(nx);
	int ny1 = static_cast < int >(ny);
	check_region(area, IntSize(nx1, ny1));
	int xlen = 0, ylen = 0;
	EMUtil::get_region_dims(area, nx1, &xlen, ny1, &ylen);

	dict["nx"] = xlen;
	dict["ny"] = ylen;
	dict["nz"] = 1;

	if (depth_type == PNG_CHAR_DEPTH) {
		dict["datatype"] = EMUtil::EM_UCHAR;
	}
	else if (depth_type == PNG_SHORT_DEPTH) {
		dict["datatype"] = EMUtil::EM_USHORT;
	}
	else {
		throw ImageReadException(filename, "unsupported PNG bit depth");
	}

	EXITFUNC;
	return 0;
}

int PngIO::write_header(const Dict & dict, int image_index, const Region*,
						EMUtil::EMDataType, bool)
{
	ENTERFUNC;

	//single image format, index can only be zero
	if(image_index != 0) {
		throw ImageWriteException(filename, "MRC file does not support stack.");
	}
	check_write_access(rw_mode, image_index);

	nx = (png_uint_32) (int) dict["nx"];
	ny = (png_uint_32) (int) dict["ny"];
	int nz = dict["nz"];
	if (nz != 1) {
		LOGERR("Only support 2D PNG file write");
		return 1;
	}

	int bit_depth = 0;
	EMUtil::EMDataType datatype = (EMUtil::EMDataType) (int) dict["datatype"];

	if (datatype == EMUtil::EM_UCHAR) {
		depth_type = PNG_CHAR_DEPTH;
		bit_depth = CHAR_BIT;
	}
	else {
		if (datatype != EMUtil::EM_USHORT) {
			LOGWARN("Don't support data type '%s' in PNG. Convert to '%s'.",
					EMUtil::get_datatype_string(datatype),
					EMUtil::get_datatype_string(EMUtil::EM_USHORT));
		}
		depth_type = PNG_SHORT_DEPTH;
		bit_depth = sizeof(unsigned short) * CHAR_BIT;
	}

	png_set_IHDR(png_ptr, info_ptr, nx, ny, bit_depth, PNG_COLOR_TYPE_GRAY,
				 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

	png_write_info(png_ptr, info_ptr);

	if (depth_type == PNG_SHORT_DEPTH) {
		png_set_swap(png_ptr);
	}

	if(dict.has_key("render_min")) rendermin=(float)dict["render_min"];
	else rendermin=0;
	if(dict.has_key("render_max")) rendermax=(float)dict["render_max"];
	else rendermax=0;
	EXITFUNC;
	return 0;
}

int PngIO::read_data(float *data, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	//single image format, index can only be zero
	image_index = 0;
	check_read_access(image_index, data);

	int nx1 = static_cast < int >(nx);
	int ny1 = static_cast < int >(ny);

	check_region(area, IntSize(nx1, ny1));

	png_init_io(png_ptr, png_file);
	png_set_sig_bytes(png_ptr, PNG_BYTES_TO_CHECK);

	int xlen = 0, ylen = 0, x0 = 0, y0 = 0;
	EMUtil::get_region_dims(area, nx1, &xlen, ny1, &ylen);
	EMUtil::get_region_origins(area, &x0, &y0);

	png_uint_32 rowbytes = png_get_rowbytes(png_ptr, info_ptr);
	unsigned char *cdata = new unsigned char[rowbytes];
	unsigned short *sdata = (unsigned short *) cdata;

	int k = 0;
	for (int i = y0; i < y0 + ylen; i++) {
		for (int pass = 0; pass < number_passes; pass++) {
			png_read_rows(png_ptr, (png_byte **) & cdata, 0, 1);
		}

		if (depth_type == PNG_CHAR_DEPTH) {
			for (int x = x0; x < x0 + xlen; x++) {
				data[k] = static_cast < float >(cdata[x]);
				k++;
			}
		}
		else if (depth_type == PNG_SHORT_DEPTH) {
			for (int x = x0; x < x0 + xlen; x++) {
				data[k] = static_cast < float >(sdata[x]);
				k++;
			}
		}
	}

	//Util::flip_image(data, nx, ny);

	if( cdata )
	{
		delete[]cdata;
		cdata = 0;
	}

	png_read_end(png_ptr, end_info);
	png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
	EXITFUNC;
	return 0;
}

int PngIO::write_data(float *data, int image_index, const Region*,
					  EMUtil::EMDataType, bool)
{
	ENTERFUNC;

	//single image format, index can only be zero
	image_index = 0;
	check_write_access(rw_mode, image_index, 1, data);

	// If we didn't get any parameters in 'render_min' or 'render_max', we need to find some good ones
	if (!rendermin && !rendermax) getRenderMinMax(data, nx, ny, rendermin, rendermax);
	
	if (depth_type == PNG_CHAR_DEPTH) {
		unsigned char *cdata = new unsigned char[nx];

		for (unsigned int y = 0; y < ny; y++) {
			for (unsigned int x = 0; x < nx; x++) {
				if(data[y * nx + x] <= rendermin){
					cdata[x] = 0;
				}
				else if(data[y * nx + x] >= rendermax) {
					cdata[x] = UCHAR_MAX;
				}
				else {
					cdata[x] = (unsigned char)((data[y * nx + x] - rendermin) / (rendermax - rendermin) * 256);
				}
			}
			png_write_row(png_ptr, (png_byte *) cdata);
		}

		if( cdata )
		{
			delete[]cdata;
			cdata = 0;
		}
	}
	else if (depth_type == PNG_SHORT_DEPTH) {
		unsigned short *sdata = new unsigned short[nx];

		for (unsigned int y = 0; y < ny; y++) {
			for (unsigned int x = 0; x < nx; x++) {
				if(data[y * nx + x] <= rendermin){
					sdata[x] = 0;
				}
				else if(data[y * nx + x] >= rendermax) {
					sdata[x] = USHRT_MAX;
				}
				else {
					sdata[x] = (unsigned short)((data[y * nx + x] - rendermin) / (rendermax - rendermin) * 65536);
				}
			}

			png_write_row(png_ptr, (png_byte *) sdata);
		}

		if( sdata )
		{
			delete[]sdata;
			sdata = 0;
		}
	}

	png_write_end(png_ptr, info_ptr);
	png_destroy_write_struct(&png_ptr, &info_ptr);

	EXITFUNC;
	return 0;
}

void PngIO::flush()
{
	png_write_flush(png_ptr);
}

bool PngIO::is_complex_mode()
{
	return false;
}

bool PngIO::is_image_big_endian()
{
	return true;
}


#endif	//EM_PNG
