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

#include "pgmio.h"
#include "geometry.h"
#include "util.h"
#include "portable_fileio.h"

#ifdef _WIN32
	//MS Visual Studio.NET does not supply isnan()
	//they have _isnan() in <cfloat>
	#include <cfloat>
#endif

using namespace EMAN;

const char *PgmIO::MAGIC_BINARY = "P5";
const char *PgmIO::MAGIC_ASCII = "P2";

PgmIO::PgmIO(const string & file, IOMode rw)
:	filename(file), rw_mode(rw), pgm_file(0), is_big_endian(true),
	initialized(false), nx(0), ny(0), maxval(0), minval(0),
	file_offset(0), rendermin(0), rendermax(0)
{}

PgmIO::~PgmIO()
{
	if (pgm_file) {
		fclose(pgm_file);
		pgm_file = 0;
	}
}

//anonymous namespace make this function local for this file
namespace
{
	int read_int_and_space(FILE * in)
	{
		char buf[32];
		int c = 0;

		int i = 0;
		while (!isspace(c = getc(in))) {
			buf[i] = static_cast < char >(c);
			i++;
		}

		return atoi(buf);
	}
}

void PgmIO::init()
{
	ENTERFUNC;

	if (initialized) {
		return;
	}

	initialized = true;

	bool is_new_file = false;
	pgm_file = sfopen(filename, rw_mode, &is_new_file, true);

	if (!is_new_file) {
		const int bufsz = 1024;
		char buf[bufsz];

		buf[0] = static_cast < char >(getc(pgm_file));
		buf[1] = static_cast < char >(getc(pgm_file));
		buf[2] = '\0';
		getc(pgm_file);

		if (!is_valid(&buf)) {
			throw ImageReadException(filename, "invalid PGM file");
		}

		char c = '\0';

		while ((c = static_cast < char >(getc(pgm_file))) == '#') {
			fgets(buf, bufsz, pgm_file);
		}
		ungetc(c, pgm_file);

		nx = read_int_and_space(pgm_file);
		ny = read_int_and_space(pgm_file);
		maxval = read_int_and_space(pgm_file);

		if (nx <= 0 || ny <= 0) {
			throw ImageReadException(filename, "file size < 0");
		}

		file_offset = portable_ftell(pgm_file);
	}
	EXITFUNC;
}

bool PgmIO::is_valid(const void *first_block)
{
	ENTERFUNC;
	bool result = false;
	if (first_block) {
		result = Util::check_file_by_magic(first_block, MAGIC_BINARY);
	}
	EXITFUNC;
	return result;
}

int PgmIO::read_header(Dict & dict, int image_index, const Region * area, bool)
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

	check_read_access(image_index);
	check_region(area, IntSize(nx, ny));
	int xlen = 0, ylen = 0;
	EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen);

	dict["nx"] = xlen;
	dict["ny"] = ylen;
	dict["nz"] = 1;

	dict["PGM.max_gray"] = maxval;
	dict["PGM.min_gray"] = minval;

	EXITFUNC;
	return 0;
}

int PgmIO::write_header(const Dict & dict, int image_index, const Region*,
						EMUtil::EMDataType, bool)
{
	ENTERFUNC;
	int err = 0;

	//single image format, index can only be zero
	if(image_index == -1) {
		image_index = 0;
	}
	if(image_index != 0) {
		throw ImageWriteException(filename, "PGM file does not support stack.");
	}
	check_write_access(rw_mode, image_index);

	int nz = dict["nz"];
	if ((int)nz != 1) {
		LOGERR("Cannot write 3D image as PGM. Your image nz = %d", nz);
		err = 1;
		throw ImageWriteException("N/A", "Cannot write 3D image as PGM.");
	}
	else {
		nx = dict["nx"];
		ny = dict["ny"];

		if(dict.has_key("min_grey")) minval = dict["min_gray"];
		if(dict.has_key("max_grey")) maxval = dict["max_gray"];

		//if we didn't get any good values from attributes, assign to 255 by default
#ifdef _WIN32
		if (maxval<=minval || _isnan(minval) || _isnan(maxval)) {
#else
			if (maxval<=minval || std::isnan(minval) || std::isnan(maxval)) {
#endif	//_WIN32
			maxval = 255;
		}

		if(dict.has_key("render_min")) rendermin=(float)dict["render_min"];	// float value representing black in the output
		if(dict.has_key("render_max")) rendermax=(float)dict["render_max"];	// float value representign white in the output

		fprintf(pgm_file, "%s\n%d %d\n%d\n", MAGIC_BINARY, nx, ny, maxval);
	}

	EXITFUNC;
	return err;
}

int PgmIO::read_data(float *data, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	//single image format, index can only be zero
	image_index = 0;
	check_read_access(image_index, data);
	check_region(area, IntSize(nx, ny));

	portable_fseek(pgm_file, file_offset, SEEK_SET);

	unsigned char *cdata = (unsigned char *) (data);
	size_t mode_size = sizeof(unsigned char);

	EMUtil::process_region_io(cdata, pgm_file, READ_ONLY, image_index,
							  mode_size, nx, ny, 1, area, true);

	int xlen = 0, ylen = 0;
	EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen);

	for (int k = xlen * ylen - 1; k >= 0; k--) {
		data[k] = static_cast < float >(cdata[k]);
	}

	EXITFUNC;
	return 0;
}

int PgmIO::write_data(float *data, int image_index, const Region* area,
					  EMUtil::EMDataType, bool)
{
	ENTERFUNC;

	//single image format, index can only be zero
	image_index = 0;
/*	if(area && (area->size[0]!=nx || area->size[1]!=ny)) {
		throw ImageWriteException("N/A", "No region writing for PGM images");
	}
*/

	check_write_access(rw_mode, image_index, 1, data);
	check_region(area, IntSize(nx, ny));

	// If we didn't get any parameters in 'render_min' or 'render_max', we need to find some good ones
	EMUtil::getRenderMinMax(data, nx, ny, rendermin, rendermax);

	unsigned char *cdata=(unsigned char *)malloc(nx*ny);	//cdata is the normalized data

	int old_add = 0;
	int new_add = 0;
	for( int j=0; j<ny; ++j ) {
		for( int i=0; i<nx; ++i) {
			old_add = j*nx+i;
			new_add = (ny-1-j)*nx + i;
			if( data[old_add]<rendermin ) {
				cdata[new_add] = 0;
			}
			else if( data[old_add]>rendermax )
			{
				cdata[new_add] = 255;
			}
			else {
				cdata[new_add] = (unsigned char)((data[old_add]-rendermin)/(rendermax-rendermin)*256.0);
			}
		}
	}

	size_t mode_size = sizeof(unsigned char);
	//fwrite(cdata, nx, ny, pgm_file);
	EMUtil::process_region_io(cdata, pgm_file, WRITE_ONLY, image_index,
							  mode_size, nx, ny, 1, area);

	free(cdata);
	EXITFUNC;
	return 0;
}

void PgmIO::flush()
{
	fflush(pgm_file);
}


bool PgmIO::is_complex_mode()
{
	return false;
}

bool PgmIO::is_image_big_endian()
{
	init();
	return is_big_endian;
}

