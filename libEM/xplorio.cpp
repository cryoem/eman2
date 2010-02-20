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
#include "xplorio.h"
#include "util.h"
#include "emassert.h"
#include "portable_fileio.h"
#include "geometry.h"

#ifdef WIN32
#include <time.h>
#include <windows.h>
#define MAXPATHLEN (MAX_PATH*4)
#else
#include <sys/param.h>
#endif


using namespace EMAN;

const string XplorIO::SECTION_MODE = "ZYX";
const int XplorIO::NFLOAT_PER_LINE = 6;
const int XplorIO::INTEGER_SIZE = 8;
const int XplorIO::FLOAT_SIZE = 12;
const char * XplorIO::OUTFORMAT = "%12.5E";

XplorIO::XplorIO(const string & file, IOMode rw)
:	filename(file), rw_mode(rw), xplor_file(0), initialized(false)
{
	is_big_endian = ByteOrder::is_host_big_endian();
	is_new_file = false;
	nlines_in_header = 0;

	nx = 0;
	ny = 0;
	nz = 0;

	apix_x = 0;
	apix_y = 0;
	apix_z = 0;

	cell_alpha = 0;
	cell_beta = 0;
	cell_gama = 0;
}

XplorIO::~XplorIO()
{
	if (xplor_file) {
		fclose(xplor_file);
		xplor_file = 0;
	}
}

void XplorIO::init()
{
	if (initialized) {
		return;
	}

	ENTERFUNC;
	initialized = true;
	xplor_file = sfopen(filename, rw_mode, &is_new_file);

	if (!is_new_file) {
		char first_block[1024];
		fread(&first_block, sizeof(char), sizeof(first_block), xplor_file);
		if (!is_valid(&first_block)) {
			throw ImageReadException(filename, "invalid XPLOR");
		}
		portable_fseek(xplor_file, 0, SEEK_SET);
		char line[1024];
		int i = 1;
		int ntitle = 0;

		int xmin = 0;
		int xmax = 0;
		int ymin = 0;
		int ymax = 0;
		int zmin = 0;
		int zmax = 0;

		float cellx = 0;
		float celly = 0;
		float cellz = 0;

		while(fgets(line, sizeof(line), xplor_file)) {
			line[strlen(line)-1] = '\0';
			if (i == 2) {
				ntitle = atoi(line);
			}
			else if (i == (ntitle+3)) {
				if (sscanf(line, "%8d%8d%8d%8d%8d%8d%8d%8d%8d", &nx, &xmin, &xmax,
						   &ny, &ymin, &ymax, &nz, &zmin, &zmax) != 9) {
					throw ImageReadException(filename, "invalid XPLOR");
				}
			}
			else if (i == (ntitle+4)) {
				if(sscanf(line, "%f %f %f %f %f %f",
						  &cellx, &celly, &cellz, &cell_alpha, &cell_beta, &cell_gama) != 6) {
					throw ImageReadException(filename, "invalid XPLOR");
				}
			}
			else if (i == (ntitle+5)) {
				break;
			}

			i++;
		}
		nlines_in_header = i;
		apix_x = cellx / nx;
		apix_y = celly / ny;
		apix_z = cellz / nz;
	}

	EXITFUNC;
}

bool XplorIO::is_valid(const void *first_block)
{
	ENTERFUNC;
	if (!first_block) {
		return false;
	}
	char *buf = (char *)(first_block);
	string line1 = Util::get_line_from_string(&buf);
	bool result = true;

	if (line1.size() != 0) {
		result = false;
	}
	else {
		string line2 = Util::get_line_from_string(&buf);
		int ntitle = 0;

		if ((int)line2.size() != INTEGER_SIZE) {
			result = false;
		}
		else {
			ntitle = atoi(line2.c_str());
			if (ntitle < 0 || ntitle > 50) {
				result = false;
			}

			else {
				for (int i = 0; i < ntitle+2; i++) {
					Util::get_line_from_string(&buf);
				}

				string modeline = Util::get_line_from_string(&buf);
				if (modeline != SECTION_MODE) {
					result = false;
				}
			}
		}
	}

	EXITFUNC;
	return result;
}

int XplorIO::read_header(Dict &dict, int image_index, const Region *area, bool)
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
	check_region(area, FloatSize(nx, ny, nz), is_new_file);

	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen, nz, &zlen);

	dict["nx"] = xlen;
	dict["ny"] = ylen;
	dict["nz"] = zlen;

	dict["apix_x"] = apix_x;
	dict["apix_y"] = apix_y;
	dict["apix_z"] = apix_z;

	dict["XPLOR.alpha"] = cell_alpha;
	dict["XPLOR.beta"] = cell_beta;
	dict["XPLOR.gama"] = cell_gama;

	EXITFUNC;
	return 0;
}

int XplorIO::write_header(const Dict & dict, int image_index, const Region* area,
						  EMUtil::EMDataType, bool)
{
	ENTERFUNC;
	//single image format, index can only be zero
	if(image_index == -1) {
		image_index = 0;
	}
	if(image_index != 0) {
		throw ImageWriteException(filename, "XPLOR file does not support stack.");
	}
	check_write_access(rw_mode, image_index);
	if (area) {
		check_region(area, FloatSize(nx, ny, nz), is_new_file);
		EXITFUNC;
		return 0;
	}

	nx = dict["nx"];
	ny = dict["ny"];
	nz = dict["nz"];
	float pixel = dict["pixel"];

	nlines_in_header = 0;
	time_t t0 = time(0);
	struct tm *t = localtime(&t0);
	rewind(xplor_file);

	fprintf(xplor_file, "\n");
	fprintf(xplor_file, "%8d\n", 1);
	fprintf(xplor_file, "\"%s\" written by EMAN at %s", filename.c_str(), asctime(t));

	int z0 = -nz / 2;
	int z1 = (nz - 1) / 2;

	if (2 * nz - 1 == nx && 2 * nz - 1 == ny) {
		z0 = 0;
		z1 = nz - 1;
	}

	fprintf(xplor_file, "%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
			nx, -nx / 2, nx % 2 ? nx / 2 : nx / 2 - 1, ny, -ny / 2,
			ny % 2 ? ny / 2 : ny / 2 - 1, nz, z0, z1);

	char fformat[256];
	sprintf(fformat, "%s%s%s%s%s%s\n",
			OUTFORMAT, OUTFORMAT,OUTFORMAT, OUTFORMAT, OUTFORMAT,OUTFORMAT);

	fprintf(xplor_file, fformat,
			nx * pixel, ny * pixel, nz * pixel, 90.0, 90.0, 90.0);
	fprintf(xplor_file, "ZYX\n");
	nlines_in_header = 5;
	flush();

	EXITFUNC;
	return 0;
}

int XplorIO::read_data(float *data, int image_index, const Region *area, bool)
{
	ENTERFUNC;

	//single image format, index can only be zero
	image_index = 0;
	check_read_access(image_index);
	FloatSize max_size = FloatSize(nx, ny, nz);
	check_region(area, max_size, is_new_file);

	// Had to put this here because this is the only function that calls
	// EMUtil::process_ascii_region_io - this function does not allow regions that
	// are outside the image dimensions. This is opposed to those functions which
	// call EMUtil::process_region_io, which can handle it. David Woolford, April 23 2009

	if (area != 0 && !area->is_region_in_box(max_size)) {
		char desc[1024];
		sprintf(desc, "Region box %s is outside image area (%d,%d,%d)",
				area->get_string().c_str(), (int)max_size[0],
				(int)max_size[1], (int)max_size[2]);
		throw ImageReadException("", desc);
	}

	Assert(nlines_in_header > 0);
	rewind(xplor_file);
	EMUtil::jump_lines(xplor_file, nlines_in_header);

	EMUtil::process_ascii_region_io(data, xplor_file, ImageIO::READ_ONLY, image_index,
									FLOAT_SIZE, nx, ny, nz, area, true,
									NFLOAT_PER_LINE, OUTFORMAT);

	EXITFUNC;
	return 0;
}


#if 0
int XplorIO::read_data(float *data, int, const Region *, bool)
{
	ENTERFUNC;
	int step = NFLOAT_PER_LINE;
	char line[1024];
	int nxy = nx * ny;
	int nlines = nxy / step;
	int nleft = nxy - nlines * step;

	for (int k = 0; k < nz; k++) {
		fgets(line, sizeof(line), xplor_file);
		int kk = 0;
		sscanf(line, "%d", &kk);
		if (kk != (k+1)) {
			LOGERR("section index = %d. It should be %d\n", kk, (k+1));
		}

		int k2 = k * nxy;

		for (int i = 0; i < nlines; i++) {
			fgets(line, sizeof(line), xplor_file);
			int i2 = k2 + i * step;
			sscanf(line, "%f %f %f %f %f %f",
				   &data[i2], &data[i2+1], &data[i2+2],
				   &data[i2+3], &data[i2+4], &data[i2+5]);
		}

		if (nleft > 0) {
			int i2 = k2 + nlines * step;
			fgets(line, sizeof(line), xplor_file);
			char *pline = line;
			for (int j = 0; j < nleft; j++) {
				sscanf(pline, "%f", &data[i2+j]);
				pline += FLOAT_SIZE;
			}
		}
	}


	EXITFUNC;
	return 0;
}
#endif


int XplorIO::write_data(float *data, int image_index, const Region* area,
						EMUtil::EMDataType, bool)
{

	ENTERFUNC;
	//single image format, index can only be zero
	image_index = 0;
	check_write_access(rw_mode, image_index, 1, data);
	check_region(area, FloatSize(nx,ny,nz), is_new_file);

	if (!is_new_file) {
		rewind(xplor_file);
		EMUtil::jump_lines(xplor_file, nlines_in_header);
	}

	int nsecs = nx * ny;
	int step = NFLOAT_PER_LINE;

	if (!area) {
		for (int k = 0; k < nz; ++k) {
			fprintf(xplor_file, "%8d\n", (k+1));

			for (int i = 0; i < nsecs - step; i += step) {
				for (int j = 0; j < step; j++) {
					fprintf(xplor_file, OUTFORMAT, data[k * nsecs + i + j]);
				}
				fprintf(xplor_file, "\n");
			}

			for (int l = (nsecs - 1) / step * step; l < nsecs; l++) {
				fprintf(xplor_file, OUTFORMAT, data[k * nsecs + l]);
			}

			fprintf(xplor_file, "\n");
		}

		// not sure what this is doing. so comment out.
		//fprintf(xplor_file, "%8d\n", -9999);
	}
	else {
		EMUtil::process_ascii_region_io(data, xplor_file, ImageIO::WRITE_ONLY,
										image_index, FLOAT_SIZE, nx, ny, nz,
										area, true, NFLOAT_PER_LINE, OUTFORMAT);

	}

	EXITFUNC;

	return 0;
}

void XplorIO::flush()
{
	fflush(xplor_file);
}

bool XplorIO::is_complex_mode()
{
	return false;
}

bool XplorIO::is_image_big_endian()
{
	init();
	return is_big_endian;
}
