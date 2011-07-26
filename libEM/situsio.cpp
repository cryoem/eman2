/**
 * $Id$
 */

/*
 * Author: Grant Tang, 06/07/2011 (gtang@bcm.edu)
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

#include "situsio.h"
#include "portable_fileio.h"
#include "util.h"

using namespace EMAN;

const int SitusIO::SITUS_HEADER_LINES=2;
const int SitusIO::FLOAT_SIZE = 12;
const int SitusIO::NFLOAT_PER_LINE = 10;
const char * SitusIO::OUTFORMAT = "%12.6f";
const int SitusIO::LINE_LENGTH = 1024;

SitusIO::SitusIO(const string & situsname, IOMode rw) :
		filename(situsname), rw_mode(rw), situsfile(0),
		initialized(false), is_new_file(false),
		apix(0.0f), origx(0.0f), origy(0.0f), origz(0.0f),
		nx(0), ny(0), nz(0)
{
}

SitusIO::~SitusIO()
{
	if (situsfile) {
		fclose(situsfile);
		situsfile = 0;
	}
}

void SitusIO::init()
{
	ENTERFUNC;
	if (initialized) {
		return;
	}

	initialized = true;
	situsfile = sfopen(filename, rw_mode, &is_new_file);

	if (!is_new_file) {
		char first_block[1024];
		fread(&first_block, sizeof(char), sizeof(first_block), situsfile);
		if (!is_valid(&first_block)) {
			throw ImageReadException(filename, "invalid SITUS file");
		}

		char * buf = (char *)first_block;
		string line1 = Util::get_line_from_string(&buf);

		sscanf(line1.c_str(), "%f %f %f %f %d %d %d", &apix, &origx, &origy, &origz, &nx, &ny, &nz);
	}

	EXITFUNC;
}

int SitusIO::read_data(float* data, int image_index, const EMAN::Region*, bool)
{
	ENTERFUNC;

	image_index = 0;	//single image format

	portable_fseek(situsfile, 0, SEEK_SET);
	EMUtil::jump_lines(situsfile, SITUS_HEADER_LINES);	//skip header lines

	int number_lines = nx*ny*nz / NFLOAT_PER_LINE + 1;

	size_t index = 0;	//linear index for image data array
	int nitems_in_line = 0;
	for (int i=0; i<number_lines; ++i) {
		char line[LINE_LENGTH];
		if (!fgets(line, sizeof(line), situsfile)) {
			printf("read situs file failed\n");
		}

		nitems_in_line = (int) (strlen(line) / FLOAT_SIZE);
		char * pline = line;
		for (int j=0; j<nitems_in_line; ++j) {
			sscanf(pline, "%f", &data[index]);
			pline += FLOAT_SIZE;
			++index;
		}
	}

	EXITFUNC;
	return 0;
}

int SitusIO::read_header(EMAN::Dict& dict, int, const EMAN::Region*, bool)
{
	ENTERFUNC;
	init();

	dict["nx"] = nx;
	dict["ny"] = ny;
	dict["nz"] = nz;

	dict["apix_x"] = apix;
	dict["apix_y"] = apix;
	dict["apix_z"] = apix;

	dict["origin_x"] = origx;
	dict["origin_y"] = origy;
	dict["origin_z"] = origz;

	EXITFUNC;
	return 0;
}

int SitusIO::write_header(const EMAN::Dict& dict, int, const EMAN::Region*, EMAN::EMUtil::EMDataType, bool)
{
	ENTERFUNC;
	init();

	apix = (float)dict["apix_x"];
	origx = (float)dict["origin_x"];
	origy = (float)dict["origin_y"];
	origz = (float)dict["origin_z"];
	nx = (int)dict["nx"];
	ny = (int)dict["ny"];
	nz = (int)dict["nz"];

	char headerline[LINE_LENGTH];
	sprintf(headerline, "%.6f %.6f %.6f %.6f %d %d %d", apix, origx, origy, origz, nx, ny, nz);

	if(!fputs(headerline, situsfile)) {
		printf("Write situs header failed\n");
	}

	if(!fputs("\n\n", situsfile)) {
		printf("Write situs header failed\n");
	}

	EXITFUNC;
	return 0;
}

int SitusIO::write_data(float* data, int, const EMAN::Region*, EMAN::EMUtil::EMDataType, bool)
{
	ENTERFUNC;

	for (size_t index=0; index<(size_t)nx*ny*nz; ++index) {
		fprintf(situsfile, OUTFORMAT, data[index]);
		if((index+1)%NFLOAT_PER_LINE == 0) {
			fputs("\n", situsfile);
		}
	}

	EXITFUNC;
	return 0;
}

bool SitusIO::is_valid(const void *first_block)
{
	ENTERFUNC;
	if (!first_block) {
		return false;
	}

	char *buf = (char *)(first_block);
	string line1 = Util::get_line_from_string(&buf);

	if(line1.size()==0) return false;

	float apix, origx, origy, origz;
	int nx, ny, nz;

	if(sscanf(line1.c_str(), "%f %f %f %f %d %d %d", &apix, &origx, &origy, &origz, &nx, &ny, &nz) != 7) return false;

	if(apix<0.01 || apix>100) return false;
	if(nx<=0 || ny<0 || nz<0) return false;

	EXITFUNC;
	return true;
}

void SitusIO::flush()
{
	fflush(situsfile);
}

bool SitusIO::is_image_big_endian()
{
	return ByteOrder::is_host_big_endian();
}

bool SitusIO::is_complex_mode()
{
	return false;
}
