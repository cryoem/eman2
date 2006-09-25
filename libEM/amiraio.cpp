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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * 
 * */
 
#include "amiraio.h"
#include "util.h"

#ifndef WIN32
	#include <sys/param.h>
#else
	#include <windows.h>
	#define MAXPATHLEN (MAX_PATH*4)
#endif	//WIN32

using namespace EMAN;

const char *AmiraIO::MAGIC = "# AmiraMesh";

AmiraIO::AmiraIO(const string & file, IOMode rw)
:	filename(file), rw_mode(rw), amira_file(0), initialized(false)
{
	is_big_endian = true;
	nx = 0;
	ny = 0;
	nz = 0;
}

AmiraIO::~AmiraIO()
{
	if (amira_file) {
		fclose(amira_file);
		amira_file = 0;
	}
}

void AmiraIO::init()
{
	ENTERFUNC;
	if (initialized) {
		return;
	}
	
	initialized = true;
	bool is_new_file = false;
	amira_file = sfopen(filename, rw_mode, &is_new_file, true);

	if (!is_new_file) {
		char buf[MAXPATHLEN];
		if (!fgets(buf, MAXPATHLEN, amira_file)) {
			throw ImageReadException(filename, "Amira Header");
		}
		if (!is_valid(buf)) {
			throw ImageReadException(filename, "invalid Amira Mesh file");
		}
	}
	EXITFUNC;
}

bool AmiraIO::is_valid(const void *first_block)
{
	ENTERFUNC;
	bool result = false;
	if (!first_block) {
		result = false;
	}
	else {
		result = Util::check_file_by_magic(first_block, MAGIC);
	}
	EXITFUNC;
	return result;
}

int AmiraIO::read_header(Dict &, int, const Region *, bool)
{
	ENTERFUNC;
	LOGWARN("Amira read header is not supported.");
	EXITFUNC;
	return 1;

}

int AmiraIO::write_header(const Dict & dict, int image_index, const Region*, EMUtil::EMDataType, bool)
{
	ENTERFUNC;
	int err = 0;
	
	check_write_access(rw_mode, image_index, 1);
	
	nx = dict["nx"];
	ny = dict["ny"];
	nz = dict["nz"];
	
	float xorigin = 0.0f;		
	if(dict.has_key("origin_row")) xorigin = dict["origin_row"];
	float yorigin = 0.0f;
	if(dict.has_key("origin_col")) yorigin = dict["origin_col"];
	float zorigin = 0.0f;
	if(dict.has_key("origin_sec")) zorigin = dict["origin_sec"];
	float pixel = 0.0f;
	if(dict.has_key("pixel")) pixel = dict["pixel"];
			
	rewind(amira_file);
	
	if (fprintf(amira_file, "# AmiraMesh 3D BINARY 2.0\n\n") <= 0) {
		LOGERR("cannot write to AmiraMesh file '%s'", filename.c_str());
		err = 1;
	}
	else {
		fprintf(amira_file, "# Dimensions in x-, y-, and z-direction\n");
		fprintf(amira_file, "define Lattice %d %d %d\n\n", nx, ny, nz);
		fprintf(amira_file, "Parameters {\n\tCoordType \"uniform\",\n\t");
		fprintf(amira_file, "# BoundingBox is xmin xmax ymin ymax zmin zmax\n\t");
				
		fprintf(amira_file, "BoundingBox %f %f %f %f %f %f\n}\n\n",
				xorigin, xorigin + pixel * (nx - 1),
				yorigin, yorigin + pixel * (ny - 1), zorigin, zorigin + pixel * (nz - 1));
				
		fprintf(amira_file, "Lattice { float ScalarField } = @1\n\n@1\n");
	}

	EXITFUNC;
	return err;
}

int AmiraIO::read_data(float *, int, const Region *, bool)
{
	ENTERFUNC;
	LOGWARN("Amira read data is not supported.");
	EXITFUNC;
	return 1;
}

int AmiraIO::write_data(float *data, int image_index, const Region*, EMUtil::EMDataType, bool)
{
	ENTERFUNC;
	
	check_write_access(rw_mode, image_index, 1, data);
	ByteOrder::become_big_endian(data, nx * ny * nz);

	if (fwrite(data, nx * nz, ny * sizeof(float), amira_file) != ny * sizeof(float)) {
		throw ImageWriteException(filename, "incomplete file write in AmiraMesh file");
	}
	
	EXITFUNC;
	return 0;
}

void AmiraIO::flush()
{
	fflush(amira_file);
}

bool AmiraIO::is_complex_mode()
{
	return false;
}

bool AmiraIO::is_image_big_endian()
{
	init();
	return is_big_endian;
}

