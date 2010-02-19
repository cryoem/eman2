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
#include <climits>
#include "vtkio.h"
#include "util.h"
#include "geometry.h"
#include "portable_fileio.h"

using namespace EMAN;

const char *VtkIO::MAGIC = "# vtk DataFile Version";

VtkIO::VtkIO(const string & vtk_filename, IOMode rw)
:	filename(vtk_filename), rw_mode(rw), vtk_file(0), initialized(false)
{
	is_big_endian = ByteOrder::is_host_big_endian();
	is_new_file = false;

	datatype = DATATYPE_UNKNOWN;
	filetype = VTK_UNKNOWN;
	nx = 0;
	ny = 0;
	nz = 0;
	originx = 0;
	originy = 0;
	originz = 0;
	spacingx = 0;
	spacingy = 0;
	spacingz = 0;
	file_offset = 0;
}

VtkIO::~VtkIO()
{
	if (vtk_file) {
		fclose(vtk_file);
		vtk_file = 0;
	}
}

static int samestr(const char *s1, const char *s2)
{
	return (strncmp(s1, s2, strlen(s2)) == 0);
}


void VtkIO::init()
{
	if (initialized) {
		return;
	}
	ENTERFUNC;
	initialized = true;

	vtk_file = sfopen(filename, rw_mode, &is_new_file);

	if (!is_new_file) {
		char buf[1024];
		int bufsz = sizeof(buf);
		if (fgets(buf, bufsz, vtk_file) == 0) {
			throw ImageReadException(filename, "first block");
		}

		if (!is_valid(&buf)) {
			throw ImageReadException(filename, "invalid VTK");
		}

		if (fgets(buf, bufsz, vtk_file) == 0) {
			throw ImageReadException(filename, "read VTK file failed");
		}

		if (fgets(buf, bufsz, vtk_file)) {
			if (samestr(buf, "ASCII")) {
				filetype = VTK_ASCII;
			}
			else if (samestr(buf, "BINARY")) {
				filetype = VTK_BINARY;
			}
		}
		else {
			throw ImageReadException(filename, "read VTK file failed");
		}

		if (fgets(buf, bufsz, vtk_file)) {
			if (samestr(buf, "DATASET")) {
				char dataset_name[128];
				sscanf(buf, "DATASET %s", dataset_name);
				DatasetType ds_type = get_datasettype_from_name(dataset_name);
				read_dataset(ds_type);
			}
		}
		else {
			throw ImageReadException(filename, "read VTK file failed");
		}

		while (fgets(buf, bufsz, vtk_file)) {
			if (samestr(buf, "SCALARS")) {
				char datatypestr[32];
				char scalartype[32];
				sscanf(buf, "SCALARS %s %s", scalartype, datatypestr);

				datatype = get_datatype_from_name(datatypestr);
				if (datatype != UNSIGNED_SHORT && datatype != FLOAT) {
					string desc = "unknown data type: " + string(datatypestr);
					throw ImageReadException(filename, desc);
				}
			}
			else if (samestr(buf, "LOOKUP_TABLE")) {
				char tablename[128];
				sscanf(buf, "LOOKUP_TABLE %s", tablename);
				if (!samestr(tablename, "default")) {
					throw ImageReadException(filename, "only default LOOKUP_TABLE supported");
				}
				else {
					break;
				}
			}
		}
#if 0
		if (filetype == VTK_BINARY) {
			throw ImageReadException(filename, "binary VTK is not supported");
		}
#endif
		file_offset = portable_ftell(vtk_file);
	}
	EXITFUNC;
}


bool VtkIO::is_valid(const void *first_block)
{
	ENTERFUNC;
	bool result = false;
	if (first_block) {
		result = Util::check_file_by_magic(first_block, MAGIC);
	}
	EXITFUNC;
	return result;
}

int VtkIO::read_header(Dict & dict, int image_index, const Region * area, bool)
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
	check_region(area, IntSize(nx, ny, nz));

	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen, nz, &zlen);

	dict["nx"] = xlen;
	dict["ny"] = ylen;
	dict["nz"] = zlen;

	dict["datatype"] = to_em_datatype(datatype);

	dict["apix_x"] = spacingx;
	dict["apix_y"] = spacingy;
	dict["apix_z"] = spacingz;

	dict["origin_x"] = originx;
	dict["origin_y"] = originy;
	dict["origin_z"] = originz;

	EXITFUNC;
	return 0;
}

int VtkIO::write_header(const Dict & dict, int image_index, const Region*,
						EMUtil::EMDataType, bool)
{
	ENTERFUNC;

	//single image format, index can only be zero
	if(image_index != 0) {
		throw ImageWriteException(filename, "VTK file does not support stack.");
	}
	check_write_access(rw_mode, image_index);

	nx = dict["nx"];
	ny = dict["ny"];
	nz = dict["nz"];

	originx = dict["origin_x"];
	originy = dict["origin_y"];
	originz = dict["origin_z"];

	spacingx = dict["apix_x"];
	spacingy = dict["apix_y"];
	spacingz = dict["apix_z"];

	fprintf(vtk_file, "# vtk DataFile Version 2.0\n");
	fprintf(vtk_file, "EMAN\n");
	fprintf(vtk_file, "BINARY\n");
	fprintf(vtk_file, "DATASET STRUCTURED_POINTS\n");
	fprintf(vtk_file, "DIMENSIONS %0d %0d %0d\nORIGIN %f %f %f\nSPACING %f %f %f\n",
			nx, ny, nz, originx, originy, originz, spacingx, spacingy, spacingz);


	fprintf(vtk_file, "POINT_DATA %0d\nSCALARS density float 1\nLOOKUP_TABLE default\n",
			nx * ny * nz);
	EXITFUNC;
	return 0;
}

int VtkIO::read_data(float *data, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	//single image format, index can only be zero
	image_index = 0;
	check_read_access(image_index, data);

	if (area) {
		LOGWARN("read VTK region is not supported yet. Read whole image instead.");
	}

	portable_fseek(vtk_file, file_offset, SEEK_SET);

	int xlen = 0, ylen = 0, zlen = 0;
	int x0 = 0, y0 = 0, z0 = 0;
	EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen, nz, &zlen);
	EMUtil::get_region_origins(area, &x0, &y0, &z0, nz, image_index);

	if (filetype == VTK_ASCII) {

		int bufsz = nx * get_mode_size(datatype) * CHAR_BIT;
		char *buf = new char[bufsz];
		int i = 0;

		while (fgets(buf, bufsz, vtk_file)) {
			size_t bufslen = strlen(buf) - 1;
			char numstr[32];
			int k = 0;
			for (size_t j = 0; j < bufslen; j++) {
				if (!isspace(buf[j])) {
					numstr[k++] = buf[j];
				}
				else {
					numstr[k] = '\0';
					data[i++] = (float)atoi(numstr);
					k = 0;
				}
			}
		}
		if( buf )
		{
			delete[]buf;
			buf = 0;
		}
	}
	else if (filetype == VTK_BINARY) {
		int nxy = nx * ny;
		int row_size = nx * get_mode_size(datatype);

		for (int i = 0; i < nz; i++) {
			int i2 = i * nxy;
			for (int j = 0; j < ny; j++) {
				fread(&data[i2 + j * nx], row_size, 1, vtk_file);
			}
		}

		if (!ByteOrder::is_host_big_endian()) {
			ByteOrder::swap_bytes(data, nx * ny * nz);
		}
	}

	EXITFUNC;
	return 0;
}

int VtkIO::write_data(float *data, int image_index, const Region* ,
					  EMUtil::EMDataType, bool)
{
	ENTERFUNC;

	//single image format, index can only be zero
	image_index = 0;
	check_write_access(rw_mode, image_index, 1, data);

	bool swapped = false;
	if (!ByteOrder::is_host_big_endian()) {
		ByteOrder::swap_bytes(data, nx * ny * nz);
		swapped = true;
	}

	fwrite(data, nx * nz, ny * sizeof(float), vtk_file);

	if (swapped) {
		ByteOrder::swap_bytes(data, nx * ny * nz);
	}
	EXITFUNC;
	return 0;
}

void VtkIO::flush()
{
	fflush(vtk_file);
}

bool VtkIO::is_complex_mode()
{
	return false;
}

bool VtkIO::is_image_big_endian()
{
	return true;
}

int VtkIO::to_em_datatype(int vtk_datatype)
{
	DataType d = static_cast < DataType > (vtk_datatype);
	switch (d) {
	case UNSIGNED_SHORT:
		return EMUtil::EM_USHORT;
	case FLOAT:
		return EMUtil::EM_FLOAT;
	default:
		break;
	}
	return EMUtil::EM_UNKNOWN;
}


int VtkIO::get_mode_size(DataType d)
{
	switch (d) {
	case UNSIGNED_CHAR:
	case CHAR:
		return sizeof(char);
	case UNSIGNED_SHORT:
	case SHORT:
		return sizeof(short);
	case UNSIGNED_INT:
	case INT:
		return sizeof(int);
	case UNSIGNED_LONG:
	case LONG:
		return sizeof(long);
	case FLOAT:
		return sizeof(float);
	case DOUBLE:
		return sizeof(double);
	default:
		LOGERR("don't support this data type '%d'", d);
		break;
	}
	return 0;
}

VtkIO::DataType VtkIO::get_datatype_from_name(const string& datatype_name)
{
	static bool initialized = false;
	static map < string, VtkIO::DataType > datatypes;

	if (!initialized) {
		datatypes["bit"] = BIT;

		datatypes["unsigned_char"] = UNSIGNED_CHAR;
		datatypes["char"] = CHAR;

		datatypes["unsigned_short"] = UNSIGNED_SHORT;
		datatypes["short"] = SHORT;

		datatypes["unsigned_int"] = UNSIGNED_INT;
		datatypes["int"] = INT;

		datatypes["unsigned_long"] = UNSIGNED_LONG;
		datatypes["long"] = LONG;

		datatypes["float"] = FLOAT;
		datatypes["double"] = DOUBLE;
		initialized = true;
	}

	DataType result = DATATYPE_UNKNOWN;

	if (datatypes.find(datatype_name) != datatypes.end()) {
		result = datatypes[datatype_name];
	}
	return result;
}

VtkIO::DatasetType VtkIO::get_datasettype_from_name(const string& dataset_name)
{

	static bool initialized = false;
	static map < string, DatasetType > types;

	if (!initialized) {
		types["STRUCTURED_POINTS"] = STRUCTURED_POINTS;
		types["STRUCTURED_GRID"] = STRUCTURED_GRID;
		types["RECTILINEAR_GRID"] = RECTILINEAR_GRID;
		types["UNSTRUCTURED_GRID"] = UNSTRUCTURED_GRID;
		types["POLYDATA"] = POLYDATA;
	}

	DatasetType result = DATASET_UNKNOWN;
	if (types.find(dataset_name) != types.end()) {
		result = types[dataset_name];
	}
	return result;
}

void VtkIO::read_dataset(DatasetType dstype)
{
	char buf[1024];
	int bufsz = sizeof(buf);

	if (dstype == STRUCTURED_POINTS) {
		int nlines = 3;
		int i = 0;
		while (i < nlines && fgets(buf, bufsz, vtk_file)) {
			if (samestr(buf, "DIMENSIONS")) {
				sscanf(buf, "DIMENSIONS %d %d %d", &nx, &ny, &nz);
			}
			else if (samestr(buf, "ORIGIN")) {
				sscanf(buf, "ORIGIN %f %f %f", &originx, &originy, &originz);
			}
			else if (samestr(buf, "SPACING") || samestr(buf, "ASPECT_RATIO")) {
				if (samestr(buf, "SPACING")) {
					sscanf(buf, "SPACING %f %f %f", &spacingx, &spacingy, &spacingz);
				}
				else {
					sscanf(buf, "ASPECT_RATIO %f %f %f", &spacingx, &spacingy, &spacingz);
				}

				if (spacingx != spacingy || spacingx != spacingz || spacingy != spacingz) {
					throw ImageReadException(filename,
											 "not support non-uniform spacing VTK so far\n");
				}
			}
			i++;
		}

		if (i != nlines) {
			throw ImageReadException(filename, "read VTK file failed");
		}
	}
	else {
		throw ImageReadException(filename, "only STRUCTURED_POINTS is supported so far");
	}
}

