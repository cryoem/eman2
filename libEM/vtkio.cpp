/**
 * $Id$
 */
#include "vtkio.h"
#include "emutil.h"
#include "log.h"
#include "util.h"
#include "geometry.h"
#include "portable_fileio.h"
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

using namespace EMAN;

const char *VtkIO::MAGIC = "# vtk DataFile Version";

VtkIO::VtkIO(string vtk_filename, IOMode rw)
    :  filename(vtk_filename), rw_mode(rw), vtk_file(0), initialized(false)
{
    is_big_endian = ByteOrder::is_machine_big_endian();

    datatype = DATA_UNKNOWN;
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


int VtkIO::init()
{
    static int err = 0;
    if (initialized) {
	return err;
    }
    Log::logger()->log("VtkIO::init()");
    initialized = true;

    vtk_file = sfopen(filename, rw_mode);
    if (!vtk_file) {
	err = 1;
	return 1;
    }

    char buf[1024];
    int bufsz = sizeof(buf);
    if (fgets(buf, bufsz, vtk_file) == 0) {
	Log::logger()->error("read VTK file '%s' failed", filename.c_str());
	return 1;
    }

    if (!is_valid(&buf)) {
	Log::logger()->error("not a valid VTK file");
	return 1;
    }

    while (fgets(buf, bufsz, vtk_file)) {
	if (samestr(buf, "LOOKUP_TABLE default")) {
	    break;
	}
	else if (samestr(buf, "ASCII")) {
	    filetype = VTK_ASCII;
	}
	else if (samestr(buf, "BINARY")) {
	    filetype = VTK_BINARY;
	}
	else if (samestr(buf, "DIMENSIONS")) {
	    sscanf(buf, "DIMENSIONS %d %d %d", &nx, &ny, &nz);
	}
	else if (samestr(buf, "ORIGIN")) {
	    sscanf(buf, "ORIGIN %f %f %f", &originx, &originy, &originz);
	}
	else if (samestr(buf, "SPACING")) {
	    sscanf(buf, "SPACING %f %f %f", &spacingx, &spacingy, &spacingz);
	    if (spacingx != spacingy || spacingx != spacingz || spacingy != spacingz) {
		Log::logger()->error("cannot handle non-uniform spacing VTK file\n");
		return 1;
	    }
	}
	else if (samestr(buf, "SCALARS")) {
	    char datatypestr[32];
	    char scalartype[32];
	    sscanf(buf, "SCALARS %s %s", scalartype, datatypestr);

	    if (samestr(datatypestr, "unsigned_short")) {
		datatype = UNSIGNED_SHORT;
	    }
	    else if (samestr(datatypestr, "float")) {
		datatype = FLOAT;
	    }
	    else {
		Log::logger()->error("unknown data type: %s", datatypestr);
		return 1;
	    }
	}
    }

    if (filetype == VTK_BINARY) {
	Log::logger()->error("binary VTK is not supported");
	return 1;
    }

    file_offset = portable_ftell(vtk_file);

    return 0;
}


bool VtkIO::is_valid(const void *first_block)
{
    Log::logger()->log("VtkIO::is_valid()");
    return Util::check_file_by_magic(first_block, MAGIC);
}

int VtkIO::read_header(Dict & dict, int image_index, const Region * area, bool)
{
    Log::logger()->log("VtkIO::read_header() from file '%s'", filename.c_str());

    if (check_read_access(image_index) != 0) {
	return 1;
    }
    if (check_region(area, Size(nx, ny, nz)) != 0) {
	return 1;
    }

    int xlen = 0, ylen = 0, zlen = 0;
    EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen, nz, &zlen);

    dict["nx"] = EMObject(xlen);
    dict["ny"] = EMObject(ylen);
    dict["nz"] = EMObject(zlen);

    dict["datatype"] = EMObject(to_em_datatype(datatype));

    dict["spacing_row"] = EMObject(spacingx);
    dict["spacing_col"] = EMObject(spacingy);
    dict["spacing_sec"] = EMObject(spacingz);

    dict["origin_row"] = EMObject(originx);
    dict["origin_col"] = EMObject(originy);
    dict["origin_sec"] = EMObject(originz);

    return 0;
}

int VtkIO::write_header(const Dict & dict, int image_index)
{
    Log::logger()->log("VtkIO::write_header() to file '%s'", filename.c_str());
    if (check_write_access(rw_mode, image_index) != 0) {
	return 1;
    }

    nx = dict["nx"].get_int();
    ny = dict["ny"].get_int();
    nz = dict["nz"].get_int();

    originx = dict["origin_row"].get_float();
    originy = dict["origin_col"].get_float();
    originz = dict["origin_sec"].get_float();

    spacingx = dict["spacing_row"].get_float();
    spacingy = dict["spacing_col"].get_float();
    spacingz = dict["spacing_sec"].get_float();

    fprintf(vtk_file, "# vtk DataFile Version 2.0\n");
    fprintf(vtk_file, "EMAN\n");
    fprintf(vtk_file, "BINARY\n");
    fprintf(vtk_file, "DATASET STRUCTURED_POINTS\n");
    fprintf(vtk_file, "DIMENSIONS %0d %0d %0d\nORIGIN %f %f %f\nSPACING %f %f %f\n",
	    nx, ny, nz, originx, originy, originz, spacingx, spacingy, spacingz);
    
    
    fprintf(vtk_file, "POINT_DATA %0d\nSCALARS density float 1\nLOOKUP_TABLE default\n",
	    nx * ny * nz);
    return 0;
}

int VtkIO::read_data(float *data, int image_index, const Region * area, bool )
{
    Log::logger()->log("VtkIO::read_data() from file '%s'", filename.c_str());
    if (check_read_access(image_index, true, data) != 0) {
	return 1;
    }

    if (area) {
	Log::logger()->warn("read VTK region is not supported yet. Read whole image instead.");
    }

    portable_fseek(vtk_file, file_offset, SEEK_SET);

    int xlen = 0, ylen = 0, zlen = 0;
    int x0 = 0, y0 = 0, z0 = 0;
    EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen, nz, &zlen);
    EMUtil::get_region_origins(area, &x0, &y0, &z0, nz, image_index);

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
		data[i++] = atoi(numstr);
		k = 0;
	    }
	}
    }

    delete [] buf;
    buf = 0;

    return 0;
}

int VtkIO::write_data(float *data, int image_index)
{
    Log::logger()->log("VtkIO::write_data() to file '%s'", filename.c_str());
    if (check_write_access(rw_mode, image_index, true, data) != 0) {
	return 1;
    }

    bool swapped = false;
    if (!ByteOrder::is_machine_big_endian()) {
	ByteOrder::swap_bytes(data, nx * ny * nz);
	swapped = true;
    }

    fwrite(data, nx * nz, ny * sizeof(float), vtk_file);

    if (swapped) {
	ByteOrder::swap_bytes(data, nx * ny * nz);
    }

    return 0;
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
    DataType d = static_cast<DataType>(vtk_datatype);
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

int VtkIO::get_nimg()
{
    if (init() != 0) {
	return 0;
    }
    return 1;
}

int VtkIO::get_mode_size(DataType d)
{
    switch (d) {
    case UNSIGNED_SHORT:
	return sizeof(short);
    case FLOAT:
	return sizeof(float);
    default:
	break;
    }
    return 0;
}
