/**
 * $Id$
 */
#include "amiraio.h"
#include "log.h"
#include "util.h"
#include <sys/param.h>

using namespace EMAN;

const char *AmiraIO::MAGIC = "# AmiraMesh";

AmiraIO::AmiraIO(string file, IOMode rw)
    :  filename(file), rw_mode(rw), amira_file(0), initialized(false)
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

int AmiraIO::init()
{
    static int err = 0;
    if (initialized) {
	return err;
    }
    Log::logger()->log("AmiraIO::init()");
    initialized = true;

    bool is_new_file = false;
    amira_file = sfopen(filename, rw_mode, &is_new_file, true);

    if (!amira_file) {
	err = 1;
	return err;
    }

    if (!is_new_file) {
	char buf[MAXPATHLEN];
	if (!fgets(buf, MAXPATHLEN, amira_file)) {
	    Log::logger()->error("cannot read from file '%s'", filename.c_str());
	    err = 1;
	    return err;
	}
	if (!is_valid(buf)) {
	    Log::logger()->error("'%s' is not a valid Amira Mesh file");
	    err = 1;
	    return err;
	}
    }

    return 0;
}

bool AmiraIO::is_valid(const void *first_block)
{
    Log::logger()->log("AmiraIO::is_valid()");
    if (!first_block) {
	return false;
    }
    return Util::check_file_by_magic(first_block, MAGIC);
}

int AmiraIO::read_header(Dict & , int , const Region * , bool )
{
    Log::logger()->log("AmiraIO::read_header() from file '%s'", filename.c_str());
    Log::logger()->warn("Amira read header is not supported.");
    return 1;

}

int AmiraIO::write_header(const Dict & dict, int image_index)
{
    Log::logger()->log("AmiraIO::write_header() to file '%s'", filename.c_str());
    if (check_write_access(rw_mode, image_index) != 0) {
	return 1;
    }

    nx = dict["nx"].get_int();
    ny = dict["ny"].get_int();
    nz = dict["nz"].get_int();

    float xorigin = dict["origin_row"].get_float();
    float yorigin = dict["origin_col"].get_float();
    float zorigin = dict["origin_sec"].get_float();
    float pixel = dict["pixel"].get_float();

    if (fprintf(amira_file, "# AmiraMesh 3D BINARY 2.0\n\n") <= 0) {
	Log::logger()->error("cannot write to AmiraMesh file '%s'", filename.c_str());
	return 1;
    }

    fprintf(amira_file, "# Dimensions in x-, y-, and z-direction\n");
    fprintf(amira_file, "define Lattice %d %d %d\n\n", nx, ny, nz);
    fprintf(amira_file, "Parameters {\n\tCoordType \"uniform\",\n\t");
    fprintf(amira_file, "# BoundingBox is xmin xmax ymin ymax zmin zmax\n\t");

    fprintf(amira_file, "BoundingBox %f %f %f %f %f %f\n}\n\n",
	    xorigin, xorigin + pixel * (nx - 1),
	    yorigin, yorigin + pixel * (ny - 1),
	    zorigin, zorigin + pixel * (nz - 1));

    fprintf(amira_file, "Lattice { float ScalarField } = @1\n\n@1\n");
    return 0;
}

int AmiraIO::read_data(float *, int , const Region * , bool )
{
    Log::logger()->log("AmiraIO::read_data() from file '%s'", filename.c_str());
    Log::logger()->warn("Amira read data is not supported.");
    return 1;
}

int AmiraIO::write_data(float *data, int image_index)
{
    Log::logger()->log("AmiraIO::write_data() to file '%s'", filename.c_str());
    if (check_write_access(rw_mode, image_index, true, data) != 0) {
	return 1;
    }

    ByteOrder::become_big_endian(data, nx * ny * nz);

    if (fwrite(data, nx * nz, ny * sizeof(float), amira_file) != ny * sizeof(float)) {
	Log::logger()->error("incomplete file write in AmiraMesh file");
	return 1;
    }

    return 0;
}

bool AmiraIO::is_complex_mode()
{
    return false;
}

bool AmiraIO::is_image_big_endian()
{
    return is_big_endian;
}

int AmiraIO::get_nimg()
{
    if (init() != 0) {
	return 0;
    }

    return 1;
}
