/**
 * $Id$
 */
#ifdef EM_PNG

#include "pngio.h"
#include "log.h"
#include "geometry.h"
#include "emutil.h"
#include "util.h"
#include <assert.h>

using namespace EMAN;

PngIO::PngIO(string file, IOMode rw)
    :  filename(file), rw_mode(rw), png_file(0), initialized(false)
{
    png_ptr = 0;
    info_ptr = 0;
    end_info = 0;

    nx = 0;
    ny = 0;
    depth_type = PNG_INVALID_DEPTH;
    number_passes = 0;
}

PngIO::~PngIO()
{
    if (png_file) {
	fclose(png_file);
	png_file = 0;
    }

    png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
    png_ptr = 0;
    info_ptr = 0;
    end_info = 0;

}

int PngIO::init()
{
    static int err = 0;
    if (initialized) {
	return err;
    }
    Log::logger()->log("PngIO::init()");
    initialized = true;

    bool is_new_file = false;
    png_file = sfopen(filename, rw_mode, &is_new_file, true);

    if (!png_file) {
	err = 1;
	return err;
    }

    if (!is_new_file) {
	png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);
    }
    else {
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);
    }

    if (!png_ptr) {
	Log::logger()->error("cannot initialize libpng data structure");
	err = 1;
	return 1;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
	Log::logger()->error("cannot create png info data structure");
	err = 1;
	return 1;
    }


    end_info = png_create_info_struct(png_ptr);
    if (!end_info) {
	Log::logger()->error("cannot create png end info structure");
	err = 1;
	return 1;
    }

    if (setjmp(png_ptr->jmpbuf)) {
	Log::logger()->error("an error occurs within png");
	err = 1;
	return 1;
    }

    png_init_io(png_ptr, png_file);

    if (!is_new_file) {
	unsigned char header[PNG_BYTES_TO_CHECK];
	fread(header, sizeof(unsigned char), PNG_BYTES_TO_CHECK, png_file);
	if (!is_valid(header)) {
	    Log::logger()->error("Image '%s' not in png format.", filename.c_str());
	    err = 1;
	    return 1;
	}

	png_set_sig_bytes(png_ptr, PNG_BYTES_TO_CHECK);

	png_read_info(png_ptr, info_ptr);

	nx = png_get_image_width(png_ptr, info_ptr);
	ny = png_get_image_height(png_ptr, info_ptr);
	int bit_depth = png_get_bit_depth(png_ptr, info_ptr);
	int color_type = png_get_color_type(png_ptr, info_ptr);

	if (nx == 0 || ny == 0) {
	    Log::logger()->error("not a valid PNG file. width = %d, height = %d", nx, ny);
	    err = 1;
	    return 1;
	}

	if (bit_depth == CHAR_BIT) {
	    depth_type = PNG_CHAR_DEPTH;
	}
	else if (bit_depth == CHAR_BIT * sizeof(short)) {
	    depth_type = PNG_SHORT_DEPTH;
	}
	else {
	    depth_type = PNG_INVALID_DEPTH;
	    Log::logger()->error("sorry, I don't know how to handle png with depth = %d bit",
				 bit_depth);
	    err = 1;
	    return 1;
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

    return 0;
}

bool PngIO::is_valid(const void *first_block)
{
    Log::logger()->log("PngIO::is_valid()");
    if (!first_block) {
	return false;
    }

    if (png_sig_cmp((png_byte *) first_block, (png_size_t) 0, PNG_BYTES_TO_CHECK) == 0) {
	return true;
    }
    return false;
}

int PngIO::read_header(Dict & dict, int image_index, const Region * area, bool is_3d)
{
    Log::logger()->log("PngIO::read_header() from file '%s'", filename.c_str());

    if (check_read_access(image_index) != 0) {
	return 1;
    }
    assert(is_3d == false);
    int nx1 = static_cast<int>(nx);
    int ny1 = static_cast<int>(ny);
    if (check_region(area, Size(nx1, ny1)) != 0) {
	return 1;
    }

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
	Log::logger()->error("invalid PNG bit depth. don't know how to handle this png type");
    }

    return 0;
}

int PngIO::write_header(const Dict & dict, int image_index)
{
    Log::logger()->log("PngIO::write_header() to file '%s'", filename.c_str());
    if (check_write_access(rw_mode, image_index) != 0) {
	return 1;
    }

    nx = (png_uint_32) dict["nx"].get_int();
    ny = (png_uint_32) dict["ny"].get_int();
    int nz = dict["nz"].get_int();

    assert(nz == 1);

    int bit_depth = 0;
    EMUtil::EMDataType datatype = (EMUtil::EMDataType) dict["datatype"].get_int();

    if (datatype == EMUtil::EM_UCHAR) {
	depth_type = PNG_CHAR_DEPTH;
	bit_depth = CHAR_BIT;
    }
    else {
	if (datatype != EMUtil::EM_USHORT) {
	    Log::logger()->warn("Don't support data type '%s' in PNG. Convert to '%s'.",
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

    return 0;
}

int PngIO::read_data(float *data, int image_index, const Region * area, bool is_3d)
{
    Log::logger()->log("PngIO::read_data() from file '%s'", filename.c_str());

    if (check_read_access(image_index, true, data) != 0) {
	return 1;
    }
    assert(!is_3d);
    int nx1 = static_cast<int>(nx);
    int ny1 = static_cast<int>(ny);

    if (check_region(area, Size(nx1, ny1)) != 0) {
	return 1;
    }

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
		data[k] = static_cast<float>(cdata[x]);
		k++;
	    }
	}
	else if (depth_type == PNG_SHORT_DEPTH) {
	    for (int x = x0; x < x0 + xlen; x++) {
		data[k] = static_cast<float>(sdata[x]);
		k++;
	    }
	}
    }

    //Util::flip_image(data, nx, ny);

    delete [] cdata;
    cdata = 0;

    png_read_end(png_ptr, end_info);

    return 0;
}

int PngIO::write_data(float *data, int image_index)
{
    Log::logger()->log("PngIO::write_data() to file '%s'", filename.c_str());
    if (check_write_access(rw_mode, image_index, true, data) != 0) {
	return 1;
    }

    if (depth_type == PNG_CHAR_DEPTH) {
	unsigned char *cdata = new unsigned char[nx];

	for (unsigned int y = 0; y < ny; y++) {
	    for (unsigned int x = 0; x < nx; x++) {
		cdata[x] = static_cast<unsigned char>(data[y * nx + x]);
	    }
	    png_write_row(png_ptr, (png_byte *) cdata);
	}

	delete [] cdata;
	cdata = 0;
    }
    else if (depth_type == PNG_SHORT_DEPTH) {
	unsigned short *sdata = new unsigned short[nx];

	for (unsigned int y = 0; y < ny; y++) {
	    for (unsigned int x = 0; x < nx; x++) {
		sdata[x] = static_cast<unsigned short>(data[y * nx + x]);
	    }
	    png_write_row(png_ptr, (png_byte *) sdata);
	}

	delete [] sdata;
	sdata = 0;
    }

    png_write_end(png_ptr, info_ptr);

    return 0;
}

bool PngIO::is_complex_mode()
{
    return false;
}

bool PngIO::is_image_big_endian()
{
    return true;
}

int PngIO::get_nimg()
{
    if (init() != 0) {
	return 0;
    }

    return 1;
}

#endif
