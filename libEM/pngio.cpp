/**
 * $Id$
 */
#ifdef EM_PNG

#include "pngio.h"
#include "log.h"
#include "geometry.h"
#include "emutil.h"
#include "util.h"


using namespace EMAN;

PngIO::PngIO(const string & file, IOMode rw)
:	filename(file), rw_mode(rw), png_file(0), initialized(false)
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
	
	check_read_access(image_index);

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

int PngIO::write_header(const Dict & dict, int image_index, const Region* area,
						EMUtil::EMDataType, bool)
{
	ENTERFUNC;

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
	EXITFUNC;
	return 0;
}

int PngIO::read_data(float *data, int image_index, const Region * area, bool)
{
	ENTERFUNC;

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

	delete[]cdata;
	cdata = 0;

	png_read_end(png_ptr, end_info);
	EXITFUNC;
	return 0;
}

int PngIO::write_data(float *data, int image_index, const Region* area,
					  EMUtil::EMDataType, bool)
{
	ENTERFUNC;

	check_write_access(rw_mode, image_index, 1, data);

	if (depth_type == PNG_CHAR_DEPTH) {
		unsigned char *cdata = new unsigned char[nx];

		for (unsigned int y = 0; y < ny; y++) {
			for (unsigned int x = 0; x < nx; x++) {
				cdata[x] = static_cast < unsigned char >(data[y * nx + x]);
			}
			png_write_row(png_ptr, (png_byte *) cdata);
		}

		delete[]cdata;
		cdata = 0;
	}
	else if (depth_type == PNG_SHORT_DEPTH) {
		unsigned short *sdata = new unsigned short[nx];

		for (unsigned int y = 0; y < ny; y++) {
			for (unsigned int x = 0; x < nx; x++) {
				sdata[x] = static_cast < unsigned short >(data[y * nx + x]);
			}
			png_write_row(png_ptr, (png_byte *) sdata);
		}

		delete[]sdata;
		sdata = 0;
	}

	png_write_end(png_ptr, info_ptr);
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


#endif
