/**
 * $Id$
 */
#include "imageio.h"
#include "geometry.h"

using namespace EMAN;

ImageIO::~ImageIO()
{
}

int ImageIO::read_ctf(Ctf &, int)
{
	return 1;
}

void ImageIO::write_ctf(const Ctf &, int)
{

}

void ImageIO::check_region(const Region * area, const FloatSize & max_size,
						   bool is_new_file)
{
	if (area) {
		if (is_new_file) {
			throw ImageReadException("", "file must exist before accessing its region");
		}
		int img_ndim = max_size.get_ndim();
		int area_ndim = area->get_ndim();
		
		if (area_ndim > img_ndim) {
			char desc[256];
			sprintf(desc, "Image is %dD. Cannot read %dD region", img_ndim, area_ndim);
			throw ImageReadException("", desc);
		}
		
		if (!area->is_region_in_box(max_size)) {
			char desc[1024];
			sprintf(desc, "Region box %s is outside image area (%d,%d,%d)",
					area->get_string().c_str(), (int)max_size[0],
					(int)max_size[1], (int)max_size[2]);
			throw ImageReadException("", desc);
		}
	}
}

void ImageIO::check_region(const Region * area, const IntSize & max_size,
						   bool is_new_file)
{
	check_region(area, FloatSize(max_size[0], max_size[1], max_size[2]),
				 is_new_file);
}

void ImageIO::check_read_access(int image_index)
{
	init();

	int nimg = get_nimg();
	if (image_index < 0 || image_index >= nimg) {
		throw OutofRangeException(0, nimg-1, image_index, "image index");
	}
}

void ImageIO::check_read_access(int image_index, const float *data)
{
	check_read_access(image_index);
	if (!data) {
		throw NullPointerException("image data is NULL");
	}
}

void ImageIO::check_write_access(IOMode iomode, int image_index, int max_nimg)
{
	init();

	if (iomode == READ_ONLY) {
		throw ImageWriteException("", "File is not openned to write");
	}

	if ((image_index < -1) || (max_nimg > 0 && image_index >= max_nimg)) {
		throw OutofRangeException(-1, max_nimg - 1, image_index, "image index");
	}
}

void ImageIO::check_write_access(IOMode iomode, int image_index,
								 int max_nimg, const float *data)
{
	check_write_access(iomode, image_index, max_nimg);
	if (!data) {
		throw NullPointerException("image data is NULL");
	}
}

FILE *ImageIO::sfopen(const string & filename, IOMode mode,
					  bool * is_new, bool overwrite)
{
	FILE *f = 0;
	if (mode == READ_ONLY) {
		f = fopen(filename.c_str(), "rb");
	}
	else if (mode == READ_WRITE) {
		if (overwrite) {
			f = fopen(filename.c_str(), "wb");
			if (is_new) {
				*is_new = true;
			}
		}
		else {
			f = fopen(filename.c_str(), "r+b");
			if (!f) {
				FILE *f1 = fopen(filename.c_str(), "wb");
				if (!f1) {
					throw FileAccessException(filename);
				}
				else {
					if (is_new) {
						*is_new = true;
					}
					fclose(f1);
					f1 = 0;
					f = fopen(filename.c_str(), "r+b");
				}
			}
		}
	}
	else if (mode == WRITE_ONLY) {
		f = fopen(filename.c_str(), "wb");
		if (is_new) {
			*is_new = true;
		}
	}
	
	if (!f) {
		throw FileAccessException(filename);
	}
	return f;
}


int ImageIO::get_nimg()
{
	init();
	return 1;
}
