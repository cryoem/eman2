/**
 * $Id$
 */
#include "imageio.h"
#include "log.h"
#include "geometry.h"


using namespace EMAN;

ImageIO::~ImageIO()
{
}

int ImageIO::read_ctf(Ctf &, int)
{
	return 1;
}

int ImageIO::write_ctf(const Ctf &, int)
{
	return 1;
}

void ImageIO::check_region(const Region * area, const IntSize & max_size)
{
	int img_ndim = max_size.get_ndim();
	if (area) {
		int area_ndim = area->get_ndim();
		
		if (area_ndim > img_ndim) {
			char desc[256];
			sprintf(desc, "Image is %dD. Cannot read %dD region", img_ndim, area_ndim);
			throw ImageReadException("", desc);
		}
#if 0
		if (!area->inside_region(max_size)) {
			LOGERR("Region box %s is outside image area->", area->get_string().c_str());			
		}
#endif

	}
}

void ImageIO::check_read_access(int image_index)
{
	init();

	int nimg = get_nimg();
	if (image_index < 0 || image_index >= nimg) {
		throw OutofRangeException(0, nimg-1, image_index, "image index");
	}
}

void ImageIO::check_read_access(int image_index, float *data)
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
								 int max_nimg, float *data)
{
	check_write_access(iomode, image_index, max_nimg);
	if (!data) {
		if (!data) {
			throw NullPointerException("image data is NULL");
		}
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
