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

int ImageIO::check_region(const Region * area, const IntSize & max_size)
{
	int img_ndim = max_size.get_ndim();
	if (area) {
		int area_ndim = area->get_ndim();

		if (area_ndim > img_ndim) {
			LOGERR("Image is %dD. Cannot read %dD region", img_ndim, area_ndim);
			return 1;
		}
#if 0
		if (!area->inside_region(max_size)) {
			LOGERR("Region box %s is outside image area->",
				   area->get_string().c_str());
			return 1;
		}
#endif

	}
	return 0;
}

int ImageIO::check_read_access(int image_index, bool check_data, float *data)
{
	if (init() != 0) {
		return 1;
	}

	if (check_data) {
		if (!data) {
			throw NullPointerException("image data is NULL");
		}
	}

	int nimg = get_nimg();
	if (image_index < 0 || image_index >= nimg) {
		throw OutofRangeException(0, nimg-1, image_index, "image index");
	}

	return 0;
}

int ImageIO::check_write_access(IOMode iomode, int image_index, int max_nimg,
								bool check_data, float *data)
{
	if (init() != 0) {
		return 1;
	}

	if (iomode == READ_ONLY) {
		LOGERR("cannot write header. File is not openned to write");
		return 1;
	}

	if ((image_index < -1) || (max_nimg > 0 && image_index >= max_nimg)) {
		throw OutofRangeException(-1, max_nimg - 1, image_index, "image index");
	}
	
	if (check_data) {
		if (!data) {
			throw NullPointerException("image data is NULL");
		}
	}

	return 0;
}

FILE *ImageIO::sfopen(string filename, IOMode mode, bool * is_new, bool overwrite)
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
					LOGERR("cannot create file '%s'", filename.c_str());
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
		LOGERR("cannot access file '%s'", filename.c_str());
	}
	return f;
}
