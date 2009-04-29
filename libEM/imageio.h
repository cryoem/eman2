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

#ifndef eman__imageio_h__
#define eman__imageio_h__ 1

#include "byteorder.h"
#include "emutil.h"

using std::vector;
using std::string;

namespace EMAN
{
	class Region;
	class FloatSize;
	class IntSize;
	class Ctf;

	/** ImageIO classes are designed for reading/writing various
	 * electron micrography image formats, including MRC, IMAGIC,
	 * SPIDER, PIF, etc.
	 *
	 * ImageIO class is the base class for all image io classes.
	 * Each subclass defines the IO routines for reading/writing a
	 * type of image format. For example, MrcIO is for reading/writing
	 * MRC images.
	 *
	 * A subclass should implement functions declared in
	 * DEFINE_IMAGEIO_FUNC macro.
	 *
	 * Image header I/O is separated from image data I/O.
	 *
	 * Some image formats (e.g., MRC, DM3, Gatan2, TIFF, PNG, EM, ICOS) only store a
	 * single 2D or 3D images. For them, image_index should always be
	 * 0. To read/write part of the image, use a region.
	 *
	 * Some image formats (e.g., Imagic) can store either a
	 * stack of 2D images or a single 3D image. From the physical
	 * image file itself, there is no difference between them. Only at
	 * reading time, they may be treated as either a stack of 2D
	 * images, or a single 3D image. A boolean hint should be given for these
	 * formats at reading time.
	 *
	 * Some image formats (e.g. HDF, PIF) can store a stack of images. Each
	 * image can be 2D or 3D.
	 *
	 * For image formats storing multiple images, valid image_index = [0, n].
	 *
	 * For image formats storing multiple images, the image append and
	 * insertion policy is:
	 *
	 *   - it should support appending image to existing file.
	 *
	 *   - it should support appending image to new file.
	 *
	 *   - it should support insert image in existing file with gap
	 *     between EOF and the new image. The gap should be zeroed.
	 *
	 *   - it should support insert image in new file with image index != 0.
	 *
	 *   - insert image in existing file overwriting existing image.
	 *
	 *   - The gap from insertion or appending should be filled with zero.
	 *
	 * Image region writing follows the following princeples:
	 *
	 *   - The file must exist before writing a region to it.
	 *
	 *   - The image header usually won't be changed except the
	 *     statistics fields and timestamp.
	 *
	 *   - The region must be inside the original image.
	 *
	 *   - If the new data are in different endian from the endian of
	 *     the original image, swap the new data.
	 *
	 *   - If the new data are of different data type from the data
     *     type of the original image, convert the new data. This may
     *     lose information when converting from longer data type to
     *     shorter data type.
	 *
	 * Each ImageIO subclass XYZ must define a static function to determine
	 * whether a given image is a valid XYZ image or not. It looks like:
	 *     static bool is_valid(const void *first_block);
	 * 'first_block' is the first block of binary data in that image file.
	 *
	 * The typical way to use an ImageIO instance is:
	 * a) To read:
	 *    ImageIO *imageio = EMUtil::get_imageio(filename, ImageIO::READ_ONLY);
	 *    int err = imageio->read_header(dict, img_index, region, is_3d);
	 *    err = imageio->read_ctf(ctf, img_index);
	 *    err = imageio->read_data(data, img_index, region, is_3d);
	 *
	 * b) To write:
	 *    similar to read.
	 */
	class ImageIO
	{
	  public:
		enum IOMode
		{ READ_ONLY = 1, READ_WRITE = 2, WRITE_ONLY = 3 };
	  public:
		virtual ~ ImageIO();

		/** Read the header from an image.
		 *
		 * @param dict A keyed-dictionary to store the header information.
		 * @param image_index The index of the image to read.
		 * @param area Define an image-region to read.
		 * @param is_3d Whether to treat the image as a single 3D or a
		 *   set of 2Ds. This is a hint for certain image formats which
		 *   has no difference between 3D image and set of 2Ds.
		 * @return 0 if OK; 1 if error.
		 */
		virtual int read_header(Dict & dict, int image_index = 0,
								const Region * area = 0, bool is_3d = false) = 0;

		/** Write a header to an image.
		 *
		 * @param dict A keyed-dictionary storing the header information.
		 * @param image_index The index of the image to write.
		 * @param area The region to write data to.
		 * @param filestoragetype The image data type used in the output file.
		 * @param use_host_endian Whether to use the host machine
		 *        endian to write out or not. If false, use the
		 *        endian opposite to the host machine's endian.
		 * @return 0 if OK; 1 if error.
		 */
		virtual int write_header(const Dict & dict,
								 int image_index = 0,
								 const Region * area = 0,
								 EMUtil::EMDataType filestoragetype = EMUtil::EM_FLOAT,
								 bool use_host_endian = true) = 0;

		/** Read the data from an image.
		 *
		 * @param data An array to store the data. It should be
		 *        created outside of this function.
		 * @param image_index The index of the image to read.
		 * @param area Define an image-region to read.
		 * @param is_3d Whether to treat the image as a single 3D or a
		 *        set of 2Ds. This is a hint for certain image formats which
		 *        has no difference between 3D image and set of 2Ds.
		 * @return 0 if OK; 1 if error.
		 */
		virtual int read_data(float *data, int image_index = 0,
							  const Region * area = 0, bool is_3d = false) = 0;

		/** Write data to an image.
		 *
		 * @param data An array storing the data.
		 * @param image_index The index of the image to write.
		 * @param area The region to write data to.
		 * @param filestoragetype The image data type used in the output file.
		 * @param use_host_endian Whether to use the host machine
		 *        endian to write out or not. If false, use the
		 *        endian opposite to the host machine's endian.
		 * @return 0 if OK; 1 if error.
		 */
		virtual int write_data(float *data,
							   int image_index = 0,
							   const Region * area = 0,
							   EMUtil::EMDataType filestoragetype = EMUtil::EM_FLOAT,
							   bool use_host_endian = true) = 0;

		/** Read CTF data from this image.
		 *
		 * @param ctf Used to store the CTF data.
		 * @param image_index The index of the image to read.
		 * @return 0 if OK; 1 if error.
		 */
		virtual int read_ctf(Ctf & ctf, int image_index = 0);

		/** Write CTF data to this image.
		 *
		 * @param ctf Ctf instance storing the CTF data.
		 * @param image_index The index of the image to write.
		 * @return 0 if OK; 1 if error.
		 */
		virtual void write_ctf(const Ctf & ctf, int image_index = 0);

		/** Flush the IO buffer.
		 */
		virtual void flush() = 0;

		/** Return the number of images in this image file. */
		virtual int get_nimg();

		/** Is this an complex image or not. */
		virtual bool is_complex_mode() = 0;

		/** Is this image in big endian or not. */
		virtual bool is_image_big_endian() = 0;

		/** Is this image format only storing 1 image or not.
		 * Some image formats like MRC only store 1 image in a file,
		 * so this function returns 'true' for them.
		 * Other image formats like IMAGIC/HDF5 may store mutliple
		 * images, so this function returns 'false' for them.
		 */
		virtual bool is_single_image_format() const
		{
			return true;
		}

		/** Convert data of this image into host endian format.
		 *
		 * @param data An array of data. It can be any type, short,
		 *        int, float, double, etc.
		 * @param n Array size.
		 */
		template < class T > void become_host_endian(T * data, size_t n = 1)
		{
			if (is_image_big_endian() != ByteOrder::is_host_big_endian()) {
				ByteOrder::swap_bytes(data, n);
			}
		}

	protected:
		/** Do some initialization beforing doing the read/write.
		 */
		virtual void init() = 0;

		/** Validate 'image_index' in file reading.
		 * @param image_index  The 'image_index'th image. Valid value
		 * = [0, nimg-1].
		 * @exception OutofRangeException If image_index is out of range.
		 */
		void check_read_access(int image_index);

		/** Validate 'image_index' and 'data' in file reading.
		 * @param image_index  The 'image_index'th image. Valid value
		 * = [0, nimg-1].
		 * @param data The array used to store the image data in reading.
		 * @exception NullPointerException If 'data' is NULL.
		 * @exception OutofRangeException If image_index is out of range.
		 */
		void check_read_access(int image_index, const float *data);

		/** Validate rw_mode and image_index in file writing.
		 * @param rw_mode File Read/Write mode.
		 * @param image_index The 'image_index'th image. Valid value =
		 * [0, max_nimg].
		 * @param max_nimg Maximum number of images in the file. If
		 * its value <= 0, don't check image_index againt max_nimg.
		 * @exception ImageWriteException Image is not opened for writing.
		 * @exception OutofRangeException If image_index is out of range.
		 */
		void check_write_access(IOMode rw_mode, int image_index, int max_nimg = 0);

		/** Validate rw_mode, image_index, and data pointer in file writing.
		 * @param rw_mode File Read/Write mode.
		 * @param image_index The 'image_index'th image. Valid value =
		 * [0, max_nimg].
		 * @param max_nimg Maximum number of images in the file. If
		 * its value <= 0, don't check image_index againt max_nimg.
		 * @param data The data array to be writting to the image file.
		 * @exception ImageWriteException Image is not opened for writing.
		 * @exception OutofRangeException If image_index is out of range.
		 */
		void check_write_access(IOMode rw_mode, int image_index, int max_nimg,
								const float *data);

		/** Validate image I/O region.
		 *
		 * @param area The image I/O region.
		 * @param max_size The upper limit of the region's size. The region
		 * must be within 'max_size'.
		 * @param is_new_file Whether it is on a new file or not.
		 * @param inbounds_only if true verifies that the region is inside the image, otherwise no check is performed
		 * @exception ImageReadException Any image reading problem.
		 */
		void check_region(const Region * area, const FloatSize & max_size,
						  bool is_new_file = false, bool inbounds_only=true);

		/** Validate image I/O region.
		 *
		 * @param area The image I/O region.
		 * @param max_size The upper limit of the region's size. The region
		 * must be within 'max_size'.
		 * @param is_new_file Whether it is on a new file or not.
		 * @param inbounds_only if true verifies that the region is inside the image, otherwise no check is performed
		 * @exception ImageReadException Any image reading problem.
		 */
		void check_region(const Region * area, const IntSize & max_size,
						  bool is_new_file = false, bool inbounds_only=true);

		/** Run fopen safely.
		 *
		 * @param filename The file name to be opened.
		 * @param mode File open mode.
		 * @param is_new Is this a new file?
		 * @param overwrite If the file exists, should we truncate it?
		 * @exception FileAccessException The file has access error.
		 * @return The opened file pointer.
		 */
		FILE *sfopen(const string & filename, IOMode mode,
					 bool * is_new = 0, bool overwrite = false);

		/** Calculate the min and max pixel value acceptedfor image nomalization,
		 * if we did not get them from image attribute dictionary, or they are not
		 * valid values
		 * rendermin = mean - 3*sigma
		 * rendermax = mean + 3*sigma
		 *
		 * @param[in] data 2D image's data array
		 * @param[in] nx x dimension size
		 * @param[in] ny y dimension size
		 * @param[out] rendermin the minmal value for normalization
		 * @param[out] rendermax the maximum value for normalization
		 * @param[in] nz z dimension size
		 * */
		void getRenderMinMax(float * data, const int nx, const int ny, float& rendermin, float& rendermax, const int nz = 1);
	};

	/** DEFINE_IMAGEIO_FUNC declares the functions that needs to
	 * be implemented by any subclass of ImageIO.
	 */
#define DEFINE_IMAGEIO_FUNC \
		int read_header(Dict & dict, int image_index = 0, const Region* area = 0, bool is_3d = false); \
		int write_header(const Dict & dict, int image_index = 0, const Region * area = 0, EMUtil::EMDataType filestoragetype = EMUtil::EM_FLOAT, bool use_host_endian = true); \
		int read_data(float* data, int image_index = 0, const Region* area = 0, bool is_3d = false); \
		int write_data(float* data, int image_index = 0, const Region * area = 0, EMUtil::EMDataType filestoragetype = EMUtil::EM_FLOAT, bool use_host_endian = true); \
		void flush(); \
		bool is_complex_mode(); \
		bool is_image_big_endian(); \
		void init()

}


#endif	//eman__imageio_h__
