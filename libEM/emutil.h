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

#ifndef eman__emutil__h__
#define eman__emutil__h__ 1

#include "emobject.h"
#include "emassert.h"
#include <string.h>

using std::string;
using std::vector;

// Defining EMDelete using templates
// Use EMDelete instead of delete as this will be very clean.
// This should perhaps be moved to somewhere else later to make it
// more widely accessible. (C. Yang 04/27/06)
template <class T>
inline void EMDeletePtr(T & x)
{
#ifdef _WIN32
	if(x != NULL) {
		delete x;
		x = NULL;
	}
#else
	{Assert(x != NULL);}
	delete x;
	x = NULL;
#endif
}

template <class T>
inline void EMDeleteArray(T & x)
{
#ifdef _WIN32
	if(x != NULL) {
		delete x;
		x = NULL;
	}
#else
	{Assert(x != NULL);}
     delete [] x;
     x = NULL;
#endif
}

namespace EMAN
{
	class Region;
	class ImageIO;

	class EMUtil
	{
	public:
		/** Image pixel data type used in EMAN.
		 * EM_U means "EM unsigned". for example, EM_USHORT means
		 * EM unsigned short.
		 */
		enum EMDataType
		{
			EM_UNKNOWN,
			EM_CHAR,
			EM_UCHAR,
			EM_SHORT,
			EM_USHORT,
			EM_INT,
			EM_UINT,
			EM_FLOAT,
			EM_DOUBLE,
			EM_SHORT_COMPLEX,
			EM_USHORT_COMPLEX,
			EM_FLOAT_COMPLEX
		};

		/** Image format types.
		 */
		enum ImageType
		{
			IMAGE_UNKNOWN,
			IMAGE_MRC,
			IMAGE_SPIDER,
			IMAGE_SINGLE_SPIDER,
			IMAGE_IMAGIC,
			IMAGE_HDF,
			IMAGE_DM3,
			IMAGE_TIFF,
			IMAGE_PGM,
			IMAGE_LST,
			IMAGE_PIF,
			IMAGE_VTK,
			IMAGE_PNG,
			IMAGE_SAL,
			IMAGE_ICOS,
			IMAGE_EMIM,
			IMAGE_GATAN2,
			IMAGE_AMIRA,
			IMAGE_XPLOR,
			IMAGE_EM,
			IMAGE_V4L,
			IMAGE_JPEG,
			IMAGE_FITS,
			IMAGE_LSTFAST
		};

		static EMData *vertical_acf(const EMData * image, int maxdy);

		static EMData *make_image_median(const vector < EMData * >&image_list);

		/** Get an image's format type from its filename extension.
		 * @param file_ext File extension.
		 * @return image format type.
		 */
		static ImageType get_image_ext_type(const string & file_ext);

		/** Get an image's format type by processing the first 1K of the image.
		 * @param filename Image file name.
		 * @return image format type.
		 */
		static ImageType get_image_type(const string & filename);

		/** Ask whether or not the given filename is a valid EM image filename
		 * This is the same thing as checking whether or not the return value of EMUtil.get_image_ext_type
		 * is IMAGE_UNKNOWN
		 * @param filename Image file name.
		 * @return whether or not it is a valid filename
		 */
		static bool is_valid_filename(const string & filename);

		/** Get the number of images in an image file.
		 * @param filename Image file name.
		 * @return Number of images in the given file.
		 */
		static int get_image_count(const string & filename);

		/** Get an ImageIO object. It may be a newly created
		 * object. Or an object stored in the cache.
		 * @param filename Image file name.
		 * @param rw_mode ImageIO read/write mode.
		 * @param image_type Image format type.
		 * @return An ImageIO object.
		 */
		static ImageIO *get_imageio(const string & filename, int rw_mode,
									ImageType image_type = IMAGE_UNKNOWN);

		/** Give each image type a meaningful name.
		 * @param type Image format type.
		 * @return A name for that type.
		 */
		static const char *get_imagetype_name(EMUtil::ImageType type);

		/** Give each data type a meaningful name
		 * @param type the EMDataType
		 * @return a name for that data type
		 * */
		static const char *get_datatype_string(EMDataType type);

		/** Get a region's dimensions.
		 * @param area The region area.
		 * @param nx Image x size.
		 * @param area_x The pointer used to return the region x size.
		 * @param ny Image y size.
		 * @param area_y The pointer used to return the region y size.
		 * @param nz Image z size.
		 * @param area_z The pointer used to return the region z size.
		 */
		static void get_region_dims(const Region * area, int nx, int *area_x, int ny,
									int *area_y, int nz = 1, int *area_z = 0);

		/** Get a region's original locations.
		 * @param area The region area.
		 * @param p_x0 The pointer used to return the region x origin.
		 * @param p_y0 The pointer used to return the region y origin.
		 * @param p_z0 The pointer used to return the region z origin.
		 * @param nz  Image z size.
		 * @param image_index Image index.
		 */
		static void get_region_origins(const Region * area, int *p_x0, int *p_y0,
									   int *p_z0 = 0, int nz = 1, int image_index = 0);

		/** Process image region IO. It eithers read a region from an
		 * image file. Or write a region to an image file.
		 * Works for regions that are outside the image data dimension area.(David Woolford, April 23 2009)
		 * @param cdata Data array.
		 * @param file The image file pointer.
		 * @param rw_mode Read/write mode. It is either READ_ONLY or WRITE_ONLY.
		 * @param image_index Image index.
		 * @param mode_size Pixel size.
		 * @param nx Image x size.
		 * @param ny Image y size.
		 * @param nz Image z size.
		 * @param area The region to read/write.
		 * @param need_flip Do we need flip the image?
		 * @param imgtype The Image type of the processed file.
		 * @param pre_row File size needed to be skipped before each row.
		 * @param post_row File size needed to be skipped after each row.
		 * @exception ImageReadException If the read has some error.
		 * @exception ImageWriteException If the write has some error.
		 */
		static void process_region_io(void *cdata, FILE * file, int rw_mode,
									  int image_index, size_t mode_size, int nx,
									  int ny, int nz = 1, const Region * area = 0,
									  bool need_flip = false, ImageType imgtype=IMAGE_UNKNOWN,
									  int pre_row = 0, int post_row = 0);


		/**
		 * Works for regions that are outside the image data dimension area.
		 * The only function that calls this is in xplorio.cpp - that function
		 * throws if the region is invalid.
		 */
		static void process_ascii_region_io(float *data, FILE * file, int rw_mode,
											int image_index, size_t mode_size,
											int nx, int ny, int nz,
											const Region * area, bool has_index_line,
											int nitems_per_line, const char *outformat);


		/** Dump a Dict object.
		 * @param dict A Dict object.
		 */
		static void dump_dict(const Dict & dict);

		/** Check whether two EMData images are of the same size.
		 * @param image1 The first EMData image.
		 * @param image2 The second EMData image.
		 * @return Whether two EMData images are of the same size.
		 */
		static bool is_same_size(const EMData * image1, const EMData * image2);

		/** Check whether two EMData images have the same CTF parameters.
		 * @param image1 The first EMData image.
		 * @param image2 The second EMData image.
		 * @return whether two EMData images have the same CTF.
		 */
		static bool is_same_ctf(const EMData * image1, const EMData * image2);


		static bool is_complex_type(EMDataType datatype);

		static void jump_lines(FILE * file, int nlines);

        static vector<string> get_euler_names(const string & euler_type);

		/** Get an attribute from a stack of image, returned as a vector
		 *
		 * @param file_name the image file name
		 * @param attr_name The header attribute name.
		 * @return the vector of attribute value
		 * @exception NotExistingObjectException when access an non-existing attribute
		 * @exception InvalidCallException when call this function for a non-stack image */
		static vector<EMObject> get_all_attributes(const string & file_name, const string & attr_name);

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
		static void getRenderMinMax(float * data, const int nx, const int ny, float& rendermin, float& rendermax, const int nz = 1);
		

		static bool cuda_available() {
#ifdef EMAN2_USING_CUDA
			return true;
#else
			return false;
#endif
		}

		inline static void* em_malloc(const size_t size) {
			return malloc(size);
		}

		inline static void* em_calloc(const size_t nmemb,const size_t size) {
			return calloc(nmemb,size);
		}

		inline static void* em_realloc(void* data,const size_t new_size) {
			return realloc(data, new_size);
		}
		inline static void em_memset(void* data, const int value, const size_t size) {
			memset(data, value, size);
		}
		inline static void em_free(void*data) {
			free(data);
		}

		inline static void em_memcpy(void* dst,const void* const src,const size_t size) {
			memcpy(dst,src,size);
		}
	  private:
		static ImageType fast_get_image_type(const string & filename,
											 const void *first_block,
											 off_t file_size);

		static void jump_lines_by_items(FILE * file, int nitems, int nitems_per_line);



		static void process_numbers_io(FILE * file, int rw_mode,
									   int nitems_per_line,
									   size_t mode_size, int start, int end,
									   float *data, int *p_i,
									   const char *outformat);

		static void exclude_numbers_io(FILE * file, int rw_mode,
									   int nitems_per_line,
									   size_t mode_size, int start, int end,
									   float * data, int *p_i,
									   const char *outformat);

		static void process_lines_io(FILE * file, int rw_mode,
									 int nitems_per_line, size_t mode_size,
									 int nitems, float *data, int *p_i,
									 const char *outformat);


	};

	struct ImageScore {
		int index;
		float score;
		ImageScore(int i=0, float s=0) : index(i), score(s) {}
	};

	class ImageSort {
	public:
		ImageSort(int n);
		~ImageSort();

		void sort();

		void set(int i, float score);
		int get_index(int i) const;
		float get_score(int i) const;

		int size() const;
	private:
		ImageScore* image_scores;
		int n;

	};
}

#endif	//eman__emutil__h__
