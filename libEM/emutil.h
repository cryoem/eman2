/**
 * $Id$
 */
#ifndef eman__emutil__h__
#define eman__emutil__h__ 1

#include "emobject.h"
#include "imageio.h"
#include <string>
#include <vector>

using std::string;
using std::vector;

namespace EMAN
{
	class Region;

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
			IMAGE_EM
		};

		static EMData *vertical_acf(const EMData * image, int maxdy);

		static EMData *make_image_median(const vector < EMData * >&image_list);

		/** Get an image's format type from its filename extension.
		 * @params file_ext File extension.
		 * @return image format type.
		 */
		static ImageType get_image_ext_type(string file_ext);
		
		/** Get an image's format type by processing the first 1K of the image.
		 * @params filename Image file name.
		 * @return image format type.
		 */
		static ImageType get_image_type(string filename);

		/** Get the number of images in an image file.
		 * @params filename Image file name.
		 * @return Number of images in the given file.
		 */
		static int get_image_count(string filename);

		/** Get an ImageIO object. It may be a newly created
		 * object. Or an object stored in the cache.
		 * @params filename Image file name.
		 * @params rw_mode ImageIO read/write mode.
		 * @params image_type Image format type.
		 * @return An ImageIO object.
		 */
		static ImageIO *get_imageio(string filename, int rw_mode, 
									ImageType image_type = IMAGE_UNKNOWN);

		/** Give each image type a meaningful name.
		 * @params type Image format type.
		 * @return A name for that type.
		 */
		static const char *get_imagetype_name(EMUtil::ImageType type);
		
		static const char *get_datatype_string(EMDataType type);

		/** Get a region's dimensions.
		 * @params area The region area.
		 * @params nx Image x size.
		 * @params area_x The pointer used to return the region x size.
		 * @params ny Image y size.
		 * @params area_y The pointer used to return the region y size.
		 * @params nz Image z size.
		 * @params area_z The pointer used to return the region z size.
		 */
		static void get_region_dims(const Region * area, int nx, int *area_x, int ny,
									int *area_y, int nz = 1, int *area_z = 0);

		/** Get a region's original locations.
		 * @params area The region area.
		 * @params p_x0 The pointer used to return the region x origin.
		 * @params p_y0 The pointer used to return the region y origin.
		 * @params p_z0 The pointer used to return the region z origin.
		 * @params nz  Image z size.
		 * @params image_index Image index.
		 */
		static void get_region_origins(const Region * area, int *p_x0, int *p_y0, int *p_z0 = 0,
									   int nz = 1, int image_index = 0);

		/** Process image region IO. It eithers read a region from an
		 * image file. Or write a region to an image file.
		 * @params cdata Data array.
		 * @params file The image file pointer.
		 * @params rw_mode Read/write mode. It is either READ_ONLY or WRITE_ONLY.
		 * @params image_index Image index.
		 * @params mode_size Pixel size.
		 * @params nx Image x size.
		 * @params ny Image y size.
		 * @params nz Image z size.
		 * @params area The region to read/write.
		 * @params need_flip Do we need flip the image?
		 * @params pre_row File size needed to be skipped before each row.
		 * @params post_row File size needed to be skipped after each row.
		 * @exception ImageReadException If the read has some error.
		 * @exception ImageWriteException If the write has some error.
		 */
		static void process_region_io(void *cdata, FILE * file, ImageIO::IOMode rw_mode,
									  int image_index, size_t mode_size, int nx, 
									  int ny, int nz = 1, const Region * area = 0, 
									  bool need_flip = false, ImageType imgtype=IMAGE_UNKNOWN,
									  int pre_row = 0, int post_row = 0);

		/** Dump a Dict object.
		 * @params dict A Dict object.
		 */
		static void dump_dict(const Dict & dict);

		/** Check whether two EMData images are of the same size.
		 * @params image1 The first EMData image.
		 * @params image2 The second EMData image.
		 * @return Whether two EMData images are of the same size.
		 */
		static bool is_same_size(const EMData * image1, const EMData * image2);

		/** Check whether two EMData images have the same CTF.
		 * parameters.
		 * @params image1 The first EMData image.
		 * @params image2 The second EMData image.
		 * @return whether two EMData images have the same CTF.
		 */
		static bool is_same_ctf(const EMData * image1, const EMData * image2);

	  private:
		static ImageType fast_get_image_type(string filename, const void *first_block,
											 off_t file_size);
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

#endif
