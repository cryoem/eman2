/**
 * $Id$
 */
#ifndef eman__imageio_h__
#define eman__imageio_h__ 1

#include "emobject.h"
#include "byteorder.h"

#include <vector>
#include <string>
#include <stdio.h>

using std::vector;
using std::string;

namespace EMAN
{
	class Region;
	class Size;
	class Ctf;

	/** ImageIO is the  base class for any image io class.
     * A subclass should implement functions declared in
     * DEFINE_IMAGEIO_FUNC macro.
     *
     * Image header I/O is separated from image data I/O.
     * valid image_index = [0, n].
	 *
	 * The typical way to use an ImageIO instance is:
     */
	class ImageIO
	{
	  public:
		enum IOMode
		{ READ_ONLY = 1, READ_WRITE = 2 };
	  public:
		virtual ~ ImageIO();

		/** read the header from an image.
		 */
		virtual int read_header(Dict & dict, int image_index = 0,
								const Region * area = 0, bool is_3d = false) = 0;

		virtual int write_header(const Dict & dict, int image_index = 0,
								 bool use_host_endian = true) = 0;

		virtual int read_data(float *data, int image_index = 0,
							  const Region * area = 0, bool is_3d = false) = 0;

		virtual int write_data(float *data, int image_index = 0, bool use_host_endian = true) = 0;

		virtual int read_ctf(Ctf & ctf, int image_index = 0);
		virtual int write_ctf(const Ctf & ctf, int image_index = 0);
		virtual int get_nimg() = 0;

		virtual bool is_complex_mode() = 0;
		virtual bool is_image_big_endian() = 0;

		template < class T > void become_host_endian(T * data, int n = 1) {
			if (is_image_big_endian() != ByteOrder::is_host_big_endian()) {
				ByteOrder::swap_bytes(data, n);
			}
		}

	  protected:
		virtual int init() = 0;
		int check_read_access(int image_index, bool check_data = false, float *data = 0);
		int check_write_access(IOMode rw_mode, int image_index, bool check_data = false,
							   float *data = 0);

		int check_region(const Region * area, const Size & max_size);

		FILE *sfopen(string filename, IOMode mode, bool * is_new = 0, bool overwrite = false);
	};

	/** macro DEFINE_IMAGEIO_FUNC declares the functions that needs to
     * be implemented by subclass of ImageIO.
     */
#define DEFINE_IMAGEIO_FUNC \
        int read_header(Dict & dict, int image_index = 0, const Region* area = 0, bool is_3d = false); \
	int write_header(const Dict & dict, int image_index = 0, bool use_host_endian = true); \
	int read_data(float* data, int image_index = 0, const Region* area = 0, bool is_3d = false); \
	int write_data(float* data, int image_index = 0, bool use_host_endian = true); \
        bool is_complex_mode(); \
        bool is_image_big_endian(); \
        int get_nimg(); \
        int init()

}


#endif
