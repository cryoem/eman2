/**
 * $Id$
 */
#ifndef eman__hdfio2_h__
#define eman__hdfio2_h__ 1

#ifdef EM_HDF5

#include "imageio.h"

#include <hdf5.h>
#include <vector>

using std::vector;

namespace EMAN
{
	/** HDF5 (hiearchical data format version 5) is supported in
	 * HdfIO. This is a revised HDF5 format file.
	 *
	 * A HDF5 file may contains multiple 2D or 3D images.
	 */
	class HdfIO2:public ImageIO
	{
	  public:
		HdfIO2(const string & filename, IOMode rw_mode = READ_ONLY);
		~HdfIO2();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);

		bool is_single_image_format() const
		{
			return false;
		}
        
		int erase_header(int image_index);
		int write_attr(hid_t loc,const char *name,EMObject obj);
		EMObject read_attr(hid_t attr);
		int get_nimg();

		int get_num_dataset();
		vector < int >get_image_indices();
		hid_t file;
		hid_t group;
		hid_t simple_space;
		string filename;
		IOMode rw_mode;
		bool initialized;
	};
}

#endif

#endif
