/**
 * $Id$
 */
#ifndef eman__hdfio_h__
#define eman__hdfio_h__ 1

#ifdef EM_HDF5

#include "imageio.h"

#include <hdf5.h>
#include <vector>

using std::vector;

namespace EMAN
{
	/** HDF5 (hiearchical data format version 5) is supported in
	 * HdfIO.
	 *
	 * A HDF5 file may contains multiple 2D or 3D images.
	 */
	class HdfIO:public ImageIO
	{
	  public:
		HdfIO(const string & filename, IOMode rw_mode = READ_ONLY);
		~HdfIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);

		bool is_single_image_format() const
		{
			return false;
		}
        
		int get_nimg();

		int get_num_dataset();
		vector < int >get_image_indices();
        
		string filename;
		IOMode rw_mode;
		bool initialized;
		bool is_new_file;
		
       
	};
}

#endif

#endif
