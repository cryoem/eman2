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
	 * 
	 * Attribute name must be within 128 charaters, including string terminator '\0'. 
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
        
        /* If this version of init() returns -1, then we have an 
         * old-style HDF5 file 
         * 
         * @return 0 for new HDF5 file, -1 for old-style HDF5 file*/
		int init_test();
			
		/* Return the number of images in this HDF file
		 * 
		 * @return the number of images in this file*/
		int get_nimg();

		int get_num_dataset();
		vector < int >get_image_indices();

	  private:
		hid_t file;
		hid_t group;
		hid_t accprop;
		hid_t simple_space;
		string filename;
		IOMode rw_mode;
		bool initialized;
		
		/* Erases any existing attributes from the image group
		 * prior to writing a new header. For a new image there
		 * won't be any, so this should be harmless. 
		 * 
		 * @param image_index 
		 * @return 0 for success*/
		int erase_header(int image_index);
		
		/* Write an attribute with specified name to a given open object.
		 * The attribute is opened and closed. returns 0 on success 
		 * 
		 * @param loc Object (dataset, group, or named datatype) id for this attibute to be attached to
		 * @param name name for this attibute
		 * @param obj value for this attribute
		 * @return 0 on success */
		int write_attr(hid_t loc, const char *name, EMObject obj);
		
		/* Read an already opened attribute and returns the results 
		 * as an EMObject. The attribute is not closed. 
		 * 
		 * @param attr Identifier of an attribute to read
		 * @return attribute value as an EMObject */
		EMObject read_attr(hid_t attr);
	};
}

#endif

#endif
