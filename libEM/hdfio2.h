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

#ifndef eman__hdfio2_h__
#define eman__hdfio2_h__ 1

#ifdef EM_HDF5

#define H5_USE_16_API

#include "imageio.h"

#include <hdf5.h>
#undef __STDINT_H
#ifndef __STDC_LIMIT_MACROS
	#define __STDC_LIMIT_MACROS 1
#endif
#ifndef __STDC_LIMIT_MACROS
	#define __STDC_CONSTANT_MACROS 1
#endif
#ifndef _WIN32
	#include <stdint.h>
#endif
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
	 * 
	 * After you make change to this class, please check the HDF5 file created 
	 * by EMAN2 with the h5check program from:
	 * ftp:://ftp.hdfgroup.org/HDF5/special_tools/h5check/
	 * to verify the HDF5 file is compliant with the HDF5 File Format Specification.
	 */
	class HdfIO2:public ImageIO
	{
	  public:
		explicit HdfIO2(const string & filename, IOMode rw_mode = READ_ONLY);
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

		/** Return the file id
		 * For single attribute read/write*/
		hid_t get_fileid() const {return file;}

	  private:
		hsize_t nx, ny, nz;
		bool is_exist;	//boolean to tell if the image (group) already exist(to be overwrite)

		hid_t file;
		hid_t group;
		hid_t accprop;
		hid_t simple_space;
		string filename;
		IOMode rw_mode;
		bool initialized;
		
		Dict meta_attr_dict;	//this is used for the meta attributes stored in /MDF/images

		/* Erases any existing attributes from the image group
		 * prior to writing a new header. For a new image there
		 * won't be any, so this should be harmless. 
		 * 
		 * @param image_index 
		 * @return 0 for success*/
		int erase_header(int image_index);

        // render_min and render_max
	    float rendermin;
	    float rendermax;
	};
}

#endif

#endif
