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

#ifndef eman__hdfio_h__
#define eman__hdfio_h__ 1

#ifdef EM_HDF5

#define H5_USE_16_API

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
	 * 
	 * After you make change to this class, please check the HDF5 file created 
	 * by EMAN2 with the h5check program from:
	 * ftp:://ftp.hdfgroup.org/HDF5/special_tools/h5check/
	 * to verify the HDF5 file is compliant with the HDF5 File Format Specification. 
	 */
	class HdfIO:public ImageIO
	{
	  public:
		enum DataType
		{ INT, FLOAT, STRING };

	  public:
		explicit HdfIO(const string & filename, IOMode rw_mode = READ_ONLY);
		~HdfIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);

		int read_ctf(Ctf & ctf, int image_index = 0);
		void write_ctf(const Ctf & ctf, int image_index = 0);

        int read_euler_angles(Dict & euler_angles, int image_index = 0);
        void write_euler_angles(const Dict & euler_angles, int image_index = 0);
        
		int read_array_attr(int image_index, const string & attr_name, void *value);
		int write_array_attr(int image_index, const string & attr_name,
							 int nitems, void *data, DataType type);

		bool is_single_image_format() const
		{
			return false;
		}
        
		int get_nimg();

        
		int read_int_attr(int image_index, const string & attr_name);
		float read_float_attr(int image_index, const string & attr_name);
		string read_string_attr(int image_index, const string & attr_name);
		int read_global_int_attr(const string & attr_name);
		float read_global_float_attr(const string & attr_name);

		int read_mapinfo_attr(int image_index, const string & attr_name);

		int write_int_attr(int image_index, const string & attr_name, int value);
		int write_float_attr(int image_index, const string & attr_name, float value);
		int write_string_attr(int image_index, const string & attr_name,
							  const string & value);

		int write_float_attr_from_dict(int image_index, const string & attr_name,
									   const Dict & dict);

		int write_global_int_attr(const string & attr_name, int value);

		int write_mapinfo_attr(int image_index, const string & attr_name, int value);

		int delete_attr(int image_index, const string & attr_name);

		int get_num_dataset();
		vector < int >get_image_indices();
        
	  private:
		enum Nametype
		{ ROOT_GROUP, CTFIT, NUMDATASET, COMPOUND_DATA_MAGIC, EULER };
        
		void create_cur_dataset(int image_index, int nx, int ny, int nz);

		int *read_dims(int image_index, int *p_ndim);

        int read_compound_dict(Nametype compound_type,
                               Dict & values, int image_index);
        
        void write_compound_dict(Nametype compound_type,
                                 const Dict & values, int image_index);
        

		static const char *HDF5_SIGNATURE;
        

		string filename;
		IOMode rw_mode;
		bool initialized;
		bool is_new_file;
		
		hid_t file;
		hid_t group;
		hid_t cur_dataset;
		int cur_image_index;

		herr_t(*old_func) (void *);
		void *old_client_data;

		static hid_t mapinfo_type;

		vector < int >image_indices;

		void hdf_err_off();
		void hdf_err_on();

		int delete_attr(const string & attr_name);

		void set_dataset(int image_index);
		int create_compound_attr(int image_index, const string & attr_name);
		void close_cur_dataset();
		static string get_item_name(Nametype type);
		void increase_num_dataset();

		static void create_enum_types();
		string get_compound_name(int id, const string & name);

		float read_float_attr(const string & attr_name);
		int write_int_attr(const string & attr_name, int value);
		int write_float_attr(const string & attr_name, float value);

		int get_hdf_dims(int image_index, int *p_nx, int *p_ny, int *p_nz);

		static herr_t file_info(hid_t loc_id, const char *name, void *opdata);
		int create_region_space(hid_t * p_dataspace_id, hid_t * p_memspace_id,
								const Region * area, int nx, int ny, int nz,
								int image_index);

       
	};
}

#endif

#endif
