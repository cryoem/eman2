/**
 * $Id$
 */
#ifndef eman__hdfio_h__
#define eman__hdfio_h__ 1

#ifdef EM_HDF5

#include "imageio.h"
#include <stdio.h>
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
		enum DataType
		{ INT, FLOAT, STRING };

	  public:
		HdfIO(string filename, IOMode rw_mode = READ_ONLY);
		~HdfIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);

		int read_ctf(Ctf & ctf, int image_index = 0);
		int write_ctf(const Ctf & ctf, int image_index = 0);

		int read_array_attr(int image_index, string attr_name, void *value);
		int write_array_attr(int image_index, string attr_name, int nitems, void *data,
							 DataType type);

		bool is_single_image_format() const
		{
			return false;
		}
		int get_nimg();
		
	  private:
		hid_t create_dataset(int image_index, int nx, int ny, int nz);

		int read_int_attr(int image_index, string attr_name);
		float read_float_attr(int image_index, string attr_name);
		string read_string_attr(int image_index, string attr_name);
		int read_global_int_attr(string attr_name);
		float read_global_float_attr(string attr_name);

		int *read_dims(int image_index, int *p_ndim);
		int read_euler_attr(int image_index, string attr_name);
		int read_mapinfo_attr(int image_index, string attr_name);

		int write_int_attr(int image_index, string attr_name, int value);
		int write_float_attr(int image_index, string attr_name, float value);
		int write_string_attr(int image_index, string attr_name, string value);

		int write_float_attr_from_dict(int image_index, string attr_name, const Dict & dict);

		int write_global_int_attr(string attr_name, int value);

		int write_euler_attr(int image_index, string attr_name, int value);
		int write_mapinfo_attr(int image_index, string attr_name, int value);

		int delete_attr(int image_index, string attr_name);

		int get_num_dataset();
		vector < int >get_image_indices();

	  private:
		enum Nametype
		{ ROOT_GROUP, CTFIT, NUMDATASET, COMPOUND_DATA_MAGIC };
		static const char *HDF5_SIGNATURE;

	  private:
		string filename;
		IOMode rw_mode;
		bool initialized;

		hid_t file;
		hid_t group;
		hid_t cur_dataset;
		int cur_image_index;

		herr_t(*old_func) (void *);
		void *old_client_data;

		static hid_t euler_type;
		static hid_t mapinfo_type;

		vector < int >image_indices;

	  private:
		void hdf_err_off();
		void hdf_err_on();

		int delete_attr(string attr_name);

		void set_dataset(int image_index);
		int create_compound_attr(int image_index, string attr_name);
		void close_dataset(hid_t dataset);
		static string get_item_name(Nametype type);
		void increase_num_dataset();

		static void create_enum_types();
		string get_compound_name(int id, string name);

		float read_float_attr(string attr_name);
		int write_int_attr(string attr_name, int value);
		int write_float_attr(string attr_name, float value);

		int get_hdf_dims(int image_index, int *p_nx, int *p_ny, int *p_nz);

		static herr_t file_info(hid_t loc_id, const char *name, void *opdata);
	};
}

#endif

#endif
