/**
 * $Id$
 */
#ifdef EM_HDF5

#include "hdfio.h"
#include "log.h"
#include "emobject.h"
#include "util.h"
#include "emutil.h"
#include "geometry.h"
#include "ctf.h"

#include <string.h>
#include <stdlib.h>
#include <sys/param.h>
#include <assert.h>

using namespace EMAN;

hid_t HdfIO::euler_type = -1;
hid_t HdfIO::mapinfo_type = -1;
const char *HdfIO::HDF5_SIGNATURE = "\211HDF\r\n\032\n";


HdfIO::HdfIO(string hdf_filename, IOMode rw)
:	filename(hdf_filename), rw_mode(rw)
{
	initialized = false;
	file = -1;
	group = -1;
	cur_dataset = -1;
	cur_image_index = -1;
	old_func = 0;
	old_client_data = 0;
}

HdfIO::~HdfIO()
{
	close_dataset(cur_dataset);
	H5Gclose(group);
	H5Fclose(file);
}

int HdfIO::init()
{
	static int err = 0;
	if (initialized) {
		return err;
	}
	LOGDEBUG("HdfIO::init()");
	initialized = true;

	bool is_new_file = false;
	FILE *tmp_file = sfopen(filename, rw_mode, &is_new_file);
	if (!tmp_file) {
		err = 1;
		return err;
	}

	if (!is_new_file) {
		char buf[128];
		if (fread(buf, sizeof(buf), 1, tmp_file) != 1) {
			LOGERR("cannot read file '%s'", filename.c_str());
			err = 1;
		}
		else {
			if (!is_valid(&buf)) {
				LOGERR("'%s' is not a valid HDF5 file", filename.c_str());
				err = 1;
			}
		}
	}

	fclose(tmp_file);
	tmp_file = 0;

	if (err) {
		return err;
	}

	if (rw_mode == READ_ONLY) {
		file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	}
	else if (rw_mode == READ_WRITE) {
		hdf_err_off();
		file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
		hdf_err_on();
		if (file < 0) {
			file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
			H5Fclose(file);
			file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
		}
	}

	string root_group_str = get_item_name(ROOT_GROUP);
	group = H5Gopen(file, root_group_str.c_str());
	cur_dataset = -1;

	H5Giterate(file, root_group_str.c_str(), NULL, file_info, &image_indices);
	create_enum_types();

	return err;
}


bool HdfIO::is_valid(const void *first_block)
{
	LOGDEBUG("HdfIO::is_valid()");
	bool valid = false;

	if (first_block) {
		const char *data = static_cast < const char *>(first_block);
		size_t hdf5_sig_len = strlen(HDF5_SIGNATURE);
		char *signature = new char[hdf5_sig_len + 1];
		strncpy(signature, data, hdf5_sig_len);
		signature[hdf5_sig_len] = '\0';

		if (strcmp(signature, HDF5_SIGNATURE) == 0) {
			valid = true;
		}
		delete[]signature;
		signature = 0;
	}
	return valid;
}

int HdfIO::read_header(Dict & dict, int image_index, const Region * area, bool)
{
	LOGDEBUG("HdfIO::read_header()");
	if (check_read_access(image_index) != 0) {
		return 1;
	}
	int nx = 0, ny = 0, nz = 0;
	if (get_hdf_dims(image_index, &nx, &ny, &nz) != 0) {
		return 1;
	}
	if (check_region(area, Size(nx, ny, nz)) != 0) {
		return 1;
	}

	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen, nz, &zlen);

	dict["nx"] = xlen;
	dict["ny"] = ylen;
	dict["nz"] = zlen;

	dict["apix_x"] = read_float_attr(image_index, "apix_x");
	dict["apix_y"] = read_float_attr(image_index, "apix_y");
	dict["apix_z"] = read_float_attr(image_index, "apix_z");

	dict["origin_row"] = read_float_attr(image_index, "origin_row");
	dict["origin_col"] = read_float_attr(image_index, "origin_col");
	dict["origin_sec"] = read_float_attr(image_index, "origin_sec");

	dict["minimum"] = read_float_attr(image_index, "minimum");
	dict["maximum"] = read_float_attr(image_index, "maximum");
	dict["mean"] = read_float_attr(image_index, "mean");
	dict["sigma"] = read_float_attr(image_index, "sigma");

	dict["micrograph_id"] = read_string_attr(image_index, "micrograph_id");

	dict["particle_center_x"] = read_float_attr(image_index, "particle_center_x");
	dict["particle_center_y"] = read_float_attr(image_index, "particle_center_y");

	dict["center_x"] = read_float_attr(image_index, "center_x");
	dict["center_y"] = read_float_attr(image_index, "center_y");

	dict["score"] = read_float_attr(image_index, "score");
	dict["good"] = read_int_attr(image_index, "good");
	dict["orientation_convention"] = (int) read_euler_attr(image_index, "orientation_convention");

	dict["datatype"] = EMUtil::EM_FLOAT;

	return 0;
}


int HdfIO::read_data(float *data, int image_index, const Region * area, bool)
{
	LOGDEBUG("HdfIO::read_data() from file '%s'", filename.c_str());
	if (check_read_access(image_index, true, data) != 0) {
		return 1;
	}

	set_dataset(image_index);

	if (cur_dataset < 0) {
		char cur_dataset_name[32];
		sprintf(cur_dataset_name, "%d", image_index);
		cur_dataset = H5Dopen(file, cur_dataset_name);

		if (cur_dataset < 0) {
			LOGERR("hdf file has no image with id = %d", image_index);
			return 1;
		}
	}

	hid_t datatype = H5Dget_type(cur_dataset);
	H5T_class_t t_class = H5Tget_class(datatype);

	if (t_class != H5T_FLOAT) {
		LOGERR("unknown data type '%d'. Can only read FLOAT data in HDF",
							 (int) t_class);
		H5Tclose(datatype);
		return 1;
	}

	int nx = 0, ny = 0, nz = 0;
	if (get_hdf_dims(image_index, &nx, &ny, &nz) != 0) {
		return 1;
	}
	if (check_region(area, Size(nx, ny, nz)) != 0) {
		return 1;
	}

	int err = 0;
	if (!area) {
		err = H5Dread(cur_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	}
	else {
		hssize_t offset[3];
		hsize_t count[3];

		int x0 = 0, y0 = 0, z0 = 0;
		int xlen = 0, ylen = 0, zlen = 0;

		EMUtil::get_region_origins(area, &x0, &y0, &z0, nz, image_index);
		EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen, nz, &zlen);

		offset[0] = static_cast < hssize_t > (x0);
		offset[1] = static_cast < hssize_t > (y0);
		offset[2] = static_cast < hssize_t > (z0);

		count[0] = static_cast < hsize_t > (xlen);
		count[1] = static_cast < hsize_t > (ylen);
		count[2] = static_cast < hsize_t > (zlen);

		hid_t dataspace = H5Dget_space(cur_dataset);

		err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
		if (err < 0) {
			LOGERR("failed to create dataspace hyperslab when reading '%s'",
								 filename.c_str());
		}
		else {
			hid_t memspace = H5Screate_simple(3, count, NULL);
			hssize_t offset_out[3];
			offset_out[0] = 0;
			offset_out[1] = 0;
			offset_out[2] = 0;

			err = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count, NULL);
			if (err < 0) {
				LOGERR("failed to create memspace hyperslab when reading '%s'",
									 filename.c_str());
			}
			else {
				err =
					H5Dread(cur_dataset, H5T_NATIVE_FLOAT, memspace, dataspace, H5P_DEFAULT, data);
			}

			H5Sclose(dataspace);
			H5Sclose(memspace);
		}
	}

	if (err < 0) {
		LOGERR("reading %dth hdf image failed", image_index);
		return 1;
	}

	H5Tclose(datatype);
	return 0;
}


int HdfIO::write_header(const Dict & dict, int image_index, bool)
{
	LOGDEBUG("HdfIO::write_header()");
	if (check_write_access(rw_mode, image_index) != 0) {
		return 1;
	}

	int nx = dict["nx"];
	int ny = dict["ny"];
	int nz = dict["nz"];

	cur_dataset = create_dataset(image_index, nx, ny, nz);
	if (cur_dataset < 0) {
		LOGERR("create dataset failed in hdf write header on file '%s'",
							 filename.c_str());
		return 1;
	}

	write_int_attr(image_index, "nx", nx);
	write_int_attr(image_index, "ny", ny);
	write_int_attr(image_index, "nz", nz);

	write_float_attr_from_dict(image_index, "apix_x", dict);
	write_float_attr_from_dict(image_index, "apix_y", dict);
	write_float_attr_from_dict(image_index, "apix_z", dict);

	write_float_attr_from_dict(image_index, "origin_row", dict);
	write_float_attr_from_dict(image_index, "origin_col", dict);
	write_float_attr_from_dict(image_index, "origin_sec", dict);

	write_float_attr_from_dict(image_index, "minimum", dict);
	write_float_attr_from_dict(image_index, "maximum", dict);
	write_float_attr_from_dict(image_index, "mean", dict);
	write_float_attr_from_dict(image_index, "sigma", dict);

	if (dict.has_key("micrograph_id")) {
		write_string_attr(image_index, "micrograph_id", (const char *) dict["micrograph_id"]);
	}

	write_float_attr_from_dict(image_index, "particle_center_x", dict);
	write_float_attr_from_dict(image_index, "particle_center_y", dict);

	write_float_attr_from_dict(image_index, "center_x", dict);
	write_float_attr_from_dict(image_index, "center_y", dict);
	write_float_attr_from_dict(image_index, "score", dict);

	if (dict.has_key("good")) {
		write_int_attr(image_index, "good", dict["good"]);
	}

	if (dict.has_key("orientation_convention")) {
		write_euler_attr(image_index, "orientation_convention", dict["orientation_convention"]);
	}

	return 0;
}

int HdfIO::write_data(float *data, int image_index, bool)
{
	LOGDEBUG("HdfIO::write_data()");
	if (check_write_access(rw_mode, image_index, true, data) != 0) {
		return 1;
	}

	int nx = read_int_attr(image_index, "nx");
	int ny = read_int_attr(image_index, "ny");
	int nz = read_int_attr(image_index, "nz");

	cur_dataset = create_dataset(image_index, nx, ny, nz);

	int err = 1;
	if (cur_dataset >= 0) {
		err = H5Dwrite(cur_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
		if (err >= 0) {
			increase_num_dataset();
			image_indices.push_back(image_index);
			err = 0;
		}
	}

	return err;
}


int *HdfIO::read_dims(int image_index, int *p_ndim)
{
	set_dataset(image_index);

	if (cur_dataset < 0) {
		char cur_dataset_name[32];
		sprintf(cur_dataset_name, "%d", image_index);
		cur_dataset = H5Dopen(file, cur_dataset_name);

		if (cur_dataset < 0) {
			fprintf(stderr, "Error in reading data dimensions in hdf: %d", image_index);
			return 0;
		}
	}

	hid_t dataspace = H5Dget_space(cur_dataset);
	int rank = H5Sget_simple_extent_ndims(dataspace);
	hsize_t *dims = new hsize_t[rank];
	H5Sget_simple_extent_dims(dataspace, dims, NULL);

	int *dims1 = new int[rank];
	for (int i = 0; i < rank; i++) {
		dims1[i] = static_cast < int >(dims[i]);
	}

	H5Sclose(dataspace);
	delete[]dims;
	dims = 0;

	(*p_ndim) = rank;
	return dims1;
}

int HdfIO::read_global_int_attr(string attr_name)
{
	int value = 0;
	hid_t attr = H5Aopen_name(group, attr_name.c_str());
	if (attr < 0) {
		//LOGWARN("no such hdf attribute '%s'", attr_name.c_str());
	}
	else {
		H5Aread(attr, H5T_NATIVE_INT, &value);
		H5Aclose(attr);
	}

	return value;
}

float HdfIO::read_global_float_attr(string attr_name)
{
	float value = 0;
	hid_t attr = H5Aopen_name(group, attr_name.c_str());
	if (attr >= 0) {
		H5Aread(attr, H5T_NATIVE_INT, &value);
		H5Aclose(attr);
	}
	return value;
}

int HdfIO::get_nimg()
{
	if (init() != 0) {
		return 0;
	}

	hdf_err_off();
	int n = read_global_int_attr(get_item_name(NUMDATASET));
	hdf_err_on();
	return n;
}


void HdfIO::increase_num_dataset()
{
	int n = get_nimg();
	n++;
	write_global_int_attr(get_item_name(NUMDATASET), n);
}


void HdfIO::set_dataset(int image_index)
{
	int need_update = 0;

	if (cur_image_index >= 0) {
		if (image_index != cur_image_index) {
			need_update = 1;
		}
	}
	else {
		cur_image_index = image_index;
		need_update = 1;
	}

	if (need_update) {
		char cur_dataset_name[32];
		sprintf(cur_dataset_name, "%d", image_index);
		hdf_err_off();
		close_dataset(cur_dataset);
		cur_dataset = H5Dopen(file, cur_dataset_name);
		hdf_err_on();
	}
}



int HdfIO::read_int_attr(int image_index, string attr_name)
{
	set_dataset(image_index);

	int value = 0;
	hid_t attr = H5Aopen_name(cur_dataset, attr_name.c_str());
	if (attr >= 0) {
		H5Aread(attr, H5T_NATIVE_INT, &value);
		H5Aclose(attr);
	}
	else {
		//LOGWARN("no such hdf attribute '%s'", attr_name.c_str());
	}

	return value;
}


float HdfIO::read_float_attr(int image_index, string attr_name)
{
	set_dataset(image_index);
	return read_float_attr(attr_name);
}


float HdfIO::read_float_attr(string attr_name)
{
	float value = 0;
	hid_t attr = H5Aopen_name(cur_dataset, attr_name.c_str());
	if (attr >= 0) {
		H5Aread(attr, H5T_NATIVE_FLOAT, &value);
		H5Aclose(attr);
	}
	else {
		//LOGWARN("no such hdf attribute '%s'", attr_name.c_str());
	}

	return value;
}



string HdfIO::read_string_attr(int image_index, string attr_name)
{
	set_dataset(image_index);
	hid_t attr = H5Aopen_name(cur_dataset, attr_name.c_str());

	string value = "";
	if (attr < 0) {
		//LOGWARN("no such hdf attribute '%s'", attr_name.c_str());
	}
	else {
		char *tmp_value = new char[MAXPATHLEN];
		hid_t datatype = H5Tcopy(H5T_C_S1);
		H5Tset_size(datatype, MAXPATHLEN);
		H5Aread(attr, datatype, tmp_value);

		H5Tclose(datatype);
		H5Aclose(attr);

		value = tmp_value;

		delete[]tmp_value;
		tmp_value = 0;
	}

	return value;
}

int HdfIO::read_array_attr(int image_index, string attr_name, void *value)
{
	set_dataset(image_index);
	int err = 0;

	hid_t attr = H5Aopen_name(cur_dataset, attr_name.c_str());
	if (attr < 0) {
		LOGERR("no such hdf attribute '%s'", attr_name.c_str());
		err = 1;
	}
	else {
		hid_t datatype = H5Aget_type(attr);
		H5Aread(attr, datatype, value);
		H5Tclose(datatype);
		H5Aclose(attr);
	}
	return err;
}

int HdfIO::read_euler_attr(int image_index, string)
{
	set_dataset(image_index);
	/*
	   Euler::EulerType val;
	   hid_t attr = H5Aopen_name(cur_dataset, attr_name.c_str());
	   if (attr < 0) {
	   LOGERR(no such hdf attribute '%s'", attr_name.c_str());
	   return attr;
	   }
	   H5Aread(attr, euler_type, &val);
	   H5Aclose(attr);

	   return static_cast<int>(val);
	 */
	return 1;
}


int HdfIO::read_mapinfo_attr(int image_index, string attr_name)
{
	set_dataset(image_index);

	MapInfoType val = ICOS_UNKNOWN;
	hid_t attr = H5Aopen_name(cur_dataset, attr_name.c_str());
	if (attr < 0) {
		LOGERR("no such hdf attribute '%s'", attr_name.c_str());
	}
	else {
		H5Aread(attr, mapinfo_type, &val);
		H5Aclose(attr);
	}
	return static_cast < int >(val);
}

int HdfIO::write_int_attr(int image_index, string attr_name, int value)
{
	set_dataset(image_index);
	return write_int_attr(attr_name, value);
}

int HdfIO::write_int_attr(string attr_name, int value)
{
	int err = -1;
	delete_attr(attr_name);

	hid_t dataspace = H5Screate(H5S_SCALAR);
	hid_t attr = H5Acreate(cur_dataset, attr_name.c_str(), H5T_NATIVE_INT, dataspace, H5P_DEFAULT);

	if (attr >= 0) {
		err = H5Awrite(attr, H5T_NATIVE_INT, &value);
	}

	H5Aclose(attr);
	H5Sclose(dataspace);

	if (err < 0) {
		return 1;
	}

	return 0;
}

int HdfIO::write_float_attr_from_dict(int image_index, string attr_name, const Dict & dict)
{
	if (dict.has_key(attr_name)) {
		return write_float_attr(image_index, attr_name, dict[attr_name]);
	}
	return 0;
}

int HdfIO::write_float_attr(int image_index, string attr_name, float value)
{
	set_dataset(image_index);
	return write_float_attr(attr_name, value);
}

int HdfIO::write_float_attr(string attr_name, float value)
{
	int err = -1;
	delete_attr(attr_name);
	hid_t dataspace = H5Screate(H5S_SCALAR);
	hid_t attr =
		H5Acreate(cur_dataset, attr_name.c_str(), H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT);

	if (attr >= 0) {
		err = H5Awrite(attr, H5T_NATIVE_FLOAT, &value);
	}

	H5Aclose(attr);
	H5Sclose(dataspace);

	if (err < 0) {
		return 1;
	}

	return 0;
}

int HdfIO::write_string_attr(int image_index, string attr_name, string value)
{
	set_dataset(image_index);

	int err = -1;
	delete_attr(attr_name);

	hid_t datatype = H5Tcopy(H5T_C_S1);
	H5Tset_size(datatype, value.size() + 1);
	hid_t dataspace = H5Screate(H5S_SCALAR);

	hid_t attr = H5Acreate(cur_dataset, attr_name.c_str(), datatype, dataspace, H5P_DEFAULT);
	if (attr >= 0) {
		err = H5Awrite(attr, datatype, (void *) value.c_str());
	}

	H5Aclose(attr);
	H5Sclose(dataspace);
	H5Tclose(datatype);

	if (err < 0) {
		return 1;
	}

	return 0;
}

int HdfIO::write_array_attr(int image_index, string attr_name,
							int nitems, void *data, DataType type)
{
	assert(nitems > 0);
	assert(data != 0);
	set_dataset(image_index);
	delete_attr(attr_name);

	if (type != INT && type != FLOAT) {
		fprintf(stderr, "can only write INTEGER and FLOAT array");
		return 1;
	}

	int err = -1;
	hsize_t dim[1];
	int perm[1];
	dim[0] = nitems;
	hsize_t sdim[] = { 1 };

	hid_t datatype = -1;

	if (type == INT) {
		datatype = H5Tarray_create(H5T_NATIVE_INT, 1, dim, perm);
	}
	else if (type == FLOAT) {
		datatype = H5Tarray_create(H5T_NATIVE_FLOAT, 1, dim, perm);
	}

	hid_t dataspace = H5Screate_simple(1, sdim, NULL);
	hid_t attr = H5Acreate(cur_dataset, attr_name.c_str(), datatype, dataspace, H5P_DEFAULT);
	if (attr >= 0) {
		err = H5Awrite(attr, datatype, data);
	}

	H5Tclose(datatype);
	H5Sclose(dataspace);
	H5Aclose(attr);

	if (err < 0) {
		return 1;
	}

	return 0;
}


int HdfIO::write_global_int_attr(string attr_name, int value)
{
	hid_t tmp_dataset = cur_dataset;
	cur_dataset = group;
	int err = write_int_attr(attr_name, value);
	cur_dataset = tmp_dataset;
	return err;
}


int HdfIO::write_euler_attr(int image_index, string attr_name, int value)
{
	set_dataset(image_index);
	delete_attr(attr_name);

	hsize_t dim[] = { 1 };
	hid_t dataspace = H5Screate_simple(1, dim, NULL);
	hid_t attr = H5Acreate(cur_dataset, attr_name.c_str(), euler_type, dataspace, H5P_DEFAULT);
	H5Awrite(attr, euler_type, &value);
	H5Sclose(dataspace);
	H5Aclose(attr);
	return 0;
}


int HdfIO::write_mapinfo_attr(int image_index, string attr_name, int value)
{
	set_dataset(image_index);
	delete_attr(attr_name);

	hsize_t dim[] = { 1 };
	hid_t dataspace = H5Screate_simple(1, dim, NULL);
	hid_t attr = H5Acreate(cur_dataset, attr_name.c_str(), mapinfo_type, dataspace, H5P_DEFAULT);
	H5Awrite(attr, mapinfo_type, &value);
	H5Sclose(dataspace);
	H5Aclose(attr);
	return 0;
}


int HdfIO::delete_attr(int image_index, string attr_name)
{
	set_dataset(image_index);

	hdf_err_off();
	int err = H5Adelete(cur_dataset, attr_name.c_str());
	hdf_err_on();

	if (err >= 0) {
		return 0;
	}
	else {
		return 1;
	}
}

int HdfIO::delete_attr(string attr_name)
{
	hdf_err_off();
	int err = H5Adelete(cur_dataset, attr_name.c_str());
	hdf_err_on();

	if (err >= 0) {
		return 0;
	}
	else {
		return 1;
	}
}

string HdfIO::get_compound_name(int id, string name)
{
	string magic = get_item_name(COMPOUND_DATA_MAGIC);
	char id_str[32];
	sprintf(id_str, "%d", id);
	string compound_name = magic + id_str + name;
	return compound_name;
}


int HdfIO::create_compound_attr(int image_index, string attr_name)
{
	string cur_dataset_name = get_compound_name(image_index, attr_name);
	cur_image_index = -1;

	hsize_t dims[1];
	dims[0] = 1;
	hid_t datatype = H5Tcopy(H5T_NATIVE_INT);
	hid_t dataspace = H5Screate_simple(1, dims, NULL);

	close_dataset(cur_dataset);
	cur_dataset = H5Dcreate(file, cur_dataset_name.c_str(), datatype, dataspace, H5P_DEFAULT);

	return 0;
}

herr_t attr_info(hid_t dataset, const char *name, void *opdata)
{
	hid_t attr = H5Aopen_name(dataset, name);
	float value = 0;
	Dict *dict = (Dict *) opdata;

	if (attr >= 0) {
		hid_t atype = H5Aget_type(attr);
		if (H5Tget_class(atype) == H5T_FLOAT) {
			H5Aread(attr, H5T_NATIVE_FLOAT, &value);
			(*dict)[name] = value;
		}
		else {
			LOGERR("can only handle float CTF parameters in HDF");
			exit(1);
		}
		H5Aclose(attr);
	}


	return 0;
}

// assume all ctf fields are floats.
int HdfIO::read_ctf(Ctf & ctf, int image_index)
{
	LOGDEBUG("HdfIO::read_ctfit()");
	if (init() != 0) {
		return 1;
	}

	int err = 0;
	set_dataset(image_index);

	string cur_dataset_name = get_compound_name(image_index, get_item_name(CTFIT));
	cur_dataset = H5Dopen(file, cur_dataset_name.c_str());

	if (cur_dataset < 0) {
		err = 1;
	}
	else {
		Dict dict;
		err = H5Aiterate(cur_dataset, 0, attr_info, &dict);
		if (err >= 0) {
			ctf.from_dict(dict);
		}
		else {
			err = 1;
		}
	}

	cur_image_index = -1;
	return err;
}

// assume all ctf fields are floats.
int HdfIO::write_ctf(const Ctf & ctf, int image_index)
{
	LOGDEBUG("HdfIO::write_ctfit()");
	if (init() != 0) {
		return 1;
	}

	set_dataset(image_index);

	string attr_name = get_item_name(CTFIT);
	string cur_dataset_name = get_compound_name(image_index, attr_name);

	cur_dataset = H5Dopen(file, cur_dataset_name.c_str());

	if (cur_dataset < 0) {
		create_compound_attr(image_index, attr_name);
	}


	Dict dict = ctf.to_dict();
	vector < string > keys = dict.keys();
	for (size_t i = 0; i < keys.size(); i++) {
		float v = dict[keys[i]];
		write_float_attr(keys[i].c_str(), v);
	}

	cur_image_index = -1;

	return 0;
}

bool HdfIO::is_complex_mode()
{
	return false;
}

// always big endian
bool HdfIO::is_image_big_endian()
{
	return true;
}


string HdfIO::get_item_name(Nametype type)
{
	switch (type) {
	case ROOT_GROUP:
		return "/";
	case CTFIT:
		return "ctfit";
	case NUMDATASET:
		return "num_dataset";
	case COMPOUND_DATA_MAGIC:
		return "compound.";
	}


	return "unknown";
}

void HdfIO::hdf_err_off()
{
	H5Eget_auto(&old_func, &old_client_data);
	H5Eset_auto(0, 0);
}

void HdfIO::hdf_err_on()
{
	H5Eset_auto(old_func, old_client_data);
}

void HdfIO::close_dataset(hid_t dataset)
{
	hdf_err_off();
	if (dataset >= 0) {
		H5Dclose(dataset);
	}
	hdf_err_on();
}

// EMAN,IMAGIC,SPIN,QUAT,MATRIX,SGIROT,SPIDER,MRC

void HdfIO::create_enum_types()
{
	static int enum_types_created = 0;

	if (!enum_types_created) {
#if 0
		Euler::EulerType e;
		euler_type = H5Tcreate(H5T_ENUM, sizeof(Euler::EulerType));
		H5Tenum_insert(euler_type, "EMAN", (e = Euler::EMAN, &e));
		H5Tenum_insert(euler_type, "IMAGIC", (e = Euler::IMAGIC, &e));
		H5Tenum_insert(euler_type, "SPIN", (e = Euler::SPIN, &e));
		H5Tenum_insert(euler_type, "QUAT", (e = Euler::QUAT, &e));
		H5Tenum_insert(euler_type, "MATRIX", (e = Euler::MATRIX, &e));
		H5Tenum_insert(euler_type, "SGIROT", (e = Euler::SGIROT, &e));
		H5Tenum_insert(euler_type, "SPIDER", (e = Euler::SPIDER, &e));
		H5Tenum_insert(euler_type, "MRC", (e = Euler::MRC, &e));

		MapInfoType m;
		mapinfo_type = H5Tcreate(H5T_ENUM, sizeof(MapInfoType));

		H5Tenum_insert(mapinfo_type, "NORMAL", (m = NORMAL, &m));
		H5Tenum_insert(mapinfo_type, "ICOS2F_FIRST_OCTANT", (m = ICOS2F_FIRST_OCTANT, &m));
		H5Tenum_insert(mapinfo_type, "ICOS2F_FULL", (m = ICOS2F_FULL, &m));
		H5Tenum_insert(mapinfo_type, "ICOS2F_HALF", (m = ICOS2F_HALF, &m));
		H5Tenum_insert(mapinfo_type, "ICOS3F_HALF", (m = ICOS3F_HALF, &m));
		H5Tenum_insert(mapinfo_type, "ICOS3F_FULL", (m = ICOS3F_FULL, &m));
		H5Tenum_insert(mapinfo_type, "ICOS5F_HALF", (m = ICOS5F_HALF, &m));
		H5Tenum_insert(mapinfo_type, "ICOS5F_FULL", (m = ICOS5F_FULL, &m));
#endif
		enum_types_created = 1;
	}
}

vector < int >HdfIO::get_image_indices()
{
	return image_indices;
}

int HdfIO::get_hdf_dims(int image_index, int *p_nx, int *p_ny, int *p_nz)
{
	*p_nx = read_int_attr(image_index, "nx");
	*p_ny = read_int_attr(image_index, "ny");
	*p_nz = read_int_attr(image_index, "nz");

	if (*p_nx == 0 || *p_ny == 0 || *p_nz == 0) {
		int ndim = 0;
		int *dims = read_dims(image_index, &ndim);

		if (ndim != 2 && ndim != 3) {
			LOGERR("only handle 2D/3D HDF5. Your file is %dD.", ndim);
			return 1;
		}
		else {
			*p_nx = dims[0];
			*p_ny = dims[1];
			if (ndim == 3) {
				*p_nz = dims[2];
			}
			else {
				*p_nz = 1;
			}
		}
	}
	return 0;
}

herr_t HdfIO::file_info(hid_t loc_id, const char *name, void *opdata)
{
	loc_id = loc_id;
	vector < int >*image_indices = static_cast < vector < int >*>(opdata);
	string magic = HdfIO::get_item_name(HdfIO::COMPOUND_DATA_MAGIC);
	const char *magic_str = magic.c_str();

	if (strncmp(name, magic_str, strlen(magic_str)) != 0) {
		int id = atoi(name);
		image_indices->push_back(id);
	}

	return 0;
}

hid_t HdfIO::create_dataset(int image_index, int nx, int ny, int nz)
{
	int ndim = 0;
	int dims[3];

	if (nz == 1) {
		ndim = 2;
	}
	else {
		ndim = 3;
	}
	dims[0] = nx;
	dims[1] = ny;
	dims[2] = nz;

	char tmp_dataset_name[32];
	sprintf(tmp_dataset_name, "%d", image_index);

	hdf_err_off();
	hid_t tmp_dataset = H5Dopen(file, tmp_dataset_name);
	hdf_err_on();

	if (tmp_dataset >= 0) {
		int ndim1 = 0;
		int *dims1 = read_dims(image_index, &ndim1);
		assert(ndim == ndim1);

		for (int i = 0; i < ndim; i++) {
			assert(dims[i] == dims1[i]);
		}
	}
	else {
		hsize_t *sdims = new hsize_t[ndim];
		for (int i = 0; i < ndim; i++) {
			sdims[i] = dims[i];
		}

		hid_t datatype = H5Tcopy(H5T_NATIVE_FLOAT);
		hid_t dataspace = H5Screate_simple(ndim, sdims, NULL);

		tmp_dataset = H5Dcreate(file, tmp_dataset_name, datatype, dataspace, H5P_DEFAULT);

		H5Tclose(datatype);
		H5Sclose(dataspace);

		delete[]sdims;
		sdims = 0;
	}

	return tmp_dataset;
}



#endif
