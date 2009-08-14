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

#ifdef EM_HDF5

#include "hdfio.h"
#include "geometry.h"
#include "ctf.h"
#include "emassert.h"
#include "transform.h"
#include <iostream>
#include <cstring>

#ifndef _WIN32
	#include <sys/param.h>
#else
	#include <windows.h>
	#define  MAXPATHLEN (MAX_PATH * 4)
#endif	//_WIN32

using namespace EMAN;

hid_t HdfIO::mapinfo_type = -1;
const char *HdfIO::HDF5_SIGNATURE = "\211HDF\r\n\032\n";

herr_t attr_info(hid_t dataset, const char *name, void *opdata)
{
	hid_t attr = H5Aopen_name(dataset, name);
	float value_float = 0.0f;
	int value_int = 0;
	string value_string = "";
	char * tmp_string = new char[1024];
	Dict *dict = (Dict *) opdata;

	if (attr >= 0) {
		hid_t atype = H5Aget_type(attr);
		if (H5Tget_class(atype) == H5T_FLOAT) {
			H5Aread(attr, atype, &value_float);
			(*dict)[name] = value_float;
		}
		else if(H5Tget_class(atype) == H5T_INTEGER) {
			H5Aread(attr, atype, &value_int);
			(*dict)[name] = value_int;
		}
		else if(H5Tget_class(atype) == H5T_STRING) {
			H5Aread(attr, atype, tmp_string);
			value_string = tmp_string;
			(*dict)[name] = value_string;
		}
		else if(H5Tget_class(atype) == H5T_ENUM || H5Tget_class(atype) == H5T_ARRAY) {
			//for old-style hdf file created by EMAN1, the euler convention is enum
			//skip those attribute to make EMAN2 read the content of the image
		}
		else {
			LOGERR("can only handle float/int/string parameters in HDF attr_info()");
			exit(1);
		}
		H5Tclose(atype);
		H5Aclose(attr);
	}

	return 0;
}

HdfIO::HdfIO(const string & hdf_filename, IOMode rw)
:	filename(hdf_filename), rw_mode(rw)
{
	initialized = false;
	is_new_file = false;
	file = -1;
	group = -1;
	cur_dataset = -1;
	cur_image_index = -1;
	old_func = 0;
	old_client_data = 0;
//	printf("HDf open\n");
}

HdfIO::~HdfIO()
{
	close_cur_dataset();
    if (group >= 0) {
        H5Gclose(group);
    }
    if (file >= 0) {
		H5Fflush(file,H5F_SCOPE_GLOBAL);	// If there were no resource leaks, this wouldn't be necessary...
		H5Fclose(file);
    }
//printf("HDf close\n");
}

void HdfIO::init()
{
	ENTERFUNC;
	if (initialized) {
		return;
	}

	initialized = true;

	FILE *tmp_file = sfopen(filename, rw_mode, &is_new_file);

	if (!is_new_file) {
		char buf[128];
		if (fread(buf, sizeof(buf), 1, tmp_file) != 1) {
			fclose(tmp_file);
			throw ImageReadException(filename, "read HDF5 first block");
		}
		else {
			if (!is_valid(buf)) {
				fclose(tmp_file);
				throw ImageReadException(filename, "invalid HDF5 file");
			}
		}
	}

	fclose(tmp_file);
	tmp_file = 0;

	if (rw_mode == READ_ONLY) {
		file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	}
	else {
		hdf_err_off();
		file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
		hdf_err_on();
		if (file < 0) {
			file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
			H5Fclose(file);
			file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
		}
	}

	if (file < 0) {
		throw FileAccessException(filename);
	}

	string root_group_str = get_item_name(ROOT_GROUP);
	group = H5Gopen(file, root_group_str.c_str());
	cur_dataset = -1;

	H5Giterate(file, root_group_str.c_str(), NULL, file_info, &image_indices);
	create_enum_types();
	EXITFUNC;

}


bool HdfIO::is_valid(const void *first_block)
{
	ENTERFUNC;

	if (first_block) {
		char signature[8] = { 137,72,68,70,13,10,26,10 };
		if (strncmp((const char *)first_block,signature,8)==0) return true;
		// const char* f=(const char *)first_block;
		// printf("bad hdf signature %d %d %d %d %d %d %d %d",f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7]);
		return false;
	}
	EXITFUNC;
	return false;
}

int HdfIO::read_header(Dict & dict, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	check_read_access(image_index);

	int nx = 0, ny = 0, nz = 0;
	if (get_hdf_dims(image_index, &nx, &ny, &nz) != 0) {
		throw ImageReadException(filename, "invalid image dimensions");
	}

	check_region(area, IntSize(nx, ny, nz));

	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen, nz, &zlen);

	set_dataset(image_index);
	int err = 0;
	if (cur_dataset < 0) {
		err = 1;
	}
	else {
		err = H5Aiterate(cur_dataset, 0, attr_info, &dict);
		if (err < 0) {
			err = 1;
		}
	}

	if(dict.has_key("ctf")) {
		Ctf * ctf_ = new EMAN1Ctf();
		ctf_->from_string((string)dict["ctf"]);
		dict.erase("ctf");
		dict["ctf"] = ctf_;
		delete ctf_;
	}

	dict["nx"] = xlen;
	dict["ny"] = ylen;
	dict["nz"] = zlen;

    Dict euler_angles;
    read_euler_angles(euler_angles, image_index);

    vector<string> euler_names = euler_angles.keys();
    for (size_t i = 0; i < euler_names.size(); i++) {
        dict[euler_names[i]] = euler_angles[euler_names[i]];
    }

	dict["datatype"] = EMUtil::EM_FLOAT;
	EXITFUNC;
	return 0;
}


int HdfIO::read_data(float *data, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	check_read_access(image_index, data);
	set_dataset(image_index);
#if 0
	if (cur_dataset < 0) {
		char cur_dataset_name[32];
		sprintf(cur_dataset_name, "%d", image_index);
		cur_dataset = H5Dopen(file, cur_dataset_name);

		if (cur_dataset < 0) {
			char desc[256];
			sprintf(desc, "no image with id = %d", image_index);
			throw ImageReadException(filename, desc);
		}
	}
#endif
	hid_t datatype = H5Dget_type(cur_dataset);
	H5T_class_t t_class = H5Tget_class(datatype);
	H5Tclose(datatype);

	if (t_class != H5T_FLOAT) {
		char desc[256];
		sprintf(desc, "unsupported HDF5 data type '%d'", (int) t_class);
		throw ImageReadException(filename, desc);
	}

	int nx = 0, ny = 0, nz = 0;
	if (get_hdf_dims(image_index, &nx, &ny, &nz) != 0) {
		throw ImageReadException(filename, "invalid image dimensions");
	}

	check_region(area, FloatSize(nx, ny, nz));

	int err = 0;
	if (!area) {
		err = H5Dread(cur_dataset, H5T_NATIVE_FLOAT, H5S_ALL,
					  H5S_ALL, H5P_DEFAULT, data);
	}
	else {
		hid_t dataspace_id = 0;
		hid_t memspace_id = 0;

		err = create_region_space(&dataspace_id, &memspace_id, area,
								  nx, ny, nz, image_index);
		if (err == 0) {
			err = H5Dread(cur_dataset, H5T_NATIVE_FLOAT, memspace_id,
						  dataspace_id, H5P_DEFAULT, data);
		}

		H5Sclose(dataspace_id);
		H5Sclose(memspace_id);
		if (err < 0) {
			throw ImageReadException(filename,
									 "creating memory space or file space id failed");
		}
	}

	if (err < 0) {
		char desc[256];
		sprintf(desc, "reading %dth HDF5 image failed", image_index);
		throw ImageReadException(filename, desc);
	}

	EXITFUNC;
	return 0;
}


int HdfIO::write_header(const Dict & dict, int image_index, const Region* area,
						EMUtil::EMDataType, bool)
{
	ENTERFUNC;
	check_write_access(rw_mode, image_index);

	if (area) {
		int nx0 = 0;
		int ny0 = 0;
		int nz0 = 0;

		if (get_hdf_dims(image_index, &nx0, &ny0, &nz0) != 0) {
			throw ImageReadException(filename, "invalid image dimensions");
		}

		check_region(area, FloatSize(nx0, ny0, nz0), is_new_file);

		EXITFUNC;
		return 0;
	}

	int nx = dict["nx"];
	int ny = dict["ny"];
	int nz = dict["nz"];

	create_cur_dataset(image_index, nx, ny, nz);
	Assert(cur_dataset >= 0);

	vector<string> keys = dict.keys();
	vector<string>::const_iterator iter;
	for (iter = keys.begin(); iter != keys.end(); iter++) {
		//handle special case for euler anglers
		if(*iter == "orientation_convention") {
			string eulerstr = (const char*) dict["orientation_convention"];
        	write_string_attr(image_index, "orientation_convention", eulerstr);

        	vector<string> euler_names = EMUtil::get_euler_names(eulerstr);
        	Dict euler_dict;
        	for (size_t i = 0; i < euler_names.size(); i++) {
            	euler_dict[euler_names[i]] = dict[euler_names[i]];
        	}

        	write_euler_angles(euler_dict, image_index);
		}
		//micrograph_id should save as string
		else if(*iter == "micrograph_id") {
			write_string_attr(image_index, "micrograph_id", (const char *) dict["micrograph_id"]);
		}
		//handle normal attributes
		else {
			EMObject attr_val = dict[*iter];
			//string val_type = EMObject::get_object_type_name(attr_val.get_type());
			EMObject::ObjectType t = attr_val.get_type();
			switch(t) {

				case EMObject::INT:
					write_int_attr(image_index, *iter, attr_val);
					break;
				case EMObject::FLOAT:
				case EMObject::DOUBLE:
					write_float_attr(image_index, *iter, attr_val);
					break;
				case EMObject::STRING:
					write_string_attr(image_index, *iter, attr_val.to_str());
					break;
//				case EMObject::EMDATA:
				case EMObject::XYDATA:
				case EMObject::FLOATARRAY:
				case EMObject::STRINGARRAY:
					throw NotExistingObjectException("EMObject", "unsupported type");
					break;
				case EMObject::UNKNOWN:
					throw NotExistingObjectException("EMObject", "unsupported type");
					break;
				default:
					throw NotExistingObjectException("EMObject", "unsupported type");
			}
		}
	}


	//flush();
	close_cur_dataset();

	EXITFUNC;
	return 0;
}

int HdfIO::write_data(float *data, int image_index, const Region* area,
					  EMUtil::EMDataType, bool)
{
	ENTERFUNC;

	check_write_access(rw_mode, image_index, 0, data);

	int nx = read_int_attr(image_index, "nx");
	int ny = read_int_attr(image_index, "ny");
	int nz = read_int_attr(image_index, "nz");


	check_region(area, FloatSize(nx, ny, nz), is_new_file);
	create_cur_dataset(image_index, nx, ny, nz);
	Assert(cur_dataset >= 0);

	int err = 0;

	if (!area) {
		err = H5Dwrite(cur_dataset, H5T_NATIVE_FLOAT, H5S_ALL,
					   H5S_ALL, H5P_DEFAULT, data);
		//if (err >= 0) {
		//increase_num_dataset();
		//image_indices.push_back(image_index);
		//}
	}
	else {
		hid_t dataspace_id = 0;
		hid_t memspace_id = 0;

		err = create_region_space(&dataspace_id, &memspace_id, area,
								  nx, ny, nz, image_index);
		if (err == 0) {
			err = H5Dwrite(cur_dataset, H5T_NATIVE_FLOAT, memspace_id,
						   dataspace_id, H5P_DEFAULT, data);
		}

		H5Sclose(dataspace_id);
		H5Sclose(memspace_id);
		if (err < 0) {
			throw ImageReadException(filename,
									 "creating memory space or file space id failed");
		}
	}

	if (err < 0) {
		throw ImageWriteException(filename, "HDF data write failed");
	}

	//flush();
	close_cur_dataset();

	EXITFUNC;
	return 0;
}

void HdfIO::flush()
{
	if (cur_dataset > 0) {
		H5Fflush(cur_dataset, H5F_SCOPE_LOCAL);
	}
}

int *HdfIO::read_dims(int image_index, int *p_ndim)
{
	set_dataset(image_index);
#if 0
	if (cur_dataset < 0) {
		char cur_dataset_name[32];
		sprintf(cur_dataset_name, "%d", image_index);
		cur_dataset = H5Dopen(file, cur_dataset_name);

		if (cur_dataset < 0) {
			throw ImageReadException(filename, "reading data dimensions");
		}
	}
#endif
	hid_t dataspace = H5Dget_space(cur_dataset);
	int rank = H5Sget_simple_extent_ndims(dataspace);
	hsize_t *dims = new hsize_t[rank];
	H5Sget_simple_extent_dims(dataspace, dims, NULL);

	int *dims1 = new int[rank];
	for (int i = 0; i < rank; i++) {
		dims1[i] = static_cast < int >(dims[i]);
	}

	H5Sclose(dataspace);
	if( dims )
	{
		delete[]dims;
		dims = 0;
	}

	(*p_ndim) = rank;
	return dims1;
}

int HdfIO::read_global_int_attr(const string & attr_name)
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

float HdfIO::read_global_float_attr(const string & attr_name)
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
	init();
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

	if (need_update || cur_dataset < 0) {
		char cur_dataset_name[32];
		sprintf(cur_dataset_name, "%d", image_index);
		hdf_err_off();
		close_cur_dataset();
		cur_dataset = H5Dopen(file, cur_dataset_name);
		if (cur_dataset < 0) {
			throw ImageReadException(filename, "open data set failed");
		}
		hdf_err_on();
		cur_image_index = image_index;
	}

	Assert(cur_dataset >= 0);
}



int HdfIO::read_int_attr(int image_index, const string & attr_name)
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


float HdfIO::read_float_attr(int image_index, const string & attr_name)
{
	set_dataset(image_index);
	return read_float_attr(attr_name);
}


float HdfIO::read_float_attr(const string & attr_name)
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



string HdfIO::read_string_attr(int image_index, const string & attr_name)
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

		if( tmp_value )
		{
			delete[]tmp_value;
			tmp_value = 0;
		}
	}

	return value;
}

int HdfIO::read_array_attr(int image_index, const string & attr_name, void *value)
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
#if 0
float HdfIO::read_euler_attr(int image_index, const string &attr_name)
{
    string full_attr_name = string(EULER_MAGIC) + attr_name;
    return read_float_attr(image_index, full_attr_name);
}
#endif

int HdfIO::read_mapinfo_attr(int image_index, const string & attr_name)
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

int HdfIO::write_int_attr(int image_index, const string & attr_name, int value)
{
	set_dataset(image_index);
	return write_int_attr(attr_name, value);
}


int HdfIO::write_int_attr(const string & attr_name, int value)
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

int HdfIO::write_float_attr_from_dict(int image_index, const string & attr_name,
									  const Dict & dict)
{
	if (dict.has_key(attr_name)) {
		return write_float_attr(image_index, attr_name, dict[attr_name]);
	}
	return 0;
}

int HdfIO::write_float_attr(int image_index, const string & attr_name, float value)
{
	set_dataset(image_index);
	return write_float_attr(attr_name, value);
}

int HdfIO::write_float_attr(const string & attr_name, float value)
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

int HdfIO::write_string_attr(int image_index, const string & attr_name,
							 const string & value)
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

int HdfIO::write_array_attr(int image_index, const string & attr_name,
							int nitems, void *data, DataType type)
{
	if (nitems <= 0) {
		return 1;
	}
	if (!data) {
		throw NullPointerException("array data is NULL");
	}

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


int HdfIO::write_global_int_attr(const string & attr_name, int value)
{
	hid_t tmp_dataset = cur_dataset;
	cur_dataset = group;
	int err = write_int_attr(attr_name, value);
	cur_dataset = tmp_dataset;
	return err;
}

#if 0
int HdfIO::write_euler_attr(int image_index, const string & attr_name, float value)
{
    string full_attr_name = string(EULER_MAGIC) + attr_name;
    return write_float_attr(image_index, full_attr_name, value);
}
#endif

int HdfIO::write_mapinfo_attr(int image_index, const string & attr_name, int value)
{
	set_dataset(image_index);
	delete_attr(attr_name);

	hsize_t dim[] = { 1 };
	hid_t dataspace = H5Screate_simple(1, dim, NULL);
	hid_t attr = H5Acreate(cur_dataset, attr_name.c_str(), mapinfo_type, dataspace, H5P_DEFAULT);
	H5Awrite(attr, mapinfo_type, &value);
	H5Aclose(attr);
	H5Sclose(dataspace);
	return 0;
}


int HdfIO::delete_attr(int image_index, const string & attr_name)
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

int HdfIO::delete_attr(const string & attr_name)
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

string HdfIO::get_compound_name(int id, const string & name)
{
	string magic = get_item_name(COMPOUND_DATA_MAGIC);
	char id_str[32];
	sprintf(id_str, "%d", id);
	string compound_name = magic + "." + id_str + "." + name;
	return compound_name;
}


int HdfIO::create_compound_attr(int image_index, const string & attr_name)
{
	string cur_dataset_name = get_compound_name(image_index, attr_name);
	cur_image_index = -1;

	hsize_t dims[1];
	dims[0] = 1;
	hid_t datatype = H5Tcopy(H5T_NATIVE_INT);
	hid_t dataspace = H5Screate_simple(1, dims, NULL);

	close_cur_dataset();
	cur_dataset = H5Dcreate(file, cur_dataset_name.c_str(), datatype, dataspace, H5P_DEFAULT);

	H5Tclose(datatype);
	H5Sclose(dataspace);
	return 0;
}

int HdfIO::read_ctf(Ctf & ctf, int image_index)
{
    Dict ctf_dict;
    int err = read_compound_dict(CTFIT, ctf_dict, image_index);
    if (!err) {
        ctf.from_dict(ctf_dict);
    }

	return err;
}

void HdfIO::write_ctf(const Ctf & ctf, int image_index)
{
    Dict ctf_dict = ctf.to_dict();
    write_compound_dict(CTFIT, ctf_dict, image_index);
}

int HdfIO::read_euler_angles(Dict & euler_angles, int image_index)
{
    int err = read_compound_dict(EULER, euler_angles, image_index);
    return err;
}


void HdfIO::write_euler_angles(const Dict & euler_angles, int image_index)
{
    write_compound_dict(EULER, euler_angles, image_index);
}

int HdfIO::read_compound_dict(Nametype compound_type,
                              Dict & values, int image_index)
{
    ENTERFUNC;
	init();

	int err = 0;

	hid_t cur_dataset_orig = cur_dataset;
	string cur_dataset_name = get_compound_name(image_index, get_item_name(compound_type));

    hdf_err_off();
    cur_dataset = H5Dopen(file, cur_dataset_name.c_str());
    hdf_err_on();

	if (cur_dataset < 0) {
		err = 1;
	}
	else {
		err = H5Aiterate(cur_dataset, 0, attr_info, &values);
		if (err < 0) {
			err = 1;
		}
	}

	H5Dclose(cur_dataset);
    cur_dataset = cur_dataset_orig;

	cur_image_index = -1;
	EXITFUNC;
	return err;
}


void HdfIO::write_compound_dict(Nametype compound_type,
                                const Dict & values, int image_index)
{
	ENTERFUNC;
	init();

	//set_dataset(image_index);
	hid_t cur_dataset_orig = cur_dataset;

	string attr_name = get_item_name(compound_type);
	string cur_dataset_name = get_compound_name(image_index, attr_name);

    hdf_err_off();
	cur_dataset = H5Dopen(file, cur_dataset_name.c_str());
    hdf_err_on();

	if (cur_dataset < 0) {
		create_compound_attr(image_index, attr_name);
	}
    else {
        // remove all existing attributes first
        Dict attr_dict;
        H5Aiterate(cur_dataset, 0, attr_info, &attr_dict);
        vector <string> attr_keys = attr_dict.keys();
        for (size_t i = 0; i < attr_keys.size(); i++) {
            H5Adelete(cur_dataset, attr_keys[i].c_str());
        }
    }

    // writing new attributes

	vector < string > keys = values.keys();
	for (size_t i = 0; i < keys.size(); i++) {
		float v = values[keys[i]];
		write_float_attr(keys[i].c_str(), v);
	}

	H5Dclose(cur_dataset);
	cur_dataset = cur_dataset_orig;

	cur_image_index = -1;
	EXITFUNC;
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
		return "compound";
    case EULER:
        return "euler_angles";
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

void HdfIO::close_cur_dataset()
{
	hdf_err_off();
	if (cur_dataset >= 0) {
		H5Dclose(cur_dataset);
		cur_dataset = -1;
		cur_image_index = -1;
	}
	hdf_err_on();
}

void HdfIO::create_enum_types()
{
	static int enum_types_created = 0;

	if (!enum_types_created) {
#if 0
        Transform3D::EulerType e;
		euler_type = H5Tcreate(H5T_ENUM, sizeof(Transform3D::EulerType));

		H5Tenum_insert(euler_type, "EMAN", (e = Transform3D::EMAN, &e));
		H5Tenum_insert(euler_type, "IMAGIC", (e = Transform3D::IMAGIC, &e));
		H5Tenum_insert(euler_type, "SPIN", (e = Transform3D::SPIN, &e));
		H5Tenum_insert(euler_type, "QUATERNION", (e = Transform3D::QUATERNION, &e));
		H5Tenum_insert(euler_type, "SGIROT", (e = Transform3D::SGIROT, &e));
		H5Tenum_insert(euler_type, "SPIDER", (e = Transform3D::SPIDER, &e));
		H5Tenum_insert(euler_type, "MRC", (e = Transform3D::MRC, &e));

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
			if( dims )
			{
				delete [] dims;
				dims = 0;
			}
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
		if( dims )
		{
			delete [] dims;
			dims = 0;
		}
	}
	return 0;
}

herr_t HdfIO::file_info(hid_t, const char *name, void *opdata)
{
//	loc_id = loc_id;
	vector < int >*image_indices = static_cast < vector < int >*>(opdata);
	string magic = HdfIO::get_item_name(HdfIO::COMPOUND_DATA_MAGIC);
	const char *magic_str = magic.c_str();

	if (strncmp(name, magic_str, strlen(magic_str)) != 0) {
		int id = atoi(name);
		image_indices->push_back(id);
	}

	return 0;
}

void HdfIO::create_cur_dataset(int image_index, int nx, int ny, int nz)
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
	cur_image_index = image_index;

	hdf_err_off();
	cur_dataset = H5Dopen(file, tmp_dataset_name);
	hdf_err_on();


	if (cur_dataset >= 0) {
		int ndim1 = 0;
		int *dims1 = read_dims(image_index, &ndim1);
		Assert(ndim == ndim1);

		for (int i = 0; i < ndim; i++) {
			Assert(dims[i] == dims1[i]);
		}
		if( dims1 )
		{
			delete [] dims1;
			dims1 = 0;
		}
	}
	else {
		hsize_t *sdims = new hsize_t[ndim];
		for (int i = 0; i < ndim; i++) {
			sdims[i] = dims[i];
		}

		hid_t datatype = H5Tcopy(H5T_NATIVE_FLOAT);
		hid_t dataspace = H5Screate_simple(ndim, sdims, NULL);

		cur_dataset = H5Dcreate(file, tmp_dataset_name, datatype, dataspace, H5P_DEFAULT);

		H5Tclose(datatype);
		H5Sclose(dataspace);

		if( sdims )
		{
			delete[]sdims;
			sdims = 0;
		}

		if (cur_dataset < 0) {
			throw ImageWriteException(filename, "create dataset failed");
		}
		else {
			increase_num_dataset();
			image_indices.push_back(image_index);
		}
	}
}

int HdfIO::create_region_space(hid_t * p_dataspace_id, hid_t * p_memspace_id,
							   const Region * area, int nx, int ny, int nz,
							   int image_index)
{
	Assert(p_dataspace_id);
	Assert(p_memspace_id);
	Assert(area);

	if (!p_dataspace_id || !p_memspace_id || !area) {
		return -1;
	}

#if H5_VERS_MINOR > 6 || (H5_VERS_MINOR == 6 && H5_VERS_RELEASE >= 4)
	hsize_t offset[3];
#else
	hssize_t offset[3];
#endif

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

	*p_dataspace_id = H5Dget_space(cur_dataset);

	int err = H5Sselect_hyperslab(*p_dataspace_id, H5S_SELECT_SET,
								  offset, NULL, count, NULL);
	if (err >= 0) {
		*p_memspace_id = H5Screate_simple(3, count, NULL);

#if H5_VERS_MINOR > 6 || (H5_VERS_MINOR == 6 && H5_VERS_RELEASE >= 4)
		hsize_t offset_out[3];
#else
		hssize_t offset_out[3];
#endif

		offset_out[0] = 0;
		offset_out[1] = 0;
		offset_out[2] = 0;

		err = H5Sselect_hyperslab(*p_memspace_id, H5S_SELECT_SET,
								  offset_out, NULL, count, NULL);
		if (err >= 0) {
			err = 0;
		}
	}

	return err;
}

int HdfIO::get_num_dataset()
{
    hdf_err_off();
    int n = read_global_int_attr(get_item_name(NUMDATASET));
    hdf_err_on();
    return n;
}


#endif	//EM_HDF5
