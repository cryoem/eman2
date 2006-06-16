/**
 * $Id$
 */
#ifdef EM_HDF5

#include "hdfio.h"
#include "geometry.h"
#include "ctf.h"
#include "Assert.h"
#include "transform.h"
#include <iostream>

#ifndef WIN32
	#include <sys/param.h>
#else
	#define  MAXPATHLEN (MAX_PATH * 4)
#endif	//WIN32

using namespace EMAN;

hid_t HdfIO::mapinfo_type = -1;
const char *HdfIO::HDF5_SIGNATURE = "\211HDF\r\n\032\n";

HdfIO::HdfIO(const string & hdf_filename, IOMode rw)
:	filename(hdf_filename), rw_mode(rw)
{
	initialized = false;
	is_new_file = false;
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

	H5Eset_auto(0, 0);	// Turn off console error logging.

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
		file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
		if (file < 0) file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
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
				case EMObject::EMDATA:
				case EMObject::XYDATA:
				case EMObject::FLOATARRAY:
				case EMObject::STRINGARRAY:
					throw NotExistingObjectException("EMObject", "unsupported type");
					break;
				case EMObject::UNKNOWN:
					std::cout << "Type::UNKNOWN the error attribute is: " << *iter << std::endl;
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

int HdfIO::get_nimg()
{
	init();
	hdf_err_off();
	int n = read_global_int_attr(get_item_name(NUMDATASET));
	hdf_err_on();
	return n;
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



#endif	//EM_HDF5
