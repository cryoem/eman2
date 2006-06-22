/**
 * $Id$
 */
#ifdef EM_HDF5

#include "hdfio2.h"
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

HdfIO2::HdfIO2(const string & hdf_filename, IOMode rw)
:	filename(hdf_filename), rw_mode(rw)
{
	initialized = false;
	file=-1;
	group=-1;
	hsize_t dims=1;
	simple_space=H5Screate_simple(1,&dims,NULL);
}

HdfIO2::~HdfIO2()
{
	H5Sclose(simple_space);
    if (group >= 0) {
        H5Gclose(group);
    }
    if (file >= 0) {
		H5Fflush(file,H5F_SCOPE_GLOBAL);	// If there were no resource leaks, this wouldn't be necessary...
		H5Fclose(file);
    }
//printf("HDf close\n");
}

EMObject HdfIO2::read_attr(hid_t loc,const char *name) {
	hid_t attr = H5Aopen_name(loc,name);
	hid_t type = H5Aget_type(attr);
	hid_t spc = H5Aget_space(attr);
	H5T_class_t cls = H5Tget_class(type);
	size_t sz = H5Tget_size(type);						// storage size, arrays handled in the 'space'
	hssize_t pts = H5Sget_simple_extent_npoints(spc);	// number of points > 1 if an array of floats

	EMObject ret(0);
	int i;
	float f,*fa;
	double d;
	char *s;
	vector <float> fv(pts);

	switch (cls) {
	case H5T_INTEGER:
		H5Aread(attr,H5T_NATIVE_INT,&i);
		ret=EMObject(i);
		break;
	case H5T_FLOAT:
		if (sz==4) {
			if (pts==1) {
				H5Aread(attr,H5T_NATIVE_FLOAT,&f);
				ret=EMObject(f);
			}
			else {
				fa=(float *)malloc(pts*sizeof(float));
				H5Aread(attr,H5T_NATIVE_FLOAT,fa);
				for (i=0; i<pts; i++) fv[i]=fa[i];
				free(fa);
				ret=EMObject(fv);
			}
		}
		else if (sz==8) {
			H5Aread(attr,H5T_NATIVE_DOUBLE,&d);
			ret=EMObject(d);
		}
		break;
	case H5T_STRING:
		s=(char *)malloc(sz+1);
		H5Aread(attr,H5T_NATIVE_CHAR,s);
		ret=EMObject(s);
		free(s);
		break;
	default:
		LOGERR("Unhandled HDF5 metadata %d", cls);
	}

	H5Sclose(spc);
	H5Tclose(type);
	H5Aclose(attr);

	return ret;
}

int HdfIO2::write_attr(hid_t loc,const char *name,EMObject obj) {
	hid_t type=0;
	hid_t spc=0;
	hsize_t dims=1;
	vector <float> fv;
	switch(obj.get_type()) {
	case EMObject::INT: type=H5Tcopy(H5T_NATIVE_INT); spc=H5Scopy(simple_space); break;
	case EMObject::FLOAT: type=H5Tcopy(H5T_NATIVE_FLOAT); spc=H5Scopy(simple_space); break;
	case EMObject::DOUBLE: type=H5Tcopy(H5T_NATIVE_DOUBLE); spc=H5Scopy(simple_space); break;
	case EMObject::STRING: 
		type=H5Tcopy(H5T_NATIVE_CHAR); 
		H5Tset_size(type,strlen((const char *)obj)+1);
		spc=H5Scopy(simple_space); 
		break;
	case EMObject::FLOATARRAY:
		type=H5Tcopy(H5T_NATIVE_FLOAT);
		fv=obj;
		dims=fv.size();
		simple_space=H5Screate_simple(1,&dims,NULL);
		break;
	case EMObject::STRINGARRAY:
	case EMObject::EMDATA:
	case EMObject::XYDATA:
		return -1;
		break;
	case EMObject::UNKNOWN:
		break;
	}

	H5Adelete(loc,name);
	hid_t attr = H5Acreate(loc,name,type,spc,H5P_DEFAULT);

	unsigned int i;
	float f,*fa;
	double d;
	const char *s;
	switch(obj.get_type()) {
	case EMObject::INT:
		i=(int)obj;
		H5Awrite(attr,H5T_NATIVE_INT,&i);
		break;
	case EMObject::FLOAT:
		f=(float)obj;
		H5Awrite(attr,H5T_NATIVE_FLOAT,&f);
		break;
	case EMObject::DOUBLE:
		d=(double)obj;
		H5Awrite(attr,H5T_NATIVE_DOUBLE,&d);
		break;
	case EMObject::STRING: 
		s=(const char *)obj;
		H5Awrite(attr,H5T_NATIVE_CHAR,s);
		break;
	case EMObject::FLOATARRAY:
		fa=(float *)malloc(fv.size()*sizeof(float));
		for (i=0; i<fv.size(); i++) fa[i]=fv[i];
		H5Awrite(attr,H5T_NATIVE_FLOAT,fa);
		free(fa);
		break;
//	case EMObject::STRINGARRAY:
//	case EMObject::EMDATA:
//	case EMObject::XYDATA:
//		return -1;
//		break;
	default:
		LOGERR("Unhandled HDF5 metadata '%s'", name);
		
	}


	H5Tclose(type);
	H5Sclose(spc);
}

void HdfIO2::init()
{
	ENTERFUNC;
	if (initialized) {
		return;
	}

	H5Eset_auto(0, 0);	// Turn off console error logging.

	if (rw_mode == READ_ONLY) {
		file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
		if (file<0) throw FileAccessException(filename);
	}
	else {
		file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
		if (file < 0) {
			file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
			if (file < 0) throw FileAccessException(filename);
		}
	}
	
	group=H5Gopen(file,"/MDF");
	if (group<0) {
		if (rw_mode == READ_ONLY) throw ImageReadException(filename,"HDF5 file has no image data (no /TEM group)");
		group=H5Gcreate(file,"/MDF",64);		// create the group for Macromolecular data
		if (group<0) throw ImageWriteException(filename,"Unable to add image group (/TEM) to HDF5 file");
		H5Gclose(group);
		group=H5Gcreate(file,"/MDF/images",4096);		// create the group for images/volumes
		if (group<0) throw ImageWriteException(filename,"Unable to add image group (/TEM/images) to HDF5 file");
		write_attr(group,"imageid_max",EMObject(-1));
	}
	initialized = true;
	EXITFUNC;

}

bool HdfIO2::is_valid(const void *first_block)
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

int HdfIO2::read_header(Dict & dict, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	EXITFUNC;
	return 0;
}


int HdfIO2::read_data(float *data, int image_index, const Region * area, bool)
{
	ENTERFUNC;
	
	EXITFUNC;
	return 0;
}


int HdfIO2::write_header(const Dict & dict, int image_index, const Region* area,
						EMUtil::EMDataType, bool)
{
	ENTERFUNC;
	// If image_index<0 append, and make sure the max value in the file is correct
	int nimg = read_attr(group,"imageid_max");
	unsigned int i;
	if (image_index<0) image_index=nimg+1;
	if (image_index>nimg) {
		write_attr(group,(const char *)"imageid_max",EMObject(image_index));
	}

	char ipath[50];
	sprintf(ipath,"/MDF/images/%d",image_index);
	hid_t igrp=H5Gopen(file,ipath);
	hid_t ds;
	if (igrp<0) {
		// Need to create a new image group
		igrp=H5Gcreate(file,ipath,64);		// The image is a group, with attributes on the group
		if (igrp<0) throw ImageWriteException(filename,"Unable to add /MDF/images/# to HDF5 file");

		// Now create the actual image dataset
		sprintf(ipath,"/MDF/images/%d/image",image_index);
		hsize_t dims[3]= { (int)dict["nx"],(int)dict["ny"],(int)dict["nz"] };
		hid_t space;
		if ((int)dict["nz"]==1) space=H5Screate_simple(2,dims,NULL);
		else space=H5Screate_simple(3,dims,NULL);
		ds=H5Dcreate(file,ipath, H5T_NATIVE_FLOAT, space, H5P_DEFAULT );
		H5Dclose(ds);
		H5Sclose(space);
	}
	
	// Write the attributes to the group
	vector <string> keys=dict.keys();

	for (i=0; i<keys.size(); i++) {
		string s("EMAN.");
		s+=keys[i];
		write_attr(igrp,s.c_str(),keys[i]); 
	}

	H5Gclose(igrp);
	EXITFUNC;
	return 0;
}

int HdfIO2::write_data(float *data, int image_index, const Region* area,
					  EMUtil::EMDataType, bool)
{
	ENTERFUNC;
	if (image_index<0)
		image_index = read_attr(group,"imageid_max");

	char ipath[50];
	sprintf(ipath,"/MDF/images/%d/image",image_index);
	hid_t ds=H5Dopen(file,ipath);
	if (ds<0) throw ImageWriteException(filename,"Image dataset does not exist");
	
	
	EXITFUNC;
	return 0;
}

int HdfIO2::get_nimg()
{
	init();
	return (int)read_attr(group,"imageid_max");
}

void HdfIO2::flush()
{
	return;
}

bool HdfIO2::is_complex_mode()
{
	return false;
}

// always big endian
bool HdfIO2::is_image_big_endian()
{
	return true;
}



#endif	//EM_HDF5
