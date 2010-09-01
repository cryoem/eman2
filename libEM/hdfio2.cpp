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

//#define DEBUGHDF	1

#include "hdfio2.h"
#include "geometry.h"
#include "ctf.h"
#include "emassert.h"
#include "transform.h"
#include "ctf.h"
#include <iostream>
#include <cstring>

#ifndef WIN32
	#include <sys/param.h>
#else
	#define  MAXPATHLEN (MAX_PATH * 4)
#endif	//WIN32

using namespace EMAN;

static const int ATTR_NAME_LEN = 128;

HdfIO2::HdfIO2(const string & hdf_filename, IOMode rw)
:	nx(1), ny(1), nz(1), is_exist(false),
	file(-1), group(-1), filename(hdf_filename),
	rw_mode(rw), initialized(false)
{
	accprop=H5Pcreate(H5P_FILE_ACCESS);

	//STDIO file driver has 2G size limit on 32 bit Linux system
	H5Pset_fapl_sec2( accprop );
	//H5Pset_fapl_stdio( accprop );

//	H5Pset_fapl_core( accprop, 1048576, 0  );
//	H5Pset_cache(accprop)
	hsize_t dims=1;
	simple_space=H5Screate_simple(1,&dims,NULL);
}

HdfIO2::~HdfIO2()
{
	H5Sclose(simple_space);
	H5Pclose(accprop);
    if (group >= 0) {
        H5Gclose(group);
    }
    if (file >= 0) {
		H5Fflush(file,H5F_SCOPE_GLOBAL);	// If there were no resource leaks, this wouldn't be necessary...
		H5Fclose(file);
    }
#ifdef DEBUGHDF
	printf("HDf close\n");
#endif
}

// This reads an already opened attribute and returns the results as an EMObject
// The attribute is not closed
EMObject HdfIO2::read_attr(hid_t attr) {
	hid_t type = H5Aget_type(attr);
	hid_t spc = H5Aget_space(attr);
	H5T_class_t cls = H5Tget_class(type);
	size_t sz = H5Tget_size(type);						// storage size, arrays handled in the 'space'
	hssize_t pts = H5Sget_simple_extent_npoints(spc);	// number of points > 1 if an array of floats or integers

	EMObject ret(0);
	char c;
	int i;
//	unsigned int ui;
	float f,*fa;
	int * ia;
//	unsigned int * uia;
	double d;
	char *s;
	vector <float> fv((size_t)pts);
	vector <int> iv((size_t)pts);
//	vector <unsigned int> uiv(pts);

	float *matrix;
	Transform* t;
	Ctf* ctf;
// 	int r, c, k=0;

	switch (cls) {
	case H5T_INTEGER:
		if(sz==1) {
			H5Aread(attr,H5T_NATIVE_CHAR,&c);
			bool b = false;
			if(c=='T') {
				b = true;
			}
			else if(c=='F') {
				b = false;
			}
			ret = EMObject(b);
		}
		else if(sz==4) {
			if(pts==1) {
				H5Aread(attr,H5T_NATIVE_INT,&i);
				ret=EMObject(i);
			}
			else {
				ia=(int *)malloc((size_t)pts*sizeof(int));
				H5Aread(attr,H5T_NATIVE_INT,ia);
				for (i=0; i<pts; i++) iv[i]=ia[i];
				free(ia);
				ret=EMObject(iv);
			}
		}
		break;
// 	case H5T_UNSIGNED_INTEGER:
// 		if(pts==1) {
// 			H5Aread(attr,H5T_NATIVE_UINT,&ui);
// 			ret=EMObject(ui);
// 		}
// 		else {
// 			uia=(unsigned int *)malloc(pts*sizeof(unsigned int));
// 			H5Aread(attr,H5T_NATIVE_UINT,uia);
// 			for (i=0; i<pts; i++) uiv[i]=uia[i];
// 			free(uia);
// 			ret=EMObject(uiv);
// 		}
// 		break;
	case H5T_FLOAT:
		if (sz==4) {
			if (pts==1) {
				H5Aread(attr,H5T_NATIVE_FLOAT,&f);
				ret=EMObject(f);
			}
			else {
				fa=(float *)malloc((size_t)pts*sizeof(float));
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
		H5Aread(attr,type,s);
//		H5Aread(attr,H5T_NATIVE_CHAR,s);
		if(s[0] == 'O' && isdigit(s[1])) {
			ctf = new EMAN1Ctf();
			ctf->from_string(string(s));
			ret = EMObject(ctf);
			delete ctf;
		}
		else if(s[0] == 'E' && isdigit(s[1])) {
			ctf = new EMAN2Ctf();
			ctf->from_string(string(s));
			ret = EMObject(ctf);
			delete ctf;
		}
		else {
			ret=EMObject(s);
		}
		free(s);
		break;
	case H5T_COMPOUND:
		matrix = (float*)malloc(12*sizeof(float));
		H5Aread(attr, type, matrix);
//		ret.create_transform3d_by_array(trans3d);
		t = new Transform(matrix);
		ret = EMObject(t);
		free(matrix);
		delete t; t=0;

// 		trans3d = (float*)malloc(16*sizeof(float));	//16 float for a Transform3D object
// 		H5Aread(attr, type, trans3d);
// //		ret.create_transform3d_by_array(trans3d);
// 		trans = new Transform3D();
// 		for(r=0; r<4; ++r) {
// 			for(c=0; c<4; ++c) {
// 				trans->set(r, c, trans3d[k]);
// 				++k;
// 			}
// 		}
// 		ret = EMObject(trans);
// 		free(trans3d);
		break;
	default:
		LOGERR("Unhandled HDF5 metadata %d", cls);
	}

	H5Sclose(spc);
	H5Tclose(type);

	return ret;
}

// This writes an attribute with specified name to a given open object
// The attribute is opened and closed. returns 0 on success
int HdfIO2::write_attr(hid_t loc,const char *name,EMObject obj) {
	hid_t type=0;
	hid_t spc=0;
	hsize_t dims=1;
	vector <float> fv;
	vector <int> iv;
	switch(obj.get_type())
	{
	case EMObject::BOOL:
		type=H5Tcopy(H5T_NATIVE_CHAR);
		spc=H5Scopy(simple_space);
		break;
	case EMObject::INT:
		type=H5Tcopy(H5T_NATIVE_INT);
		spc=H5Scopy(simple_space);
		break;
	case EMObject::UNSIGNEDINT:
		type=H5Tcopy(H5T_NATIVE_UINT);
		spc=H5Scopy(simple_space);
		break;
	case EMObject::FLOAT:
		type=H5Tcopy(H5T_NATIVE_FLOAT);
		spc=H5Scopy(simple_space);
		break;
	case EMObject::DOUBLE:
		type=H5Tcopy(H5T_NATIVE_DOUBLE);
		spc=H5Scopy(simple_space);
		break;
	case EMObject::STRING:
	case EMObject::CTF:
		type=H5Tcopy(H5T_C_S1);
		H5Tset_size(type,strlen((const char *)obj)+1);
		spc=H5Screate(H5S_SCALAR);
		break;
	case EMObject::FLOATARRAY:
		type=H5Tcopy(H5T_NATIVE_FLOAT);
		fv=obj;
		dims=fv.size();
		spc=H5Screate_simple(1,&dims,NULL);
		break;
	case EMObject::INTARRAY:
		type=H5Tcopy(H5T_NATIVE_INT);
		iv=obj;
		dims=iv.size();
		spc=H5Screate_simple(1,&dims,NULL);
		break;
	case EMObject::TRANSFORM:
		type = H5Tcreate(H5T_COMPOUND, 12 * sizeof(float)); //Transform is a 3x4 matrix
		H5Tinsert(type, "00", 0, H5T_NATIVE_FLOAT);
		H5Tinsert(type, "01", 1*sizeof(float), H5T_NATIVE_FLOAT);
		H5Tinsert(type, "02", 2*sizeof(float), H5T_NATIVE_FLOAT);
		H5Tinsert(type, "03", 3*sizeof(float), H5T_NATIVE_FLOAT);
		H5Tinsert(type, "10", 4*sizeof(float), H5T_NATIVE_FLOAT);
		H5Tinsert(type, "11", 5*sizeof(float), H5T_NATIVE_FLOAT);
		H5Tinsert(type, "12", 6*sizeof(float), H5T_NATIVE_FLOAT);
		H5Tinsert(type, "13", 7*sizeof(float), H5T_NATIVE_FLOAT);
		H5Tinsert(type, "20", 8*sizeof(float), H5T_NATIVE_FLOAT);
		H5Tinsert(type, "21", 9*sizeof(float), H5T_NATIVE_FLOAT);
		H5Tinsert(type, "22", 10*sizeof(float), H5T_NATIVE_FLOAT);
		H5Tinsert(type, "23", 11*sizeof(float), H5T_NATIVE_FLOAT);
		H5Tpack(type);

		dims = 1;	//one compound type
		spc = H5Screate_simple(1, &dims, NULL);
		break;
	case EMObject::STRINGARRAY:
	case EMObject::EMDATA:
	case EMObject::XYDATA:
	case EMObject::FLOAT_POINTER:
	case EMObject::INT_POINTER:
	case EMObject::VOID_POINTER:
		return -1;
		break;
	case EMObject::UNKNOWN:
		break;
	}

    //we need this delete attribute call here, even we called erase_header()
    //at the biginning of write_header(), since the  "imageid_max" need be updated correctly.
	if( H5Adelete(loc,name) < 0 ) {
#ifdef DEBUGHDF
		LOGERR("Attribute %s deletion error in write_attr().\n", name);
#endif
	}
	else {
#ifdef DEBUGHDF
		printf("delete attribute %s successfully in write_attr().\n", name);
#endif
	}
	hid_t attr = H5Acreate(loc,name,type,spc,H5P_DEFAULT);

	bool b;
	char c;
	int i;
	float f,*fa;
	int * ia;
	unsigned int ui;
	double d;
	const char *s;
	Transform * tp;
	switch(obj.get_type()) {
	case EMObject::BOOL:
		b = (bool)obj;
		if(b) {
			c = 'T';
		} else {
			c = 'F';
		}
		H5Awrite(attr,type,&c);
		break;
	case EMObject::INT:
		i=(int)obj;
		H5Awrite(attr,type,&i);
		break;
	case EMObject::UNSIGNEDINT:
		ui=(unsigned int)obj;
		H5Awrite(attr,type,&ui);
		break;
	case EMObject::FLOAT:
		f=(float)obj;
		H5Awrite(attr,type,&f);
		break;
	case EMObject::DOUBLE:
		d=(double)obj;
		H5Awrite(attr,type,&d);
		break;
	case EMObject::STRING:
	case EMObject::CTF:
		s=(const char *)obj;
		H5Awrite(attr,type,s);
		break;
	case EMObject::FLOATARRAY:
		fa=(float *)malloc(fv.size()*sizeof(float));
		for (ui=0; ui<fv.size(); ui++) fa[ui]=fv[ui];
		H5Awrite(attr,type,fa);
		free(fa);
		break;
	case EMObject::INTARRAY:
		ia=(int *)malloc(iv.size()*sizeof(int));
		for (ui=0; ui<iv.size(); ui++) ia[ui]=iv[ui];
		H5Awrite(attr,type,ia);
		free(ia);
		break;
	case EMObject::TRANSFORM:
	{
		tp = (Transform *)obj;
		fa = (float *)malloc(12*sizeof(float));
		int r, c, k=0;
		for(r=0; r<3; ++r) {
			for(c=0; c<4; ++c) {
				fa[k] = tp->at(r,c);
				++k;
			}
		}
		H5Awrite(attr,type,fa);
		free(fa);
	}
		break;
//	case EMObject::STRINGARRAY:
//	case EMObject::EMDATA:
//	case EMObject::XYDATA:
//		return -1;
//		break;
	default:
		LOGERR("Unhandled HDF5 metadata '%s'", name);

	}

	herr_t ret1 = H5Tclose(type);
	herr_t ret2 = H5Sclose(spc);
	herr_t ret3 = H5Aclose(attr);
	if(ret1>=0 && ret2>=0 && ret3>=0) {
		return 0;
	}
	else {
		LOGERR("close error in write_attr()\n");
		return -1;
	}
}

// Initializes the file for read-only or read-write access
// Data is stored under /MDF/images
// An attribute named imageid_max stores the number of the highest
// numbered image in the file.
// A group is then made for each individual image, all metadata for the
// individual images is currently associated with the GROUP, not the dataset
// dataset-specific data could also be associated with the dataset in
// future. At the moment, there is only a single dataset in each group.
void HdfIO2::init()
{
	ENTERFUNC;
	if (initialized) {
		return;
	}
#ifdef DEBUGHDF
	printf("init\n");
#endif

	H5Eset_auto(0, 0);	// Turn off console error logging.

	if (rw_mode == READ_ONLY) {
		file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, accprop);
		if (file<0) throw FileAccessException(filename);
	}
	else {
		file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, accprop);
		if (file < 0) {
			file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, accprop);
			if (file < 0) {
				throw FileAccessException(filename);
			}
			else {
#ifdef DEBUGHDF
				printf("File truncated or new file created\n");
#endif
			}
		}
	}

	group=H5Gopen(file,"/MDF/images");
	if (group<0) {
		if (rw_mode == READ_ONLY) throw ImageReadException(filename,"HDF5 file has no image data (no /MDF group)");
		group=H5Gcreate(file,"/MDF",64);		// create the group for Macromolecular data
		if (group<0) throw ImageWriteException(filename,"Unable to add image group (/MDF) to HDF5 file");
		H5Gclose(group);
		group=H5Gcreate(file,"/MDF/images",4096);		// create the group for images/volumes
		if (group<0) throw ImageWriteException(filename,"Unable to add image group (/MDF/images) to HDF5 file");
		write_attr(group,"imageid_max",EMObject(-1));
	}
	initialized = true;
	EXITFUNC;
}


// If this version of init() returns -1, then we have an old-style HDF5 file
int HdfIO2::init_test()
{
	ENTERFUNC;
	if (initialized) {
		return 1;
	}
#ifdef DEBUGHDF
	printf("init_test\n");
#endif

	H5Eset_auto(0, 0);	// Turn off console error logging.

	hid_t fileid = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5Pcreate(H5P_FILE_ACCESS));
	hid_t groupid = H5Gopen(fileid, "/");
	hid_t attid = H5Aopen_name(groupid, "num_dataset");

	if (attid < 0) {
		H5Gclose(groupid);
		H5Fclose(fileid);
		init();
		EXITFUNC;
		return 0;
	}
	else {
		H5Aclose(attid);
		H5Gclose(groupid);
		H5Fclose(fileid);
		EXITFUNC;
		return -1;
	}
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

// Reads all of the attributes from the /MDF/images/<imgno> group
int HdfIO2::read_header(Dict & dict, int image_index, const Region * area, bool)
{
	ENTERFUNC;
	init();
#ifdef DEBUGHDF
	printf("read_head %d\n", image_index);
#endif
	int i;
	// Each image is in a group for later expansion. Open the group
	char ipath[50];
	sprintf(ipath,"/MDF/images/%d", image_index);
	hid_t igrp=H5Gopen(file, ipath);

	int nattr=H5Aget_num_attrs(igrp);

	char name[ATTR_NAME_LEN];
	for (i=0; i<nattr; i++) {
		hid_t attr=H5Aopen_idx(igrp, i);
		H5Aget_name(attr,127,name);
		if (strncmp(name,"EMAN.", 5)!=0) {
			H5Aclose(attr);
			continue;
		}
		EMObject val=read_attr(attr);
		dict[name+5]=val;
		H5Aclose(attr);
	}

	if(dict.has_key("ctf")) {
		string ctfString = (string)dict["ctf"];
		if(ctfString.substr(0, 1) == "O") {
			Ctf * ctf_ = new EMAN1Ctf();
			ctf_->from_string(ctfString);
			dict.erase("ctf");
			dict["ctf"] = ctf_;
			delete ctf_;
		}
		else if(ctfString.substr(0, 1) == "E") {
			Ctf * ctf_ = new EMAN2Ctf();
			ctf_->from_string(ctfString);
			dict.erase("ctf");
			dict["ctf"] = ctf_;
			delete ctf_;
		}
	}

	if(area) {
		check_region(area, IntSize(dict["nx"], dict["ny"], dict["nz"]), false, false);

		dict["nx"] = area->get_width();
		dict["ny"] = area->get_height();
		dict["nz"] = area->get_depth();

		if( dict.has_key("apix_x") && dict.has_key("apix_y") && dict.has_key("apix_z") )
		{
			if( dict.has_key("origin_x") && dict.has_key("origin_y") && dict.has_key("origin_z") )
			{
				float xorigin = dict["origin_x"];
				float yorigin = dict["origin_y"];
				float zorigin = dict["origin_z"];

				float apix_x = dict["apix_x"];
				float apix_y = dict["apix_y"];
				float apix_z = dict["apix_z"];

				dict["origin_x"] = xorigin + apix_x * area->origin[0];
				dict["origin_y"] = yorigin + apix_y * area->origin[1];
				dict["origin_z"] = zorigin + apix_z * area->origin[2];
			}
		}
	}

	H5Gclose(igrp);

	//Get the data type from data set, HDF5 file header attribute 'datatype' may be wrong
	sprintf(ipath,"/MDF/images/%d/image",image_index);
	hid_t ds=H5Dopen(file,ipath);

	if(ds>0) {	//ds>0 means successfully open the dataset
		hid_t dt = H5Dget_type(ds);

		switch(H5Tget_size(dt)) {
		case 4:
			dict["datatype"] = (int)EMUtil::EM_FLOAT;
			break;
		case 2:
			dict["datatype"] = (int)EMUtil::EM_USHORT;
			break;
		case 1:
			dict["datatype"] = (int)EMUtil::EM_UCHAR;
			break;
		default:
			throw ImageReadException(filename, "EMAN does not support this data type.");
		}

		H5Tclose(dt);
	}

	H5Dclose(ds);

	EXITFUNC;
	return 0;
}

// This erases any existing attributes from the image group
// prior to writing a new header. For a new image there
// won't be any, so this should be harmless.
int HdfIO2::erase_header(int image_index)
{
	ENTERFUNC;

	if(image_index < 0) return 0; //image_index<0 for appending image, no need for erasing

	init();
#ifdef DEBUGHDF
	printf("erase_head %d\n",image_index);
#endif
	int i;
	// Each image is in a group for later expansion. Open the group
	char ipath[50];
	sprintf(ipath,"/MDF/images/%d", image_index);
	hid_t igrp=H5Gopen(file, ipath);

	int nattr=H5Aget_num_attrs(igrp);

	char name[ATTR_NAME_LEN];
	for (i=0; i<nattr; i++) {
		hid_t attr = H5Aopen_idx(igrp, 0); //use 0 as index here, since the H5Adelete() will change the index
		H5Aget_name(attr,127,name);
		H5Aclose(attr);
		if( H5Adelete(igrp,name) < 0 ) {
			LOGERR("attribute %s deletion error in erase_header().\n", name);
		}
	}

	H5Gclose(igrp);
	EXITFUNC;
	return 0;
}


int HdfIO2::read_data(float *data, int image_index, const Region *area, bool)
{
	ENTERFUNC;
#ifdef DEBUGHDF
	printf("read_data %d\n",image_index);
#endif

	char ipath[50];
	sprintf(ipath,"/MDF/images/%d/image",image_index);
	hid_t ds=H5Dopen(file,ipath);
	if (ds<0) throw ImageWriteException(filename,"Image does not exist");
	hid_t spc=H5Dget_space(ds);
	hid_t dt = H5Dget_type(ds);

	hsize_t dims_out[3];
	hsize_t rank = H5Sget_simple_extent_ndims(spc);

	H5Sget_simple_extent_dims(spc, dims_out, NULL);
	if(rank == 1) {
		nx = dims_out[0];
		ny = 1;
		nz = 1;
	}
	else if(rank == 2) {
		nx = dims_out[1];
		ny = dims_out[0];
		nz = 1;
	}
	else if(rank == 3) {
		nx = dims_out[2];
		ny = dims_out[1];
		nz = dims_out[0];
	}

	if (area) {
		hid_t memoryspace = 0;

 		/*Get the file dataspace - the region we want to read in the file*/
		int x0 = 0, y0 = 0, z0 = 0;		//the coordinates for up left corner, trim to be within image bound
		int x1 = 0, y1 = 0, z1 = 0;		//the coordinates for down right corner, trim to be within image bound
		int nx1 = 1, ny1 = 1, nz1 = 1;	//dimensions of the sub-region, actual region read form file
 		if(rank == 3) {
			hsize_t     doffset[3];             /* hyperslab offset in the file */
			doffset[2] = (hsize_t)(area->x_origin() < 0 ? 0 : area->x_origin());
			doffset[1] = (hsize_t)(area->y_origin() < 0 ? 0 : area->y_origin());
			doffset[0] = (hsize_t)(area->z_origin() < 0 ? 0 : area->z_origin());
			x0 = (int)doffset[0];
			y0 = (int)doffset[1];
			z0 = (int)doffset[2];

			z1 = (int)(area->x_origin() + area->get_width());
			z1 = (int)(z1 > static_cast<int>(nx) ? nx : z1);

			y1 = (int)(area->y_origin() + area->get_height());
			y1 = (int)(y1 > static_cast<int>(ny) ? ny : y1);
			if(y1 <= 0) {
				y1 = 1;
			}

			x1 = (int)(area->z_origin() + area->get_depth());
			x1 = (int)(x1 > static_cast<int>(nz) ? nz : x1);
			if(x1 <= 0) {
				x1 = 1;
			}

			if(x1 < x0 || y1< y0 || z1 < z0) return 0; //out of bounds, this is fine, nothing happens

			hsize_t     dcount[3];              /* size of the hyperslab in the file */
			dcount[0] = x1 - doffset[0];
			dcount[1] = y1 - doffset[1];
			dcount[2] = z1 - doffset[2];

			H5Sselect_hyperslab (spc, H5S_SELECT_SET, (const hsize_t*)doffset, NULL, (const hsize_t*)dcount, NULL);

			/*Define memory dataspace - the memory we will created for the region*/
			hsize_t     dims[3];              /* size of the region in the memory */
			dims[0] = dcount[2]?dcount[2]:1;
			dims[1]	= dcount[1]?dcount[1]:1;
			dims[2] = dcount[0]?dcount[0]:1;
			nx1 = (int)dims[0];
			ny1 = (int)dims[1];
			nz1 = (int)dims[2];

			memoryspace = H5Screate_simple(3, dims, NULL);
 		}
 		else if(rank == 2) {
 			hsize_t     doffset[2];             /* hyperslab offset in the file */
			doffset[1] = (hsize_t)(area->x_origin() < 0 ? 0 : area->x_origin());
			doffset[0] = (hsize_t)(area->y_origin() < 0 ? 0 : area->y_origin());
			x0 = (int)doffset[0];
			y0 = (int)doffset[1];
			z0 = 1;

			y1 = (int)(area->x_origin() + area->get_width());
			y1 = (int)(y1 > static_cast<int>(nx) ? nx : y1);

			x1 = (int)(area->y_origin() + area->get_height());
			x1 = (int)(x1 > static_cast<int>(ny) ? ny : x1);
			if(x1 <= 0) {
				x1 = 1;
			}

			z1 = 1;

			if(x1 < x0 || y1< y0) return 0; //out of bounds, this is fine, nothing happens

			hsize_t     dcount[2];              /* size of the hyperslab in the file */
			dcount[0] = x1 - doffset[0];
			dcount[1] = y1 - doffset[1];

			H5Sselect_none(spc);
			H5Sselect_hyperslab (spc, H5S_SELECT_SET, (const hsize_t*)doffset, NULL, (const hsize_t*)dcount, NULL);

			/*Define memory dataspace - the memory we will created for the region*/
			hsize_t     dims[2];              /* size of the region in the memory */
			dims[0] = (hsize_t)(dcount[1]?dcount[1]:1);
			dims[1]	= (hsize_t)(dcount[0]?dcount[0]:1);
			nx1 = (int)dims[0];
			ny1 = (int)dims[1];
			nz1 = 1;

			memoryspace = H5Screate_simple(2, dims, NULL);
 		}

 		if( (area->x_origin()>=0) && (area->y_origin()>=0) && (area->z_origin()>=0) && ((hsize_t)(area->x_origin() + area->get_width())<=nx) && ((hsize_t)(area->y_origin() + area->get_height())<=ny) && ((hsize_t)(area->z_origin() + area->get_depth())<=nz) ){	//the region is in boundary
 			H5Dread(ds,H5T_NATIVE_FLOAT,memoryspace,spc,H5P_DEFAULT,data);
 		}
 		else {	//the region are partial out of boundary
 			/* When the requested region has some part out of image boundary,
 			 * we need read the sub-area which is within image,
 			 * and fill the out of boundary part with zero.
 			 * We actually read the sub-region from HDF by hyperslab I/O,
 			 * then copy it back to the pre-allocated region.*/
 			float * subdata = new float[nx1*ny1*nz1];


 			H5Dread(ds,H5T_NATIVE_FLOAT,memoryspace,spc,H5P_DEFAULT,subdata);

 			int xd0=0, yd0=0, zd0=0;	//The coordinates of the top-left corner sub-region in region
 			size_t clipped_row_size = 0;
 			if(rank == 3) {
				xd0 = (int) (area->x_origin() < 0 ? -area->x_origin() : 0);
				yd0 = (int) (area->y_origin() < 0 ? -area->y_origin() : 0);
				zd0 = (int) (area->z_origin() < 0 ? -area->z_origin() : 0);
 				clipped_row_size = (z1-z0)* sizeof(float);
 			}
 			else if(rank == 2) {
 				xd0 = (int) (area->x_origin() < 0 ? -area->x_origin() : 0);
 				yd0 = (int) (area->y_origin() < 0 ? -area->y_origin() : 0);
 				clipped_row_size = (y1-y0)* sizeof(float);
 			}

 			int src_secsize = nx1 * ny1;
 			int dst_secsize = (int)(area->get_width())*(int)(area->get_height());

 			float * src = subdata;
 			float * dst = data + zd0*dst_secsize + yd0*(int)(area->get_width()) + xd0;

 			int src_gap = src_secsize - (y1-y0) * nx1;
 			int dst_gap = dst_secsize - (y1-y0) * (int)(area->get_width());

 			for(int i = 0; i<nz1; ++i) {
 				for(int j = 0; j<ny1; ++j) {
 					EMUtil::em_memcpy(dst, src, clipped_row_size);

 					src += nx1;
 					dst += (int)(area->get_width());
 				}
 				src += src_gap;
 				dst += dst_gap;
 			}

 			delete [] subdata;
 		}
 		H5Sclose(memoryspace);
	} else {
		hsize_t size = nx*ny*nz;
		size_t i=0;
		size_t j=0;
		unsigned short *usdata = (unsigned short *) data;
		unsigned char *cdata = (unsigned char *) data;
		switch(H5Tget_size(dt)) {
		case 4:
			H5Dread(ds,H5T_NATIVE_FLOAT,spc,spc,H5P_DEFAULT,data);
			break;
		case 2:
			H5Dread(ds,H5T_NATIVE_USHORT,spc,spc,H5P_DEFAULT,usdata);
			for (i = 0; i < size; ++i) {
				j = size - 1 - i;
				data[j] = static_cast < float >(usdata[j]);
			}
			break;
		case 1:
			H5Dread(ds,H5T_NATIVE_UCHAR,spc,spc,H5P_DEFAULT,cdata);
			for (i = 0; i < size; ++i) {
				j = size - 1 - i;
				data[j] = static_cast < float >(cdata[j]);
			}
			break;
		default:
			throw ImageReadException(filename, "EMAN does not support this data type.");
		}
	}

	H5Tclose(dt);
	H5Sclose(spc);
	H5Dclose(ds);
	EXITFUNC;
	return 0;
}


// Writes all attributes in 'dict' to the image group
// Creation of the image dataset is also handled here
int HdfIO2::write_header(const Dict & dict, int image_index, const Region* area,
						EMUtil::EMDataType dt, bool)
{
#ifdef DEBUGHDF
	printf("write_head %d\n",image_index);
#endif
	ENTERFUNC;
	init();

	nx = (int)dict["nx"];
	ny = (int)dict["ny"];
	nz = (int)dict["nz"];

	if(image_index<0) {
		image_index = get_nimg();
	}

	// If image_index<0 append, and make sure the max value in the file is correct
	// though this is normally handled by EMData.write_image()
	hid_t attr=H5Aopen_name(group,"imageid_max");
	int nimg = read_attr(attr);
	H5Aclose(attr);

	unsigned int i;
	if (image_index<0) image_index=nimg+1;
	if (image_index>nimg) {
		write_attr(group,(const char *)"imageid_max",EMObject(image_index));
	}

	// Each image is in a group for later expansion. Open the group, create if necessary
	char ipath[50];
	sprintf(ipath,"/MDF/images/%d",image_index);
	hid_t igrp=H5Gopen(file,ipath);

	if (igrp<0) {	//group not existed
		is_exist = false;
		// Need to create a new image group
		igrp=H5Gcreate(file,ipath,64);		// The image is a group, with attributes on the group
		if (igrp<0) throw ImageWriteException(filename,"Unable to add /MDF/images/# to HDF5 file");
	}
	//if group already existed, erase the header and unlink the existing dataset first
	else {
		is_exist = true;
		int nattr=H5Aget_num_attrs(igrp);
		char name[ATTR_NAME_LEN];
		Dict dict2;
		for (int i=0; i<nattr; i++) {
			hid_t attr=H5Aopen_idx(igrp, i);
			H5Aget_name(attr,127,name);
			if (strncmp(name,"EMAN.", 5)!=0) {
				H5Aclose(attr);
				continue;
			}
			EMObject val=read_attr(attr);
			dict2[name+5]=val;
			H5Aclose(attr);

			if(!dict2.has_key("datatype")) {//by default, HDF5 is written as float
				dict2["datatype"] = (int)EMUtil::EM_FLOAT;
			}
		}

		erase_header(image_index);

		//change the size or data type of a image,
		//the existing data set is invalid, unlink it
		if( (int)dict["nx"]*(int)dict["ny"]*(int)dict["nz"] !=
			(int)dict2["nx"]*(int)dict2["ny"]*(int)dict2["nz"] ||
			dict["datatype"] != dict2["datatype"] ) {
			sprintf(ipath,"/MDF/images/%d/image",image_index);
			H5Gunlink(igrp, ipath);
		}
	}

	if(area) {
		check_region(area, IntSize(dict["nx"], dict["ny"], dict["nz"]), false, true);
	}

	// Write the attributes to the group
	vector <string> keys=dict.keys();

	for (i=0; i<keys.size(); i++) {
		string s("EMAN.");
		s+=keys[i];
		write_attr(igrp,s.c_str(),dict[keys[i]]);
	}

	H5Gclose(igrp);
	EXITFUNC;
	return 0;
}

// Writes the actual image data to the corresponding dataset (already created)
int HdfIO2::write_data(float *data, int image_index, const Region* area,
					  EMUtil::EMDataType dt, bool)
{
	ENTERFUNC;

#ifdef DEBUGHDF
	printf("write_data %d\n",image_index);
#endif

	if (image_index<0) {
		hid_t attr=H5Aopen_name(group,"imageid_max");
		image_index = read_attr(attr);
		H5Aclose(attr);
	}

	hid_t spc;	//dataspace
	hid_t ds;	//dataset
	char ipath[50];
	sprintf(ipath,"/MDF/images/%d/image",image_index);

	// Now create the actual image dataset
	if (nz==1)  {
		hsize_t dims[2]= { ny,nx };
		spc=H5Screate_simple(2,dims,NULL);
	}
	else {
		hsize_t dims[3]= { nz, ny, nx };
		spc=H5Screate_simple(3,dims,NULL);
	}

	ds=H5Dopen(file,ipath);
	if(ds<0) {//new dataset
		switch(dt) {
		case EMUtil::EM_FLOAT:
			ds=H5Dcreate(file,ipath, H5T_NATIVE_FLOAT, spc, H5P_DEFAULT );
			break;
		case EMUtil::EM_USHORT:
			ds=H5Dcreate(file,ipath, H5T_NATIVE_USHORT, spc, H5P_DEFAULT );
			break;
		case EMUtil::EM_UCHAR:
			ds=H5Dcreate(file,ipath, H5T_NATIVE_UCHAR, spc, H5P_DEFAULT );
			break;
		default:
			throw ImageWriteException(filename,"HDF5 does not support this data format");
		}
	}

	//convert data to unsigned short, unsigned char...
	size_t size = nx*ny*nz;
	unsigned char *cdata = 0;
	unsigned short *usdata = 0;
	float rendermin = 0.0f;
	float rendermax = 0.0f;
	EMUtil::getRenderMinMax(data, nx, ny, rendermin, rendermax, nz);

	if(area) {
		hsize_t doffset[3];		/*hyperslab offset in the file*/
		doffset[0] = (hsize_t)(area->x_origin());
		doffset[1] = (hsize_t)(area->y_origin());
		doffset[2] = (hsize_t)(area->z_origin());

		hsize_t dcount[3];		/*size of the hyperslab in the file*/
		dcount[0] = (hsize_t)(area->get_width());
		dcount[1] = (hsize_t)(area->get_height()?area->get_height():1);
		dcount[2] = (hsize_t)(area->get_depth()?area->get_depth():1);

		H5Sselect_hyperslab(spc, H5S_SELECT_SET, (const hsize_t*)doffset, NULL, (const hsize_t*)dcount, NULL);

		/*Create memory space with size of the region.*/
		hsize_t dims[3];	/*size of the region in the memory*/
		dims[0] = (hsize_t)(area->get_width());
		dims[1] = (hsize_t)(area->get_height()?area->get_height():1);
		dims[2] = (hsize_t)(area->get_depth()?area->get_depth():1);

		hid_t memoryspace = H5Screate_simple(3, dims, NULL);
		switch(dt) {
		case EMUtil::EM_FLOAT:
			H5Dwrite(ds, H5T_NATIVE_FLOAT, memoryspace, spc, H5P_DEFAULT, data);
			break;
		default:
			throw ImageWriteException(filename,"HDF5 does not support regional writing for this data format");
		}
		H5Sclose(memoryspace);
	}
	else {
		switch(dt) {
		case EMUtil::EM_FLOAT:
			H5Dwrite(ds,H5T_NATIVE_FLOAT,spc,spc,H5P_DEFAULT,data);
			break;
		case EMUtil::EM_USHORT:
			usdata = new unsigned short[size];
			for (size_t i = 0; i < size; ++i) {
				if(data[i] <= rendermin) {
					usdata[i] = 0;
				}
				else if(data[i] >= rendermax) {
					usdata[i] = USHRT_MAX;
				}
				else {
					usdata[i]=(unsigned short)((data[i]-rendermin)/(rendermax-rendermin)*USHRT_MAX);
				}
			}
			H5Dwrite(ds,H5T_NATIVE_USHORT,spc,spc,H5P_DEFAULT,usdata);
			if(usdata) {delete [] usdata; usdata=0;}
			break;
		case EMUtil::EM_UCHAR:
			cdata = new unsigned char[size];
			for (size_t i = 0; i < size; ++i) {
				if(data[i] <= rendermin) {
					cdata[i] = 0;
				}
				else if(data[i] >= rendermax){
					cdata[i] = UCHAR_MAX;
				}
				else {
					cdata[i]=(unsigned char)((data[i]-rendermin)/(rendermax-rendermin)*UCHAR_MAX);
				}
			}
			H5Dwrite(ds,H5T_NATIVE_UCHAR,spc,spc,H5P_DEFAULT,cdata);
			if(cdata) {delete [] cdata; cdata=0;}
			break;
		default:
			throw ImageWriteException(filename,"HDF5 does not support this data format");
		}
	}

	H5Sclose(spc);
	H5Dclose(ds);
	EXITFUNC;
	return 0;
}

int HdfIO2::get_nimg()
{
	init();
	hid_t attr=H5Aopen_name(group,"imageid_max");
	int n = read_attr(attr);
	H5Aclose(attr);

	return n+1;
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
