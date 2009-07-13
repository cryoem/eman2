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

#include "amiraio.h"
#include "util.h"

#ifndef WIN32
	#include <sys/param.h>
#else
	#include <windows.h>
	#define MAXPATHLEN (MAX_PATH*4)
#endif	//WIN32

#include <cstdio>

using namespace EMAN;

const char *AmiraIO::MAGIC = "# AmiraMesh";

AmiraIO::AmiraIO(const string & file, IOMode rw)
:	filename(file), rw_mode(rw), amira_file(0),
	is_big_endian(true), initialized(false), dt(EMUtil::EM_UNKNOWN),
	 nx(0), ny(0), nz(0),
	 pixel(0), xorigin(0), yorigin(0), zorigin(0)
{
}

AmiraIO::~AmiraIO()
{
	if (amira_file) {
		fclose(amira_file);
		amira_file = 0;
	}
}

void AmiraIO::init()
{
	ENTERFUNC;

	if (initialized) {
		return;
	}

	initialized = true;
	bool is_new_file = false;

	amira_file = sfopen(filename, rw_mode, &is_new_file, true);

	if (!is_new_file) {
		char buf[MAXPATHLEN];
		if (!fgets(buf, MAXPATHLEN, amira_file)) {
			throw ImageReadException(filename, "Amira Header");
		}

		if (!is_valid(buf)) {
			throw ImageReadException(filename, "invalid Amira Mesh file");
		}

		if(strstr(buf,"BINARY-LITTLE-ENDIAN")!=0) {
			is_big_endian = false;
		}
		else if(strstr(buf,"BINARY")!=0) {
			is_big_endian = true;
		}
		else if(strstr(buf,"2.0")!=0) {
			is_big_endian = true; //AMIRA convention(<=2.0): BIGENDIAN only
		}
	}
	EXITFUNC;
}

bool AmiraIO::is_valid(const void *first_block)
{
	ENTERFUNC;
	bool result = false;
	if (!first_block) {
		result = false;
	}
	else {
		result = Util::check_file_by_magic(first_block, MAGIC);
	}
	EXITFUNC;
	return result;
}

int AmiraIO::read_header(Dict & dict, int, const Region *, bool)
{
	ENTERFUNC;
	init();

	char ll[MAXPATHLEN+10];
	char datatype[16] = "";

	do {
		fgets(ll,MAXPATHLEN,amira_file);
//		printf("%s", ll);
		if(char* s=strstr(ll,"define Lattice ")) {
			if(sscanf(s+15,"%d %d %d",&nx, &ny, &nz) == 3) {
				dict["nx"] = nx;
				dict["ny"] = ny;
				dict["nz"] = nz;
//				printf("nx = %d, ny = %d, nz = %d\n", nx, ny, nz);
			};
		}
		else if(char* s=strstr(ll,"BoundingBoxXY")) {		//amiramesh 2.1
			float bx0, bx1, by0, by1;
			if (sscanf(s+13,"%f %f %f %f",&bx0, &bx1, &by0, &by1) == 4 ) {
				pixel = (bx1-bx0)/(nx-1);
				xorigin = bx0;
				yorigin = by0;
				dict["apix_x"] = pixel;
				dict["apix_y"] = pixel;
				dict["apix_z"] = pixel;
				dict["origin_x"] = xorigin;
				dict["origin_y"] = yorigin;
//				printf("pixel=%g: xorigin=%g yorigin=%g \n", pixel, xorigin, yorigin);
			}
		}
		else if(char* s=strstr(ll,"BoundingBox")) {		//amiramesh 2.0
			float bx0, bx1, by0, by1, bz0, bz1;
			if (sscanf(s+11,"%f %f %f %f %f %f",&bx0, &bx1, &by0, &by1, &bz0, &bz1) == 6 ) {
				pixel = (bx1-bx0)/(nx-1);
				xorigin = bx0;
				yorigin = by0;
				zorigin = bz0;
				dict["apix_x"] = pixel;
				dict["apix_y"] = pixel;
				dict["apix_z"] = pixel;
				dict["origin_x"] = xorigin;
				dict["origin_y"] = yorigin;
				dict["origin_z"] = zorigin;
//				printf("pixel=%g: xorigin=%g yorigin=%g zorigin=%g\n", pixel, xorigin, yorigin, zorigin);
			}
		}
		else if(char* s=strstr(ll,"Lattice { ")) {
			sscanf(s+10,"%s ", datatype);
			if(!strncmp(datatype, "float", 5)) {
				dt = EMUtil::EM_FLOAT;
				dict["datatype"] = EMUtil::EM_FLOAT;
			}
			else if(!strncmp(datatype, "short", 5)) {
				dt = EMUtil::EM_SHORT;
				dict["datatype"] = EMUtil::EM_SHORT;
			}
			else if(!strncmp(datatype, "byte", 4)) {
				dt = EMUtil::EM_CHAR;
				dict["datatype"] = EMUtil::EM_CHAR;
			}
			else {
				fprintf(stderr,"AmiraIO::read_header: data type \"%s\" is not supported yet\n", datatype);
				return -1;
			}
		}
	} while (! ( ll[0]=='@' && ll[1]=='1') );

	EXITFUNC;
	return 0;

}

int AmiraIO::write_header(const Dict & dict, int image_index, const Region*, EMUtil::EMDataType, bool)
{
	ENTERFUNC;
	int err = 0;

	check_write_access(rw_mode, image_index, 1);

	nx = dict["nx"];
	ny = dict["ny"];
	nz = dict["nz"];

	float xorigin = 0.0f;
	if(dict.has_key("origin_x")) xorigin = dict["origin_x"];
	float yorigin = 0.0f;
	if(dict.has_key("origin_y")) yorigin = dict["origin_y"];
	float zorigin = 0.0f;
	if(dict.has_key("origin_z")) zorigin = dict["origin_z"];
	float pixel = 0.0f;
	if(dict.has_key("apix_x")) pixel = dict["apix_x"];

	rewind(amira_file);

	string line1;
	if(ByteOrder::is_host_big_endian()) {
		line1= "# AmiraMesh 3D BINARY 2.1\n\n";
	}
	else {
		line1= "# AmiraMesh BINARY-LITTLE-ENDIAN 2.1\n\n";
	}

	string type = "float";

	if (fprintf(amira_file,line1.c_str()) <= 0) {
		LOGERR("cannot write to AmiraMesh file '%s'", filename.c_str());
		err = 1;
	}
	else {
		fprintf(amira_file,"define Lattice %d %d %d\n\n",nx,ny,nz);
		fprintf(amira_file,"Parameters {\n");
		fprintf(amira_file, "\tContent \"%dx%dx%d %s, uniform coordinates\",\n", nx,ny,nz,type.c_str());
		fprintf(amira_file, "\tCoordType \"uniform\",\n");
		fprintf(amira_file,"\tBoundingBox %.2f %.2f %.2f %.2f %.2f %.2f\n}\n\n",xorigin,xorigin+pixel*(nx-1),yorigin,yorigin+pixel*(ny-1),zorigin,zorigin+pixel*(nz-1));

		fprintf(amira_file, "Lattice { float ScalarField } @1\n\n# Data section follows\n@1\n");
	}

	EXITFUNC;
	return err;
}

int AmiraIO::read_data(float * rdata, int, const Region *, bool)
{
	ENTERFUNC;

	size_t size = nx*ny*nz;
	switch(dt) {
	case EMUtil::EM_FLOAT:
		fread(rdata,nx*nz,ny*sizeof(float),amira_file);
		if( (is_big_endian && ByteOrder::is_host_big_endian()) && !(is_big_endian || ByteOrder::is_host_big_endian()) ) {
			char tmpdata;
			char *cdata=(char*)rdata;
			for(size_t i=0;i<size;++i){
				//swap 0-3
				tmpdata=cdata[4*i+3];
				cdata[4*i+3]=cdata[4*i];
				cdata[4*i] = tmpdata;
				//swap 1-2
				tmpdata=cdata[4*i+2];
				cdata[4*i+2]=cdata[4*i+1];
				cdata[4*i+1] = tmpdata;
			}
		}
		break;
	case EMUtil::EM_SHORT:
	{
		short *datashort = (short*)malloc(sizeof(short)*nx*ny*nz);
		fread(datashort,nx*nz,ny*sizeof(short),amira_file);
		if( (is_big_endian && ByteOrder::is_host_big_endian()) && !(is_big_endian || ByteOrder::is_host_big_endian()) ) {
			char tmpdata;
			char *cdata=(char*)datashort;
			for(size_t i=0;i<size;i++){
				//swap 0-1
				tmpdata=cdata[2*i+1];
				cdata[2*i+1]=cdata[2*i];
				cdata[2*i] = tmpdata;
			}
			for(size_t i=0; i<size; i++) rdata[i] = float(datashort[i]);
			free(datashort);
		}
	}
		break;
	case EMUtil::EM_CHAR:
	{
		char *databyte = (char*)malloc(sizeof(char)*nx*ny*nz);
		fread(databyte,nx*nz,ny*sizeof(char),amira_file);
		for(size_t i=0; i<size; i++) rdata[i] = float(databyte[i]);
		free(databyte);
	}
		break;
	default:
		fprintf(stderr,"AmiraIO::read_data: data type is not supported yet\n");
		return -1;
	}

	EXITFUNC;
	return 0;
}

int AmiraIO::write_data(float *data, int image_index, const Region*, EMUtil::EMDataType, bool)
{
	ENTERFUNC;

	check_write_access(rw_mode, image_index, 1, data);
//	ByteOrder::become_big_endian(data, nx * ny * nz);

	if (fwrite(data, nx * nz, ny * sizeof(float), amira_file) != ny * sizeof(float)) {
		throw ImageWriteException(filename, "incomplete file write in AmiraMesh file");
	}

	EXITFUNC;
	return 0;
}

void AmiraIO::flush()
{
	fflush(amira_file);
}

bool AmiraIO::is_complex_mode()
{
	return false;
}

bool AmiraIO::is_image_big_endian()
{
	init();
	return is_big_endian;
}

