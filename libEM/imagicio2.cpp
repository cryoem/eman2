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

#include <cstring>
#include "imagicio.h"
#include "portable_fileio.h"
#include "util.h"
#include "geometry.h"
#include "ctf.h"
#include "emassert.h"
#include "transform.h"

#ifdef _WIN32
	#include <ctime>
#endif	//_WIN32

using namespace EMAN;

const char *ImagicIO::HED_EXT = "hed";
const char *ImagicIO::IMG_EXT = "img";
const char *ImagicIO::REAL_TYPE_MAGIC = "REAL";
const char *ImagicIO::CTF_MAGIC = "!-";

ImagicIO2::ImagicIO2(string file, IOMode rw)
:	filename(file), rw_mode(rw), hed_file(0), img_file(0), initialized(false)
{
	hed_filename = Util::change_filename_ext(filename, HED_EXT);
	img_filename = Util::change_filename_ext(filename, IMG_EXT);

	is_big_endian = ByteOrder::is_host_big_endian();
	is_new_hed = false;
	is_new_img = false;
	memset(&imagich, 0, sizeof(Imagic4D));
	imagich.count = -1;
	datatype = IMAGIC_UNKNOWN_TYPE;
	nz = 0;
}

ImagicIO2::~ImagicIO2()
{
	if (hed_file) {
		fclose(hed_file);
		hed_file = 0;
	}

	if (img_file) {
		fclose(img_file);
		img_file = 0;
	}
}

void ImagicIO2::init()
{
	ENTERFUNC;

	if (initialized) {
		return;
	}

	initialized = true;

	is_new_hed = false;
	is_new_img = false;

	hed_file = sfopen(hed_filename, rw_mode, &is_new_hed);
	img_file = sfopen(img_filename, rw_mode, &is_new_img);

	if (is_new_hed != is_new_img) {
		LOGWARN("IMAGIC header file and data file should both exist or both not exist");
	}

	if (!is_new_hed) {
		if (fread(&imagich, sizeof(ImagicHeader), 1, hed_file) != 1) {
			throw ImageReadException(hed_filename, "IMAGIC header");
		}

		if (!is_valid(&imagich)) {
			throw ImageReadException(hed_filename, "invalid IMAGIC file");
		}

		datatype = get_datatype_from_name(imagich.type);

		if (datatype != IMAGIC_USHORT && datatype != IMAGIC_FLOAT) {
			LOGERR("unsupported imagic data type: %s", imagich.type);
			throw ImageReadException(hed_filename, "unsupported imagic data type");
		}

		is_big_endian = ByteOrder::is_data_big_endian(&imagich.ny);
		make_header_host_endian(imagich);
		rewind(hed_file);
	}
	EXITFUNC;
}

bool ImagicIO2::is_valid(const void *first_block)
{
	ENTERFUNC;

	if (!first_block) {
		return false;
	}

	const int *data = static_cast < const int *>(first_block);
	int count = data[1];
	int headrec = data[3];
	int month = data[4];
	int hour = data[7];
	int nx = data[13];
	int ny = data[12];
	int realtype = data[68];
	
	bool data_big_endian = ByteOrder::is_data_big_endian(&headrec);

	if (data_big_endian != ByteOrder::is_host_big_endian()) {
		ByteOrder::swap_bytes(&count);
		ByteOrder::swap_bytes(&headrec);
		ByteOrder::swap_bytes(&month);
		ByteOrder::swap_bytes(&hour);
		ByteOrder::swap_bytes(&nx);
		ByteOrder::swap_bytes(&ny);
		ByteOrder::swap_bytes(&realtype);
	}

	const int max_dim = 1 << 20;
	bool result = false;

	// this field realtype is unique to new Imagic-5 format
	if(realtype != VAX_VMS && realtype != LINUX_WINDOWS && realtype != SGI_IBM) {
		EXITFUNC;
		return result;
	}
	
	if (headrec == 1 &&
		count >= 0 && count < max_dim &&
		nx > 0 && nx < max_dim &&
		ny > 0 && ny < max_dim && 
		month >= 0 && month <=12 && 
		hour >= 0 && hour <= 24) {
		result = true;
	}

	EXITFUNC;
	return result;
}

int ImagicIO2::read_header(Dict & dict, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	check_read_access(image_index);

	int nimg = 1;

	Imagic4D hed;
	if (image_index == 0) {
		hed = imagich;
	}
	else {
		memset(&hed, 0, sizeof(Imagic4D));
		portable_fseek(hed_file, sizeof(Imagic5H) * image_index, SEEK_SET);
		fread(&hed, sizeof(Imagic4D), 1, hed_file);
		make_header_host_endian(hed);
	}
	check_region(area, FloatSize(hed.nx, hed.ny, nimg), is_new_hed, false);

    datatype = get_datatype_from_name(imagich.type);

	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, hed.nx, &xlen, hed.ny, &ylen, hed.izlp, &zlen);

	dict["nx"] = xlen;
	dict["ny"] = ylen;
	dict["nz"] = zlen;

	dict["IMAGIC.imgnum"] = hed.imgnum;
	dict["IMAGIC.count"] = hed.count;
	dict["IMAGIC.error"] = hed.error;
	
	dict["IMAGIC.month"] = hed.month;
	dict["IMAGIC.day"] = hed.mday;
	dict["IMAGIC.year"] = hed.year;
	dict["IMAGIC.hour"] = hed.hour;
	dict["IMAGIC.minute"] = hed.minute;
	dict["IMAGIC.sec"] = hed.sec;

	dict["IMAGIC.reals"] = hed.reals;
	dict["IMAGIC.izold"] = hed.izold;
	
	dict["datatype"] = to_em_datatype(hed.type);

	dict["IMAGIC.ixold"] = hed.ixold;
	dict["IMAGIC.iyold"] = hed.iyold;
	
	dict["mean"] = hed.avdens;
	dict["sigma"] = hed.sigma;
	dict["variance"] = hed.varia;
	dict["IMAGIC.oldav"] = hed.oldav;
	dict["maximum"] = hed.max;
	dict["minimum"] = hed.min;
	
	dict["IMAGIC.cellx"] = hed.cellx;
	dict["IMAGIC.celly"] = hed.celly;
	dict["IMAGIC.cellz"] = hed.cellz;
	dict["IMAGIC.cell_alpha"] = hed.clpha; 
	dict["IMAGIC.cell_beta"] = hed.cbeta;
	dict["IMAGIC.cell_gamma"] = hed.cgamma;

	dict["IMAGIC.label"] = hed.label; 
	dict["IMAGIC.mapc"] = hed.mapc;
	dict["IMAGIC.mapr"] = hed.mapr;
	dict["IMAGIC.maps"] = hed.maps;
	
	dict["IMAGIC.ispg"] = hed.ispg;
	dict["IMAGIC.nxstart"] = hed.nxstart;
	dict["IMAGIC.nystart"] = hed.nystart;
	dict["IMAGIC.nzstart"] = hed.nzstart;
	
	dict["IMAGIC.nxintv"] = hed.nxintv;
	dict["IMAGIC.nyintv"] = hed.nyintv;
	dict["IMAGIC.nzintv"] = hed.nzintv;
	
	dict["IMAGIC.i4lp"] = hed.i4lp;	//number of 3D volumes in 4D data
	
	
	char tmp[5] = { hed.type[0],hed.type[1],hed.type[2],hed.type[3],0 };
	dict["IMAGIC.type"] = tmp;
	

	dict["IMAGIC.oldav"] = hed.oldav;
	
	dict["ptcl_repr"] = hed.mrc2;			// raw images represented by this image

	dict["orientation_convention"] = "EMAN";
    const float alt = hed.mrc1[1]*180.0f/M_PI;
    const float az = hed.mrc1[2]*180.0f/M_PI;
    const float phi = hed.mrc1[0]*180.0f/M_PI;
	dict["euler_alt"] = alt;
	dict["euler_az"] = az;
	dict["euler_phi"] = phi;
	Transform *trans = new Transform();
	trans->set_rotation(Dict("type", "eman", "alt", alt, "az", az, "phi", phi));
	if( hed.count==0 ) {
		dict["xform.projection"] = trans;
	}
	else {
		dict["xform.projection"] = trans;
		dict["xform.align3d"] = trans;
	}
	Ctf * ctf_ = read_ctf(hed);
	if( ctf_ != 0) {
		dict["ctf"] = ctf_;
	}
	if(trans) {delete trans; trans=0;}
	if(ctf_) {delete ctf_; ctf_=0;}
	EXITFUNC;
	return 0;
}

ImagicIO2::DataType ImagicIO2::get_datatype_from_name(const char *name)
{
	DataType t = IMAGIC_UNKNOWN_TYPE;

	if (strncmp(name, "PACK",4) == 0) {
		t = IMAGIC_CHAR;
	}
	else if (strncmp(name, "INTG",4) == 0) {
		t = IMAGIC_SHORT;
	}
	else if (strncmp(name, REAL_TYPE_MAGIC,4) == 0) {
		t = IMAGIC_FLOAT;
	}
	else if (strncmp(name, "COMP",4) == 0) {
		t = IMAGIC_FLOAT_COMPLEX;
	}
	else if (strncmp(name, "RECO",4) == 0) {
		t = IMAGIC_FFT_FLOAT_COMPLEX;
	}
	return t;
}

int ImagicIO2::to_em_datatype(DataType t)
{
	switch (t) {
	case IMAGIC_CHAR:
		return EMUtil::EM_CHAR;
	case IMAGIC_SHORT:
		return EMUtil::EM_SHORT;
	case IMAGIC_FLOAT:
		return EMUtil::EM_FLOAT;
	case IMAGIC_FLOAT_COMPLEX:
	case IMAGIC_FFT_FLOAT_COMPLEX:
		return EMUtil::EM_FLOAT_COMPLEX;
	default:
		break;
	}

	return EMUtil::EM_UNKNOWN;
}

Ctf * ImagicIO2::read_ctf(const Imagic5H& hed) const
{
	ENTERFUNC;

	Ctf * ctf_ = 0;
	size_t n = strlen(CTF_MAGIC);

	if (strncmp(imagich.label, CTF_MAGIC, n) == 0) {
		ctf_ = new EMAN1Ctf();
		string header_label(hed.label);
		// Note: this block was making things crash because it assumed the following if statement
		// was true - I added the if statement (d.woolford)
		if (header_label.size() > 2) {
			string sctf = "O" + header_label.substr(2);
			ctf_->from_string(sctf);
		}
	}

	EXITFUNC;
	return ctf_;
}

void ImagicIO2::write_ctf(const Ctf * const ctf, int image_index)
{
	ENTERFUNC;
	init();

	size_t n = strlen(CTF_MAGIC);
	strcpy(imagich.label, CTF_MAGIC);
	string ctf_ = ctf->to_string().substr(1);
	strncpy(&imagich.label[n], ctf_.c_str(), sizeof(imagich.label) - n);

	rewind(hed_file);
	if (fwrite(&imagich, sizeof(ImagicHeader), 1, hed_file) != 1) {
		throw ImageWriteException(hed_filename, "Imagic Header");
	}

	EXITFUNC;
}