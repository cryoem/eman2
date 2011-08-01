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
#include <climits>
#include "imagicio2.h"
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

const char *ImagicIO2::HED_EXT = "hed";
const char *ImagicIO2::IMG_EXT = "img";
const char *ImagicIO2::REAL_TYPE_MAGIC = "REAL";
const char *ImagicIO2::CTF_MAGIC = "!-";

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
		if (fread(&imagich, sizeof(Imagic4D), 1, hed_file) != 1) {
			throw ImageReadException(hed_filename, "IMAGIC4D header");
		}

//		if (!is_valid(&imagich)) {
//			throw ImageReadException(hed_filename, "invalid IMAGIC file");
//		}

		datatype = get_datatype_from_name(imagich.type);

		if (datatype != IMAGIC_SHORT && datatype != IMAGIC_FLOAT) {
			LOGERR("unsupported imagic data type: %s", imagich.type);
			throw ImageReadException(hed_filename, "unsupported imagic data type");
		}

		is_big_endian = ByteOrder::is_data_big_endian(&imagich.ny);
		make_header_host_endian(imagich);
		rewind(hed_file);
	}

	EXITFUNC;
}

int ImagicIO2::init_test()
{
	ENTERFUNC;

	if (initialized) {
		return 1;
	}

	FILE *in = fopen(hed_filename.c_str(), "rb");
	if (!in) {
		throw FileAccessException(filename);
	}

	char first_block[1024];
	size_t n = fread(first_block, sizeof(char), sizeof(first_block), in);

	if (n == 0) {
		LOGERR("file '%s' is an empty file", filename.c_str());
		fclose(in);
		return -1;
	}
	fclose(in);

	const int *data = reinterpret_cast <const int *>(first_block);
	int nx = data[13];
	int ny = data[12];
	int izold = data[11];

	if(izold==nx*ny) {
		EXITFUNC;
		return -1;	//old style IMAGIC file
	}
	else {
		EXITFUNC;
		return 0;	//new IMAGIC4D file
	}
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
	int hour = data[7];
	int minute = data[8];
	int second = data[9];
	int rsize = data[10];
	int nx = data[13];
	int ny = data[12];
	int nz = data[60];
	int realtype = data[68];

	bool data_big_endian = ByteOrder::is_data_big_endian(&headrec);

	if (data_big_endian != ByteOrder::is_host_big_endian()) {
		ByteOrder::swap_bytes(&count);
		ByteOrder::swap_bytes(&headrec);
		ByteOrder::swap_bytes(&hour);
		ByteOrder::swap_bytes(&rsize);
		ByteOrder::swap_bytes(&nx);
		ByteOrder::swap_bytes(&ny);
		ByteOrder::swap_bytes(&nz);
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
		nz > 0 && nz < max_dim &&
		hour >= 0 && hour < 24 &&
		minute >=0 && minute <60 &&
		second >=0 && second <60) {
		result = true;
	}

	EXITFUNC;
	return result;
}

int ImagicIO2::read_header(Dict & dict, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	check_read_access(image_index);

	Imagic4D hed;
	if (image_index == 0) {
		hed = imagich;
	}
	else {
		memset(&hed, 0, sizeof(Imagic4D));
		portable_fseek(hed_file, sizeof(Imagic4D) * image_index, SEEK_SET);
		fread(&hed, sizeof(Imagic4D), 1, hed_file);
		make_header_host_endian(hed);
	}

	int nz = hed.izlp ? hed.izlp : 1;
	check_region(area, FloatSize(hed.nx, hed.ny, nz), is_new_hed, false);

    datatype = get_datatype_from_name(imagich.type);

	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, hed.nx, &xlen, hed.ny, &ylen, nz, &zlen);

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

	dict["IMAGIC.rsize"] = hed.rsize;
	dict["IMAGIC.izold"] = hed.izold;
	
	dict["datatype"] = to_em_datatype(datatype);
	dict["IMAGIC.type"] = hed.type+'\0';

	dict["IMAGIC.ixold"] = hed.ixold;
	dict["IMAGIC.iyold"] = hed.iyold;
	
	dict["mean"] = hed.avdens;
	dict["sigma"] = hed.sigma;

	dict["maximum"] = hed.densmax;
	dict["minimum"] = hed.densmin;
	
	dict["IMAGIC.complex"] = hed.complex;
	dict["IMAGIC.defocus1"] = hed.defocus1;
	dict["IMAGIC.defocus2"] = hed.defocus2;
	dict["IMAGIC.defangle"] = hed.defangle;
	dict["IMAGIC.sinostart"] = hed.sinostart;
	dict["IMAGIC.sinoend"] = hed.sinoend;

	dict["IMAGIC.label"] = hed.label+'\0';

	dict["IMAGIC.ccc3d"] = hed.ccc3d;
	dict["IMAGIC.ref3d"] = hed.ref3d;
	dict["IMAGIC.mident"] = hed.mident;
	dict["IMAGIC.ezshift"] = hed.ezshift;
	
	dict["IMAGIC.ealpha"] = hed.ealpha;
	dict["IMAGIC.ebeta"] = hed.ebeta;
	dict["IMAGIC.egamma"] = hed.egamma;
	
	dict["IMAGIC.nalisum"] = hed.nalisum;
	dict["IMAGIC.pgroup"] = hed.pgroup;
	
	dict["IMAGIC.i4lp"] = hed.i4lp;	//number of objects (1D, 2D, or 3D) in 4D data
	
	dict["IMAGIC.alpha"] = hed.alpha;
	dict["IMAGIC.beta"] = hed.beta;
	dict["IMAGIC.gamma"] = hed.gamma;

	dict["IMAGIC.IMAVERS"] = hed.imavers;
	dict["IMAGIC.REALTYPE"] = hed.realtype;

	dict["IMAGIC.ANGLE"] = hed.angle;
	dict["IMAGIC.VOLTAGE"] = hed.voltage;
	dict["IMAGIC.SPABERR"] = hed.spaberr;
	dict["IMAGIC.PCOHER"] = hed.pcoher;
	dict["IMAGIC.CCC"] = hed.ccc;
	dict["IMAGIC.ERRAR"] = hed.errar;
	dict["IMAGIC.ERR3D"] = hed.err3d;
	dict["IMAGIC.REF"] = hed.ref;
	dict["IMAGIC.CLASSNO"] = hed.ref;
	dict["IMAGIC.LOCOLD"] = hed.locold;
	dict["IMAGIC.REPQUAL"] = hed.repqual;
	dict["IMAGIC.ZSHIFT"] = hed.zshift;
	dict["IMAGIC.XSHIFT"] = hed.xshift;
	dict["IMAGIC.YSHIFT"] = hed.yshift;
	dict["IMAGIC.NUMCLS"] = hed.numcls;
	dict["IMAGIC.OVQUAL"] = hed.ovqual;
	dict["IMAGIC.EANGLE"] = hed.eangle;
	dict["IMAGIC.EXSHIFT"] = hed.exshift;
	dict["IMAGIC.EYSHIFT"] = hed.eyshift;
	dict["IMAGIC.CMTOTVAR"] = hed.cmtotvar;
	dict["IMAGIC.INFORMAT"] = hed.informat;
	dict["IMAGIC.NUMEIGEN"] = hed.numeigen;
	dict["IMAGIC.NIACTIVE"] = hed.niactive;
	dict["IMAGIC.RESOLX"] = hed.resolx;
	dict["IMAGIC.RESOLY"] = hed.resoly;
	dict["IMAGIC.RESOLZ"] = hed.resolz;
	if(hed.errar==-1.0) {
		dict["IMAGIC.FABOSA1"] = hed.alpha2;
		dict["IMAGIC.FABOSA2"] = hed.beta2;
		dict["IMAGIC.FABOSA3"] = hed.gamma2;
	}
	else {
		dict["IMAGIC.ALPHA2"] = hed.alpha2;
		dict["IMAGIC.BETA2"] = hed.beta2;
		dict["IMAGIC.GAMMA2"] = hed.gamma2;
	}
	dict["IMAGIC.NMETRIC"] = hed.nmetric;
	dict["IMAGIC.ACTMSA"] = hed.actmsa;

	vector<float> v_coosmsa(hed.coosmsa, hed.coosmsa+69);
	dict["IMAGIC.COOSMSA"] = v_coosmsa;

	dict["IMAGIC.EIGVAL"] = hed.coosmsa[19];
	dict["IMAGIC.HISTORY"] = hed.history+'\0';

	dict["orientation_convention"] = "IMAGIC";
	const float alpha = hed.alpha;
	const float beta = hed.beta;
	const float gamma = hed.gamma;
	dict["euler_alpha"] = alpha;
	dict["euler_beta"] = beta;
	dict["euler_gamma"] = gamma;
	Transform *trans = new Transform();
	trans->set_rotation(Dict("type", "imagic", "alpha", alpha, "beta", beta, "gamma", gamma));
	if( nz<=1 ) {
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

ImagicIO2::DataType ImagicIO2::get_datatype_from_name(const char *name) const
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

int ImagicIO2::to_em_datatype(DataType t) const
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

Ctf * ImagicIO2::read_ctf(const Imagic4D& hed) const
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
	if (fwrite(&imagich, sizeof(Imagic4D), 1, hed_file) != 1) {
		throw ImageWriteException(hed_filename, "Imagic Header");
	}

	EXITFUNC;
}

bool ImagicIO2::is_complex_mode()
{
	init();
	if (datatype == IMAGIC_FLOAT_COMPLEX || datatype == IMAGIC_FFT_FLOAT_COMPLEX) {
		return true;
	}
	return false;
}

void ImagicIO2::flush()
{
	fflush(img_file);
	fflush(hed_file);
}

int ImagicIO2::write_header(EMAN::Dict const& dict, int image_index,
		EMAN::Region const* area, EMUtil::EMDataType, bool use_host_endian)
{
	ENTERFUNC;

	if(image_index<0) {
		image_index = get_nimg();
	}
	check_write_access(rw_mode, image_index);

	if (area) {
		check_region(area, FloatSize(imagich.nx, imagich.ny, imagich.izlp),
					 is_new_hed);
		EXITFUNC;
		return 0;
	}

	int nx = dict["nx"];
	int ny = dict["ny"];
	int nz = dict["nz"];
	int nimg=0;		//# images currently in file

	if (!is_new_hed) {
        datatype = get_datatype_from_name(imagich.type);

		if (imagich.nx != nx || imagich.ny != ny || imagich.izlp != nz) {
			char desc[256];
			sprintf(desc, "new IMAGIC size %dx%dx%d is not equal to existing size %dx%dx%d",
					nx, ny, nz, imagich.nx, imagich.ny, imagich.izlp);
			throw ImageWriteException(filename, desc);
		}

        if (datatype!=IMAGIC_FLOAT) {
			throw ImageWriteException(filename, "Attempted write to non REAL Imagic file");
		}

        rewind(hed_file);
		nimg=image_index+1;
	}
	else {
		nimg = 1;	//new file writing
	}

	Imagic4D new_hed;
	memset(&new_hed, 0, sizeof(Imagic4D));

	time_t cur_time = time(0);
	struct tm *tm = localtime(&cur_time);

	new_hed.error = 0;
	new_hed.headrec = 1;	//always 1 a the moment

	new_hed.mday = tm->tm_mday;
	new_hed.month = tm->tm_mon;
	new_hed.year = tm->tm_year + 1900;
	new_hed.hour = tm->tm_hour;
	new_hed.minute = tm->tm_min;
	new_hed.sec = tm->tm_sec;

	size_t img_size = nx*ny*nz;
	if(img_size > (size_t)INT_MAX) {
		new_hed.rsize = -1;
	}
	else {
		new_hed.rsize = (int)img_size;
	}

	new_hed.nx = nx;
	new_hed.ny = ny;
	new_hed.izlp = nz;

	strncpy(new_hed.type, REAL_TYPE_MAGIC,4);
	new_hed.avdens = (float)dict["mean"];
	new_hed.sigma = (float)dict["sigma"];
	new_hed.densmax = (float)dict["maximum"];
	new_hed.densmin = (float)dict["minimum"];

	new_hed.ixold = 0;
	new_hed.iyold = 0;

	string new_label = dict.has_key("IMAGIC.label") ? (string) dict["IMAGIC.label"] : "";
	sprintf(new_hed.label, new_label.c_str() );

	new_hed.i4lp = nimg;

	Transform * t = 0;
	if(nz<=1 && dict.has_key("xform.projection")) {
		t = dict["xform.projection"];
	}
	else if(nz>1 && dict.has_key("xform.align3d")) {
		t = dict["xform.align3d"];
	}

	if(t) {
		Dict d = t->get_rotation("imagic");
		new_hed.alpha = d["alpha"];
		new_hed.beta = d["beta"];
		new_hed.gamma = d["gamma"];
		delete t;
		t=0;
	}
	else {
		if(dict.has_key("euler_alpha")) new_hed.alpha = dict["euler_alpha"];
		if(dict.has_key("euler_beta")) new_hed.beta = dict["euler_beta"];
		if(dict.has_key("euler_gamma")) new_hed.gamma = dict["euler_gamma"];
	}

	new_hed.realtype = generate_machine_stamp();

	new_hed.resolx = dict["apix_x"];
	new_hed.resoly = dict["apix_y"];
	new_hed.resolz = dict["apix_z"];

	if ( (is_big_endian != ByteOrder::is_host_big_endian()) || !use_host_endian)  swap_header(new_hed);

	// overwrite existing header if necessary
	if (image_index>=0 && image_index<nimg) {
		portable_fseek(hed_file, sizeof(Imagic4D)*image_index, SEEK_SET);
		new_hed.imgnum=image_index+1;
		if (is_big_endian != ByteOrder::is_host_big_endian())
				ByteOrder::swap_bytes((int *) &new_hed.imgnum,1);
		fwrite(&new_hed, sizeof(Imagic4D),1,hed_file);
	}

	// update the 1st header with total # images
	int ifol = nimg-1;
	if (is_big_endian != ByteOrder::is_host_big_endian()) {
		ByteOrder::swap_bytes((int *) &nimg,1);
		ByteOrder::swap_bytes((int *) &ifol,1);
	}
	portable_fseek(hed_file, sizeof(int), SEEK_SET);
	fwrite(&ifol, sizeof(int), 1, hed_file);
	portable_fseek(hed_file, 61*sizeof(int), SEEK_SET);	//I4LP(62) is the number of "objects" in file
	fwrite(&nimg, sizeof(int), 1, hed_file);

	// header in machine order
	if ( (is_big_endian != ByteOrder::is_host_big_endian()) || !use_host_endian)  swap_header(new_hed);
	imagich=new_hed;
	imagich.count=nimg;
	is_new_hed = false;

	if( dict.has_key("ctf") ) {
		Ctf * ctf_ = dict["ctf"];
		write_ctf(ctf_);
		if(ctf_) {delete ctf_; ctf_=0;}
	}

	EXITFUNC;
	return 0;
}

int ImagicIO2::generate_machine_stamp() const
{
	int machinestamp;

#ifdef __sgi
	machinestamp = SGI_IBM;
#elif defined __OPENVMS__
	machinestamp = VAX_VMS;
#else
	machinestamp = LINUX_WINDOWS;
#endif

	return machinestamp;
}

void ImagicIO2::make_header_host_endian(Imagic4D& hed) const
{
	if (is_big_endian != ByteOrder::is_host_big_endian()) {
		swap_header(hed);
	}
}

void ImagicIO2::swap_header(Imagic4D & hed) const
{
	ByteOrder::swap_bytes((int *) &hed, NUM_4BYTES_PRE_IYLP);
	ByteOrder::swap_bytes(&hed.ixold, NUM_4BYTES_AFTER_IXOLD);
	ByteOrder::swap_bytes((int *) &hed.ccc3d, NUM_4BYTES_AFTER_NAME);
}

int ImagicIO2::write_data(float* data, int image_index, const Region * area, EMAN::EMUtil::EMDataType, bool use_host_endian)
{
	ENTERFUNC;

	check_write_access(rw_mode, image_index, 0, data);
	check_region(area, FloatSize(imagich.nx, imagich.ny, imagich.izlp), is_new_hed);

	if (image_index == -1) {
		portable_fseek(img_file, 0, SEEK_END);
	}
	else {
		size_t img_size = imagich.nx * imagich.ny * imagich.izlp * sizeof(float);
		portable_fseek(img_file, img_size*image_index, SEEK_SET);
	}

	if(is_new_img) {
		if(!use_host_endian) {
			ByteOrder::swap_bytes(data, imagich.nx * imagich.ny * imagich.izlp);
		}
	}
	else if (is_big_endian != ByteOrder::is_host_big_endian()) {
		ByteOrder::swap_bytes(data, imagich.nx * imagich.ny * imagich.izlp);
	}

	EMUtil::process_region_io(data, img_file, WRITE_ONLY, 0,
							  sizeof(float), imagich.nx, imagich.ny,
							  imagich.izlp, area, true);

	EXITFUNC;
	return 0;
}

int ImagicIO2::read_data(float* data, int image_index, EMAN::Region const* area, bool)
{
	ENTERFUNC;

	check_read_access(image_index, data);
	Assert(datatype != IMAGIC_UNKNOWN_TYPE);

	int nx = imagich.ny;
	int ny = imagich.nx;
	int nz = imagich.izlp ? imagich.izlp : 1;
	size_t img_size = (size_t)nx*ny*nz;

	check_region(area, FloatSize(nx, ny, nz), is_new_hed, false);

	portable_fseek(img_file, img_size*image_index*sizeof(float), SEEK_SET);

	short *sdata = (short *) data;
	unsigned char *cdata = (unsigned char *) data;
	size_t mode_size = get_datatype_size(datatype);

	/**The image_index option does not work in EMUtil::process_region_io(), so I set the
	 * file pointer to the right place before actually read data*/
	EMUtil::process_region_io(cdata, img_file, READ_ONLY, 0, mode_size, nx, ny, nz, area, true);

	if (datatype == IMAGIC_FLOAT) {
		become_host_endian(data, img_size);
	}
	else if (datatype == IMAGIC_SHORT) {
		become_host_endian(sdata, img_size);

		for (ptrdiff_t j = img_size - 1; j >= 0; --j) {
			data[j] = static_cast < float >(sdata[j]);
		}
	}
	else {
		throw ImageReadException(filename, "unknown imagic data type");
	}

	EXITFUNC;
	return 0;
}

int ImagicIO2::get_nimg()
{
	init();
	return imagich.count + 1;
}

bool ImagicIO2::is_image_big_endian()
{
	init();
	return is_big_endian;
}

size_t ImagicIO2::get_datatype_size(DataType t) const
{
	size_t s = 0;
	switch (t) {
	case IMAGIC_CHAR:
		s = sizeof(unsigned char);
		break;
	case IMAGIC_SHORT:
		s = sizeof(unsigned short);
		break;
	case IMAGIC_FLOAT:
	case IMAGIC_FLOAT_COMPLEX:
	case IMAGIC_FFT_FLOAT_COMPLEX:
		s = sizeof(float);
		break;
	default:
		s = 0;
	}

	return s;
}
