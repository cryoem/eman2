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

ImagicIO::ImagicIO(string file, IOMode rw)
:	filename(file), rw_mode(rw), hed_file(0), img_file(0), initialized(false)
{
	hed_filename = Util::change_filename_ext(filename, HED_EXT);
	img_filename = Util::change_filename_ext(filename, IMG_EXT);

	is_big_endian = ByteOrder::is_host_big_endian();
	is_new_hed = false;
	is_new_img = false;
	memset(&imagich, 0, sizeof(ImagicHeader));
	imagich.count = -1;
	datatype = IMAGIC_UNKNOWN_TYPE;
	nz = 0;
}

ImagicIO::~ImagicIO()
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

void ImagicIO::init()
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

bool ImagicIO::is_valid(const void *first_block)
{
	ENTERFUNC;

	if (!first_block) {
		return false;
	}

	const int *data = static_cast < const int *>(first_block);
	int count = data[1];
	int headrec = data[3];
	int month = data[5];
	int hour = data[7];
	int nx = data[13];
	int ny = data[12];

	bool data_big_endian = ByteOrder::is_data_big_endian(&headrec);

	if (data_big_endian != ByteOrder::is_host_big_endian()) {
		ByteOrder::swap_bytes(&count);
		ByteOrder::swap_bytes(&headrec);
		ByteOrder::swap_bytes(&month);
		ByteOrder::swap_bytes(&hour);
		ByteOrder::swap_bytes(&nx);
		ByteOrder::swap_bytes(&ny);
	}

	const int max_dim = 1 << 20;
	bool result = false;

	if (headrec == 1 &&
		count >= 0 && count < max_dim &&
		nx > 0 && nx < max_dim &&
		ny > 0 && ny < max_dim && month >= 0 && month <= 12 && hour >= 0 && hour <= 24) {
		result = true;
	}

	EXITFUNC;
	return result;
}

int ImagicIO::read_header(Dict & dict, int image_index, const Region * area, bool is_3d)
{
	ENTERFUNC;

	check_read_access(image_index);

	int nimg = 1;

	if (is_3d) {
		nimg = imagich.count + 1;

		if (nimg <= 1) {
			LOGWARN("this is not a 3D IMAGIC. Read as a 2D");
		}
	}

	ImagicHeader hed;
	if (image_index == 0) {
		hed = imagich;
	}
	else {
		memset(&hed, 0, sizeof(ImagicHeader));
		portable_fseek(hed_file, sizeof(ImagicHeader) * image_index, SEEK_SET);
		fread(&hed, sizeof(ImagicHeader), 1, hed_file);
		make_header_host_endian(hed);
	}
	check_region(area, FloatSize(hed.nx, hed.ny, nimg), is_new_hed,false);

    datatype = get_datatype_from_name(imagich.type);

	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, hed.nx, &xlen, hed.ny, &ylen, nimg, &zlen);

	dict["nx"] = xlen;
	dict["ny"] = ylen;
	dict["nz"] = zlen;

	dict["datatype"] = to_em_datatype(datatype);

	dict["minimum"] = hed.min;
	dict["maximum"] = hed.max;
	dict["mean"] = hed.avdens;
	dict["sigma"] = hed.sigma;

	dict["IMAGIC.imgnum"] = hed.imgnum;
	dict["IMAGIC.count"] = hed.count;
	dict["IMAGIC.error"] = hed.error;

	dict["IMAGIC.headrec"] = hed.headrec;
	dict["IMAGIC.mday"] = hed.mday;
	dict["IMAGIC.month"] = hed.month;

	dict["IMAGIC.year"] = hed.year;
	dict["IMAGIC.hour"] = hed.hour;
	dict["IMAGIC.minute"] = hed.minute;

	dict["IMAGIC.sec"] = hed.sec;
	dict["IMAGIC.reals"] = hed.reals;
	dict["IMAGIC.pixels"] = hed.pixels;

	char tmp[5] = { hed.type[0],hed.type[1],hed.type[2],hed.type[3],0 };
	dict["IMAGIC.type"] = tmp;
	dict["IMAGIC.ixold"] = hed.ixold;
	dict["IMAGIC.iyold"] = hed.iyold;

	dict["IMAGIC.oldav"] = hed.oldav;
	dict["IMAGIC.label"] = hed.label;
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
	EXITFUNC;
	return 0;
}

int ImagicIO::write_header(const Dict & dict, int image_index,
						   const Region * area, EMUtil::EMDataType, bool use_host_endian)
{
	ENTERFUNC;

	check_write_access(rw_mode, image_index);
	nz = dict["nz"];
	if (nz > 1 && image_index != 0) {
		throw ImageWriteException(filename, "to write 3D IMAGIC image, image index must be 0");
	}

	if (area) {
		check_region(area, FloatSize(imagich.nx, imagich.ny, imagich.count+1),
					 is_new_hed);
		EXITFUNC;
		return 0;
	}

	int nx = dict["nx"];
	int ny = dict["ny"];
	int nimg=0;		//# images currently in file


	if (!is_new_hed) {
        datatype = get_datatype_from_name(imagich.type);

		if (imagich.nx != nx || imagich.ny != ny) {
			char desc[256];
			sprintf(desc, "new IMAGIC size %dx%d is not equal to existing size %dx%d",
					nx, ny, imagich.nx, imagich.ny);
			throw ImageWriteException(filename, desc);
		}

        if (datatype!=IMAGIC_FLOAT) {
			throw ImageWriteException(filename, "Attempted write to non REAL Imagic file");
		}

        rewind(hed_file);
		nimg=imagich.count+1;
	}

	ImagicHeader new_hed;
	memset(&new_hed, 0, sizeof(ImagicHeader));

	time_t cur_time = time(0);
	struct tm *tm = localtime(&cur_time);

	new_hed.error = 0;
	new_hed.headrec = 1;

	new_hed.mday = tm->tm_mday;
	new_hed.month = tm->tm_mon;
	new_hed.year = tm->tm_year + 1900;
	new_hed.hour = tm->tm_hour;
	new_hed.minute = tm->tm_min;
	new_hed.sec = tm->tm_sec;

	new_hed.reals = nx * ny;
	new_hed.pixels = nx * ny;
	new_hed.ny = ny;
	new_hed.nx = nx;

	new_hed.ixold = 0;
	new_hed.iyold = 0;
	new_hed.oldav = 0;

	new_hed.min = (float)dict["minimum"];
	new_hed.max = (float)dict["maximum"];
	new_hed.avdens = (float)dict["mean"];
	new_hed.sigma = (float)dict["sigma"];

	if(nz<=1 && dict.has_key("xform.projection")) {
		Transform * t = dict["xform.projection"];
		Dict d = t->get_rotation("eman");
		new_hed.mrc1[1] = (float)d["alt"]*M_PI/180.0f;
		new_hed.mrc1[2] = (float)d["az"]*M_PI/180.0f;
		new_hed.mrc1[0] = (float)d["phi"]*M_PI/180.0f;
	}
	else if(nz>1 && dict.has_key("xform.align3d")) {
		Transform * t = dict["xform.align3d"];
		Dict d = t->get_rotation("eman");
		new_hed.mrc1[1] = (float)d["alt"]*M_PI/180.0f;
		new_hed.mrc1[2] = (float)d["az"]*M_PI/180.0f;
		new_hed.mrc1[0] = (float)d["phi"]*M_PI/180.0f;
	}
	else {
		if(dict.has_key("euler_alt")) new_hed.mrc1[1] = (float)dict["euler_alt"]*M_PI/180.0f;
		if(dict.has_key("euler_az")) new_hed.mrc1[2] = (float)dict["euler_az"]*M_PI/180.0f;
		if(dict.has_key("euler_phi")) new_hed.mrc1[0] = (float)dict["euler_phi"]*M_PI/180.0f;
	}

	if(dict.has_key("ptcl_repr")) new_hed.mrc2 = (int)dict["ptcl_repr"];

	string new_label = dict.has_key("IMAGIC.label") ? (string) dict["IMAGIC.label"] : "";
	sprintf(new_hed.label, new_label.c_str() );

	new_hed.lbuf = nx;
	new_hed.inn = 1;
	new_hed.iblp = ny;
	new_hed.ifb = 0;
	new_hed.lbw = 0;
	new_hed.lbr = -1;
	new_hed.lastlr = -1;
	new_hed.lastlw = 1;
	new_hed.num = 8;
	new_hed.nhalf = nx / 2;
	new_hed.ibsd = nx * 2;
	new_hed.ihfl = 7;
	new_hed.lcbr = -1;
	new_hed.lcbw = 1;
	new_hed.imstr = -1;
	new_hed.imstw = -1;
	new_hed.istart = 1;
	new_hed.iend = nx;
	new_hed.leff = nx;
	new_hed.linbuf = nx * 2;
	new_hed.ntotbuf = -1;
	new_hed.icstart = 1;
	new_hed.icend = nx / 2;
	strncpy(new_hed.type, REAL_TYPE_MAGIC,4);


	// header in file order
	if ( (is_big_endian != ByteOrder::is_host_big_endian()) || !use_host_endian)  swap_header(new_hed);

	// overwrite existing header if necessary
	if (image_index>=0 && image_index<nimg) {
		portable_fseek(hed_file, sizeof(ImagicHeader)*image_index, SEEK_SET);
		new_hed.imgnum=image_index+1;
		if (is_big_endian != ByteOrder::is_host_big_endian())
				ByteOrder::swap_bytes((int *) &new_hed.imgnum,1);
		fwrite(&new_hed, sizeof(ImagicHeader),1,hed_file);
	}

	// How many images does the file need when we're done ?
	int required_len;
	if (nz>1) required_len=nz;
	else {
		if (image_index<0) required_len=nimg+1;
		else if (image_index+1>nimg) required_len=image_index+1;
		else required_len=nimg;
	}

	// Extend the file to the necessary length
	portable_fseek(hed_file, 0, SEEK_END);
	while (nimg<required_len) {
		nimg++;
		new_hed.imgnum=nimg;
		if (is_big_endian != ByteOrder::is_host_big_endian())
				ByteOrder::swap_bytes((int *) &new_hed.imgnum,1);
		fwrite(&new_hed, sizeof(ImagicHeader),1,hed_file);
	}

	// update the 1st header with total # images
	portable_fseek(hed_file, sizeof(int), SEEK_SET);
	nimg--;
	if (is_big_endian != ByteOrder::is_host_big_endian())
			ByteOrder::swap_bytes((int *) &nimg,1);
	fwrite(&nimg, sizeof(int), 1, hed_file);

	// header in machine order
	if ( (is_big_endian != ByteOrder::is_host_big_endian()) || !use_host_endian)  swap_header(new_hed);
	imagich=new_hed;
	imagich.count=nimg;
	is_new_hed = false;

	if( dict.has_key("ctf") ) {
		Ctf * ctf_ = dict["ctf"];
		write_ctf(ctf_);
	}

	EXITFUNC;
	return 0;
}

int ImagicIO::read_data(float *data, int image_index, const Region * area, bool is_3d)
{
	ENTERFUNC;

	check_read_access(image_index, data);
    Assert(datatype != IMAGIC_UNKNOWN_TYPE);
	int nimg = 1;
	if (is_3d) {
		nimg = imagich.count + 1;
	}

	if (is_3d && imagich.count < 1) {
		LOGWARN("this is not a 3D IMAGIC. Read as a 2D");
		is_3d = false;
	}

	check_region(area, FloatSize(imagich.nx, imagich.ny, nimg), is_new_hed,false);

	rewind(img_file);

	unsigned short *sdata = (unsigned short *) data;
	unsigned char *cdata = (unsigned char *) data;
	size_t mode_size = get_datatype_size(datatype);
	EMUtil::process_region_io(cdata, img_file, READ_ONLY, image_index, mode_size,
							  imagich.nx, imagich.ny, nimg, area, true);

#if 0
	int row_size = imagich.nx * mode_size;
	int sec_size = imagich.nx * imagich.ny * mode_size;

	for (int k = 0; k < nimg; k++) {
		for (int i = imagich.ny - 1; i >= 0; i--) {
			if (fread(&cdata[k * sec_size + i * row_size], row_size, 1, img_file) != 1) {
				LOGERR("incomplete data read: %d/%d blocks on file '%s'",
									 i, imagich.ny, filename.c_str());
				return 1;
			}
		}
	}
#endif

	int img_size = imagich.nx * imagich.ny * nimg;

	if (datatype == IMAGIC_FLOAT) {
		become_host_endian(data, img_size);
	}
	else if (datatype == IMAGIC_USHORT) {
		become_host_endian((unsigned short *) cdata, img_size);

		for (int j = img_size - 1; j >= 0; j--) {
			data[j] = static_cast < float >(sdata[j]);
		}
	}
	else {
		throw ImageReadException(filename, "unknown imagic data type");
	}

	EXITFUNC;
	return 0;
}

int ImagicIO::write_data(float *data, int image_index, const Region* area,
						 EMUtil::EMDataType, bool use_host_endian)
{
	ENTERFUNC;

	check_write_access(rw_mode, image_index, 0, data);
	check_region(area, FloatSize(imagich.nx, imagich.ny, imagich.count+1),
				 is_new_hed);

	if (nz == 1) {
		if (image_index == -1) {
			portable_fseek(img_file, 0, SEEK_END);
		}
		else {
			size_t sec_size = imagich.nx * imagich.ny * sizeof(float);
			portable_fseek(img_file, ((off_t) sec_size) * image_index, SEEK_SET);
		}
	}

	if(is_new_img) {
		if(!use_host_endian) {
			ByteOrder::swap_bytes(data, imagich.nx * imagich.ny * nz);
		}
	}
	else if (is_big_endian != ByteOrder::is_host_big_endian()) {
		ByteOrder::swap_bytes(data, imagich.nx * imagich.ny * nz);
	}
#if 0
	int n_pad_imgs = 0;
	int old_num_imgs = imagich.count + 1;
	if (image_index > old_num_imgs) {
		n_pad_imgs = image_index - old_num_imgs;
	}
#endif

	// New way to write data which includes region writing.
	// If it is tested to be OK, remove the old code in the
	// #if 0  ... #endif block.
	EMUtil::process_region_io(data, img_file, WRITE_ONLY, 0,
							  sizeof(float), imagich.nx, imagich.ny,
							  nz, area, true);


#if 0
	size_t row_size = imagich.nx * sizeof(float);
	int nxy = imagich.nx * imagich.ny;

	for (int i = 0; i < nz; i++) {
		for (int j = imagich.ny - 1; j >= 0; j--) {
			fwrite(&data[i * nxy + j * imagich.nx], row_size, 1, img_file);
		}

		if (!is_new_img && (is_big_endian != ByteOrder::is_host_big_endian())) {
			ByteOrder::swap_bytes(data, imagich.nx * imagich.ny);
		}
	}
#endif
	EXITFUNC;
	return 0;
}

void ImagicIO::flush()
{
	fflush(img_file);
	fflush(hed_file);
}

Ctf * ImagicIO::read_ctf(const ImagicHeader& hed) const
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

void ImagicIO::write_ctf(const Ctf * const ctf, int image_index)
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

bool ImagicIO::is_complex_mode()
{
	init();
	if (datatype == IMAGIC_FLOAT_COMPLEX || datatype == IMAGIC_FFT_FLOAT_COMPLEX) {
		return true;
	}
	return false;
}

bool ImagicIO::is_image_big_endian()
{
	init();
	return is_big_endian;
}

int ImagicIO::get_nimg()
{
	init();
	return (imagich.count + 1);
}

ImagicIO::DataType ImagicIO::get_datatype_from_name(const char *name)
{
	DataType t = IMAGIC_UNKNOWN_TYPE;

	if (strncmp(name, "PACK",4) == 0) {
		t = IMAGIC_UCHAR;
	}
	else if (strncmp(name, "INTG",4) == 0) {
		t = IMAGIC_USHORT;
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

size_t ImagicIO::get_datatype_size(DataType t)
{
	size_t s = 0;
	switch (t) {
	case IMAGIC_UCHAR:
		s = sizeof(unsigned char);
		break;
	case IMAGIC_USHORT:
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

int ImagicIO::to_em_datatype(DataType t)
{
	switch (t) {
	case IMAGIC_UCHAR:
		return EMUtil::EM_UCHAR;
	case IMAGIC_USHORT:
		return EMUtil::EM_USHORT;
	case IMAGIC_FLOAT:
		return EMUtil::EM_FLOAT;
	case IMAGIC_FLOAT_COMPLEX:
		return EMUtil::EM_FLOAT_COMPLEX;
	default:
		break;
	}

	return EMUtil::EM_UNKNOWN;
}

void ImagicIO::make_header_host_endian(ImagicHeader & hed)
{
	if (is_big_endian != ByteOrder::is_host_big_endian()) {
		swap_header(hed);
	}
}


void ImagicIO::swap_header(ImagicHeader & hed)
{
	ByteOrder::swap_bytes((int *) &hed, NUM_4BYTES_PRE_IXOLD);
	ByteOrder::swap_bytes(&hed.ixold, NUM_4BYTES_AFTER_IXOLD);
	ByteOrder::swap_bytes((int *) &hed.space, NUM_4BYTES_AFTER_SPACE);
}
