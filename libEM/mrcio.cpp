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

#include "mrcio.h"
#include "portable_fileio.h"
#include "geometry.h"
#include "util.h"
#include "ctf.h"
#include "transform.h"

using namespace EMAN;

const char *MrcIO::CTF_MAGIC = "!-";
const char *MrcIO::SHORT_CTF_MAGIC = "!$";

MrcIO::MrcIO(const string & mrc_filename, IOMode rw)
:	filename(mrc_filename), rw_mode(rw), mrcfile(0), mode_size(0)
{
	memset(&mrch, 0, sizeof(MrcHeader));
	is_ri = 0;
	is_big_endian = ByteOrder::is_host_big_endian();
	is_new_file = false;
	initialized = false;
}

MrcIO::~MrcIO()
{
	if (mrcfile) {
		fclose(mrcfile);
		mrcfile = 0;
	}
}

void MrcIO::init()
{
	ENTERFUNC;

	if (initialized) {
		return;
	}

	initialized = true;
	mrcfile = sfopen(filename, rw_mode, &is_new_file);

	if (!is_new_file) {
		if (fread(&mrch, sizeof(MrcHeader), 1, mrcfile) != 1) {
			throw ImageReadException(filename, "MRC header");
		}

		if (!is_valid(&mrch)) {
			throw ImageReadException(filename, "invalid MRC");
		}

		is_big_endian = ByteOrder::is_data_big_endian(&mrch.nz);
		if (is_big_endian != ByteOrder::is_host_big_endian()) {
			swap_header(mrch);
		}
		//become_host_endian((int *) &mrch, NUM_4BYTES_PRE_MAP);
		//become_host_endian((int *) &mrch.machinestamp, NUM_4BYTES_AFTER_MAP);
		mode_size = get_mode_size(mrch.mode);
		if(is_complex_mode()) {
			is_ri = 1;
		}

		if (mrch.nxstart != 0 || mrch.nystart != 0 || mrch.nzstart != 0) {
			LOGWARN("nx/ny/nz start not zero");
		}

		if (is_complex_mode()) {
			mrch.nx *= 2;
		}

		if (mrch.xlen == 0) {
			mrch.xlen = 1.0;
		}

		if (mrch.ylen == 0) {
			mrch.ylen = 1.0;
		}

		if (mrch.zlen == 0) {
			mrch.zlen = 1.0;
		}
	}
	EXITFUNC;
}


bool MrcIO::is_image_big_endian()
{
	init();
	return is_big_endian;
}

bool MrcIO::is_valid(const void *first_block, off_t file_size)
{
	ENTERFUNC;

	if (!first_block) {
		return false;
	}

	const int *data = static_cast < const int *>(first_block);
	int nx = data[0];
	int ny = data[1];
	int nz = data[2];
	int mrcmode = data[3];
	int nsymbt = data[23];	//this field specify the extra bytes for symmetry information

	bool data_big_endian = ByteOrder::is_data_big_endian(&nz);

	if (data_big_endian != ByteOrder::is_host_big_endian()) {
		ByteOrder::swap_bytes(&nx);
		ByteOrder::swap_bytes(&ny);
		ByteOrder::swap_bytes(&nz);
		ByteOrder::swap_bytes(&mrcmode);
	}

	if (mrcmode == MRC_SHORT_COMPLEX || mrcmode == MRC_FLOAT_COMPLEX) {
		nx *= 2;
	}

	const int max_dim = 1 << 20;

	if ((mrcmode >= MRC_UCHAR && mrcmode < MRC_UNKNOWN) &&
		(nx > 1 && nx < max_dim) && (ny > 0 && ny < max_dim) && (nz > 0 && nz < max_dim)) {
#ifndef SPIDERMRC // Spider MRC files don't satisfy the following test
		if (file_size > 0) {
			off_t file_size1 = (off_t)nx * (off_t)ny * (off_t)nz * (off_t)get_mode_size(mrcmode) + (off_t)sizeof(MrcHeader) + nsymbt;
			if (file_size == file_size1) {
				return true;
			}
			return false;
		}
		else {
			return true;
		}
#endif // SPIDERMRC
		return true;
	}
	EXITFUNC;
	return false;
}

int MrcIO::read_header(Dict & dict, int image_index, const Region * area, bool )
{
	ENTERFUNC;

	init();
	
	//single image format, index can only be zero
	if(image_index == -1) {
		image_index = 0;
	}

	if(image_index != 0) {
		throw ImageReadException(filename, "no stack allowed for MRC image. For take 2D slice out of 3D image, read the 3D image first, then use get_clip().");
	}
//	check_read_access(image_index);
	check_region(area, FloatSize(mrch.nx, mrch.ny, mrch.nz), is_new_file,false);

	dict["apix_x"] = mrch.xlen / mrch.mx;
	dict["apix_y"] = mrch.ylen / mrch.my;
	dict["apix_z"] = mrch.zlen / mrch.mz;

	dict["minimum"] = mrch.amin;
	dict["maximum"] = mrch.amax;
	dict["mean"] = mrch.amean;
	dict["datatype"] = to_em_datatype(mrch.mode);

	if (is_complex_mode()) {
		dict["is_complex"] = 1;
		dict["is_complex_ri"] = 1;
	}

	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, mrch.nx, &xlen, mrch.ny, &ylen, mrch.nz, &zlen);

	dict["nx"] = xlen;
	dict["ny"] = ylen;
	dict["nz"] = zlen;

	if (area) {
		dict["origin_x"] = mrch.xorigin + mrch.xlen * area->origin[0];
		dict["origin_y"] = mrch.yorigin + mrch.xlen * area->origin[1];

		if (area->get_ndim() == 3 && mrch.nz > 1) {
			dict["origin_z"] = mrch.zorigin + mrch.xlen * area->origin[2];
		}
		else {
			dict["origin_z"] = mrch.zorigin;
		}
	}
	else {
		dict["origin_x"] = mrch.xorigin;
		dict["origin_y"] = mrch.yorigin;
		dict["origin_z"] = mrch.zorigin;
	}

	dict["MRC.nxstart"] = mrch.nxstart;
	dict["MRC.nystart"] = mrch.nystart;
	dict["MRC.nzstart"] = mrch.nzstart;

	dict["MRC.mx"] = mrch.mx;
	dict["MRC.my"] = mrch.my;
	dict["MRC.mz"] = mrch.mz;

	dict["MRC.nx"] = mrch.nx;
	dict["MRC.ny"] = mrch.ny;
	dict["MRC.nz"] = mrch.nz;

	dict["MRC.xlen"] = mrch.xlen;
	dict["MRC.ylen"] = mrch.ylen;
	dict["MRC.zlen"] = mrch.zlen;

	dict["MRC.alpha"] = mrch.alpha;
	dict["MRC.beta"] = mrch.beta;
	dict["MRC.gamma"] = mrch.gamma;

	dict["MRC.mapc"] = mrch.mapc;
	dict["MRC.mapr"] = mrch.mapr;
	dict["MRC.maps"] = mrch.maps;

	dict["MRC.ispg"] = mrch.ispg;
	dict["MRC.nsymbt"] = mrch.nsymbt;
	dict["MRC.machinestamp"] = mrch.machinestamp;

	dict["MRC.rms"] = mrch.rms;
	dict["MRC.nlabels"] = mrch.nlabels;
	for (int i = 0; i < mrch.nlabels; i++) {
		char label[32];
		sprintf(label, "MRC.label%d", i);
		dict[string(label)] = mrch.labels[i];
	}

	EMAN1Ctf ctf_;
	if(read_ctf(ctf_) == 0) {
		vector<float> vctf = ctf_.to_vector();
		dict["ctf"] = vctf;
	}

	Dict dic;
	dic.put("type", "imagic");
	dic.put("alpha", mrch.alpha);
	dic.put("beta", mrch.beta);
	dic.put("gamma", mrch.gamma);
	dic.put("tx", mrch.xorigin);
	dic.put("ty", mrch.yorigin);
	dic.put("tz", mrch.zorigin);
	Transform * trans = new Transform(dic);
	if(zlen<=1) {
		dict["xform.projection"] = trans;
	}
	else {
		dict["xform.align3d"] = trans;
	}

	if(trans) {delete trans; trans=0;}
	EXITFUNC;
	return 0;
}

int MrcIO::write_header(const Dict & dict, int image_index, const Region* area,
						EMUtil::EMDataType filestoragetype, bool use_host_endian)
{
	ENTERFUNC;

	//single image format, index can only be zero
	image_index = 0;
	check_write_access(rw_mode, image_index, 1);
	if (area) {
		check_region(area, FloatSize(mrch.nx, mrch.ny, mrch.nz), is_new_file);
		EXITFUNC;
		return 0;
	}

	int new_mode = to_mrcmode(filestoragetype, (int) dict["is_complex"]);
	int nx = dict["nx"];
	int ny = dict["ny"];
	int nz = dict["nz"];
	is_ri =  dict["is_complex_ri"];

	bool opposite_endian = false;

	if (!is_new_file) {
		if (is_big_endian != ByteOrder::is_host_big_endian()) {
			opposite_endian = true;
		}
#if 0
		if (new_mode != mrch.mode) {
			LOGERR("cannot write to different mode file %s", filename.c_str());
			return 1;
		}
#endif
		portable_fseek(mrcfile, 0, SEEK_SET);
	}
	else {
		mrch.alpha = mrch.beta = mrch.gamma = 90.0f;
		mrch.mapc = 1;
		mrch.mapr = 2;
		mrch.maps = 3;
		mrch.nxstart = mrch.nystart = mrch.nzstart = 0;
	}

	if(nz<=1 && dict.has_key("xform.projection")) {
		Transform * t = dict["xform.projection"];
		Dict d = t->get_params("imagic");
		mrch.alpha = d["alpha"];
		mrch.beta = d["beta"];
		mrch.gamma = d["gamma"];
		mrch.xorigin = d["tx"];
		mrch.yorigin = d["ty"];
		mrch.zorigin = d["tz"];
		if(t) {delete t; t=0;}
	}
	else if(nz>1 && dict.has_key("xform.align3d")) {
		Transform * t = dict["xform.align3d"];
		Dict d = t->get_params("imagic");
		mrch.alpha = d["alpha"];
		mrch.beta = d["beta"];
		mrch.gamma = d["gamma"];
		mrch.xorigin = d["tx"];
		mrch.yorigin = d["ty"];
		mrch.zorigin = d["tz"];
		if(t) {delete t; t=0;}
	}

	if(dict.has_key("origin_x") && dict.has_key("origin_y") && dict.has_key("origin_z")){
		mrch.xorigin = (float)dict["origin_x"];
		mrch.yorigin = (float)dict["origin_y"];

		if (is_new_file) {
			mrch.zorigin = (float)dict["origin_z"];
		}
		else {
			mrch.zorigin = (float) dict["origin_z"] - (float) dict["apix_z"] * image_index;
		}
	}

	if (dict.has_key("MRC.nlabels")) {
		mrch.nlabels = dict["MRC.nlabels"];
	}

	for (int i = 0; i < MRC_NUM_LABELS; i++) {
		char label[32];
		sprintf(label, "MRC.label%d", i);
		if (dict.has_key(label)) {
			sprintf(&mrch.labels[i][0], "%s", (const char *) dict[label]);
			mrch.nlabels = i + 1;
		}
	}

	if (mrch.nlabels < (MRC_NUM_LABELS - 1)) {
		sprintf(&mrch.labels[mrch.nlabels][0], "EMAN %s", Util::get_time_label().c_str());
		mrch.nlabels++;
	}

	mrch.labels[mrch.nlabels][0] = '\0';
	mrch.mode = new_mode;

	if (is_complex_mode()) {
		mrch.nx = nx / 2;
	}
	else {
		mrch.nx = nx;
	}
	mrch.ny = ny;

	if (is_new_file) {
		mrch.nz = nz;
	}
	else if (image_index >= mrch.nz) {
		mrch.nz = image_index + 1;
	}

	mrch.ispg = 0;
	mrch.nsymbt = 0;
	mrch.amin = dict["minimum"];
	mrch.amax = dict["maximum"];
	mrch.amean = dict["mean"];

	/** the folowing lines are commented out.
	 * To make EMAN2 consistent with IMOD. Especially "header" command in IMOD. */
//	if(dict.has_key("MRC.mx")) {
//		mrch.mx = dict["MRC.mx"];
//	}
//	else {
		mrch.mx = nx;
//	}
//	if(dict.has_key("MRC.my")) {
//		mrch.my = dict["MRC.my"];
//	}
//	else {
		mrch.my = ny;
//	}
//	if(dict.has_key("MRC.mz")) {
//		mrch.mz = dict["MRC.mz"];
//	}
//	else {
		mrch.mz = nz;
//	}

	mrch.xlen = mrch.mx * (float) dict["apix_x"];
	mrch.ylen = mrch.my * (float) dict["apix_y"];
	mrch.zlen = mrch.mz * (float) dict["apix_z"];

	if(dict.has_key("MRC.nxstart")) {
		mrch.nxstart = dict["MRC.nxstart"];
	}
	else {
		mrch.nxstart = -nx / 2;
	}
	if(dict.has_key("MRC.nystart")) {
		mrch.nystart = dict["MRC.nystart"];
	}
	else {
		mrch.nystart = -ny / 2;
	}
	if(dict.has_key("MRC.nzstart")) {
		mrch.nzstart = dict["MRC.nzstart"];
	}
	else {
		mrch.nzstart = -nz / 2;
	}

	sprintf(mrch.map, "MAP ");
	mrch.machinestamp = generate_machine_stamp();
	if(dict.has_key("MRC.rms")) {
		mrch.rms = (float)dict["MRC.rms"];
	}

	MrcHeader mrch2 = mrch;

	if (opposite_endian || !use_host_endian) {
		swap_header(mrch2);
	}

	if (fwrite(&mrch2, sizeof(MrcHeader), 1, mrcfile) != 1) {
		throw ImageWriteException(filename, "MRC header");
	}

	mode_size = get_mode_size(mrch.mode);
	is_new_file = false;

	if( dict.has_key("ctf") ) {
		vector<float> vctf = dict["ctf"];
		EMAN1Ctf ctf_;
		ctf_.from_vector(vctf);
		write_ctf(ctf_);
	}

	EXITFUNC;
	return 0;
}

int MrcIO::read_data(float *rdata, int image_index, const Region * area, bool )
{
	ENTERFUNC;

	//single image format, index can only be zero
	image_index = 0;
	check_read_access(image_index, rdata);

	if (area && is_complex_mode()) {
		LOGERR("Error: cannot read a region of a complex image.");
		return 1;
	}

	check_region(area, FloatSize(mrch.nx, mrch.ny, mrch.nz), is_new_file, false);

	unsigned char *cdata = (unsigned char *) rdata;
	short *sdata = (short *) rdata;
	unsigned short *usdata = (unsigned short *) rdata;

	portable_fseek(mrcfile, sizeof(MrcHeader)+mrch.nsymbt, SEEK_SET);

	EMUtil::process_region_io(cdata, mrcfile, READ_ONLY,
							  image_index, mode_size,
							  mrch.nx, mrch.ny, mrch.nz, area);

	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, mrch.nx, &xlen, mrch.ny, &ylen, mrch.nz, &zlen);

	size_t size = xlen * ylen * zlen;

	if (mrch.mode != MRC_UCHAR) {
		if (mode_size == sizeof(short)) {
			become_host_endian < short >(sdata, size);
		}
		else if (mode_size == sizeof(float)) {
			become_host_endian < float >(rdata, size);
		}
	}

	if (mrch.mode == MRC_UCHAR) {
		for (size_t i = 0; i < size; ++i) {
			size_t j = size - 1 - i;
			//rdata[i] = static_cast<float>(cdata[i]/100.0f - 1.28f);
			rdata[j] = static_cast < float >(cdata[j]);
		}
	}
	else if (mrch.mode == MRC_SHORT ) {
		for (size_t i = 0; i < size; ++i) {
			size_t j = size - 1 - i;
			rdata[j] = static_cast < float >(sdata[j]);
		}
	}
	else if (mrch.mode == MRC_USHORT) {
		for (size_t i = 0; i < size; ++i) {
			size_t j = size - 1 - i;
			rdata[j] = static_cast < float >(usdata[j]);
		}
	}

	if (is_complex_mode()) {
		if(!is_ri) Util::ap2ri(rdata, size);
		Util::flip_complex_phase(rdata, size);
		Util::rotate_phase_origin(rdata, xlen, ylen, zlen);
	}
	EXITFUNC;
	return 0;
}

int MrcIO::write_data(float *data, int image_index, const Region* area,
					  EMUtil::EMDataType, bool use_host_endian)
{
	ENTERFUNC;
	//single image format, index can only be zero
	image_index = 0;
	check_write_access(rw_mode, image_index, 1, data);
	check_region(area, FloatSize(mrch.nx, mrch.ny, mrch.nz), is_new_file);

	int nx = mrch.nx;
	int ny = mrch.ny;
	int nz = mrch.nz;
	size_t size = nx * ny * nz;

	if (is_complex_mode()) {
		nx *= 2;
		if (!is_ri) {
			Util::ap2ri(data, size);
			is_ri = 1;
		}
		Util::flip_complex_phase(data, size);
		Util::rotate_phase_origin(data, nx, ny, nz);
	}

	portable_fseek(mrcfile, sizeof(MrcHeader), SEEK_SET);

	if ( (is_big_endian != ByteOrder::is_host_big_endian()) || !use_host_endian) {
		if (mrch.mode != MRC_UCHAR) {
			if (mode_size == sizeof(short)) {
				ByteOrder::swap_bytes((short*) data, size);
			}
			else if (mode_size == sizeof(float)) {
				ByteOrder::swap_bytes((float*) data, size);
			}
		}
	}
	mode_size = get_mode_size(mrch.mode);

//	int xlen = 0, ylen = 0, zlen = 0;
//	EMUtil::get_region_dims(area, nx, &xlen, mrch.ny, &ylen, mrch.nz, &zlen);
//	int size = xlen * ylen * zlen;
	void * ptr_data = data;

	float rendermin = 0.0f;
	float rendermax = 0.0f;
	getRenderMinMax(data, nx, ny, rendermin, rendermax);

	unsigned char *cdata = 0;
	short *sdata = 0;
	unsigned short *usdata = 0;
	if (mrch.mode == MRC_UCHAR) {
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
		ptr_data = cdata;
	}
	else if (mrch.mode == MRC_SHORT || mrch.mode == MRC_SHORT_COMPLEX) {
		sdata = new short[size];
		for (size_t i = 0; i < size; ++i) {
			if(data[i] <= rendermin) {
				sdata[i] = SHRT_MIN;
			}
			else if(data[i] >= rendermax) {
				sdata[i] = SHRT_MAX;
			}
			else {
				sdata[i]=(short)(((data[i]-rendermin)/(rendermax-rendermin))*(SHRT_MAX-SHRT_MIN) - SHRT_MAX);
			}
		}
		ptr_data = sdata;
	}
	else if (mrch.mode == MRC_USHORT) {
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
		ptr_data = usdata;
	}

	// New way to write data which includes region writing.
	// If it is tested to be OK, remove the old code in the
	// #if 0  ... #endif block.
	EMUtil::process_region_io(ptr_data, mrcfile, WRITE_ONLY, image_index,
							  mode_size, nx, mrch.ny, mrch.nz, area);

	if(cdata) {delete [] cdata; cdata=0;}
	if(sdata) {delete [] sdata; sdata=0;}
	if(usdata) {delete [] usdata; usdata=0;}

#if 0
	int row_size = nx * get_mode_size(mrch.mode);
	int sec_size = nx * ny;

	unsigned char *cbuf = new unsigned char[row_size];
	unsigned short *sbuf = (unsigned short *) cbuf;

	for (int i = 0; i < nz; i++) {
		int i2 = i * sec_size;
		for (int j = 0; j < ny; j++) {
			int k = i2 + j * nx;
			void *pbuf = 0;

			switch (mrch.mode) {
			case MRC_UCHAR:
				for (int l = 0; l < nx; l++) {
					cbuf[l] = static_cast < unsigned char >(data[k + l]);
				}
				pbuf = cbuf;
				fwrite(cbuf, row_size, 1, mrcfile);
				break;

			case MRC_SHORT:
			case MRC_SHORT_COMPLEX:
				for (int l = 0; l < nx; l++) {
					sbuf[l] = static_cast < short >(data[k + l]);
				}
				pbuf = sbuf;
				fwrite(sbuf, row_size, 1, mrcfile);
				break;

			case MRC_USHORT:
				for (int l = 0; l < nx; l++) {
					sbuf[l] = static_cast < unsigned short >(data[k + l]);
				}
				pbuf = sbuf;
				fwrite(sbuf, row_size, 1, mrcfile);
				break;

			case MRC_FLOAT:
			case MRC_FLOAT_COMPLEX:
				pbuf = &data[k];
				break;
			}
			if (pbuf) {
				fwrite(pbuf, row_size, 1, mrcfile);
			}
		}
	}

	if(cbuf)
	{
		delete[]cbuf;
		cbuf = 0;
	}
#endif

	EXITFUNC;
	return 0;
}


bool MrcIO::is_complex_mode()
{
	init();
	if (mrch.mode == MRC_SHORT_COMPLEX || mrch.mode == MRC_FLOAT_COMPLEX) {
		return true;
	}
	return false;
}


int MrcIO::read_ctf(Ctf & ctf, int)
{
	ENTERFUNC;
	init();
	size_t n = strlen(CTF_MAGIC);

	int err = 1;
	if (strncmp(&mrch.labels[0][0], CTF_MAGIC, n) == 0) {
		err = ctf.from_string(string(&mrch.labels[0][n]));
	}
	EXITFUNC;
	return err;
}

void MrcIO::write_ctf(const Ctf & ctf, int)
{
	ENTERFUNC;
	init();

	string ctf_str = ctf.to_string();
	sprintf(&mrch.labels[0][0], "%s%s", CTF_MAGIC, ctf_str.c_str());
	rewind(mrcfile);

	if (fwrite(&mrch, sizeof(MrcHeader), 1, mrcfile) != 1) {
		throw ImageWriteException(filename, "write CTF info to header failed");
	}
	EXITFUNC;
}

void MrcIO::flush()
{
	fflush(mrcfile);
}


int MrcIO::get_mode_size(int mm)
{
	MrcIO::MrcMode m = static_cast < MrcMode > (mm);

	int msize = 0;
	switch (m) {
	case MRC_UCHAR:
		msize = sizeof(char);
		break;
	case MRC_SHORT:
	case MRC_USHORT:
	case MRC_SHORT_COMPLEX:
		msize = sizeof(short);
		break;
	case MRC_FLOAT:
	case MRC_FLOAT_COMPLEX:
		msize = sizeof(float);
		break;
	default:
		msize = 0;
	}

	return msize;
}

int MrcIO::to_em_datatype(int m)
{
	EMUtil::EMDataType e = EMUtil::EM_UNKNOWN;

	switch (m) {
	case MRC_UCHAR:
		e = EMUtil::EM_UCHAR;
		break;
	case MRC_SHORT:
		e = EMUtil::EM_SHORT;
		break;
	case MRC_USHORT:
		e = EMUtil::EM_USHORT;
		break;
	case MRC_SHORT_COMPLEX:
		e = EMUtil::EM_SHORT_COMPLEX;
		break;
	case MRC_FLOAT:
		e = EMUtil::EM_FLOAT;
		break;
	case MRC_FLOAT_COMPLEX:
		e = EMUtil::EM_FLOAT_COMPLEX;
		break;
	default:
		e = EMUtil::EM_UNKNOWN;
	}
	return e;
}


int MrcIO::to_mrcmode(int e, int is_complex)
{
	MrcMode m = MRC_UNKNOWN;
	EMUtil::EMDataType em_type = static_cast < EMUtil::EMDataType > (e);

	switch (em_type) {
	case EMUtil::EM_UCHAR:
		m = MRC_UCHAR;
		break;
	case EMUtil::EM_USHORT:
		if (is_complex) {
			m = MRC_SHORT_COMPLEX;
		}
		else {
			m = MRC_USHORT;
		}
		break;
	case EMUtil::EM_SHORT:
		if (is_complex) {
			m = MRC_SHORT_COMPLEX;
		}
		else {
			m = MRC_SHORT;
		}
		break;
	case EMUtil::EM_SHORT_COMPLEX:
	case EMUtil::EM_USHORT_COMPLEX:
		m = MRC_SHORT_COMPLEX;
		break;
	case EMUtil::EM_CHAR:
	case EMUtil::EM_INT:
	case EMUtil::EM_UINT:
	case EMUtil::EM_FLOAT:
		if (is_complex) {
			m = MRC_FLOAT_COMPLEX;
		}
		else {
			m = MRC_FLOAT;
		}
		break;
	case EMUtil::EM_FLOAT_COMPLEX:
		m = MRC_FLOAT_COMPLEX;
		break;
	default:
		m = MRC_FLOAT;
	}

	return m;
}



int MrcIO::generate_machine_stamp()
{
	int stamp = 0;
	char *p = (char *) (&stamp);

	if (ByteOrder::is_host_big_endian()) {
		p[0] = 0x11;
		p[1] = 0x11;
		p[2] = 0;
		p[3] = 0;
	}
	else {
		p[0] = 0x44;
		p[1] = 0x41;
		p[2] = 0;
		p[3] = 0;
	}
	return stamp;
}

void MrcIO::swap_header(MrcHeader& mrch)
{
	ByteOrder::swap_bytes((int *) &mrch, NUM_4BYTES_PRE_MAP);
	ByteOrder::swap_bytes((int *) &mrch.machinestamp, NUM_4BYTES_AFTER_MAP);
}
