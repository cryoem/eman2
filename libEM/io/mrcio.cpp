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

#include <sys/stat.h>

#include "mrcio.h"
#include "portable_fileio.h"
#include "geometry.h"
#include "util.h"
#include "ctf.h"
#include "transform.h"

using namespace EMAN;

const char *MrcIO::CTF_MAGIC = "!-";
const char *MrcIO::SHORT_CTF_MAGIC = "!$";

MrcIO::MrcIO(const string & fname, IOMode rw)
:	ImageIO(fname, rw), mode_size(0),
		isFEI(false), is_ri(0), is_new_file(false),
		is_transpose(false), is_stack(false), stack_size(1),
		is_8_bit_packed(false), use_given_dimensions(true)
{
	memset(&mrch, 0, sizeof(MrcHeader));
	is_big_endian = ByteOrder::is_host_big_endian();
}

MrcIO::~MrcIO()
{
	if (file) {
		fclose(file);
		file = NULL;
	}
}

void MrcIO::init()
{
	ENTERFUNC;

	if (initialized) {
		return;
	}

	setbuf (stdout, NULL);

	IOMode rwmode = (rw_mode == WRITE_ONLY ? READ_WRITE : rw_mode);

	int error_type;
	struct stat status;

	error_type = stat(filename.c_str(), & status);
	is_new_file = (error_type != 0);

	initialized = true;

//	mrcfile = sfopen(filename, rwmode, &is_new_file);
	file = sfopen(filename, rwmode, NULL);

	string ext = Util::get_filename_ext(filename);

	if (ext != "") {
		isFEI    = (ext == "raw"  || ext == "RAW");
		is_stack = (ext == "mrcs" || ext == "MRCS");
	}

	if (! is_new_file) {
		if (fread(&mrch, sizeof(MrcHeader), 1, file) != 1) {
			throw ImageReadException(filename, "MRC header");
		}

		bool do_swap, have_err;

		check_swap((const int *) (& mrch), filename.c_str(), true,
					   do_swap, have_err);

		if (have_err  ||  ! is_valid(&mrch)) {
			throw ImageReadException(filename, "invalid MRC");
		}

		if (do_swap) {
			swap_header(mrch);

			is_big_endian = ! ByteOrder::is_host_big_endian();
		}
		else {
			is_big_endian =   ByteOrder::is_host_big_endian();
		}

		int max_labels = Util::get_min(mrch.nlabels, (int) MRC_NUM_LABELS);

		for (int ilabel = 0; ilabel < max_labels; ilabel++) {
			Util::replace_non_ascii(mrch.labels[ilabel], MRC_LABEL_SIZE);
		}

		// become_host_endian((int *) &mrch, NUM_4BYTES_PRE_MAP);
		// become_host_endian((int *) &mrch.machinestamp, NUM_4BYTES_AFTER_MAP);
		mode_size = get_mode_size(mrch.mode);

		if (is_complex_mode()) {
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

		if (mrch.nlabels > 0) {
			if (string(mrch.labels[0],3) == "Fei") {
				isFEI = true;
			}
		}

		is_transpose = (mrch.mapc == 2 && mrch.mapr == 1);

		if (is_stack) {
			stack_size = mrch.nz;
			mrch.nz = 1;
		}

		// Stuff added for 8 bit packed mode,
		// with 2 4-bit values packed into each 8-bit byte:
		float ny_to_nx_ratio;

		if (mrch.nx > 0) {
			ny_to_nx_ratio = (float) mrch.ny / (float) mrch.nx;
		}
		else {
			ny_to_nx_ratio = 1.0;
		}

		bool have_packed_label    = (mrch.nlabels > 0  &&
			  strstr(mrch.labels[0],   "4 bits packed") != NULL);
		bool have_packed_filename =
			 (strstr(filename.c_str(), "4_bits_packed") != NULL);

		bool ny_twice_nx = (fabs(ny_to_nx_ratio - 2.0f) <= 0.1f);
		bool double_nx   = (ny_twice_nx  &&  mrch.mode == MRC_UCHAR  &&
			 (have_packed_label  ||  have_packed_filename));

		use_given_dimensions = (! double_nx);

		is_8_bit_packed = (mrch.mode == MRC_UHEX || (mrch.imod_flags & 16)  ||  double_nx);

		if (getenv("DISALLOW_PACKED_FORMAT")) {
			use_given_dimensions = true;
			is_8_bit_packed = false;
		}

		if (is_8_bit_packed) {
			mrch.mode = MRC_UHEX;

			if (! use_given_dimensions) {
				mrch.nx = mrch.nx * 2;
			}
		}
	}

	EXITFUNC;
}

bool MrcIO::is_image_big_endian()
{
	init();

	return is_big_endian;
}

void MrcIO::check_swap(const int * data, const char * filnam, bool show_errors,
							  bool & do_swap, bool & have_err)
{
	int nx      = data[0];
	int ny      = data[1];
	int nz      = data[2];
	int mrcmode = data[3];
	int mapc    = data[16];	// Which axis corresponds to Columns  (X=1, Y=2, Z=3)
	int mapr    = data[17];	// Which axis corresponds to Rows     (X=1, Y=2, Z=3)
	int maps    = data[18];	// Which axis corresponds to Sections (X=1, Y=2, Z=3)
	int nsymbt  = data[23];	// the extra bytes for symmetry information
	int mach    = data[44]; // Machine endian stamp

	int nxw, nyw, nzw, modew, mapcw, maprw, mapsw, nsymw, machw;

	nxw   = nx;
	nyw   = ny;
	nzw   = nz;
	modew = mrcmode;
	mapcw = mapc;
	maprw = mapr;
	mapsw = maps;
	nsymw = nsymbt;
	machw = mach;

	ByteOrder::swap_bytes(&nxw);
	ByteOrder::swap_bytes(&nyw);
	ByteOrder::swap_bytes(&nzw);
	ByteOrder::swap_bytes(&modew);
	ByteOrder::swap_bytes(&mapcw);
	ByteOrder::swap_bytes(&maprw);
	ByteOrder::swap_bytes(&mapsw);
	ByteOrder::swap_bytes(&nsymw);
	ByteOrder::swap_bytes(&machw);

	const int max_dim = 1 << 20;

	int actual_stamp = MrcIO::generate_machine_stamp();

	bool debug = (getenv("DEBUG_MRC_SANITY"));

	string errmsg;

	have_err = false;

	if (mach == actual_stamp) {
		do_swap = false;
	}
	else if (machw == actual_stamp) {
		do_swap = true;
	}
	else {
		if (mrcmode == 0) {
			if (1 <= mapr  &&  mapr <= 3  &&
				 1 <= mapc  &&  mapc <= 3  &&
				 1 <= maps  &&  maps <= 3) {
				do_swap = false;
			}
			else if (1 <= maprw  &&  maprw <= 3  &&
				      1 <= mapcw  &&  mapcw <= 3  &&
				      1 <= mapsw  &&  mapsw <= 3) {
				do_swap = true;
			}
			else {
				double ave_xyz  = ((double) nx  + (double) ny  + (double) nz)  / 3.0;
				double ave_xyzw = ((double) nxw + (double) nyw + (double) nzw) / 3.0;

				if      (nx  > 0  &&  ny  > 0  &&  nz  > 0  &&  ave_xyz  <= max_dim) {
					do_swap = false;
				}
				else if (nxw > 0  &&  nyw > 0  &&  nzw > 0  &&  ave_xyzw <= max_dim) {
					do_swap = true;
				}
				else {
					have_err = true;
					do_swap  = false;
					errmsg = "MRC image dimensions nonpositive or too large.";
				}
			}
		}
		else {
			if (mrcmode > 0  &&  mrcmode < 128) {
				do_swap = false;
			}
			else if (modew > 0  &&  modew < 128) {
				do_swap = true;
			}
			else {
				have_err = true;
				do_swap  = false;
				errmsg = "MRC mode is not from 0 to 127.";
			}
		}
	}

	if (debug  ||  (have_err  &&  show_errors)) {
		if (filnam) {
			printf ("file name = '%s'.\n", filnam);
		}

		printf ("stamp: mach, file, swapd = %8.0x %8.0x %8.0x\n",
				  actual_stamp, mach, machw);
		printf ("stamp: mach, file, swapd = %d %d %d\n",
				  actual_stamp, mach, machw);
		printf ("nx, ny, nz, mode         = %d %d %d %d\n", nx,    ny,    nz,   mrcmode);
		printf ("nx, ny, nz, mode swapped = %d %d %d %d\n", nxw,   nyw,   nzw,  modew);
		printf ("mapc, mapr, maps         = %d %d %d\n",    mapc,  mapr,  maps);
		printf ("mapc, mapr, maps swapped = %d %d %d\n",    mapcw, maprw, mapsw);
		printf ("%s\n", errmsg.c_str());
	}
}

bool MrcIO::is_valid(const void * first_block, off_t file_size)
{
	ENTERFUNC;

	if (! first_block) { // Null pointer - no data
		return false;
	}

	const int * data = (const int *) first_block;

	int nx      = data[0];
	int ny      = data[1];
	int nz      = data[2];
	int mrcmode = data[3];
	int nsymbt  = data[23];	// the extra bytes for symmetry information

	const int max_dim = 1 << 20;

	bool do_swap, have_err;

	check_swap(data, NULL, false, do_swap, have_err);

	if (have_err) {
		return false;
	}

	if (do_swap) {
		ByteOrder::swap_bytes(&nx);
		ByteOrder::swap_bytes(&ny);
		ByteOrder::swap_bytes(&nz);
		ByteOrder::swap_bytes(&mrcmode);
		ByteOrder::swap_bytes(&nsymbt);
	}

	if (mrcmode == MRC_SHORT_COMPLEX || mrcmode == MRC_FLOAT_COMPLEX) {
		nx *= 2;
	}

	if ((mrcmode >= MRC_UCHAR &&
		(mrcmode < MRC_UNKNOWN || mrcmode == MRC_UHEX)) &&
		(nx > 1 && nx < max_dim) && (ny > 0 && ny < max_dim)
		 && (nz > 0 && nz < max_dim)) {

//#ifndef SPIDERMRC // Spider MRC files don't satisfy the following test

		if (file_size > 0) {
			off_t file_size1 = (off_t)nx * (off_t)ny * (off_t)nz *
				(off_t)get_mode_size(mrcmode) +
				(off_t)sizeof(MrcHeader) + nsymbt;

			if (file_size == file_size1) {
				return true;
			}

//			return false;

         // when size doesn't match, print error message instead of make it fail

			LOGWARN("image size check fails, still try to read it...");
		}
		else {
			return true;
		}

//#endif // SPIDERMRC

		return true;
	}

	EXITFUNC;

	return false;
}

int MrcIO::read_header(Dict & dict, int image_index, const Region * area, bool is_3d)
{
	init();

	if (isFEI) {
		return read_fei_header(dict, image_index, area, is_3d);
	}
	else {
		return read_mrc_header(dict, image_index, area, is_3d);
	}
}

int MrcIO::read_mrc_header(Dict & dict, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	if (image_index < 0) {
		image_index = 0;
	}

	if (image_index != 0  &&  ! is_stack) {
		throw ImageReadException(filename,
			"no stack allowed for MRC image. For take 2D slice out of 3D image, "
			"read the 3D image first, then use get_clip().");
	}

	check_region(area, FloatSize(mrch.nx, mrch.ny, mrch.nz), is_new_file, false);

	int xlen = 0, ylen = 0, zlen = 0;

	EMUtil::get_region_dims(area, mrch.nx, & xlen, mrch.ny, & ylen, mrch.nz, & zlen);

	dict["nx"] = xlen;
	dict["ny"] = ylen;
	dict["nz"] = zlen;
	dict["MRC.nx"] = mrch.nx;
	dict["MRC.ny"] = mrch.ny;
	dict["MRC.nz"] = mrch.nz;

	dict["datatype"] = to_em_datatype(mrch.mode);

	dict["MRC.nxstart"] = mrch.nxstart;
	dict["MRC.nystart"] = mrch.nystart;
	dict["MRC.nzstart"] = mrch.nzstart;

	dict["MRC.mx"] = mrch.mx;
	dict["MRC.my"] = mrch.my;
	dict["MRC.mz"] = mrch.mz;

	dict["MRC.mode"] = mrch.mode;

	dict["MRC.xlen"] = mrch.xlen;
	dict["MRC.ylen"] = mrch.ylen;
	dict["MRC.zlen"] = mrch.zlen;

	dict["MRC.alpha"] = mrch.alpha;
	dict["MRC.beta"] = mrch.beta;
	dict["MRC.gamma"] = mrch.gamma;

	dict["MRC.mapc"] = mrch.mapc;
	dict["MRC.mapr"] = mrch.mapr;
	dict["MRC.maps"] = mrch.maps;

	dict["MRC.minimum"] = mrch.amin;
	dict["MRC.maximum"] = mrch.amax;
	dict["MRC.mean"] = mrch.amean;
//	dict["minimum"] = mrch.amin;
//	dict["maximum"] = mrch.amax;
//	dict["mean"] = mrch.amean;

	dict["MRC.ispg"] = mrch.ispg;
	dict["MRC.nsymbt"] = mrch.nsymbt;

	dict["MRC.exttyp"] = string(mrch.exttyp, 4);
	dict["MRC.nversion"] = mrch.nversion;

	dict["MRC.imod_stamp"] = string(mrch.imod_stamp, 4);
	dict["MRC.imod_flags"] = mrch.imod_flags;
	dict["MRC.map"] = string(mrch.map, 4);		/* constant string "MAP "  */

	dict["IMODMRC.dvid"] = mrch.dvid;
	dict["IMODMRC.numintegers"] = mrch.numintegers;
	dict["IMODMRC.numfloats"] = mrch.numfloats;
	dict["IMODMRC.sub"] = mrch.sub;
	dict["IMODMRC.zfac"] = mrch.zfac;

	dict["IMODMRC.min2"] = mrch.min2;
	dict["IMODMRC.max2"] = mrch.max2;
	dict["IMODMRC.min3"] = mrch.min3;
	dict["IMODMRC.max3"] = mrch.max3;

	dict["IMODMRC.idtype"] = mrch.idtype;
	dict["IMODMRC.lens"] = mrch.lens;
	dict["IMODMRC.nd1"] = mrch.nd1;
	dict["IMODMRC.nd2"] = mrch.nd2;
	dict["IMODMRC.vd1"] = mrch.vd1;
	dict["IMODMRC.vd2"] = mrch.vd2;

	float apx = mrch.xlen / mrch.mx;
	float apy = mrch.ylen / mrch.my;
	float apz = mrch.zlen / mrch.mz;

	if (apx > 1000 || apx < 0.01) {
		dict["apix_x"] = 1.0f;
	}
	else {
		dict["apix_x"] = apx;
	}

	if (apy > 1000 || apy < 0.01) {
		dict["apix_y"] = 1.0f;
	}
	else {
		dict["apix_y"] = apy;
	}

	if (apz > 1000 || apz < 0.01) {
		dict["apix_z"] = 1.0f;
	}
	else {
		dict["apix_z"] = apz;
	}

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

	if (is_complex_mode()) {
		dict["is_complex"] = 1;
		dict["is_complex_ri"] = 1;
	}

	dict["MRC.machinestamp"] = mrch.machinestamp;

	dict["MRC.rms"] = mrch.rms;
	dict["sigma"] = mrch.rms;
	dict["MRC.nlabels"] = mrch.nlabels;

	char label[32];

	for (int i = 0; i < mrch.nlabels; i++) {
		sprintf(label, "MRC.label%d", i);
		dict[string(label)] = string(mrch.labels[i], MRC_LABEL_SIZE);
	}

	EMAN1Ctf ctf_;

	if (read_ctf(ctf_) == 0) {
		vector<float> vctf = ctf_.to_vector();
		dict["ctf"] = vctf;
	}

	if (is_transpose) {
		dict["nx"] = ylen;
		dict["ny"] = xlen;
		dict["MRC.nx"] = mrch.ny;
		dict["MRC.ny"] = mrch.nx;
		dict["MRC.mx"] = mrch.my;
		dict["MRC.my"] = mrch.mx;
		dict["apix_x"] = mrch.ylen / mrch.my;
		dict["apix_y"] = mrch.xlen / mrch.mx;
		dict["origin_x"] = mrch.yorigin;
		dict["origin_y"] = mrch.xorigin;
		dict["MRC.nxstart"] = mrch.nystart;
		dict["MRC.nystart"] = mrch.nxstart;
	}

	EXITFUNC;

	return 0;
}

int MrcIO::read_fei_header(Dict & dict, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	/* Read extended image header by specified image index */
	FeiMrcExtHeader feiexth;

	portable_fseek(file, sizeof(MrcHeader)+sizeof(FeiMrcExtHeader)*image_index, SEEK_SET);

	if (fread(&feiexth, sizeof(FeiMrcExtHeader), 1, file) != 1) {
		throw ImageReadException(filename, "FEI MRC extended header");
	}

	EXITFUNC;

	return 0;
}

int MrcIO::write_header(const Dict & dict, int image_index, const Region* area,
						EMUtil::EMDataType filestoragetype, bool use_host_endian)
{
	ENTERFUNC;

	string ext = Util::get_filename_ext(filename);
	is_stack = (ext == "mrcs"  ||  ext == "MRCS");

	bool append = (image_index == -1);

	if (image_index == -1) {
		image_index = 0;
	}

	if (image_index != 0  &&  ! is_stack) {
		throw ImageWriteException(filename, "MRC file does not support stack.");
	}

	int max_images = 0;

	if (! is_stack) {
		max_images = 1;
	}

	check_write_access(rw_mode, image_index, max_images);

	if (area) {
		check_region(area, FloatSize(mrch.nx, mrch.ny, mrch.nz), is_new_file);

		EXITFUNC;

		return 0;
	}

	int new_mode = to_mrcmode(filestoragetype, (int) dict["is_complex"]);
	int nx = dict["nx"];
	int ny = dict["ny"];
	int nz = dict["nz"];

	bool got_one_image = (nz > 1);

	is_ri =  dict["is_complex_ri"];

	bool opposite_endian = false;

	if (! is_new_file) {
		if (is_big_endian != ByteOrder::is_host_big_endian()) {
			opposite_endian = true;
		}

		portable_fseek(file, 0, SEEK_SET);
	}
	else {
		mrch.alpha = mrch.beta = mrch.gamma = 90.0f;
		mrch.mapc = 1;
		mrch.mapr = 2;
		mrch.maps = 3;
		mrch.nxstart = mrch.nystart = mrch.nzstart = 0;
	}

	if (nz <= 1 && dict.has_key("xform.projection") &&
					 ! dict.has_key("UCSF.chimera")) {
		Transform * t = dict["xform.projection"];
		Dict d = t->get_params("eman");
//		mrch.alpha   = d["alpha"];
//		mrch.beta    = d["beta"];
//		mrch.gamma   = d["gamma"];
		mrch.xorigin = d["tx"];
		mrch.yorigin = d["ty"];
		mrch.zorigin = d["tz"];

		if (t) {delete t; t = NULL;}
	}
	else if (nz > 1 && dict.has_key("xform.align3d") &&
						  ! dict.has_key("UCSF.chimera")) {
		Transform * t = dict["xform.align3d"];
		Dict d = t->get_params("eman");
//		mrch.alpha   = d["alpha"];
//		mrch.beta    = d["beta"];
//		mrch.gamma   = d["gamma"];
		mrch.xorigin = d["tx"];
		mrch.yorigin = d["ty"];
		mrch.zorigin = d["tz"];

		if (t) {delete t; t = NULL;}
	}

	if (dict.has_key("origin_x") && dict.has_key("origin_y") && dict.has_key("origin_z")){
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

	char label[32];

	for (int i = 0; i < MRC_NUM_LABELS; i++) {
		sprintf(label, "MRC.label%d", i);

		if (dict.has_key(label)) {
			strncpy(mrch.labels[i], (const char *) dict[label], MRC_LABEL_SIZE);
			mrch.labels[i][MRC_LABEL_SIZE-1] = 0;

			Util::replace_non_ascii(mrch.labels[i], MRC_LABEL_SIZE);

			mrch.nlabels = i + 1;
		}
	}

	if (mrch.nlabels < (MRC_NUM_LABELS - 1)) {
		strcpy(mrch.labels[mrch.nlabels], "EMAN ");
		strcat(mrch.labels[mrch.nlabels], Util::get_time_label().c_str());
		mrch.nlabels++;
	}

	mrch.mode = new_mode;

	if (is_complex_mode()) {
		mrch.nx = nx / 2;
	}
	else {
		mrch.nx = nx;
	}

	mrch.ny = ny;

	if (is_stack) {
		if (got_one_image) {
			stack_size = nz;
			image_index = 0;
		}
		else if (is_new_file) {
			stack_size = 1;
			image_index = stack_size - 1;
		}
		else if (append) {
			stack_size++;
			image_index = stack_size - 1;
		}
		else if (image_index >= stack_size) {
			stack_size = image_index + 1;
		}

		nz = stack_size;
		mrch.nz = nz;
	}
	else {
		if (is_new_file) {
			mrch.nz = nz;
		}
		else if (append) {
			nz++;
			mrch.nz = nz;
		}
		else if (image_index >= nz) {
			nz = image_index + 1;
			mrch.nz = nz;
		}
		else {
			mrch.nz = nz;
		}

		stack_size = 1;
	}

	mrch.ispg = dict.has_key("MRC.ispg") ? (int)dict["MRC.ispg"] : 0;
	mrch.nsymbt = 0;
	mrch.amin = dict["minimum"];
	mrch.amax = dict["maximum"];
	mrch.amean = dict["mean"];
	mrch.rms = dict["sigma"];

	mrch.mx = nx;
	mrch.my = ny;
	mrch.mz = nz;

	mrch.xlen = mrch.mx * (float) dict["apix_x"];
	mrch.ylen = mrch.my * (float) dict["apix_y"];
	mrch.zlen = mrch.mz * (float) dict["apix_z"];

	if (dict.has_key("MRC.nxstart")) {
		mrch.nxstart = dict["MRC.nxstart"];
	}
	else {
		mrch.nxstart = -nx / 2;
	}

	if (dict.has_key("MRC.nystart")) {
		mrch.nystart = dict["MRC.nystart"];
	}
	else {
		mrch.nystart = -ny / 2;
	}

	if (dict.has_key("MRC.nzstart")) {
		mrch.nzstart = dict["MRC.nzstart"];
	}
	else {
		mrch.nzstart = -nz / 2;
	}

	strncpy(mrch.map, "MAP ", 4);
	mrch.machinestamp = generate_machine_stamp();

	MrcHeader mrch2 = mrch;

	if (opposite_endian || !use_host_endian) {
		swap_header(mrch2);
	}

	if (fwrite(&mrch2, sizeof(MrcHeader), 1, file) != 1) {
		throw ImageWriteException(filename, "MRC header");
	}

	mode_size = get_mode_size(mrch.mode);
	is_new_file = false;

	if (is_stack  &&  ! got_one_image) {
		mrch.nz = 1;
	}

	EMUtil::getRenderLimits(dict, rendermin, rendermax, renderbits);

	// Do not write ctf to mrc header in EMAN2

//	if( dict.has_key("ctf") ) {
//		vector<float> vctf = dict["ctf"];
//		EMAN1Ctf ctf_;
//		ctf_.from_vector(vctf);
//		write_ctf(ctf_);
//	}

	EXITFUNC;

	return 0;
}

int MrcIO::read_data(float *rdata, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	if (! (isFEI || is_stack)) {
		// single image format, index can only be zero
		image_index = 0;
	}

	if (is_transpose && area != 0) {
		printf("Warning: This image dimension is in (y,x,z), "
				"region I/O not supported, return the whole image instead.");
	}

	check_read_access(image_index, rdata);

	if (area && is_complex_mode()) {
		LOGERR("Error: cannot read a region of a complex image.");

		return 1;
	}

	signed char *    scdata = (signed char *)    rdata;
	unsigned char *  cdata  = (unsigned char *)  rdata;
	short *          sdata  = (short *)          rdata;
	unsigned short * usdata = (unsigned short *) rdata;

	size_t size = 0;
	int xlen = 0, ylen = 0, zlen = 0;

	if (isFEI) {	// FEI extended MRC
		check_region(area, FloatSize(feimrch.nx, feimrch.ny, feimrch.nz), is_new_file, false);
		portable_fseek(file, sizeof(MrcHeader)+feimrch.next, SEEK_SET);

		EMUtil::process_region_io(cdata, file, READ_ONLY,
								  image_index, mode_size,
								  feimrch.nx, feimrch.ny, feimrch.nz, area);

		EMUtil::get_region_dims(area, feimrch.nx, &xlen, feimrch.ny, &ylen, feimrch.nz, &zlen);

		size = (size_t)xlen * ylen * zlen;
	}
	else {	// regular MRC
		check_region(area, FloatSize(mrch.nx, mrch.ny, mrch.nz), is_new_file, false);
		portable_fseek(file, sizeof(MrcHeader)+mrch.nsymbt, SEEK_SET);

		size_t modesize;

		if (mrch.mode == MRC_UHEX) {
			// Have MRC packed 8 bit format with 2 4-bit values in each 8-bit byte,
			// so the mode size is effectively half a byte, signalled by this value:
			modesize = 11111111;
		}
		else {
			modesize = mode_size;
		}

		EMUtil::process_region_io(cdata, file, READ_ONLY,
								  image_index, modesize,
								  mrch.nx, mrch.ny, mrch.nz, area);

		EMUtil::get_region_dims(area, mrch.nx, &xlen, mrch.ny, &ylen, mrch.nz, &zlen);

		size = (size_t)xlen * ylen * zlen;
	}

	if (mrch.mode != MRC_UCHAR  &&  mrch.mode != MRC_CHAR  &&
	    mrch.mode != MRC_UHEX) {

		if (mode_size == sizeof(short)) {
			become_host_endian < short >(sdata, size);
		}
		else if (mode_size == sizeof(float)) {
			become_host_endian < float >(rdata, size);
		}
	}

	if (mrch.mode == MRC_UHEX) {
		size_t num_pairs = size / 2;
		size_t num_pts   = num_pairs * 2;

		size_t ipt = num_pts;

		for (size_t ipair = num_pairs; ipair >= 1; ipair--) {
			unsigned int v = (unsigned int) cdata[ipair-1];
			ipt--;
			rdata[ipt] = (float)(v >> 4); // v / 16;
			ipt--;
			rdata[ipt] = (float)(v & 15); // v % 16;
		}
	}
	else if (mrch.mode == MRC_UCHAR) {
		for (size_t i = 0; i < size; ++i) {
			size_t j = size - 1 - i;
			// rdata[i] = static_cast<float>(cdata[i]/100.0f - 1.28f);
			rdata[j] = static_cast < float >(cdata[j]);
		}
	}
	else if (mrch.mode == MRC_CHAR) {
		for (size_t i = 0; i < size; ++i) {
			size_t j = size - 1 - i;
			// rdata[i] = static_cast<float>(cdata[i]/100.0f - 1.28f);
			rdata[j] = static_cast < float >(scdata[j]);
		}
	}
	else if (mrch.mode == MRC_SHORT) {
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

	if (is_transpose) {
		transpose(rdata, xlen, ylen, zlen);
	}

	if (is_complex_mode()) {
		if (! is_ri) {
			Util::ap2ri(rdata, size);
		}

		Util::flip_complex_phase(rdata, size);
		Util::rotate_phase_origin(rdata, xlen, ylen, zlen);
	}

	EXITFUNC;

	return 0;
}

template<class T>
void MrcIO::update_stats(const vector<T> &data)
{
	const auto [min, max] = std::minmax_element(std::begin(data), std::end(data));
	mrch.amin  = *min;
	mrch.amax  = *max;

	auto size = data.size();
	mrch.amean = (size > 0 ? std::accumulate(std::begin(data), std::end(data), 0.0, std::plus<T>()) / (double) size : 0.0);

	double square_sum = std::accumulate(std::begin(data), std::end(data), 0.0, [this](auto x1, auto x2){
		auto d1 = x1 - mrch.amean;
		auto d2 = x2 - mrch.amean;

		return d1 * d1 + d2 * d2;
	});

	mrch.rms = (float) (size > 1 ? std::sqrt(square_sum / (double) (size-1)) : 0.0);

	portable_fseek(file, 0, SEEK_SET);

	if (fwrite(& mrch, sizeof(MrcHeader), 1, file) != 1) {
		throw ImageWriteException(filename, "Error writing MRC header to update statistics.");
	}

	portable_fseek(file, sizeof(MrcHeader), SEEK_SET);
}

template<class T>
auto MrcIO::write_compressed(float *data, size_t size, int image_index, const Region* area) {
	void * ptr_data = data;

	if constexpr (!std::is_same<T, float>::value) {
		auto [rendered_data, count] = getRenderedDataAndRendertrunc<T>(data, size);

		update_stats<T>(rendered_data);
		ptr_data = rendered_data.data();
	}

	EMUtil::process_region_io(ptr_data, file, WRITE_ONLY, image_index,
							  mode_size, mrch.nx, mrch.ny, mrch.nz, area);
}

int MrcIO::write_data(float *data, int image_index, const Region* area,
					  EMUtil::EMDataType dt, bool use_host_endian)
{
	ENTERFUNC;

	if (is_stack  &&  (image_index == -1))
		image_index = stack_size - 1;

	int max_images = 0;

	if (! is_stack) {
		max_images  = 1;
		image_index = 0;
	}

	check_write_access(rw_mode, image_index, max_images, data);
	check_region(area, FloatSize(mrch.nx, mrch.ny, mrch.nz), is_new_file);

	int nx, ny, nz;

	if (! area) {
		nx = mrch.nx;
		ny = mrch.ny;
		nz = mrch.nz;
	}
	else {
		nx = (int)(area->get_width());
		ny = (int)(area->get_height());
		nz = (int)(area->get_depth());
	}

	if (is_stack  &&  ! (nz > 1)) {
		nz = 1;
		mrch.nz = 1;
	}

	size_t size = (size_t)nx * ny * nz;

	if (is_complex_mode()) {
		nx *= 2;

		if (! is_ri) {
			Util::ap2ri(data, size);
			is_ri = 1;
		}

		Util::flip_complex_phase(data, size);
		Util::rotate_phase_origin(data, nx, ny, nz);
	}

	portable_fseek(file, sizeof(MrcHeader), SEEK_SET);

	if(dt == EMUtil::EM_COMPRESSED) {
		if (renderbits <= 0)       mrch.mode = MRC_FLOAT;
		else if (renderbits <= 8)  mrch.mode = MRC_UCHAR;
		else if (renderbits <= 16) mrch.mode = MRC_USHORT;
	}

	if ((is_big_endian != ByteOrder::is_host_big_endian()) || ! use_host_endian) {
		if (mrch.mode != MRC_UCHAR  &&  mrch.mode != MRC_CHAR) {
			if (mode_size == sizeof(short))
				ByteOrder::swap_bytes((short*) data, size);
			else if (mode_size == sizeof(float))
				ByteOrder::swap_bytes((float*) data, size);
		}
	}

	mode_size = get_mode_size(mrch.mode);

	int truebits=0;
	switch(mrch.mode) {
		case MRC_UCHAR:	truebits=8; break;
		case MRC_SHORT: truebits=16; break;
		case MRC_FLOAT: truebits=32; break;
		case MRC_SHORT_COMPLEX: truebits=32; break;
		case MRC_FLOAT_COMPLEX: truebits=25; break;
		case MRC_USHORT: truebits=16; break;
		case MRC_UCHAR3: truebits=8; break;
		case MRC_CHAR: truebits=8; break;
		case MRC_UHEX: truebits=4; break;
		case MRC_UNKNOWN: truebits=0; break;
	}
	if (renderbits==0 || renderbits>truebits) renderbits=truebits;
	
	EMUtil::getRenderMinMax(data, nx, ny, rendermin, rendermax, renderbits, nz);

	if (mrch.mode == MRC_UCHAR)
		write_compressed<unsigned char>(data, size, image_index, area);
	else if (mrch.mode == MRC_CHAR)
		write_compressed<char>(data, size, image_index, area);
	else if (mrch.mode == MRC_SHORT || mrch.mode == MRC_SHORT_COMPLEX)
		write_compressed<short>(data, size, image_index, area);
	else if (mrch.mode == MRC_USHORT)
		write_compressed<unsigned short>(data, size, image_index, area);
	else
		write_compressed<float>(data, size, image_index, area);

	EXITFUNC;
	return 0;
}

bool MrcIO::is_complex_mode()
{
	init();

	return mrch.mode == MRC_SHORT_COMPLEX || mrch.mode == MRC_FLOAT_COMPLEX;
}

int MrcIO::read_ctf(Ctf & ctf, int)
{
	ENTERFUNC;

	init();

	size_t n = strlen(CTF_MAGIC);
	int err = 1;

	if (strncmp(mrch.labels[0], CTF_MAGIC, n) == 0) {
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

	strcpy (mrch.labels[0], CTF_MAGIC);
	strncat(mrch.labels[0], ctf_str.c_str(),
			  MRC_LABEL_SIZE - strlen(CTF_MAGIC) - 1);

	rewind(file);

	if (fwrite(&mrch, sizeof(MrcHeader), 1, file) != 1) {
		throw ImageWriteException(filename, "write CTF info to header failed");
	}

	EXITFUNC;
}

void MrcIO::flush()
{
	fflush(file);
}

int MrcIO::get_mode_size(int mm)
{
	MrcIO::MrcMode m = static_cast < MrcMode > (mm);

	int msize = 0;
	switch (m) {
	case MRC_CHAR:
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
	case MRC_CHAR:
		e = EMUtil::EM_CHAR;
		break;
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
	case EMUtil::EM_CHAR:
		m = MRC_CHAR;

		break;
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

int MrcIO::get_nimg()
{
	init();

	return stack_size;
}

int MrcIO::transpose(float *data, int xlen, int ylen, int zlen) const
{
	float * tmp = new float[xlen*ylen];

	for(size_t z=0; z<(size_t)zlen; ++z) {
		for(size_t y=0; y<(size_t)ylen; ++y) {
			for(size_t x=0; x<(size_t)xlen; ++x) {
				tmp[x*ylen+y] = data[z*xlen*ylen+y*xlen+x];
			}
		}

		std::copy(tmp, tmp+xlen*ylen, data+z*xlen*ylen);
	}

	delete [] tmp;

	return 0;
}
