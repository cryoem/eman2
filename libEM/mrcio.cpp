/**
 * $Id$
 */
#include "mrcio.h"
#include "byteorder.h"
#include "portable_fileio.h"
#include "geometry.h"
#include "log.h"
#include "util.h"
#include "emutil.h"
#include "ctf.h"

#include <math.h>
#include <string.h>


using namespace EMAN;

const char *MrcIO::CTF_MAGIC = "!-";
const char *MrcIO::SHORT_CTF_MAGIC = "!$";

MrcIO::MrcIO(string mrc_filename, IOMode rw)
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

	bool data_big_endian = ByteOrder::is_data_big_endian(&nz);

	if (data_big_endian != ByteOrder::is_host_big_endian()) {
		ByteOrder::swap_bytes(&nx);
		ByteOrder::swap_bytes(&ny);
		ByteOrder::swap_bytes(&nz);
		ByteOrder::swap_bytes(&mrcmode);
	}

	if (mrcmode == MRC_USHORT_COMPLEX || mrcmode == MRC_FLOAT_COMPLEX) {
		nx *= 2;
	}

	const int max_dim = 1 << 20;

	if ((mrcmode >= MRC_UCHAR && mrcmode < MRC_UNKNOWN) &&
		(nx > 1 && nx < max_dim) && (ny > 0 && ny < max_dim) && (nz > 0 && nz < max_dim)) {
		if (file_size > 0) {
			off_t file_size1 = nx * ny * nz * get_mode_size(mrcmode) + sizeof(MrcHeader);
			if (file_size == file_size1) {
				return true;
			}
		}
		else {
			return true;
		}
	}
	EXITFUNC;
	return false;
}

int MrcIO::read_header(Dict & dict, int image_index, const Region * area, bool )
{
	ENTERFUNC;

	check_read_access(image_index);
	check_region(area, IntSize(mrch.nx, mrch.ny, mrch.nz));

	dict["apix_x"] = mrch.xlen / (mrch.nx - 1);
	dict["apix_y"] = mrch.ylen / (mrch.ny - 1);

	int nz2 = mrch.nz - 1;
	if (mrch.nz == 1) {
		nz2 = 1;
	}

	dict["apix_z"] = mrch.zlen / nz2;

	dict["minimum"] = mrch.amin;
	dict["maximum"] = mrch.amax;
	dict["mean"] = mrch.amean;
	dict["datatype"] = to_em_datatype(mrch.mode);

	if (is_complex_mode()) {
		dict["is_complex"] = 1;
		dict["is_ri"] = 1;
	}

	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, mrch.nx, &xlen, mrch.ny, &ylen, mrch.nz, &zlen);

	dict["nx"] = xlen;
	dict["ny"] = ylen;
	dict["nz"] = zlen;

	if (area) {
		dict["origin_row"] = mrch.xorigin + mrch.xlen * area->origin[0];
		dict["origin_col"] = mrch.yorigin + mrch.xlen * area->origin[1];

		if (area->get_ndim() == 3 && mrch.nz > 1) {
			dict["origin_sec"] = mrch.zorigin + mrch.xlen * area->origin[2];
		}
		else {
			dict["origin_sec"] = mrch.zorigin;
		}
	}
	else {
		dict["origin_row"] = mrch.xorigin;
		dict["origin_col"] = mrch.yorigin;
		dict["origin_sec"] = mrch.zorigin;
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
	EXITFUNC;
	return 0;
}

int MrcIO::write_header(const Dict & dict, int image_index, const Region* , bool)
{
	ENTERFUNC;

	check_write_access(rw_mode, image_index, 1);

	int new_mode = to_mrcmode(dict["datatype"], (int) dict["is_complex"]);
	int nx = dict["nx"];
	int ny = dict["ny"];
	int nz = dict["nz"];
	is_ri =  dict["is_ri"];

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

	mrch.mx = nx - 1;
	mrch.my = ny - 1;
	mrch.mz = nz - 1;

	mrch.xlen = (nx - 1) * (float) dict["apix_x"];
	mrch.ylen = (ny - 1) * (float) dict["apix_y"];
	mrch.zlen = (nz - 1) * (float) dict["apix_z"];

	mrch.nxstart = -nx / 2;
	mrch.nystart = -ny / 2;
	mrch.nzstart = -nz / 2;

	mrch.xorigin = dict["origin_row"];
	mrch.yorigin = dict["origin_col"];

	if (is_new_file) {
		mrch.zorigin = dict["origin_sec"];
	}
	else {
		mrch.zorigin = (float) dict["origin_sec"] - (float) dict["apix_z"] * image_index;
	}

	sprintf(mrch.map, "MAP");
	mrch.machinestamp = generate_machine_stamp();

	MrcHeader mrch2 = mrch;
	
	if (opposite_endian) {
		swap_header(mrch2);
	}
	
	if (fwrite(&mrch2, sizeof(MrcHeader), 1, mrcfile) != 1) {
		throw ImageWriteException(filename, "MRC header");
	}

	mode_size = get_mode_size(mrch.mode);
	is_new_file = false;
	EXITFUNC;
	return 0;
}

int MrcIO::read_data(float *rdata, int image_index, const Region * area, bool )
{
	ENTERFUNC;

	check_read_access(image_index, rdata);

	if (area && is_complex_mode()) {
		LOGERR("Error: cannot read a region of a complex image.");
		return 1;
	}

	check_region(area, IntSize(mrch.nx, mrch.ny, mrch.nz));

	unsigned char *cdata = (unsigned char *) rdata;
	short *sdata = (short *) rdata;

	portable_fseek(mrcfile, sizeof(MrcHeader), SEEK_SET);

	EMUtil::process_region_io(cdata, mrcfile, READ_ONLY,
							  image_index, mode_size,
							  mrch.nx, mrch.ny, mrch.nz, area);

	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, mrch.nx, &xlen, mrch.ny, &ylen, mrch.nz, &zlen);

	int size = xlen * ylen * zlen;

	if (mrch.mode != MRC_UCHAR) {
		if (mode_size == sizeof(short)) {
			become_host_endian < short >(sdata, size);
		}
		else if (mode_size == sizeof(float)) {
			become_host_endian < float >(rdata, size);
		}
	}

	if (mrch.mode == MRC_UCHAR) {
		for (int i = size - 1; i >= 0; i--) {
			//rdata[i] = static_cast<float>(cdata[i]/100.0f - 1.28f);
			rdata[i] = static_cast < float >(cdata[i]);
		}
	}
	else if (mrch.mode == MRC_USHORT) {
		for (int i = size - 1; i >= 0; i--) {
			rdata[i] = static_cast < float >(sdata[i]);		
		}
	}

	if (is_complex_mode()) {
		Util::ap2ri(rdata, size);
		Util::flip_complex_phase(rdata, size);
	}
	EXITFUNC;
	return 0;
}

int MrcIO::write_data(float *data, int image_index, const Region* area, bool)
{
	ENTERFUNC;

	check_write_access(rw_mode, image_index, 1, data);

	int nx = mrch.nx;
	int ny = mrch.ny;
	int nz = mrch.nz;

	if (is_complex_mode()) {
		nx *= 2;
		size_t size = nx * ny * nz;
		if (!is_ri) {
			Util::ap2ri(data, size);
			is_ri = 1;
		}
		Util::flip_complex_phase(data, size);
	}

	portable_fseek(mrcfile, sizeof(MrcHeader), SEEK_SET);


	if (is_big_endian != ByteOrder::is_host_big_endian()) {
		if (mrch.mode != MRC_UCHAR) {
			size_t size = nz * nx * ny;
			if (mode_size == sizeof(short)) {
				ByteOrder::swap_bytes((short*) data, size);
			}
			else if (mode_size == sizeof(float)) {
				ByteOrder::swap_bytes((float*) data, size);
			}
		}
	}
	mode_size = get_mode_size(mrch.mode);
	
	unsigned char *cdata = (unsigned char *) data;
	short *sdata = (short *) data;
	
	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, nx, &xlen, mrch.ny, &ylen, mrch.nz, &zlen);
	int size = xlen * ylen * zlen;
	void * ptr_data = data;

	if (mrch.mode == MRC_UCHAR) {
		for (int i = 0; i < size; i++) {
			cdata[i] = static_cast < unsigned char >(data[i]);
		}
		ptr_data = cdata;
	}
	else if (mrch.mode == MRC_USHORT ||
			 mrch.mode == MRC_USHORT_COMPLEX) {
		for (int i = 0; i < size; i++) {
			sdata[i] = static_cast < unsigned short >(data[i]);
		}
		ptr_data = sdata;
	}
	
	EMUtil::process_region_io(ptr_data, mrcfile, WRITE_ONLY, image_index,
							  mode_size, nx, mrch.ny, mrch.nz, area);
	
	
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

			case MRC_USHORT:
			case MRC_USHORT_COMPLEX:
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

	delete[]cbuf;
	cbuf = 0;
#endif
	
	if (is_complex_mode()) {
		size_t size = nx * ny * nz;
		if (!is_ri) {
			Util::ap2ri(data, size);
			is_ri = 1;
		}
		Util::flip_complex_phase(data, size);
	}


	EXITFUNC;
	return 0;
}


bool MrcIO::is_complex_mode()
{
	init();
	if (mrch.mode == MRC_USHORT_COMPLEX || mrch.mode == MRC_FLOAT_COMPLEX) {
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

int MrcIO::write_ctf(const Ctf & ctf, int)
{
	ENTERFUNC;
	init();
	
	string ctf_str = ctf.to_string();
	sprintf(&mrch.labels[0][0], "%s%s", CTF_MAGIC, ctf_str.c_str());
	rewind(mrcfile);

	if (fwrite(&mrch, sizeof(MrcHeader), 1, mrcfile) != 1) {
		LOGERR("cannot write MRC header to file '%s'", filename.c_str());
		return 1;
	}
	EXITFUNC;
	return 0;
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
	case MRC_USHORT:
	case MRC_USHORT_COMPLEX:
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
	case MRC_USHORT:
		e = EMUtil::EM_USHORT;
		break;
	case MRC_USHORT_COMPLEX:
		e = EMUtil::EM_USHORT_COMPLEX;
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
			m = MRC_USHORT_COMPLEX;
		}
		else {
			m = MRC_USHORT;
		}
		break;
	case EMUtil::EM_SHORT_COMPLEX:
	case EMUtil::EM_USHORT_COMPLEX:
		m = MRC_USHORT_COMPLEX;
		break;
	case EMUtil::EM_CHAR:
	case EMUtil::EM_SHORT:
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
		p[0] = 0x44;
		p[1] = 0x44;
		p[2] = 0;
		p[3] = 0;
	}
	else {
		p[0] = 0x11;
		p[1] = 0x11;
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
