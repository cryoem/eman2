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
#include <assert.h>

using namespace EMAN;

const char *MrcIO::CTF_MAGIC = "!-";
const char *MrcIO::SHORT_CTF_MAGIC = "!$";

MrcIO::MrcIO(string mrc_filename, IOMode rw)
    :  filename(mrc_filename), rw_mode(rw), mrcfile(0), mode_size(0)
{
    memset(&mrch, 0, sizeof(MrcHeader));
    is_ri = false;
    is_big_endian = ByteOrder::is_machine_big_endian();
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

int MrcIO::init()
{
    static int err = 0;

    if (initialized) {
	return err;
    }
    initialized = true;

    mrcfile = sfopen(filename, rw_mode, &is_new_file);

    if (!mrcfile) {
	err = 1;
	return err;
    }

    if (!is_new_file) {
	if (fread(&mrch, sizeof(MrcHeader), 1, mrcfile) != 1) {
	    Log::logger()->error("read header failed on file: '%s'", filename.c_str());
	    err = 1;
	    return err;
	}

	if (!is_valid(&mrch)) {
	    Log::logger()->error("'%s' is not a valid MRC file", filename.c_str());
	    err = 1;
	    return err;
	}

	is_big_endian = ByteOrder::is_data_big_endian(&mrch.nz);
	become_platform_endian((int *) &mrch, NUM_4BYTES_PRE_MAP);
	become_platform_endian((int *) &mrch.machinestamp, NUM_4BYTES_AFTER_MAP);
	mode_size = get_mode_size(mrch.mode);

	if (mrch.nxstart != 0 || mrch.nystart != 0 || mrch.nzstart != 0) {
	    Log::logger()->warn("nx/ny/nz start not zero");
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

    return 0;
}


bool MrcIO::is_image_big_endian()
{
    init();
    return is_big_endian;
}

bool MrcIO::is_valid(const void *first_block, off_t file_size)
{
    Log::logger()->log("MrcIO::is_valid");

    if (!first_block) {
	return false;
    }

    const int *data = static_cast<const int *>(first_block);
    int nx = data[0];
    int ny = data[1];
    int nz = data[2];
    int mrcmode = data[3];

    bool data_big_endian = ByteOrder::is_data_big_endian(&nz);

    if (data_big_endian != ByteOrder::is_machine_big_endian()) {
	ByteOrder::swap_bytes(&nx);
	ByteOrder::swap_bytes(&ny);
	ByteOrder::swap_bytes(&nz);
	ByteOrder::swap_bytes(&mrcmode);
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

    return false;
}

int MrcIO::read_header(Dict & dict, int image_index, const Region * area, bool is_3d)
{
    Log::logger()->log("MrcIO::read_header() from file '%s'", filename.c_str());

    if (check_read_access(image_index) != 0) {
	return 1;
    }

    if (check_region(area, Size(mrch.nx, mrch.ny, mrch.nz)) != 0) {
	return 1;
    }
    if (area && area->get_ndim() > 2 && image_index > 0) {
	Log::logger()->warn("when reading 3D region in MRC, image index must be 0");
	image_index = 0;
    }

    assert(is_3d == false);

    dict["spacing_row"] = mrch.xlen / (mrch.nx - 1);
    dict["spacing_col"] = mrch.ylen / (mrch.ny - 1);

    int nz2 = mrch.nz - 1;
    if (mrch.nz == 1) {
	nz2 = 1;
    }
    
    dict["spacing_sec"] = mrch.zlen / nz2;

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
	dict["origin_row"] = mrch.xorigin + mrch.xlen * area->origin.x;
	dict["origin_col"] = mrch.yorigin + mrch.xlen * area->origin.y;

	if (area->get_ndim() == 3 && mrch.nz > 1) {
	    dict["origin_sec"] = mrch.zorigin + mrch.xlen * area->origin.z;
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
	dict[string(label)] = mrch.labels[0];
    }
    
    return 0;
}

int MrcIO::write_header(const Dict & dict, int image_index)
{
    Log::logger()->log("MrcIO::write_header() to file '%s'", filename.c_str());
    if (check_write_access(rw_mode, image_index) != 0) {
	return 1;
    }
    
    int new_mode = to_mrcmode(dict["datatype"].get_int(), dict["is_complex"].get_int());
    int nx = dict["nx"].get_int();
    int ny = dict["ny"].get_int();
    int nz = dict["nz"].get_int();
    is_ri = (bool) dict["is_ri"].get_int();
    
    if (!is_new_file) {
	if (is_big_endian != ByteOrder::is_machine_big_endian()) {
	    Log::logger()->error("cannot write to existing file '%s' which is in opposite byte order",
				 filename.c_str());
	    return 1;
	}
	
	if (nx != mrch.nx || ny != mrch.ny) {
	    Log::logger()->error("cannot write to different size file %s\n", filename.c_str());
	    return 1;
	}

	if (new_mode != mrch.mode) {
	    Log::logger()->error("cannot write to different mode file %sn", filename.c_str());
	    return 1;
	}

	if (mrch.nlabels < (MRC_NUM_LABELS - 1)) {
	    sprintf(&mrch.labels[mrch.nlabels][0], "EMAN %s", Util::get_time_label().c_str());
	    mrch.labels[mrch.nlabels][MRC_LABEL_SIZE - 1] = 0;
	    mrch.nlabels++;
	}

	portable_fseek(mrcfile, 0, SEEK_SET);
    }
    else {
	mrch.alpha = mrch.beta = mrch.gamma = 90.0;
	mrch.mapc = 1;
	mrch.mapr = 2;
	mrch.maps = 3;
	mrch.nxstart = mrch.nystart = mrch.nzstart = 0;
	mrch.nlabels = 1;
	sprintf(&mrch.labels[0][0], "EMAN %s", Util::get_time_label().c_str());
	mrch.labels[0][MRC_LABEL_SIZE - 1] = '\0';
	mrch.labels[1][0] = '\0';
    }

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
    mrch.amin = dict["minimum"].get_float();
    mrch.amax = dict["maximum"].get_float();
    mrch.amean = dict["mean"].get_float();

    mrch.mx = nx - 1;
    mrch.my = ny - 1;
    mrch.mz = nz - 1;

    mrch.xlen = (nx - 1) * dict["spacing_row"].get_float();
    mrch.ylen = (ny - 1) * dict["spacing_col"].get_float();
    mrch.zlen = (nz - 1) * dict["spacing_sec"].get_float();

    mrch.nxstart = -nx / 2;
    mrch.nystart = -ny / 2;
    mrch.nzstart = -nz / 2;
    
    mrch.xorigin = dict["origin_row"].get_float();
    mrch.yorigin = dict["origin_col"].get_float();

    if (is_new_file) {
	mrch.zorigin = dict["origin_sec"].get_float();
    }
    else {
	mrch.zorigin = dict["origin_sec"].get_float() -
	    dict["spacing_sec"].get_float() * image_index;
    }

    sprintf(mrch.map, "MAP");
    mrch.machinestamp = Util::generate_machine_stamp();
    /*
      if (mrcfile != stdout) {
      Util::file_lock_wait(mrcfile);
      }
    */

    if (fwrite(&mrch, sizeof(MrcHeader), 1, mrcfile) != 1) {
	Log::logger()->error("cannot write mrc header to file '%s'", filename.c_str());
	return 1;
    }

    mode_size = get_mode_size(mrch.mode);

    return 0;
}

int MrcIO::read_data(float *rdata, int image_index, const Region * area, bool is_3d)
{
    Log::logger()->log("MrcIO::read_data() from file '%s'", filename.c_str());

    if (check_read_access(image_index, true, rdata) != 0) {
	return 1;
    }

    assert(!is_3d);

    if (area && is_complex_mode()) {
	Log::logger()->error("Error: cannot read a region of a complex image.");
	return 1;
    }
    if (area && area->get_ndim() > 2 && image_index > 0) {
	Log::logger()->warn("when reading 3D region in MRC, image index must be 0");
	image_index = 0;
    }

    if (check_region(area, Size(mrch.nx, mrch.ny, mrch.nz)) != 0) {
	return 1;
    }

    unsigned char *cdata = (unsigned char *) rdata;
    short *sdata = (short *) rdata;

    portable_fseek(mrcfile, sizeof(MrcHeader) + mrch.nsymbt, SEEK_SET);

    int err = EMUtil::get_region_data(cdata, mrcfile, image_index, mode_size,
				      mrch.nx, mrch.ny, mrch.nz, area);
    if (err) {
	return 1;
    }

    int xlen = 0, ylen = 0, zlen = 0;
    EMUtil::get_region_dims(area, mrch.nx, &xlen, mrch.ny, &ylen, mrch.nz, &zlen);

    int size = xlen * ylen * zlen;

    if (mrch.mode != MRC_UCHAR) {
	if (mode_size == sizeof(short)) {
	    become_platform_endian<short>(sdata, size);
	}
	else if (mode_size == sizeof(float)) {
	    become_platform_endian<float>(rdata, size);
	}
    }

    for (int i = size - 1; i >= 0; i--) {
	switch (mrch.mode) {
	case MRC_UCHAR:
	    //rdata[i] = static_cast<float>(cdata[i]/100.0 - 1.28);
	    rdata[i] = static_cast<float>(cdata[i]);
	    break;
	case MRC_USHORT:
	    rdata[i] = static_cast<float>(sdata[i]);
	    break;
	}
    }

    if (is_complex_mode()) {	
	Util::ap2ri(rdata, size);
	Util::flip_complex_phase(rdata, size);
    }

    return 0;
}

int MrcIO::write_data(float *data, int image_index)
{
    Log::logger()->log("MrcIO::write_data() to file '%s'", filename.c_str());
    if (check_write_access(rw_mode, image_index, true, data) != 0) {
	return 1;
    }

    int nx = mrch.nx;
    int ny = mrch.ny;
    int nz = mrch.nz;

    Util::file_lock_wait(mrcfile);

    if (is_complex_mode()) {
	nx *= 2;
	int size = nx * ny * nz;
	if (!is_ri) {
	    Util::ap2ri(data, size);
	    is_ri = true;
	}
	Util::flip_complex_phase(data, size);
    }

    int nz1 = nz;
    if (image_index > 0) {
	portable_fseek(mrcfile, nx * ny * image_index * mode_size, SEEK_CUR);
	nz1 = 1;
    }

    int row_size = nx * get_mode_size(mrch.mode);
    int sec_size = nx * ny;

    unsigned char *cbuf = new unsigned char[row_size];
    unsigned short *sbuf = (unsigned short *) cbuf;
    int l = 0;

    for (int i = 0; i < nz1; i++) {
	for (int j = 0; j < ny; j++) {
	    int k = i * sec_size + j * nx;
	    switch (mrch.mode) {
	    case MRC_UCHAR:
		for (l = 0; l < nx; l++) {
		    cbuf[l] = static_cast<unsigned char>(data[k + l]);
		}
		fwrite(cbuf, row_size, 1, mrcfile);
		break;

	    case MRC_USHORT:
	    case MRC_USHORT_COMPLEX:
		for (l = 0; l < nx; l++) {
		    sbuf[l] = static_cast<unsigned short>(data[k + l]);
		}
		fwrite(sbuf, row_size, 1, mrcfile);
		break;

	    case MRC_FLOAT:
	    case MRC_FLOAT_COMPLEX:
		fwrite(&data[k], row_size, 1, mrcfile);
		break;
	    }
	}
    }
    
    delete [] cbuf;
    cbuf = 0;

    if (is_complex_mode()) {
	int size = nx * ny * nz1;
	if (!is_ri) {
	    Util::ap2ri(data, size);
	    is_ri = true;
	}
	Util::flip_complex_phase(data, size);
    }

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
    Log::logger()->log("MrcIO::read_ctfit()");
    if (init() != 0) {
	return 1;
    }

    size_t n = strlen(CTF_MAGIC);

    int err = 1;
    if (strncmp(&mrch.labels[0][0], CTF_MAGIC, n) == 0) {
	err = ctf.from_string(string(&mrch.labels[0][n]));
    }

    return err;
}

int MrcIO::write_ctf(const Ctf & ctf, int )
{
    Log::logger()->log("MrcIO::write_ctfit()");
    if (init() != 0) {
	return 1;
    }

    string ctf_str = ctf.to_string();

    mrch.nlabels = static_cast<int>((ctf_str.size() + sizeof(CTF_MAGIC)) / MRC_LABEL_SIZE + 1);
    sprintf(&mrch.labels[0][0], "%s%s", CTF_MAGIC, ctf_str.c_str());

    rewind(mrcfile);
    if (fwrite(&mrch, sizeof(MrcHeader), 1, mrcfile) != 1) {
	Log::logger()->error("cannot write MRC header to file '%s'", filename.c_str());
	return 1;
    }

    return 0;
}

int MrcIO::get_nimg()
{
    if (init() != 0) {
	return 0;
    }

    return 1;
}

int MrcIO::get_mode_size(int mm)
{
    MrcIO::MrcMode m = static_cast<MrcMode>(mm);

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


int MrcIO::to_mrcmode(int e, bool is_complex)
{
    Log::logger()->log("MrcIO::to_mrcmode");
    
    MrcMode m = MRC_UNKNOWN;
    EMUtil::EMDataType em_type = static_cast<EMUtil::EMDataType>(e);

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
	Log::logger()->error("unknown MRC mode: %s", EMUtil::get_datatype_string(em_type));
	m = MRC_UNKNOWN;
    }

    
    return m;
}
