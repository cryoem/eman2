/**
 * $Id$
 */
#include "imagicio.h"
#include "log.h"
#include "portable_fileio.h"
#include "util.h"
#include "emutil.h"
#include "geometry.h"
#include "ctf.h"
#include <assert.h>

using namespace EMAN;

const char *ImagicIO::HED_EXT = ".hed";
const char *ImagicIO::IMG_EXT = ".img";
const char *ImagicIO::REAL_TYPE_MAGIC = "REAL";
const char *ImagicIO::CTF_MAGIC = "!-";

ImagicIO::ImagicIO(string file, IOMode rw)
    :  filename(file), rw_mode(rw), hed_file(0), img_file(0), initialized(false)
{
    hed_filename = Util::get_filename_by_ext(filename, HED_EXT);
    img_filename = Util::get_filename_by_ext(filename, IMG_EXT);

    is_big_endian = ByteOrder::is_machine_big_endian();
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

int ImagicIO::init()
{
    static int err = 0;
    if (initialized) {
	return err;
    }
    Log::logger()->log("ImagicIO::init()");
    initialized = true;

    is_new_hed = false;
    is_new_img = false;

    hed_file = sfopen(hed_filename, rw_mode, &is_new_hed);
    img_file = sfopen(img_filename, rw_mode, &is_new_img);

    if (!hed_file) {
	err = 1;
	return err;
    }

    if (is_new_hed != is_new_img) {
	Log::logger()->warn("IMAGIC header file and data file should both exist or both not exist");
    }

    if (!is_new_hed) {
	if (fread(&imagich, sizeof(ImagicHeader), 1, hed_file) != 1) {
	    Log::logger()->error("cannot read IMAGIC header from file '%s'", hed_filename.c_str());
	    err = 1;
	    return err;
	}

	if (!is_valid(&imagich)) {
	    Log::logger()->error("file '%s' is not a valid IMAGIC file", hed_filename.c_str());
	    err = 1;
	    return err;
	}

	datatype = get_datatype_from_name(imagich.type);

	if (datatype != IMAGIC_USHORT && datatype != IMAGIC_FLOAT) {
	    Log::logger()->error("unsupported imagic data type: %s", imagich.type);
	    err = 1;
	    return err;
	}

	is_big_endian = ByteOrder::is_data_big_endian(&imagich.nx);
	make_header_right_endian(imagich);
	rewind(hed_file);
    }

    return 0;
}

bool ImagicIO::is_valid(const void *first_block)
{
    Log::logger()->log("ImagicIO::is_valid()");
    if (!first_block) {
	return false;
    }

    const int *data = static_cast<const int *>(first_block);
    int count = data[1];
    int headrec = data[3];
    int month = data[5];
    int hour = data[7];
    int nx = data[13];
    int ny = data[12];

    bool data_big_endian = ByteOrder::is_data_big_endian(&headrec);

    if (data_big_endian != ByteOrder::is_machine_big_endian()) {
	ByteOrder::swap_bytes(&count);
	ByteOrder::swap_bytes(&headrec);
	ByteOrder::swap_bytes(&month);
	ByteOrder::swap_bytes(&hour);
	ByteOrder::swap_bytes(&nx);
	ByteOrder::swap_bytes(&ny);
    }

    const int max_dim = 1 << 20;

    if (headrec == 1 &&
	count >= 0 && count < max_dim &&
	nx > 0 && nx < max_dim &&
	ny > 0 && ny < max_dim && month >= 0 && month <= 12 && hour >= 0 && hour <= 24) {
	return true;
    }
    return false;
}

int ImagicIO::read_header(Dict & dict, int image_index, const Region * area, bool is_3d)
{
    Log::logger()->log("ImagicIO::read_header() from file '%s'", hed_filename.c_str());

    if (check_read_access(image_index) != 0) {
	return 1;
    }

    int nimg = 1;

    if (is_3d) {
	nimg = imagich.count + 1;
	if (nimg <= 1) {
	    Log::logger()->warn("this is not a 3D IMAGIC. Read as a 2D");
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
	make_header_right_endian(hed);
    }

    if (check_region(area, Size(hed.nx, hed.ny, nimg)) != 0) {
	return 1;
    }

    int xlen = 0, ylen = 0, zlen = 0;
    EMUtil::get_region_dims(area, hed.nx, &xlen, hed.ny, &ylen, nimg, &zlen);

    dict["nx"] = EMObject(xlen);
    dict["ny"] = EMObject(ylen);
    dict["nz"] = EMObject(zlen);

    dict["datatype"] = EMObject(to_em_datatype(datatype));

    dict["minimum"] = EMObject(hed.min);
    dict["maximum"] = EMObject(hed.max);
    dict["mean"] = EMObject(hed.avdens);
    dict["sigma"] = EMObject(hed.sigma);

    dict["IMAGIC.imgnum"] = EMObject(hed.imgnum);
    dict["IMAGIC.count"] = EMObject(hed.count);
    dict["IMAGIC.error"] = EMObject(hed.error);

    dict["IMAGIC.headrec"] = EMObject(hed.headrec);
    dict["IMAGIC.mday"] = EMObject(hed.mday);
    dict["IMAGIC.month"] = EMObject(hed.month);

    dict["IMAGIC.year"] = EMObject(hed.year);
    dict["IMAGIC.hour"] = EMObject(hed.hour);
    dict["IMAGIC.minute"] = EMObject(hed.minute);

    dict["IMAGIC.sec"] = EMObject(hed.sec);
    dict["IMAGIC.reals"] = EMObject(hed.reals);
    dict["IMAGIC.pixels"] = EMObject(hed.pixels);

    dict["IMAGIC.type"] = EMObject(hed.type);
    dict["IMAGIC.ixold"] = EMObject(hed.ixold);
    dict["IMAGIC.iyold"] = EMObject(hed.iyold);

    dict["IMAGIC.oldav"] = EMObject(hed.oldav);
    dict["IMAGIC.label"] = EMObject(hed.label);
    dict["IMAGIC.mrc2"] = EMObject(hed.mrc2);

    return 0;
}

int ImagicIO::write_header(const Dict & dict, int image_index)
{
    Log::logger()->log("ImagicIO::write_header() to file '%s'", hed_filename.c_str());
    if (check_write_access(rw_mode, image_index) != 0) {
	return 1;
    }

    nz = dict["nz"].get_int();
    int n_new_img = nz;
    if (n_new_img > 1 && image_index != 0) {
	Log::logger()->error("to write 3D IMAGIC image, image index must be 0");
	return 1;
    }

    int nx = dict["nx"].get_int();
    int ny = dict["ny"].get_int();

    if (!is_new_hed) {
	make_header_right_endian(imagich);
	if (imagich.nx != nx || imagich.ny != ny) {
	    Log::logger()->error("new IMAGIC size %dx%d is not equal to existing size %dx%d",
				 nx, ny, imagich.nx, imagich.ny);
	    return 1;
	}
	rewind(hed_file);
    }

    ImagicHeader new_hed;
    memset(&new_hed, 0, sizeof(ImagicHeader));

    if (image_index == -1) {
	portable_fseek(hed_file, 0, SEEK_END);
	new_hed.imgnum = imagich.count + 2;
    }
    else {
	new_hed.imgnum = image_index + 1;
	portable_fseek(hed_file, sizeof(ImagicHeader) * image_index, SEEK_SET);
    }

    time_t cur_time = time(0);
    struct tm *tm = localtime(&cur_time);

    new_hed.count = 0;
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

    new_hed.min = dict["minimum"].get_float();
    new_hed.max = dict["maximum"].get_float();
    new_hed.avdens = dict["mean"].get_float();
    new_hed.sigma = dict["sigma"].get_float();
    /*
      new_hed.mrc1[0] = euler.phi();
      new_hed.mrc1[1] = euler.alt();
      new_hed.mrc1[2] = euler.az();
    */
    new_hed.mrc2 = n_new_img;

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

    if (!is_new_hed) {
	if (is_big_endian != ByteOrder::is_machine_big_endian()) {
	    swap_header(new_hed);
	}
    }

    strcpy(new_hed.type, REAL_TYPE_MAGIC);

    int n_pad_heds = 0;
    int old_num_imgs = imagich.count + 1;
    if (image_index > old_num_imgs) {
	n_pad_heds = image_index - old_num_imgs;
    }

    for (int i = 0; i < n_pad_heds + n_new_img; i++) {
	fwrite(&new_hed, sizeof(ImagicHeader), 1, hed_file);
    }

    if (is_new_hed && imagich.nx == 0) {
	imagich = new_hed;
	imagich.count = -1;
    }

    imagich.count += (n_pad_heds + n_new_img);

    if (is_big_endian != ByteOrder::is_machine_big_endian()) {
	ByteOrder::swap_bytes(&imagich.count);
    }
    portable_fseek(hed_file, sizeof(int), SEEK_SET);
    fwrite(&imagich.count, sizeof(int), 1, hed_file);

    if (!is_new_hed) {
	if (is_big_endian != ByteOrder::is_machine_big_endian()) {
	    swap_header(new_hed);
	}
    }

    return 0;
}

int ImagicIO::read_data(float *data, int image_index, const Region * area, bool is_3d)
{
    Log::logger()->log("ImagicIO::read_data() from file '%s'", img_filename.c_str());

    if (check_read_access(image_index, true, data) != 0) {
	return 1;
    }

    int nimg = 1;
    if (is_3d) {
	nimg = imagich.count + 1;
    }

    if (is_3d && imagich.count < 1) {
	Log::logger()->warn("this is not a 3D IMAGIC. Read as a 2D");
	is_3d = false;
    }
    if (check_region(area, Size(imagich.nx, imagich.ny, nimg)) != 0) {
	return 1;
    }

    rewind(img_file);

    unsigned short *sdata = (unsigned short *) data;
    unsigned char *cdata = (unsigned char *) data;
    size_t mode_size = get_datatype_size(datatype);

    int err = EMUtil::get_region_data(cdata, img_file, image_index, mode_size,
				      imagich.nx, imagich.ny, nimg, area, true);
    if (err) {
	return 1;
    }

#if 0
    int row_size = imagich.nx * mode_size;
    int sec_size = imagich.nx * imagich.ny * mode_size;

    for (int k = 0; k < nimg; k++) {
	for (int i = imagich.ny - 1; i >= 0; i--) {
	    if (fread(&cdata[k * sec_size + i * row_size], row_size, 1, img_file) != 1) {
		Log::logger()->error("incomplete data read: %d/%d blocks on file '%s'",
				     i, imagich.ny, filename.c_str());
		return 1;
	    }
	}
    }
#endif

    int img_size = imagich.nx * imagich.ny * nimg;

    if (datatype == IMAGIC_FLOAT) {
	become_platform_endian(data, img_size);
    }
    else if (datatype == IMAGIC_USHORT) {
	become_platform_endian((unsigned short *) cdata, img_size);

	for (int j = img_size - 1; j >= 0; j--) {
	    data[j] = static_cast<float>(sdata[j]);
	}
    }
    else {
	Log::logger()->error("unknown imagic data type");
	return 1;
    }

    return 0;
}

int ImagicIO::write_data(float *data, int image_index)
{
    Log::logger()->log("ImagicIO::write_data() to file '%s'", img_filename.c_str());
    if (check_write_access(rw_mode, image_index, true, data) != 0) {
	return 1;
    }

    if (nz == 1) {
	if (image_index == -1) {
	    portable_fseek(img_file, 0, SEEK_END);
	}
	else {
	    size_t sec_size = imagich.nx * imagich.ny * sizeof(float);
	    portable_fseek(img_file, ((off_t) sec_size) * image_index, SEEK_SET);
	}
    }

    if (!is_new_img && (is_big_endian != ByteOrder::is_machine_big_endian())) {
	ByteOrder::swap_bytes(data, imagich.nx * imagich.ny * nz);
    }

    int nimg = 0;
    if (nz == 1) {
	if (image_index > imagich.count) {
	    nimg = image_index - imagich.count;
	}
	nimg++;
    }
    else {
	nimg = nz;
    }

    int row_size = imagich.nx * sizeof(float);
    int sec_dim = imagich.nx * imagich.ny;

    for (int i = 0; i < nimg; i++) {
	for (int j = imagich.ny - 1; j >= 0; j--) {
	    fwrite(&data[i * sec_dim + j * imagich.nx], row_size, 1, img_file);
	}
    }

    if (!is_new_img && (is_big_endian != ByteOrder::is_machine_big_endian())) {
	ByteOrder::swap_bytes(data, imagich.nx * imagich.ny);
    }

    return 0;
}

int ImagicIO::read_ctf(Ctf & ctf, int )
{
    Log::logger()->log("ImagicIO::read_ctfit()");
    if (init() != 0) {
	return 1;
    }
    size_t n = strlen(CTF_MAGIC);
    int err = 1;
    if (strncmp(imagich.label, CTF_MAGIC, n) == 0) {
	err = ctf.from_string(string(&imagich.label[n]));
    }

    return err;
}

int ImagicIO::write_ctf(const Ctf & ctf, int )
{
    Log::logger()->log("ImagicIO::write_ctfit()");
    if (init() != 0) {
	return 1;
    }

    size_t n = strlen(CTF_MAGIC);
    strcpy(imagich.label, CTF_MAGIC);
    strncpy(&imagich.label[n], ctf.to_string().c_str(), sizeof(imagich.label) - n);

    rewind(hed_file);
    if (fwrite(&imagich, sizeof(ImagicHeader), 1, hed_file) != 1) {
	Log::logger()->error("cannot write Imagic header to file '%s'", hed_filename.c_str());
	return 1;
    }

    return 0;
}

bool ImagicIO::is_complex_mode()
{
    if (datatype == IMAGIC_FLOAT_COMPLEX || datatype == IMAGIC_FFT_FLOAT_COMPLEX) {
	return true;
    }
    return false;
}

bool ImagicIO::is_image_big_endian()
{
    return is_big_endian;
}

int ImagicIO::get_nimg()
{
    if (init() != 0) {
	return 0;
    }

    return (imagich.count + 1);
}

ImagicIO::DataType ImagicIO::get_datatype_from_name(const char *name)
{
    DataType t = IMAGIC_UNKNOWN_TYPE;

    if (strcmp(name, "PACK") == 0) {
	t = IMAGIC_UCHAR;
    }
    else if (strcmp(name, "INTG") == 0) {
	t = IMAGIC_USHORT;
    }
    else if (strcmp(name, REAL_TYPE_MAGIC) == 0) {
	t = IMAGIC_FLOAT;
    }
    else if (strcmp(name, "COMP") == 0) {
	t = IMAGIC_FLOAT_COMPLEX;
    }
    else if (strcmp(name, "RECO") == 0) {
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

void ImagicIO::make_header_right_endian(ImagicHeader & hed)
{
    become_platform_endian((int *) &hed, NUM_4BYTES_PRE_IXOLD);
    become_platform_endian(&hed.ixold, NUM_4BYTES_AFTER_IXOLD);
    become_platform_endian((int *) &hed.space, NUM_4BYTES_AFTER_SPACE);
}


void ImagicIO::swap_header(ImagicHeader & hed)
{
    ByteOrder::swap_bytes((int *) &hed, NUM_4BYTES_PRE_IXOLD);
    ByteOrder::swap_bytes(&hed.ixold, NUM_4BYTES_AFTER_IXOLD);
    ByteOrder::swap_bytes((int *) &hed.space, NUM_4BYTES_AFTER_SPACE);
}
