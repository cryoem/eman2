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

#include <time.h>

using namespace EMAN;

const char *ImagicIO::HED_EXT = ".hed";
const char *ImagicIO::IMG_EXT = ".img";
const char *ImagicIO::REAL_TYPE_MAGIC = "REAL";
const char *ImagicIO::CTF_MAGIC = "!-";

ImagicIO::ImagicIO(string file, IOMode rw)
:	filename(file), rw_mode(rw), hed_file(0), img_file(0), initialized(false)
{
	hed_filename = Util::get_filename_by_ext(filename, HED_EXT);
	img_filename = Util::get_filename_by_ext(filename, IMG_EXT);

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

		is_big_endian = ByteOrder::is_data_big_endian(&imagich.nx);
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

	check_region(area, IntSize(hed.nx, hed.ny, nimg));

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

	dict["IMAGIC.type"] = hed.type;
	dict["IMAGIC.ixold"] = hed.ixold;
	dict["IMAGIC.iyold"] = hed.iyold;

	dict["IMAGIC.oldav"] = hed.oldav;
	dict["IMAGIC.label"] = hed.label;
	dict["IMAGIC.mrc2"] = hed.mrc2;
	EXITFUNC;
	return 0;
}

int ImagicIO::write_header(const Dict & dict, int image_index, const Region* , bool)
{
	ENTERFUNC;

	check_write_access(rw_mode, image_index);

	nz = dict["nz"];
	int n_new_img = nz;
	if (n_new_img > 1 && image_index != 0) {
		throw ImageWriteException(filename, "to write 3D IMAGIC image, image index must be 0");
	}

	int nx = dict["nx"];
	int ny = dict["ny"];

	if (!is_new_hed) {
		make_header_host_endian(imagich);
		if (imagich.nx != nx || imagich.ny != ny) {
			char desc[256];
			sprintf(desc, "new IMAGIC size %dx%d is not equal to existing size %dx%d",
					nx, ny, imagich.nx, imagich.ny);
			throw ImageWriteException(filename, desc);
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

		if (image_index > (imagich.count + 1)) {
			portable_fseek(hed_file, 0, SEEK_END);
		}
		else {
			portable_fseek(hed_file, sizeof(ImagicHeader) * image_index, SEEK_SET);
		}
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

	new_hed.min = dict["minimum"];
	new_hed.max = dict["maximum"];
	new_hed.avdens = dict["mean"];
	new_hed.sigma = dict["sigma"];
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
		if (is_big_endian != ByteOrder::is_host_big_endian()) {
			swap_header(new_hed);
		}
	}
	
	strcpy(new_hed.type, REAL_TYPE_MAGIC);

	int n_pad_heds = 0;
	int old_num_imgs = imagich.count + 1;
	if (image_index > old_num_imgs) {
		n_pad_heds = image_index - old_num_imgs;
	}

	if (image_index < old_num_imgs && n_new_img == 1) {
		n_new_img=0;
	}

	int imgnum_bak = new_hed.imgnum;

	for (int i = 0; i < n_pad_heds; i++) {
		new_hed.imgnum = i + 1 + old_num_imgs;
		fwrite(&new_hed, sizeof(ImagicHeader), 1, hed_file);
	}
	for (int i = 0; i < n_new_img; i++) {
		new_hed.imgnum = i + 1 + image_index;
		fwrite(&new_hed, sizeof(ImagicHeader), 1, hed_file);
	}
	
	new_hed.imgnum = imgnum_bak;

	if (is_new_hed && imagich.nx == 0) {
		imagich = new_hed;
		imagich.count = -1;
	}

	if (old_num_imgs < (image_index+1)) {
		imagich.count += (n_pad_heds + n_new_img);
		
		if (is_big_endian != ByteOrder::is_host_big_endian()) {
			ByteOrder::swap_bytes(&imagich.count);
		}
		
		portable_fseek(hed_file, sizeof(int), SEEK_SET);
		fwrite(&imagich.count, sizeof(int), 1, hed_file);
	}
	
	if (!is_new_hed) {
		if (is_big_endian != ByteOrder::is_host_big_endian()) {
			swap_header(new_hed);
		}
	}

	EXITFUNC;
	return 0;
}

int ImagicIO::read_data(float *data, int image_index, const Region * area, bool is_3d)
{
	ENTERFUNC;

	check_read_access(image_index, data);

	int nimg = 1;
	if (is_3d) {
		nimg = imagich.count + 1;
	}

	if (is_3d && imagich.count < 1) {
		LOGWARN("this is not a 3D IMAGIC. Read as a 2D");
		is_3d = false;
	}
	
	check_region(area, IntSize(imagich.nx, imagich.ny, nimg));

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

int ImagicIO::write_data(float *data, int image_index, const Region* , bool)
{
	ENTERFUNC;

	check_write_access(rw_mode, image_index, 0, data);
	
	if (nz == 1) {
		if (image_index == -1) {
			portable_fseek(img_file, 0, SEEK_END);
		}
		else {
			size_t sec_size = imagich.nx * imagich.ny * sizeof(float);
			portable_fseek(img_file, ((off_t) sec_size) * image_index, SEEK_SET);
		}
	}

	if (!is_new_img && (is_big_endian != ByteOrder::is_host_big_endian())) {
		ByteOrder::swap_bytes(data, imagich.nx * imagich.ny * nz);
	}
#if 0
	int n_pad_imgs = 0;
	int old_num_imgs = imagich.count + 1;	
	if (image_index > old_num_imgs) {
		n_pad_imgs = image_index - old_num_imgs;
	}
#endif
	
	size_t row_size = imagich.nx * sizeof(float);
	int sec_dim = imagich.nx * imagich.ny;

	for (int i = 0; i < nz; i++) {
		for (int j = imagich.ny - 1; j >= 0; j--) {
			fwrite(&data[i * sec_dim + j * imagich.nx], row_size, 1, img_file);
		}
	}

	if (!is_new_img && (is_big_endian != ByteOrder::is_host_big_endian())) {
		ByteOrder::swap_bytes(data, imagich.nx * imagich.ny);
	}
	EXITFUNC;
	return 0;
}

void ImagicIO::flush()
{	
	fflush(img_file);
	fflush(hed_file);
}

int ImagicIO::read_ctf(Ctf & ctf, int)
{
	ENTERFUNC;
	init();

	size_t n = strlen(CTF_MAGIC);
	int err = 1;
	if (strncmp(imagich.label, CTF_MAGIC, n) == 0) {
		err = ctf.from_string(string(&imagich.label[n]));
	}
	EXITFUNC;
	return err;
}

int ImagicIO::write_ctf(const Ctf & ctf, int)
{
	ENTERFUNC;
	init();

	size_t n = strlen(CTF_MAGIC);
	strcpy(imagich.label, CTF_MAGIC);
	strncpy(&imagich.label[n], ctf.to_string().c_str(), sizeof(imagich.label) - n);

	rewind(hed_file);
	if (fwrite(&imagich, sizeof(ImagicHeader), 1, hed_file) != 1) {
		throw ImageWriteException(hed_filename, "Imagic Header");
	}
	
	EXITFUNC;
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
	init();
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
