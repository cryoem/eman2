/**
 * $Id$
 */

#include "emdata.h"
#include "log.h"
#include "transform.h"
#include "io.h"
#include "ctf.h"
#include "filter.h"
#include "aligner.h"
#include "cmp.h"
#include "emfft.h"
#include "projector.h"

#include <float.h>
#include <math.h>
#include <assert.h>
#include <algorithm>

using namespace EMAN;

EMData::EMData()
{
	rdata = 0;
	supp = 0;
	ctf = 0;
	parent = 0;

	rfp = 0;
	flags = 0;
	// used to replace cube 'pixel'
	attr_dict["apix_x"] = 1.0f;
	attr_dict["apix_y"] = 1.0f;
	attr_dict["apix_z"] = 1.0f;

	attr_dict["is_complex"] = 0;
	attr_dict["is_ri"] = 0;

	nx = 0;
	ny = 0;
	nz = 0;
}

EMData::~EMData()
{
	if (rdata) {
		free(rdata);
		rdata = 0;
	}

	if (supp) {
		free(supp);
		supp = 0;
	}

	if (ctf) {
		delete ctf;
		ctf = 0;
	}

	if (rfp) {
		delete rfp;
		rfp = 0;
	}
}

void EMData::read_image(string filename, int img_index, bool nodata,
						const Region * r, bool is_3d)
{
	LOGDEBUG("\n\nEMData::read_image() on file '%s'", filename.c_str());

	ImageIO *imageio = EMUtil::get_imageio(filename, ImageIO::READ_ONLY);
	
	if (!imageio) {
		throw UnknownImageFormat(filename, __FILE__, __LINE__, "cannot create an image io");
	}
	else {
		int err = imageio->read_header(attr_dict, img_index, r, is_3d);
		if (err) {
			throw ImageHeaderReadError(filename, __FILE__, __LINE__,
									   "imageio read header failed");
		}
		else {
			if (imageio->is_complex_mode()) {
				set_complex(true);
			}
			if ((int) attr_dict["is_ri"] == 1) {
				set_ri(true);
			}

			nx = attr_dict["nx"];
			ny = attr_dict["ny"];
			nz = attr_dict["nz"];

			if (!ctf) {
				ctf = new SimpleCtf();
			}
			err = imageio->read_ctf(*ctf, img_index);
			if (err) {
				delete ctf;
				ctf = 0;
			}

			if (!nodata) {
				set_size(nx, ny, nz);
				int err = imageio->read_data(rdata, img_index, r, is_3d);
				if (err) {
					throw ImageDataReadError(filename, __FILE__, __LINE__,
											 "imageio read data failed");
				}
				else {
					flags |= EMDATA_NEEDUPD;
				}
			}
		}
	}

	if (nodata) {
		flags |= EMDATA_NODATA;
	}
}


void EMData::write_image(string filename, int img_index, EMUtil::ImageType imgtype,
						 bool header_only, bool use_host_endian) 
{
	LOGDEBUG("\n\nEMData::write_image() on file '%s'", filename.c_str());

	if (imgtype == EMUtil::IMAGE_UNKNOWN) {
		char *ext = strrchr(filename.c_str(), '.');
		if (ext) {
			ext++;
			imgtype = EMUtil::get_image_ext_type(ext);
		}
	}

	ImageIO *imageio = EMUtil::get_imageio(filename, ImageIO::READ_WRITE, imgtype);
	if (!imageio) {
		throw UnknownImageFormat(filename, __FILE__, __LINE__, "cannot create an image io");
	}
	else {
		update_stat();
		if (img_index < 0) {
			img_index = imageio->get_nimg();
		}
		int err = imageio->write_header(attr_dict, img_index, use_host_endian);
		if (err) {
			throw ImageHeaderWriteError(filename, __FILE__, __LINE__,
										"imageio write header failed");
		}
		else {
			if (ctf) {
				imageio->write_ctf(*ctf, img_index);
			}

			if (!header_only) {
				err = imageio->write_data(rdata, img_index, use_host_endian);
				if (err) {
					throw ImageDataWriteError(filename, __FILE__, __LINE__,
											  "imageio write data failed");
				}
			}
		}
	}
}

void EMData::append_image(string filename, EMUtil::ImageType imgtype, bool header_only)
{
	return write_image(filename, -1, imgtype, header_only);
}


void EMData::filter(string filtername, const Dict & params)
{
	Filter *f = Factory < Filter >::get(filtername, params);
	if (f) {
		f->process(this);
	}
}

float EMData::cmp(string cmpname, const Dict & params)
{
	float result = 0;

	Cmp *c = Factory < Cmp >::get(cmpname, params);
	if (c) {
		result = c->cmp(this);
	}
	return result;
}

EMData *EMData::align(string aligner_name, const Dict & params, string cmp_name)
{
	EMData *result = 0;

	Aligner *a = Factory < Aligner >::get(aligner_name, params);
	if (a) {
		if (cmp_name == "") {
			result = a->align(this);
		}
		else {
			result = a->align(this, cmp_name);
		}
	}

	return result;
}

EMData *EMData::project(string projector_name, const Dict & params)
{
	EMData *result = 0;
	Projector *p = Factory < Projector >::get(projector_name, params);
	if (p) {
		result = p->project3d(this);
	}

	return result;
}


EMData *EMData::copy(bool with_parent)
{
	EMData *ret = new EMData();

	ret->attr_dict = attr_dict;
	ret->set_size(nx, ny, nz);
	float *data = ret->get_data();
	memcpy(data, rdata, nx * ny * nz * sizeof(float));
	ret->done_data();
	ret->update();

	if (ctf) {
		ret->ctf = new SimpleCtf();
		ret->ctf->copy_from(ctf);
	}

	if (with_parent) {
		ret->parent = this;
	}
	else {
		ret->parent = 0;
	}

	ret->rfp = 0;

	ret->flags = flags & (EMDATA_COMPLEX | EMDATA_RI | EMDATA_HASCTF);

	ret->all_translation = all_translation;
	ret->all_rotation = all_rotation;

	ret->name = name;
	ret->path = path;
	ret->pathnum = pathnum;

	return ret;
}

EMData *EMData::copy_head()
{
	EMData *ret = new EMData();
	ret->attr_dict = attr_dict;
	ret->set_size(nx, ny, nz);

	if (ctf) {
		ret->ctf = new SimpleCtf();
		ret->ctf->copy_from(ctf);
	}

	ret->parent = this;
	ret->rfp = 0;

	ret->flags = flags & (EMDATA_COMPLEX | EMDATA_RI | EMDATA_HASCTF);

	ret->all_translation = all_translation;
	ret->all_rotation = all_rotation;

	ret->name = name;
	ret->path = path;
	ret->pathnum = pathnum;

	ret->update();

	return ret;
}


EMData *EMData::get_clip(const Region & area)
{
	if (get_ndim() != area.get_ndim()) {
		LOGERR("cannot get %dD clip out of %dD image", get_ndim(), area.get_ndim());
		return 0;
	}

	EMData *result = new EMData();
	int zsize = area.size.z;
	if (zsize == 0 || nz <= 1) {
		zsize = 1;
	}

	result->set_size(area.size.x, area.size.y, zsize);

	int x0 = (int) area.origin.x;
	x0 = x0 < 0 ? 0 : x0;

	int y0 = (int) area.origin.y;
	y0 = y0 < 0 ? 0 : y0;

	int z0 = (int) area.origin.z;
	z0 = z0 < 0 ? 0 : z0;
								
	int x1 = (int) (area.origin.x + area.size.x);
	x1 = x1 > nx ? nx : x1;
										
	int y1 = (int) (area.origin.y + area.size.y);
	y1 = y1 > ny ? ny : y1;

	int z1 = (int) (area.origin.z + area.size.z);
	z1 = z1 > nz ? nz : z1;
	if (z1 <= 0) {
		z1 = 1;
	}
	
	int xd0 = (int) (area.origin.x < 0 ? -area.origin.x : 0);
	int yd0 = (int) (area.origin.y < 0 ? -area.origin.y : 0);
	int zd0 = (int) (area.origin.z < 0 ? -area.origin.z : 0);

	int clipped_row_size = x1 * sizeof(float);
	int src_secsize = nx * ny;
	int dst_secsize = area.size.x * area.size.y;

	float *src = rdata + z0 * src_secsize + y0 * nx + x0;
	float *dst = result->get_data();
	dst += zd0 * dst_secsize + yd0 * area.size.x + xd0;

	int src_gap = src_secsize - (y1-y0) * nx;
	int dst_gap = dst_secsize - (y1-y0) * area.size.x;
	
	for (int i = z0; i < z1; i++) {
		for (int j = y0; j < y1; j++) {
			memcpy(dst, src, clipped_row_size);
			src += nx;
			dst += area.size.x;
		}
		src += src_gap;
		dst += dst_gap;
	}

	done_data();
	result->done_data();

	result->attr_dict["apix_x"] = attr_dict["apix_x"];
	result->attr_dict["apix_y"] = attr_dict["apix_y"];
	result->attr_dict["apix_z"] = attr_dict["apix_z"];

	float xorigin = attr_dict["origin_row"];
	float yorigin = attr_dict["origin_col"];
	float zorigin = attr_dict["origin_sec"];

	float apix_x = attr_dict["apix_x"];
	float apix_y = attr_dict["apix_y"];
	float apix_z = attr_dict["apix_z"];

	result->set_xyz_origin(xorigin + apix_x * area.origin.x,
						   yorigin + apix_y * area.origin.y,
						   zorigin + apix_z * area.origin.z);

	result->update();
	result->set_parent(0);
	result->set_name(name);
	result->set_path(path);
	result->set_pathnum(pathnum);

	return result;
}

void EMData::insert_clip(EMData * block, const IntPoint &origin)
{
	int nx1 = block->get_xsize();
	int ny1 = block->get_ysize();
	int nz1 = block->get_zsize();

	Region area(nx1, ny1, nz1, origin.x, origin.y, origin.z);
#if 0
	if (area.inside_region(Size(nx, ny, nz))) {
		LOGERR("outside of destination image not supported.");
		return;
	}
#endif
	int x0 = origin.x;
	int y0 = origin.y;
	int y1 = origin.y + ny1;
	int z0 = origin.z;
	int z1 = origin.z + nz1;

	size_t inserted_row_size = nx1 * sizeof(float);
	int nxy = nx * ny;
	int inserted_nxy = nx1 * ny1;

	float *inserted_data = block->get_data();

	float *src = inserted_data;
	float *dst = rdata + y0 * nx + x0;

	for (int i = z0; i < z1; i++) {
		src += inserted_nxy;
		dst += nxy;

		for (int j = y0; j < y1; j++) {
			memcpy(dst, src, inserted_row_size);
			src += nx1;
			dst += nx;;
		}
	}

	flags |= EMDATA_NEEDUPD;
}

EMData *EMData::get_rotated_clip(FloatPoint &center, Rotation &orient, Size &size, float scale)
{
	return NULL;
}

void EMData::insert_scaled_sum(EMData *block, const FloatPoint &center, float scale, float mult) 
{

if (get_ndim()==3) {
	// Start by determining the region to operate on
	int xs=(int)floor(block->get_xsize()*scale/2.0);
	int ys=(int)floor(block->get_ysize()*scale/2.0);
	int zs=(int)floor(block->get_zsize()*scale/2.0);
	int x0=(int)center.x-xs;
	int x1=(int)center.x+xs;
	int y0=(int)center.y-ys;
	int y1=(int)center.y+ys;
	int z0=(int)center.z-zs;
	int z1=(int)center.z+zs;
	
	if (x1<0||y1<0||z1<0||x0>get_xsize()||y0>get_ysize()||z0>get_zsize()) return;	// object is completely outside the target volume
	
	// make sure we stay inside the volume
	if (x0<0) x0=0;
	if (y0<0) y0=0;
	if (z0<0) z0=0;
	if (x1>get_xsize()) x1=get_xsize();
	if (y1>get_ysize()) y1=get_ysize();
	if (z1>get_zsize()) z1=get_zsize();
	
	float bx=block->get_xsize()/2.0;
	float by=block->get_ysize()/2.0;
	float bz=block->get_zsize()/2.0;
	
	for (int x=x0; x<x1; x++) {
		for (int y=y0; y<y1; y++) {
			for (int z=z0; z<z1; z++) {
				rdata[x + y * nx + z * nx * ny] += 
					block->get_value_at_interp((x-center.x)/scale+bx,(y-center.y)/scale+by,(z-center.z)/scale+bz);
			}
		}
	}
	update();
}
else if (get_ndim()==2) {
	// Start by determining the region to operate on
	int xs=(int)floor(block->get_xsize()*scale/2.0);
	int ys=(int)floor(block->get_ysize()*scale/2.0);
	int x0=(int)center.x-xs;
	int x1=(int)center.x+xs;
	int y0=(int)center.y-ys;
	int y1=(int)center.y+ys;
	
	if (x1<0||y1<0||x0>get_xsize()||y0>get_ysize()) return;	// object is completely outside the target volume
	
	// make sure we stay inside the volume
	if (x0<0) x0=0;
	if (y0<0) y0=0;
	if (x1>get_xsize()) x1=get_xsize();
	if (y1>get_ysize()) y1=get_ysize();
	
	float bx=block->get_xsize()/2.0;
	float by=block->get_ysize()/2.0;
	
	for (int x=x0; x<x1; x++) {
		for (int y=y0; y<y1; y++) {
			rdata[x + y * nx] += 
				block->get_value_at_interp((x-center.x)/scale+bx,(y-center.y)/scale+by);
		}
	}
	update();
}
else LOGERR("insert_scaled_sum supports only 2D and 3D data");

}

EMData *EMData::get_top_half() const
{
	EMData *half = new EMData();
	half->attr_dict = attr_dict;
	half->set_size(nx, ny, nz / 2);

	float *half_data = half->get_data();
	memcpy(half_data, &rdata[nz / 2 * nx * ny], sizeof(float) * nx * ny * nz / 2);
	half->done_data();

	float apix_z = attr_dict["apix_z"];
	float origin_sec = attr_dict["origin_sec"];
	origin_sec += apix_z * nz / 2;
	half->attr_dict["origin_sec"] = origin_sec;
	half->update();

	return half;
}



EMData *EMData::do_fft()
{
	if (flags & EMDATA_COMPLEX) {
		return this;
	}

	int nx2 = nx + 2;

	EMData *dat = copy_head();
	dat->set_size(nx2, ny, nz);


	float *d = dat->get_data();
	get_data();


	EMfft::real_to_complex_nd(rdata, d, nx, ny, nz);

	if (nz == 1) {
		int l = ny / 2 * nx2;

		for (int i = 0; i < ny / 2; i++) {
			int inx2 = i * nx2;
			for (int j = 0; j < nx2; j++) {
				int k = j + inx2;
				float f = d[k];
				d[k] = d[k + l];
				d[k + l] = f;
			}
		}
	}
	else if (ny != 1) {
		char *t = new char[nx2 * sizeof(float)];

		int k = nx2 * ny * (nz + 1) / 2;
		int l = nx2 * ny * (nz - 1) / 2;
		size_t jj = nx2 * sizeof(float);
		int ii = 0;

		for (int j = 0; j < nz / 2; j++) {
			for (int i = 0; i < ny; i++) {
				memcpy(t, d + ii, jj);

				if (i < ny / 2) {
					memcpy(d + ii, d + ii + k, jj);
					memcpy(d + ii + k, t, jj);
				}
				else {
					memcpy(d + ii, d + ii + l, jj);
					memcpy(d + ii + l, t, jj);
				}
				ii += nx2;
			}
		}
		delete[]t;
		t = 0;
	}


	float scale = 1.0f / (nx * ny * nz);
	dat->mult(scale);

	dat->done_data();
	dat->set_complex(true);
	dat->set_ri(true);

	int i = flags;
	done_data();
	flags = i & ~EMDATA_BUSY;

	return dat;
}

EMData *EMData::do_ift()
{
	if (!is_complex()) {
		LOGERR("cannot run IFT on real-space image");
		return 0;
	}

	if (!is_ri()) {
		LOGWARN("run IFT on AP data, only RI should be used. ");
	}

	EMData *dat = copy_head();
#if 0
	if (strlen(name) > 74)
		name[75] = 0;
#endif

	dat->set_size(nx, ny, nz);
	get_data();
	ap2ri();

	float *d = dat->get_data();
	int ndim = get_ndim();

	if (ndim >= 2) {
		memcpy((char *) d, (char *) rdata, nx * ny * nz * sizeof(float));
	}


	if (ndim == 1) {
		EMfft::complex_to_real_nd(rdata, d, nx - 2, ny, nz);
	}
	else if (ndim == 2) {
		int l = ny / 2 * nx;
		for (int i = 0; i < ny / 2; i++) {
			for (int j = 0; j < nx; j++) {
				int k = j + i * nx;
				float f = d[k];
				d[k] = d[k + l];
				d[k + l] = f;
			}
		}
	}
	else {
		char *t = new char[(nx + 2) * sizeof(float)];
		int k = nx * ny * (nz + 1) / 2;
		int l = nx * ny * (nz - 1) / 2;
		size_t jj = nx * sizeof(float);
		int ii = 0;

		for (int j = 0; j < nz / 2; j++) {
			for (int i = 0; i < ny; i++) {
				memcpy(t, d + ii, jj);

				if (i < ny / 2) {
					memcpy(d + ii, d + ii + k, jj);
					memcpy(d + ii + k, t, jj);
				}
				else {
					memcpy(d + ii, d + ii + l, jj);
					memcpy(d + ii + l, t, jj);
				}

				ii += nx;
			}
		}

		delete[]t;
		t = 0;
	}

	if (ndim >= 2) {
		EMfft::complex_to_real_nd(d, d, nx - 2, ny, nz);

		size_t row_size = (nx - 2) * sizeof(float);
		for (int i = 1; i < ny * nz; i++) {
			memcpy((char *) &d[i * (nx - 2)], (char *) &d[i * nx], row_size);
		}
	}

	dat->done_data();
#if 1
	dat->set_size(nx - 2, ny, nz);
#endif
	dat->update();
	dat->set_complex(false);
	dat->set_ri(false);

	int i = flags;
	done_data();
	flags = i & ~EMDATA_BUSY;

	return dat;
}


FloatPoint EMData::normalize_slice(EMData * slice, float alt, float az, float phi)
{
	if (!is_complex() || !slice->is_complex() || !parent) {
		LOGERR("normalize slice only works on complex images");
		return FloatPoint();
	}

	slice->ap2ri();
	get_data();

	float *norm = parent->get_data();
	float *dat = slice->get_data();
	float mx[6];
	mx[0] = (cos(phi) * cos(az) - cos(alt) * sin(az) * sin(phi));
	mx[1] = -(sin(phi) * cos(az) + cos(alt) * sin(az) * cos(phi));
	mx[2] = (cos(phi) * sin(az) + cos(alt) * cos(az) * sin(phi));
	mx[3] = (-sin(phi) * sin(az) + cos(alt) * cos(az) * cos(phi));
	mx[4] = sin(alt) * sin(phi);
	mx[5] = sin(alt) * cos(phi);

	float r = 0;
	float rn = 0;
	float pr = 0;
	float prn = 0;
	int nxy = nx * ny;

	for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nx / 2; x++) {
			float rad = hypot(x, y - ny / 2);

			if (rad < ny / 2 - 1) {
				float xx = x * mx[0] + (y - ny / 2) * mx[1];
				float yy = x * mx[2] + (y - ny / 2) * mx[3];
				float zz = x * mx[4] + (y - ny / 2) * mx[5];
				float cc = 1.0f;
				if (xx < 0) {
					xx = -xx;
					yy = -yy;
					zz = -zz;
					cc = -1.0f;
				}

				yy += ny / 2;
				zz += nz / 2;

				int x0 = 2 * (int) floor(xx + 0.5f);
				int y0 = (int) floor(yy + 0.5f);
				int z0 = (int) floor(zz + 0.5f);
				int i = x0 + y0 * nx + z0 * nxy;

				if (rdata[i] != 0 && rdata[i + 1] != 0) {
					float dt0 = hypot(rdata[i], rdata[i + 1]);
					float dt1 = hypot(dat[x * 2 + y * nx], dat[x * 2 + 1 + y * nx]);
					r += norm[i] * dt1;
					rn += dt0;

					float p1 = 0;
					float p2 = 0;

					if (rdata[i + 1] != 0 || rdata[i] != 0) {
						p1 = atan2(cc * rdata[i + 1], rdata[i]);
					}

					if (dat[x * 2 + 1 + y * nx] != 0 || dat[x * 2 + y * nx] != 0) {
						p2 = atan2(dat[x * 2 + 1 + y * nx], dat[x * 2 + y * nx]);
					}

					if (rad > 3.0f) {
						pr += Util::angle_sub_2pi(p1, p2) * dt0 * dt1 * norm[i];
						prn += dt0 * dt1 * norm[i];
					}
				}
			}
		}
	}

	float phaseres = (prn == 0) ? 0 : pr / prn;

	if (rn != 0) {
		r = r / rn;
	}
	else {
		r = 1.0f;
	}

	done_data();
	parent->done_data();
	slice->done_data();
	slice->update();

	return FloatPoint(r, phaseres);
}


void EMData::calc_hist(vector < float >&hist, float histmin, float histmax, bool add)
{
	static int prime[] = { 1, 3, 7, 11, 17, 23, 37, 59, 127, 253, 511 };

	if (histmin == histmax) {
		histmin = get_attr("minimum");
		histmax = get_attr("maximum");
	}

	if (!add) {
		for (size_t i = 0; i < hist.size(); i++) {
			hist[i] = 0;
		}
	}

	int p0 = 0;
	int p1 = 0;
	size_t size = nx * ny * nz;
	if (size < 300000) {
		p0 = 0;
		p1 = 0;
	}
	else if (size < 2000000) {
		p0 = 2;
		p1 = 3;
	}
	else if (size < 8000000) {
		p0 = 4;
		p1 = 6;
	}
	else {
		p0 = 7;
		p1 = 9;
	}

	if (!(flags & EMDATA_COMPLEX) && p0 > 0) {
		p0++;
		p1++;
	}

	int di = 0;
	float norm = 0;
	size_t n = hist.size();

	for (int k = p0; k <= p1; k++) {
		if (flags & EMDATA_COMPLEX) {
			di = prime[k] * 2;
		}
		else {
			di = prime[k];
		}

		norm += size / (float) di;
		float w = n / (histmax - histmin);

		for (int i = size - di; i >= 0; i -= di) {
			int j = Util::round((rdata[i] - histmin) * w);
			if (j >= 0 && j < (int) n) {
				hist[j] += 1;
			}
		}
	}

	for (size_t i = 0; i < hist.size(); i++) {
		hist[i] = hist[i] / norm;
	}
}



EMData *EMData::little_big_dot(EMData * with, bool do_sigma)
{
	EMData *ret = copy_head();
	ret->to_zero();

	int nx2 = with->get_xsize();
	int ny2 = with->get_ysize();
	float em = with->get_edge_mean();

	float *data = get_data();
	float *with_data = with->get_data();
	float *ret_data = ret->get_data();

	float sum2 = (Util::square((float)with->get_attr("sigma")) +
				  Util::square((float)with->get_attr("mean")));

	if (do_sigma) {
		for (int j = ny2 / 2; j < ny - ny2 / 2; j++) {
			for (int i = nx2 / 2; i < nx - nx2 / 2; i++) {
				float sum = 0;
				float sum1 = 0;
				float summ = 0;
				int k = 0;

				for (int jj = j - ny2 / 2; jj < j + ny2 / 2; jj++) {
					for (int ii = i - nx2 / 2; ii < i + nx2 / 2; ii++) {
						int l = ii + jj * nx;
						sum1 += Util::square(data[l]);
						summ += data[l];
						sum += data[l] * with_data[k];
						k++;
					}
				}
				float tmp_f1 = (sum1 / 2.0f - sum) / (nx2 * ny2);
				float tmp_f2 = Util::square((float)with->get_attr("mean") -
											summ / (nx2 * ny2));
				ret_data[i + j * nx] = sum2 + tmp_f1 - tmp_f2;
			}
		}
	}
	else {
		for (int j = ny2 / 2; j < ny - ny2 / 2; j++) {
			for (int i = nx2 / 2; i < nx - nx2 / 2; i++) {
				float eml = 0;
				float dot = 0;
				float dot2 = 0;

				for (int ii = i - nx2 / 2; ii < i + nx2 / 2; ii++) {
					eml += data[ii + (j - ny2 / 2) * nx] + data[ii + (j + ny2 / 2 - 1) * nx];
				}

				for (int jj = j - ny2 / 2; jj < j + ny2 / 2; jj++) {
					eml += data[i - nx2 / 2 + jj * nx] + data[i + nx2 / 2 - 1 + jj * nx];
				}

				eml /= (nx2 + ny2) * 2.0f;
				int k = 0;

				for (int jj = j - ny2 / 2; jj < j + ny2 / 2; jj++) {
					for (int ii = i - nx2 / 2; ii < i + nx2 / 2; ii++) {
						dot += (data[ii + jj * nx] - eml) * (with_data[k] - em);
						dot2 += Util::square(data[ii + jj * nx] - eml);
						k++;
					}
				}

				dot2 = sqrt(dot2);

				if (dot2 == 0) {
					ret_data[i + j * nx] = 0;
				}
				else {
					ret_data[i + j * nx] = dot / (nx2 * ny2 * dot2 * (float)with->get_attr("sigma"));
				}
			}
		}
	}

	done_data();
	with->done_data();
	ret->done_data();

	return ret;
}

void EMData::render_amp8(unsigned char *data, int x0, int y0, int ixsize, int iysize,
						 int bpl, float scale, int mingray, int maxgray,
						 float render_min, float render_max)
{
	if (is_complex()) {
		ri2ap();
	}

	if (render_max <= render_min) {
		render_max = render_min + 0.01f;
	}

	float rm = render_min;
	float inv_scale = 1.0f / scale;
	int ysize = iysize;
	int xsize = ixsize;

	int ymin = 0;
	if (iysize * inv_scale > ny) {
		ymin = (int) (iysize - ny / inv_scale);
	}

	float gs = (maxgray - mingray) / (render_max - render_min);
	if (render_max < render_min) {
		gs = 0;
		rm = FLT_MAX;
	}

	int dsx = -1;
	int dsy = 0;
	int remx = 0;
	int remy = 0;
	const int scale_n = 100000;

	int addi = 0;
	int addr = 0;
	if (inv_scale == floor(inv_scale)) {
		dsx = (int) inv_scale;
		dsy = (int) (inv_scale * nx);
	}
	else {
		addi = (int) floor(inv_scale);
		addr = (int) (scale_n * (inv_scale - floor(inv_scale)));
	}

	int xmin = 0;
	if (x0 < 0) {
		xmin = (int) (-x0 / inv_scale);
		xsize -= (int) floor(x0 / inv_scale);
		x0 = 0;
	}

	if ((xsize - xmin) * inv_scale > (nx - x0)) {
		xsize = (int) ((nx - x0) / inv_scale + xmin);
	}
	int ymax = ysize - 1;
	if (y0 < 0) {
		ymax = (int) (ysize + y0 / inv_scale - 1);
		ymin += (int) floor(y0 / inv_scale);
		y0 = 0;
	}

	if (xmin < 0) {
		xmin = 0;
	}

	if (ymin < 0) {
		ymin = 0;
	}

	if (xsize > ixsize) {
		xsize = ixsize;
	}

	if (ymax > iysize) {
		ymax = iysize;
	}

	int lmax = nx * ny - 1;

	if (is_complex()) {
		if (dsx != -1) {
			int l = y0 * nx;
			for (int j = ymax; j >= ymin; j--) {
				int ll = x0;
				for (int i = xmin; i < xsize; i++) {
					if (l + ll > lmax || ll >= nx - 2) {
						break;
					}
					int k = 0;
					if (ll >= nx / 2) {
						if (l >= (ny - inv_scale) * nx) {
							k = 2 * (ll - nx / 2) + 2;
						}
						else {
							k = 2 * (ll - nx / 2) + l + 2 + nx;
						}
					}
					else {
						k = nx * ny - (l + 2 * ll) - 2;
					}
					float t = rdata[k];
					if (t <= rm) {
						k = mingray;
					}
					else if (t >= render_max) {
						k = maxgray;
					}
					else {
						k = (int) (gs * (t - render_min));
						k += mingray;
					}
					data[i + j * bpl] = static_cast < unsigned char >(k);
					ll += dsx;
				}
				l += dsy;
			}
		}
		else {
			remy = 10;
			int l = y0 * nx;
			for (int j = ymax; j >= ymin; j--) {
				int br = l;
				remx = 10;
				int ll = x0;
				for (int i = xmin; i < xsize - 1; i++) {
					if (l + ll > lmax || ll >= nx - 2) {
						break;
					}
					int k = 0;
					if (ll >= nx / 2) {
						if (l >= (ny * nx - nx)) {
							k = 2 * (ll - nx / 2) + 2;
						}
						else {
							k = 2 * (ll - nx / 2) + l + 2 + nx;
						}
					}
					else {
						k = nx * ny - (l + 2 * ll) - 2;
					}

					float t = rdata[k];
					if (t <= rm)
						k = mingray;
					else if (t >= render_max) {
						k = maxgray;
					}
					else {
						k = (int) (gs * (t - render_min));
						k += mingray;
					}
					data[i + j * bpl] = static_cast < unsigned char >(k);
					ll += addi;
					remx += addr;
					if (remx > scale_n) {
						remx -= scale_n;
						ll++;
					}
				}
				l = br + addi * nx;
				remy += addr;
				if (remy > scale_n) {
					remy -= scale_n;
					l += nx;
				}
			}
		}
	}
	else {
		if (dsx != -1) {
			int l = x0 + y0 * nx;
			for (int j = ymax; j >= ymin; j--) {
				int br = l;
				for (int i = xmin; i < xsize; i++) {
					if (l > lmax) {
						break;
					}
					int k = 0;
					float t = rdata[l];
					if (t <= rm) {
						k = mingray;
					}
					else if (t >= render_max) {
						k = maxgray;
					}
					else {
						k = (int) (gs * (t - render_min));
						k += mingray;
					}
					data[i + j * bpl] = static_cast < unsigned char >(k);
					l += dsx;
				}
				l = br + dsy;
			}
		}
		else {
			remy = 10;
			int l = x0 + y0 * nx;
			for (int j = ymax; j >= ymin; j--) {
				int br = l;
				remx = 10;
				for (int i = xmin; i < xsize; i++) {
					if (l > lmax)
						break;
					int k = 0;
					float t = rdata[l];
					if (t <= rm) {
						k = mingray;
					}
					else if (t >= render_max) {
						k = maxgray;
					}
					else {
						k = (int) (gs * (t - render_min));
						k += mingray;
					}
					data[i + j * bpl] = static_cast < unsigned char >(k);
					l += addi;
					remx += addr;
					if (remx > scale_n) {
						remx -= scale_n;
						l++;
					}
				}
				l = br + addi * nx;
				remy += addr;
				if (remy > scale_n) {
					remy -= scale_n;
					l += nx;
				}
			}
		}
	}

	done_data();
}


void EMData::render_amp24(unsigned char *data, int x0, int y0, int ixsize, int iysize,
						  int bpl, float scale, int mingray, int maxgray,
						  float render_min, float render_max, void *ref,
						  void cmap(void *, int coord, unsigned char *tri))
{
	if (is_complex()) {
		ri2ap();
	}

	if (render_max <= render_min) {
		render_max = render_min + 0.01f;
	}

	float rm = render_min;
	float inv_scale = 1.0f / scale;
	int ysize = iysize;
	int xsize = ixsize;
	const int scale_n = 100000;

	int ymin = 0;
	if (iysize * inv_scale > ny) {
		ymin = (int) (iysize - ny / inv_scale);
	}
	float gs = (maxgray - mingray) / (render_max - render_min);
	if (render_max < render_min) {
		gs = 0;
		rm = FLT_MAX;
	}
	int dsx = -1;
	int dsy = 0;
	if (inv_scale == floor(inv_scale)) {
		dsx = (int) inv_scale;
		dsy = (int) (inv_scale * nx);
	}
	int addi = 0;
	int addr = 0;

	if (dsx == -1) {
		addi = (int) floor(inv_scale);
		addr = (int) (scale_n * (inv_scale - floor(inv_scale)));
	}

	int remx = 0;
	int remy = 0;
	int xmin = 0;
	if (x0 < 0) {
		xmin = (int) (-x0 / inv_scale);
		xsize -= (int) floor(x0 / inv_scale);
		x0 = 0;
	}

	if ((xsize - xmin) * inv_scale > (nx - x0)) {
		xsize = (int) ((nx - x0) / inv_scale + xmin);
	}
	int ymax = ysize - 1;
	if (y0 < 0) {
		ymax = (int) (ysize + y0 / inv_scale - 1);
		ymin += (int) floor(y0 / inv_scale);
		y0 = 0;
	}


	if (xmin < 0) {
		xmin = 0;
	}

	if (ymin < 0) {
		ymin = 0;
	}
	if (xsize > ixsize) {
		xsize = ixsize;
	}
	if (ymax > iysize) {
		ymax = iysize;
	}

	int lmax = nx * ny - 1;
	unsigned char tri[3];

	if (is_complex()) {
		if (dsx != -1) {
			int l = y0 * nx;
			for (int j = ymax; j >= ymin; j--) {
				int ll = x0;
				for (int i = xmin; i < xsize; i++, ll += dsx) {
					if (l + ll > lmax || ll >= nx - 2) {
						break;
					}
					int kk = 0;
					if (ll >= nx / 2) {
						if (l >= (ny - inv_scale) * nx) {
							kk = 2 * (ll - nx / 2) + 2;
						}
						else {
							kk = 2 * (ll - nx / 2) + l + 2 + nx;
						}
					}
					else {
						kk = nx * ny - (l + 2 * ll) - 2;
					}
					int k = 0;
					float t = rdata[kk];
					if (t <= rm) {
						k = mingray;
					}
					else if (t >= render_max) {
						k = maxgray;
					}
					else {
						k = (int) (gs * (t - render_min));
						k += mingray;
					}
					tri[0] = static_cast < unsigned char >(k);
					cmap(ref, kk, tri);
					data[i * 3 + j * bpl] = tri[0];
					data[i * 3 + 1 + j * bpl] = tri[1];
					data[i * 3 + 2 + j * bpl] = tri[2];
				}
				l += dsy;
			}
		}
		else {
			remy = 10;
			for (int j = ymax, l = y0 * nx; j >= ymin; j--) {
				int br = l;
				remx = 10;
				for (int i = xmin, ll = x0; i < xsize - 1; i++) {
					if (l + ll > lmax || ll >= nx - 2) {
						break;
					}
					int kk = 0;
					if (ll >= nx / 2) {
						if (l >= (ny * nx - nx)) {
							kk = 2 * (ll - nx / 2) + 2;
						}
						else {
							kk = 2 * (ll - nx / 2) + l + 2 + nx;
						}
					}
					else {
						kk = nx * ny - (l + 2 * ll) - 2;
					}
					int k = 0;
					float t = rdata[kk];
					if (t <= rm) {
						k = mingray;
					}
					else if (t >= render_max) {
						k = maxgray;
					}
					else {
						k = (int) (gs * (t - render_min));
						k += mingray;
					}
					tri[0] = static_cast < unsigned char >(k);
					cmap(ref, kk, tri);
					data[i * 3 + j * bpl] = tri[0];
					data[i * 3 + 1 + j * bpl] = tri[1];
					data[i * 3 + 2 + j * bpl] = tri[2];
					ll += addi;
					remx += addr;
					if (remx > scale_n) {
						remx -= scale_n;
						ll++;
					}
				}
				l = br + addi * nx;
				remy += addr;
				if (remy > scale_n) {
					remy -= scale_n;
					l += nx;
				}
			}
		}
	}
	else {
		if (dsx != -1) {
			for (int j = ymax, l = x0 + y0 * nx; j >= ymin; j--) {
				int br = l;
				for (int i = xmin; i < xsize; i++, l += dsx) {
					if (l > lmax) {
						break;
					}
					float t = rdata[l];
					int k = 0;
					if (t <= rm) {
						k = mingray;
					}
					else if (t >= render_max) {
						k = maxgray;
					}
					else {
						k = (int) (gs * (t - render_min));
						k += mingray;
					}
					tri[0] = static_cast < unsigned char >(k);
					cmap(ref, l, tri);
					data[i * 3 + j * bpl] = tri[0];
					data[i * 3 + 1 + j * bpl] = tri[1];
					data[i * 3 + 2 + j * bpl] = tri[2];
				}
				l = br + dsy;
			}
		}
		else {
			remy = 10;
			for (int j = ymax, l = x0 + y0 * nx; j >= ymin; j--) {
				int br = l;
				remx = 10;
				for (int i = xmin; i < xsize; i++) {
					if (l > lmax) {
						break;
					}
					float t = rdata[l];
					int k = 0;
					if (t <= rm) {
						k = mingray;
					}
					else if (t >= render_max) {
						k = maxgray;
					}
					else {
						k = (int) (gs * (t - render_min));
						k += mingray;
					}
					tri[0] = static_cast < unsigned char >(k);
					cmap(ref, l, tri);
					data[i * 3 + j * bpl] = tri[0];
					data[i * 3 + 1 + j * bpl] = tri[1];
					data[i * 3 + 2 + j * bpl] = tri[2];
					l += addi;
					remx += addr;
					if (remx > scale_n) {
						remx -= scale_n;
						l++;
					}
				}
				l = br + addi * nx;
				remy += addr;
				if (remy > scale_n) {
					remy -= scale_n;
					l += nx;
				}
			}
		}
	}

	done_data();
}



void EMData::calc_az_dist(int n, float a0, float da, float *d, float rmin, float rmax)
{
	if (nz > 1) {
		throw ImageDimError(__FILE__, __LINE__, 
							"Can only calculate AZ distribution on 2D images");
	}

	get_data();
	float *yc = new float[n];

	for (int i = 0; i < n; i++) {
		d[i] = 0;
		yc[i] = 0.00001f;
	}

	if (is_complex()) {
		int c = 0;
		for (int y = 0; y < ny; y++) {
			for (int x = 0; x < nx; x += 2, c += 2) {
				float x1 = x / 2.0f;
				float y1 = y - ny / 2.0f;
				float r = hypot(x1, y1);

				if (r >= rmin && r <= rmax) {
					float a = 0;

					if (y != ny / 2 || x != 0) {
						a = (atan2(y1, x1) - a0) / da;
					}

					int i = static_cast < int >(floor(a));
					a -= i;

					if (i == 0) {
						d[0] += rdata[c] * (1.0f - a);
						yc[0] += (1.0f - a);
					}
					else if (i == n - 1) {
						d[n - 1] += rdata[c] * a;
						yc[n - 1] += a;
					}
					else if (i > 0 && i < (n - 1)) {
						float h = 0;
						if (is_ri()) {
							h = hypot(rdata[c], rdata[c + 1]);
						}
						else {
							h = rdata[c];
						}

						d[i] += h * (1.0f - a);
						yc[i] += (1.0f - a);
						d[i + 1] += h * a;
						yc[i + 1] += a;
					}
				}
			}
		}
	}
	else {
		int c = 0;
		for (int y = 0; y < ny; y++) {
			for (int x = 0; x < nx; x++, c++) {
				float y1 = y - (ny - 1) / 2.0f;
				float x1 = x - (nx - 1) / 2.0f;
				float r = hypot(x1, y1);

				if (r >= rmin && r <= rmax) {
					float a = 0;
					if (x1 != 0 || y1 != 0) {
						a = atan2(y1, x1);
						if (a < 0) {
							a += M_PI * 2.0f;
						}
					}

					a = (a - a0) / da;
					int i = static_cast < int >(floor(a));
					a -= i;

					if (i == 0) {
						d[0] += rdata[c] * (1.0f - a);
						yc[0] += (1.0f - a);
					}
					else if (i == n - 1) {
						d[n - 1] += rdata[c] * a;
						yc[n - 1] += (a);
					}
					else if (i > 0 && i < (n - 1)) {
						d[i] += rdata[c] * (1.0f - a);
						yc[i] += (1.0f - a);
						d[i + 1] += rdata[c] * a;
						yc[i + 1] += a;
					}
				}
			}
		}
	}


	for (int i = 0; i < n; i++) {
		d[i] /= yc[i];
	}

	done_data();
	delete[]yc;
	yc = 0;
}


void EMData::ri2ap()
{
	if (!is_complex() || !is_ri()) {
		return;
	}

	int size = nx * ny * nz;
	for (int i = 0; i < size; i += 2) {
		float f = hypot(rdata[i], rdata[i + 1]);
		if (rdata[i] == 0 && rdata[i + 1] == 0) {
			rdata[i + 1] = 0;
		}
		else {
			rdata[i + 1] = atan2(rdata[i + 1], rdata[i]);
		}
		rdata[i] = f;
	}

	set_ri(false);
	flags |= EMDATA_NEEDUPD;
}

void EMData::ap2ri()
{
	if (!is_complex() || is_ri()) {
		return;
	}

	Util::ap2ri(rdata, nx * ny * nz);
	set_ri(true);
	flags |= EMDATA_NEEDUPD;
}


vector < float >EMData::calc_fourier_shell_correlation(EMData * with)
{
	if (!with) {
		return vector < float >();
	}

	if (!EMUtil::is_same_size(this, with)) {
		LOGERR("size mismatch in calc_fourier_shell_correlation");
		return vector < float >();
	}

	EMData *f1 = this;
	if (!is_complex()) {
		f1 = do_fft();
	}

	f1->ap2ri();

	EMData *f2 = with;
	if (!with->is_complex()) {
		f2 = with->do_fft();
	}
	f2->ap2ri();

	float *d1 = f1->get_data();
	float *d2 = f2->get_data();

	int f1_nx = f1->get_xsize();

	float *ret = new float[ny / 2];
	float *n1 = new float[ny / 2];
	float *n2 = new float[ny / 2];
	for (int i = 0; i < ny / 2; i++) {
		ret[i] = 0;
		n1[i] = 0;
		n2[i] = 0;
	}

	int nsec = nx * ny;

	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < f1_nx; i += 2) {
				int r = Util::round(Util::hypot3(i / 2, j - ny / 2, k - nz / 2));
				if (r >= 1 && r < ny / 2) {
					int ii = i + j * f1_nx + k * nsec;
					ret[r] += d1[ii] * d2[ii] + d1[ii + 1] * d2[ii + 1];
					n1[r] += d1[ii] * d1[ii] + d1[ii + 1] * d1[ii + 1];
					n2[r] += d2[ii] * d2[ii] + d2[ii + 1] * d2[ii + 1];
				}
			}
		}
	}

	vector < float >result(ny / 2);

	result[0] = 1;
	for (int i = 1; i < ny / 2; i++) {
		result[i] = ret[i] / (sqrt(n1[i] * n2[i]));
	}

	delete[]ret;
	ret = 0;

	delete[]n1;
	n1 = 0;
	delete[]n2;
	n2 = 0;

	return result;
}


void EMData::add(float f)
{
	if (f != 0) {
		flags |= EMDATA_NEEDUPD;
		int size = nx * ny * nz;
		for (int i = 0; i < size; i++) {
			rdata[i] += f;
		}
	}
}

void EMData::add(const EMData & em) 
{
	if (nx != em.get_xsize() || ny != em.get_ysize() || nz != em.get_zsize()) {
		throw NotSameImageSizeError(__FILE__, __LINE__,
									"cannot add 2 images with different sizes");
	}
	else {
		flags |= EMDATA_NEEDUPD;
		const float *src_data = em.get_data();
		int size = nx * ny * nz;

		for (int i = 0; i < size; i++) {
			rdata[i] += src_data[i];
		}
	}
}

void EMData::sub(float f)
{
	if (f != 0) {
		flags |= EMDATA_NEEDUPD;
		size_t size = nx * ny * nz;
		for (size_t i = 0; i < size; i++) {
			rdata[i] -= f;
		}
	}
}

void EMData::sub(const EMData & em) 
{
	if (nx != em.get_xsize() || ny != em.get_ysize() || nz != em.get_zsize()) {
		throw NotSameImageSizeError(__FILE__, __LINE__,
									"cannot substract 2 images with different sizes");
	}
	else {
		flags |= EMDATA_NEEDUPD;
		const float *src_data = em.get_data();
		size_t size = nx * ny * nz;
		
		for (size_t i = 0; i < size; i++) {
			rdata[i] -= src_data[i];
		}
	}
}



void EMData::mult(float f)
{
	if (is_complex()) {
		ap2ri();
	}
	if (f != 1) {
		flags |= EMDATA_NEEDUPD;
		size_t size = nx * ny * nz;
		for (size_t i = 0; i < size; i++) {
			rdata[i] *= f;
		}
	}
}

void EMData::mult(const EMData & em) 
{
	if (nx != em.get_xsize() || ny != em.get_ysize() || nz != em.get_zsize()) {
		throw NotSameImageSizeError(__FILE__, __LINE__,
									"cannot multiply 2 images with different sizes");
	}
	else {
		flags |= EMDATA_NEEDUPD;
		const float *src_data = em.get_data();
		size_t size = nx * ny * nz;

		for (size_t i = 0; i < size; i++) {
			rdata[i] *= src_data[i];
		}
	}
}

void EMData::div(float f)
{
	if (f != 0) {
		flags |= EMDATA_NEEDUPD;
		size_t size = nx * ny * nz;
		for (size_t i = 0; i < size; i++) {
			rdata[i] *= f;
		}
	}
}

void EMData::div(const EMData & em) 
{
	if (nx != em.get_xsize() || ny != em.get_ysize() || nz != em.get_zsize()) {
		throw NotSameImageSizeError(__FILE__, __LINE__,
									"cannot divide 2 images with different sizes");
	}
	else {
		flags |= EMDATA_NEEDUPD;
		const float *src_data = em.get_data();
		size_t size = nx * ny * nz;

		for (size_t i = 0; i < size; i++) {
			rdata[i] /= src_data[i];
		}
	}
}


void EMData::update_stat()
{
	if (!(flags & EMDATA_NEEDUPD)) {
		return;
	}

	float max = -FLT_MAX;
	float min = -max;

	double sum = 0;
	double square_sum = 0;

	int step = 1;
	if (is_complex() && !is_ri()) {
		step = 2;
	}

	int i = 0;
	double n_nonzero = 0;

	for (int j = 0; j < nz; j++) {
		for (int k = 0; k < ny; k++) {
			for (int l = 0; l < nx; l += step) {
				double v = rdata[i];
				if (v > max) {
					max = v;
				}

				if (v < min) {
					min = v;
				}
				if (v != 0) {
					n_nonzero++;
				}

				sum += v;
				square_sum += v * v;
				i += step;
			}
		}
	}

	size_t size = nx * ny * nz;
	float mean = sum * step / size;

	if (n_nonzero == 0) {
		n_nonzero = 1;
	}

	float mean_nonzero = sum * step / n_nonzero;
	float tmp1 = square_sum * step / size - mean * mean;

	if (tmp1 < 0) {
		tmp1 = 0;
	}

	float sigma = sqrt(tmp1);
	float sigma_nonzero = square_sum * step / n_nonzero - mean_nonzero * mean_nonzero;

	double kurtosis_sum = 0;
	double skewness_sum = 0;
	
	for (size_t k = 0; k < size; k++) {
		float t = (rdata[k] - mean) / sigma;
		kurtosis_sum += pow(t, 4.0f);
		skewness_sum +=  pow(t, 3.0f);
	}

	float kurtosis = (float)(kurtosis_sum / size - 3.0);
	float skewness = (float)(skewness_sum / size);
	
	attr_dict["minimum"] = min;
	attr_dict["maximum"] = max;
	attr_dict["mean"] = mean;
	attr_dict["sigma"] = sigma;
	attr_dict["kurtosis"] = kurtosis;
	attr_dict["skewness"] = skewness;
	
	attr_dict["square_sum"] = square_sum;
	attr_dict["mean_nonzero"] = mean_nonzero;
	attr_dict["sigma_nonzero"] = sigma_nonzero;
	attr_dict["is_complex"] = (int) is_complex();
	attr_dict["is_ri"] = (int) is_ri();

	flags &= ~EMDATA_NEEDUPD;
}


void EMData::set_size(int x, int y, int z)
{
	assert(x > 0 && y > 0 && z > 0);
	int old_nx = nx;
	nx = x;
	ny = y;
	nz = z;

	rdata = static_cast < float *>(realloc(rdata, x * y * z * sizeof(float)));
	update();

	attr_dict["nx"] = x;
	attr_dict["ny"] = y;
	attr_dict["nz"] = z;

	if (old_nx == 0) {
		memset(rdata, 0, x * y * z * sizeof(float));
	}

	if (supp) {
		free(supp);
		supp = 0;
	}
}

float *EMData::get_data() const
{
	return rdata;
}

void EMData::done_data()
{
	flags |= EMDATA_NEEDUPD;
}

Dict EMData::get_attr_dict()
{
	update_stat();
	return Dict(attr_dict);
}


void EMData::dump_data(string filename)
{
	string headerfile = filename + ".head";
	string datafile = filename + +".data";

	FILE *hfile = fopen(headerfile.c_str(), "wb");
	if (!hfile) {
		LOGERR("cannot open dump header file: %s", headerfile.c_str());
		return;
	}

	vector < string > keys = attr_dict.keys();
	for (size_t i = 0; i < keys.size(); i++) {
		fprintf(hfile, "%s = %s\n", keys[i].c_str(), attr_dict[keys[i]].to_str().c_str());
	}
	fclose(hfile);
	hfile = 0;

	FILE *dfile = fopen(datafile.c_str(), "wb");
	if (!dfile) {
		LOGERR("cannot open dump data file: %s", datafile.c_str());
		return;
	}

	size_t row_size = nx * sizeof(float);
	int nxy = nx * ny;

	for (int i = 0; i < nz; i++) {
		for (int j = 0; j < ny; j++) {
			fwrite(&rdata[i * nxy + j * nx], row_size, 1, dfile);
		}
	}
	fclose(dfile);
	dfile = 0;
}


void EMData::set_ctf(Ctf * new_ctf)
{
	if (!ctf) {
		ctf = new SimpleCtf();
	}

	ctf->copy_from(new_ctf);
}


vector < EMData * >EMData::read_images(string filename, vector < int >img_indices,
									   bool header_only) 
{
	int total_img = EMUtil::get_image_count(filename);
	size_t num_img = img_indices.size();

	for (size_t i = 0; i < num_img; i++) {
		assert(img_indices[i] >= 0 && img_indices[i] < total_img);
	}

	size_t n = num_img == 0 ? total_img : num_img;

	vector < EMData * >v;
	for (size_t j = 0; j < n; j++) {
		EMData *d = new EMData();
		size_t k = num_img == 0 ? j : img_indices[j];
		try {
			d->read_image(filename, k, header_only);
		}
		catch(...) {
			delete d;
			d = 0;
			throw;
		}

		v.push_back(d);
	}
	return v;
}


vector < EMData * >EMData::read_images_ext(string filename, int img_index_start,
										   int img_index_end, bool header_only,
										   string ext) 
{
	assert(img_index_end >= img_index_start);

	if (img_index_start < 0) {
		img_index_start = 0;
	}

	string new_filename = filename.insert(filename.rfind("."), ext);
	int num_img = EMUtil::get_image_count(new_filename);

	assert(img_index_start < num_img);

	if (img_index_end >= num_img) {
		img_index_end = num_img - 1;
	}

	vector < EMData * >v;

	for (int i = img_index_start; i < img_index_end; i++) {
		EMData *d = new EMData();
		try {
			d->read_image(new_filename, i, header_only);
		}
		catch(...) {
			delete d;
			d = 0;
			throw;
		}
		v.push_back(d);
	}
	return v;
}

EMData & EMData::operator+=(float n)
{
	add(n);
	update_stat();
	return *this;
}

EMData & EMData::operator-=(float n)
{
	*this += (-n);
	return *this;
}

EMData & EMData::operator*=(float n)
{
	mult(n);
	update_stat();
	return *this;
}

EMData & EMData::operator/=(float n)
{
	if (n == 0) {
		LOGERR("divided by zero");
		return *this;
	}
	*this *= (1.0f / n);
	return *this;
}


EMData & EMData::operator+=(const EMData & em)
{
	add(em);
	update_stat();	
	return *this;
}

EMData & EMData::operator-=(const EMData & em)
{
	sub(em);
	update_stat();	
	return *this;
}

EMData & EMData::operator*=(const EMData & em)
{
	mult(em);
	update_stat();
	return *this;
}

EMData & EMData::operator/=(const EMData & em)
{
	div(em);
	update_stat();	
	return *this;
}


EMData EMAN::operator+(const EMData & em, float n)
{
	EMData r = em;
	r += n;
	return r;
}

EMData EMAN::operator-(const EMData & em, float n)
{
	EMData r = em;
	r -= n;
	return r;
}

EMData EMAN::operator*(const EMData & em, float n)
{
	EMData r = em;
	r *= n;
	return r;
}

EMData EMAN::operator/(const EMData & em, float n)
{
	EMData r = em;
	r /= n;
	return r;
}


EMData EMAN::operator+(float n, const EMData & em)
{
	EMData r = em;
	r += n;
	return r;
}

EMData EMAN::operator-(float n, const EMData & em)
{
	EMData r = em;
	r -= n;
	return r;
}

EMData EMAN::operator*(float n, const EMData & em)
{
	EMData r = em;
	r *= n;
	return r;
}

EMData EMAN::operator/(float n, const EMData & em)
{
	EMData r = em;
	r /= n;
	return r;
}


EMData EMAN::operator+(const EMData & a, const EMData & b)
{
	EMData r = a;
	r += b;
	return r;
}

EMData EMAN::operator-(const EMData & a, const EMData & b)
{
	EMData r = a;
	r -= b;
	return r;
}

EMData EMAN::operator*(const EMData & a, const EMData & b)
{
	EMData r = a;
	r *= b;
	return r;
}

EMData EMAN::operator/(const EMData & a, const EMData & b)
{
	EMData r = a;
	r /= b;
	return r;
}

double EMData::dot_rotate_translate(EMData * with, float dx, float dy, float da)
{
	if (!EMUtil::is_same_size(this, with)) {
		LOGERR("dot_rotate_translate requires same size image");
		return 0;
	}

	if (nz != 1) {
		LOGERR("dot_rotate_translate can only handle 2D images");
		return 0;
	}

	float *this_data = 0;
	if (parent) {
		this_data = parent->get_data();
	}
	else {
		this_data = get_data();
	}

	float *with_data = with->get_data();
	float mx0 = cos(da);
	float mx1 = sin(da);
	float y = -ny / 2.0f;
	float my0 = mx0 * (-nx / 2.0f - 1.0f) + nx / 2.0f - dx;
	float my1 = -mx1 * (-nx / 2.0f - 1.0f) + ny / 2.0f - dy;
	double result = 0;

	for (int j = 0; j < ny; j++) {
		float x2 = my0 + mx1 * y;
		float y2 = my1 + mx0 * y;

		int ii = Util::fast_floor(x2);
		int jj = Util::fast_floor(y2);
		float t = x2 - ii;
		float u = y2 - jj;

		for (int i = 0; i < nx; i++) {
			t += mx0;
			u -= mx1;

			if (t >= 1.0f) {
				ii++;
				t -= 1.0f;
			}

			if (u >= 1.0f) {
				jj++;
				u -= 1.0f;
			}

			if (t < 0) {
				ii--;
				t += 1.0f;
			}

			if (u < 0) {
				jj--;
				u += 1.0f;
			}

			if (ii >= 0 && ii <= nx - 2 && jj >= 0 && jj <= ny - 2) {
				int k0 = ii + jj * nx;
				int k1 = k0 + 1;
				int k2 = k0 + nx + 1;
				int k3 = k0 + nx;

				float tt = 1 - t;
				float uu = 1 - u;

				result += (this_data[k0] * tt * uu + this_data[k1] * t * uu +
						   this_data[k2] * t * u + this_data[k3] * tt * u) * with_data[i + j * nx];
			}
		}
		y += 1.0f;
	}

	return result;
}

void EMData::rotate_x(int dx)
{
	float *tmp = new float[nx];
	size_t float_size = sizeof(float);

	for (int y = 0; y < ny; y++) {
		int y_nx = y * nx;
		for (int x = 0; x < nx; x++) {
			tmp[x] = rdata[y_nx + (x + dx) % nx];
		}
		memcpy(&rdata[y * nx], tmp, nx * float_size);
	}

	done_data();
	delete[]tmp;
	tmp = 0;
}

void EMData::set_xyz_origin(float origin_x, float origin_y, float origin_z)
{
	attr_dict["origin_row"] = origin_x;
	attr_dict["origin_col"] = origin_y;
	attr_dict["origin_sec"] = origin_z;
}

float *EMData::setup4slice(bool redo)
{
	if (!is_complex() || nz == 1) {
		return 0;
	}

	if (supp) {
		if (redo) {
			free(supp);
		}
		else {
			return supp;
		}
	}

	const int SUPP_ROW_SIZE = 8;
	const int SUPP_ROW_OFFSET = 4;
	const int supp_size = SUPP_ROW_SIZE + SUPP_ROW_OFFSET;

	supp = (float *) calloc(supp_size, ny * nz * sizeof(float));
	int nxy = nx * ny;
	int supp_xy = supp_size * ny;

	for (int z = 0; z < nz; z++) {
		int cur_z1 = z * nxy;
		int cur_z2 = z * supp_xy;

		for (int y = 0; y < ny; y++) {
			int cur_y1 = y * nx;
			int cur_y2 = y * supp_size;

			for (int x = 0; x < SUPP_ROW_SIZE; x++) {
				int k = (x + SUPP_ROW_OFFSET) + cur_y2 + cur_z2;
				supp[k] = rdata[x + cur_y1 + cur_z1];
			}
		}
	}

	for (int z = 1, zz = nz - 1; z < nz; z++, zz--) {
		int cur_z1 = zz * nxy;
		int cur_z2 = z * supp_xy;

		for (int y = 1, yy = ny - 1; y < ny; y++, yy--) {
			supp[y * 12 + cur_z2] = rdata[4 + yy * nx + cur_z1];
			supp[1 + y * 12 + cur_z2] = -rdata[5 + yy * nx + cur_z1];
			supp[2 + y * 12 + cur_z2] = rdata[2 + yy * nx + cur_z1];
			supp[3 + y * 12 + cur_z2] = -rdata[3 + yy * nx + cur_z1];
		}
	}

	return supp;
}

void EMData::to_zero()
{
	memset(rdata, 0, nx * ny * nz * sizeof(float));
	done_data();
}

void EMData::scale(float s)
{
	Transform t;
	t.set_scale_instance(Vec3 < float >(s, 1, 1));
	rotate_translate(t);
}

void EMData::translate(float dx, float dy, float dz)
{
	translate(Vec3 < float >(dx, dy, dz));
}


void EMData::translate(const Vec3 < float >&translation)
{
	if (nz != 1) {
		LOGERR("3D translation not supported yet");
		return;
	}

	float *parent_data = get_data();
	float *this_data = get_data();

	int x0 = nx - 1;
	int x1 = -1;
	int x2 = -1;

	if (translation[0] < 0) {
		x0 = 0;
		x1 = nx;
		x2 = 1;
	}

	int y0 = ny - 1;
	int y1 = -1;
	int y2 = -1;

	if (translation[1] < 0) {
		y0 = 0;
		y1 = ny;
		y2 = 1;
	}

	for (int x = x0; x != x1; x += x2) {
		for (int y = y0; y != y1; y += y2) {
			int xp = static_cast < int >(x - translation[0]);
			int yp = static_cast < int >(y - translation[1]);

			if (xp < 0 || yp < 0 || xp >= nx || yp >= ny) {
				this_data[x + y * nx] = 0;
			}
			else {
				this_data[x + y * nx] = parent_data[xp + yp * nx];
			}
		}
	}

	done_data();
	all_translation += translation;
}

/** slow; need to fix later
 */
void EMData::rotate(float alt, float az, float phi)
{
	Transform t(Rotation(alt, az, phi, Rotation::EMAN));
	rotate_translate(t);
}

void EMData::rotate(const Rotation & r)
{
	Transform t(r);
	rotate_translate(t);
}

void EMData::rotate_translate(float alt, float az, float phi, float dx, float dy, float dz)
{
	rotate_translate(Rotation(alt, az, phi, Rotation::EMAN), Vec3 < float >(dx, dy, dz));
}

void EMData::rotate_translate(const Rotation & rotation, const Vec3 < float >&translation)
{
	float *src_data = 0;
	float *des_data = 0;

	if (parent) {
		src_data = parent->get_data();
		des_data = get_data();
	}
	else {
		src_data = get_data();
		des_data = (float *) malloc(nx * ny * nz * sizeof(float));
	}

	if (nz == 1) {
		float mx[2];

		mx[0] = cos(rotation.eman_alt());
		mx[1] = sin(rotation.eman_alt());
		float dxx = mx[0];
		float dyx = -mx[1];
		float dxy = mx[1];
		float dyy = mx[0];

		int k = 0;
		float x = -mx[1] * ny / 2.0f + nx / 2.0f - translation[0];
		float y = -mx[0] * ny / 2.0f + ny / 2.0f - translation[1];

		for (int j = 0; j < ny; j++) {
			float x2 = x - mx[0] * nx / 2.0f;
			float y2 = y + mx[1] * nx / 2.0f;

			for (int i = 0; i < nx; i++) {
				if (x2 < 0 || x2 > nx - 2.0f || y2 < 0 || y2 > ny - 2.0f) {
					des_data[k] = 0;
				}
				else {
					int ii = (int) floor(x2 + 0.5f);
					int jj = (int) floor(y2 + 0.5f);
					des_data[k] = src_data[ii + jj * nx];
				}
				k++;
				x2 += dxx;
				y2 += dyx;
			}

			x += dxy;
			y += dyy;
		}
	}
	else {
		Matrix3f mx = rotation.get_matrix3();

		float x4 = (mx[0][0] * (-nx / 2.0f) + mx[0][1] * (-ny / 2.0f) +
					mx[0][2] * (-nz / 2.0f)) + nx / 2.0f - translation[0];
		float y4 = (mx[1][0] * (-nx / 2.0f) + mx[1][1] * (-ny / 2.0f) +
					mx[1][2] * (-nz / 2.0f)) + ny / 2.0f - translation[1];
		float z4 = (mx[2][0] * (-nx / 2.0f) + mx[2][1] * (-ny / 2.0f) +
					mx[2][2] * (-nz / 2.0f)) + nz / 2.0f - translation[2];

		int xy = nx * ny;
		int mr = 0;

		if (nx < ny && nx < nz) {
			mr = (nx - 2) * (nx - 2) / 4;
		}
		else if (ny < nz) {
			mr = (ny - 2) * (ny - 2) / 4;
		}
		else {
			mr = (nz - 2) * (nz - 2) / 4;
		}

		int l = 0;
		for (int k = -nz / 2; k < nz / 2; k++) {
			float x3 = x4;
			float y3 = y4;
			float z3 = z4;

			for (int j = -ny / 2; j < ny / 2; j++) {
				float x2 = x3;
				float y2 = y3;
				float z2 = z3;
				for (int i = -nx / 2; i < nx / 2; i++) {
					if (i * i + j * j + k * k >= mr) {
						des_data[l] = 0;
					}
					else {
						float x = floor(x2 + 0.5f);
						float y = floor(y2 + 0.5f);
						float z = floor(z2 + 0.5f);

						int ii = (int) (x + y * nx + z * xy);
						des_data[l] = src_data[ii];
					}
					l++;
					x2 += mx[0][0];
					y2 += mx[1][0];
					z2 += mx[2][0];
				}

				x3 += mx[0][1];
				y3 += mx[1][1];
				z3 += mx[2][1];
			}

			x4 += mx[0][2];
			y4 += mx[1][2];
			z4 += mx[2][2];
		}
	}

	if (parent) {
		parent->done_data();
	}
	else {
		free(rdata);
		rdata = des_data;
	}

	done_data();

	attr_dict["origin_row"] = attr_dict["origin_row"];
	attr_dict["origin_col"] = attr_dict["origin_col"];
	attr_dict["origin_sec"] = attr_dict["origin_sec"];

	update();

	all_rotation += Vec3 < float >(rotation.eman_alt(), rotation.eman_az(), rotation.eman_phi());
	all_translation += translation;
}



void EMData::rotate_translate(const Transform & xform)
{
	float scale = xform.get_scale()[0];
	Vec3 < float >dcenter = xform.get_center();
	float dxc = dcenter[0];
	float dyc = dcenter[1];
	float dzc = dcenter[2];

	Vec3 < float >translation = xform.get_post_translate();
	Rotation rotation = xform.get_rotation();

	int nx2 = nx;
	int ny2 = ny;
	float inv_scale = 1.0f;

	if (scale != 0) {
		inv_scale = 1.0f / scale;
	}

	float *src_data = 0;
	float *des_data = 0;

	if (parent) {
		src_data = parent->get_data();
		des_data = get_data();
		nx2 = parent->get_xsize();
		ny2 = parent->get_ysize();
	}
	else {
		src_data = get_data();
		des_data = (float *) malloc(nx * ny * nz * sizeof(float));
	}

	if (nz == 1) {
		float mx0 = inv_scale * cos(rotation.eman_alt());
		float mx1 = inv_scale * sin(rotation.eman_alt());

		float x2c = nx / 2.0f - dxc - translation[0];
		float y2c = ny / 2.0f - dyc - translation[1];
		float y = -ny / 2.0f + dxc;

		for (int j = 0; j < ny; j++, y += 1.0f) {
			float x = -nx / 2.0f + dyc;

			for (int i = 0; i < nx; i++, x += 1.0f) {
				float x2 = (mx0 * x + mx1 * y) + x2c;
				float y2 = (-mx1 * x + mx0 * y) + y2c;

				if (x2 < 0 || x2 > nx2 - 1 || y2 < 0 || y2 > ny2 - 1) {
					des_data[i + j * nx] = 0;
				}
				else {
					int ii = Util::fast_floor(x2);
					int jj = Util::fast_floor(y2);
					int k0 = ii + jj * nx;
					int k1 = k0 + 1;
					int k2 = k0 + nx + 1;
					int k3 = k0 + nx;

					if (ii == nx2 - 1) {
						k1--;
						k2--;
					}
					if (jj == ny2 - 1) {
						k2 -= nx2;
						k3 -= nx2;
					}

					float t = x2 - ii;
					float u = y2 - jj;
					float tt = 1 - t;
					float uu = 1 - u;

					float p0 = src_data[k0] * tt * uu;
					float p1 = src_data[k1] * t * uu;
					float p3 = src_data[k3] * tt * u;
					float p2 = src_data[k2] * t * u;

					des_data[i + j * nx] = p0 + p1 + p2 + p3;
				}
			}
		}
	}

	else if (nx == (nx / 2 * 2 + 1) && nx == ny && (2 * nz - 1) == nx) {
		// make sure this is right
		Matrix3f mx = rotation.get_matrix3();
		mx *= inv_scale;
		int nxy = nx * ny;
		int l = 0;

		for (int k = 0; k < nz; k++) {
			for (int j = -ny / 2; j < ny - ny / 2; j++) {
				for (int i = -nx / 2; i < nx - nx / 2; i++, l++) {
					float x2 = mx[0][0] * i + mx[0][1] * j + mx[0][2] * k + nx / 2;
					float y2 = mx[1][0] * i + mx[1][1] * j + mx[1][2] * k + ny / 2;
					float z2 = mx[2][0] * i + mx[2][1] * j + mx[2][2] * k + 0 / 2;	// 0/2?

					if (x2 >= 0 && y2 >= 0 && z2 > -(nz - 1) && x2 < nx && y2 < ny && z2 < nz - 1) {
						if (z2 < 0) {
							x2 = nx - 1 - x2;
							z2 = -z2;
						}

						int x = Util::fast_floor(x2);
						int y = Util::fast_floor(y2);
						int z = Util::fast_floor(z2);

						float t = x2 - x;
						float u = y2 - y;
						float v = z2 - z;

						int ii = (int) (x + y * nx + z * nxy);

						des_data[l] += Util::trilinear_interpolate(src_data[ii], src_data[ii + 1],
																   src_data[ii + nx],
																   src_data[ii + nx + 1],
																   src_data[ii + nx * ny],
																   src_data[ii + nxy + 1],
																   src_data[ii + nxy + nx],
																   src_data[ii + nxy + nx + 1], t,
																   u, v);
					}
				}
			}
		}
	}
	else {
		Matrix3f mx = rotation.get_matrix3();
		mx *= inv_scale;

		float x4 = (mx[0][0] * (-nx / 2.0f + dxc) + mx[0][1] * (-ny / 2.0f + dyc) +
					mx[0][2] * (-nz / 2.0f + dzc)) + nx / 2.0f - dxc - translation[0];

		float y4 = (mx[1][3] * (-nx / 2.0f + dxc) + mx[1][4] * (-ny / 2.0f + dyc) +
					mx[1][5] * (-nz / 2.0f + dzc)) + ny / 2.0f - dyc - translation[1];

		float z4 = (mx[2][6] * (-nx / 2.0f + dxc) + mx[2][7] * (-ny / 2.0f + dyc) +
					mx[2][8] * (-nz / 2.0f + dzc)) + nz / 2.0f - dzc - translation[2];

		int nxy = nx * ny;
		int l = 0;

		for (int k = -nz / 2; k < nz / 2; k++) {
			float x3 = x4;
			float y3 = y4;
			float z3 = z4;

			for (int j = -ny / 2; j < ny / 2; j++) {
				float x2 = x3;
				float y2 = y3;
				float z2 = z3;

				for (int i = -nx / 2; i < nx / 2; i++, l++) {

					if (x2 < 0 || y2 < 0 || z2 < 0 || x2 >= nx - 1 || y2 >= ny - 1 || z2 >= nz - 1) {
						des_data[l] = 0;
					}
					else {
						int x = Util::fast_floor(x2);
						int y = Util::fast_floor(y2);
						int z = Util::fast_floor(z2);

						float t = x2 - x;
						float u = y2 - y;
						float v = z2 - z;
						int ii = (int) (x + y * nx + z * nxy);

						des_data[l] = Util::trilinear_interpolate(src_data[ii], src_data[ii + 1],
																  src_data[ii + nx],
																  src_data[ii + nx + 1],
																  src_data[ii + nx * ny],
																  src_data[ii + nxy + 1],
																  src_data[ii + nxy + nx],
																  src_data[ii + nxy + nx + 1], t, u,
																  v);
					}

					x2 += mx[0][0];
					y2 += mx[1][1];
					z2 += mx[2][0];
				}

				x3 += mx[0][1];
				y3 += mx[1][1];
				z3 += mx[2][1];
			}

			x4 += mx[0][2];
			y4 += mx[1][2];
			z4 += mx[2][2];
		}


		if (parent) {
			parent->done_data();
		}
		else {
			free(rdata);
			rdata = des_data;
		}

		scale_pixel(inv_scale);

		attr_dict["origin_row"] = (float) attr_dict["origin_row"] * inv_scale;
		attr_dict["origin_col"] = (float) attr_dict["origin_col"] * inv_scale;
		attr_dict["origin_sec"] = (float) attr_dict["origin_sec"] * inv_scale;

		done_data();
		update();
	}


	all_rotation += Vec3 < float >(rotation.eman_alt(), rotation.eman_az(), rotation.eman_phi());
	all_translation += translation;
}


void EMData::rotate_180()
{
	if (nx != ny) {
		throw NotSquareImageError(__FILE__, __LINE__, 
								  "Rot180 on non-square image poorly defined.");
	}

	float *d = get_data();

	for (int x = 1; x < nx; x++) {
		int y = 0;
		for (y = 1; y < ny; y++) {
			if (x == nx / 2 && y == ny / 2) {
				break;
			}
			int i = x + y * nx;
			int k = nx - x + (ny - y) * nx;

			float t = d[i];
			d[i] = d[k];
			d[k] = t;
		}
		if (x == nx / 2 && y == ny / 2) {
			break;
		}
	}

	done_data();
	all_rotation[0] += M_PI;
}

EMData *EMData::do_radon()
{
	if (nx != ny || nz != 1) {
		LOGERR("Radon transform only works on square 2D images");
		return 0;
	}

	EMData *result = new EMData();
	result->set_size(nx, ny, 1);
	result->to_zero();
	float *result_data = result->get_data();

	EMData *this_copy = this;
	if (!parent) {
		this_copy = copy();
	}

	for (int i = 0; i < nx; i++) {
		this_copy->rotate(M_PI * 2.0f * i / nx, 0, 0);

		float *copy_data = this_copy->get_data();

		for (int y = 0; y < nx; y++) {
			for (int x = 0; x < nx; x++) {
				if (Util::square(x - nx / 2) + Util::square(y - nx / 2) <= nx * nx / 4) {
					result_data[i + y * nx] += copy_data[x + y * nx];
				}
			}
		}

		this_copy->done_data();
	}

	result->done_data();

	if (!parent) {
		delete this_copy;
		this_copy = 0;
	}

	return result;

}


void EMData::mean_shrink(int shrink_factor)
{
	if (shrink_factor <= 1) {
		throw ShrinkFactorError(shrink_factor, __FILE__, __LINE__,
								 "mean shrink: shrink factor must > 1");
	}

	if ((nx % shrink_factor != 0) || (ny % shrink_factor != 0) ||
		(nz > 1 && (nz % shrink_factor != 0))) {
		throw ShrinkFactorError(shrink_factor, __FILE__, __LINE__,
								 "Image size not divisible by shrink factor");
	}

	int shrinked_nx = nx / shrink_factor;
	int shrinked_ny = ny / shrink_factor;
	int shrinked_nz = 1;

	int threed_shrink_factor = shrink_factor * shrink_factor;
	int z_shrink_factor = 1;

	if (nz > 1) {
		shrinked_nz = nz / shrink_factor;
		threed_shrink_factor *= shrink_factor;
		z_shrink_factor = shrink_factor;
	}

	float *data = get_data();
	int nxy = nx * ny;
	int shrinked_nxy = shrinked_nx * shrinked_ny;

	for (int k = 0; k < shrinked_nz; k++) {
		int k_min = k * shrink_factor;
		int k_max = k * shrink_factor + z_shrink_factor;
		int cur_k = k * shrinked_nxy;

		for (int j = 0; j < shrinked_ny; j++) {
			int j_min = j * shrink_factor;
			int j_max = j * shrink_factor + shrink_factor;
			int cur_j = j * shrinked_nx + cur_k;

			for (int i = 0; i < shrinked_nx; i++) {
				int i_min = i * shrink_factor;
				int i_max = i * shrink_factor + shrink_factor;

				float sum = 0;
				for (int kk = k_min; kk < k_max; kk++) {
					int cur_kk = kk * nxy;

					for (int jj = j_min; jj < j_max; jj++) {
						int cur_jj = jj * nx + cur_kk;
						for (int ii = i_min; ii < i_max; ii++) {
							sum += data[ii + cur_jj];
						}
					}
				}
				data[i + cur_j] = sum / threed_shrink_factor;
			}
		}
	}

	done_data();
	set_size(shrinked_nx, shrinked_ny, shrinked_nz);
	scale_pixel(shrink_factor);
}

void EMData::median_shrink(int shrink_factor)
{
	if (shrink_factor <= 1) {
		throw ShrinkFactorError(shrink_factor, __FILE__, __LINE__,
								 "median shrink: shrink factor must > 1");
	}

	if ((nx % shrink_factor != 0) || (ny % shrink_factor != 0) ||
		(nz > 1 && (nz % shrink_factor != 0))) {
		throw ShrinkFactorError(shrink_factor, __FILE__, __LINE__,
								 "Image size not divisible by shrink factor");
	}

	int threed_shrink_factor = shrink_factor * shrink_factor;
	int size = nx * ny;
	int nx_old = nx;
	int ny_old = ny;

	int shrinked_nx = nx / shrink_factor;
	int shrinked_ny = ny / shrink_factor;
	int shrinked_nz = 1;

	int z_shrink_factor = 1;

	if (nz > 1) {
		threed_shrink_factor *= shrink_factor;
		size *= nz;
		shrinked_nz = nz / shrink_factor;
		z_shrink_factor = shrink_factor;
	}

	float *mbuf = new float[threed_shrink_factor];
	float *data_copy = new float[size];

	memcpy(data_copy, get_data(), size * sizeof(float));
	set_size(shrinked_nx, shrinked_ny, shrinked_nz);
	scale_pixel(shrink_factor);

	get_data();

	int nxy_old = nx_old * ny_old;
	int nxy_new = nx * ny;

	for (int l = 0; l < nz; l++) {
		int l_min = l * shrink_factor;
		int l_max = l * shrink_factor + z_shrink_factor;
		int cur_l = l * nxy_new;

		for (int j = 0; j < ny; j++) {
			int j_min = j * shrink_factor;
			int j_max = (j + 1) * shrink_factor;
			int cur_j = j * nx + cur_l;

			for (int i = 0; i < nx; i++) {
				int i_min = i * shrink_factor;
				int i_max = (i + 1) * shrink_factor;

				int k = 0;
				for (int l2 = l_min; l2 < l_max; l2++) {
					int cur_l2 = l2 * nxy_old;

					for (int j2 = j_min; j2 < j_max; j2++) {
						int cur_j2 = j2 * nx_old + cur_l2;

						for (int i2 = i_min; i2 < i_max; i2++) {
							mbuf[k] = data_copy[i2 + cur_j2];
							k++;
						}
					}
				}


				for (k = 0; k < threed_shrink_factor / 2 + 1; k++) {
					for (int i2 = k + 1; i2 < threed_shrink_factor; i2++) {
						if (mbuf[i2] < mbuf[k]) {
							float f = mbuf[i2];
							mbuf[i2] = mbuf[k];
							mbuf[k] = f;
						}
					}
				}

				rdata[i + cur_j] = mbuf[threed_shrink_factor / 2];
			}
		}
	}

	done_data();

	delete[]data_copy;
	data_copy = 0;

	delete[]mbuf;
	mbuf = 0;
}

IntPoint EMData::calc_min_location() const
{
	int di = 1;
	if (is_complex() && !is_ri()) {
		di = 2;
	}

	float min = FLT_MAX;
	int min_x = 0;
	int min_y = 0;
	int min_z = 0;
	int nxy = nx * ny;

	for (int j = 0; j < nz; j++) {
		int cur_z = j * nxy;

		for (int k = 0; k < ny; k++) {
			int cur_y = k * nx + cur_z;

			for (int l = 0; l < nx; l += di) {
				float t = rdata[l + cur_y];
				if (t < min) {
					min_x = l;
					min_y = k;
					min_z = j;
					min = t;
				}
			}
		}
	}

	return IntPoint(min_x, min_y, min_z);
}

IntPoint EMData::calc_max_location() const
{
	int di = 1;
	if (is_complex() && !is_ri()) {
		di = 2;
	}

	float max = -FLT_MAX;
	int max_x = 0;
	int max_y = 0;
	int max_z = 0;
	int nxy = nx * ny;

	for (int j = 0; j < nz; j++) {
		int cur_z = j * nxy;

		for (int k = 0; k < ny; k++) {
			int cur_y = k * nx + cur_z;

			for (int l = 0; l < nx; l += di) {
				float t = rdata[l + cur_y];
				if (t > max) {
					max_x = l;
					max_y = k;
					max_z = j;
					max = t;
				}
			}
		}
	}

	return IntPoint(max_x, max_y, max_z);
}



int EMData::calc_min_index() const
{
	IntPoint min_location = calc_min_location();
	int i = min_location.x + min_location.y * nx + min_location.z * nx * ny;
	return i;
}


int EMData::calc_max_index() const
{
	IntPoint max_location = calc_max_location();
	int i = max_location.x + max_location.y * nx + max_location.z * nx * ny;
	return i;
}

EMData *EMData::calc_ccfx(EMData * with, int y0, int y1, bool no_sum)
{
	if (!with) {
		LOGERR("Empty 'with' image. calc_ccfx requires 'with'.");
		return 0;
	}

	if (!EMUtil::is_same_size(this, with)) {
		LOGERR("calc_ccfx must have 2 images of the same size");
		return 0;
	}
	if (nz > 1) {
		LOGERR("calc_ccfx only works on 2D images");
		return 0;
	}

	EMData *cf = new EMData();
	if (no_sum) {
		cf->set_size(nx, y1 - y0 + 1, 1);
	}
	else {
		cf->set_size(nx, 1, 1);
	}

	cf->set_name("CCFx");
	cf->set_path("/tmp/eman.ccf");

	if (y1 <= y0) {
		y1 = ny;
	}

	if (y0 >= y1) {
		y0 = 0;
	}

	if (y0 < 0) {
		y0 = 0;
	}

	if (y1 > ny) {
		y1 = ny;
	}

	if (no_sum) {
		float *cfd = cf->get_data();
		float *with_data = with->get_data();

		for (int y = y0; y < y1; y++) {
			int cur_y = y * nx;

			for (int x = 0; x < nx; x++) {
				float dot = 0;
				for (int i = 0; i < nx; i++) {
					int k1 = (i + x) % nx + cur_y;
					dot += rdata[i + cur_y] * with_data[k1];
				}
				cfd[x + (y - y0) * nx] = dot;
			}
		}

		cf->done_data();
		return cf;
	}
	else {
		float *f1 = (float *) calloc(nx * sizeof(float), 1);
		float *f2 = (float *) calloc(nx * sizeof(float), 1);

		float *cfd = cf->get_data();
		float *d1 = get_data();
		float *d2 = with->get_data();
		size_t row_size = nx * sizeof(float);

		if (!is_complex_x()) {
			for (int j = 0; j < ny; j++) {
				EMfft::real_to_complex_1d(d1 + j * nx, f1, nx);
				memcpy(d1 + j * nx, f1, row_size);
			}

			set_complex_x(true);
		}
		if (!with->is_complex_x()) {
			for (int j = 0; j < with->get_ysize(); j++) {
				EMfft::real_to_complex_1d(d2 + j * nx, f2, nx);
				memcpy(d2 + j * nx, f2, row_size);
			}

			with->set_complex_x(true);
		}

		for (int j = y0; j < y1; j++) {
			float *f1a = d1 + j * nx;
			float *f2a = d2 + j * nx;

			f1[0] = f1a[0] * f2a[0];
			f1[nx / 2] = f1a[nx / 2] * f2a[nx / 2];

			for (int i = 1; i < nx / 2; i++) {
				float re1 = f1a[i];
				float re2 = f2a[i];
				float im1 = f1a[nx - i];
				float im2 = f2a[nx - i];

				f1[i] = re1 * re2 + im1 * im2;
				f1[nx - i] = im1 * re2 - re1 * im2;
			}

			EMfft::complex_to_real_1d(f1, f2, nx);

			if (no_sum) {
				for (int i = 0; i < nx; i++) {
					cfd[i + nx * (j - y0)] = f2[i];
				}
			}
			else {
				for (int i = 0; i < nx; i++) {
					cfd[i] += f2[i];
				}
			}
		}

		free(f1);
		f1 = 0;
		free(f2);
		f2 = 0;
	}

	cf->done_data();
	done_data();
	with->done_data();


	return cf;
}

#if 0
void EMData::calc_rcf(EMData * with, vector < float >&sum_array)
{
	int array_size = sum_array.size();
	float da = 2 * M_PI / array_size;
	float *dat = new float[array_size + 2];
	float *dat2 = new float[array_size + 2];
	int nx2 = nx * 9 / 20;

	float lim = 0;
	if (fabs(translation[0]) < fabs(translation[1])) {
		lim = fabs(translation[1]);
	}
	else {
		lim = fabs(translation[0]);
	}

	nx2 -= static_cast < int >(floor(lim));

	for (int i = 0; i < array_size; i++) {
		sum_array[i] = 0;
	}

	float sigma = attr_dict["sigma"];
	float with_sigma = with->get_attr_dict().get("sigma");

	for (int i = 8; i < nx2; i += 6) {
		calc_az_dist(array_size, 0, da, dat, i, i + 6);
		with->calc_az_dist(array_size, 0, da, dat2, i, i + 6);

		EMfft::real_to_complex_1d(dat, dat, array_size);
		EMfft::real_to_complex_1d(dat2, dat2, array_size);

		for (int j = 0; j < array_size + 2; j += 2) {
			float max = dat[j] * dat2[j] + dat[j + 1] * dat2[j + 1];
			float max2 = dat[j + 1] * dat2[j] - dat2[j + 1] * dat[j];
			dat[j] = max;
			dat[j + 1] = max2;
		}

		EMfft::complex_to_real_1d(dat, dat, array_size);
		float norm = array_size * array_size * (4.0f * sigma) * (4.0f * with_sigma);

		for (int j = 0; j < array_size; j++) {
			sum_array[j] += dat[j] * (float) i / norm;
		}
	}

	delete[]dat;
	dat = 0;

	delete[]dat2;
	dat2 = 0;
}

#endif


void EMData::to_one()
{
	if (is_complex()) {
		set_ri(true);
	}

	for (int i = 0; i < nx * ny * nz; i++) {
		rdata[i] = 1.0f;
	}

	update();
}




EMData *EMData::calc_ccf(EMData * with, bool tocorner, EMData * filter)
{
	EMData *f1 = 0;

	if (filter) {
		f1 = do_fft()->copy();
	}
	else {
		f1 = do_fft();
	}

	if (!f1) {
		return 0;
	}

	EMData *with_fft = 0;
	if (with) {
		with_fft = with->do_fft();
		if (!with_fft) {
			return 0;
		}
	}

	f1->ap2ri();
	EMData *cf = 0;

	if (with) {
		cf = with_fft->copy(false);
		cf->ap2ri();
	}
	else {
		cf = f1->copy(false);
	}

	if (filter) {
		if (!EMUtil::is_same_size(filter, cf)) {
			LOGERR("calc_ccf: cannot apply improperly sized filter!");
			return 0;
		}

		cf->mult(*filter);
		f1->mult(*filter);
	}

	if (with && !EMUtil::is_same_size(f1, cf)) {
		LOGERR("calc_ccf: Images must be the same size.");
		return 0;
	}

	float *rdata1 = f1->get_data();
	float *rdata2 = cf->get_data();

	int cf_size = cf->get_xsize() * cf->get_ysize() * cf->get_zsize();

	if (with == this) {
		for (int i = 0; i < cf_size; i += 2) {
			rdata2[i] = rdata1[i] * rdata2[i] + rdata1[i + 1] * rdata2[i + 1];
		}
	}
	else if (with) {
		for (int i = 0; i < cf_size; i += 2) {
			rdata2[i] = rdata1[i] * rdata2[i] + rdata1[i + 1] * rdata2[i + 1];
			rdata2[i + 1] = rdata1[i + 1] * rdata2[i] - rdata1[i] * rdata2[i + 1];
		}
	}
	else {
		for (int i = 0; i < cf_size; i += 2) {
			rdata2[i] = rdata1[i] * rdata2[i] - rdata1[i + 1] * rdata2[i + 1];
			rdata2[i + 1] = rdata1[i + 1] * rdata2[i] + rdata1[i] * rdata2[i + 1];
		}
	}

	f1->done_data();

	if (tocorner) {
		cf->filter("Phase180");
	}

	EMData *f2 = cf->do_ift();

	delete cf;
	cf = 0;

	if (filter) {
		delete f1;
		f1 = 0;
	}

	f2->set_name("CCF");
	f2->set_path("/tmp/eman.ccf");

	return f2;
}

EMData *EMData::make_rotational_footprint(bool premasked, bool unwrap)
{
	static EMData *filt = 0;

	if (!(flags & EMDATA_NEWRFP) && rfp) {
		return rfp;
	}

	if (nx & 1) {
		return 0;
	}

	if (rfp) {
		delete rfp;
		rfp = 0;
	}

	if (!filt) {
		filt = new EMData();
		filt->set_complex(true);
	}

	int cs = (((nx * 7 / 4) & 0xfffff8) - nx) / 2;

	EMData *tmp2 = 0;
	Region r1;
	if (nz == 1) {
		r1 = Region(-cs, -cs, nx + 2 * cs, ny + 2 * cs);
	}
	else {
		r1 = Region(-cs, -cs, -cs, nx + 2 * cs, ny + 2 * cs, nz + 2 * cs);
	}
	tmp2 = get_clip(r1);

	if (!premasked) {
		tmp2->filter("SharpMask", Dict("outer_radius", nx / 2, "value", 0));
	}

	if (filt->get_xsize() != tmp2->get_xsize() + 2 || filt->get_ysize() != tmp2->get_ysize() ||
		filt->get_zsize() != tmp2->get_zsize()) {
		filt->set_size(tmp2->get_xsize() + 2, tmp2->get_ysize(), tmp2->get_zsize());
		filt->to_one();

		filt->filter("GaussLowpass", Dict("lowpass", 3));
		filt->filter("GaussHighpass", Dict("highpass", 9999));
	}

	EMData *tmp = tmp2->calc_mutual_correlation(tmp2, true, filt);
	delete tmp2;
	tmp2 = 0;

	Region r2;
	if (nz == 1) {
		r2 = Region(cs - nx / 4, cs - ny / 4, nx * 3 / 2, ny * 3 / 2);
	}
	else {
		r2 = Region(cs - nx / 4, cs - ny / 4, cs - nz / 4, nx * 3 / 2, ny * 3 / 2, nz * 3 / 2);
	}
	tmp2 = tmp->get_clip(r2);
	rfp = tmp2;
	flags &= ~EMDATA_NEWRFP;
	delete tmp;
	tmp = 0;

	if (nz == 1) {
		if (!unwrap) {
			tmp2->filter("SharpMask", Dict("outer_radius", -1, "value", 0));
			rfp = 0;
			return tmp2;
		}

		rfp = tmp2->unwrap();
		delete tmp2;
		tmp2 = 0;
	}

	return rfp;
}


EMData *EMData::calc_mutual_correlation(EMData * with, bool tocorner, EMData * filter)
{
	if (with && !EMUtil::is_same_size(this, with)) {
		LOGERR("to calculate mutual correlation, Images must be the same size.");
		return 0;
	}

	EMData *this_fft = 0;
	if (filter) {
		this_fft = do_fft()->copy();
	}
	else {
		this_fft = do_fft();
	}

	if (!this_fft) {
		return 0;
	}

	this_fft->ap2ri();
	EMData *cf = 0;

	if (with) {
		EMData *cf1 = with->do_fft();
		if (!cf1) {
			return 0;
		}
		cf = cf1->copy(false);
		cf->ap2ri();
	}
	else {
		cf = this_fft->copy(false);
	}

	if (filter) {
		if (!EMUtil::is_same_size(filter, cf)) {
			LOGERR("cannot apply improperly sized filter");
			return 0;
		}

		cf->mult(*filter);
		this_fft->mult(*filter);
	}

	float *rdata1 = this_fft->get_data();
	float *rdata2 = cf->get_data();
	int this_fft_size = this_fft->get_xsize() * this_fft->get_ysize() * this_fft->get_zsize();

	if (with == this) {
		for (int i = 0; i < this_fft_size; i += 2) {
			rdata2[i] = sqrt(rdata1[i] * rdata2[i] + rdata1[i + 1] * rdata2[i + 1]);
			rdata2[i + 1] = 0;
		}

		this_fft->done_data();
		cf->done_data();
	}
	else {
		for (int i = 0; i < this_fft_size; i += 2) {
			rdata2[i] = (rdata1[i] * rdata2[i] + rdata1[i + 1] * rdata2[i + 1]);
			rdata2[i + 1] = (rdata1[i + 1] * rdata2[i] - rdata1[i] * rdata2[i + 1]);
		}

		this_fft->done_data();
		cf->done_data();
		rdata1 = cf->get_data();

		for (int i = 0; i < this_fft_size; i += 2) {
			float t = Util::square(rdata1[i]) + Util::square(rdata1[i + 1]);
			if (t != 0) {
				t = pow(t, (float) 0.25);
				rdata1[i] /= t;
				rdata1[i + 1] /= t;
			}
		}
		cf->done_data();
	}

	if (tocorner) {
		cf->filter("Phase180");
	}

	EMData *f2 = cf->do_ift();

	delete cf;
	cf = 0;

	if (filter) {
		delete this_fft;
		this_fft = 0;
	}

	f2->set_name("MCF");
	f2->set_path("/tmp/eman.mcf");

	return f2;
}

EMData *EMData::unwrap(int r1, int r2, int xs, int dx, int dy, bool do360)
{
	int p = 1;
	if (do360) {
		p = 2;
	}

	EMData *ret = new EMData();

	if (xs < 1) {
		xs = (int) floor(p * M_PI * ny / 4);
		xs = Util::calc_best_fft_size(xs);
	}

	if (r1 < 0) {
		r1 = 4;
	}

	int rr = ny / 2 - 2 - (int) floor(hypot(dx, dy));
	if (r2 <= r1 || r2 > rr) {
		r2 = rr;
	}

	ret->set_size(xs, r2 - r1, 1);
	float *d = get_data();
	float *dd = ret->get_data();

	for (int x = 0; x < xs; x++) {
		float si = sin(x * M_PI * p / xs);
		float co = cos(x * M_PI * p / xs);

		for (int y = 0; y < r2 - r1; y++) {
			float xx = (y + r1) * co + nx / 2 + dx;
			float yy = (y + r1) * si + ny / 2 + dy;
			float t = xx - floor(xx);
			float u = yy - floor(yy);
			int k = (int) floor(xx) + (int) (floor(yy)) * nx;
			dd[x + y * xs] =
				Util::bilinear_interpolate(d[k], d[k + 1], d[k + nx + 1], d[k + nx], t,
										   u) * (y + r1);
		}
	}
	done_data();
	ret->done_data();

	return ret;
}


void EMData::add_incoherent(EMData * obj)
{
	if (!obj) {
		throw NullEMDataObjectError(__FILE__, __LINE__);
	}

	if (!obj->is_complex() || !is_complex()) {
		LOGERR("add incoherent can only work on complex images");
		return;
	}

	if (!EMUtil::is_same_size(this, obj)) {
		throw NotSameImageSizeError(__FILE__, __LINE__,
									"add incoherent can only add same size image");
	}

	ri2ap();
	obj->ri2ap();

	float *dest = get_data();
	float *src = obj->get_data();
	int size = nx * ny * nz;
	for (int j = 0; j < size; j += 2) {
		dest[j] = hypot(src[j], dest[j]);
		dest[j + 1] = 0;
	}

	obj->done_data();
	done_data();
	update();
}


vector < float >EMData::calc_radial_dist(int n, float x0, float dx)
{
	float *yc = new float[n];
	float *d = new float[n];
	for (int i = 0; i < n; i++) {
		yc[i] = 0;
		d[i] = 0;
	}

	int c = 0;
	float half_nz = 0;
	if (nz > 1) {
		half_nz = nz / 2.0f;
	}

	int step = 1;
	if (is_complex()) {
		step = 2;
	}

	for (int z = 0; z < nz; z++) {
		for (int y = 0; y < ny; y++) {
			for (int x = 0; x < nx; x += step, c += step) {
				float r = 0;
				if (is_complex()) {
					r = (Util::hypot3(x / 2.0f, y - ny / 2.0f, z - half_nz) - x0) / dx;
				}
				else {
					r = (Util::hypot3(x - nx / 2.0f, y - ny / 2.0f, z - half_nz) - x0) / dx;
				}

				int i = (int) floor(r);
				r -= i;
				if (i == 0) {
					d[0] += Util::square(rdata[c]) * (1 - r);
					yc[0] += (1 - r);
				}
				else if (i == n - 1) {
					d[n - 1] += Util::square(rdata[c]) * r;
					yc[n - 1] += r;
				}
				else if (i > 0 && i < n - 1) {
					float h = 0;
					if (is_complex()) {
						if (is_ri()) {
							h = Util::square(rdata[c]) + Util::square(rdata[c + 1]);
						}
						else {
							h = rdata[c] * rdata[c];
						}
					}
					else {
						h = rdata[c];
					}

					d[i] += h * (1 - r);
					yc[i] += (1 - r);
					d[i + 1] += h * r;
					yc[i + 1] += r;
				}
			}
		}
	}

	for (int i = 0; i < n; i++) {
		if (yc[i] != 0) {
			d[i] /= yc[i];
		}
	}

	for (int i = 1; i < n - 1; i++) {
		if (yc[i] < 0.1) {
			d[i] = (d[i - 1] + d[i + 1]) / 2.0f;
		}
	}

	delete[]yc;
	yc = 0;

	vector < float >dv(n);
	for (int i = 0; i < n; i++) {
		dv[i] = d[i];
	}

	delete[]d;
	d = 0;

	return dv;
}

vector < float >EMData::calc_radial_dist(int n, float x0, float dx, float acen, float awid)
{
	if (nz > 1) {
		LOGERR("calc radial dist: can only handle 2D images.");
		return vector < float >();
	}

	float *yc = new float[n];
	float *yc2 = new float[n];
	float *d = new float[n];
	float *d2 = new float[n];

	for (int i = 0; i < n; i++) {
		yc[i] = 0;
		yc2[i] = 0;
		d[i] = 0;
		d2[i] = 0;
	}

	int step = 1;
	if (is_complex()) {
		step = 2;
	}

	int c = 0;
	for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nx; x += step, c += step) {
			float a = 0;
			if (y != ny / 2.0f || x != 0) {
				if (is_complex()) {
					a = atan2(y - ny / 2.0f, x / 2.0f);
				}
				else {
					a = atan2(y - ny / 2.0f, x - nx / 2.0f);
				}
			}

			if (fabs(Util::angle_sub_pi(a, acen)) <= awid) {
				float r = 0;
				if (is_complex()) {
					r = (hypot(x / 2.0f, y - ny / 2.0f) - x0) / dx;
				}
				else {
					r = (hypot(x - (nx - 1) / 2.0f, y - (ny - 1) / 2.0f) - x0) / dx;
				}

				int i = (int) floor(r);
				r -= i;

				if (i == 0) {
					d[i] += Util::square(rdata[c]) * (1 - r);
					yc[i] += 1 - r;
				}
				else if (i == n - 1) {
					d[i] += Util::square(rdata[c]) * r;
					yc[i] += r;
				}
				else if (i > 0 && i < n - 1) {
					float h = 0;
					if (is_complex()) {
						if (is_ri()) {
							h = Util::square(rdata[c]) + Util::square(rdata[c + 1]);
						}
						else {
							h = rdata[c] * rdata[c];
						}
					}
					else {
						h = rdata[c];
					}

					d[i] += h * (1 - r);
					yc[i] += 1 - r;
					d[i + 1] += h * r;
					yc[i + 1] += r;
				}
			}
		}
	}

	const float scale1 = 0.57f;
	const float scale2 = 0.105f;

	for (int i = 0; i < n; i++) {
		d2[i] = d[i];
		yc2[i] = yc[i];
	}

	for (int i = 1; i < n; i++) {
		d2[i - 1] += d[i] * scale1;
		yc2[i - 1] += yc[i] * scale1;
	}

	for (int i = 2; i < n; i++) {
		d2[i - 2] += d[i] * scale2;
		yc2[i - 2] += yc[i] * scale2;
	}

	for (int i = 0; i < n - 1; i++) {
		d2[i + 1] += d[i] * scale1;
		yc2[i + 1] += yc[i] * scale1;
	}

	for (int i = 0; i < n - 2; i++) {
		d2[i + 2] += d[i] * scale2;
		yc2[i + 2] += yc[i] * scale2;
	}

	for (int i = 0; i < n; i++) {
		if (yc2[i] != 0) {
			d[i] = d2[i] / yc2[i];
		}
		else {
			d[i] = 0;
		}
	}

	vector < float >dv(n);
	for (int i = 0; i < n; i++) {
		dv[i] = d[i];
	}

	delete[]yc;
	yc = 0;

	delete[]yc2;
	yc2 = 0;

	delete[]d2;
	d2 = 0;

	delete[]d;
	d = 0;

	return dv;
}


float EMData::calc_dist(EMData * second_img, int y_index) const
{
	if (second_img->get_xsize() != nx || ny != 1) {
		return -1;
	}

	if (y_index > second_img->get_ysize() || y_index < 0) {
		return -1;
	}

	float ret = 0;
	float *d1 = get_data();
	float *d2 = second_img->get_data() + second_img->get_xsize() * y_index;

	for (int i = 0; i < nx; i++) {
		ret += Util::square(d1[i] - d2[i]);
	}
	return sqrt(ret);
}



EMData *EMData::calc_flcf(EMData * with, int radius, string mask_filter)
{
	if (!with) {
		LOGERR("input image is NULL");
		return 0;
	}

	Dict filter_dict;
	if (mask_filter == "SharpMask") {
		filter_dict["value"] = 0;
	}

	EMData *img1 = this->copy(false);
	EMData *img2 = with->copy(false);

	int img1_nx = img1->get_xsize();
	int img1_ny = img1->get_ysize();
	int img1_nz = img1->get_zsize();
	int img1_size = img1_nx * img1_ny * img1_nz;

	float img1min = img1->get_attr("minimum");
	img1->add(-img1min);

	float img2min = img2->get_attr("minimum");
	img2->add(-img2min);

	filter_dict["outer_radius"] = radius;

	EMData *img1_copy = img1->copy(false);
	img1_copy->to_one();
	img1_copy->filter(mask_filter, filter_dict);
	img1_copy->filter("Phase180");

	int num = 0;
	float *img1_copy_data = img1_copy->get_data();

	for (int i = 0; i < img1_size; i++) {
		if (img1_copy_data[i] == 1) {
			num++;
		}
	}

	img2->filter(mask_filter, filter_dict);

	float *img2_data = img2->get_data();
	double lsum = 0;
	double sumsq = 0;

	for (int i = 0; i < img1_size; i++) {
		lsum += img2_data[i];
		sumsq += img2_data[i] * img2_data[i];
	}

	float sq = ((num * sumsq - lsum * lsum) / (num * num));
	if (sq < 0) {
		LOGERR("sigma < 0");
		return 0;
	}

	float mean = lsum / num;
	float sigma = sqrt(sq);
	float th = 0.00001f;

	if (sq > th) {
		for (int i = 0; i < img1_size; i++) {
			img2_data[i] = (img2_data[i] - mean) / sigma;
		}
	}
	else {
		for (int i = 0; i < img1_size; i++) {
			img2_data[i] -= mean;
		}
	}

	img2->done_data();

	EMData *img2_copy = img2->copy(false);
	delete img2;
	img2 = 0;

	img2_copy->filter(mask_filter, filter_dict);
	img2_copy->filter("Phase180");

	delete img1_copy;
	img1_copy = 0;

	EMData *img1_copy2 = img1->copy(false);

	img1_copy2->filter("Square");

	EMData *ccf = img1->calc_ccf(img2_copy);
	delete img2_copy;
	img2_copy = 0;

	ccf->mult(img1_size);

	EMData *conv1 = img1->convolute(img1_copy2);
	delete img1;
	img1 = 0;

	conv1->mult(img1_size);
	conv1->mult(1.0f / num);

	EMData *conv2 = img1_copy2->convolute(img1_copy2);
	delete img1_copy2;
	img1_copy2 = 0;

	conv2->mult(img1_size);
	conv1->filter("Square");
	conv1->mult(1.0f / (num * num));

	EMData *conv2_copy = conv2->copy(false);
	delete conv2;
	conv2 = 0;

	conv2_copy->sub(*conv1);
	delete conv1;
	conv1 = 0;

	conv2_copy->mult(1.0f / num);
	conv2_copy->filter("Sqrt");

	EMData *ccf_copy = ccf->copy(false);
	delete ccf;
	ccf = 0;

	ccf_copy->mult(1.0f / num);

	float *lcfd = ccf_copy->get_data();
	float *vdd = conv2_copy->get_data();

	for (int i = 0; i < img1_size; i++) {
		if (vdd[i] > 0) {
			lcfd[i] /= vdd[i];
		}
	}
	delete conv2_copy;
	conv2_copy = 0;

	ccf_copy->done_data();
	EMData *lcf = ccf_copy->copy(false);
	delete ccf_copy;
	ccf_copy = 0;

	return lcf;
}

EMData *EMData::convolute(EMData * with)
{
	EMData *f1 = do_fft();
	if (!f1) {
		return 0;
	}

	f1->ap2ri();

	EMData *cf = 0;
	if (with) {
		EMData *cf1 = with->do_fft();
		if (!cf1) {
			return 0;
		}
		cf = cf1->copy(false);
		cf->ap2ri();
	}
	else {
		cf = f1->copy(false);
	}

	if (with && !EMUtil::is_same_size(f1, cf)) {
		LOGERR("Convolution: Images must be the same size.");
		return 0;
	}

	float *rdata1 = f1->get_data();
	float *rdata2 = cf->get_data();
	int cf_size = cf->get_xsize() * cf->get_ysize() * cf->get_zsize();

	for (int i = 0; i < cf_size; i += 2) {
		rdata2[i] = rdata1[i] * rdata2[i] - rdata1[i + 1] * rdata2[i + 1];
		rdata2[i + 1] = rdata1[i + 1] * rdata2[i] + rdata1[i] * rdata2[i + 1];
	}

	cf->done_data();
	EMData *f2 = cf->do_ift();

	delete cf;
	cf = 0;

	return f2;

}


void EMData::common_lines(EMData * image1, EMData * image2,
						  int mode, int steps, bool horizontal)
{
	if (!image1 || !image2) {
		throw NullEMDataObjectError(__FILE__, __LINE__, 
									"common lines: cannot handle NULL image");
	}

	if (mode < 0 || mode > 2) {
		throw InvalidModeError(0, 2, mode, __FILE__, __LINE__, "invalid mode");
	}

	if (!image1->is_complex()) {
		image1 = image1->do_fft();
	}
	if (!image2->is_complex()) {
		image2 = image2->do_fft();
	}

	image1->ap2ri();
	image2->ap2ri();

	if (!EMUtil::is_same_size(image1, image2)) {
		throw NotSameImageSizeError(__FILE__, __LINE__, 
									"common lines: input images must be the same size");
	}

	int image2_nx = image2->get_xsize();
	int image2_ny = image2->get_ysize();

	int rmax = image2_ny / 4 - 1;
	int array_size = steps * rmax * 2;
	float *im1 = new float[array_size];
	float *im2 = new float[array_size];
	for (int i = 0; i < array_size; i++) {
		im1[i] = 0;
		im2[i] = 0;
	}

	set_size(steps * 2, steps * 2, 1);

	float *image1_data = image1->get_data();
	float *image2_data = image2->get_data();

	float da = M_PI / steps;
	float a = -M_PI / 2.0f + da / 2.0f;
	int jmax = 0;

	for (int i = 0; i < steps * 2; i += 2, a += da) {
		float s1 = 0;
		float s2 = 0;
		int i2 = i * rmax;
		int j = 0;

		for (float r = 3.0f; r < rmax - 3.0f; j += 2, r += 1.0f) {
			float x = r * cos(a);
			float y = r * sin(a);

			if (x < 0) {
				x = -x;
				y = -y;
				LOGERR("CCL ERROR %d, %f !\n", i, -x);
			}

			int k = (int) (floor(x) * 2 + floor(y + image2_ny / 2) * image2_nx);
			int l = i2 + j;
			float x2 = x - floor(x);
			float y2 = y - floor(y);

			im1[l] = Util::bilinear_interpolate(image1_data[k],
												image1_data[k + 2],
												image1_data[k + 2 + image2_nx],
												image1_data[k + image2_nx], x2, y2);

			im2[l] = Util::bilinear_interpolate(image2_data[k],
												image2_data[k + 2],
												image2_data[k + 2 + image2_nx],
												image2_data[k + image2_nx], x2, y2);

			k++;

			im1[l + 1] = Util::bilinear_interpolate(image1_data[k],
													image1_data[k + 2],
													image1_data[k + 2 + image2_nx],
													image1_data[k + image2_nx], x2, y2);

			im2[l + 1] = Util::bilinear_interpolate(image2_data[k],
													image2_data[k + 2],
													image2_data[k + 2 + image2_nx],
													image2_data[k + image2_nx], x2, y2);

			s1 += Util::square_sum(im1[l], im1[l + 1]);
			s2 += Util::square_sum(im2[l], im2[l + 1]);
		}

		jmax = j - 1;
		float sqrt_s1 = sqrt(s1);
		float sqrt_s2 = sqrt(s2);

		int l = 0;
		for (float r = 1; r < rmax; r += 1.0f) {
			int i3 = i2 + l;
			im1[i3] /= sqrt_s1;
			im1[i3 + 1] /= sqrt_s1;
			im2[i3] /= sqrt_s2;
			im2[i3 + 1] /= sqrt_s2;
			l += 2;
		}
	}

	switch (mode) {
	case 0:
		for (int m1 = 0; m1 < 2; m1++) {
			for (int m2 = 0; m2 < 2; m2++) {

				if (m1 == 0 && m2 == 0) {
					for (int i = 0; i < steps; i++) {
						int i2 = i * rmax * 2;
						for (int j = 0; j < steps; j++) {
							int l = i + j * steps * 2;
							int j2 = j * rmax * 2;
							rdata[l] = 0;
							for (int k = 0; k < jmax; k++) {
								rdata[l] += im1[i2 + k] * im2[j2 + k];
							}
						}
					}
				}
				else {
					int steps2 = steps * m2 + steps * steps * 2 * m1;

					for (int i = 0; i < steps; i++) {
						int i2 = i * rmax * 2;
						for (int j = 0; j < steps; j++) {
							int j2 = j * rmax * 2;
							int l = i + j * steps * 2 + steps2;
							rdata[l] = 0;

							for (int k = 0; k < jmax; k += 2) {
								i2 += k;
								j2 += k;
								rdata[l] += im1[i2] * im2[j2];
								rdata[l] += -im1[i2 + 1] * im2[j2 + 1];
							}
						}
					}
				}
			}
		}

		break;
	case 1:
		for (int m1 = 0; m1 < 2; m1++) {
			for (int m2 = 0; m2 < 2; m2++) {
				int steps2 = steps * m2 + steps * steps * 2 * m1;
				int p1_sign = 1;
				if (m1 != m2) {
					p1_sign = -1;
				}

				for (int i = 0; i < steps; i++) {
					int i2 = i * rmax * 2;

					for (int j = 0; j < steps; j++) {
						int j2 = j * rmax * 2;

						int l = i + j * steps * 2 + steps2;
						rdata[l] = 0;
						float a = 0;

						for (int k = 0; k < jmax; k += 2) {
							i2 += k;
							j2 += k;

							float a1 = hypot(im1[i2], im1[i2 + 1]);
							float p1 = atan2(im1[i2 + 1], im1[i2]);
							float p2 = atan2(im2[j2 + 1], im2[j2]);

							rdata[l] += Util::angle_sub_2pi(p1_sign * p1, p2) * a1;
							a += a1;
						}

						rdata[l] /= (float)(a * M_PI / 180.0f);
					}
				}
			}
		}

		break;
	case 2:
		for (int m1 = 0; m1 < 2; m1++) {
			for (int m2 = 0; m2 < 2; m2++) {
				int steps2 = steps * m2 + steps * steps * 2 * m1;

				for (int i = 0; i < steps; i++) {
					int i2 = i * rmax * 2;

					for (int j = 0; j < steps; j++) {
						int j2 = j * rmax * 2;
						int l = i + j * steps * 2 + steps2;
						rdata[l] = 0;

						for (int k = 0; k < jmax; k += 2) {
							i2 += k;
							j2 += k;
							rdata[l] += hypot(im1[i2], im1[i2 + 1]) * hypot(im2[j2], im2[j2 + 1]);
						}
					}
				}
			}
		}

		break;
	default:
		break;
	}

	if (horizontal) {
		float *tmp_array = new float[ny];
		for (int i = 1; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				tmp_array[j] = get_value_at(i, j);
			}
			for (int j = 0; j < ny; j++) {
				set_value_at(i, j, 0, tmp_array[(j + i) % ny]);
			}
		}
		delete[]tmp_array;
		tmp_array = 0;
	}

	delete[]im1;
	im1 = 0;

	delete im2;
	im2 = 0;


	image1->done_data();
	image2->done_data();
	done_data();
	update();
}



void EMData::common_lines_real(EMData * image1, EMData * image2,
							   int steps, bool horiz)
{
	if (!image1 || !image2) {
		throw NullEMDataObjectError(__FILE__, __LINE__, 
									"common lines real: cannot handle NULL image");
	}

	if (!EMUtil::is_same_size(image1, image2)) {
		throw NotSameImageSizeError(__FILE__, __LINE__, 
									"common lines real: input images must be the same size");
	}

	int steps2 = steps * 2;
	int image_ny = image1->get_ysize();
	EMData *image1_copy = image1->copy();
	EMData *image2_copy = image2->copy();

	float *im1 = new float[steps2 * image_ny];
	float *im2 = new float[steps2 * image_ny];

	EMData *images[] = { image1_copy, image2_copy };
	float *ims[] = { im1, im2 };

	for (int m = 0; m < 2; m++) {
		float *im = ims[m];
		float a = M_PI / steps2;

		for (int i = 0; i < steps2; i++) {
			images[i]->rotate(-a, 0, 0);
			float *data = images[i]->get_data();

			for (int j = 0; j < image_ny; j++) {
				float sum = 0;
				for (int k = 0; k < image_ny; k++) {
					sum += data[j * image_ny + k];
				}
				im[i * image_ny + j] = sum;
			}

			float sum1 = 0;
			float sum2 = 0;
			for (int j = 0; j < image_ny; j++) {
				int l = i * image_ny + j;
				sum1 += im[l];
				sum2 += im[l] * im[l];
			}

			float mean = sum1 / image_ny;
			float sigma = sqrt(sum2 / image_ny - sum1 * sum1);

			for (int j = 0; j < image_ny; j++) {
				int l = i * image_ny + j;
				im[l] = (im[l] - mean) / sigma;
			}

			images[i]->done_data();
			a += M_PI / steps;
		}
	}

	set_size(steps2, steps2, 1);
	float *data1 = get_data();

	if (horiz) {
		for (int i = 0; i < steps2; i++) {
			for (int j = 0, l = i; j < steps2; j++, l++) {
				if (l == steps2) {
					l = 0;
				}

				float sum = 0;
				for (int k = 0; k < image_ny; k++) {
					sum += im1[i * image_ny + k] * im2[l * image_ny + k];
				}
				data1[i + j * steps2] = sum;
			}
		}
	}
	else {
		for (int i = 0; i < steps2; i++) {
			for (int j = 0; j < steps2; j++) {
				float sum = 0;
				for (int k = 0; k < image_ny; k++) {
					sum += im1[i * image_ny + k] * im2[j * image_ny + k];
				}
				data1[i + j * steps2] = sum;
			}
		}
	}

	done_data();

	delete image1_copy;
	image1_copy = 0;

	delete image2_copy;
	image2_copy = 0;

	delete[]im1;
	im1 = 0;

	delete[]im2;
	im2 = 0;
}


void EMData::cut_slice(EMData * map, float dz, Rotation * ort,
					   bool interpolate, float dx, float dy)
{
	if (!map) {
		throw NullEMDataObjectError(__FILE__, __LINE__, "input map is empty");
	}

	Rotation r(0, 0, 0, Rotation::EMAN);
	if (!ort) {
		ort = &r;
	}

	float *sdata = map->get_data();
	float *ddata = get_data();

	Matrix3f mx = ort->get_matrix3();

	int map_nx = map->get_xsize();
	int map_ny = map->get_ysize();
	int map_nz = map->get_zsize();
	int map_nxy = map_nx * map_ny;

	float mdz0 = dz * mx[0][2] + map_nx / 2;
	float mdz1 = dz * mx[1][2] + map_ny / 2;
	float mdz2 = dz * mx[2][2] + map_nz / 2;

	for (int y = 0; y < ny; y++) {
		int y2 = (int) (y - ny / 2 - dy);
		float my2_0 = y2 * mx[0][1] + mdz0;
		float my2_1 = y2 * mx[1][1] + mdz1;
		float my2_2 = y2 * mx[2][1] + mdz2;

		for (int x = 0; x < nx; x++) {
			int x2 = (int) (x - nx / 2 - dx);
			float xx = x2 * mx[0][0] + my2_0;
			float yy = x2 * mx[1][0] + my2_1;
			float zz = x2 * mx[2][0] + my2_2;
			int l = x + y * nx;

			if (xx < 0 || yy < 0 || zz < 0 || xx > map_nx - 2 || yy > map_ny - 2 || zz > map_nz - 2) {
				ddata[l] = 0;
			}
			else {
				float t = xx - floor(xx);
				float u = yy - floor(yy);
				float v = zz - floor(zz);

				if (interpolate) {
					int k = (int) (floor(xx) + floor(yy) * map_nx + floor(zz) * map_nxy);

					ddata[l] = Util::trilinear_interpolate(sdata[k],
														   sdata[k + 1],
														   sdata[k + map_nx],
														   sdata[k + map_nx + 1],
														   sdata[k + map_nxy],
														   sdata[k + map_nxy + 1],
														   sdata[k + map_nx + map_nxy],
														   sdata[k + map_nx + map_nxy + 1],
														   t, u, v);
				}
				else {
					int k = Util::round(xx) + Util::round(yy) * map_nx + Util::round(zz) * map_nxy;
					ddata[l] = sdata[k];
				}
			}
		}
	}

	done_data();
	map->done_data();
}


void EMData::uncut_slice(EMData * map, float dz, Rotation * ort, float dx, float dy)
{
	if (!map) {
		throw NullEMDataObjectError(__FILE__, __LINE__,
									"input map is empty");
	}

	Rotation r(0, 0, 0, Rotation::EMAN);
	if (!ort) {
		ort = &r;
	}

	float *ddata = map->get_data();
	float *sdata = get_data();

	Matrix3f mx = ort->get_matrix3();

	int map_nx = map->get_xsize();
	int map_ny = map->get_ysize();
	int map_nz = map->get_zsize();
	int map_nxy = map_nx * map_ny;

	float mdz0 = dz * mx[0][2] + map_nx / 2;
	float mdz1 = dz * mx[1][2] + map_ny / 2;
	float mdz2 = dz * mx[2][2] + map_nz / 2;

	for (int y = 0; y < ny; y++) {
		int y2 = (int) (y - ny / 2 - dy);

		float my2_0 = y2 * mx[0][1] + mdz0;
		float my2_1 = y2 * mx[1][1] + mdz1;
		float my2_2 = y2 * mx[2][1] + mdz2;

		for (int x = 0; x < nx; x++) {
			int x2 = (int) (x - nx / 2 - dx);

			float xx = x2 * mx[0][0] + my2_0;
			float yy = x2 * mx[1][0] + my2_1;
			float zz = x2 * mx[2][0] + my2_2;

			if (xx >= 0 && yy >= 0 && zz >= 0 && xx <= map_nx - 2 && yy <= map_ny - 2
				&& zz <= map_nz - 2) {
				int k = Util::round(xx) + Util::round(yy) * map_nx + Util::round(zz) * map_nxy;
				ddata[k] = sdata[x + y * nx];
			}
		}
	}

	done_data();
	map->done_data();
}

vector<Pixel> EMData::calc_highest_locations(float threshold)
{
	vector<Pixel> result;
	
	int di = 1;
	if (is_complex() && !is_ri()) {
		di = 2;
	}

	int nxy = nx * ny;

	for (int j = 0; j < nz; j++) {
		int cur_z = j * nxy;

		for (int k = 0; k < ny; k++) {
			int cur_y = k * nx + cur_z;

			for (int l = 0; l < nx; l += di) {
				float v = rdata[l + cur_y];
				if (v > threshold) {
					result.push_back(Pixel(l, k, j, v));
				}
			}
		}
	}

	std::sort(result.begin(), result.end());
	
	return result;
}


float EMData::get_edge_mean() const
{
	int di = 0;
	double edge_sum = 0;
	float edge_mean = 0;
	int nxy = nx * ny;

	if (nz == 1) {
		for (int i = 0, j = (ny - 1) * nx; i < nx; i++, j++) {
			edge_sum += rdata[i] + rdata[j];
		}
		for (int i = 0, j = nx - 1; i < nxy; i += nx, j += nx) {
			edge_sum += rdata[i] + rdata[j];
		}
		edge_mean = edge_sum / (nx * 2 + ny * 2);
	}
	else {
		if (nx == ny && nx == nz * 2 - 1) {
			for (int j = (nxy * (nz - 1)); j < nxy * nz; j++, di++) {
				edge_sum += rdata[j];
			}
		}
		else {
			for (int i = 0, j = (nxy * (nz - 1)); i < nxy; i++, j++, di++) {
				edge_sum += rdata[i] + rdata[j];
			}
		}

		int nxy2 = nx * (ny - 1);
		for (int k = 1; k < nz - 1; k++) {
			int k2 = k * nxy;
			int k3 = k2 + nxy2;
			for (int i = 0; i < nx; i++, di++) {
				edge_sum += rdata[i + k2] + rdata[i + k3];
			}
		}
		for (int k = 1; k < nz - 1; k++) {
			int k2 = k * nxy;
			int k3 = nx - 1 + k2;
			for (int i = 1; i < ny - 1; i++, di++) {
				edge_sum += rdata[i * nx + k2] + rdata[i * nx + k3];
			}
		}

		edge_mean = edge_sum / (di * 2);
	}

	return edge_mean;
}

float EMData::get_circle_mean()
{
	static bool busy = false;
	static EMData *mask = 0;

	while (busy);
	busy = true;

	if (!mask || !EMUtil::is_same_size(this, mask)) {
		if (!mask) {
			mask = new EMData();
		}
		mask->set_size(nx, ny, nz);

		float radius = ny / 2 - 2;
		mask->filter("SharpMask", Dict("inner_radius", radius - 1, "outer_radius", radius + 1));

		int n = 0;
		float *d = mask->get_data();

		for (int i = 0; i < nx * ny * nz; i++) {
			if (d[i]) {
				n++;
			}
		}
		mask->done_data();
		mask->mult(1.0f / n);
	}

	float result = dot(mask);
	busy = false;

	return result;
}



float EMData::dot(EMData * with, bool evenonly)
{
	DotCmp dot_cmp;
	Dict cmp_params;
	cmp_params["with"] = with;
	cmp_params["evenonly"] = (int) evenonly;
	float r = dot_cmp.cmp(this);
	return r;
}


void EMData::setup_insert_slice(int size)
{
	const float scale = 1.0e-10;
	set_size(size + 2, size, size);
	set_complex(true);
	set_ri(true);
	to_zero();

	get_data();
	for (int i = 0; i < nx * ny * nz; i += 2) {
		float f = Util::get_frand(0, 2.0f * M_PI);
		rdata[i] = scale * sin(f);
		rdata[i + 1] = scale * cos(f);
	}

	done_data();

	if (!parent) {
		parent = new EMData();
	}
	parent->set_size(size + 2, size, size);
}


EMData *EMData::get_row(int row_index) const
{
	EMData *ret = new EMData();
	ret->set_size(nx, 1, 1);
	memcpy(ret->get_data(), get_data() + nx * row_index, nx * sizeof(float));
	ret->done_data();
	return ret;
}


void EMData::set_row(const EMData * d, int row_index)
{
	float *dst = get_data();
	float *src = d->get_data();
	memcpy(dst + nx * row_index, src, nx * sizeof(float));
	done_data();
}

EMData *EMData::get_col(int col_index) const
{
	EMData *ret = new EMData();
	ret->set_size(ny, 1, 1);
	float *dst = ret->get_data();
	float *src = get_data();

	for (int i = 0; i < ny; i++) {
		dst[i] = src[i * nx + col_index];
	}

	ret->done_data();
	return ret;
}


void EMData::set_col(const EMData * d, int n)
{
	float *dst = get_data();
	float *src = d->get_data();

	for (int i = 0; i < ny; i++) {
		dst[i * nx + n] = src[i];
	}

	done_data();
}

float EMData::calc_density_center()
{
	float center = get_attr("mean");
	float sigma = get_attr("sigma");
	float ds = sigma / 2;
	size_t size = nx * ny * nz;
	float *d = get_data();
	float sigma1 = sigma / 20;
	float sigma2 = sigma / 1000;

	while (ds > sigma1) {
		double sum = 0;
		int norm = 0;

		for (size_t i = 0; i < size; i++) {
			if (fabs(d[i] - center) < ds) {
				sum += d[i];
				norm++;
			}
		}
		if (!norm) {
			break;
		}
		double mean = sum / norm;
		if (fabs(mean - center) < sigma2) {
			ds *= 0.5f;
		}
		center = mean;
	}

	return center;
}

float EMData::calc_sigma_diff()
{
	float *d = get_data();
	float mean = get_attr("mean");
	float sigma = get_attr("sigma");
	
	double sum_up = 0;
	double sum_down = 0;
	int nup = 0;
	int ndown = 0;

	size_t size = nx * ny * nz;

	for (size_t i = 0; i < size; i++) {
		if (d[i] > mean) {
			sum_up += Util::square(d[i] - mean);
			nup++;
		}
		else {
			sum_down += Util::square(mean - d[i]);
			ndown++;
		}
	}

	float sigup = sqrt(sum_up / nup);
	float sigdown = sqrt(sum_down / ndown);
	float sig_diff = fabs(sigup - sigdown) / sigma;

	return sig_diff;

}

int EMData::get_ndim() const
{
	if (nz == 1) {
		if (ny == 1) {
			return 1;
		}
		else {
			return 2;
		}
	}

	return 3;
}


void EMData::apply_radial_func(float x0, float step, vector < float >array, bool interp)
{
	if (is_complex()) {
		return;
	}

	int n = static_cast < int >(array.size());

	ap2ri();
	get_data();

	size_t ndims = get_ndim();

	if (ndims == 2) {
		int k = 0;
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i += 2, k += 2) {
				float r = hypot(i / 2.0f, (j - ny / 2.0f));
				r = (r - x0) / step;

				int l = 0;
				if (interp) {
					l = (int) floor(r);
				}
				else {
					l = (int) floor(r + 1);
				}

				r -= l;

				float f = 0;
				if (l >= n - 2) {
					f = array[n - 1];
				}
				else {
					if (interp) {
						f = (array[l] * (1.0f - r) + array[l + 1] * r);
					}
					else {
						f = array[l];
					}
				}

				rdata[k] *= f;
				rdata[k + 1] *= f;
			}
		}
	}
	else if (ndims == 3) {
		int k = 0;
		for (int m = 0; m < nz; m++) {
			float mnz = (m - nz / 2.0f) * (m - nz / 2.0f);
			for (int j = 0; j < ny; j++) {
				float jny = (j - ny / 2.0f) * (j - ny / 2.0f);
				for (int i = 0; i < nx; i += 2, k += 2) {
					float r = sqrt((i * i / 4.0f) + jny + mnz);
					r = (r - x0) / step;

					int l = 0;
					if (interp) {
						l = (int) floor(r);
					}
					else {
						l = (int) floor(r + 1);
					}

					r -= l;

					float f = 0;
					if (l >= n - 2) {
						f = array[n - 1];
					}
					else {
						if (interp) {
							f = (array[l] * (1.0f - r) + array[l + 1] * r);
						}
						else {
							f = array[l];
						}
					}

					rdata[k] *= f;
					rdata[k + 1] *= f;
				}
			}
		}

	}

	done_data();
	update();
}
