/**
 * $Id$
 */

#include "emdata.h"
#include "log.h"
#include "all_imageio.h"
#include "ctf.h"
#include "processor.h"
#include "aligner.h"
#include "cmp.h"
#include "emfft.h"
#include "projector.h"
#include "byteorder.h"
#include "emconstants.h"

#include "gsl_sf_bessel.h"
#include "gsl_errno.h"
#include <complex>
#include <algorithm>
#include <boost/array.hpp>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#ifdef WIN32
	#define M_PI 3.14159265358979323846f
#endif	//WIN32

using namespace EMAN;
using namespace std;
using namespace boost;

int EMData::totalalloc=0;		// mainly used for debugging/memory leak purposes

EMData::EMData()
{
	ENTERFUNC;

	rdata = 0;
	supp = 0;
	ctf = 0;

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
	xoff = yoff = zoff = 0;

	EMData::totalalloc++;
#ifdef MEMDEBUG
	printf("EMDATA+  %4d    %p\n",EMData::totalalloc,this);
#endif
	EXITFUNC;
}

EMData::~EMData()
{
	ENTERFUNC;
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

	EMData::totalalloc--;
#ifdef MEMDEBUG
	printf("EMDATA-  %4d    %p\n",EMData::totalalloc,this);
#endif
	EXITFUNC;
}

void EMData::read_image(const string & filename, int img_index, bool nodata,
						const Region * region, bool is_3d)
{
	ENTERFUNC;

	ImageIO *imageio = EMUtil::get_imageio(filename, ImageIO::READ_ONLY);

	if (!imageio) {
		throw ImageFormatException("cannot create an image io");
	}
	else {
		int err = imageio->read_header(attr_dict, img_index, region, is_3d);
		if (err) {
			throw ImageReadException(filename, "imageio read header failed");
		}
		else {
			if (imageio->is_complex_mode()) {
				set_complex(true);
				set_fftpad(true);
			}
			if (attr_dict.has_key("is_fftodd")) {
				if (1 == int(attr_dict["is_fftodd"]))
					set_fftodd(true);
			}
			if ((int) attr_dict["is_ri"] == 1) {
				set_ri(true);
			}

			save_byteorder_to_dict(imageio);

			nx = attr_dict["nx"];
			ny = attr_dict["ny"];
			nz = attr_dict["nz"];

			if (!ctf) {
				ctf = new SimpleCtf();
			}
			err = imageio->read_ctf(*ctf, img_index);
			if (err) {
				if( ctf )
				{
					delete ctf;
					ctf = 0;
				}
				flags &= ~EMDATA_HASCTFF;
			}
			else {
				flags |= EMDATA_HASCTFF;
			}

			if (!nodata) {
				set_size(nx, ny, nz);
				int err = imageio->read_data(rdata, img_index, region, is_3d);
				if (err) {
					throw ImageReadException(filename, "imageio read data failed");
				}
				else {
					flags |= EMDATA_NEEDUPD;
				}
			}
		}
	}

#ifndef IMAGEIO_CACHE
	if( imageio )
	{
		delete imageio;
		imageio = 0;
	}
#endif

	EXITFUNC;
}


void EMData::write_image(const string & filename, int img_index,
						 EMUtil::ImageType imgtype,
						 bool header_only, const Region * region,
						 EMUtil::EMDataType filestoragetype,
						 bool use_host_endian)
{
	ENTERFUNC;

	if (is_complex() && is_shuffled())
		fft_shuffle();

	if (imgtype == EMUtil::IMAGE_UNKNOWN) {
		char *ext = strrchr(filename.c_str(), '.');
		if (ext) {
			ext++;
			imgtype = EMUtil::get_image_ext_type(ext);
		}
	}

	ImageIO::IOMode rwmode = ImageIO::READ_WRITE;

	if (Util::is_file_exist(filename)) {
		LOGVAR("file exists");
		if (!header_only && region == 0) {
			ImageIO * tmp_imageio = EMUtil::get_imageio(filename, ImageIO::READ_ONLY,
														imgtype);
			if (tmp_imageio->is_single_image_format()) {
				rwmode = ImageIO::WRITE_ONLY;
			}
#ifndef IMAGEIO_CACHE
			if( tmp_imageio )
			{
				delete tmp_imageio;
				tmp_imageio = 0;
			}
#endif
		}
	}

	LOGVAR("getimageio %d",rwmode);
	ImageIO *imageio = EMUtil::get_imageio(filename, rwmode, imgtype);
	if (!imageio) {
		throw ImageFormatException("cannot create an image io");
	}
	else {
		update_stat();
		if (img_index < 0) {
			img_index = imageio->get_nimg();
		}
		LOGVAR("header write %d",img_index);
		int err = imageio->write_header(attr_dict, img_index, region, filestoragetype,
										use_host_endian);
		if (err) {
			throw ImageWriteException(filename, "imageio write header failed");
		}
		else {
			if (ctf) {
				imageio->write_ctf(*ctf, img_index);
			}

			if (!header_only) {
				if (imgtype == EMUtil::IMAGE_LST) {

					const char *reffile = attr_dict["LST.reffile"];
					if (strcmp(reffile, "") == 0) {
						reffile = path.c_str();
					}
					int refn = attr_dict["LST.refn"];
					if (refn < 0) {
						refn = pathnum;
					}

					const char *comment = attr_dict["LST.comment"];
					char *lstdata = new char[1024];
					sprintf(lstdata, "%d\t%s", refn, reffile);
					if (comment != "") {
						sprintf(lstdata+strlen(lstdata), "\t%s\n", comment);
					}
					else {
						strcat(lstdata, "\n");
					}
					err = imageio->write_data((float*)lstdata, img_index,
											  region, filestoragetype, use_host_endian);
					if( lstdata )
					{
						delete [] lstdata;
						lstdata = 0;
					}
				}
				else {
					err = imageio->write_data(rdata, img_index, region, filestoragetype,
											  use_host_endian);
				}
				if (err) {
					imageio->flush();
					throw ImageWriteException(filename, "imageio write data failed");
				}
			}
		}
	}
	imageio->flush();

#ifndef IMAGEIO_CACHE
	if( imageio )
	{
		delete imageio;
		imageio = 0;
	}
#endif



	EXITFUNC;
}

void EMData::append_image(const string & filename,
						  EMUtil::ImageType imgtype, bool header_only)
{
	ENTERFUNC;
	write_image(filename, -1, imgtype, header_only, 0);
	EXITFUNC;
}

void EMData::write_lst(const string & filename, const string & reffile,
					   int refn, const string & comment)
{
	ENTERFUNC;
	attr_dict["LST.reffile"] = reffile;
	attr_dict["LST.refn"] = refn;
	attr_dict["LST.comment"] = comment;
	write_image(filename, -1, EMUtil::IMAGE_LST, false);
	EXITFUNC;
}


void EMData::process(const string & processorname, const Dict & params)
{
	ENTERFUNC;
	Processor *f = Factory < Processor >::get(processorname, params);
	if (f) {
		f->process(this);
		if( f )
		{
			delete f;
			f = 0;
		}
	}
	EXITFUNC;
}

float EMData::cmp(const string & cmpname, EMData * with, const Dict & params)
{
	ENTERFUNC;
	float result = 0;
	Cmp *c = Factory < Cmp >::get(cmpname, params);
	if (c) {
		result = c->cmp(this, with);
		if( c )
		{
			delete c;
			c = 0;
		}
	}

	EXITFUNC;
	return result;
}

EMData *EMData::align(const string & aligner_name, EMData * to_img,
					  const Dict & params, const string & cmp_name, const Dict& cmp_params)
{
	ENTERFUNC;
	EMData *result = 0;
	Aligner *a = Factory < Aligner >::get(aligner_name, params);
	if (a) {
		if (cmp_name == "") {
			result = a->align(this, to_img);
		}
		else {
			result = a->align(this, to_img, cmp_name, cmp_params);
		}
		if( a )
		{
			delete a;
			a = 0;
		}
	}

	EXITFUNC;
	return result;
}

EMData *EMData::project(const string & projector_name, const Dict & params)
{
	ENTERFUNC;
	EMData *result = 0;
	Projector *p = Factory < Projector >::get(projector_name, params);
	if (p) {
		result = p->project3d(this);
		if( p )
		{
			delete p;
			p = 0;
		}
	}

	EXITFUNC;
	return result;
}

EMData *EMData::copy() const
{
	ENTERFUNC;
	EMData *ret = new EMData();

	ret->set_size(nx, ny, nz);
	float *data = ret->get_data();
	memcpy(data, rdata, nx * ny * nz * sizeof(float));
	ret->done_data();


	if (ctf) {
		ret->ctf = new SimpleCtf();
		ret->ctf->copy_from(ctf);
	}

	ret->rfp = 0;
	ret->flags = flags & (EMDATA_COMPLEX 
							| EMDATA_RI 
							| EMDATA_PAD 
							| EMDATA_SHUFFLE
							| EMDATA_FLIP
							| EMDATA_FH
							| EMDATA_FFTODD);

	ret->all_translation = all_translation;

	ret->path = path;
	ret->pathnum = pathnum;
	ret->attr_dict = attr_dict;
	ret->update();

	EXITFUNC;
	return ret;
}

EMData *EMData::copy_head() const
{
	ENTERFUNC;
	EMData *ret = new EMData();
	ret->attr_dict = attr_dict;
	ret->set_size(nx, ny, nz);

	if (ctf) {
		ret->ctf = new SimpleCtf();
		ret->ctf->copy_from(ctf);
	}

	ret->rfp = 0;

	ret->flags = flags & (EMDATA_COMPLEX | EMDATA_RI | EMDATA_PAD | EMDATA_FFTODD);

	ret->all_translation = all_translation;

	ret->path = path;
	ret->pathnum = pathnum;

	ret->update();

	EXITFUNC;
	return ret;
}

EMData *EMData::get_rotated_clip(const Transform3D &xform,
								 const IntSize &size, float scale)
{
	EMData *result = new EMData();
	result->set_size(size[0],size[1],size[2]);

	for (int z=-size[2]/2; z<size[2]/2; z++) {
		for (int y=-size[1]/2; y<size[1]/2; y++) {
			for (int x=-size[0]/2; x<size[0]/2; x++) {
				Vec3f xv=Vec3f((float)x,(float)y,(float)z)*xform;
				float v = 0;

				if (xv[0]<0||xv[1]<0||xv[2]<0||xv[0]>nx-2||xv[1]>ny-2||xv[2]>nz-2) v=0.;
				else v=sget_value_at_interp(xv[0],xv[1],xv[2]);
				result->set_value_at(x+size[0]/2,y+size[1]/2,z+size[2]/2,v);
			}
		}
	}
	result->update();

	return result;
}

EMData *EMData::get_clip(const Region & area)
{
	ENTERFUNC;
	if (get_ndim() != area.get_ndim()) {
		LOGERR("cannot get %dD clip out of %dD image", get_ndim(), area.get_ndim());
		return 0;
	}

	EMData *result = new EMData();
	int zsize = (int)area.size[2];
	if (zsize == 0 || nz <= 1) {
		zsize = 1;
	}
	int ysize = (ny<=1 && (int)area.size[1]==0 ? 1 : (int)area.size[1]);

	result->set_size((int)area.size[0], ysize, zsize);

	int x0 = (int) area.origin[0];
	x0 = x0 < 0 ? 0 : x0;

	int y0 = (int) area.origin[1];
	y0 = y0 < 0 ? 0 : y0;

	int z0 = (int) area.origin[2];
	z0 = z0 < 0 ? 0 : z0;

	int x1 = (int) (area.origin[0] + area.size[0]);
	x1 = x1 > nx ? nx : x1;

	int y1 = (int) (area.origin[1] + area.size[1]);
	y1 = y1 > ny ? ny : y1;

	int z1 = (int) (area.origin[2] + area.size[2]);
	z1 = z1 > nz ? nz : z1;
	if (z1 <= 0) {
		z1 = 1;
	}

	int xd0 = (int) (area.origin[0] < 0 ? -area.origin[0] : 0);
	int yd0 = (int) (area.origin[1] < 0 ? -area.origin[1] : 0);
	int zd0 = (int) (area.origin[2] < 0 ? -area.origin[2] : 0);

	size_t clipped_row_size = (x1-x0) * sizeof(float);
	int src_secsize = nx * ny;
	int dst_secsize = (int)(area.size[0] * area.size[1]);

	float *src = rdata + z0 * src_secsize + y0 * nx + x0;
	float *dst = result->get_data();
	dst += zd0 * dst_secsize + yd0 * (int)area.size[0] + xd0;

	int src_gap = src_secsize - (y1-y0) * nx;
	int dst_gap = dst_secsize - (y1-y0) * (int)area.size[0];

	for (int i = z0; i < z1; i++) {
		for (int j = y0; j < y1; j++) {
			memcpy(dst, src, clipped_row_size);
			src += nx;
			dst += (int)area.size[0];
		}
		src += src_gap;
		dst += dst_gap;
	}

	done_data();
	result->done_data();

	if( attr_dict.has_key("apix_x") && attr_dict.has_key("apix_y") &&
		attr_dict.has_key("apix_z") )
	{
		result->attr_dict["apix_x"] = attr_dict["apix_x"];
		result->attr_dict["apix_y"] = attr_dict["apix_y"];
		result->attr_dict["apix_z"] = attr_dict["apix_z"];

		if( attr_dict.has_key("origin_row") && attr_dict.has_key("origin_col") &&
		    attr_dict.has_key("origin_sec") )
		{
			float xorigin = attr_dict["origin_row"];
			float yorigin = attr_dict["origin_col"];
			float zorigin = attr_dict["origin_sec"];

			float apix_x = attr_dict["apix_x"];
			float apix_y = attr_dict["apix_y"];
			float apix_z = attr_dict["apix_z"];

			result->set_xyz_origin(xorigin + apix_x * area.origin[0],
							   	   yorigin + apix_y * area.origin[1],
							       zorigin + apix_z * area.origin[2]);
		}
	}

	result->update();

	result->set_path(path);
	result->set_pathnum(pathnum);

	EXITFUNC;
	return result;
}

void EMData::insert_clip(EMData * block, const IntPoint &origin)
{
	ENTERFUNC;
	int nx1 = block->get_xsize();
	int ny1 = block->get_ysize();
	int nz1 = block->get_zsize();

	Region area(origin[0], origin[1], origin[2],nx1, ny1, nz1);

	if (area.inside_region((float)nx, (float)ny, (float)nz)) {
		throw ImageFormatException("outside of destination image not supported.");
	}

	int x0 = origin[0];
	int y0 = origin[1];
	int y1 = origin[1] + ny1;
	int z0 = origin[2];
	int z1 = origin[2] + nz1;

	size_t inserted_row_size = nx1 * sizeof(float);
	float *inserted_data = block->get_data();
	float *src = inserted_data;
	float *dst = rdata + z0 * nx * ny + y0 * nx + x0;

	for (int i = z0; i < z1; i++) {

		for (int j = y0; j < y1; j++) {
			memcpy(dst, src, inserted_row_size);
			src += nx1;
			dst += nx;
		}
		dst += nx * (ny - ny1);
	}

	flags |= EMDATA_NEEDUPD;
	EXITFUNC;
}

void EMData::insert_scaled_sum(EMData *block, const FloatPoint &center,
							   float scale, float)
{
	ENTERFUNC;
if (get_ndim()==3) {
	// Start by determining the region to operate on
	int xs=(int)floor(block->get_xsize()*scale/2.0);
	int ys=(int)floor(block->get_ysize()*scale/2.0);
	int zs=(int)floor(block->get_zsize()*scale/2.0);
	int x0=(int)center[0]-xs;
	int x1=(int)center[0]+xs;
	int y0=(int)center[1]-ys;
	int y1=(int)center[1]+ys;
	int z0=(int)center[2]-zs;
	int z1=(int)center[2]+zs;

	if (x1<0||y1<0||z1<0||x0>get_xsize()||y0>get_ysize()||z0>get_zsize()) return;	// object is completely outside the target volume

	// make sure we stay inside the volume
	if (x0<0) x0=0;
	if (y0<0) y0=0;
	if (z0<0) z0=0;
	if (x1>get_xsize()) x1=get_xsize();
	if (y1>get_ysize()) y1=get_ysize();
	if (z1>get_zsize()) z1=get_zsize();

	float bx=block->get_xsize()/2.0f;
	float by=block->get_ysize()/2.0f;
	float bz=block->get_zsize()/2.0f;

	for (int x=x0; x<x1; x++) {
		for (int y=y0; y<y1; y++) {
			for (int z=z0; z<z1; z++) {
				rdata[x + y * nx + z * nx * ny] +=
					block->sget_value_at_interp((x-center[0])/scale+bx,(y-center[1])/scale+by,(z-center[2])/scale+bz);
			}
		}
	}
	update();
}
else if (get_ndim()==2) {
	// Start by determining the region to operate on
	int xs=(int)floor(block->get_xsize()*scale/2.0);
	int ys=(int)floor(block->get_ysize()*scale/2.0);
	int x0=(int)center[0]-xs;
	int x1=(int)center[0]+xs;
	int y0=(int)center[1]-ys;
	int y1=(int)center[1]+ys;

	if (x1<0||y1<0||x0>get_xsize()||y0>get_ysize()) return;	// object is completely outside the target volume

	// make sure we stay inside the volume
	if (x0<0) x0=0;
	if (y0<0) y0=0;
	if (x1>get_xsize()) x1=get_xsize();
	if (y1>get_ysize()) y1=get_ysize();

	float bx=block->get_xsize()/2.0f;
	float by=block->get_ysize()/2.0f;

	for (int x=x0; x<x1; x++) {
		for (int y=y0; y<y1; y++) {
			rdata[x + y * nx] +=
				block->sget_value_at_interp((x-center[0])/scale+bx,(y-center[1])/scale+by);
		}
	}
	update();
}
else {
	LOGERR("insert_scaled_sum supports only 2D and 3D data");
	throw ImageDimensionException("2D/3D only");
}

	EXITFUNC;
}

EMData *EMData::get_top_half() const
{
	ENTERFUNC;

	if (get_ndim() != 3) {
		throw ImageDimensionException("3D only");
	}

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

	EXITFUNC;
	return half;
}

void
EMData::onelinenn(int j, int n, int n2, MCArray3D& x,
		          MIArray3D& nr, MCArray2D& bi, const Transform3D& tf) {
	int jp = (j >= 0) ? j+1 : n+j+1;
	// loop over x
	for (int i = 0; i <= n2; i++) {
        if (((i*i+j*j) < n*n/4) && !((0 == i) && (j < 0))) {
			float xnew = i*tf[0][0] + j*tf[1][0];
			float ynew = i*tf[0][1] + j*tf[1][1];
			float znew = i*tf[0][2] + j*tf[1][2];
			complex<float> btq;
			if (xnew < 0.) {
				xnew = -xnew;
				ynew = -ynew;
				znew = -znew;
				btq = conj(bi[i][jp]);
			} else {
				btq = bi[i][jp];
			}
			int ixn = int(xnew + 0.5 + n) - n;
			int iyn = int(ynew + 0.5 + n) - n;
			int izn = int(znew + 0.5 + n) - n;
			if ((ixn <= n2) && (iyn >= -n2) && (iyn <= n2)
				            && (izn >= -n2) && (izn <= n2)) {
				if (ixn >= 0) {
					int iza, iya;
					if (izn >= 0) {
						iza = izn + 1;
					} else {
						iza = n + izn + 1;
					}
					if (iyn >= 0) {
						iya = iyn + 1;
					} else {
						iya = n + iyn + 1;
					}
					x[ixn][iya][iza] += btq;
					nr[ixn][iya][iza]++;
				} else {
					int izt, iyt;
					if (izn > 0) {
						izt = n - izn + 1;
					} else {
						izt = -izn + 1;
					}
					if (iyn > 0) {
						iyt = n - iyn + 1;
					} else {
						iyt = -iyn + 1;
					}
					x[-ixn][iyt][izt] += conj(btq);
					nr[-ixn][iyt][izt]++;
				}
			}

		}
	}
}

void
EMData::nn(MIArray3D& nr, EMData* myfft, const Transform3D& tf) {
	ENTERFUNC;
	int nxc = attr_dict["nxc"]; // # of complex elements along x
	// let's treat nr, bi, and local data as matrices
	MCArray3D x = get_3dcview(0,1,1);
	MCArray2D bi = myfft->get_2dcview(0,1);
	// loop over frequencies in y
	for (int iy = -ny/2 + 1; iy <= ny/2; iy++) {
		onelinenn(iy, ny, nxc, x, nr, bi, tf);
	}
	EXITFUNC;
}

void
EMData::symplane0(MIArray3D& w) {
	ENTERFUNC;
	int nxc = attr_dict["nxc"];
	int n = nxc*2;
	// let's treat the local data as a matrix
	MCArray3D x = get_3dcview(0,1,1);
	for (int iza = 2; iza <= nxc; iza++) {
		for (int iya = 2; iya <= nxc; iya++) {
			x[0][iya][iza] += conj(x[0][n-iya+2][n-iza+2]);
			w[0][iya][iza] += w[0][n-iya+2][n-iza+2];
			x[0][n-iya+2][n-iza+2] = conj(x[0][iya][iza]);
			w[0][n-iya+2][n-iza+2] = w[0][iya][iza];
			x[0][n-iya+2][iza] += conj(x[0][iya][n-iza+2]);
			w[0][n-iya+2][iza] += w[0][iya][n-iza+2];
			x[0][iya][n-iza+2] = conj(x[0][n-iya+2][iza]);
			w[0][iya][n-iza+2] = w[0][n-iya+2][iza];
		}
	}
	for (int iya = 2; iya <= nxc; iya++) {
		x[0][iya][1] += conj(x[0][n-iya+2][1]);
		w[0][iya][1] += w[0][n-iya+2][1];
		x[0][n-iya+2][1] = conj(x[0][iya][1]);
		w[0][n-iya+2][1] = w[0][iya][1];
	}
	for (int iza = 2; iza <= nxc; iza++) {
		x[0][1][iza] += conj(x[0][1][n-iza+2]);
		w[0][1][iza] += w[0][1][n-iza+2];
		x[0][1][n-iza+2] = conj(x[0][1][iza]);
		w[0][1][n-iza+2] = w[0][1][iza];
	}
	EXITFUNC;
}

EMData* EMData::window_center(int l) {
	ENTERFUNC;
	// sanity checks
	int n = nx;
	if (is_complex()) {
		LOGERR("Need real-space data for window_center()");
		throw ImageFormatException(
			"Complex input image; real-space expected.");
	}
	if (is_fftpadded()) {
		// image has been fft-padded, compute the real-space size
		n -= (2 - int(is_fftodd()));
	}
	int corner = (n-l)/2;
	int ndim = get_ndim();
	EMData* ret;
	switch (ndim) {
		case 3:
			if ((n != ny) || (n != nz)) {
				LOGERR("Need the real-space image to be cubic.");
				throw ImageFormatException(
						"Need cubic real-space image.");
			}
			ret = get_clip(Region(corner, corner, corner, l, l, l));
			break;
		case 2:
			if (n != ny) {
				LOGERR("Need the real-space image to be square.");
				throw ImageFormatException(
						"Need square real-space image.");
			}
			ret = get_clip(Region(corner, corner, l, l));
			break;
		case 1:
			ret = get_clip(Region(corner, l));
			break;
		default:
			throw ImageDimensionException(
					"window_center only supports 1-d, 2-d, and 3-d images");
	}
	return ret;
	EXITFUNC;
}

void EMData::postift_depad_corner_inplace() {
	ENTERFUNC;
	int npad = attr_dict["npad"];
	if (0 == npad) npad = 1;
	int offset = is_fftodd() ? 1 : 2;
	int nxold = (nx - offset)/npad;
#ifdef _WIN32
	int nyold = _MAX(ny/npad, 1);
	int nzold = _MAX(nz/npad, 1);
#else
	int nyold = std::max<int>(ny/npad, 1);
	int nzold = std::max<int>(nz/npad, 1);
#endif	//_WIN32
	int bytes = nxold*sizeof(float);
	MArray3D src = get_3dview();
	float* dest = get_data();
	for (int iz=0; iz < nzold; iz++) {
		for (int iy = 0; iy < nyold; iy++) {
			memmove(dest, &src[0][iy][iz], bytes);
			dest += nxold;
		}
	}
	set_size(nxold, nyold, nzold);
	set_fftpad(false);
	done_data();
	update();
	set_complex(false);
	if(ny==1 && nz==1) {
		set_complex_x(false);
	}
	EXITFUNC;
}

void EMData::center_origin()
{
	ENTERFUNC;
	if (is_complex()) {
		LOGERR("Real image expected. Input image is complex.");
		throw ImageFormatException("Real image expected. Input image is complex.");
	}
	MArray3D dat = get_3dview();
	for (int iz = 0; iz < nz; iz++) {
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				// next line multiplies by +/- 1
				dat[ix][iy][iz] *= -2*((ix+iy+iz)%2) + 1;
			}
		}
	}
	done_data();
	update();
	EXITFUNC;
}

void EMData::center_origin_fft()
{
	ENTERFUNC;
	if (!is_complex()) {
		LOGERR("complex image expected. Input image is real image.");
		throw ImageFormatException("complex image expected. Input image is real image.");
	}

	if (!is_ri()) {
		LOGWARN("Only RI should be used. ");
	}
	MCArray3D dat = get_3dcview();
	// iz in [1,nz], iy in [1,ny], ix in [0,nx/2]
	boost::array<MCArray3D::index,3> bases = {{0, 1, 1}};
	dat.reindex(bases);
	int xmax = (is_fftodd())
		? (nx-1)/2 + 1
		: (nx-2)/2;
	for (int iz = 1; iz <= nz; iz++) {
		for (int iy = 1; iy <= ny; iy++) {
			for (int ix = 0; ix <= xmax; ix++) {
				// next line multiplies by +/- 1
				dat[ix][iy][iz] *= static_cast<float>(-2*((ix+iy+iz)%2) + 1);
			}
		}
	}
	done_data();
	update();
	EXITFUNC;
}

// #G2#
EMData* EMData::zeropad_ntimes(int npad) {
	ENTERFUNC;
	if (is_complex()) 
		throw ImageFormatException("Zero padding complex images not supported");
	EMData* newimg = copy_head();
	int nxpad = npad*nx;
	int nypad = npad*ny;
	int nzpad = npad*nz;
	if (1 == ny) {
		// 1-d image, don't want to pad along y or z
		// Also, assuming that we can't have an image sized as nx=5, ny=1, nz=5.
		nypad = ny;
		nzpad = nz;
	} else if (nz == 1) {
		// 2-d image, don't want to pad along z
		nzpad = nz;
	}
	newimg->set_size(nxpad,nypad,nzpad);
	newimg->to_zero();
	size_t bytes = nx*sizeof(float);
	MArray3D dest = newimg->get_3dview();
	MArray3D src = this->get_3dview();
	int xstart = (nx != 1) ? (nxpad - nx)/2 + nx%2 : 0;
	int ystart = (ny != 1) ? (nypad - ny)/2 + ny%2 : 0;
	int zstart = (nz != 1) ? (nzpad - nz)/2 + nz%2 : 0;
	for (int iz = 0; iz < nz; iz++) {
		for (int iy = 0; iy < ny; iy++) {
			memcpy(&dest[xstart][iy+ystart][iz+zstart], &src[0][iy][iz], bytes);
		}
	}
	newimg->done_data();
	return newimg;
	EXITFUNC;
}

/** #G2#
Purpose: Create a new [npad-times zero-padded] fft-extended real image.
Method: Pad with zeros npad-times (npad may be 1, which is the default) and extend for fft,
return new real image.
Input: f real n-dimensional image
npad specify number of times to extend the image with zeros (default npad = 1, meaning no
padding)
Output: real image that may have been zero-padded and has been extended along x for fft.
 */
EMData* EMData::pad_fft(int npad) {
	ENTERFUNC;
	if (is_complex()) 
		throw ImageFormatException("Padding of complex images not supported");
	EMData* newimg = copy_head();
	newimg->to_zero();
	if (is_fftpadded() == false) {
		int nxpad = npad*nx;
		int nypad = npad*ny;
		int nzpad = npad*nz;
		if (1 == ny) {
			// 1-d image, don't want to pad along y or z
			// Also, assuming that we can't have an image sized as nx=5, ny=1, nz=5.
			nypad = ny;
			nzpad = nz;
		} else if (nz == 1) {
			// 2-d image, don't want to pad along z
			nzpad = nz;
		}
		size_t bytes;
		size_t offset;
		// Not currently padded, so we want to pad for ffts
		offset = 2 - nxpad%2;
		bytes = nx*sizeof(float);
		newimg->set_size(nxpad+offset, nypad, nzpad);
		newimg->to_zero();
		newimg->set_fftpad(true);
		newimg->set_attr("npad", npad);
		if (offset == 1)
			newimg->set_fftodd(true);
		MArray3D dest = newimg->get_3dview();
		MArray3D src = get_3dview();
		for (int iz = 0; iz < nz; iz++) {
			for (int iy = 0; iy < ny; iy++) {
				memcpy(&dest[0][iy][iz], &src[0][iy][iz], bytes);
			}
		}
	} else {
		// Image already padded, so we want to remove the padding
		// (Note: The npad passed in is ignored in favor of the one
		//  stored in the image.)
		npad = get_attr("npad");
		if (0 == npad) npad = 1;
		int nxold = (nx - 2 + is_fftodd())/npad; // using the value of is_fftodd() <- FIXME
	#ifdef _WIN32
		int nyold = _MAX(ny/npad, 1);
		int nzold = _MAX(nz/npad, 1);
	#else
		int nyold = std::max<int>(ny/npad, 1);
		int nzold = std::max<int>(nz/npad, 1);
	#endif	//_WIN32
		int bytes = nxold*sizeof(float);
		newimg->set_size(nxold, nyold, nzold);
		newimg->to_zero();
		newimg->set_fftpad(false);
		MArray3D dest = newimg->get_3dview();
		MArray3D src = this->get_3dview();
		for (int iz = 0; iz < nzold; iz++) {
			for (int iy = 0; iy < nyold; iy++) {
				memcpy(&dest[0][iy][iz], &src[0][iy][iz], bytes);
			}
		}
	}
	newimg->done_data();
	return newimg;
}


EMData *EMData::FH2F(int Size, float OverSamplekB)  // PRB
{
	int nx=get_xsize();
	int ny=get_ysize();
	int nz=get_zsize();
	float ScalFactor=4.1f;
	int Center = (int) floor((Size+1.0)/2.0 +.1);
	int CenterM= Center-1;
	int CountMax = (Center+1)*Center/2;

	int   *PermMatTr           = new int[CountMax];
	float *RValsSorted         = new float[CountMax];
	float *weightofkValsSorted = new float[CountMax];
	int *SizeReturned = new int[1];
	Util::Radialize(PermMatTr, RValsSorted,weightofkValsSorted,Size, SizeReturned);
	int RIntMax= SizeReturned[0];  // replaces CountMax; the latter should now never be used.
//	kVec2Use = (0:1/OverSamplek:RValsSorted(RIntMax)+1/OverSamplek); %   in pixels  (otherwise need *2*pi/Size)

	int mMax = (int) floor( ScalFactor*RValsSorted[RIntMax-1]+10.0);

	int    kIntMax  = 2+ (int) floor( RValsSorted[RIntMax-1]*OverSamplekB);
	float *kVec2Use = new float[kIntMax];
	for (int kk=0; kk<kIntMax; kk++){
		kVec2Use[kk]= ((float) kk)/OverSamplekB;}



#ifdef DEBUG
	printf("nx=%d, ny=%d, nz=%d Center=%d mMax=%d CountMax=%d kIntMax=%d Centerm1=%d  Size=%d\n\n",
	    nx,ny,nz, Center, mMax, CountMax, kIntMax,  CenterM, Size);
#endif

	MArray2D rhoOfkmB = get_2dview();

//     check mMax's are equal
//     check kIntMax's are equal

	if ( (nx==2*(mMax+1)) && (ny==kIntMax) &&(nz==1) ) {



	EMData* tempCopy = copy();
	tempCopy->set_size(2*(mMax+1),RIntMax);
	tempCopy->to_zero();
	MArray2D rhoOfkandm = tempCopy->get_2dview();
//	float rhoOfkandm = new MArray2D
//	printf("rhoOfkandm \n");
	for (int mr=0; mr <2*(mMax+1); mr++){
		float *Row= new float[kIntMax];
		float *RowOut= new float[RIntMax];
		for (int ii=0; ii<kIntMax; ii++){ Row[ii]=rhoOfkmB[mr][ii];}
		Util::spline_mat(kVec2Use, Row, kIntMax,  RValsSorted, RowOut, RIntMax ); 
		for (int ii=0; ii<RIntMax; ii++){
			rhoOfkandm[mr][ii] = RowOut[ii];
//			printf("%3.3f  ",RowOut[ii]);
		}
//		printf(" \n");
//		rhoOfkandm(m+1,:) = spline(kVec2Use,rhoOfkmBReIm(m+1,1:kIntMax),kIntMax,RValsSorted);
	}
	tempCopy->done_data();

//          So far so good PRB ....

	EMData* outCopy = tempCopy -> copy();
	outCopy->set_size(2*Size,Size,1);
	outCopy->to_zero();
	MArray2D ImBWfftRm = outCopy->get_2dview();

	int Count =0, kInt, kIntm1;
	complex <float> ImfTemp;
	float kValue, thetak;
	
	for (int jkx=0; jkx <Center; jkx++) { // These index the outputted picture
		for (int jky=0; jky<=jkx; jky++){
			kInt = PermMatTr[Count];
			kIntm1= kInt-1;
			Count++;
			float fjkx = float(jkx);
			float fjky = float(jky);

			kValue =sqrt(fjkx*fjkx +  fjky*fjky )  ;
//        		mMaxR= floor(ScalFactor*kValue +10);

 //                   How many copies

			thetak = atan2(fjky,fjkx);
			ImfTemp = rhoOfkandm[0][kIntm1];
        		for (int mm= 1; mm <mMax;mm++) {  // The index for m
				complex <float> fact(0,-mm*thetak);
				complex <float> expfact= exp(fact);
				complex <float> tempRho(rhoOfkandm[2*mm][kIntm1],rhoOfkandm[2*mm+1][kIntm1]);
				ImfTemp +=   expfact * tempRho + float(1-2*(mm%2))  *conj(expfact*tempRho);//pow(float(-1),mm)
        		}
 			ImBWfftRm[2*(CenterM+jkx)][CenterM+jky]   = ImfTemp.real();
			ImBWfftRm[2*(CenterM+jkx)+1][CenterM+jky] = ImfTemp.imag();
//			printf("jkx=%d, jky=%d; %f + %f i \n",jkx,jky,ImfTemp.real(), ImfTemp.imag());

			if (jky>0) {
				thetak = atan2(-fjky,fjkx);
				ImfTemp = rhoOfkandm[0][kIntm1];
				for (int mm= 1; mm<mMax; mm++) { // The index for m
					complex <float> fact(0,-mm*thetak);
					complex <float> expfact= exp(fact);
					complex <float> tempRho(rhoOfkandm[2*mm][kIntm1],rhoOfkandm[2*mm+1][kIntm1]);
					ImfTemp +=   expfact * tempRho +  float(1-2*(mm%2))  *conj(expfact*tempRho);
				}
				ImBWfftRm[2*(CenterM+jkx)][CenterM-jky]   = ImfTemp.real();
				ImBWfftRm[2*(CenterM+jkx)+1][CenterM-jky] = ImfTemp.imag();
			}

			if (jkx>0) {
            			thetak = atan2(fjky,-fjkx);
				ImfTemp = rhoOfkandm[0][kIntm1];
				for (int mm= 1; mm<mMax; mm++) { // The index for m
					complex <float> fact(0,-mm*thetak);
					complex <float> expfact= exp(fact);
					complex <float> tempRho(rhoOfkandm[2*mm][kIntm1],rhoOfkandm[2*mm+1][kIntm1]);
					ImfTemp +=   expfact * tempRho +  float(1-2*(mm%2)) *conj(expfact*tempRho);
				}
				ImBWfftRm[2*(CenterM-jkx)  ][CenterM+jky] = ImfTemp.real();
				ImBWfftRm[2*(CenterM-jkx)+1][CenterM+jky] = ImfTemp.imag();
			}

 			if (jkx>0 && jky>0) {
				thetak = atan2(-fjky,-fjkx);
				ImfTemp = rhoOfkandm[0][kIntm1];
				for (int mm= 1; mm<mMax; mm++) {  // The index for m
					complex <float> fact(0,-mm*thetak);
					complex <float> expfact= exp(fact);
					complex <float> tempRho(rhoOfkandm[2*mm][kIntm1],rhoOfkandm[2*mm+1][kIntm1]);
					ImfTemp +=   expfact * tempRho +  float(1-2*(mm%2)) *conj(expfact*tempRho);
				}
				ImBWfftRm[2*(CenterM-jkx)  ][CenterM-jky] = ImfTemp.real();
				ImBWfftRm[2*(CenterM-jkx)+1][CenterM-jky] = ImfTemp.imag();
			}

			if (jky< jkx) {
				thetak = atan2(fjkx,fjky);
				ImfTemp = rhoOfkandm[0][kIntm1];
				for (int mm= 1; mm<mMax; mm++){ // The index for m
					complex <float> fact(0,-mm*thetak);
					complex <float> expfact= exp(fact);
					complex <float> tempRho(rhoOfkandm[2*mm][kIntm1],rhoOfkandm[2*mm+1][kIntm1]);
					ImfTemp +=   expfact * tempRho +  float(1-2*(mm%2)) *conj(expfact*tempRho);
				}
				ImBWfftRm[2*(CenterM+jky)  ][CenterM+jkx] = ImfTemp.real();
				ImBWfftRm[2*(CenterM+jky)+1][CenterM+jkx] = ImfTemp.imag();

				if (jky>0){
					thetak = atan2(fjkx,-fjky);
					ImfTemp = rhoOfkandm[0][kIntm1];
					for (int mm= 1; mm <mMax; mm++) { // The index for m
						complex <float> fact(0,-mm*thetak);
						complex <float> expfact= exp(fact);
						complex <float> tempRho(rhoOfkandm[2*mm][kIntm1],rhoOfkandm[2*mm+1][kIntm1]);
						ImfTemp +=  expfact * tempRho +  float(1-2*(mm%2)) *conj(expfact*tempRho);
					}
					ImBWfftRm[2*(CenterM-jky)  ][CenterM+jkx] = ImfTemp.real();
					ImBWfftRm[2*(CenterM-jky)+1][CenterM+jkx] = ImfTemp.imag();
				}

				 if (jkx>0) {
					 thetak = atan2(-fjkx,fjky);
					 ImfTemp = rhoOfkandm[0][kIntm1];
					for (int mm= 1; mm <mMax; mm++) { // The index for m
						complex <float> fact(0,-mm*thetak);
						complex <float> expfact= exp(fact);
						complex <float> tempRho(rhoOfkandm[2*mm][kIntm1],rhoOfkandm[2*mm+1][kIntm1]);
						ImfTemp +=  expfact * tempRho +  float(1-2*(mm%2)) *conj(expfact*tempRho);
 					}
					ImBWfftRm[2*(CenterM+jky)  ][CenterM-jkx] = ImfTemp.real();
					ImBWfftRm[2*(CenterM+jky)+1][CenterM-jkx] = ImfTemp.imag();
 				}

	 			if (jkx>0 && jky>0) {
					thetak = atan2(-fjkx,-fjky);
					ImfTemp = rhoOfkandm[0][kIntm1];
					for (int mm= 1; mm <mMax; mm++) { // The index for m
						complex <float> fact(0,-mm*thetak);
						complex <float> expfact= exp(fact);
						complex <float> tempRho(rhoOfkandm[2*mm][kIntm1],rhoOfkandm[2*mm+1][kIntm1]);
						ImfTemp +=  expfact * tempRho +  float(1-2*(mm%2)) *conj(expfact*tempRho);
					}
					ImBWfftRm[2*(CenterM-jky)  ][CenterM-jkx] = ImfTemp.real();
					ImBWfftRm[2*(CenterM-jky)+1][CenterM-jkx] = ImfTemp.imag();
 				}
 			} // ends jky <jkx


		} // ends jky
	} // ends jkx
	outCopy->done_data();
	outCopy->set_complex(true);
	if(outCopy->get_ysize()==1 && outCopy->get_zsize()==1) {
		outCopy->set_complex_x(true);
	}
	outCopy->set_ri(true);
	outCopy->set_FH(false);
	outCopy->set_fftodd(true);
	outCopy->set_shuffled(true);
	return outCopy;
	} else {
		LOGERR("can't be an FH image not this size");
		throw ImageFormatException("something strange about this image: not a FH");

	}
}  // ends FH2F


EMData *EMData::real2FH(float OverSamplekB) // PRB
{
	int nx=get_xsize();
	int ny=get_ysize();
	int nz=get_zsize();
	int Center = (int) floor( (nx+1.0)/2.0 +.01);
#ifdef DEBUG
	printf("nx=%d, ny=%d, nz=%d Center=%d\n", nx,ny,nz, Center);
#endif	//DEBUG
	float ScalFactor=4.1f;
	gsl_set_error_handler_off();

	if ( (nz==1) && (nx==ny) && (!is_complex())  && (Center*2)==(nx+1)){
#ifdef DEBUG
		printf("entered if \n");fflush(stdout);
#endif	//DEBUG
		MArray2D ImBW = this ->get_2dview();
		int Size=nx;
		int iMax = (int) floor( (Size-1.0)/2 +.01);
		int CountMax = (iMax+2)*(iMax+1)/2;
		int *PermMatTr  = new int[CountMax];
		float *RValsSorted  = new float[CountMax];
		float *weightofkValsSorted = new float[CountMax];
		int *SizeReturned = new int[1];
		Util::Radialize(PermMatTr, RValsSorted,weightofkValsSorted,Size, SizeReturned);
	  	int RIntMax= SizeReturned[0];

		int mMax = (int) floor( ScalFactor*RValsSorted[RIntMax-1]+10.0);

		int kIntMax=2+ (int) floor( RValsSorted[RIntMax-1]*OverSamplekB);
		float *kVec2Use= new float[kIntMax];
		for (int kk=0; kk<kIntMax; kk++){
			kVec2Use[kk]= ((float) kk)/OverSamplekB;}

		float *krVec= new float[kIntMax*RIntMax];
		int Count=0;
		for (int jk=0; jk<kIntMax; jk++ ){
			for (int jR=0; jR<RIntMax; jR++ ){
				krVec[Count]=2.0f*M_PI*RValsSorted[jR]
					*kVec2Use[jk]/( (float) Size);
				Count++;
//				printf("krVec[%d]=%f \n",Count,krVec[Count-1]);fflush(stdout);
		}} // end building up krVec
		float krVecMin= kVec2Use[1]*RValsSorted[1];
		float krVecMax = krVec[kIntMax*RIntMax-1]+krVecMin;
		int Number2Use = (int) floor(OverSamplekB*krVecMax+1.0);
		float *krVec2Use      = new float[Number2Use+1];
		float *sampledBesselJ = new float[Number2Use+1];
#ifdef DEBUG
		printf("Size=%d, iMax=%d, SizeReturned=%d, RIntMax=%d, \n"
		      "mMax=%d, kIntMax=%d, krVecMin=%f, krVecMax=%f,  Number2Use=%d  \n\n",
			Size, iMax, SizeReturned[0], RIntMax, mMax, kIntMax,
			       krVecMin,krVecMax,Number2Use);fflush(stdout);
#endif	//DEBUG
		for (int jkr=0; jkr<= Number2Use; jkr++) {
			krVec2Use[jkr] =((float)jkr)*krVecMax/
			            ((float)Number2Use);
//			printf("krVec2Use[%d]=%f \n",jkr+1,krVec2Use[jkr]);fflush(stdout);
		}


		EMData* FH = copy(); // glibc detected ** malloc(); memory corruption
//		printf("finished O \n");fflush(stdout);
		FH->set_size(2*(mMax+1),kIntMax);
		FH->to_zero();
		MArray2D rhoOfkmB = FH->get_2dview();

		int CenterM= Center-1; // to convert from Matlab to C++
		complex <float> *rhoOfRandmTemp = new complex <float>[RIntMax];
		complex <float> rhoTemp;
		int PCount=0;


		for (int m=0; m <=mMax; m++){
		//    if m==mMax, tic, end
			complex <float> tempF(0.0f,-1.0f);
			complex <float> overallFactor = pow(tempF,m);  //(-i)^m ;  % I dropped off the 2 pi
			complex <float> mI(0.0f,static_cast<float>(m));
			for (int ii=0; ii< RIntMax; ii++){ rhoOfRandmTemp[ii]=0;}
			for (int jx=0; jx <Center ; jx++) {
				for (int jy=0; jy <=jx; jy++){
					float fjx=float(jx);
					float fjy= float(jy);
          				Count = (jx*jx+jx)/2 +1 +jy;
					PCount = PermMatTr[Count-1];
//					printf("PCount=%d, Count=%d \n", PCount, Count);
  				        rhoTemp =  complex <float> (ImBW[CenterM+jx][CenterM+jy]) *exp(mI* complex <float> (atan2(+fjy,+fjx)))
				         +   complex <float> (ImBW[CenterM+jx][CenterM-jy]) * exp(mI*complex <float>(atan2(-fjy,+fjx)))
				         +   complex <float> (ImBW[CenterM-jx][CenterM+jy]) * exp(mI*complex <float>(atan2(+fjy,-fjx)))
				         +   complex <float> (ImBW[CenterM-jx][CenterM-jy]) * exp(mI*complex <float>(atan2(-fjy,-fjx)))
			               	 +   complex <float> (ImBW[CenterM+jy][CenterM+jx]) * exp(mI*complex <float>(atan2(+fjx,+fjy)))
					 +   complex <float> (ImBW[CenterM+jy][CenterM-jx]) * exp(mI*complex <float>(atan2(-fjx,+fjy)))
					 +   complex <float> (ImBW[CenterM-jy][CenterM+jx]) * exp(mI*complex <float>(atan2(+fjx,-fjy)))
					 +   complex <float> (ImBW[CenterM-jy][CenterM-jx]) * exp(mI*complex <float>(atan2(-fjx,-fjy)));
            				if (((jx+jy)==0)&&(m>0) ){
						rhoTemp=0;}
//			printf("m=%d, jx=%d, jy=%d, rhoTemp= %f+ %f i\n", m,jx,jy,(rhoTemp.real()), (rhoTemp.imag()) );fflush(stdout);
//			{" %f,%f %f,%f %f,%f %f,%f \n",
//			       ImBW[CenterM+jx][CenterM+jy] ,ImBW[CenterM+jx][CenterM-jy]  , ImBW[CenterM-jx][CenterM+jy] ,ImBW[CenterM-jx][CenterM-jy],
//			       ImBW[CenterM+jy][CenterM+jx] ,ImBW[CenterM+jy][CenterM-jx]  , ImBW[CenterM-jy][CenterM+jx] ,ImBW[CenterM-jy][CenterM-jx]);
            				rhoOfRandmTemp[PCount-1] +=
				            rhoTemp/((float)pow(2.,(int)( (jx==0)  +(jy==0)+ (jy==jx))));

			}} // end walk through lattice
//			printf("\n m=%d rhoOfRandmTemp" ,m  );fflush(stdout);
//			for (int ss=0; ss< RIntMax; ss++){
//				printf(" %3.1f+ %3.1fi \t",(rhoOfRandmTemp[ss].real()), (rhoOfRandmTemp[ss].imag())   );fflush(stdout);}

// calculate product
			float tempp;
//			printf("\n m=%d sampledBesselJ" ,m  );fflush(stdout);
			for (int st=0; st<= Number2Use; st++){
				tempp=krVec2Use[st];
				sampledBesselJ[st] = static_cast<float>(gsl_sf_bessel_Jn(m,tempp));
//				printf(" %3.2f  \t",sampledBesselJ[st]   );fflush(stdout);
			} // good so far
//			sampledBesselJ  = BesselJ(m,krVec2Use);
			float *tempMB = new float [kIntMax*RIntMax];
			Util::spline_mat(krVec2Use, sampledBesselJ, Number2Use+1,krVec,tempMB,kIntMax*RIntMax ); 
//			printf("\n tempMB m=%d y2" ,m  );fflush(stdout);
			complex <float> *rowV = new complex <float> [kIntMax];

//			for (int st=0; st< kIntMax*RIntMax; st++){printf(" %3.2f  \t",tempMB[st]   );fflush(stdout);} // good so far

//   tempMB,krVec is in blocks of RIntMax
//			printf("\n rowV m=%d \t" ,m  );fflush(stdout);
			for (int st=0; st < kIntMax; st++) {
					rowV[st]=0;
					for (int sv=0; sv < RIntMax; sv++) {
						rowV[st]+=  rhoOfRandmTemp[sv] *tempMB[sv+st*RIntMax];
					}
					 rowV[st] *= overallFactor;
//					printf(" %1.3f +%1.3fi \t" , rowV[st].real(), rowV[st].imag() );fflush(stdout);
			}
			for (int st=0; st < kIntMax; st++) {
					rhoOfkmB[2*m  ][st] = rowV[st].real();
					rhoOfkmB[2*m+1][st] = rowV[st].imag();
			}
// 			rowV = overallFactor*rhoOfRandmTemp*tempMBB;
//			rhoOfkmB(m+1,1:kIntMax) = rowV ;

//			if m==mMax, toc, end

// %'final interpolation'
// %     rhoOfkm(m+1,:) = spline(kVec2Use,rowV,RValsSorted); ;


		} // ends m loop
		done_data();
		FH-> done_data();
		FH->set_complex(true);
		if(FH->get_ysize()==1 && FH->get_zsize()==1) {
			FH->set_complex_x(true);
		}
	    	FH->set_ri(true);
	    	FH->set_FH(true);
	    	FH->set_fftodd(true);
		return FH;
	} else {
		LOGERR("2D real square odd image expected.");
		throw ImageFormatException("2D real square odd image expected.");
	}
}


EMData *EMData::do_fft()
{
	ENTERFUNC;

	if ( is_complex() ) {
		LOGERR("real image expected. Input image is complex image.");
		throw ImageFormatException("real image expected. Input image is complex image.");
	}

	int nxreal = nx;
	int offset = 2 - nx%2;
	int nx2 = nx + offset;
	EMData* dat = copy_head();
	dat->set_size(nx2, ny, nz);
	dat->to_zero();
	if (offset == 1)
		dat->set_fftodd(true);

	float *d = dat->get_data();
	EMfft::real_to_complex_nd(rdata, d, nxreal, ny, nz);

#if 0 // Remove funky reordering
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
		if( t )
		{
			delete[]t;
			t = 0;
		}
	}
#endif // 0

	dat->done_data();
	dat->set_complex(true);
	if(dat->get_ysize()==1 && dat->get_zsize()==1) {
		dat->set_complex_x(true);
	}
	dat->set_ri(true);

	int i = flags;
	done_data();
	flags = i & ~EMDATA_BUSY;

	EXITFUNC;
	return dat;
}

EMData *EMData::do_fft_inplace()
{
	ENTERFUNC;

	if ( is_complex() ) {
		LOGERR("real image expected. Input image is complex image.");
		throw ImageFormatException("real image expected. Input image is complex image.");
	}
	
	size_t offset;
	int nxreal;
	if (!is_fftpadded()) {
		// need to extend the matrix along x
		// meaning nx is the un-fftpadded size
		nxreal = nx;
		MArray3D arr = get_3dview();
		offset = 2 - nx%2;
		if (1 == offset) set_fftodd(true);
		int nxnew = nx + offset;
		set_size(nxnew, ny, nz);
		// now need to relocate the data in rdata
		MArray3D dest = get_3dview();
		array<std::size_t,3> dims = {{nxreal, ny, nz}};
		MArray3D src(get_data(), dims, boost::fortran_storage_order());
		for (int iz = nz-1; iz >= 0; iz--)
			for (int iy = ny-1; iy >= 0; iy--)
				for (int ix = nxreal-1; ix >= 0; ix--)
					dest[ix][iy][iz] = src[ix][iy][iz];
		set_fftpad(true);
	} else {
		offset = is_fftodd() ? 1 : 2;
		nxreal = nx - offset;
	}
	EMfft::real_to_complex_nd(rdata, rdata, nxreal, ny, nz);
	done_data();
	set_complex(true);
	if(ny==1 && nz==1) {
		set_complex_x(true);
	}
	set_ri(true);
	int i = flags;
	done_data();
	flags = i & ~EMDATA_BUSY;
	EXITFUNC;
	return this;
}

EMData *EMData::do_ift()
{
	ENTERFUNC;

	if (!is_complex()) {
		LOGERR("complex image expected. Input image is real image.");
		throw ImageFormatException("complex image expected. Input image is real image.");
	}

	if (!is_ri()) {
		LOGWARN("run IFT on AP data, only RI should be used. ");
	}

	EMData* dat = copy_head();
	dat->set_size(nx, ny, nz);
	ap2ri();

	float *d = dat->get_data();
	int ndim = get_ndim();

	if (ndim >= 2) {
		memcpy((char *) d, (char *) rdata, nx * ny * nz * sizeof(float));
	}


	int offset = is_fftodd() ? 1 : 2;
	if (ndim == 1) {
		EMfft::complex_to_real_nd(rdata, d, nx - offset, ny, nz);
	}
#if 0 // remove funky reordering
	else if (ndim == 2 ) {
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
		char *t = new char[(nx + offset) * sizeof(float)];
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

		if( t )
		{
			delete[]t;
			t = 0;
		}
	}
#endif // 0

	if (ndim >= 2) {
		EMfft::complex_to_real_nd(d, d, nx - offset, ny, nz);

		size_t row_size = (nx - offset) * sizeof(float);
		for (int i = 1; i < ny * nz; i++) {
			memmove((char *) &d[i * (nx - offset)], (char *) &d[i * nx], row_size);
		}
	}


	// SCALE the inverse FFT
	float scale = 1.0f / ((nx - offset) * ny * nz);
	dat->mult(scale);
	dat->done_data();
#if 1
	dat->set_size(nx - offset, ny, nz);
#endif
	dat->update();
	dat->set_complex(false);
	if(dat->get_ysize()==1 && dat->get_zsize()==1) {
		dat->set_complex_x(false);
	}
	dat->set_ri(false);

	int i = flags;
	done_data();
	flags = i & ~EMDATA_BUSY;

	EXITFUNC;
	return dat;
}

EMData *EMData::do_ift_inplace()
{
	ENTERFUNC;

	if (!is_complex()) {
		LOGERR("complex image expected. Input image is real image.");
		throw ImageFormatException("complex image expected. Input image is real image.");
	}

	if (!is_ri()) {
		LOGWARN("run IFT on AP data, only RI should be used. ");
	}

	ap2ri();

	int ndim = get_ndim();


	// turn off data shuffling if we're doing an in-place transform
	int offset = is_fftodd() ? 1 : 2;
	if (ndim == 1) {
		EMfft::complex_to_real_nd(rdata, rdata, nx - offset, ny, nz);
	}

	if (ndim >= 2) {
		EMfft::complex_to_real_nd(rdata, rdata, nx - offset, ny, nz);
	}

	// SCALE the inverse FFT
	float scale = 1.0f / ((nx - offset) * ny * nz);
	mult(scale);
	done_data();
	update();
	set_complex(false);
	if(ny==1 && nz==1) {
		set_complex_x(false);
	}
	set_ri(false);

	int i = flags;
	done_data();
	flags = i & ~EMDATA_BUSY;

	EXITFUNC;
	return this;
}

EMData* EMData::get_fft_amplitude2D()
{
	ENTERFUNC;

//	int ndim = get_ndim();
	if (!is_complex()) {
		LOGERR("complex image expected. Input image is real image.");
		throw ImageFormatException("complex image expected. Input image is a real image.");
	}
	if (nz>1) {
		LOGERR("2D image expected. Input image is 3D");
		throw ImageFormatException("2D odd square complex image"
			" expected Input image is 3D.");
	}


	int nx2 = nx/2;

	EMData *dat = copy_head();

	dat->set_size(nx2, ny, nz);
	dat->to_zero();
	MArray2D rnewdata = this-> get_2dview();
	MArray2D d        = dat -> get_2dview(); // pointer to an emData object

	float temp=0;

	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx2; i++) {
			temp  = rnewdata[2*i  ][j] *rnewdata[2*i  ][j] ;
			temp += rnewdata[2*i+1][j] *rnewdata[2*i+1][j] ;
//			printf("i= %d, j=%d, temp=%f",i,j,temp);
			d[i][j] = sqrt(temp);
		}
	}


	done_data();

	dat->done_data();
	dat->update();
	dat->set_complex(false);
	dat->set_ri(false);

	EXITFUNC;
	return dat;
}


EMData* EMData::get_fft_amplitude()
{
	ENTERFUNC;

	if (!is_complex()) {
		LOGERR("complex image expected. Input image is real image.");
		throw ImageFormatException("complex image expected. Input image is a real image.");
	}

	ri2ap();

	int nx2 = nx - 2;
	EMData *dat = copy_head();
	dat->set_size(nx2, ny, nz);
	dat->to_zero();

	float *d = dat->get_data();

	int ndim = get_ndim();

	if (ndim == 3) {
		for (int k = 1; k < nz; k++) {
			for (int j = 1; j < ny; j++) {
				for (int i = 0; i < nx2/2; i++) {
					d[k*nx2*ny+j*nx2+nx2/2+i] = rdata[k*nx*ny+j*nx+2*i];
					d[(nz-k)*nx2*ny+(ny-j)*nx2+nx2/2-i] = rdata[k*nx*ny+j*nx+2*i];
				}
			}
		}
	}
	else if (ndim == 2) {
		for (int j = 1; j < ny; j++) {
			for (int i = 0; i < nx2/2; i++) {
				d[j*nx2+nx2/2+i] = rdata[j*nx+2*i];
				d[(ny-j)*nx2+nx2/2-i] = rdata[j*nx+2*i];
			}
		}
	}

	done_data();

	dat->done_data();
	dat->update();
	dat->set_complex(false);
	if(dat->get_ysize()==1 && dat->get_zsize()==1) {
		dat->set_complex_x(false);
	}
	dat->set_ri(false);

	EXITFUNC;
	return dat;
}

EMData* EMData::get_fft_phase()
{
	ENTERFUNC;

	if (!is_complex()) {
		LOGERR("complex image expected. Input image is real image.");
		throw ImageFormatException("complex image expected. Input image is a real image.");
	}

	ri2ap();

	int nx2 = nx - 2;
	EMData *dat = copy_head();
	dat->set_size(nx2, ny, nz);
	dat->to_zero();

	float *d = dat->get_data();

	int ndim = get_ndim();
	if (ndim == 3) {
		for (int k = 1; k < nz; k++) {
			for (int j = 1; j < ny; j++) {
				for (int i = 0; i < nx2/2; i++) {
					d[k*nx2*ny+j*nx2+nx2/2+i] = rdata[k*nx*ny+j*nx+2*i+1];
					d[(nz-k)*nx2*ny+(ny-j)*nx2+nx2/2-i] = -rdata[k*nx*ny+j*nx+2*i+1];
				}
			}
		}
	}
	else if (ndim == 2) {
		for (int j = 1; j < ny; j++) {
			for (int i = 0; i < nx2/2; i++) {
				d[j*nx2+nx2/2+i] = rdata[j*nx+2*i+1];
				d[(ny-j)*nx2+nx2/2-i] = -rdata[j*nx+2*i+1];
			}
		}
	}
	done_data();

	dat->done_data();
	dat->update();
	dat->set_complex(false);
	if(dat->get_ysize()==1 && dat->get_zsize()==1) {
		dat->set_complex_x(false);
	}
	dat->set_ri(false);

	EXITFUNC;
	return dat;
}

vector < float > EMData::calc_hist(int hist_size, float histmin, float histmax)
{
	ENTERFUNC;

	static size_t prime[] = { 1, 3, 7, 11, 17, 23, 37, 59, 127, 253, 511 };

	if (histmin == histmax) {
		histmin = get_attr("minimum");
		histmax = get_attr("maximum");
	}

	vector <float> hist(256, 0.0);

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

	size_t di = 0;
	float norm = 0;
	size_t n = hist.size();

	for (int k = p0; k <= p1; ++k) {
		if (flags & EMDATA_COMPLEX) {
			di = prime[k] * 2;
		}
		else {
			di = prime[k];
		}

		norm += (float)size / (float) di;
		float w = (float)n / (histmax - histmin);

		for(size_t i=0; i<=size-di; i += di) {
			int j = Util::round((rdata[i] - histmin) * w);
			if (j >= 0 && j < (int) n) {
				hist[j] += 1;
			}
		}
	}

	for (size_t i = 0; i < hist.size(); ++i) {
		if (norm != 0) {
			hist[i] = hist[i] / norm;
		}
	}
	
	return hist;
	
	EXITFUNC;
}



EMData *EMData::little_big_dot(EMData * with, bool do_sigma)
{
	ENTERFUNC;

	if (get_ndim() > 2) {
		throw ImageDimensionException("1D/2D only");
	}

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

	EXITFUNC;
	return ret;
}

std::string EMData::render_amp8(int x0, int y0, int ixsize, int iysize,
						 int bpl, float scale, int mingray, int maxgray,
						 float render_min, float render_max,int asrgb)
{
	ENTERFUNC;

	if (get_ndim() != 2) {
		throw ImageDimensionException("2D only");
	}

	if (is_complex()) {
		ri2ap();
	}

	if (render_max <= render_min) {
		render_max = render_min + 0.01f;
	}

	if (asrgb) asrgb=3;
	else asrgb=1;

	std::string ret=std::string();
	ret.resize(iysize*bpl);
	ret.assign(iysize*bpl,char(mingray));
	unsigned char *data=(unsigned char *)ret.data();

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

	if (xmin < 0) xmin = 0;
	if (ymin < 0) ymin = 0;
	if (xsize > ixsize) xsize = ixsize;
	if (ymax > iysize) ymax = iysize;

	int lmax = nx * ny - 1;

	if (is_complex()) {
		if (dsx != -1) {
			int l = y0 * nx;
			for (int j = ymax; j >= ymin; j--) {
				int ll = x0;
				for (int i = xmin; i < xsize; i++) {
					if (l + ll > lmax || ll >= nx - 2) break;

					int k = 0;
					if (ll >= nx / 2) {
						if (l >= (ny - inv_scale) * nx) k = 2 * (ll - nx / 2) + 2;
						else k = 2 * (ll - nx / 2) + l + 2 + nx;
					}
					else k = nx * ny - (l + 2 * ll) - 2;
					float t = rdata[k];
					if (t <= rm)  k = mingray;
					else if (t >= render_max) k = maxgray;
					else {
						k = (int) (gs * (t - render_min));
						k += mingray;
					}
					data[i * asrgb + j * bpl] = static_cast < unsigned char >(k);
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
						if (l >= (ny * nx - nx)) k = 2 * (ll - nx / 2) + 2;
						else k = 2 * (ll - nx / 2) + l + 2 + nx;
					}
					else k = nx * ny - (l + 2 * ll) - 2;

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
					data[i * asrgb + j * bpl] = static_cast < unsigned char >(k);
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
					if (t <= rm) k = mingray;
					else if (t >= render_max) k = maxgray;
					else {
						k = (int) (gs * (t - render_min));
						k += mingray;
					}
					data[i * asrgb + j * bpl] = static_cast < unsigned char >(k);
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
					if (l > lmax) break;
					int k = 0;
					float t = rdata[l];
					if (t <= rm) k = mingray;
					else if (t >= render_max) k = maxgray;
					else {
						k = (int) (gs * (t - render_min));
						k += mingray;
					}
					data[i * asrgb + j * bpl] = static_cast < unsigned char >(k);
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

	// this replicates r -> g,b
	if (asrgb==3) {
		for (int j=ymin*bpl; j<=ymax*bpl; j+=bpl) {
			for (int i=xmin; i<xsize*3; i+=3) {
				data[i+j+1]=data[i+j+2]=data[i+j];
			}
		}
	}

	done_data();
	EXITFUNC;

    //	return PyString_FromStringAndSize((const char*) data,iysize*bpl);
	return ret;
}


void EMData::render_amp24( int x0, int y0, int ixsize, int iysize,
						  int bpl, float scale, int mingray, int maxgray,
						  float render_min, float render_max, void *ref,
						  void cmap(void *, int coord, unsigned char *tri))
{
	ENTERFUNC;

	if (get_ndim() != 2) {
		throw ImageDimensionException("2D only");
	}

	if (is_complex()) {
		ri2ap();
	}

	if (render_max <= render_min) {
		render_max = render_min + 0.01f;
	}

	std::string ret=std::string();
	ret.resize(iysize*bpl);
	unsigned char *data=(unsigned char *)ret.data();

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
	EXITFUNC;
}

vector<float> EMData::calc_az_dist(int n, float a0, float da, float rmin, float rmax)
{
	ENTERFUNC;

	if (get_ndim() > 2) {
		throw ImageDimensionException("no 3D image");
	}

	float *yc = new float[n];

	vector<float>	vd(n);
	for (int i = 0; i < n; i++) {
		yc[i] = 0.00001f;
	}

	if (is_complex()) {
		int c = 0;
		for (int y = 0; y < ny; y++) {
			for (int x = 0; x < nx; x += 2, c += 2) {
				float x1 = x / 2.0f;
				float y1 = y - ny / 2.0f;
				float r = (float)hypot(x1, y1);

				if (r >= rmin && r <= rmax) {
					float a = 0;

					if (y != ny / 2 || x != 0) {
						a = (atan2(y1, x1) - a0) / da;
					}

					int i = static_cast < int >(floor(a));
					a -= i;

					if (i == 0) {
						vd[0] += rdata[c] * (1.0f - a);
						yc[0] += (1.0f - a);
					}
					else if (i == n - 1) {
						vd[n - 1] += rdata[c] * a;
						yc[n - 1] += a;
					}
					else if (i > 0 && i < (n - 1)) {
						float h = 0;
						if (is_ri()) {
							h = (float)hypot(rdata[c], rdata[c + 1]);
						}
						else {
							h = rdata[c];
						}

						vd[i] += h * (1.0f - a);
						yc[i] += (1.0f - a);
						vd[i + 1] += h * a;
						yc[i + 1] += a;
					}
				}
			}
		}
	}
	else {
		int c = 0;
		float half_nx = (nx - 1) / 2.0f;
		float half_ny = (ny - 1) / 2.0f;

		for (int y = 0; y < ny; y++) {
			for (int x = 0; x < nx; x++, c++) {
				float y1 = y - half_ny;
				float x1 = x - half_nx;
				float r = (float)hypot(x1, y1);

				if (r >= rmin && r <= rmax) {
					float a = 0;
					if (x1 != 0 || y1 != 0) {
						a = atan2(y1, x1);
						if (a < 0) {
							a += M_PI * 2;
						}
					}

					a = (a - a0) / da;
					int i = static_cast < int >(floor(a));
					a -= i;

					if (i == 0) {
						vd[0] += rdata[c] * (1.0f - a);
						yc[0] += (1.0f - a);
					}
					else if (i == n - 1) {
						vd[n - 1] += rdata[c] * a;
						yc[n - 1] += (a);
					}
					else if (i > 0 && i < (n - 1)) {
						vd[i] += rdata[c] * (1.0f - a);
						yc[i] += (1.0f - a);
						vd[i + 1] += rdata[c] * a;
						yc[i + 1] += a;
					}
				}
			}
		}
	}


	for (int i = 0; i < n; i++) {
		vd[i] /= yc[i];
	}

	done_data();
	if( yc )
	{
		delete[]yc;
		yc = 0;
	}
	
	return vd;

	EXITFUNC;
}


void EMData::ri2ap()
{
	ENTERFUNC;

	if (!is_complex() || !is_ri()) {
		return;
	}

	int size = nx * ny * nz;
	for (int i = 0; i < size; i += 2) {
		float f = (float)hypot(rdata[i], rdata[i + 1]);
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
	EXITFUNC;
}

void EMData::ap2ri()
{
	ENTERFUNC;

	if (!is_complex() || is_ri()) {
		return;
	}

	Util::ap2ri(rdata, nx * ny * nz);
	set_ri(true);
	flags |= EMDATA_NEEDUPD;
	EXITFUNC;
}


vector < float >EMData::calc_fourier_shell_correlation(EMData * with, float w)
{
	ENTERFUNC;

/*
 ******************************************************
 *DISCLAIMER
 * 08/16/05 P.A.Penczek
 * The University of Texas
 * Pawel.A.Penczek@uth.tmc.edu
 * Please do not modify the content of calc_fourier_shell_correlation
 ******************************************************/
/*
Fourier Ring/Shell Correlation
Purpose: Calculate CCF in Fourier space as a function of spatial frequency
         between a pair of 2-3D images. 
Method: Calculate FFT (if needed), calculate FSC.
Input:  f - real or complex 2-3D image
	g - real or complex 2-3D image
        w - float ring width
Output: 2D 3xk real image.
        k - length of FSC curve, depends on dimensions of the image and ring width
	1 column - FSC,
	2 column - normalized frequency [0,0.5]
	3 column - currently n /error of the FSC = 1/sqrt(n), where n is the number of Fourier
	           coefficients within given shell.
*/  
	int needfree=0, nx, ny, nz, nx2, ny2, nz2, ix, iy, iz, kz, ky;
	float  dx2, dy2, dz2, argx, argy, argz;

	if (!with) {
		throw NullPointerException("NULL input image");
	}

	
	EMData *f = this;
	EMData *g = with;
	
	nx  = f->get_xsize();
	ny  = f->get_ysize();
	nz  = f->get_zsize();

	if (ny==0 && nz==0) {
		throw ImageFormatException( "Cannot calculate FSC for 1D images");
	}

	if (!equalsize(f, g)) {
		LOGERR("FSC requires congruent images");
		throw ImageDimensionException("FSC requires congruent images");
	}

	if (f->is_complex()) nx = (nx - 2 + f->is_fftodd()); // nx is the real-space size of the input image
	int lsd2 = (nx + 2 - nx%2) ; // Extended x-dimension of the complex image

//  Process f if real
	EMData* fpimage = NULL;
	if(f->is_complex()) fpimage = f;
	else {fpimage= norm_pad_ft(f, false, false); needfree|=1;} // Extend and do the FFT if f is real


//  Process g if real
	EMData* gpimage = NULL;
	if(g->is_complex()) gpimage = g; 
	else {gpimage= norm_pad_ft(g, false, false); needfree|=2;} // Extend and do the FFT if f is real


	float *d1 = fpimage->get_data();
	float *d2 = gpimage->get_data();

	nx2=nx/2; ny2 = ny/2; nz2 = nz/2;
	dx2 = 1.0f/float(nx2)/float(nx2); 
	dy2 = 1.0f/float(ny2)/float(ny2);

#ifdef _WIN32
	dz2 = 1.0f / _MAX(float(nz2),1.0f)/_MAX(float(nz2),1.0f);
	
	int inc = Util::round(float( _MAX( _MAX(nx2,ny2),nz2) )/w );
#else
	dz2 = 1.0f/std::max(float(nz2),1.0f)/std::max(float(nz2),1.0f);
	
	int inc = Util::round(float(std::max(std::max(nx2,ny2),nz2))/w);
#endif	//_WIN32

	double *ret = new double[inc+1];
	double *n1  = new double[inc+1];
	double *n2  = new double[inc+1];
	float  *lr  = new float[inc+1];
	for (int i = 0; i <= inc; i++) {
		ret[i] = 0; n1[i] = 0; n2[i] = 0; lr[i]=0;
	}

	for ( iz = 0; iz <= nz-1; iz++) {
		if(iz>nz2) kz=iz-nz; else kz=iz; argz = float(kz*kz)*dz2;
		for ( iy = 0; iy <= ny-1; iy++) {
			if(iy>ny2) ky=iy-ny; else ky=iy; argy = argz + float(ky*ky)*dy2;
			for ( ix = 0; ix <= lsd2-1; ix+=2) {
			// Skip Friedel related values
			   if(ix>0 || (kz>=0 && (ky>=0 || kz!=0))) {
				argx = 0.5f*sqrt(argy + float(ix*ix)*0.25f*dx2);
				int r = Util::round(inc*2*argx);
				if(r <= inc) {
					int ii = ix + (iy  + iz * ny)* lsd2;
					ret[r] += d1[ii] * double(d2[ii]) + d1[ii + 1] * double(d2[ii + 1]);
					n1[r]  += d1[ii] * double(d1[ii]) + d1[ii + 1] * double(d1[ii + 1]);
					n2[r]  += d2[ii] * double(d2[ii]) + d2[ii + 1] * double(d2[ii + 1]);
					lr[r]  +=1;
				}
			   }
			}
		}
	}

	vector < float >result((inc+1)*3);

	for (int i = 0; i <= inc; i++) {
		if(lr[i]>0) {
			result[i]     = float(i)/float(2*inc);
			result[i+inc+1]           = float(ret[i] / (sqrt(n1[i] * n2[i])));
			result[i+2*(inc+1)] = lr[i]  /*1.0f/sqrt(float(lr[i]))*/;}
		else {
			result[i]           = 0.0f;
			result[i+inc+1]     = 0.0f;
			result[i+2*(inc+1)] = 0.0f;}
	}

	if( ret )
	{
		delete[]ret;
		ret = 0;
	}

	if( n1 )
	{
		delete[]n1;
		n1 = 0;
	}
	if( n2 )
	{
		delete[]n2;
		n2 = 0;
	}

	if (needfree&1)
	{
		if( fpimage )
		{
			delete fpimage;
			fpimage = 0;
		}
	}
	if (needfree&2)
	{
		if( gpimage )
		{
			delete gpimage;
			gpimage = 0;
		}
	}

	EXITFUNC;
	return result;
}


void EMData::add(float f,int keepzero)
{
	ENTERFUNC;

	if( is_real() )
	{
		if (f != 0) {
			flags |= EMDATA_NEEDUPD;
			size_t size = nx * ny * nz;
			if (keepzero) {
				for (size_t i = 0; i < size; i++) {
					if (rdata[i]) rdata[i] += f;
				}
			}
			else {
				for (size_t i = 0; i < size; i++) {
					rdata[i] += f;
				}
			}
		}
	}
	else if( is_complex() )
	{
		if( f!=0 )
		{
			flags |= EMDATA_NEEDUPD;
			size_t size = nx*ny*nz; //size of data
			if( keepzero )
			{
				for(size_t i=0; i<size; i+=2)
				{
					if (rdata[i]) rdata[i] += f;
				}
			}
			else
			{
				for(size_t i=0; i<size; i+=2)
				{
					rdata[i] += f;
				}
			}
		}
	}
	else
	{
		throw ImageFormatException("This image is neither a real nor a complex image.");
	}

	EXITFUNC;
}

//for add operation, real and complex image is the same
void EMData::add(const EMData & image)
{
	ENTERFUNC;

	if (nx != image.get_xsize() || ny != image.get_ysize() || nz != image.get_zsize()) {
		throw ImageFormatException( "images not same sizes");
	}
	else if( (is_real()^image.is_real()) == true )
	{
		throw ImageFormatException( "not support add between real image and complex image");
	}
	else {
		flags |= EMDATA_NEEDUPD;
		const float *src_data = image.get_data();
		int size = nx * ny * nz;

		for (int i = 0; i < size; i++) {
			rdata[i] += src_data[i];
		}
	}
	EXITFUNC;
}

void EMData::sub(float f)
{
	ENTERFUNC;

	if( is_real() )
	{
		if (f != 0) {
			flags |= EMDATA_NEEDUPD;
			size_t size = nx * ny * nz;
			for (size_t i = 0; i < size; i++) {
				rdata[i] -= f;
			}
		}
	}
	else if( is_complex() )
	{
		if( f != 0 )
		{
			flags |= EMDATA_NEEDUPD;
			size_t size = nx * ny * nz;
			for( size_t i=0; i<size; i+=2 )
			{
				rdata[i] -= f;
			}
		}
	}
	else
	{
		throw ImageFormatException("This image is neither a real nor a complex image.");
	}

	EXITFUNC;
}

//for sub operation, real and complex image is the same
void EMData::sub(const EMData & em)
{
	ENTERFUNC;

	if (nx != em.get_xsize() || ny != em.get_ysize() || nz != em.get_zsize()) {
		throw ImageFormatException("images not same sizes");
	}
	else if( (is_real()^em.is_real()) == true )
	{
		throw ImageFormatException( "not support sub between real image and complex image");
	}
	else {
		flags |= EMDATA_NEEDUPD;
		const float *src_data = em.get_data();
		size_t size = nx * ny * nz;

		for (size_t i = 0; i < size; i++) {
			rdata[i] -= src_data[i];
		}
	}
	EXITFUNC;
}



void EMData::mult(float f)
{
	ENTERFUNC;

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
	EXITFUNC;
}

void EMData::mult(const EMData & em)
{
	ENTERFUNC;

	if (nx != em.get_xsize() || ny != em.get_ysize() || nz != em.get_zsize()) {
		throw ImageFormatException( "images not same sizes");
	}
	else if( (is_real()^em.is_real()) == true )
	{
		throw ImageFormatException( "not support multiply between real image and complex image");
	}
	else
	{
		flags |= EMDATA_NEEDUPD;
		const float *src_data = em.get_data();
		size_t size = nx * ny * nz;
		if( is_real() )
		{
			for (size_t i = 0; i < size; i++) {
				rdata[i] *= src_data[i];
			}
		}
		else
		{
			typedef complex<float> comp;
			for( size_t i = 0; i < size; i+=2 )
			{
				comp c_src( src_data[i], src_data[i+1] );
				comp c_rdat( rdata[i], rdata[i+1] );
				comp c_result = c_src * c_rdat;
				rdata[i] = c_result.real();
				rdata[i+1] = c_result.imag();
			}
		}
	}

	EXITFUNC;
}

void EMData::div(float f)
{
	ENTERFUNC;

	if (is_complex()) {
		ap2ri();
	}
	if (f != 0) {
		flags |= EMDATA_NEEDUPD;
		size_t size = nx * ny * nz;
		for (size_t i = 0; i < size; i++) {
			rdata[i] /= f;
		}
	}
	EXITFUNC;
}

void EMData::div(const EMData & em)
{
	ENTERFUNC;

	if (nx != em.get_xsize() || ny != em.get_ysize() || nz != em.get_zsize()) {
		throw ImageFormatException( "images not same sizes");
	}
	else if( (is_real()^em.is_real()) == true )
	{
		throw ImageFormatException( "not support division between real image and complex image");
	}
	else {
		flags |= EMDATA_NEEDUPD;
		const float *src_data = em.get_data();
		size_t size = nx * ny * nz;

		if( is_real() )
		{
			for (size_t i = 0; i < size; i++) {
				rdata[i] /= src_data[i];
			}
		}
		else
		{
			typedef complex<float> comp;
			for( size_t i = 0; i < size; i+=2 )
			{
				comp c_src( src_data[i], src_data[i+1] );
				comp c_rdat( rdata[i], rdata[i+1] );
				comp c_result = c_rdat / c_src;
				rdata[i] = c_result.real();
				rdata[i+1] = c_result.imag();
			}
		}
	}

	EXITFUNC;
}

void EMData::cconj() {
	ENTERFUNC;
	if (!is_complex() || !is_ri())
		throw ImageFormatException("EMData::conj requires a complex, ri image");
	int nxreal = nx -2 + int(is_fftodd());
	int nxhalf = nxreal/2;
	for (int iz = 0; iz < nz; iz++)
		for (int iy = 0; iy < ny; iy++)
			for (int ix = 0; ix <= nxhalf; ix++)
				cmplx(ix,iy,iz) = conj(cmplx(ix,iy,iz));
	EXITFUNC;
}


void EMData::update_stat()
{
	ENTERFUNC;

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

	int n_nonzero = 0;

	for (int i = 0; i < nx*ny*nz; i += step) {
		float v = rdata[i];
	#ifdef _WIN32
		max = _cpp_max(max,v);
		min = _cpp_min(min,v);
	#else
		max=std::max<float>(max,v);
		min=std::min<float>(min,v);
	#endif	//_WIN32
		sum += v;
		square_sum += v * (double)(v);
		if (v != 0) n_nonzero++;
	}

	int n = nx * ny * nz / step;
	double mean = sum / n;

#ifdef _WIN32
	float sigma = (float)sqrt( _cpp_max(0.0,(square_sum - sum*sum / n)/(n-1)));
	n_nonzero = _cpp_max(1,n_nonzero);
	double sigma_nonzero = sqrt( _cpp_max(0,(square_sum  - sum*sum/n_nonzero)/(n_nonzero-1)));
#else
	float sigma = (float)sqrt(std::max<double>(0.0,(square_sum - sum*sum / n)/(n-1)));
	n_nonzero = std::max<int>(1,n_nonzero);
	double sigma_nonzero = sqrt(std::max<double>(0,(square_sum  - sum*sum/n_nonzero)/(n_nonzero-1)));
#endif	//_WIN32
	double mean_nonzero = sum / n_nonzero; // previous version overcounted! G2
	
	attr_dict["minimum"] = min;
	attr_dict["maximum"] = max;
	attr_dict["mean"] = (float)(mean);
	attr_dict["sigma"] = (float)(sigma);
	attr_dict["square_sum"] = (float)(square_sum);
	attr_dict["mean_nonzero"] = (float)(mean_nonzero);
	attr_dict["sigma_nonzero"] = (float)(sigma_nonzero);
	attr_dict["is_complex"] = (int) is_complex();
	attr_dict["is_ri"] = (int) is_ri();

	flags &= ~EMDATA_NEEDUPD;
	EXITFUNC;
}

MArray2D EMData::get_2dview() const
{
	const int ndims = 2;
	if (get_ndim() != ndims) {
		throw ImageDimensionException("2D only");
	}
	boost::array<std::size_t,ndims> dims = {{nx, ny}};
	MArray2D marray(rdata, dims, boost::fortran_storage_order());
	return marray;
}

MArray3D EMData::get_3dview() const
{
	const int ndims = 3;
	boost::array<std::size_t,ndims> dims = {{nx, ny, nz}};
	MArray3D marray(rdata, dims, boost::fortran_storage_order());
	return marray;
}

MCArray2D EMData::get_2dcview() const
{
	const int ndims = 2;
	if (get_ndim() != ndims) {
		throw ImageDimensionException("2D only");
	}
	boost::array<std::size_t,ndims> dims = {{nx/2, ny}};
	complex<float>* cdata = reinterpret_cast<complex<float>*>(rdata);
	MCArray2D marray(cdata, dims, boost::fortran_storage_order());
	return marray;
}

MCArray3D EMData::get_3dcview() const
{
	const int ndims = 3;
	boost::array<std::size_t,ndims> dims = {{nx/2, ny, nz}};
	complex<float>* cdata = reinterpret_cast<complex<float>*>(rdata);
	MCArray3D marray(cdata, dims, boost::fortran_storage_order());
	return marray;
}

MCArray3D* EMData::get_3dcviewptr() const
{
	const int ndims = 3;
	boost::array<std::size_t,ndims> dims = {{nx/2, ny, nz}};
	complex<float>* cdata = reinterpret_cast<complex<float>*>(rdata);
	MCArray3D* marray = new MCArray3D(cdata, dims,
									  boost::fortran_storage_order());
	return marray;
}

MArray2D EMData::get_2dview(int x0, int y0) const
{
	const int ndims = 2;
	if (get_ndim() != ndims) {
		throw ImageDimensionException("2D only");
	}
	boost::array<std::size_t,ndims> dims = {{nx, ny}};
	MArray2D marray(rdata, dims, boost::fortran_storage_order());
	boost::array<std::size_t,ndims> bases={{x0, y0}};
	marray.reindex(bases);
	return marray;
}

MArray3D EMData::get_3dview(int x0, int y0, int z0) const
{
	const int ndims = 3;
	boost::array<std::size_t,ndims> dims = {{nx, ny, nz}};
	MArray3D marray(rdata, dims, boost::fortran_storage_order());
	boost::array<std::size_t,ndims> bases={{x0, y0, z0}};
	marray.reindex(bases);
	return marray;
}

MCArray2D EMData::get_2dcview(int x0, int y0) const
{
	const int ndims = 2;
	if (get_ndim() != ndims) {
		throw ImageDimensionException("2D only");
	}
	boost::array<std::size_t,ndims> dims = {{nx/2, ny}};
	complex<float>* cdata = reinterpret_cast<complex<float>*>(rdata);
	MCArray2D marray(cdata, dims, boost::fortran_storage_order());
	boost::array<std::size_t,ndims> bases={{x0, y0}};
	marray.reindex(bases);
	return marray;
}

MCArray3D EMData::get_3dcview(int x0, int y0, int z0) const
{
	const int ndims = 3;
	boost::array<std::size_t,ndims> dims = {{nx/2, ny, nz}};
	complex<float>* cdata = reinterpret_cast<complex<float>*>(rdata);
	MCArray3D marray(cdata, dims, boost::fortran_storage_order());
	boost::array<std::size_t,ndims> bases={{x0, y0, z0}};
	marray.reindex(bases);
	return marray;
}


EMData *EMData::get_row(int row_index) const
{
	ENTERFUNC;

	if (get_ndim() > 2) {
		throw ImageDimensionException("1D/2D image only");
	}

	EMData *ret = new EMData();
	ret->set_size(nx, 1, 1);
	memcpy(ret->get_data(), get_data() + nx * row_index, nx * sizeof(float));
	ret->done_data();
	EXITFUNC;
	return ret;
}


void EMData::set_row(const EMData * d, int row_index)
{
	ENTERFUNC;

	if (get_ndim() > 2) {
		throw ImageDimensionException("1D/2D image only");
	}
	if (d->get_ndim() != 1) {
		throw ImageDimensionException("1D image only");
	}

	float *dst = get_data();
	float *src = d->get_data();
	memcpy(dst + nx * row_index, src, nx * sizeof(float));
	done_data();
	EXITFUNC;
}

EMData *EMData::get_col(int col_index) const
{
	ENTERFUNC;

	if (get_ndim() != 2) {
		throw ImageDimensionException("2D image only");
	}

	EMData *ret = new EMData();
	ret->set_size(ny, 1, 1);
	float *dst = ret->get_data();
	float *src = get_data();

	for (int i = 0; i < ny; i++) {
		dst[i] = src[i * nx + col_index];
	}

	ret->done_data();
	EXITFUNC;
	return ret;
}


void EMData::set_col(const EMData * d, int n)
{
	ENTERFUNC;

	if (get_ndim() != 2) {
		throw ImageDimensionException("2D image only");
	}
	if (d->get_ndim() != 1) {
		throw ImageDimensionException("1D image only");
	}

	float *dst = get_data();
	float *src = d->get_data();

	for (int i = 0; i < ny; i++) {
		dst[i * nx + n] = src[i];
	}

	done_data();
	EXITFUNC;
}


EMObject EMData::get_attr(const string & key)
{
	ENTERFUNC;

	update_stat();

	float mean = attr_dict["mean"];
	float sigma = attr_dict["sigma"];
	size_t size = nx * ny * nz;

	if (key == "kurtosis") {
		double kurtosis_sum = 0;

		for (size_t k = 0; k < size; k++) {
			float t = (rdata[k] - mean) / sigma;
			float tt = t * t;
			kurtosis_sum += tt * tt;
		}

		float kurtosis = (float)(kurtosis_sum / size - 3.0);
		attr_dict["kurtosis"] = kurtosis;
		return attr_dict["kurtosis"];
	}
	else if (key == "skewness") {
		double skewness_sum = 0;
		for (size_t k = 0; k < size; k++) {
			float t = (rdata[k] - mean) / sigma;
			skewness_sum +=  t * t * t;
		}
		float skewness = (float)(skewness_sum / size);
		attr_dict["skewness"] = skewness;
		return attr_dict["skewness"];
	}


	EXITFUNC;
	return attr_dict[key];
}

void EMData::set_size(int x, int y, int z)
{
	ENTERFUNC;

	if (x <= 0) {
		throw InvalidValueException(x, "x size <= 0");
	}
	else if (y <= 0) {
		throw InvalidValueException(y, "y size <= 0");
	}
	else if (z <= 0) {
		throw InvalidValueException(z, "z size <= 0");
	}

	int old_nx = nx;
	nx = x;
	ny = y;
	nz = z;

	size_t size = (size_t)x * (size_t)y * (size_t)z * sizeof(float);
	rdata = static_cast < float *>(realloc(rdata, size));
	update();

	attr_dict["nx"] = x;
	attr_dict["ny"] = y;
	attr_dict["nz"] = z;

	if (old_nx == 0) {
		memset(rdata, 0, size);
	}

	if (supp) {
		free(supp);
		supp = 0;
	}
	EXITFUNC;
}

float *EMData::get_data() const
{
	//flags |= EMDATA_BUSY;
	return rdata;
}

void EMData::done_data()
{
	flags &= (~EMDATA_BUSY);
	flags |= EMDATA_NEEDUPD;
}

Dict EMData::get_attr_dict()
{
	update_stat();
	return Dict(attr_dict);
}

void EMData::set_attr_dict(const Dict & new_dict)
{
	attr_dict = new_dict;
}


void EMData::set_ctf(Ctf * new_ctf)
{
	ENTERFUNC;

	if (!ctf) {
		ctf = new SimpleCtf();
	}

	ctf->copy_from(new_ctf);
	EXITFUNC;
}


vector < EMData * >EMData::read_images(const string & filename, vector < int >img_indices,
									   bool header_only)
{
	ENTERFUNC;

	int total_img = EMUtil::get_image_count(filename);
	size_t num_img = img_indices.size();

	for (size_t i = 0; i < num_img; i++) {
		if (img_indices[i] < 0 && img_indices[i] >= total_img) {
			throw OutofRangeException(0, total_img, img_indices[i], "image index");
		}
	}

	size_t n = num_img == 0 ? total_img : num_img;

	vector < EMData * >v;
	for (size_t j = 0; j < n; j++) {
		EMData *d = new EMData();
		size_t k = num_img == 0 ? j : img_indices[j];
		try {
			d->read_image(filename, (int)k, header_only);
		}
		catch(E2Exception &e) {
			if( d )
			{
				delete d;
				d = 0;
			}
			throw(e);
		}

		v.push_back(d);
	}
	EXITFUNC;
	return v;
}


vector < EMData * >EMData::read_images_ext(const string & filename, int img_index_start,
										   int img_index_end, bool header_only,
										   const string & ext)
{
	ENTERFUNC;

	if (img_index_end < img_index_start) {
		throw InvalidValueException(img_index_end, "image index end < image index start");
	}
	string new_filename = filename;
	new_filename = new_filename.insert(new_filename.rfind("."), ext);
	int num_img = EMUtil::get_image_count(new_filename);

	if (img_index_start < 0 || img_index_start >= num_img) {
		throw OutofRangeException(0, num_img-1, img_index_start, "image index start");
	}

	if (img_index_end >= num_img) {
		img_index_end = num_img - 1;
	}

	vector < EMData * >v;

	for (int i = img_index_start; i < img_index_end; i++) {
		EMData *d = new EMData();
		try {
			d->read_image(new_filename, i, header_only);
		}
		catch(E2Exception &e) {
			if( d )
			{
				delete d;
				d = 0;
			}
			throw(e);
		}
		v.push_back(d);
	}
	EXITFUNC;
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

EMData * EMData::power(int n)
{
	if( n<0 ) {
		throw InvalidValueException(n, "the power of negative integer not supported.");
	}
	
	EMData * r = this->copy();
	if( n == 0 ) {
		r->to_one();
	}
	else if( n>1 ) {
		for( int i=1; i<n; i++ ) {
			*r *= *this;
		}
	}

	r->update_stat();
	return r;
}

EMData * EMAN::operator+(const EMData & em, float n)
{
	EMData * r = em.copy();
	r->add(n);
	return r;
}

EMData * EMAN::operator-(const EMData & em, float n)
{
	EMData* r = em.copy();
	r->sub(n);
	return r;
}

EMData * EMAN::operator*(const EMData & em, float n)
{
	EMData* r = em.copy();
	r ->mult(n);
	return r;
}

EMData * EMAN::operator/(const EMData & em, float n)
{
	EMData * r = em.copy();
	r->div(n);
	return r;
}


EMData * EMAN::operator+(float n, const EMData & em)
{
	EMData * r = em.copy();
	r->add(n);
	return r;
}

EMData * EMAN::operator-(float n, const EMData & em)
{
	EMData * r = em.copy();
	r->sub(n);
	return r;
}

EMData * EMAN::operator*(float n, const EMData & em)
{
	EMData * r = em.copy();
	r->mult(n);
	return r;
}

EMData * EMAN::operator/(float n, const EMData & em)
{
	EMData * r = em.copy();
	r->div(n);
	return r;
}


EMData * EMAN::operator+(const EMData & a, const EMData & b)
{
	EMData * r = a.copy();
	r->add(b);
	return r;
}

EMData * EMAN::operator-(const EMData & a, const EMData & b)
{
	EMData * r = a.copy();
	r->sub(b);
	return r;
}

EMData * EMAN::operator*(const EMData & a, const EMData & b)
{
	EMData * r = a.copy();
	r->mult(b);
	return r;
}

EMData * EMAN::operator/(const EMData & a, const EMData & b)
{
	EMData * r = a.copy();
	r->div(b);
	return r;
}

double EMData::dot_rotate_translate(EMData * with, float dx, float dy, float da)
{
	ENTERFUNC;

	if (!EMUtil::is_same_size(this, with)) {
		LOGERR("images not same size");
		throw ImageFormatException("images not same size");
	}

	if (get_ndim() == 3) {
		LOGERR("1D/2D Images only");
		throw ImageDimensionException("1D/2D only");
	}

	float *this_data = 0;

	this_data = get_data();

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

	EXITFUNC;
	return result;
}

void EMData::rotate_x(int dx)
{
	ENTERFUNC;

	if (get_ndim() > 2) {
		throw ImageDimensionException("no 3D image");
	}

	float *tmp = new float[nx];
	size_t row_size = nx * sizeof(float);

	for (int y = 0; y < ny; y++) {
		int y_nx = y * nx;
		for (int x = 0; x < nx; x++) {
			tmp[x] = rdata[y_nx + (x + dx) % nx];
		}
		memcpy(&rdata[y_nx], tmp, row_size);
	}

	done_data();
	if( tmp )
	{
		delete[]tmp;
		tmp = 0;
	}
	EXITFUNC;
}

void EMData::set_xyz_origin(float origin_x, float origin_y, float origin_z)
{
	attr_dict["origin_row"] = origin_x;
	attr_dict["origin_col"] = origin_y;
	attr_dict["origin_sec"] = origin_z;
}

float *EMData::setup4slice(bool redo)
{
	ENTERFUNC;

	if (!is_complex()) {
		throw ImageFormatException("complex image only");
	}

	if (get_ndim() != 3) {
		throw ImageDimensionException("3D only");
	}

	if (supp) {
		if (redo) {
			free(supp);
			supp = 0;
		}
		else {
			EXITFUNC;
			return supp;
		}
	}

	const int SUPP_ROW_SIZE = 8;
	const int SUPP_ROW_OFFSET = 4;
	const int supp_size = SUPP_ROW_SIZE + SUPP_ROW_OFFSET;

	supp = (float *) calloc(supp_size * ny * nz, sizeof(float));
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

	EXITFUNC;
	return supp;
}

void EMData::to_zero()
{
	ENTERFUNC;

	memset(rdata, 0, nx * ny * nz * sizeof(float));
	done_data();
	EXITFUNC;
}

void EMData::scale(float s)
{
	ENTERFUNC;
	Transform3D t;
	t.set_scale(s);
	rotate_translate(t);
	EXITFUNC;
}

void EMData::translate(int dx, int dy, int dz)
{
	ENTERFUNC;
	translate(Vec3i(dx, dy, dz));
	EXITFUNC;
}

void EMData::translate(float dx, float dy, float dz)
{
	ENTERFUNC;
	int dx_ = Util::round(dx);
	int dy_ = Util::round(dy);
	int dz_ = Util::round(dz);
	if( ( (dx-dx_) == 0 ) && ( (dy-dy_) == 0 ) && ( (dz-dz_) == 0 )) {
		translate(dx_, dy_, dz_);
	}
	else {
		translate(Vec3f(dx, dy, dz));
	}
	EXITFUNC;
}

void EMData::translate(const Vec3i &translation)
{
	ENTERFUNC;

	//if traslation is 0, do nothing
	if( translation[0] == 0 && translation[1] == 0 && translation[2] == 0) {
		EXITFUNC;
		return;
	}

	float *this_data = get_data();
	int data_size = sizeof(float)*get_xsize()*get_ysize()*get_zsize();
	float *tmp_data = (float *)malloc(data_size);
	memcpy(tmp_data, this_data, data_size);

	int x0, x1, x2;
	if( translation[0] < 0 ) {
		x0 = 0;
		x1 = nx;
		x2 = 1;
	}
	else {
		x0 = nx-1;
		x1 = -1;
		x2 = -1;
	}

	int y0, y1, y2;
	if( translation[1] < 0 ) {
		y0 = 0;
		y1 = ny;
		y2 = 1;
	}
	else {
		y0 = ny-1;
		y1 = -1;
		y2 = -1;
	}

	int z0, z1, z2;
	if( translation[2] < 0 ) {
		z0 = 0;
		z1 = nz;
		z2 = 1;
	}
	else {
		z0 = nz-1;
		z1 = -1;
		z2 = -1;
	}

	int xp, yp, zp;
	int tx = translation[0];
	int ty = translation[1];
	int tz = translation[2];
	for (int y = y0; y != y1; y += y2) {
		for (int x = x0; x != x1; x += x2) {
			for (int z = z0; z != z1; z+=z2) {
				xp = x - tx;
				yp = y - ty;
				zp = z - tz;
				if (xp < 0 || yp < 0 || zp<0 || xp >= nx || yp >= ny || zp >= nz) {
					this_data[x + y * nx + z * nx * ny] = 0;
				}
				else {
					this_data[x + y * nx + z * nx * ny] = tmp_data[xp + yp * nx + zp * nx * ny];
				}
			}
		}
	}

	if( tmp_data ) {
		free(tmp_data);
		tmp_data = 0;
	}

	done_data();
	all_translation += translation;

	EXITFUNC;
}

void EMData::translate(const Vec3f &translation)
{
	ENTERFUNC;

	//if traslation is 0, do nothing
	if( translation[0] == 0.0f && translation[1] == 0.0f && translation[2] == 0.0f ) {
		EXITFUNC;
		return;
	}

	float *this_data = get_data();
	EMData *tmp_emdata = copy();

	int x0, x1, x2;
	if( translation[0] < 0 ) {
		x0 = 0;
		x1 = nx;
		x2 = 1;
	}
	else {
		x0 = nx-1;
		x1 = -1;
		x2 = -1;
	}

	int y0, y1, y2;
	if( translation[1] < 0 ) {
		y0 = 0;
		y1 = ny;
		y2 = 1;
	}
	else {
		y0 = ny-1;
		y1 = -1;
		y2 = -1;
	}

	int z0, z1, z2;
	if( translation[2] < 0 ) {
		z0 = 0;
		z1 = nz;
		z2 = 1;
	}
	else {
		z0 = nz-1;
		z1 = -1;
		z2 = -1;
	}

	if( nz == 1 ) 	//2D translation
	{
		int tx = Util::round(translation[0]);
		int ty = Util::round(translation[1]);
		float ftx = translation[0];
		float fty = translation[1];
		int xp, yp;
		for (int y = y0; y != y1; y += y2) {
			for (int x = x0; x != x1; x += x2) {
				xp = x - tx;
				yp = y - ty;

				if (xp < 0 || yp < 0 || xp >= nx || yp >= ny) {
					this_data[x + y * nx] = 0;
				}
				else {
					float fxp = static_cast<float>(x) - ftx;
					float fyp = static_cast<float>(y) - fty;
					this_data[x + y * nx] = tmp_emdata->sget_value_at_interp(fxp, fyp);
				}
			}
		}
	}
	else 	//3D translation
	{
		int tx = Util::round(translation[0]);
		int ty = Util::round(translation[1]);
		int tz = Util::round(translation[2]);
		float ftx = translation[0];
		float fty = translation[1];
		float ftz = translation[2];
		int xp, yp, zp;
		for (int z = z0; z != z1; z += z2) {
			for (int y = y0; y != y1; y += y2) {
				for (int x = x0; x != x1; x += x2) {
					xp = x - tx;
					yp = y - ty;
					zp = z - tz;
					if (xp < 0 || yp < 0 || zp<0 || xp >= nx || yp >= ny || zp >= nz) {
						this_data[x + y * nx] = 0;
					}
					else {
						float fxp = static_cast<float>(x) - ftx;
						float fyp = static_cast<float>(y) - fty;
						float fzp = static_cast<float>(z) - ftz;
						this_data[x + y * nx + z * nx * ny] = tmp_emdata->sget_value_at_interp(fxp, fyp, fzp);
					}
				}
			}
		}

	}

	if( tmp_emdata ) {
		delete tmp_emdata;
		tmp_emdata = 0;
	}
	done_data();
	update();
	all_translation += translation;
	EXITFUNC;
}

void EMData::rotate(float az, float alt, float phi)
{
	Transform3D t(az, alt, phi);
	rotate_translate(t);
}

void EMData::rotate(const Transform3D & t)
{
	rotate_translate(t);
}

void EMData::rotate_translate(float az, float alt, float phi, float dx, float dy, float dz)
{
	Transform3D t(Vec3f(dx, dy, dz),  az, alt, phi);
	rotate_translate(t);
}


void EMData::rotate_translate(float az, float alt, float phi, float dx, float dy,
							  float dz, float pdx, float pdy, float pdz)
{
	Transform3D t(Vec3f(dx, dy, dz), Vec3f(pdx,pdy,pdz), az, alt, phi);
	rotate_translate(t);
}



void EMData::rotate_translate(const Transform3D & xform)
{
	ENTERFUNC;

	float scale = xform.get_scale();
	Vec3f dcenter = xform.get_center();
	Vec3f translation = xform.get_posttrans();
	Dict rotation = xform.get_rotation(Transform3D::EMAN);

	int nx2 = nx;
	int ny2 = ny;
	float inv_scale = 1.0f;

	if (scale != 0) {
		inv_scale = 1.0f / scale;
	}

	float *src_data = 0;
	float *des_data = 0;

	src_data = get_data();
	des_data = (float *) malloc(nx * ny * nz * sizeof(float));

	if (nz == 1) {
		float mx0 = inv_scale * cos((float)rotation["phi"]);
		float mx1 = inv_scale * sin((float)rotation["phi"]);

		float x2c = nx / 2.0f - dcenter[0] - translation[0];
		float y2c = ny / 2.0f - dcenter[1] - translation[1];
		float y = -ny / 2.0f + dcenter[0];

		for (int j = 0; j < ny; j++, y += 1.0f) {
			float x = -nx / 2.0f + dcenter[1];

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
		Transform3D mx = xform;
		mx.set_scale(inv_scale);
		int nxy = nx * ny;
		int l = 0;

		for (int k = 0; k < nz; k++) {
			for (int j = -ny / 2; j < ny - ny / 2; j++) {
				for (int i = -nx / 2; i < nx - nx / 2; i++, l++) {
					float x2 = mx[0][0] * i + mx[0][1] * j + mx[0][2] * k + nx / 2;
					float y2 = mx[1][0] * i + mx[1][1] * j + mx[1][2] * k + ny / 2;
					float z2 = mx[2][0] * i + mx[2][1] * j + mx[2][2] * k + 0 / 2;

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
		Transform3D mx = xform;
		mx.set_scale(inv_scale);

		Vec3f dcenter2 = Vec3f((float)nx,(float)ny,(float)nz)/(-2.0f) + dcenter;
		Vec3f v4 = dcenter2 * mx  - dcenter2 - translation; // verify this

		int nxy = nx * ny;
		int l = 0;

		for (int k = 0; k < nz; k++) {
			Vec3f v3 = v4;

			for (int j = 0; j < ny; j++) {
				Vec3f v2 = v3;

				for (int i = 0; i < nx; i++, l++) {

					if (v2[0] < 0 || v2[1] < 0 || v2[2] < 0 ||
						v2[0] >= nx - 1 || v2[1] >= ny - 1 || v2[2] >= nz - 1) {
						des_data[l] = 0;
					}
					else {
						int x = Util::fast_floor(v2[0]);
						int y = Util::fast_floor(v2[1]);
						int z = Util::fast_floor(v2[2]);
						Vec3f tuv = v2 - Vec3f((float)x,(float)y,(float)z);
						int ii = x + y * nx + z * nxy;

						des_data[l] = Util::trilinear_interpolate(src_data[ii],
																  src_data[ii + 1],
																  src_data[ii + nx],
																  src_data[ii + nx + 1],
																  src_data[ii + nx * ny],
																  src_data[ii + nxy + 1],
																  src_data[ii + nxy + nx],
																  src_data[ii + nxy + nx + 1],
																  tuv[0],
																  tuv[1],
																  tuv[2]
																  );

					}

					v2 += mx.get_matrix3_col(0);
				}
				v3 += mx.get_matrix3_col(1);
			}
			v4 += mx.get_matrix3_col(2); //  or should it be row?   PRB April 2005
		}

	}

	if( rdata )
	{
		free(rdata);
		rdata = 0;
	}
	rdata = des_data;

	scale_pixel(inv_scale);

	attr_dict["origin_row"] = (float) attr_dict["origin_row"] * inv_scale;
	attr_dict["origin_col"] = (float) attr_dict["origin_col"] * inv_scale;
	attr_dict["origin_sec"] = (float) attr_dict["origin_sec"] * inv_scale;

	done_data();
	update();


	all_translation += translation;
	EXITFUNC;
}


void EMData::rotate_180()
{
	ENTERFUNC;

	if (nx != ny) {
		throw ImageFormatException("non-square image");
	}

	if (get_ndim() != 2) {
		throw ImageDimensionException("2D only");
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
	EXITFUNC;
}

EMData *EMData::do_radon()
{
	ENTERFUNC;

	if (get_ndim() != 2) {
		throw ImageDimensionException("2D only");
	}

	if (nx != ny) {
		throw ImageFormatException("square image only");
	}

	EMData *result = new EMData();
	result->set_size(nx, ny, 1);
	result->to_zero();
	float *result_data = result->get_data();

	EMData *this_copy = this;
	this_copy = copy();

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

	if( this_copy )
	{
		delete this_copy;
		this_copy = 0;
	}

	EXITFUNC;
	return result;
}


void EMData::mean_shrink(float shrink_factor0)
{
	ENTERFUNC;
	int shrink_factor = int(shrink_factor0);

	if (shrink_factor0 <= 1.0F || ((shrink_factor0 != shrink_factor) && (shrink_factor0 != 1.5F) ) ) {
		throw InvalidValueException(shrink_factor0,
									"mean shrink: shrink factor must be >1 integer or 1.5");
	}

/*	if ((nx % shrink_factor != 0) || (ny % shrink_factor != 0) ||
		(nz > 1 && (nz % shrink_factor != 0))) {
		throw InvalidValueException(shrink_factor,
									"Image size not divisible by shrink factor");
	}*/

	// here handle the special averaging by 1.5 for 2D case
	if (shrink_factor0==1.5 ) {
		if (nz > 1 ) throw InvalidValueException(shrink_factor0, "mean shrink: only support 2D images for shrink factor = 1.5");

		int shrinked_nx = (int(nx / 1.5)+1)/2*2;	// make sure the output size is even
		int shrinked_ny = (int(ny / 1.5)+1)/2*2;
		int nx0 = nx, ny0 = ny;	// the original size

		EMData* orig = copy();
		set_size(shrinked_nx, shrinked_ny, 1);	// now nx = shrinked_nx, ny = shrinked_ny
		to_zero();

		float *data = get_data(), *data0 = orig->get_data();

		for (int j = 0; j < ny; j++) {
			int jj = int(j * 1.5);
			float jw0 = 1.0F, jw1 = 0.5F;	// 3x3 -> 2x2, so each new pixel should have 2.25 of the old pixels
			if ( j%2 ) {
				jw0 = 0.5F;
				jw1 = 1.0F;
			}
			for (int i = 0; i < nx; i++) {
				int ii = int(i * 1.5);
				float iw0 = 1.0F, iw1 = 0.5F;
				float w = 0.0F;

				if ( i%2 ) {
					iw0 = 0.5F;
					iw1 = 1.0F;
				}
				if ( jj < ny0 ) {
					if ( ii < nx0 ) {
						data[j * nx + i] = data0[ jj * nx0 + ii ] * jw0 * iw0 ;
						w += jw0 * iw0 ;
						if ( ii+1 < nx0 ) {
							data[j * nx + i] += data0[ jj * nx0 + ii + 1] * jw0 * iw1;
							w += jw0 * iw1;
						}
					}
					if ( jj +1 < ny0 ) {
						if ( ii < nx0 ) {
							data[j * nx + i] += data0[ (jj+1) * nx0 + ii ] * jw1 * iw0;
							w += jw1 * iw0;
							if ( ii+1 < nx0 ) {
								data[j * nx + i] += data0[ (jj+1) * nx0 + ii + 1] * jw1 * iw1;
								w += jw1 * iw1;
							}
						}
					}
				}
				if ( w>0 ) data[j * nx + i] /= w;
			}
		}
		orig->done_data();
		if( orig )
		{
			delete orig;
			orig = 0;
		}
		done_data();
		update();

		return;
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
	scale_pixel((float)shrink_factor);
	EXITFUNC;
}

void EMData::median_shrink(int shrink_factor)
{
	ENTERFUNC;

	if (shrink_factor <= 1) {
		throw InvalidValueException(shrink_factor,
									"median shrink: shrink factor must > 1");
	}

	if ((nx % shrink_factor != 0) || (ny % shrink_factor != 0) ||
		(nz > 1 && (nz % shrink_factor != 0))) {
		throw InvalidValueException(shrink_factor,
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
	scale_pixel((float)shrink_factor);

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

	if( data_copy )
	{
		delete[]data_copy;
		data_copy = 0;
	}

	if( mbuf )
	{
		delete[]mbuf;
		mbuf = 0;
	}
	EXITFUNC;
}

// NOTE : x axis is from 0 to 0.5  (Nyquist), and thus properly handles non-square images
// complex only
void EMData::apply_radial_func(float x0, float step, vector < float >array, bool interp)
{
	ENTERFUNC;

	if (!is_complex()) throw ImageFormatException("apply_radial_func requires a complex image");

	int n = static_cast < int >(array.size());

//	printf("%f %f %f\n",array[0],array[25],array[50]);

	ap2ri();

	size_t ndims = get_ndim();

	if (ndims == 2) {
		int k = 0;
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i += 2, k += 2) {
				float r;
				if (j<ny/2) r = (float)hypot(i/(float)(nx*2), j/(float)ny);
				else r = (float)hypot(i/(float)(nx*2), (ny-j)/(float)ny);
				r = (r - x0) / step;

				int l = 0;
				if (interp) {
					l = (int) floor(r);
				}
				else {
					l = (int) floor(r + 1);
				}


				float f = 0;
				if (l >= n - 2) {
					f = array[n - 1];
				}
				else {
					if (interp) {
						r -= l;
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
			float mnz;
			if (m<nz/2) mnz=m*m/(float)(nz*nz);
			else { mnz=(nz-m)/(float)nz; mnz*=mnz; }

			for (int j = 0; j < ny; j++) {
				float jny;
				if (j<ny/2) jny= j*j/(float)(ny*ny);
				else { jny=(ny-j)/(float)ny; jny*=jny; }

				for (int i = 0; i < nx; i += 2, k += 2) {
					float r = sqrt((i * i / (nx*nx*4.0f)) + jny + mnz);
					r = (r - x0) / step;

					int l = 0;
					if (interp) {
						l = (int) floor(r);
					}
					else {
						l = (int) floor(r + 1);
					}


					float f = 0;
					if (l >= n - 2) {
						f = array[n - 1];
					}
					else {
						if (interp) {
							r -= l;
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
	EXITFUNC;
}

float EMData::calc_center_density()
{
	ENTERFUNC;

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
		float mean = (float)(sum / norm);
		if (fabs(mean - center) < sigma2) {
			ds *= 0.5f;
		}
		center = mean;
	}
	EXITFUNC;

	return center;
}

float EMData::calc_sigma_diff()
{
	ENTERFUNC;

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

	float sigup = sqrt((float)sum_up / nup);
	float sigdown = sqrt((float)sum_down / ndown);
	float sig_diff = fabs(sigup - sigdown) / sigma;


	EXITFUNC;
	return sig_diff;

}


IntPoint EMData::calc_min_location() const
{
	ENTERFUNC;

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
	ENTERFUNC;

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

	EXITFUNC;
	return IntPoint(max_x, max_y, max_z);
}



int EMData::calc_min_index() const
{
	IntPoint min_location = calc_min_location();
	int i = min_location[0] + min_location[1] * nx + min_location[2] * nx * ny;
	return i;
}


int EMData::calc_max_index() const
{
	IntPoint max_location = calc_max_location();
	int i = max_location[0] + max_location[1] * nx + max_location[2] * nx * ny;
	return i;
}

EMData *EMData::calc_ccfx(EMData * with, int y0, int y1, bool no_sum)
{
	ENTERFUNC;

	if (!with) {
		LOGERR("NULL 'with' image. ");
		throw NullPointerException("NULL input image");
	}

	if (!EMUtil::is_same_size(this, with)) {
		LOGERR("images not same size: (%d,%d,%d) != (%d,%d,%d)",
			   nx, ny, nz,
			   with->get_xsize(), with->get_ysize(), with->get_zsize());
		throw ImageFormatException("images not same size");
	}
	if (get_ndim() > 2) {
		LOGERR("2D images only");
		throw ImageDimensionException("2D images only");
	}

	EMData *cf = new EMData();
	if (no_sum) {
		cf->set_size(nx, y1 - y0 + 1, 1);
	}
	else {
		cf->set_size(nx, 1, 1);
	}

	cf->set_attr("label", "CCFx");
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
		float *f1 = (float *) calloc(nx, sizeof(float));
		float *f2 = (float *) calloc(nx, sizeof(float));

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

		if( f1 )
		{
			free(f1);
			f1 = 0;
		}
		if( f2 )
		{
			free(f2);
			f2 = 0;
		}
	}

	cf->done_data();
	done_data();
	with->done_data();


	EXITFUNC;
	return cf;
}

#if 0
void EMData::calc_rcf(EMData * with, vector < float >&sum_array)
{
	ENTERFUNC;

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

	vector<float> vdata, vdata2;
	for (int i = 8; i < nx2; i += 6) {
		vdata = calc_az_dist(array_size, 0, da, i, i + 6);
		vdata2 = with->calc_az_dist(array_size, 0, da, i, i + 6);
		assert(vdata.size() <= array_size + 2);
		assert(cdata2.size() <= array_size + 2);
		std::copy(vdata.begin(), vdata.end(), dat);
		std::copy(vdata2.begin(), vdata2.end(), dat2);

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

	if( dat )
	{
		delete[]dat;
		dat = 0;
	}

	if( dat2 )
	{
		delete[]dat2;
		dat2 = 0;
	}
	EXITFUNC;
}

#endif


void EMData::to_one()
{
	ENTERFUNC;

	if (is_complex()) {
		set_ri(true);
	}

	for (int i = 0; i < nx * ny * nz; i++) {
		rdata[i] = 1.0f;
	}

	update();
	EXITFUNC;
}




EMData *EMData::calc_ccf(EMData * with, fp_flag fpflag) {
	if (with==this) return self_correlation(this,fpflag);
	else {
		if (with) return correlation(this, with, fpflag);
		else return autocorrelation(this, fpflag);
	}
}


EMData *EMData::make_rotational_footprint(bool premasked, bool unwrap)
{
	ENTERFUNC;

	static EMData obj_filt;
	EMData* filt = &obj_filt;
	filt->set_complex(true);

	if (rfp) {
		return rfp;
	}

	if (nx & 1) {
		LOGERR("even image xsize only");
		throw ImageFormatException("even image xsize only");
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
		tmp2->process("eman1.mask.sharp", Dict("outer_radius", nx / 2, "value", 0));
	}

	if (filt->get_xsize() != tmp2->get_xsize() + 2 || filt->get_ysize() != tmp2->get_ysize() ||
		filt->get_zsize() != tmp2->get_zsize()) {
		filt->set_size(tmp2->get_xsize() + 2, tmp2->get_ysize(), tmp2->get_zsize());
		filt->to_one();

		filt->process("eman1.filter.highpass.gaussian", Dict("highpass", 1.5/nx));
	}

	EMData *tmp = tmp2->calc_mutual_correlation(tmp2, true, filt);
	if( tmp2 )
	{
		delete tmp2;
		tmp2 = 0;
	}

	Region r2;
	if (nz == 1) {
		r2 = Region(cs - nx / 4, cs - ny / 4, nx * 3 / 2, ny * 3 / 2);
	}
	else {
		r2 = Region(cs - nx / 4, cs - ny / 4, cs - nz / 4, nx * 3 / 2, ny * 3 / 2, nz * 3 / 2);
	}
	tmp2 = tmp->get_clip(r2);
	rfp = tmp2;

	if( tmp )
	{
		delete tmp;
		tmp = 0;
	}

	EMData * result = rfp;

	if (nz == 1) {
		if (!unwrap) {
			tmp2->process("eman1.mask.sharp", Dict("outer_radius", -1, "value", 0));
			rfp = 0;
			result = tmp2;
		}
		else {
			rfp = tmp2->unwrap();
			if( tmp2 )
			{
				delete tmp2;
				tmp2 = 0;
			}
			result = rfp;
		}
	}

	EXITFUNC;
	return result;
}


EMData *EMData::calc_mutual_correlation(EMData * with, bool tocorner, EMData * filter)
{
	ENTERFUNC;

	if (with && !EMUtil::is_same_size(this, with)) {
		LOGERR("images not same size");
		throw ImageFormatException( "images not same size");
	}

	EMData *this_fft = 0;
	this_fft = do_fft();

	if (!this_fft) {
		LOGERR("FFT returns NULL image");
		throw NullPointerException("FFT returns NULL image");
	}

	this_fft->ap2ri();
	EMData *cf = 0;

	if (with) {
		cf = with->do_fft();
		if (!cf) {
			LOGERR("FFT returns NULL image");
			throw NullPointerException("FFT returns NULL image");
		}
		cf->ap2ri();
	}
	else {
		cf = this_fft->copy();
	}

	if (filter) {
		if (!EMUtil::is_same_size(filter, cf)) {
			LOGERR("improperly sized filter");
			throw ImageFormatException("improperly sized filter");
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
		cf->process("eman1.xform.phaseorigin");
	}

	EMData *f2 = cf->do_ift();

	if( cf )
	{
		delete cf;
		cf = 0;
	}

	if( this_fft )
	{
		delete this_fft;
		this_fft = 0;
	}

	f2->set_attr("label", "MCF");
	f2->set_path("/tmp/eman.mcf");

	EXITFUNC;
	return f2;
}

EMData *EMData::unwrap(int r1, int r2, int xs, int dx, int dy, bool do360)
{
	ENTERFUNC;

	if (get_ndim() != 2) {
		throw ImageDimensionException("2D image only");
	}

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

	EXITFUNC;
	return ret;
}


void EMData::add_incoherent(EMData * obj)
{
	ENTERFUNC;

	if (!obj) {
		LOGERR("NULL image");
		throw NullPointerException("NULL image");
	}

	if (!obj->is_complex() || !is_complex()) {
		throw ImageFormatException("complex images only");
	}

	if (!EMUtil::is_same_size(this, obj)) {
		throw ImageFormatException("images not same size");
	}

	ri2ap();
	obj->ri2ap();

	float *dest = get_data();
	float *src = obj->get_data();
	int size = nx * ny * nz;
	for (int j = 0; j < size; j += 2) {
		dest[j] = (float) hypot(src[j], dest[j]);
		dest[j + 1] = 0;
	}

	obj->done_data();
	done_data();
	update();
	EXITFUNC;
}


vector < float >EMData::calc_radial_dist(int n, float x0, float dx)
{
	ENTERFUNC;

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
					r = (Util::hypot3(x / 2.0f, (float)(y<ny/2?y:ny-y), (float)(z<half_nz?z:nz-z)) - x0) / dx;
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

	if( yc )
	{
		delete[]yc;
		yc = 0;
	}

	vector < float >dv(n);
	for (int i = 0; i < n; i++) {
		dv[i] = d[i];
	}

	if( d )
	{
		delete[]d;
		d = 0;
	}

	EXITFUNC;
	return dv;
}

vector < float >EMData::calc_radial_dist(int n, float x0, float dx, float acen, float awid)
{
	ENTERFUNC;

	if (nz > 1) {
		LOGERR("2D images only.");
		throw ImageDimensionException("2D images only");
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
					a = atan2(y<ny/2.0f?y:ny-y, x / 2.0f);
				}
				else {
					a = atan2(y - ny / 2.0f, x - nx / 2.0f);
				}
			}

			if (fabs(Util::angle_sub_pi(a, acen)) <= awid) {
				float r = 0;
				if (is_complex()) {
					r = ((float)hypot(x / 2.0f, y - ny / 2.0f) - x0) / dx;
				}
				else {
					r = ((float)hypot(x - (nx - 1) / 2.0f, y - (ny - 1) / 2.0f) - x0) / dx;
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

	if( yc )
	{
		delete[]yc;
		yc = 0;
	}

	if( yc2 )
	{
		delete[]yc2;
		yc2 = 0;
	}

	if( d2 )
	{
		delete[]d2;
		d2 = 0;
	}

	if( d )
	{
		delete[]d;
		d = 0;
	}

	EXITFUNC;
	return dv;
}

EMData* EMData::rotavg()
{
	ENTERFUNC;

	if (nz > 1) {
		LOGERR("2D images only.");
		throw ImageDimensionException("2D images only");
	}
	MArray2D arr = get_2dview(-nx/2, -ny/2);
#ifdef _WIN32
	int rmax = _MIN(nx/2 + nx%2, ny/2 + ny%2);
#else
	int rmax = std::min(nx/2 + nx%2, ny/2 + ny%2);
#endif	//_WIN32
	EMData* ret = new EMData();
	ret->set_size(rmax+1, 1, 1);
	ret->to_zero();
	float* retarr = ret->get_data();
	vector<float> count(rmax+1);
	for (int j = -ny/2; j < ny/2 + ny%2; j++) {
		if (abs(j) > rmax) continue;
		for (int i = -nx/2; i < nx/2 + nx%2; i++) {
			float r = sqrt(float(j*j) + float(i*i));
			int ir = int(r);
			if (ir >= rmax) continue;
			float frac = r - float(ir);
			retarr[ir] += arr[i][j]*(1.0f - frac);
			retarr[ir+1] += arr[i][j]*frac;
			count[ir] += 1.0f - frac;
			count[ir+1] += frac;
		}
	}
	for (int ir = 0; ir <= rmax; ir++) {
	#ifdef _WIN32
		retarr[ir] /= _MAX(count[ir],1.0f);
	#else
		retarr[ir] /= std::max(count[ir],1.0f);
	#endif	//_WIN32
	}

	ret->update();
	ret->done_data();
	EXITFUNC;
	return ret;
}

float EMData::calc_dist(EMData * second_img, int y_index) const
{
	ENTERFUNC;

	if (get_ndim() != 1) {
		throw ImageDimensionException("'this' image is 1D only");
	}

	if (second_img->get_xsize() != nx || ny != 1) {
		throw ImageFormatException("image xsize not same");
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
	EXITFUNC;
	return sqrt(ret);
}



EMData *EMData::calc_flcf(EMData * with, int radius, const string & mask_filter)
{
	ENTERFUNC;

	if (!with) {
		LOGERR("input image is NULL");
		throw NullPointerException("input image is NULL");
	}

	Dict filter_dict;
	if (mask_filter == "eman1.mask.sharp") {
		filter_dict["value"] = 0;
	}

	EMData *img1 = this->copy();
	EMData *img2 = with->copy();

	int img1_nx = img1->get_xsize();
	int img1_ny = img1->get_ysize();
	int img1_nz = img1->get_zsize();
	int img1_size = img1_nx * img1_ny * img1_nz;

	float img1min = img1->get_attr("minimum");
	img1->add(-img1min);

	float img2min = img2->get_attr("minimum");
	img2->add(-img2min);

	filter_dict["outer_radius"] = radius;

	EMData *img1_copy = img1->copy();
	img1_copy->to_one();
	img1_copy->process(mask_filter, filter_dict);
	img1_copy->process("eman1.xform.phaseorigin");

	int num = 0;
	float *img1_copy_data = img1_copy->get_data();

	for (int i = 0; i < img1_size; i++) {
		if (img1_copy_data[i] == 1) {
			num++;
		}
	}

	img2->process(mask_filter, filter_dict);

	float *img2_data = img2->get_data();
	double lsum = 0;
	double sumsq = 0;

	for (int i = 0; i < img1_size; i++) {
		lsum += img2_data[i];
		sumsq += img2_data[i] * img2_data[i];
	}

	float sq = (float)((num * sumsq - lsum * lsum) / (num * num));
	if (sq < 0) {
		LOGERR("sigma < 0");
		throw ImageFormatException("image sigma < 0");
	}

	float mean = (float)lsum / num;
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

	EMData *img2_copy = img2->copy();
	if( img2 )
	{
		delete img2;
		img2 = 0;
	}

	img2_copy->process(mask_filter, filter_dict);
	img2_copy->process("eman1.xform.phaseorigin");

	if( img1_copy )
	{
		delete img1_copy;
		img1_copy = 0;
	}

	EMData *img1_copy2 = img1->copy();

	img1_copy2->process("eman1.math.squared");

	EMData *ccf = img1->calc_ccf(img2_copy);
	if( img2_copy )
	{
		delete img2_copy;
		img2_copy = 0;
	}

	ccf->mult(img1_size);

	EMData *conv1 = img1->convolute(img1_copy2);
	if( img1 )
	{
		delete img1;
		img1 = 0;
	}

	conv1->mult(img1_size);
	conv1->mult(1.0f / num);

	EMData *conv2 = img1_copy2->convolute(img1_copy2);
	if( img1_copy2 )
	{
		delete img1_copy2;
		img1_copy2 = 0;
	}

	conv2->mult(img1_size);
	conv1->process("eman1.math.squared");
	conv1->mult(1.0f / (num * num));

	EMData *conv2_copy = conv2->copy();
	if( conv2 )
	{
		delete conv2;
		conv2 = 0;
	}

	conv2_copy->sub(*conv1);
	if( conv1 )
	{
		delete conv1;
		conv1 = 0;
	}

	conv2_copy->mult(1.0f / num);
	conv2_copy->process("eman1.math.sqrt");

	EMData *ccf_copy = ccf->copy();
	if( ccf )
	{
		delete ccf;
		ccf = 0;
	}

	ccf_copy->mult(1.0f / num);

	float *lcfd = ccf_copy->get_data();
	float *vdd = conv2_copy->get_data();

	for (int i = 0; i < img1_size; i++) {
		if (vdd[i] > 0) {
			lcfd[i] /= vdd[i];
		}
	}
	if( conv2_copy )
	{
		delete conv2_copy;
		conv2_copy = 0;
	}

	ccf_copy->done_data();
	EMData *lcf = ccf_copy->copy();
	if( ccf_copy )
	{
		delete ccf_copy;
		ccf_copy = 0;
	}

	EXITFUNC;
	return lcf;
}

EMData *EMData::convolute(EMData * with)
{
	ENTERFUNC;

	EMData *f1 = do_fft();
	if (!f1) {
		LOGERR("FFT returns NULL image");
		throw NullPointerException("FFT returns NULL image");
	}

	f1->ap2ri();

	EMData *cf = 0;
	if (with) {
		cf = with->do_fft();
		if (!cf) {
			LOGERR("FFT returns NULL image");
			throw NullPointerException("FFT returns NULL image");
		}
		cf->ap2ri();
	}
	else {
		cf = f1->copy();
	}

	if (with && !EMUtil::is_same_size(f1, cf)) {
		LOGERR("images not same size");
		throw ImageFormatException("images not same size");
	}

	float *rdata1 = f1->get_data();
	float *rdata2 = cf->get_data();
	int cf_size = cf->get_xsize() * cf->get_ysize() * cf->get_zsize();

	float re,im;
	for (int i = 0; i < cf_size; i += 2) {
		re = rdata1[i] * rdata2[i] - rdata1[i + 1] * rdata2[i + 1];
		im = rdata1[i + 1] * rdata2[i] + rdata1[i] * rdata2[i + 1];
		rdata2[i]=re;
		rdata2[i+1]=im;
	}

	cf->done_data();
	EMData *f2 = cf->do_ift();

	if( cf )
	{
		delete cf;
		cf = 0;
	}

	if( f1 )
	{
		delete f1;
		f1=0;
	}

	EXITFUNC;
	return f2;
}


void EMData::common_lines(EMData * image1, EMData * image2,
						  int mode, int steps, bool horizontal)
{
	ENTERFUNC;

	if (!image1 || !image2) {
		throw NullPointerException("NULL image");
	}

	if (mode < 0 || mode > 2) {
		throw OutofRangeException(0, 2, mode, "invalid mode");
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
		throw ImageFormatException("images not same sizes");
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

							float a1 = (float) hypot(im1[i2], im1[i2 + 1]);
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
							rdata[l] += (float) (hypot(im1[i2], im1[i2 + 1]) * hypot(im2[j2], im2[j2 + 1]));
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
		if( tmp_array )
		{
			delete[]tmp_array;
			tmp_array = 0;
		}
	}

	if( im1 )
	{
		delete[]im1;
		im1 = 0;
	}

	if( im2 )
	{
		delete im2;
		im2 = 0;
	}


	image1->done_data();
	image2->done_data();
	if( image1 )
	{
		delete image1;
		image1 = 0;
	}
	if( image2 )
	{
		delete image2;
		image2 = 0;
	}
	done_data();
	update();
	EXITFUNC;
}



void EMData::common_lines_real(EMData * image1, EMData * image2,
							   int steps, bool horiz)
{
	ENTERFUNC;

	if (!image1 || !image2) {
		throw NullPointerException("NULL image");
	}

	if (!EMUtil::is_same_size(image1, image2)) {
		throw ImageFormatException("images not same size");
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

	if( image1_copy )
	{
		delete image1_copy;
		image1_copy = 0;
	}

	if( image2_copy )
	{
		delete image2_copy;
		image2_copy = 0;
	}

	if( im1 )
	{
		delete[]im1;
		im1 = 0;
	}

	if( im2 )
	{
		delete[]im2;
		im2 = 0;
	}
	EXITFUNC;
}


void EMData::cut_slice(const EMData * map, float dz, Transform3D * ort,
					   bool interpolate, float dx, float dy)
{
	ENTERFUNC;

	if (!map) {
		throw NullPointerException("NULL image");
	}

	Transform3D r(0, 0, 0); // EMAN by default
	if (!ort) {
		ort = &r;
	}

	float *sdata = map->get_data();
	float *ddata = get_data();

	int map_nx = map->get_xsize();
	int map_ny = map->get_ysize();
	int map_nz = map->get_zsize();
	int map_nxy = map_nx * map_ny;

	float mdz0 = dz * (*ort)[0][2] + map_nx / 2;
	float mdz1 = dz * (*ort)[1][2] + map_ny / 2;
	float mdz2 = dz * (*ort)[2][2] + map_nz / 2;

	for (int y = 0; y < ny; y++) {
		int y2 = (int) (y - ny / 2 - dy);
		float my2_0 = y2 * (*ort)[0][1] * y2 + mdz0;
		float my2_1 = y2 * (*ort)[1][1] * y2 + mdz1;
		float my2_2 = y2 * (*ort)[2][1] * y2 + mdz2;

		for (int x = 0; x < nx; x++) {
			int x2 = (int) (x - nx / 2 - dx);
			float xx = x2 * (*ort)[0][0] + my2_0;
			float yy = x2 * (*ort)[1][0] + my2_1;
			float zz = x2 * (*ort)[2][0] + my2_2;
			int l = x + y * nx;

			if (xx < 0 || yy < 0 || zz < 0 || xx > map_nx - 2 ||
				yy > map_ny - 2 || zz > map_nz - 2) {
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

	EXITFUNC;
}


void EMData::uncut_slice(EMData * map, float dz, Transform3D * ort, float dx, float dy)
{
	ENTERFUNC;

	if (!map) {
		throw NullPointerException("NULL image");
	}

	Transform3D r( 0, 0, 0); // EMAN by default
	if (!ort) {
		ort = &r;
	}

	float *ddata = map->get_data();
	float *sdata = get_data();

	int map_nx = map->get_xsize();
	int map_ny = map->get_ysize();
	int map_nz = map->get_zsize();
	int map_nxy = map_nx * map_ny;

	float mdz0 = dz * (*ort)[0][2] + map_nx / 2;
	float mdz1 = dz * (*ort)[1][2] + map_ny / 2;
	float mdz2 = dz * (*ort)[2][2] + map_nz / 2;

	for (int y = 0; y < ny; y++) {
		int y2 = (int) (y - ny / 2 - dy);

		float my2_0 = y2 * (*ort)[0][1] + mdz0;
		float my2_1 = y2 * (*ort)[1][1] + mdz1;
		float my2_2 = y2 * (*ort)[2][1] + mdz2;

		for (int x = 0; x < nx; x++) {
			int x2 = (int) (x - nx / 2 - dx);

			float xx = x2 * (*ort)[0][0] + my2_0;
			float yy = x2 * (*ort)[1][0] + my2_1;
			float zz = x2 * (*ort)[2][0] + my2_2;

			if (xx >= 0 && yy >= 0 && zz >= 0 && xx <= map_nx - 2 && yy <= map_ny - 2
				&& zz <= map_nz - 2) {
				int k = Util::round(xx) + Util::round(yy) * map_nx + Util::round(zz) * map_nxy;
				ddata[k] = sdata[x + y * nx];
			}
		}
	}

	done_data();
	map->done_data();
	EXITFUNC;
}

vector<Pixel> EMData::calc_highest_locations(float threshold)
{
	ENTERFUNC;

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

	EXITFUNC;
	return result;
}


float EMData::get_edge_mean() const
{
	ENTERFUNC;

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
		edge_mean = (float)edge_sum / (nx * 2 + ny * 2);
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

		edge_mean = (float)edge_sum / (di * 2);
	}
	EXITFUNC;

	return edge_mean;
}

float EMData::get_circle_mean()
{
	ENTERFUNC;

	static bool busy = false;
	static EMData *mask = 0;

	while (busy);
	busy = true;

	if (!mask || !EMUtil::is_same_size(this, mask)) {
		if (!mask) {
			mask = new EMData();
		}
		mask->set_size(nx, ny, nz);
		mask->to_one();

		float radius = (float)(ny / 2 - 2);
		mask->process("eman1.mask.sharp", Dict("inner_radius", radius - 1,
									   "outer_radius", radius + 1));

		int n = 0;
		float *d = mask->get_data();

		for (int i = 0; i < nx * ny * nz; i++) {
			if (d[i]) {
				n++;
			}
		}
		mask->mult(1.0f / n);
	}

	float result = dot(mask);
	busy = false;

	EXITFUNC;
	return result;
}


// just a shortcut for cmp("dot")
float EMData::dot(EMData * with)
{
	ENTERFUNC;
	if (!with) {
		throw NullPointerException("Null EMData Image");
	}
	DotCmp dot_cmp;
	float r = -dot_cmp.cmp(this, with);
	EXITFUNC;
	return r;
}

float EMData::sget_value_at(int x, int y, int z) const
{
	if (x < 0 || x >= nx || y < 0 || y >= ny || z < 0 || z >= nz) {
		return 0;
	}
	return rdata[x + y * nx + z * nx * ny];
}

float EMData::sget_value_at(int x, int y) const
{
	if (x < 0 || x >= nx || y < 0 || y >= ny) {
		return 0;
	}
	return rdata[x + y * nx];
}

float EMData::sget_value_at(size_t i) const
{
	size_t size = nx*ny;
	size *= nz;
	if (i >= size) {
		return 0;
	}
	return rdata[i];
}

float EMData::sget_value_at_interp(float xx, float yy) const
{
	int x = static_cast < int >(floor(xx));
	int y = static_cast < int >(floor(yy));

	float p1 = sget_value_at(x, y);
	float p2 = sget_value_at(x + 1, y);
	float p3 = sget_value_at(x + 1, y + 1);
	float p4 = sget_value_at(x, y + 1);

	float result = Util::bilinear_interpolate(p1, p2, p3, p4, xx - x, yy - y);
	return result;
}

float EMData::sget_value_at_interp(float xx, float yy, float zz) const
{
	int x = (int) floor(xx);
	int y = (int) floor(yy);
	int z = (int) floor(zz);

	float p1 = sget_value_at(x, y, z);
	float p2 = sget_value_at(x + 1, y, z);
	float p3 = sget_value_at(x, y + 1, z);
	float p4 = sget_value_at(x + 1, y + 1, z);

	float p5 = sget_value_at(x, y, z + 1);
	float p6 = sget_value_at(x + 1, y, z + 1);
	float p7 = sget_value_at(x, y + 1, z + 1);
	float p8 = sget_value_at(x + 1, y + 1, z + 1);

	float result = Util::trilinear_interpolate(p1, p2, p3, p4, p5, p6, p7, p8,
											   xx - x, yy - y, zz - z);

	return result;
}

void EMData::save_byteorder_to_dict(ImageIO * imageio)
{
	string image_endian = "ImageEndian";
	string host_endian = "HostEndian";

	if (imageio->is_image_big_endian()) {
		attr_dict[image_endian] = "big";
	}
	else {
		attr_dict[image_endian] = "little";
	}

	if (ByteOrder::is_host_big_endian()) {
		attr_dict[host_endian] = "big";
	}
	else {
		attr_dict[host_endian] = "little";
	}
}

void EMData::print_image(const string str, ostream& out) {
	out << "Printing EMData object: " << str << std::endl;
	MArray3D mat = get_3dview();
	int nx = get_xsize();
	int ny = get_ysize();
	int nz = get_zsize();
	for (int iz = 0; iz < nz; iz++) {
		out << "(z = " << iz << " slice)" << std::endl;
		for (int ix = 0; ix < nx; ix++) {
			for (int iy = 0; iy < ny; iy++) {
				out << setiosflags(std::ios::fixed)
					<< setiosflags(std::ios_base::scientific)
					<< std::setw(12)
					 << std::setprecision(5) << mat[ix][iy][iz] << "  ";
				if (((iy+1) % 6) == 0) {
					out << std::endl << "   ";
				}
			}
			out << std::endl;
		}
	}
}

EMData * EMData::real() //real part has half of x dimension
{
	ENTERFUNC;
	
	EMData * e = new EMData();

	if( is_real() ) // a real image, return a copy of itself
	{
		e = this->copy();
	}
	else //for a complex image
	{
		if( !is_ri() ) //complex image in amplitude/phase foramt, convert it to real/imaginary first
		{
			ap2ri();
		}
		int nx = get_xsize();
		int ny = get_ysize();
		int nz = get_zsize();
		e->set_size(nx/2, ny, nz);
		float * edata = e->get_data();
		for( int i=0; i<nx; i++ )
		{
			for( int j=0; j<ny; j++ )
			{
				for( int k=0; k<nz; k++ )
				{
					if( i%2 == 0 )
					{
						//complex data in format [real, complex, real, complex...]
						edata[i/2+j*(nx/2)+k*(nx/2)*ny] = rdata[i+j*nx+k*nx*ny];
					}
				}
			}
		}
	}
	
	e->set_complex(false);
	if(e->get_ysize()==1 && e->get_zsize()==1) {
		e->set_complex_x(false);
	}
	e->update_stat();
	return e;

	EXITFUNC;
}

EMData * EMData::imag()
{
	ENTERFUNC;

	EMData * e = new EMData();

	if( is_real() ) {	//a real image has no imaginary part, throw exception
		throw InvalidCallException("No imaginary part for a real image, this function call require a complex image.");
	}
	else {	//for complex image
		if( !is_ri() ) {	//complex image in amplitude/phase foramt, convert it to real/imaginary first
			ap2ri();
		}
		int nx = get_xsize();
		int ny = get_ysize();
		int nz = get_zsize();
		e->set_size(nx/2, ny, nz);
		float * edata = e->get_data();
		for( int i=0; i<nx; i++ ) {
			for( int j=0; j<ny; j++ ) {
				for( int k=0; k<nz; k++ ) {
					if( i%2 == 1 ) {
						//complex data in format [real, complex, real, complex...]
						edata[i/2+j*(nx/2)+k*(nx/2)*ny] = rdata[i+j*nx+k*nx*ny];
					}
				}
			}
		}
	}
	
	e->set_complex(false);
	if(e->get_ysize()==1 && e->get_zsize()==1) {
		e->set_complex_x(false);
	}
	e->update_stat();
	return e;

	EXITFUNC;
}

EMData * EMData::real2complex(const float img)
{
	ENTERFUNC;

	if( is_complex() ) {
		throw InvalidCallException("This function call only apply to real image");
	}
	
	EMData * e = new EMData();
	int nx = get_xsize();
	int ny = get_ysize();
	int nz = get_zsize();
	e->set_size(nx*2, ny, nz);
	
	MArray3D edata = e->get_3dview();
	MArray3D data  = this->get_3dview();

	for( int k=0; k<nz; k++ ) {
		for( int j=0; j<ny; j++ ) {
			for( int i=0; i<nx; i++ ) {			
				edata[i*2][j][k] = data[i][j][k];
				edata[i*2+1][j][k] = img;
			}
		}
	}
	
	e->set_complex(true);
	if(e->get_ysize()==1 && e->get_zsize()==1) {
		e->set_complex_x(true);
	}
	e->set_ri(true);
	e->update_stat();
	return e;

	EXITFUNC;
}

EMData* EMData::symvol(string symmetry) {
	ENTERFUNC;
	int nsym = Transform3D::get_nsym(symmetry); // number of symmetries
	Transform3D sym;
	int llim = -nx/2;
	int ulim = (nx/2) -1 + (nx % 2);
	// set up output volume
	EMData* svol = new EMData;
	svol->set_size(nx, ny, nz);
	svol->to_zero();
	// set up coord grid
	const int nsize = 27;
	int x[nsize], y[nsize], z[nsize];
	float f[nsize];
	for (int i = 0; i < nsize; i+=3) {
		x[i] = -1;
		x[i+1] = 0;
		x[i+2] = 1;
		int imod = (i/3) % 3;
		y[i] = imod - 1;
		y[i+1] = imod - 1;
		y[i+2] = imod - 1;
	}
	for (int i = 0; i < nsize; i++) z[i] = (i/9) - 1;
	// calc radius within which the rotation is valid
	int iradrt = (0 == nx % 2) ? nx/2 - 1 : nx / 2;
	int iradi = (iradrt-1)*(iradrt-1);
	// actual work -- loop over symmetries and symmetrize
	MArray3D q1 = get_3dview();
	MArray3D q2 = svol->get_3dview();
	for (int isym = 0; isym < nsym; isym++) {
		Transform3D rm = sym.get_sym(symmetry, isym);
		if ((1.0 == rm[0][0]) && (1.0 == rm[1][1]) && (1.0 == rm[2][2])) {
			// symmetry is the identity matrix
			for (int iz = 0; iz < nz; iz++)
				for (int iy = 0; iy < ny; iy++)
					for (int ix = 0; ix < nx; ix++)
						q2[ix][iy][iz] += q1[ix][iy][iz];
		} else {
			// symmetry is something interesting
			Vec3f qrt, qr;
			for (int iz = llim; iz <= ulim; iz++) {
				for (int iy = llim; iy <= ulim; iy++) {
					qrt[0] = rm[0][1]*iy + rm[0][2]*iz;
					qrt[1] = rm[1][1]*iy + rm[1][2]*iz;
					qrt[2] = rm[2][1]*iy + rm[2][2]*iz;
					for (int ix = llim; ix <= ulim; ix++) {
						int icrd = ix*ix + iy*iy + iz*iz;
						if (icrd <= iradi) {
							qr[0] = qrt[0] + rm[0][0]*ix;
							qr[1] = qrt[1] + rm[1][0]*ix;
							qr[2] = qrt[2] + rm[2][0]*ix;
							// iox -- integer location in -nx/2...nx/2
							int iox = static_cast<int>(floorf(qr[0]));
							int ioy = static_cast<int>(floorf(qr[1]));
							int ioz = static_cast<int>(floorf(qr[2]));
							// dx -- offset from integer array
							float dx = qr[0] - iox;
							float dy = qr[1] - ioy;
							float dz = qr[2] - ioz;
							// find intensities on 3x3x3 grid
							for (int i = 0; i < nsize; i++) {
								int jx = iox + x[i] - llim;
								int jy = ioy + y[i] - llim;
								int jz = ioz + z[i] - llim;
								f[i] = q1[jx][jy][jz];
							}
							// eval intensity at px, py, pz
							int jx = ix - llim;
							int jy = iy - llim;
							int jz = iz - llim;
							q2[jx][jy][jz] += Util::triquad(dx,dy,dz,f);
						} else {
							// rotated position is outside volume
							int jx = ix - llim;
							int jy = iy - llim;
							int jz = iz - llim;
							q2[jx][jy][jz] += q1[jx][jy][jz];
						}
					}
				}
			}
		}
	}
	// normalize
	for (int iz = 0; iz < nz; iz++)
		for (int iy = 0; iy < ny; iy++)
			for (int ix = 0; ix < nx; ix++)
				q2[ix][iy][iz] /= nsym;
	svol->done_data();
	svol->update();
	EXITFUNC;
	return svol;
}

EMData*
EMData::rot_trans2D(float ang, float delx, float dely) {
	if (1 >= ny) 
		throw ImageDimensionException("Can't rotate 1D image");
	if (1 < nz) 
		throw ImageDimensionException("Volume not currently supported");
	if (0.f == ang) {
		EMData* ret = copy();
		return ret;
	}
	update();
	float background = get_attr("mean");
	if (ang > pi) ang -= static_cast<float>(twopi);
	if (ang < -pi) ang += static_cast<float>(twopi);
	float cang = cos(ang);
	float sang = sin(ang);
	EMData* ret = copy_head();
	// center of the image
	int xc = nx/2;
	int yc = ny/2;
	// shift center for rotation (if desired)
	float shiftxc = xc + delx;
	float shiftyc = yc + dely;
	MArray2D in = get_2dview();
	MArray2D out = ret->get_2dview();
	for (int iy = 0; iy < ny; iy++) {
		float y = float(iy) - shiftyc;
		float ycang = y*cang + shiftyc;
		float ysang = -y*sang + shiftxc;
		for (int ix = 0; ix < nx; ix++) {
			out[ix][iy] = background;
			float x = float(ix) - shiftxc;
			float xold = x*cang + ysang;
			float yold = x*sang + ycang;
			int iyold = int(yold);
			float q = yold - float(iyold);
			float qcomp = 1.f - q;
			int ixold = int(xold);
			// Note: nx-2 or ny-2 below because need room for
			// (forward) interpolation
			if ((yold>=0 && iyold<=(ny-2)) && (xold>=0 && ixold<=(nx-2))) {
				// inside boundaries of input image
				float p = xold - ixold;
				float pcomp = 1.f - p;
				out[ix][iy] = q*(pcomp*in[ixold][iyold+1]
						         + p*in[ixold+1][iyold+1])
					        + qcomp*(pcomp*in[ixold][iyold]
									 + p*in[ixold+1][iyold]);
			}
		}
	}
	ret->done_data();
	ret->update();
	return ret;
}

EMData*
EMData::rot_scale_trans2D(float ang, float scale, float delx,
		                  float dely) {
	if (1 >= ny)
		throw ImageDimensionException("Can't rotate 1D image");
	if (0.f == scale) scale = 1.f; // silently fix common user error
	EMData* ret = copy_head();
	delx = fmod(delx, float(nx));
	dely = fmod(dely, float(ny));
	// center of image
	int xc = nx/2;
	int yc = ny/2;
	// shifted center for rotation
	float shiftxc = xc + delx;
	float shiftyc = yc + dely;
	// bounds if origin at center
	float ymin = -ny/2.0f;
	float xmin = -nx/2.0f;
	float ymax = -ymin;
	float xmax = -xmin;
	if (0 == nx%2) xmax--;
	if (0 == ny%2) ymax--;
	// trig
	float cang = cos(ang);
	float sang = sin(ang);
	MArray3D out = ret->get_3dview();
	for (int iz = 0; iz < nz; iz++) {
		for (int iy = 0; iy < ny; iy++) {
			float y = float(iy) - shiftyc;
		#ifdef _WIN32
			if (y < ymin) y = _MIN(y+ny,ymax);
			if (y > ymax) y = _MAX(y-ny,ymin);
		#else
			if (y < ymin) y = std::min(y+ny,ymax);
			if (y > ymax) y = std::max(y-ny,ymin);
		#endif	//_WIN32
			float ycang = y*cang/scale + yc;
			float ysang = -y*sang/scale + xc;
			for (int ix = 0; ix < nx; ix++) {
				float x = float(ix) - shiftxc;
			#ifdef _WIN32
				if (x < xmin) x = _MIN(x+nx,xmax);
				if (x > xmax) x = _MAX(x-nx,xmin);
			#else
				if (x < xmin) x = std::min(x+nx,xmax);
				if (x > xmax) x = std::max(x-nx,xmin);
			#endif	//_WIN32
				float xold = x*cang/scale + ysang;
				float yold = x*sang/scale + ycang;
				out[ix][iy][iz] =
					Util::quadri(this, xold, yold, iz);
			}
		}
	}
	return ret;
}

float EMData::getconvpt2d_kbi0(float x, float y, 
		Util::KaiserBessel::kbi0_win win, int size) {
	const int nxhalf = nx/2;
	const int nyhalf = ny/2;
	const int bd = size/2;
	float* wxarr = new float[size];
	float* wyarr = new float[size];
	float* wx = wxarr + bd; // wx[-bd] = wxarr[0]
	float* wy = wyarr + bd;
	int ixc = int(x + 0.5f*Util::sgn(x));
	int iyc = int(y + 0.5f*Util::sgn(y));
	if (abs(ixc) > nxhalf)
		throw InvalidValueException(ixc, "getconv: X value out of range");
	if (abs(iyc) > nyhalf)
		throw InvalidValueException(ixc, "getconv: Y value out of range");
	for (int i = -bd; i <= bd; i++) {
		int iyp = iyc + i;
		wy[i] = win(y - iyp);
		int ixp = ixc + i;
		wx[i] = win(x - ixp);
	}
	MArray2D imgarr = get_2dview(-nxhalf, -nyhalf);
	float conv = 0.f, wsum = 0.f;
	for (int iy = -bd; iy <= bd; iy++) {
		int iyp = iyc + iy;
		for (int ix = -bd; ix <= bd; ix++) {
			int ixp = ixc + ix;
			float wg = wx[ix]*wy[iy];
			conv += imgarr[ixp][iyp]*wg;
			wsum += wg;
		}
	}
	delete [] wxarr;
	delete [] wyarr;
	//return conv/wsum;
	return conv;
}

#if 0 // FIXME: broken
EMData* EMData::rotconvtrunc2d_kbi0(float ang, float alpha, int size) {
    // truncate anything outside r=min(nx/2,ny/x)-window
    int nx = get_xsize();
    int ny = get_ysize();
    int nxhalf = nx/2;
    int nyhalf = ny/2;
#ifdef _WIN32
	float rmax = float(_MIN(nxhalf,nyhalf)) - float(size/2);
#else
    float rmax = float(std::min(nxhalf,nyhalf)) - float(size/2);
#endif	//_WIN32
    float rmax2 = rmax*rmax;
	Util::KaiserBessel kb(alpha, size-1);
    if (1 >= ny) 
        throw ImageDimensionException("Can't rotate 1D image");
	EMData* ret = copy_head();
    float cod = cos(ang);
    float sid = sin(ang);
    MArray2D out = ret->get_2dview(-nxhalf,-nyhalf);
    MArray2D in  = get_2dview(-nxhalf,-nyhalf);
    for (int iy = -nyhalf; iy < nyhalf + ny%2; iy++) {
        float ycod = iy*cod;
        float ysid = -iy*sid;
        for (int ix = -nxhalf; ix < nxhalf + nx%2; ix++) {
            if (ix*ix + iy*iy <= rmax2) {
                float yold = ix*sid + ycod;
                float xold = ix*cod + ysid;
                out[ix][iy] =  getconvpt2d_kbi0(xold, yold, kb.get_kbi0_win());
            } else {
                out[ix][iy] = in[ix][iy];
            }
        }
    }
	ret->done_data();
    return ret;
}
#endif // 0

complex<float> EMData::extractpoint(float nuxnew, float nuynew,
		Util::KaiserBessel& kb) {
	if (2 != get_ndim())
		throw ImageDimensionException("extractpoint needs a 2-D image.");
	if (!is_complex()) 
		throw ImageFormatException("extractpoint requires a fourier image");
	int nxreal = nx - 2;
	if (nxreal != ny)
		throw ImageDimensionException("extractpoint requires ny == nx");
	int nhalf = nxreal/2; 
	int kbsize = kb.get_window_size();
	int kbmin = -kbsize/2;
	int kbmax = -kbmin;
	bool flip = (nuxnew < 0.f);
	if (flip) {
		nuxnew *= -1;
		nuynew *= -1;
	}
	// put (xnew,ynew) on a grid.  The indices will be wrong for
	// the Fourier elements in the image, but the grid sizing will
	// be correct.
	int ixn = int(Util::round(nuxnew));
	int iyn = int(Util::round(nuynew));
	// displacements of (xnew,ynew) from the grid
	float nuxdispl = nuxnew - ixn;
	float nuydispl = nuynew - iyn;
	// set up some temporary weighting arrays
	float* wy0 = new float[kbmax - kbmin + 1];
	float* wy = wy0 - kbmin; // wy[kbmin:kbmax]
	float* wx0 = new float[kbmax - kbmin + 1];
	float* wx = wx0 - kbmin;
	for (int i = kbmin; i <= kbmax; i++) {
		wy[i] = kb.i0win_tab(nuydispl - i);
		//wy[i] = (0 == i) ? 1.f : 0.f; // FIXME: remove after debugging
		wx[i] = kb.i0win_tab(nuxdispl - i);
		//wx[i] = (0 == i) ? 1.f : 0.f; // FIXME: remove after debugging
	}
	// restrict loops to non-zero elements
	int iymin = 0;
	for (int iy = kbmin; iy <= -1; iy++) {
		if (wy[iy] != 0.f) {
			iymin = iy;
			break;
		}
	}
	int iymax = 0;
	for (int iy = kbmax; iy >= 1; iy--) {
		if (wy[iy] != 0.f) {
			iymax = iy;
			break;
		}
	}
	int ixmin = 0;
	for (int ix = kbmin; ix <= -1; ix++) {
		if (wx[ix] != 0.f) {
			ixmin = ix;
			break;
		}
	}
	int ixmax = 0;
	for (int ix = kbmax; ix >= 1; ix--) {
		if (wx[ix] != 0.f) {
			ixmax = ix;
			break;
		}
	}
	double wsum = 0.f;
	complex<float> result(0.f,0.f);
	if ((ixn >= -kbmin) && (ixn <= nhalf-1-kbmax)
			&& (iyn >= -nhalf-kbmin) && (iyn <= nhalf-1-kbmax)) {
		// (xin,yin) not within window border from the edge
		for (int iy = iymin; iy <= iymax; iy++) {
			int iyp = iyn + iy;
			for (int ix = ixmin; ix <= ixmax; ix++) {
				int ixp = ixn + ix;
				float w = wx[ix]*wy[iy];
				complex<float> val = cmplx(ixp,iyp);
				result += val*w;
				wsum += w;
			}
		}
	} else {
		// points that "stick out"
		for (int iy = iymin; iy <= iymax; iy++) {
			int iyp = iyn + iy;
			for (int ix = ixmin; ix <= ixmax; ix++) {
				int ixp = ixn + ix;
				bool mirror = false;
				int ixt= ixp, iyt= iyp;
				if ((ixt > nhalf) || (ixt < -nhalf)) {
					ixt = Util::sgn(ixt)*(nxreal - abs(ixt));
					iyt *= -1;
					mirror = !mirror;
				}
				if ((iyt >= nhalf) || (iyt < -nhalf)) {
					if (ixt != 0) {
						ixt = -ixt;
						iyt = Util::sgn(iyt)*(nxreal-abs(iyt));
						mirror = !mirror;
					} else {
						iyt -= Util::sgn(iyt)*nxreal;
					}
				}
				if (ixt < 0) {
					ixt = -ixt;
					iyt = -iyt;
					mirror = !mirror;
				}
				if (iyt == nhalf) iyt = -nhalf;
				float w = wx[ix]*wy[iy];
				wsum += w;
				complex<float> val = this->cmplx(ixt,iyt);
				if (mirror) 
					result += conj(val)*w;
				else
					result += val*w;
			}
		}
	}
	if (flip) 
		result = conj(result)/static_cast<float>(wsum);
	else
		result /= static_cast<float>(wsum);
	delete [] wx0;
	delete [] wy0;
	return result;
}

void EMData::center_padded() {
	int npad = get_attr("npad");
	if (1 == npad) return;
	EMData& self = *this;
	self.set_array_offsets();
	int nxorig = nx/npad;
	int nyorig = ny/npad;
	int nxcorner = (nx - nxorig)/2;
	int nycorner = (ny - nyorig)/2;
	for (int iy = nyorig-1; iy >= 0; iy--) 
		for (int ix = nxorig-1; ix >= 0; ix--)
			std::swap(self(nxcorner+ix,nycorner+iy),self(ix,iy));
}

void EMData::fft_shuffle() {
	if (!is_complex()) 
		throw ImageFormatException("fft_shuffle requires a fourier image");
	vector<int> offsets = get_array_offsets();
	set_array_offsets(); // clear offsets before shuffling
	EMData& self = *this;
	if (0 == ny%2) {
		int offset = ny/2;
		for (int iy = 0; iy < ny/2; iy++) 
			// swap column iy and iy + offset
			for (int ix = 0; ix < nx; ix++)
				std::swap(self(ix,iy),self(ix,iy+offset));
	} else {
		//stupid algorithm; too lazy to find better
		int shifts = ny/2 + int(is_shuffled());
		for (int ix = 0; ix < nx; ix++) {
			for (int shift = 0; shift < shifts; shift++) {
				float temp = self(ix,0);
				for (int iy = 1; iy < ny; iy++) 
					self(ix,iy-1) = self(ix,iy);
				self(ix,ny-1) = temp;
			}
		}
	}
	set_shuffled(!is_shuffled()); // toggle
	set_array_offsets(offsets); // reset offsets
	done_data();
}

EMData* EMData::fouriergridrot2d(float ang, Util::KaiserBessel& kb) {
	if (2 != get_ndim())
		throw ImageDimensionException("fouriergridrot2d needs a 2-D image.");
	if (!is_complex()) 
		throw ImageFormatException("fouriergridrot2d requires a fourier image");
	if (!is_shuffled()) 
		fft_shuffle();

	int nxreal = nx - 2 + int(is_fftodd());
	if (nxreal != ny)
		throw ImageDimensionException("fouriergridrot2d requires ny == nx(real)");
	if (0 != nxreal%2)
		throw ImageDimensionException("fouriergridrot2d needs an even image.");
	int nxhalf = nxreal/2;
	float nxhalf2 = nxhalf*float(nxhalf);
	int nyhalf = ny/2;
	EMData* result = copy();
	set_array_offsets(0,-nyhalf);
	result->set_array_offsets(0,-nyhalf);
	float cang = cos(ang);
	float sang = sin(ang);
	for (int iy = -nyhalf; iy < nyhalf; iy++) {
		float ycang = iy*cang;
		float ysang = -iy*sang;
		float iy2 = iy*float(iy);
		for (int ix = 0; ix <= nxhalf; ix++) {
			float ix2 = ix*float(ix);
			if (ix2 + iy2 <= nxhalf2) {
				float nuyold = ix*sang + ycang;
				float nuxold = ix*cang + ysang;
				result->cmplx(ix,iy) = extractpoint(nuxold,nuyold,kb);
			} else {
				result->cmplx(ix,iy) = complex<float>(0.f,0.f);
			}
		}
	}
	result->set_array_offsets();
	result->fft_shuffle(); // reset to an unshuffled result
	result->done_data();
	set_array_offsets();
	fft_shuffle(); // reset to an unshuffled complex image
	return result;
}


Dict EMData::masked_stats(const EMData* mask) {
	if (is_complex())
		throw ImageFormatException(
				"Complex images not supported by EMData::masked_stats");
	float* ptr = get_data();
	float* mptr = mask->get_data();
	long double sum1 = 0.L;
	long double sum2 = 0.L;
	long nmask = 0L;
	for (long i = 0; i < nx*ny*nz; i++,ptr++,mptr++) {
		if (*mptr > 0.5f) {
			nmask++;
			sum1 += *ptr;
			sum2 += (*ptr)*(*ptr);
		}
	}
	float avg = static_cast<float>(sum1/nmask);
	float sig2 = static_cast<float>(sum2/nmask - avg*avg);
	float sig = sqrt(sig2);
	Dict mydict;
	mydict["avg"] = avg; mydict["sigma"] = sig; mydict["nmask"] = int(nmask);
	return mydict;
}

/* vim: set ts=4 noet: */
