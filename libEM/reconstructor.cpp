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

#include "reconstructor.h"
#include "plugins/reconstructor_template.h"
#include "ctf.h"
#include "emassert.h"
#include "symmetry.h"
#include <cstring>
#include <fstream>
#include <iomanip>
#include <boost/bind.hpp>
#include <boost/format.hpp>

#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_fit.h>

#ifdef EMAN2_USING_CUDA
#include "cuda/cuda_reconstructor.h"
#endif

using namespace EMAN;
using std::complex;


#include <iostream>
using std::cerr;
using std::endl;
using std::cout; // for debug

#include <algorithm>
// find, for_each

#include <iomanip>
using std::setprecision;

template < typename T > void checked_delete( T*& x )
{
    typedef char type_must_be_complete[ sizeof(T)? 1: -1 ];
    (void) sizeof(type_must_be_complete);
    delete x;
    x = NULL;
}

const string FourierReconstructor::NAME = "fourier";
const string FourierPlaneReconstructor::NAME = "fourier_plane";
const string FourierReconstructorSimple2D::NAME = "fouriersimple2D";
const string WienerFourierReconstructor::NAME = "wiener_fourier";
const string BackProjectionReconstructor::NAME = "back_projection";
const string nn4Reconstructor::NAME = "nn4";
const string nn4_rectReconstructor::NAME = "nn4_rect";
const string nnSSNR_Reconstructor::NAME = "nnSSNR";
const string nn4_ctfReconstructor::NAME = "nn4_ctf";
const string nn4_ctfwReconstructor::NAME = "nn4_ctfw";
const string nn4_ctf_rectReconstructor::NAME = "nn4_ctf_rect";
const string nnSSNR_ctfReconstructor::NAME = "nnSSNR_ctf";

template <> Factory < Reconstructor >::Factory()
{
	force_add<FourierReconstructor>();
	force_add<FourierPlaneReconstructor>();
	force_add<FourierReconstructorSimple2D>();
//	force_add(&BaldwinWoolfordReconstructor::NEW);
	force_add<WienerFourierReconstructor>();
	force_add<BackProjectionReconstructor>();
	force_add<nn4Reconstructor>();
	force_add<nn4_rectReconstructor>();
	force_add<nnSSNR_Reconstructor>();
	force_add<nn4_ctfReconstructor>();
	force_add<nn4_ctfwReconstructor>();
	force_add<nn4_ctf_rectReconstructor>();
	force_add<nnSSNR_ctfReconstructor>();
//	force_add<XYZReconstructor>();
}

class ctf_store_real
{
public:

    static void init( int winsize, const Ctf* ctf ) {
		Dict params = ctf->to_dict();

		m_winsize = winsize;

		m_voltage = params["voltage"];
		m_pixel   = params["apix"];
		m_cs      = params["cs"];
		m_ampcont = params["ampcont"];
		m_bfactor = params["bfactor"];
		m_defocus = params["defocus"];
		m_dza     = params["dfdiff"];
		m_azz     = params["dfang"];
		m_winsize2= m_winsize*m_winsize;
		m_vecsize = m_winsize2/4;
    }

    // static float get_ctf( int r2 ,int i, int j) {
	// 	float  ak = std::sqrt( r2/float(m_winsize2) )/m_pixel;
	// 	if(m_dza == 0.0f)  return Util::tf( m_defocus, ak, m_voltage, m_cs, m_ampcont, m_bfactor, 1);
	// 	else {
	// 		float az = atan2(float(j), float(i));
	// 		float dzz = m_defocus - m_dza/2.0f*sin(2*(az+m_azz*M_PI/180.0f));
	// 		return Util::tf( dzz, ak, m_voltage, m_cs, m_ampcont, m_bfactor, 1);
	// 	}
	// }
	
    static EMData* get_ctf_real() {
		return Util::ctf_img_real(m_winsize, m_winsize, 1, m_defocus, m_pixel, m_voltage, m_cs, m_ampcont, m_bfactor, m_dza, m_azz, 1);
	}

private:

	static int m_winsize, m_winsize2, m_vecsize;
	static float m_cs;
	static float m_voltage;
	static float m_pixel;
	static float m_ampcont;
	static float m_bfactor;
	static float m_defocus;
	static float m_dza;
	static float m_azz;
};

int ctf_store_real::m_winsize, ctf_store_real::m_winsize2, ctf_store_real::m_vecsize;

float ctf_store_real::m_cs, ctf_store_real::m_voltage, ctf_store_real::m_pixel;
float ctf_store_real::m_ampcont, ctf_store_real::m_bfactor;
float ctf_store_real::m_defocus, ctf_store_real::m_dza, ctf_store_real::m_azz;


void FourierReconstructorSimple2D::setup()
{
	nx = params.set_default("nx",0);

	if ( nx < 0 ) throw InvalidValueException(nx, "nx must be positive");

	bool is_fftodd = (nx % 2 == 1);

	ny = nx;
	nx += 2-is_fftodd;

	image = new EMData();
	image->set_size(nx, ny);
	image->set_complex(true);
	image->set_fftodd(is_fftodd);
	image->set_ri(true);

	tmp_data = new EMData();
	tmp_data->set_size(nx/2, nx);
}

int FourierReconstructorSimple2D::insert_slice(const EMData* const slice, const Transform & euler, const float)
{

	// Are these exceptions really necessary? (d.woolford)
	if (!slice) throw NullPointerException("EMData pointer (input image) is NULL");

	if ( slice->get_ndim() != 1 ) throw ImageDimensionException("Image dimension must be 1");

	// I could also test to make sure the image is the right dimensions...
	if (slice->is_complex()) throw ImageFormatException("The image is complex, expecting real");

	EMData* working_slice = slice->process("xform.phaseorigin.tocorner");

	// Fourier transform the slice
	working_slice->do_fft_inplace();

	float* rdata = image->get_data();
	float* norm = tmp_data->get_data();
	float* dat = working_slice->get_data();

	float g[4];
	int offset[4];
	float dt[2];
	offset[0] = 0; offset[1] = 2; offset[2] = nx; offset[3] = nx+2;

	float alt = -((float)(euler.get_rotation("2d"))["alpha"])*M_PI/180.0f;
	for (int x = 0; x < working_slice->get_xsize() / 2; x++) {

		float rx = (float) x;

		float xx = rx*cos(alt);
		float yy = rx*sin(alt);
		float cc = 1.0;

		if (xx < 0) {
			xx = -xx;
			yy = -yy;
			cc = -1.0;
		}

		yy += ny / 2;


		dt[0] = dat[2*x];
		dt[1] = cc * dat[2*x+1];

		// PHIL IS INTERESTED FROM HERE DOWN
		int x0 = (int) floor(xx);
		int y0 = (int) floor(yy);

		int i = 2*x0 + y0*nx;

		float dx = xx - x0;
		float dy = yy - y0;

		g[0] = Util::agauss(1, dx, dy, 0, EMConsts::I2G);
		g[1] = Util::agauss(1, 1 - dx, dy, 0, EMConsts::I2G);
		g[2] = Util::agauss(1, dx, 1 - dy, 0, EMConsts::I2G);
		g[3] = Util::agauss(1, 1 - dx, 1 - dy, 0, EMConsts::I2G);

		// At the extreme we can only do some much...
		if ( x0 == nx-2 ) {
			int k = i + offset[0];
			rdata[k] += g[0] * dt[0];
			rdata[k + 1] += g[0] * dt[1];
			norm[k/2] += g[0];

			k = i + offset[2];
			rdata[k] += g[2] * dt[0];
			rdata[k + 1] += g[2] * dt[1];
			norm[k/2] += g[2];
			continue;

		}
		// capture and accommodate for periodic boundary conditions in the x direction
		if ( x0 > nx-2 ) {
			int dif = x0 - (nx-2);
			x0 -= dif;
		}
		// At the extreme we can only do some much...
		if ( y0 == ny -1 ) {
			int k = i + offset[0];
			rdata[k] += g[0] * dt[0];
			rdata[k + 1] += g[0] * dt[1];
			norm[k/2] += g[0];

			k = i + offset[1];
			rdata[k] += g[1] * dt[0];
			rdata[k + 1] += g[1] * dt[1];
			norm[k/2] += g[1];
			continue;
		}
		// capture and accommodate for periodic boundary conditions in the y direction
		if ( y0 > ny-1) {
			int dif = y0 - (ny-1);
			y0 -= dif;
		}

		if (x0 >= nx - 2 || y0 >= ny - 1) continue;




		for (int j = 0; j < 4; j++)
		{
			int k = i + offset[j];
			rdata[k] += g[j] * dt[0];
			rdata[k + 1] += g[j] * dt[1];
			norm[k/2] += g[j];

		}
	}

	return 0;

}

EMData *FourierReconstructorSimple2D::finish(bool)
{
	normalize_threed();

	image->process_inplace("xform.fourierorigin.tocorner");
	image->do_ift_inplace();
	image->depad();
	image->process_inplace("xform.phaseorigin.tocenter");

	EMData *ret = image;
	image = 0;
	return  ret;
}

void ReconstructorVolumeData::normalize_threed(const bool sqrt_damp,const bool wiener)
// normalizes the 3-D Fourier volume. Also imposes appropriate complex conjugate relationships
{
	float* norm = tmp_data->get_data();
	float* rdata = image->get_data();
	
	size_t nnx=tmp_data->get_xsize();
	size_t nnxy=tmp_data->get_ysize()*nnx;
	
//	printf("%d %d %d %d %d %d\n",subnx,subny,subnz,image->get_xsize(),image->get_ysize(),image->get_zsize());

	// FIXME should throw a sensible error
	if ( 0 == norm ) throw NullPointerException("The normalization volume was null!");
	if ( 0 == rdata ) throw NullPointerException("The complex reconstruction volume was null!");

	// add_complex_at handles complex conjugate addition, but the normalization image doesn't
	// take this into account, so we need to sum the appropriate values
	// only works for whole volumes!
	if (subx0==0 && subnx==nx && suby0==0 && subny==ny && subz0==0 && subnz==nz) {
//		printf("cc gain correction\n");
		for (int z=-nz/2; z<nz/2; z++) {
			for (int y=0; y<=ny/2; y++) {
				if (z<=0 && (y==0||y==ny/2)) continue;
				// x=0 plane
				size_t idx1=(y==0?0:ny-y)*nnx+(z<=0?-z:nz-z)*nnxy;
				size_t idx2=y*nnx+(z<0?nz+z:z)*nnxy;
				if (idx1==idx2) continue;
				float snorm=norm[idx1]+norm[idx2];
				norm[idx1]=norm[idx2]=snorm;
				
				// This is the x=nx-1 plane
				idx1+=nnx-1;
				idx2+=nnx-1;
				snorm=norm[idx1]+norm[idx2];
				norm[idx1]=norm[idx2]=snorm;
			}
		}
		// special cases not handled elsewhere
		norm[0     +   0*nnx+nz/2*nnxy]*=2;
		norm[0     +ny/2*nnx+   0*nnxy]*=2;
		norm[0     +ny/2*nnx+nz/2*nnxy]*=2;
		norm[nx/2-1+   0*nnx+nz/2*nnxy]*=2;
		norm[nx/2-1+ny/2*nnx+   0*nnxy]*=2;
	}
//	else printf("Subregion, no CC plane correction\n");
	
	// The math is a bit tricky to explain. Wiener filter is basically SNR/(1+SNR)
	// In this case, data have already been effectively multiplied by SNR (one image at a time),
	// so without Wiener filter, we just divide by total SNR. With Wiener filter, we just add
	// 1.0 to total SNR, and we're done		--steve
	float wfilt=0.0;
	if (wiener) wfilt=1.0;		
		
	// actual normalization
	for (size_t i = 0; i < (size_t)subnx * subny * subnz; i += 2) {
		float d = norm[i/2];
		if (sqrt_damp) d*=sqrt(d);
		if (d == 0) {
			rdata[i] = 0;
			rdata[i + 1] = 0;
		}
		else {
//			rdata[i]=1.0/d;
//			rdata[i+1]=0.0;		// for debugging only
			rdata[i] /= d+wfilt;
			rdata[i + 1] /= d+wfilt;
		}
	}

	
	// This task should now be handled by use of add_complex_at, which adds both values when appropriate
	
	// enforce complex conj, only works on subvolumes if the complete conjugate plane is in the volume
// 	if (subx0==0 && subnx>1 && subny==ny && subnz==nz) {
// 		for (int z=0; z<=nz/2; z++) {
// 			for (int y=1; y<=ny; y++) {
// 				if (y==0 && z==0) continue;
// 				// x=0
// 				size_t i=(size_t)(y%ny)*subnx+(size_t)(z%nz)*subnx*subny;
// 				size_t i2=(size_t)(ny-y)*subnx+(size_t)((nz-z)%nz)*subnx*subny;
// 				float ar=(rdata[i]+rdata[i2])/2.0f;
// 				float ai=(rdata[i+1]-rdata[i2+1])/2.0f;
// 				rdata[i]=ar;
// 				rdata[i2]=ar;
// 				rdata[i+1]=ai;
// 				rdata[i2+1]=-ai;
// 			}
// 		}
// 	}

// 	if (subx0+subnx==nx && subnx>1 && subny==ny && subnz==nz) {
// 		for (int z=0; z<=nz/2; z++) {
// 			for (int y=1; y<=ny; y++) {
// 				if (y==0 && z==0) continue;
// 				// x=0
// 				size_t i=(size_t)(y%ny)*subnx+(size_t)(z%nz)*subnx*subny+subnx-2;
// 				size_t i2=(size_t)(ny-y)*subnx+(size_t)((nz-z)%nz)*subnx*subny+subnx-2;
// 				float ar=(rdata[i]+rdata[i2])/2.0f;
// 				float ai=(rdata[i+1]-rdata[i2+1])/2.0f;
// 				rdata[i]=ar;
// 				rdata[i2]=ar;
// 				rdata[i+1]=ai;
// 				rdata[i2+1]=-ai;
// 			}
// 		}
// 	}
}

void FourierPlaneReconstructor::load_default_settings() {}
void FourierPlaneReconstructor::free_memory() {}
void FourierPlaneReconstructor::load_inserter() {}
void FourierPlaneReconstructor::setup() {}
void FourierPlaneReconstructor::setup_seed(EMData* seed,float seed_weight) {}
void FourierPlaneReconstructor::clear() {}
EMData* FourierPlaneReconstructor::preprocess_slice( const EMData* const slice,  const Transform& t ) { return NULL; }
int FourierPlaneReconstructor::insert_slice(const EMData* const input_slice, const Transform & arg, const float weight) { return 0; }
void FourierPlaneReconstructor::do_insert_slice_work(const EMData* const input_slice, const Transform & arg,const float weight) {}
int FourierPlaneReconstructor::determine_slice_agreement(EMData*  input_slice, const Transform & arg, const float weight,bool sub) { return 0; }
void FourierPlaneReconstructor::do_compare_slice_work(EMData* input_slice, const Transform & arg,float weight) {}
bool FourierPlaneReconstructor::pixel_at(const float& xx, const float& yy, const float& zz, float *dt) { return false; }
EMData *FourierPlaneReconstructor::finish(bool doift) { return NULL; }

void FourierReconstructor::load_default_settings()
{
	inserter=0;
	image=0;
	tmp_data=0;
}

void FourierReconstructor::free_memory()
{
	if (image) { delete image; image=0; }
	if (tmp_data) { delete tmp_data; tmp_data=0; }
	if ( inserter != 0 )
	{
		delete inserter;
		inserter = 0;
	}
}

#include <sstream>

void FourierReconstructor::load_inserter()
{
//	ss
//	string mode = (string)params["mode"];
	Dict parms;
	parms["data"] = image;
	parms["norm"] = tmp_data->get_data();
	// These aren't necessary because we deal with them before calling the inserter
// 	parms["subnx"] = nx;
// 	parms["subny"] = ny;
// 	parms["subnx"] = nz;
// 	parms["subx0"] = x0;
// 	parms["suby0"] = y0;
// 	parms["subz0"] = z0;

	if ( inserter != 0 )
	{
		delete inserter;
	}

	inserter = Factory<FourierPixelInserter3D>::get((string)params["mode"], parms);
	inserter->init();
}

void FourierReconstructor::setup()
{
	// default setting behavior - does not override if the parameter is already set
	params.set_default("mode","gauss_2");
	params.set_default("verbose",(int)0);
	
	vector<int> size=params["size"];

	nx = size[0];
	ny = size[1];
	nz = size[2];
	nx2=nx/2-1;
	ny2=ny/2;
	nz2=nz/2;

#ifdef RECONDEBUG
	for (int i=0; i<125; i++) ddata[i]=dnorm[i]=0.0;
#endif

	// Adjust nx if for Fourier transform even odd issues
	bool is_fftodd = (nx % 2 == 1);
	// The Fourier transform requires one extra pixel in the x direction,
	// which is two spaces in memory, one each for the complex and the
	// real components
	nx += 2-is_fftodd;

	if (params.has_key("subvolume")) {
		vector<int> sub=params["subvolume"];
		subx0=sub[0];
		suby0=sub[1];
		subz0=sub[2];
		subnx=sub[3];
		subny=sub[4];
		subnz=sub[5];

		if (subx0<0 || suby0<0 || subz0<0 || subx0+subnx>nx || suby0+subny>ny || subz0+subnz>nz)
			throw ImageDimensionException("The subvolume cannot extend outside the reconstructed volume");

	}
	else {
		subx0=suby0=subz0=0;
		subnx=nx;
		subny=ny;
		subnz=nz;
	}


	// Odd dimension support is here atm, but not above.
	if (image) delete image;
	image = new EMData();
	image->set_size(subnx, subny, subnz);
	image->set_complex(true);
	image->set_fftodd(is_fftodd);
	image->set_ri(true);
//	printf("%d %d %d\n\n",subnx,subny,subnz);
	image->to_zero();

	if (params.has_key("subvolume")) {
		image->set_attr("subvolume_x0",subx0);
		image->set_attr("subvolume_y0",suby0);
		image->set_attr("subvolume_z0",subz0);
		image->set_attr("subvolume_full_nx",nx);
		image->set_attr("subvolume_full_ny",ny);
		image->set_attr("subvolume_full_nz",nz);
	}
	
	if (tmp_data) delete tmp_data;
	tmp_data = new EMData();
	tmp_data->set_size(subnx/2, subny, subnz);
	tmp_data->to_zero();
	tmp_data->update();

	load_inserter();

	if ( (bool) params["verbose"] )
	{
		cout << "3D Fourier dimensions are " << nx << " " << ny << " " << nz << endl;
		cout << "3D Fourier subvolume is " << subnx << " " << subny << " " << subnz << endl;
		printf ("You will require approximately %1.3g GB of memory to reconstruct this volume\n",((float)subnx*subny*subnz*sizeof(float)*1.5)/1000000000.0);
	}
}

void FourierReconstructor::setup_seed(EMData* seed,float seed_weight) {
	// default setting behavior - does not override if the parameter is already set
	params.set_default("mode","gauss_2");

	vector<int> size=params["size"];

	nx = size[0];
	ny = size[1];
	nz = size[2];
	nx2=nx/2-1;
	ny2=ny/2;
	nz2=nz/2;


	// Adjust nx if for Fourier transform even odd issues
	bool is_fftodd = (nx % 2 == 1);
	// The Fourier transform requires one extra pixel in the x direction,
	// which is two spaces in memory, one each for the complex and the
	// real components
	nx += 2-is_fftodd;

	if (params.has_key("subvolume")) {
		vector<int> sub=params["subvolume"];
		subx0=sub[0];
		suby0=sub[1];
		subz0=sub[2];
		subnx=sub[3];
		subny=sub[4];
		subnz=sub[5];

		if (subx0<0 || suby0<0 || subz0<0 || subx0+subnx>nx || suby0+subny>ny || subz0+subnz>nz)
			throw ImageDimensionException("The subvolume cannot extend outside the reconstructed volume");

	}
	else {
		subx0=suby0=subz0=0;
		subnx=nx;
		subny=ny;
		subnz=nz;
	}

	if (seed->get_xsize()!=subnx || seed->get_ysize()!=subny || seed->get_zsize()!=subnz || !seed->is_complex())
		throw ImageDimensionException("The dimensions of the seed volume do not match the reconstruction size");

	// Odd dimension support is here atm, but not above.
	image = seed;
	if (params.has_key("subvolume")) {
		image->set_attr("subvolume_x0",subx0);
		image->set_attr("subvolume_y0",suby0);
		image->set_attr("subvolume_z0",subz0);
		image->set_attr("subvolume_full_nx",nx);
		image->set_attr("subvolume_full_ny",ny);
		image->set_attr("subvolume_full_nz",nz);
	}

	if (tmp_data) delete tmp_data;
	tmp_data = new EMData();
	tmp_data->set_size(subnx/2, subny, subnz);
	tmp_data->to_value(seed_weight);

	load_inserter();

	if ( (bool) params["quiet"] == false )
	{
		cout << "Seeded direct Fourier inversion";
		cout << "3D Fourier dimensions are " << nx << " " << ny << " " << nz << endl;
		cout << "3D Fourier subvolume is " << subnx << " " << subny << " " << subnz << endl;
		cout << "You will require approximately " << setprecision(3) << (subnx*subny*subnz*sizeof(float)*1.5)/1000000000.0 << "GB of memory to reconstruct this volume" << endl;
	}
}

void FourierReconstructor::clear()
{
	bool zeroimage = true;
	bool zerotmpimg = true;
	
#ifdef EMAN2_USING_CUDA
	if(EMData::usecuda == 1) {
		if(image->getcudarwdata()) {
			to_zero_cuda(image->getcudarwdata(),image->get_xsize(),image->get_ysize(),image->get_zsize());
			zeroimage = false;
		}
		if(tmp_data->getcudarwdata()) {
			//to_zero_cuda(tmp_data->getcudarwdata(),image->get_xsize(),image->get_ysize(),image->get_zsize());
			zerotmpimg = false;
		}
	}
#endif

	if(zeroimage) image->to_zero();
	if(zerotmpimg) tmp_data->to_zero();
	
}

EMData* FourierReconstructor::preprocess_slice( const EMData* const slice,  const Transform& t )
{
#ifdef EMAN2_USING_CUDA
	if(EMData::usecuda == 1) {
		if(!slice->getcudarwdata()) slice->copy_to_cuda(); //copy slice to cuda using the const version
	}
#endif
	// Shift the image pixels so the real space origin is now located at the phase origin (at the bottom left of the image)
	EMData* return_slice = 0;
	Transform tmp(t);
	tmp.set_rotation(Dict("type","eman")); // resets the rotation to 0 implicitly, this way we only do 2d translation,scaling and mirroring

	if (tmp.is_identity()) return_slice=slice->copy();
	else return_slice = slice->process("xform",Dict("transform",&tmp));

	return_slice->process_inplace("xform.phaseorigin.tocorner");

//	printf("@@@ %d\n",(int)return_slice->get_attr("nx"));
	// Fourier transform the slice
	
#ifdef EMAN2_USING_CUDA
	if(EMData::usecuda == 1 && return_slice->getcudarwdata()) {
		return_slice->do_fft_inplace_cuda(); //a CUDA FFT inplace is quite slow as there is a lot of mem copying.
	}else{
		return_slice->do_fft_inplace();
	}
#else
	return_slice->do_fft_inplace();
#endif

//	printf("%d\n",(int)return_slice->get_attr("nx"));

	return_slice->mult((float)sqrt(1.0f/(return_slice->get_ysize())*return_slice->get_xsize()));

	// Shift the Fourier transform so that it's origin is in the center (bottom) of the image.
//	return_slice->process_inplace("xform.fourierorigin.tocenter");

	return_slice->set_attr("reconstruct_preproc",(int)1);
	return return_slice;
}


int FourierReconstructor::insert_slice(const EMData* const input_slice, const Transform & arg, const float weight)
{
	// Are these exceptions really necessary? (d.woolford)
	if (!input_slice) throw NullPointerException("EMData pointer (input image) is NULL");

#ifdef EMAN2_USING_CUDA
	if(EMData::usecuda == 1) {
		if(!input_slice->getcudarwdata()) input_slice->copy_to_cuda(); //copy slice to cuda using the const version
	}
#endif

	Transform * rotation;
/*	if ( input_slice->has_attr("xform.projection") ) {
		rotation = (Transform*) (input_slice->get_attr("xform.projection")); // assignment operator
	} else {*/
	rotation = new Transform(arg); // assignment operator
// 	}

	EMData *slice;
	if (input_slice->get_attr_default("reconstruct_preproc",(int) 0)) slice=input_slice->copy();
	else slice = preprocess_slice( input_slice, *rotation);


	// We must use only the rotational component of the transform, scaling, translation and mirroring
	// are not implemented in Fourier space, but are in preprocess_slice
	rotation->set_scale(1.0);
	rotation->set_mirror(false);
	rotation->set_trans(0,0,0);

	// Finally to the pixel wise slice insertion
	//slice->copy_to_cuda();
//	EMData *s2=slice->do_ift();
//	s2->write_image("is.hdf",-1);
	do_insert_slice_work(slice, *rotation, weight);
	
	delete rotation; rotation=0;
	delete slice;

// 	image->update();
	return 0;
}

void FourierReconstructor::do_insert_slice_work(const EMData* const input_slice, const Transform & arg,const float weight)
{
	// Reload the inserter if the mode has changed
// 	string mode = (string) params["mode"];
// 	if ( mode != inserter->get_name() )	load_inserter();

// 	int y_in = input_slice->get_ysize();
// 	int x_in = input_slice->get_xsize();
// 	// Adjust the dimensions to account for odd and even ffts
// 	if (input_slice->is_fftodd()) x_in -= 1;
// 	else x_in -= 2;

	vector<Transform> syms = Symmetry3D::get_symmetries((string)params["sym"]);

	float inx=(float)(input_slice->get_xsize());		// x/y dimensions of the input image
	float iny=(float)(input_slice->get_ysize());

#ifdef EMAN2_USING_CUDA
	if(EMData::usecuda == 1) {
		if(!image->getcudarwdata()){
			image->copy_to_cuda();
			tmp_data->copy_to_cuda();
		}
		float * m = new float[12];
		for ( vector<Transform>::const_iterator it = syms.begin(); it != syms.end(); ++it ) {
			Transform t3d = arg*(*it);
			t3d.copy_matrix_into_array(m);
			//cout << "using CUDA " << image->getcudarwdata() << endl;
			insert_slice_cuda(m,input_slice->getcudarwdata(),image->getcudarwdata(),tmp_data->getcudarwdata(),inx,iny,image->get_xsize(),image->get_ysize(),image->get_zsize(), weight);
		}
		delete m;
		return;
	}
#endif
	for ( vector<Transform>::const_iterator it = syms.begin(); it != syms.end(); ++it ) {
		Transform t3d = arg*(*it);
		for (int y = -iny/2; y < iny/2; y++) {
			for (int x = 0; x < inx/2; x++) {

				float rx = (float) x/(inx-2.0f);	// coords relative to Nyquist=.5
				float ry = (float) y/iny;

				Vec3f coord(rx,ry,0);
				coord = coord*t3d; // transpose multiplication
				float xx = coord[0]; // transformed coordinates in terms of Nyquist
				float yy = coord[1];
				float zz = coord[2];

				// Map back to real pixel coordinates in output volume
				xx=xx*(nx-2);
				yy=yy*ny;
				zz=zz*nz;

// 				if (x==10 && y==0) printf("10,0 -> %1.2f,%1.2f,%1.2f\t(%5.2f %5.2f %5.2f   %5.2f %5.2f %5.2f   %5.2f %5.2f %5.2f) %1.0f %d\n",
// 					xx,yy,zz,t3d.at(0,0),t3d.at(0,1),t3d.at(0,2),t3d.at(1,0),t3d.at(1,1),t3d.at(1,2),t3d.at(2,0),t3d.at(2,1),t3d.at(2,2),inx,nx);
// 				if (x==0 && y==10 FourierReconstructor:) printf("0,10 -> %1.2f,%1.2f,%1.2f\t(%5.2f %5.2f %5.2f   %5.2f %5.2f %5.2f   %5.2f %5.2f %5.2f)\n",
// 					xx,yy,zz,t3d.at(0,0),t3d.at(0,1),t3d.at(0,2),t3d.at(1,0),t3d.at(1,1),t3d.at(1,2),t3d.at(2,0),t3d.at(2,1),t3d.at(2,2));

				//printf("%3.1f %3.1f %3.1f\t %1.4f %1.4f\t%1.4f\n",xx,yy,zz,input_slice->get_complex_at(x,y).real(),input_slice->get_complex_at(x,y).imag(),weight);
//				if (floor(xx)==45 && floor(yy)==45 &&floor(zz)==0) printf("%d. 45 45 0\t %d %d\t %1.4f %1.4f\t%1.4f\n",(int)input_slice->get_attr("n"),x,y,input_slice->get_complex_at(x,y).real(),input_slice->get_complex_at(x,y).imag(),weight);
//				if (floor(xx)==21 && floor(yy)==21 &&floor(zz)==0) printf("%d. 21 21 0\t %d %d\t %1.4f %1.4f\t%1.4f\n",(int)input_slice->get_attr("n"),x,y,input_slice->get_complex_at(x,y).real(),input_slice->get_complex_at(x,y).imag(),weight);
				inserter->insert_pixel(xx,yy,zz,input_slice->get_complex_at(x,y),weight);
			}
		}
	}
}

int FourierReconstructor::determine_slice_agreement(EMData*  input_slice, const Transform & arg, const float weight,bool sub)
{
	// Are these exceptions really necessary? (d.woolford)
	if (!input_slice) throw NullPointerException("EMData pointer (input image) is NULL");

#ifdef EMAN2_USING_CUDA
	if(EMData::usecuda == 1) {
		if(!input_slice->getcudarwdata()) input_slice->copy_to_cuda(); //copy slice to cuda using the const version
	}
#endif

	Transform * rotation;
	rotation = new Transform(arg); // assignment operator

 	EMData *slice;
 	if (input_slice->get_attr_default("reconstruct_preproc",(int) 0)) slice=input_slice->copy();
 	else slice = preprocess_slice( input_slice, *rotation);


	// We must use only the rotational component of the transform, scaling, translation and mirroring
	// are not implemented in Fourier space, but are in preprocess_slice
	rotation->set_scale(1.0);
	rotation->set_mirror(false);
	rotation->set_trans(0,0,0);
	if (sub) do_insert_slice_work(slice, *rotation, -weight);
	// Remove the current slice first (not threadsafe, but otherwise performance would be awful)
	
	// Compare
	do_compare_slice_work(slice, *rotation,weight);

	input_slice->set_attr("reconstruct_norm",slice->get_attr("reconstruct_norm"));
	input_slice->set_attr("reconstruct_absqual",slice->get_attr("reconstruct_absqual"));
//	input_slice->set_attr("reconstruct_qual",slice->get_attr("reconstruct_qual"));
	input_slice->set_attr("reconstruct_weight",slice->get_attr("reconstruct_weight"));

	// Now put the slice back
	if (sub) do_insert_slice_work(slice, *rotation, weight);

	delete rotation;
	delete slice;

// 	image->update();
	return 0;

}

void FourierReconstructor::do_compare_slice_work(EMData* input_slice, const Transform & arg,float weight)
{

	float dt[3];	// This stores the complex and weight from the volume
	float dt2[2];	// This stores the local image complex
	float *dat = input_slice->get_data();
	vector<Transform> syms = Symmetry3D::get_symmetries((string)params["sym"]);

	float inx=(float)(input_slice->get_xsize());		// x/y dimensions of the input image
	float iny=(float)(input_slice->get_ysize());

	double dot=0;		// summed pixel*weight dot product
	double vweight=0;		// sum of weights
	double power=0;		// sum of inten*weight from volume
	double power2=0;		// sum of inten*weight from image
	bool use_cpu = true;
	
#ifdef EMAN2_USING_CUDA
	if(EMData::usecuda == 1) {
		if(!input_slice->getcudarwdata()) input_slice->copy_to_cuda();
		if(!image->getcudarwdata()){
			image->copy_to_cuda();
			tmp_data->copy_to_cuda();
		}
		for ( vector<Transform>::const_iterator it = syms.begin(); it != syms.end(); ++it ) {
			Transform t3d = arg*(*it);
			float * m = new float[12];
			t3d.copy_matrix_into_array(m);
			float4 stats = determine_slice_agreement_cuda(m,input_slice->getcudarwdata(),image->getcudarwdata(),tmp_data->getcudarwdata(),inx,iny,image->get_xsize(),image->get_ysize(),image->get_zsize(), weight);
			dot = stats.x;
			vweight = stats.y;
			power = stats.z;
			power2 = stats.w;
			//cout << "CUDA stats " << stats.x << " " << stats.y << " " << stats.z << " " << stats.w << endl;
			use_cpu = false;
		}
	}
#endif
	if(use_cpu) {
		for ( vector<Transform>::const_iterator it = syms.begin(); it != syms.end(); ++it ) {
			Transform t3d = arg*(*it);
			for (int y = -iny/2; y < iny/2; y++) {
				for (int x = 0; x <=  inx/2; x++) {
					if (x==0 && y==0) continue;		// We don't want to use the Fourier origin

					float rx = (float) x/(inx-2);	// coords relative to Nyquist=.5
					float ry = (float) y/iny;

// 					if ((rx * rx + Util::square(ry - max_input_dim "xform.projection"/ 2)) > rl)
// 					continue;

					Vec3f coord(rx,ry,0);
					coord = coord*t3d; // transpose multiplication
					float xx = coord[0]; // transformed coordinates in terms of Nyquist
					float yy = coord[1];
					float zz = coord[2];


					if (fabs(xx)>0.5 || fabs(yy)>=0.5 || fabs(zz)>=0.5) continue;

					// Map back to actual pixel coordinates in output volume
					xx=xx*(nx-2);
					yy=yy*ny;
					zz=zz*nz;


					int idx = (int)(x * 2 + inx*(y<0?iny+y:y));
					dt2[0] = dat[idx];
					dt2[1] = dat[idx+1];

					// value returned indirectly in dt
					if (!pixel_at(xx,yy,zz,dt) || dt[2]==0) continue;

//					printf("%f\t%f\t%f\t%f\t%f\n",dt[0],dt[1],dt[2],dt2[0],dt2[1]);
					dot+=(dt[0]*dt2[0]+dt[1]*dt2[1])*dt[2];
					vweight+=dt[2];
					power+=(dt[0]*dt[0]+dt[1]*dt[1])*dt[2];
					power2+=(dt2[0]*dt2[0]+dt2[1]*dt2[1])*dt[2];
				}
			}
			//cout << dot << " " << vweight << " " << power << " " << power2 << endl;
		}
	}
	
	dot/=sqrt(power*power2);		// normalize the dot product
//	input_slice->set_attr("reconstruct_norm",(float)(power2<=0?1.0:sqrt(power/power2)/(inx*iny)));
	input_slice->set_attr("reconstruct_norm",(float)(power2<=0?1.0:sqrt(power/power2)));
	input_slice->set_attr("reconstruct_absqual",(float)dot);
	float rw=weight<=0?1.0f:1.0f/weight;
	input_slice->set_attr("reconstruct_qual",(float)(dot*rw/((rw-1.0)*dot+1.0)));	// here weight is a proxy for SNR
	input_slice->set_attr("reconstruct_weight",(float)vweight/(float)(subnx*subny*subnz));
//	printf("** %g\t%g\t%g\t%g ##\n",dot,vweight,power,power2);
	//printf("** %f %f %f ##\n",(float)(power2<=0?1.0:sqrt(power/power2)/(inx*iny)),(float)dot,(float)(dot*weight/((weight-1.0)*dot+1.0)));
}

bool FourierReconstructor::pixel_at(const float& xx, const float& yy, const float& zz, float *dt)
{
	int x0 = (int) floor(xx);
	int y0 = (int) floor(yy);
	int z0 = (int) floor(zz);
	
	float *rdata=image->get_data();
	float *norm=tmp_data->get_data();
	float normsum=0,normsum2=0;

	dt[0]=dt[1]=dt[2]=0.0;

	if (nx==subnx) {			// normal full reconstruction
		if (x0<-nx2-1 || y0<-ny2-1 || z0<-nz2-1 || x0>nx2 || y0>ny2 || z0>nz2 ) return false;

		// no error checking on add_complex_fast, so we need to be careful here
		int x1=x0+1;
		int y1=y0+1;
		int z1=z0+1;
		if (x0<-nx2) x0=-nx2;
		if (x1>nx2) x1=nx2;
		if (y0<-ny2) y0=-ny2;
		if (y1>ny2) y1=ny2;
		if (z0<-nz2) z0=-nz2;
		if (z1>nz2) z1=nz2;
		
		size_t idx=0;
		float r, gg;
		for (int k = z0 ; k <= z1; k++) {
			for (int j = y0 ; j <= y1; j++) {
				for (int i = x0; i <= x1; i ++) {
					r = Util::hypot3sq((float) i - xx, j - yy, k - zz);
					idx=image->get_complex_index_fast(i,j,k);
					gg = Util::fast_exp(-r / EMConsts::I2G);
					
					dt[0]+=gg*rdata[idx];
					dt[1]+=(i<0?-1.0f:1.0f)*gg*rdata[idx+1];
					dt[2]+=norm[idx/2]*gg;
					normsum2+=gg;
					normsum+=gg*norm[idx/2];				
				}
			}
		}
		if (normsum==0) return false;
		dt[0]/=normsum;
		dt[1]/=normsum;
		dt[2]/=normsum2;
//		printf("%1.2f,%1.2f,%1.2f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n",xx,yy,zz,dt[0],dt[1],dt[2],rdata[idx],rdata[idx+1]);
		return true;
	} 
	else {					// for subvolumes, not optimized yet
		size_t idx;
		float r, gg;
		for (int k = z0 ; k <= z0 + 1; k++) {
			for (int j = y0 ; j <= y0 + 1; j++) {
				for (int i = x0; i <= x0 + 1; i ++) {
					r = Util::hypot3sq((float) i - xx, j - yy, k - zz);
					idx=image->get_complex_index(i,j,k,subx0,suby0,subz0,nx,ny,nz);
					gg = Util::fast_exp(-r / EMConsts::I2G)*norm[idx/2];
					
					dt[0]+=gg*rdata[idx];
					dt[1]+=(i<0?-1.0f:1.0f)*gg*rdata[idx+1];
					dt[2]+=norm[idx/2];
					normsum+=gg;
				}
			}
		}
		
		if (normsum==0)  return false;
		return true;
	}
}


EMData *FourierReconstructor::finish(bool doift)
{
// 	float *norm = tmp_data->get_data();
// 	float *rdata = image->get_data();
#ifdef EMAN2_USING_CUDA
	if(EMData::usecuda == 1 && image->getcudarwdata()){
		cout << "copy back from CUDA" << endl;
		image->copy_from_device();
		tmp_data->copy_from_device();
	}
#endif
	
	bool sqrtnorm=params.set_default("sqrtnorm",false);
	normalize_threed(sqrtnorm);
	
// This compares single precision sum to double precision sum near the origin
#ifdef RECONDEBUG
	for (int k=0; k<5; k++) {
		for (int j=0; j<5; j++) {
			for (int i=0; i<5; i++) {
				int idx=i*2+j*10+k*50;
				ddata[idx]/=dnorm[idx];
				ddata[idx+1]/=dnorm[idx+1];
				printf("%d %d %d   %1.4lg\t%1.4g     %1.4lg\t%1.4g\n",i,j,k,ddata[idx],image->get_value_at(i*2,j,k),ddata[idx+1],image->get_value_at(i*2+1,j,k));
			}
		}
	}
#endif
	
// 	tmp_data->write_image("density.mrc");

	// we may as well delete the tmp data now... it saves memory and the calling program might
	// need memory after it gets the return volume.
	// If this delete didn't happen now, it would happen when the deconstructor was called --david
	// no longer a good idea with the new iterative scheme -- steve
// 	if ( tmp_data != 0 )
// 	{
// 		delete tmp_data;
// 		tmp_data = 0;
// 	}

/*	image->process_inplace("xform.fourierorigin.tocorner");*/

	if (doift) {
		image->do_ift_inplace();
		image->depad();
		image->process_inplace("xform.phaseorigin.tocenter");
	}
	// If the image was padded it should be the original size, as the client would expect
	//  I blocked the rest, it is almost certainly incorrect  PAP 07/31/08
	// No, it's not incorrect. You are wrong. You have the meaning of nx mixed up. DSAW 09/23/cd
	// This should now be handled in the calling program --steve 11/03/09
// 	bool is_fftodd = (nx % 2 == 1);
// 	if ( (nx-2*(!is_fftodd)) != output_x || ny != output_y || nz != output_z )
// 	{
// 		FloatPoint origin( (nx-output_x)/2, (ny-output_y)/2, (nz-output_z)/2 );
// 		FloatSize region_size( output_x, output_y, output_z);
// 		Region clip_region( origin, region_size );
// 		image->clip_inplace( clip_region );
// 	}

	// Should be an "if (verbose)" here or something
	//print_stats(quality_scores);

	image->update();
	
	if (params.has_key("savenorm") && strlen((const char *)params["savenorm"])>0) {
		if (tmp_data->get_ysize()%2==0 && tmp_data->get_zsize()%2==0) tmp_data->process_inplace("xform.fourierorigin.tocenter");
		tmp_data->write_image((const char *)params["savenorm"]);
	}

	delete tmp_data;
	tmp_data=0;
	//Since we give up the ownership of the pointer to-be-returned,	it's caller's responsibility to delete the returned image.
	//So we wrap this function with return_value_policy< manage_new_object >() in libpyReconstructor2.cpp to hand over ownership to Python.
	EMData *ret=image;
	image=0;
	
	return ret;
}

int WienerFourierReconstructor::insert_slice(const EMData* const input_slice, const Transform & arg, const float weight)
{
	// Are these exceptions really necessary? (d.woolford)
	if (!input_slice) throw NullPointerException("EMData pointer (input image) is NULL");

	Transform * rotation;
/*	if ( input_slice->has_attr("xform.projection") ) {
		rotation = (Transform*) (input_slice->get_attr("xform.projection")); // assignment operator
	} else {*/
	rotation = new Transform(arg); // assignment operator
// 	}

	if (!input_slice->has_attr("ctf_snr_total")) 
		throw NotExistingObjectException("ctf_snr_total","No SNR information present in class-average. Must use the ctf.auto or ctfw.auto averager.");

	EMData *slice;
	if (input_slice->get_attr_default("reconstruct_preproc",(int) 0)) slice=input_slice->copy();
	else slice = preprocess_slice( input_slice, *rotation);


	// We must use only the rotational component of the transform, scaling, translation and mirroring
	// are not implemented in Fourier space, but are in preprocess_slice
	rotation->set_scale(1.0);
	rotation->set_mirror(false);
	rotation->set_trans(0,0,0);

	// Finally to the pixel wise slice insertion
	do_insert_slice_work(slice, *rotation, weight);

	delete rotation; rotation=0;
	delete slice;

// 	image->update();
	return 0;
}

void WienerFourierReconstructor::do_insert_slice_work(const EMData* const input_slice, const Transform & arg,const float inweight)
{

	vector<Transform> syms = Symmetry3D::get_symmetries((string)params["sym"]);

	float inx=(float)(input_slice->get_xsize());		// x/y dimensions of the input image
	float iny=(float)(input_slice->get_ysize());
	
	int undo_wiener=(int)input_slice->get_attr_default("ctf_wiener_filtered",0);	// indicates whether we need to undo a wiener filter before insertion
//	if (undo_wiener) throw UnexpectedBehaviorException("wiener_fourier does not yet accept already Wiener filtered class-averages. Suggest using ctf.auto averager for now.");
	
	vector<float> snr=input_slice->get_attr("ctf_snr_total");
	float sub=1.0;
	if (inweight<0) sub=-1.0;
	float weight;
	
	for ( vector<Transform>::const_iterator it = syms.begin(); it != syms.end(); ++it ) {
		Transform t3d = arg*(*it);
		for (int y = -iny/2; y < iny/2; y++) {
			for (int x = 0; x <=  inx/2; x++) {

				float rx = (float) x/(inx-2.0f);	// coords relative to Nyquist=.5
				float ry = (float) y/iny;

				// This deals with the SNR weight
				float rn = (float)hypot(rx,ry);
				if (rn>=.5) continue;		// no SNR in the corners, and we're going to mask them later anyway
				rn*=snr.size()*2.0f;
				int rni=(int)floor(rn);
				if ((unsigned int)rni>=snr.size()-1) weight=snr[snr.size()-1]*sub;
				else {
					rn-=rni;
					weight=Util::linear_interpolate(snr[rni],snr[rni+1],rn);
				}
//				if (weight>500.0) printf("%f %d %d %f %f %d %f\n",weight,x,y,rx,ry,rni);
				
				Vec3f coord(rx,ry,0);
				coord = coord*t3d; // transpose multiplication
				float xx = coord[0]; // transformed coordinates in terms of Nyquist
				float yy = coord[1];
				float zz = coord[2];

				// Map back to real pixel coordinates in output volume
				xx=xx*(nx-2);
				yy=yy*ny;
				zz=zz*nz;

//				printf("%f\n",weight);
				if (undo_wiener) inserter->insert_pixel(xx,yy,zz,(input_slice->get_complex_at(x,y))*((weight+1.0f)/weight),weight*sub);
				else inserter->insert_pixel(xx,yy,zz,input_slice->get_complex_at(x,y),weight*sub);
			}
		}
	}
}

int WienerFourierReconstructor::determine_slice_agreement(EMData*  input_slice, const Transform & arg, const float weight,bool sub)
{
	// Are these exceptions really necessary? (d.woolford)
	if (!input_slice) throw NullPointerException("EMData pointer (input image) is NULL");

	Transform * rotation;
	rotation = new Transform(arg); // assignment operator

 	EMData *slice;
 	if (input_slice->get_attr_default("reconstruct_preproc",(int) 0)) slice=input_slice->copy();
 	else slice = preprocess_slice( input_slice, *rotation);


	// We must use only the rotational component of the transform, scaling, translation and mirroring
	// are not implemented in Fourier space, but are in preprocess_slice
	rotation->set_scale(1.0);
	rotation->set_mirror(false);
	rotation->set_trans(0,0,0);

//	tmp_data->write_image("dbug.hdf",0);
	
	// Remove the current slice first (not threadsafe, but otherwise performance would be awful)
	if (sub) do_insert_slice_work(slice, *rotation, -weight);

	// Compare
	do_compare_slice_work(slice, *rotation,weight);

	input_slice->set_attr("reconstruct_norm",slice->get_attr("reconstruct_norm"));
	input_slice->set_attr("reconstruct_absqual",slice->get_attr("reconstruct_absqual"));
//	input_slice->set_attr("reconstruct_qual",slice->get_attr("reconstruct_qual"));
	input_slice->set_attr("reconstruct_weight",slice->get_attr("reconstruct_weight"));

	// Now put the slice back
	if (sub) do_insert_slice_work(slice, *rotation, weight);


	delete rotation;
	delete slice;

// 	image->update();
	return 0;

}

void WienerFourierReconstructor::do_compare_slice_work(EMData* input_slice, const Transform & arg,float weight)
{

	float dt[3];	// This stores the complex and weight from the volume
	float dt2[2];	// This stores the local image complex
	float *dat = input_slice->get_data();
	vector<Transform> syms = Symmetry3D::get_symmetries((string)params["sym"]);

	float inx=(float)(input_slice->get_xsize());		// x/y dimensions of the input image
	float iny=(float)(input_slice->get_ysize());

	double dot=0;		// summed pixel*weight dot product
	double vweight=0;		// sum of weights
	double power=0;		// sum of inten*weight from volume
	double power2=0;		// sum of inten*weight from image
	for ( vector<Transform>::const_iterator it = syms.begin(); it != syms.end(); ++it ) {
		Transform t3d = arg*(*it);
		for (int y = -iny/2; y < iny/2; y++) {
			for (int x = 0; x <=  inx/2; x++) {
				if (x==0 && y==0) continue;		// We don't want to use the Fourier origin

				float rx = (float) x/(inx-2);	// coords relative to Nyquist=.5
				float ry = (float) y/iny;

// 				if ((rx * rx + Util::square(ry - max_input_dim / 2)) > rl)
// 					continue;

				Vec3f coord(rx,ry,0);
				coord = coord*t3d; // transpose multiplication
				float xx = coord[0]; // transformed coordinates in terms of Nyquist
				float yy = coord[1];
				float zz = coord[2];


				if (fabs(xx)>0.5 || fabs(yy)>=0.5 || fabs(zz)>=0.5) continue;

				// Map back to actual pixel coordinates in output volume
				xx=xx*(nx-2);
				yy=yy*ny;
				zz=zz*nz;


				int idx = (int)(x * 2 + inx*(y<0?iny+y:y));
				dt2[0] = dat[idx];
				dt2[1] = dat[idx+1];

				// value returned indirectly in dt
				if (!pixel_at(xx,yy,zz,dt) || dt[2]<=0) continue;

//				printf("%f\t%f\t%f\t%f\t%f\n",dt[0],dt[1],dt[2],dt2[0],dt2[1]);
				dot+=(dt[0]*dt2[0]+dt[1]*dt2[1])*dt[2];
				vweight+=dt[2];
				power+=(dt[0]*dt[0]+dt[1]*dt[1])*dt[2];
				power2+=(dt2[0]*dt2[0]+dt2[1]*dt2[1])*dt[2];
			}
		}
	}

	dot/=sqrt(power*power2);		// normalize the dot product
//	input_slice->set_attr("reconstruct_norm",(float)(power2<=0?1.0:sqrt(power/power2)/(inx*iny)));
	input_slice->set_attr("reconstruct_norm",(float)(power2<=0?1.0:sqrt(power/power2)));
	input_slice->set_attr("reconstruct_absqual",(float)dot);
	float rw=weight<=0?1.0f:1.0f/weight;
	input_slice->set_attr("reconstruct_qual",(float)(dot*rw/((rw-1.0)*dot+1.0)));	// here weight is a proxy for SNR
	input_slice->set_attr("reconstruct_weight",(float)vweight/(float)(subnx*subny*subnz));
//	printf("** %g\t%g\t%g\t%g ##\n",dot,vweight,power,power2);
	//printf("** %f %f %f ##\n",(float)(power2<=0?1.0:sqrt(power/power2)/(inx*iny)),(float)dot,(float)(dot*weight/((weight-1.0)*dot+1.0)));
}

bool WienerFourierReconstructor::pixel_at(const float& xx, const float& yy, const float& zz, float *dt)
{
	int x0 = (int) floor(xx);
	int y0 = (int) floor(yy);
	int z0 = (int) floor(zz);
	
	float *rdata=image->get_data();
	float *norm=tmp_data->get_data();
	float normsum=0,normsum2=0;

	dt[0]=dt[1]=dt[2]=0.0;

	if (nx==subnx) {			// normal full reconstruction
		if (x0<-nx2-1 || y0<-ny2-1 || z0<-nz2-1 || x0>nx2 || y0>ny2 || z0>nz2 ) return false;

		// no error checking on add_complex_fast, so we need to be careful here
		int x1=x0+1;
		int y1=y0+1;
		int z1=z0+1;
		if (x0<-nx2) x0=-nx2;
		if (x1>nx2) x1=nx2;
		if (y0<-ny2) y0=-ny2;
		if (y1>ny2) y1=ny2;
		if (z0<-nz2) z0=-nz2;
		if (z1>nz2) z1=nz2;
		
		size_t idx=0;
		float r, gg;
		for (int k = z0 ; k <= z1; k++) {
			for (int j = y0 ; j <= y1; j++) {
				for (int i = x0; i <= x1; i ++) {
					r = Util::hypot3sq((float) i - xx, j - yy, k - zz);
					idx=image->get_complex_index_fast(i,j,k);
					gg = Util::fast_exp(-r / EMConsts::I2G);
					
					dt[0]+=gg*rdata[idx];
					dt[1]+=(i<0?-1.0f:1.0f)*gg*rdata[idx+1];
					dt[2]+=norm[idx/2]*gg;
					normsum2+=gg;
					normsum+=gg*norm[idx/2];				
				}
			}
		}
		if (normsum==0) return false;
		dt[0]/=normsum;
		dt[1]/=normsum;
		dt[2]/=normsum2;
//		printf("%1.2f,%1.2f,%1.2f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n",xx,yy,zz,dt[0],dt[1],dt[2],rdata[idx],rdata[idx+1]);
		return true;
	} 
	else {					// for subvolumes, not optimized yet
		size_t idx;
		float r, gg;
		for (int k = z0 ; k <= z0 + 1; k++) {
			for (int j = y0 ; j <= y0 + 1; j++) {
				for (int i = x0; i <= x0 + 1; i ++) {
					r = Util::hypot3sq((float) i - xx, j - yy, k - zz);
					idx=image->get_complex_index(i,j,k,subx0,suby0,subz0,nx,ny,nz);
					gg = Util::fast_exp(-r / EMConsts::I2G)*norm[idx/2];
					
					dt[0]+=gg*rdata[idx];
					dt[1]+=(i<0?-1.0f:1.0f)*gg*rdata[idx+1];
					dt[2]+=norm[idx/2];
					normsum+=gg;				
				}
			}
		}
		
		if (normsum==0)  return false;
		return true;
	}
}


EMData *WienerFourierReconstructor::finish(bool doift)
{

	bool sqrtnorm=params.set_default("sqrtnorm",false);
	normalize_threed(sqrtnorm,true);		// true is the wiener filter

	if (doift) {
		image->do_ift_inplace();
		image->depad();
		image->process_inplace("xform.phaseorigin.tocenter");
	}

	image->update();
	
	if (params.has_key("savenorm") && strlen((const char *)params["savenorm"])>0) {
		if (tmp_data->get_ysize()%2==0 && tmp_data->get_zsize()%2==0) tmp_data->process_inplace("xform.fourierorigin.tocenter");
		tmp_data->write_image((const char *)params["savenorm"]);
	}

	delete tmp_data;
	tmp_data=0;
	EMData *ret=image;
	image=0;
	
	return ret;
}

/*
void BaldwinWoolfordReconstructor::setup()
{
	//This is a bit of a hack - but for now it suffices
	params.set_default("mode","nearest_neighbor");
	WienerFourierReconstructor::setup();
	// Set up the Baldwin Kernel
	int P = (int)((1.0+0.25)*max_input_dim+1);
	float r = (float)(max_input_dim+1)/(float)P;
	dfreq = 0.2f;
	if (W != 0) delete [] W;
	int maskwidth = params.set_default("maskwidth",2);
	W = Util::getBaldwinGridWeights(maskwidth, (float)P, r,dfreq,0.5f,0.2f);
}

EMData* BaldwinWoolfordReconstructor::finish(bool doift)
{
	tmp_data->write_image("density.mrc");
	image->process_inplace("xform.fourierorigin.tocorner");
	image->do_ift_inplace();
	image->depad();
	image->process_inplace("xform.phaseorigin.tocenter");

	if ( (bool) params.set_default("postmultiply", false) )
	{
		cout << "POST MULTIPLYING" << endl;
	// now undo the Fourier convolution with real space division
		float* d = image->get_data();
		float N = (float) image->get_xsize()/2.0f;
		N *= N;
		size_t rnx = image->get_xsize();
		size_t rny = image->get_ysize();
		size_t rnxy = rnx*rny;
		int cx = image->get_xsize()/2;
		int cy = image->get_ysize()/2;
		int cz = image->get_zsize()/2;
		size_t idx;
		for (int k = 0; k < image->get_zsize(); ++k ){
			for (int j = 0; j < image->get_ysize(); ++j ) {
				for (int i =0; i < image->get_xsize(); ++i ) {
					float xd = (float)(i-cx); xd *= xd;
					float yd = (float)(j-cy); yd *= yd;
					float zd = (float)(k-cz); zd *= zd;
					float weight = exp((xd+yd+zd)/N);
					idx = k*rnxy + j*rnx + i;
					d[idx] *=  weight;
				}
			}
		}
	}
	image->update();
	return  image;
}

#include <iomanip>

// int BaldwinWoolfordReconstructor::insert_slice_weights(const Transform& t3d)
// {
// 	bool fftodd = image->is_fftodd();
// 	int rnx = nx-2*!fftodd;
//
// 	float y_scale = 1.0, x_scale = 1.0;
//
// 	if ( ny != rnx  )
// 	{
// 		if ( rnx > ny ) y_scale = (float) rnx / (float) ny;
// 		else x_scale = (float) ny / (float) rnx;
// 	}
//
// 	int tnx = tmp_data->get_xsize();
// 	int tny = tmp_data->get_ysize();
// 	int tnz = tmp_data->get_zsize();
//
// 	vector<Transform> syms = Symmetry3D::get_symmetries((string)params["sym"]);
// 	for ( vector<Transform>::const_iterator it = syms.begin(); it != syms.end(); ++it ) {
// 		Transform n3d = t3d*(*it);
//
// 		for (int y = 0; y < tny; y++) {
// 			for (int x = 0; x < tnx; x++) {
//
// 				float rx = (float) x;
// 				float ry = (float) y;
//
// 				if ( ny != rnx )
// 				{
// 					if ( rnx > ny ) ry *= y_scale;
// 					else rx *= x_scale;
// 				}
// // 				float xx = rx * n3d[0][0] + (ry - tny/2) * n3d[1][0];
// // 				float yy = rx * n3d[0][1] + (ry - tny/2) * n3d[1][1];
// // 				float zz = rx * n3d[0][2] + (ry - tny/2) * n3d[1][2];
//
// 				Vec3f coord(rx,(ry - tny/2),0);
// 				coord = coord*n3d; // transpose multiplication
// 				float xx = coord[0];
// 				float yy = coord[1];
// 				float zz = coord[2];
//
// 				if (xx < 0 ){
// 					xx = -xx;
// 					yy = -yy;
// 					zz = -zz;
// 				}
//
// 				yy += tny/2;
// 				zz += tnz/2;
// 				insert_density_at(xx,yy,zz);
// 			}
// 		}
// 	}
//
// 	return 0;
// }

void BaldwinWoolfordReconstructor::insert_density_at(const float& x, const float& y, const float& z)
{
	int xl = Util::fast_floor(x);
	int yl = Util::fast_floor(y);
	int zl = Util::fast_floor(z);

	// w is the windowing width
	int w = params.set_default("maskwidth",2);
	float wsquared = (float) w*w;
	float dw = 1.0f/w;
	dw *= dw;

	// w minus one - this control the number of
	// pixels/voxels to the left of the main pixel
	// that will have density
	int wmox = w-1;
	int wmoy = w-1;
	int wmoz = w-1;

	// If any coordinate is incedental with a vertex, then
	// make sure there is symmetry in density accruing.
	// i.e. the window width must be equal in both directions
	if ( ((float) xl) == x ) wmox = w;
	if ( ((float) yl) == y ) wmoy = w;
	if ( ((float) zl) == z ) wmoz = w;

	float* d = tmp_data->get_data();
	int tnx = tmp_data->get_xsize();
	int tny = tmp_data->get_ysize();
	int tnz = tmp_data->get_zsize();
	size_t tnxy = tnx*tny;

	int mode = params.set_default("mode","nearest_neighbor");

	for(int k = zl-wmoz; k <= zl+w; ++k ) {
		for(int j = yl-wmoy; j <= yl+w; ++j) {
			for( int i = xl-wmox; i <= xl+w; ++i) {
				float fac = 1.0;
				int ic = i, jc = j, kc = k;

				// Fourier space is periodic, which is enforced
				// by the next 6 if statements. These if statements
				// assume that the Fourier DC components is at
				// (0,ny/2,nz/2).
				if ( i <= 0 ) {

					if ( x != 0 && i == 0 ) fac = 1.0;
					else if ( x == 0 && i < 0) continue;
// 					if (i < 0 ) ic = -i;
					if (i < 0 ) {
						continue;
						ic = -i;
						jc = tny-jc;
						kc = tnz-kc;
					}
				}
				if ( ic >= tnx ) ic = 2*tnx-ic-1;
				if ( jc < 0 ) jc = tny+jc;
				if ( jc >= tny ) jc = jc-tny;
				if ( kc < 0 ) kc = tnz+kc;
				if ( kc >= tnz ) kc = kc-tnz;
// 				if ( ic >= tnx ) continue;
// 				if ( jc < 0 ) continue;
// 				if ( jc >= tny ) continue;
// 				if ( kc < 0 ) continue;
// 				if ( kc >= tnz ) continue;
				// This shouldn't happen
				// Debug remove later
				if ( ic < 0 ) { cout << "wo 1" << endl; }
				if ( ic >= tnx  ){ cout << "wo 2" << endl; }
				if ( jc < 0 ) { cout << "wo 3" << endl; }
				if ( jc >= tny ) { cout << "wo 4" << endl; }
				if ( kc < 0 ) { cout << "wo 5" << endl; }
				if ( kc >= tnz ) { cout << "wo 6" << endl; }


				float zd = (z-(float)k);
				float yd = (y-(float)j);
				float xd = (x-(float)i);
				zd *= zd; yd *= yd; xd *= xd;
				float distsquared = xd+yd+zd;
				// We enforce a spherical kernel
				if ( mode == 1 && distsquared > wsquared ) continue;

// 				float f = fac*exp(-dw*(distsquared));
				float f = fac*exp(-2.467f*(distsquared));
				// Debug - this error should never occur.
				if ( (kc*tnxy+jc*tnx+ic) >= tnxy*tnz ) throw OutofRangeException(0,tnxy*tnz,kc*tnxy+jc*tnx+ic, "in density insertion" );
				d[kc*tnxy+jc*tnx+ic] += f;
			}
		}
	}
}

int BaldwinWoolfordReconstructor::insert_slice(const EMData* const input_slice, const Transform & t, const float weight)
{
	Transform * rotation;
	if ( input_slice->has_attr("xform.projection") ) {
		rotation = (Transform*) (input_slice->get_attr("xform.projection")); // assignment operator
	} else {
		rotation = new Transform(t); // assignment operator
	}
	Transform tmp(*rotation);
	tmp.set_rotation(Dict("type","eman")); // resets the rotation to 0 implicitly

	Vec2f trans = tmp.get_trans_2d();
	float scale = tmp.get_scale();
	bool mirror = tmp.get_mirror();
	EMData* slice = 0;
	if (trans[0] != 0 || trans[1] != 0 || scale != 1.0 ) {
		slice = input_slice->process("xform",Dict("transform",&tmp));
	} else if ( mirror == true ) {
		slice = input_slice->process("xform.flip",Dict("axis","x"));
	}
	if ( slice == 0 ) {
		slice = input_slice->process("xform.phaseorigin.tocorner");
	} else {
		slice->process_inplace("xform.phaseorigin.tocorner");
	}

	slice->do_fft_inplace();
	slice->process_inplace("xform.fourierorigin.tocenter");
	float *dat = slice->get_data();
	float dt[2];

	bool fftodd = image->is_fftodd();
	int rnx = nx-2*!fftodd;

	float y_scale = 1.0, x_scale = 1.0;

	if ( ny != rnx  )
	{
		if ( rnx > ny ) y_scale = (float) rnx / (float) ny;
		else x_scale = (float) ny / (float) rnx;
	}

	int tnx = tmp_data->get_xsize();
	int tny = tmp_data->get_ysize();
	int tnz = tmp_data->get_zsize();

	vector<Transform> syms = Symmetry3D::get_symmetries((string)params["sym"]);
// 	float weight = params.set_default("weight",1.0f);

	rotation->set_scale(1.0); rotation->set_mirror(false); rotation->set_trans(0,0,0);
	for ( vector<Transform>::const_iterator it = syms.begin(); it != syms.end(); ++it ) {
		Transform t3d = (*rotation)*(*it);

		for (int y = 0; y < tny; y++) {
			for (int x = 0; x < tnx; x++) {
				float rx = (float) x;
				float ry = (float) y;

				if ( ny != rnx )
				{
					if ( rnx > ny ) ry *= y_scale;
					else rx *= x_scale;
				}

// 				float xx = rx * n3d[0][0] + (ry - tny/2) * n3d[1][0];
// 				float yy = rx * n3d[0][1] + (ry - tny/2) * n3d[1][1];
// 				float zz = rx * n3d[0][2] + (ry - tny/2) * n3d[1][2];

				Vec3f coord(rx,(ry - tny/2),0);
				coord = coord*t3d; // transpose multiplication
				float xx = coord[0];
				float yy = coord[1];
				float zz = coord[2];


				float cc = 1;
				if (xx < 0 ){
					xx = -xx;
					yy = -yy;
					zz = -zz;
					cc = -1;
				}

				yy += tny/2;
				zz += tnz/2;

				int idx = x * 2 + y * (slice->get_xsize());
				dt[0] = dat[idx];
				dt[1] = cc * dat[idx+1];

				insert_pixel(xx,yy,zz,dt);
			}
		}
	}

	if(rotation) {delete rotation; rotation=0;}
	delete slice;

	return 0;
}

void BaldwinWoolfordReconstructor::insert_pixel(const float& x, const float& y, const float& z, const float dt[2])
{
	int xl = Util::fast_floor(x);
	int yl = Util::fast_floor(y);
	int zl = Util::fast_floor(z);

	// w is the windowing width
	int w = params.set_default("maskwidth",2);
	float wsquared = (float) w*w;
	float dw = 1.0f/w;
	dw *= dw;

	int wmox = w-1;
	int wmoy = w-1;
	int wmoz = w-1;

	// If any coordinate is incedental with a vertex, then
	// make sure there is symmetry in density accruing.
	// i.e. the window width must be equal in both directions
	if ( ((float) xl) == x ) wmox = w;
	if ( ((float) yl) == y ) wmoy = w;
	if ( ((float) zl) == z ) wmoz = w;

	float* we = tmp_data->get_data();
	int tnx = tmp_data->get_xsize();
	int tny = tmp_data->get_ysize();
	int tnz = tmp_data->get_zsize();
	int tnxy = tnx*tny;

	int rnx = 2*tnx;
	int rnxy = 2*tnxy;

	int mode = params.set_default("mode","nearest_neighbor");

	float* d = image->get_data();
	for(int k = zl-wmoz; k <= zl+w; ++k ) {
		for(int j = yl-wmoy; j <= yl+w; ++j) {
			for( int i = xl-wmox; i <= xl+w; ++i) {
				float fac = 1.0;
				int ic = i, jc = j, kc = k;

				// Fourier space is periodic, which is enforced
				// by the next 6 if statements. These if statements
				// assume that the Fourier DC component is at
				// (0,ny/2,nz/2).
				float negfac=1.0;
				if ( i <= 0 ) {
					if ( x != 0 && i == 0 ) fac = 1.0;
					else if ( x == 0 && i < 0) continue;
					if (i < 0 ) {
						continue;
						ic = -i;
						jc = tny-jc;
						kc = tnz-kc;
						negfac=-1.0;
					}
				}
				if ( ic >= tnx ) ic = 2*tnx-ic-1;
				if ( jc < 0 ) jc = tny+jc;
				if ( jc >= tny ) jc = jc-tny;
				if ( kc < 0 ) kc = tnz+kc;
				if ( kc >= tnz ) kc = kc-tnz;
// 				if ( ic >= tnx ) continue;
// 				if ( jc < 0 ) continue;
// 				if ( jc >= tny ) continue;
// 				if ( kc < 0 ) continue;
// 				if ( kc >= tnz ) continue;

				float zd = (z-(float)k);
				float yd = (y-(float)j);
				float xd = (x-(float)i);
				zd *= zd; yd *= yd; xd *= xd;
				float distsquared = xd+yd+zd;
// 				float f = fac*exp(-dw*(distsquared));
				float f = fac*exp(-2.467f*(distsquared));
				float weight = f/we[kc*tnxy+jc*tnx+ic];
				// debug - this error should never occur
				if ( (kc*rnxy+jc*rnx+2*ic+1) >= rnxy*tnz ) throw OutofRangeException(0,rnxy*tnz,kc*rnxy+jc*rnx+2*ic+1, "in pixel insertion" );
				size_t k = kc*rnxy+jc*rnx+2*ic;

				float factor, dist,residual;
				int sizeW,sizeWmid,idx;
				switch (mode) {
					case 0:
						d[k] += weight*f*dt[0];
						d[k+1] += negfac*weight*f*dt[1];
						cout << "hello" << endl;
					break;

					case 1:
						// We enforce a spherical kernel
						if ( distsquared > wsquared ) continue;

						sizeW = (int)(1+2*w/dfreq);
						sizeWmid = sizeW/2;

						dist = sqrtf(distsquared);
						idx = (int)(sizeWmid + dist/dfreq);
						if (idx >= sizeW) throw InvalidValueException(idx, "idx was greater than or equal to sizeW");
						residual = dist/dfreq - (int)(dist/dfreq);
						if ( fabs(residual) > 1) throw InvalidValueException(residual, "Residual was too big");

						factor = (W[idx]*(1.0f-residual)+W[idx+1]*residual)*weight;

						d[k] += dt[0]*factor;
						d[k+1] += dt[1]*factor;
					break;

					default:
						throw InvalidValueException(mode, "The mode was unsupported in BaldwinWoolfordReconstructor::insert_pixel");
					break;
				}
			}
		}
	}
}

// void BaldwinWoolfordReconstructor::insert_pixel(const float& x, const float& y, const float& z, const float dt[2])
// {
// 	int xl = Util::fast_floor(x);
// 	int yl = Util::fast_floor(y);
// 	int zl = Util::fast_floor(z);
//
// 	// w is the windowing width
// 	int w = params.set_default("maskwidth",2);
// 	float dw = 1.0/w;
// 	dw *= dw;
// // 	dw = 2;
// // 	cout << w << endl;
// 	// 	int w = 3;
// 	// w minus one - this control the number of
// 	// pixels/voxels to the left of the main pixel
// 	// that will have density
// 	int wmox = w-1;
// 	int wmoy = w-1;
// 	int wmoz = w-1;
//
// 	// If any coordinate is incedental with a vertex, then
// 	// make sure there is symmetry in density accruing.
// 	// i.e. the window width must be equal in both directions
// 	if ( ((float) xl) == x ) wmox = w;
// 	if ( ((float) yl) == y ) wmoy = w;
// 	if ( ((float) zl) == z ) wmoz = w;
//
// 	float* d = tmp_data->get_data();
// 	int tnx = tmp_data->get_xsize();
// 	int tny = tmp_data->get_ysize();
// 	int tnz = tmp_data->get_zsize();
// 	int tnxy = tnx*tny;
//
// 	float weight = 1.0;
// //
// 	for(int k = zl-wmoz; k <= zl+w; ++k ) {
// 		for(int j = yl-wmoy; j <= yl+w; ++j) {
// 			for( int i = xl-wmox; i <= xl+w; ++i) {
// 				float fac = 1.0;
// 				int ic = i, jc = j, kc = k;
//
// 				// Fourier space is periodic, which is enforced
// 				// by the next 6 if statements. These if statements
// 				// assume that the Fourier DC components is at
// 				// (0,ny/2,nz/2).
// 				if ( i <= 0 ) {
// 					if ( x != 0 && i == 0 ) fac = 1.0;
// 					else if ( x == 0 && i < 0) continue;
// // 					if (i < 0 ) ic = -i;
// 					if (i < 0 ) {
// 						ic = -i;
// 						jc = tny-jc;
// 						kc = tnz-kc;
// 					}
// 				}
// 				if ( ic >= tnx ) ic = 2*tnx-ic-1;
// 				if ( jc < 0 ) jc = tny+jc;
// 				if ( jc >= tny ) jc = jc-tny;
// 				if ( kc < 0 ) kc = tnz+kc;
// 				if ( kc >= tnz ) kc = kc-tnz;
// 				// This shouldn't happen
// 				// Debug remove later
// 				if ( ic < 0 ) { cout << "wo 1" << endl; }
// 				if ( ic >= tnx  ){ cout << "wo 2" << endl; }
// 				if ( jc < 0 ) { cout << "wo 3" << endl; }
// 				if ( jc >= tny ) { cout << "wo 4" << endl; }
// 				if ( kc < 0 ) { cout << "wo 5" << endl; }
// 				if ( kc >= tnz ) { cout << "wo 6" << endl; }
//
//
// 				float zd = (z-(float)k);
// 				float yd = (y-(float)j);
// 				float xd = (x-(float)i);
// 				zd *= zd; yd *= yd; xd *= xd;
// 				// Debug - this error should never occur.
// 				if ( (kc*tnxy+jc*tnx+ic) >= tnxy*tnz ) throw OutofRangeException(0,tnxy*tnz,kc*tnxy+jc*tnx+ic, "in weight determination insertion" );
// // 				float f = fac*exp(-dw*(xd+yd+zd)*0.5);
// 				float f = exp(-2.467*(xd+yd+zd));
// 				weight += f*(d[kc*tnxy+jc*tnx+ic] - f);
// 			}
// 		}
// 	}
// 	weight = 1.0/weight;
// 	int rnx = 2*tnx;
// 	int rnxy = 2*tnxy;
// 	d = image->get_data();
// 	for(int k = zl-wmoz; k <= zl+w; ++k ) {
// 		for(int j = yl-wmoy; j <= yl+w; ++j) {
// 			for( int i = xl-wmox; i <= xl+w; ++i) {
// 				float fac = 1.0;
// 				int ic = i, jc = j, kc = k;
//
// 				// Fourier space is periodic, which is enforced
// 				// by the next 6 if statements. These if statements
// 				// assume that the Fourier DC components is at
// 				// (0,ny/2,nz/2).
// 				float negfac=1.0;
// 				if ( i <= 0 ) {
// 					if ( x != 0 && i == 0 ) fac = 1.0;
// 					else if ( x == 0 && i < 0) continue;
// 					if (i < 0 ) {
// 						continue;
// 						ic = -i;
// 						jc = tny-jc;
// 						kc = tnz-kc;
// 						negfac=-1.0;
// 					}
// 				}
// 				if ( ic >= tnx ) ic = 2*tnx-ic-1;
// 				if ( jc < 0 ) jc = tny+jc;
// 				if ( jc >= tny ) jc = jc-tny;
// 				if ( kc < 0 ) kc = tnz+kc;
// 				if ( kc >= tnz ) kc = kc-tnz;
// 				// This shouldn't happen
// 				// Debug remove later
//
//
// 				float zd = (z-(float)k);
// 				float yd = (y-(float)j);
// 				float xd = (x-(float)i);
// 				zd *= zd; yd *= yd; xd *= xd;
// // 				float f = fac*exp(-dw*(xd+yd+zd));
// 				float f = exp(-4.934*(xd+yd+zd));
// 				// Debug - this error should never occur.
// 				if ( (kc*rnxy+jc*rnx+2*ic+1) >= rnxy*tnz ) throw OutofRangeException(0,rnxy*tnz,kc*rnxy+jc*rnx+2*ic+1, "in pixel insertion" );
//
// 				d[kc*rnxy+jc*rnx+2*ic] += weight*f*dt[0];
// 				d[kc*rnxy+jc*rnx+2*ic+1] += negfac*weight*f*dt[1];
// 			}
// 		}
// 	}
// }
*/


void BackProjectionReconstructor::setup()
{
	image = new EMData();
	vector<int> size=params["size"];
	nx = size[0];
	ny = size[1];
	nz = size[2];
	image->set_size(nx, ny, nz);
}

EMData* BackProjectionReconstructor::preprocess_slice(const EMData* const slice, const Transform& t)
{

// 	EMData* return_slice = slice->process("normalize.edgemean");
// 	return_slice->process_inplace("filter.linearfourier");

	EMData* return_slice;

	return_slice = slice->process("filter.linearfourier");
//	return_slice = slice->copy();
	
// 	Transform tmp(t);
// 	tmp.set_rotation(Dict("type","eman")); // resets the rotation to 0 implicitly
// 	Vec2f trans = tmp.get_trans_2d();
// 	float scale = tmp.get_scale();
// 	bool mirror = tmp.get_mirror();
// 	if (trans[0] != 0 || trans[1] != 0 || scale != 1.0 ) {
// 		return_slice->transform(tmp);
// 	} 
// 	if ( mirror == true ) {
// 		return_slice->process_inplace("xform.flip",Dict("axis","x"));
// 	}

	return return_slice;
}

int BackProjectionReconstructor::insert_slice(const EMData* const input, const Transform &t, const float)
{
	if (!input) {
		LOGERR("try to insert NULL slice");
		return 1;
	}

	if (input->get_xsize() != input->get_ysize() || input->get_xsize() != nx) {
		LOGERR("tried to insert image that was not correction dimensions");
		return 1;
	}

// 	Transform * transform;
// 	if ( input->has_attr("xform.projection") ) {
// 		transform = (Transform*) (input->get_attr("xform.projection")); // assignment operator
// 	} else {
// 		transform = new Transform(t); // assignment operator
// 	}
	EMData* slice = preprocess_slice(input, t);

	// Clearly weight isn't a useful concept in back-projection without compensating with an exact-filter
// 	float weight = params["weight"];
// 	slice->mult(weight);

	EMData *tmp = new EMData();
	tmp->set_size(nx, ny, nz);

	float *slice_data = slice->get_data();
	float *tmp_data = tmp->get_data();

	size_t nxy = nx * ny;
	size_t nxy_size = nxy * sizeof(float);;
	for (int i = 0; i < nz; ++i) {
		memcpy(&tmp_data[nxy * i], slice_data, nxy_size);
	}
	tmp->update();

// 	transform->set_scale(1.0);
// 	transform->set_mirror(false);
// 	transform->set_trans(0,0,0);
// 	transform->invert();

	tmp->transform(t);
	image->add(*tmp);

// 	if(transform) {delete transform; transform=0;}
	delete tmp;
	delete slice;

	return 0;
}

int BackProjectionReconstructor::determine_slice_agreement(EMData*  input_slice, const Transform & arg, const float weight,bool sub)
{
	// Are these exceptions really necessary? (d.woolford)
	if (!input_slice) throw NullPointerException("EMData pointer (input image) is NULL");

	input_slice->set_attr("reconstruct_norm",1.0f);
	input_slice->set_attr("reconstruct_absqual",1.0f);
	input_slice->set_attr("reconstruct_weight",1.0f);

	return 0;

}

EMData *BackProjectionReconstructor::finish(bool)
{

	Symmetry3D* sym = Factory<Symmetry3D>::get((string)params["sym"]);
	vector<Transform> syms = sym->get_syms();

	for ( vector<Transform>::const_iterator it = syms.begin(); it != syms.end(); ++it ) {

//		it->printme();
		Transform t=*it;
		EMData *tmpcopy = image->process("xform",Dict("transform",(EMObject)&t));
		image->add(*tmpcopy);
		delete tmpcopy;
	}

	image->mult(1.0f/(float)sym->get_nsym());
	delete sym;
	if (image->get_xsize()==image->get_ysize() && image->get_ysize()==image->get_zsize()) {
		image->process_inplace("mask.sharp",Dict("outer_radius",image->get_xsize()/2-1));
	}
	else printf("No masking %d %d %d\n",image->get_xsize(),image->get_ysize(),image->get_zsize());
		
	EMData *ret = image;
	image = 0 ;
	return ret;
}

EMData* EMAN::padfft_slice( const EMData* const slice, const Transform& t, int npad )
{
	int nx = slice->get_xsize();
	int ny = slice->get_ysize();
	int padffted= slice->get_attr_default("padffted", 0);
	int ndim = (ny==1) ? 1 : 2;
	int extension = 2*padffted;  //  If 2, it means it is a Fourier file.

	if( ndim==2 && (nx-extension)!=ny )
	{
		// FIXME: What kind of exception should we throw here?
		throw std::runtime_error("Tried to padfft a 2D slice which is not square.");
	}

	EMData* padfftslice = NULL;
	if( padffted == 0) {
		// process 2D slice or 1D line -- subtract the average outside of the circle, zero-pad, fft extend, and fft
		EMData* temp = slice->average_circ_sub();

		padfftslice = temp->norm_pad( false, npad );
		checked_delete( temp );

		padfftslice->do_fft_inplace();
	} else {
		padfftslice = new EMData(*slice);
	}

	// shift the projection
	Vec2f trans = t.get_trans_2d();
	float sx = -trans[0];
	float sy = -trans[1];
	if(sx != 0.0f || sy != 0.0) padfftslice->process_inplace("filter.shift", Dict("x_shift", sx, "y_shift", sy, "z_shift", 0.0f));

	int remove = slice->get_attr_default("remove", 0);
	padfftslice->set_attr( "remove", remove );

	padfftslice->center_origin_fft();
	return padfftslice;
}

nn4Reconstructor::nn4Reconstructor()
{
	m_volume = NULL;
	m_wptr   = NULL;
}

nn4Reconstructor::nn4Reconstructor( const string& symmetry, int size, int npad )
{
	m_volume = NULL;
	m_wptr   = NULL;
	setup( symmetry, size, npad );
	load_default_settings();
	print_params();
}

nn4Reconstructor::~nn4Reconstructor()
{
	//if( m_delete_volume ) checked_delete(m_volume); 

	//if( m_delete_weight ) checked_delete( m_wptr );

	//checked_delete( m_result );
}

enum weighting_method { NONE, ESTIMATE, VORONOI };

float max2d( int kc, const vector<float>& pow_a )
{
	float max = 0.0;
	for( int i=-kc; i <= kc; ++i ) {
		for( int j=-kc; j <= kc; ++j ) {
			if( i==0 && j==0 ) continue;
			{
				int c = 2*kc+1 - std::abs(i) - std::abs(j);
				max = max + pow_a[c];
			}
		}
	}
	return max;
}

float max3d( int kc, const vector<float>& pow_a )
{
	float max = 0.0;
	for( int i=-kc; i <= kc; ++i ) {
		for( int j=-kc; j <= kc; ++j ) {
			for( int k=-kc; k <= kc; ++k ) {
				if( i==0 && j==0 && k==0 ) continue;
				// if( i!=0 )
				{
					int c = 3*kc+1 - std::abs(i) - std::abs(j) - std::abs(k);
					max = max + pow_a[c];
					// max = max + c * c;
					// max = max + c;
				}
			}
		}
	}
	return max;
}


void nn4Reconstructor::setup()
{
	int size = params["size"];
	int npad = params["npad"];


	string symmetry;
	if( params.has_key("symmetry") )  symmetry = params["symmetry"].to_str();
	else                               symmetry = "c1";

	if( params.has_key("ndim") )  m_ndim = params["ndim"];
	else                          m_ndim = 3;

	if( params.has_key( "snr" ) )  m_osnr = 1.0f/float( params["snr"] );
	else                           m_osnr = 0.0;

 	setup( symmetry, size, npad );
}

void nn4Reconstructor::setup( const string& symmetry, int size, int npad )
{
	m_weighting = ESTIMATE;
	m_wghta = 0.2f;

	m_symmetry = symmetry;
	m_npad = npad;
	m_nsym = Transform::get_nsym(m_symmetry);

	m_vnx = size;
	m_vny = size;
	m_vnz = (m_ndim==3) ? size : 1;

	m_vnxp = size*npad;
	m_vnyp = size*npad;
	m_vnzp = (m_ndim==3) ? size*npad : 1;

	m_vnxc = m_vnxp/2;
	m_vnyc = m_vnyp/2;
	m_vnzc = (m_ndim==3) ? m_vnzp/2 : 1;

	buildFFTVolume();
	buildNormVolume();
	
}


void nn4Reconstructor::buildFFTVolume() {
	int offset = 2 - m_vnxp%2;

	m_volume = params["fftvol"];

	if( m_volume->get_xsize() != m_vnxp+offset && m_volume->get_ysize() != m_vnyp && m_volume->get_zsize() != m_vnzp ) {
		m_volume->set_size(m_vnxp+offset,m_vnyp,m_vnzp);
		m_volume->to_zero();
	}
	// ----------------------------------------------------------------
	// Added by Zhengfan Yang on 03/15/07
	// Original author: please check whether my revision is correct and
	// other Reconstructor need similiar revision.
	if ( m_vnxp % 2 == 0 )  m_volume->set_fftodd(0);
	else                    m_volume->set_fftodd(1);
	// ----------------------------------------------------------------

	m_volume->set_nxc(m_vnxp/2);
	m_volume->set_complex(true);
	m_volume->set_ri(true);
	m_volume->set_fftpad(true);
	m_volume->set_attr("npad", m_npad);
	m_volume->set_array_offsets(0,1,1);
}

void nn4Reconstructor::buildNormVolume() {

	m_wptr = params["weight"];

	if( m_wptr->get_xsize() != m_vnxc+1 &&
		m_wptr->get_ysize() != m_vnyp &&
		m_wptr->get_zsize() != m_vnzp ) {
		m_wptr->set_size(m_vnxc+1,m_vnyp,m_vnzp);
		m_wptr->to_zero();
	}
	m_wptr->set_array_offsets(0,1,1);
}

void printImage( const EMData* line )
{
	Assert( line->get_zsize()==1 );


	int nx = line->get_xsize();
	int ny = line->get_ysize();
	for( int j=0; j < ny; ++j ) {
		for( int i=0; i < nx; ++i )  printf( "%10.3f ", line->get_value_at(i,j) );
		printf( "\n" );
	}
}



int nn4Reconstructor::insert_slice(const EMData* const slice, const Transform& t, const float weight) {
	// sanity checks
	if (!slice) {
		LOGERR("try to insert NULL slice");
		return 1;
	}

	int padffted= slice->get_attr_default( "padffted", 0 );
	if( m_ndim==3 ) {
		if ( padffted==0 && (slice->get_xsize()!=slice->get_ysize() || slice->get_xsize()!=m_vnx)  ) {
			// FIXME: Why doesn't this throw an exception?
			LOGERR("Tried to insert a slice that is the wrong size.");
			return 1;
		}
	} else {
		Assert( m_ndim==2 );
		if( slice->get_ysize() !=1 ) {
			LOGERR( "for 2D reconstruction, a line is excepted" );
			return 1;
		}
	}
	if( weight > 0.0f )  {

		EMData* padfft = padfft_slice( slice, t,  m_npad );

		if( m_ndim==3 ) {
			insert_padfft_slice( padfft, t, weight );
		} else {
			float alpha = padfft->get_attr( "alpha" );
			alpha = alpha/180.0f*M_PI;
			for(int i=0; i < m_vnxc+1; ++i ) {
				float xnew = i*cos(alpha);
				float ynew = -i*sin(alpha);
				float btqr = padfft->get_value_at( 2*i, 0, 0 );
				float btqi = padfft->get_value_at( 2*i+1, 0, 0 );
				if( xnew < 0.0 ) {
					xnew *= -1;
					ynew *= -1;
					btqi *= -1;
				}

				int ixn = int(xnew+0.5+m_vnxp) - m_vnxp;
				int iyn = int(ynew+0.5+m_vnyp) - m_vnyp;

				if(iyn < 0 ) iyn += m_vnyp;

				(*m_volume)( 2*ixn, iyn+1, 1 )   += btqr * weight;
				(*m_volume)( 2*ixn+1, iyn+1, 1 ) += btqi * weight;
				(*m_wptr)(ixn,iyn+1, 1) += weight;
			}

		}
		checked_delete( padfft );
	}
	return 0;
}

int nn4Reconstructor::insert_padfft_slice( EMData* padfft, const Transform& t, float weight )
{
	Assert( padfft != NULL );

	vector<Transform> tsym = t.get_sym_proj(m_symmetry);
	for (unsigned int isym=0; isym < tsym.size(); isym++)  m_volume->nn( m_wptr, padfft, tsym[isym], weight);
	
	return 0;
}


#define  tw(i,j,k)      tw[ i-1 + (j-1+(k-1)*iy)*ix ]

void circumfnn( EMData* win , int npad)
{
	float *tw = win->get_data();
	//  correct for the fall-off of NN interpolation using sinc functions
	//  mask and subtract circumference average
	int ix = win->get_xsize();
	int iy = win->get_ysize();
	int iz = win->get_zsize();
	int L2 = (ix/2)*(ix/2);
	int L2P = (ix/2-1)*(ix/2-1);

	int IP = ix/2+1;
	int JP = iy/2+1;
	int KP = iz/2+1;

	//  sinc functions tabulated for fall-off
	float* sincx = new float[IP+1];
	float* sincy = new float[JP+1];
	float* sincz = new float[KP+1];

	sincx[0] = 1.0f;
	sincy[0] = 1.0f;
	sincz[0] = 1.0f;

	float cor;
	if( npad == 1 )  cor = 1.0;
	else  cor = 4.0;

	float cdf = M_PI/(cor*ix);
	for (int i = 1; i <= IP; ++i)  sincx[i] = sin(i*cdf)/(i*cdf);
	cdf = M_PI/(cor*iy);
	for (int i = 1; i <= JP; ++i)  sincy[i] = sin(i*cdf)/(i*cdf);
	cdf = M_PI/(cor*iz);
	for (int i = 1; i <= KP; ++i)  sincz[i] = sin(i*cdf)/(i*cdf);
	for (int k = 1; k <= iz; ++k) {
		int kkp = abs(k-KP);
		for (int j = 1; j <= iy; ++j) {
			cdf = sincy[abs(j- JP)]*sincz[kkp];
			for (int i = 1; i <= ix; ++i)  tw(i,j,k) /= (sincx[abs(i-IP)]*cdf);
		}
	}

	delete[] sincx;
	delete[] sincy;
	delete[] sincz;

	float  TNR = 0.0f;
	size_t m = 0;
	for (int k = 1; k <= iz; ++k) {
		for (int j = 1; j <= iy; ++j) {
			for (int i = 1; i <= ix; ++i) {
				size_t LR = (k-KP)*(k-KP)+(j-JP)*(j-JP)+(i-IP)*(i-IP);
				if (LR >= (size_t)L2P && LR<=(size_t)L2) {
					TNR += tw(i,j,k);
					++m;
				}
			}
		}
	}

	TNR /=float(m);
	
	
	for (int k = 1; k <= iz; ++k) {
		for (int j = 1; j <= iy; ++j) {
			for (int i = 1; i <= ix; ++i) {
				size_t LR = (k-KP)*(k-KP)+(j-JP)*(j-JP)+(i-IP)*(i-IP);
				if (LR<=(size_t)L2) tw(i,j,k) -= TNR;
				else                tw(i,j,k) = 0.0f;

			}
		}
	}

}


void circumftrl( EMData* win , int npad)
{
	float *tw = win->get_data();
	//  correct for the fall-off of tri-linear interpolation using sinc^2 functions
	//  mask and subtract circumference average
	int ix = win->get_xsize();
	int iy = win->get_ysize();
	int iz = win->get_zsize();
	int L2 = (ix/2)*(ix/2);
	int L2P = (ix/2-1)*(ix/2-1);

	int IP = ix/2+1;
	int JP = iy/2+1;
	int KP = iz/2+1;

	//  sinc functions tabulated for fall-off
	float* sincx = new float[IP+1];
	float* sincy = new float[JP+1];
	float* sincz = new float[KP+1];

	sincx[0] = 1.0f;
	sincy[0] = 1.0f;
	sincz[0] = 1.0f;

	float cor;
	if( npad == 1 )  cor = 1.0;
	else  cor = 4.0;

	float cdf = M_PI/(cor*ix);
	for (int i = 1; i <= IP; ++i)  sincx[i] = pow(sin(i*cdf)/(i*cdf),2);
	cdf = M_PI/(cor*iy);
	for (int i = 1; i <= JP; ++i)  sincy[i] = pow(sin(i*cdf)/(i*cdf),2);
	cdf = M_PI/(cor*iz);
	for (int i = 1; i <= KP; ++i)  sincz[i] = pow(sin(i*cdf)/(i*cdf),2);
	for (int k = 1; k <= iz; ++k) {
		int kkp = abs(k-KP);
		for (int j = 1; j <= iy; ++j) {
			cdf = sincy[abs(j- JP)]*sincz[kkp];
			for (int i = 1; i <= ix; ++i)  tw(i,j,k) /= (sincx[abs(i-IP)]*cdf);
		}
	}

	delete[] sincx;
	delete[] sincy;
	delete[] sincz;

	float  TNR = 0.0f;
	size_t m = 0;
	for (int k = 1; k <= iz; ++k) {
		for (int j = 1; j <= iy; ++j) {
			for (int i = 1; i <= ix; ++i) {
				size_t LR = (k-KP)*(k-KP)+(j-JP)*(j-JP)+(i-IP)*(i-IP);
				if (LR >= (size_t)L2P && LR<=(size_t)L2) {
					TNR += tw(i,j,k);
					++m;
				}
			}
		}
	}

	TNR /=float(m);
	
	
	for (int k = 1; k <= iz; ++k) {
		for (int j = 1; j <= iy; ++j) {
			for (int i = 1; i <= ix; ++i) {
				size_t LR = (k-KP)*(k-KP)+(j-JP)*(j-JP)+(i-IP)*(i-IP);
				if (LR<=(size_t)L2) tw(i,j,k) -= TNR;
				else                tw(i,j,k) = 0.0f;

			}
		}
	}

}

EMData* nn4Reconstructor::finish(bool) {

	if( m_ndim == 3 ) {
		m_volume->symplane0(m_wptr);
	} else {
		for( int i=1; i <= m_vnyp; ++i ) {

			if( (*m_wptr)(0, i, 1)==0.0 ) {
				int j = m_vnyp + 1 - i;
				(*m_wptr)(0, i, 1) = (*m_wptr)(0, j, 1);
				(*m_volume)(0, i, 1) = (*m_volume)(0, j, 1);
				(*m_volume)(1, i, 1) = (*m_volume)(1, j, 1);
			}
		}
	}


	int box = 7;
	int kc = (box-1)/2;
	vector< float > pow_a( m_ndim*kc+1, 1.0 );
	for( unsigned int i=1; i < pow_a.size(); ++i ) pow_a[i] = pow_a[i-1] * exp(m_wghta);
	pow_a.back()=0.0;

	float alpha = 0.0;
	if( m_ndim==3) {
		int vol = box*box*box;
		float max = max3d( kc, pow_a );
		alpha = ( 1.0f - 1.0f/(float)vol ) / max;
	} else {
		int ara = box*box;
		float max = max2d( kc, pow_a );
		alpha = ( 1.0f - 1.0f/(float)ara ) / max;
	}
	int ix,iy,iz;
	for (iz = 1; iz <= m_vnzp; iz++) {
		for (iy = 1; iy <= m_vnyp; iy++) {
			for (ix = 0; ix <= m_vnxc; ix++) {
				if ( (*m_wptr)(ix,iy,iz) > 0) {//(*v) should be treated as complex!!
					float tmp = (-2*((ix+iy+iz)%2)+1)/((*m_wptr)(ix,iy,iz)+m_osnr);
					if( m_weighting == ESTIMATE ) {
						int cx = ix;
						int cy = (iy<=m_vnyc) ? iy - 1 : iy - 1 - m_vnyp;
						int cz = (iz<=m_vnzc) ? iz - 1 : iz - 1 - m_vnzp;
						float sum = 0.0;
						for( int ii = -kc; ii <= kc; ++ii ) {
							int nbrcx = cx + ii;
							if( nbrcx >= m_vnxc ) continue;
							for( int jj= -kc; jj <= kc; ++jj ) {
								int nbrcy = cy + jj;
								if( nbrcy <= -m_vnyc || nbrcy >= m_vnyc ) continue;

								int kcz = (m_ndim==3) ? kc : 0;
								for( int kk = -kcz; kk <= kcz; ++kk ) {
									int nbrcz = cz + kk;
									if( nbrcz <= -m_vnyc || nbrcz >= m_vnyc ) continue;
									if( nbrcx < 0 ) {
										nbrcx = -nbrcx;
							    			nbrcy = -nbrcy;
							    			nbrcz = -nbrcz;
									}
									int nbrix = nbrcx;
									int nbriy = nbrcy >= 0 ? nbrcy + 1 : nbrcy + 1 + m_vnyp;
									int nbriz = nbrcz >= 0 ? nbrcz + 1 : nbrcz + 1 + m_vnzp;
									if( (*m_wptr)( nbrix, nbriy, nbriz ) == 0 ) {
										int c = m_ndim*kc+1 - std::abs(ii) - std::abs(jj) - std::abs(kk);
										sum = sum + pow_a[c];
									}
								}
							}
						}
						float wght = 1.0f / ( 1.0f - alpha * sum );
						tmp = tmp * wght;
					}
//cout<<" mvol "<<ix<<"  "<<iy<<"  "<<iz<<"  "<<(*m_volume)(2*ix,iy,iz)<<"  "<<(*m_volume)(2*ix+1,iy,iz)<<"  "<<tmp<<"  "<<m_osnr<<endl;
					(*m_volume)(2*ix,iy,iz)   *= tmp;
					(*m_volume)(2*ix+1,iy,iz) *= tmp;
				}
			}
		}
	}

	//if(m_ndim==2) printImage( m_volume );

	// back fft
	m_volume->do_ift_inplace();

	// EMData* win = m_volume->window_center(m_vnx);
	int npad = m_volume->get_attr("npad");
	m_volume->depad();
	circumfnn( m_volume, npad );
	m_volume->set_array_offsets( 0, 0, 0 );

	return 0;
}
#undef  tw


nn4_rectReconstructor::nn4_rectReconstructor()
{
	m_volume = NULL;
	m_wptr   = NULL;
}

nn4_rectReconstructor::nn4_rectReconstructor( const string& symmetry, int size, int npad )
{
	m_volume = NULL;
	m_wptr   = NULL;
	setup( symmetry, size, npad );
	load_default_settings();
	print_params();
}

nn4_rectReconstructor::~nn4_rectReconstructor()
{
	//if( m_delete_volume ) checked_delete(m_volume);

	//if( m_delete_weight ) checked_delete( m_wptr );

	//checked_delete( m_result );
}


void nn4_rectReconstructor::setup()
{
	m_sizeofprojection = params["sizeprojection"];
	int npad = params["npad"];
	m_count=0;

	string symmetry;
	if( params.has_key("symmetry") )  symmetry = params["symmetry"].to_str();
	else                               symmetry = "c1";

	if( params.has_key("ndim") )  m_ndim = params["ndim"];             
	else 				m_ndim = 3;
    
	if( params.has_key( "snr" ) )  m_osnr = 1.0f/float( params["snr"] );
	else                           m_osnr = 0.0;

	setup( symmetry, m_sizeofprojection, npad );
}

void nn4_rectReconstructor::setup( const string& symmetry, int sizeprojection, int npad )
{
	m_weighting = ESTIMATE;
	m_wghta = 0.2f;
	m_symmetry = symmetry;
	m_npad = npad;
	m_nsym = Transform::get_nsym(m_symmetry);

	if( params.has_key("sizex") )  m_vnx = params["sizex"];
	else if(params.has_key("xratio")) 
		{
		float temp=params["xratio"];
		m_vnx=int(float(sizeprojection)*temp);
		}
	else                           m_vnx=sizeprojection;

	if( params.has_key("sizey") )  m_vny = params["sizey"];
	else if (params.has_key("yratio"))  
               {
		float temp=params["yratio"];
		 m_vny=int(float(sizeprojection)*temp);
		}
	else m_vny=sizeprojection;

	if( params.has_key("sizez") ) 
		m_vnz = params["sizez"];
	else 
		if (params.has_key("zratio"))
		{
			float temp=params["zratio"];
		 	m_vnz=int(float(sizeprojection)*temp);
		}
		else                          
			m_vnz = (m_ndim==3) ? sizeprojection : 1;
	
	m_xratio=float(m_vnx)/float(sizeprojection);	
	m_yratio=float(m_vny)/float(sizeprojection);
	m_zratio=float(m_vnz)/float(sizeprojection);

	m_vnxp = m_vnx*npad;
	m_vnyp = m_vny*npad;
	m_vnzp = (m_ndim==3) ? m_vnz*npad : 1;

	m_vnxc = m_vnxp/2;
	m_vnyc = m_vnyp/2;
	m_vnzc = (m_ndim==3) ? m_vnzp/2 : 1;

	buildFFTVolume();
	buildNormVolume();
}


void nn4_rectReconstructor::buildFFTVolume() {
	int offset = 2 - m_vnxp%2;

	m_volume = params["fftvol"];

	if( m_volume->get_xsize() != m_vnxp+offset && m_volume->get_ysize() != m_vnyp && m_volume->get_zsize() != m_vnzp ) {
		m_volume->set_size(m_vnxp+offset,m_vnyp,m_vnzp);
		m_volume->to_zero();
	}
	// ----------------------------------------------------------------
	// Added by Zhengfan Yang on 03/15/07
	// Original author: please check whether my revision is correct and
	// other Reconstructor need similiar revision.
	if ( m_vnxp % 2 == 0 )  m_volume->set_fftodd(0);
	else                    m_volume->set_fftodd(1);
	// ----------------------------------------------------------------

	m_volume->set_nxc(m_vnxp/2);
	m_volume->set_complex(true);
	m_volume->set_ri(true);
	m_volume->set_fftpad(true);
	m_volume->set_attr("npad", m_npad);
	m_volume->set_array_offsets(0,1,1);
}

void nn4_rectReconstructor::buildNormVolume() {

	m_wptr = params["weight"];

	if( m_wptr->get_xsize() != m_vnxc+1 &&
		m_wptr->get_ysize() != m_vnyp &&
		m_wptr->get_zsize() != m_vnzp ) {
		m_wptr->set_size(m_vnxc+1,m_vnyp,m_vnzp);
		m_wptr->to_zero();
	}

	m_wptr->set_array_offsets(0,1,1);
}

int nn4_rectReconstructor::insert_slice(const EMData* const slice, const Transform& t, const float weight) {
	// sanity checks


	if (!slice) {
		LOGERR("try to insert NULL slice");
		return 1;
	}

	int padffted= slice->get_attr_default( "padffted", 0 );
	if( m_ndim==3 ) {
		if ( padffted==0 && (slice->get_xsize()!=slice->get_ysize() || slice->get_xsize()!=m_sizeofprojection)  ) {
			// FIXME: Why doesn't this throw an exception?
			LOGERR("Tried to insert a slice that is the wrong size.");
			return 1;
		}
	} 
	if (m_ndim==2) {
		if( slice->get_ysize() !=1 ) {
			LOGERR( "for 2D reconstruction, a line is excepted" );
			return 1;
		}
	}
	if( weight > 0.0f )  {

		EMData* padfft = padfft_slice( slice, t,  m_npad );
		
		//Assert( mult > 0 );

		if( m_ndim==3 ) {
			insert_padfft_slice( padfft, t, weight );		
		} else {
			float ellipse_length,ellipse_step,cos_alpha,sin_alpha;
			int ellipse_length_int;
			float alpha = padfft->get_attr( "alpha" );
			alpha = alpha/180.0f*M_PI;
			int loop_range;
			float temp1,temp2;
					temp1=m_xratio*cos(alpha)*float(m_sizeofprojection*m_npad)/2;
			temp2=m_yratio*sin(alpha)*float(m_sizeofprojection*m_npad)/2;
			ellipse_length=sqrt(temp1*temp1+temp2*temp2);
			ellipse_length_int=int(ellipse_length);
			ellipse_step=0.5f*(m_sizeofprojection*m_npad)/float(ellipse_length_int);
			loop_range=ellipse_length_int;
			cos_alpha=temp1/ellipse_length;
			sin_alpha=temp2/ellipse_length;
			if(m_count%100==0) {
				std::cout<<"#############################################################"<<std::endl;
				std::cout<<"line insert start=="<<m_count<<std::endl;
				std::cout<<"ellipse length=="<<ellipse_length_int<<"ellips step=="<<ellipse_step<<std::endl;
				std::cout<<"loop_range"<<loop_range<<std::endl;
							std::cout<<"x and y ratio=="<<m_xratio<<"  "<<m_yratio<<std::endl;
				std::cout<<"cos sin of alpha=="<<cos(alpha)<<"   "<<sin(alpha)<<std::endl;
				std::cout<<"cos sin of alpha_new==="<<cos_alpha<<sin_alpha<<std::endl;
				std::cout<<"alpah dig==="<<cos_alpha<<sin_alpha<<std::endl;
				std::cout<<"prjection maximum==="<<loop_range*ellipse_step<<"ideal maximum"<<m_sizeofprojection*m_npad/2<<std::endl;
				std::cout<<"x_size=="<<m_volume->get_xsize()<<"y_size=="<<m_volume->get_ysize()<<std::endl;
				std::cout<<"#############################################################"<<std::endl;
			}
			for(int i=0; i <=loop_range; ++i ) {
				float xnew = i*cos_alpha;
				float ynew = -i*sin_alpha;
				if(m_count%100==0&&i==loop_range)
					std::cout<<"x_new=="<<xnew<<"Y_new=="<<ynew<<std::endl;
				float btqr=0,btqi=0;
				float xprj=i*ellipse_step;
				float t=xprj-int(xprj);
				btqr = (1-t)*padfft->get_value_at( 2*int(xprj), 0, 0 )+t*padfft->get_value_at( 2*(1+int(xprj)), 0, 0 );
				btqi = (1-t)*padfft->get_value_at( 2*int(xprj)+1, 0, 0 )+t*padfft->get_value_at( 2*(1+int(xprj))+1, 0, 0 );
				if( xnew < 0.0 ) {
					xnew *= -1;
					ynew *= -1;
					btqi *= -1;
				}

				int ixn = int(xnew+0.5+m_vnxp) - m_vnxp;
				int iyn = int(ynew+0.5+m_vnyp) - m_vnyp;

				if(iyn < 0 ) iyn += m_vnyp;
				if(m_count%100==0&&i==loop_range)
					std::cout<<"xnn=="<<ixn<<"ynn=="<<iyn<<std::endl;
				(*m_volume)( 2*ixn, iyn+1, 1 )   += btqr * weight;
				(*m_volume)( 2*ixn+1, iyn+1, 1 ) += btqi * weight;
				(*m_wptr)(ixn,iyn+1, 1) += weight;
			}


		}
		checked_delete( padfft );
		return 0;
	} else return 0;
}




int nn4_rectReconstructor::insert_padfft_slice( EMData* padded, const Transform& t, float weight )
{
	Assert( padded != NULL );
		
	vector<Transform> tsym = t.get_sym_proj(m_symmetry);
	for (unsigned int isym=0; isym < tsym.size(); isym++)
		m_volume->insert_rect_slice(m_wptr, padded, tsym[isym], m_sizeofprojection, m_xratio, m_yratio, m_zratio, m_npad, weight);

	return 0;

}

#define  tw(i,j,k)      tw[ i-1 + (j-1+(k-1)*iy)*ix ]
void circumfnn_rect( EMData* win , int npad)
{
	float *tw = win->get_data();
	//  correct for the fall-off
	//  mask and subtract circumfnnerence average
	int ix = win->get_xsize();
	int iy = win->get_ysize();
	int iz = win->get_zsize();

	int IP = ix/2+1;
	int JP = iy/2+1;
	int KP = iz/2+1;
	
	//  sinc functions tabulated for fall-off
	float* sincx = new float[IP+1];
	float* sincy = new float[JP+1];
	float* sincz = new float[KP+1];

	sincx[0] = 1.0f;
	sincy[0] = 1.0f;
	sincz[0] = 1.0f;

	float cdf = M_PI/float(npad*2*ix);
	for (int i = 1; i <= IP; ++i)  sincx[i] = sin(i*cdf)/(i*cdf);
	cdf = M_PI/float(npad*2*iy);
	for (int i = 1; i <= JP; ++i)  sincy[i] = sin(i*cdf)/(i*cdf);
	cdf = M_PI/float(npad*2*iz);
	for (int i = 1; i <= KP; ++i)  sincz[i] = sin(i*cdf)/(i*cdf);
	for (int k = 1; k <= iz; ++k) {
		int kkp = abs(k-KP);
		for (int j = 1; j <= iy; ++j) {
			cdf = sincy[abs(j- JP)]*sincz[kkp];
			for (int i = 1; i <= ix; ++i)  tw(i,j,k) /= (sincx[abs(i-IP)]*cdf);
		}
	}

	delete[] sincx;
	delete[] sincy;
	delete[] sincz;
	
	
	
	float dxx = 1.0f/float(0.25*ix*ix);
	float dyy = 1.0f/float(0.25*iy*iy);
	
	
	
	float LR2=(float(ix)/2-1)*(float(ix)/2-1)*dxx;

	float  TNR = 0.0f;
	size_t m = 0;
	for (int k = 1; k <= iz; ++k) {
		for (int j = 1; j <= iy; ++j) {
			for (int i = 1; i <= ix; ++i) {
				float LR = (j-JP)*(j-JP)*dyy+(i-IP)*(i-IP)*dxx;
				if (LR<=1.0f && LR >= LR2) {
					TNR += tw(i,j,k);
					++m;
				}
			}
		}
	}

	TNR /=float(m);
	

	for (int k = 1; k <= iz; ++k) {
		for (int j = 1; j <= iy; ++j) {
			for (int i = 1; i <= ix; ++i) {
				float LR = (j-JP)*(j-JP)*dyy+(i-IP)*(i-IP)*dxx;
				if (LR<=1.0f)  tw(i,j,k)=tw(i,j,k)-TNR;
			 	else 		tw(i,j,k) = 0.0f;	
			}
		}
	}

}
#undef tw 

EMData* nn4_rectReconstructor::finish(bool)
{
	
        if( m_ndim==3 ) {
		m_volume->symplane0_rect(m_wptr);
	} else {
		for( int i=1; i <= m_vnyp; ++i ) {

			if( (*m_wptr)(0, i, 1)==0.0 ) {
				int j = m_vnyp + 1 - i;
				(*m_wptr)(0, i, 1) = (*m_wptr)(0, j, 1);
				(*m_volume)(0, i, 1) = (*m_volume)(0, j, 1);
				(*m_volume)(1, i, 1) = (*m_volume)(1, j, 1);
			}
		}
	}


	int box = 7;
	int kc = (box-1)/2;
	vector< float > pow_a( m_ndim*kc+1, 1.0 );
	for( unsigned int i=1; i < pow_a.size(); ++i ) pow_a[i] = pow_a[i-1] * exp(m_wghta);
	pow_a.back()=0.0;

	float alpha = 0.0;
	if( m_ndim==3) {
		int vol = box*box*box;
		float max = max3d( kc, pow_a );
		alpha = ( 1.0f - 1.0f/(float)vol ) / max;
	} else {
		int ara = box*box;
		float max = max2d( kc, pow_a );
		alpha = ( 1.0f - 1.0f/(float)ara ) / max;
	}

	int ix,iy,iz;
	for (iz = 1; iz <= m_vnzp; iz++) {
		for (iy = 1; iy <= m_vnyp; iy++) {
			for (ix = 0; ix <= m_vnxc; ix++) {
				if ( (*m_wptr)(ix,iy,iz) > 0) {//(*v) should be treated as complex!!
					float tmp;
					tmp = (-2*((ix+iy+iz)%2)+1)/((*m_wptr)(ix,iy,iz)+m_osnr);
					
					if( m_weighting == ESTIMATE ) {
						int cx = ix;
						int cy = (iy<=m_vnyc) ? iy - 1 : iy - 1 - m_vnyp;
						int cz = (iz<=m_vnzc) ? iz - 1 : iz - 1 - m_vnzp;
						float sum = 0.0;
						for( int ii = -kc; ii <= kc; ++ii ) {
							int nbrcx = cx + ii;
							if( nbrcx >= m_vnxc ) continue;
							for( int jj= -kc; jj <= kc; ++jj ) {
								int nbrcy = cy + jj;
								if( nbrcy <= -m_vnyc || nbrcy >= m_vnyc ) continue;

								int kcz = (m_ndim==3) ? kc : 0;
								for( int kk = -kcz; kk <= kcz; ++kk ) {
									int nbrcz = cz + kk;
									if( nbrcz <= -m_vnyc || nbrcz >= m_vnyc ) continue;
									if( nbrcx < 0 ) {
										nbrcx = -nbrcx;
							    			nbrcy = -nbrcy;
							    			nbrcz = -nbrcz;
									}
									int nbrix = nbrcx;
									int nbriy = nbrcy >= 0 ? nbrcy + 1 : nbrcy + 1 + m_vnyp;
									int nbriz = nbrcz >= 0 ? nbrcz + 1 : nbrcz + 1 + m_vnzp;
									if( (*m_wptr)( nbrix, nbriy, nbriz ) == 0 ) {
										int c = m_ndim*kc+1 - std::abs(ii) - std::abs(jj) - std::abs(kk);
										sum = sum + pow_a[c];
									}
								}
							}
						}
						float wght = 1.0f / ( 1.0f - alpha * sum );
						tmp = tmp * wght;
					}
					(*m_volume)(2*ix,iy,iz)   *= tmp;
					(*m_volume)(2*ix+1,iy,iz) *= tmp;
				}
			}
		}
	}

	//if(m_ndim==2) printImage( m_volume );

	// back fft
	m_volume->do_ift_inplace();

	
 	int npad = m_volume->get_attr("npad");
 	m_volume->depad();
 	circumfnn_rect( m_volume, npad );
 	m_volume->set_array_offsets( 0, 0, 0 );

	return 0;
}


// Added By Zhengfan Yang on 03/16/07
// Beginning of the addition
// --------------------------------------------------------------------------------

nnSSNR_Reconstructor::nnSSNR_Reconstructor()
{
	m_volume = NULL;
	m_wptr   = NULL;
	m_wptr2  = NULL;
}

nnSSNR_Reconstructor::nnSSNR_Reconstructor( const string& symmetry, int size, int npad)
{
	m_volume = NULL;
	m_wptr   = NULL;
	m_wptr2  = NULL;

	setup( symmetry, size, npad );
}

nnSSNR_Reconstructor::~nnSSNR_Reconstructor()
{
	//if( m_delete_volume ) checked_delete(m_volume);

	//if( m_delete_weight ) checked_delete( m_wptr );

	//if( m_delete_weight2 ) checked_delete( m_wptr2 );

	//checked_delete( m_result );
}

void nnSSNR_Reconstructor::setup()
{
	int size = params["size"];
	int npad = params["npad"];

	string symmetry;
	if( params.has_key("symmetry") ) symmetry = params["symmetry"].to_str();
	else				 symmetry = "c1";

	setup( symmetry, size, npad );
}

void nnSSNR_Reconstructor::setup( const string& symmetry, int size, int npad )
{

	m_weighting = ESTIMATE;
	m_wghta = 0.2f;
	m_wghtb = 0.004f;

	m_symmetry = symmetry;
	m_npad = npad;
	m_nsym = Transform::get_nsym(m_symmetry);

	m_vnx = size;
	m_vny = size;
	m_vnz = size;

	m_vnxp = size*npad;
	m_vnyp = size*npad;
	m_vnzp = size*npad;

	m_vnxc = m_vnxp/2;
	m_vnyc = m_vnyp/2;
	m_vnzc = m_vnzp/2;

	buildFFTVolume();
	buildNormVolume();
	buildNorm2Volume();

}

void nnSSNR_Reconstructor::buildFFTVolume() {

	m_volume = params["fftvol"]; 
	m_volume->set_size(m_vnxp+ 2 - m_vnxp%2,m_vnyp,m_vnzp);
	m_volume->to_zero();
	if ( m_vnxp % 2 == 0 ) m_volume->set_fftodd(0);
	else                   m_volume->set_fftodd(1);

	m_volume->set_nxc(m_vnxc);
	m_volume->set_complex(true);
	m_volume->set_ri(true);
	m_volume->set_fftpad(true);
	m_volume->set_attr("npad", m_npad);
	m_volume->set_array_offsets(0,1,1);
}

void nnSSNR_Reconstructor::buildNormVolume() {

	m_wptr = params["weight"];
	m_wptr->set_size(m_vnxc+1,m_vnyp,m_vnzp);
	m_wptr->to_zero();
	m_wptr->set_array_offsets(0,1,1);
}

void nnSSNR_Reconstructor::buildNorm2Volume() {

	m_wptr2 = params["weight2"];
	m_wptr2->set_size(m_vnxc+1,m_vnyp,m_vnzp);
	m_wptr2->to_zero();
	m_wptr2->set_array_offsets(0,1,1);
}


int nnSSNR_Reconstructor::insert_slice(const EMData* const slice, const Transform& t, const float weight) {
	// sanity checks
	if (!slice) {
		LOGERR("try to insert NULL slice");
		return 1;
	}
	if( weight > 0.0f )  {
		int padffted=slice->get_attr_default( "padffted", 0 );

		if ( padffted==0 && (slice->get_xsize()!=slice->get_ysize() || slice->get_xsize()!=m_vnx)  ) {
			// FIXME: Why doesn't this throw an exception?
			LOGERR("Tried to insert a slice that has wrong size.");
			return 1;
		}

		EMData* padfft = padfft_slice( slice, t, m_npad );

		insert_padfft_slice( padfft, t, weight );

		checked_delete( padfft );
		return 0;
	} else return 0;
}

int nnSSNR_Reconstructor::insert_padfft_slice( EMData* padfft, const Transform& t, float weight )
{
	Assert( padfft != NULL );
	// insert slice for all symmetry related positions
	for (int isym=0; isym < m_nsym; isym++) {
		Transform tsym = t.get_sym(m_symmetry, isym);
		m_volume->nn_SSNR( m_wptr, m_wptr2, padfft, tsym, weight);
	}
	return 0;
}


EMData* nnSSNR_Reconstructor::finish(bool)
{
/*
  I changed the code on 05/15/2008 so it only returns variance.
  Lines commented out are marked by //#
  The version prior to the currect changes is r1.190
  PAP
*/
	int kz, ky;
	//#int iix, iiy, iiz;
	int box = 7;
	int kc = (box-1)/2;
	float alpha = 0.0;
	float argx, argy, argz;
	vector< float > pow_a( 3*kc+1, 1.0 );
	float w = params["w"];
	EMData* SSNR = params["SSNR"];
	//#EMData* vol_ssnr = new EMData();
	//#vol_ssnr->set_size(m_vnxp, m_vnyp, m_vnzp);
	//#vol_ssnr->to_zero();
	//#  new line follow
	EMData* vol_ssnr = params["vol_ssnr"];
	vol_ssnr->set_size(m_vnxp+ 2 - m_vnxp%2, m_vnyp ,m_vnzp);
	vol_ssnr->to_zero();
	if ( m_vnxp % 2 == 0 ) vol_ssnr->set_fftodd(0);
	else                   vol_ssnr->set_fftodd(1);
	vol_ssnr->set_nxc(m_vnxc);
	vol_ssnr->set_complex(true);
	vol_ssnr->set_ri(true);
	vol_ssnr->set_fftpad(false);
	//#

	float dx2 = 1.0f/float(m_vnxc)/float(m_vnxc);
	float dy2 = 1.0f/float(m_vnyc)/float(m_vnyc);
	float dz2 = 1.0f/Util::get_max(float(m_vnzc),1.0f)/Util::get_max(float(m_vnzc),1.0f);
	int   inc = Util::round(float(Util::get_max(m_vnxc,m_vnyc,m_vnzc))/w);

	SSNR->set_size(inc+1,4,1);

	float *nom    = new float[inc+1];
	float *denom  = new float[inc+1];
	int  *nn     = new int[inc+1];
	int  *ka     = new int[inc+1];
	float wght = 1.0f;
	for (int i = 0; i <= inc; i++) {
		nom[i] = 0.0f;
		denom[i] = 0.0f;
		nn[i] = 0;
		ka[i] = 0;
	}

	m_volume->symplane1(m_wptr, m_wptr2);

	if ( m_weighting == ESTIMATE ) {
		int vol = box*box*box;
		for( unsigned int i=1; i < pow_a.size(); ++i ) pow_a[i] = pow_a[i-1] * exp(m_wghta);
		pow_a[3*kc] = 0.0;
		float max = max3d( kc, pow_a );
		alpha = ( 1.0f - 1.0f/(float)vol ) / max;
	}

	for (int iz = 1; iz <= m_vnzp; iz++) {
		if ( iz-1 > m_vnzc ) kz = iz-1-m_vnzp; else kz = iz-1;
		argz = float(kz*kz)*dz2;
		for (int iy = 1; iy <= m_vnyp; iy++) {
			if ( iy-1 > m_vnyc ) ky = iy-1-m_vnyp; else ky = iy-1;
			argy = argz + float(ky*ky)*dy2;
			for (int ix = 0; ix <= m_vnxc; ix++) {
				float Kn = (*m_wptr)(ix,iy,iz);
				argx = std::sqrt(argy + float(ix*ix)*dx2);
				int r = Util::round(float(inc)*argx);
				if ( r >= 0 && Kn > 4.5f ) {
					if ( m_weighting == ESTIMATE ) {
						int cx = ix;
						int cy = (iy<=m_vnyc) ? iy - 1 : iy - 1 - m_vnyp;
						int cz = (iz<=m_vnzc) ? iz - 1 : iz - 1 - m_vnzp;

						float sum = 0.0;
						for( int ii = -kc; ii <= kc; ++ii ) {
							int nbrcx = cx + ii;
							if( nbrcx >= m_vnxc ) continue;
						    for ( int jj= -kc; jj <= kc; ++jj ) {
								int nbrcy = cy + jj;
								if( nbrcy <= -m_vnyc || nbrcy >= m_vnyc ) continue;
								for( int kk = -kc; kk <= kc; ++kk ) {
									int nbrcz = cz + jj;
									if ( nbrcz <= -m_vnyc || nbrcz >= m_vnyc ) continue;
									if( nbrcx < 0 ) {
										nbrcx = -nbrcx;
										nbrcy = -nbrcy;
										nbrcz = -nbrcz;
									}
		                        				int nbrix = nbrcx;
									int nbriy = nbrcy >= 0 ? nbrcy + 1 : nbrcy + 1 + m_vnyp;
									int nbriz = nbrcz >= 0 ? nbrcz + 1 : nbrcz + 1 + m_vnzp;
									if( (*m_wptr)( nbrix, nbriy, nbriz ) == 0 ) {
										int c = 3*kc+1 - std::abs(ii) - std::abs(jj) - std::abs(kk);
										sum = sum + pow_a[c];
									}
								}
							}
						}
						wght = 1.0f / ( 1.0f - alpha * sum );
					} // end of ( m_weighting == ESTIMATE )
					float nominator = std::norm(m_volume->cmplx(ix,iy,iz)/Kn);
					float denominator = ((*m_wptr2)(ix,iy,iz)-nominator)/(Kn-1.0f);
					// Skip Friedel related values
					if( (ix>0 || (kz>=0 && (ky>=0 || kz!=0)))) {
						if ( r <= inc ) {
							nom[r]   += nominator*wght;
							denom[r] += denominator/Kn*wght;
							nn[r]    += 2;
							ka[r]    += int(Kn);
						}
/*
						//#float  tmp = Util::get_max(nominator/denominator/Kn-1.0f,0.0f);
						//  Create SSNR as a 3D array (-n/2:n/2+n%2-1)
						iix = m_vnxc + ix; iiy = m_vnyc + ky; iiz = m_vnzc + kz;
						if( iix >= 0 && iix < m_vnxp && iiy >= 0 && iiy < m_vnyp && iiz >= 0 && iiz < m_vnzp )
							(*vol_ssnr)(iix, iiy, iiz) = tmp;
						// Friedel part
						iix = m_vnxc - ix; iiy = m_vnyc - ky; iiz = m_vnzc - kz;
						if( iix >= 0 && iix < m_vnxp && iiy >= 0 && iiy < m_vnyp && iiz >= 0 && iiz < m_vnzp )
							(*vol_ssnr)(iix, iiy, iiz) = tmp;
*/

					}
					(*vol_ssnr)(2*ix, iy-1, iz-1) = nominator*wght;//Kn;//denominator*wght;//
					//(*vol_ssnr)(2*ix, iy-1, iz-1) =  real(m_volume->cmplx(ix,iy,iz))*wght/Kn;
					//(*vol_ssnr)(2*ix+1, iy-1, iz-1) = imag(m_volume->cmplx(ix,iy,iz))*wght/Kn;
				} // end of Kn>4.5
			}
		}
	}

	for (int i = 0; i <= inc; i++)  {
		(*SSNR)(i,0,0) = nom[i];  ///(*SSNR)(i,0,0) = nom[i]/denom[i] - 1;///
		(*SSNR)(i,1,0) = denom[i];    // variance
		(*SSNR)(i,2,0) = static_cast<float>(nn[i]);
		(*SSNR)(i,3,0) = static_cast<float>(ka[i]);
	}
	vol_ssnr->update();
	
	delete[] nom;
	delete[] denom;
	delete[] nn;
	delete[] ka;

	return 0;
}

// -----------------------------------------------------------------------------------
// End of this addition

//####################################################################################
//** nn4 ctf reconstructor

nn4_ctfReconstructor::nn4_ctfReconstructor()
{
	m_volume  = NULL;
	m_wptr    = NULL;
}

nn4_ctfReconstructor::nn4_ctfReconstructor( const string& symmetry, int size, int npad, float snr, int sign )
{
	setup( symmetry, size, npad, snr, sign );
}

nn4_ctfReconstructor::~nn4_ctfReconstructor()
{
	//if( m_delete_volume ) checked_delete(m_volume);

	//if( m_delete_weight ) checked_delete( m_wptr );

	//checked_delete( m_result );
}

void nn4_ctfReconstructor::setup()
{
	if( ! params.has_key("size") ) throw std::logic_error("Error: image size is not given");

	int size = params["size"];
	int npad = params.has_key("npad") ? int(params["npad"]) : 2;
	// int sign = params.has_key("sign") ? int(params["sign"]) : 1;
	int sign = 1;
	string symmetry = params.has_key("symmetry")? params["symmetry"].to_str() : "c1";

	float snr = params["snr"];

	m_varsnr = params.has_key("varsnr") ? int(params["varsnr"]) : 0;
	setup( symmetry, size, npad, snr, sign );

}

void nn4_ctfReconstructor::setup( const string& symmetry, int size, int npad, float snr, int sign )
{
	m_weighting = ESTIMATE;
	if( params.has_key("weighting") ) {
		if( int( params["weighting"])==0 ) m_weighting = NONE;
	}



	m_wghta = 0.2f;
	m_wghtb = 0.004f;

	m_symmetry = symmetry;
	m_npad = npad;
	m_sign = sign;
	m_nsym = Transform::get_nsym(m_symmetry);

	m_snr = snr;

	m_vnx = size;
	m_vny = size;
	m_vnz = size;

	m_vnxp = size*npad;
	m_vnyp = size*npad;
	m_vnzp = size*npad;

	m_vnxc = m_vnxp/2;
	m_vnyc = m_vnyp/2;
	m_vnzc = m_vnzp/2;

	buildFFTVolume();
	buildNormVolume();
}

void nn4_ctfReconstructor::buildFFTVolume() {
	int offset = 2 - m_vnxp%2;

	m_volume = params["fftvol"];

	if( m_volume->get_xsize() != m_vnxp+offset && m_volume->get_ysize() != m_vnyp && m_volume->get_zsize() != m_vnzp ) {
		m_volume->set_size(m_vnxp+offset,m_vnyp,m_vnzp);
		m_volume->to_zero();
	}

	m_volume->set_nxc(m_vnxp/2);
	m_volume->set_complex(true);
	m_volume->set_ri(true);
	m_volume->set_fftpad(true);
	m_volume->set_attr("npad", m_npad);
	m_volume->set_array_offsets(0,1,1);
}

void nn4_ctfReconstructor::buildNormVolume()
{
	m_wptr = params["weight"];

	if( m_wptr->get_xsize() != m_vnxc+1 && m_wptr->get_ysize() != m_vnyp && m_wptr->get_zsize() != m_vnzp ) {
               m_wptr->set_size(m_vnxc+1,m_vnyp,m_vnzp);
               m_wptr->to_zero();
	}

	m_wptr->set_array_offsets(0,1,1);

}

int nn4_ctfReconstructor::insert_slice(const EMData* const slice, const Transform& t, const float weight)
{
	// sanity checks
	if (!slice) {
		LOGERR("try to insert NULL slice");
		return 1;
	}
	if( weight > 0.0f )  {
		int buffed = slice->get_attr_default( "buffed", 0 );
			if( buffed > 0 ) {
				insert_buffed_slice( slice, weight );
				return 0;
			}

		int padffted= slice->get_attr_default("padffted", 0);
		if( padffted==0 && (slice->get_xsize()!=slice->get_ysize() || slice->get_xsize()!=m_vnx)  ) {
			// FIXME: Why doesn't this throw an exception?
			LOGERR("Tried to insert a slice that is the wrong size.");
			return 1;
		}

		EMData* padfft = padfft_slice( slice, t, m_npad );

		float tmp = padfft->get_attr_default("ctf_applied", 0);
		int   ctf_applied = (int) tmp;

		// Generate 2D CTF (EMData object)
    	ctf_store_real::init( padfft->get_ysize(), padfft->get_attr( "ctf" ) );
    	EMData* ctf2d = ctf_store_real::get_ctf_real(); //This is in 2D projection plane

		int nx=ctf2d->get_xsize(),ny=ctf2d->get_ysize(),nz=ctf2d->get_zsize();
		float *ctf2d_ptr  = ctf2d->get_data();

		size_t size = (size_t)nx*ny*nz;
		if (!ctf_applied) {
			for (int i = 0; i < size; ++i) padfft->cmplx(i) *= ctf2d_ptr[i]; // Multiply padfft by CTF
		}

		for (int i = 0; i < size; ++i) ctf2d_ptr[i] *= ctf2d_ptr[i];     // Square 2D CTF
		
		insert_padfft_slice(padfft, ctf2d, t, weight);

		checked_delete( ctf2d );  
		checked_delete( padfft );

	}
	return 0;
}

int nn4_ctfReconstructor::insert_buffed_slice( const EMData* buffed, float weight )
{
	const float* bufdata = buffed->get_data();
	float* cdata = m_volume->get_data();
	float* wdata = m_wptr->get_data();

	int npoint = buffed->get_xsize()/4;
	for( int i=0; i < npoint; ++i ) {

		int pos2 = int( bufdata[4*i] );
		int pos1 = pos2 * 2;
		cdata[pos1  ] += bufdata[4*i+1]*weight;
		cdata[pos1+1] += bufdata[4*i+2]*weight;
		wdata[pos2  ] += bufdata[4*i+3]*weight;
/*
        std::cout << "pos1, pos2, ctfv1, ctfv2, ctf2: ";
        std::cout << pos1 << " " << bufdata[5*i+1] << " " << bufdata[5*i+2] << " ";
        std::cout << pos2 << " " << bufdata[5*i+4] << std::endl;
 */
	}
	return 0;
}

int nn4_ctfReconstructor::insert_padfft_slice( EMData* padfft, EMData* ctf2d2, const Transform& t, float weight)
{
	Assert( padfft != NULL );
	
	vector<float> abc_list;
	int abc_list_len = 0;	
	if (m_volume->has_attr("smear"))
	{
		abc_list = m_volume->get_attr("smear");
		abc_list_len = abc_list.size();
	}
				
	vector<Transform> tsym = t.get_sym_proj(m_symmetry);
	for (unsigned int isym=0; isym < tsym.size(); isym++) {
		if (abc_list_len == 0)
			m_volume->nn_ctf_exists(m_wptr, padfft, ctf2d2, tsym[isym], weight);
		else
			for (int i = 0; i < abc_list_len; i += 4) {
				m_volume->nn_ctf_exists(m_wptr, padfft, ctf2d2, tsym[isym] * Transform(Dict("type", "SPIDER", "phi",  abc_list[i], "theta", abc_list[i+1], "psi", abc_list[i+2])), weight * abc_list[i+3]);
			}
	}
	return 0;
}

EMData* nn4_ctfReconstructor::finish(bool)
{
	m_volume->set_array_offsets(0, 1, 1);
	m_wptr->set_array_offsets(0, 1, 1);
	m_volume->symplane0_ctf(m_wptr);

	int box = 7;
	int vol = box*box*box;
	int kc = (box-1)/2;
	vector< float > pow_a( 3*kc+1, 1.0 );
	for( unsigned int i=1; i < pow_a.size(); ++i ) pow_a[i] = pow_a[i-1] * exp(m_wghta);
	pow_a[3*kc]=0.0;


	float max = max3d( kc, pow_a );
	float alpha = ( 1.0f - 1.0f/(float)vol ) / max;
	float osnr = 1.0f/m_snr;

	// normalize
	int ix,iy,iz;
	for (iz = 1; iz <= m_vnzp; iz++) {
		for (iy = 1; iy <= m_vnyp; iy++) {
			for (ix = 0; ix <= m_vnxc; ix++) {
				if ( (*m_wptr)(ix,iy,iz) > 0.0f) {//(*v) should be treated as complex!!
					float tmp=0.0f;
					if( m_varsnr )  {
					    int iyp = (iy<=m_vnyc) ? iy - 1 : iy-m_vnyp-1;
					    int izp = (iz<=m_vnzc) ? iz - 1 : iz-m_vnzp-1;
						float freq = sqrt( (float)(ix*ix+iyp*iyp+izp*izp) );
						tmp = (-2*((ix+iy+iz)%2)+1)/((*m_wptr)(ix,iy,iz)+freq*osnr);//*m_sign;
					} else  {
						tmp = (-2*((ix+iy+iz)%2)+1)/((*m_wptr)(ix,iy,iz)+osnr);//*m_sign;
					}

			if( m_weighting == ESTIMATE ) {
				int cx = ix;
				int cy = (iy<=m_vnyc) ? iy - 1 : iy - 1 - m_vnyp;
				int cz = (iz<=m_vnzc) ? iz - 1 : iz - 1 - m_vnzp;
				float sum = 0.0;
				for( int ii = -kc; ii <= kc; ++ii ) {
					int nbrcx = cx + ii;
					if( nbrcx >= m_vnxc ) continue;
					for( int jj= -kc; jj <= kc; ++jj ) {
						int nbrcy = cy + jj;
						if( nbrcy <= -m_vnyc || nbrcy >= m_vnyc ) continue;
						for( int kk = -kc; kk <= kc; ++kk ) {
							int nbrcz = cz + jj;
							if( nbrcz <= -m_vnyc || nbrcz >= m_vnyc ) continue;
							if( nbrcx < 0 ) {
								nbrcx = -nbrcx;
								nbrcy = -nbrcy;
								nbrcz = -nbrcz;
							}

							int nbrix = nbrcx;
							int nbriy = nbrcy >= 0 ? nbrcy + 1 : nbrcy + 1 + m_vnyp;
							int nbriz = nbrcz >= 0 ? nbrcz + 1 : nbrcz + 1 + m_vnzp;
							if( (*m_wptr)( nbrix, nbriy, nbriz ) == 0.0 ) {
								int c = 3*kc+1 - std::abs(ii) - std::abs(jj) - std::abs(kk);
								sum = sum + pow_a[c];
					          		  // if(ix%20==0 && iy%20==0 && iz%20==0)
					           		 //   std::cout << boost::format( "%4d %4d %4d %4d %10.3f" ) % nbrix % nbriy % nbriz % c % sum << std::endl;
							}
						}
					}
				}
				float wght = 1.0f / ( 1.0f - alpha * sum );
/*
                        if(ix%10==0 && iy%10==0)
                        {
                            std::cout << boost::format( "%4d %4d %4d " ) % ix % iy %iz;
                            std::cout << boost::format( "%10.3f %10.3f %10.3f " )  % tmp % wght % sum;
                            std::  << boost::format( "%10.3f %10.3e " ) % pow_b[r] % alpha;
                            std::cout << std::endl;
                        }
 */
				tmp = tmp * wght;
				}
				(*m_volume)(2*ix,iy,iz) *= tmp;
				(*m_volume)(2*ix+1,iy,iz) *= tmp;
				}
			}
		}
	}

	// back fft
	m_volume->do_ift_inplace();
	int npad = m_volume->get_attr("npad");
	m_volume->depad();
	circumfnn( m_volume, npad );
	m_volume->set_array_offsets( 0, 0, 0 );

	return 0;
}


//####################################################################################
//** nn4 ctfw reconstructor

nn4_ctfwReconstructor::nn4_ctfwReconstructor()
{
	m_volume  = NULL;
	m_wptr    = NULL;
}

nn4_ctfwReconstructor::nn4_ctfwReconstructor( const string& symmetry, int size, int npad, float snr, int sign, int do_ctf )
{
	setup( symmetry, size, npad, snr, sign, do_ctf );
}

nn4_ctfwReconstructor::~nn4_ctfwReconstructor()
{
	////if( m_delete_volume ) checked_delete(m_volume);

	//if( m_delete_weight ) checked_delete( m_wptr );

	//checked_delete( m_result );
}

void nn4_ctfwReconstructor::setup()
{
	if( ! params.has_key("size") ) throw std::logic_error("Error: image size is not given");

	int size = params["size"];
	int npad = params.has_key("npad") ? int(params["npad"]) : 2;
	// int sign = params.has_key("sign") ? int(params["sign"]) : 1;
	int sign = 1;
	string symmetry = params.has_key("symmetry")? params["symmetry"].to_str() : "c1";

	float snr = params["snr"];
	int do_ctf = params["do_ctf"];

	m_varsnr = params.has_key("varsnr") ? int(params["varsnr"]) : 0;
	setup( symmetry, size, npad, snr, sign, do_ctf );

}

void nn4_ctfwReconstructor::setup( const string& symmetry, int size, int npad, float snr, int sign, int do_ctf )
{
	m_weighting = ESTIMATE;
	if( params.has_key("weighting") ) {
		if( int( params["weighting"])==0 ) m_weighting = NONE;
	}

	m_wghta = 0.2f;
	m_wghtb = 0.004f;

	m_symmetry = symmetry;
	m_npad = npad;
	m_sign = sign;
	m_nsym = Transform::get_nsym(m_symmetry);
	m_do_ctf = do_ctf;

	m_snr = snr;

	m_vnx = size;
	m_vny = size;
	m_vnz = size;

	m_vnxp = size*npad;
	m_vnyp = size*npad;
	m_vnzp = size*npad;

	m_vnxc = m_vnxp/2;
	m_vnyc = m_vnyp/2;
	m_vnzc = m_vnzp/2;

	buildFFTVolume();
	buildNormVolume();
	m_refvol = params["refvol"];

}

void nn4_ctfwReconstructor::buildFFTVolume() {
	int offset = 2 - m_vnxp%2;

	m_volume = params["fftvol"];

	if( m_volume->get_xsize() != m_vnxp+offset && m_volume->get_ysize() != m_vnyp && m_volume->get_zsize() != m_vnzp ) {
		m_volume->set_size(m_vnxp+offset,m_vnyp,m_vnzp);
		m_volume->to_zero();
	}

	m_volume->set_nxc(m_vnxp/2);
	m_volume->set_complex(true);
	m_volume->set_ri(true);
	m_volume->set_fftpad(true);
	m_volume->set_attr("npad", m_npad);
	m_volume->set_array_offsets(0,1,1);
}

void nn4_ctfwReconstructor::buildNormVolume()
{
	m_wptr = params["weight"];

	if( m_wptr->get_xsize() != m_vnxc+1 && m_wptr->get_ysize() != m_vnyp && m_wptr->get_zsize() != m_vnzp ) {
               m_wptr->set_size(m_vnxc+1,m_vnyp,m_vnzp);
               m_wptr->to_zero();
	}

	m_wptr->set_array_offsets(0,1,1);

}

int nn4_ctfwReconstructor::insert_slice(const EMData* const slice, const Transform& t, const float weight)
{
	// sanity checks
	if (!slice) {
		LOGERR("try to insert NULL slice");
		return 1;
	}
	if(weight >0.0f) {
		//cout<<"  insert_slice "<<m_do_ctf<<endl;
		/*
		int buffed = slice->get_attr_default( "buffed", 0 );
			if( buffed > 0 ) {
				insert_buffed_slice( slice, weight );
				return 0;
			}
		*/
		EMData* padfft = padfft_slice( slice, t, m_npad );;

		EMData* ctf2d = NULL;
		if( m_do_ctf == 1 ) {
			float tmp = padfft->get_attr_default("ctf_applied", 0);
			int   ctf_applied = (int) tmp;

			// Generate 2D CTF (EMData object)
			ctf_store_real::init( padfft->get_ysize(), padfft->get_attr( "ctf" ) );
			ctf2d = ctf_store_real::get_ctf_real(); //This is in 2D projection plane

			int nx=ctf2d->get_xsize(),ny=ctf2d->get_ysize(),nz=ctf2d->get_zsize();
			float *ctf2d_ptr  = ctf2d->get_data();

			size_t size = (size_t)nx*ny*nz;
			if (!ctf_applied) {
				for (int i = 0; i < size; ++i) padfft->cmplx(i) *= ctf2d_ptr[i]; // Multiply padfft by CTF
			}

			for (int i = 0; i < size; ++i) ctf2d_ptr[i] *= ctf2d_ptr[i];     // Square 2D CTF
		} else {
			int nx=padfft->get_xsize(),ny=padfft->get_ysize(),nz=padfft->get_zsize();
			//cout<<"  size of padfft "<<nx<<"   "<<ny<<"   "<<nz<<endl;
			ctf2d = new EMData();
			ctf2d->set_size(nx/2,ny,nz);
			float *ctf2d_ptr  = ctf2d->get_data();
			size_t size = (size_t)nx*ny*nz/2;
			for (int i = 0; i < size; ++i) ctf2d_ptr[i] = 1.0;
		}

		EMData* bckgnoise;
		bckgnoise = slice->get_attr("bckgnoise");

		insert_padfft_slice_weighted(padfft, ctf2d, bckgnoise, t, weight);

		checked_delete( ctf2d );  
		checked_delete( padfft );

	}
	return 0;
}

int nn4_ctfwReconstructor::insert_padfft_slice_weighted( EMData* padfft, EMData* ctf2d2, EMData* bckgnoise, const Transform& t, float weight )
{
	Assert( padfft != NULL );

	vector<float> abc_list;
	int abc_list_len = 0;	
	if (m_volume->has_attr("smear"))
	{
		abc_list = m_volume->get_attr("smear");
		abc_list_len = abc_list.size();
	}

	vector<Transform> tsym = t.get_sym_proj(m_symmetry);
	for (unsigned int isym=0; isym < tsym.size(); isym++) {
		if (abc_list_len == 0)
			m_volume->nn_ctfw(m_wptr, padfft, ctf2d2, bckgnoise, tsym[isym], weight);
		else
			for (int i = 0; i < abc_list_len; i += 4) 
				m_volume->nn_ctfw(m_wptr, padfft, ctf2d2, bckgnoise, tsym[isym] * Transform(Dict("type", "SPIDER", "phi",  abc_list[i], "theta", abc_list[i+1], "psi", abc_list[i+2])), weight * abc_list[i+3]);
	}
	return 0;
}


EMData* nn4_ctfwReconstructor::finish(bool compensate)
{
	m_volume->set_array_offsets(0, 1, 1);
	m_wptr->set_array_offsets(0, 1, 1);
	//cout <<  "  will set refvol  "  <<endl;
	//m_refvol->set_array_offsets(0, 1, 1);
	m_volume->symplane0_ctf(m_wptr);

	/*
	int box = 7;
	int vol = box*box*box;
	int kc = (box-1)/2;
	vector< float > pow_a( 3*kc+1, 1.0 );
	for( unsigned int i=1; i < pow_a.size(); ++i ) pow_a[i] = pow_a[i-1] * exp(m_wghta);
	pow_a[3*kc]=0.0;


	float max = max3d( kc, pow_a );
	float alpha = ( 1.0f - 1.0f/(float)vol ) / max;
	float osnr = 1.0f/m_snr;
	*/


	int ix,iy,iz;
	vector<float> count(m_vnyc+1, 0.0f);
	//  refvol carries fsc
	int  limitres = m_vnyc-1;
	if( (*m_refvol)(0) > 0.0f )  { // If fsc is set to zero, it will be straightforward reconstruction with snr = 1
		for (ix = 0; ix < m_vnyc; ix++) {
				//cout<<"  fsc  "<< m_vnyc-ix-1 <<"   "<<m_vnyc<<"   "<<(*m_refvol)(m_vnyc-ix-1)<<endl;
			  if( (*m_refvol)(m_vnyc-ix-1) == 0.0f )  limitres = m_vnyc-ix-2;
		}
//cout<<"   limitres   "<<limitres<<endl;

		vector<float> sigma2(m_vnyc+1, 0.0f);

		// compute sigma2
		for (iz = 1; iz <= m_vnzp; iz++) {
			int   izp = (iz<=m_vnzc) ? iz - 1 : iz-m_vnzp-1;
			float argz = float(izp*izp);
			for (iy = 1; iy <= m_vnyp; iy++) {
				int   iyp = (iy<=m_vnyc) ? iy - 1 : iy-m_vnyp-1;
				float argy = argz + float(iyp*iyp);
				for (ix = 0; ix <= m_vnxc; ix++) {
					if(ix>0 || (izp>=0 && (iyp>=0 || izp!=0))) {  //Skip Friedel related values
						float r = std::sqrt(argy + float(ix*ix));
						int  ir = int(r);
						if (ir <= limitres) {
							float frac = r - float(ir);
							float qres = 1.0f - frac;
							float temp = (*m_wptr)(ix,iy,iz);
							//cout<<" WEIGHTS "<<jx<<"  "<<jy<<"  "<<ir<<"  "<<temp<<"  "<<frac<<endl;
							//cout<<" WEIGHTS "<<ix<<"  "<<iy-1<<"  "<<iz-1<<"  "<<temp<<"  "<<endl;
							sigma2[ir]   += temp*qres;
							sigma2[ir+1] += temp*frac;
							count[ir]    += qres;
							count[ir+1]  += frac;
						}
					}
				}
			}
		}
		for (ix = 0; ix <= limitres; ix++) {
			if( count[ix] > 0.0f )  sigma2[ix] = sigma2[ix]/count[ix];
			//cout<<"  sigma2  "<< ix <<"   "<<sigma2[ix]<<endl;
		}
		float fudge = m_refvol->get_attr("fudge");
		// now counter will serve to keep fsc-derived stuff
		for (ix = 0; ix <= limitres; ix++)  count[ix] = fudge * sigma2[ix] * (1.0f - (*m_refvol)(ix))/(*m_refvol)(ix);  //fudge?
		count[limitres+1] = count[limitres];
		//for (ix = 0; ix <= limitres+1; ix++)  cout<<"  tau2  "<< ix <<"   "<<count[ix]<<endl;
	}


	// normalize
	float osnr = 1.0f;
	for (iz = 1; iz <= m_vnzp; iz++) {
		int   izp = (iz<=m_vnzc) ? iz - 1 : iz-m_vnzp-1;
		float argz = float(izp*izp);
		for (iy = 1; iy <= m_vnyp; iy++) {
			int   iyp = (iy<=m_vnyc) ? iy - 1 : iy-m_vnyp-1;
			float argy = argz + float(iyp*iyp);
			for (ix = 0; ix <= m_vnxc; ix++) {
				float r = std::sqrt(argy + float(ix*ix));
				int  ir = int(r);
				if (ir <= limitres) {
					if ( (*m_wptr)(ix,iy,iz) > 0.0f) {
						if( (*m_refvol)(0) > 0.0f && ir > -1) {
							float frac = r - float(ir);
							float qres = 1.0f - frac;
							osnr = qres*count[ir] + frac*count[ir+1];
							//if(osnr == 0.0f)  osnr = 1.0f/(0.001*(*m_wptr)(ix,iy,iz));
							//cout<<"  "<<iz<<"   "<<iy<<"   "<<"   "<<ix<<"   "<<ir<<"   "<<(*m_wptr)(ix,iy,iz)<<"   "<<osnr<<"      "<<(*m_volume)(2*ix,iy,iz)<<"      "<<(*m_volume)(2*ix+1,iy,iz)<<endl;
						}  else osnr = 0.0f;

						float tmp = ((*m_wptr)(ix,iy,iz)+osnr);

						if(tmp>0.0f) {
							tmp = (-2*((ix+iy+iz)%2)+1)/tmp;
						//cout<<" mvol "<<ix<<"  "<<iy<<"  "<<iz<<"  "<<(*m_volume)(2*ix,iy,iz)<<"  "<<(*m_volume)(2*ix+1,iy,iz)<<"  "<<tmp<<"  "<<osnr<<endl;
							(*m_volume)(2*ix,iy,iz)   *= tmp;
							(*m_volume)(2*ix+1,iy,iz) *= tmp;
						} else {
							(*m_volume)(2*ix,iy,iz)   = 0.0f;
							(*m_volume)(2*ix+1,iy,iz) = 0.0f;
						}
					}
				} else {
					(*m_volume)(2*ix,iy,iz)   = 0.0f;
					(*m_volume)(2*ix+1,iy,iz) = 0.0f;
				}
			}
		}
	}

	// back fft
	m_volume->do_ift_inplace();
	int npad = m_volume->get_attr("npad");
	m_volume->depad();
	if( compensate )  circumftrl( m_volume, npad );
	m_volume->set_array_offsets( 0, 0, 0 );

	return 0;
}

/*
// For postprocessing only multiply by +/- 1
						if(tmp>0.0f) {
							int mum = (-2*((ix+iy+iz)%2)+1);
						//cout<<" mvol "<<ix<<"  "<<iy<<"  "<<iz<<"  "<<(*m_volume)(2*ix,iy,iz)<<"  "<<(*m_volume)(2*ix+1,iy,iz)<<"  "<<tmp<<"  "<<osnr<<endl;
							(*m_volume)(2*ix,iy,iz)   *= mum;
							(*m_volume)(2*ix+1,iy,iz) *= mum;
							(*m_wptr)(ix,iy,iz) *= mum;
						} else {
							(*m_volume)(2*ix,iy,iz)   = 0.0f;
							(*m_volume)(2*ix+1,iy,iz) = 0.0f;
						}
					}
				} else {
					(*m_volume)(2*ix,iy,iz)   = 0.0f;
					(*m_volume)(2*ix+1,iy,iz) = 0.0f;
				}
			}
		}
	}

	// back fft
	m_volume->do_ift_inplace();
	int npad = m_volume->get_attr("npad");
	m_volume->depad();
	circumftrl( m_volume, npad );
	m_volume->set_array_offsets( 0, 0, 0 );

	return 0;
}
*/





/*
int nn4_ctfwReconstructor::insert_buffed_slice( const EMData* buffed, float weight )
{
	const float* bufdata = buffed->get_data();
	float* cdata = m_volume->get_data();
	float* wdata = m_wptr->get_data();

	int npoint = buffed->get_xsize()/4;
	for( int i=0; i < npoint; ++i ) {

		int pos2 = int( bufdata[4*i] );
		int pos1 = pos2 * 2;
		cdata[pos1  ] += bufdata[4*i+1]*weight;
		cdata[pos1+1] += bufdata[4*i+2]*weight;
		wdata[pos2  ] += bufdata[4*i+3]*weight;

        //std::cout << "pos1, pos2, ctfv1, ctfv2, ctf2: ";
        //std::cout << pos1 << " " << bufdata[5*i+1] << " " << bufdata[5*i+2] << " ";
        //std::cout << pos2 << " " << bufdata[5*i+4] << std::endl;
 
	}
	return 0;
}
*/
#ifdef False
EMData* nn4_ctfwReconstructor::finish(bool)
{
	m_volume->set_array_offsets(0, 1, 1);
	m_wptr->set_array_offsets(0, 1, 1);
	m_volume->symplane0_ctf(m_wptr);

	int box = 7;
	int vol = box*box*box;
	int kc = (box-1)/2;
	vector< float > pow_a( 3*kc+1, 1.0 );
	for( unsigned int i=1; i < pow_a.size(); ++i ) pow_a[i] = pow_a[i-1] * exp(m_wghta);
	pow_a[3*kc]=0.0;


	float max = max3d( kc, pow_a );
	float alpha = ( 1.0f - 1.0f/(float)vol ) / max;
	float osnr = 1.0f/m_snr;

	// normalize
	int ix,iy,iz;
	for (iz = 1; iz <= m_vnzp; iz++) {
		for (iy = 1; iy <= m_vnyp; iy++) {
			for (ix = 0; ix <= m_vnxc; ix++) {
				if ( (*m_wptr)(ix,iy,iz) > 0.0f) {//(*v) should be treated as complex!!
					float tmp=0.0f;
					if( m_varsnr )  {
					    int iyp = (iy<=m_vnyc) ? iy - 1 : iy-m_vnyp-1;
					    int izp = (iz<=m_vnzc) ? iz - 1 : iz-m_vnzp-1;
						float freq = sqrt( (float)(ix*ix+iyp*iyp+izp*izp) );
						tmp = (-2*((ix+iy+iz)%2)+1)/((*m_wptr)(ix,iy,iz)+freq*osnr);//*m_sign;
					} else  {
						tmp = (-2*((ix+iy+iz)%2)+1)/((*m_wptr)(ix,iy,iz)+osnr);//*m_sign;
					}

			if( m_weighting == ESTIMATE ) {
				int cx = ix;
				int cy = (iy<=m_vnyc) ? iy - 1 : iy - 1 - m_vnyp;
				int cz = (iz<=m_vnzc) ? iz - 1 : iz - 1 - m_vnzp;
				float sum = 0.0;
				for( int ii = -kc; ii <= kc; ++ii ) {
					int nbrcx = cx + ii;
					if( nbrcx >= m_vnxc ) continue;
					for( int jj= -kc; jj <= kc; ++jj ) {
						int nbrcy = cy + jj;
						if( nbrcy <= -m_vnyc || nbrcy >= m_vnyc ) continue;
						for( int kk = -kc; kk <= kc; ++kk ) {
							int nbrcz = cz + jj;
							if( nbrcz <= -m_vnyc || nbrcz >= m_vnyc ) continue;
							if( nbrcx < 0 ) {
								nbrcx = -nbrcx;
								nbrcy = -nbrcy;
								nbrcz = -nbrcz;
							}

							int nbrix = nbrcx;
							int nbriy = nbrcy >= 0 ? nbrcy + 1 : nbrcy + 1 + m_vnyp;
							int nbriz = nbrcz >= 0 ? nbrcz + 1 : nbrcz + 1 + m_vnzp;
							if( (*m_wptr)( nbrix, nbriy, nbriz ) == 0.0 ) {
								int c = 3*kc+1 - std::abs(ii) - std::abs(jj) - std::abs(kk);
								sum = sum + pow_a[c];
					          		  // if(ix%20==0 && iy%20==0 && iz%20==0)
					           		 //   std::cout << boost::format( "%4d %4d %4d %4d %10.3f" ) % nbrix % nbriy % nbriz % c % sum << std::endl;
							}
						}
					}
				}
				float wght = 1.0f / ( 1.0f - alpha * sum );
/*
                        if(ix%10==0 && iy%10==0)
                        {
                            std::cout << boost::format( "%4d %4d %4d " ) % ix % iy %iz;
                            std::cout << boost::format( "%10.3f %10.3f %10.3f " )  % tmp % wght % sum;
                            std::  << boost::format( "%10.3f %10.3e " ) % pow_b[r] % alpha;
                            std::cout << std::endl;
                        }
 */
				tmp = tmp * wght;
				}
				(*m_volume)(2*ix,iy,iz) *= tmp;
				(*m_volume)(2*ix+1,iy,iz) *= tmp;
				}
			}
		}
	}

	// back fft
	m_volume->do_ift_inplace();
	int npad = m_volume->get_attr("npad");
	m_volume->depad();
	circumfnn( m_volume, npad );
	m_volume->set_array_offsets( 0, 0, 0 );

	return 0;
}


{
	m_volume->set_array_offsets(0, 1, 1);
	m_wptr->set_array_offsets(0, 1, 1);
	m_volume->symplane0_ctf(m_wptr);

	int box = 7;
	int vol = box*box*box;
	int kc = (box-1)/2;
	vector< float > pow_a( 3*kc+1, 1.0 );
	for( unsigned int i=1; i < pow_a.size(); ++i ) pow_a[i] = pow_a[i-1] * exp(m_wghta);
	pow_a[3*kc]=0.0;


	float max = max3d( kc, pow_a );
	//float alpha = ( 1.0f - 1.0f/(float)vol ) / max;
	float osnr = 1.0f/m_snr;
 
    vector<float> sigma2(m_vnyc+1, 0.0f);
    vector<float> count(m_vnyc+1, 0.0f);

	int ix,iy,iz;
	// compute sigma2
	for (iz = 1; iz <= m_vnzp; iz++) {
		int   izp = (iz<=m_vnzc) ? iz - 1 : iz-m_vnzp-1;
		float argz = float(izp*izp);
		for (iy = 1; iy <= m_vnyp; iy++) {
			int   iyp = (iy<=m_vnyc) ? iy - 1 : iy-m_vnyp-1;
            float argy = argz + float(iyp*iyp);
			for (ix = 0; ix <= m_vnxc; ix++) {
			    if(ix>0 || (izp>=0 && (iyp>=0 || izp!=0))) {  //Skip Friedel related values
                    float r = std::sqrt(argy + float(ix*ix));
                    int  ir = int(r);
                    if (ir <= m_vnyc) {
                        float frac = r - float(ir);
                        float qres = 1.0f - frac;
                        float temp = (*m_wptr)(ix,iy,iz);
                        //cout<<" WEIGHTS "<<jx<<"  "<<jy<<"  "<<ir<<"  "<<temp<<"  "<<frac<<endl;
                        //cout<<" WEIGHTS "<<ix<<"  "<<iy-1<<"  "<<iz-1<<"  "<<temp<<"  "<<endl;
                        sigma2[ir]   += temp*qres;
                        sigma2[ir+1] += temp*frac;
                        count[ir]    += qres;
                        count[ir+1]  += frac;
                    }
                }
            }
        }
    }
    for (ix = 0; ix <= m_vnyc+1; ix++) {
        if( sigma2[ix] > 0.0f )  sigma2[ix] = count[ix]/sigma2[ix];
        cout<<"  1/sigma2  "<< ix <<"   "<<sigma2[ix]<<endl;
    }
    // now counter will serve to keep fsc-derived stuff
	//  refvol carries fsc
	float fudge = m_refvol->get_attr("fudge");
    for (ix = 0; ix <= m_vnyc+1; ix++)
		  count[ix] = Util::get_max(0.0f, Util::get_min( 0.999f, (*m_refvol)(ix) ) );
    for (ix = 0; ix <= m_vnyc+1; ix++)  count[ix] = count[ix]/(1.0f - count[ix]) * sigma2[ix];
    for (ix = 0; ix <= m_vnyc+1; ix++)  {
        if ( count[ix] >0.0f) count[ix] = fudge/count[ix];  //fudge?
    }
for (ix = 0; ix <= m_vnyc+1; ix++)  cout<<"  tau2  "<< ix <<"   "<<count[ix]<<"  m_wptr  "<<(*m_wptr)(ix,1,1)<<endl;

	// normalize
	for (iz = 1; iz <= m_vnzp; iz++) {
		int   izp = (iz<=m_vnzc) ? iz - 1 : iz-m_vnzp-1;
		float argz = float(izp*izp);
		for (iy = 1; iy <= m_vnyp; iy++) {
			int   iyp = (iy<=m_vnyc) ? iy - 1 : iy-m_vnyp-1;
            float argy = argz + float(iyp*iyp);
			for (ix = 0; ix <= m_vnxc; ix++) {
                    float r = std::sqrt(argy + float(ix*ix));
                    int  ir = int(r);
			//cout<<"  m_wptr  "<<(*m_wptr)(ix,iy,iz)<<endl;
                    if (ir <= m_vnyc) {
                        float frac = r - float(ir);
                        float qres = 1.0f - frac;
                        osnr = qres*count[ir] + frac*count[ir+1];
                        if(osnr == 0.0f)  osnr = 1.0f/(0.001*(*m_wptr)(ix,iy,iz));
                        //cout<<"  "<<iz<<"   "<<iy<<"   "<<"   "<<ix<<"   "<<(*m_wptr)(ix,iy,iz)<<"   "<<osnr<<"      "<<(*m_volume)(2*ix,iy,iz)<<"      "<<(*m_volume)(2*ix+1,iy,iz)<<endl;
 					    float tmp=((*m_wptr)(ix,iy,iz)+osnr);
					    //if( m_varsnr )  tmp = (-2*((ix+iy+iz)%2)+1)/((*m_wptr)(ix,iy,iz)+freq*osnr)*m_sign;
					    //else {
					    //cout<<"  "<<iz<<"  "<<iy<<"  "<<"  "<<ix<<"  "<<iz<<"  "<<"  "<<(*m_wptr)(ix,iy,iz)<<"  "<<osnr<<"  "<<endl;
					    if(tmp>0.0f) {
					        tmp = (-2*((ix+iy+iz)%2)+1)/tmp;


				/*if ( (*m_wptr)(ix,iy,iz) > 0.0f) {//(*v) should be treated as complex!!
					float tmp=0.0f;
					if( m_varsnr )  {
					    int iyp = (iy<=m_vnyc) ? iy - 1 : iy-m_vnyp-1;
					    int izp = (iz<=m_vnzc) ? iz - 1 : iz-m_vnzp-1;
						float freq = sqrt( (float)(ix*ix+iyp*iyp+izp*izp) );
						tmp = (-2*((ix+iy+iz)%2)+1)/((*m_wptr)(ix,iy,iz)+freq*osnr);//   *m_sign;
					} else  {
						tmp = (-2*((ix+iy+iz)%2)+1)/((*m_wptr)(ix,iy,iz)+osnr);//   *m_sign;
					}
				*/

			/*if( m_weighting == ESTIMATE ) {
				int cx = ix;
				int cy = (iy<=m_vnyc) ? iy - 1 : iy - 1 - m_vnyp;
				int cz = (iz<=m_vnzc) ? iz - 1 : iz - 1 - m_vnzp;
				float sum = 0.0;
				for( int ii = -kc; ii <= kc; ++ii ) {
					int nbrcx = cx + ii;
					if( nbrcx >= m_vnxc ) continue;
					for( int jj= -kc; jj <= kc; ++jj ) {
						int nbrcy = cy + jj;
						if( nbrcy <= -m_vnyc || nbrcy >= m_vnyc ) continue;
						for( int kk = -kc; kk <= kc; ++kk ) {
							int nbrcz = cz + jj;
							if( nbrcz <= -m_vnyc || nbrcz >= m_vnyc ) continue;
							if( nbrcx < 0 ) {
								nbrcx = -nbrcx;
								nbrcy = -nbrcy;
								nbrcz = -nbrcz;
							}

							int nbrix = nbrcx;
							int nbriy = nbrcy >= 0 ? nbrcy + 1 : nbrcy + 1 + m_vnyp;
							int nbriz = nbrcz >= 0 ? nbrcz + 1 : nbrcz + 1 + m_vnzp;
							if( (*m_wptr)( nbrix, nbriy, nbriz ) == 0.0 ) {
								int c = 3*kc+1 - std::abs(ii) - std::abs(jj) - std::abs(kk);
								sum = sum + pow_a[c];
					          		  // if(ix%20==0 && iy%20==0 && iz%20==0)
					           		 //   std::cout << boost::format( "%4d %4d %4d %4d %10.3f" ) % nbrix % nbriy % nbriz % c % sum << std::endl;
							}
						}
					}
				}
				float wght = 1.0f / ( 1.0f - alpha * sum );
/
                        if(ix%10==0 && iy%10==0)
                        {
                            std::cout << boost::format( "%4d %4d %4d " ) % ix % iy %iz;
                            std::cout << boost::format( "%10.3f %10.3f %10.3f " )  % tmp % wght % sum;
                            std::  << boost::format( "%10.3f %10.3e " ) % pow_b[r] % alpha;
                            std::cout << std::endl;
                        }
 /
				tmp = tmp * wght;
				}*/
				(*m_volume)(2*ix,iy,iz) *= tmp;
				(*m_volume)(2*ix+1,iy,iz) *= tmp;
				} else {
				(*m_volume)(2*ix,iy,iz)   = 0.0f;
				(*m_volume)(2*ix+1,iy,iz) = 0.0f;
				}
				}
			}
		}
	}

	// back fft
	m_volume->do_ift_inplace();
	int npad = m_volume->get_attr("npad");
	m_volume->depad();
	circumfnn( m_volume, npad );
	m_volume->set_array_offsets( 0, 0, 0 );

	return 0;
}
#endif


//####################################################################################
//** nn4 ctf rect reconstructor

nn4_ctf_rectReconstructor::nn4_ctf_rectReconstructor()
{
	m_volume  = NULL;
	m_wptr    = NULL;
}

nn4_ctf_rectReconstructor::nn4_ctf_rectReconstructor( const string& symmetry, int size, int npad, float snr, int sign )
{
	setup( symmetry, size, npad, snr, sign );
}

nn4_ctf_rectReconstructor::~nn4_ctf_rectReconstructor()
{
	//if( m_delete_volume ) checked_delete(m_volume);

	//if( m_delete_weight ) checked_delete( m_wptr );

	//checked_delete( m_result );
}

void nn4_ctf_rectReconstructor::setup()
{
	if( ! params.has_key("sizeprojection") ) throw std::logic_error("Error: projection size is not given");
	m_sizeofprojection = params["sizeprojection"];
	int npad = params.has_key("npad") ? int(params["npad"]) : 4;
	// int sign = params.has_key("sign") ? int(params["sign"]) : 1;
	int sign = 1;
	string symmetry = params.has_key("symmetry")? params["symmetry"].to_str() : "c1";

	float snr = params["snr"];

	m_varsnr = params.has_key("varsnr") ? int(params["varsnr"]) : 0;
	setup( symmetry, m_sizeofprojection, npad, snr, sign );

}

void nn4_ctf_rectReconstructor::setup( const string& symmetry, int sizeprojection, int npad, float snr, int sign )
{
	m_weighting = ESTIMATE;
	if( params.has_key("weighting") ) {
		int tmp = int( params["weighting"] );
		if( tmp==0 ) m_weighting = NONE;
	}

	m_wghta = 0.2f;
	m_wghtb = 0.004f;

	m_symmetry = symmetry;
	m_npad = npad;
	m_sign = sign;
	m_nsym = Transform::get_nsym(m_symmetry);

	m_snr = snr;
	if (params.has_key("sizex"))  m_vnx = params["sizex"];
	else if (params.has_key("xratio")) {
		float temp=params["xratio"];
		m_vnx=int(float(sizeprojection)*temp);
	} else  m_vnx=sizeprojection;

	if (params.has_key("sizey"))  m_vny = params["sizey"];
	else if (params.has_key("yratio")) {
		float temp=params["yratio"];
		m_vny=int(float(sizeprojection)*temp);
	}
	else m_vny=sizeprojection;
	
	if (params.has_key("sizez"))  m_vnz = params["sizez"];
	else if (params.has_key("zratio")) {
		float temp=params["zratio"];
		m_vnz=int(float(sizeprojection)*temp);
	}
	else m_vnz=sizeprojection;

	
	m_xratio=float(m_vnx)/float(sizeprojection);	
	m_yratio=float(m_vny)/float(sizeprojection);
	m_zratio=float(m_vnz)/float(sizeprojection);

	//std::cout<<"xratio=="<<m_xratio<<"yratio=="<<m_yratio<<std::endl;
	//std::cout<<"sx=="<<m_vnx<<"sy=="<<m_vny<<"sz=="<<m_vnz<<std::endl;

	m_vnxp = m_vnx*npad;
	m_vnyp = m_vny*npad;
	m_vnzp = m_vnz*npad;

	m_vnxc = m_vnxp/2;
	m_vnyc = m_vnyp/2;
	m_vnzc = m_vnzp/2;

	buildFFTVolume();
	buildNormVolume();
}

void nn4_ctf_rectReconstructor::buildFFTVolume() {
	int offset = 2 - m_vnxp%2;

	m_volume = params["fftvol"];

	if( m_volume->get_xsize() != m_vnxp+offset && m_volume->get_ysize() != m_vnyp && m_volume->get_zsize() != m_vnzp ) {
		m_volume->set_size(m_vnxp+offset,m_vnyp,m_vnzp);
		m_volume->to_zero();
	}

	m_volume->set_nxc(m_vnxp/2);
	m_volume->set_complex(true);
	m_volume->set_ri(true);
	m_volume->set_fftpad(true);
	m_volume->set_attr("npad", m_npad);
	m_volume->set_array_offsets(0,1,1);
}

void nn4_ctf_rectReconstructor::buildNormVolume()
{
	m_wptr = params["weight"];

	if( m_wptr->get_xsize() != m_vnxc+1 && m_wptr->get_ysize() != m_vnyp && m_wptr->get_zsize() != m_vnzp ) {
               m_wptr->set_size(m_vnxc+1,m_vnyp,m_vnzp);
               m_wptr->to_zero();
	}

	m_wptr->set_array_offsets(0,1,1);

}

int nn4_ctf_rectReconstructor::insert_slice(const EMData* const slice, const Transform& t, const float weight)
{
	// sanity checks
	if (!slice) {
		LOGERR("try to insert NULL slice");
		return 1;
	}
	if(weight >0.0f )  {
		int buffed = slice->get_attr_default( "buffed", 0 );
		if( buffed > 0 ) {
			insert_buffed_slice( slice, weight );
			return 0;
		}

		int padffted= slice->get_attr_default("padffted", 0);
		//if( padffted==0 && (slice->get_xsize()!=slice->get_ysize() || slice->get_xsize()!=m_vnx)  )
		if( padffted==0 && (slice->get_xsize()!=slice->get_ysize())  )
		{
			// FIXME: Why doesn't this throw an exception?
			LOGERR("Tried to insert a slice that is the wrong size.");
			return 1;
		}

		EMData* padfft = padfft_slice( slice, t, m_npad );

		insert_padfft_slice( padfft, t, weight );

		checked_delete( padfft );
	}
	return 0;
}

int nn4_ctf_rectReconstructor::insert_buffed_slice( const EMData* buffed, float weight )
{
	const float* bufdata = buffed->get_data();
	float* cdata = m_volume->get_data();
	float* wdata = m_wptr->get_data();

	int npoint = buffed->get_xsize()/4;
	for( int i=0; i < npoint; ++i ) {

		int pos2 = int( bufdata[4*i] );
		int pos1 = pos2 * 2;
		cdata[pos1  ] += bufdata[4*i+1]*weight;
		cdata[pos1+1] += bufdata[4*i+2]*weight;
		wdata[pos2  ] += bufdata[4*i+3]*weight;
/*
        std::cout << "pos1, pos2, ctfv1, ctfv2, ctf2: ";
        std::cout << pos1 << " " << bufdata[5*i+1] << " " << bufdata[5*i+2] << " ";
        std::cout << pos2 << " " << bufdata[5*i+4] << std::endl;
 */
	}
	return 0;
}


int nn4_ctf_rectReconstructor::insert_padfft_slice( EMData* padfft, const Transform& t, float weight )
{
	float tmp = padfft->get_attr("ctf_applied");
	int   ctf_applied = (int) tmp;
	vector<Transform> tsym = t.get_sym_proj(m_symmetry);
	for (unsigned int isym=0; isym < tsym.size(); isym++) {
		if(ctf_applied) m_volume->insert_rect_slice_ctf_applied(m_wptr, padfft, tsym[isym], m_sizeofprojection, m_xratio,m_yratio, m_zratio, m_npad, weight);			
		else            m_volume->insert_rect_slice_ctf(m_wptr, padfft, tsym[isym], m_sizeofprojection, m_xratio, m_yratio, m_zratio, m_npad, weight);	
	}

	return 0;
}

EMData* nn4_ctf_rectReconstructor::finish(bool)
{
	m_volume->set_array_offsets(0, 1, 1);
	m_wptr->set_array_offsets(0, 1, 1);
	m_volume->symplane0_rect(m_wptr);

	int box = 7;
	int vol = box*box*box;
	int kc = (box-1)/2;
	vector< float > pow_a( 3*kc+1, 1.0 );
	for( unsigned int i=1; i < pow_a.size(); ++i ) pow_a[i] = pow_a[i-1] * exp(m_wghta);
	pow_a[3*kc]=0.0;


	float max = max3d( kc, pow_a );
	float alpha = ( 1.0f - 1.0f/(float)vol ) / max;
	float osnr = 1.0f/m_snr;

	// normalize
	int ix,iy,iz;
	for (iz = 1; iz <= m_vnzp; iz++) {
		for (iy = 1; iy <= m_vnyp; iy++) {
			for (ix = 0; ix <= m_vnxc; ix++) {
				if ( (*m_wptr)(ix,iy,iz) > 0.0f) {//(*v) should be treated as complex!!
                    float tmp=0.0f;
                    if( m_varsnr ) {
                        int iyp = (iy<=m_vnyc) ? iy - 1 : iy-m_vnyp-1;
                        int izp = (iz<=m_vnzc) ? iz - 1 : iz-m_vnzp-1;
			            float freq = sqrt( (float)(ix*ix/(m_xratio*m_xratio)+iyp*iyp/(m_zratio*m_yratio)+izp*izp) );
                        tmp = (-2*((ix+iy+iz)%2)+1)/((*m_wptr)(ix,iy,iz)+freq*osnr)*m_sign;
                    } else {
                        tmp = (-2*((ix+iy+iz)%2)+1)/((*m_wptr)(ix,iy,iz)+osnr)*m_sign;
                    }

					if( m_weighting == ESTIMATE ) {
						int cx = ix;
						int cy = (iy<=m_vnyc) ? iy - 1 : iy - 1 - m_vnyp;
						int cz = (iz<=m_vnzc) ? iz - 1 : iz - 1 - m_vnzp;
						float sum = 0.0;
						for( int ii = -kc; ii <= kc; ++ii ) {
							int nbrcx = cx + ii;
							if( nbrcx >= m_vnxc ) continue;
							for( int jj= -kc; jj <= kc; ++jj ) {
								int nbrcy = cy + jj;
								if( nbrcy <= -m_vnyc || nbrcy >= m_vnyc ) continue;
								for( int kk = -kc; kk <= kc; ++kk ) {
									int nbrcz = cz + jj;
									if( nbrcz <= -m_vnyc || nbrcz >= m_vnyc ) continue;
									if( nbrcx < 0 ) {
										nbrcx = -nbrcx;
										nbrcy = -nbrcy;
										nbrcz = -nbrcz;
									}

									int nbrix = nbrcx;
									int nbriy = nbrcy >= 0 ? nbrcy + 1 : nbrcy + 1 + m_vnyp;
									int nbriz = nbrcz >= 0 ? nbrcz + 1 : nbrcz + 1 + m_vnzp;
									if( (*m_wptr)( nbrix, nbriy, nbriz ) == 0.0 ) {
										int c = 3*kc+1 - std::abs(ii) - std::abs(jj) - std::abs(kk);
										sum = sum + pow_a[c];
							          		  // if(ix%20==0 && iy%20==0 && iz%20==0)
							           		 //   std::cout << boost::format( "%4d %4d %4d %4d %10.3f" ) % nbrix % nbriy % nbriz % c % sum << std::endl;
									}
								}
							}
						}
						float wght = 1.0f / ( 1.0f - alpha * sum );
/*
                        if(ix%10==0 && iy%10==0)
                        {
                            std::cout << boost::format( "%4d %4d %4d " ) % ix % iy %iz;
                            std::cout << boost::format( "%10.3f %10.3f %10.3f " )  % tmp % wght % sum;
                            std::  << boost::format( "%10.3f %10.3e " ) % pow_b[r] % alpha;
                            std::cout << std::endl;
                        }
 */
						tmp = tmp * wght;
					}
					(*m_volume)(2*ix,iy,iz) *= tmp;
					(*m_volume)(2*ix+1,iy,iz) *= tmp;
				}
			}
		}
	}

	// back fft
	m_volume->do_ift_inplace();
	int npad = m_volume->get_attr("npad");
	m_volume->depad();
	circumfnn_rect( m_volume, npad );
	m_volume->set_array_offsets( 0, 0, 0 );
	return 0;
}



// Added By Zhengfan Yang on 04/11/07
// Beginning of the addition
// --------------------------------------------------------------------------------

nnSSNR_ctfReconstructor::nnSSNR_ctfReconstructor()
{
	m_volume  = NULL;
	m_wptr    = NULL;
	m_wptr2   = NULL;
	m_wptr3   = NULL;
}

nnSSNR_ctfReconstructor::nnSSNR_ctfReconstructor( const string& symmetry, int size, int npad, float snr, int sign)
{
	m_volume  = NULL;
	m_wptr    = NULL;
	m_wptr2   = NULL;
	m_wptr3   = NULL;

	setup( symmetry, size, npad, snr, sign );
}

nnSSNR_ctfReconstructor::~nnSSNR_ctfReconstructor()
{

	//if( m_delete_volume )  checked_delete(m_volume);
	//if( m_delete_weight )  checked_delete( m_wptr );
	//if( m_delete_weight2 ) checked_delete( m_wptr2 );
	//if( m_delete_weight3 ) checked_delete( m_wptr3 );
	//checked_delete( m_result );
}

void nnSSNR_ctfReconstructor::setup()
{
	int  size = params["size"];
	int  npad = params["npad"];
	int  sign = params["sign"];
	float snr = params["snr"];
	string symmetry;
	if( params.has_key("symmetry") )  symmetry = params["symmetry"].to_str();
	else                              symmetry = "c1";
	
	setup( symmetry, size, npad, snr, sign );
}
void nnSSNR_ctfReconstructor::setup( const string& symmetry, int size, int npad, float snr, int sign )
{

	m_weighting = ESTIMATE;
	m_wghta     = 0.2f;
	m_wghtb     = 0.004f;
	wiener      = 1;

	m_symmetry  = symmetry;
	m_npad      = npad;
	m_nsym      = Transform::get_nsym(m_symmetry);

	m_sign      = sign;
	m_snr	    = snr;

	m_vnx	    = size;
	m_vny	    = size;
	m_vnz	    = size;

	m_vnxp      = size*npad;
	m_vnyp      = size*npad;
	m_vnzp      = size*npad;

	m_vnxc      = m_vnxp/2;
	m_vnyc      = m_vnyp/2;
	m_vnzc      = m_vnzp/2;

	buildFFTVolume();
	buildNormVolume();
	buildNorm2Volume();
	buildNorm3Volume();
}

void nnSSNR_ctfReconstructor::buildFFTVolume() {

	int offset = 2 - m_vnxp%2;
	m_volume = params["fftvol"];

	m_volume->set_size(m_vnxp+offset,m_vnyp,m_vnzp);
	m_volume->to_zero();

	m_volume->set_fftodd(m_vnxp % 2);

	m_volume->set_nxc(m_vnxc);
	m_volume->set_complex(true);
	m_volume->set_ri(true); //(real, imaginary) instead of polar coordinate
	m_volume->set_fftpad(true);
	m_volume->set_attr("npad", m_npad);
	m_volume->set_array_offsets(0,1,1);
}


void nnSSNR_ctfReconstructor::buildNormVolume()
{
	m_wptr = params["weight"];
	m_wptr->set_size(m_vnxc+1,m_vnyp,m_vnzp);
	m_wptr->to_zero();
	m_wptr->set_array_offsets(0,1,1);
}

void nnSSNR_ctfReconstructor::buildNorm2Volume() {

	m_wptr2 = params["weight2"];
	m_wptr2->set_size(m_vnxc+1,m_vnyp,m_vnzp);
	m_wptr2->to_zero();
	m_wptr2->set_array_offsets(0,1,1);
}

void nnSSNR_ctfReconstructor::buildNorm3Volume() {

	m_wptr3 = params["weight3"];
	m_wptr3->set_size(m_vnxc+1,m_vnyp,m_vnzp);
	m_wptr3->to_zero();
	m_wptr3->set_array_offsets(0,1,1);
}

int nnSSNR_ctfReconstructor::insert_slice(const EMData *const  slice, const Transform& t, const float weight) {
	// sanity checks
	if (!slice) {
		LOGERR("try to insert NULL slice");
		return 1;
	}
	if( weight > 0.0f ) {
		int padffted= slice->get_attr_default("padffted", 0);
		if ( padffted==0 && (slice->get_xsize()!=slice->get_ysize() || slice->get_xsize()!=m_vnx)  ) {
			// FIXME: Why doesn't this throw an exception?
			LOGERR("Tried to insert a slice that is the wrong size.");
			return 1;
		}
		EMData* padfft = padfft_slice( slice, t, m_npad );

		insert_padfft_slice( padfft, t, weight );

		checked_delete( padfft );
	}
	return 0;
}
int nnSSNR_ctfReconstructor::insert_padfft_slice( EMData* padfft, const Transform& t, float weight )
{

	// insert slice for all symmetry related positions
	for (int isym=0; isym < m_nsym; isym++) {
		Transform tsym = t.get_sym(m_symmetry, isym);
		m_volume->nn_SSNR_ctf(m_wptr, m_wptr2, m_wptr3, padfft, tsym, weight);
	}
	return 0;
}

EMData* nnSSNR_ctfReconstructor::finish(bool)
{
/*
  I changed the code on 05/15/2008 so it only returns variance.
  Lines commented out are marked by //#
  The version prior to the currect changes is r1.190
  PAP
*/
	/***
	    m_volume ctf*(P^2D->3D(F^3D))
	    m_wptr   ctf^2
	    m_wptr2  |P^2D->3D(F^3D)|^2
	    m_wptr3  Kn
	    nominator = sum_rot [ wght*signal ]
	    denominator  = sum_rot[ wght*variance ]
						      ***/
	int kz, ky;
	int box = 7;
	int kc  = (box-1)/2;
	float alpha = 0.0;
	float argx, argy, argz;
	vector< float > pow_a( 3*kc+1, 1.0 );
	float w = params["w"];
	float dx2 = 1.0f/float(m_vnxc)/float(m_vnxc);
	float dy2 = 1.0f/float(m_vnyc)/float(m_vnyc);
	float dz2 = 1.0f/Util::get_max(float(m_vnzc),1.0f)/Util::get_max(float(m_vnzc),1.0f);
	int inc = Util::round(float(Util::get_max(m_vnxc,m_vnyc,m_vnzc))/w);

	EMData* SSNR = params["SSNR"];
	SSNR->set_size(inc+1,4,1);
	//#EMData* vol_ssnr = new EMData();
	//#vol_ssnr->set_size(m_vnxp, m_vnyp, m_vnzp);
	//#vol_ssnr->to_zero();
	//#  new linea follow
	EMData* vol_ssnr = params["vol_ssnr"];
	vol_ssnr->set_size(m_vnxp+ 2 - m_vnxp%2, m_vnyp ,m_vnzp);
	vol_ssnr->to_zero();
	if ( m_vnxp % 2 == 0 ) vol_ssnr->set_fftodd(0);
	else                   vol_ssnr->set_fftodd(1);
	vol_ssnr->set_nxc(m_vnxc);
	vol_ssnr->set_complex(true);
	vol_ssnr->set_ri(true);
	vol_ssnr->set_fftpad(false);
	//#
	float *nom    = new float[inc+1];
	float *denom  = new float[inc+1];
	int  *ka     = new int[inc+1];
	int  *nn     = new int[inc+1];
	float wght = 1.f;
	for (int i = 0; i <= inc; i++) {
		nom[i]   = 0.0f;
		denom[i] = 0.0f;
		nn[i]	 = 0;
		ka[i]	 = 0;
	}
	m_volume->symplane2(m_wptr, m_wptr2, m_wptr3);
	if ( m_weighting == ESTIMATE ) {
		int vol = box*box*box;
		for( unsigned int i=1; i < pow_a.size(); ++i ) pow_a[i] = pow_a[i-1] * exp(m_wghta);
		pow_a[3*kc] = 0.0;
		float max = max3d( kc, pow_a );
		alpha = ( 1.0f - 1.0f/(float)vol ) / max;
	}
	for (int iz = 1; iz <= m_vnzp; iz++) {
		if ( iz-1 > m_vnzc ) kz = iz-1-m_vnzp; else kz = iz-1;
		argz = float(kz*kz)*dz2;
		for (int iy = 1; iy <= m_vnyp; iy++) {
			if ( iy-1 > m_vnyc ) ky = iy-1-m_vnyp; else ky = iy-1;
			argy = argz + float(ky*ky)*dy2;
			for (int ix = 0; ix <= m_vnxc; ix++) {
				float Kn = (*m_wptr3)(ix,iy,iz);
				argx = std::sqrt(argy + float(ix*ix)*dx2);
				int r = Util::round(float(inc)*argx);
				if ( r >= 0 && Kn > 4.5f ) {
					if ( m_weighting == ESTIMATE ) {
						int cx = ix;
						int cy = (iy<=m_vnyc) ? iy - 1 : iy - 1 - m_vnyp;
						int cz = (iz<=m_vnzc) ? iz - 1 : iz - 1 - m_vnzp;
						float sum = 0.0;
						for( int ii = -kc; ii <= kc; ++ii ) {
							int nbrcx = cx + ii;
							if( nbrcx >= m_vnxc ) continue;
							for ( int jj= -kc; jj <= kc; ++jj ) {
								int nbrcy = cy + jj;
								if( nbrcy <= -m_vnyc || nbrcy >= m_vnyc ) continue;
								for( int kk = -kc; kk <= kc; ++kk ) {
									int nbrcz = cz + jj;
									if ( nbrcz <= -m_vnyc || nbrcz >= m_vnyc ) continue;
									if( nbrcx < 0 ) {
										nbrcx = -nbrcx;
										nbrcy = -nbrcy;
										nbrcz = -nbrcz;
									}
									int nbrix = nbrcx;
									int nbriy = nbrcy >= 0 ? nbrcy + 1 : nbrcy + 1 + m_vnyp;
									int nbriz = nbrcz >= 0 ? nbrcz + 1 : nbrcz + 1 + m_vnzp;
									if( (*m_wptr)( nbrix, nbriy, nbriz ) == 0 ) {
										int c = 3*kc+1 - std::abs(ii) - std::abs(jj) - std::abs(kk);
										sum = sum + pow_a[c];
									}
								}
							}
						}
// 						int r = std::abs(cx) + std::abs(cy) + std::abs(cz);
						wght = 1.0f / ( 1.0f - alpha * sum );
					} // end of ( m_weighting == ESTIMATE )
					float nominator   = std::norm(m_volume->cmplx(ix,iy,iz))/(*m_wptr)(ix,iy,iz);
					float denominator = ((*m_wptr2)(ix,iy,iz)-std::norm(m_volume->cmplx(ix,iy,iz))/(*m_wptr)(ix,iy,iz))/(Kn-1.0f);
					// Skip Friedel related values
					if( (ix>0 || (kz>=0 && (ky>=0 || kz!=0)))) {
						if ( r <= inc ) {
							nom[r]   += nominator*wght;
							denom[r] += denominator/(*m_wptr)(ix,iy,iz)*wght;
							nn[r]	 += 2;
							ka[r]	 += int(Kn);
						}
/*
						float  tmp = Util::get_max(nominator/denominator/(*m_wptr)(ix,iy,iz)-1.0f,0.0f);
						//  Create SSNR as a 3D array (-n/2:n/2+n%2-1)
						int iix = m_vnxc + ix; int iiy = m_vnyc + ky; int iiz = m_vnzc + kz;
						if( iix >= 0 && iix < m_vnxp && iiy >= 0 && iiy < m_vnyp && iiz >= 0 && iiz < m_vnzp )
							(*vol_ssnr)(iix, iiy, iiz) = tmp;
						// Friedel part
						iix = m_vnxc - ix; iiy = m_vnyc - ky; iiz = m_vnzc - kz;
						if( iix >= 0 && iix < m_vnxp && iiy >= 0 && iiy < m_vnyp && iiz >= 0 && iiz < m_vnzp )
							(*vol_ssnr)(iix, iiy, iiz) = tmp;
*/
					}
					(*vol_ssnr)(2*ix, iy-1, iz-1) = denominator*wght;
				} // end of Kn>4.5 or whatever
			}
		}
	}
	for (int i = 0; i <= inc; i++) {
		(*SSNR)(i,0,0) = nom[i];
		(*SSNR)(i,1,0) = denom[i];
		(*SSNR)(i,2,0) = static_cast<float>(nn[i]);
		(*SSNR)(i,3,0) = static_cast<float>(ka[i]);
	}
	vol_ssnr->update();

	delete[] nom;
	delete[] denom;
	delete[] nn;
	delete[] ka;

	return 0;
}

// -----------------------------------------------------------------------------------
// End of this addition


void EMAN::dump_reconstructors()
{
	dump_factory < Reconstructor > ();
}

map<string, vector<string> > EMAN::dump_reconstructors_list()
{
	return dump_factory_list < Reconstructor > ();
}





using std::ofstream;
using std::ifstream;


newfile_store::newfile_store( const string& filename, int npad, bool ctf )
    : m_bin_file( filename + ".bin" ),
      m_txt_file( filename + ".txt" )
{
	m_npad = npad;
	m_ctf = ctf;
}

newfile_store::~newfile_store( )
{
}

void newfile_store::add_image( EMData* emdata, const Transform& tf )
{
	if( m_bin_of == NULL ) {
	    m_bin_of = shared_ptr<ofstream>( new ofstream(m_bin_file.c_str(), std::ios::out|std::ios::binary) );
	    m_txt_of = shared_ptr<ofstream>( new ofstream(m_txt_file.c_str()) );
	}


	EMData* padfft = padfft_slice( emdata, tf, m_npad );

	int nx = padfft->get_xsize();
	int ny = padfft->get_ysize();
	int n2 = ny / 2;
	int n = ny;

	float voltage=0.0f, pixel=0.0f, Cs=0.0f, ampcont=0.0f, bfactor=0.0f, defocus=0.0f;

	if( m_ctf ) {
		Ctf* ctf = emdata->get_attr( "ctf" );
		Dict params = ctf->to_dict();
		voltage = params["voltage"];
		pixel	= params["apix"];
		Cs	= params["cs"];
		ampcont = params["ampcont"];
		bfactor = params["bfactor"];
		defocus = params["defocus"];
		if(ctf) {delete ctf; ctf=0;}
	}

	vector<point_t> points;
	for( int j=-ny/2+1; j <= ny/2; j++ ) {
		int jp = (j>=0) ? j+1 : ny+j+1;
		for( int i=0; i <= n2; ++i ) {
			int r2 = i*i + j*j;
			if( (r2<ny*ny/4) && !( (i==0) && (j<0) ) ) {
				float ctf;
				if( m_ctf ) {
					float ak = std::sqrt( r2/float(ny*ny) )/pixel;
					ctf = Util::tf( defocus, ak, voltage, Cs, ampcont, bfactor, 1);
				} else {
					ctf = 1.0;
				}

				float xnew = i*tf[0][0] + j*tf[1][0];
				float ynew = i*tf[0][1] + j*tf[1][1];
				float znew = i*tf[0][2] + j*tf[1][2];
				std::complex<float> btq;
				if (xnew < 0.) {
					xnew = -xnew;
					ynew = -ynew;
					znew = -znew;
					btq = conj(padfft->cmplx(i,jp-1));
				} else {
					btq = padfft->cmplx(i,jp-1);
				}

				int ixn = int(xnew + 0.5 + n) - n;
				int iyn = int(ynew + 0.5 + n) - n;
				int izn = int(znew + 0.5 + n) - n;
				if ((ixn <= n2) && (iyn >= -n2) && (iyn <= n2) && (izn >= -n2) && (izn <= n2)) {
					int ixf, iyf, izf;
					if (ixn >= 0) {
						int iza, iya;
						if (izn >= 0)
						    iza = izn + 1;
						else
						    iza = n + izn + 1;

						if (iyn >= 0)
						    iya = iyn + 1;
						else
						    iya = n + iyn + 1;

						ixf = ixn;
						iyf = iya;
						izf = iza;
					} else {
						int izt, iyt;
						if (izn > 0)
						    izt = n - izn + 1;
						else
						    izt = -izn + 1;

						if (iyn > 0)
						    iyt = n - iyn + 1;
						else
						    iyt = -iyn + 1;

						ixf = -ixn;
						iyf = iyt;
						izf = izt;
					}


					int pos2 = ixf + (iyf-1)*nx/2 + (izf-1)*ny*nx/2;
					float ctfv1 = btq.real() * ctf;
					float ctfv2 = btq.imag() * ctf;
					float ctf2 = ctf*ctf;

					point_t p;
					p.pos2 = pos2;
					p.real = ctfv1;
					p.imag = ctfv2;
					p.ctf2 = ctf2;

					points.push_back( p );
				}
			}
		}
	}


	int npoint = points.size();
	std::istream::off_type offset = (m_offsets.size()==0) ? 0 : m_offsets.back();
	offset += npoint*sizeof(point_t);
	m_offsets.push_back( offset );

	*m_txt_of << m_offsets.back() << std::endl;
	m_bin_of->write( (char*)(&points[0]), sizeof(point_t)*npoint );
	checked_delete( padfft );
}

void newfile_store::get_image( int id, EMData* buf )
{
	if( m_offsets.size()==0 ) {
		ifstream is( m_txt_file.c_str() );
		std::istream::off_type off;
		while( is >> off ) {
		    m_offsets.push_back( off );
		}

		m_bin_if = shared_ptr<std::ifstream>( new ifstream(m_bin_file.c_str(), std::ios::in|std::ios::binary) );
	}

	Assert( m_bin_if != NULL );

	std::istream::off_type offset = (id==0) ? 0 : m_offsets[id-1];
	Assert( offset >= 0 );
	m_bin_if->seekg( offset, std::ios::beg );


	if( m_bin_if->bad() || m_bin_if->fail() || m_bin_if->eof() ) {
		std::cout << "bad or fail or eof while fetching id, offset: " << id << " " << offset << std::endl;
		throw std::logic_error( "bad happen" );
	}

	int bufsize = (m_offsets[id] - offset)/sizeof(float);
	if( buf->get_xsize() != bufsize ) {
		buf->set_size( bufsize, 1, 1 );
	}

	char* data = (char*)(buf->get_data());
	m_bin_if->read( data, sizeof(float)*bufsize );
	buf->update();
}

void newfile_store::read( int nprj )
{
	if( m_offsets.size()==0 ) {
	    ifstream is( m_txt_file.c_str() );
	    std::istream::off_type off;
	    while( is >> off ) {
		m_offsets.push_back( off );
	    }
	}

	if( m_bin_if==NULL ) {
	    m_bin_if = shared_ptr< ifstream>( new ifstream(m_bin_file.c_str(), std::ios::in|std::ios::binary) );
	}


	int npoint = m_offsets[0]/sizeof(point_t);
	std::ios::off_type prjsize = m_offsets[0];

	try {
	    m_points.resize(nprj * npoint);
	}
	catch( std::exception& e ) {
	    std::cout << "Error: " << e.what() << std::endl;
	}

	int ip = 0;
	for( int i=0; i < nprj; ++i ) {
		m_bin_if->read( (char*)(&m_points[ip]), prjsize );
		if( m_bin_if->bad() || m_bin_if->fail() || m_bin_if->eof() )  {
			std::cout << "Error: file hander bad or fail or eof" << std::endl;
			return;
		}
		ip += npoint;
	}
}

void newfile_store::add_tovol( EMData* fftvol, EMData* wgtvol, const vector<int>& mults, int pbegin, int pend )
{
	float* vdata = fftvol->get_data();
	float* wdata = wgtvol->get_data();

	int npoint = m_offsets[0]/sizeof(point_t);
//    Assert( int(mults.size())==nprj );
	Assert( int(m_points.size())== (pend - pbegin)*npoint );

	for( int iprj=pbegin; iprj < pend; ++iprj ) {
		int m = mults[iprj];
		if( m==0 ) continue;

		int ipt = (iprj-pbegin)*npoint;
		for( int i=0; i < npoint; ++i )  {
		    int pos2 = m_points[ipt].pos2;
		    int pos1 = pos2*2;

		    wdata[pos2]   += m_points[ipt].ctf2*m;
		    vdata[pos1]   += m_points[ipt].real*m;
		    vdata[pos1+1] += m_points[ipt].imag*m;
		    ++ipt;
		}
	}
}

void newfile_store::restart()
{
    m_bin_if = shared_ptr< ifstream>( new ifstream(m_bin_file.c_str(), std::ios::in|std::ios::binary) );
}

file_store::file_store(const string& filename, int npad, int write, bool ctf)
    : m_bin_file(filename + ".bin"),
      m_txt_file(filename + ".txt")
{
    m_ctf = ctf;
    m_prev = -1;
    m_npad = npad;
    m_write = write;
}

file_store::~file_store()
{
}

void file_store::add_image( EMData* emdata, const Transform& tf )
{

    EMData* padfft = padfft_slice( emdata, tf, m_npad );

    float* data = padfft->get_data();

    if( m_write && m_bin_ohandle == NULL )
    {
        m_bin_ohandle = shared_ptr< ofstream >( new ofstream(m_bin_file.c_str(), std::ios::out | std::ios::binary) );
        m_txt_ohandle = shared_ptr< ofstream >( new ofstream(m_txt_file.c_str() ) );
        if( m_ctf )
		 *m_txt_ohandle << "Cs pixel voltage ctf_applied amp_contrast defocus ";

	*m_txt_ohandle << "phi theta psi" << std::endl;
    }

    m_x_out = padfft->get_xsize();
    m_y_out = padfft->get_ysize();
    m_z_out = padfft->get_zsize();
    m_totsize = m_x_out*m_y_out*m_z_out;

    if( m_ctf )
    {
        Ctf* ctf = padfft->get_attr( "ctf" );
        Dict ctf_params = ctf->to_dict();

        m_ctf_applied = padfft->get_attr( "ctf_applied" );

        m_Cs = ctf_params["cs"];
        m_pixel = ctf_params["apix"];
        m_voltage = ctf_params["voltage"];
        m_amp_contrast = ctf_params["ampcont"];
        m_defocuses.push_back( ctf_params["defocus"] );
        if(ctf) {delete ctf; ctf=0;}
    }

    Dict params = tf.get_rotation( "spider" );
    float phi = params.get( "phi" );
    float tht = params.get( "theta" );
    float psi = params.get( "psi" );


    m_phis.push_back( phi );
    m_thetas.push_back( tht );
    m_psis.push_back( psi );

    if( m_write )
    {
        m_bin_ohandle->write( (char*)data, sizeof(float)*m_totsize );

        if( m_ctf )
        {
            *m_txt_ohandle << m_Cs << " ";
            *m_txt_ohandle << m_pixel << " ";
            *m_txt_ohandle << m_voltage << " ";
            *m_txt_ohandle << m_ctf_applied << " ";
            *m_txt_ohandle << m_amp_contrast << " ";
            *m_txt_ohandle << m_defocuses.back() << " ";
        }
        *m_txt_ohandle << m_phis.back() << " ";
        *m_txt_ohandle << m_thetas.back() << " ";
        *m_txt_ohandle << m_psis.back() << " ";
        *m_txt_ohandle << m_x_out << " ";
        *m_txt_ohandle << m_y_out << " ";
        *m_txt_ohandle << m_z_out << " ";
        *m_txt_ohandle << m_totsize << std::endl;
    }

    checked_delete(padfft);

}

void file_store::get_image( int id, EMData* padfft )
{

    if( m_phis.size() == 0 ) {
        ifstream m_txt_ifs( m_txt_file.c_str() );

		if( !m_txt_ifs ) std::cerr << "Error: file " << m_txt_file << " does not exist" << std::endl;


		string line;
		std::getline( m_txt_ifs, line );

		float first, defocus, phi, theta, psi;



		while( m_txt_ifs >> first ) {

			if( m_ctf ) {
				m_Cs = first;
				m_txt_ifs >> m_pixel >> m_voltage;
				m_txt_ifs >> m_ctf_applied >> m_amp_contrast;
				m_txt_ifs >> defocus >> phi >> theta >> psi;
				m_defocuses.push_back( defocus );
			} else {
				phi = first;
				m_txt_ifs >> theta >> psi;
			}

			m_txt_ifs >> m_x_out >> m_y_out >> m_z_out >> m_totsize;
			m_phis.push_back( phi );
			m_thetas.push_back( theta );
			m_psis.push_back( psi );
		}
	}

	Assert( m_ihandle != NULL );

	std::istream::off_type offset = id*sizeof(float)*m_totsize;
	Assert( offset >= 0 );

	if( offset > 0 ) m_ihandle->seekg(offset, std::ios::beg);


	if( m_ihandle->bad() )
	{
		std::cout << "bad while fetching id, offset: " << id << " " << offset << std::endl;
		throw std::logic_error( "bad happen" );
	}

	if( m_ihandle->fail() )
	{
		std::cout << "fail while fetching id, offset, curoff: " << id << " " << offset << std::endl;
		throw std::logic_error( "fail happen" );
	}

	if( m_ihandle->eof() )
	{
		std::cout << "eof while fetching id, offset: " << id << " " << offset << std::endl;
		throw std::logic_error( "eof happen" );
	}

	if( padfft->get_xsize() != m_x_out ||
		padfft->get_ysize() != m_y_out ||
		padfft->get_zsize() != m_z_out )
	{
		padfft->set_size(m_x_out, m_y_out, m_z_out);
	}

	char* data = (char*)(padfft->get_data());
	m_ihandle->read( data, sizeof(float)*m_totsize );
	padfft->update();

	if( m_ctf ) {
		padfft->set_attr( "Cs", m_Cs );
		padfft->set_attr( "Pixel_size", m_pixel );
		padfft->set_attr( "voltage", m_voltage );
		padfft->set_attr( "ctf_applied", m_ctf_applied );
		padfft->set_attr( "amp_contrast", m_amp_contrast );
		padfft->set_attr( "defocus", m_defocuses[id] );
	}

	padfft->set_attr( "padffted", 1 );
	padfft->set_attr( "phi", m_phis[id] );
	padfft->set_attr( "theta", m_thetas[id] );
	padfft->set_attr( "psi", m_psis[id] );

}

void file_store::restart( )
{
    if( m_ihandle == NULL ) m_ihandle = shared_ptr< ifstream >( new ifstream(m_bin_file.c_str(), std::ios::in | std::ios::binary) );


    if( m_ihandle->bad() || m_ihandle->fail() || m_ihandle->eof() ) m_ihandle->open( m_bin_file.c_str(), std::ios::binary );

    m_ihandle->seekg( 0, std::ios::beg );
}

/* vim: set ts=4 noet: */
