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

#include "emdata.h"
#include "emfft.h"

#include <cstring>
#include <cstdio>

#include  "gsl_sf_result.h"
#include  "gsl_sf_bessel.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <utility>
#include <cmath>
#include "util.h"

//#ifdef EMAN2_USING_CUDA
//#include "cuda/cuda_processor.h"
//#endif

using namespace EMAN;
using namespace std;
typedef vector< pair<float,int> > vp;

EMData *EMData::do_fft() const
{
	ENTERFUNC;
#ifdef FFT_CACHING
	if (fftcache!=0) {
		return fftcache->copy();
	}
#endif //FFT_CACHING

	if (is_complex() ) { // ming add 08/17/2010
#ifdef NATIVE_FFT
		LOGERR(" NATIVE_FFT does not support complex to complex.");  // PAP
		throw ImageFormatException("real image expected. Input image is complex image.");
#else
		EMData *temp_in=copy();
		EMData *dat= copy_head();
		int offset;
		if(is_fftpadded()) {
			offset = is_fftodd() ? 1 : 2;
		}
		else offset=0;
		//printf("offset=%d\n",offset);
		EMfft::complex_to_complex_nd(temp_in->get_data(),dat->get_data(),nx-offset,ny,nz);

		if(dat->get_ysize()==1 && dat->get_zsize()==1) dat->set_complex_x(true);

		dat->update();
		delete temp_in;
		EXITFUNC;
		return dat;
#endif // NATIVE_FFT
	} else {
		int nxreal = nx;
		int offset = 2 - nx%2;
		int nx2 = nx + offset;
		EMData* dat = copy_head();
		dat->set_size(nx2, ny, nz);
		//dat->to_zero();  // do not need it, real_to_complex will do it right anyway
		if (offset == 1) dat->set_fftodd(true);
		else             dat->set_fftodd(false);

		float *d = dat->get_data();
		//std::cout<<" do_fft "<<rdata[5]<<"  "<<d[5]<<std::endl;
		EMfft::real_to_complex_nd(get_data(), d, nxreal, ny, nz);

		dat->update();
		dat->set_fftpad(true);
		dat->set_complex(true);
		if(dat->get_ysize()==1 && dat->get_zsize()==1) dat->set_complex_x(true);
		dat->set_ri(true);

		EXITFUNC;
#ifdef FFT_CACHING
//		printf("%p %d\n",this,nxyz);
		if (nxyz<80000000) {
			fftcache=dat->copy();
		}
#endif //FFT_CACHING
		return dat;
	}
}

void EMData::do_fft_inplace()
{
	ENTERFUNC;

	if ( is_complex() ) {
		LOGERR("real image expected. Input image is complex image.");
		throw ImageFormatException("real image expected. Input image is complex image.");
	}

	size_t offset;
	int nxreal;
	get_data(); // Required call if GPU caching is being used. Otherwise harmless
	if (!is_fftpadded()) {
		// need to extend the matrix along x
		// meaning nx is the un-fftpadded size
		nxreal = nx;
		offset = 2 - nx%2;
		if (1 == offset) set_fftodd(true);
		else             set_fftodd(false);
		int nxnew = nx + offset;
		set_size(nxnew, ny, nz);
		for (int iz = nz-1; iz >= 0; iz--) {
			for (int iy = ny-1; iy >= 0; iy--) {
				for (int ix = nxreal-1; ix >= 0; ix--) {
					size_t oldxpos = ix + (iy + iz*ny)*(size_t)nxreal;
					size_t newxpos = ix + (iy + iz*ny)*(size_t)nxnew;
					(*this)(newxpos) = (*this)(oldxpos);
				}
			}
		}
		set_fftpad(true);
	} else {
		offset = is_fftodd() ? 1 : 2;
		nxreal = nx - offset;
	}
	EMfft::real_to_complex_nd(rdata, rdata, nxreal, ny, nz);

	set_complex(true);
	if(ny==1 && nz==1)  set_complex_x(true);
	set_ri(true);

	update();

	EXITFUNC;
	//return this;
}

#ifdef EMAN2_USING_CUDA

#include "cuda/cuda_emfft.h"

EMData *EMData::do_fft_cuda()
{
	ENTERFUNC;

	if ( is_complex() ) {
		LOGERR("real image expected. Input image is complex image.");
		throw ImageFormatException("real image expected. Input image is complex image.");
	}

	int offset = 2 - nx%2;
	EMData* dat = new EMData(0,0,nx+offset,ny,nz,attr_dict);
	if(!dat->rw_alloc()) throw UnexpectedBehaviorException("Bad alloc");
	//cout << "Doing CUDA FFT " << cudarwdata << endl;
	if(cudarwdata == 0){copy_to_cuda();}
	cuda_dd_fft_real_to_complex_nd(cudarwdata, dat->cudarwdata, nx, ny,nz, 1);

	if (offset == 1) dat->set_fftodd(true);
	else             dat->set_fftodd(false);

	dat->set_fftpad(true);
	dat->set_complex(true);
	if(dat->get_ysize()==1 && dat->get_zsize()==1) dat->set_complex_x(true);
	dat->set_ri(true);
	dat->update();

	EXITFUNC;
	return dat;
}

void EMData::do_fft_inplace_cuda()
{
	ENTERFUNC;

	if ( is_complex() ) {
		LOGERR("real image expected. Input image is complex image.");
		throw ImageFormatException("real image expected. Input image is complex image.");
	}

	int offset = 2 - nx%2;
	float* tempcudadata = 0;
	cudaError_t error = cudaMalloc((void**)&tempcudadata,(nx + offset)*ny*nz*sizeof(float));
	if( error != cudaSuccess) throw ImageFormatException("Couldn't allocate memory.");
	
	//cout << "Doing CUDA FFT inplace" << cudarwdata << endl;
	if(cudarwdata == 0){copy_to_cuda();}
	cuda_dd_fft_real_to_complex_nd(cudarwdata, tempcudadata, nx, ny,nz, 1);
	// this section is a bit slight of hand it actually does the FFT out of place but this avoids and EMData object creation and detruction...
	cudaError_t ferror = cudaFree(cudarwdata);
	if ( ferror != cudaSuccess) throw UnexpectedBehaviorException( "CudaFree failed:" + string(cudaGetErrorString(error)));
	cudarwdata = tempcudadata;
	num_bytes = (nx + offset)*ny*nz*sizeof(float);

	if (offset == 1) set_fftodd(true);
	else             set_fftodd(false);

	nx = nx + offset; // don't want to call set_size b/c that will delete my cudadata, remember what I am doing is a bit slignt of hand....
	set_fftpad(true);
	set_complex(true);
	if(get_ysize()==1 && get_zsize()==1) set_complex_x(true);
	set_ri(true);
	update();

	EXITFUNC;
//	return this;
}

EMData *EMData::do_ift_cuda()
{
	ENTERFUNC;

	if (!is_complex()) {
		LOGERR("complex image expected. Input image is real image.");
		throw ImageFormatException("complex image expected. Input image is real image.");
	}

	if (!is_ri()) {
		throw ImageFormatException("complex ri expected. Got amplitude/phase.");
	}

	int offset = is_fftodd() ? 1 : 2;
	EMData* dat = new EMData(0,0,nx-offset,ny,nz,attr_dict);
	if(!dat->rw_alloc()) throw UnexpectedBehaviorException("Bad alloc");
	
	if(cudarwdata == 0){copy_to_cuda();}

	
	int ndim = get_ndim();
	if ( ndim == 1 ) {
		cuda_dd_fft_complex_to_real_nd(cudarwdata,dat->cudarwdata, nx-offset,1,1,1);
	} else if (ndim == 2) {
		cuda_dd_fft_complex_to_real_nd(cudarwdata,dat->cudarwdata, ny,nx-offset,1,1);
	} else if (ndim == 3) {
		cuda_dd_fft_complex_to_real_nd(cudarwdata,dat->cudarwdata, nz,ny,nx-offset,1);
	} else throw ImageDimensionException("No cuda FFT support of images with dimensions exceeding 3");
	
	// SCALE the inverse FFT
	float scale = 1.0f/static_cast<float>((dat->get_size()));
	dat->mult(scale); 

	dat->set_fftpad(false);
	dat->set_fftodd(false);
	dat->set_complex(false);
	if(dat->get_ysize()==1 && dat->get_zsize()==1)  dat->set_complex_x(false);
	dat->set_ri(false);
//	dat->gpu_update();
	dat->update(); 
	
	EXITFUNC;
	return dat;
}

/*
   FFT in place does not depad, hence this routine is of limited use b/c mem operations on the device are quite SLOW, JFF
   use
*/

void EMData::do_ift_inplace_cuda()
{
	ENTERFUNC;

	if (!is_complex()) {
		LOGERR("complex image expected. Input image is real image.");
		throw ImageFormatException("complex image expected. Input image is real image.");
	}

	if (!is_ri()) {
		LOGWARN("run IFT on AP data, only RI should be used. ");
	}

	int offset = is_fftodd() ? 1 : 2;
	
	if(cudarwdata == 0){copy_to_cuda();}
	
	int ndim = get_ndim();
	if ( ndim == 1 ) {
		cuda_dd_fft_complex_to_real_nd(cudarwdata,cudarwdata, nx-offset,1,1,1);
	} else if (ndim == 2) {
		cuda_dd_fft_complex_to_real_nd(cudarwdata,cudarwdata, ny,nx-offset,1,1);
	} else if (ndim == 3) {
		cuda_dd_fft_complex_to_real_nd(cudarwdata,cudarwdata, nz,ny,nx-offset,1);
	} else throw ImageDimensionException("No cuda FFT support of images with dimensions exceeding 3");
#if defined	FFTW2 || defined FFTW3 //native fft and ACML already done normalization
	// SCALE the inverse FFT
	int nxo = nx - offset;
	float scale = 1.0f / (nxo * ny * nz);
	mult(scale); //if we are just doing a CCF, this is a waste!
#endif //FFTW2 || FFTW3

	set_fftpad(true);
	set_complex(false);

	if(ny==1 && nz==1) set_complex_x(false);
	set_ri(false);
	update();
	
	EXITFUNC;
//	return this;
}

#endif //EMAN2_USING_CUDA

EMData *EMData::do_ift()
{
	ENTERFUNC;

	if (!is_complex()) {
		LOGERR("complex image expected. Input image is real image.");
		throw ImageFormatException("complex image expected. Input image is real image.");
	}

	if (!is_ri()) {
		LOGWARN("run IFT on AP data, only RI should be used. Converting.");
	}

	get_data(); // Required call if GPU caching is being used. Otherwise harmless
	EMData* dat = copy_head();
	dat->set_size(nx, ny, nz);
	ap2ri();

	float *d = dat->get_data();
	int ndim = get_ndim();

	/* Do inplace IFT on a image copy, because the complex to real transform of
	 * nd will destroy its input array even for out-of-place transforms.
	 */
	memcpy((char *) d, (char *) rdata, (size_t)nx * ny * nz * sizeof(float));

	int offset = is_fftodd() ? 1 : 2;
	//cout << "Sending offset " << offset << " " << nx-offset << endl;
	if (ndim == 1) {
		EMfft::complex_to_real_nd(d, d, nx - offset, ny, nz);
	} else {
		EMfft::complex_to_real_nd(d, d, nx - offset, ny, nz);

		size_t row_size = (nx - offset) * sizeof(float);
		for (size_t i = 1; i < (size_t)ny * nz; i++) {
			memmove((char *) &d[i * (nx - offset)], (char *) &d[i * nx], row_size);
		}
	}

	dat->set_size(nx - offset, ny, nz);	//remove the padding
#if defined	FFTW2 || defined FFTW3 //native fft and ACML already done normalization
	// SCALE the inverse FFT
	float scale = 1.0f / ((nx - offset) * ny * nz);
	dat->mult(scale);
#endif	//FFTW2 || FFTW3
	dat->set_fftodd(false);
	dat->set_fftpad(false);
	dat->set_complex(false);
	if(dat->get_ysize()==1 && dat->get_zsize()==1)  dat->set_complex_x(false);
	dat->set_ri(false);
	dat->update();


	EXITFUNC;
	return dat;
}

/*
   FFT in place does not depad, return real x-extended image (needs to be depadded before use as PAP does in CCF routines)
   use
*/
void EMData::do_ift_inplace()
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

	int offset = is_fftodd() ? 1 : 2;
	float* data = get_data();
	EMfft::complex_to_real_nd(data, data, nx - offset, ny, nz);

#if defined	FFTW2 || defined FFTW3 	//native fft and ACML already done normalization
	// SCALE the inverse FFT
	int nxo = nx - offset;
	float scale = 1.0f / ((size_t)nxo * ny * nz);
	mult(scale);
#endif //FFTW2 || FFTW3

	set_fftpad(true);
	set_complex(false);
	if(ny==1 && nz==1) set_complex_x(false);
	set_ri(false);
	update();

	EXITFUNC;
//	return this;
}
#undef rdata


std::string EMData::render_ap24(int x0, int y0, int ixsize, int iysize,
						 int bpl, float scale, int mingray, int maxgray,
						 float render_min, float render_max,float gamma,int flags)
{
	ENTERFUNC;

	int asrgb;
	int hist=(flags&2)/2;
	int invy=(flags&4)?1:0;

	if (!is_complex()) throw ImageDimensionException("complex only");

	if (get_ndim() != 2) {
		throw ImageDimensionException("2D only");
	}

	if (is_complex()) ri2ap();

	if (render_max <= render_min) {
		render_max = render_min + 0.01f;
	}

	if (gamma<=0) gamma=1.0;

	// Calculating a full floating point gamma for
	// each pixel in the image slows rendering unacceptably
	// however, applying a gamma-mapping to an 8 bit colorspace
	// has unaccepable accuracy. So, we oversample the 8 bit colorspace
	// as a 12 bit colorspace and apply the gamma mapping to that
	// This should produce good accuracy for gamma values
	// larger than 0.5 (and a high upper limit)
	static int smg0=0,smg1=0;	// while this destroys threadsafety in the rendering process
	static float sgam=0;		// it is necessary for speed when rendering large numbers of small images
	static unsigned char gammamap[4096];
	if (gamma!=1.0 && (smg0!=mingray || smg1!=maxgray || sgam!=gamma)) {
		for (int i=0; i<4096; i++) {
			if (mingray<maxgray) gammamap[i]=(unsigned char)(mingray+(maxgray-mingray+0.999)*pow(((float)i/4096.0f),gamma));
			else gammamap[4095-i]=(unsigned char)(mingray+(maxgray-mingray+0.999)*pow(((float)i/4096.0f),gamma));
		}
	}
	smg0=mingray;	// so we don't recompute the map unless something changes
	smg1=maxgray;
	sgam=gamma;

	if (flags&8) asrgb=4;
	else if (flags&1) asrgb=3;
	else throw ImageDimensionException("must set flag 1 or 8");

	std::string ret=std::string();
//	ret.resize(iysize*bpl);
	ret.assign(iysize*bpl+hist*1024,char(mingray));
	unsigned char *data=(unsigned char *)ret.data();
	unsigned int *histd=(unsigned int *)(data+iysize*bpl);
	if (hist) {
		for (int i=0; i<256; i++) histd[i]=0;
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
	float gs2 = 4095.999f / (render_max - render_min);
//	float gs2 = 1.0 / (render_max - render_min);
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

	int mid=nx*ny/2;
	float* image_data = get_data();
	if (dsx != -1) {
		int l = y0 * nx;
		for (int j = ymax; j >= ymin; j--) {
			int ll = x0;
			for (int i = xmin; i < xsize; i++) {
				if (l + ll > lmax || ll >= nx - 2) break;

				int k = 0;
				unsigned char p;
				int ph;
				if (ll >= nx / 2) {
					if (l >= (ny - inv_scale) * nx) k = 2 * (ll - nx / 2) + 2;
					else k = 2 * (ll - nx / 2) + l + 2 + nx;
					if (k>=mid) k-=mid;		// These 2 lines handle the Fourier origin being in the corner, not the middle
					else k+=mid;
					ph = (int)(image_data[k+1]*768/(2.0*M_PI))+384;	// complex phase as integer 0-767
				}
				else {
					k = nx * ny - (l + 2 * ll) - 2;
					ph = (int)(-image_data[k+1]*768/(2.0*M_PI))+384;	// complex phase as integer 0-767
					if (k>=mid) k-=mid;		// These 2 lines handle the Fourier origin being in the corner, not the middle
					else k+=mid;
				}
				float t = image_data[k];
				if (t <= rm)  p = mingray;
				else if (t >= render_max) p = maxgray;
				else if (gamma!=1.0) {
					k=(int)(gs2 * (t-render_min));		// map float value to 0-4096 range
					p = gammamap[k];					// apply gamma using precomputed gamma map
				}
				else {
					p = (unsigned char) (gs * (t - render_min));
					p += mingray;
				}
				if (ph<256) {
					data[i * asrgb + j * bpl] = p*(255-ph)/256;
					data[i * asrgb + j * bpl+1] = p*ph/256;
					data[i * asrgb + j * bpl+2] = 0;
				}
				else if (ph<512) {
					data[i * asrgb + j * bpl+1] = p*(511-ph)/256;
					data[i * asrgb + j * bpl+2] = p*(ph-256)/256;
					data[i * asrgb + j * bpl] = 0;
				}
				else {
					data[i * asrgb + j * bpl+2] = p*(767-ph)/256;
					data[i * asrgb + j * bpl] = p*(ph-512)/256;
					data[i * asrgb + j * bpl+1] = 0;
				}
				if (hist) histd[p]++;
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
				unsigned char p;
				int ph;
				if (ll >= nx / 2) {
					if (l >= (ny * nx - nx)) k = 2 * (ll - nx / 2) + 2;
					else k = 2 * (ll - nx / 2) + l + 2 + nx;
					if (k>=mid) k-=mid;		// These 2 lines handle the Fourier origin being in the corner, not the middle
					else k+=mid;
					ph = (int)(image_data[k+1]*768/(2.0*M_PI))+384;	// complex phase as integer 0-767
				}
				else {
					k = nx * ny - (l + 2 * ll) - 2;
					if (k>=mid) k-=mid;		// These 2 lines handle the Fourier origin being in the corner, not the middle
					else k+=mid;
					ph = (int)(-image_data[k+1]*768/(2.0*M_PI))+384;	// complex phase as integer 0-767
				}

				float t = image_data[k];
				if (t <= rm)
					p = mingray;
				else if (t >= render_max) {
					p = maxgray;
				}
				else if (gamma!=1.0) {
					k=(int)(gs2 * (t-render_min));		// map float value to 0-4096 range
					p = gammamap[k];					// apply gamma using precomputed gamma map
				}
				else {
					p = (unsigned char) (gs * (t - render_min));
					p += mingray;
				}
				if (ph<256) {
					data[i * asrgb + j * bpl] = p*(255-ph)/256;
					data[i * asrgb + j * bpl+1] = p*ph/256;
					data[i * asrgb + j * bpl+2] = 0;
				}
				else if (ph<512) {
					data[i * asrgb + j * bpl+1] = p*(511-ph)/256;
					data[i * asrgb + j * bpl+2] = p*(ph-256)/256;
					data[i * asrgb + j * bpl] = 0;
				}
				else {
					data[i * asrgb + j * bpl+2] = p*(767-ph)/256;
					data[i * asrgb + j * bpl] = p*(ph-512)/256;
					data[i * asrgb + j * bpl+1] = 0;
				}
				if (hist) histd[p]++;
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

	// this replicates r -> g,b
	if (asrgb==4) {
		for (int j=ymin*bpl; j<=ymax*bpl; j+=bpl) {
			for (int i=xmin; i<xsize*4; i+=4) {
				data[i+j+3]=255;
			}
		}
	}

	EXITFUNC;

	// ok, ok, not the most efficient place to do this, but it works
	if (invy) {
		int x,y;
		char swp;
		for (y=0; y<iysize/2; y++) {
			for (x=0; x<ixsize; x++) {
				swp=ret[y*bpl+x];
				ret[y*bpl+x]=ret[(iysize-y-1)*bpl+x];
				ret[(iysize-y-1)*bpl+x]=swp;
			}
		}
	}

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
	if ( iysize * inv_scale > ny) {
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
	float* image_data = get_data();
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
					float t = image_data[kk];
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
					float t = image_data[kk];
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
					float t = image_data[l];
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
					float t = image_data[l];
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

	EXITFUNC;
}

void EMData::ap2ri()
{
	ENTERFUNC;

	if (!is_complex() || is_ri()) {
		return;
	}
	
	Util::ap2ri(get_data(), (size_t)nx * ny * nz);
	set_ri(true);
	update();
	EXITFUNC;
}

void EMData::ri2inten()
{
	ENTERFUNC;

	if (!is_complex()) return;
	if (!is_ri()) ap2ri();

	float * data = get_data();
	size_t size = (size_t)nx * ny * nz;
	for (size_t i = 0; i < size; i += 2) {
		data[i]=data[i]*data[i]+data[i+1]*data[i+1];
		data[i+1]=0;
	}

	set_attr("is_intensity", int(1));
	update();
	EXITFUNC;
}


void EMData::ri2ap()
{
	ENTERFUNC;

	if (!is_complex() || !is_ri()) {
		return;
	}

	float * data = get_data();

	size_t size = (size_t)nx * ny * nz;
	for (size_t i = 0; i < size; i += 2) {
#ifdef	_WIN32
		float f = (float)_hypot(data[i], data[i + 1]);
#else
		float f = (float)hypot(data[i], data[i + 1]);
#endif
		if (data[i] == 0 && data[i + 1] == 0) {
			data[i + 1] = 0;
		}
		else {
			data[i + 1] = atan2(data[i + 1], data[i]);
		}
		data[i] = f;
	}

	set_ri(false);
	update();
	EXITFUNC;
}


float calc_bessel(const int n, const float& x) {
	gsl_sf_result result;
//	int success =
	gsl_sf_bessel_Jn_e(n,(double)x, &result);
	return (float)result.val;
}

EMData*   EMData::bispecRotTransInvN(int N, int NK)
{

	int EndP = this -> get_xsize(); // length(fTrueVec);
	int Mid  = (int) ((1+EndP)/2);
	int End = 2*Mid-1;

        int CountxyMax = End*End;

	int   *SortfkInds       = new    int[CountxyMax];
	int   *kVecX            = new    int[CountxyMax];
	int   *kVecY            = new    int[CountxyMax];
	float *fkVecR           = new  float[CountxyMax];
	float *fkVecI           = new  float[CountxyMax];
	float *absD1fkVec       = new  float[CountxyMax];
	float *absD1fkVecSorted = new  float[CountxyMax];

	float *jxjyatan2         = new  float[End*End]; //  jxjyatan2[jy*End + jx]  = atan2(jy+1-Mid , jx +1 -Mid);

	EMData * ThisCopy = new EMData(End,End);

	for (int jx=0; jx <End ; jx++) {
		for (int jy=0; jy <End ; jy++) {
			float ValNow = this -> get_value_at(jx,jy);
			ThisCopy -> set_value_at(jx,jy,ValNow);
//		cout<< " jxM= " << jx+1<<" jyM= " << jy+1<< "ValNow" << ValNow << endl; //    Works
	}}


	EMData* fk = ThisCopy -> do_fft();
	fk          ->process_inplace("xform.fourierorigin.tocenter");

//	EMData* fk
	EMData* fkRCopy = new EMData(End,End);
	EMData* fkICopy = new EMData(End,End);
	EMData* fkCopy  = new EMData(End,End);


	for  (int jCount= 0; jCount<End*End; jCount++) {
//		jCount = jy*End + jx;
		int jx             = jCount%End ;
		int jy             = (jCount-jx)/End ;
		jxjyatan2[jCount]  = atan2((float)(jy+1-Mid) , (float)(jx +1-Mid));
	}


	for (int kEx= 0; kEx<2*Mid; kEx=kEx+2) { // kEx twice the value of the Fourier
						// x variable: EMAN index for real, imag
		int kx    = kEx/2;		// kx  is  the value of the Fourier variable
		int kIx   = kx+Mid-1; // This is the value of the index for a matlab image (-1)
		int kCx   =  -kx ;
		int kCIx  = kCx+ Mid-1 ;
		for (int kEy= 0 ; kEy<End; kEy++) { // This is the value of the EMAN index
    		int kIy              =  kEy       ; //  This is the value of the index for a matlab image (-1)
			int ky               =  kEy+1-Mid; // (kEy+ Mid-1)%End - Mid+1 ;  // This is the actual value of the Fourier variable
			float realVal        =  fk -> get_value_at(kEx  ,kEy) ;
			float imagVal        =  fk -> get_value_at(kEx+1,kEy) ;
			float absVal         =  ::sqrt(realVal*realVal+imagVal*imagVal);
			float fkAng 	      =  atan2(imagVal,realVal);

			float NewRealVal   ;
			float NewImagVal   ;
			float AngMatlab    ;

			if (kIx==Mid-1) {
//				AngMatlab = -fkAng - 2.*M_PI*(kIy+ 1-Mid)*(Mid)/End;
//			cout<< "i= " << i << " kIx= " << kIx << " kIy=" << kIy << " fkVecR[i] =" << fkVecR[i]<< " fkVecI[i]="  << fkVecI[i] <<"  angle[i]= "  << AngMatlab << endl;
			}

			if (kIx>Mid-1){
//			cout<< "i= " << i << " kIx= " << kIx << " kIy=" << kIy << " fkVecR[i] =" << fkVecR[i]<< " fkVecI[i]="  << fkVecI[i] <<"  angle[i]= "  << AngMatlab << endl;
			}

			AngMatlab = fkAng - 2.0f*M_PI*(kx +ky)*(Mid)/End;
			NewRealVal  =   absVal*cos(AngMatlab);
			NewImagVal  =   absVal*sin(AngMatlab);


			fkVecR[kIy+kIx *End] =  NewRealVal ;
			fkVecR[kIy+kCIx*End] =  NewRealVal ;
			fkVecI[kIy+kIx *End] =  NewImagVal ;
			fkVecI[kIy+kCIx*End] = -NewImagVal ;
        	absD1fkVec[kIy + kIx  *End] = absVal;
        	absD1fkVec[kIy + kCIx *End] = absVal;
			kVecX[kIy+kIx  *End] =  kx      ;
        	kVecX[kIy+kCIx *End] =  kCx    ;
			kVecY[kIy+kIx  *End] =  ky     ;
			kVecY[kIy+kCIx *End] =  ky     ;
//			printf("kx=%d,ky=%d,tempVal =%f+ i %4.2f \n",kx,ky,realVal,imagVal );
//			cout << "kx = " << kx << "; ky = "<< ky << "; val is" << realVal<<"+ i "<<imagVal<< endl;

//			cout << "kIMx = "<< kIx+1 << "; kIMy = "<< kIy+1 <<"; fkAng*9/ 2pi is " << fkAng*9/2/M_PI<<  endl;
//			cout << "kIMx = "<< kIx+1 << "; kIMy = "<< kIy+1 <<"; absval is " << absVal<<  "; realval is " << NewRealVal<< "; imagval is " << NewImagVal<< endl;
			fkCopy  -> set_value_at(kIx ,kIy, absVal);
			fkCopy  -> set_value_at(kCIx,kIy, absVal);
			fkRCopy -> set_value_at(kIx, kIy, NewRealVal);
			fkRCopy -> set_value_at(kCIx,kIy, NewRealVal);
			fkICopy -> set_value_at(kIx, kIy, NewImagVal);
			fkICopy -> set_value_at(kCIx,kIy,-NewImagVal);

		}
	}
#ifdef _WIN32
	_unlink("fkCopy.???");
	_unlink("fk?Copy.???");
#else
	int rslt;
	rslt = system("rm -f fkCopy.???");  rslt++;
	rslt = system("rm -f fk?Copy.???"); rslt++;
#endif	//_WIN32
	fkCopy  -> write_image("fkCopy.img");
	fkRCopy -> write_image("fkRCopy.img");
	fkICopy -> write_image("fkICopy.img");

	cout << "Starting the sort "<< endl;

	vector< pair<float, int> > absInds;
	for(int i  = 0; i < CountxyMax; ++i ) {
		pair<float,int> p;
		p = make_pair(absD1fkVec[i],i);	// p = make_pair(rand(),i);
		absInds.push_back( p);
	}

	std::sort(absInds.begin(),absInds.end());

	for(int i  = 0; i < CountxyMax; ++i ) {
		pair<float,int> p   ;
		p = absInds[i]         ;
		absD1fkVecSorted[CountxyMax-1-i] =  p.first ;
		SortfkInds[CountxyMax-1-i]       =  p.second ;
	}

	cout << "Ending the sort "<< endl;

// float AngsMat[] ={2.8448, -0.3677,-0.2801,-1.0494,-1.7836,-2.5179, 2.9959, 3.0835,-0.1290,-0.8876,2.1829, 2.2705,1.5011,0.7669,0.0327,-0.7366,-0.6489,2.4215,-1.6029,1.4676,1.5552,0.7859,0.0517,-0.6825,-1.4518,-1.3642,1.7063,-1.7845,1.2859,1.3736,0.6043,-0.1299,-0.8642,-1.6335,-1.5459,1.5247,-1.6546,1.4159,1.5036,0.7342,0,-0.7342,-1.5036,-1.4159,1.6546,-1.5247,1.5459,1.6335,0.8642,0.1299,-0.6043,-1.3736,-1.286,1.7846,-1.7063,1.3642,1.4519,0.6825,-0.0517,-0.7859,-1.5553,-1.4676,1.6029,-2.4216,0.649,0.7366,-0.0327,-0.767,-1.5012,-2.2705,-2.1829,0.8877,0.1291,-3.0836,-2.9959,2.5179,1.7837,1.0495,0.2801,0.3677,-2.8449};


 	for(int i  = 0; i < CountxyMax; ++i ) {  // creates a new fkVec
//		int Si  = SortfkInds[i];
//		int kIx = (int)  Si/End;  kIx = (int)  i/End; // i = kIx*End+kIy
//		int kIy = Si  - kIx*End;  kIy = i  - kIx*End;
//		int iC = (End-1-kIx)*End + (End-1-kIy);
//		if (i<30) { cout<< "i= " << i << " kIx= " << kIx << " kIy=" << kIy << " valAft=" << absD1fkVecSorted[i]<< " valBef="  <<     absD1fkVec[Si] << "  SortfkInds = " << Si <<endl; }// This worked
//		cout<< "i= " << i << " kIx= " << kIx << " kIy=" << kIy << " fkVecR[i] =" << fkVecR[i]<< " fkVecI[i]="  << fkVecI[i] <<"  angle[i]= "  << fkAng << endl;
 	}
	cout<< "Ratio of Last Amplitude to First Amplitude= " << absD1fkVecSorted[NK] /absD1fkVecSorted[0]  << endl;

//	pause;

// 	for(int i  = 0; i < NK; ++i ) { // Prints out the new fkVec ,  CountxyMax
//		int Si= SortfkInds[i];
//		int kIx = (int)  Si/End; // i = kIx*End+kIy
//		int kIy = Si  - kIx*End;
//		cout << " kIxM= " << kIx+1 << " kIyM=" << kIy+1 << " fkVecAbs=" << ::sqrt(fkVecR[Si]*fkVecR[Si] +  fkVecI[Si]* fkVecI[Si]) << " fkVecAbs=" << absD1fkVecSorted[i] << " kx= " << kVecX[Si] <<  " ky=" << kVecY[Si] <<  endl;
// 	}

//       angEMAN+angMat+angDiff    =0  mod 2 pi

// 	angDiff=  2*pi*((-4):4)*(Mid)/End; angEMAN+angMat+angDiff= integer*2 *pi
//		[  absD1fkVecSorted, SortfkInds] =sort( absD1fkVec,'descend') ;
//	Util::sort_mat(&absD1fkVec[0],&absD1fkVec[Countxy],&SortfkInds[0],&SortfkInds[Countxy]);


//      Let radial sampling be 0:0.5:(Mid-1)

 //	int NK=  min(12,CountxyMax) ;



	cout << "NK = " << NK << endl;
	float frR= 3.0/4.0;
	int LradRange= (int) (floor(Mid/frR)) ;

        float *radRange = new float[LradRange]; //= 0:.75:(Mid-1);
	radRange[0]=0;
	for (int irad=1; irad < LradRange; irad++){
			radRange[irad] = radRange[irad-1] + frR; }



         // should equal to (2*Mid-1)
	cout << "Starting the calculation of invariants for N= " << N << endl;

/*	int NMax=5;            */

	EMData* RotTransInv = new EMData();
	RotTransInv -> set_size(LradRange,LradRange);


//	float  *RotTransInv       = new float[LradRange*LradRange ] ;
//	float  *RotTransInvN      = new float[LradRange*LradRange*(NMax+1) ] ;

//	for (int N=0 ; N<NMax; N++) {

	for (int jr1=0; jr1 < LradRange ; jr1++ ) { // LradRange
		float r1= radRange[jr1];
//		cout << "Pre jr2 "<< endl;
		for (int jr2=0;  jr2<LradRange;  jr2++ ) { //LradRange
			float r2= radRange[jr2];
			float RotTransInvTemp=0;
			for (int jCountkxy =0; jCountkxy<NK; jCountkxy++){
				int Countkxy =SortfkInds[jCountkxy] ;   // 1: CountxyMax
				int kx = kVecX[Countkxy] ;
				int ky = kVecY[Countkxy] ;
				float k2 = (float)(kx*kx+ky*ky);
				if (k2==0) { continue;}
				float phiK =0;
				if (k2>0) phiK= jxjyatan2[ (ky+Mid-1)*End + kx+Mid-1];  phiK= atan2((float)ky,(float)kx);

				float fkR     = fkVecR[Countkxy] ;
				float fkI     = fkVecI[Countkxy]  ;
/*				printf("jCountkxy=%d, Countkxy=%d,absD1fkVec(Countkxy)=%f,\t\t kx=%d, ky=%d \n", jCountkxy, Countkxy, absD1fkVec[Countkxy], kx, ky);*/

				for (int jCountqxy =0; jCountqxy<NK; jCountqxy++){
					int Countqxy =SortfkInds[jCountqxy] ;   // Countqxy is the index for absD1fkVec
					int qx   = kVecX[Countqxy] ;
					int qy   = kVecY[Countqxy] ;
					int q2   = qx*qx+qy*qy;
					if (q2==0) {continue;}
					float phiQ =0;
					if (q2>0) phiQ = jxjyatan2[ (qy+Mid-1)*End + qx+Mid-1];   phiQ=atan2((float)qy,(float)qx);
					float fqR     = fkVecR[Countqxy]  ;
					float fqI     = fkVecI[Countqxy]  ;
					int kCx  = (-kx-qx);
					int kCy  = (-ky-qy);
					int kCIx = ((kCx+Mid+2*End)%End);// labels of the image in C
					int kCIy = ((kCy+Mid+2*End)%End);
					kCx  = kCIx-Mid; // correct
					kCy  = kCIy-Mid; // correct
					int CountCxy = kCIx*End+kCIy;
					float fCR     = fkVecR[CountCxy];
					float fCI     = fkVecI[CountCxy];
					if (jr1+jr2==-1) {
					printf("jCountqxy=%d , Countqxy=%d, absD1fkVec(Countqxy)=%f,qx=%d, qy=%d \n", jCountqxy, Countqxy, absD1fkVec[Countqxy],qx, qy);
					printf(" CountCxy=%d,absD1fkVec[CountCxy]=%f,  kCx=%d,     kCy=%d \n",CountCxy, absD1fkVec[CountCxy], kCx, kCy );
					}
					for (int p=0; p<NK; p++){
//						printf("p=%d, SortfkInds[p]=%d, CountCxy =%d \n", p,SortfkInds[p], CountCxy);
						if (SortfkInds[p]==CountCxy){
							float Arg1 = 2.0f*M_PI*r1*::sqrt((float) q2)/End;
							float Arg2 = 2.0f*M_PI*r2*::sqrt((float) k2)/End;
//							printf("Arg1=%4.2f, Arg2=%4.2f,  \n",Arg1, Arg2 );
//							if (Arg1+ Arg2<15) {
								float bispectemp  = (fkR*(fqR*fCR -fqI*fCI) -fkI*(fqI*fCR  +fqR*fCI))
								* cos(N*(phiK-phiQ+M_PI));
								bispectemp  -= (fkR*(fqR*fCI + fqI*fCR) +fkI*(fqR*fCR - fqI*fCI))
								* sin(N*(phiK-phiQ+M_PI));
								float bess1 = calc_bessel(N, Arg1 );
								float bess2 = calc_bessel(N, Arg2 );
//			printf("fkr=%4.2f, fqr=%4.2f, bess1=%4.2f,bess2=%4.2f \n",fkR, fqR, bess1, bess2);
/*			printf("p =%d, SortfkInds[p]=%d, CountCxy=%d, Arg1 =%4.2f, bess1=%4.2f,  \n",
				p, SortfkInds[p],CountCxy, Arg1, bess1);*/
								RotTransInvTemp   = RotTransInvTemp  + bispectemp  * bess1*bess2 ;
//							}
						}
					}
				} // jCountqxy
			} // jCountkxy
			RotTransInv -> set_value_at(jr1,jr2, RotTransInvTemp)   ;
/*		RotTransInvN[jr1 + LradRange*jr2+LradRange*LradRange*N] = RotTransInvTemp  ;*/
		} //jr2
	} //jr1
// }//N

	return  RotTransInv ;


}



/*
// find example
#include <iostream>
#include <algorithm>
#include <vector>
using namespace std;

int main () {
  int myints[] = { 10, 20, 30 ,40 };
  int * p;

  // pointer to array element:
  p = find(myints,myints+4,30);
  ++p;
  cout << "The element following 30 is " << *p << endl;

  vector<int> myvector (myints,myints+4);
  vector<int>::iterator it;

  // iterator to vector element:
  it = find (myvector.begin(), myvector.end(), 30);
  ++it;
  cout << "The element following 30 is " << *it << endl;

  return 0;
}*/

EMData*   EMData::bispecRotTransInvDirect(int type)
{

	int EndP = this -> get_xsize(); // length(fTrueVec);
	int Mid  = (int) ((1+EndP)/2);
	int End = 2*Mid-1;

        int CountxyMax = End*End;

//	int   *SortfkInds       = new    int[CountxyMax];
	int   *kVecX            = new    int[CountxyMax];
	int   *kVecY            = new    int[CountxyMax];
	float *fkVecR           = new  float[CountxyMax];
	float *fkVecI           = new  float[CountxyMax];
	float *absD1fkVec       = new  float[CountxyMax];
//	float *absD1fkVecSorted = new  float[CountxyMax];


	float *jxjyatan2         = new  float[End*End];


	EMData * ThisCopy = new EMData(End,End);

	for (int jx=0; jx <End ; jx++) {  // create jxjyatan2
		for (int jy=0; jy <End ; jy++) {
			float ValNow = this -> get_value_at(jx,jy);
			ThisCopy -> set_value_at(jx,jy,ValNow);
			jxjyatan2[jy*End + jx]  = atan2((float)(jy+1-Mid) , (float)(jx +1 -Mid));
//		cout<< " jxM= " << jx+1<<" jyM= " << jy+1<< "ValNow" << ValNow << endl; //    Works
	}}


	EMData* fk = ThisCopy -> do_fft();
	fk          ->process_inplace("xform.fourierorigin.tocenter");

//	Create kVecX , kVecy etc

	for (int kEx= 0; kEx<2*Mid; kEx=kEx+2) { // kEx twice the value of the Fourier
						// x variable: EMAN index for real, imag
		int kx    = kEx/2;		// kx  is  the value of the Fourier variable
	        int kIx   = kx+Mid-1; // This is the value of the index for a matlab image (-1)
		int kCx   = -kx ;
		int kCIx  = kCx+ Mid-1 ;
		for (int kEy= 0 ; kEy<End; kEy++) { // This is the value of the EMAN index
    		 	int kIy              =  kEy       ; //  This is the value of the index for a matlab image (-1)
			int ky               =  kEy+1-Mid; // (kEy+ Mid-1)%End - Mid+1 ;  // This is the actual value of the Fourier variable
			float realVal        =  fk -> get_value_at(kEx  ,kEy) ;
			float imagVal        =  fk -> get_value_at(kEx+1,kEy) ;
			float absVal         =  ::sqrt(realVal*realVal+imagVal*imagVal);
			float fkAng 	     =  atan2(imagVal,realVal);

			float NewRealVal   ;
			float NewImagVal   ;
			float AngMatlab    ;

			if (kIx==Mid-1) {
//				AngMatlab = -fkAng - 2.*M_PI*(kIy+ 1-Mid)*(Mid)/End;
			}

			if (kIx>Mid-1){
//			cout<< "i= " << i << " kIx= " << kIx << " kIy=" << kIy << " fkVecR[i] =" << fkVecR[i]<< " fkVecI[i]="  << fkVecI[i] <<"  angle[i]= "  << AngMatlab << endl;
			}

			AngMatlab = fkAng - 2.0f*M_PI*(kx +ky)*(Mid)/End;
			NewRealVal  =   absVal*cos(AngMatlab);
			NewImagVal  =   absVal*sin(AngMatlab);


			fkVecR[ kIy +kIx *End] =  NewRealVal ;
			fkVecR[(End-1-kIy)+kCIx*End] =  NewRealVal ;
			fkVecI[ kIy +kIx *End] =  NewImagVal ;
			fkVecI[(End-1-kIy)+kCIx*End] = -NewImagVal ;
        		absD1fkVec[(End-1-kIy) + kIx  *End] = absVal;
        		absD1fkVec[(End-1-kIy) + kCIx *End] = absVal;
			kVecX[kIy+kIx  *End] =  kx      ;
        		kVecX[kIy+kCIx *End] =  kCx    ;
			kVecY[kIy+kIx  *End] =  ky     ;
			kVecY[kIy+kCIx *End] =  ky     ;

 //			cout << " kIxM= " << kIx+1 << " kIy=" << kIy+1 << " fkVecR[i] =" << NewRealVal << " fkVecI[i]="  << NewImagVal <<"  angle[i]= "  << AngMatlab << " Total Index" << kIy+kIx *End << endl;

//			printf("kx=%d,ky=%d,tempVal =%f+ i %4.2f \n",kx,ky,realVal,imagVal );
//			cout << "kx = " << kx << "; ky = "<< ky << "; val is" << realVal<<"+ i "<<imagVal<< endl;

//			cout << "kIMx = "<< kIx+1 << "; kIMy = "<< kIy+1 <<"; fkAng*9/ 2pi is " << fkAng*9/2/M_PI<<  endl;
//			cout << "kIMx = "<< kIx+1 << "; kIMy = "<< kIy+1 <<"; absval is " << absVal<<  "; realval is " << NewRealVal<< "; imagval is " << NewImagVal<< endl;
		}
	}


//	for (int TotalInd = 0 ;  TotalInd < CountxyMax ; TotalInd++){
//	        int kx     = kVecX[TotalInd]; // This is the value of the index for a matlab image (-1)
//	        int kIx    = kx+Mid-1; // This is the value of the index for a matlab image (-1)
//	        int ky     = kVecY[TotalInd];
//	        int kIy    = ky+Mid-1; // This is the value of the index for a matlab image (-1)
		//float fkR  = fkVecR[kIy+kIx *End]  ;
		//float fkI  = fkVecI[kIy+kIx *End]  ;
//	}

	float frR= 3.0/4.0;
	frR= 1;
	int LradRange= (int) (1+floor(Mid/frR -.1)) ;

        float *radRange = new float[LradRange]; //= 0:.75:(Mid-1);
	for (int irad=0; irad < LradRange; irad++){
			radRange[irad] =  frR*irad;
//			cout << " irad = " << irad << " radRange[irad]= " <<  radRange[irad] <<  " LradRange= " << LradRange << endl;
	}

	cout << "Starting the calculation of invariants" << endl;


	if (type==0) {
		int LthetaRange  = 59;
		float ftR        = (2.0f*M_PI/LthetaRange );
		float *thetaRange = new float[LthetaRange]; //= 0:.75:(Mid-1);

		for (int ith=0; ith < LthetaRange; ith++){
				thetaRange[ith] =  ftR*ith; }

		int TotalVol = LradRange*LradRange*LthetaRange;

		float *RotTransInv   = new  float[TotalVol];
		float *WeightInv     = new  float[TotalVol];

		for (int jW=0; jW<TotalVol; jW++) {
			RotTransInv[jW] = 0;
			WeightInv[jW]   = 0;
		}

		for (int jW=0; jW<TotalVol; jW++) {
			RotTransInv[jW] = 0;
			WeightInv[jW]   = 0;
		}
	//	float  *RotTransInv       = new float[LradRange*LradRange ] ;
	//	float  *RotTransInvN      = new float[LradRange*LradRange*(NMax+1) ] ;

		for (int Countkxy =0; Countkxy<CountxyMax; Countkxy++){  // Main Section for type 0
			int kx = kVecX[Countkxy] ;
			int ky = kVecY[Countkxy] ;
			float k2 = ::sqrt((float)(kx*kx+ky*ky));
			float phiK =0;
			if (k2>0)    phiK = jxjyatan2[ (ky+Mid-1)*End + kx+Mid-1]; //  phiK=atan2(ky,kx);
			float fkR     = fkVecR[(ky+Mid-1) + (kx+Mid-1) *End] ;
			float fkI     = fkVecI[(ky+Mid-1) + (kx+Mid-1) *End]  ;
	//		printf("Countkxy=%d,\t kx=%d, ky=%d, fkR=%3.2f,fkI=%3.2f \n", Countkxy, kx, ky, fkR, fkI);

			if ((k2==0)|| (k2>Mid) ) { continue;}

			for (int Countqxy =0; Countqxy<CountxyMax; Countqxy++){   // This is the innermost loop
				int qx   = kVecX[Countqxy] ;
				int qy   = kVecY[Countqxy] ;
				float q2   = ::sqrt((float)(qx*qx+qy*qy));
				if ((q2==0)|| (q2>Mid) ) {continue;}
				float phiQ =0;
				if (q2>0)   phiQ = jxjyatan2[ (qy+Mid-1)*End + qx+Mid-1]; // phiQ=atan2(qy,qx);
				float fqR     = fkVecR[(qy+Mid-1) + (qx+Mid-1) *End] ;
				float fqI     = fkVecI[(qy+Mid-1) + (qx+Mid-1) *End]  ;
				int kCx  = (-kx-qx);
				int kCy  = (-ky-qy);
				int kCIx = ((kCx+Mid+2*End)%End);// labels of the image in C
				int kCIy = ((kCy+Mid+2*End)%End);
				kCx  = ((kCIx+End-1)%End)+1-Mid; // correct
				kCy  = ((kCIy+End-1)%End)+1-Mid ; // correct

//				float C2   = ::sqrt((float)(kCx*kCx+ kCy*kCy));
				int CountCxy  = (kCx+Mid-1)*End+(kCy+Mid-1);
				float fCR     = fkVecR[CountCxy];
				float fCI     = fkVecI[CountCxy];
	/*			if (Countkxy==1) {
					printf(" Countqxy=%d, absD1fkVec(Countqxy)=%f,qx=%d, qy=%d \n", Countqxy, absD1fkVec[Countqxy],qx, qy);
					printf(" CountCxy=%d, absD1fkVec[CountCxy]=%f,kCx=%d,kCy=%d \n",CountCxy, absD1fkVec[CountCxy], kCx, kCy );
				}*/
//				float   phiC = jxjyatan2[ (kCy+Mid-1)*End + kCx+Mid-1];
				float   phiQK = (4*M_PI+phiQ-phiK);
				while (phiQK> (2*M_PI)) phiQK -= (2*M_PI);



				float bispectemp  = (fkR*(fqR*fCR -fqI*fCI) -fkI*(fqI*fCR  +fqR*fCI));

				if  ((q2<k2) )  continue;
//				if  ((q2<k2) || (C2<k2) || (C2<q2))  continue;

	//				printf(" CountCxy=%d, absD1fkVec[CountCxy]=%f,kCx=%d,kCy=%d \n",CountCxy, absD1fkVec[CountCxy], kCx, kCy );

	//                      up to here, matched perfectly with Matlab

				int k2IndLo  = 0; while ((k2>=radRange[k2IndLo+1]) && (k2IndLo+1 < LradRange ) ) k2IndLo +=1;
				int k2IndHi = k2IndLo;
				float k2Lo= radRange[k2IndLo];
				if (k2IndLo+1< LradRange) {
					k2IndHi   = k2IndLo+1;
				}
//				float k2Hi= radRange[k2IndHi];

				float kCof =k2-k2Lo;

				int q2IndLo  = 0; while ((q2>=radRange[q2IndLo+1]) && (q2IndLo+1 < LradRange ) ) q2IndLo +=1;
				int q2IndHi=q2IndLo;
				float q2Lo= radRange[q2IndLo];
				if (q2IndLo+1 < LradRange)  {
					q2IndHi   = q2IndLo+1 ;
				}
				float qCof = q2-q2Lo;

				if ((qCof<0) || (qCof >1) ) {
					cout<< "Weird! qCof="<< qCof <<  " q2="<< q2 << " q2IndLo="<< q2IndLo << endl ;
					int x    ;
					cin >> x ;
				}

				int thetaIndLo = 0; while ((phiQK>=thetaRange[thetaIndLo+1])&& (thetaIndLo+1<LthetaRange)) thetaIndLo +=1;
				int thetaIndHi = thetaIndLo;

				float thetaLo  = thetaRange[thetaIndLo];
				float thetaHi = thetaLo;
				float thetaCof = 0;

				if (thetaIndLo+1< LthetaRange) {
					thetaIndHi = thetaIndLo +1;
				}else{
					thetaIndHi=0;
				}

				thetaHi    = thetaRange[thetaIndHi];

				if (thetaHi==thetaLo) {
					thetaCof =0 ;
				} else {
					thetaCof   = (phiQK-thetaLo)/(thetaHi-thetaLo);
				}

				if ((thetaCof>2*M_PI)  ) {
					cout<< "Weird! thetaCof="<< thetaCof <<endl ;
					thetaCof=0;
				}


	// 			if ((thetaIndLo>=58) || (k2IndLo >= LradRange-1) || (q2IndLo >= LradRange-1) ) {


				for (int jk =1; jk<=2; jk++){
				for (int jq =1; jq<=2; jq++){
				for (int jtheta =1; jtheta<=2; jtheta++){

					float Weight = (kCof+(1-2*kCof)*(jk==1))*(qCof+(1-2*qCof)*(jq==1))
							* (thetaCof+(1-2*thetaCof)*(jtheta==1));


					int k2Ind      =  k2IndLo*(jk==1)      +   k2IndHi*(jk==2);
					int q2Ind      =  q2IndLo*(jq==1)      +   q2IndHi*(jq==2);
					int thetaInd   =  thetaIndLo*(jtheta==1)  + thetaIndHi*(jtheta ==2);
					int TotalInd   = thetaInd*LradRange*LradRange+q2Ind*LradRange+k2Ind;
	/*				if (TotalInd+1 >=  LthetaRange*LradRange*LradRange) {
						cout << "Weird!!! TotalInd="<< TotalInd << " IndMax" << LthetaRange*LradRange*LradRange << " LradRange=" << LradRange << endl;
						cout << "k2Ind= "<< k2Ind  << " q2Ind="<< q2Ind  << " thetaInd="<< thetaInd  << " q2IndLo="<< q2IndLo  << " q2IndHi="<< q2IndHi  <<  endl;
						cout << "k2=" << k2 << "q2=" << q2 << " phiQK=" << phiQK*180.0/M_PI<< endl;
					}*/

					RotTransInv[TotalInd] += Weight*bispectemp;
					WeightInv[TotalInd]   +=  Weight;
	//				cout << "k2Ind= "<< k2Ind  << " q2Ind="<< q2Ind  << "Weight=" << Weight << endl;
				}}}
			} // Countqxy
		} // Countkxy

		cout << "Finished Main Section " << endl;

	/*		RotTransInvN[jr1 + LradRange*jr2+LradRange*LradRange*N] = RotTransInvTemp  ;*/

		cout << " LradRange " <<LradRange <<" LthetaRange " << LthetaRange << endl;
		EMData *RotTransInvF  = new  EMData(LradRange,LradRange,LthetaRange);
		EMData *WeightImage   = new  EMData(LradRange,LradRange,LthetaRange);

	// 	cout << "FFFFFFF" << endl;
	//
	// 	RotTransInvF -> set_size(LradRange,LradRange,LthetaRange);
	//
	// 	cout << "GGGG" << endl;

		for (int jtheta =0; jtheta < LthetaRange; jtheta++){	// write out to RotTransInvF
		for (int jq =0; jq<LradRange; jq++){ // LradRange
		for (int jk =0; jk<LradRange ; jk++){// LradRange
	//		cout << "Hi There" << endl;
			int TotalInd   = jtheta*LradRange*LradRange+jq*LradRange+jk;
			float Weight = WeightInv[TotalInd];
			WeightImage    -> set_value_at(jk,jq,jtheta,Weight);
			RotTransInvF   -> set_value_at(jk,jq,jtheta,0);
			if (Weight <= 0) continue;
			RotTransInvF -> set_value_at(jk,jq,jtheta,RotTransInv[TotalInd] / Weight);//  include /Weight
			int newjtheta = (LthetaRange- jtheta)%LthetaRange;
			RotTransInvF -> set_value_at(jq,jk,newjtheta,RotTransInv[TotalInd]/Weight );//  include /Weight
				}
			}
		}

		cout << " Almost Done " << endl;
#ifdef	_WIN32
		_unlink("WeightImage.???");
#else
		int rslt;
		rslt = system("rm -f WeightImage.???"); rslt++;
#endif	//_WIN32
		WeightImage  -> write_image("WeightImage.img");

		return  RotTransInvF ;
	} // Finish type 0

	if (type==1) {
		int TotalVol = LradRange*LradRange;

		float *RotTransInv   = new  float[TotalVol];
		float *WeightInv     = new  float[TotalVol];

		for (int jW=0; jW<TotalVol; jW++) {
			RotTransInv[jW] = 0;
			WeightInv[jW]   = 0;
		}


		for (int Countkxy =0; Countkxy<CountxyMax; Countkxy++){
			int kx = kVecX[Countkxy] ;
			int ky = kVecY[Countkxy] ;
			float k2 = ::sqrt((float)(kx*kx+ky*ky));
			float fkR     = fkVecR[(ky+Mid-1) + (kx+Mid-1) *End] ;
			float fkI     = fkVecI[(ky+Mid-1) + (kx+Mid-1) *End]  ;
	//		printf("Countkxy=%d,\t kx=%d, ky=%d, fkR=%3.2f,fkI=%3.2f \n", Countkxy, kx, ky, fkR, fkI);

			if ((k2==0)|| (k2>Mid) ) { continue;}

			for (int Countqxy =0; Countqxy<CountxyMax; Countqxy++){   // This is the innermost loop

//                      up to here, matched perfectly with Matlab
				int qx   = kVecX[Countqxy] ;
				int qy   = kVecY[Countqxy] ;
				float q2   = ::sqrt((float)(qx*qx+qy*qy));
				if ((q2==0)|| (q2>Mid) ) {continue;}
				if  ((q2<k2) )   continue;

				float fqR     = fkVecR[(qy+Mid-1) + (qx+Mid-1) *End] ;
				float fqI     = fkVecI[(qy+Mid-1) + (qx+Mid-1) *End]  ;

				int kCx  = (-kx-qx);
				int kCy  = (-ky-qy);
				int kCIx = ((kCx+Mid+2*End)%End);// labels of the image in C
				int kCIy = ((kCy+Mid+2*End)%End);
				kCx  = ((kCIx+End-1)%End)+1-Mid; // correct
				kCy  = ((kCIy+End-1)%End)+1-Mid ; // correct

//				float C2   = ::sqrt((float)(kCx*kCx+ kCy*kCy));
				int CountCxy  = (kCx+Mid-1)*End+(kCy+Mid-1);
				float fCR     = fkVecR[CountCxy];
				float fCI     = fkVecI[CountCxy];


				float bispectemp  = (fkR*(fqR*fCR -fqI*fCI) -fkI*(fqI*fCR  +fqR*fCI));


				int k2IndLo  = 0; while ((k2>=radRange[k2IndLo+1]) && (k2IndLo+1 < LradRange ) ) k2IndLo +=1;
				int k2IndHi = k2IndLo;
				float k2Lo= radRange[k2IndLo];
				if (k2IndLo+1< LradRange) {
					k2IndHi   = k2IndLo+1;
				}
//				float k2Hi= radRange[k2IndHi];

				float kCof =k2-k2Lo;

				int q2IndLo  = 0; while ((q2>=radRange[q2IndLo+1]) && (q2IndLo+1 < LradRange ) ) q2IndLo +=1;
				int q2IndHi=q2IndLo;
				float q2Lo= radRange[q2IndLo];
				if (q2IndLo+1 < LradRange)  {
					q2IndHi   = q2IndLo+1 ;
				}
				float qCof = q2-q2Lo;


				for (int jk =1; jk<=2; jk++){
				for (int jq =1; jq<=2; jq++){

					float Weight = (kCof+(1-2*kCof)*(jk==1))*(qCof+(1-2*qCof)*(jq==1));

					int k2Ind      =  k2IndLo*(jk==1)      +   k2IndHi*(jk==2);
					int q2Ind      =  q2IndLo*(jq==1)      +   q2IndHi*(jq==2);
					int TotalInd   = q2Ind*LradRange+k2Ind;
					RotTransInv[TotalInd] += Weight*bispectemp;
					WeightInv[TotalInd]   +=  Weight;
	//				cout << "k2Ind= "<< k2Ind  << " q2Ind="<< q2Ind  << "Weight=" << Weight << endl;
				}}
			} // Countqxy
		} // Countkxy

//		cout << "Finished Main Section " << endl;
//		cout << " LradRange " <<LradRange <<  endl;


		EMData *RotTransInvF  = new  EMData(LradRange,LradRange);
		EMData *WeightImage   = new  EMData(LradRange,LradRange);

		for (int jk =0; jk<LradRange ; jk++){// LradRange
		for (int jq =jk; jq<LradRange; jq++){ // LradRange
			int TotalInd      = jq*LradRange+jk;
			int TotalIndBar   = jq*LradRange+jk;
			float Weight = WeightInv[TotalInd] + WeightInv[TotalIndBar];
			if (Weight <=0) continue;
			WeightImage    -> set_value_at(jk,jq,Weight);
			WeightImage    -> set_value_at(jq,jk,Weight);
#ifdef _WIN32
			float ValNow  = pow( (RotTransInv[TotalInd] + RotTransInv[TotalIndBar]) / Weight, 1.0f/3.0f )  ;
#else
			float ValNow  = cbrt( (RotTransInv[TotalInd] + RotTransInv[TotalIndBar]) / Weight )  ;
#endif	//_WIN32
			RotTransInvF -> set_value_at(jk,jq,ValNow);//  include /Weight
 			RotTransInvF -> set_value_at(jq,jk,ValNow );//  include /Weight
		}}

#ifdef	_WIN32
		_unlink("WeightImage.???");
#else
		int rslt;
		rslt = system("rm -f WeightImage.???"); rslt++;
#endif	//_WIN32
		WeightImage  -> write_image("WeightImage.img");

		return  RotTransInvF ;
	}
	return 0;
}


void EMData::insert_clip(const EMData * const block, const IntPoint &origin) {
	int nx1 = block->get_xsize();
	int ny1 = block->get_ysize();
	int nz1 = block->get_zsize();

	Region area(origin[0], origin[1], origin[2],nx1, ny1, nz1);

	//make sure the block fits in EMData 
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

	if (x1 < x0 || y1 < y0 || z1 < z0) return; // out of bounds, this is fine, nothing happens

	size_t clipped_row_size = (x1-x0) * sizeof(float);
	size_t src_secsize =  (size_t)(nx1 * ny1);
	size_t dst_secsize = (size_t)(nx * ny);

/*
#ifdef EMAN2_USING_CUDA
	if(block->cudarwdata){
		// this is VERY slow.....
		float *cudasrc = block->cudarwdata + zd0 * src_secsize + yd0 * nx1 + xd0;
		if(!cudarwdata) rw_alloc();
		float *cudadst = cudarwdata + z0 * dst_secsize + y0 * nx + x0;
		for (int i = z0; i < z1; i++) {
			for (int j = y0; j < y1; j++) {
				//printf("%x %x %d\n", cudadst, cudasrc, clipped_row_size);
				cudaMemcpy(cudadst,cudasrc,clipped_row_size,cudaMemcpyDeviceToDevice);
				cudasrc += nx1;
				cudadst += nx;
			}
			cudasrc += src_gap;
			cudadst += dst_gap;
		}
		return;
	}
#endif
*/
	float *src = block->get_data() + zd0 * src_secsize + yd0 * nx1 + xd0;
	float *dst = get_data() + z0 * dst_secsize + y0 * nx + x0;
	
	size_t src_gap = src_secsize - (y1-y0) * nx1;
	size_t dst_gap = dst_secsize - (y1-y0) * nx;
	
	for (int i = z0; i < z1; i++) {
		for (int j = y0; j < y1; j++) {
			EMUtil::em_memcpy(dst, src, clipped_row_size);
			src += nx1;
			dst += nx;
		}
		src += src_gap;
		dst += dst_gap;
	}
	
#ifdef EMAN2_USING_CUDA	
	if(block->cudarwdata){
		copy_to_cuda(); // copy back to device as padding is faster on the host
	}
#endif

	update();
	EXITFUNC;
}


void EMData::insert_scaled_sum(EMData *block, const FloatPoint &center,
						   float scale, float)
{
	ENTERFUNC;
	float * data = get_data();
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
		if (x1>=get_xsize()) x1=get_xsize()-1;
		if (y1>=get_ysize()) y1=get_ysize()-1;
		if (z1>=get_zsize()) z1=get_zsize()-1;

		float bx=block->get_xsize()/2.0f;
		float by=block->get_ysize()/2.0f;
		float bz=block->get_zsize()/2.0f;

		size_t idx;
		for (int x=x0; x<=x1; x++) {
			for (int y=y0; y<=y1; y++) {
				for (int z=z0; z<=z1; z++) {
					idx = x + y * nx + (size_t)z * nx * ny;
					data[idx] +=
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
		if (x1>=get_xsize()) x1=get_xsize()-1;
		if (y1>=get_ysize()) y1=get_ysize()-1;

		float bx=block->get_xsize()/2.0f;
		float by=block->get_ysize()/2.0f;

		for (int x=x0; x<=x1; x++) {
			for (int y=y0; y<=y1; y++) {
				data[x + y * nx] +=
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
// 			else if ( m == 0 )
// 			{
// 				if ( n_f == -ny/2 )
// 				{
// 					t2++;
// // 					continue;
// 					for (int y = 0; y < return_slice->get_ysize(); ++y) {
// 						for (int x = 0; x < return_slice->get_xsize(); ++x) {
// 							double cur_val = return_slice->get_value_at(x,y);
// 							return_slice->set_value_at(x,y,cur_val+dat[idx]*std::pow(-1.0f,y));
// 						}
// 					}
// 					if (phase > 0.01 ) cout << "foo 2 " << phase << " " << amp << " " << dat[idx] << endl;
// 				}
// 				else
// 				{
// 					if ( n_f < 1 ) continue;
// 					t3++;
// 					for (int y = 0; y < return_slice->get_ysize(); ++y) {
// 						for (int x = 0; x < return_slice->get_xsize(); ++x) {
// 							double cur_val = return_slice->get_value_at(x,y);
// 							return_slice->set_value_at(x,y,cur_val+2*amp*cos(ndash*y+phase));
// 						}
// 					}
// 				}
// 			}
// 			else if ( n_f == -ny/2 )
// 			{
// 				if ( m == ((nx-2)/2) )
// 				{
// 					t4++;
// 					for (int y = 0; y < return_slice->get_ysize(); ++y) {
// 						for (int x = 0; x < return_slice->get_xsize(); ++x) {
// 							double cur_val = return_slice->get_value_at(x,y);
// 							return_slice->set_value_at(x,y,cur_val+dat[idx]*std::pow(-1.0f,x+y));
// 						}
// 					}
// 					if (phase > 0.01 ) cout << "foo 4 " << phase << " " << amp << " " << dat[idx] << endl;
// 				}
// 				else
// 				{
// 					t5++;
// 					for (int y = 0; y < return_slice->get_ysize(); ++y) {
// 						for (int x = 0; x < return_slice->get_xsize(); ++x) {
// 							double cur_val = return_slice->get_value_at(x,y);
// 							return_slice->set_value_at(x,y,cur_val+2*amp*cos(mdash*x+phase));
// 						}
// 					}
// 				}
// 			}
// 			else if ( n_f == 0 )
// 			{
// 				if ( m == ((nx-2)/2) )
// 				{
// 					t6++;
// 					for (int y = 0; y < return_slice->get_ysize(); ++y) {
// 						for (int x = 0; x < return_slice->get_xsize(); ++x) {
// 							double cur_val = return_slice->get_value_at(x,y);
// 							return_slice->set_value_at(x,y,cur_val+dat[idx]*std::pow(-1.0f,x));
// 						}
// 					}
// 					if (phase > 0.01 ) cout << "foo 3 " << phase << " " << amp << " " << dat[idx] << endl;
// 				}
// 				else
// 				{
// 					t7++;
// 					for (int y = 0; y < return_slice->get_ysize(); ++y) {
// 						for (int x = 0; x < return_slice->get_xsize(); ++x) {
// 							double cur_val = return_slice->get_value_at(x,y);
// 							return_slice->set_value_at(x,y,cur_val+2*amp*cos(mdash*x+M_PI*y + phase));
// 						}
// 					}
// 				}
// 			}
// 			else if ( m == ((nx-2)/2) )
// 			{
// 				if ( n_f < 1 ) continue;
// 				t8++;
// 				for (int y = 0; y < return_slice->get_ysize(); ++y) {
// 					for (int x = 0; x < return_slice->get_xsize(); ++x) {
// 						double cur_val = return_slice->get_value_at(x,y);
// 						return_slice->set_value_at(x,y,cur_val+2*amp*cos(ndash*y+M_PI*x+phase));
// 					}
// 				}
// 			}
// }
