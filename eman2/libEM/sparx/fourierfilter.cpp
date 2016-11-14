/*
 * Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
 * Copyright (c) 2000-2006 The University of Texas - Houston Medical School
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
 */

#include "emdata.h"
#include "processor.h"
#include <algorithm>
#include <cstdlib>

using namespace EMAN;
using namespace std;

/*
 ******************************************************
 *DISCLAIMER
 * 03/29/05 P.A.Penczek
 * The University of Texas
 * Pawel.A.Penczek@uth.tmc.edu
 * Please do not modify the content of this document
 * without a written permission of the author.
 ******************************************************/
/*
   Fourier filter
Purpose: Apply selected Fourier space filters to 1-2-3D image.
Method: Calculate FFT (if needed), apply the filter, if the input was real, do the iFFT, otherwise exit. .
Real input -> padded workspace -> Reallocate workspace to real output.
Complex input -> complex output.
Input: f real or complex 1-2-3D image
Output: 1-2-3D filtered image (real or complex).
 */
EMData* Processor::EMFourierFilterFunc(EMData * fimage, Dict params, bool doInPlace)
{
	int    nx, ny, nz, nyp2, nzp2, ix, iy, iz, jx, jy, jz;
	float  dx, dy, dz, omega=0, omegaL=0, omegaH=0;
	float  center=0, gamma=0, argx, argy, argz;
	float  aa, eps, ord=0, cnst=0, aL, aH, cnstL=0, cnstH=0;
	bool   complex_input;
	vector<float> table;
	int undoctf=0;
	float voltage=100.0f, ak=0.0f, cs=2.0f, ps=1.0f, b_factor=0.0f, wgh=0.1f, sign=-1.0f, dza = 0.0f, azz = 0.0f;
	float  tf=0.0f;
	if (!fimage)  return NULL;
	const int ndim = fimage->get_ndim();
	// Set amount of Fourier padding
	// dopad should be a bool, but EMObject Dict's can't store bools.
	int dopad = params["dopad"];
	int npad;
	if (0 == dopad) {
		// no padding
		npad = 1;
	} else if (1 == dopad) {
		// 2x padding (hard-wired)
		npad = 2;
	} else if (2 == dopad) {
        npad = 4;
	} else {
		// invalid entry
		LOGERR("The dopad parameter must be 0 (false) or 1 (true)");
		return NULL; // FIXME: replace w/ exception throw
	}

	// If the input image is already a Fourier image, then we want to
	// have this routine return a Fourier image
	complex_input = fimage->is_complex();
	if ( complex_input && 1 == dopad ) {
		// Cannot pad Fourier input image
		LOGERR("Cannot pad Fourier input image");
		return NULL; // FIXME: replace w/ exception throw
	}

	Util::KaiserBessel* kbptr = 0;


	nx  = fimage->get_xsize();
	ny  = fimage->get_ysize();
	nz  = fimage->get_zsize();
		// We manifestly assume no zero-padding here, just the
		// necessary extension along x for the fft
	if (complex_input) nx = (nx - 2 + fimage->is_fftodd());

	const int nxp = npad*nx;
	const int nyp = (ny > 1) ? npad*ny : 1;
	const int nzp = (nz > 1) ? npad*nz : 1;

	int lsd2 = (nxp + 2 - nxp%2) / 2; // Extended x-dimension of the complex image
	int lsd3 = lsd2 - 1;

	//  Perform padding (if necessary) and fft, if the image is not already an fft image
	EMData* fp = NULL; // workspace image
	if (complex_input) {
		if (doInPlace) {
			// it's okay to change the original image
			fp = fimage;
		} else {
			// fimage must remain pristine
			fp = fimage->copy();
		}
	} else {
		if (doInPlace) {
			if (npad>1) {
				LOGERR("Cannot pad with inplace filter");
				return NULL;	// FIXME, exception
			}
			fp=fimage;
			fp->do_fft_inplace();
		} else {
			fp = fimage->norm_pad( false, npad, 1);
			fp->do_fft_inplace();
		}
	}
	fp->set_array_offsets(1,1,1);

	//  And the filter type is:
	int filter_type = params["filter_type"];

	nyp2 = nyp/2; nzp2 = nzp/2;
	dx = 1.0f/float(nxp);
#ifdef _WIN32
	dy = 1.0f/_cpp_max(float(nyp),1.0f);
	dz = 1.0f/_cpp_max(float(nzp),1.0f);
#else
	dy = 1.0f/std::max(float(nyp),1.0f);
	dz = 1.0f/std::max(float(nzp),1.0f);
#endif	//_WIN32
	float dx2 = dx*dx, dy2 = dy*dy, dz2 = dz*dz;

	vector<float>::size_type tsize;
	float sz[3];
	float szmax;
	vector<float>::size_type maxsize;
	float xshift=0.0, yshift=0.0, zshift=0.0;

	// For the given type of filter set up any necessary parameters for the
	// filter calculation.  FIXME: Need parameter bounds checking!
	switch (filter_type) {
		case TOP_HAT_LOW_PASS:
		case TOP_HAT_HIGH_PASS:
			omega = params["cutoff_abs"];
			omega = 1.0f/omega/omega;
			break;
		case TOP_HAT_BAND_PASS:
			omegaL = params["low_cutoff_frequency"];
			omegaH = params["high_cutoff_frequency"];
			omegaL = 1.0f/omegaL/omegaL;
			omegaH = 1.0f/omegaH/omegaH;
			break;
		case TOP_HOMOMORPHIC:
			omegaL = params["low_cutoff_frequency"];
			omegaH = params["high_cutoff_frequency"];
			gamma  = params["value_at_zero_frequency"];
			omegaL = 1.0f/omegaL/omegaL;
			omegaH = 1.0f/omegaH/omegaH;
			break;
		case GAUSS_LOW_PASS:
		case GAUSS_HIGH_PASS:
		case GAUSS_INVERSE:
			omega = params["cutoff_abs"];
			omega = 0.5f/omega/omega;
			break;
		case GAUSS_BAND_PASS:
			omega = params["cutoff_abs"];
			center = params["center"];
			omega = 0.5f/omega/omega;
			break;
		case GAUSS_HOMOMORPHIC:
			omega = params["cutoff_abs"];
			gamma = params["value_at_zero_frequency"];
			omega = 0.5f/omega/omega;
			gamma = 1.0f-gamma;
			break;
		case BUTTERWORTH_LOW_PASS:
		case BUTTERWORTH_HIGH_PASS:
			omegaL = params["low_cutoff_frequency"];
			omegaH = params["high_cutoff_frequency"];
			eps = 0.882f;
			aa = 10.624f;
			ord = 2.0f*log10(eps/sqrt(aa*aa-1.0f))/log10(omegaL/omegaH);
			omegaL = omegaL/pow(eps,2.0f/ord);
			break;
		case BUTTERWORTH_HOMOMORPHIC:
			omegaL = params["low_cutoff_frequency"];
			omegaH = params["high_cutoff_frequency"];
			gamma  = params["value_at_zero_frequency"];
			eps = 0.882f;
			aa  = 10.624f;
			ord = 2.0f*log10(eps/sqrt(pow(aa,2)-1.0f))/log10(omegaL/omegaH);
			omegaL = omegaL/pow(eps,2.0f/ord);
			gamma = 1.0f-gamma;
			break;
		case SHIFT:
			xshift = params["x_shift"];
			yshift = params["y_shift"];
			zshift = params["z_shift"];
			//origin_type = params["origin_type"];
			break;
		case TANH_LOW_PASS:
		case TANH_HIGH_PASS:
			omega = params["cutoff_abs"];
			aa = params["fall_off"];
			cnst = float(pihalf/aa/omega);
			break;
		case TANH_HOMOMORPHIC:
			omega = params["cutoff_abs"];
			aa = params["fall_off"];
			gamma = params["value_at_zero_frequency"];
			cnst = float(pihalf/aa/omega);
			gamma=1.0f-gamma;
			break;
		case TANH_BAND_PASS:
			omegaL = params["low_cutoff_frequency"];
			aL = params["Low_fall_off"];
			omegaH = params["high_cutoff_frequency"];
			aH = params["high_fall_off"];
			cnstL = float(pihalf/aL/(omegaH-omegaL));
			cnstH = float(pihalf/aH/(omegaH-omegaL));
			break;
		case CTF_:
			dz       = params["defocus"];
			cs       = params["Cs"];
			voltage  = params["voltage"];
			ps       = params["Pixel_size"];
			b_factor = params["B_factor"];
			wgh      = params["amp_contrast"];
			sign     = params["sign"];
			undoctf  = params["undo"];
			ix       = params["binary"];
			if(ix == 1) {undoctf = 2;  b_factor=0.0;} //ignore B-factor for the binary CTF
			dza = params["dza"];
			azz = params["azz"];
			break;
		case KAISER_I0:
		case KAISER_SINH:
		case KAISER_I0_INVERSE:
		case KAISER_SINH_INVERSE:
			{
				float alpha = params["alpha"];
				int       K = params["K"];
				float     r = params["r"];
				float     v = params["v"];
				int       N = params["N"];
			        kbptr = new Util::KaiserBessel(alpha, K, r, v, N);
				break;
			}//without this bracket, compiler on water will complain about crosses initialization
		case RADIAL_TABLE:
			table = params["table"];
			tsize = table.size();
			sz[0] = static_cast<float>(lsd2);
			sz[1] = static_cast<float>(nyp2);
			sz[2] = static_cast<float>(nzp2);
			szmax = *max_element(&sz[0],&sz[3]);
			// for 2d, sqrt(2) = 1.414, relax a little bit to 1.6
			// for 3d, sqrt(3) = 1.732, relax a little bit to 1.9
			if (nzp > 1) {maxsize = vector<float>::size_type(1.9*szmax);} else {maxsize = vector<float>::size_type(1.6*szmax);}
			for (vector<float>::size_type i = tsize+1; i < maxsize; i++) table.push_back(0.f);
			break;
		default:
			LOGERR("Unknown Fourier Filter type");
			return NULL; // FIXME: replace w/ exception throw
	}
	// Perform actual calculation
	//  Gaussian bandpass is the only one with center for frequencies
	switch (filter_type) {
		case GAUSS_BAND_PASS:
			for ( iz = 1; iz <= nzp; iz++) {
				jz=iz-1; if (jz>nzp2) jz=jz-nzp; argz = float(jz*jz)*dz2;
				for ( iy = 1; iy <= nyp; iy++) {
					jy=iy-1; if (jy>nyp2) jy=jy-nyp; argy = argz + float(jy*jy)*dy2;
					for ( ix = 1; ix <= lsd2; ix++) {
						jx=ix-1; argx = argy + float(jx*jx)*dx2;
						fp->cmplx(ix,iy,iz) *= exp(-pow(sqrt(argx)-center,2)*omega);
					}
				}
			}
			break;
		case TOP_HAT_LOW_PASS:
			for ( iz = 1; iz <= nzp; iz++) {
				jz=iz-1; if (jz>nzp2) jz=jz-nzp; argz = float(jz*jz)*dz2;
				for ( iy = 1; iy <= nyp; iy++) {
					jy=iy-1; if (jy>nyp2) jy=jy-nyp; argy = argz + float(jy*jy)*dy2;
					for ( ix = 1; ix <= lsd2; ix++) {
						jx=ix-1; argx = argy + float(jx*jx)*dx2;
						if (argx*omega>1.0f) fp->cmplx(ix,iy,iz) = 0.0f;
					}
				}
			}
			break;
		case TOP_HAT_HIGH_PASS:
			for ( iz = 1; iz <= nzp; iz++) {
				jz=iz-1; if (jz>nzp2) jz=jz-nzp; argz = float(jz*jz)*dz2;
				for ( iy = 1; iy <= nyp; iy++) {
					jy=iy-1; if (jy>nyp2) jy=jy-nyp; argy = argz + float(jy*jy)*dy2;
					for ( ix = 1; ix <= lsd2; ix++) {
						jx=ix-1; argx = argy + float(jx*jx)*dx2;
						if (argx*omega<=1.0f) fp->cmplx(ix,iy,iz) = 0.0f;
					}
				}
			}				break;
		case TOP_HAT_BAND_PASS:
			for ( iz = 1; iz <= nzp; iz++) {
				jz=iz-1; if (jz>nzp2) jz=jz-nzp; argz = float(jz*jz)*dz2;
				for ( iy = 1; iy <= nyp; iy++) {
					jy=iy-1; if (jy>nyp2) jy=jy-nyp; argy = argz + float(jy*jy)*dy2;
					for ( ix = 1; ix <= lsd2; ix++) {
						jx=ix-1; argx = argy + float(jx*jx)*dx2;
						if (argx*omegaL<1.0f || argx*omegaH>=1.0f) fp->cmplx(ix,iy,iz) = 0.0f;
					}
				}
			}
			break;
		case TOP_HOMOMORPHIC:
			for ( iz = 1; iz <= nzp; iz++) {
				jz=iz-1; if (jz>nzp2) jz=jz-nzp; argz = float(jz*jz)*dz2;
				for ( iy = 1; iy <= nyp; iy++) {
					jy=iy-1; if (jy>nyp2) jy=jy-nyp; argy = argz + float(jy*jy)*dy2;
					for ( ix = 1; ix <= lsd2; ix++) {
						jx=ix-1; argx = argy + float(jx*jx)*dx2;
						if (argx*omegaH>1.0f)      fp->cmplx(ix,iy,iz)  = 0.0f;
						else if (argx*omegaL<=1.0f) fp->cmplx(ix,iy,iz) *= gamma;
					}
				}
			}
			break;
		case GAUSS_LOW_PASS :
			for ( iz = 1; iz <= nzp; iz++) {
				jz=iz-1; if (jz>nzp2) jz=jz-nzp; argz = float(jz*jz)*dz2;
				for ( iy = 1; iy <= nyp; iy++) {
					jy=iy-1; if (jy>nyp2) jy=jy-nyp; argy = argz + float(jy*jy)*dy2;
					for ( ix = 1; ix <= lsd2; ix++) {
						jx=ix-1; argx = argy + float(jx*jx)*dx2;
						fp->cmplx(ix,iy,iz) *= exp(-argx*omega);
					}
				}
			}
			break;
		case GAUSS_HIGH_PASS:
			for ( iz = 1; iz <= nzp; iz++) {
				jz=iz-1; if (jz>nzp2) jz=jz-nzp; argz = float(jz*jz)*dz2;
				for ( iy = 1; iy <= nyp; iy++) {
					jy=iy-1; if (jy>nyp2) jy=jy-nyp; argy = argz + float(jy*jy)*dy2;
					for ( ix = 1; ix <= lsd2; ix++) {
						jx=ix-1; argx = argy + float(jx*jx)*dx2;
						fp->cmplx(ix,iy,iz) *= 1.0f-exp(-argx*omega);
					}
				}
			}
			break;
		case GAUSS_HOMOMORPHIC:
			for ( iz = 1; iz <= nzp; iz++) {
				jz=iz-1; if (jz>nzp2) jz=jz-nzp; argz = float(jz*jz)*dz2;
				for ( iy = 1; iy <= nyp; iy++) {
					jy=iy-1; if (jy>nyp2) jy=jy-nyp; argy = argz + float(jy*jy)*dy2;
					for ( ix = 1; ix <= lsd2; ix++) {
						jx=ix-1; argx = argy + float(jx*jx)*dx2;
						fp->cmplx(ix,iy,iz) *= 1.0f-gamma*exp(-argx*omega);
					}
				}
			}
			break;
		case GAUSS_INVERSE :
			for ( iz = 1; iz <= nzp; iz++) {
				jz=iz-1; if (jz>nzp2) jz=jz-nzp; argz = float(jz*jz)*dz2;
				for ( iy = 1; iy <= nyp; iy++) {
					jy=iy-1; if (jy>nyp2) jy=jy-nyp; argy = argz + float(jy*jy)*dy2;
					for ( ix = 1; ix <= lsd2; ix++) {
						jx=ix-1; argx = argy + float(jx*jx)*dx2;
						fp->cmplx(ix,iy,iz) *= exp(argx*omega);
					}
				}
			}
			break;
		case KAISER_I0:   // K-B filter
			for ( iz = 1; iz <= nzp; iz++) {
				jz=iz-1; if (jz>nzp2) jz=jz-nzp;
				float nuz = jz*dz;
				for ( iy = 1; iy <= nyp; iy++) {
					jy=iy-1; if (jy>nyp2) jy=jy-nyp;
					float nuy = jy*dy;
					for ( ix = 1; ix <= lsd2; ix++) {
						jx=ix-1;
						float nux = jx*dx;
						//if (!kbptr)
						//	throw
						//		NullPointerException("kbptr null!");
						switch (ndim) {
							case 3:
								fp->cmplx(ix,iy,iz) *= kbptr->i0win(nux)*kbptr->i0win(nuy)*kbptr->i0win(nuz);
								break;
							case 2:
								fp->cmplx(ix,iy,iz) *= kbptr->i0win(nux)*kbptr->i0win(nuy);
								break;
							case 1:
								fp->cmplx(ix,iy,iz)*= kbptr->i0win(nux);
								break;
						}
					}
				}
			}
			break;
		case KAISER_SINH:   //  Sinh filter
			for ( iz = 1; iz <= nzp; iz++) {
				jz=iz-1; if (jz>nzp2) jz=jz-nzp;
				for ( iy = 1; iy <= nyp; iy++) {
					jy=iy-1; if(jy>nyp2) jy=jy-nyp;
					for ( ix = 1; ix <= lsd2; ix++) {
						jx=ix-1;
						//if (!kbptr)
						//	throw
						//		NullPointerException("kbptr null!");
						switch (ndim) {
							case 3:
								fp->cmplx(ix,iy,iz)*= kbptr->sinhwin((float)jx)*kbptr->sinhwin((float)jy)*kbptr->sinhwin((float)jz);
								break;
							case 2:
								fp->cmplx(ix,iy,iz)*= kbptr->sinhwin((float)jx)*kbptr->sinhwin((float)jy);
								break;
							case 1:
								fp->cmplx(ix,iy,iz)*= kbptr->sinhwin((float)jx);
								//float argu = kbptr->sinhwin((float) jx);
								//cout << jx<<"  "<< nux<<"  "<<argu<<endl;
								break;
						}
					}
				}
			}
			break;
		case KAISER_I0_INVERSE:   // 1./(K-B filter)
			for ( iz = 1; iz <= nzp; iz++) {
				jz=iz-1; if (jz>nzp2) jz=jz-nzp;
				float nuz = jz*dz;
				for ( iy = 1; iy <= nyp; iy++) {
					jy=iy-1; if(jy>nyp2) jy=jy-nyp;
					float nuy = jy*dy;
					for ( ix = 1; ix <= lsd2; ix++) {
						jx=ix-1;
						float nux = jx*dx;
					//if (!kbptr)
					//	throw
					//		NullPointerException("kbptr null!");
						switch (ndim) {
							case 3:
								fp->cmplx(ix,iy,iz) /= (kbptr->i0win(nux)*kbptr->i0win(nuy)*kbptr->i0win(nuz));
								break;
							case 2:
								fp->cmplx(ix,iy,iz) /= (kbptr->i0win(nux)*kbptr->i0win(nuy));
								break;
							case 1:
								fp->cmplx(ix,iy,iz) /= kbptr->i0win(nux);
								break;
						}
					}
				}
			}
			break;
		case KAISER_SINH_INVERSE:  // 1./sinh
			for ( iz = 1; iz <= nzp; iz++) {
				jz=iz-1; if (jz>nzp2) jz=jz-nzp;
				for ( iy = 1; iy <= nyp; iy++) {
					jy=iy-1; if (jy>nyp2) jy=jy-nyp;
					for ( ix = 1; ix <= lsd2; ix++) {
						jx=ix-1;
						//if (!kbptr)
						//	throw
						//		NullPointerException("kbptr null!");
						switch (ndim) {
							case 3:
								fp->cmplx(ix,iy,iz) /= (kbptr->sinhwin((float)jx)*kbptr->sinhwin((float)jy)*kbptr->sinhwin((float)jz));
								break;
							case 2:
								fp->cmplx(ix,iy,iz) /= (kbptr->sinhwin((float)jx)*kbptr->sinhwin((float)jy));
								break;
							case 1:
								fp->cmplx(ix,iy,iz) /= kbptr->sinhwin((float)jx);
								//float argu = kbptr->sinhwin((float) jx);
								//cout << jx<<"  "<< nux<<"  "<<argu<<endl;
								break;
						}
					}
				}
			}
			break;
		case BUTTERWORTH_LOW_PASS:
			for ( iz = 1; iz <= nzp; iz++) {
				jz=iz-1; if (jz>nzp2) jz=jz-nzp; argz = float(jz*jz)*dz2;
				for ( iy = 1; iy <= nyp; iy++) {
					jy=iy-1; if (jy>nyp2) jy=jy-nyp; argy = argz + float(jy*jy)*dy2;
					for ( ix = 1; ix <= lsd2; ix++) {
						jx=ix-1; argx = argy + float(jx*jx)*dx2;
						fp->cmplx(ix,iy,iz) *= sqrt(1.0f/(1.0f+pow(sqrt(argx)/omegaL,ord)));
					}
				}
			}
			break;
		case BUTTERWORTH_HIGH_PASS:
			for ( iz = 1; iz <= nzp; iz++) {
				jz=iz-1; if (jz>nzp2) jz=jz-nzp; argz = float(jz*jz)*dz2;
				for ( iy = 1; iy <= nyp; iy++) {
					jy=iy-1; if (jy>nyp2) jy=jy-nyp; argy = argz + float(jy*jy)*dy2;
					for ( ix = 1; ix <= lsd2; ix++) {
						jx=ix-1; argx = argy + float(jx*jx)*dx2;
						fp->cmplx(ix,iy,iz) *= 	1.0f-sqrt(1.0f/(1.0f+pow(sqrt(argx)/omegaL,ord)));
					}
				}
			}
			break;
		case BUTTERWORTH_HOMOMORPHIC:
			for ( iz = 1; iz <= nzp; iz++) {
				jz=iz-1; if (jz>nzp2) jz=jz-nzp; argz = float(jz*jz)*dz2;
				for ( iy = 1; iy <= nyp; iy++) {
					jy=iy-1; if (jy>nyp2) jy=jy-nyp; argy = argz + float(jy*jy)*dy2;
					for ( ix = 1; ix <= lsd2; ix++) {
						jx=ix-1; argx = argy + float(jx*jx)*dx2;
						fp->cmplx(ix,iy,iz) *= 	1.0f-gamma*sqrt(1.0f/(1.0f+pow(sqrt(argx)/omegaL,ord)));
					}
				}
			}
			break;
		case SHIFT:
			//if (origin_type) {
				for ( iz = 1; iz <= nzp; iz++) {
					jz=iz-1; if (jz>nzp2) jz=jz-nzp;
					for ( iy = 1; iy <= nyp; iy++) {
						jy=iy-1; if (jy>nyp2) jy=jy-nyp;
						for ( ix = 1; ix <= lsd2; ix++) {
							jx=ix-1;
							fp->cmplx(ix,iy,iz) *= 	exp(-float(twopi)*iimag*(xshift*jx/nx + yshift*jy/ny+ zshift*jz/nz));
						}
					}
				}
			/*} else {
				for ( iz = 1; iz <= nzp; iz++) {
					jz=iz-1; if (jz>nzp2) jz=jz-nzp;
					if  (iz>nzp2) { kz=iz-nzp2; } else { kz=iz+nzp2; }
					for ( iy = 1; iy <= nyp; iy++) {
						jy=iy-1; if (jy>nyp2) jy=jy-nyp;
						if  (iy>nyp2) { ky=iy-nyp2; } else { ky=iy+nyp2; }
						for ( ix = 1; ix <= lsd2; ix++) {
							jx=ix-1;
							fp->cmplx(ix,ky,kz) *= 	exp(-float(twopi)*iimag*(xshift*jx/nx + yshift*jy/ny+ zshift*jz/nz));
						}
					}
				}
			}*/
			break;
		case TANH_LOW_PASS:
			for ( iz = 1; iz <= nzp; iz++) {
				jz=iz-1; if (jz>nzp2) jz=jz-nzp; argz = float(jz*jz)*dz2;
				for ( iy = 1; iy <= nyp; iy++) {
					jy=iy-1; if (jy>nyp2) jy=jy-nyp; argy = argz + float(jy*jy)*dy2;
					for ( ix = 1; ix <= lsd2; ix++) {
						jx=ix-1; argx = sqrt(argy + float(jx*jx)*dx2);
						fp->cmplx(ix,iy,iz) *= 	0.5f*(tanh(cnst*(argx+omega))-tanh(cnst*(argx-omega)));
					}
				}
			}
			break;
		case TANH_HIGH_PASS:
			for ( iz = 1; iz <= nzp; iz++) {
				jz=iz-1; if (jz>nzp2) jz=jz-nzp; argz = float(jz*jz)*dz2;
				for ( iy = 1; iy <= nyp; iy++) {
					jy=iy-1; if (jy>nyp2) jy=jy-nyp; argy = argz + float(jy*jy)*dy2;
					for ( ix = 1; ix <= lsd2; ix++) {
						jx=ix-1; sqrt(argx = argy + float(jx*jx)*dx2);
						fp->cmplx(ix,iy,iz) *= 	1.0f-0.5f*(tanh(cnst*(argx+omega))-tanh(cnst*(argx-omega)));
					}
				}
			}
			break;
		case TANH_HOMOMORPHIC:
			for ( iz = 1; iz <= nzp; iz++) {
				jz=iz-1; if (jz>nzp2) jz=jz-nzp; argz = float(jz*jz)*dz2;
				for ( iy = 1; iy <= nyp; iy++) {
					jy=iy-1; if (jy>nyp2) jy=jy-nyp; argy = argz + float(jy*jy)*dy2;
					for ( ix = 1; ix <= lsd2; ix++) {
						jx=ix-1; argx = sqrt(argy + float(jx*jx)*dx2);
						fp->cmplx(ix,iy,iz) *= 1.0f-gamma*0.5f*(tanh(cnst*(argx+omega))-tanh(cnst*(argx-omega)));
					}
				}
			}
			break;
		case TANH_BAND_PASS:
			for ( iz = 1; iz <= nzp; iz++) {
				jz=iz-1; if (jz>nzp2) jz=jz-nzp; argz = float(jz*jz)*dz2;
				for ( iy = 1; iy <= nyp; iy++) {
					jy=iy-1; if (jy>nyp2) jy=jy-nyp; argy = argz + float(jy*jy)*dy2;
					for ( ix = 1; ix <= lsd2; ix++) {
						jx=ix-1; argx = sqrt(argy + float(jx*jx)*dx2);
						fp->cmplx(ix,iy,iz) *= 0.5f*(tanh(cnstH*(argx+omegaH))-tanh(cnstH*(argx-omegaH))-tanh(cnstL*(argx+omegaL))+tanh(cnstL*(argx-omegaL)));
					}
				}
			}
			break;
		case RADIAL_TABLE:
			for ( iz = 1; iz <= nzp; iz++) {
				jz=iz-1; if (jz>nzp2) jz=jz-nzp; argz = float(jz*jz)*dz2;
				for ( iy = 1; iy <= nyp; iy++) {
					jy=iy-1; if (jy>nyp2) jy=jy-nyp; argy = argz + float(jy*jy)*dy2;
					for ( ix = 1; ix <= lsd2; ix++) {
						jx=ix-1; argx = argy + float(jx*jx)*dx2;
						float rf = sqrt( argx )*nxp;
						int  ir = int(rf);
						float df = rf - float(ir);
						float f = table[ir] + df * (table[ir+1] - table[ir]); // (1-df)*table[ir]+df*table[ir+1];
						fp->cmplx(ix,iy,iz) *= f;
					}
				}
			}
			break;
		case CTF_:
			for ( iz = 1; iz <= nzp; iz++) {
				jz=iz-1; if (jz>nzp2) jz=jz-nzp;
				for ( iy = 1; iy <= nyp; iy++) {
					jy=iy-1; if (jy>nyp2) jy=jy-nyp;
					for ( ix = 1; ix <= lsd2; ix++) {
						jx=ix-1;
						if(ny>1 && nz<=1 ) {
							//  astigmatism makes sense only on 2D
							ak = sqrt(static_cast<float>(jx)/lsd3*static_cast<float>(jx)/lsd3 +
		        					static_cast<float>(jy)/nyp2*static_cast<float>(jy)/nyp2)/ps/2.0f;
							if(dza == 0.0f)  tf = Util::tf(dz, ak, voltage, cs, wgh, b_factor, sign);
							else {
								float az = atan2(static_cast<float>(jy)/nyp2, static_cast<float>(jx)/lsd3);
								float dzz = dz - dza/2.0f*sin(2*(az+azz*M_PI/180.0f));
								tf = Util::tf(dzz, ak, voltage, cs, wgh, b_factor, sign);
							}
						}  else if(ny<=1) {
							ak=sqrt(static_cast<float>(jx)/lsd3*static_cast<float>(jx)/lsd3)/ps/2.0f;
							tf = Util::tf(dz, ak, voltage, cs, wgh, b_factor, sign);
						}  else if(nz>1)  {
							ak=sqrt(static_cast<float>(jx)/lsd3*static_cast<float>(jx)/lsd3 +
								static_cast<float>(jy)/nyp2*static_cast<float>(jy)/nyp2 +
								static_cast<float>(jz)/nzp2*static_cast<float>(jz)/nzp2)/ps/2.0f;
							tf  =Util::tf(dz, ak, voltage, cs, wgh, b_factor, sign);
						}
						switch (undoctf) {
						case 0:
						    fp->cmplx(ix,iy,iz) *= tf;
						    break;
						case 1:
						    if( tf>0 && tf <  1e-5 ) tf =  1e-5f;
						    if( tf<0 && tf > -1e-5 ) tf = -1e-5f;
						    fp->cmplx(ix,iy,iz) /= tf;
						    break;
						case 2:
						    if(tf < 0.0f) fp->cmplx(ix,iy,iz) *= -1.0f;
						    break;
						}
					}
				}
			}
			break;
	}
	delete kbptr; kbptr = 0;
	if (!complex_input) {
		fp->do_ift_inplace();
		fp->depad();
	}

		// Return a complex (Fourier) filtered image
		// Note: fp and fimage are the _same_ object if doInPlace
		// is true, so in that case fimage has been filtered.
		// We always return an image (pointer), but if the function
		// was called with doInPlace == true then the calling function
		// will probably ignore the return value.

		// ELSE Return a real-space filtered image
		//
		// On 12/15/2006 Wei Zhang comment:
		// If input is reald and doInPlace == true, we might need delete fp after copy its
		// data back to fimage, since fp is allocated inside this function and is ignored
		// by caller if doInPlace == true. As a reminder, the caller is EMFourierFuncInPlace
		//
	fp->set_array_offsets(0,0,0);
	fp->update();
	if (doInPlace && !complex_input) {
		// copy workspace data into the original image
		float* orig = fimage->get_data();
		float* work = fp->get_data();
		for (size_t i = 0; i < (size_t)nx*ny*nz; i++) orig[i] = work[i];
		fimage->update();
	}
	EXITFUNC;
	return fp;
}
