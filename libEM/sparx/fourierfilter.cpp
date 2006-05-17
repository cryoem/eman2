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
	float lambda,ak,cs,ps,b_factor,wgh,sign;
	if (!fimage) {
		return NULL;
	}
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

	nx  = fimage->get_xsize();
	ny  = fimage->get_ysize();
	nz  = fimage->get_zsize();
		// We manifestly assume no zero-padding here, just the 
			omega = params["cutoff_frequency"];
		// necessary extension along x for the fft
	if (fimage->is_complex()) nx = (nx - 2 + fimage->is_fftodd()); 

	const int nxp = npad*nx;
	const int nyp = (ny > 1) ? npad*ny : 1;
	const int nzp = (nz > 1) ? npad*nz : 1;

	int lsd2 = (nxp + 2 - nxp%2) / 2; // Extended x-dimension of the complex image

	Util::KaiserBessel* kbptr = 0;

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
		fp = fimage->pad_fft(npad); 
		fp->do_fft_inplace();
	}
	fp->set_array_offsets(1,1,1);

	//  And the filter type is:
	int filter_type = params["filter_type"];

	nyp2 = nyp/2; nzp2 = nzp/2;
	dx = 1.0f/float(nxp); 
#ifdef _WIN32
	dy = 1.0f/_MAX(float(nyp),1.0f);
	dz = 1.0f/_MAX(float(nzp),1.0f);	
#else
	dy = 1.0f/std::max(float(nyp),1.0f);
	dz = 1.0f/std::max(float(nzp),1.0f);
#endif	//_WIN32
	float dx2 = dx*dx, dy2 = dy*dy, dz2 = dz*dz;

	vector<float>::size_type tsize;
	float sz[3];
	float szmax;
	vector<float>::size_type maxsize;
	float xshift, yshift, zshift;

	// For the given type of filter set up any necessary parameters for the
	// filter calculation.  FIXME: Need parameter bounds checking!
	switch (filter_type) {
		case TOP_HAT_LOW_PASS: 
		case TOP_HAT_HIGH_PASS:
			omega = params["cutoff_frequency"];
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
			gamma = params["value_at_zero_frequency"];
			omegaL = 1.0f/omegaL/omegaL; 
			omegaH = 1.0f/omegaH/omegaH;
			break;
		case GAUSS_LOW_PASS: 
		case GAUSS_HIGH_PASS: 
		case GAUSS_INVERSE:
			omega = params["sigma"]; 
			omega = 0.5f/omega/omega;
			break;
		case GAUSS_BAND_PASS: 
			omega = params["sigma"]; 
			center = params["center"];
			omega = 0.5f/omega/omega;
			break;
		case GAUSS_HOMOMORPHIC: 
			omega = params["sigma"]; 
			gamma = params["value_at_zero_frequency"];
			omega = 0.5f/omega/omega;
			gamma=1.0f-gamma;
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
			gamma = params["value_at_zero_frequency"];
			eps = 0.882f;
			aa = 10.624f;
			ord = 2.0f*log10(eps/sqrt(pow(aa,2)-1.0f))/log10(omegaL/omegaH);
			omegaL = omegaL/pow(eps,2.0f/ord);
			gamma=1.0f-gamma;
			break;
		case SHIFT:
			xshift = params["x_shift"];
			yshift = params["y_shift"];
			zshift = params["z_shift"];
			break;
		case TANH_LOW_PASS: 
		case TANH_HIGH_PASS: 
			omega = params["cutoff_frequency"];
			aa = params["fall_off"];
			cnst = float(pihalf/aa/omega);
			break;
		case TANH_HOMOMORPHIC: 
			omega = params["cutoff_frequency"];
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
		     dz=params["defocus"];
		     cs=params["cs"];
		     lambda=params["voltage"];  	       
		     ps=params["ps"];
		     b_factor=params["b_factor"];
		     wgh=params["wgh"];
		     sign=params["sign"];
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
			// for 2d, sqrt(2)/2 ~ 0.75
			maxsize = vector<float>::size_type(0.75*szmax);
			// for 3d, sqrt(3)/2 ~ 0.9
			if (nzp > 1)
				maxsize = 
					vector<float>::size_type(0.9*szmax);
			for (vector<float>::size_type i = tsize+1; i < maxsize; i++) 
				table.push_back(0.f);
			break;
		default: 
			LOGERR("Unknown Fourier Filter type"); 
			return NULL; // FIXME: replace w/ exception throw
	}
	// Perform actual calculation
	//  Gaussian bandpass is the only one with center for frequencies
	if(filter_type == GAUSS_BAND_PASS) {
		for ( iz = 1; iz <= nzp; iz++) {
			jz=iz-1; if(jz>nzp2) jz=jz-nzp; 
			argz = (float(jz)-center)*(float(jz)-center)*dz2;
			for ( iy = 1; iy <= nyp; iy++) {
				jy=iy-1; if(jy>nyp2) jy=jy-nyp; 
				argy = argz + (float(jy)-center)*(float(jy)-center)*dy2;
				for ( ix = 1; ix <= lsd2; ix++) {
					jx=ix-1; argx = argy + (float(jx)-center)*(float(jx)-center)/float((nxp-1)*(nxp-1));
					// RHS of filter calculation for Gaussian bandpass
					fp->cmplx(ix,iy,iz) *= exp(-0.125f*argx*omega);
				}
			}
		}
	} else {
		for ( iz = 1; iz <= nzp; iz++) {
			jz=iz-1; if(jz>nzp2) jz=jz-nzp; argz = float(jz*jz)*dz2;
			float nuz = jz*dz;
			for ( iy = 1; iy <= nyp; iy++) {
				jy=iy-1;
				if(jy>nyp2) jy=jy-nyp; argy = argz + float(jy*jy)*dy2;
				float nuy = jy*dy;
				for ( ix = 1; ix <= lsd2; ix++) {
					jx=ix-1; argx = argy + float(jx*jx)*dx2;
					float nux = jx*dx;
					// RHS of filter calculation 
					switch (filter_type) {
						case TOP_HAT_LOW_PASS:
							if(argx*omega>1.0f) fp->cmplx(ix,iy,iz) = 0; break;
						case TOP_HAT_HIGH_PASS:
							if(argx*omega<=1.0f) fp->cmplx(ix,iy,iz) = 0; break;
						case TOP_HAT_BAND_PASS: 
							if(argx*omegaL>1.0f && argx*omegaH<=1.0f) 
								fp->cmplx(ix,iy,iz) = 0; break;
						case TOP_HOMOMORPHIC: 
							if(argx*omegaH>1.0f) fp->cmplx(ix,iy,iz) = 0; 
							else if(argx*omegaL<=1.0f) fp->cmplx(ix,iy,iz) *= gamma;  
							break;
						case GAUSS_LOW_PASS :
							fp->cmplx(ix,iy,iz) *= exp(-argx*omega); break;
						case GAUSS_HIGH_PASS:
							fp->cmplx(ix,iy,iz) *= 1.0f-exp(-argx*omega); break;
						case GAUSS_HOMOMORPHIC: 
							fp->cmplx(ix,iy,iz) *= 1.0f-gamma*exp(-argx*omega); break;
						case GAUSS_INVERSE :
							fp->cmplx(ix,iy,iz) *= exp(argx*omega); break;
						case KAISER_I0:   // K-B filter
							if (!kbptr) 
								throw 
									NullPointerException("kbptr null!");
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
							break;
						case KAISER_SINH:   //  Sinh filter
							if (!kbptr) 
								throw 
									NullPointerException("kbptr null!");
							switch (ndim) {
								case 3:
									fp->cmplx(ix,iy,iz)*= kbptr->sinhwin(jx)*kbptr->sinhwin(jy)*kbptr->sinhwin(jz);
									break;
								case 2:
									fp->cmplx(ix,iy,iz)*= kbptr->sinhwin(jx)*kbptr->sinhwin(jy);
									break;
								case 1:
									fp->cmplx(ix,iy,iz)*= kbptr->sinhwin(jx);
									//float argu = kbptr->sinhwin((float) jx);
									//cout << jx<<"  "<< nux<<"  "<<argu<<endl;
									break;
							}
							break;
						case KAISER_I0_INVERSE:   // 1./(K-B filter)
							if (!kbptr) 
								throw 
									NullPointerException("kbptr null!");
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
							break;
						case KAISER_SINH_INVERSE:  // 1./sinh
							if (!kbptr) 
								throw 
									NullPointerException("kbptr null!");
							switch (ndim) {
								case 3:
									fp->cmplx(ix,iy,iz) /= (kbptr->sinhwin(jx)*kbptr->sinhwin(jy)*kbptr->sinhwin(jz));
									break;
								case 2:
									fp->cmplx(ix,iy,iz) /= (kbptr->sinhwin(jx)*kbptr->sinhwin(jy));
									break;
								case 1:
									fp->cmplx(ix,iy,iz) /= kbptr->sinhwin(jx);
									//float argu = kbptr->sinhwin((float) jx);
									//cout << jx<<"  "<< nux<<"  "<<argu<<endl;
									break;
							}
							break;
						case BUTTERWORTH_LOW_PASS: 
							fp->cmplx(ix,iy,iz) 
								*= sqrt(1.0f/(1.0f+pow(sqrt(argx)/omegaL,ord))); break;
						case BUTTERWORTH_HIGH_PASS: 
							fp->cmplx(ix,iy,iz) *= 
								1.0f-sqrt(1.0f/(1.0f+pow(sqrt(argx)/omegaL,ord))); break;
						case BUTTERWORTH_HOMOMORPHIC: 
							fp->cmplx(ix,iy,iz) *= 
								1.0f-gamma*sqrt(1.0f/(1.0f+pow(sqrt(argx)/omegaL,ord))); break;
						case SHIFT:
							fp->cmplx(ix,iy,iz) *= 
								exp(-float(twopi)*iimag*(xshift*jx/nx + yshift*jy/ny+ zshift*jz/nz));
							break;
						case TANH_LOW_PASS : 
							argx= sqrt(argx); 
							fp->cmplx(ix,iy,iz) *= 
								0.5f*(tanh(cnst*(argx+omega)))-tanh(cnst*(argx-omega)); break;
						case TANH_HIGH_PASS: 
							argx= sqrt(argx); 
							fp->cmplx(ix,iy,iz) *= 
								1.0f-0.5f*(tanh(cnst*(argx+omega)))-tanh(cnst*(argx-omega)); break;
						case TANH_HOMOMORPHIC: 
							argx= sqrt(argx); 
							fp->cmplx(ix,iy,iz) *= 
								1.0f-gamma*0.5f*(tanh(cnst*(argx+omega)))-tanh(cnst*(argx-omega)); break;
						case TANH_BAND_PASS: 
							argx= sqrt(argx); 
							fp->cmplx(ix,iy,iz) *= 
								0.5f*(tanh(cnstH*(argx+omegaH)))-tanh(cnstH*(argx-omegaH)-tanh(cnstL*(argx+omegaL)))+tanh(cnstL*(argx-omegaL)); break;
						case RADIAL_TABLE:
							float rf = sqrt( static_cast<float>(ix) * static_cast<float>(ix) + 
											 static_cast<float>(iy) * static_cast<float>(iy) + 
											 static_cast<float>(iz) * static_cast<float>(iz) );
							int ir = int(rf);
							float df = rf - float(ir);
							float f = (1-df)*table[ir]+df*table[ir+1];
							fp->cmplx(ix,iy,iz) *= f;
							break;
						case CTF_:
						  if(ny>1 && nz<=1 )						  					
						   { ak=sqrt(jx/lsd2*jx/lsd2+jy/nyp2*jy/nyp2)/ps/2.0f;}
						  else	if(ny<=1)
						  	     { ak=sqrt(jx/lsd2*jx/lsd2)/ps/2.0f;}
						  else  if(nz>1)
						     { ak=sqrt(jx/lsd2*jx/lsd2+jy/nyp2*jy/nyp2+jz/nzp2*jz/nzp2)/ps/2.0f;}          
						       				
						  fp->cmplx(ix,iy,iz) *= Util::tf(dz,ak,12.398f/sqrt(lambda *(1022.f+lambda)),cs*1.0e-7f,atan(wgh/(1.0-wgh)),b_factor,sign);			  
						 
						  break;      
					}
				}
			}
		}
	}
	if (!complex_input) {
		fp->do_ift_inplace();
		fp->postift_depad_corner_inplace();
	}

		// Return a complex (Fourier) filtered image
		// Note: fp and fimage are the _same_ object if doInPlace
		// is true, so in that case fimage has been filtered.
		// We always return an image (pointer), but if the function 
		// was called with doInPlace == true then the calling function
		// will probably ignore the return value.

		// ELSE Return a real-space filtered image
	fp->set_array_offsets(0,0,0);
	fp->done_data();
	if (doInPlace && !complex_input) {
		// copy workspace data into the original image
		float* orig = fimage->get_data();
		float* work = fp->get_data();
		for (int i = 0; i < nx*ny*nz; i++) orig[i] = work[i];
		fimage->done_data();
	}
	return fp;
	EXITFUNC;
}

/* vim: set ts=4 noet: */
