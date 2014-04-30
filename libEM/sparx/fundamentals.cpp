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

using namespace EMAN;
using std::complex;

namespace EMAN {

EMData* 
periodogram(EMData* f) {
	// These are actual dimensions
	int nx  = f->get_xsize();
	int ny  = f->get_ysize();
	int nz  = f->get_zsize();

 		// We manifestly assume no zero-padding here, just the
		// necessary extension along x for the fft

	if (f->is_complex()) nx = (nx - 2 + f->is_fftodd()); // nx is the real-space size of the input image
	int lsd2 = (nx + 2 - nx%2) / 2; // Extended x-dimension of the complex image

//  Process f if real
	EMData* fp = NULL;
	if(f->is_complex()) fp = f->copy(); // we need to make a full copy so that we don't damage the original
	else {
		
		fp = f->norm_pad(false, 1); // Extend and do the FFT if f is real
		fp->do_fft_inplace();


	}
	fp->set_array_offsets(1,1,1);

	//  Periodogram: fp:=|fp|**2
	for (int iz = 1; iz <= nz; iz++) {
		for (int iy = 1; iy <= ny; iy++) {
			for (int ix = 1; ix <= lsd2; ix++) {
				float fpr = real(fp->cmplx(ix,iy,iz));
				float fpi = imag(fp->cmplx(ix,iy,iz));
				fp->cmplx(ix,iy,iz) = fpr*fpr + fpi*fpi;
			}
		}
	}
	//  Create power as a 3D array (-n/2:n/2+n%2-1)
	int nyt, nzt;
	int nx2 = nx/2;
	int ny2 = ny/2; if(ny2 == 0) nyt =0; else nyt=ny;
	int nz2 = nz/2; if(nz2 == 0) nzt =0; else nzt=nz;
	int nx2p = nx2+nx%2;
	int ny2p = ny2+ny%2;
	int nz2p = nz2+nz%2;
	EMData& power = *(new EMData()); // output image
	power.set_size(nx, ny, nz);
	power.set_array_offsets(-nx2,-ny2,-nz2);
//If instead of preservation of the norm one would prefer to have peak of a PW of a single sine wave equal one
//                             multiply power by the scale below, or the other way around.
	float scale = 4.0f/float (nx*nx)/float (ny*ny)/float (nz*nz);
	for (int iz = 1; iz <= nz; iz++) {
		int jz=iz-1; 
		if(jz>=nz2p) jz=jz-nzt;
		for (int iy = 1; iy <= ny; iy++) {
			int jy=iy-1; 
			if(jy>=ny2p) jy=jy-nyt;
			for (int ix = 1; ix <= lsd2; ix++) {
				int jx=ix-1;
				if(jx>=nx2p) jx=jx-nx;
				power(jx,jy,jz) = real(fp->cmplx(ix,iy,iz)) * scale;
			}
		}
	}

//  Create the Friedel related half
	int  nzb, nze, nyb, nye, nxb, nxe;
	nxb =-nx2+(nx+1)%2;
	nxe = nx2-(nx+1)%2;
	if(ny2 == 0) {nyb =0; nye = 0;} else {nyb =-ny2+(ny+1)%2; nye = ny2-(ny+1)%2;}
	if(nz2 == 0) {nzb =0; nze = 0;} else {nzb =-nz2+(nz+1)%2; nze = nz2-(nz+1)%2;}
	for (int iz = nzb; iz <= nze; iz++) {
		for (int iy = nyb; iy <= nye; iy++) {
			for (int ix = 1; ix <= nxe; ix++) { // Note this loop begins with 1 - FFT should create correct Friedel related 0 plane
				power(-ix,-iy,-iz) = power(ix,iy,iz);
			}
		}
	}
	if(ny2 != 0)  {
		if(nz2 != 0)  {
			if(nz%2 == 0) {  //if nz even, fix the first slice
				for (int iy = nyb; iy <= nye; iy++) {
					for (int ix = nxb; ix <= -1; ix++) {
						power(ix,iy,-nz2) = power(-ix,-iy,-nz2);
					}
				}
				if(ny%2 == 0) {  //if ny even, fix the first line
					for (int ix = nxb; ix <= -1; ix++) {
						power(ix,-ny2,-nz2) = power(-ix,-ny2,-nz2);
					}
				}
			}
		}
		if(ny%2 == 0) {  //if ny even, fix the first column
			for (int iz = nzb; iz <= nze; iz++) {
				for (int ix = nxb; ix <= -1; ix++) {
					power(ix,-ny2,-iz) = power(-ix,-ny2,iz);
				}
			}
		}
		
	}

	if( fp ) {
		delete fp; // avoid a memory leak!
		fp = 0;
	}
	//power[0][0][0]=power[1][0][0];  //Steve requested the original origin.
	
	int sz[3];
	sz[0] = nx;
	sz[1] = ny;
	sz[2] = nz;
	int max_size = *std::max_element(&sz[0],&sz[3]);
	// set the pixel size for the power spectrum, only ration of the frequency pixel size is considered 	
	power.set_attr("apix_x", float(max_size)/nx);
	if(ny2 > 0) power.set_attr("apix_y", float(max_size)/ny);
	if(nz2 > 0) power.set_attr("apix_z", float(max_size)/nz);
	
	power.update();
	power.set_array_offsets(0,0,0);
	return &power;
//OVER AND OUT
}


/*
fourierproduct
Purpose: Calculate various fundamental image processing operations 
         based on Fourier products for 1-2-3D images. 
Method: Calculate FFT (if needed),
	Real f input -> padded fp workspace
	Complex f input -> set fp to f
	Real g input -> padded gp workspace
	Complex g input -> set gp to g
	No g input -> set gp tp fp
	 -> real output -> if f real delete workspace fp and if g real delete workspace gp.
Input: f real or complex 1-2-3D image, g - real or complex 1-2-3D image.
Output: 1-2-3D real image with the result
        (correlation F*conjg(G), convolution F*G, 
		autocorrelation |F|^2, self-correlation |F|)
        and with the origin at the image center int[n/2]+1.
*/
	EMData* 
	fourierproduct(EMData* f, EMData* g, fp_flag flag, fp_type ptype, bool center) {
		//std::complex<float> phase_mult;
		// Not only does the value of "flag" determine how we handle
		// periodicity, but it also determines whether or not we should
		// normalize the results.  Here's some convenience bools:
		bool donorm = (0 == flag%2) ? true : false;
		// the 2x padding is hardcoded for now
		int  npad  = (flag >= 3) ? 2 : 1;  // amount of padding used
		// g may be NULL.  If so, have g point to the same object as f.  In that
		// case we need to be careful later on not to try to delete g's workspace
		// as well as f's workspace, since they will be the same.
		bool  gexists = true;
		if (!g) { g = f; gexists = false; }
		if ( f->is_complex() || g->is_complex() ) {
			// Fourier input only allowed for circulant
			if (CIRCULANT != flag) {
				LOGERR("Cannot perform normalization or padding on Fourier type.");
				throw InvalidValueException(flag, "Cannot perform normalization or padding on Fourier type.");
			}
		}
		// These are actual dimensions of f (and real-space sizes for ny and nz)
		int nx  = f->get_xsize();
		int ny  = f->get_ysize();
		int nz  = f->get_zsize();
		// We manifestly assume no zero-padding here, just the 
		// necessary extension along x for the fft
		if (!f->is_real()) nx = (nx - 2 + (f->is_fftodd() ? 1 : 0)); 

		// these are padded dimensions
		const int nxp = npad*nx;
		const int nyp = (ny > 1) ? npad*ny : 1; // don't pad y for 1-d image
		const int nzp = (nz > 1) ? npad*nz : 1; // don't pad z for 2-d image

		// now one half of the padded, fft-extended size along x
		const int lsd2 = (nxp + 2 - nxp%2) / 2; 

		EMData* fp = NULL;
		if (f->is_complex()) { 
			// If f is already a fourier object then fp is a copy of f.
			// (The fp workspace is modified, so we copy f to keep f pristine.)
			fp=f->copy();
		} else {
			//  [normalize] [pad] compute fft
			fp = f->norm_pad(donorm, npad);
			if (donorm) {
				fp->div( sqrtf(nx) * sqrtf(ny) * sqrtf(nz) );
			}
			fp->do_fft_inplace();
		}
		// The [padded] fft-extended version of g is gp.
		EMData* gp = NULL;
		if(f==g) {
			// g is an alias for f, so gp should be an alias for fp
			gp=fp;
		} else if (g->is_complex()) {
			// g is already a Fourier object, so gp is just an alias for g
			// (The gp workspace is not modified, so we don't need a copy.)
			gp = g;
		} else {
			// normal case: g is real and different from f, so compute gp
			gp = g->norm_pad(donorm, npad);
			if (donorm) {
				gp->div( sqrtf(nx) * sqrtf(ny) * sqrtf(nz) );
			}
			gp->do_fft_inplace();
		}
		// Get complex matrix views of fp and gp; matrices start from 1 (not 0)
		fp->set_array_offsets(1,1,1);
		gp->set_array_offsets(1,1,1);

		// If the center flag is true, put the center of the correlation in the middle
		// If it is false, put it in (0,0), this approach saves time, but it is diffcult to manage the result
		if (center) {
			//  Multiply two functions (the real work of this routine)
			int itmp = nx/2;
			//float sx  = float(-twopi*float(itmp)/float(nxp));
			float sxn = 2*float(itmp)/float(nxp);
			float sx = -M_PI*sxn;
			itmp = ny/2;
			//float sy  = float(-twopi*float(itmp)/float(nyp));
			float syn = 2*float(itmp)/float(nyp);
			float sy = -M_PI*syn;
			itmp = nz/2;
			//float sz  = float(-twopi*float(itmp)/float(nzp));
			float szn = 2*float(itmp)/float(nzp);
			float sz = -M_PI*szn;
			if ( (flag > 2) || (nx%2==0 && (ny%2==0 || ny==1) && (nz%2==0 || nz==1)) ) {  // padded images have always even size:  if ( (padded) || (even) ) ...
				switch (ptype) {
					case AUTOCORRELATION:
					// fpmat := |fpmat|^2
					// Note nxp are padded dimensions
						for (int iz = 1; iz <= nzp; iz++) {
							for (int iy = 1; iy <= nyp; iy++) {
								for (int ix = 1; ix <= lsd2; ix++) {
									float fpr = real(fp->cmplx(ix,iy,iz));
									float fpi = imag(fp->cmplx(ix,iy,iz));
									fp->cmplx(ix,iy,iz) = complex<float>(fpr*fpr+fpi*fpi, 0.0f);
								}
							}
						}
						break;
					case SELF_CORRELATION:
					// fpmat:=|fpmat|
					// Note nxp are padded dimensions
						for (int iz = 1; iz <= nzp; iz++) {
							for (int iy = 1; iy <= nyp; iy++) {
								for (int ix = 1; ix <= lsd2; ix++) {
									fp->cmplx(ix,iy,iz) = complex<float>(abs(fp->cmplx(ix,iy,iz)), 0.0f);
								}
							}
						}
						break;
					case CORRELATION:
					// fpmat:=fpmat*conjg(gpmat)
					// Note nxp are padded dimensions
						for (int iz = 1; iz <= nzp; iz++) {
							for (int iy = 1; iy <= nyp; iy++) {
								for (int ix = 1; ix <= lsd2; ix++) {
									fp->cmplx(ix,iy,iz) *= conj(gp->cmplx(ix,iy,iz));
								}
							}
						}
					break;
					case CONVOLUTION:
					// fpmat:=fpmat*gpmat
					// Note nxp are padded dimensions
						for (int iz = 1; iz <= nzp; iz++) {
							for (int iy = 1; iy <= nyp; iy++) {
								for (int ix = 1; ix <= lsd2; ix++) {
									fp->cmplx(ix,iy,iz) *= gp->cmplx(ix,iy,iz);
								}
							}
						}
						break;
					default:
						LOGERR("Illegal option in Fourier Product");
						throw InvalidValueException(ptype, "Illegal option in Fourier Product");
				}					
				for (int iz = 1; iz <= nzp; iz++) {
					for (int iy = 1; iy <= nyp; iy++) {
						for (int ix = (iz+iy+1)%2+1; ix <= lsd2; ix+=2) {
							fp->cmplx(ix,iy,iz) = -fp->cmplx(ix,iy,iz);
						}
					}
				}
			} else {
				switch (ptype) {
					case AUTOCORRELATION:
					// fpmat := |fpmat|^2
					// Note nxp are padded dimensions
						for (int iz = 1; iz <= nzp; iz++) {
						int jz=iz-1; if(jz>nzp/2) jz=jz-nzp; float argz=sz*jz;
							for (int iy = 1; iy <= nyp; iy++) {
							int jy=iy-1; if(jy>nyp/2) jy=jy-nyp; float argy=sy*jy+argz;
								for (int ix = 1; ix <= lsd2; ix++) {
									int jx=ix-1; float arg=sx*jx+argy;
									float fpr = real(fp->cmplx(ix,iy,iz));
									float fpi = imag(fp->cmplx(ix,iy,iz));
									fp->cmplx(ix,iy,iz)= (fpr*fpr + fpi*fpi) *std::complex<float>(cos(arg),sin(arg));
								}
							}
						}
						break;
					case SELF_CORRELATION:
					// fpmat:=|fpmat|
					// Note nxp are padded dimensions
						for (int iz = 1; iz <= nzp; iz++) {
						int jz=iz-1; if(jz>nzp/2) jz=jz-nzp; float argz=sz*jz;
							for (int iy = 1; iy <= nyp; iy++) {
							int jy=iy-1; if(jy>nyp/2) jy=jy-nyp; float argy=sy*jy+argz;
								for (int ix = 1; ix <= lsd2; ix++) {
									int jx=ix-1; float arg=sx*jx+argy;
									fp->cmplx(ix,iy,iz) = abs(fp->cmplx(ix,iy,iz)) *std::complex<float>(cos(arg),sin(arg));
								}
							}
						}
						break;
					case CORRELATION:
					// fpmat:=fpmat*conjg(gpmat)
					// Note nxp are padded dimensions
						for (int iz = 1; iz <= nzp; iz++) {
						int jz=iz-1; if(jz>nzp/2) jz=jz-nzp; float argz=sz*jz;
							for (int iy = 1; iy <= nyp; iy++) {
							int jy=iy-1; if(jy>nyp/2) jy=jy-nyp; float argy=sy*jy+argz;
								for (int ix = 1; ix <= lsd2; ix++) {
									int jx=ix-1; float arg=sx*jx+argy;
									fp->cmplx(ix,iy,iz) *= conj(gp->cmplx(ix,iy,iz)) *std::complex<float>(cos(arg),sin(arg));
								}
							}
						}
					break;
					case CONVOLUTION:
					// fpmat:=fpmat*gpmat
					// Note nxp are padded dimensions
						if(npad == 1) {
							sx -= 4*(nx%2)/float(nx);
							sy -= 4*(ny%2)/float(ny);
							sz -= 4*(nz%2)/float(nz);
						}
						for (int iz = 1; iz <= nzp; iz++) {
							int jz=iz-1; if(jz>nzp/2) jz=jz-nzp; float argz=sz*jz;
							for (int iy = 1; iy <= nyp; iy++) {
								int jy=iy-1; if(jy>nyp/2) jy=jy-nyp; float argy=sy*jy+argz;
								for (int ix = 1; ix <= lsd2; ix++) {
									int jx=ix-1; float arg=sx*jx+argy;
									fp->cmplx(ix,iy,iz) *= gp->cmplx(ix,iy,iz) *std::complex<float>(cos(arg),sin(arg));
								}
							}
						}
						break;
					default:
						LOGERR("Illegal option in Fourier Product");
						throw InvalidValueException(ptype, "Illegal option in Fourier Product");
				}
			}
		} else {
			// If the center flag is false, then just do basic multiplication
			// Here I aterd the method of complex calculation. This method is much faster than the previous one.
			switch (ptype) {
				case AUTOCORRELATION:
					for (int iz = 1; iz <= nzp; iz++) {
						for (int iy = 1; iy <= nyp; iy++) {
							for (int ix = 1; ix <= lsd2; ix++) {
								float fpr = real(fp->cmplx(ix,iy,iz));
								float fpi = imag(fp->cmplx(ix,iy,iz));
								fp->cmplx(ix,iy,iz) = complex<float>(fpr*fpr+fpi*fpi, 0.0f);
							}
						}
					}
					break;
				case SELF_CORRELATION:
					for (int iz = 1; iz <= nzp; iz++) {
						for (int iy = 1; iy <= nyp; iy++) {
							for (int ix = 1; ix <= lsd2; ix++) {
								fp->cmplx(ix,iy,iz) = complex<float>(abs(fp->cmplx(ix,iy,iz)), 0.0f);
							}
						}
					}
					break;
				case CORRELATION:
					//phase_mult = 1;
					for (int iz = 1; iz <= nzp; iz++) {
						for (int iy = 1; iy <= nyp; iy++) {
							for (int ix = 1; ix <= lsd2; ix++) {
								fp->cmplx(ix,iy,iz)*= conj(gp->cmplx(ix,iy,iz));
							}
						}
					}
					break;
				case CONVOLUTION:
					if(npad == 1) {
						float sx = -M_PI*2*(nx%2)/float(nx);
						float sy = -M_PI*2*(ny%2)/float(ny);
						float sz = -M_PI*2*(nz%2)/float(nz);
						for (int iz = 1; iz <= nzp; iz++) {
							int jz=iz-1; if(jz>nzp/2) jz=jz-nzp; float argz=sz*jz;
							for (int iy = 1; iy <= nyp; iy++) {
								int jy=iy-1; if(jy>nyp/2) jy=jy-nyp; float argy=sy*jy+argz;
								for (int ix = 1; ix <= lsd2; ix++) {
									int jx=ix-1; float arg=sx*jx+argy;
									fp->cmplx(ix,iy,iz) *= gp->cmplx(ix,iy,iz) *std::complex<float>(cos(arg),sin(arg));
								}
							}
						}
					} else {
						for (int iz = 1; iz <= nzp; iz++) {
							for (int iy = 1; iy <= nyp; iy++) {
								for (int ix = 1; ix <= lsd2; ix++) {
									fp->cmplx(ix,iy,iz)*= gp->cmplx(ix,iy,iz);
								}
							}
						}
					}
					break;
				default:
					LOGERR("Illegal option in Fourier Product");
					throw InvalidValueException(ptype, "Illegal option in Fourier Product");
			}
		}
		// Now done w/ gp, so let's get rid of it (if it's not an alias of fp or simply g was complex on input);
		if (gexists && (f != g) && (!g->is_complex())) {
			if( gp ) {
				delete gp;
				gp = 0;
			}
		}
		// back transform
		fp->do_ift_inplace();
		if(center && npad ==2)  fp->depad();
		else                    fp->depad_corner();

		//vector<int> saved_offsets = fp->get_array_offsets();  I do not know what the meaning of it was, did not work anyway PAP
		fp->set_array_offsets(1,1,1);

		// Lag normalization
		if(flag>4)  {
			int nxc=nx/2+1, nyc=ny/2+1, nzc=nz/2+1;
			for (int iz = 1; iz <= nz; iz++) {
				const float lagz = float(nz) / (nz-abs(iz-nzc));
				for (int iy = 1; iy <= ny; iy++) {
					const float lagyz = lagz * (float(ny) / (ny-abs(iy-nyc)));
					for (int ix = 1; ix <= nx; ix++) {
						(*fp)(ix,iy,iz) *= lagyz * (float(nx) / (nx-abs(ix-nxc)));
					}
				}
			}
		}
		//OVER AND OUT
		//fp->set_array_offsets(saved_offsets);  This was strange and did not work, PAP
		fp->set_array_offsets(0,0,0);
		fp->update();
		return fp;
	}
}

