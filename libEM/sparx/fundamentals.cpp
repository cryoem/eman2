#include "emdata.h"

using namespace EMAN;

namespace EMAN {
/*
******************************************************
*DISCLAIMER
* 03/29/05 P.A.Penczek
* The University of Texas
* Pawel.A.Penczek@uth.tmc.edu
* Please do not modify the content of this document
* without a written permission of the author.
*/
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
	fourierproduct(EMData* f, EMData* g, fp_flag flag, fp_type ptype) {
		int normfact;
		// Not only does the value of "flag" determine how we handle
		// periodicity, but it also determines whether or not we should
		// normalize the results.  Here's some convenience bools:
		bool donorm = (0 == flag%2) ? true : false;
		bool dopad  = (flag >= 3) ? true : false;
		// g may be NULL.  If so, have g point to the same object as f.  In that
		// case we need to be careful later on not to try to delete g's workspace
		// as well as f's workspace, since they will be the same.
		bool  gexists = true;
		if (!g) { g = f; gexists = false; }
		if (!equalsize(f, g)) {
			LOGERR("Fourier product requires congruent images");
			throw ImageDimensionException("Fourier product requires congruent images");
		}
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

		int npad = (flag>2) ? 2 : 1; // amount of padding used
		// the 2x padding is hardcoded for now
		// these are padded dimensions
		const int nxp = npad*nx;
		const int nyp = (ny > 1) ? npad*ny : 1; // don't pad y for 1-d image
		const int nzp = (nz > 1) ? npad*nz : 1; // don't pad z for 2-d image

		// now one half of the padded, fft-extended size along x
		const int lsd2 = (nxp + 2 - nxp%2) / 2; 
		// The [padded] fft-extended fourier version of f is fp.

		EMData* fp = NULL;
		if (f->is_complex()) { 
			// If f is already a fourier object then fp is a copy of f.
			// (The fp workspace is modified, so we copy f to keep f pristine.)
			fp=f->copy();
		} else {
			//  [normalize] [pad] compute fft
			fp = norm_pad_ft(f, donorm, dopad); 
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
			gp = norm_pad_ft(g, donorm, dopad);
		}
		// Get complex matrix views of fp and gp; matrices start from 1 (not 0)
		fp->set_array_offsets(1,1,1);
		gp->set_array_offsets(1,1,1);
		//  Multiply two functions (the real work of this routine)
		int itmp= nx/2;
		float sx  = float(-twopi*float(itmp)/float(nxp));
		itmp= ny/2;
		float sy  = float(-twopi*float(itmp)/float(nyp));
		itmp= nz/2;
		float sz  = float(-twopi*float(itmp)/float(nzp));
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
							fp->cmplx(ix,iy,iz)=
								(fpr*fpr + fpi*fpi)
								*std::complex<float>(cos(arg),sin(arg));
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
							fp->cmplx(ix,iy,iz) = 
								abs(fp->cmplx(ix,iy,iz))
								*std::complex<float>(cos(arg),sin(arg));
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
							fp->cmplx(ix,iy,iz) *=
								conj(gp->cmplx(ix,iy,iz));
								//*std::complex<float>(cos(arg),sin(arg));
						}
					}
				}
				break;
			case CONVOLUTION:
				// fpmat:=fpmat*gpmat
				// Note nxp are padded dimensions
				for (int iz = 1; iz <= nzp; iz++) {
					int jz=iz-1; if(jz>nzp/2) jz=jz-nzp; float argz=sz*jz;
					for (int iy = 1; iy <= nyp; iy++) {
						int jy=iy-1; if(jy>nyp/2) jy=jy-nyp; float argy=sy*jy+argz;
						for (int ix = 1; ix <= lsd2; ix++) {
							int jx=ix-1; float arg=sx*jx+argy;
							fp->cmplx(ix,iy,iz) *=
								gp->cmplx(ix,iy,iz)
								*std::complex<float>(cos(arg),sin(arg));
						}
					}
				}
				break;
			default:
				LOGERR("Illegal option in Fourier Product");
				throw InvalidValueException(ptype, "Illegal option in Fourier Product");
		}
		// Now done w/ gp, so let's get rid of it (if it's not an alias of fp);
		if (gexists && (f != g)) 
		{
			if( gp )
			{
				delete gp;
				gp = 0;
			}
		}
		// back transform
		fp->do_ift_inplace();
		fp->postift_depad_corner_inplace();

		vector<int> saved_offsets = fp->get_array_offsets();
		fp->set_array_offsets(1,1,1);

		normfact = (nxp/nx)*(nyp/ny)*(nzp/nz);  // Normalization factor for the padded operations
		if(normfact>1) {
			for (int iz = 1; iz <= nz; iz++) {
				for (int iy = 1; iy <= ny; iy++) {
					for (int ix = 1; ix <= nx; ix++) {
						(*fp)(ix,iy,iz) *= normfact;
					}
				}
			}
		}
		// Lag normalization
		if(flag>4)  {
			normfact = nx*ny*nz;  // Normalization factor
			int nxc=nx/2+1, nyc=ny/2+1, nzc=nz/2+1;
			for (int iz = 1; iz <= nz; iz++) {
				float lagz=float(normfact/(nz-abs(iz-nzc)));
				for (int iy = 1; iy <= ny; iy++) {
					float lagyz=lagz/(ny-abs(iy-nyc));
					for (int ix = 1; ix <= nx; ix++) {
						(*fp)(ix,iy,iz) *= lagyz/(nx-abs(ix-nxc));
					}
				}
			}	
		}
		//OVER AND OUT
		fp->set_array_offsets(saved_offsets);
		fp->done_data();
		return fp;  
	}

	bool equalsize(EMData* f, EMData* g) {
		if (g == f) return true;
		return ((f->get_xsize() == g->get_xsize()) &&
				(f->get_ysize() == g->get_ysize()) &&
				(f->get_zsize() == g->get_zsize()));
	}
}

/* vim: set ts=4 noet: */
