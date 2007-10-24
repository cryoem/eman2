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
using std::vector;
using std::complex;

void EMData::center_origin()
{
	ENTERFUNC;
	if (is_complex()) {
		LOGERR("Real image expected. Input image is complex.");
		throw ImageFormatException("Real image expected. Input image is complex.");
	}
	for (int iz = 0; iz < nz; iz++) {
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				// next line multiplies by +/- 1
				(*this)(ix,iy,iz) *= -2*((ix+iy+iz)%2) + 1;
			}
		}
	}
	update();
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
	vector<int> saved_offsets = get_array_offsets();
	// iz in [1,nz], iy in [1,ny], ix in [0,nx/2], but nx comes in as extended and is the same for odd
	//                                                 and even, so we can ignore the difference...
	//                         in short, as nx is extended, it should be  ix in [0,(nx-2)/2],  corrected PAP 05/20
	set_array_offsets(0,1,1);
	int nxc = nx/2;

	if (is_fftodd()) {
		for (int iz = 1; iz <= nz; iz++) {
			for (int iy = 1; iy <= ny; iy++) {
				for (int ix = 0; ix < nxc; ix++) {
					cmplx(ix,iy,iz) *= float(-2*((ix+iy+iz)%2) + 1);
					float temp = float(iz-1+iy-1+ix)/float(ny)*M_PI; 		
					complex<float> temp2 = complex<float>(cos(temp), -sin(temp));
					cmplx(ix,iy,iz) *= temp2;
				}
			}
		}
	} else {
		for (int iz = 1; iz <= nz; iz++) {
			for (int iy = 1; iy <= ny; iy++) {
				for (int ix = 0; ix < nxc; ix++) {
					// next line multiplies by +/- 1
					cmplx(ix,iy,iz) *= float(-2*((ix+iy+iz)%2) + 1);
				}
			}
		}
	}
	set_array_offsets(saved_offsets);
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
	
	// Commented by Zhengfan Yang on 03/14/07
	// This is incorrect for npad=1 and nx (or ny, nz) is odd
	// int xstart = (nx != 1) ? (nxpad - nx)/2 + nx%2 : 0;
	// int ystart = (ny != 1) ? (nypad - ny)/2 + ny%2 : 0;
	// int zstart = (nz != 1) ? (nzpad - nz)/2 + nz%2 : 0;
	// Should have better way to write it, but temporarily we'll use the following
	int xstart = ( nx == 1 || npad == 1 ) ? 0 : (nxpad - nx)/2 + nx%2;
	int ystart = ( ny == 1 || npad == 1 ) ? 0 : (nypad - ny)/2 + ny%2;
	int zstart = ( nz == 1 || npad == 1 ) ? 0 : (nzpad - nz)/2 + nz%2;
	
	for (int iz = 0; iz < nz; iz++) for (int iy = 0; iy < ny; iy++) memcpy( &(*newimg)(xstart,iy+ystart,iz+zstart), &(*this)(0,iy,iz), bytes);
	newimg->update();
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
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
	EMData* newimg = copy_head();
	//newimg->to_zero();  makes no sense, no size set
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
		if(npad > 1)  newimg->to_zero();  // If only extension for FFT, no need to zero...
		newimg->set_fftpad(true);
		newimg->set_attr("npad", npad);
		if (offset == 1) newimg->set_fftodd(true);
		else             newimg->set_fftodd(false);
		for (int iz = 0; iz < nz; iz++) {
			for (int iy = 0; iy < ny; iy++) {
				memcpy(&(*newimg)(0,iy,iz), &(*this)(0,iy,iz), bytes);
			}
		}
	} else {
		// Image already padded, so we want to remove the padding
		// (Note: The npad passed in is ignored in favor of the one
		//  stored in the image.)
		npad = get_attr("npad");
		if (0 == npad) npad = 1;
		int nxold = (nx - 2 + int(is_fftodd()))/npad; 
	#ifdef _WIN32
		int nyold = _cpp_max(ny/npad, 1);
		int nzold = _cpp_max(nz/npad, 1);
	#else
		int nyold = std::max<int>(ny/npad, 1);
		int nzold = std::max<int>(nz/npad, 1);
	#endif	//_WIN32
		int bytes = nxold*sizeof(float);
		newimg->set_size(nxold, nyold, nzold);
		newimg->set_fftpad(false);
		for (int iz = 0; iz < nzold; iz++) {
			for (int iy = 0; iy < nyold; iy++) {
				memcpy(&(*newimg)(0,iy,iz), &(*this)(0,iy,iz), bytes);
			}
		}
	}
	newimg->update();
	set_array_offsets(saved_offsets);
	return newimg;
}


void EMData::postift_depad_corner_inplace() {
	ENTERFUNC;
 	if(is_fftpadded() == true) {
		vector<int> saved_offsets = get_array_offsets();
		set_array_offsets(0,0,0);
		int npad = attr_dict["npad"];
		if (0 == npad) npad = 1;
		int offset = is_fftodd() ? 1 : 2;
		int nxold = (nx - offset)/npad;
#ifdef _WIN32
		int nyold = _cpp_max(ny/npad, 1);
		int nzold = _cpp_max(nz/npad, 1);
#else
		int nyold = std::max<int>(ny/npad, 1);
		int nzold = std::max<int>(nz/npad, 1);
#endif  //_WIN32
		int bytes = nxold*sizeof(float);
		float* dest = get_data();
		for (int iz=0; iz < nzold; iz++) {
			for (int iy = 0; iy < nyold; iy++) {
				memmove(dest, &(*this)(0,iy,iz), bytes);
				dest += nxold;
			}
		}
		set_attr("npad", 1);
		set_size(nxold, nyold, nzold);
		set_fftpad(false);
		update();
		set_complex(false);
		if(ny==1 && nz==1) set_complex_x(false);
		set_array_offsets(saved_offsets);
	}
	EXITFUNC;
}


#define  fint(i,j,k)  fint[(i-1) + ((j-1) + (k-1)*ny)*lsd]
#define  fout(i,j,k)  fout[(i-1) + ((j-1) + (k-1)*nyn)*lsdn]
EMData *EMData::FourInterpol(int nxn, int nyni, int nzni, bool RetReal) {

	int nyn, nzn, lsd, lsdn, inx, iny, inz;
	int i, j, k;
	if (is_complex()) 
		throw ImageFormatException("Input image has to be real");

	if(ny > 1) {
		nyn = nyni;
		if(nz > 1) {
			nzn = nzni;
		}  else {
			nzn = 1;
		}
	} else {
		nyn = 1; nzn = 1;
	}
	if(nxn<nx || nyn<ny || nzn<nz)	throw ImageDimensionException("Cannot reduce the image size");
	lsd = nx + 2 - nx%2;
	lsdn = nxn + 2 - nxn%2;
//  do out of place ft
	EMData *temp_ft = do_fft();
	EMData *ret = this->copy();
	ret->set_size(lsdn, nyn, nzn);
	ret->to_zero();
	float *fout = ret->get_data();
	float *fint = temp_ft->get_data();
//  TO KEEP EXACT VALUES ON THE ORIGINAL GRID ONE SHOULD USE
//  SQ2     = 2.0. HOWEVER, TOTAL ENERGY WILL NOT BE CONSERVED
	float  sq2 = 1.0f/std::sqrt(2.0f);
	float  anorm = (float) nxn* (float) nyn* (float) nzn/(float) nx/ (float) ny/ (float) nz;
	for (i = 0; i < lsd*ny*nz; i++)  fint[i] *= anorm;
	inx = nxn-nx; iny = nyn - ny; inz = nzn - nz;
	for (k=1; k<=nz/2+1; k++) for (j=1; j<=ny/2+1; j++) for (i=1; i<=lsd; i++) fout(i,j,k)=fint(i,j,k);
	if(nyn>1) {
	//cout << "  " <<nxn<<"  " <<nyn<<" A " <<nzn<<endl;
		for (k=1; k<=nz/2+1; k++) for (j=ny/2+2+iny; j<=nyn; j++) for (i=1; i<=lsd; i++) fout(i,j,k)=fint(i,j-iny,k);
		if(nzn>1) {
			for (k=nz/2+2+inz; k<=nzn; k++) {
				for (j=1; j<=ny/2+1; j++) {
					for (i=1; i<=lsd; i++) {
						fout(i,j,k)=fint(i,j,k-inz);
					}
				}
				for (j=ny/2+2+iny; j<=nyn; j++) {
					for (i=1; i<=lsd; i++) {
						fout(i,j,k)=fint(i,j-iny,k-inz);
					}
				}
			}
		}
	}
//       WEIGHTING FACTOR USED FOR EVEN NSAM. REQUIRED SINCE ADDING ZERO FOR
//       INTERPOLATION WILL INTRODUCE A COMPLEX CONJUGATE FOR NSAM/2'TH
//       ELEMENT.
	if(nx%2 == 0 && inx !=0) {
		for (k=1; k<=nzn; k++) {
			for (j=1; j<=nyn; j++) {
				fout(nx+1,j,k) *= sq2;
				fout(nx+2,j,k) *= sq2;
			}
		}
		if(nyn>1) {
			for (k=1; k<=nzn; k++) {
			  for (i=1; i<=lsd; i++) {
			    fout(i,ny/2+1+iny,k) = sq2*fout(i,ny/2+1,k);
			    fout(i,ny/2+1,k) *= sq2;
			  }
			}
			if(nzn>1) {
				for (j=1; j<=nyn; j++) {
					for (i=1; i<=lsd; i++) {
						fout(i,j,nz/2+1+inz) = sq2*fout(i,j,nz/2+1);
						fout(i,j,nz/2+1) *= sq2;
					}
				}
			}
		}
	}
	ret->set_complex(true);
/*
//  For padding from odd to even dimension additional shift by 1 pixel is necessary.
	float  xshift = 0.f, yshift = 0.f, zshift = 0.f;
	int nyn2, nzn2;
	if(nxn > nx && nx%2 == 1)  xshift = 1.0f;
	if(ny > 1) {
		if(nyn > ny && ny%2 == 1)  yshift = 1.0f;
		nyn2 = nyn/2;
		if(nz > 1) {
			if(nzn > nz && nz%2 == 1)  zshift = 1.0f;
			nzn2 = nzn/2;
		}  else {
			nzn2 = 0;
		}
	} else {
		nyn2 = 0; nzn2 = 0;
	}
	if(xshift == 1.0 || yshift == 1.0 || zshift == 1.0)  {
		ret->set_array_offsets(1,1,1);
		int  lsdn2 = lsd/2;
		for (int iz = 1; iz <= nzn; iz++) {
			int jz=iz-1; if(jz>nzn2) jz=jz-nzn;
			for (int iy = 1; iy <= nyn; iy++) {
				int jy=iy-1; if(jy>nyn2) jy=jy-nyn;
				for (int ix = 1; ix <= lsdn2; ix++) {
					int jx=ix-1;
					ret->cmplx(ix,iy,iz) *= 
					exp(-float(twopi)*iimag*(xshift*jx/nxn + yshift*jy/nyn+ zshift*jz/nzn));
				}
			}
		}
		ret->set_array_offsets(0,0,0);
	}*/
	ret->set_ri(1);
	ret->set_fftpad(true);
	ret->set_attr("npad", 1);
	if (nxn%2 == 1) {ret->set_fftodd(true);} else {ret->set_fftodd(false);}
	if(RetReal) {
		ret->do_ift_inplace();
		ret->postift_depad_corner_inplace();
	}
	ret->update();
	
	/*Dict d1 = temp_ft->get_attr_dict();
	Dict d2 = ret->get_attr_dict();
	printf("-----------------Attribute Dict for temp_ft--------------\n");
	EMUtil::dump_dict(d1);
	printf("-----------------Attribute Dict for ret--------------\n");
	EMUtil::dump_dict(d2);*/
	delete temp_ft;
	temp_ft = 0;
	return ret;
}

EMData *EMData::FourTruncate(int nxn, int nyni, int nzni, bool RetReal) {

	int nyn, nzn, lsd, lsdn, inx, iny, inz;
	int i, j, k;
	if (is_complex()) 
		throw ImageFormatException("Input image has to be real");

	if(ny > 1) {
		nyn = nyni;
		if(nz > 1) {
			nzn = nzni;
		}  else {
			nzn = 1;
		}
	} else {
		nyn = 1; nzn = 1;
	}
	if(nxn>nx || nyn>ny || nzn>nz)	throw ImageDimensionException("Cannot increase the image size");
	lsd = nx + 2 - nx%2;
	lsdn = nxn + 2 - nxn%2;
//  do out of place ft
	EMData *temp_ft = do_fft();
	EMData *ret = this->copy();
	ret->set_size(lsdn, nyn, nzn);
	float *fout = ret->get_data();
	float *fint = temp_ft->get_data();
//  TO KEEP EXACT VALUES ON THE ORIGINAL GRID ONE SHOULD USE
//  SQ2     = 2.0. HOWEVER, TOTAL ENERGY WILL NOT BE CONSERVED
	//float  sq2 = std::sqrt(2.0f);
	float  anorm = (float) nxn* (float) nyn* (float) nzn/(float) nx/ (float) ny/ (float) nz;
	for (i = 0; i < lsd*ny*nz; i++)  fint[i] *= anorm;
	inx = nx - nxn;  iny = ny - nyn;  inz = nz - nzn;
	for (k=1; k<=nzn/2+1; k++) for (j=1; j<=nyn/2+1; j++) for (i=1; i<=lsdn; i++) fout(i,j,k)=fint(i,j,k);
	if(nyn>1) {
		for (k=1; k<=nzn/2+1; k++) for (j=nyn/2+2; j<=nyn; j++) for (i=1; i<=lsdn; i++) fout(i,j,k)=fint(i,j+iny,k);
		if(nzn>1) {
			for (k=nzn/2+2; k<=nzn; k++) {
				for (j=1; j<=nyn/2+1; j++) {
					for (i=1; i<=lsdn; i++) {
						fout(i,j,k)=fint(i,j,k+inz);
					}
				}
				for (j=nyn/2+2; j<=nyn; j++) {
					for (i=1; i<=lsdn; i++) {
						fout(i,j,k)=fint(i,j+iny,k+inz);
					}
				}
			}
		}
	}
//       WEIGHTING FACTOR USED FOR EVEN NSAM. REQUIRED SINCE ADDING ZERO FOR
//       INTERPOLATION WILL INTRODUCE A COMPLEX CONJUGATE FOR NSAM/2'TH
//       ELEMENT.
	/*
	if(nxn%2 == 0 && inx !=0) {
		for (k=1; k<=nzn; k++) {
			for (j=1; j<=nyn; j++) {
				fout(nxn+1,j,k) *= sq2;
				fout(nxn+2,j,k) *= sq2;
			}
		}
		if(nyn>1) {
			for (k=1; k<=nzn; k++) {
			  for (i=1; i<=lsdn; i++) {
			    fout(i,nyn/2+1+iny,k) = sq2*fout(i,nyn/2+1,k);
			    fout(i,nyn/2+1,k) *= sq2;
			  }
			}
			if(nzn>1) {
				for (j=1; j<=nyn; j++) {
					for (i=1; i<=lsdn; i++) {
						fout(i,j,nzn/2+1+inz) = sq2*fout(i,j,nzn/2+1);
						fout(i,j,nzn/2+1) *= sq2;
					}
				}
			}
		}
	}*/
	ret->set_complex(true);
	ret->set_ri(1);
	ret->set_fftpad(true);
	ret->set_attr("npad", 1);
	if (nxn%2 == 1) {ret->set_fftodd(true);} else {ret->set_fftodd(false);}
	if(RetReal) {
		ret->do_ift_inplace();
		ret->postift_depad_corner_inplace();
	}
	ret->update();
	
	/*Dict d1 = temp_ft->get_attr_dict();
	Dict d2 = ret->get_attr_dict();
	printf("-----------------Attribute Dict for temp_ft--------------\n");
	EMUtil::dump_dict(d1);
	printf("-----------------Attribute Dict for ret--------------\n");
	EMUtil::dump_dict(d2);*/
	delete temp_ft;
	temp_ft = 0;
	return ret;
}

EMData *EMData::FourInterpol_i(int nxn, int nyni, int nzni, bool RetReal) {

	int nyn, nzn, lsd, lsdn, inx, iny, inz;
	int i, j, k;

	if(ny > 1) {
		nyn = nyni;
		if(nz > 1) {
			nzn = nzni;
		}  else {
			nzn = 1;
		}
	} else {
		nyn = 1; nzn = 1;
	}
	if(nxn<nx || nyn<ny || nzn<nz)	throw ImageDimensionException("Cannot reduce the image size");
	lsd = nx-2 + 2 - nx%2;
	lsdn = nxn + 2 - nxn%2;
//  do out of place ft
	EMData *temp_ft = this->copy();
	EMData *ret = this->copy();
	ret->set_size(lsdn, nyn, nzn);
	ret->to_zero();
	float *fout = ret->get_data();
	float *fint = temp_ft->get_data();
//  TO KEEP EXACT VALUES ON THE ORIGINAL GRID ONE SHOULD USE
//  SQ2     = 2.0. HOWEVER, TOTAL ENERGY WILL NOT BE CONSERVED
	float  sq2 = 1.0f/std::sqrt(2.0f);
	float  anorm = (float) nxn* (float) nyn* (float) nzn/(float) nx/ (float) ny/ (float) nz;
	for (i = 0; i < lsd*ny*nz; i++)  fint[i] *= anorm;
	inx = nxn-(nx-2); iny = nyn - ny; inz = nzn - nz;
	for (k=1; k<=nz/2+1; k++) for (j=1; j<=ny/2+1; j++) for (i=1; i<=lsd; i++) fout(i,j,k)=fint(i,j,k);
	if(nyn>1) {
	//cout << "  " <<nxn<<"  " <<nyn<<" A " <<nzn<<endl;
		for (k=1; k<=nz/2+1; k++) for (j=ny/2+2+iny; j<=nyn; j++) for (i=1; i<=lsd; i++) fout(i,j,k)=fint(i,j-iny,k);
		if(nzn>1) {
			for (k=nz/2+2+inz; k<=nzn; k++) {
				for (j=1; j<=ny/2+1; j++) {
					for (i=1; i<=lsd; i++) {
						fout(i,j,k)=fint(i,j,k-inz);
					}
				}
				for (j=ny/2+2+iny; j<=nyn; j++) {
					for (i=1; i<=lsd; i++) {
						fout(i,j,k)=fint(i,j-iny,k-inz);
					}
				}
			}
		}
	}
//       WEIGHTING FACTOR USED FOR EVEN NSAM. REQUIRED SINCE ADDING ZERO FOR
//       INTERPOLATION WILL INTRODUCE A COMPLEX CONJUGATE FOR NSAM/2'TH
//       ELEMENT.
	if(nx%2 == 0 && inx !=0) {
		for (k=1; k<=nzn; k++) {
			for (j=1; j<=nyn; j++) {
				fout(nx-2+1,j,k) *= sq2;
				fout(nx-2+2,j,k) *= sq2;
			}
		}
		if(nyn>1) {
			for (k=1; k<=nzn; k++) {
			  for (i=1; i<=lsd; i++) {
			    fout(i,ny/2+1+iny,k) = sq2*fout(i,ny/2+1,k);
			    fout(i,ny/2+1,k) *= sq2;
			  }
			}
			if(nzn>1) {
				for (j=1; j<=nyn; j++) {
					for (i=1; i<=lsd; i++) {
						fout(i,j,nz/2+1+inz) = sq2*fout(i,j,nz/2+1);
						fout(i,j,nz/2+1) *= sq2;
					}
				}
			}
		}
	}
	ret->set_complex(true);
/*
//  For padding from odd to even dimension additional shift by 1 pixel is necessary.
	float  xshift = 0.f, yshift = 0.f, zshift = 0.f;
	int nyn2, nzn2;
	if(nxn > nx && nx%2 == 1)  xshift = 1.0f;
	if(ny > 1) {
		if(nyn > ny && ny%2 == 1)  yshift = 1.0f;
		nyn2 = nyn/2;
		if(nz > 1) {
			if(nzn > nz && nz%2 == 1)  zshift = 1.0f;
			nzn2 = nzn/2;
		}  else {
			nzn2 = 0;
		}
	} else {
		nyn2 = 0; nzn2 = 0;
	}
	if(xshift == 1.0 || yshift == 1.0 || zshift == 1.0)  {
		ret->set_array_offsets(1,1,1);
		int  lsdn2 = lsd/2;
		for (int iz = 1; iz <= nzn; iz++) {
			int jz=iz-1; if(jz>nzn2) jz=jz-nzn;
			for (int iy = 1; iy <= nyn; iy++) {
				int jy=iy-1; if(jy>nyn2) jy=jy-nyn;
				for (int ix = 1; ix <= lsdn2; ix++) {
					int jx=ix-1;
					ret->cmplx(ix,iy,iz) *= 
					exp(-float(twopi)*iimag*(xshift*jx/nxn + yshift*jy/nyn+ zshift*jz/nzn));
				}
			}
		}
		ret->set_array_offsets(0,0,0);
	}*/
	ret->set_ri(1);
	ret->set_fftpad(true);
	ret->set_attr("npad", 1);
	if (nxn%2 == 1) {ret->set_fftodd(true);} else {ret->set_fftodd(false);}
	if(RetReal) {
		ret->do_ift_inplace();
		ret->postift_depad_corner_inplace();
	}
	ret->update();
	
	/*Dict d1 = temp_ft->get_attr_dict();
	Dict d2 = ret->get_attr_dict();
	printf("-----------------Attribute Dict for temp_ft--------------\n");
	EMUtil::dump_dict(d1);
	printf("-----------------Attribute Dict for ret--------------\n");
	EMUtil::dump_dict(d2);*/
	delete temp_ft;
	temp_ft = 0;
	return ret;
}

EMData *EMData::Four_ds(int nxn, int nyni, int nzni, bool RetReal) {

	int nyn, nzn, lsd, lsdn, inx, iny, inz;
	int i, j;

	if(ny > 1) {
		nyn = nyni;
		if(nz > 1) {
			nzn = nzni;
		}  else {
			nzn = 1;
		}
	} else {
		nyn = 1; nzn = 1;
	}
	lsd = nx-2 + 2 - nx%2;
	lsdn = nxn + 2 - nxn%2;
//  do out of place ft
	EMData *temp_ft = this->copy();
	EMData *ret = this->copy();
	ret->set_size(lsdn, nyn, nzn);
	ret->to_zero();
	float *fout = ret->get_data();
	float *fint = temp_ft->get_data();
//  TO KEEP EXACT VALUES ON THE ORIGINAL GRID ONE SHOULD USE
//  SQ2     = 2.0. HOWEVER, TOTAL ENERGY WILL NOT BE CONSERVED
//	float  sq2 = 1.0f/std::sqrt(2.0f);
	float  anorm = (float) nxn* (float) nyn* (float) nzn/(float) nx/ (float) ny/ (float) nz;
	for (i = 0; i < lsd*ny*nz; i++)  fint[i] *= anorm;
	inx = nxn-(nx-2); iny = nyn - ny; inz = nzn - nz;
	for (j=1; j<=nyn; j++)
		for (i=1; i<=lsdn; i++) 
			fout(i,j,1)=fint((i-1)/2*4+2-i%2,j*2-1,1);
	ret->set_complex(true);
	ret->set_ri(1);
	//ret->set_fftpad(true);
	//ret->set_attr("npad", 1);
	if (nxn%2 == 1) {ret->set_fftodd(true);} else {ret->set_fftodd(false);}
	if(RetReal) {
		ret->do_ift_inplace();
		ret->postift_depad_corner_inplace();
	}
	ret->update();
	
	delete temp_ft;
	temp_ft = 0;
	return ret;
}

EMData *EMData::Four_shuf_ds_cen_us(int nxn, int nyni, int nzni, bool RetReal) {

	int nyn, nzn, lsd, lsdn, inx, iny, inz;
	int i, j;

	nyn = nyni;
	nzn = 1;
	lsd = nx;
	lsdn = nxn + 2 - nxn%2;

	EMData *temp_ft = this->copy();
	EMData *ret = this->copy();
	ret->set_size(lsdn, nyn, nzn);
	ret->to_zero();
	float *fout = ret->get_data();
	float *fint = temp_ft->get_data();
//  TO KEEP EXACT VALUES ON THE ORIGINAL GRID ONE SHOULD USE
//  SQ2     = 2.0. HOWEVER, TOTAL ENERGY WILL NOT BE CONSERVED
	float  sq2 = 1.0f/std::sqrt(2.0f);

	for (i = 0; i < lsd*ny*nz; i++)  fint[i] *= 4;

	inx = nxn-(nx-2); iny = nyn - ny; inz = nzn - nz;
	for (j=1; j<=ny/4; j++)
		for (i=1; i<=(nx-2)/2+2; i++) {
			int g = (i-1)/2+1;
			if ((g+j)%2 == 0) {
				fout(i,j,1)=fint(g*4-2-i%2,j*2-1+ny/2,1);
			} else {
				fout(i,j,1)=-fint(g*4-2-i%2,j*2-1+ny/2,1);			
			}
		}
			
	for (j=ny/4+1; j<=ny/4+1; j++)
		for (i=1; i<=(nx-2)/2+2; i++) {
			int g = (i-1)/2+1;
			if ((g+j)%2 == 0) {
				fout(i,j,1)=fint(g*4-2-i%2,j*2-1-ny/2,1);
			} else {			
				fout(i,j,1)=-fint(g*4-2-i%2,j*2-1-ny/2,1);
			}			
		}
		
	for (j=ny/4+2; j<=ny/2; j++)
		for (i=1; i<=(nx-2)/2+2; i++) {
			int g = (i-1)/2+1;
			if ((g+j)%2 == 0) {
				fout(i,j+ny/2,1)=fint(g*4-2-i%2,j*2-1-ny/2,1);
			} else {			
				fout(i,j+ny/2,1)=-fint(g*4-2-i%2,j*2-1-ny/2,1);
			}			
		}
	
	if (nx%2 == 0) {
		for (j=1; j<=nyn; j++) {
			fout((nx-2)/2+1,j,1) *= sq2;
			fout((nx-2)/2+2,j,1) *= sq2;
		}
		for (i=1; i<=lsd/2+1; i++) {
			fout(i,ny/4+1+ny/2,1) = sq2*fout(i,ny/4+1,1);
			fout(i,ny/4+1,1) *= sq2;
		}
	}
	
	ret->set_complex(true);
	ret->set_ri(1);

	if (nxn%2 == 1) {ret->set_fftodd(true);} else {ret->set_fftodd(false);}
	if(RetReal) {
		ret->do_ift_inplace();
		ret->postift_depad_corner_inplace();
	}
	ret->update();
	
	delete temp_ft;
	temp_ft = 0;
	return ret;
}

#undef fint
#undef fout

#define  fint(jx,jy,jz)  fint[jx + (jy + jz*ny)*nox]
EMData *EMData::filter_by_image(EMData* image, bool RetReal) {


	bool   complex_input = this->is_complex();
	nx  = this->get_xsize();
	ny  = this->get_ysize();
	nz  = this->get_zsize();
	int nox;
	if (complex_input) nox = (nx - 2 + this->is_fftodd()); else nox = nx;

	int lsd2 = (nox + 2 - nox%2) / 2; // Extended x-dimension of the complex image

	EMData* fp = NULL; // output image
	if(complex_input) {
		// fimage must remain pristine
		fp = this->copy();
	} else {
		fp = this->pad_fft(1); 
		fp->do_fft_inplace();
	}
	fp->set_array_offsets(1,1,1);
	int nx2 = nox/2;
	int ny2 = ny/2;
	int nz2 = nz/2;
	float *fint = image->get_data();
	for ( int iz = 1; iz <= nz; iz++) {
		int jz=nz2-iz+1; if(jz<0) jz += nz;
		for ( int iy = 1; iy <= ny; iy++) {
			int jy=ny2-iy+1; if(jy<0) jy += ny;
			for ( int ix = 1; ix <= lsd2; ix++) {
				int jx = nx2-ix+1;
				fp->cmplx(ix,iy,iz) *= fint(jx,jy,jz);
			}
		}
	}

	fp->set_ri(1);
	fp->set_fftpad(true);
	fp->set_attr("npad", 1);
	if (nx%2 == 1) {fp->set_fftodd(true);} else {fp->set_fftodd(false);}
	if(RetReal) {
		fp->do_ift_inplace();
		fp->postift_depad_corner_inplace();
	}
	fp->set_array_offsets(0,0,0);
	fp->update();

	return fp;
}
#undef fint


namespace EMAN {
	/* #G2#
	Purpose: Create a new [normalized] [zero-padded] fft image. 
	Method: Normalize, pad with zeros, extend for fft, create new fft image,
	return new fft image.  
	Input: f real n-dimensional image
	flag specify normalize, pad, and/or extend 
	Output: fft image of normalized, zero-padded input image
	 */
	EMData* norm_pad_ft(EMData* fimage, bool donorm, bool dopad, int npad) {
		int nx = fimage->get_xsize();
		int ny = fimage->get_ysize();
		int nz = fimage->get_zsize();
		float mean = 0., stddev = 1.;
		if(donorm) { // Normalization requested
			mean = fimage->get_attr("mean");
			stddev = fimage->get_attr("sigma");
		}
		//  Padding requested
		if (dopad) {
			// sanity check
			if (npad < 2) npad = 2;
		} else {
			npad = 1;
		}
		EMData& fpimage = *fimage->pad_fft(npad);
		//  Perform the actual normalization (only on the
		//  non-zero section of the image)
		if(donorm) { // Normalization requested
			// indexing starts at 1
			vector<int> saved_offsets = fpimage.get_array_offsets();
			fpimage.set_array_offsets(1,1,1);
			for (int iz = 1; iz <= nz; iz++) for (int iy = 1; iy <= ny; iy++) for (int ix = 1; ix <= nx; ix++) fpimage(ix,iy,iz) = (fpimage(ix,iy,iz)-mean)/stddev;
			fpimage.set_array_offsets(saved_offsets); // reset
		}
		//mean = fpimage.get_attr("mean");  // I do not know why these were needed and I doubt they would work
		                               //  properly for fft extended image
		//stddev = fpimage.get_attr("sigma");
		fpimage.do_fft_inplace();
		return &fpimage;
	}

}

/* vim: set ts=4 noet: */
