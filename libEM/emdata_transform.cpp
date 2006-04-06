/**
 * $Id$
 */
#include "emdata.h"
#include "emfft.h"

using namespace EMAN;

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
	vector<int> saved_offsets = get_array_offsets();
	// iz in [1,nz], iy in [1,ny], ix in [0,nx/2]
	set_array_offsets(0,1,1);
	int xmax = (is_fftodd())
		? (nx-1)/2 + 1
		: (nx-2)/2;
	for (int iz = 1; iz <= nz; iz++) {
		for (int iy = 1; iy <= ny; iy++) {
			for (int ix = 0; ix <= xmax; ix++) {
				// next line multiplies by +/- 1
				cmplx(ix,iy,iz) *= static_cast<float>(-2*((ix+iy+iz)%2) + 1);
			}
		}
	}
	set_array_offsets(saved_offsets);
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
	int xstart = (nx != 1) ? (nxpad - nx)/2 + nx%2 : 0;
	int ystart = (ny != 1) ? (nypad - ny)/2 + ny%2 : 0;
	int zstart = (nz != 1) ? (nzpad - nz)/2 + nz%2 : 0;
	for (int iz = 0; iz < nz; iz++) {
		for (int iy = 0; iy < ny; iy++) {
			memcpy(&(*newimg)(xstart,iy+ystart,iz+zstart),
				   &(*this)(0,iy,iz), bytes);
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
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
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
		if (offset == 1) newimg->set_fftodd(true);
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
		for (int iz = 0; iz < nzold; iz++) {
			for (int iy = 0; iy < nyold; iy++) {
				memcpy(&(*newimg)(0,iy,iz), &(*this)(0,iy,iz), bytes);
			}
		}
	}
	newimg->done_data();
	set_array_offsets(saved_offsets);
	return newimg;
}


void EMData::postift_depad_corner_inplace() {
	ENTERFUNC;
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
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
	float* dest = get_data();
	for (int iz=0; iz < nzold; iz++) {
		for (int iy = 0; iy < nyold; iy++) {
			memmove(dest, &(*this)(0,iy,iz), bytes);
			dest += nxold;
		}
	}
	set_size(nxold, nyold, nzold);
	set_fftpad(false);
	update();
	set_complex(false);
	if(ny==1 && nz==1) {
		set_complex_x(false);
	}
	set_array_offsets(saved_offsets);
	EXITFUNC;
}


#define  fint(i,j,k)  fint[(i-1) + ((j-1) + (k-1)*ny)*lsd]
#define  fout(i,j,k)  fout[(i-1) + ((j-1) + (k-1)*nyn)*lsdn]
EMData *EMData::FourInterpol(int nxn, int nyni, int nzni) {

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
	lsd = nx+ 2 -nx%2;
	lsdn = nxn+ 2 -nxn%2;
//  do out of place ft
        EMData *temp_ft = do_fft();
	EMData *ret = new EMData();
	ret->set_size(lsdn, nyn, nzn);
	ret->to_zero();
	float *fout = ret->get_data();
	float *fint = temp_ft->get_data();
//  TO KEEP THE EXACT VALUES ON THE PREVIOUS GRID ONE SHOULD USE
//  SQ2     = 2.0. HOWEVER, TOTAL ENERGY WILL NOT BE CONSERVED
	float  sq2 = 1.0f/sqrt(2.0f);
	float  anorm = (float) nxn* (float) nyn* (float) nzn/(float) nx/ (float) ny/ (float) nz;
	//for (i = 0; i < lsd*ny*nz; i++)  fout[i] = fint[i];
	for (i = 0; i < lsd*ny*nz; i++)  fint[i] *= anorm;
	inx = nxn-nx; iny = nyn - ny; inz = nzn - nz;
	for (k=1; k<=nz/2+1; k++) {
	  for (j=1; j<=ny/2+1; j++) {
	    for (i=1; i<=lsd; i++) {
	      fout(i,j,k)=fint(i,j,k);
	    }
	  }
	}
	if(nyn>1) {
	//cout << "  " <<nxn<<"  " <<nyn<<" A " <<nzn<<endl;
	 for (k=1; k<=nz/2+1; k++) {
	   for (j=ny/2+2+iny; j<=nyn; j++) {
	     for (i=1; i<=lsd; i++) {
	       fout(i,j,k)=fint(i,j-iny,k);
	     }
	   }
	 }
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
	ret->set_ri(1);
	ret->set_fftpad(true);
	ret->set_attr("npad", 1);
	if (nxn%2 == 1) {ret->set_fftodd(true);}else{ret->set_fftodd(false);}
	ret->do_ift_inplace();
	ret->postift_depad_corner_inplace();
	ret->done_data();
	
	/*Dict d1 = temp_ft->get_attr_dict();
	Dict d2 = ret->get_attr_dict();
	printf("-----------------Attribute Dict for temp_ft--------------\n");
	EMUtil::dump_dict(d1);
	printf("-----------------Attribute Dict for ret--------------\n");
	EMUtil::dump_dict(d2);*/
	
	return ret;
}
#undef fint
#undef fout


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
	if (offset == 1)  dat->set_fftodd(true);

	float *d = dat->get_data();
	EMfft::real_to_complex_nd(rdata, d, nxreal, ny, nz);

	dat->done_data();
	dat->set_complex(true);
	if(dat->get_ysize()==1 && dat->get_zsize()==1) {
		dat->set_complex_x(true);
	}
	dat->set_ri(true);


	done_data();

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
		offset = 2 - nx%2;
		if (1 == offset) set_fftodd(true);
		int nxnew = nx + offset;
		set_size(nxnew, ny, nz);
		for (int iz = nz-1; iz >= 0; iz--) {
			for (int iy = ny-1; iy >= 0; iy--) {
				for (int ix = nxreal-1; ix >= 0; ix--) {
					size_t oldxpos = ix + (iy + iz*ny)*nxreal;
					size_t newxpos = ix + (iy + iz*ny)*nx;
					(*this)(newxpos) = (*this)(oldxpos);
				}
			}
		}
		/*// zero out padding   SHOULD NOT BE NECCESSARY PAP 01/28/06
		for (int iz=0; iz < nz; iz++)
			for (int iy=0; iy < ny; iy++)
				for (int ix=nxreal; ix < nx; ix++)
					(*this)(ix,iy,iz) = 0.f;*/
		set_fftpad(true);
	} else {
		offset = is_fftodd() ? 1 : 2;
		nxreal = nx - offset;
	}
	EMfft::real_to_complex_nd(rdata, rdata, nxreal, ny, nz);

	set_complex(true);
	if(ny==1 && nz==1) {
		set_complex_x(true);
	}
	set_ri(true);

	done_data();

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

	if (ndim >= 2) {
		EMfft::complex_to_real_nd(d, d, nx - offset, ny, nz);

		size_t row_size = (nx - offset) * sizeof(float);
		for (int i = 1; i < ny * nz; i++) {
			memmove((char *) &d[i * (nx - offset)], (char *) &d[i * nx], row_size);
		}
	}

#if defined	FFTW2 || defined FFTW3	//native fft and ACML already done normalization
	// SCALE the inverse FFT
	float scale = 1.0f / ((nx - offset) * ny * nz);
	dat->mult(scale);
#endif	//NATIVE_FFT || ACML

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


	done_data();

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

	int offset = is_fftodd() ? 1 : 2;
	EMfft::complex_to_real_nd(rdata, rdata, nx - offset, ny, nz);

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


	done_data();

	EXITFUNC;
	return this;
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
	update();
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
	update();
	EXITFUNC;
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

	update();
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

