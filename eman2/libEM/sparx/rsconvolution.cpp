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

#include <algorithm>

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

using namespace EMAN;
using std::swap;

namespace {
	// K(i,j,k)*f(a-i, b-j, c-k) <-- internal, so no boundary condition issues
	inline float mult_internal(EMData& K, EMData& f, 
						       int kzmin, int kzmax, int kymin, int kymax, 
							   int kxmin, int kxmax, 
							   int iz, int iy, int ix) {
		float sum = 0.f;
		for (int kz = kzmin; kz <= kzmax; kz++) {
			for (int ky = kymin; ky <= kymax; ky++) {
				for (int kx = kxmin; kx <= kxmax; kx++) {
					float Kp = K(kx,ky,kz);
					float fp = f(ix-kx,iy-ky,iz-kz);
					sum += Kp*fp;
				}
			}
		}
		return sum;
	}
	// K(i,j,k)*f(a-i, b-j, c-k) <-- Circulant boundary conditions
	inline float mult_circ(EMData& K, EMData& f, int kzmin, int kzmax,
			               int kymin, int kymax, int kxmin, int kxmax,
						   int nzf, int nyf, int nxf, int iz, int iy, int ix) {
		float sum = 0.f;
		for (int kz = kzmin; kz <= kzmax; kz++) {
			int jz = (iz - kz) % nzf;
			if (jz < 0) jz += nzf;
			for (int ky = kymin; ky <= kymax; ky++) {
				int jy = (iy - ky) % nyf; 
				if (jy < 0) jy += nyf;
				for (int kx = kxmin; kx <= kxmax; kx++) {
					int jx = (ix - kx) % nxf; 
					if (jx < 0) jx += nxf;
					float Kp = K(kx,ky,kz);
					float fp = f(jx,jy,jz);
					sum += Kp*fp;
				}
			}
		}
		return sum;
	}
	// In the future we may want to add other boundary conditions here
	

	// This function selects the k-th smallest element in the array table[0..n-1]
	// i.e. k = 1 gets the smallest element
	//      k = n gets the largest element
	float select_kth_smallest(float *table, int n, int k) {

		int i,j,left,middle,right;
		float temp;
		bool flag = 0;

		left = 0;
		right = n-1;
		k = k-1;  // This is necessary because the index of array begins from 0
		while (flag == 0) {
			if ( left+1 < right ) {
				middle = (left+right)/2;
				swap(table[middle],table[left+1]);
				if ( table[left+1] > table [right] )
					swap(table[left+1], table[right]);
				if ( table[left] > table[right] )
					swap(table[left], table[right]);
				if ( table[left+1] > table[left] )
					swap(table[left+1], table[left]);
				i = left+1;
				j = right;
				temp = table[left];
				do {
					i++; 
					while (table[i] < temp) i++;
					j--;
					while (table[j] > temp) j--;
					if (j >= i) 
						swap(table[i], table[j]);
				} while (j >= i);
				table[left] = table[j];
				table[j] = temp;
				if (j >= k) right = j-1;
				if (j <= k) left = i;
			} else {
				if ( right == left+1 && table[right] < table[left] ) 
					swap(table[left], table[right]);
				flag = 1;
			}
		}
		return table[k];
	}
	
	inline float median(EMData& f, int nxk, int nyk, int nzk, kernel_shape myshape, int iz, int iy, int ix) {
		size_t index = 0;
		int dimension = 3;
		float median_value = 0.f;
		float *table = 0;

		int nxf = (&f)->get_xsize();
		int nyf = (&f)->get_ysize();
		int nzf = (&f)->get_zsize();

		int nxk2 = (nxk-1)/2;
		int nyk2 = (nyk-1)/2;
		int nzk2 = (nzk-1)/2;

		int kzmin = iz-nzk2;
		int kzmax = iz+nzk2;
		int kymin = iy-nyk2;
		int kymax = iy+nyk2;
		int kxmin = ix-nxk2;
		int kxmax = ix+nxk2;

		if ( nzf == 1 ) {
			dimension--;
			if ( nyf == 1 )  dimension--; 
		}

		switch (myshape) {
		case BLOCK:
			switch (dimension) {
			case 1: 
				table = (float*)malloc(nxk*sizeof(float));
				break;
			case 2: table = (float*)malloc(nxk*nyk*sizeof(float));
				break;
			case 3: table = (float*)malloc(nxk*nyk*nzk*sizeof(float));
			 	break;
			}	
			for (int kz = kzmin; kz <= kzmax; kz++) {
				int jz = kz < 0 ? kz+nzf : kz % nzf;
				for (int ky = kymin; ky <= kymax; ky++) {
					int jy = ky < 0 ? ky+nyf : ky % nyf; 
					for (int kx = kxmin; kx <= kxmax; kx++) {
						int jx = kx < 0 ? kx+nxf : kx % nxf; 
						table[index] = f(jx,jy,jz);
						index++;
					}
				}
			}
			break;
		case CIRCULAR:
			switch (dimension) {
			case 1: 
				table = (float*)malloc(nxk*sizeof(float));
				break;
			case 2: table = (float*)malloc(nxk*nxk*sizeof(float));
				break;
			case 3: table = (float*)malloc(nxk*nxk*nxk*sizeof(float));
			 	break;
			}	
			for (int kz = kzmin; kz <= kzmax; kz++) {
				int jz = kz < 0 ? kz+nzf : kz % nzf;
				for (int ky = kymin; ky <= kymax; ky++) {
					int jy = ky < 0 ? ky+nyf : ky % nyf; 
					for (int kx = kxmin; kx <= kxmax; kx++) {
						int jx = kx < 0 ? kx+nxf : kx % nxf; 
						if ( (kz-iz)*(kz-iz)+(ky-iy)*(ky-iy)+(kx-ix)*(kx-ix) <= nxk2*nxk2 ) {
							table[index] = f(jx,jy,jz);
							index++;
						}
					}
				}
			}
			break;
		case CROSS:
			if ( nzf != 1 )  {
				table = (float*)malloc((nxk+nyk+nzk-2)*sizeof(float));
				for (int kz = kzmin; kz <= kzmax; kz++) {
					int jz = kz < 0 ? kz+nzf : kz % nzf;
					if ( kz != iz ) { table[index] = f(ix,iy,jz); index++; }
				}
				for (int ky = kymin; ky <= kymax; ky++) {
					int jy = ky < 0 ? ky+nyf : ky % nyf; 
					if ( ky != iy ) { table[index] = f(ix,jy,iz); index++; }
				}
				for (int kx = kxmin; kx <= kxmax; kx++) {
					int jx = kx < 0 ? kx+nxf : kx % nxf; 
					table[index] = f(jx,iy,iz);
					index++;
				}
			} else if  ( nyf != 1 ) {
				table = (float*)malloc((nxk+nyk-1)*sizeof(float));
				for (int ky = kymin; ky <= kymax; ky++) {
					int jy = ky < 0 ? ky+nyf : ky % nyf; 
					if ( ky != iy ) { table[index] = f(ix,jy,iz); index++; }
				}
				for (int kx = kxmin; kx <= kxmax; kx++) {
					int jx = kx < 0 ? kx+nxf : kx % nxf; 
					table[index] = f(jx,iy,iz);
					index++;
				}
			} else {
				table = (float*)malloc(nxk*sizeof(float));
				for (int kx = kxmin; kx <= kxmax; kx++) {
					int jx = kx < 0 ? kx+nxf : kx % nxf; 
					table[index] = f(jx,iy,iz);
					index++;
				}
			}
			break;
		default: throw ImageDimensionException("Illegal Kernal Shape!");
		}
		median_value=select_kth_smallest(table, index, (index+1)/2);
		free((void *)table);
		return median_value;
	}
}

namespace EMAN {

    EMData* rsconvolution(EMData* f, EMData* K) {//Does not work properly in 3D, corners are not done, PAP 07/16/09
		// Kernel should be the smaller image
		int nxf=f->get_xsize(); int nyf=f->get_ysize(); int nzf=f->get_zsize();
		int nxK=K->get_xsize(); int nyK=K->get_ysize(); int nzK=K->get_zsize();
		if ((nxf<nxK)&&(nyf<nyK)&&(nzf<nzK)) {
			// whoops, f smaller than K
			swap(f,K); swap(nxf,nxK); swap(nyf,nyK); swap(nzf,nzK);
		} else if ((nxK<=nxf)&&(nyK<=nyf)&&(nzK<=nzf)) {
			// that's what it should be, so do nothing
			;
		} else {
			// incommensurate sizes
			throw ImageDimensionException("input images are incommensurate");
		}
		// Kernel needs to be _odd_ in size
		if ((nxK % 2 != 1) || (nyK % 2 != 1) || (nzK % 2 != 1))
			throw ImageDimensionException("Real-space convolution kernel"
				" must have odd nx,ny,nz (so the center is well-defined).");
		EMData* result = new EMData();
		result->set_size(nxf, nyf, nzf);
		result->to_zero();
		// kernel corners, need to check for degenerate case
		int kxmin = -nxK/2; int kymin = -nyK/2; int kzmin = -nzK/2;
		int kxmax = (1 == nxK % 2) ? -kxmin : -kxmin - 1;
		int kymax = (1 == nyK % 2) ? -kymin : -kymin - 1;
		int kzmax = (1 == nzK % 2) ? -kzmin : -kzmin - 1;
		vector<int> K_saved_offsets = K->get_array_offsets();
		K->set_array_offsets(kxmin,kymin,kzmin);
		// interior boundaries, need to check for degenerate cases
		int izmin = 0, izmax = 0, iymin = 0, iymax = 0, ixmin = 0, ixmax = 0;
		if (1 != nzf) {
			izmin = -kzmin;
			izmax = nzf - 1 - kzmax;
		}
		if (1 != nyf) {
			iymin = -kymin;
			iymax = nyf - 1 - kymax;
		}
		if (1 != nxf) {
			ixmin = -kxmin;
			ixmax = nxf - 1 - kxmax;
		}
		// interior (no boundary condition issues here)
		for (int iz = izmin; iz <= izmax; iz++) {
			for (int iy = iymin; iy <= iymax; iy++) {
				for (int ix = ixmin; ix <= ixmax; ix++) {
					(*result)(ix,iy,iz) =
						mult_internal(*K, *f, 
								      kzmin, kzmax, kymin, kymax, kxmin, kxmax,
									  iz, iy, ix);
				}
			}
		}
		// corners
		// corner sizes, with checking for degenerate cases
		int sz = (1 == nzK) ? 1 : -kzmin + kzmax;
		int sy = (1 == nyK) ? 1 : -kymin + kymax;
		int sx = (1 == nxK) ? 1 : -kxmin + kxmax;
		// corner starting locations, with checking for degenerate cases
		int zstart = (0 == izmin) ? 0 : izmin - 1;
		int ystart = (0 == iymin) ? 0 : iymin - 1;
		int xstart = (0 == ixmin) ? 0 : ixmin - 1;
		// corners
		for (int cz = 0; cz < sz; cz++) {
			int iz = (zstart - cz) % nzf;
			if (iz < 0) iz += nzf;
			for (int cy = 0; cy < sy; cy++) {
				int iy = (ystart - cy) % nyf;
				if (iy < 0) iy += nyf;
				for (int cx=0; cx < sx; cx++) {
					int ix = (xstart - cx) % nxf;
					if (ix < 0) ix += nxf;
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, 
								 kymax, kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}
		// remaining stripes -- should use a more elegant (non-3D-specific) method here
		// ix < ixmin
		for (int ix = 0; ix < ixmin; ix++) {
			for (int iy = iymin; iy <= iymax; iy++) {
				for (int iz = izmin; iz <= izmax; iz++) {
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, kymax, 
								 kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}
		// ix > ixmax
		for (int ix = ixmax+1; ix < nxf; ix++) {
			for (int iy = iymin; iy <= iymax; iy++) {
				for (int iz = izmin; iz <= izmax; iz++) {
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, kymax, 
								 kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}
		// iy < iymin
		for (int iy = 0; iy < iymin; iy++) {
			for (int ix = ixmin; ix <= ixmax; ix++) {
				for (int iz = izmin; iz <= izmax; iz++) {
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, kymax, 
								 kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}
		// iy > iymax
		for (int iy = iymax+1; iy < nyf; iy++) {
			for (int ix = ixmin; ix <= ixmax; ix++) {
				for (int iz = izmin; iz <= izmax; iz++) {
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, kymax, 
								 kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}
		// iz < izmin
		for (int iz = 0; iz < izmin; iz++) {
			for (int ix = ixmin; ix <= ixmax; ix++) {
				for (int iy = iymin; iy <= iymax; iy++) {
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, kymax, 
								 kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}
		// iz > izmax
		for (int iz = izmax+1; iz < nzf; iz++) {
			for (int ix = ixmin; ix <= ixmax; ix++) {
				for (int iy = iymin; iy <= iymax; iy++) {
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, kymax, 
								 kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}
		
		
		// ix < ixmin, iy < iymin
		for (int ix = 0; ix < ixmin; ix++) {
			for (int iy = 0; iy < iymin; iy++) {
				for (int iz = izmin; iz <= izmax; iz++) {
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, kymax, 
								 kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}
		
		// ix < ixmin, iy > iymax
		for (int ix = 0; ix < ixmin; ix++) {
			for (int iy = iymax+1; iy < nyf; iy++) {
				for (int iz = izmin; iz <= izmax; iz++) {
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, kymax, 
								 kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}

        // ix > ixmax, iy < iymin
		for (int ix = ixmax+1; ix < nxf; ix++) {
			for (int iy = 0; iy < iymin; iy++) {
				for (int iz = izmin; iz <= izmax; iz++) {
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, kymax, 
								 kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}
		
		// ix > ixmax, iy > iymax
		for (int ix = ixmax+1; ix < nxf; ix++) {
			for (int iy = iymax+1; iy < nyf; iy++) {
				for (int iz = izmin; iz <= izmax; iz++) {
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, kymax, 
								 kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}


		
        // ix < ixmin, iz < izmin
		for (int ix = 0; ix < ixmin; ix++) {
			for (int iy = iymin; iy <= iymax; iy++) {
				for (int iz = 0; iz < izmin; iz++) {
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, kymax, 
								 kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}
		
		 // ix < ixmin, iz > izmax
		for (int ix = 0; ix < ixmin; ix++) {
			for (int iy = iymin; iy <= iymax; iy++) {
				for (int iz = izmax+1; iz < nzf; iz++) {
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, kymax, 
								 kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}


         // ix > ixmin, iz < izmin
		for (int ix = ixmax+1; ix < nxf; ix++) {
			for (int iy = iymin; iy <= iymax; iy++) {
				for (int iz = 0; iz < izmin; iz++) {
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, kymax, 
								 kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}
		
		 // ix > ixmin, iz > izmax
		for (int ix = ixmax+1; ix < nxf; ix++) {
			for (int iy = iymin; iy <= iymax; iy++) {
				for (int iz = izmax+1; iz < nzf; iz++) {
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, kymax, 
								 kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}

		

       // iy < iymin, iz < izmin
	   
	   for (int iz = 0; iz < izmin; iz++) {
			for (int ix = ixmin; ix <= ixmax; ix++) {
				for (int iy = 0; iy < iymin; iy++) {
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, kymax, 
								 kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}


       // iy < iymin, iz > izmax
	   
	   for (int iz = izmax+1; iz < nzf; iz++) {
			for (int ix = ixmin; ix <= ixmax; ix++) {
				for (int iy = 0; iy < iymin; iy++) {
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, kymax, 
								 kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}
		
		
		// iy > iymax, iz < izmin
	   
	   for (int iz = 0; iz < izmin; iz++) {
			for (int ix = ixmin; ix <= ixmax; ix++) {
				for (int iy = iymax+1; iy < nyf; iy++) {
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, kymax, 
								 kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}


       // iy > iymax, iz > izmax
	   
	   for (int iz = izmax+1; iz < nzf; iz++) {
			for (int ix = ixmin; ix <= ixmax; ix++) {
				for (int iy = iymax+1; iy < nyf; iy++) {
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, kymax, 
								 kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}

		
		K->set_array_offsets(K_saved_offsets);
		result->update();
		return result;
	}

    EMData* filt_median_(EMData* f, int nxk, int nyk, int nzk, kernel_shape myshape) {
		
 		int nxf = f->get_xsize();
		int nyf = f->get_ysize(); 
		int nzf = f->get_zsize();
		
		if ( nxk > nxf || nyk > nyf || nzk > nzf ) {
			// Kernel should be smaller than the size of image
			throw ImageDimensionException("Kernel should be smaller than the size of image.");
		}	

		if ( nxk % 2 != 1 || nyk % 2 != 1 || nzk % 2 != 1 ) {
			// Kernel needs to be odd in size
			throw ImageDimensionException("Real-space kernel must have odd size so that the center is well-defined.");
		}

		if ( myshape == CIRCULAR ) {
			// For CIRCULAR kernal, size must be same on all dimensions
			if ( (nzf != 1 && ( nxk != nyk || nxk != nzk )) || (nzf == 1 && nyf != 1 && nxk != nyk) ) {
				throw ImageDimensionException("For CIRCULAR kernal, size must be same on all dimensions.");
			}
		}

		EMData* result = new EMData();
		result->set_size(nxf, nyf, nzf);
		result->to_zero();

		for (int iz = 0; iz <= nzf-1; iz++) {
			for (int iy = 0; iy <= nyf-1; iy++) {
				for (int ix = 0; ix <= nxf-1; ix++) {
					(*result)(ix,iy,iz) = median (*f, nxk, nyk, nzk, myshape, iz, iy, ix);					
				}
			}
		}
		
		return result;
	}

    EMData* filt_dilation_(EMData* f, EMData* K, morph_type mydilation) {

 		int nxf = f->get_xsize();
		int nyf = f->get_ysize(); 
		int nzf = f->get_zsize();

		int nxk = K->get_xsize();
		int nyk = K->get_ysize();
		int nzk = K->get_zsize();
	
		if ( nxf < nxk && nyf < nyk && nzf < nzk ) {
			// whoops, f smaller than K
			swap(f,K); swap(nxf,nxk); swap(nyf,nyk); swap(nzf,nzk);
		} else if ( nxk > nxf || nyk > nyf || nzk > nzf ) {
			// Incommensurate sizes
			throw ImageDimensionException("Two input images are incommensurate.");
		}

		if ( nxk % 2 != 1 || nyk % 2 != 1 || nzk % 2 != 1 ) {
			// Kernel needs to be odd in size
			throw ImageDimensionException("Kernel should have odd nx,ny,nz so that the center is well-defined.");
		}

		int nxk2 = (nxk-1)/2;
		int nyk2 = (nyk-1)/2;
		int nzk2 = (nzk-1)/2;

		if ( mydilation == BINARY ) {
			// Check whether two images are truly binary.
	 		for (int iz = 0; iz <= nzf-1; iz++) {
				for (int iy = 0; iy <= nyf-1; iy++) {
					for (int ix = 0; ix <= nxf-1; ix++) {
						int fxyz=(int)(*f)(ix,iy,iz);
						if ( fxyz != 0 && fxyz != 1 ) {
							throw ImageDimensionException("One of the two images is not binary.");
						}
					}
				}
			}
	 		for (int iz = 0; iz <= nzk-1; iz++) {
				for (int iy = 0; iy <= nyk-1; iy++) {
					for (int ix = 0; ix <= nxk-1; ix++) {
						int kxyz=(int)(*K)(ix,iy,iz);
						if ( kxyz != 0 && kxyz != 1 ) {
							throw ImageDimensionException("One of the two images is not binary.");
						}
					}
				}
			}
		}

		EMData* result = new EMData();
		result->set_size(nxf, nyf, nzf);
		result->to_zero();

		for (int iz = 0; iz <= nzf-1; iz++) {
			for (int iy = 0; iy <= nyf-1; iy++) {
				for (int ix = 0; ix <= nxf-1; ix++) {
//					int kzmin = iz-nzk2 < 0     ?   0   : iz-nzk2 ;
//					int kzmax = iz+nzk2 > nzf-1 ? nzf-1 : iz+nzk2 ;
//					int kymin = iy-nyk2 < 0     ?   0   : iy-nyk2 ;
//					int kymax = iy+nyk2 > nyf-1 ? nyf-1 : iy+nyk2 ;
//					int kxmin = ix-nxk2 < 0     ?   0   : ix-nxk2 ;
//					int kxmax = ix+nxk2 > nxf-1 ? nxf-1 : ix+nxk2 ;
					if ( mydilation == BINARY ) {
						int fxyz = (int)(*f)(ix,iy,iz);
						if ( fxyz == 1 ) {
							for (int jz = -nzk2; jz <= nzk2; jz++) {
								for (int jy = -nyk2; jy <= nyk2; jy++) {
									for (int jx= -nxk2; jx <= nxk2; jx++) {
										if ( (int)(*K)(jx+nxk2,jy+nyk2,jz+nzk2) == 1 ) {
											int fz = iz+jz;
											int fy = iy+jy;
											int fx = ix+jx;
											if ( fz >= 0 && fz <= nzf-1 && fy >= 0 && fy <= nyf-1 && fx >= 0 && fx <= nxf-1 )
												(*result)(fx,fy,fz) = 1;
											}
										}
									}
								}
							}
					} else if ( mydilation == GRAYLEVEL ) {
							float pmax = (*f)(ix,iy,iz)+(*K)(nxk2,nyk2,nzk2); 
							for (int jz = -nzk2; jz <= nzk2; jz++) {
								for (int jy = -nyk2; jy <= nyk2; jy++) {
									for (int jx = -nxk2; jx <= nxk2; jx++) {
										int fz = iz+jz;
										int fy = iy+jy;
										int fx = ix+jx;
										if ( fz >= 0 && fz <= nzf-1 && fy >= 0 && fy <= nyf-1 && fx >= 0 && fx <= nxf-1 ) {
											float kxyz = (*K)(jx+nxk2,jy+nyk2,jz+nzk2);
											float fxyz = (*f)(fx,fy,fz);											
											if ( kxyz+fxyz > pmax )  pmax = kxyz+fxyz;
										}
									}
								}
							}
							(*result)(ix,iy,iz) = pmax;
					} else {
						throw ImageDimensionException("Illegal dilation type!");
					}
				}
			}
		}		
		return result;
    }

    EMData* filt_erosion_(EMData* f, EMData* K, morph_type myerosion) {

 		int nxf = f->get_xsize();
		int nyf = f->get_ysize(); 
		int nzf = f->get_zsize();

		int nxk = K->get_xsize();
		int nyk = K->get_ysize();
		int nzk = K->get_zsize();
	
		if ( nxf < nxk && nyf < nyk && nzf < nzk ) {
			// whoops, f smaller than K
			swap(f,K); swap(nxf,nxk); swap(nyf,nyk); swap(nzf,nzk);
		} else if ( nxk > nxf || nyk > nyf || nzk > nzf ) {
			// Incommensurate sizes
			throw ImageDimensionException("Two input images are incommensurate.");
		}

		if ( nxk % 2 != 1 || nyk % 2 != 1 || nzk % 2 != 1 ) {
			// Kernel needs to be odd in size
			throw ImageDimensionException("Kernel should have odd nx,ny,nz so that the center is well-defined.");
		}

		int nxk2 = (nxk-1)/2;
		int nyk2 = (nyk-1)/2;
		int nzk2 = (nzk-1)/2;

		if ( myerosion == BINARY ) {
			// Check whether two images are truly binary.
	 		for (int iz = 0; iz <= nzf-1; iz++) {
				for (int iy = 0; iy <= nyf-1; iy++) {
					for (int ix = 0; ix <= nxf-1; ix++) {
						int fxyz=(int)(*f)(ix,iy,iz);
						if ( fxyz != 0 && fxyz != 1 ) {
							throw ImageDimensionException("One of the two images is not binary.");
						}
					}
				}
			}
	 		for (int iz = 0; iz <= nzk-1; iz++) {
				for (int iy = 0; iy <= nyk-1; iy++) {
					for (int ix = 0; ix <= nxk-1; ix++) {
						int kxyz=(int)(*K)(ix,iy,iz);
						if ( kxyz != 0 && kxyz != 1 ) {
							throw ImageDimensionException("One of the two images is not binary.");
						}
					}
				}
			}
		}

		EMData* result = new EMData();
		result->set_size(nxf, nyf, nzf);
		result->to_one();

		for (int iz = 0; iz <= nzf-1; iz++) {
			for (int iy = 0; iy <= nyf-1; iy++) {
				for (int ix = 0; ix <= nxf-1; ix++) {
					if ( myerosion == BINARY ) {
						int fxyz = (int)(*f)(ix,iy,iz);
						if ( fxyz == 0 ) {
							for (int jz = -nzk2; jz <= nzk2; jz++) {
								for (int jy = -nyk2; jy <= nyk2; jy++) {
									for (int jx= -nxk2; jx <= nxk2; jx++) {
										if ( (int)(*K)(jx+nxk2,jy+nyk2,jz+nzk2) == 1 ) {
											int fz = iz+jz;
											int fy = iy+jy;
											int fx = ix+jx;
											if ( fz >= 0 && fz <= nzf-1 && fy >= 0 && fy <= nyf-1 && fx >= 0 && fx <= nxf-1 )
												(*result)(fx,fy,fz) = 0;
											}
										}
									}
								}
							}
					} else if ( myerosion == GRAYLEVEL ) {
							float pmin = (*f)(ix,iy,iz)-(*K)(nxk2,nyk2,nzk2); 
							for (int jz = -nzk2; jz <= nzk2; jz++) {
								for (int jy = -nyk2; jy <= nyk2; jy++) {
									for (int jx = -nxk2; jx <= nxk2; jx++) {
										int fz = iz+jz;
										int fy = iy+jy;
										int fx = ix+jx;
										if ( fz >= 0 && fz <= nzf-1 && fy >= 0 && fy <= nyf-1 && fx >= 0 && fx <= nxf-1 ) {
											float kxyz = (*K)(jx+nxk2,jy+nyk2,jz+nzk2);
											float fxyz = (*f)(fx,fy,fz);											
											if ( fxyz-kxyz < pmin )  pmin = fxyz-kxyz;
										}
									}
								}
							}
							(*result)(ix,iy,iz) = pmin;
					} else {
						throw ImageDimensionException("Illegal dilation type!");
					}
				}
			}
		}		
		return result;
    }

}

/*
namespace {
	// K(i,j,k)*f(a-i, b-j, c-k) <-- internal, so no boundary condition issues
	inline float kmlt_internal(EMData& K, EMData& f, 
						       int kzmin, int kzmax, int kymin, int kymax, 
							   int kxmin, int kxmax, 
							   int iz, int iy, int ix) {
		float sum = 0.f;
		for (int kz = kzmin; kz <= kzmax; kz++) {
			for (int ky = kymin; ky <= kymax; ky++) {
				for (int kx = kxmin; kx <= kxmax; kx++) {
					float Kp = K(kx,ky,kz);
					float fp = f(ix-kx,iy-ky,iz-kz);
					sum += Kp*fp;
				}
			}
		}
		return sum;
	}
	// K(i,j,k)*f(a-i, b-j, c-k) <-- Circulant boundary conditions
	inline float kmlt_circ(EMData& K, EMData& f, int kzmin, int kzmax,
			               int kymin, int kymax, int kxmin, int kxmax,
						   int nzf, int nyf, int nxf, int iz, int iy, int ix) {
		float sum = 0.f;
		for (int kz = kzmin; kz <= kzmax; kz++) {
			int jz = (iz - kz) % nzf;
			if (jz < 0) jz += nzf;
			for (int ky = kymin; ky <= kymax; ky++) {
				int jy = (iy - ky) % nyf; 
				if (jy < 0) jy += nyf;
				for (int kx = kxmin; kx <= kxmax; kx++) {
					int jx = (ix - kx) % nxf; 
					if (jx < 0) jx += nxf;
					float Kp = K(kx,ky,kz);
					float fp = f(jx,jy,jz);
					sum += Kp*fp;
				}
			}
		}
		return sum;
	}
	// In the future we may want to add other boundary conditions here
}
namespace EMAN {

    EMData* rscp(EMData* f) {
		// Kernel should be the smaller image
		int nxf=f->get_xsize(); int nyf=f->get_ysize(); int nzf=f->get_zsize();
		const int npad = 2;
		const int m = Util::get_min(nxf,nyf,nzf);
		const int n = m*npad;

		const int K = 6;  //params["kb_K"];
		const float alpha = 1.75;  //params["kb_alpha"];
		Util::KaiserBessel kb(alpha, K, m/2,K/(2.*n),n);

                int nxK = K/2+1; nyK=nxK; nzK=nxK;

		EMData* result = new EMData();
		result->set_size(nxf, nyf, nzf);
		result->to_zero();
		// kernel corners, need to check for degenerate case
		int kxmin = -nxK/2; int kymin = -nyK/2; int kzmin = -nzK/2;
		int kxmax = (1 == nxK % 2) ? -kxmin : -kxmin - 1;
		int kymax = (1 == nyK % 2) ? -kymin : -kymin - 1;
		int kzmax = (1 == nzK % 2) ? -kzmin : -kzmin - 1;
		// interior boundaries, need to check for degenerate cases
		int izmin = 0, izmax = 0, iymin = 0, iymax = 0, ixmin = 0, ixmax = 0;
		if (1 != nzf) {
			izmin = -kzmin;
			izmax = nzf - 1 - kzmax;
		}
		if (1 != nyf) {
			iymin = -kymin;
			iymax = nyf - 1 - kymax;
		}
		if (1 != nxf) {
			ixmin = -kxmin;
			ixmax = nxf - 1 - kxmax;
		}
		// interior (no boundary condition issues here)
		for (int iz = izmin; iz <= izmax; iz++) {
			for (int iy = iymin; iy <= iymax; iy++) {
				for (int ix = ixmin; ix <= ixmax; ix++) {
					(*result)(ix,iy,iz) =
						mult_internal(*K, *f, 
								      kzmin, kzmax, kymin, kymax, kxmin, kxmax,
									  iz, iy, ix);
				}
			}
		}
		//   INITIALLY SKIP IT / corners  
		// corner sizes, with checking for degenerate cases
		int sz = (1 == nzK) ? 1 : -kzmin + kzmax;
		int sy = (1 == nyK) ? 1 : -kymin + kymax;
		int sx = (1 == nxK) ? 1 : -kxmin + kxmax;
		// corner starting locations, with checking for degenerate cases
		int zstart = (0 == izmin) ? 0 : izmin - 1;
		int ystart = (0 == iymin) ? 0 : iymin - 1;
		int xstart = (0 == ixmin) ? 0 : ixmin - 1;
		// corners
		for (int cz = 0; cz < sz; cz++) {
			int iz = (zstart - cz) % nzf;
			if (iz < 0) iz += nzf;
			for (int cy = 0; cy < sy; cy++) {
				int iy = (ystart - cy) % nyf;
				if (iy < 0) iy += nyf;
				for (int cx=0; cx < sx; cx++) {
					int ix = (xstart - cx) % nxf;
					if (ix < 0) ix += nxf;
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, 
								 kymax, kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}
		// remaining stripes -- should use a more elegant (non-3D-specific) method here
		// ix < ixmin
		for (int ix = 0; ix < ixmin; ix++) {
			for (int iy = iymin; iy <= iymax; iy++) {
				for (int iz = izmin; iz <= izmax; iz++) {
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, kymax, 
								 kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}
		// ix > ixmax
		for (int ix = ixmax+1; ix < nxf; ix++) {
			for (int iy = iymin; iy <= iymax; iy++) {
				for (int iz = izmin; iz <= izmax; iz++) {
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, kymax, 
								 kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}
		// iy < iymin
		for (int iy = 0; iy < iymin; iy++) {
			for (int ix = ixmin; ix <= ixmax; ix++) {
				for (int iz = izmin; iz <= izmax; iz++) {
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, kymax, 
								 kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}
		// iy > iymax
		for (int iy = iymax+1; iy < nyf; iy++) {
			for (int ix = ixmin; ix <= ixmax; ix++) {
				for (int iz = izmin; iz <= izmax; iz++) {
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, kymax, 
								 kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}
		// iz < izmin
		for (int iz = 0; iz < izmin; iz++) {
			for (int ix = ixmin; ix <= ixmax; ix++) {
				for (int iy = iymin; iy <= iymax; iy++) {
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, kymax, 
								 kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}
		// iz > izmax
		for (int iz = izmax+1; iz < nzf; iz++) {
			for (int ix = ixmin; ix <= ixmax; ix++) {
				for (int iy = iymin; iy <= iymax; iy++) {
					(*result)(ix,iy,iz) =
						mult_circ(*K, *f, kzmin, kzmax, kymin, kymax, 
								 kxmin, kxmax,
								 nzf, nyf, nxf, iz, iy, ix);
				}
			}
		}
		//K->set_array_offsets(K_saved_offsets);
		result->update();
		return result;
	}

}
*/        
/* vim: set ts=4 noet: */
