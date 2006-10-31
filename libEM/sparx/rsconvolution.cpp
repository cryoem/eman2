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
}

namespace EMAN {

    EMData* rsconvolution(EMData* f, EMData* K) {
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
		K->set_array_offsets(K_saved_offsets);
		result->done_data();
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
		result->done_data();
		return result;
	}

}
*/        
/* vim: set ts=4 noet: */
