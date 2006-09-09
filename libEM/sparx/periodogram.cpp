/*
� * Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include "emdata.h"

using namespace std;
using namespace EMAN;

namespace EMAN {
/*
******************************************************
*DISCLAIMER
* 05/25/05 P.A.Penczek
* The University of Texas
* Pawel.A.Penczek@uth.tmc.edu
* Please do not modify the content of this document
* without a written permission of the author.
*/
/*
Periodogram
Purpose: Calculate a periodogram (a squared modulus of the FT) of a 1-2-3D image. 
Method: Calculate FFT (if needed), squared modulus, shift the origin to n/2+1,
        create the Friedel-related part.
	Real input -> padded workspace -> real output -> delete workspace.
	Complex input -> copy input to workspace -> real output -> delete workspace.
Input: f real or complex 1-2-3D image
Output: 1-2-3D real image with periodogram.  |F(0,0,0)|^2 is at the image center
        int[n/2]+1.
*/  
/*
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
	else fp= norm_pad_ft(f, false, false); // Extend and do the FFT if f is real

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
		
	if( fp )
	{
		delete fp; // avoid a memory leak!
		fp = 0;
	}
	//power[0][0][0]=power[1][0][0];  //Steve requested the original origin.

	power.done_data();
	power.set_array_offsets(0,0,0);
	return &power;
//OVER AND OUT
}
*/

}
/* vim: set ts=4 noet nospell: */
