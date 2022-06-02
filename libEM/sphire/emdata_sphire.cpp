/*
 * Copyright (c) 2019- Baylor College of Medicine
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

#include <stack>
#include "ctf.h"
#include "emdata.h"
#include <iostream>
#include <cmath>
#include <cstring>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <vector>
using std::vector;
using std::cout;
using namespace EMAN;
using namespace std;

#define    QUADPI      		    3.141592653589793238462643383279502884197
#define    PI2                  QUADPI/2.0
#define    TWOPI                2*QUADPI

#define deg_rad  QUADPI/180.0
#define rad_deg  180.0/QUADPI


EMData* EMData::rotavg_sphire() {

	ENTERFUNC;

	int rmax;
	EMData* ret = new EMData();
	vector<int> saved_offsets = get_array_offsets();

	vector<float> count;

	if (ny<2 && nz <2) {
		LOGERR("No 1D images.");
		throw ImageDimensionException("No 1D images!");
	}


	if( this->is_complex() )  {
		//  We will assume square image for the time being
		rmax = ny/2;
		ret->set_size(rmax+1, 1, 1);
		ret->to_zero();
		count.resize(rmax+1);
		set_array_offsets(1,1,1);
		int nz2 = nz/2;
		int ny2 = ny/2;
		int nx2 = nx/2;
		int jx, jy, jz;
		float argy, argz;
			for ( int iz = 1; iz <= nz; iz++) {
				jz=iz-1; if (jz>nz2) jz=jz-nz; argz = float(jz*jz);
				for ( int iy = 1; iy <= ny; iy++) {
					jy=iy-1; if (jy>ny2) jy=jy-ny; argy = argz + float(jy*jy);
					for ( int ix = 1; ix <= nx2; ix++) {
					jx=ix-1;
					float r = std::sqrt(argy + float(jx*jx));
					int  ir = int(r);
					if (ir >= rmax) continue;
					float frac = r - float(ir);
					float qres = 1.0f - frac;
					float temp = std::real(cmplx(ix,iy,iz));
					// cout<<"  "<<jx<<"  "<<jy<<"  "<<ir<<"  "<<temp<<"  "<<frac<<endl;
					(*ret)(ir)   += temp*qres;
					(*ret)(ir+1) += temp*frac;
					count[ir]    += qres;
					count[ir+1]  += frac;
					}
				}
			}

	} else {

		float apix[3];
		apix[0] = get_attr_default("apix_x",1.0);
		apix[1] = get_attr_default("apix_y",1.0);
		apix[2] = get_attr_default("apix_z",1.0);
		float min_apix = *std::min_element(&apix[0],&apix[3]);

		//here,only the relative value of apix_x, apix_y, apix_z are considered
		float apix_x = apix[0]/min_apix;
		float apix_y = apix[1]/min_apix;
		float apix_z = 1.0;

		if( nz > 1)   apix_z=apix[2]/min_apix;

		float apix_x2 = apix_x*apix_x;
		float apix_y2 = apix_y*apix_y;
		float apix_z2 = apix_z*apix_z;

		set_array_offsets(-nx/2,-ny/2,-nz/2);


		//int rmax = std::min(nx/2 + nx%2, ny/2 + ny%2);
		if ( nz == 1 )  rmax = std::min(nx/2 + nx%2, ny/2 + ny%2);
		else            rmax = std::min(nx/2 + nx%2, std::min(ny/2 + ny%2, nz/2 + nz%2));

		float rmax_ratio = 0.0f;
		if      (rmax == nx/2 + nx%2 ) rmax_ratio = apix_x;
		else if (rmax == ny/2 + ny%2)  rmax_ratio = apix_y;
		else                           rmax_ratio = apix_z;

		ret->set_size(rmax+1, 1, 1);
		ret->to_zero();
		count.resize(rmax+1);
		for (int k = -nz/2; k < nz/2 + nz%2; k++) {
			if (abs( k*apix_z) > rmax*rmax_ratio ) continue;
			for (int j = -ny/2; j < ny/2 + ny%2; j++) {
				if (abs( j*apix_y ) > rmax*rmax_ratio) continue;
				for (int i = -nx/2; i < nx/2 + nx%2; i++) {
				    float r = std::sqrt(float(k*k*apix_z2) + float(j*j*apix_y2) + float(i*i*apix_x2))/rmax_ratio;
				    int ir = int(r);
				    if (ir >= rmax) continue;
				    float frac = r - float(ir);
				    (*ret)(ir) += (*this)(i,j,k)*(1.0f - frac);
				    (*ret)(ir+1) += (*this)(i,j,k)*frac;
				    count[ir] += 1.0f - frac;
				    count[ir+1] += frac;

				}
			}
		}
	}
	for (int ir = 0; ir <= rmax; ir++) {
		(*ret)(ir) /= std::max(count[ir],1.0f);
	}
	set_array_offsets(saved_offsets);
	ret->update();
	EXITFUNC;
	return ret;
	}


EMData* EMData::rotavg_i_sphire() {

	int rmax;
	ENTERFUNC;
	if ( ny == 1 && nz == 1 ) {
		LOGERR("Input image must be 2-D or 3-D!");
		throw ImageDimensionException("Input image must be 2-D or 3-D!");
	}

	EMData* avg1D  = new EMData();
	EMData* result = new EMData();

	result->set_size(nx,ny,nz);
	result->to_zero();
	result->set_array_offsets(-nx/2, -ny/2, -nz/2);

	if ( nz == 1 ) {
		rmax = std::min(nx/2 + nx%2, ny/2 + ny%2);
	} else {
		rmax = std::min(nx/2 + nx%2, std::min(ny/2 + ny%2, nz/2 + nz%2));
	}

	avg1D = rotavg_sphire();
	float padded_value = 0.0, r;
	int i, j, k, ir;
	size_t number_of_pixels = 0;
	for ( k = -nz/2; k < nz/2 + nz%2; k++) {
		if (abs(k) > rmax) continue;
		for ( j = -ny/2; j < ny/2 + ny%2; j++) {
			if (abs(j) > rmax) continue;
			for (i = -nx/2; i < nx/2 + nx%2; i++) {
				r = std::sqrt(float(k*k) + float(j*j) + float(i*i));
				ir = int(r);
				if (ir > rmax || ir < rmax-2 ) continue ;
				else {
	      				padded_value += (*avg1D)(ir) ;
	      				number_of_pixels++ ;
				}
			}
		}
	}
	padded_value /= number_of_pixels;
	for ( k = -nz/2; k < nz/2 + nz%2; k++) {
		for ( j = -ny/2; j < ny/2 + ny%2; j++) {
			for ( i = -nx/2; i < nx/2 + nx%2; i++)  {
				r = std::sqrt(float(k*k) + float(j*j) + float(i*i));
				ir = int(r);
				if (ir >= rmax) (*result)(i,j,k) = padded_value ;
				else            (*result)(i,j,k) = (*avg1D)(ir)+((*avg1D)(ir+1)-(*avg1D)(ir))*(r - float(ir));

			}
		}
	}
	result->update();
	result->set_array_offsets(0,0,0);
	EXITFUNC;
	return result;
}
