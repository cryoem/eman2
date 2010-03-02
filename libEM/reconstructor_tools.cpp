/**
 * $Id$
 */

/*
 * Author: David Woolford, 07/25/2007 (woolford@bcm.edu)
 * Copyright (c) 2000-2007 Baylor College of Medicine
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

#include <cstring>
#include "reconstructor_tools.h"

using namespace EMAN;

const string FourierInserter3DMode1::NAME = "nearest_neighbor";
const string FourierInserter3DMode2::NAME = "gauss_2";
const string FourierInserter3DMode3::NAME = "gauss_3";
//const string FourierInserter3DMode4::NAME = "gauss_4";
const string FourierInserter3DMode5::NAME = "gauss_5";
const string FourierInserter3DMode6::NAME = "gauss_5_slow";
const string FourierInserter3DMode7::NAME = "hypergeom_5";
const string FourierInserter3DMode8::NAME = "experimental";

template <> Factory < FourierPixelInserter3D >::Factory()
{
	force_add<FourierInserter3DMode1>();
	force_add<FourierInserter3DMode2>();
	force_add<FourierInserter3DMode3>();
//	force_add<FourierInserter3DMode4>();
	force_add<FourierInserter3DMode5>();
	force_add<FourierInserter3DMode6>();
	force_add<FourierInserter3DMode7>();
//	force_add(&FourierInserter3DMode8::NEW);
}



void FourierPixelInserter3D::init()
{
	if ( params.has_key("data") )
	{
		data = params["data"];
		if ( data == 0 )
			throw NotExistingObjectException("data", "error the data pointer was 0 in FourierPixelInserter3D::init");
	}
	else throw NotExistingObjectException("data", "the data pointer was not defined in FourierPixelInserter3D::init");

	if ( params.has_key("norm"))
	{
		norm = params["norm"];
		if ( norm == 0 )
			throw NotExistingObjectException("norm", "error the norm pointer was 0 in FourierPixelInserter3D::init");
	}
	else throw NotExistingObjectException("norm", "the norm pointer was not defined in FourierPixelInserter3D::init");

	nx=data->get_xsize();
	ny=data->get_ysize();
	nz=data->get_zsize();
	nxyz=nx*ny*nz;
	nx2=nx/2-1;
	ny2=ny/2;
	nz2=nz/2;
	
	if (data->has_attr("subvolume_x0") && data->has_attr("subvolume_full_nx")) {
		subx0=data->get_attr("subvolume_x0");
		suby0=data->get_attr("subvolume_y0");
		subz0=data->get_attr("subvolume_z0");
		fullnx=data->get_attr("subvolume_full_nx");
		fullny=data->get_attr("subvolume_full_ny");
		fullnz=data->get_attr("subvolume_full_nz");
	}
	else {
		subx0=suby0=subz0=-1;
	}
}

bool FourierInserter3DMode1::insert_pixel(const float& xx, const float& yy, const float& zz, const std::complex<float> dt, const float& weight)
{
	int x0 = (int) floor(xx + 0.5f);
	int y0 = (int) floor(yy + 0.5f);
	int z0 = (int) floor(zz + 0.5f);

	size_t off;
	if (subx0<0) off=data->add_complex_at(x0,y0,z0,dt*weight);
	else off=data->add_complex_at(x0,y0,z0,subx0,suby0,subz0,fullnx,fullny,fullnz,dt*weight);
	if (off!=nxyz) norm[off/2]+=weight;
	else return false;
	
	return true;
}

bool FourierInserter3DMode2::insert_pixel(const float& xx, const float& yy, const float& zz, const std::complex<float> dt,const float& weight)
{
	int x0 = (int) floor(xx);
	int y0 = (int) floor(yy);
	int z0 = (int) floor(zz);
	
	if (subx0<0) {			// normal full reconstruction
		if (x0<-nx2-1 || y0<-ny2-1 || z0<-nz2-1 || x0>nx2 || y0>ny2 || z0>nz2 ) return false;

		// no error checking on add_complex_fast, so we need to be careful here
		int x1=x0+1;
		int y1=y0+1;
		int z1=z0+1;
		if (x0<-nx2) x0=-nx2;
		if (x1>nx2) x1=nx2;
		if (y0<-ny2) y0=-ny2;
		if (y1>ny2) y1=ny2;
		if (z0<-nz2) z0=-nz2;
		if (z1>nz2) z1=nz2;
		
//		float h=2.0/((1.0+pow(Util::hypot3sq(xx,yy,zz),.5))*EMConsts::I2G);
		float h=1.0/EMConsts::I2G;
		size_t idx;
		float r, gg;
		int pc=0;
		for (int k = z0 ; k <= z1; k++) {
			for (int j = y0 ; j <= y1; j++) {
				for (int i = x0; i <= x1; i ++) {
					r = Util::hypot3sq((float) i - xx, j - yy, k - zz);
//					gg=weight;
					gg = Util::fast_exp(-r *h)*weight;
//					gg = Util::fast_exp(-r / EMConsts::I2G)*weight;
//					gg = sqrt(Util::fast_exp(-r / EMConsts::I2G))*weight;
					
					size_t off;
					off=data->add_complex_at_fast(i,j,k,dt*gg);
//					off=data->add_complex_at(i,j,k,dt*gg);
					norm[off/2]+=gg;
				}
			}
		}
		return true;
	} 
	else {					// for subvolumes, not optimized yet
		size_t idx;
		float r, gg;
		int pc=0;
		for (int k = z0 ; k <= z0 + 1; k++) {
			for (int j = y0 ; j <= y0 + 1; j++) {
				for (int i = x0; i <= x0 + 1; i ++) {
					r = Util::hypot3sq((float) i - xx, j - yy, k - zz);
					gg = Util::fast_exp(-r / EMConsts::I2G)*weight;

					size_t off;
					if (subx0<0) off=data->add_complex_at(i,j,k,dt*gg);
					else off=data->add_complex_at(i,j,k,subx0,suby0,subz0,fullnx,fullny,fullnz,dt*gg);
					if (off!=nxyz) { norm[off/2]+=gg; pc+=1; }
				}
			}
		}
		
		if (pc>0)  return true;
		return false;
	}
}


bool FourierInserter3DMode3::insert_pixel(const float& xx, const float& yy, const float& zz, const std::complex<float> dt,const float& weight)
{
	int x0 = (int) floor(xx-.5);
	int y0 = (int) floor(yy-.5);
	int z0 = (int) floor(zz-.5);
	
	if (subx0<0) {			// normal full reconstruction
		if (x0<-nx2-2 || y0<-ny2-2 || z0<-nz2-2 || x0>nx2+1 || y0>ny2+1 || z0>nz2+1 ) return false;

		// no error checking on add_complex_fast, so we need to be careful here
		int x1=x0+2;
		int y1=y0+2;
		int z1=z0+2;
		if (x0<-nx2) x0=-nx2;
		if (x1>nx2) x1=nx2;
		if (y0<-ny2) y0=-ny2;
		if (y1>ny2) y1=ny2;
		if (z0<-nz2) z0=-nz2;
		if (z1>nz2) z1=nz2;
		
//		float h=2.0/((1.0+pow(Util::hypot3sq(xx,yy,zz),.5))*EMConsts::I2G);
//		float h=2.0/EMConsts::I3G;
		float h=32.0/((8.0+Util::hypot3(xx,yy,zz))*EMConsts::I3G);
//		float w=weight;
		float w=weight/(1.0+6.0*Util::fast_exp(-h)+12*Util::fast_exp(-h*2.0)+8*Util::fast_exp(-h*3.0));	// approx normalization so higer radii aren't upweighted relative to lower due to wider Gaussian
		size_t idx;
		float r, gg;
		int pc=0;
		for (int k = z0 ; k <= z1; k++) {
			for (int j = y0 ; j <= y1; j++) {
				for (int i = x0; i <= x1; i ++) {
					r = Util::hypot3sq((float) i - xx, (float)j - yy, (float)k - zz);
//					gg=weight;
					gg = Util::fast_exp(-r *h)*w;
//					gg = Util::fast_exp(-r / EMConsts::I2G)*weight;
//					gg = sqrt(Util::fast_exp(-r / EMConsts::I2G))*weight;
					
					size_t off;
					off=data->add_complex_at_fast(i,j,k,dt*gg);
					norm[off/2]+=gg;
				}
			}
		}
		return true;
	}
	printf("region writing not supported in mode 3\n");
	return false;
}


bool FourierInserter3DMode5::insert_pixel(const float& xx, const float& yy, const float& zz, const std::complex<float> dt,const float& weight)
{
	int x0 = (int) floor(xx-1.5);
	int y0 = (int) floor(yy-1.5);
	int z0 = (int) floor(zz-1.5);
	
	if (subx0<0) {			// normal full reconstruction
		if (x0<-nx2-4 || y0<-ny2-4 || z0<-nz2-4 || x0>nx2+3 || y0>ny2+3 || z0>nz2+3 ) return false;

		// no error checking on add_complex_fast, so we need to be careful here
		int x1=x0+4;
		int y1=y0+4;
		int z1=z0+4;
		if (x0<-nx2) x0=-nx2;
		if (x1>nx2) x1=nx2;
		if (y0<-ny2) y0=-ny2;
		if (y1>ny2) y1=ny2;
		if (z0<-nz2) z0=-nz2;
		if (z1>nz2) z1=nz2;
		
//		float h=2.0/((1.0+pow(Util::hypot3sq(xx,yy,zz),.5))*EMConsts::I2G);
//		float h=2.0/EMConsts::I3G;
		float h=32.0/((8.0+Util::hypot3(xx,yy,zz))*EMConsts::I3G); 
		float w=weight/(1.0+6.0*Util::fast_exp(-h)+12*Util::fast_exp(-h*2.0)+8*Util::fast_exp(-h*3.0)+
			6.0*Util::fast_exp(-h*4.0)+24.0*Util::fast_exp(-h*5.0)+24.0*Util::fast_exp(-h*6.0)+12.0*Util::fast_exp(-h*8.0)+
			24.0*Util::fast_exp(-h*9.0)+8.0*Util::fast_exp(-h*12.0));	// approx normalization so higer radii aren't upweighted relative to lower due to wider Gaussian
		size_t idx;
		float r, gg;
		int pc=0;
		for (int k = z0 ; k <= z1; k++) {
			for (int j = y0 ; j <= y1; j++) {
				for (int i = x0; i <= x1; i ++) {
					r = Util::hypot3sq((float) i - xx, j - yy, k - zz);
//					gg=weight;
					gg = Util::fast_exp(-r *h)*w;
//					gg = Util::fast_exp(-r / EMConsts::I2G)*weight;
//					gg = sqrt(Util::fast_exp(-r / EMConsts::I2G))*weight;
					
					size_t off;
					off=data->add_complex_at_fast(i,j,k,dt*gg);
					norm[off/2]+=gg;
				}
			}
		}
		return true;
	}
	printf("region writing not supported in mode 3\n");
	return false;
}


bool FourierInserter3DMode6::insert_pixel(const float& xx, const float& yy, const float& zz, const std::complex<float> dt,const float& weight)
{
	int x0 = 2 * (int) floor(xx + 0.5f);
	int y0 = (int) floor(yy + 0.5f);
	int z0 = (int) floor(zz + 0.5f);

	if (x0<-6||y0<-3||z0<-3||x0>nx+6||y0>ny+3||z0>nz+3) return false;

	size_t idx;
	float r, gg;
	for (int k = z0 - 1; k <= z0 + 1; k++) {
		for (int j = y0 - 1; j <= y0 + 1; j++) {
			for (int i = x0 -2; i <= x0 + 2; i += 2) {
				if (k<0 || j<0 || i<0 || k>=nz || j>=ny || i>=nx) continue;
				r = Util::hypot3sq((float) i / 2 - xx, j - yy, k - zz);
				gg = exp(-r / EMConsts::I5G)*weight;

				size_t off;
				if (subx0<0) off=data->add_complex_at(i,j,k,dt*gg);
				else off=data->add_complex_at(i,j,k,subx0,suby0,subz0,fullnx,fullny,fullnz,dt*gg);
				if (off!=nxyz) norm[off/2]+=gg;

			}
		}
	}
	return true;
 }


bool FourierInserter3DMode7::insert_pixel(const float& xx, const float& yy, const float& zz, const std::complex<float> dt,const float& weight)
{
	int x0 = 2 * (int) floor(xx + 0.5f);
	int y0 = (int) floor(yy + 0.5f);
	int z0 = (int) floor(zz + 0.5f);

	if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
		return false;

	int l = x0 - 4;
	if (x0 <= 2)
		l = 0;

	size_t ii;
	float r, gg;
	for (int k = z0 - 2; k <= z0 + 2; k++) {
		for (int j = y0 - 2; j <= y0 + 2; j++) {
			for (int i = l; i <= x0 + 4; i += 2) {
				r =	sqrt(Util::hypot3sq((float) i / 2 - xx, (float) j - yy, (float) k - zz));
				gg = Interp::hyperg(r)*weight;

				size_t off;
				if (subx0<0) off=data->add_complex_at(i,j,k,dt*gg);
				else off=data->add_complex_at(i,j,k,subx0,suby0,subz0,fullnx,fullny,fullnz,dt*gg);
				if (off!=nxyz) norm[off/2]+=gg;

			}
		}
	}

	if (x0 <= 2) {
		float xx_b = -xx;
		float yy_b = -(yy - ny / 2) + ny / 2;
		float zz_b = -(zz - nz / 2) + nz / 2;
		x0 = 2 * (int) floor(xx_b + 0.5f);
		y0 = (int) floor(yy_b + 0.5f);
		z0 = (int) floor(zz_b + 0.5f);

		if (y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
			return false;

		size_t ii;
		float r, gg;
		for (int k = z0 - 2; k <= z0 + 2; k++) {
			for (int j = y0 - 2; j <= y0 + 2; j++) {
				for (int i = 0; i <= x0 + 4; i += 2) {
					r = sqrt(Util::hypot3sq((float) i / 2 - xx_b, (float) j - yy_b,
								   (float) k - zz_b));
					gg = Interp::hyperg(r)*weight;

					size_t off;
					if (subx0<0) off=data->add_complex_at(i,j,k,dt*gg);
					else off=data->add_complex_at(i,j,k,subx0,suby0,subz0,fullnx,fullny,fullnz,dt*gg);
					if (off!=nxyz) norm[off/2]+=gg;

				}
			}
		}
	}

	return true;
}


void FourierInserter3DMode8::init()
{
	FourierPixelInserter3D::init();
	int P = (int)((1.0+0.25)*nx+1);
	float r = (float)(nx+1)/(float)P;
	mFreqCutoff = 2;
	mDFreq = 0.2f;
	if (W != 0) delete [] W;
	W = Util::getBaldwinGridWeights(mFreqCutoff, (float)P, r,mDFreq,0.5f,0.2f);

}
bool FourierInserter3DMode8::insert_pixel(const float& qx, const float& qy, const float& qz, const std::complex<float> fq,const float& weight)
{
	int x0 = (int) floor(qx);
	int y0 = (int) floor(qy);
	int z0 = (int) floor(qz);

	int sizeW = (int)(1+2*mFreqCutoff/mDFreq);
	int sizeWmid = sizeW/2;

// 	for (int z = z0-mFreqCutoff; z < z0+mFreqCutoff; ++z){
// 		for (int y = y0-mFreqCutoff; y < y0+mFreqCutoff; ++y){
// 			for (int x = x0-mFreqCutoff; x < x0+mFreqCutoff; ++x){
// 				if ( x < 0 || x >= nx ) continue;
// 				if ( y < 0 || y >= ny ) continue;
// 				if ( z < 0 || z >= nz ) continue;
// 				float dist = (float)((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0));
// 				dist = sqrtf(dist);
// 				// We enforce a spherical kernel
// 				if ( dist > mFreqCutoff ) continue;
// 				int idx = (int)(sizeWmid + dist/mDFreq);
// 				if (idx >= sizeW) throw;
// 				float residual = dist/mDFreq - (int)(dist/mDFreq);
// 				if ( fabs(residual) > 1) throw;
// 
// 				float factor = W[idx]*(1.0f-residual)+W[idx+1]*residual*weight;
// 
// 				size_t k = z*nx*ny + y*nx + 2*x;
// 
// // 				float c = Util::agauss(1, x-x0,y-y0,z-z0, EMConsts::I2G);
// 				rdata[k] += fq[0]*factor;
// 				rdata[k+1] += fq[1]*factor;
// 
// 
// 				norm[k/2] += weight;
// 
// 			}
// 		}
// 	}

	return true;
}
