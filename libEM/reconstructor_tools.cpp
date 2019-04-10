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
#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include "reconstructor_tools.h"

using namespace EMAN;

const string FourierInserter3DMode1::NAME = "nearest_neighbor";
const string FourierInserter3DMode2::NAME = "gauss_2";
const string FourierInserter3DMode3::NAME = "gauss_3";
//const string FourierInserter3DMode4::NAME = "gauss_4";
const string FourierInserter3DMode5::NAME = "gauss_5";
const string FourierInserter3DMode6::NAME = "gauss_var";
const string FourierInserter3DMode7::NAME = "gridding_5";
const string FourierInserter3DMode11::NAME = "gridding_7";
const string FourierInserter3DMode8::NAME = "experimental";
const string FourierInserter3DMode9::NAME = "kaiser_bessel";
const string FourierInserter3DMode10::NAME = "kaiser_bessel_derived";

template <> Factory < FourierPixelInserter3D >::Factory()
{
	force_add<FourierInserter3DMode1>();
	force_add<FourierInserter3DMode2>();
	force_add<FourierInserter3DMode3>();
//	force_add<FourierInserter3DMode4>();
	force_add<FourierInserter3DMode5>();
	force_add<FourierInserter3DMode6>();
	force_add<FourierInserter3DMode7>();
	force_add<FourierInserter3DMode8>();
//	force_add(&FourierInserter3DMode8::NEW);
	force_add<FourierInserter3DMode9>();
	force_add<FourierInserter3DMode10>();
	force_add<FourierInserter3DMode11>();
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
	nxyz=(size_t)nx*ny*nz;
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
	if (static_cast<int>(off)!=nxyz) norm[off/2]+=weight;
	else return false;

	return true;
}

bool FourierInserter3DMode2::insert_pixel(const float& xx, const float& yy, const float& zz, const std::complex<float> dt,const float& weight)
{
	int x0 = (int) floor(xx);
	int y0 = (int) floor(yy);
	int z0 = (int) floor(zz);

	// note that subnx differs in the inserters. In the reconstructors it subx0 is 0 for the full volume. Here it is -1
	if (subx0<0) {			// normal full reconstruction
		if (x0<-nx2-1 || y0<-ny2-1 || z0<-nz2-1 || x0>nx2 || y0>ny2 || z0>nz2 ) return false;

		int x1=x0+1;
		int y1=y0+1;
		int z1=z0+1;
// 		if (x0<-nx2) x0=-nx2;
// 		if (x1>nx2) x1=nx2;
// 		if (y0<-ny2) y0=-ny2;
// 		if (y1>ny2) y1=ny2;
// 		if (z0<-nz2) z0=-nz2;
// 		if (z1>nz2) z1=nz2;

//		float h=2.0/((1.0+pow(Util::hypot3sq(xx,yy,zz),.5))*EMConsts::I2G);
		float h=1.0f/EMConsts::I2G;
		//size_t idx;
		float r, gg;
//		int pc=0;
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
					if (off!=nxyz) norm[off/2]+=gg;
				}
			}
		}
		return true;
	}
	else {					// for subvolumes, not optimized yet
		//size_t idx;
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
					if (static_cast<int>(off)!=nxyz) { norm[off/2]+=gg; pc+=1; }
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
		float h=32.0f/((8.0f+Util::hypot3(xx,yy,zz))*EMConsts::I3G);
//		float w=weight;
		float w=weight/(1.0f+6.0f*Util::fast_exp(-h)+12*Util::fast_exp(-h*2.0f)+8*Util::fast_exp(-h*3.0f));	// approx normalization so higer radii aren't upweighted relative to lower due to wider Gaussian
		//size_t idx;
		float r, gg;
//		int pc=0;
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
	int x0 = (int) floor(xx-2.5);
	int y0 = (int) floor(yy-2.5);
	int z0 = (int) floor(zz-2.5);

	if (subx0<0) {			// normal full reconstruction
		if (x0<-nx2-4 || y0<-ny2-4 || z0<-nz2-4 || x0>nx2+3 || y0>ny2+3 || z0>nz2+3 ) return false;

		// no error checking on add_complex_fast, so we need to be careful here
		int x1=x0+5;
		int y1=y0+5;
		int z1=z0+5;
		if (x0<-nx2) x0=-nx2;
		if (x1>nx2) x1=nx2;
		if (y0<-ny2) y0=-ny2;
		if (y1>ny2) y1=ny2;
		if (z0<-nz2) z0=-nz2;
		if (z1>nz2) z1=nz2;

//		float h=2.0/((1.0+pow(Util::hypot3sq(xx,yy,zz),.5))*EMConsts::I2G);
		float h=1.0f/EMConsts::I5G;
		float w=weight;

		// Not sure exactly what this was doing? Using wider Gaussian at high radius?
// 		float h=32.0f/((8.0f+Util::hypot3(xx,yy,zz))*EMConsts::I3G);
// 		float w=weight/(1.0f+6.0f*Util::fast_exp(-h)+12*Util::fast_exp(-h*2.0f)+8*Util::fast_exp(-h*3.0f)+
// 			6.0f*Util::fast_exp(-h*4.0f)+24.0f*Util::fast_exp(-h*5.0f)+24.0f*Util::fast_exp(-h*6.0f)+12.0f*Util::fast_exp(-h*8.0f)+
// 			24.0f*Util::fast_exp(-h*9.0f)+8.0f*Util::fast_exp(-h*12.0f));	// approx normalization so higer radii aren't upweighted relative to lower due to wider Gaussian
		//size_t idx;
		float r, gg;
//		int pc=0;
		for (int k = z0 ; k <= z1; k++) {
			for (int j = y0 ; j <= y1; j++) {
				for (int i = x0; i <= x1; i ++) {
					r = Util::hypot3sq((float) i - xx, j - yy, k - zz);
//					gg=weight;
					gg = Util::fast_exp(-r *h);
//					gg = Util::fast_exp(-r / EMConsts::I2G)*weight;
//					gg = sqrt(Util::fast_exp(-r / EMConsts::I2G))*weight;

					size_t off;
					off=data->add_complex_at_fast(i,j,k,dt*gg*w);
					norm[off/2]+=gg*w;		// This would use a Gaussian WEIGHT with square kernel
//					norm[off/2]+=w;			// This would use a Gaussian KERNEL rather than WEIGHT 

#ifdef RECONDEBUG
					std::complex<double> v1=dt*gg*w,v2=gg*w;

					if (k<5 && j<5&& i<5&& k>=0 && j>=0 && i>=0) {
						int idx=i*2+j*10+k*50;
						ddata[idx]+=v1.real();
						ddata[idx+1]+=v1.imag();
						dnorm[idx]+=v2.real();
						dnorm[idx+1]+=v2.imag();
					}
#endif
				}
			}
		}
		return true;
	}
	printf("region writing not supported in mode 5\n");
	return false;
}


bool FourierInserter3DMode6::insert_pixel(const float& xx, const float& yy, const float& zz, const std::complex<float> dt,const float& weight)
{
	int x0 = (int) floor(xx-2.5);
	int y0 = (int) floor(yy-2.5);
	int z0 = (int) floor(zz-2.5);

	if (subx0<0) {			// normal full reconstruction
		if (x0<-nx2-4 || y0<-ny2-4 || z0<-nz2-4 || x0>nx2+3 || y0>ny2+3 || z0>nz2+3 ) return false;

		// no error checking on add_complex_fast, so we need to be careful here
		int x1=x0+5;
		int y1=y0+5;
		int z1=z0+5;
		if (x0<-nx2) x0=-nx2;
		if (x1>nx2) x1=nx2;
		if (y0<-ny2) y0=-ny2;
		if (y1>ny2) y1=ny2;
		if (z0<-nz2) z0=-nz2;
		if (z1>nz2) z1=nz2;

//		float h=2.0/((1.0+pow(Util::hypot3sq(xx,yy,zz),.5))*EMConsts::I2G);
//		float h=1.0f/EMConsts::I5G;
//		float h=1.0f/(Util::hypot3sq(xx/nx2,yy/ny2,zz/nz2)*EMConsts::I5G*2.0+.1);		// gaussian kernel is used as a weight not a kernel. We increase radius of integration with resolution. Changed away from this on 11/10/17
		float h=1.0f/((Util::hypot3sq(xx,yy,zz)/4000)+.15);		// gaussian kernel is used as a weight not a kernel. We increase radius of integration with resolution. Changed away from this on 11/10/17
		h=h<0.1?0.1:h;	// new formula has h from 0.2 to 5
//		if (yy==0 &&zz==0) printf("%0.0f\t%0.0f\t%0.0f\t%g\n",xx,yy,zz,h);
//		printf("%1.0f\t%1.0f\t%1.0f\t%1.4f\n",xx,yy,zz,h);
		float w=weight;

// 		float h=32.0f/((8.0f+Util::hypot3(xx,yy,zz))*EMConsts::I3G);
// 		float w=weight/(1.0f+6.0f*Util::fast_exp(-h)+12*Util::fast_exp(-h*2.0f)+8*Util::fast_exp(-h*3.0f)+
// 			6.0f*Util::fast_exp(-h*4.0f)+24.0f*Util::fast_exp(-h*5.0f)+24.0f*Util::fast_exp(-h*6.0f)+12.0f*Util::fast_exp(-h*8.0f)+
// 			24.0f*Util::fast_exp(-h*9.0f)+8.0f*Util::fast_exp(-h*12.0f));	// approx normalization so higer radii aren't upweighted relative to lower due to wider Gaussian
		//size_t idx;
		float r, gg;
//		int pc=0;
		for (int k = z0 ; k <= z1; k++) {
			for (int j = y0 ; j <= y1; j++) {
				for (int i = x0; i <= x1; i ++) {
					r = Util::hypot3sq((float) i - xx, j - yy, k - zz);
//					gg=weight;
					gg = Util::fast_exp(-r *h);
//					if (gg<.00001) continue;		// skip tiny weights for speed
//					gg = Util::fast_exp(-r / EMConsts::I2G)*weight;
//					gg = sqrt(Util::fast_exp(-r / EMConsts::I2G))*weight;

					size_t off;
					off=data->add_complex_at_fast(i,j,k,dt*gg*w);
//					if (i==16 && j==7 && k==4) printf("%g\t%g\t%g\t%g\t%g\n",dt.real(),dt.imag(),gg*w,gg,w);
					norm[off/2]+=gg*w;		// This would use a Gaussian WEIGHT with square kernel
//					norm[off/2]+=w;			// This would use a Gaussian KERNEL rather than WEIGHT 

#ifdef RECONDEBUG
					std::complex<double> v1=dt*gg*w,v2=gg*w;

					if (k<5 && j<5&& i<5&& k>=0 && j>=0 && i>=0) {
						int idx=i*2+j*10+k*50;
						ddata[idx]+=v1.real();
						ddata[idx+1]+=v1.imag();
						dnorm[idx]+=v2.real();
						dnorm[idx+1]+=v2.imag();
					}
#endif
				}
			}
		}
		return true;
	}
	printf("region writing not supported in mode 5\n");
	return false;
}

// Kernel determined by examples/kernel_opt
const float FourierInserter3DMode7::kernel[9][9][9] = {
0.3293228,0.2999721,0.2230405,0.1263839,0.0412195,-0.0116440,-0.0287601,-0.0214417,0.0000000,
0.2999721,0.2731577,0.2028874,0.1146381,0.0369511,-0.0111675,-0.0266030,-0.0197382,0.0000000,
0.2230405,0.2028874,0.1501083,0.0839287,0.0258589,-0.0098285,-0.0208839,-0.0152399,0.0000000,
0.1263839,0.1146381,0.0839287,0.0455768,0.0122079,-0.0078769,-0.0135036,-0.0094886,0.0000000,
0.0412195,0.0369511,0.0258589,0.0122079,0.0007022,-0.0056645,-0.0066442,-0.0042390,0.0000000,
-0.0116440,-0.0111675,-0.0098285,-0.0078769,-0.0056645,-0.0035595,-0.0018569,-0.0007100,0.0000000,
-0.0287601,-0.0266030,-0.0208839,-0.0135036,-0.0066442,-0.0018569,0.0004320,0.0008101,0.0000000,
-0.0214417,-0.0197382,-0.0152398,-0.0094886,-0.0042390,-0.0007100,0.0008101,0.0008532,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.2999721,0.2731577,0.2028874,0.1146381,0.0369511,-0.0111675,-0.0266030,-0.0197382,0.0000000,
0.2731577,0.2486666,0.1844966,0.1039437,0.0330965,-0.0106899,-0.0246017,-0.0181663,0.0000000,
0.2028874,0.1844966,0.1363422,0.0759894,0.0230832,-0.0093549,-0.0192968,-0.0140161,0.0000000,
0.1146381,0.1039437,0.0759894,0.0410967,0.0107719,-0.0074293,-0.0124545,-0.0087120,0.0000000,
0.0369511,0.0330965,0.0230832,0.0107719,0.0004174,-0.0052781,-0.0061010,-0.0038740,0.0000000,
-0.0111675,-0.0106899,-0.0093549,-0.0074294,-0.0052781,-0.0032679,-0.0016753,-0.0006270,0.0000000,
-0.0266030,-0.0246017,-0.0192968,-0.0124545,-0.0061010,-0.0016753,0.0004302,0.0007650,0.0000000,
-0.0197382,-0.0181663,-0.0140161,-0.0087120,-0.0038740,-0.0006270,0.0007650,0.0007953,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.2230405,0.2028874,0.1501083,0.0839287,0.0258589,-0.0098285,-0.0208839,-0.0152399,0.0000000,
0.2028874,0.1844966,0.1363422,0.0759894,0.0230832,-0.0093549,-0.0192968,-0.0140161,0.0000000,
0.1501083,0.1363422,0.1003225,0.0552536,0.0158836,-0.0080490,-0.0150930,-0.0107868,0.0000000,
0.0839287,0.0759894,0.0552535,0.0294215,0.0070640,-0.0062153,-0.0096798,-0.0066651,0.0000000,
0.0258589,0.0230832,0.0158836,0.0070640,-0.0002933,-0.0042462,-0.0046694,-0.0029154,0.0000000,
-0.0098285,-0.0093549,-0.0080490,-0.0062153,-0.0042462,-0.0024995,-0.0012019,-0.0004126,0.0000000,
-0.0208839,-0.0192968,-0.0150930,-0.0096798,-0.0046694,-0.0012019,0.0004197,0.0006424,0.0000000,
-0.0152398,-0.0140161,-0.0107868,-0.0066651,-0.0029154,-0.0004126,0.0006424,0.0006405,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.1263839,0.1146381,0.0839287,0.0455768,0.0122079,-0.0078769,-0.0135036,-0.0094886,0.0000000,
0.1146381,0.1039437,0.0759894,0.0410967,0.0107719,-0.0074294,-0.0124545,-0.0087120,0.0000000,
0.0839287,0.0759894,0.0552535,0.0294215,0.0070640,-0.0062153,-0.0096798,-0.0066651,0.0000000,
0.0455767,0.0410967,0.0294215,0.0149533,0.0025724,-0.0045680,-0.0061197,-0.0040604,0.0000000,
0.0122078,0.0107719,0.0070640,0.0025724,-0.0010797,-0.0028931,-0.0028473,-0.0017048,0.0000000,
-0.0078770,-0.0074294,-0.0062153,-0.0045680,-0.0028931,-0.0015217,-0.0006148,-0.0001526,0.0000000,
-0.0135036,-0.0124545,-0.0096798,-0.0061197,-0.0028473,-0.0006148,0.0003888,0.0004754,0.0000000,
-0.0094886,-0.0087120,-0.0066651,-0.0040604,-0.0017048,-0.0001527,0.0004754,0.0004373,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.0412195,0.0369511,0.0258589,0.0122079,0.0007022,-0.0056645,-0.0066442,-0.0042390,0.0000000,
0.0369511,0.0330965,0.0230832,0.0107719,0.0004174,-0.0052781,-0.0061010,-0.0038740,0.0000000,
0.0258589,0.0230832,0.0158836,0.0070640,-0.0002933,-0.0042462,-0.0046694,-0.0029154,0.0000000,
0.0122078,0.0107719,0.0070640,0.0025724,-0.0010797,-0.0028931,-0.0028473,-0.0017048,0.0000000,
0.0007022,0.0004174,-0.0002933,-0.0010797,-0.0015774,-0.0015967,-0.0011988,-0.0006268,0.0000000,
-0.0056645,-0.0052781,-0.0042463,-0.0028931,-0.0015967,-0.0006377,-0.0001117,0.0000590,0.0000000,
-0.0066442,-0.0061010,-0.0046694,-0.0028473,-0.0011988,-0.0001117,0.0003294,0.0003044,0.0000000,
-0.0042390,-0.0038740,-0.0029154,-0.0017048,-0.0006268,0.0000590,0.0003044,0.0002424,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
-0.0116440,-0.0111675,-0.0098285,-0.0078769,-0.0056645,-0.0035595,-0.0018569,-0.0007100,0.0000000,
-0.0111675,-0.0106899,-0.0093549,-0.0074294,-0.0052781,-0.0032679,-0.0016753,-0.0006270,0.0000000,
-0.0098285,-0.0093549,-0.0080490,-0.0062153,-0.0042463,-0.0024995,-0.0012019,-0.0004126,0.0000000,
-0.0078770,-0.0074294,-0.0062153,-0.0045680,-0.0028931,-0.0015217,-0.0006148,-0.0001526,0.0000000,
-0.0056645,-0.0052781,-0.0042463,-0.0028931,-0.0015967,-0.0006377,-0.0001117,0.0000590,0.0000000,
-0.0035595,-0.0032679,-0.0024995,-0.0015217,-0.0006377,-0.0000555,0.0001796,0.0001647,0.0000000,
-0.0018569,-0.0016753,-0.0012019,-0.0006148,-0.0001117,0.0001796,0.0002446,0.0001630,0.0000000,
-0.0007100,-0.0006270,-0.0004126,-0.0001527,0.0000590,0.0001647,0.0001630,0.0000978,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
-0.0287600,-0.0266030,-0.0208839,-0.0135036,-0.0066442,-0.0018569,0.0004320,0.0008101,0.0000000,
-0.0266030,-0.0246017,-0.0192968,-0.0124545,-0.0061010,-0.0016753,0.0004302,0.0007650,0.0000000,
-0.0208839,-0.0192968,-0.0150930,-0.0096798,-0.0046694,-0.0012019,0.0004197,0.0006424,0.0000000,
-0.0135036,-0.0124545,-0.0096798,-0.0061197,-0.0028473,-0.0006148,0.0003888,0.0004754,0.0000000,
-0.0066442,-0.0061010,-0.0046694,-0.0028473,-0.0011988,-0.0001117,0.0003294,0.0003044,0.0000000,
-0.0018569,-0.0016753,-0.0012019,-0.0006148,-0.0001117,0.0001796,0.0002446,0.0001630,0.0000000,
0.0004320,0.0004302,0.0004197,0.0003888,0.0003294,0.0002446,0.0001503,0.0000679,0.0000000,
0.0008101,0.0007650,0.0006424,0.0004754,0.0003044,0.0001630,0.0000679,0.0000181,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
-0.0214417,-0.0197382,-0.0152398,-0.0094886,-0.0042390,-0.0007100,0.0008101,0.0008532,0.0000000,
-0.0197382,-0.0181663,-0.0140161,-0.0087120,-0.0038740,-0.0006270,0.0007650,0.0007953,0.0000000,
-0.0152398,-0.0140161,-0.0107868,-0.0066651,-0.0029154,-0.0004126,0.0006424,0.0006405,0.0000000,
-0.0094886,-0.0087120,-0.0066651,-0.0040604,-0.0017048,-0.0001527,0.0004754,0.0004373,0.0000000,
-0.0042390,-0.0038740,-0.0029154,-0.0017048,-0.0006268,0.0000590,0.0003044,0.0002424,0.0000000,
-0.0007100,-0.0006270,-0.0004126,-0.0001527,0.0000590,0.0001647,0.0001630,0.0000978,0.0000000,
0.0008101,0.0007650,0.0006424,0.0004754,0.0003044,0.0001630,0.0000679,0.0000181,0.0000000,
0.0008532,0.0007953,0.0006405,0.0004373,0.0002424,0.0000978,0.0000181,-0.0000082,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
};

bool FourierInserter3DMode7::insert_pixel(const float& xx, const float& yy, const float& zz, const std::complex<float> dt,const float& weight)
{
	int x0 = Util::fast_floor(xx-1.5);
	int y0 = Util::fast_floor(yy-1.5);
	int z0 = Util::fast_floor(zz-1.5);

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

		float w=weight;

		float r, gg;
		for (int k = z0 ; k <= z1; k++) {
			for (int j = y0 ; j <= y1; j++) {
				for (int i = x0; i <= x1; i ++) {
					gg=FourierInserter3DMode7::kernel[abs(Util::fast_floor((i-xx)*3.0f+0.5))][abs(Util::fast_floor((j-yy)*3.0f+0.5))][abs(Util::fast_floor((k-zz)*3.0f+0.5))]; 

					size_t off;
					off=data->add_complex_at_fast(i,j,k,dt*gg*w);
					norm[off/2]+=w;
//					if (i==67&&j==19&&k==1) printf("%1.1f  %1.1f  %1.1f\t%d %d %d\t%d %d %d\t%f\t%f\t%f\t%f\t%f\t%f\n",xx,yy,zz,i,j,k,abs(Util::fast_floor((i-xx)*3.0f+0.5)),abs(Util::fast_floor((j-yy)*3.0f+0.5)),abs(Util::fast_floor((k-zz)*3.0f+0.5)),gg,w,dt.real(),dt.imag(),norm[off/2],data->get_value_at(off));

				}
			}
		}
		return true;
	}
	printf("region writing not supported in mode \n");
	return false;
}


void FourierInserter3DMode8::init()
{
	FourierPixelInserter3D::init();
// 	int P = (int)((1.0+0.25)*nx+1);
// 	float r = (float)(nx+1)/(float)P;
// 	mFreqCutoff = 2;
// 	mDFreq = 0.2f;
// 	if (W != 0) delete [] W;
// 	W = Util::getBaldwinGridWeights(mFreqCutoff, (float)P, r,mDFreq,0.5f,0.2f);

}

bool FourierInserter3DMode8::insert_pixel(const float& xx, const float& yy, const float& zz, const std::complex<float> dt,const float& weight)
{
	static FILE *out400=NULL,*out862=NULL,*out962,*out872,*out4093=NULL;
	
	int x0 = (int) floor(xx);
	int y0 = (int) floor(yy);
	int z0 = (int) floor(zz);

	if (out400==NULL) {
		out400=fopen("pxl4_0_0.txt","w");
		out862=fopen("pxl8_6_2.txt","w");
		out962=fopen("pxl9_6_2.txt","w");
		out872=fopen("pxl8_7_2.txt","w");
		out4093=fopen("pxl40_9_3.txt","w");
	}
		
	// note that subnx differs in the inserters. In the reconstructors it subx0 is 0 for the full volume. Here it is -1
	if (subx0<0) {			// normal full reconstruction
		if (x0<-nx2-1 || y0<-ny2-1 || z0<-nz2-1 || x0>nx2 || y0>ny2 || z0>nz2 ) return false;

		int x1=x0+1;
		int y1=y0+1;
		int z1=z0+1;
// 		if (x0<-nx2) x0=-nx2;
// 		if (x1>nx2) x1=nx2;
// 		if (y0<-ny2) y0=-ny2;
// 		if (y1>ny2) y1=ny2;
// 		if (z0<-nz2) z0=-nz2;
// 		if (z1>nz2) z1=nz2;

//		float h=2.0/((1.0+pow(Util::hypot3sq(xx,yy,zz),.5))*EMConsts::I2G);
		float h=1.0f/EMConsts::I2G;
		//size_t idx;
		float r, gg;
//		int pc=0;
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
					if (off!=nxyz) norm[off/2]+=gg;
					
					if (i==4&&j==0&&k==0) { fprintf(out400,"%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\n",dt.real(),dt.imag(),gg,std::abs(dt),std::arg(dt)); fflush(out400); }
					if (i==8&&j==6&&k==2) { fprintf(out862,"%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\n",dt.real(),dt.imag(),gg,std::abs(dt),std::arg(dt)); fflush(out862); }
					if (i==9&&j==6&&k==2) { fprintf(out962,"%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\n",dt.real(),dt.imag(),gg,std::abs(dt),std::arg(dt)); fflush(out962); }
					if (i==8&&j==7&&k==2) { fprintf(out872,"%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\n",dt.real(),dt.imag(),gg,std::abs(dt),std::arg(dt)); fflush(out872); }
					if (i==40&&j==9&&k==3) { fprintf(out4093,"%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\n",dt.real(),dt.imag(),gg,std::abs(dt),std::arg(dt)); fflush(out4093); }
				}
			}
		}
		return true;
	}
	else {					// for subvolumes, not optimized yet
		//size_t idx;
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
					if (static_cast<int>(off)!=nxyz) { norm[off/2]+=gg; pc+=1; }
				}
			}
		}

		if (pc>0)  return true;
		return false;
	}
//	int x0 = (int) floor(qx);
//	int y0 = (int) floor(qy);
//	int z0 = (int) floor(qz);

//	int sizeW = (int)(1+2*mFreqCutoff/mDFreq);
//	int sizeWmid = sizeW/2;

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

bool FourierInserter3DMode9::insert_pixel(const float& xx, const float& yy, const float& zz, const std::complex<float> dt,const float& weight)
{
	int N = 4;		// kernel width

	int x0 = (int) floor(xx-N/2);
	int y0 = (int) floor(yy-N/2);
	int z0 = (int) floor(zz-N/2);

	if (subx0<0) {		// normal full reconstruction
		if (x0<-nx2-7 || y0<-ny2-7 || z0<-nz2-7 || x0>nx2+6 || y0>ny2+6 || z0>nz2+6 ) return false;

		// no error checking on add_complex_fast, so we need to be careful here
		int x1=x0+N;
		int y1=y0+N;
		int z1=z0+N;

		if (x0<-nx2) x0=-nx2;
		if (x1>nx2) x1=nx2;
		if (y0<-ny2) y0=-ny2;
		if (y1>ny2) y1=ny2;
		if (z0<-nz2) z0=-nz2;
		if (z1>nz2) z1=nz2;

		float w=weight;
		float a=15.0;
		float r, kb;

		for (int k = z0 ; k <= z1; k++) {
			for (int j = y0 ; j <= y1; j++) {
				for (int i = x0; i <= x1; i ++) {
					r = Util::hypot3sq((float) i - xx, j - yy, k - zz);
					kb = gsl_sf_bessel_i0_scaled(M_PI * a * sqrt(1.0f - Util::square((r/(nx2-1))-1))) /
					     gsl_sf_bessel_i0_scaled(M_PI * a);
					size_t off;
					off = data->add_complex_at_fast(i,j,k,dt*kb*w);
					norm[off/2]+=w;
				}
			}
		}
		return true;
	}
	printf("region writing not supported in mode 9\n");
	return false;
}

// imprecise KBD kernel/window
bool FourierInserter3DMode10::insert_pixel(const float& xx, const float& yy, const float& zz, const std::complex<float> dt,const float& weight)
{
	const int N = 8;		// kernel width

	int x0 = (int) floor(xx-N/2);
	int y0 = (int) floor(yy-N/2);
	int z0 = (int) floor(zz-N/2);

	if (subx0<0) {		// normal full reconstruction
		if (x0<-nx2-7 || y0<-ny2-7 || z0<-nz2-7 || x0>nx2+6 || y0>ny2+6 || z0>nz2+6 ) return false;

		// no error checking on add_complex_fast, so we need to be careful here
		int x1=x0+N;
		int y1=y0+N;
		int z1=z0+N;

		if (x0<-nx2) x0=-nx2;
		if (x1>nx2) x1=nx2;
		if (y0<-ny2) y0=-ny2;
		if (y1>ny2) y1=ny2;
		if (z0<-nz2) z0=-nz2;
		if (z1>nz2) z1=nz2;

		float w=weight;
		float ws [ N/2 + 1 ];
		float alpha = 32.0;
		float wm = 0.0;

		// compute 1D weights... not exactly correct, but somewhat close.
		for ( int p = 0; p <= N/2; p++) {
			double tmp = gsl_sf_bessel_i0_scaled(M_PI * alpha * sqrt(1.0f - Util::square((((N/2)+p)/(N-1))-1))) / gsl_sf_bessel_i0_scaled(M_PI * alpha);
			ws[p] = (float) tmp;
			wm += (float) tmp;
		}

		float r, kb, dn;
		for (int k = z0 ; k <= z1; k++) {
			for (int j = y0 ; j <= y1; j++) {
				for (int i = x0; i <= x1; i ++) {
					r = Util::hypot3sq((float) i - xx, j - yy, k - zz);
					kb = 0.0;
					//quasi radial...true cumulative radial weights are much more time consuming to code.
					for (int p = 0; p <= Util::get_min(r,(float) N/2); p++) {
						kb += ws[p];
					}
					dn = sqrt(kb/wm);
					size_t off = data->add_complex_at_fast(i,j,k,dt*dn*w);
					norm[off/2]+=w;
				}
			}
		}
		return true;
	}
	printf("region writing not supported in mode 10\n");
	return false;
}

const float FourierInserter3DMode11::kernel[12][12][12] = {
0.3756136,0.3403310,0.2474435,0.1296879,0.0245121,-0.0418448,-0.0631182,-0.0511850,-0.0262455,-0.0054893,0.0036277,0.0000000,
0.3403310,0.3082430,0.2237843,0.1167696,0.0212886,-0.0388006,-0.0578490,-0.0467286,-0.0238678,-0.0049178,0.0033657,0.0000000,
0.2474435,0.2237843,0.1615623,0.0828749,0.0129467,-0.0306454,-0.0438698,-0.0349361,-0.0175894,-0.0034180,0.0026655,0.0000000,
0.1296879,0.1167696,0.0828749,0.0402484,0.0028031,-0.0198868,-0.0258289,-0.0198082,-0.0095752,-0.0015321,0.0017471,0.0000000,
0.0245121,0.0212886,0.0129467,0.0028031,-0.0054675,-0.0095100,-0.0091333,-0.0059726,-0.0023192,0.0001220,0.0008690,0.0000000,
-0.0418448,-0.0388006,-0.0306454,-0.0198868,-0.0095100,-0.0018223,0.0022627,0.0032331,0.0023989,0.0011163,0.0002267,0.0000000,
-0.0631182,-0.0578490,-0.0438698,-0.0258289,-0.0091333,0.0022627,0.0071368,0.0068533,0.0041017,0.0013575,-0.0001089,0.0000000,
-0.0511850,-0.0467286,-0.0349361,-0.0198082,-0.0059726,0.0032331,0.0068533,0.0061573,0.0035159,0.0010550,-0.0001880,0.0000000,
-0.0262455,-0.0238678,-0.0175894,-0.0095752,-0.0023192,0.0023989,0.0041017,0.0035159,0.0019496,0.0005548,-0.0001273,0.0000000,
-0.0054893,-0.0049178,-0.0034180,-0.0015321,0.0001220,0.0011163,0.0013575,0.0010550,0.0005548,0.0001483,-0.0000403,0.0000000,
0.0036277,0.0033657,0.0026655,0.0017471,0.0008690,0.0002267,-0.0001089,-0.0001880,-0.0001273,-0.0000403,0.0000092,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.3403310,0.3082430,0.2237843,0.1167696,0.0212886,-0.0388006,-0.0578490,-0.0467286,-0.0238678,-0.0049178,0.0033657,0.0000000,
0.3082430,0.2790696,0.2022997,0.1050786,0.0184299,-0.0359603,-0.0530014,-0.0426434,-0.0216941,-0.0043990,0.0031230,0.0000000,
0.2237843,0.2022996,0.1458102,0.0744136,0.0110401,-0.0283535,-0.0401427,-0.0318353,-0.0159554,-0.0030383,0.0024746,0.0000000,
0.1167696,0.1050786,0.0744136,0.0358774,0.0020773,-0.0183257,-0.0235546,-0.0179765,-0.0086344,-0.0013296,0.0016238,0.0000000,
0.0212885,0.0184299,0.0110400,0.0020773,-0.0051866,-0.0086670,-0.0082158,-0.0053132,-0.0020141,0.0001649,0.0008100,0.0000000,
-0.0388006,-0.0359603,-0.0283535,-0.0183257,-0.0086670,-0.0015308,0.0022350,0.0030945,0.0022783,0.0010565,0.0002142,0.0000000,
-0.0578490,-0.0530014,-0.0401427,-0.0235546,-0.0082158,0.0022350,0.0066778,0.0063753,0.0038095,0.0012622,-0.0000978,0.0000000,
-0.0467286,-0.0426434,-0.0318353,-0.0179765,-0.0053133,0.0030945,0.0063753,0.0057001,0.0032491,0.0009749,-0.0001724,0.0000000,
-0.0238678,-0.0216940,-0.0159554,-0.0086344,-0.0020141,0.0022783,0.0038095,0.0032491,0.0017978,0.0005110,-0.0001172,0.0000000,
-0.0049178,-0.0043990,-0.0030383,-0.0013296,0.0001649,0.0010565,0.0012622,0.0009749,0.0005110,0.0001362,-0.0000371,0.0000000,
0.0033657,0.0031230,0.0024746,0.0016238,0.0008100,0.0002142,-0.0000978,-0.0001724,-0.0001172,-0.0000371,0.0000085,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.2474435,0.2237843,0.1615623,0.0828749,0.0129467,-0.0306454,-0.0438698,-0.0349361,-0.0175894,-0.0034180,0.0026655,0.0000000,
0.2237843,0.2022997,0.1458102,0.0744136,0.0110401,-0.0283535,-0.0401427,-0.0318353,-0.0159554,-0.0030383,0.0024746,0.0000000,
0.1615622,0.1458101,0.1044314,0.0522464,0.0061329,-0.0222221,-0.0302624,-0.0236374,-0.0116456,-0.0020445,0.0019643,0.0000000,
0.0828748,0.0744135,0.0522463,0.0244677,0.0002465,-0.0141589,-0.0175347,-0.0131427,-0.0061593,-0.0008027,0.0012941,0.0000000,
0.0129467,0.0110400,0.0061328,0.0002465,-0.0044008,-0.0064289,-0.0058001,-0.0035858,-0.0012204,0.0002714,0.0006519,0.0000000,
-0.0306454,-0.0283535,-0.0222221,-0.0141589,-0.0064289,-0.0007715,0.0021423,0.0027101,0.0019478,0.0008934,0.0001802,0.0000000,
-0.0438697,-0.0401427,-0.0302624,-0.0175348,-0.0058001,0.0021423,0.0054435,0.0050956,0.0030280,0.0010071,-0.0000689,0.0000000,
-0.0349361,-0.0318353,-0.0236374,-0.0131427,-0.0035858,0.0027101,0.0050956,0.0044802,0.0025380,0.0007613,-0.0001310,0.0000000,
-0.0175894,-0.0159554,-0.0116456,-0.0061593,-0.0012204,0.0019478,0.0030280,0.0025380,0.0013935,0.0003943,-0.0000903,0.0000000,
-0.0034180,-0.0030383,-0.0020445,-0.0008027,0.0002714,0.0008934,0.0010071,0.0007613,0.0003943,0.0001039,-0.0000288,0.0000000,
0.0026655,0.0024746,0.0019643,0.0012941,0.0006519,0.0001802,-0.0000689,-0.0001310,-0.0000903,-0.0000288,0.0000065,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.1296879,0.1167696,0.0828749,0.0402484,0.0028031,-0.0198868,-0.0258289,-0.0198082,-0.0095752,-0.0015321,0.0017471,0.0000000,
0.1167696,0.1050786,0.0744136,0.0358774,0.0020773,-0.0183257,-0.0235546,-0.0179765,-0.0086344,-0.0013296,0.0016238,0.0000000,
0.0828748,0.0744135,0.0522463,0.0244677,0.0002465,-0.0141589,-0.0175347,-0.0131427,-0.0061593,-0.0008027,0.0012941,0.0000000,
0.0402484,0.0358773,0.0244676,0.0102938,-0.0018371,-0.0087088,-0.0098086,-0.0069817,-0.0030272,-0.0001546,0.0008601,0.0000000,
0.0028031,0.0020773,0.0002465,-0.0018371,-0.0032679,-0.0035382,-0.0027388,-0.0014218,-0.0002427,0.0003867,0.0004426,0.0000000,
-0.0198868,-0.0183257,-0.0141589,-0.0087088,-0.0035382,0.0001652,0.0019638,0.0021631,0.0014894,0.0006692,0.0001335,0.0000000,
-0.0258289,-0.0235545,-0.0175347,-0.0098086,-0.0027388,0.0019638,0.0038001,0.0034089,0.0020002,0.0006705,-0.0000328,0.0000000,
-0.0198081,-0.0179764,-0.0131427,-0.0069817,-0.0014218,0.0021631,0.0034089,0.0028848,0.0016103,0.0004821,-0.0000781,0.0000000,
-0.0095751,-0.0086344,-0.0061593,-0.0030272,-0.0002427,0.0014894,0.0020002,0.0016103,0.0008674,0.0002423,-0.0000557,0.0000000,
-0.0015321,-0.0013296,-0.0008027,-0.0001546,0.0003867,0.0006692,0.0006705,0.0004821,0.0002423,0.0000617,-0.0000183,0.0000000,
0.0017471,0.0016238,0.0012941,0.0008601,0.0004426,0.0001335,-0.0000328,-0.0000781,-0.0000557,-0.0000183,0.0000038,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.0245121,0.0212886,0.0129467,0.0028031,-0.0054675,-0.0095100,-0.0091333,-0.0059726,-0.0023192,0.0001220,0.0008690,0.0000000,
0.0212885,0.0184299,0.0110400,0.0020773,-0.0051866,-0.0086670,-0.0082158,-0.0053132,-0.0020141,0.0001649,0.0008100,0.0000000,
0.0129467,0.0110400,0.0061328,0.0002465,-0.0044008,-0.0064289,-0.0058001,-0.0035858,-0.0012204,0.0002714,0.0006519,0.0000000,
0.0028031,0.0020773,0.0002465,-0.0018371,-0.0032679,-0.0035382,-0.0027388,-0.0014218,-0.0002427,0.0003867,0.0004426,0.0000000,
-0.0054675,-0.0051866,-0.0044008,-0.0032679,-0.0020109,-0.0008640,-0.0000115,0.0004597,0.0005759,0.0004525,0.0002390,0.0000000,
-0.0095100,-0.0086670,-0.0064289,-0.0035382,-0.0008640,0.0009487,0.0016881,0.0015623,0.0010061,0.0004364,0.0000852,0.0000000,
-0.0091333,-0.0082158,-0.0058001,-0.0027388,-0.0000115,0.0016881,0.0021848,0.0017821,0.0010130,0.0003454,-0.0000014,0.0000000,
-0.0059726,-0.0053132,-0.0035858,-0.0014218,0.0004596,0.0015623,0.0017821,0.0013691,0.0007327,0.0002175,-0.0000298,0.0000000,
-0.0023192,-0.0020141,-0.0012204,-0.0002427,0.0005759,0.0010061,0.0010130,0.0007327,0.0003725,0.0000992,-0.0000240,0.0000000,
0.0001220,0.0001649,0.0002714,0.0003867,0.0004525,0.0004364,0.0003454,0.0002175,0.0000992,0.0000219,-0.0000088,0.0000000,
0.0008690,0.0008101,0.0006519,0.0004426,0.0002390,0.0000852,-0.0000014,-0.0000298,-0.0000240,-0.0000088,0.0000011,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
-0.0418448,-0.0388006,-0.0306454,-0.0198868,-0.0095100,-0.0018223,0.0022627,0.0032331,0.0023989,0.0011163,0.0002267,0.0000000,
-0.0388006,-0.0359603,-0.0283535,-0.0183257,-0.0086670,-0.0015308,0.0022350,0.0030945,0.0022783,0.0010565,0.0002142,0.0000000,
-0.0306454,-0.0283535,-0.0222221,-0.0141589,-0.0064289,-0.0007715,0.0021423,0.0027101,0.0019478,0.0008934,0.0001802,0.0000000,
-0.0198868,-0.0183257,-0.0141589,-0.0087088,-0.0035382,0.0001652,0.0019638,0.0021631,0.0014894,0.0006692,0.0001335,0.0000000,
-0.0095100,-0.0086670,-0.0064289,-0.0035382,-0.0008640,0.0009487,0.0016881,0.0015623,0.0010061,0.0004364,0.0000852,0.0000000,
-0.0018223,-0.0015308,-0.0007715,0.0001652,0.0009487,0.0013519,0.0013308,0.0010089,0.0005874,0.0002396,0.0000447,0.0000000,
0.0022627,0.0022350,0.0021423,0.0019638,0.0016881,0.0013308,0.0009379,0.0005713,0.0002847,0.0001029,0.0000171,0.0000000,
0.0032331,0.0030945,0.0027101,0.0021631,0.0015623,0.0010089,0.0005713,0.0002745,0.0001051,0.0000275,0.0000025,0.0000000,
0.0023990,0.0022783,0.0019478,0.0014894,0.0010061,0.0005874,0.0002847,0.0001051,0.0000226,-0.0000019,-0.0000026,0.0000000,
0.0011163,0.0010565,0.0008934,0.0006692,0.0004364,0.0002396,0.0001029,0.0000275,-0.0000019,-0.0000062,-0.0000025,0.0000000,
0.0002267,0.0002142,0.0001802,0.0001335,0.0000852,0.0000447,0.0000171,0.0000025,-0.0000026,-0.0000025,-0.0000010,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
-0.0631182,-0.0578490,-0.0438698,-0.0258289,-0.0091333,0.0022627,0.0071368,0.0068533,0.0041017,0.0013575,-0.0001089,0.0000000,
-0.0578490,-0.0530014,-0.0401427,-0.0235546,-0.0082158,0.0022350,0.0066778,0.0063753,0.0038095,0.0012622,-0.0000978,0.0000000,
-0.0438697,-0.0401427,-0.0302624,-0.0175348,-0.0058001,0.0021423,0.0054435,0.0050956,0.0030280,0.0010071,-0.0000689,0.0000000,
-0.0258289,-0.0235545,-0.0175347,-0.0098086,-0.0027388,0.0019638,0.0038001,0.0034089,0.0020002,0.0006705,-0.0000328,0.0000000,
-0.0091333,-0.0082158,-0.0058001,-0.0027388,-0.0000115,0.0016881,0.0021848,0.0017821,0.0010130,0.0003454,-0.0000014,0.0000000,
0.0022627,0.0022350,0.0021423,0.0019638,0.0016881,0.0013308,0.0009379,0.0005713,0.0002847,0.0001029,0.0000171,0.0000000,
0.0071368,0.0066778,0.0054435,0.0038001,0.0021848,0.0009379,0.0002000,-0.0000882,-0.0001034,-0.0000291,0.0000213,0.0000000,
0.0068532,0.0063753,0.0050956,0.0034089,0.0017821,0.0005713,-0.0000882,-0.0002787,-0.0002039,-0.0000665,0.0000155,0.0000000,
0.0041017,0.0038095,0.0030280,0.0020002,0.0010130,0.0002847,-0.0001034,-0.0002039,-0.0001440,-0.0000496,0.0000066,0.0000000,
0.0013574,0.0012622,0.0010071,0.0006705,0.0003454,0.0001029,-0.0000291,-0.0000665,-0.0000496,-0.0000193,0.0000001,0.0000000,
-0.0001089,-0.0000978,-0.0000689,-0.0000328,-0.0000014,0.0000171,0.0000213,0.0000155,0.0000066,0.0000001,-0.0000022,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
-0.0511850,-0.0467286,-0.0349361,-0.0198082,-0.0059726,0.0032331,0.0068533,0.0061573,0.0035159,0.0010550,-0.0001880,0.0000000,
-0.0467286,-0.0426434,-0.0318353,-0.0179765,-0.0053133,0.0030945,0.0063753,0.0057001,0.0032491,0.0009749,-0.0001724,0.0000000,
-0.0349361,-0.0318353,-0.0236374,-0.0131427,-0.0035858,0.0027101,0.0050956,0.0044802,0.0025380,0.0007613,-0.0001310,0.0000000,
-0.0198081,-0.0179764,-0.0131427,-0.0069817,-0.0014218,0.0021631,0.0034089,0.0028848,0.0016103,0.0004821,-0.0000781,0.0000000,
-0.0059726,-0.0053132,-0.0035858,-0.0014218,0.0004596,0.0015623,0.0017821,0.0013691,0.0007327,0.0002175,-0.0000298,0.0000000,
0.0032331,0.0030945,0.0027101,0.0021631,0.0015623,0.0010089,0.0005713,0.0002745,0.0001051,0.0000275,0.0000025,0.0000000,
0.0068532,0.0063753,0.0050956,0.0034089,0.0017821,0.0005713,-0.0000882,-0.0002787,-0.0002039,-0.0000665,0.0000155,0.0000000,
0.0061573,0.0057001,0.0044802,0.0028848,0.0013691,0.0002745,-0.0002787,-0.0003841,-0.0002511,-0.0000812,0.0000141,0.0000000,
0.0035159,0.0032491,0.0025380,0.0016103,0.0007327,0.0001051,-0.0002039,-0.0002511,-0.0001613,-0.0000537,0.0000066,0.0000000,
0.0010550,0.0009749,0.0007613,0.0004821,0.0002175,0.0000275,-0.0000665,-0.0000812,-0.0000537,-0.0000199,0.0000000,0.0000000,
-0.0001880,-0.0001724,-0.0001310,-0.0000781,-0.0000298,0.0000025,0.0000155,0.0000141,0.0000066,0.0000000,-0.0000024,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
-0.0262455,-0.0238678,-0.0175894,-0.0095752,-0.0023192,0.0023989,0.0041017,0.0035159,0.0019496,0.0005548,-0.0001273,0.0000000,
-0.0238678,-0.0216940,-0.0159554,-0.0086344,-0.0020141,0.0022783,0.0038095,0.0032491,0.0017978,0.0005110,-0.0001172,0.0000000,
-0.0175894,-0.0159554,-0.0116456,-0.0061593,-0.0012204,0.0019478,0.0030280,0.0025380,0.0013935,0.0003943,-0.0000903,0.0000000,
-0.0095751,-0.0086344,-0.0061593,-0.0030272,-0.0002427,0.0014894,0.0020002,0.0016103,0.0008674,0.0002423,-0.0000557,0.0000000,
-0.0023192,-0.0020141,-0.0012204,-0.0002427,0.0005759,0.0010061,0.0010130,0.0007327,0.0003725,0.0000992,-0.0000240,0.0000000,
0.0023990,0.0022783,0.0019478,0.0014894,0.0010061,0.0005874,0.0002847,0.0001051,0.0000226,-0.0000019,-0.0000026,0.0000000,
0.0041017,0.0038095,0.0030280,0.0020002,0.0010130,0.0002847,-0.0001034,-0.0002039,-0.0001440,-0.0000496,0.0000066,0.0000000,
0.0035159,0.0032491,0.0025380,0.0016102,0.0007327,0.0001051,-0.0002039,-0.0002511,-0.0001613,-0.0000537,0.0000066,0.0000000,
0.0019496,0.0017978,0.0013935,0.0008674,0.0003725,0.0000226,-0.0001440,-0.0001613,-0.0001020,-0.0000349,0.0000026,0.0000000,
0.0005548,0.0005110,0.0003943,0.0002423,0.0000992,-0.0000019,-0.0000496,-0.0000537,-0.0000349,-0.0000136,-0.0000009,0.0000000,
-0.0001273,-0.0001172,-0.0000903,-0.0000557,-0.0000240,-0.0000026,0.0000066,0.0000066,0.0000026,-0.0000009,-0.0000020,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
-0.0054893,-0.0049178,-0.0034180,-0.0015321,0.0001220,0.0011163,0.0013575,0.0010550,0.0005548,0.0001483,-0.0000403,0.0000000,
-0.0049178,-0.0043990,-0.0030383,-0.0013296,0.0001649,0.0010565,0.0012622,0.0009749,0.0005110,0.0001362,-0.0000371,0.0000000,
-0.0034180,-0.0030383,-0.0020445,-0.0008027,0.0002714,0.0008934,0.0010071,0.0007613,0.0003943,0.0001039,-0.0000288,0.0000000,
-0.0015321,-0.0013296,-0.0008027,-0.0001546,0.0003867,0.0006692,0.0006705,0.0004821,0.0002423,0.0000617,-0.0000183,0.0000000,
0.0001220,0.0001649,0.0002714,0.0003867,0.0004525,0.0004364,0.0003454,0.0002175,0.0000992,0.0000219,-0.0000088,0.0000000,
0.0011163,0.0010565,0.0008934,0.0006692,0.0004364,0.0002396,0.0001029,0.0000275,-0.0000019,-0.0000062,-0.0000025,0.0000000,
0.0013574,0.0012622,0.0010071,0.0006705,0.0003454,0.0001029,-0.0000291,-0.0000665,-0.0000496,-0.0000193,0.0000001,0.0000000,
0.0010550,0.0009749,0.0007613,0.0004821,0.0002175,0.0000275,-0.0000665,-0.0000812,-0.0000537,-0.0000199,0.0000000,0.0000000,
0.0005548,0.0005110,0.0003943,0.0002423,0.0000992,-0.0000019,-0.0000496,-0.0000537,-0.0000349,-0.0000136,-0.0000009,0.0000000,
0.0001483,0.0001362,0.0001039,0.0000617,0.0000219,-0.0000062,-0.0000193,-0.0000199,-0.0000136,-0.0000063,-0.0000015,0.0000000,
-0.0000403,-0.0000371,-0.0000288,-0.0000183,-0.0000088,-0.0000025,0.0000001,0.0000000,-0.0000009,-0.0000015,-0.0000012,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.0036277,0.0033657,0.0026655,0.0017471,0.0008690,0.0002267,-0.0001089,-0.0001880,-0.0001273,-0.0000403,0.0000092,0.0000000,
0.0033657,0.0031230,0.0024746,0.0016238,0.0008100,0.0002142,-0.0000978,-0.0001724,-0.0001172,-0.0000371,0.0000085,0.0000000,
0.0026655,0.0024746,0.0019643,0.0012941,0.0006519,0.0001802,-0.0000689,-0.0001310,-0.0000903,-0.0000288,0.0000065,0.0000000,
0.0017471,0.0016238,0.0012941,0.0008601,0.0004426,0.0001335,-0.0000328,-0.0000781,-0.0000557,-0.0000183,0.0000038,0.0000000,
0.0008690,0.0008101,0.0006519,0.0004426,0.0002390,0.0000852,-0.0000014,-0.0000298,-0.0000240,-0.0000088,0.0000011,0.0000000,
0.0002267,0.0002142,0.0001802,0.0001335,0.0000852,0.0000447,0.0000171,0.0000025,-0.0000026,-0.0000025,-0.0000010,0.0000000,
-0.0001089,-0.0000978,-0.0000689,-0.0000328,-0.0000014,0.0000171,0.0000213,0.0000155,0.0000066,0.0000001,-0.0000022,0.0000000,
-0.0001880,-0.0001724,-0.0001310,-0.0000781,-0.0000298,0.0000025,0.0000155,0.0000141,0.0000066,0.0000000,-0.0000024,0.0000000,
-0.0001273,-0.0001172,-0.0000903,-0.0000557,-0.0000240,-0.0000026,0.0000066,0.0000066,0.0000026,-0.0000009,-0.0000020,0.0000000,
-0.0000403,-0.0000371,-0.0000288,-0.0000183,-0.0000088,-0.0000025,0.0000001,0.0000000,-0.0000009,-0.0000015,-0.0000012,0.0000000,
0.0000092,0.0000085,0.0000065,0.0000038,0.0000011,-0.0000010,-0.0000022,-0.0000024,-0.0000020,-0.0000012,-0.0000005,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
};

bool FourierInserter3DMode11::insert_pixel(const float& xx, const float& yy, const float& zz, const std::complex<float> dt,const float& weight)
{
	int x0 = Util::fast_floor(xx-2.5);
	int y0 = Util::fast_floor(yy-2.5);
	int z0 = Util::fast_floor(zz-2.5);

	if (subx0<0) {			// normal full reconstruction
		if (x0<-nx2-6 || y0<-ny2-6 || z0<-nz2-6 || x0>nx2+5 || y0>ny2+5 || z0>nz2+5 ) return false;

		// no error checking on add_complex_fast, so we need to be careful here
		int x1=x0+6;
		int y1=y0+6;
		int z1=z0+6;
		if (x0<-nx2) x0=-nx2;
		if (x1>nx2) x1=nx2;
		if (y0<-ny2) y0=-ny2;
		if (y1>ny2) y1=ny2;
		if (z0<-nz2) z0=-nz2;
		if (z1>nz2) z1=nz2;

		float w=weight;

		float r, gg;
		for (int k = z0 ; k <= z1; k++) {
			for (int j = y0 ; j <= y1; j++) {
				for (int i = x0; i <= x1; i ++) {
					gg=FourierInserter3DMode11::kernel[abs(Util::fast_floor((i-xx)*3.0f+0.5))][abs(Util::fast_floor((j-yy)*3.0f+0.5))][abs(Util::fast_floor((k-zz)*3.0f+0.5))]; 

					size_t off;
					off=data->add_complex_at_fast(i,j,k,dt*gg*w);
					norm[off/2]+=w;
//					if (i==67&&j==19&&k==1) printf("%1.1f  %1.1f  %1.1f\t%d %d %d\t%d %d %d\t%f\t%f\t%f\t%f\t%f\t%f\n",xx,yy,zz,i,j,k,abs(Util::fast_floor((i-xx)*3.0f+0.5)),abs(Util::fast_floor((j-yy)*3.0f+0.5)),abs(Util::fast_floor((k-zz)*3.0f+0.5)),gg,w,dt.real(),dt.imag(),norm[off/2],data->get_value_at(off));

				}
			}
		}
		return true;
	}
	printf("region writing not supported in mode \n");
	return false;
}
