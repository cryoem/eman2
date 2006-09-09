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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include "emdata.h"

using namespace EMAN;

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

		//vector<int> saved_offsets = fp->get_array_offsets();  I do not know what the meaning of it was, did not work anyway PAP
		fp->set_array_offsets(1,1,1);

		normfact = (nxp/nx)*(nyp/ny)*(nzp/nz);  // Normalization factor for the padded operations
		if(normfact>1) {
		 for (int iz = 1; iz <= nz; iz++) for (int iy = 1; iy <= ny; iy++) for (int ix = 1; ix <= nx; ix++) (*fp)(ix,iy,iz) *= normfact;
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
		//fp->set_array_offsets(saved_offsets);  This was strange and did not work, PAP
		fp->set_array_offsets(0,0,0);
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

#ifdef NEVER
 Voronoi code from Phani

test/voronoi/VORONOI/C/voronoi/final.cpp
gcc -g -c final.cpp
this will call as a binary genwtsc.c compiled and exe called getwts
in directory test/voronoi/VORONOI/C
gcc -g -ansi -I../F2C/libf2c/ -c stridpack.c
gcc -g -ansi genwts.c stridpack.o -L../F2C/libf2c/ -lf2c -lm -o getwts



********************final.cpp
#include <numeric>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <stdio.h> 
#include <vector> 
#include <math.h>

#include <cstdlib>

 
using namespace std;


#define		quadpi 	 	3.141592653589793238462643383279502884197
#define 	twopi 		2*quadpi
#define		dgr_to_rad	quadpi/180
#define		rad_to_deg	180/quadpi

struct	tmpstruct{
 	double theta1,phi1;
	int key1;
 	};

bool cmp1(tmpstruct tmp1,tmpstruct tmp2)
{
 return(tmp1.theta1 < tmp2.theta1);
}

bool cmp2(tmpstruct tmp1,tmpstruct tmp2)
{
 return(tmp1.key1 < tmp2.key1);
}

void hsortd(vector <double> &theta,vector <double> &phi,vector <int> &key,int option)
{
vector<tmpstruct> tmp(theta.size());
for(int i = 0;i<theta.size();i++)
	{
		tmp[i].theta1 = theta[i];
		tmp[i].phi1 = phi[i];
		tmp[i].key1 = key[i];
	}
if (option == 1) sort(tmp.begin(),tmp.end(),cmp1);
if (option == 2) sort(tmp.begin(),tmp.end(),cmp2);

for(int i = 0;i<theta.size();i++)
	{
		theta[i] = tmp[i].theta1;
		phi[i] 	 = tmp[i].phi1;
		key[i]   = tmp[i].key1;
	}
}



vector <double> xyz_to_ang(vector<double> x,vector<double> y,vector<double> z)
{
 int len = x.size();
 
 vector<double> theta(len,0),phi(len,0);
 double temp;
 //FILE *fp = fopen("phi.txt","w");
 for(int i = 0;i<len;i++)
    {
     if (fabs(z[i])>0.99999999) continue;
     theta[i]  = acos(z[i]);
     temp       = sin(theta[i]); temp = 1/temp;
     phi[i]  = atan2(y[i]*temp,x[i]*temp);
     theta[i] *= rad_to_deg;
     phi[i]   *= rad_to_deg; 
     if (theta[i]>90)
         {
	  theta[i] = 180 - theta[i];
	  phi[i] += 180;
	 }
     phi[i] = fmod(phi[i]+360,360.0);
   //  fprintf(fp,"%d  %lf  %lf\n",i+1,phi[i],theta[i]/rad_to_deg);
    }
    
 vector <double> angles;
 angles.assign(theta.begin(),theta.end());
 angles.insert(angles.end(),phi.begin(),phi.end());        
 
// fclose(fp);
 
 return (angles); 
}




void ang_to_xyz(vector<double> key_theta_phi,double *x,double *y,double *z,int *key)
{
 // x is theta; Y is phi;
 int len = key_theta_phi.size()/3;
 
 x = new double[len];
 y = new double[len];
 z = new double[len];
 key = new int[len]; 
 
  
 double costheta,sintheta,cosphi,sinphi;
 for(int i = 0;i<len;i++) 
 	{
		key[i] = (int)key_theta_phi[i];
		cosphi = cos(key_theta_phi[(2*len)+i]*dgr_to_rad);
		sinphi = sin(key_theta_phi[(2*len)+i]*dgr_to_rad);
		if(fabs(key_theta_phi[len+i]-90)< 1.0e-5){
			x[i] = cosphi;
			y[i] = sinphi;
			z[i] = 0.0;
		}
		else{
			costheta = cos(key_theta_phi[len+i]*dgr_to_rad);
			sintheta = sin(key_theta_phi[len+i]*dgr_to_rad);
			x[i] = cosphi*sintheta;
			y[i] = sinphi*sintheta;
			z[i] = costheta;
		}
	}	 
}
//-----------------------------------------------------------------------------
vector<double> disorder2(vector<double> theta,vector<double> phi)
{
  int len = phi.size();
  vector <double> key(len,1.0);
  for(int i = 1;i<=len;i++)
  	key[i] = i+(0.0);
 
  std::random_shuffle(key.begin(),key.end());
  double temp;
  for(int i = 1;i<=len;i++)
  	{
		std::swap(phi[i],phi[(int)key[i]]);
		std::swap(theta[i],theta[(int)key[i]]);
 	}

  vector<double> key_theta_phi;
  key_theta_phi.assign(key.begin(),key.end());
  key_theta_phi.insert(key_theta_phi.end(),theta.begin(),theta.end());
  key_theta_phi.insert(key_theta_phi.end(),phi.begin(),phi.end());

  return (key_theta_phi); 
}
//-----------------------------------------------------------------------------
void flip23(double *x,double *y,double *z,int *key,int k,int len)
{
 int i = k;
 while(i==k)
 	i = rand()%len;
 
 std::swap(key[i],key[k]);
 std::swap(x[i],x[k]);
 std::swap(y[i],y[k]);
 std::swap(z[i],z[k]); 
} 
//-----------------------------------------------------------------------------


void vrdg(double *dm3,int n,int nang2)
{ 
 vector <double> dm,x,y,z;
 vector <int> key;
 double tmp[3];
 FILE *fp;
 
 for(int i =0;i<nang2*n;i++)
 { 
  x.push_back(0.0);
  y.push_back(0.0);
  z.push_back(0.0);
  key.push_back(0);
  if(i%n ==0)
  	{dm.push_back(0.0);dm.push_back(0.0);}	
 }
 // angles for reverse gridding
double  phis=quadpi/nang2;
 
//fp = fopen("dmc.txt","w"); 
 for (int i=1;i<=nang2;i++){ 
	dm[2*i -2]	=	cos( phis*(i-1) );
	dm[2*i -1]	=	sin( phis*(i-1) );
//	fprintf(fp,"%lf	%lf\n",dm[2*i-2],dm[2*i-1]); 
//	cout<<dm[2*i-2]<<"  "<<dm[2*i-1]<<endl;
 }	
  
// fclose(fp);
// fp = fopen("XYZ_c.txt","w");
 for(int k=1;k<=n;k++){
 	for(int in=1;in<=nang2;in++)
		{
			x[(in-1)+(k-1)*nang2] = dm[2*in - 2] * dm3[9*k - 9] + dm[2*in - 1] * dm3[9*k - 6];
			y[(in-1)+(k-1)*nang2] = dm[2*in - 2] * dm3[9*k - 8] + dm[2*in - 1] * dm3[9*k - 5];
			z[(in-1)+(k-1)*nang2] = dm[2*in - 2] * dm3[9*k - 7] + dm[2*in - 1] * dm3[9*k - 4];
//			fprintf(fp,"%d	%20.18lf	%20.18lf	%20.18lf\n",(in-1)+(k-1)*nang2,x[(in-1)+(k-1)*nang2],y[(in-1)+(k-1)*nang2],z[(in-1)+(k-1)*nang2]); 
		}
	}
// fclose(fp);
 
 vector <double> ang = xyz_to_ang(x,y,z);
 x.clear();
 y.clear();
 z.clear();
 
// fp = fopen("thetaphi_c.txt","w");
 vector <double> phi,theta;
 for(int i = 0;i<ang.size()/2;i++){
 	theta.push_back(ang[i]);
	phi.push_back(ang[n*nang2 + i]);
	//cout<<theta[i]<<"  "<<phi[i]<<endl;
//	fprintf(fp,"%d	%30.28lf	%30.28lf\n",i,theta[i],phi[i]); 	
 }
// fclose(fp);
 
 ang.clear();   

  
 //  sort by theta

 for(int i=1;i<=n*nang2;i++)
 	key[i-1]=i;

 hsortd(theta,phi,key,1);

 
 fp = fopen("theta_phi.txt","w");  
 for(int i = 1;i<=nang2*n;i++){	
	fprintf(fp,"%30.28lf   %30.28lf\n",theta[i-1],phi[i-1]);
}
 fclose(fp);
 
 fp = fopen("parameters.txt","w");
 //fprintf(fp,"%d	%d	%d	%d",n*nang2,1,n*nang2 +1,n*nang2);
 fprintf(fp,"%d	%d\n",nang2,n); 
 fclose(fp);
 
 
// cout<<" at line 228"<<endl;	
  
int ik = system ("../getwts");
 	//vector<double> weight = voronoi(theta,phi,1,n*nang2+1,n*nang2,true);
	 //voronoidiag(theta,phi); 
 ik = remove("theta_phi.txt");
 ik = remove("parameters.txt"); 
  
 fp = fopen("myweights.txt","r");
 vector<double> weight(key.size(),0.0); 
 int temp1;
 double temp2;
 for(int i = 1;i<=nang2*n;i++){
 	fscanf(fp,"%lf\n",&temp2);
	weight[i-1] = temp2;
}
 fclose(fp);
 ik = remove("myweights.txt");
  
 hsortd(weight,weight,key,2);

// cout<<" at line 250"<<endl;


 fp = fopen("/home/phani/test/grj15/finalweights.txt","w");
 for(int k = 1,count = 0;k<= n;k++){
 	for(int in = 1;in<=nang2;in++){
		fprintf(fp,"%d  %30.28lf\n",count+1,weight[count]);
		count++;
	}
 }	
 fclose(fp); 
}	 

//------------------------------------------------------------------------------------------------------------

int main(){
	FILE *fp;
	int n,nang2;
	nang2 = 236;n = 84;

	double *dm3 = new double[n*9];
	fp = fopen("fort.67","r");
	for(int i = 0;i<9*n;i++){
		if(feof(fp)){ break;cout<<i+1- 9*n;}
		fscanf(fp,"%lf",dm3+i);		
	}
	fclose(fp);
     	//cout<<" at line 160"<<endl;	
	vrdg(&dm3[0],n,nang2); 
	delete []dm3;
	return 0;
}


 
************************   genwtsc.c

#include<stdio.h>
/*#include<stdlib.h>*/
#include<math.h> 
#include<f2c.h>


#define		mymax(x,y)		(((x)>(y))?(x):(y))
#define 	mymin(x,y)		(((x)<(y))?(x):(y))
#define		sign(x,y)		(((((y)>0)?(1):(-1))*(y!=0))*(x))
  

#define		quadpi	 	 	3.141592653589793238462643383279502884197
#define 	twopi 			2*quadpi
#define		dgr_to_rad		quadpi/180
#define		deg_to_rad		quadpi/180
#define		rad_to_deg		180/quadpi
#define		rad_to_dgr		180/quadpi
#define		TRUE			1
#define		FALSE			0





#define theta(i)		theta	[i-1]
#define phi(i)			phi	[i-1]
#define weight(i)		weight	[i-1]
#define lband(i)		lband	[i-1]
#define aat(i)			aat	[i-1]
#define ts(i)			ts	[i-1]
#define	thetast(i)		thetast	[i-1]  



int  trmsh3_(int    *, double *, double *, double *, 
             double *, int    *, int    *, int    *,
             int    *, int    *, int    *, int    *,
             int    *, int    *, double *, int    *);
double areav_(int    *, int *, double *, double *, 
              double *, int *, int    *, int    *, 
              int *);
	      

void voronoidiag(double *theta,double *phi,double* weight,int nang2,int nn);
void voronoi(double *theta,double *phi,double *weight,int lenw,int low,int medium,int nt,int last);	      
void angstep(double* thetast,int len);
void flip23(double *x,double *y,double *z,int *key,int k,int len);
void ang_to_xyz(double *x,double *y,double *z,int len);
void disorder2(double *x,double *y,int *key,int len);
void swap(double *x,double *y);
double dot_product(double *x,double *y);
 


int voronoidiag_(doublereal *theta, doublereal *phi,doublereal *weight, integer *nang2, integer *nn);
int angstep_(doublereal*,integer* );
 





		/*************************************************/
		/*                   Main			 */
		/*************************************************/

int main(){

	int nang2,n,i;
	double *theta,*phi,*weight;
	
	FILE *fp;
	fp = fopen("/home/phani/test/voronoi/VORONOI/C/voronoi/parameters.txt","r");
 	fscanf(fp,"%d	%d\n",&nang2,&n);
	fclose(fp);
		 
	weight  = 	(double*) calloc(n*nang2,sizeof(double));
	theta  = 	(double*) calloc(n*nang2,sizeof(double));
	phi  = 	(double*) calloc(n*nang2,sizeof(double));


	if ((theta==NULL) || (phi==NULL) || (weight==NULL)){
	       fprintf(stderr,"memory allocation failure!\n");
	       exit(1);
	}
	
	fp = fopen("/home/phani/test/voronoi/VORONOI/C/voronoi/theta_phi.txt","r");
 	for(i =1;i<=n*nang2;i++){
		fscanf(fp, "%lf %lf\n", theta+i-1, phi+i-1);
		weight[i-1]= 0.0;
	}
	fclose(fp);
 	
	voronoidiag(theta,phi,weight,nang2,n);
	/*i = voronoidiag_(theta, phi, weight, &nang2, &n);*/
	
	/*voronoi(theta,phi,weight,n*nang2,1,(n*nang2)+1,n*nang2,1);*/
	 
	
	fp = fopen("/home/phani/test/voronoi/VORONOI/C/voronoi/myweights.txt","w");
	for(i=1;i<=n*nang2;i++)
		fprintf(fp,"%30.28lf\n",weight[i-1]);
	fclose(fp);

	free(weight);
	free(theta);
	free(phi);
	
	printf("\n out of getwts \n");
	
	return(0);
}





		/*************************************************/
		/*************************************************/

/***************************************************************************************************/ 

void voronoidiag(double *theta,double *phi,double* weight,int nang2,int nn){
	
	int	*lband;
	double	aat,*ts; 
	double	aa,acum,area;
	int	last;
  	
	int n		=	nang2*nn;
	int numth 	= 	1;/*call getthreads(numth)*/ 
	int nbt		=	mymax((int)(sqrt((n/500.0))) , 3);
	
	int i,it,l,k;
	int nband,lb,low,medium,lhigh,lbw,lenw;
	int count;
	char ch;
/**********************************************************************/

	printf(" in voronoidiag");
	
	lband 	=	(int*)calloc(nbt,sizeof(int));
	ts 	=	(double*)calloc(nbt,sizeof(double));
	/*aat 	=	(double*)calloc(numth,sizeof(int));*/
	
	printf(" in voronoidiag");
	if(lband == NULL || ts == NULL ){
       		fprintf(stderr,"memory allocation failure!\n");
		exit(1);
	}
/*	printf("\n entering main loop in voronoidiag");*/
	nband=nbt;
	count = 0;
	printf(" in voronoidiag, %d",nbt);
	
	
	while(nband>0){
					
		angstep(ts,nband);
		for(i=1;i<=nband;i++)
			printf(" %lf ",ts[i-1]);
		
 
		l=1;
		for(i=1;i<=n;i++){
			if(theta(i)>ts(l)){
				lband(l)=i;
				l=l+1;
				if(l>nband)  exit(1);
			}
		}
		
			
		printf("\n @ line 155  %d ",l);

		lband(l)=n+1;
 		acum=0.0;

		for(it=l;it>=1;it-=numth){
			/*printf("\n @ line 161");*/
			 
			count++;
			if (count>10)	goto label3;
/*				printf("\n count - %d, itertion  - %d",count,i);*/
			if(i==l) last	=	 TRUE;
			else	 last	=	 FALSE;

			if(l==1){
				lb=1;
				low=1;
				medium=n+1;
				lhigh=n-lb+1;
				lbw=1;					
			} 				
			else if(i==1){
				lb=1;
				low=1;
				medium=lband(1);
				lhigh=lband(2)-1;
				lbw=1;
			}
			else if(i==l){
				if(l==2)	lb=1;
 				else		lb=lband(l-2);
				low=lband(l-1)-lb+1;
				medium=lband(l)-lb+1;
				lhigh=n-lb+1;
				lbw=lband(i-1);
			}
			else{
				if(i==2)	lb=1;
				else		lb=lband(i-2);
				low=lband(i-1)-lb+1;
				medium=lband(i)-lb+1;
				lhigh=lband(i+1)-1-lb+1;
				lbw=lband(i-1);
			}
			lenw=medium-low;
			
			printf(" /n lb - %d , lbw - %d, low - %d, medium  - %d, lhigh - %d",lb,lbw,low,medium,lhigh);
			goto label4;
/*			lb = 1;lbw = 1;lenw = n;low = 1;medium = n+1;lhigh = n; last = 1;*/ 
 			voronoi(&phi(lb),&theta(lb),&weight(lbw),lenw,low,medium,lhigh,last); 
			
/*c  test the areas in all sections*/

			printf("\n %d  nband",nband);
			if(nband>1){ 
				if(i==1)	area=quadpi*2.0*(1.0-cos(ts(1)*dgr_to_rad));
				else		area=quadpi*2.0*(cos(ts(i-1)*dgr_to_rad)-cos(ts(i)*dgr_to_rad));

				aa = 0.0;
				for(k = lbw;k<=lbw+lenw-1;k++)
					aa = aa+weight(k);
					/*aa=sum(weight(lbw:lbw+lenw-1));*/

				acum=acum+aa;
		   		aat=aa/area;
			}


			if(fabs(aat-1.0)>0.02){
				nband=mymax(0,mymin( (int)(((float)nband) * 0.75) ,nband-1) );
				goto  label2;
			}
		}
		acum=acum/quadpi/2.0; 
		exit(1);
label2:		
	 	 
		continue;
	}/*while*/
label3:
	
	/*free(aat);*/
	free(ts);
	free(lband);
label4:
i = 2;

}


/***************************************************************************************************/

void voronoi(double *theta,double *phi,double *weight,int lenw,int low,int medium,int nt,int last){
/*
lenw 	= 	n*nang2;	//length of weight
low  	= 	1;
nt 	= 	n*nang2;
medium	= 	n*nang2 +1;
last	= 	1;//TRUE
*/
/*phi(nt),theta(nt),weight(lenw)*/
/* list(nt6),lptr(nt6),lend(n),iwk(n),good(n),key(n),indx(n),lcnt(n)*/
  	
	
	
	FILE *fp;
     	double areav;
      	/*double  *phi, *theta,*weight;*/
      	int *list, *lptr, *lend, *iwk, *key,*lcnt,*indx,*good;
      
	int len;
	int nt6, n, ier,nout,lnew,mdup,nd;
	int i,k,mt,status;
	

     	double *ds,*x,*y,*z;
	double tol=1.0e-8;
	double tmp1[3],tmp2[3],a;

      	/*ds(n),x(n),y(n),z(n)*/
      	/*logical  last*/
	/*
	fp = fopen("parameters.txt","r");
	fscanf(fp,"%d	%d	%d	%d",&lenw,&low,&medium,&nt);
	fclose(fp);
	*/
	
	if(last){
		if(medium>nt)  n=nt+nt;
		else   n=nt+nt-medium+1;
	} 
	else{
		n=nt;
	}
	 
	nt6 = n*6;
 /*	printf("\n in voronoi \n %d  %d, %d  %d,%d",lenw,low,medium,nt,last);		*/

	list = (int*)calloc(nt6,sizeof(int));  
	lptr = (int*)calloc(nt6,sizeof(int));  
	lend = (int*)calloc(n  ,sizeof(int));  
	iwk  = (int*)calloc(n  ,sizeof(int));  
	good = (int*)calloc(n  ,sizeof(int));  
	key  = (int*)calloc(n  ,sizeof(int));  
	indx = (int*)calloc(n  ,sizeof(int)); 
	lcnt = (int*)calloc(n  ,sizeof(int)); 
	
	ds	= 	(double*) calloc(n,sizeof(double)); 
	x	= 	(double*) calloc(n,sizeof(double)); 
	y	= 	(double*) calloc(n,sizeof(double)); 
	z	= 	(double*) calloc(n,sizeof(double)); 
	/*weight  = 	(double*) calloc(lenw,sizeof(double)); */
	
	/*printf("\n @line 111");*/
	 
       if (list == NULL ||
       lptr == NULL ||
       lend == NULL ||
       iwk  == NULL ||
       good == NULL ||
       key  == NULL ||
       indx == NULL ||
       lcnt == NULL ||
       x    == NULL ||
       y    == NULL ||
       z    == NULL ||
       ds   == NULL) {
       /*printf("memory allocation failure!\n");*/
       fprintf(stderr,"memory allocation failure!\n");
       exit(1);
   }
	
/*	printf("\n @line 111");*/
	
/* 	fp = fopen("theta_phi.txt","r");
 	for(i =1;i<=nt;i++)
		fscanf(fp, "%lf %lf\n", x+i-1, y+i-1);
	fclose(fp);
*/
	
	
	for(i = 1;i<=nt;i++){
		x[i-1] = theta[i-1];
		y[i-1] = phi[i-1];
	}
/*	printf("\n @line 342");*/
	if (last) { 
		for(i=nt+1;i<=n;i++){			
			x[i-1]=180.0-x[2*nt-i];
			y[i-1]=180.0+y[2*nt-i];
		}
	}
	
 	/*printf("\n @line 350");*/
	disorder2(x,y,key,n);       
	ang_to_xyz(x,y,z,n);        
label1: 
	for(k=1;k<=2;k++){
		for(i=k+1;i<=3;i++){
			tmp1[0]= x[i-1]-x[k-1];tmp1[1]= y[i-1]-y[k-1];tmp1[2] = z[i-1]-z[k-1];
			/*tmp2[0]= x[k-1];tmp2[1]= y[k-1];tmp2[2] = z[k-1];
			tmp2[0] = dot_product(tmp1,tmp1);*/
			if( dot_product(tmp1,tmp1) < tol){
				flip23(x,y,z,key,k-1,n);
				goto label1;
		      	} 
		}
	}
  
/*	printf("\n @line 366");*/
        status = trmsh3_(&n,&tol,x,y,z,&nout,list,lptr,lend,&lnew,indx,lcnt,iwk,good,ds,&ier);
	                  
	if (status != 0) {
		printf(" error in trmsh3\n");
		/*return status;*/
	}
	
	printf("nout = %d\n", nout);
	
	mdup=n-nout;
	printf("number of duplicated=%d\n", mdup);
	if (ier == -2) {
		printf("*** Error in TRMESH:the first three nodes are collinear***\n");
		exit(1);
	}
	else if (ier > 0) {
		printf("*** Error in TRMESH:  duplicate nodes encountered ***\n");
		exit(1);
	}
 	
 	nd=0;
	for (k=1;k<=n;k++){
		if (indx[k-1]>0){
			nd++;
			good[nd-1]	=	k;
		}
	}
	 
	 
	 
     	for(i = 1;i<=nout;i++) {
		k=good[i-1];
 		if (key[k-1] >= low && key[k-1]<medium){
			a = areav_(&i,&nout,x,y,z,list,lptr,lend,&ier);
    			if (ier != 0){
				printf("    *** error in areav:  ier = %d ***\n", ier);
				weight[key[k-1]-low] =-1.0;
			}
			else {
				weight[key[k-1]-low]	=	a/lcnt[i-1];
			}
		}
	}
	 		
		
	for(i = 1;i<=n;i++){
		mt=-indx[i-1];
		if (mt>0){
			k=good[mt-1];
			if ( key[i-1]>=low && key[i-1]<medium){
				if(key[k-1]>=low && key[k-1]<medium){
					weight[key[i-1]-low]=weight[key[k-1]-low];
				}
				else{
					a = areav_(&mt,&nout,x,y,z,list,lptr,lend,&ier);
					if (ier != 0){
						printf("    *** error in areav:  ier = %d ***\n", ier);
						weight[key[i-1]-low] =-1.0;
					}
					else {
						weight[key[i-1]-low]	=	a/lcnt[i-1];
					}
				}
			}
		}
      	}
 	/* 
	 
	fp = fopen("myweights.txt","w");
	for(i=1;i<=lenw;i++)
		fprintf(fp,"%d  %30.28lf\n",i,weight[i-1]);
	fclose(fp);
	*/
	
			
 	free(list);
	free(lend);
	free(iwk);
	free(good);
	free(key); 

	free(indx);
	free(lcnt);
	free(ds);
	free(x);
	free(y);
	free(z);
/*	printf("\n exit from voronoi @456");*/
	/*free(weight);	*/
}

/***************************************************************************************************/

void ang_to_xyz(double *x,double *y,double *z,int len){
 /* x is theta; Y is phi;*/
   
 double costheta,sintheta,cosphi,sinphi;
 int i;
 for(i = 0;i<len;i++) 
 	{
		cosphi = cos(y[i]*dgr_to_rad);
		sinphi = sin(y[i]*dgr_to_rad);
		if(fabs(x[i]-90)< 1.0e-5){
			x[i] = cosphi;
			y[i] = sinphi;
			z[i] = 0.0;
		}
		else{
			costheta = cos(x[i]*dgr_to_rad);
			sintheta = sin(x[i]*dgr_to_rad);
			x[i] = cosphi*sintheta;
			y[i] = sinphi*sintheta;
			z[i] = costheta;
		}
	}	 
}

/***************************************************************************************************/


void disorder2(double *x,double *y,int *key,int len){
	
	int k,it,i;
	double tmp;
	for(i=0;i<len;i++)
		key[i]	=	i+1;
	
	
	for(i = 0;i<len;i++){
		k = rand()%len;
		it=key[k];
		key[k]=key[i];
		key[i]=it;
		
		tmp=x[k];x[k]=x[i];x[i] = tmp; 
		tmp=y[k];y[k]=y[i];y[i] = tmp; 
		/*tmp=z[k];z[k]=z[i];z[i] = tmp;*/
	}
}

/***************************************************************************************************/
void flip23(double *x,double *y,double *z,int *key,int k,int len)
{
 int i = k;
 while(i==k)
 	i = rand()%len;
 
 int tmp = key[i];key[i] = key[k];key[k]  = tmp;
 swap(x+i,x+k);
 swap(y+i,y+k);
 swap(z+i,z+k); 
} 


/***************************************************************************************************/
void swap(double *x,double *y)
{
 double tmp;
 tmp = *x;
 *x = *y;
 *y = tmp;
}

double dot_product(double *x,double *y)
{
 return(x[0]*y[0] +x[1]*y[1] +x[2]*y[2]);
}
/***************************************************************************************************/ 
 
void angstep(double* thetast,int len){
	double t1,t2,tmp;
	int i;
 	if(len>1){
		t1=0;
		for(i=1;i<=len-1;i++){
			tmp=cos(t1)-1/((float)len);	
			t2=acos(sign(mymin(1.0,fabs(tmp)),tmp));
			thetast(i)=t2 * rad_to_deg;
			t1=t2;
		} 
	}
	thetast(len)=90.0;
}


/***************************************************************************************************/


/* Subroutine */ int voronoidiag_(double *theta, doublereal *phi, 
	double *weight, integer *nang2, integer *nn)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
    real r__1;
    doublereal d__1;

    /* Builtin functions */
	/*    double sqrt(double);*/
    integer i_nint(real *);
    /* Subroutine */ int s_stop(char *, ftnlen);
    /*double cos(double);*/

    /* Local variables */
    static integer i__, l, n;
    static double aa;
    static integer lb, it;
    static double *ts;
    static integer nbc,angtmp;
    static double *aat;
    static integer nbt, lbw, low;
    static double area;
    static double acum;
    static logical last;
    static integer lenw, *lband, nband, lhigh, numth, medium;

/* 	INTEGER, DIMENSION(:), ALLOCATABLE :: LBAND */
/* 	REAL, DIMENSION(:), ALLOCATABLE :: TS,AAT */

/* Logical unit number and output file name: */
    n = *nang2 * *nn;
    
 	       
    /*Parameter adjustments */
    --weight;
    --phi;
    --theta;

 
    /* Function Body */
    n = *nang2 * *nn;
    numth = 1;/* 	CALL GETTHREADS(NUMTH) */
    

 /* Computing MAX */
    r__1 = sqrt(((double)n /500.0));
    i__1 = i_nint(&r__1);
    nbt = max(i__1,3);
/*    nbt	= 1;*/
    
    
    
    lband  = 	(int*) calloc(nbt,sizeof(int));
    aat	   = 	(double*) calloc(numth,sizeof(double));
    ts	   = 	(double*) calloc(nbt,sizeof(double));
    
    
/* ********************************************************************** */
/* 	ALLOCATE(AAT(NUMTH),TS(NBT),LBAND(NBT),STAT=ISTAT) */
/* 	IF (ISTAT.NE.0)  THEN */
/* 	  PRINT  *,' ALLOCATION FAILED, NBAND' */
/* 	  STOP */
/* 	ENDIF */
    nband = nbt;
    while(nband > 0) {
/* 	PRINT  *,' NBAND=',NBAND */
 	printf("\n nband -%d",nband);

/*  TS contains angles selected such that the number of nodes in each band should */
/*    be the same */
/*  TS(1) - limit of the first band, contains angles from 0 to TS(1) */
/*  TS(NBAND) - is always 90 */

 	angstep(ts,nband);

	l = 1;
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (theta[i__] > ts[l - 1]) {
		lband[l - 1] = i__;
		++l;
		if (l > nband) {
		    s_stop("", (ftnlen)0);
		}
	    }
	}
	lband[l - 1] = n + 1;
/* 	PRINT  *,(LBAND(I),I=1,L) */
 
	acum = 0.f;
 
	i__1 = -numth;
	for (it = l; i__1 < 0 ? it >= 1 : it <= 1; it += i__1) {
 	    i__3 = 1, i__4 = it - numth + 1;
	    i__2 = max(i__3,i__4);
	    for (i__ = it; i__ >= i__2; --i__) {
		if (i__ == l) {
		    last = TRUE_;
		} else {
		    last = FALSE_;
		}
		if (l == 1) {
		    lb = 1;
		    low = 1;
		    medium = n + 1;
		    lhigh = n - lb + 1;
		    lbw = 1;
		} else if (i__ == 1) {
		    lb = 1;
		    low = 1;
		    medium = lband[0];
		    lhigh = lband[1] - 1;
		    lbw = 1;
		} else if (i__ == l) {
		    if (l == 2) {
			lb = 1;
		    } else {
			lb = lband[l - 3];
		    }
		    low = lband[l - 2] - lb + 1;
		    medium = lband[l - 1] - lb + 1;
		    lhigh = n - lb + 1;
		    lbw = lband[i__ - 2];
		} else {
		    if (i__ == 2) {
			lb = 1;
		    } else {
			lb = lband[i__ - 3];
		    }
		    low = lband[i__ - 2] - lb + 1;
		    medium = lband[i__ - 1] - lb + 1;
		    lhigh = lband[i__] - 1 - lb + 1;
		    lbw = lband[i__ - 2];
		}
		lenw = medium - low;
/* 	print  *,lbw,lenw,LB,LOW,MEDIUM,LHIGH */
/* 	print  *,LOW+LB-1,MEDIUM+LB-1 */
/* 	PRINT  *,'  SECTION NUMBER',I */
 	       
	printf("\n lbw = %d" ,lbw);
	voronoi(&theta[lb],&phi[lb],&weight[lbw],lenw,low,medium,lhigh,last);
	
	       
	       
	       
	       
/*  Test the areas in all sections */
		if (nband > 1) {
		    if (i__ == 1) {
			area = (1.f - cos(ts[0] * .017453292519943295f)) * 
				6.2831853071795862;
		    } else {
			area = (cos(ts[i__ - 2] * dgr_to_rad) - cos(
				ts[i__ - 1] * dgr_to_rad)) * 
				2*quadpi;
		    }
/* 	  AA=SUM(WEIGHT(1:LBW+LENW-1))/AREA */
		    aa = 0.;
		    i__3 = lbw + lenw - 1;
		    for (nbc = lbw; nbc <= i__3; ++nbc) {
			aa += weight[nbc];
		    }
		    acum += aa;
/* Computing MAX */
		    i__3 = 1, i__4 = it - numth + 1;
		    aat[i__ - max(i__3,i__4)] = aa / area;
		}
/*  end of parallel loop */
	    }
/* Computing MAX */
	    i__3 = 1, i__4 = it - numth + 1;
	    i__2 = max(i__3,i__4);
	    for (i__ = it; i__ >= i__2; --i__) {
/*         WRITE (*,121) I,AAT(I-MAX(1,IT-NUMTH+1)+1) */
/* Computing MAX */
		i__3 = 1, i__4 = it - numth + 1;
		if ((d__1 = aat[i__ - max(i__3,i__4)] - 1.f, abs(d__1)) > .02)
			 {
 		    r__1 = (real) nband * .75f;
		    i__5 = i_nint(&r__1), i__6 = nband - 1;
		    i__3 = 0, i__4 = min(i__5,i__6);
		    nband = max(i__3,i__4);
		    goto L128;
		}
	    }

/*  End of loop over Bands */
	}
	acum = acum / quadpi / 2.f;
 
 	exit(1);
/*  Do not check the sum - if all individual have small errors, the sum will have a small error */
/* 	  IF(DABS(ACUM-1.0)>1.0D-3) THEN */
/* 	    NBAND=MAX0(1,MIN0(NINT(REAL(NBAND)*0.9),NBAND-1)) */
/* 	    GOTO  128 */
/* 	  ELSE */
/* 	   EXIT */
/* 	  ENDIF */
L128:
	;
    }
/*      WRITE (LOUT,120) SUM(WEIGHT(1:ND2))/(16.D0*ATAN(1.D0)), */
/*     &		MAXVAL(WEIGHT(1:ND2)) */
/*    DEALLOCATE(AAT,TS,LBAND) */

free(aat);
free(lband);
free(ts);    
	++weight;
	++theta;
	++phi;
    return 0;
} /* voronoidiag_ */

/****************************************************************************************************************************/
/* Subroutine */ int angstep_(double  *thetast, integer *n)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3;

     /* Local variables */
    static integer i__;
    static double t1, t2, tmp;

/*  Here the equation follows from the integral Ssin(theta)dtheta=1/n */
/*  i.e., equal area bands */
    /* Parameter adjustments */
    --thetast;

    /* Function Body */
    if (*n > 1) {
	t1 = 0.f;
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    tmp = cos(t1) - 1.f / (double) (*n);
/* Computing MIN */
	    r__2 = 1.f, r__3 = dabs(tmp);
	    r__1 = dmin(r__2,r__3);
	    t2 = acos(sign(r__1, tmp));
/* 	 PRINT  *,T1*RAD_TO_DGR,T2*RAD_TO_DGR,COS(T1)-COS(T2) */
	    thetast[i__] = t2 * rad_to_dgr;
	    t1 = t2;
	}
    }
    thetast[*n] = 90.f;
    return 0;
} /* angstep_ */

/****************************************************************************************************************************/

#undef	thetast


#undef	theta
#undef	phi
#undef	weight

#undef	aat
#undef	ls
#undef	ts

**********************   stripack.c
/* stridpack.f -- translated by f2c (version 19991025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    doublereal y;
} stcom_;

#define stcom_1 stcom_

/* Table of constant values */

static integer c__1 = 1;

/* ____________________STRID______________________________________ */
/* Subroutine */ int trmsh3_(n0, tol, x, y, z__, n, list, lptr, lend, lnew, 
	indx, lcnt, near__, next, dist, ier)
integer *n0;
doublereal *tol, *x, *y, *z__;
integer *n, *list, *lptr, *lend, *lnew, *indx, *lcnt, *near__, *next;
doublereal *dist;
integer *ier;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    extern logical left_();
    static doublereal d__;
    static integer i__, j;
    static doublereal d1, d2;
    static integer nexti;
    static doublereal d3;
    static integer i0;
    extern /* Subroutine */ int addnod_();
    static integer lp, kt, ku, lpl, nku;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   01/20/03 */

/*   This is an alternative to TRMESH with the inclusion of */
/* an efficient means of removing duplicate or nearly dupli- */
/* cate nodes. */

/*   This subroutine creates a Delaunay triangulation of a */
/* set of N arbitrarily distributed points, referred to as */
/* nodes, on the surface of the unit sphere.  Refer to Sub- */
/* routine TRMESH for definitions and a list of additional */
/* subroutines.  This routine is an alternative to TRMESH */
/* with the inclusion of an efficient means of removing dup- */
/* licate or nearly duplicate nodes. */

/*   The algorithm has expected time complexity O(N*log(N)) */
/* for random nodal distributions. */


/* On input: */

/*       N0 = Number of nodes, possibly including duplicates. */
/*            N0 .GE. 3. */

/*       TOL = Tolerance defining a pair of duplicate nodes: */
/*             bound on the deviation from 1 of the cosine of */
/*             the angle between the nodes.  Note that */
/*             |1-cos(A)| is approximately A*A/2. */

/* The above parameters are not altered by this routine. */

/*       X,Y,Z = Arrays of length at least N0 containing the */
/*               Cartesian coordinates of nodes.  (X(K),Y(K), */
/*               Z(K)) is referred to as node K, and K is re- */
/*               ferred to as a nodal index.  It is required */
/*               that X(K)**2 + Y(K)**2 + Z(K)**2 = 1 for all */
/*               K.  The first three nodes must not be col- */
/*               linear (lie on a common great circle). */

/*       LIST,LPTR = Arrays of length at least 6*N0-12. */

/*       LEND = Array of length at least N0. */

/*       INDX = Array of length at least N0. */

/*       LCNT = Array of length at least N0 (length N is */
/*              sufficient). */

/*       NEAR,NEXT,DIST = Work space arrays of length at */
/*                        least N0.  The space is used to */
/*                        efficiently determine the nearest */
/*                        triangulation node to each un- */
/*                        processed node for use by ADDNOD. */

/* On output: */

/*       N = Number of nodes in the triangulation.  3 .LE. N */
/*           .LE. N0, or N = 0 if IER < 0. */

/*       X,Y,Z = Arrays containing the Cartesian coordinates */
/*               of the triangulation nodes in the first N */
/*               locations.  The original array elements are */
/*               shifted down as necessary to eliminate dup- */
/*               licate nodes. */

/*       LIST = Set of nodal indexes which, along with LPTR, */
/*              LEND, and LNEW, define the triangulation as a */
/*              set of N adjacency lists -- counterclockwise- */
/*              ordered sequences of neighboring nodes such */
/*              that the first and last neighbors of a bound- */
/*              ary node are boundary nodes (the first neigh- */
/*              bor of an interior node is arbitrary).  In */
/*              order to distinguish between interior and */
/*              boundary nodes, the last neighbor of each */
/*              boundary node is represented by the negative */
/*              of its index. */

/*       LPTR = Set of pointers (LIST indexes) in one-to-one */
/*              correspondence with the elements of LIST. */
/*              LIST(LPTR(I)) indexes the node which follows */
/*              LIST(I) in cyclical counterclockwise order */
/*              (the first neighbor follows the last neigh- */
/*              bor). */

/*       LEND = Set of pointers to adjacency lists.  LEND(K) */
/*              points to the last neighbor of node K for */
/*              K = 1,...,N.  Thus, LIST(LEND(K)) < 0 if and */
/*              only if K is a boundary node. */

/*       LNEW = Pointer to the first empty location in LIST */
/*              and LPTR (list length plus one).  LIST, LPTR, */
/*              LEND, and LNEW are not altered if IER < 0, */
/*              and are incomplete if IER > 0. */

/*       INDX = Array of output (triangulation) nodal indexes */
/*              associated with input nodes.  For I = 1 to */
/*              N0, INDX(I) is the index (for X, Y, and Z) of */
/*              the triangulation node with the same (or */
/*              nearly the same) coordinates as input node I. */

/*       LCNT = Array of integer weights (counts) associated */
/*              with the triangulation nodes.  For I = 1 to */
/*              N, LCNT(I) is the number of occurrences of */
/*              node I in the input node set, and thus the */
/*              number of duplicates is LCNT(I)-1. */

/*       NEAR,NEXT,DIST = Garbage. */

/*       IER = Error indicator: */
/*             IER =  0 if no errors were encountered. */
/*             IER = -1 if N0 < 3 on input. */
/*             IER = -2 if the first three nodes are */
/*                      collinear. */
/*             IER = -3 if Subroutine ADDNOD returns an error */
/*                      flag.  This should not occur. */

/* Modules required by TRMSH3:  ADDNOD, BDYADD, COVSPH, */
/*                                INSERT, INTADD, JRAND, */
/*                                LEFT, LSTPTR, STORE, SWAP, */
/*                                SWPTST, TRFIND */

/* Intrinsic function called by TRMSH3:  ABS */

/* *********************************************************** */


/* Local parameters: */

/* D =        (Negative cosine of) distance from node KT to */
/*              node I */
/* D1,D2,D3 = Distances from node KU to nodes 1, 2, and 3, */
/*              respectively */
/* I,J =      Nodal indexes */
/* I0 =       Index of the node preceding I in a sequence of */
/*              unprocessed nodes:  I = NEXT(I0) */
/* KT =       Index of a triangulation node */
/* KU =       Index of an unprocessed node and DO-loop index */
/* LP =       LIST index (pointer) of a neighbor of KT */
/* LPL =      Pointer to the last neighbor of KT */
/* NEXTI =    NEXT(I) */
/* NKU =      NEAR(KU) */

    /* Parameter adjustments */
    --dist;
    --next;
    --near__;
    --indx;
    --lend;
    --z__;
    --y;
    --x;
    --list;
    --lptr;
    --lcnt;

    /* Function Body */
    if (*n0 < 3) {
	*n = 0;
	*ier = -1;
	return 0;
    }

/* Store the first triangle in the linked list. */

    if (! left_(&x[1], &y[1], &z__[1], &x[2], &y[2], &z__[2], &x[3], &y[3], &
	    z__[3])) {

/*   The first triangle is (3,2,1) = (2,1,3) = (1,3,2). */

	list[1] = 3;
	lptr[1] = 2;
	list[2] = -2;
	lptr[2] = 1;
	lend[1] = 2;

	list[3] = 1;
	lptr[3] = 4;
	list[4] = -3;
	lptr[4] = 3;
	lend[2] = 4;

	list[5] = 2;
	lptr[5] = 6;
	list[6] = -1;
	lptr[6] = 5;
	lend[3] = 6;

    } else if (! left_(&x[2], &y[2], &z__[2], &x[1], &y[1], &z__[1], &x[3], &
	    y[3], &z__[3])) {

/*   The first triangle is (1,2,3):  3 Strictly Left 1->2, */
/*     i.e., node 3 lies in the left hemisphere defined by */
/*     arc 1->2. */

	list[1] = 2;
	lptr[1] = 2;
	list[2] = -3;
	lptr[2] = 1;
	lend[1] = 2;

	list[3] = 3;
	lptr[3] = 4;
	list[4] = -1;
	lptr[4] = 3;
	lend[2] = 4;

	list[5] = 1;
	lptr[5] = 6;
	list[6] = -2;
	lptr[6] = 5;
	lend[3] = 6;

    } else {

/*   The first three nodes are collinear. */

	*n = 0;
	*ier = -2;
	return 0;
    }

/* Initialize LNEW, INDX, and LCNT, and test for N = 3. */

    *lnew = 7;
    indx[1] = 1;
    indx[2] = 2;
    indx[3] = 3;
    lcnt[1] = 1;
    lcnt[2] = 1;
    lcnt[3] = 1;
    if (*n0 == 3) {
	*n = 3;
	*ier = 0;
	return 0;
    }

/* A nearest-node data structure (NEAR, NEXT, and DIST) is */
/*   used to obtain an expected-time (N*log(N)) incremental */
/*   algorithm by enabling constant search time for locating */
/*   each new node in the triangulation. */

/* For each unprocessed node KU, NEAR(KU) is the index of the */
/*   triangulation node closest to KU (used as the starting */
/*   point for the search in Subroutine TRFIND) and DIST(KU) */
/*   is an increasing function of the arc length (angular */
/*   distance) between nodes KU and NEAR(KU):  -Cos(a) for */
/*   arc length a. */

/* Since it is necessary to efficiently find the subset of */
/*   unprocessed nodes associated with each triangulation */
/*   node J (those that have J as their NEAR entries), the */
/*   subsets are stored in NEAR and NEXT as follows:  for */
/*   each node J in the triangulation, I = NEAR(J) is the */
/*   first unprocessed node in J's set (with I = 0 if the */
/*   set is empty), L = NEXT(I) (if I > 0) is the second, */
/*   NEXT(L) (if L > 0) is the third, etc.  The nodes in each */
/*   set are initially ordered by increasing indexes (which */
/*   maximizes efficiency) but that ordering is not main- */
/*   tained as the data structure is updated. */

/* Initialize the data structure for the single triangle. */

    near__[1] = 0;
    near__[2] = 0;
    near__[3] = 0;
    for (ku = *n0; ku >= 4; --ku) {
	d1 = -(x[ku] * x[1] + y[ku] * y[1] + z__[ku] * z__[1]);
	d2 = -(x[ku] * x[2] + y[ku] * y[2] + z__[ku] * z__[2]);
	d3 = -(x[ku] * x[3] + y[ku] * y[3] + z__[ku] * z__[3]);
	if (d1 <= d2 && d1 <= d3) {
	    near__[ku] = 1;
	    dist[ku] = d1;
	    next[ku] = near__[1];
	    near__[1] = ku;
	} else if (d2 <= d1 && d2 <= d3) {
	    near__[ku] = 2;
	    dist[ku] = d2;
	    next[ku] = near__[2];
	    near__[2] = ku;
	} else {
	    near__[ku] = 3;
	    dist[ku] = d3;
	    next[ku] = near__[3];
	    near__[3] = ku;
	}
/* L1: */
    }

/* Loop on unprocessed nodes KU.  KT is the number of nodes */
/*   in the triangulation, and NKU = NEAR(KU). */

    kt = 3;
    i__1 = *n0;
    for (ku = 4; ku <= i__1; ++ku) {
	nku = near__[ku];

/* Remove KU from the set of unprocessed nodes associated */
/*   with NEAR(KU). */

	i__ = nku;
	if (near__[i__] == ku) {
	    near__[i__] = next[ku];
	} else {
	    i__ = near__[i__];
L2:
	    i0 = i__;
	    i__ = next[i0];
	    if (i__ != ku) {
		goto L2;
	    }
	    next[i0] = next[ku];
	}
	near__[ku] = 0;

/* Bypass duplicate nodes. */

	if (dist[ku] <= *tol - 1.) {
	    indx[ku] = -nku;
	    ++lcnt[nku];
	    goto L6;
	}

/* Add a new triangulation node KT with LCNT(KT) = 1. */

	++kt;
	x[kt] = x[ku];
	y[kt] = y[ku];
	z__[kt] = z__[ku];
	indx[ku] = kt;
	lcnt[kt] = 1;
	addnod_(&nku, &kt, &x[1], &y[1], &z__[1], &list[1], &lptr[1], &lend[1]
		, lnew, ier);
	if (*ier != 0) {
	    *n = 0;
	    *ier = -3;
	    return 0;
	}

/* Loop on neighbors J of node KT. */

	lpl = lend[kt];
	lp = lpl;
L3:
	lp = lptr[lp];
	j = (i__2 = list[lp], abs(i__2));

/* Loop on elements I in the sequence of unprocessed nodes */
/*   associated with J:  KT is a candidate for replacing J */
/*   as the nearest triangulation node to I.  The next value */
/*   of I in the sequence, NEXT(I), must be saved before I */
/*   is moved because it is altered by adding I to KT's set. */

	i__ = near__[j];
L4:
	if (i__ == 0) {
	    goto L5;
	}
	nexti = next[i__];

/* Test for the distance from I to KT less than the distance */
/*   from I to J. */

	d__ = -(x[i__] * x[kt] + y[i__] * y[kt] + z__[i__] * z__[kt]);
	if (d__ < dist[i__]) {

/* Replace J by KT as the nearest triangulation node to I: */
/*   update NEAR(I) and DIST(I), and remove I from J's set */
/*   of unprocessed nodes and add it to KT's set. */

	    near__[i__] = kt;
	    dist[i__] = d__;
	    if (i__ == near__[j]) {
		near__[j] = nexti;
	    } else {
		next[i0] = nexti;
	    }
	    next[i__] = near__[kt];
	    near__[kt] = i__;
	} else {
	    i0 = i__;
	}

/* Bottom of loop on I. */

	i__ = nexti;
	goto L4;

/* Bottom of loop on neighbors J. */

L5:
	if (lp != lpl) {
	    goto L3;
	}
L6:
	;
    }
    *n = kt;
    *ier = 0;
    return 0;
} /* trmsh3_ */

/* stripack.dbl sent by Robert on 06/03/03 */
/* Subroutine */ int addnod_(nst, k, x, y, z__, list, lptr, lend, lnew, ier)
integer *nst, *k;
doublereal *x, *y, *z__;
integer *list, *lptr, *lend, *lnew, *ier;
{
    /* Initialized data */

    static doublereal tol = 0.;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int swap_();
    static integer lpo1s, l;
    static doublereal p[3], b1, b2, b3;
    static integer i1, i2, i3, kk;
    extern /* Subroutine */ int bdyadd_();
    static integer lp;
    extern /* Subroutine */ int intadd_(), trfind_(), covsph_();
    static integer in1, io1, io2, km1;
    extern integer lstptr_();
    extern logical swptst_();
    static integer lpf, ist, lpo1;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   01/08/03 */

/*   This subroutine adds node K to a triangulation of the */
/* convex hull of nodes 1,...,K-1, producing a triangulation */
/* of the convex hull of nodes 1,...,K. */

/*   The algorithm consists of the following steps:  node K */
/* is located relative to the triangulation (TRFIND), its */
/* index is added to the data structure (INTADD or BDYADD), */
/* and a sequence of swaps (SWPTST and SWAP) are applied to */
/* the arcs opposite K so that all arcs incident on node K */
/* and opposite node K are locally optimal (satisfy the cir- */
/* cumcircle test).  Thus, if a Delaunay triangulation is */
/* input, a Delaunay triangulation will result. */


/* On input: */

/*       NST = Index of a node at which TRFIND begins its */
/*             search.  Search time depends on the proximity */
/*             of this node to K.  If NST < 1, the search is */
/*             begun at node K-1. */

/*       K = Nodal index (index for X, Y, Z, and LEND) of the */
/*           new node to be added.  K .GE. 4. */

/*       X,Y,Z = Arrays of length .GE. K containing Car- */
/*               tesian coordinates of the nodes. */
/*               (X(I),Y(I),Z(I)) defines node I for */
/*               I = 1,...,K. */

/* The above parameters are not altered by this routine. */

/*       LIST,LPTR,LEND,LNEW = Data structure associated with */
/*                             the triangulation of nodes 1 */
/*                             to K-1.  The array lengths are */
/*                             assumed to be large enough to */
/*                             add node K.  Refer to Subrou- */
/*                             tine TRMESH. */

/* On output: */

/*       LIST,LPTR,LEND,LNEW = Data structure updated with */
/*                             the addition of node K as the */
/*                             last entry unless IER .NE. 0 */
/*                             and IER .NE. -3, in which case */
/*                             the arrays are not altered. */

/*       IER = Error indicator: */
/*             IER =  0 if no errors were encountered. */
/*             IER = -1 if K is outside its valid range */
/*                      on input. */
/*             IER = -2 if all nodes (including K) are col- */
/*                      linear (lie on a common geodesic). */
/*             IER =  L if nodes L and K coincide for some */
/*                      L < K.  Refer to TOL below. */

/* Modules required by ADDNOD:  BDYADD, COVSPH, INSERT, */
/*                                INTADD, JRAND, LSTPTR, */
/*                                STORE, SWAP, SWPTST, */
/*                                TRFIND */

/* Intrinsic function called by ADDNOD:  ABS */

/* *********************************************************** */


/* Local parameters: */

/* B1,B2,B3 = Unnormalized barycentric coordinates returned */
/*              by TRFIND. */
/* I1,I2,I3 = Vertex indexes of a triangle containing K */
/* IN1 =      Vertex opposite K:  first neighbor of IO2 */
/*              that precedes IO1.  IN1,IO1,IO2 are in */
/*              counterclockwise order. */
/* IO1,IO2 =  Adjacent neighbors of K defining an arc to */
/*              be tested for a swap */
/* IST =      Index of node at which TRFIND begins its search */
/* KK =       Local copy of K */
/* KM1 =      K-1 */
/* L =        Vertex index (I1, I2, or I3) returned in IER */
/*              if node K coincides with a vertex */
/* LP =       LIST pointer */
/* LPF =      LIST pointer to the first neighbor of K */
/* LPO1 =     LIST pointer to IO1 */
/* LPO1S =    Saved value of LPO1 */
/* P =        Cartesian coordinates of node K */
/* TOL =      Tolerance defining coincident nodes:  bound on */
/*              the deviation from 1 of the cosine of the */
/*              angle between the nodes. */
/*              Note that |1-cos(A)| is approximately A*A/2. */

    /* Parameter adjustments */
    --lend;
    --z__;
    --y;
    --x;
    --list;
    --lptr;

    /* Function Body */

    kk = *k;
    if (kk < 4) {
	goto L3;
    }

/* Initialization: */

    km1 = kk - 1;
    ist = *nst;
    if (ist < 1) {
	ist = km1;
    }
    p[0] = x[kk];
    p[1] = y[kk];
    p[2] = z__[kk];

/* Find a triangle (I1,I2,I3) containing K or the rightmost */
/*   (I1) and leftmost (I2) visible boundary nodes as viewed */
/*   from node K. */

    trfind_(&ist, p, &km1, &x[1], &y[1], &z__[1], &list[1], &lptr[1], &lend[1]
	    , &b1, &b2, &b3, &i1, &i2, &i3);

/*   Test for collinear or (nearly) duplicate nodes. */

    if (i1 == 0) {
	goto L4;
    }
    l = i1;
    if (p[0] * x[l] + p[1] * y[l] + p[2] * z__[l] >= 1. - tol) {
	goto L5;
    }
    l = i2;
    if (p[0] * x[l] + p[1] * y[l] + p[2] * z__[l] >= 1. - tol) {
	goto L5;
    }
    if (i3 != 0) {
	l = i3;
	if (p[0] * x[l] + p[1] * y[l] + p[2] * z__[l] >= 1. - tol) {
	    goto L5;
	}
	intadd_(&kk, &i1, &i2, &i3, &list[1], &lptr[1], &lend[1], lnew);
    } else {
	if (i1 != i2) {
	    bdyadd_(&kk, &i1, &i2, &list[1], &lptr[1], &lend[1], lnew);
	} else {
	    covsph_(&kk, &i1, &list[1], &lptr[1], &lend[1], lnew);
	}
    }
    *ier = 0;

/* Initialize variables for optimization of the */
/*   triangulation. */

    lp = lend[kk];
    lpf = lptr[lp];
    io2 = list[lpf];
    lpo1 = lptr[lpf];
    io1 = (i__1 = list[lpo1], abs(i__1));

/* Begin loop:  find the node opposite K. */

L1:
    lp = lstptr_(&lend[io1], &io2, &list[1], &lptr[1]);
    if (list[lp] < 0) {
	goto L2;
    }
    lp = lptr[lp];
    in1 = (i__1 = list[lp], abs(i__1));

/* Swap test:  if a swap occurs, two new arcs are */
/*             opposite K and must be tested. */

    lpo1s = lpo1;
    if (! swptst_(&in1, &kk, &io1, &io2, &x[1], &y[1], &z__[1])) {
	goto L2;
    }
    swap_(&in1, &kk, &io1, &io2, &list[1], &lptr[1], &lend[1], &lpo1);
    if (lpo1 == 0) {

/*   A swap is not possible because KK and IN1 are already */
/*     adjacent.  This error in SWPTST only occurs in the */
/*     neutral case and when there are nearly duplicate */
/*     nodes. */

	lpo1 = lpo1s;
	goto L2;
    }
    io1 = in1;
    goto L1;

/* No swap occurred.  Test for termination and reset */
/*   IO2 and IO1. */

L2:
    if (lpo1 == lpf || list[lpo1] < 0) {
	return 0;
    }
    io2 = io1;
    lpo1 = lptr[lpo1];
    io1 = (i__1 = list[lpo1], abs(i__1));
    goto L1;

/* KK < 4. */

L3:
    *ier = -1;
    return 0;

/* All nodes are collinear. */

L4:
    *ier = -2;
    return 0;

/* Nodes L and K coincide. */

L5:
    *ier = l;
    return 0;
} /* addnod_ */

doublereal angle_(v1, v2, v3)
doublereal *v1, *v2, *v3;
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double sqrt(), acos();

    /* Local variables */
    extern logical left_();
    static doublereal a;
    static integer i__;
    static doublereal ca, s21, s23, u21[3], u23[3];


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   06/03/03 */

/*   Given a sequence of three nodes (V1,V2,V3) on the sur- */
/* face of the unit sphere, this function returns the */
/* interior angle at V2 -- the dihedral angle between the */
/* plane defined by V2 and V3 (and the origin) and the plane */
/* defined by V2 and V1 or, equivalently, the angle between */
/* the normals V2 X V3 and V2 X V1.  Note that the angle is */
/* in the range 0 to Pi if V3 Left V1->V2, Pi to 2*Pi other- */
/* wise.  The surface area of a spherical polygon with CCW- */
/* ordered vertices V1, V2, ..., Vm is Asum - (m-2)*Pi, where */
/* Asum is the sum of the m interior angles computed from the */
/* sequences (Vm,V1,V2), (V1,V2,V3), (V2,V3,V4), ..., */
/* (Vm-1,Vm,V1). */


/* On input: */

/*       V1,V2,V3 = Arrays of length 3 containing the Carte- */
/*                  sian coordinates of unit vectors.  These */
/*                  vectors, if nonzero, are implicitly */
/*                  scaled to have length 1. */

/* Input parameters are not altered by this function. */

/* On output: */

/*       ANGLE = Angle defined above, or 0 if V2 X V1 = 0 or */
/*               V2 X V3 = 0. */

/* Module required by ANGLE:  LEFT */

/* Intrinsic functions called by ANGLE:  ACOS, SQRT */

/* *********************************************************** */


/* Local parameters: */

/* A =       Interior angle at V2 */
/* CA =      cos(A) */
/* I =       DO-loop index and index for U21 and U23 */
/* S21,S23 = Sum of squared components of U21 and U23 */
/* U21,U23 = Unit normal vectors to the planes defined by */
/*             pairs of triangle vertices */


/* Compute cross products U21 = V2 X V1 and U23 = V2 X V3. */

    /* Parameter adjustments */
    --v3;
    --v2;
    --v1;

    /* Function Body */
    u21[0] = v2[2] * v1[3] - v2[3] * v1[2];
    u21[1] = v2[3] * v1[1] - v2[1] * v1[3];
    u21[2] = v2[1] * v1[2] - v2[2] * v1[1];

    u23[0] = v2[2] * v3[3] - v2[3] * v3[2];
    u23[1] = v2[3] * v3[1] - v2[1] * v3[3];
    u23[2] = v2[1] * v3[2] - v2[2] * v3[1];

/* Normalize U21 and U23 to unit vectors. */

    s21 = 0.;
    s23 = 0.;
    for (i__ = 1; i__ <= 3; ++i__) {
	s21 += u21[i__ - 1] * u21[i__ - 1];
	s23 += u23[i__ - 1] * u23[i__ - 1];
/* L1: */
    }

/* Test for a degenerate triangle associated with collinear */
/*   vertices. */

    if (s21 == 0. || s23 == 0.) {
	ret_val = 0.;
	return ret_val;
    }
    s21 = sqrt(s21);
    s23 = sqrt(s23);
    for (i__ = 1; i__ <= 3; ++i__) {
	u21[i__ - 1] /= s21;
	u23[i__ - 1] /= s23;
/* L2: */
    }

/* Compute the angle A between normals: */

/*   CA = cos(A) = <U21,U23> */

    ca = u21[0] * u23[0] + u21[1] * u23[1] + u21[2] * u23[2];
    if (ca < -1.) {
	ca = -1.;
    }
    if (ca > 1.) {
	ca = 1.;
    }
    a = acos(ca);

/* Adjust A to the interior angle:  A > Pi iff */
/*   V3 Right V1->V2. */

    if (! left_(&v1[1], &v1[2], &v1[3], &v2[1], &v2[2], &v2[3], &v3[1], &v3[2]
	    , &v3[3])) {
	a = acos(-1.) * 2. - a;
    }
    ret_val = a;
    return ret_val;
} /* angle_ */

doublereal areas_(v1, v2, v3)
doublereal *v1, *v2, *v3;
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double sqrt(), acos();

    /* Local variables */
    static integer i__;
    static doublereal a1, a2, a3, s12, s31, s23, u12[3], u23[3], u31[3], ca1, 
	    ca2, ca3;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   06/22/98 */

/*   This function returns the area of a spherical triangle */
/* on the unit sphere. */


/* On input: */

/*       V1,V2,V3 = Arrays of length 3 containing the Carte- */
/*                  sian coordinates of unit vectors (the */
/*                  three triangle vertices in any order). */
/*                  These vectors, if nonzero, are implicitly */
/*                  scaled to have length 1. */

/* Input parameters are not altered by this function. */

/* On output: */

/*       AREAS = Area of the spherical triangle defined by */
/*               V1, V2, and V3 in the range 0 to 2*PI (the */
/*               area of a hemisphere).  AREAS = 0 (or 2*PI) */
/*               if and only if V1, V2, and V3 lie in (or */
/*               close to) a plane containing the origin. */

/* Modules required by AREAS:  None */

/* Intrinsic functions called by AREAS:  ACOS, SQRT */

/* *********************************************************** */


/* Local parameters: */

/* A1,A2,A3 =    Interior angles of the spherical triangle */
/* CA1,CA2,CA3 = cos(A1), cos(A2), and cos(A3), respectively */
/* I =           DO-loop index and index for Uij */
/* S12,S23,S31 = Sum of squared components of U12, U23, U31 */
/* U12,U23,U31 = Unit normal vectors to the planes defined by */
/*                 pairs of triangle vertices */


/* Compute cross products Uij = Vi X Vj. */

    /* Parameter adjustments */
    --v3;
    --v2;
    --v1;

    /* Function Body */
    u12[0] = v1[2] * v2[3] - v1[3] * v2[2];
    u12[1] = v1[3] * v2[1] - v1[1] * v2[3];
    u12[2] = v1[1] * v2[2] - v1[2] * v2[1];

    u23[0] = v2[2] * v3[3] - v2[3] * v3[2];
    u23[1] = v2[3] * v3[1] - v2[1] * v3[3];
    u23[2] = v2[1] * v3[2] - v2[2] * v3[1];

    u31[0] = v3[2] * v1[3] - v3[3] * v1[2];
    u31[1] = v3[3] * v1[1] - v3[1] * v1[3];
    u31[2] = v3[1] * v1[2] - v3[2] * v1[1];

/* Normalize Uij to unit vectors. */

    s12 = 0.;
    s23 = 0.;
    s31 = 0.;
    for (i__ = 1; i__ <= 3; ++i__) {
	s12 += u12[i__ - 1] * u12[i__ - 1];
	s23 += u23[i__ - 1] * u23[i__ - 1];
	s31 += u31[i__ - 1] * u31[i__ - 1];
/* L2: */
    }

/* Test for a degenerate triangle associated with collinear */
/*   vertices. */

    if (s12 == 0. || s23 == 0. || s31 == 0.) {
	ret_val = 0.;
	return ret_val;
    }
    s12 = sqrt(s12);
    s23 = sqrt(s23);
    s31 = sqrt(s31);
    for (i__ = 1; i__ <= 3; ++i__) {
	u12[i__ - 1] /= s12;
	u23[i__ - 1] /= s23;
	u31[i__ - 1] /= s31;
/* L3: */
    }

/* Compute interior angles Ai as the dihedral angles between */
/*   planes: */
/*           CA1 = cos(A1) = -<U12,U31> */
/*           CA2 = cos(A2) = -<U23,U12> */
/*           CA3 = cos(A3) = -<U31,U23> */

    ca1 = -u12[0] * u31[0] - u12[1] * u31[1] - u12[2] * u31[2];
    ca2 = -u23[0] * u12[0] - u23[1] * u12[1] - u23[2] * u12[2];
    ca3 = -u31[0] * u23[0] - u31[1] * u23[1] - u31[2] * u23[2];
    if (ca1 < -1.) {
	ca1 = -1.;
    }
    if (ca1 > 1.) {
	ca1 = 1.;
    }
    if (ca2 < -1.) {
	ca2 = -1.;
    }
    if (ca2 > 1.) {
	ca2 = 1.;
    }
    if (ca3 < -1.) {
	ca3 = -1.;
    }
    if (ca3 > 1.) {
	ca3 = 1.;
    }
    a1 = acos(ca1);
    a2 = acos(ca2);
    a3 = acos(ca3);

/* Compute AREAS = A1 + A2 + A3 - PI. */

    ret_val = a1 + a2 + a3 - acos(-1.);
    if (ret_val < 0.) {
	ret_val = 0.;
    }
    return ret_val;
} /* areas_ */

doublereal areav_(k, n, x, y, z__, list, lptr, lend, ier)
integer *k, *n;
doublereal *x, *y, *z__;
integer *list, *lptr, *lend, *ier;
{
    /* Initialized data */

    static doublereal amax = 6.28;

    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer ierr;
    static doublereal asum, a;
    extern doublereal areas_();
    static doublereal c0[3], c2[3], c3[3];
    static logical first;
    static integer n1, n2, n3;
    static doublereal v1[3], v2[3], v3[3];
    static integer lp;
    extern /* Subroutine */ int circum_();
    static integer lpl;


/* *********************************************************** */

/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   10/25/02 */

/*   Given a Delaunay triangulation and the index K of an */
/* interior node, this subroutine returns the (surface) area */
/* of the Voronoi region associated with node K.  The Voronoi */
/* region is the polygon whose vertices are the circumcenters */
/* of the triangles that contain node K, where a triangle */
/* circumcenter is the point (unit vector) lying at the same */
/* angular distance from the three vertices and contained in */
/* the same hemisphere as the vertices. */


/* On input: */

/*       K = Nodal index in the range 1 to N. */

/*       N = Number of nodes in the triangulation.  N > 3. */

/*       X,Y,Z = Arrays of length N containing the Cartesian */
/*               coordinates of the nodes (unit vectors). */

/*       LIST,LPTR,LEND = Data structure defining the trian- */
/*                        gulation.  Refer to Subroutine */
/*                        TRMESH. */

/* Input parameters are not altered by this function. */

/* On output: */

/*       AREAV = Area of Voronoi region K unless IER > 0, */
/*               in which case AREAV = 0. */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered. */
/*             IER = 1 if K or N is outside its valid range */
/*                     on input. */
/*             IER = 2 if K indexes a boundary node. */
/*             IER = 3 if an error flag is returned by CIRCUM */
/*                     (null triangle). */
/*             IER = 4 if AREAS returns a value greater than */
/*                     AMAX (defined below). */

/* Modules required by AREAV:  AREAS, CIRCUM */

/* *********************************************************** */


/* Maximum valid triangle area is less than 2*Pi: */

    /* Parameter adjustments */
    --lend;
    --z__;
    --y;
    --x;
    --list;
    --lptr;

    /* Function Body */

/* Test for invalid input. */

    if (*k < 1 || *k > *n || *n <= 3) {
	goto L11;
    }

/* Initialization:  Set N3 to the last neighbor of N1 = K. */
/*   FIRST = TRUE only for the first triangle. */
/*   The Voronoi region area is accumulated in ASUM. */

    n1 = *k;
    v1[0] = x[n1];
    v1[1] = y[n1];
    v1[2] = z__[n1];
    lpl = lend[n1];
    n3 = list[lpl];
    if (n3 < 0) {
	goto L12;
    }
    lp = lpl;
    first = TRUE_;
    asum = 0.;

/* Loop on triangles (N1,N2,N3) containing N1 = K. */

L1:
    n2 = n3;
    lp = lptr[lp];
    n3 = list[lp];
    v2[0] = x[n2];
    v2[1] = y[n2];
    v2[2] = z__[n2];
    v3[0] = x[n3];
    v3[1] = y[n3];
    v3[2] = z__[n3];
    if (first) {

/* First triangle:  compute the circumcenter C3 and save a */
/*   copy in C0. */

	circum_(v1, v2, v3, c3, &ierr);
	if (ierr != 0) {
	    goto L13;
	}
	c0[0] = c3[0];
	c0[1] = c3[1];
	c0[2] = c3[2];
	first = FALSE_;
    } else {

/* Set C2 to C3, compute the new circumcenter C3, and compute */
/*   the area A of triangle (V1,C2,C3). */

	c2[0] = c3[0];
	c2[1] = c3[1];
	c2[2] = c3[2];
	circum_(v1, v2, v3, c3, &ierr);
	if (ierr != 0) {
	    goto L13;
	}
	a = areas_(v1, c2, c3);
	if (a > amax) {
	    goto L14;
	}
	asum += a;
    }

/* Bottom on loop on neighbors of K. */

    if (lp != lpl) {
	goto L1;
    }

/* Compute the area of triangle (V1,C3,C0). */

    a = areas_(v1, c3, c0);
    if (a > amax) {
	goto L14;
    }
    asum += a;

/* No error encountered. */

    *ier = 0;
    ret_val = asum;
    return ret_val;

/* Invalid input. */

L11:
    *ier = 1;
    ret_val = 0.;
    return ret_val;

/* K indexes a boundary node. */

L12:
    *ier = 2;
    ret_val = 0.;
    return ret_val;

/* Error in CIRCUM. */

L13:
    *ier = 3;
    ret_val = 0.;
    return ret_val;

/* AREAS value larger than AMAX. */

L14:
    *ier = 4;
    ret_val = 0.;
    return ret_val;
} /* areav_ */

doublereal areav_new__(k, n, x, y, z__, list, lptr, lend, ier)
integer *k, *n;
doublereal *x, *y, *z__;
integer *list, *lptr, *lend, *ier;
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double acos();

    /* Local variables */
    static integer ierr;
    static doublereal asum;
    static integer m;
    extern doublereal angle_();
    static real areav;
    static doublereal c1[3], c2[3], c3[3];
    static integer n1, n2, n3;
    static doublereal v1[3], v2[3], v3[3];
    static integer lp;
    extern /* Subroutine */ int circum_();
    static doublereal c1s[3], c2s[3];
    static integer lpl;


/* *********************************************************** */

/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   06/03/03 */

/*   Given a Delaunay triangulation and the index K of an */
/* interior node, this subroutine returns the (surface) area */
/* of the Voronoi region associated with node K.  The Voronoi */
/* region is the polygon whose vertices are the circumcenters */
/* of the triangles that contain node K, where a triangle */
/* circumcenter is the point (unit vector) lying at the same */
/* angular distance from the three vertices and contained in */
/* the same hemisphere as the vertices.  The Voronoi region */
/* area is computed as Asum-(m-2)*Pi, where m is the number */
/* of Voronoi vertices (neighbors of K) and Asum is the sum */
/* of interior angles at the vertices. */


/* On input: */

/*       K = Nodal index in the range 1 to N. */

/*       N = Number of nodes in the triangulation.  N > 3. */

/*       X,Y,Z = Arrays of length N containing the Cartesian */
/*               coordinates of the nodes (unit vectors). */

/*       LIST,LPTR,LEND = Data structure defining the trian- */
/*                        gulation.  Refer to Subroutine */
/*                        TRMESH. */

/* Input parameters are not altered by this function. */

/* On output: */

/*       AREAV = Area of Voronoi region K unless IER > 0, */
/*               in which case AREAV = 0. */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered. */
/*             IER = 1 if K or N is outside its valid range */
/*                     on input. */
/*             IER = 2 if K indexes a boundary node. */
/*             IER = 3 if an error flag is returned by CIRCUM */
/*                     (null triangle). */

/* Modules required by AREAV:  ANGLE, CIRCUM */

/* Intrinsic functions called by AREAV:  ACOS, DBLE */

/* *********************************************************** */


/* Test for invalid input. */

    /* Parameter adjustments */
    --lend;
    --z__;
    --y;
    --x;
    --list;
    --lptr;

    /* Function Body */
    if (*k < 1 || *k > *n || *n <= 3) {
	goto L11;
    }

/* Initialization:  Set N3 to the last neighbor of N1 = K. */
/*   The number of neighbors and the sum of interior angles */
/*   are accumulated in M and ASUM, respectively. */

    n1 = *k;
    v1[0] = x[n1];
    v1[1] = y[n1];
    v1[2] = z__[n1];
    lpl = lend[n1];
    n3 = list[lpl];
    if (n3 < 0) {
	goto L12;
    }
    lp = lpl;
    m = 0;
    asum = 0.;

/* Loop on triangles (N1,N2,N3) containing N1 = K. */

L1:
    ++m;
    n2 = n3;
    lp = lptr[lp];
    n3 = list[lp];
    v2[0] = x[n2];
    v2[1] = y[n2];
    v2[2] = z__[n2];
    v3[0] = x[n3];
    v3[1] = y[n3];
    v3[2] = z__[n3];
    if (m == 1) {

/* First triangle:  compute the circumcenter C2 and save a */
/*   copy in C1S. */

	circum_(v1, v2, v3, c2, &ierr);
	if (ierr != 0) {
	    goto L13;
	}
	c1s[0] = c2[0];
	c1s[1] = c2[1];
	c1s[2] = c2[2];
    } else if (m == 2) {

/* Second triangle:  compute the circumcenter C3 and save a */
/*   copy in C2S. */

	circum_(v1, v2, v3, c3, &ierr);
	if (ierr != 0) {
	    goto L13;
	}
	c2s[0] = c3[0];
	c2s[1] = c3[1];
	c2s[2] = c3[2];
    } else {

/* Set C1 to C2, set C2 to C3, compute the new circumcenter */
/*   C3, and compute the interior angle at C2 from the */
/*   sequence of vertices (C1,C2,C3). */

	c1[0] = c2[0];
	c1[1] = c2[1];
	c1[2] = c2[2];
	c2[0] = c3[0];
	c2[1] = c3[1];
	c2[2] = c3[2];
	circum_(v1, v2, v3, c3, &ierr);
	if (ierr != 0) {
	    goto L13;
	}
	asum += angle_(c1, c2, c3);
    }

/* Bottom on loop on neighbors of K. */

    if (lp != lpl) {
	goto L1;
    }

/* C3 is the last vertex.  Compute its interior angle from */
/*   the sequence (C2,C3,C1S). */

    asum += angle_(c2, c3, c1s);

/* Compute the interior angle at C1S from */
/*   the sequence (C3,C1S,C2S). */

    asum += angle_(c3, c1s, c2s);

/* No error encountered. */

    *ier = 0;
    ret_val = asum - (doublereal) (m - 2) * acos(-1.);
    return ret_val;

/* Invalid input. */

L11:
    *ier = 1;
    areav = (float)0.;
    return ret_val;

/* K indexes a boundary node. */

L12:
    *ier = 2;
    areav = (float)0.;
    return ret_val;

/* Error in CIRCUM. */

L13:
    *ier = 3;
    areav = (float)0.;
    return ret_val;
} /* areav_new__ */

/* Subroutine */ int bdyadd_(kk, i1, i2, list, lptr, lend, lnew)
integer *kk, *i1, *i2, *list, *lptr, *lend, *lnew;
{
    static integer lsav, nsav, next, k, n1, n2, lp;
    extern /* Subroutine */ int insert_();


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/11/96 */

/*   This subroutine adds a boundary node to a triangulation */
/* of a set of KK-1 points on the unit sphere.  The data */
/* structure is updated with the insertion of node KK, but no */
/* optimization is performed. */

/*   This routine is identical to the similarly named routine */
/* in TRIPACK. */


/* On input: */

/*       KK = Index of a node to be connected to the sequence */
/*            of all visible boundary nodes.  KK .GE. 1 and */
/*            KK must not be equal to I1 or I2. */

/*       I1 = First (rightmost as viewed from KK) boundary */
/*            node in the triangulation that is visible from */
/*            node KK (the line segment KK-I1 intersects no */
/*            arcs. */

/*       I2 = Last (leftmost) boundary node that is visible */
/*            from node KK.  I1 and I2 may be determined by */
/*            Subroutine TRFIND. */

/* The above parameters are not altered by this routine. */

/*       LIST,LPTR,LEND,LNEW = Triangulation data structure */
/*                             created by Subroutine TRMESH. */
/*                             Nodes I1 and I2 must be in- */
/*                             cluded in the triangulation. */

/* On output: */

/*       LIST,LPTR,LEND,LNEW = Data structure updated with */
/*                             the addition of node KK.  Node */
/*                             KK is connected to I1, I2, and */
/*                             all boundary nodes in between. */

/* Module required by BDYADD:  INSERT */

/* *********************************************************** */


/* Local parameters: */

/* K =     Local copy of KK */
/* LP =    LIST pointer */
/* LSAV =  LIST pointer */
/* N1,N2 = Local copies of I1 and I2, respectively */
/* NEXT =  Boundary node visible from K */
/* NSAV =  Boundary node visible from K */

    /* Parameter adjustments */
    --lend;
    --lptr;
    --list;

    /* Function Body */
    k = *kk;
    n1 = *i1;
    n2 = *i2;

/* Add K as the last neighbor of N1. */

    lp = lend[n1];
    lsav = lptr[lp];
    lptr[lp] = *lnew;
    list[*lnew] = -k;
    lptr[*lnew] = lsav;
    lend[n1] = *lnew;
    ++(*lnew);
    next = -list[lp];
    list[lp] = next;
    nsav = next;

/* Loop on the remaining boundary nodes between N1 and N2, */
/*   adding K as the first neighbor. */

L1:
    lp = lend[next];
    insert_(&k, &lp, &list[1], &lptr[1], lnew);
    if (next == n2) {
	goto L2;
    }
    next = -list[lp];
    list[lp] = next;
    goto L1;

/* Add the boundary nodes between N1 and N2 as neighbors */
/*   of node K. */

L2:
    lsav = *lnew;
    list[*lnew] = n1;
    lptr[*lnew] = *lnew + 1;
    ++(*lnew);
    next = nsav;

L3:
    if (next == n2) {
	goto L4;
    }
    list[*lnew] = next;
    lptr[*lnew] = *lnew + 1;
    ++(*lnew);
    lp = lend[next];
    next = list[lp];
    goto L3;

L4:
    list[*lnew] = -n2;
    lptr[*lnew] = lsav;
    lend[k] = *lnew;
    ++(*lnew);
    return 0;
} /* bdyadd_ */

/* Subroutine */ int bnodes_(n, list, lptr, lend, nodes, nb, na, nt)
integer *n, *list, *lptr, *lend, *nodes, *nb, *na, *nt;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k, n0, lp, nn, nst;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   06/26/96 */

/*   Given a triangulation of N nodes on the unit sphere */
/* created by Subroutine TRMESH, this subroutine returns an */
/* array containing the indexes (if any) of the counterclock- */
/* wise-ordered sequence of boundary nodes -- the nodes on */
/* the boundary of the convex hull of the set of nodes.  (The */
/* boundary is empty if the nodes do not lie in a single */
/* hemisphere.)  The numbers of boundary nodes, arcs, and */
/* triangles are also returned. */


/* On input: */

/*       N = Number of nodes in the triangulation.  N .GE. 3. */

/*       LIST,LPTR,LEND = Data structure defining the trian- */
/*                        gulation.  Refer to Subroutine */
/*                        TRMESH. */

/* The above parameters are not altered by this routine. */

/*       NODES = Integer array of length at least NB */
/*               (NB .LE. N). */

/* On output: */

/*       NODES = Ordered sequence of boundary node indexes */
/*               in the range 1 to N (in the first NB loca- */
/*               tions). */

/*       NB = Number of boundary nodes. */

/*       NA,NT = Number of arcs and triangles, respectively, */
/*               in the triangulation. */

/* Modules required by BNODES:  None */

/* *********************************************************** */


/* Local parameters: */

/* K =   NODES index */
/* LP =  LIST pointer */
/* N0 =  Boundary node to be added to NODES */
/* NN =  Local copy of N */
/* NST = First element of nodes (arbitrarily chosen to be */
/*         the one with smallest index) */

    /* Parameter adjustments */
    --lend;
    --list;
    --lptr;
    --nodes;

    /* Function Body */
    nn = *n;

/* Search for a boundary node. */

    i__1 = nn;
    for (nst = 1; nst <= i__1; ++nst) {
	lp = lend[nst];
	if (list[lp] < 0) {
	    goto L2;
	}
/* L1: */
    }

/* The triangulation contains no boundary nodes. */

    *nb = 0;
    *na = (nn - 2) * 3;
    *nt = nn - 2 << 1;
    return 0;

/* NST is the first boundary node encountered.  Initialize */
/*   for traversal of the boundary. */

L2:
    nodes[1] = nst;
    k = 1;
    n0 = nst;

/* Traverse the boundary in counterclockwise order. */

L3:
    lp = lend[n0];
    lp = lptr[lp];
    n0 = list[lp];
    if (n0 == nst) {
	goto L4;
    }
    ++k;
    nodes[k] = n0;
    goto L3;

/* Store the counts. */

L4:
    *nb = k;
    *nt = (*n << 1) - *nb - 2;
    *na = *nt + *n - 1;
    return 0;
} /* bnodes_ */

/* Subroutine */ int circle_(k, xc, yc, ier)
integer *k;
doublereal *xc, *yc;
integer *ier;
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double atan(), cos(), sin();

    /* Local variables */
    static doublereal a, c__;
    static integer i__;
    static doublereal s;
    static integer k2, k3;
    static doublereal x0, y0;
    static integer kk, np1;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   04/06/90 */

/*   This subroutine computes the coordinates of a sequence */
/* of N equally spaced points on the unit circle centered at */
/* (0,0).  An N-sided polygonal approximation to the circle */
/* may be plotted by connecting (XC(I),YC(I)) to (XC(I+1), */
/* YC(I+1)) for I = 1,...,N, where XC(N+1) = XC(1) and */
/* YC(N+1) = YC(1).  A reasonable value for N in this case */
/* is 2*PI*R, where R is the radius of the circle in device */
/* coordinates. */


/* On input: */

/*       K = Number of points in each quadrant, defining N as */
/*           4K.  K .GE. 1. */

/*       XC,YC = Arrays of length at least N+1 = 4K+1. */

/* K is not altered by this routine. */

/* On output: */

/*       XC,YC = Cartesian coordinates of the points on the */
/*               unit circle in the first N+1 locations. */
/*               XC(I) = cos(A*(I-1)), YC(I) = sin(A*(I-1)), */
/*               where A = 2*PI/N.  Note that XC(N+1) = XC(1) */
/*               and YC(N+1) = YC(1). */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered. */
/*             IER = 1 if K < 1 on input. */

/* Modules required by CIRCLE:  None */

/* Intrinsic functions called by CIRCLE:  ATAN, COS, DBLE, */
/*                                          SIN */

/* *********************************************************** */


/* Local parameters: */

/* I =     DO-loop index and index for XC and YC */
/* KK =    Local copy of K */
/* K2 =    K*2 */
/* K3 =    K*3 */
/* NP1 =   N+1 = 4*K + 1 */
/* A =     Angular separation between adjacent points */
/* C,S =   Cos(A) and sin(A), respectively, defining a */
/*           rotation through angle A */
/* X0,Y0 = Cartesian coordinates of a point on the unit */
/*           circle in the first quadrant */

    /* Parameter adjustments */
    --yc;
    --xc;

    /* Function Body */
    kk = *k;
    k2 = kk << 1;
    k3 = kk * 3;
    np1 = (kk << 2) + 1;

/* Test for invalid input, compute A, C, and S, and */
/*   initialize (X0,Y0) to (1,0). */

    if (kk < 1) {
	goto L2;
    }
    a = atan(1.) * 2. / (doublereal) kk;
    c__ = cos(a);
    s = sin(a);
    x0 = 1.;
    y0 = 0.;

/* Loop on points (X0,Y0) in the first quadrant, storing */
/*   the point and its reflections about the x axis, the */
/*   y axis, and the line y = -x. */

    i__1 = kk;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xc[i__] = x0;
	yc[i__] = y0;
	xc[i__ + kk] = -y0;
	yc[i__ + kk] = x0;
	xc[i__ + k2] = -x0;
	yc[i__ + k2] = -y0;
	xc[i__ + k3] = y0;
	yc[i__ + k3] = -x0;

/*   Rotate (X0,Y0) counterclockwise through angle A. */

	x0 = c__ * x0 - s * y0;
	y0 = s * x0 + c__ * y0;
/* L1: */
    }

/* Store the coordinates of the first point as the last */
/*   point. */

    xc[np1] = xc[1];
    yc[np1] = yc[1];
    *ier = 0;
    return 0;

/* K < 1. */

L2:
    *ier = 1;
    return 0;
} /* circle_ */

/* Subroutine */ int circum_(v1, v2, v3, c__, ier)
doublereal *v1, *v2, *v3, *c__;
integer *ier;
{
    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer i__;
    static doublereal cnorm, e1[3], e2[3], cu[3];


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   10/27/02 */

/*   This subroutine returns the circumcenter of a spherical */
/* triangle on the unit sphere:  the point on the sphere sur- */
/* face that is equally distant from the three triangle */
/* vertices and lies in the same hemisphere, where distance */
/* is taken to be arc-length on the sphere surface. */


/* On input: */

/*       V1,V2,V3 = Arrays of length 3 containing the Carte- */
/*                  sian coordinates of the three triangle */
/*                  vertices (unit vectors) in CCW order. */

/* The above parameters are not altered by this routine. */

/*       C = Array of length 3. */

/* On output: */

/*       C = Cartesian coordinates of the circumcenter unless */
/*           IER > 0, in which case C is not defined.  C = */
/*           (V2-V1) X (V3-V1) normalized to a unit vector. */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered. */
/*             IER = 1 if V1, V2, and V3 lie on a common */
/*                     line:  (V2-V1) X (V3-V1) = 0. */
/*             (The vertices are not tested for validity.) */

/* Modules required by CIRCUM:  None */

/* Intrinsic function called by CIRCUM:  SQRT */

/* *********************************************************** */


/* Local parameters: */

/* CNORM = Norm of CU:  used to compute C */
/* CU =    Scalar multiple of C:  E1 X E2 */
/* E1,E2 = Edges of the underlying planar triangle: */
/*           V2-V1 and V3-V1, respectively */
/* I =     DO-loop index */

    /* Parameter adjustments */
    --c__;
    --v3;
    --v2;
    --v1;

    /* Function Body */
    for (i__ = 1; i__ <= 3; ++i__) {
	e1[i__ - 1] = v2[i__] - v1[i__];
	e2[i__ - 1] = v3[i__] - v1[i__];
/* L1: */
    }

/* Compute CU = E1 X E2 and CNORM**2. */

    cu[0] = e1[1] * e2[2] - e1[2] * e2[1];
    cu[1] = e1[2] * e2[0] - e1[0] * e2[2];
    cu[2] = e1[0] * e2[1] - e1[1] * e2[0];
    cnorm = cu[0] * cu[0] + cu[1] * cu[1] + cu[2] * cu[2];

/* The vertices lie on a common line if and only if CU is */
/*   the zero vector. */

    if (cnorm != 0.) {

/*   No error:  compute C. */

	cnorm = sqrt(cnorm);
	for (i__ = 1; i__ <= 3; ++i__) {
	    c__[i__] = cu[i__ - 1] / cnorm;
/* L2: */
	}

/* If the vertices are nearly identical, the problem is */
/*   ill-conditioned and it is possible for the computed */
/*   value of C to be 180 degrees off:  <C,V1> near -1 */
/*   when it should be positive. */

	if (c__[1] * v1[1] + c__[2] * v1[2] + c__[3] * v1[3] < -.5) {
	    c__[1] = -c__[1];
	    c__[2] = -c__[2];
	    c__[3] = -c__[3];
	}
	*ier = 0;
    } else {

/*   CU = 0. */

	*ier = 1;
    }
    return 0;
} /* circum_ */

/* Subroutine */ int covsph_(kk, n0, list, lptr, lend, lnew)
integer *kk, *n0, *list, *lptr, *lend, *lnew;
{
    static integer lsav, next, k, lp;
    extern /* Subroutine */ int insert_();
    static integer nst;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/17/96 */

/*   This subroutine connects an exterior node KK to all */
/* boundary nodes of a triangulation of KK-1 points on the */
/* unit sphere, producing a triangulation that covers the */
/* sphere.  The data structure is updated with the addition */
/* of node KK, but no optimization is performed.  All boun- */
/* dary nodes must be visible from node KK. */


/* On input: */

/*       KK = Index of the node to be connected to the set of */
/*            all boundary nodes.  KK .GE. 4. */

/*       N0 = Index of a boundary node (in the range 1 to */
/*            KK-1).  N0 may be determined by Subroutine */
/*            TRFIND. */

/* The above parameters are not altered by this routine. */

/*       LIST,LPTR,LEND,LNEW = Triangulation data structure */
/*                             created by Subroutine TRMESH. */
/*                             Node N0 must be included in */
/*                             the triangulation. */

/* On output: */

/*       LIST,LPTR,LEND,LNEW = Data structure updated with */
/*                             the addition of node KK as the */
/*                             last entry.  The updated */
/*                             triangulation contains no */
/*                             boundary nodes. */

/* Module required by COVSPH:  INSERT */

/* *********************************************************** */


/* Local parameters: */

/* K =     Local copy of KK */
/* LP =    LIST pointer */
/* LSAV =  LIST pointer */
/* NEXT =  Boundary node visible from K */
/* NST =   Local copy of N0 */

    /* Parameter adjustments */
    --lend;
    --lptr;
    --list;

    /* Function Body */
    k = *kk;
    nst = *n0;

/* Traverse the boundary in clockwise order, inserting K as */
/*   the first neighbor of each boundary node, and converting */
/*   the boundary node to an interior node. */

    next = nst;
L1:
    lp = lend[next];
    insert_(&k, &lp, &list[1], &lptr[1], lnew);
    next = -list[lp];
    list[lp] = next;
    if (next != nst) {
	goto L1;
    }

/* Traverse the boundary again, adding each node to K's */
/*   adjacency list. */

    lsav = *lnew;
L2:
    lp = lend[next];
    list[*lnew] = next;
    lptr[*lnew] = *lnew + 1;
    ++(*lnew);
    next = list[lp];
    if (next != nst) {
	goto L2;
    }

    lptr[*lnew - 1] = lsav;
    lend[k] = *lnew - 1;
    return 0;
} /* covsph_ */

/* Subroutine */ int crlist_(n, ncol, x, y, z__, list, lend, lptr, lnew, ltri,
	 listc, nb, xc, yc, zc, rc, ier)
integer *n, *ncol;
doublereal *x, *y, *z__;
integer *list, *lend, *lptr, *lnew, *ltri, *listc, *nb;
doublereal *xc, *yc, *zc, *rc;
integer *ier;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double acos();

    /* Local variables */
    static integer ierr;
    static doublereal c__[3], t;
    static integer i1, i2, i3, i4, n0, n1, n2, n3, n4;
    static doublereal v1[3], v2[3], v3[3];
    static integer lp, kt, nn, nt;
    extern /* Subroutine */ int circum_();
    static integer nm2, kt1, kt2;
    extern integer lstptr_();
    extern logical swptst_();
    static integer kt11, kt12, kt21, kt22, lpl, lpn;
    static logical swp;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   03/05/03 */

/*   Given a Delaunay triangulation of nodes on the surface */
/* of the unit sphere, this subroutine returns the set of */
/* triangle circumcenters corresponding to Voronoi vertices, */
/* along with the circumradii and a list of triangle indexes */
/* LISTC stored in one-to-one correspondence with LIST/LPTR */
/* entries. */

/*   A triangle circumcenter is the point (unit vector) lying */
/* at the same angular distance from the three vertices and */
/* contained in the same hemisphere as the vertices.  (Note */
/* that the negative of a circumcenter is also equidistant */
/* from the vertices.)  If the triangulation covers the sur- */
/* face, the Voronoi vertices are the circumcenters of the */
/* triangles in the Delaunay triangulation.  LPTR, LEND, and */
/* LNEW are not altered in this case. */

/*   On the other hand, if the nodes are contained in a sin- */
/* gle hemisphere, the triangulation is implicitly extended */
/* to the entire surface by adding pseudo-arcs (of length */
/* greater than 180 degrees) between boundary nodes forming */
/* pseudo-triangles whose 'circumcenters' are included in the */
/* list.  This extension to the triangulation actually con- */
/* sists of a triangulation of the set of boundary nodes in */
/* which the swap test is reversed (a non-empty circumcircle */
/* test).  The negative circumcenters are stored as the */
/* pseudo-triangle 'circumcenters'.  LISTC, LPTR, LEND, and */
/* LNEW contain a data structure corresponding to the ex- */
/* tended triangulation (Voronoi diagram), but LIST is not */
/* altered in this case.  Thus, if it is necessary to retain */
/* the original (unextended) triangulation data structure, */
/* copies of LPTR and LNEW must be saved before calling this */
/* routine. */


/* On input: */

/*       N = Number of nodes in the triangulation.  N .GE. 3. */
/*           Note that, if N = 3, there are only two Voronoi */
/*           vertices separated by 180 degrees, and the */
/*           Voronoi regions are not well defined. */

/*       NCOL = Number of columns reserved for LTRI.  This */
/*              must be at least NB-2, where NB is the number */
/*              of boundary nodes. */

/*       X,Y,Z = Arrays of length N containing the Cartesian */
/*               coordinates of the nodes (unit vectors). */

/*       LIST = Integer array containing the set of adjacency */
/*              lists.  Refer to Subroutine TRMESH. */

/*       LEND = Set of pointers to ends of adjacency lists. */
/*              Refer to Subroutine TRMESH. */

/* The above parameters are not altered by this routine. */

/*       LPTR = Array of pointers associated with LIST.  Re- */
/*              fer to Subroutine TRMESH. */

/*       LNEW = Pointer to the first empty location in LIST */
/*              and LPTR (list length plus one). */

/*       LTRI = Integer work space array dimensioned 6 by */
/*              NCOL, or unused dummy parameter if NB = 0. */

/*       LISTC = Integer array of length at least 3*NT, where */
/*               NT = 2*N-4 is the number of triangles in the */
/*               triangulation (after extending it to cover */
/*               the entire surface if necessary). */

/*       XC,YC,ZC,RC = Arrays of length NT = 2*N-4. */

/* On output: */

/*       LPTR = Array of pointers associated with LISTC: */
/*              updated for the addition of pseudo-triangles */
/*              if the original triangulation contains */
/*              boundary nodes (NB > 0). */

/*       LNEW = Pointer to the first empty location in LISTC */
/*              and LPTR (list length plus one).  LNEW is not */
/*              altered if NB = 0. */

/*       LTRI = Triangle list whose first NB-2 columns con- */
/*              tain the indexes of a clockwise-ordered */
/*              sequence of vertices (first three rows) */
/*              followed by the LTRI column indexes of the */
/*              triangles opposite the vertices (or 0 */
/*              denoting the exterior region) in the last */
/*              three rows.  This array is not generally of */
/*              any use. */

/*       LISTC = Array containing triangle indexes (indexes */
/*               to XC, YC, ZC, and RC) stored in 1-1 corres- */
/*               pondence with LIST/LPTR entries (or entries */
/*               that would be stored in LIST for the */
/*               extended triangulation):  the index of tri- */
/*               angle (N1,N2,N3) is stored in LISTC(K), */
/*               LISTC(L), and LISTC(M), where LIST(K), */
/*               LIST(L), and LIST(M) are the indexes of N2 */
/*               as a neighbor of N1, N3 as a neighbor of N2, */
/*               and N1 as a neighbor of N3.  The Voronoi */
/*               region associated with a node is defined by */
/*               the CCW-ordered sequence of circumcenters in */
/*               one-to-one correspondence with its adjacency */
/*               list (in the extended triangulation). */

/*       NB = Number of boundary nodes unless IER = 1. */

/*       XC,YC,ZC = Arrays containing the Cartesian coordi- */
/*                  nates of the triangle circumcenters */
/*                  (Voronoi vertices).  XC(I)**2 + YC(I)**2 */
/*                  + ZC(I)**2 = 1.  The first NB-2 entries */
/*                  correspond to pseudo-triangles if NB > 0. */

/*       RC = Array containing circumradii (the arc lengths */
/*            or angles between the circumcenters and associ- */
/*            ated triangle vertices) in 1-1 correspondence */
/*            with circumcenters. */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered. */
/*             IER = 1 if N < 3. */
/*             IER = 2 if NCOL < NB-2. */
/*             IER = 3 if a triangle is degenerate (has ver- */
/*                     tices lying on a common geodesic). */

/* Modules required by CRLIST:  CIRCUM, LSTPTR, SWPTST */

/* Intrinsic functions called by CRLIST:  ABS, ACOS */

/* *********************************************************** */


/* Local parameters: */

/* C =         Circumcenter returned by Subroutine CIRCUM */
/* I1,I2,I3 =  Permutation of (1,2,3):  LTRI row indexes */
/* I4 =        LTRI row index in the range 1 to 3 */
/* IERR =      Error flag for calls to CIRCUM */
/* KT =        Triangle index */
/* KT1,KT2 =   Indexes of a pair of adjacent pseudo-triangles */
/* KT11,KT12 = Indexes of the pseudo-triangles opposite N1 */
/*               and N2 as vertices of KT1 */
/* KT21,KT22 = Indexes of the pseudo-triangles opposite N1 */
/*               and N2 as vertices of KT2 */
/* LP,LPN =    LIST pointers */
/* LPL =       LIST pointer of the last neighbor of N1 */
/* N0 =        Index of the first boundary node (initial */
/*               value of N1) in the loop on boundary nodes */
/*               used to store the pseudo-triangle indexes */
/*               in LISTC */
/* N1,N2,N3 =  Nodal indexes defining a triangle (CCW order) */
/*               or pseudo-triangle (clockwise order) */
/* N4 =        Index of the node opposite N2 -> N1 */
/* NM2 =       N-2 */
/* NN =        Local copy of N */
/* NT =        Number of pseudo-triangles:  NB-2 */
/* SWP =       Logical variable set to TRUE in each optimiza- */
/*               tion loop (loop on pseudo-arcs) iff a swap */
/*               is performed */
/* V1,V2,V3 =  Vertices of triangle KT = (N1,N2,N3) sent to */
/*               Subroutine CIRCUM */

    /* Parameter adjustments */
    --lend;
    --z__;
    --y;
    --x;
    ltri -= 7;
    --list;
    --lptr;
    --listc;
    --xc;
    --yc;
    --zc;
    --rc;

    /* Function Body */
    nn = *n;
    *nb = 0;
    nt = 0;
    if (nn < 3) {
	goto L21;
    }

/* Search for a boundary node N1. */

    i__1 = nn;
    for (n1 = 1; n1 <= i__1; ++n1) {
	lp = lend[n1];
	if (list[lp] < 0) {
	    goto L2;
	}
/* L1: */
    }

/* The triangulation already covers the sphere. */

    goto L9;

/* There are NB .GE. 3 boundary nodes.  Add NB-2 pseudo- */
/*   triangles (N1,N2,N3) by connecting N3 to the NB-3 */
/*   boundary nodes to which it is not already adjacent. */

/*   Set N3 and N2 to the first and last neighbors, */
/*     respectively, of N1. */

L2:
    n2 = -list[lp];
    lp = lptr[lp];
    n3 = list[lp];

/*   Loop on boundary arcs N1 -> N2 in clockwise order, */
/*     storing triangles (N1,N2,N3) in column NT of LTRI */
/*     along with the indexes of the triangles opposite */
/*     the vertices. */

L3:
    ++nt;
    if (nt <= *ncol) {
	ltri[nt * 6 + 1] = n1;
	ltri[nt * 6 + 2] = n2;
	ltri[nt * 6 + 3] = n3;
	ltri[nt * 6 + 4] = nt + 1;
	ltri[nt * 6 + 5] = nt - 1;
	ltri[nt * 6 + 6] = 0;
    }
    n1 = n2;
    lp = lend[n1];
    n2 = -list[lp];
    if (n2 != n3) {
	goto L3;
    }

    *nb = nt + 2;
    if (*ncol < nt) {
	goto L22;
    }
    ltri[nt * 6 + 4] = 0;
    if (nt == 1) {
	goto L7;
    }

/* Optimize the exterior triangulation (set of pseudo- */
/*   triangles) by applying swaps to the pseudo-arcs N1-N2 */
/*   (pairs of adjacent pseudo-triangles KT1 and KT2 > KT1). */
/*   The loop on pseudo-arcs is repeated until no swaps are */
/*   performed. */

L4:
    swp = FALSE_;
    i__1 = nt - 1;
    for (kt1 = 1; kt1 <= i__1; ++kt1) {
	for (i3 = 1; i3 <= 3; ++i3) {
	    kt2 = ltri[i3 + 3 + kt1 * 6];
	    if (kt2 <= kt1) {
		goto L5;
	    }

/*   The LTRI row indexes (I1,I2,I3) of triangle KT1 = */
/*     (N1,N2,N3) are a cyclical permutation of (1,2,3). */

	    if (i3 == 1) {
		i1 = 2;
		i2 = 3;
	    } else if (i3 == 2) {
		i1 = 3;
		i2 = 1;
	    } else {
		i1 = 1;
		i2 = 2;
	    }
	    n1 = ltri[i1 + kt1 * 6];
	    n2 = ltri[i2 + kt1 * 6];
	    n3 = ltri[i3 + kt1 * 6];

/*   KT2 = (N2,N1,N4) for N4 = LTRI(I,KT2), where */
/*     LTRI(I+3,KT2) = KT1. */

	    if (ltri[kt2 * 6 + 4] == kt1) {
		i4 = 1;
	    } else if (ltri[kt2 * 6 + 5] == kt1) {
		i4 = 2;
	    } else {
		i4 = 3;
	    }
	    n4 = ltri[i4 + kt2 * 6];

/*   The empty circumcircle test is reversed for the pseudo- */
/*     triangles.  The reversal is implicit in the clockwise */
/*     ordering of the vertices. */

	    if (! swptst_(&n1, &n2, &n3, &n4, &x[1], &y[1], &z__[1])) {
		goto L5;
	    }

/*   Swap arc N1-N2 for N3-N4.  KTij is the triangle opposite */
/*     Nj as a vertex of KTi. */

	    swp = TRUE_;
	    kt11 = ltri[i1 + 3 + kt1 * 6];
	    kt12 = ltri[i2 + 3 + kt1 * 6];
	    if (i4 == 1) {
		i2 = 2;
		i1 = 3;
	    } else if (i4 == 2) {
		i2 = 3;
		i1 = 1;
	    } else {
		i2 = 1;
		i1 = 2;
	    }
	    kt21 = ltri[i1 + 3 + kt2 * 6];
	    kt22 = ltri[i2 + 3 + kt2 * 6];
	    ltri[kt1 * 6 + 1] = n4;
	    ltri[kt1 * 6 + 2] = n3;
	    ltri[kt1 * 6 + 3] = n1;
	    ltri[kt1 * 6 + 4] = kt12;
	    ltri[kt1 * 6 + 5] = kt22;
	    ltri[kt1 * 6 + 6] = kt2;
	    ltri[kt2 * 6 + 1] = n3;
	    ltri[kt2 * 6 + 2] = n4;
	    ltri[kt2 * 6 + 3] = n2;
	    ltri[kt2 * 6 + 4] = kt21;
	    ltri[kt2 * 6 + 5] = kt11;
	    ltri[kt2 * 6 + 6] = kt1;

/*   Correct the KT11 and KT22 entries that changed. */

	    if (kt11 != 0) {
		i4 = 4;
		if (ltri[kt11 * 6 + 4] != kt1) {
		    i4 = 5;
		    if (ltri[kt11 * 6 + 5] != kt1) {
			i4 = 6;
		    }
		}
		ltri[i4 + kt11 * 6] = kt2;
	    }
	    if (kt22 != 0) {
		i4 = 4;
		if (ltri[kt22 * 6 + 4] != kt2) {
		    i4 = 5;
		    if (ltri[kt22 * 6 + 5] != kt2) {
			i4 = 6;
		    }
		}
		ltri[i4 + kt22 * 6] = kt1;
	    }
L5:
	    ;
	}
/* L6: */
    }
    if (swp) {
	goto L4;
    }

/* Compute and store the negative circumcenters and radii of */
/*   the pseudo-triangles in the first NT positions. */

L7:
    i__1 = nt;
    for (kt = 1; kt <= i__1; ++kt) {
	n1 = ltri[kt * 6 + 1];
	n2 = ltri[kt * 6 + 2];
	n3 = ltri[kt * 6 + 3];
	v1[0] = x[n1];
	v1[1] = y[n1];
	v1[2] = z__[n1];
	v2[0] = x[n2];
	v2[1] = y[n2];
	v2[2] = z__[n2];
	v3[0] = x[n3];
	v3[1] = y[n3];
	v3[2] = z__[n3];
	circum_(v2, v1, v3, c__, &ierr);
	if (ierr != 0) {
	    goto L23;
	}

/*   Store the negative circumcenter and radius (computed */
/*     from <V1,C>). */

	xc[kt] = -c__[0];
	yc[kt] = -c__[1];
	zc[kt] = -c__[2];
	t = -(v1[0] * c__[0] + v1[1] * c__[1] + v1[2] * c__[2]);
	if (t < -1.) {
	    t = -1.;
	}
	if (t > 1.) {
	    t = 1.;
	}
	rc[kt] = acos(t);
/* L8: */
    }

/* Compute and store the circumcenters and radii of the */
/*   actual triangles in positions KT = NT+1, NT+2, ... */
/*   Also, store the triangle indexes KT in the appropriate */
/*   LISTC positions. */

L9:
    kt = nt;

/*   Loop on nodes N1. */

    nm2 = nn - 2;
    i__1 = nm2;
    for (n1 = 1; n1 <= i__1; ++n1) {
	lpl = lend[n1];
	lp = lpl;
	n3 = list[lp];

/*   Loop on adjacent neighbors N2,N3 of N1 for which N2 > N1 */
/*     and N3 > N1. */

L10:
	lp = lptr[lp];
	n2 = n3;
	n3 = (i__2 = list[lp], abs(i__2));
	if (n2 <= n1 || n3 <= n1) {
	    goto L11;
	}
	++kt;

/*   Compute the circumcenter C of triangle KT = (N1,N2,N3). */

	v1[0] = x[n1];
	v1[1] = y[n1];
	v1[2] = z__[n1];
	v2[0] = x[n2];
	v2[1] = y[n2];
	v2[2] = z__[n2];
	v3[0] = x[n3];
	v3[1] = y[n3];
	v3[2] = z__[n3];
	circum_(v1, v2, v3, c__, &ierr);
	if (ierr != 0) {
	    goto L23;
	}

/*   Store the circumcenter, radius and triangle index. */

	xc[kt] = c__[0];
	yc[kt] = c__[1];
	zc[kt] = c__[2];
	t = v1[0] * c__[0] + v1[1] * c__[1] + v1[2] * c__[2];
	if (t < -1.) {
	    t = -1.;
	}
	if (t > 1.) {
	    t = 1.;
	}
	rc[kt] = acos(t);

/*   Store KT in LISTC(LPN), where Abs(LIST(LPN)) is the */
/*     index of N2 as a neighbor of N1, N3 as a neighbor */
/*     of N2, and N1 as a neighbor of N3. */

	lpn = lstptr_(&lpl, &n2, &list[1], &lptr[1]);
	listc[lpn] = kt;
	lpn = lstptr_(&lend[n2], &n3, &list[1], &lptr[1]);
	listc[lpn] = kt;
	lpn = lstptr_(&lend[n3], &n1, &list[1], &lptr[1]);
	listc[lpn] = kt;
L11:
	if (lp != lpl) {
	    goto L10;
	}
/* L12: */
    }
    if (nt == 0) {
	goto L20;
    }

/* Store the first NT triangle indexes in LISTC. */

/*   Find a boundary triangle KT1 = (N1,N2,N3) with a */
/*     boundary arc opposite N3. */

    kt1 = 0;
L13:
    ++kt1;
    if (ltri[kt1 * 6 + 4] == 0) {
	i1 = 2;
	i2 = 3;
	i3 = 1;
	goto L14;
    } else if (ltri[kt1 * 6 + 5] == 0) {
	i1 = 3;
	i2 = 1;
	i3 = 2;
	goto L14;
    } else if (ltri[kt1 * 6 + 6] == 0) {
	i1 = 1;
	i2 = 2;
	i3 = 3;
	goto L14;
    }
    goto L13;
L14:
    n1 = ltri[i1 + kt1 * 6];
    n0 = n1;

/*   Loop on boundary nodes N1 in CCW order, storing the */
/*     indexes of the clockwise-ordered sequence of triangles */
/*     that contain N1.  The first triangle overwrites the */
/*     last neighbor position, and the remaining triangles, */
/*     if any, are appended to N1's adjacency list. */

/*   A pointer to the first neighbor of N1 is saved in LPN. */

L15:
    lp = lend[n1];
    lpn = lptr[lp];
    listc[lp] = kt1;

/*   Loop on triangles KT2 containing N1. */

L16:
    kt2 = ltri[i2 + 3 + kt1 * 6];
    if (kt2 != 0) {

/*   Append KT2 to N1's triangle list. */

	lptr[lp] = *lnew;
	lp = *lnew;
	listc[lp] = kt2;
	++(*lnew);

/*   Set KT1 to KT2 and update (I1,I2,I3) such that */
/*     LTRI(I1,KT1) = N1. */

	kt1 = kt2;
	if (ltri[kt1 * 6 + 1] == n1) {
	    i1 = 1;
	    i2 = 2;
	    i3 = 3;
	} else if (ltri[kt1 * 6 + 2] == n1) {
	    i1 = 2;
	    i2 = 3;
	    i3 = 1;
	} else {
	    i1 = 3;
	    i2 = 1;
	    i3 = 2;
	}
	goto L16;
    }

/*   Store the saved first-triangle pointer in LPTR(LP), set */
/*     N1 to the next boundary node, test for termination, */
/*     and permute the indexes:  the last triangle containing */
/*     a boundary node is the first triangle containing the */
/*     next boundary node. */

    lptr[lp] = lpn;
    n1 = ltri[i3 + kt1 * 6];
    if (n1 != n0) {
	i4 = i3;
	i3 = i2;
	i2 = i1;
	i1 = i4;
	goto L15;
    }

/* No errors encountered. */

L20:
    *ier = 0;
    return 0;

/* N < 3. */

L21:
    *ier = 1;
    return 0;

/* Insufficient space reserved for LTRI. */

L22:
    *ier = 2;
    return 0;

/* Error flag returned by CIRCUM: KT indexes a null triangle. */

L23:
    *ier = 3;
    return 0;
} /* crlist_ */

/* Subroutine */ int delarc_(n, io1, io2, list, lptr, lend, lnew, ier)
integer *n, *io1, *io2, *list, *lptr, *lend, *lnew, *ier;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int delnb_();
    static integer n1, n2, n3, lp;
    extern integer lstptr_();
    static integer lph, lpl;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/17/96 */

/*   This subroutine deletes a boundary arc from a triangula- */
/* tion.  It may be used to remove a null triangle from the */
/* convex hull boundary.  Note, however, that if the union of */
/* triangles is rendered nonconvex, Subroutines DELNOD, EDGE, */
/* and TRFIND (and hence ADDNOD) may fail.  Also, Function */
/* NEARND should not be called following an arc deletion. */

/*   This routine is identical to the similarly named routine */
/* in TRIPACK. */


/* On input: */

/*       N = Number of nodes in the triangulation.  N .GE. 4. */

/*       IO1,IO2 = Indexes (in the range 1 to N) of a pair of */
/*                 adjacent boundary nodes defining the arc */
/*                 to be removed. */

/* The above parameters are not altered by this routine. */

/*       LIST,LPTR,LEND,LNEW = Triangulation data structure */
/*                             created by Subroutine TRMESH. */

/* On output: */

/*       LIST,LPTR,LEND,LNEW = Data structure updated with */
/*                             the removal of arc IO1-IO2 */
/*                             unless IER > 0. */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered. */
/*             IER = 1 if N, IO1, or IO2 is outside its valid */
/*                     range, or IO1 = IO2. */
/*             IER = 2 if IO1-IO2 is not a boundary arc. */
/*             IER = 3 if the node opposite IO1-IO2 is al- */
/*                     ready a boundary node, and thus IO1 */
/*                     or IO2 has only two neighbors or a */
/*                     deletion would result in two triangu- */
/*                     lations sharing a single node. */
/*             IER = 4 if one of the nodes is a neighbor of */
/*                     the other, but not vice versa, imply- */
/*                     ing an invalid triangulation data */
/*                     structure. */

/* Module required by DELARC:  DELNB, LSTPTR */

/* Intrinsic function called by DELARC:  ABS */

/* *********************************************************** */


/* Local parameters: */

/* LP =       LIST pointer */
/* LPH =      LIST pointer or flag returned by DELNB */
/* LPL =      Pointer to the last neighbor of N1, N2, or N3 */
/* N1,N2,N3 = Nodal indexes of a triangle such that N1->N2 */
/*              is the directed boundary edge associated */
/*              with IO1-IO2 */

    /* Parameter adjustments */
    --lend;
    --list;
    --lptr;

    /* Function Body */
    n1 = *io1;
    n2 = *io2;

/* Test for errors, and set N1->N2 to the directed boundary */
/*   edge associated with IO1-IO2:  (N1,N2,N3) is a triangle */
/*   for some N3. */

    if (*n < 4 || n1 < 1 || n1 > *n || n2 < 1 || n2 > *n || n1 == n2) {
	*ier = 1;
	return 0;
    }

    lpl = lend[n2];
    if (-list[lpl] != n1) {
	n1 = n2;
	n2 = *io1;
	lpl = lend[n2];
	if (-list[lpl] != n1) {
	    *ier = 2;
	    return 0;
	}
    }

/* Set N3 to the node opposite N1->N2 (the second neighbor */
/*   of N1), and test for error 3 (N3 already a boundary */
/*   node). */

    lpl = lend[n1];
    lp = lptr[lpl];
    lp = lptr[lp];
    n3 = (i__1 = list[lp], abs(i__1));
    lpl = lend[n3];
    if (list[lpl] <= 0) {
	*ier = 3;
	return 0;
    }

/* Delete N2 as a neighbor of N1, making N3 the first */
/*   neighbor, and test for error 4 (N2 not a neighbor */
/*   of N1).  Note that previously computed pointers may */
/*   no longer be valid following the call to DELNB. */

    delnb_(&n1, &n2, n, &list[1], &lptr[1], &lend[1], lnew, &lph);
    if (lph < 0) {
	*ier = 4;
	return 0;
    }

/* Delete N1 as a neighbor of N2, making N3 the new last */
/*   neighbor. */

    delnb_(&n2, &n1, n, &list[1], &lptr[1], &lend[1], lnew, &lph);

/* Make N3 a boundary node with first neighbor N2 and last */
/*   neighbor N1. */

    lp = lstptr_(&lend[n3], &n1, &list[1], &lptr[1]);
    lend[n3] = lp;
    list[lp] = -n1;

/* No errors encountered. */

    *ier = 0;
    return 0;
} /* delarc_ */

/* Subroutine */ int delnb_(n0, nb, n, list, lptr, lend, lnew, lph)
integer *n0, *nb, *n, *list, *lptr, *lend, *lnew, *lph;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, lp, nn, lpb, lpl, lpp, lnw;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/29/98 */

/*   This subroutine deletes a neighbor NB from the adjacency */
/* list of node N0 (but N0 is not deleted from the adjacency */
/* list of NB) and, if NB is a boundary node, makes N0 a */
/* boundary node.  For pointer (LIST index) LPH to NB as a */
/* neighbor of N0, the empty LIST,LPTR location LPH is filled */
/* in with the values at LNEW-1, pointer LNEW-1 (in LPTR and */
/* possibly in LEND) is changed to LPH, and LNEW is decremen- */
/* ted.  This requires a search of LEND and LPTR entailing an */
/* expected operation count of O(N). */

/*   This routine is identical to the similarly named routine */
/* in TRIPACK. */


/* On input: */

/*       N0,NB = Indexes, in the range 1 to N, of a pair of */
/*               nodes such that NB is a neighbor of N0. */
/*               (N0 need not be a neighbor of NB.) */

/*       N = Number of nodes in the triangulation.  N .GE. 3. */

/* The above parameters are not altered by this routine. */

/*       LIST,LPTR,LEND,LNEW = Data structure defining the */
/*                             triangulation. */

/* On output: */

/*       LIST,LPTR,LEND,LNEW = Data structure updated with */
/*                             the removal of NB from the ad- */
/*                             jacency list of N0 unless */
/*                             LPH < 0. */

/*       LPH = List pointer to the hole (NB as a neighbor of */
/*             N0) filled in by the values at LNEW-1 or error */
/*             indicator: */
/*             LPH > 0 if no errors were encountered. */
/*             LPH = -1 if N0, NB, or N is outside its valid */
/*                      range. */
/*             LPH = -2 if NB is not a neighbor of N0. */

/* Modules required by DELNB:  None */

/* Intrinsic function called by DELNB:  ABS */

/* *********************************************************** */


/* Local parameters: */

/* I =   DO-loop index */
/* LNW = LNEW-1 (output value of LNEW) */
/* LP =  LIST pointer of the last neighbor of NB */
/* LPB = Pointer to NB as a neighbor of N0 */
/* LPL = Pointer to the last neighbor of N0 */
/* LPP = Pointer to the neighbor of N0 that precedes NB */
/* NN =  Local copy of N */

    /* Parameter adjustments */
    --lend;
    --list;
    --lptr;

    /* Function Body */
    nn = *n;

/* Test for error 1. */

    if (*n0 < 1 || *n0 > nn || *nb < 1 || *nb > nn || nn < 3) {
	*lph = -1;
	return 0;
    }

/*   Find pointers to neighbors of N0: */

/*     LPL points to the last neighbor, */
/*     LPP points to the neighbor NP preceding NB, and */
/*     LPB points to NB. */

    lpl = lend[*n0];
    lpp = lpl;
    lpb = lptr[lpp];
L1:
    if (list[lpb] == *nb) {
	goto L2;
    }
    lpp = lpb;
    lpb = lptr[lpp];
    if (lpb != lpl) {
	goto L1;
    }

/*   Test for error 2 (NB not found). */

    if ((i__1 = list[lpb], abs(i__1)) != *nb) {
	*lph = -2;
	return 0;
    }

/*   NB is the last neighbor of N0.  Make NP the new last */
/*     neighbor and, if NB is a boundary node, then make N0 */
/*     a boundary node. */

    lend[*n0] = lpp;
    lp = lend[*nb];
    if (list[lp] < 0) {
	list[lpp] = -list[lpp];
    }
    goto L3;

/*   NB is not the last neighbor of N0.  If NB is a boundary */
/*     node and N0 is not, then make N0 a boundary node with */
/*     last neighbor NP. */

L2:
    lp = lend[*nb];
    if (list[lp] < 0 && list[lpl] > 0) {
	lend[*n0] = lpp;
	list[lpp] = -list[lpp];
    }

/*   Update LPTR so that the neighbor following NB now fol- */
/*     lows NP, and fill in the hole at location LPB. */

L3:
    lptr[lpp] = lptr[lpb];
    lnw = *lnew - 1;
    list[lpb] = list[lnw];
    lptr[lpb] = lptr[lnw];
    for (i__ = nn; i__ >= 1; --i__) {
	if (lend[i__] == lnw) {
	    lend[i__] = lpb;
	    goto L5;
	}
/* L4: */
    }

L5:
    i__1 = lnw - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (lptr[i__] == lnw) {
	    lptr[i__] = lpb;
	}
/* L6: */
    }

/* No errors encountered. */

    *lnew = lnw;
    *lph = lpb;
    return 0;
} /* delnb_ */

/* Subroutine */ int delnod_(k, n, x, y, z__, list, lptr, lend, lnew, lwk, 
	iwk, ier)
integer *k, *n;
doublereal *x, *y, *z__;
integer *list, *lptr, *lend, *lnew, *lwk, *iwk, *ier;
{
    /* Format strings */
    static char fmt_100[] = "(//5x,\002*** Error in OPTIM (called from \002\
,\002DELNOD):  NIT = \002,i4,\002, IER = \002,i1,\002 ***\002/)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    extern logical left_();
    static logical bdry;
    static integer ierr, lwkl;
    extern /* Subroutine */ int swap_();
    static integer i__, j;
    extern /* Subroutine */ int delnb_();
    extern integer nbcnt_();
    extern /* Subroutine */ int optim_();
    static integer nfrst, n1, n2;
    static doublereal x1, x2, y1, y2, z1, z2;
    static integer nl, lp, nn, nr;
    static doublereal xl, yl, zl, xr, yr, zr;
    extern integer lstptr_();
    static integer nnb, lp21, lpf, lph, lpl, lpn, iwl, nit, lnw, lpl2;

    /* Fortran I/O blocks */
    static cilist io___194 = { 0, 6, 0, fmt_100, 0 };



/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   11/30/99 */

/*   This subroutine deletes node K (along with all arcs */
/* incident on node K) from a triangulation of N nodes on the */
/* unit sphere, and inserts arcs as necessary to produce a */
/* triangulation of the remaining N-1 nodes.  If a Delaunay */
/* triangulation is input, a Delaunay triangulation will */
/* result, and thus, DELNOD reverses the effect of a call to */
/* Subroutine ADDNOD. */


/* On input: */

/*       K = Index (for X, Y, and Z) of the node to be */
/*           deleted.  1 .LE. K .LE. N. */

/* K is not altered by this routine. */

/*       N = Number of nodes in the triangulation on input. */
/*           N .GE. 4.  Note that N will be decremented */
/*           following the deletion. */

/*       X,Y,Z = Arrays of length N containing the Cartesian */
/*               coordinates of the nodes in the triangula- */
/*               tion. */

/*       LIST,LPTR,LEND,LNEW = Data structure defining the */
/*                             triangulation.  Refer to Sub- */
/*                             routine TRMESH. */

/*       LWK = Number of columns reserved for IWK.  LWK must */
/*             be at least NNB-3, where NNB is the number of */
/*             neighbors of node K, including an extra */
/*             pseudo-node if K is a boundary node. */

/*       IWK = Integer work array dimensioned 2 by LWK (or */
/*             array of length .GE. 2*LWK). */

/* On output: */

/*       N = Number of nodes in the triangulation on output. */
/*           The input value is decremented unless 1 .LE. IER */
/*           .LE. 4. */

/*       X,Y,Z = Updated arrays containing nodal coordinates */
/*               (with elements K+1,...,N+1 shifted up one */
/*               position, thus overwriting element K) unless */
/*               1 .LE. IER .LE. 4. */

/*       LIST,LPTR,LEND,LNEW = Updated triangulation data */
/*                             structure reflecting the dele- */
/*                             tion unless 1 .LE. IER .LE. 4. */
/*                             Note that the data structure */
/*                             may have been altered if IER > */
/*                             3. */

/*       LWK = Number of IWK columns required unless IER = 1 */
/*             or IER = 3. */

/*       IWK = Indexes of the endpoints of the new arcs added */
/*             unless LWK = 0 or 1 .LE. IER .LE. 4.  (Arcs */
/*             are associated with columns, or pairs of */
/*             adjacent elements if IWK is declared as a */
/*             singly-subscripted array.) */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered. */
/*             IER = 1 if K or N is outside its valid range */
/*                     or LWK < 0 on input. */
/*             IER = 2 if more space is required in IWK. */
/*                     Refer to LWK. */
/*             IER = 3 if the triangulation data structure is */
/*                     invalid on input. */
/*             IER = 4 if K indexes an interior node with */
/*                     four or more neighbors, none of which */
/*                     can be swapped out due to collineari- */
/*                     ty, and K cannot therefore be deleted. */
/*             IER = 5 if an error flag (other than IER = 1) */
/*                     was returned by OPTIM.  An error */
/*                     message is written to the standard */
/*                     output unit in this case. */
/*             IER = 6 if error flag 1 was returned by OPTIM. */
/*                     This is not necessarily an error, but */
/*                     the arcs may not be optimal. */

/*   Note that the deletion may result in all remaining nodes */
/* being collinear.  This situation is not flagged. */

/* Modules required by DELNOD:  DELNB, LEFT, LSTPTR, NBCNT, */
/*                                OPTIM, SWAP, SWPTST */

/* Intrinsic function called by DELNOD:  ABS */

/* *********************************************************** */


/* Local parameters: */

/* BDRY =    Logical variable with value TRUE iff N1 is a */
/*             boundary node */
/* I,J =     DO-loop indexes */
/* IERR =    Error flag returned by OPTIM */
/* IWL =     Number of IWK columns containing arcs */
/* LNW =     Local copy of LNEW */
/* LP =      LIST pointer */
/* LP21 =    LIST pointer returned by SWAP */
/* LPF,LPL = Pointers to the first and last neighbors of N1 */
/* LPH =     Pointer (or flag) returned by DELNB */
/* LPL2 =    Pointer to the last neighbor of N2 */
/* LPN =     Pointer to a neighbor of N1 */
/* LWKL =    Input value of LWK */
/* N1 =      Local copy of K */
/* N2 =      Neighbor of N1 */
/* NFRST =   First neighbor of N1:  LIST(LPF) */
/* NIT =     Number of iterations in OPTIM */
/* NR,NL =   Neighbors of N1 preceding (to the right of) and */
/*             following (to the left of) N2, respectively */
/* NN =      Number of nodes in the triangulation */
/* NNB =     Number of neighbors of N1 (including a pseudo- */
/*             node representing the boundary if N1 is a */
/*             boundary node) */
/* X1,Y1,Z1 = Coordinates of N1 */
/* X2,Y2,Z2 = Coordinates of N2 */
/* XL,YL,ZL = Coordinates of NL */
/* XR,YR,ZR = Coordinates of NR */


/* Set N1 to K and NNB to the number of neighbors of N1 (plus */
/*   one if N1 is a boundary node), and test for errors.  LPF */
/*   and LPL are LIST indexes of the first and last neighbors */
/*   of N1, IWL is the number of IWK columns containing arcs, */
/*   and BDRY is TRUE iff N1 is a boundary node. */

    /* Parameter adjustments */
    iwk -= 3;
    --lend;
    --lptr;
    --list;
    --z__;
    --y;
    --x;

    /* Function Body */
    n1 = *k;
    nn = *n;
    if (n1 < 1 || n1 > nn || nn < 4 || *lwk < 0) {
	goto L21;
    }
    lpl = lend[n1];
    lpf = lptr[lpl];
    nnb = nbcnt_(&lpl, &lptr[1]);
    bdry = list[lpl] < 0;
    if (bdry) {
	++nnb;
    }
    if (nnb < 3) {
	goto L23;
    }
    lwkl = *lwk;
    *lwk = nnb - 3;
    if (lwkl < *lwk) {
	goto L22;
    }
    iwl = 0;
    if (nnb == 3) {
	goto L3;
    }

/* Initialize for loop on arcs N1-N2 for neighbors N2 of N1, */
/*   beginning with the second neighbor.  NR and NL are the */
/*   neighbors preceding and following N2, respectively, and */
/*   LP indexes NL.  The loop is exited when all possible */
/*   swaps have been applied to arcs incident on N1. */

    x1 = x[n1];
    y1 = y[n1];
    z1 = z__[n1];
    nfrst = list[lpf];
    nr = nfrst;
    xr = x[nr];
    yr = y[nr];
    zr = z__[nr];
    lp = lptr[lpf];
    n2 = list[lp];
    x2 = x[n2];
    y2 = y[n2];
    z2 = z__[n2];
    lp = lptr[lp];

/* Top of loop:  set NL to the neighbor following N2. */

L1:
    nl = (i__1 = list[lp], abs(i__1));
    if (nl == nfrst && bdry) {
	goto L3;
    }
    xl = x[nl];
    yl = y[nl];
    zl = z__[nl];

/*   Test for a convex quadrilateral.  To avoid an incorrect */
/*     test caused by collinearity, use the fact that if N1 */
/*     is a boundary node, then N1 LEFT NR->NL and if N2 is */
/*     a boundary node, then N2 LEFT NL->NR. */

    lpl2 = lend[n2];
    if (! ((bdry || left_(&xr, &yr, &zr, &xl, &yl, &zl, &x1, &y1, &z1)) && (
	    list[lpl2] < 0 || left_(&xl, &yl, &zl, &xr, &yr, &zr, &x2, &y2, &
	    z2)))) {

/*   Nonconvex quadrilateral -- no swap is possible. */

	nr = n2;
	xr = x2;
	yr = y2;
	zr = z2;
	goto L2;
    }

/*   The quadrilateral defined by adjacent triangles */
/*     (N1,N2,NL) and (N2,N1,NR) is convex.  Swap in */
/*     NL-NR and store it in IWK unless NL and NR are */
/*     already adjacent, in which case the swap is not */
/*     possible.  Indexes larger than N1 must be decremented */
/*     since N1 will be deleted from X, Y, and Z. */

    swap_(&nl, &nr, &n1, &n2, &list[1], &lptr[1], &lend[1], &lp21);
    if (lp21 == 0) {
	nr = n2;
	xr = x2;
	yr = y2;
	zr = z2;
	goto L2;
    }
    ++iwl;
    if (nl <= n1) {
	iwk[(iwl << 1) + 1] = nl;
    } else {
	iwk[(iwl << 1) + 1] = nl - 1;
    }
    if (nr <= n1) {
	iwk[(iwl << 1) + 2] = nr;
    } else {
	iwk[(iwl << 1) + 2] = nr - 1;
    }

/*   Recompute the LIST indexes and NFRST, and decrement NNB. */

    lpl = lend[n1];
    --nnb;
    if (nnb == 3) {
	goto L3;
    }
    lpf = lptr[lpl];
    nfrst = list[lpf];
    lp = lstptr_(&lpl, &nl, &list[1], &lptr[1]);
    if (nr == nfrst) {
	goto L2;
    }

/*   NR is not the first neighbor of N1. */
/*     Back up and test N1-NR for a swap again:  Set N2 to */
/*     NR and NR to the previous neighbor of N1 -- the */
/*     neighbor of NR which follows N1.  LP21 points to NL */
/*     as a neighbor of NR. */

    n2 = nr;
    x2 = xr;
    y2 = yr;
    z2 = zr;
    lp21 = lptr[lp21];
    lp21 = lptr[lp21];
    nr = (i__1 = list[lp21], abs(i__1));
    xr = x[nr];
    yr = y[nr];
    zr = z__[nr];
    goto L1;

/*   Bottom of loop -- test for termination of loop. */

L2:
    if (n2 == nfrst) {
	goto L3;
    }
    n2 = nl;
    x2 = xl;
    y2 = yl;
    z2 = zl;
    lp = lptr[lp];
    goto L1;

/* Delete N1 and all its incident arcs.  If N1 is an interior */
/*   node and either NNB > 3 or NNB = 3 and N2 LEFT NR->NL, */
/*   then N1 must be separated from its neighbors by a plane */
/*   containing the origin -- its removal reverses the effect */
/*   of a call to COVSPH, and all its neighbors become */
/*   boundary nodes.  This is achieved by treating it as if */
/*   it were a boundary node (setting BDRY to TRUE, changing */
/*   a sign in LIST, and incrementing NNB). */

L3:
    if (! bdry) {
	if (nnb > 3) {
	    bdry = TRUE_;
	} else {
	    lpf = lptr[lpl];
	    nr = list[lpf];
	    lp = lptr[lpf];
	    n2 = list[lp];
	    nl = list[lpl];
	    bdry = left_(&x[nr], &y[nr], &z__[nr], &x[nl], &y[nl], &z__[nl], &
		    x[n2], &y[n2], &z__[n2]);
	}
	if (bdry) {

/*   IF a boundary node already exists, then N1 and its */
/*     neighbors cannot be converted to boundary nodes. */
/*     (They must be collinear.)  This is a problem if */
/*     NNB > 3. */

	    i__1 = nn;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (list[lend[i__]] < 0) {
		    bdry = FALSE_;
		    goto L5;
		}
/* L4: */
	    }
	    list[lpl] = -list[lpl];
	    ++nnb;
	}
    }
L5:
    if (! bdry && nnb > 3) {
	goto L24;
    }

/* Initialize for loop on neighbors.  LPL points to the last */
/*   neighbor of N1.  LNEW is stored in local variable LNW. */

    lp = lpl;
    lnw = *lnew;

/* Loop on neighbors N2 of N1, beginning with the first. */

L6:
    lp = lptr[lp];
    n2 = (i__1 = list[lp], abs(i__1));
    delnb_(&n2, &n1, n, &list[1], &lptr[1], &lend[1], &lnw, &lph);
    if (lph < 0) {
	goto L23;
    }

/*   LP and LPL may require alteration. */

    if (lpl == lnw) {
	lpl = lph;
    }
    if (lp == lnw) {
	lp = lph;
    }
    if (lp != lpl) {
	goto L6;
    }

/* Delete N1 from X, Y, Z, and LEND, and remove its adjacency */
/*   list from LIST and LPTR.  LIST entries (nodal indexes) */
/*   which are larger than N1 must be decremented. */

    --nn;
    if (n1 > nn) {
	goto L9;
    }
    i__1 = nn;
    for (i__ = n1; i__ <= i__1; ++i__) {
	x[i__] = x[i__ + 1];
	y[i__] = y[i__ + 1];
	z__[i__] = z__[i__ + 1];
	lend[i__] = lend[i__ + 1];
/* L7: */
    }

    i__1 = lnw - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (list[i__] > n1) {
	    --list[i__];
	}
	if (list[i__] < -n1) {
	    ++list[i__];
	}
/* L8: */
    }

/*   For LPN = first to last neighbors of N1, delete the */
/*     preceding neighbor (indexed by LP). */

/*   Each empty LIST,LPTR location LP is filled in with the */
/*     values at LNW-1, and LNW is decremented.  All pointers */
/*     (including those in LPTR and LEND) with value LNW-1 */
/*     must be changed to LP. */

/*  LPL points to the last neighbor of N1. */

L9:
    if (bdry) {
	--nnb;
    }
    lpn = lpl;
    i__1 = nnb;
    for (j = 1; j <= i__1; ++j) {
	--lnw;
	lp = lpn;
	lpn = lptr[lp];
	list[lp] = list[lnw];
	lptr[lp] = lptr[lnw];
	if (lptr[lpn] == lnw) {
	    lptr[lpn] = lp;
	}
	if (lpn == lnw) {
	    lpn = lp;
	}
	for (i__ = nn; i__ >= 1; --i__) {
	    if (lend[i__] == lnw) {
		lend[i__] = lp;
		goto L11;
	    }
/* L10: */
	}

L11:
	for (i__ = lnw - 1; i__ >= 1; --i__) {
	    if (lptr[i__] == lnw) {
		lptr[i__] = lp;
	    }
/* L12: */
	}
/* L13: */
    }

/* Update N and LNEW, and optimize the patch of triangles */
/*   containing K (on input) by applying swaps to the arcs */
/*   in IWK. */

    *n = nn;
    *lnew = lnw;
    if (iwl > 0) {
	nit = iwl << 2;
	optim_(&x[1], &y[1], &z__[1], &iwl, &list[1], &lptr[1], &lend[1], &
		nit, &iwk[3], &ierr);
	if (ierr != 0 && ierr != 1) {
	    goto L25;
	}
	if (ierr == 1) {
	    goto L26;
	}
    }

/* Successful termination. */

    *ier = 0;
    return 0;

/* Invalid input parameter. */

L21:
    *ier = 1;
    return 0;

/* Insufficient space reserved for IWK. */

L22:
    *ier = 2;
    return 0;

/* Invalid triangulation data structure.  NNB < 3 on input or */
/*   N2 is a neighbor of N1 but N1 is not a neighbor of N2. */

L23:
    *ier = 3;
    return 0;

/* N1 is interior but NNB could not be reduced to 3. */

L24:
    *ier = 4;
    return 0;

/* Error flag (other than 1) returned by OPTIM. */

L25:
    *ier = 5;
    s_wsfe(&io___194);
    do_fio(&c__1, (char *)&nit, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ierr, (ftnlen)sizeof(integer));
    e_wsfe();
    return 0;

/* Error flag 1 returned by OPTIM. */

L26:
    *ier = 6;
    return 0;
} /* delnod_ */

/* Subroutine */ int drwarc_(lun, p, q, tol, nseg)
integer *lun;
doublereal *p, *q, *tol;
integer *nseg;
{
    /* Format strings */
    static char fmt_100[] = "(2f12.6,\002 moveto\002,2f12.6,\002 lineto\002)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt();
    integer s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    static doublereal enrm;
    static integer i__, k;
    static doublereal s, p1[3], p2[3], u1, u2, v1, v2;
    static integer na;
    static doublereal dp[3], du, dv, pm[3], um, vm, err;

    /* Fortran I/O blocks */
    static cilist io___213 = { 1, 0, 0, fmt_100, 0 };
    static cilist io___214 = { 1, 0, 0, fmt_100, 0 };



/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   03/04/03 */

/*   Given unit vectors P and Q corresponding to northern */
/* hemisphere points (with positive third components), this */
/* subroutine draws a polygonal line which approximates the */
/* projection of arc P-Q onto the plane containing the */
/* equator. */

/*   The line segment is drawn by writing a sequence of */
/* 'moveto' and 'lineto' Postscript commands to unit LUN.  It */
/* is assumed that an open file is attached to the unit, */
/* header comments have been written to the file, a window- */
/* to-viewport mapping has been established, etc. */

/* On input: */

/*       LUN = Logical unit number in the range 0 to 99. */

/*       P,Q = Arrays of length 3 containing the endpoints of */
/*             the arc to be drawn. */

/*       TOL = Maximum distance in world coordinates between */
/*             the projected arc and polygonal line. */

/* Input parameters are not altered by this routine. */

/* On output: */

/*       NSEG = Number of line segments in the polygonal */
/*              approximation to the projected arc.  This is */
/*              a decreasing function of TOL.  NSEG = 0 and */
/*              no drawing is performed if P = Q or P = -Q */
/*              or an error is encountered in writing to unit */
/*              LUN. */

/* STRIPACK modules required by DRWARC:  None */

/* Intrinsic functions called by DRWARC:  ABS, DBLE, SQRT */

/* *********************************************************** */


/* Local parameters: */

/* DP =    (Q-P)/NSEG */
/* DU,DV = Components of the projection Q'-P' of arc P->Q */
/*           onto the projection plane */
/* ENRM =  Euclidean norm (or squared norm) of Q'-P' or PM */
/* ERR =   Orthogonal distance from the projected midpoint */
/*           PM' to the line defined by P' and Q': */
/*           |Q'-P' X PM'-P'|/|Q'-P'| */
/* I,K =   DO-loop indexes */
/* NA =    Number of arcs (segments) in the partition of P-Q */
/* P1,P2 = Pairs of adjacent points in a uniform partition of */
/*           arc P-Q into NSEG segments; obtained by normal- */
/*           izing PM values */
/* PM =    Midpoint of arc P-Q or a point P + k*DP in a */
/*           uniform partition of the line segment P-Q into */
/*           NSEG segments */
/* S =     Scale factor 1/NA */
/* U1,V1 = Components of P' */
/* U2,V2 = Components of Q' */
/* UM,VM = Components of the midpoint PM' */


/* Compute the midpoint PM of arc P-Q. */

    /* Parameter adjustments */
    --q;
    --p;

    /* Function Body */
    enrm = 0.;
    for (i__ = 1; i__ <= 3; ++i__) {
	pm[i__ - 1] = p[i__] + q[i__];
	enrm += pm[i__ - 1] * pm[i__ - 1];
/* L1: */
    }
    if (enrm == 0.) {
	goto L5;
    }
    enrm = sqrt(enrm);
    pm[0] /= enrm;
    pm[1] /= enrm;
    pm[2] /= enrm;

/* Project P, Q, and PM to P' = (U1,V1), Q' = (U2,V2), and */
/*   PM' = (UM,VM), respectively. */

    u1 = p[1];
    v1 = p[2];
    u2 = q[1];
    v2 = q[2];
    um = pm[0];
    vm = pm[1];

/* Compute the orthogonal distance ERR from PM' to the line */
/*   defined by P' and Q'.  This is the maximum deviation */
/*   between the projected arc and the line segment.  It is */
/*   undefined if P' = Q'. */

    du = u2 - u1;
    dv = v2 - v1;
    enrm = du * du + dv * dv;
    if (enrm == 0.) {
	goto L5;
    }
    err = (d__1 = du * (vm - v1) - (um - u1) * dv, abs(d__1)) / sqrt(enrm);

/* Compute the number of arcs into which P-Q will be parti- */
/*   tioned (the number of line segments to be drawn): */
/*   NA = ERR/TOL. */

    na = (integer) (err / *tol + 1.);

/* Initialize for loop on arcs P1-P2, where the intermediate */
/*   points are obtained by normalizing PM = P + k*DP for */
/*   DP = (Q-P)/NA */

    s = 1. / (doublereal) na;
    for (i__ = 1; i__ <= 3; ++i__) {
	dp[i__ - 1] = s * (q[i__] - p[i__]);
	pm[i__ - 1] = p[i__];
	p1[i__ - 1] = p[i__];
/* L2: */
    }

/* Loop on arcs P1-P2, drawing the line segments associated */
/*   with the projected endpoints. */

    i__1 = na - 1;
    for (k = 1; k <= i__1; ++k) {
	enrm = 0.;
	for (i__ = 1; i__ <= 3; ++i__) {
	    pm[i__ - 1] += dp[i__ - 1];
	    enrm += pm[i__ - 1] * pm[i__ - 1];
/* L3: */
	}
	if (enrm == 0.) {
	    goto L5;
	}
	enrm = sqrt(enrm);
	p2[0] = pm[0] / enrm;
	p2[1] = pm[1] / enrm;
	p2[2] = pm[2] / enrm;
	io___213.ciunit = *lun;
	i__2 = s_wsfe(&io___213);
	if (i__2 != 0) {
	    goto L5;
	}
	i__2 = do_fio(&c__1, (char *)&p1[0], (ftnlen)sizeof(doublereal));
	if (i__2 != 0) {
	    goto L5;
	}
	i__2 = do_fio(&c__1, (char *)&p1[1], (ftnlen)sizeof(doublereal));
	if (i__2 != 0) {
	    goto L5;
	}
	i__2 = do_fio(&c__1, (char *)&p2[0], (ftnlen)sizeof(doublereal));
	if (i__2 != 0) {
	    goto L5;
	}
	i__2 = do_fio(&c__1, (char *)&p2[1], (ftnlen)sizeof(doublereal));
	if (i__2 != 0) {
	    goto L5;
	}
	i__2 = e_wsfe();
	if (i__2 != 0) {
	    goto L5;
	}
	p1[0] = p2[0];
	p1[1] = p2[1];
	p1[2] = p2[2];
/* L4: */
    }
    io___214.ciunit = *lun;
    i__1 = s_wsfe(&io___214);
    if (i__1 != 0) {
	goto L5;
    }
    i__1 = do_fio(&c__1, (char *)&p1[0], (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L5;
    }
    i__1 = do_fio(&c__1, (char *)&p1[1], (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L5;
    }
    i__1 = do_fio(&c__1, (char *)&q[1], (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L5;
    }
    i__1 = do_fio(&c__1, (char *)&q[2], (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L5;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L5;
    }

/* No error encountered. */

    *nseg = na;
    return 0;

/* Invalid input value of P or Q. */

L5:
    *nseg = 0;
    return 0;
} /* drwarc_ */

/* Subroutine */ int edge_(in1, in2, x, y, z__, lwk, iwk, list, lptr, lend, 
	ier)
integer *in1, *in2;
doublereal *x, *y, *z__;
integer *lwk, *iwk, *list, *lptr, *lend, *ier;
{
    /* Format strings */
    static char fmt_130[] = "(//5x,\002*** Error in EDGE:  Invalid triangul\
a\002,\002tion or null triangles on boundary\002/9x,\002IN1 =\002,i4,\002, I\
N2=\002,i4/)";
    static char fmt_140[] = "(//5x,\002*** Error in OPTIM (called from EDG\
E):\002,\002  NIT = \002,i4,\002, IER = \002,i1,\002 ***\002/)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    extern logical left_();
    static integer ierr;
    extern /* Subroutine */ int swap_();
    static integer next, iwcp1, n1lst, i__, iwend;
    extern /* Subroutine */ int optim_();
    static integer n0, n1, n2;
    static doublereal x0, x1, x2, y0, y1, y2, z0, z1, z2;
    static integer n1frst, nl, lp, nr;
    static doublereal dp12;
    static integer lp21, iwc, iwf, lft, lpl, iwl, nit;
    static doublereal dp1l, dp2l, dp1r, dp2r;

    /* Fortran I/O blocks */
    static cilist io___249 = { 0, 6, 0, fmt_130, 0 };
    static cilist io___250 = { 0, 6, 0, fmt_140, 0 };



/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/30/98 */

/*   Given a triangulation of N nodes and a pair of nodal */
/* indexes IN1 and IN2, this routine swaps arcs as necessary */
/* to force IN1 and IN2 to be adjacent.  Only arcs which */
/* intersect IN1-IN2 are swapped out.  If a Delaunay triangu- */
/* lation is input, the resulting triangulation is as close */
/* as possible to a Delaunay triangulation in the sense that */
/* all arcs other than IN1-IN2 are locally optimal. */

/*   A sequence of calls to EDGE may be used to force the */
/* presence of a set of edges defining the boundary of a non- */
/* convex and/or multiply connected region, or to introduce */
/* barriers into the triangulation.  Note that Subroutine */
/* GETNP will not necessarily return closest nodes if the */
/* triangulation has been constrained by a call to EDGE. */
/* However, this is appropriate in some applications, such */
/* as triangle-based interpolation on a nonconvex domain. */


/* On input: */

/*       IN1,IN2 = Indexes (of X, Y, and Z) in the range 1 to */
/*                 N defining a pair of nodes to be connected */
/*                 by an arc. */

/*       X,Y,Z = Arrays of length N containing the Cartesian */
/*               coordinates of the nodes. */

/* The above parameters are not altered by this routine. */

/*       LWK = Number of columns reserved for IWK.  This must */
/*             be at least NI -- the number of arcs that */
/*             intersect IN1-IN2.  (NI is bounded by N-3.) */

/*       IWK = Integer work array of length at least 2*LWK. */

/*       LIST,LPTR,LEND = Data structure defining the trian- */
/*                        gulation.  Refer to Subroutine */
/*                        TRMESH. */

/* On output: */

/*       LWK = Number of arcs which intersect IN1-IN2 (but */
/*             not more than the input value of LWK) unless */
/*             IER = 1 or IER = 3.  LWK = 0 if and only if */
/*             IN1 and IN2 were adjacent (or LWK=0) on input. */

/*       IWK = Array containing the indexes of the endpoints */
/*             of the new arcs other than IN1-IN2 unless */
/*             IER > 0 or LWK = 0.  New arcs to the left of */
/*             IN1->IN2 are stored in the first K-1 columns */
/*             (left portion of IWK), column K contains */
/*             zeros, and new arcs to the right of IN1->IN2 */
/*             occupy columns K+1,...,LWK.  (K can be deter- */
/*             mined by searching IWK for the zeros.) */

/*       LIST,LPTR,LEND = Data structure updated if necessary */
/*                        to reflect the presence of an arc */
/*                        connecting IN1 and IN2 unless IER > */
/*                        0.  The data structure has been */
/*                        altered if IER >= 4. */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered. */
/*             IER = 1 if IN1 < 1, IN2 < 1, IN1 = IN2, */
/*                     or LWK < 0 on input. */
/*             IER = 2 if more space is required in IWK. */
/*                     Refer to LWK. */
/*             IER = 3 if IN1 and IN2 could not be connected */
/*                     due to either an invalid data struc- */
/*                     ture or collinear nodes (and floating */
/*                     point error). */
/*             IER = 4 if an error flag other than IER = 1 */
/*                     was returned by OPTIM. */
/*             IER = 5 if error flag 1 was returned by OPTIM. */
/*                     This is not necessarily an error, but */
/*                     the arcs other than IN1-IN2 may not */
/*                     be optimal. */

/*   An error message is written to the standard output unit */
/* in the case of IER = 3 or IER = 4. */

/* Modules required by EDGE:  LEFT, LSTPTR, OPTIM, SWAP, */
/*                              SWPTST */

/* Intrinsic function called by EDGE:  ABS */

/* *********************************************************** */


/* Local parameters: */

/* DPij =     Dot product <Ni,Nj> */
/* I =        DO-loop index and column index for IWK */
/* IERR =     Error flag returned by Subroutine OPTIM */
/* IWC =      IWK index between IWF and IWL -- NL->NR is */
/*              stored in IWK(1,IWC)->IWK(2,IWC) */
/* IWCP1 =    IWC + 1 */
/* IWEND =    Input or output value of LWK */
/* IWF =      IWK (column) index of the first (leftmost) arc */
/*              which intersects IN1->IN2 */
/* IWL =      IWK (column) index of the last (rightmost) are */
/*              which intersects IN1->IN2 */
/* LFT =      Flag used to determine if a swap results in the */
/*              new arc intersecting IN1-IN2 -- LFT = 0 iff */
/*              N0 = IN1, LFT = -1 implies N0 LEFT IN1->IN2, */
/*              and LFT = 1 implies N0 LEFT IN2->IN1 */
/* LP =       List pointer (index for LIST and LPTR) */
/* LP21 =     Unused parameter returned by SWAP */
/* LPL =      Pointer to the last neighbor of IN1 or NL */
/* N0 =       Neighbor of N1 or node opposite NR->NL */
/* N1,N2 =    Local copies of IN1 and IN2 */
/* N1FRST =   First neighbor of IN1 */
/* N1LST =    (Signed) last neighbor of IN1 */
/* NEXT =     Node opposite NL->NR */
/* NIT =      Flag or number of iterations employed by OPTIM */
/* NL,NR =    Endpoints of an arc which intersects IN1-IN2 */
/*              with NL LEFT IN1->IN2 */
/* X0,Y0,Z0 = Coordinates of N0 */
/* X1,Y1,Z1 = Coordinates of IN1 */
/* X2,Y2,Z2 = Coordinates of IN2 */


/* Store IN1, IN2, and LWK in local variables and test for */
/*   errors. */

    /* Parameter adjustments */
    --lend;
    --lptr;
    --list;
    iwk -= 3;
    --z__;
    --y;
    --x;

    /* Function Body */
    n1 = *in1;
    n2 = *in2;
    iwend = *lwk;
    if (n1 < 1 || n2 < 1 || n1 == n2 || iwend < 0) {
	goto L31;
    }

/* Test for N2 as a neighbor of N1.  LPL points to the last */
/*   neighbor of N1. */

    lpl = lend[n1];
    n0 = (i__1 = list[lpl], abs(i__1));
    lp = lpl;
L1:
    if (n0 == n2) {
	goto L30;
    }
    lp = lptr[lp];
    n0 = list[lp];
    if (lp != lpl) {
	goto L1;
    }

/* Initialize parameters. */

    iwl = 0;
    nit = 0;

/* Store the coordinates of N1 and N2. */

L2:
    x1 = x[n1];
    y1 = y[n1];
    z1 = z__[n1];
    x2 = x[n2];
    y2 = y[n2];
    z2 = z__[n2];

/* Set NR and NL to adjacent neighbors of N1 such that */
/*   NR LEFT N2->N1 and NL LEFT N1->N2, */
/*   (NR Forward N1->N2 or NL Forward N1->N2), and */
/*   (NR Forward N2->N1 or NL Forward N2->N1). */

/*   Initialization:  Set N1FRST and N1LST to the first and */
/*     (signed) last neighbors of N1, respectively, and */
/*     initialize NL to N1FRST. */

    lpl = lend[n1];
    n1lst = list[lpl];
    lp = lptr[lpl];
    n1frst = list[lp];
    nl = n1frst;
    if (n1lst < 0) {
	goto L4;
    }

/*   N1 is an interior node.  Set NL to the first candidate */
/*     for NR (NL LEFT N2->N1). */

L3:
    if (left_(&x2, &y2, &z2, &x1, &y1, &z1, &x[nl], &y[nl], &z__[nl])) {
	goto L4;
    }
    lp = lptr[lp];
    nl = list[lp];
    if (nl != n1frst) {
	goto L3;
    }

/*   All neighbors of N1 are strictly left of N1->N2. */

    goto L5;

/*   NL = LIST(LP) LEFT N2->N1.  Set NR to NL and NL to the */
/*     following neighbor of N1. */

L4:
    nr = nl;
    lp = lptr[lp];
    nl = (i__1 = list[lp], abs(i__1));
    if (left_(&x1, &y1, &z1, &x2, &y2, &z2, &x[nl], &y[nl], &z__[nl])) {

/*   NL LEFT N1->N2 and NR LEFT N2->N1.  The Forward tests */
/*     are employed to avoid an error associated with */
/*     collinear nodes. */

	dp12 = x1 * x2 + y1 * y2 + z1 * z2;
	dp1l = x1 * x[nl] + y1 * y[nl] + z1 * z__[nl];
	dp2l = x2 * x[nl] + y2 * y[nl] + z2 * z__[nl];
	dp1r = x1 * x[nr] + y1 * y[nr] + z1 * z__[nr];
	dp2r = x2 * x[nr] + y2 * y[nr] + z2 * z__[nr];
	if ((dp2l - dp12 * dp1l >= 0. || dp2r - dp12 * dp1r >= 0.) && (dp1l - 
		dp12 * dp2l >= 0. || dp1r - dp12 * dp2r >= 0.)) {
	    goto L6;
	}

/*   NL-NR does not intersect N1-N2.  However, there is */
/*     another candidate for the first arc if NL lies on */
/*     the line N1-N2. */

	if (! left_(&x2, &y2, &z2, &x1, &y1, &z1, &x[nl], &y[nl], &z__[nl])) {
	    goto L5;
	}
    }

/*   Bottom of loop. */

    if (nl != n1frst) {
	goto L4;
    }

/* Either the triangulation is invalid or N1-N2 lies on the */
/*   convex hull boundary and an edge NR->NL (opposite N1 and */
/*   intersecting N1-N2) was not found due to floating point */
/*   error.  Try interchanging N1 and N2 -- NIT > 0 iff this */
/*   has already been done. */

L5:
    if (nit > 0) {
	goto L33;
    }
    nit = 1;
    n1 = n2;
    n2 = *in1;
    goto L2;

/* Store the ordered sequence of intersecting edges NL->NR in */
/*   IWK(1,IWL)->IWK(2,IWL). */

L6:
    ++iwl;
    if (iwl > iwend) {
	goto L32;
    }
    iwk[(iwl << 1) + 1] = nl;
    iwk[(iwl << 1) + 2] = nr;

/*   Set NEXT to the neighbor of NL which follows NR. */

    lpl = lend[nl];
    lp = lptr[lpl];

/*   Find NR as a neighbor of NL.  The search begins with */
/*     the first neighbor. */

L7:
    if (list[lp] == nr) {
	goto L8;
    }
    lp = lptr[lp];
    if (lp != lpl) {
	goto L7;
    }

/*   NR must be the last neighbor, and NL->NR cannot be a */
/*     boundary edge. */

    if (list[lp] != nr) {
	goto L33;
    }

/*   Set NEXT to the neighbor following NR, and test for */
/*     termination of the store loop. */

L8:
    lp = lptr[lp];
    next = (i__1 = list[lp], abs(i__1));
    if (next == n2) {
	goto L9;
    }

/*   Set NL or NR to NEXT. */

    if (left_(&x1, &y1, &z1, &x2, &y2, &z2, &x[next], &y[next], &z__[next])) {
	nl = next;
    } else {
	nr = next;
    }
    goto L6;

/* IWL is the number of arcs which intersect N1-N2. */
/*   Store LWK. */

L9:
    *lwk = iwl;
    iwend = iwl;

/* Initialize for edge swapping loop -- all possible swaps */
/*   are applied (even if the new arc again intersects */
/*   N1-N2), arcs to the left of N1->N2 are stored in the */
/*   left portion of IWK, and arcs to the right are stored in */
/*   the right portion.  IWF and IWL index the first and last */
/*   intersecting arcs. */

    iwf = 1;

/* Top of loop -- set N0 to N1 and NL->NR to the first edge. */
/*   IWC points to the arc currently being processed.  LFT */
/*   .LE. 0 iff N0 LEFT N1->N2. */

L10:
    lft = 0;
    n0 = n1;
    x0 = x1;
    y0 = y1;
    z0 = z1;
    nl = iwk[(iwf << 1) + 1];
    nr = iwk[(iwf << 1) + 2];
    iwc = iwf;

/*   Set NEXT to the node opposite NL->NR unless IWC is the */
/*     last arc. */

L11:
    if (iwc == iwl) {
	goto L21;
    }
    iwcp1 = iwc + 1;
    next = iwk[(iwcp1 << 1) + 1];
    if (next != nl) {
	goto L16;
    }
    next = iwk[(iwcp1 << 1) + 2];

/*   NEXT RIGHT N1->N2 and IWC .LT. IWL.  Test for a possible */
/*     swap. */

    if (! left_(&x0, &y0, &z0, &x[nr], &y[nr], &z__[nr], &x[next], &y[next], &
	    z__[next])) {
	goto L14;
    }
    if (lft >= 0) {
	goto L12;
    }
    if (! left_(&x[nl], &y[nl], &z__[nl], &x0, &y0, &z0, &x[next], &y[next], &
	    z__[next])) {
	goto L14;
    }

/*   Replace NL->NR with N0->NEXT. */

    swap_(&next, &n0, &nl, &nr, &list[1], &lptr[1], &lend[1], &lp21);
    iwk[(iwc << 1) + 1] = n0;
    iwk[(iwc << 1) + 2] = next;
    goto L15;

/*   Swap NL-NR for N0-NEXT, shift columns IWC+1,...,IWL to */
/*     the left, and store N0-NEXT in the right portion of */
/*     IWK. */

L12:
    swap_(&next, &n0, &nl, &nr, &list[1], &lptr[1], &lend[1], &lp21);
    i__1 = iwl;
    for (i__ = iwcp1; i__ <= i__1; ++i__) {
	iwk[(i__ - 1 << 1) + 1] = iwk[(i__ << 1) + 1];
	iwk[(i__ - 1 << 1) + 2] = iwk[(i__ << 1) + 2];
/* L13: */
    }
    iwk[(iwl << 1) + 1] = n0;
    iwk[(iwl << 1) + 2] = next;
    --iwl;
    nr = next;
    goto L11;

/*   A swap is not possible.  Set N0 to NR. */

L14:
    n0 = nr;
    x0 = x[n0];
    y0 = y[n0];
    z0 = z__[n0];
    lft = 1;

/*   Advance to the next arc. */

L15:
    nr = next;
    ++iwc;
    goto L11;

/*   NEXT LEFT N1->N2, NEXT .NE. N2, and IWC .LT. IWL. */
/*     Test for a possible swap. */

L16:
    if (! left_(&x[nl], &y[nl], &z__[nl], &x0, &y0, &z0, &x[next], &y[next], &
	    z__[next])) {
	goto L19;
    }
    if (lft <= 0) {
	goto L17;
    }
    if (! left_(&x0, &y0, &z0, &x[nr], &y[nr], &z__[nr], &x[next], &y[next], &
	    z__[next])) {
	goto L19;
    }

/*   Replace NL->NR with NEXT->N0. */

    swap_(&next, &n0, &nl, &nr, &list[1], &lptr[1], &lend[1], &lp21);
    iwk[(iwc << 1) + 1] = next;
    iwk[(iwc << 1) + 2] = n0;
    goto L20;

/*   Swap NL-NR for N0-NEXT, shift columns IWF,...,IWC-1 to */
/*     the right, and store N0-NEXT in the left portion of */
/*     IWK. */

L17:
    swap_(&next, &n0, &nl, &nr, &list[1], &lptr[1], &lend[1], &lp21);
    i__1 = iwf;
    for (i__ = iwc - 1; i__ >= i__1; --i__) {
	iwk[(i__ + 1 << 1) + 1] = iwk[(i__ << 1) + 1];
	iwk[(i__ + 1 << 1) + 2] = iwk[(i__ << 1) + 2];
/* L18: */
    }
    iwk[(iwf << 1) + 1] = n0;
    iwk[(iwf << 1) + 2] = next;
    ++iwf;
    goto L20;

/*   A swap is not possible.  Set N0 to NL. */

L19:
    n0 = nl;
    x0 = x[n0];
    y0 = y[n0];
    z0 = z__[n0];
    lft = -1;

/*   Advance to the next arc. */

L20:
    nl = next;
    ++iwc;
    goto L11;

/*   N2 is opposite NL->NR (IWC = IWL). */

L21:
    if (n0 == n1) {
	goto L24;
    }
    if (lft < 0) {
	goto L22;
    }

/*   N0 RIGHT N1->N2.  Test for a possible swap. */

    if (! left_(&x0, &y0, &z0, &x[nr], &y[nr], &z__[nr], &x2, &y2, &z2)) {
	goto L10;
    }

/*   Swap NL-NR for N0-N2 and store N0-N2 in the right */
/*     portion of IWK. */

    swap_(&n2, &n0, &nl, &nr, &list[1], &lptr[1], &lend[1], &lp21);
    iwk[(iwl << 1) + 1] = n0;
    iwk[(iwl << 1) + 2] = n2;
    --iwl;
    goto L10;

/*   N0 LEFT N1->N2.  Test for a possible swap. */

L22:
    if (! left_(&x[nl], &y[nl], &z__[nl], &x0, &y0, &z0, &x2, &y2, &z2)) {
	goto L10;
    }

/*   Swap NL-NR for N0-N2, shift columns IWF,...,IWL-1 to the */
/*     right, and store N0-N2 in the left portion of IWK. */

    swap_(&n2, &n0, &nl, &nr, &list[1], &lptr[1], &lend[1], &lp21);
    i__ = iwl;
L23:
    iwk[(i__ << 1) + 1] = iwk[(i__ - 1 << 1) + 1];
    iwk[(i__ << 1) + 2] = iwk[(i__ - 1 << 1) + 2];
    --i__;
    if (i__ > iwf) {
	goto L23;
    }
    iwk[(iwf << 1) + 1] = n0;
    iwk[(iwf << 1) + 2] = n2;
    ++iwf;
    goto L10;

/* IWF = IWC = IWL.  Swap out the last arc for N1-N2 and */
/*   store zeros in IWK. */

L24:
    swap_(&n2, &n1, &nl, &nr, &list[1], &lptr[1], &lend[1], &lp21);
    iwk[(iwc << 1) + 1] = 0;
    iwk[(iwc << 1) + 2] = 0;

/* Optimization procedure -- */

    *ier = 0;
    if (iwc > 1) {

/*   Optimize the set of new arcs to the left of IN1->IN2. */

	nit = iwc - 1 << 2;
	i__1 = iwc - 1;
	optim_(&x[1], &y[1], &z__[1], &i__1, &list[1], &lptr[1], &lend[1], &
		nit, &iwk[3], &ierr);
	if (ierr != 0 && ierr != 1) {
	    goto L34;
	}
	if (ierr == 1) {
	    *ier = 5;
	}
    }
    if (iwc < iwend) {

/*   Optimize the set of new arcs to the right of IN1->IN2. */

	nit = iwend - iwc << 2;
	i__1 = iwend - iwc;
	optim_(&x[1], &y[1], &z__[1], &i__1, &list[1], &lptr[1], &lend[1], &
		nit, &iwk[(iwc + 1 << 1) + 1], &ierr);
	if (ierr != 0 && ierr != 1) {
	    goto L34;
	}
	if (ierr == 1) {
	    goto L35;
	}
    }
    if (*ier == 5) {
	goto L35;
    }

/* Successful termination (IER = 0). */

    return 0;

/* IN1 and IN2 were adjacent on input. */

L30:
    *ier = 0;
    return 0;

/* Invalid input parameter. */

L31:
    *ier = 1;
    return 0;

/* Insufficient space reserved for IWK. */

L32:
    *ier = 2;
    return 0;

/* Invalid triangulation data structure or collinear nodes */
/*   on convex hull boundary. */

L33:
    *ier = 3;
    s_wsfe(&io___249);
    do_fio(&c__1, (char *)&(*in1), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*in2), (ftnlen)sizeof(integer));
    e_wsfe();
    return 0;

/* Error flag (other than 1) returned by OPTIM. */

L34:
    *ier = 4;
    s_wsfe(&io___250);
    do_fio(&c__1, (char *)&nit, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ierr, (ftnlen)sizeof(integer));
    e_wsfe();
    return 0;

/* Error flag 1 returned by OPTIM. */

L35:
    *ier = 5;
    return 0;
} /* edge_ */

/* Subroutine */ int getnp_(x, y, z__, list, lptr, lend, l, npts, df, ier)
doublereal *x, *y, *z__;
integer *list, *lptr, *lend, *l, *npts;
doublereal *df;
integer *ier;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, n1;
    static doublereal x1, y1, z1;
    static integer nb, ni, lp, np, lm1;
    static doublereal dnb, dnp;
    static integer lpl;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/28/98 */

/*   Given a Delaunay triangulation of N nodes on the unit */
/* sphere and an array NPTS containing the indexes of L-1 */
/* nodes ordered by angular distance from NPTS(1), this sub- */
/* routine sets NPTS(L) to the index of the next node in the */
/* sequence -- the node, other than NPTS(1),...,NPTS(L-1), */
/* that is closest to NPTS(1).  Thus, the ordered sequence */
/* of K closest nodes to N1 (including N1) may be determined */
/* by K-1 calls to GETNP with NPTS(1) = N1 and L = 2,3,...,K */
/* for K .GE. 2. */

/*   The algorithm uses the property of a Delaunay triangula- */
/* tion that the K-th closest node to N1 is a neighbor of one */
/* of the K-1 closest nodes to N1. */


/* On input: */

/*       X,Y,Z = Arrays of length N containing the Cartesian */
/*               coordinates of the nodes. */

/*       LIST,LPTR,LEND = Triangulation data structure.  Re- */
/*                        fer to Subroutine TRMESH. */

/*       L = Number of nodes in the sequence on output.  2 */
/*           .LE. L .LE. N. */

/* The above parameters are not altered by this routine. */

/*       NPTS = Array of length .GE. L containing the indexes */
/*              of the L-1 closest nodes to NPTS(1) in the */
/*              first L-1 locations. */

/* On output: */

/*       NPTS = Array updated with the index of the L-th */
/*              closest node to NPTS(1) in position L unless */
/*              IER = 1. */

/*       DF = Value of an increasing function (negative cos- */
/*            ine) of the angular distance between NPTS(1) */
/*            and NPTS(L) unless IER = 1. */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered. */
/*             IER = 1 if L < 2. */

/* Modules required by GETNP:  None */

/* Intrinsic function called by GETNP:  ABS */

/* *********************************************************** */


/* Local parameters: */

/* DNB,DNP =  Negative cosines of the angular distances from */
/*              N1 to NB and to NP, respectively */
/* I =        NPTS index and DO-loop index */
/* LM1 =      L-1 */
/* LP =       LIST pointer of a neighbor of NI */
/* LPL =      Pointer to the last neighbor of NI */
/* N1 =       NPTS(1) */
/* NB =       Neighbor of NI and candidate for NP */
/* NI =       NPTS(I) */
/* NP =       Candidate for NPTS(L) */
/* X1,Y1,Z1 = Coordinates of N1 */

    /* Parameter adjustments */
    --x;
    --y;
    --z__;
    --list;
    --lptr;
    --lend;
    --npts;

    /* Function Body */
    lm1 = *l - 1;
    if (lm1 < 1) {
	goto L6;
    }
    *ier = 0;

/* Store N1 = NPTS(1) and mark the elements of NPTS. */

    n1 = npts[1];
    x1 = x[n1];
    y1 = y[n1];
    z1 = z__[n1];
    i__1 = lm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ni = npts[i__];
	lend[ni] = -lend[ni];
/* L1: */
    }

/* Candidates for NP = NPTS(L) are the unmarked neighbors */
/*   of nodes in NPTS.  DNP is initially greater than -cos(PI) */
/*   (the maximum distance). */

    dnp = 2.;

/* Loop on nodes NI in NPTS. */

    i__1 = lm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ni = npts[i__];
	lpl = -lend[ni];
	lp = lpl;

/* Loop on neighbors NB of NI. */

L2:
	nb = (i__2 = list[lp], abs(i__2));
	if (lend[nb] < 0) {
	    goto L3;
	}

/* NB is an unmarked neighbor of NI.  Replace NP if NB is */
/*   closer to N1. */

	dnb = -(x[nb] * x1 + y[nb] * y1 + z__[nb] * z1);
	if (dnb >= dnp) {
	    goto L3;
	}
	np = nb;
	dnp = dnb;
L3:
	lp = lptr[lp];
	if (lp != lpl) {
	    goto L2;
	}
/* L4: */
    }
    npts[*l] = np;
    *df = dnp;

/* Unmark the elements of NPTS. */

    i__1 = lm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ni = npts[i__];
	lend[ni] = -lend[ni];
/* L5: */
    }
    return 0;

/* L is outside its valid range. */

L6:
    *ier = 1;
    return 0;
} /* getnp_ */

/* Subroutine */ int insert_(k, lp, list, lptr, lnew)
integer *k, *lp, *list, *lptr, *lnew;
{
    static integer lsav;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/17/96 */

/*   This subroutine inserts K as a neighbor of N1 following */
/* N2, where LP is the LIST pointer of N2 as a neighbor of */
/* N1.  Note that, if N2 is the last neighbor of N1, K will */
/* become the first neighbor (even if N1 is a boundary node). */

/*   This routine is identical to the similarly named routine */
/* in TRIPACK. */


/* On input: */

/*       K = Index of the node to be inserted. */

/*       LP = LIST pointer of N2 as a neighbor of N1. */

/* The above parameters are not altered by this routine. */

/*       LIST,LPTR,LNEW = Data structure defining the trian- */
/*                        gulation.  Refer to Subroutine */
/*                        TRMESH. */

/* On output: */

/*       LIST,LPTR,LNEW = Data structure updated with the */
/*                        addition of node K. */

/* Modules required by INSERT:  None */

/* *********************************************************** */


    /* Parameter adjustments */
    --lptr;
    --list;

    /* Function Body */
    lsav = lptr[*lp];
    lptr[*lp] = *lnew;
    list[*lnew] = *k;
    lptr[*lnew] = lsav;
    ++(*lnew);
    return 0;
} /* insert_ */

logical inside_(p, lv, xv, yv, zv, nv, listv, ier)
doublereal *p;
integer *lv;
doublereal *xv, *yv, *zv;
integer *nv, *listv, *ier;
{
    /* Initialized data */

    static doublereal eps = .001;

    /* System generated locals */
    integer i__1;
    logical ret_val;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static logical even;
    static integer ierr;
    static logical pinr, qinr;
    static doublereal qnrm, vnrm, b[3], d__;
    static integer k, n;
    static doublereal q[3];
    static integer i1, i2, k0;
    static doublereal v1[3], v2[3], cn[3], bp, bq;
    static integer ni;
    static doublereal pn[3], qn[3], vn[3];
    extern /* Subroutine */ int intrsc_();
    static integer imx;
    static logical lft1, lft2;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   12/27/93 */

/*   This function locates a point P relative to a polygonal */
/* region R on the surface of the unit sphere, returning */
/* INSIDE = TRUE if and only if P is contained in R.  R is */
/* defined by a cyclically ordered sequence of vertices which */
/* form a positively-oriented simple closed curve.  Adjacent */
/* vertices need not be distinct but the curve must not be */
/* self-intersecting.  Also, while polygon edges are by defi- */
/* nition restricted to a single hemisphere, R is not so */
/* restricted.  Its interior is the region to the left as the */
/* vertices are traversed in order. */

/*   The algorithm consists of selecting a point Q in R and */
/* then finding all points at which the great circle defined */
/* by P and Q intersects the boundary of R.  P lies inside R */
/* if and only if there is an even number of intersection */
/* points between Q and P.  Q is taken to be a point immedi- */
/* ately to the left of a directed boundary edge -- the first */
/* one that results in no consistency-check failures. */

/*   If P is close to the polygon boundary, the problem is */
/* ill-conditioned and the decision may be incorrect.  Also, */
/* an incorrect decision may result from a poor choice of Q */
/* (if, for example, a boundary edge lies on the great cir- */
/* cle defined by P and Q).  A more reliable result could be */
/* obtained by a sequence of calls to INSIDE with the ver- */
/* tices cyclically permuted before each call (to alter the */
/* choice of Q). */


/* On input: */

/*       P = Array of length 3 containing the Cartesian */
/*           coordinates of the point (unit vector) to be */
/*           located. */

/*       LV = Length of arrays XV, YV, and ZV. */

/*       XV,YV,ZV = Arrays of length LV containing the Carte- */
/*                  sian coordinates of unit vectors (points */
/*                  on the unit sphere).  These values are */
/*                  not tested for validity. */

/*       NV = Number of vertices in the polygon.  3 .LE. NV */
/*            .LE. LV. */

/*       LISTV = Array of length NV containing the indexes */
/*               (for XV, YV, and ZV) of a cyclically-ordered */
/*               (and CCW-ordered) sequence of vertices that */
/*               define R.  The last vertex (indexed by */
/*               LISTV(NV)) is followed by the first (indexed */
/*               by LISTV(1)).  LISTV entries must be in the */
/*               range 1 to LV. */

/* Input parameters are not altered by this function. */

/* On output: */

/*       INSIDE = TRUE if and only if P lies inside R unless */
/*                IER .NE. 0, in which case the value is not */
/*                altered. */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered. */
/*             IER = 1 if LV or NV is outside its valid */
/*                     range. */
/*             IER = 2 if a LISTV entry is outside its valid */
/*                     range. */
/*             IER = 3 if the polygon boundary was found to */
/*                     be self-intersecting.  This error will */
/*                     not necessarily be detected. */
/*             IER = 4 if every choice of Q (one for each */
/*                     boundary edge) led to failure of some */
/*                     internal consistency check.  The most */
/*                     likely cause of this error is invalid */
/*                     input:  P = (0,0,0), a null or self- */
/*                     intersecting polygon, etc. */

/* Module required by INSIDE:  INTRSC */

/* Intrinsic function called by INSIDE:  SQRT */

/* *********************************************************** */


/* Local parameters: */

/* B =         Intersection point between the boundary and */
/*               the great circle defined by P and Q */
/* BP,BQ =     <B,P> and <B,Q>, respectively, maximized over */
/*               intersection points B that lie between P and */
/*               Q (on the shorter arc) -- used to find the */
/*               closest intersection points to P and Q */
/* CN =        Q X P = normal to the plane of P and Q */
/* D =         Dot product <B,P> or <B,Q> */
/* EPS =       Parameter used to define Q as the point whose */
/*               orthogonal distance to (the midpoint of) */
/*               boundary edge V1->V2 is approximately EPS/ */
/*               (2*Cos(A/2)), where <V1,V2> = Cos(A). */
/* EVEN =      TRUE iff an even number of intersection points */
/*               lie between P and Q (on the shorter arc) */
/* I1,I2 =     Indexes (LISTV elements) of a pair of adjacent */
/*               boundary vertices (endpoints of a boundary */
/*               edge) */
/* IERR =      Error flag for calls to INTRSC (not tested) */
/* IMX =       Local copy of LV and maximum value of I1 and */
/*               I2 */
/* K =         DO-loop index and LISTV index */
/* K0 =        LISTV index of the first endpoint of the */
/*               boundary edge used to compute Q */
/* LFT1,LFT2 = Logical variables associated with I1 and I2 in */
/*               the boundary traversal:  TRUE iff the vertex */
/*               is strictly to the left of Q->P (<V,CN> > 0) */
/* N =         Local copy of NV */
/* NI =        Number of intersections (between the boundary */
/*               curve and the great circle P-Q) encountered */
/* PINR =      TRUE iff P is to the left of the directed */
/*               boundary edge associated with the closest */
/*               intersection point to P that lies between P */
/*               and Q (a left-to-right intersection as */
/*               viewed from Q), or there is no intersection */
/*               between P and Q (on the shorter arc) */
/* PN,QN =     P X CN and CN X Q, respectively:  used to */
/*               locate intersections B relative to arc Q->P */
/* Q =         (V1 + V2 + EPS*VN/VNRM)/QNRM, where V1->V2 is */
/*               the boundary edge indexed by LISTV(K0) -> */
/*               LISTV(K0+1) */
/* QINR =      TRUE iff Q is to the left of the directed */
/*               boundary edge associated with the closest */
/*               intersection point to Q that lies between P */
/*               and Q (a right-to-left intersection as */
/*               viewed from Q), or there is no intersection */
/*               between P and Q (on the shorter arc) */
/* QNRM =      Euclidean norm of V1+V2+EPS*VN/VNRM used to */
/*               compute (normalize) Q */
/* V1,V2 =     Vertices indexed by I1 and I2 in the boundary */
/*               traversal */
/* VN =        V1 X V2, where V1->V2 is the boundary edge */
/*               indexed by LISTV(K0) -> LISTV(K0+1) */
/* VNRM =      Euclidean norm of VN */

    /* Parameter adjustments */
    --p;
    --zv;
    --yv;
    --xv;
    --listv;

    /* Function Body */

/* Store local parameters, test for error 1, and initialize */
/*   K0. */

    imx = *lv;
    n = *nv;
    if (n < 3 || n > imx) {
	goto L11;
    }
    k0 = 0;
    i1 = listv[1];
    if (i1 < 1 || i1 > imx) {
	goto L12;
    }

/* Increment K0 and set Q to a point immediately to the left */
/*   of the midpoint of edge V1->V2 = LISTV(K0)->LISTV(K0+1): */
/*   Q = (V1 + V2 + EPS*VN/VNRM)/QNRM, where VN = V1 X V2. */

L1:
    ++k0;
    if (k0 > n) {
	goto L14;
    }
    i1 = listv[k0];
    if (k0 < n) {
	i2 = listv[k0 + 1];
    } else {
	i2 = listv[1];
    }
    if (i2 < 1 || i2 > imx) {
	goto L12;
    }
    vn[0] = yv[i1] * zv[i2] - zv[i1] * yv[i2];
    vn[1] = zv[i1] * xv[i2] - xv[i1] * zv[i2];
    vn[2] = xv[i1] * yv[i2] - yv[i1] * xv[i2];
    vnrm = sqrt(vn[0] * vn[0] + vn[1] * vn[1] + vn[2] * vn[2]);
    if (vnrm == 0.) {
	goto L1;
    }
    q[0] = xv[i1] + xv[i2] + eps * vn[0] / vnrm;
    q[1] = yv[i1] + yv[i2] + eps * vn[1] / vnrm;
    q[2] = zv[i1] + zv[i2] + eps * vn[2] / vnrm;
    qnrm = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
    q[0] /= qnrm;
    q[1] /= qnrm;
    q[2] /= qnrm;

/* Compute CN = Q X P, PN = P X CN, and QN = CN X Q. */

    cn[0] = q[1] * p[3] - q[2] * p[2];
    cn[1] = q[2] * p[1] - q[0] * p[3];
    cn[2] = q[0] * p[2] - q[1] * p[1];
    if (cn[0] == 0. && cn[1] == 0. && cn[2] == 0.) {
	goto L1;
    }
    pn[0] = p[2] * cn[2] - p[3] * cn[1];
    pn[1] = p[3] * cn[0] - p[1] * cn[2];
    pn[2] = p[1] * cn[1] - p[2] * cn[0];
    qn[0] = cn[1] * q[2] - cn[2] * q[1];
    qn[1] = cn[2] * q[0] - cn[0] * q[2];
    qn[2] = cn[0] * q[1] - cn[1] * q[0];

/* Initialize parameters for the boundary traversal. */

    ni = 0;
    even = TRUE_;
    bp = -2.;
    bq = -2.;
    pinr = TRUE_;
    qinr = TRUE_;
    i2 = listv[n];
    if (i2 < 1 || i2 > imx) {
	goto L12;
    }
    lft2 = cn[0] * xv[i2] + cn[1] * yv[i2] + cn[2] * zv[i2] > 0.;

/* Loop on boundary arcs I1->I2. */

    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	i1 = i2;
	lft1 = lft2;
	i2 = listv[k];
	if (i2 < 1 || i2 > imx) {
	    goto L12;
	}
	lft2 = cn[0] * xv[i2] + cn[1] * yv[i2] + cn[2] * zv[i2] > 0.;
	if (lft1 == lft2) {
	    goto L2;
	}

/*   I1 and I2 are on opposite sides of Q->P.  Compute the */
/*     point of intersection B. */

	++ni;
	v1[0] = xv[i1];
	v1[1] = yv[i1];
	v1[2] = zv[i1];
	v2[0] = xv[i2];
	v2[1] = yv[i2];
	v2[2] = zv[i2];
	intrsc_(v1, v2, cn, b, &ierr);

/*   B is between Q and P (on the shorter arc) iff */
/*     B Forward Q->P and B Forward P->Q       iff */
/*     <B,QN> > 0 and <B,PN> > 0. */

	if (b[0] * qn[0] + b[1] * qn[1] + b[2] * qn[2] > 0. && b[0] * pn[0] + 
		b[1] * pn[1] + b[2] * pn[2] > 0.) {

/*   Update EVEN, BQ, QINR, BP, and PINR. */

	    even = ! even;
	    d__ = b[0] * q[0] + b[1] * q[1] + b[2] * q[2];
	    if (d__ > bq) {
		bq = d__;
		qinr = lft2;
	    }
	    d__ = b[0] * p[1] + b[1] * p[2] + b[2] * p[3];
	    if (d__ > bp) {
		bp = d__;
		pinr = lft1;
	    }
	}
L2:
	;
    }

/* Test for consistency:  NI must be even and QINR must be */
/*   TRUE. */

    if (ni != ni / 2 << 1 || ! qinr) {
	goto L1;
    }

/* Test for error 3:  different values of PINR and EVEN. */

    if (pinr != even) {
	goto L13;
    }

/* No error encountered. */

    *ier = 0;
    ret_val = even;
    return ret_val;

/* LV or NV is outside its valid range. */

L11:
    *ier = 1;
    return ret_val;

/* A LISTV entry is outside its valid range. */

L12:
    *ier = 2;
    return ret_val;

/* The polygon boundary is self-intersecting. */

L13:
    *ier = 3;
    return ret_val;

/* Consistency tests failed for all values of Q. */

L14:
    *ier = 4;
    return ret_val;
} /* inside_ */

/* Subroutine */ int intadd_(kk, i1, i2, i3, list, lptr, lend, lnew)
integer *kk, *i1, *i2, *i3, *list, *lptr, *lend, *lnew;
{
    static integer k, n1, n2, n3, lp;
    extern /* Subroutine */ int insert_();
    extern integer lstptr_();


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/17/96 */

/*   This subroutine adds an interior node to a triangulation */
/* of a set of points on the unit sphere.  The data structure */
/* is updated with the insertion of node KK into the triangle */
/* whose vertices are I1, I2, and I3.  No optimization of the */
/* triangulation is performed. */

/*   This routine is identical to the similarly named routine */
/* in TRIPACK. */


/* On input: */

/*       KK = Index of the node to be inserted.  KK .GE. 1 */
/*            and KK must not be equal to I1, I2, or I3. */

/*       I1,I2,I3 = Indexes of the counterclockwise-ordered */
/*                  sequence of vertices of a triangle which */
/*                  contains node KK. */

/* The above parameters are not altered by this routine. */

/*       LIST,LPTR,LEND,LNEW = Data structure defining the */
/*                             triangulation.  Refer to Sub- */
/*                             routine TRMESH.  Triangle */
/*                             (I1,I2,I3) must be included */
/*                             in the triangulation. */

/* On output: */

/*       LIST,LPTR,LEND,LNEW = Data structure updated with */
/*                             the addition of node KK.  KK */
/*                             will be connected to nodes I1, */
/*                             I2, and I3. */

/* Modules required by INTADD:  INSERT, LSTPTR */

/* *********************************************************** */


/* Local parameters: */

/* K =        Local copy of KK */
/* LP =       LIST pointer */
/* N1,N2,N3 = Local copies of I1, I2, and I3 */

    /* Parameter adjustments */
    --lend;
    --lptr;
    --list;

    /* Function Body */
    k = *kk;

/* Initialization. */

    n1 = *i1;
    n2 = *i2;
    n3 = *i3;

/* Add K as a neighbor of I1, I2, and I3. */

    lp = lstptr_(&lend[n1], &n2, &list[1], &lptr[1]);
    insert_(&k, &lp, &list[1], &lptr[1], lnew);
    lp = lstptr_(&lend[n2], &n3, &list[1], &lptr[1]);
    insert_(&k, &lp, &list[1], &lptr[1], lnew);
    lp = lstptr_(&lend[n3], &n1, &list[1], &lptr[1]);
    insert_(&k, &lp, &list[1], &lptr[1], lnew);

/* Add I1, I2, and I3 as neighbors of K. */

    list[*lnew] = n1;
    list[*lnew + 1] = n2;
    list[*lnew + 2] = n3;
    lptr[*lnew] = *lnew + 1;
    lptr[*lnew + 1] = *lnew + 2;
    lptr[*lnew + 2] = *lnew;
    lend[k] = *lnew + 2;
    *lnew += 3;
    return 0;
} /* intadd_ */

/* Subroutine */ int intrsc_(p1, p2, cn, p, ier)
doublereal *p1, *p2, *cn, *p;
integer *ier;
{
    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer i__;
    static doublereal t, d1, d2, pp[3], ppn;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/19/90 */

/*   Given a great circle C and points P1 and P2 defining an */
/* arc A on the surface of the unit sphere, where A is the */
/* shorter of the two portions of the great circle C12 assoc- */
/* iated with P1 and P2, this subroutine returns the point */
/* of intersection P between C and C12 that is closer to A. */
/* Thus, if P1 and P2 lie in opposite hemispheres defined by */
/* C, P is the point of intersection of C with A. */


/* On input: */

/*       P1,P2 = Arrays of length 3 containing the Cartesian */
/*               coordinates of unit vectors. */

/*       CN = Array of length 3 containing the Cartesian */
/*            coordinates of a nonzero vector which defines C */
/*            as the intersection of the plane whose normal */
/*            is CN with the unit sphere.  Thus, if C is to */
/*            be the great circle defined by P and Q, CN */
/*            should be P X Q. */

/* The above parameters are not altered by this routine. */

/*       P = Array of length 3. */

/* On output: */

/*       P = Point of intersection defined above unless IER */
/*           .NE. 0, in which case P is not altered. */

/*       IER = Error indicator. */
/*             IER = 0 if no errors were encountered. */
/*             IER = 1 if <CN,P1> = <CN,P2>.  This occurs */
/*                     iff P1 = P2 or CN = 0 or there are */
/*                     two intersection points at the same */
/*                     distance from A. */
/*             IER = 2 if P2 = -P1 and the definition of A is */
/*                     therefore ambiguous. */

/* Modules required by INTRSC:  None */

/* Intrinsic function called by INTRSC:  SQRT */

/* *********************************************************** */


/* Local parameters: */

/* D1 =  <CN,P1> */
/* D2 =  <CN,P2> */
/* I =   DO-loop index */
/* PP =  P1 + T*(P2-P1) = Parametric representation of the */
/*         line defined by P1 and P2 */
/* PPN = Norm of PP */
/* T =   D1/(D1-D2) = Parameter value chosen so that PP lies */
/*         in the plane of C */

    /* Parameter adjustments */
    --p;
    --cn;
    --p2;
    --p1;

    /* Function Body */
    d1 = cn[1] * p1[1] + cn[2] * p1[2] + cn[3] * p1[3];
    d2 = cn[1] * p2[1] + cn[2] * p2[2] + cn[3] * p2[3];

    if (d1 == d2) {
	*ier = 1;
	return 0;
    }

/* Solve for T such that <PP,CN> = 0 and compute PP and PPN. */

    t = d1 / (d1 - d2);
    ppn = 0.;
    for (i__ = 1; i__ <= 3; ++i__) {
	pp[i__ - 1] = p1[i__] + t * (p2[i__] - p1[i__]);
	ppn += pp[i__ - 1] * pp[i__ - 1];
/* L1: */
    }

/* PPN = 0 iff PP = 0 iff P2 = -P1 (and T = .5). */

    if (ppn == 0.) {
	*ier = 2;
	return 0;
    }
    ppn = sqrt(ppn);

/* Compute P = PP/PPN. */

    for (i__ = 1; i__ <= 3; ++i__) {
	p[i__] = pp[i__ - 1] / ppn;
/* L2: */
    }
    *ier = 0;
    return 0;
} /* intrsc_ */

integer jrand_(n, ix, iy, iz)
integer *n, *ix, *iy, *iz;
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static real u, x;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/28/98 */

/*   This function returns a uniformly distributed pseudo- */
/* random integer in the range 1 to N. */


/* On input: */

/*       N = Maximum value to be returned. */

/* N is not altered by this function. */

/*       IX,IY,IZ = Integer seeds initialized to values in */
/*                  the range 1 to 30,000 before the first */
/*                  call to JRAND, and not altered between */
/*                  subsequent calls (unless a sequence of */
/*                  random numbers is to be repeated by */
/*                  reinitializing the seeds). */

/* On output: */

/*       IX,IY,IZ = Updated integer seeds. */

/*       JRAND = Random integer in the range 1 to N. */

/* Reference:  B. A. Wichmann and I. D. Hill, "An Efficient */
/*             and Portable Pseudo-random Number Generator", */
/*             Applied Statistics, Vol. 31, No. 2, 1982, */
/*             pp. 188-190. */

/* Modules required by JRAND:  None */

/* Intrinsic functions called by JRAND:  INT, MOD, REAL */

/* *********************************************************** */


/* Local parameters: */

/* U = Pseudo-random number uniformly distributed in the */
/*     interval (0,1). */
/* X = Pseudo-random number in the range 0 to 3 whose frac- */
/*       tional part is U. */

    *ix = *ix * 171 % 30269;
    *iy = *iy * 172 % 30307;
    *iz = *iz * 170 % 30323;
    x = (real) (*ix) / (float)30269. + (real) (*iy) / (float)30307. + (real) (
	    *iz) / (float)30323.;
    u = x - (integer) x;
    ret_val = (real) (*n) * u + (float)1.;
    return ret_val;
} /* jrand_ */

logical left_(x1, y1, z1, x2, y2, z2, x0, y0, z0)
doublereal *x1, *y1, *z1, *x2, *y2, *z2, *x0, *y0, *z0;
{
    /* System generated locals */
    logical ret_val;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/15/96 */

/*   This function determines whether node N0 is in the */
/* (closed) left hemisphere defined by the plane containing */
/* N1, N2, and the origin, where left is defined relative to */
/* an observer at N1 facing N2. */


/* On input: */

/*       X1,Y1,Z1 = Coordinates of N1. */

/*       X2,Y2,Z2 = Coordinates of N2. */

/*       X0,Y0,Z0 = Coordinates of N0. */

/* Input parameters are not altered by this function. */

/* On output: */

/*       LEFT = TRUE if and only if N0 is in the closed */
/*              left hemisphere. */

/* Modules required by LEFT:  None */

/* *********************************************************** */

/* LEFT = TRUE iff <N0,N1 X N2> = det(N0,N1,N2) .GE. 0. */

    ret_val = *x0 * (*y1 * *z2 - *y2 * *z1) - *y0 * (*x1 * *z2 - *x2 * *z1) + 
	    *z0 * (*x1 * *y2 - *x2 * *y1) >= 0.;
    return ret_val;
} /* left_ */

integer lstptr_(lpl, nb, list, lptr)
integer *lpl, *nb, *list, *lptr;
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer nd, lp;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/15/96 */

/*   This function returns the index (LIST pointer) of NB in */
/* the adjacency list for N0, where LPL = LEND(N0). */

/*   This function is identical to the similarly named */
/* function in TRIPACK. */


/* On input: */

/*       LPL = LEND(N0) */

/*       NB = Index of the node whose pointer is to be re- */
/*            turned.  NB must be connected to N0. */

/*       LIST,LPTR = Data structure defining the triangula- */
/*                   tion.  Refer to Subroutine TRMESH. */

/* Input parameters are not altered by this function. */

/* On output: */

/*       LSTPTR = Pointer such that LIST(LSTPTR) = NB or */
/*                LIST(LSTPTR) = -NB, unless NB is not a */
/*                neighbor of N0, in which case LSTPTR = LPL. */

/* Modules required by LSTPTR:  None */

/* *********************************************************** */


/* Local parameters: */

/* LP = LIST pointer */
/* ND = Nodal index */

    /* Parameter adjustments */
    --lptr;
    --list;

    /* Function Body */
    lp = lptr[*lpl];
L1:
    nd = list[lp];
    if (nd == *nb) {
	goto L2;
    }
    lp = lptr[lp];
    if (lp != *lpl) {
	goto L1;
    }

L2:
    ret_val = lp;
    return ret_val;
} /* lstptr_ */

integer nbcnt_(lpl, lptr)
integer *lpl, *lptr;
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer k, lp;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/15/96 */

/*   This function returns the number of neighbors of a node */
/* N0 in a triangulation created by Subroutine TRMESH. */

/*   This function is identical to the similarly named */
/* function in TRIPACK. */


/* On input: */

/*       LPL = LIST pointer to the last neighbor of N0 -- */
/*             LPL = LEND(N0). */

/*       LPTR = Array of pointers associated with LIST. */

/* Input parameters are not altered by this function. */

/* On output: */

/*       NBCNT = Number of neighbors of N0. */

/* Modules required by NBCNT:  None */

/* *********************************************************** */


/* Local parameters: */

/* K =  Counter for computing the number of neighbors */
/* LP = LIST pointer */

    /* Parameter adjustments */
    --lptr;

    /* Function Body */
    lp = *lpl;
    k = 1;

L1:
    lp = lptr[lp];
    if (lp == *lpl) {
	goto L2;
    }
    ++k;
    goto L1;

L2:
    ret_val = k;
    return ret_val;
} /* nbcnt_ */

integer nearnd_(p, ist, n, x, y, z__, list, lptr, lend, al)
doublereal *p;
integer *ist, *n;
doublereal *x, *y, *z__;
integer *list, *lptr, *lend;
doublereal *al;
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Builtin functions */
    double acos();

    /* Local variables */
    static integer l;
    static doublereal b1, b2, b3;
    static integer i1, i2, i3, listp[25], n1, n2, lptrp[25], n3, lp, nn, nr;
    extern /* Subroutine */ int trfind_();
    static doublereal ds1;
    static integer lp1, lp2;
    static doublereal dx1, dx2, dx3, dy1, dy2, dy3, dz1, dz2, dz3;
    extern integer lstptr_();
    static integer lpl;
    static doublereal dsr;
    static integer nst;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/28/98 */

/*   Given a point P on the surface of the unit sphere and a */
/* Delaunay triangulation created by Subroutine TRMESH, this */
/* function returns the index of the nearest triangulation */
/* node to P. */

/*   The algorithm consists of implicitly adding P to the */
/* triangulation, finding the nearest neighbor to P, and */
/* implicitly deleting P from the triangulation.  Thus, it */
/* is based on the fact that, if P is a node in a Delaunay */
/* triangulation, the nearest node to P is a neighbor of P. */


/* On input: */

/*       P = Array of length 3 containing the Cartesian coor- */
/*           dinates of the point P to be located relative to */
/*           the triangulation.  It is assumed without a test */
/*           that P(1)**2 + P(2)**2 + P(3)**2 = 1. */

/*       IST = Index of a node at which TRFIND begins the */
/*             search.  Search time depends on the proximity */
/*             of this node to P. */

/*       N = Number of nodes in the triangulation.  N .GE. 3. */

/*       X,Y,Z = Arrays of length N containing the Cartesian */
/*               coordinates of the nodes. */

/*       LIST,LPTR,LEND = Data structure defining the trian- */
/*                        gulation.  Refer to TRMESH. */

/* Input parameters are not altered by this function. */

/* On output: */

/*       NEARND = Nodal index of the nearest node to P, or 0 */
/*                if N < 3 or the triangulation data struc- */
/*                ture is invalid. */

/*       AL = Arc length (angular distance in radians) be- */
/*            tween P and NEARND unless NEARND = 0. */

/*       Note that the number of candidates for NEARND */
/*       (neighbors of P) is limited to LMAX defined in */
/*       the PARAMETER statement below. */

/* Modules required by NEARND:  JRAND, LSTPTR, TRFIND, STORE */

/* Intrinsic functions called by NEARND:  ABS, ACOS */

/* *********************************************************** */


/* Local parameters: */

/* B1,B2,B3 =  Unnormalized barycentric coordinates returned */
/*               by TRFIND */
/* DS1 =       (Negative cosine of the) distance from P to N1 */
/* DSR =       (Negative cosine of the) distance from P to NR */
/* DX1,..DZ3 = Components of vectors used by the swap test */
/* I1,I2,I3 =  Nodal indexes of a triangle containing P, or */
/*               the rightmost (I1) and leftmost (I2) visible */
/*               boundary nodes as viewed from P */
/* L =         Length of LISTP/LPTRP and number of neighbors */
/*               of P */
/* LMAX =      Maximum value of L */
/* LISTP =     Indexes of the neighbors of P */
/* LPTRP =     Array of pointers in 1-1 correspondence with */
/*               LISTP elements */
/* LP =        LIST pointer to a neighbor of N1 and LISTP */
/*               pointer */
/* LP1,LP2 =   LISTP indexes (pointers) */
/* LPL =       Pointer to the last neighbor of N1 */
/* N1 =        Index of a node visible from P */
/* N2 =        Index of an endpoint of an arc opposite P */
/* N3 =        Index of the node opposite N1->N2 */
/* NN =        Local copy of N */
/* NR =        Index of a candidate for the nearest node to P */
/* NST =       Index of the node at which TRFIND begins the */
/*               search */


/* Store local parameters and test for N invalid. */

    /* Parameter adjustments */
    --p;
    --lend;
    --z__;
    --y;
    --x;
    --list;
    --lptr;

    /* Function Body */
    nn = *n;
    if (nn < 3) {
	goto L6;
    }
    nst = *ist;
    if (nst < 1 || nst > nn) {
	nst = 1;
    }

/* Find a triangle (I1,I2,I3) containing P, or the rightmost */
/*   (I1) and leftmost (I2) visible boundary nodes as viewed */
/*   from P. */

    trfind_(&nst, &p[1], n, &x[1], &y[1], &z__[1], &list[1], &lptr[1], &lend[
	    1], &b1, &b2, &b3, &i1, &i2, &i3);

/* Test for collinear nodes. */

    if (i1 == 0) {
	goto L6;
    }

/* Store the linked list of 'neighbors' of P in LISTP and */
/*   LPTRP.  I1 is the first neighbor, and 0 is stored as */
/*   the last neighbor if P is not contained in a triangle. */
/*   L is the length of LISTP and LPTRP, and is limited to */
/*   LMAX. */

    if (i3 != 0) {
	listp[0] = i1;
	lptrp[0] = 2;
	listp[1] = i2;
	lptrp[1] = 3;
	listp[2] = i3;
	lptrp[2] = 1;
	l = 3;
    } else {
	n1 = i1;
	l = 1;
	lp1 = 2;
	listp[l - 1] = n1;
	lptrp[l - 1] = lp1;

/*   Loop on the ordered sequence of visible boundary nodes */
/*     N1 from I1 to I2. */

L1:
	lpl = lend[n1];
	n1 = -list[lpl];
	l = lp1;
	lp1 = l + 1;
	listp[l - 1] = n1;
	lptrp[l - 1] = lp1;
	if (n1 != i2 && lp1 < 25) {
	    goto L1;
	}
	l = lp1;
	listp[l - 1] = 0;
	lptrp[l - 1] = 1;
    }

/* Initialize variables for a loop on arcs N1-N2 opposite P */
/*   in which new 'neighbors' are 'swapped' in.  N1 follows */
/*   N2 as a neighbor of P, and LP1 and LP2 are the LISTP */
/*   indexes of N1 and N2. */

    lp2 = 1;
    n2 = i1;
    lp1 = lptrp[0];
    n1 = listp[lp1 - 1];

/* Begin loop:  find the node N3 opposite N1->N2. */

L2:
    lp = lstptr_(&lend[n1], &n2, &list[1], &lptr[1]);
    if (list[lp] < 0) {
	goto L3;
    }
    lp = lptr[lp];
    n3 = (i__1 = list[lp], abs(i__1));

/* Swap test:  Exit the loop if L = LMAX. */

    if (l == 25) {
	goto L4;
    }
    dx1 = x[n1] - p[1];
    dy1 = y[n1] - p[2];
    dz1 = z__[n1] - p[3];

    dx2 = x[n2] - p[1];
    dy2 = y[n2] - p[2];
    dz2 = z__[n2] - p[3];

    dx3 = x[n3] - p[1];
    dy3 = y[n3] - p[2];
    dz3 = z__[n3] - p[3];
    if (dx3 * (dy2 * dz1 - dy1 * dz2) - dy3 * (dx2 * dz1 - dx1 * dz2) + dz3 * 
	    (dx2 * dy1 - dx1 * dy2) <= 0.) {
	goto L3;
    }

/* Swap:  Insert N3 following N2 in the adjacency list for P. */
/*        The two new arcs opposite P must be tested. */

    ++l;
    lptrp[lp2 - 1] = l;
    listp[l - 1] = n3;
    lptrp[l - 1] = lp1;
    lp1 = l;
    n1 = n3;
    goto L2;

/* No swap:  Advance to the next arc and test for termination */
/*           on N1 = I1 (LP1 = 1) or N1 followed by 0. */

L3:
    if (lp1 == 1) {
	goto L4;
    }
    lp2 = lp1;
    n2 = n1;
    lp1 = lptrp[lp1 - 1];
    n1 = listp[lp1 - 1];
    if (n1 == 0) {
	goto L4;
    }
    goto L2;

/* Set NR and DSR to the index of the nearest node to P and */
/*   an increasing function (negative cosine) of its distance */
/*   from P, respectively. */

L4:
    nr = i1;
    dsr = -(x[nr] * p[1] + y[nr] * p[2] + z__[nr] * p[3]);
    i__1 = l;
    for (lp = 2; lp <= i__1; ++lp) {
	n1 = listp[lp - 1];
	if (n1 == 0) {
	    goto L5;
	}
	ds1 = -(x[n1] * p[1] + y[n1] * p[2] + z__[n1] * p[3]);
	if (ds1 < dsr) {
	    nr = n1;
	    dsr = ds1;
	}
L5:
	;
    }
    dsr = -dsr;
    if (dsr > 1.) {
	dsr = 1.;
    }
    *al = acos(dsr);
    ret_val = nr;
    return ret_val;

/* Invalid input. */

L6:
    ret_val = 0;
    return ret_val;
} /* nearnd_ */

/* Subroutine */ int optim_(x, y, z__, na, list, lptr, lend, nit, iwk, ier)
doublereal *x, *y, *z__;
integer *na, *list, *lptr, *lend, *nit, *iwk, *ier;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer iter;
    extern /* Subroutine */ int swap_();
    static integer i__, maxit, n1, n2, lp, io1, io2;
    extern logical swptst_();
    static integer nna, lp21, lpl, lpp;
    static logical swp;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/30/98 */

/*   Given a set of NA triangulation arcs, this subroutine */
/* optimizes the portion of the triangulation consisting of */
/* the quadrilaterals (pairs of adjacent triangles) which */
/* have the arcs as diagonals by applying the circumcircle */
/* test and appropriate swaps to the arcs. */

/*   An iteration consists of applying the swap test and */
/* swaps to all NA arcs in the order in which they are */
/* stored.  The iteration is repeated until no swap occurs */
/* or NIT iterations have been performed.  The bound on the */
/* number of iterations may be necessary to prevent an */
/* infinite loop caused by cycling (reversing the effect of a */
/* previous swap) due to floating point inaccuracy when four */
/* or more nodes are nearly cocircular. */


/* On input: */

/*       X,Y,Z = Arrays containing the nodal coordinates. */

/*       NA = Number of arcs in the set.  NA .GE. 0. */

/* The above parameters are not altered by this routine. */

/*       LIST,LPTR,LEND = Data structure defining the trian- */
/*                        gulation.  Refer to Subroutine */
/*                        TRMESH. */

/*       NIT = Maximum number of iterations to be performed. */
/*             NIT = 4*NA should be sufficient.  NIT .GE. 1. */

/*       IWK = Integer array dimensioned 2 by NA containing */
/*             the nodal indexes of the arc endpoints (pairs */
/*             of endpoints are stored in columns). */

/* On output: */

/*       LIST,LPTR,LEND = Updated triangulation data struc- */
/*                        ture reflecting the swaps. */

/*       NIT = Number of iterations performed. */

/*       IWK = Endpoint indexes of the new set of arcs */
/*             reflecting the swaps. */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered. */
/*             IER = 1 if a swap occurred on the last of */
/*                     MAXIT iterations, where MAXIT is the */
/*                     value of NIT on input.  The new set */
/*                     of arcs is not necessarily optimal */
/*                     in this case. */
/*             IER = 2 if NA < 0 or NIT < 1 on input. */
/*             IER = 3 if IWK(2,I) is not a neighbor of */
/*                     IWK(1,I) for some I in the range 1 */
/*                     to NA.  A swap may have occurred in */
/*                     this case. */
/*             IER = 4 if a zero pointer was returned by */
/*                     Subroutine SWAP. */

/* Modules required by OPTIM:  LSTPTR, SWAP, SWPTST */

/* Intrinsic function called by OPTIM:  ABS */

/* *********************************************************** */


/* Local parameters: */

/* I =       Column index for IWK */
/* IO1,IO2 = Nodal indexes of the endpoints of an arc in IWK */
/* ITER =    Iteration count */
/* LP =      LIST pointer */
/* LP21 =    Parameter returned by SWAP (not used) */
/* LPL =     Pointer to the last neighbor of IO1 */
/* LPP =     Pointer to the node preceding IO2 as a neighbor */
/*             of IO1 */
/* MAXIT =   Input value of NIT */
/* N1,N2 =   Nodes opposite IO1->IO2 and IO2->IO1, */
/*             respectively */
/* NNA =     Local copy of NA */
/* SWP =     Flag set to TRUE iff a swap occurs in the */
/*             optimization loop */

    /* Parameter adjustments */
    --x;
    --y;
    --z__;
    iwk -= 3;
    --list;
    --lptr;
    --lend;

    /* Function Body */
    nna = *na;
    maxit = *nit;
    if (nna < 0 || maxit < 1) {
	goto L7;
    }

/* Initialize iteration count ITER and test for NA = 0. */

    iter = 0;
    if (nna == 0) {
	goto L5;
    }

/* Top of loop -- */
/*   SWP = TRUE iff a swap occurred in the current iteration. */

L1:
    if (iter == maxit) {
	goto L6;
    }
    ++iter;
    swp = FALSE_;

/*   Inner loop on arcs IO1-IO2 -- */

    i__1 = nna;
    for (i__ = 1; i__ <= i__1; ++i__) {
	io1 = iwk[(i__ << 1) + 1];
	io2 = iwk[(i__ << 1) + 2];

/*   Set N1 and N2 to the nodes opposite IO1->IO2 and */
/*     IO2->IO1, respectively.  Determine the following: */

/*     LPL = pointer to the last neighbor of IO1, */
/*     LP = pointer to IO2 as a neighbor of IO1, and */
/*     LPP = pointer to the node N2 preceding IO2. */

	lpl = lend[io1];
	lpp = lpl;
	lp = lptr[lpp];
L2:
	if (list[lp] == io2) {
	    goto L3;
	}
	lpp = lp;
	lp = lptr[lpp];
	if (lp != lpl) {
	    goto L2;
	}

/*   IO2 should be the last neighbor of IO1.  Test for no */
/*     arc and bypass the swap test if IO1 is a boundary */
/*     node. */

	if ((i__2 = list[lp], abs(i__2)) != io2) {
	    goto L8;
	}
	if (list[lp] < 0) {
	    goto L4;
	}

/*   Store N1 and N2, or bypass the swap test if IO1 is a */
/*     boundary node and IO2 is its first neighbor. */

L3:
	n2 = list[lpp];
	if (n2 < 0) {
	    goto L4;
	}
	lp = lptr[lp];
	n1 = (i__2 = list[lp], abs(i__2));

/*   Test IO1-IO2 for a swap, and update IWK if necessary. */

	if (! swptst_(&n1, &n2, &io1, &io2, &x[1], &y[1], &z__[1])) {
	    goto L4;
	}
	swap_(&n1, &n2, &io1, &io2, &list[1], &lptr[1], &lend[1], &lp21);
	if (lp21 == 0) {
	    goto L9;
	}
	swp = TRUE_;
	iwk[(i__ << 1) + 1] = n1;
	iwk[(i__ << 1) + 2] = n2;
L4:
	;
    }
    if (swp) {
	goto L1;
    }

/* Successful termination. */

L5:
    *nit = iter;
    *ier = 0;
    return 0;

/* MAXIT iterations performed without convergence. */

L6:
    *nit = maxit;
    *ier = 1;
    return 0;

/* Invalid input parameter. */

L7:
    *nit = 0;
    *ier = 2;
    return 0;

/* IO2 is not a neighbor of IO1. */

L8:
    *nit = iter;
    *ier = 3;
    return 0;

/* Zero pointer returned by SWAP. */

L9:
    *nit = iter;
    *ier = 4;
    return 0;
} /* optim_ */

/* Subroutine */ int projct_(px, py, pz, ox, oy, oz, ex, ey, ez, vx, vy, vz, 
	init, x, y, z__, ier)
doublereal *px, *py, *pz, *ox, *oy, *oz, *ex, *ey, *ez, *vx, *vy, *vz;
logical *init;
doublereal *x, *y, *z__;
integer *ier;
{
    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal s, sc, xe, ye, ze, xh, yh, zh, xv, yv, zv, xw, yw, zw, 
	    oes, xoe, yoe, zoe, xep, yep, zep;


/* *********************************************************** */

/*                        From PLTPACK, SCRPLOT, and STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/18/90 */

/*   Given a projection plane and associated coordinate sys- */
/* tem defined by an origin O, eye position E, and up-vector */
/* V, this subroutine applies a perspective depth transform- */
/* ation T to a point P = (PX,PY,PZ), returning the point */
/* T(P) = (X,Y,Z), where X and Y are the projection plane */
/* coordinates of the point that lies in the projection */
/* plane and on the line defined by P and E, and Z is the */
/* depth associated with P. */

/*   The projection plane is defined to be the plane that */
/* contains O and has normal defined by O and E. */

/*   The depth Z is defined in such a way that Z < 1, T maps */
/* lines to lines (and planes to planes), and if two distinct */
/* points have the same projection plane coordinates, then */
/* the one closer to E has a smaller depth.  (Z increases */
/* monotonically with orthogonal distance from P to the plane */
/* that is parallel to the projection plane and contains E.) */
/* This depth value facilitates depth sorting and depth buf- */
/* fer methods. */


/* On input: */

/*       PX,PY,PZ = Cartesian coordinates of the point P to */
/*                  be mapped onto the projection plane.  The */
/*                  half line that contains P and has end- */
/*                  point at E must intersect the plane. */

/*       OX,OY,OZ = Coordinates of O (the origin of a coordi- */
/*                  nate system in the projection plane).  A */
/*                  reasonable value for O is a point near */
/*                  the center of an object or scene to be */
/*                  viewed. */

/*       EX,EY,EZ = Coordinates of the eye-position E defin- */
/*                  ing the normal to the plane and the line */
/*                  of sight for the projection.  E must not */
/*                  coincide with O or P, and the angle be- */
/*                  tween the vectors O-E and P-E must be */
/*                  less than 90 degrees.  Note that E and P */
/*                  may lie on opposite sides of the projec- */
/*                  tion plane. */

/*       VX,VY,VZ = Coordinates of a point V which defines */
/*                  the positive Y axis of an X-Y coordinate */
/*                  system in the projection plane as the */
/*                  half-line containing O and the projection */
/*                  of O+V onto the plane.  The positive X */
/*                  axis has direction defined by the cross */
/*                  product V X (E-O). */

/* The above parameters are not altered by this routine. */

/*       INIT = Logical switch which must be set to TRUE on */
/*              the first call and when the values of O, E, */
/*              or V have been altered since a previous call. */
/*              If INIT = FALSE, it is assumed that only the */
/*              coordinates of P have changed since a previ- */
/*              ous call.  Previously stored quantities are */
/*              used for increased efficiency in this case. */

/* On output: */

/*       INIT = Switch with value reset to FALSE if IER = 0. */

/*       X,Y = Projection plane coordinates of the point */
/*             that lies in the projection plane and on the */
/*             line defined by E and P.  X and Y are not */
/*             altered if IER .NE. 0. */

/*       Z = Depth value defined above unless IER .NE. 0. */

/*       IER = Error indicator. */
/*             IER = 0 if no errors were encountered. */
/*             IER = 1 if the inner product of O-E with P-E */
/*                     is not positive, implying that E is */
/*                     too close to the plane. */
/*             IER = 2 if O, E, and O+V are collinear.  See */
/*                     the description of VX,VY,VZ. */

/* Modules required by PROJCT:  None */

/* Intrinsic function called by PROJCT:  SQRT */

/* *********************************************************** */


/* Local parameters: */

/* OES =         Norm squared of OE -- inner product (OE,OE) */
/* S =           Scale factor for computing projections */
/* SC =          Scale factor for normalizing VN and HN */
/* XE,YE,ZE =    Local copies of EX, EY, EZ */
/* XEP,YEP,ZEP = Components of the vector EP from E to P */
/* XH,YH,ZH =    Components of a unit vector HN defining the */
/*                 positive X-axis in the plane */
/* XOE,YOE,ZOE = Components of the vector OE from O to E */
/* XV,YV,ZV =    Components of a unit vector VN defining the */
/*                 positive Y-axis in the plane */
/* XW,YW,ZW =    Components of the vector W from O to the */
/*                 projection of P onto the plane */

    if (*init) {

/* Compute parameters defining the transformation: */
/*   17 adds, 27 multiplies, 3 divides, 2 compares, and */
/*   2 square roots. */

/* Set the coordinates of E to local variables, compute */
/*   OE = E-O and OES, and test for OE = 0. */

	xe = *ex;
	ye = *ey;
	ze = *ez;
	xoe = xe - *ox;
	yoe = ye - *oy;
	zoe = ze - *oz;
	oes = xoe * xoe + yoe * yoe + zoe * zoe;
	if (oes == 0.) {
	    goto L1;
	}

/* Compute S = (OE,V)/OES and VN = V - S*OE. */

	s = (xoe * *vx + yoe * *vy + zoe * *vz) / oes;
	xv = *vx - s * xoe;
	yv = *vy - s * yoe;
	zv = *vz - s * zoe;

/* Normalize VN to a unit vector. */

	sc = xv * xv + yv * yv + zv * zv;
	if (sc == 0.) {
	    goto L2;
	}
	sc = 1. / sqrt(sc);
	xv = sc * xv;
	yv = sc * yv;
	zv = sc * zv;

/* Compute HN = VN X OE (normalized). */

	xh = yv * zoe - yoe * zv;
	yh = xoe * zv - xv * zoe;
	zh = xv * yoe - xoe * yv;
	sc = sqrt(xh * xh + yh * yh + zh * zh);
	if (sc == 0.) {
	    goto L2;
	}
	sc = 1. / sc;
	xh = sc * xh;
	yh = sc * yh;
	zh = sc * zh;
    }

/* Apply the transformation:  13 adds, 12 multiplies, */
/*                            1 divide, and 1 compare. */

/* Compute EP = P-E, S = OES/(OE,EP), and W = OE - S*EP. */

    xep = *px - xe;
    yep = *py - ye;
    zep = *pz - ze;
    s = xoe * xep + yoe * yep + zoe * zep;
    if (s >= 0.) {
	goto L1;
    }
    s = oes / s;
    xw = xoe - s * xep;
    yw = yoe - s * yep;
    zw = zoe - s * zep;

/* Map W into X = (W,HN), Y = (W,VN), compute Z = 1+S, and */
/*   reset INIT. */

    *x = xw * xh + yw * yh + zw * zh;
    *y = xw * xv + yw * yv + zw * zv;
    *z__ = s + 1.;
    *init = FALSE_;
    *ier = 0;
    return 0;

/* (OE,EP) .GE. 0. */

L1:
    *ier = 1;
    return 0;

/* O, E, and O+V are collinear. */

L2:
    *ier = 2;
    return 0;
} /* projct_ */

/* Subroutine */ int scoord_(px, py, pz, plat, plon, pnrm)
doublereal *px, *py, *pz, *plat, *plon, *pnrm;
{
    /* Builtin functions */
    double sqrt(), atan2(), asin();


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   08/27/90 */

/*   This subroutine converts a point P from Cartesian coor- */
/* dinates to spherical coordinates. */


/* On input: */

/*       PX,PY,PZ = Cartesian coordinates of P. */

/* Input parameters are not altered by this routine. */

/* On output: */

/*       PLAT = Latitude of P in the range -PI/2 to PI/2, or */
/*              0 if PNRM = 0.  PLAT should be scaled by */
/*              180/PI to obtain the value in degrees. */

/*       PLON = Longitude of P in the range -PI to PI, or 0 */
/*              if P lies on the Z-axis.  PLON should be */
/*              scaled by 180/PI to obtain the value in */
/*              degrees. */

/*       PNRM = Magnitude (Euclidean norm) of P. */

/* Modules required by SCOORD:  None */

/* Intrinsic functions called by SCOORD:  ASIN, ATAN2, SQRT */

/* *********************************************************** */

    *pnrm = sqrt(*px * *px + *py * *py + *pz * *pz);
    if (*px != 0. || *py != 0.) {
	*plon = atan2(*py, *px);
    } else {
	*plon = 0.;
    }
    if (*pnrm != 0.) {
	*plat = asin(*pz / *pnrm);
    } else {
	*plat = 0.;
    }
    return 0;
} /* scoord_ */

doublereal store_(x)
doublereal *x;
{
    /* System generated locals */
    doublereal ret_val;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   05/09/92 */

/*   This function forces its argument X to be stored in a */
/* memory location, thus providing a means of determining */
/* floating point number characteristics (such as the machine */
/* precision) when it is necessary to avoid computation in */
/* high precision registers. */


/* On input: */

/*       X = Value to be stored. */

/* X is not altered by this function. */

/* On output: */

/*       STORE = Value of X after it has been stored and */
/*               possibly truncated or rounded to the single */
/*               precision word length. */

/* Modules required by STORE:  None */

/* *********************************************************** */

    stcom_1.y = *x;
    ret_val = stcom_1.y;
    return ret_val;
} /* store_ */

/* Subroutine */ int swap_(in1, in2, io1, io2, list, lptr, lend, lp21)
integer *in1, *in2, *io1, *io2, *list, *lptr, *lend, *lp21;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer lpsav, lp;
    extern integer lstptr_();
    static integer lph;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   06/22/98 */

/*   Given a triangulation of a set of points on the unit */
/* sphere, this subroutine replaces a diagonal arc in a */
/* strictly convex quadrilateral (defined by a pair of adja- */
/* cent triangles) with the other diagonal.  Equivalently, a */
/* pair of adjacent triangles is replaced by another pair */
/* having the same union. */


/* On input: */

/*       IN1,IN2,IO1,IO2 = Nodal indexes of the vertices of */
/*                         the quadrilateral.  IO1-IO2 is re- */
/*                         placed by IN1-IN2.  (IO1,IO2,IN1) */
/*                         and (IO2,IO1,IN2) must be trian- */
/*                         gles on input. */

/* The above parameters are not altered by this routine. */

/*       LIST,LPTR,LEND = Data structure defining the trian- */
/*                        gulation.  Refer to Subroutine */
/*                        TRMESH. */

/* On output: */

/*       LIST,LPTR,LEND = Data structure updated with the */
/*                        swap -- triangles (IO1,IO2,IN1) and */
/*                        (IO2,IO1,IN2) are replaced by */
/*                        (IN1,IN2,IO2) and (IN2,IN1,IO1) */
/*                        unless LP21 = 0. */

/*       LP21 = Index of IN1 as a neighbor of IN2 after the */
/*              swap is performed unless IN1 and IN2 are */
/*              adjacent on input, in which case LP21 = 0. */

/* Module required by SWAP:  LSTPTR */

/* Intrinsic function called by SWAP:  ABS */

/* *********************************************************** */


/* Local parameters: */

/* LP,LPH,LPSAV = LIST pointers */


/* Test for IN1 and IN2 adjacent. */

    /* Parameter adjustments */
    --lend;
    --lptr;
    --list;

    /* Function Body */
    lp = lstptr_(&lend[*in1], in2, &list[1], &lptr[1]);
    if ((i__1 = list[lp], abs(i__1)) == *in2) {
	*lp21 = 0;
	return 0;
    }

/* Delete IO2 as a neighbor of IO1. */

    lp = lstptr_(&lend[*io1], in2, &list[1], &lptr[1]);
    lph = lptr[lp];
    lptr[lp] = lptr[lph];

/* If IO2 is the last neighbor of IO1, make IN2 the */
/*   last neighbor. */

    if (lend[*io1] == lph) {
	lend[*io1] = lp;
    }

/* Insert IN2 as a neighbor of IN1 following IO1 */
/*   using the hole created above. */

    lp = lstptr_(&lend[*in1], io1, &list[1], &lptr[1]);
    lpsav = lptr[lp];
    lptr[lp] = lph;
    list[lph] = *in2;
    lptr[lph] = lpsav;

/* Delete IO1 as a neighbor of IO2. */

    lp = lstptr_(&lend[*io2], in1, &list[1], &lptr[1]);
    lph = lptr[lp];
    lptr[lp] = lptr[lph];

/* If IO1 is the last neighbor of IO2, make IN1 the */
/*   last neighbor. */

    if (lend[*io2] == lph) {
	lend[*io2] = lp;
    }

/* Insert IN1 as a neighbor of IN2 following IO2. */

    lp = lstptr_(&lend[*in2], io2, &list[1], &lptr[1]);
    lpsav = lptr[lp];
    lptr[lp] = lph;
    list[lph] = *in1;
    lptr[lph] = lpsav;
    *lp21 = lph;
    return 0;
} /* swap_ */

logical swptst_(n1, n2, n3, n4, x, y, z__)
integer *n1, *n2, *n3, *n4;
doublereal *x, *y, *z__;
{
    /* System generated locals */
    logical ret_val;

    /* Local variables */
    static doublereal x4, y4, z4, dx1, dx2, dx3, dy1, dy2, dy3, dz1, dz2, dz3;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   03/29/91 */

/*   This function decides whether or not to replace a */
/* diagonal arc in a quadrilateral with the other diagonal. */
/* The decision will be to swap (SWPTST = TRUE) if and only */
/* if N4 lies above the plane (in the half-space not contain- */
/* ing the origin) defined by (N1,N2,N3), or equivalently, if */
/* the projection of N4 onto this plane is interior to the */
/* circumcircle of (N1,N2,N3).  The decision will be for no */
/* swap if the quadrilateral is not strictly convex. */


/* On input: */

/*       N1,N2,N3,N4 = Indexes of the four nodes defining the */
/*                     quadrilateral with N1 adjacent to N2, */
/*                     and (N1,N2,N3) in counterclockwise */
/*                     order.  The arc connecting N1 to N2 */
/*                     should be replaced by an arc connec- */
/*                     ting N3 to N4 if SWPTST = TRUE.  Refer */
/*                     to Subroutine SWAP. */

/*       X,Y,Z = Arrays of length N containing the Cartesian */
/*               coordinates of the nodes.  (X(I),Y(I),Z(I)) */
/*               define node I for I = N1, N2, N3, and N4. */

/* Input parameters are not altered by this routine. */

/* On output: */

/*       SWPTST = TRUE if and only if the arc connecting N1 */
/*                and N2 should be swapped for an arc con- */
/*                necting N3 and N4. */

/* Modules required by SWPTST:  None */

/* *********************************************************** */


/* Local parameters: */

/* DX1,DY1,DZ1 = Coordinates of N4->N1 */
/* DX2,DY2,DZ2 = Coordinates of N4->N2 */
/* DX3,DY3,DZ3 = Coordinates of N4->N3 */
/* X4,Y4,Z4 =    Coordinates of N4 */

    /* Parameter adjustments */
    --z__;
    --y;
    --x;

    /* Function Body */
    x4 = x[*n4];
    y4 = y[*n4];
    z4 = z__[*n4];
    dx1 = x[*n1] - x4;
    dx2 = x[*n2] - x4;
    dx3 = x[*n3] - x4;
    dy1 = y[*n1] - y4;
    dy2 = y[*n2] - y4;
    dy3 = y[*n3] - y4;
    dz1 = z__[*n1] - z4;
    dz2 = z__[*n2] - z4;
    dz3 = z__[*n3] - z4;

/* N4 lies above the plane of (N1,N2,N3) iff N3 lies above */
/*   the plane of (N2,N1,N4) iff Det(N3-N4,N2-N4,N1-N4) = */
/*   (N3-N4,N2-N4 X N1-N4) > 0. */

    ret_val = dx3 * (dy2 * dz1 - dy1 * dz2) - dy3 * (dx2 * dz1 - dx1 * dz2) + 
	    dz3 * (dx2 * dy1 - dx1 * dy2) > 0.;
    return ret_val;
} /* swptst_ */

/* Subroutine */ int trans_(n, rlat, rlon, x, y, z__)
integer *n;
doublereal *rlat, *rlon, *x, *y, *z__;
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(), sin();

    /* Local variables */
    static integer i__;
    static doublereal theta;
    static integer nn;
    static doublereal cosphi, phi;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   04/08/90 */

/*   This subroutine transforms spherical coordinates into */
/* Cartesian coordinates on the unit sphere for input to */
/* Subroutine TRMESH.  Storage for X and Y may coincide with */
/* storage for RLAT and RLON if the latter need not be saved. */


/* On input: */

/*       N = Number of nodes (points on the unit sphere) */
/*           whose coordinates are to be transformed. */

/*       RLAT = Array of length N containing latitudinal */
/*              coordinates of the nodes in radians. */

/*       RLON = Array of length N containing longitudinal */
/*              coordinates of the nodes in radians. */

/* The above parameters are not altered by this routine. */

/*       X,Y,Z = Arrays of length at least N. */

/* On output: */

/*       X,Y,Z = Cartesian coordinates in the range -1 to 1. */
/*               X(I)**2 + Y(I)**2 + Z(I)**2 = 1 for I = 1 */
/*               to N. */

/* Modules required by TRANS:  None */

/* Intrinsic functions called by TRANS:  COS, SIN */

/* *********************************************************** */


/* Local parameters: */

/* COSPHI = cos(PHI) */
/* I =      DO-loop index */
/* NN =     Local copy of N */
/* PHI =    Latitude */
/* THETA =  Longitude */

    /* Parameter adjustments */
    --z__;
    --y;
    --x;
    --rlon;
    --rlat;

    /* Function Body */
    nn = *n;
    i__1 = nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	phi = rlat[i__];
	theta = rlon[i__];
	cosphi = cos(phi);
	x[i__] = cosphi * cos(theta);
	y[i__] = cosphi * sin(theta);
	z__[i__] = sin(phi);
/* L1: */
    }
    return 0;
} /* trans_ */

/* Subroutine */ int trfind_(nst, p, n, x, y, z__, list, lptr, lend, b1, b2, 
	b3, i1, i2, i3)
integer *nst;
doublereal *p;
integer *n;
doublereal *x, *y, *z__;
integer *list, *lptr, *lend;
doublereal *b1, *b2, *b3;
integer *i1, *i2, *i3;
{
    /* Initialized data */

    static integer ix = 1;
    static integer iy = 2;
    static integer iz = 3;

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer next;
    static doublereal q[3];
    extern integer jrand_();
    extern doublereal store_();
    static integer n0, n1, n2, n3, n4, nf;
    static doublereal s12;
    static integer nl, lp;
    static doublereal xp, yp, zp;
    static integer n1s, n2s;
    extern integer lstptr_();
    static doublereal eps, tol, ptn1, ptn2;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   11/30/99 */

/*   This subroutine locates a point P relative to a triangu- */
/* lation created by Subroutine TRMESH.  If P is contained in */
/* a triangle, the three vertex indexes and barycentric coor- */
/* dinates are returned.  Otherwise, the indexes of the */
/* visible boundary nodes are returned. */


/* On input: */

/*       NST = Index of a node at which TRFIND begins its */
/*             search.  Search time depends on the proximity */
/*             of this node to P. */

/*       P = Array of length 3 containing the x, y, and z */
/*           coordinates (in that order) of the point P to be */
/*           located. */

/*       N = Number of nodes in the triangulation.  N .GE. 3. */

/*       X,Y,Z = Arrays of length N containing the Cartesian */
/*               coordinates of the triangulation nodes (unit */
/*               vectors).  (X(I),Y(I),Z(I)) defines node I */
/*               for I = 1 to N. */

/*       LIST,LPTR,LEND = Data structure defining the trian- */
/*                        gulation.  Refer to Subroutine */
/*                        TRMESH. */

/* Input parameters are not altered by this routine. */

/* On output: */

/*       B1,B2,B3 = Unnormalized barycentric coordinates of */
/*                  the central projection of P onto the un- */
/*                  derlying planar triangle if P is in the */
/*                  convex hull of the nodes.  These parame- */
/*                  ters are not altered if I1 = 0. */

/*       I1,I2,I3 = Counterclockwise-ordered vertex indexes */
/*                  of a triangle containing P if P is con- */
/*                  tained in a triangle.  If P is not in the */
/*                  convex hull of the nodes, I1 and I2 are */
/*                  the rightmost and leftmost (boundary) */
/*                  nodes that are visible from P, and */
/*                  I3 = 0.  (If all boundary nodes are vis- */
/*                  ible from P, then I1 and I2 coincide.) */
/*                  I1 = I2 = I3 = 0 if P and all of the */
/*                  nodes are coplanar (lie on a common great */
/*                  circle. */

/* Modules required by TRFIND:  JRAND, LSTPTR, STORE */

/* Intrinsic function called by TRFIND:  ABS */

/* *********************************************************** */


    /* Parameter adjustments */
    --p;
    --lend;
    --z__;
    --y;
    --x;
    --list;
    --lptr;

    /* Function Body */

/* Local parameters: */

/* EPS =      Machine precision */
/* IX,IY,IZ = Integer seeds for JRAND */
/* LP =       LIST pointer */
/* N0,N1,N2 = Nodes in counterclockwise order defining a */
/*              cone (with vertex N0) containing P, or end- */
/*              points of a boundary edge such that P Right */
/*              N1->N2 */
/* N1S,N2S =  Initially-determined values of N1 and N2 */
/* N3,N4 =    Nodes opposite N1->N2 and N2->N1, respectively */
/* NEXT =     Candidate for I1 or I2 when P is exterior */
/* NF,NL =    First and last neighbors of N0, or first */
/*              (rightmost) and last (leftmost) nodes */
/*              visible from P when P is exterior to the */
/*              triangulation */
/* PTN1 =     Scalar product <P,N1> */
/* PTN2 =     Scalar product <P,N2> */
/* Q =        (N2 X N1) X N2  or  N1 X (N2 X N1) -- used in */
/*              the boundary traversal when P is exterior */
/* S12 =      Scalar product <N1,N2> */
/* TOL =      Tolerance (multiple of EPS) defining an upper */
/*              bound on the magnitude of a negative bary- */
/*              centric coordinate (B1 or B2) for P in a */
/*              triangle -- used to avoid an infinite number */
/*              of restarts with 0 <= B3 < EPS and B1 < 0 or */
/*              B2 < 0 but small in magnitude */
/* XP,YP,ZP = Local variables containing P(1), P(2), and P(3) */
/* X0,Y0,Z0 = Dummy arguments for DET */
/* X1,Y1,Z1 = Dummy arguments for DET */
/* X2,Y2,Z2 = Dummy arguments for DET */

/* Statement function: */

/* DET(X1,...,Z0) .GE. 0 if and only if (X0,Y0,Z0) is in the */
/*                       (closed) left hemisphere defined by */
/*                       the plane containing (0,0,0), */
/*                       (X1,Y1,Z1), and (X2,Y2,Z2), where */
/*                       left is defined relative to an ob- */
/*                       server at (X1,Y1,Z1) facing */
/*                       (X2,Y2,Z2). */


/* Initialize variables. */

    xp = p[1];
    yp = p[2];
    zp = p[3];
    n0 = *nst;
    if (n0 < 1 || n0 > *n) {
	n0 = jrand_(n, &ix, &iy, &iz);
    }

/* Compute the relative machine precision EPS and TOL. */

    eps = 1.;
L1:
    eps /= 2.;
    d__1 = eps + 1.;
    if (store_(&d__1) > 1.) {
	goto L1;
    }
    eps *= 2.;
    tol = eps * 4.;

/* Set NF and NL to the first and last neighbors of N0, and */
/*   initialize N1 = NF. */

L2:
    lp = lend[n0];
    nl = list[lp];
    lp = lptr[lp];
    nf = list[lp];
    n1 = nf;

/* Find a pair of adjacent neighbors N1,N2 of N0 that define */
/*   a wedge containing P:  P LEFT N0->N1 and P RIGHT N0->N2. */

    if (nl > 0) {

/*   N0 is an interior node.  Find N1. */

L3:
	if (xp * (y[n0] * z__[n1] - y[n1] * z__[n0]) - yp * (x[n0] * z__[n1] 
		- x[n1] * z__[n0]) + zp * (x[n0] * y[n1] - x[n1] * y[n0]) < 
		0.) {
	    lp = lptr[lp];
	    n1 = list[lp];
	    if (n1 == nl) {
		goto L6;
	    }
	    goto L3;
	}
    } else {

/*   N0 is a boundary node.  Test for P exterior. */

	nl = -nl;
	if (xp * (y[n0] * z__[nf] - y[nf] * z__[n0]) - yp * (x[n0] * z__[nf] 
		- x[nf] * z__[n0]) + zp * (x[n0] * y[nf] - x[nf] * y[n0]) < 
		0.) {

/*   P is to the right of the boundary edge N0->NF. */

	    n1 = n0;
	    n2 = nf;
	    goto L9;
	}
	if (xp * (y[nl] * z__[n0] - y[n0] * z__[nl]) - yp * (x[nl] * z__[n0] 
		- x[n0] * z__[nl]) + zp * (x[nl] * y[n0] - x[n0] * y[nl]) < 
		0.) {

/*   P is to the right of the boundary edge NL->N0. */

	    n1 = nl;
	    n2 = n0;
	    goto L9;
	}
    }

/* P is to the left of arcs N0->N1 and NL->N0.  Set N2 to the */
/*   next neighbor of N0 (following N1). */

L4:
    lp = lptr[lp];
    n2 = (i__1 = list[lp], abs(i__1));
    if (xp * (y[n0] * z__[n2] - y[n2] * z__[n0]) - yp * (x[n0] * z__[n2] - x[
	    n2] * z__[n0]) + zp * (x[n0] * y[n2] - x[n2] * y[n0]) < 0.) {
	goto L7;
    }
    n1 = n2;
    if (n1 != nl) {
	goto L4;
    }
    if (xp * (y[n0] * z__[nf] - y[nf] * z__[n0]) - yp * (x[n0] * z__[nf] - x[
	    nf] * z__[n0]) + zp * (x[n0] * y[nf] - x[nf] * y[n0]) < 0.) {
	goto L6;
    }

/* P is left of or on arcs N0->NB for all neighbors NB */
/*   of N0.  Test for P = +/-N0. */

    d__2 = (d__1 = x[n0] * xp + y[n0] * yp + z__[n0] * zp, abs(d__1));
    if (store_(&d__2) < 1. - eps * 4.) {

/*   All points are collinear iff P Left NB->N0 for all */
/*     neighbors NB of N0.  Search the neighbors of N0. */
/*     Note:  N1 = NL and LP points to NL. */

L5:
	if (xp * (y[n1] * z__[n0] - y[n0] * z__[n1]) - yp * (x[n1] * z__[n0] 
		- x[n0] * z__[n1]) + zp * (x[n1] * y[n0] - x[n0] * y[n1]) >= 
		0.) {
	    lp = lptr[lp];
	    n1 = (i__1 = list[lp], abs(i__1));
	    if (n1 == nl) {
		goto L14;
	    }
	    goto L5;
	}
    }

/* P is to the right of N1->N0, or P = +/-N0.  Set N0 to N1 */
/*   and start over. */

    n0 = n1;
    goto L2;

/* P is between arcs N0->N1 and N0->NF. */

L6:
    n2 = nf;

/* P is contained in a wedge defined by geodesics N0-N1 and */
/*   N0-N2, where N1 is adjacent to N2.  Save N1 and N2 to */
/*   test for cycling. */

L7:
    n3 = n0;
    n1s = n1;
    n2s = n2;

/* Top of edge-hopping loop: */

L8:
    *b3 = xp * (y[n1] * z__[n2] - y[n2] * z__[n1]) - yp * (x[n1] * z__[n2] - 
	    x[n2] * z__[n1]) + zp * (x[n1] * y[n2] - x[n2] * y[n1]);
    if (*b3 < 0.) {

/*   Set N4 to the first neighbor of N2 following N1 (the */
/*     node opposite N2->N1) unless N1->N2 is a boundary arc. */

	lp = lstptr_(&lend[n2], &n1, &list[1], &lptr[1]);
	if (list[lp] < 0) {
	    goto L9;
	}
	lp = lptr[lp];
	n4 = (i__1 = list[lp], abs(i__1));

/*   Define a new arc N1->N2 which intersects the geodesic */
/*     N0-P. */

	if (xp * (y[n0] * z__[n4] - y[n4] * z__[n0]) - yp * (x[n0] * z__[n4] 
		- x[n4] * z__[n0]) + zp * (x[n0] * y[n4] - x[n4] * y[n0]) < 
		0.) {
	    n3 = n2;
	    n2 = n4;
	    n1s = n1;
	    if (n2 != n2s && n2 != n0) {
		goto L8;
	    }
	} else {
	    n3 = n1;
	    n1 = n4;
	    n2s = n2;
	    if (n1 != n1s && n1 != n0) {
		goto L8;
	    }
	}

/*   The starting node N0 or edge N1-N2 was encountered */
/*     again, implying a cycle (infinite loop).  Restart */
/*     with N0 randomly selected. */

	n0 = jrand_(n, &ix, &iy, &iz);
	goto L2;
    }

/* P is in (N1,N2,N3) unless N0, N1, N2, and P are collinear */
/*   or P is close to -N0. */

    if (*b3 >= eps) {

/*   B3 .NE. 0. */

	*b1 = xp * (y[n2] * z__[n3] - y[n3] * z__[n2]) - yp * (x[n2] * z__[n3]
		 - x[n3] * z__[n2]) + zp * (x[n2] * y[n3] - x[n3] * y[n2]);
	*b2 = xp * (y[n3] * z__[n1] - y[n1] * z__[n3]) - yp * (x[n3] * z__[n1]
		 - x[n1] * z__[n3]) + zp * (x[n3] * y[n1] - x[n1] * y[n3]);
	if (*b1 < -tol || *b2 < -tol) {

/*   Restart with N0 randomly selected. */

	    n0 = jrand_(n, &ix, &iy, &iz);
	    goto L2;
	}
    } else {

/*   B3 = 0 and thus P lies on N1->N2. Compute */
/*     B1 = Det(P,N2 X N1,N2) and B2 = Det(P,N1,N2 X N1). */

	*b3 = 0.;
	s12 = x[n1] * x[n2] + y[n1] * y[n2] + z__[n1] * z__[n2];
	ptn1 = xp * x[n1] + yp * y[n1] + zp * z__[n1];
	ptn2 = xp * x[n2] + yp * y[n2] + zp * z__[n2];
	*b1 = ptn1 - s12 * ptn2;
	*b2 = ptn2 - s12 * ptn1;
	if (*b1 < -tol || *b2 < -tol) {

/*   Restart with N0 randomly selected. */

	    n0 = jrand_(n, &ix, &iy, &iz);
	    goto L2;
	}
    }

/* P is in (N1,N2,N3). */

    *i1 = n1;
    *i2 = n2;
    *i3 = n3;
    if (*b1 < (float)0.) {
	*b1 = (float)0.;
    }
    if (*b2 < (float)0.) {
	*b2 = (float)0.;
    }
    return 0;

/* P Right N1->N2, where N1->N2 is a boundary edge. */
/*   Save N1 and N2, and set NL = 0 to indicate that */
/*   NL has not yet been found. */

L9:
    n1s = n1;
    n2s = n2;
    nl = 0;

/*           Counterclockwise Boundary Traversal: */

L10:
    lp = lend[n2];
    lp = lptr[lp];
    next = list[lp];
    if (xp * (y[n2] * z__[next] - y[next] * z__[n2]) - yp * (x[n2] * z__[next]
	     - x[next] * z__[n2]) + zp * (x[n2] * y[next] - x[next] * y[n2]) 
	    >= 0.) {

/*   N2 is the rightmost visible node if P Forward N2->N1 */
/*     or NEXT Forward N2->N1.  Set Q to (N2 X N1) X N2. */

	s12 = x[n1] * x[n2] + y[n1] * y[n2] + z__[n1] * z__[n2];
	q[0] = x[n1] - s12 * x[n2];
	q[1] = y[n1] - s12 * y[n2];
	q[2] = z__[n1] - s12 * z__[n2];
	if (xp * q[0] + yp * q[1] + zp * q[2] >= 0.) {
	    goto L11;
	}
	if (x[next] * q[0] + y[next] * q[1] + z__[next] * q[2] >= 0.) {
	    goto L11;
	}

/*   N1, N2, NEXT, and P are nearly collinear, and N2 is */
/*     the leftmost visible node. */

	nl = n2;
    }

/* Bottom of counterclockwise loop: */

    n1 = n2;
    n2 = next;
    if (n2 != n1s) {
	goto L10;
    }

/* All boundary nodes are visible from P. */

    *i1 = n1s;
    *i2 = n1s;
    *i3 = 0;
    return 0;

/* N2 is the rightmost visible node. */

L11:
    nf = n2;
    if (nl == 0) {

/* Restore initial values of N1 and N2, and begin the search */
/*   for the leftmost visible node. */

	n2 = n2s;
	n1 = n1s;

/*           Clockwise Boundary Traversal: */

L12:
	lp = lend[n1];
	next = -list[lp];
	if (xp * (y[next] * z__[n1] - y[n1] * z__[next]) - yp * (x[next] * 
		z__[n1] - x[n1] * z__[next]) + zp * (x[next] * y[n1] - x[n1] *
		 y[next]) >= 0.) {

/*   N1 is the leftmost visible node if P or NEXT is */
/*     forward of N1->N2.  Compute Q = N1 X (N2 X N1). */

	    s12 = x[n1] * x[n2] + y[n1] * y[n2] + z__[n1] * z__[n2];
	    q[0] = x[n2] - s12 * x[n1];
	    q[1] = y[n2] - s12 * y[n1];
	    q[2] = z__[n2] - s12 * z__[n1];
	    if (xp * q[0] + yp * q[1] + zp * q[2] >= 0.) {
		goto L13;
	    }
	    if (x[next] * q[0] + y[next] * q[1] + z__[next] * q[2] >= 0.) {
		goto L13;
	    }

/*   P, NEXT, N1, and N2 are nearly collinear and N1 is the */
/*     rightmost visible node. */

	    nf = n1;
	}

/* Bottom of clockwise loop: */

	n2 = n1;
	n1 = next;
	if (n1 != n1s) {
	    goto L12;
	}

/* All boundary nodes are visible from P. */

	*i1 = n1;
	*i2 = n1;
	*i3 = 0;
	return 0;

/* N1 is the leftmost visible node. */

L13:
	nl = n1;
    }

/* NF and NL have been found. */

    *i1 = nf;
    *i2 = nl;
    *i3 = 0;
    return 0;

/* All points are collinear (coplanar). */

L14:
    *i1 = 0;
    *i2 = 0;
    *i3 = 0;
    return 0;
} /* trfind_ */

/* Subroutine */ int trlist_(n, list, lptr, lend, nrow, nt, ltri, ier)
integer *n, *list, *lptr, *lend, *nrow, *nt, *ltri, *ier;
{
    /* System generated locals */
    integer ltri_dim1, ltri_offset, i__1, i__2;

    /* Local variables */
    static logical arcs;
    static integer lpln1, i__, j, i1, i2, i3, n1, n2, n3, ka, kn, lp, kt, nm2,
	     lp2, lpl, isv;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/20/96 */

/*   This subroutine converts a triangulation data structure */
/* from the linked list created by Subroutine TRMESH to a */
/* triangle list. */

/* On input: */

/*       N = Number of nodes in the triangulation.  N .GE. 3. */

/*       LIST,LPTR,LEND = Linked list data structure defin- */
/*                        ing the triangulation.  Refer to */
/*                        Subroutine TRMESH. */

/*       NROW = Number of rows (entries per triangle) re- */
/*              served for the triangle list LTRI.  The value */
/*              must be 6 if only the vertex indexes and */
/*              neighboring triangle indexes are to be */
/*              stored, or 9 if arc indexes are also to be */
/*              assigned and stored.  Refer to LTRI. */

/* The above parameters are not altered by this routine. */

/*       LTRI = Integer array of length at least NROW*NT, */
/*              where NT is at most 2N-4.  (A sufficient */
/*              length is 12N if NROW=6 or 18N if NROW=9.) */

/* On output: */

/*       NT = Number of triangles in the triangulation unless */
/*            IER .NE. 0, in which case NT = 0.  NT = 2N-NB-2 */
/*            if NB .GE. 3 or 2N-4 if NB = 0, where NB is the */
/*            number of boundary nodes. */

/*       LTRI = NROW by NT array whose J-th column contains */
/*              the vertex nodal indexes (first three rows), */
/*              neighboring triangle indexes (second three */
/*              rows), and, if NROW = 9, arc indexes (last */
/*              three rows) associated with triangle J for */
/*              J = 1,...,NT.  The vertices are ordered */
/*              counterclockwise with the first vertex taken */
/*              to be the one with smallest index.  Thus, */
/*              LTRI(2,J) and LTRI(3,J) are larger than */
/*              LTRI(1,J) and index adjacent neighbors of */
/*              node LTRI(1,J).  For I = 1,2,3, LTRI(I+3,J) */
/*              and LTRI(I+6,J) index the triangle and arc, */
/*              respectively, which are opposite (not shared */
/*              by) node LTRI(I,J), with LTRI(I+3,J) = 0 if */
/*              LTRI(I+6,J) indexes a boundary arc.  Vertex */
/*              indexes range from 1 to N, triangle indexes */
/*              from 0 to NT, and, if included, arc indexes */
/*              from 1 to NA, where NA = 3N-NB-3 if NB .GE. 3 */
/*              or 3N-6 if NB = 0.  The triangles are or- */
/*              dered on first (smallest) vertex indexes. */

/*       IER = Error indicator. */
/*             IER = 0 if no errors were encountered. */
/*             IER = 1 if N or NROW is outside its valid */
/*                     range on input. */
/*             IER = 2 if the triangulation data structure */
/*                     (LIST,LPTR,LEND) is invalid.  Note, */
/*                     however, that these arrays are not */
/*                     completely tested for validity. */

/* Modules required by TRLIST:  None */

/* Intrinsic function called by TRLIST:  ABS */

/* *********************************************************** */


/* Local parameters: */

/* ARCS =     Logical variable with value TRUE iff are */
/*              indexes are to be stored */
/* I,J =      LTRI row indexes (1 to 3) associated with */
/*              triangles KT and KN, respectively */
/* I1,I2,I3 = Nodal indexes of triangle KN */
/* ISV =      Variable used to permute indexes I1,I2,I3 */
/* KA =       Arc index and number of currently stored arcs */
/* KN =       Index of the triangle that shares arc I1-I2 */
/*              with KT */
/* KT =       Triangle index and number of currently stored */
/*              triangles */
/* LP =       LIST pointer */
/* LP2 =      Pointer to N2 as a neighbor of N1 */
/* LPL =      Pointer to the last neighbor of I1 */
/* LPLN1 =    Pointer to the last neighbor of N1 */
/* N1,N2,N3 = Nodal indexes of triangle KT */
/* NM2 =      N-2 */


/* Test for invalid input parameters. */

    /* Parameter adjustments */
    --lend;
    --list;
    --lptr;
    ltri_dim1 = *nrow;
    ltri_offset = 1 + ltri_dim1 * 1;
    ltri -= ltri_offset;

    /* Function Body */
    if (*n < 3 || *nrow != 6 && *nrow != 9) {
	goto L11;
    }

/* Initialize parameters for loop on triangles KT = (N1,N2, */
/*   N3), where N1 < N2 and N1 < N3. */

/*   ARCS = TRUE iff arc indexes are to be stored. */
/*   KA,KT = Numbers of currently stored arcs and triangles. */
/*   NM2 = Upper bound on candidates for N1. */

    arcs = *nrow == 9;
    ka = 0;
    kt = 0;
    nm2 = *n - 2;

/* Loop on nodes N1. */

    i__1 = nm2;
    for (n1 = 1; n1 <= i__1; ++n1) {

/* Loop on pairs of adjacent neighbors (N2,N3).  LPLN1 points */
/*   to the last neighbor of N1, and LP2 points to N2. */

	lpln1 = lend[n1];
	lp2 = lpln1;
L1:
	lp2 = lptr[lp2];
	n2 = list[lp2];
	lp = lptr[lp2];
	n3 = (i__2 = list[lp], abs(i__2));
	if (n2 < n1 || n3 < n1) {
	    goto L8;
	}

/* Add a new triangle KT = (N1,N2,N3). */

	++kt;
	ltri[kt * ltri_dim1 + 1] = n1;
	ltri[kt * ltri_dim1 + 2] = n2;
	ltri[kt * ltri_dim1 + 3] = n3;

/* Loop on triangle sides (I2,I1) with neighboring triangles */
/*   KN = (I1,I2,I3). */

	for (i__ = 1; i__ <= 3; ++i__) {
	    if (i__ == 1) {
		i1 = n3;
		i2 = n2;
	    } else if (i__ == 2) {
		i1 = n1;
		i2 = n3;
	    } else {
		i1 = n2;
		i2 = n1;
	    }

/* Set I3 to the neighbor of I1 that follows I2 unless */
/*   I2->I1 is a boundary arc. */

	    lpl = lend[i1];
	    lp = lptr[lpl];
L2:
	    if (list[lp] == i2) {
		goto L3;
	    }
	    lp = lptr[lp];
	    if (lp != lpl) {
		goto L2;
	    }

/*   I2 is the last neighbor of I1 unless the data structure */
/*     is invalid.  Bypass the search for a neighboring */
/*     triangle if I2->I1 is a boundary arc. */

	    if ((i__2 = list[lp], abs(i__2)) != i2) {
		goto L12;
	    }
	    kn = 0;
	    if (list[lp] < 0) {
		goto L6;
	    }

/*   I2->I1 is not a boundary arc, and LP points to I2 as */
/*     a neighbor of I1. */

L3:
	    lp = lptr[lp];
	    i3 = (i__2 = list[lp], abs(i__2));

/* Find J such that LTRI(J,KN) = I3 (not used if KN > KT), */
/*   and permute the vertex indexes of KN so that I1 is */
/*   smallest. */

	    if (i1 < i2 && i1 < i3) {
		j = 3;
	    } else if (i2 < i3) {
		j = 2;
		isv = i1;
		i1 = i2;
		i2 = i3;
		i3 = isv;
	    } else {
		j = 1;
		isv = i1;
		i1 = i3;
		i3 = i2;
		i2 = isv;
	    }

/* Test for KN > KT (triangle index not yet assigned). */

	    if (i1 > n1) {
		goto L7;
	    }

/* Find KN, if it exists, by searching the triangle list in */
/*   reverse order. */

	    for (kn = kt - 1; kn >= 1; --kn) {
		if (ltri[kn * ltri_dim1 + 1] == i1 && ltri[kn * ltri_dim1 + 2]
			 == i2 && ltri[kn * ltri_dim1 + 3] == i3) {
		    goto L5;
		}
/* L4: */
	    }
	    goto L7;

/* Store KT as a neighbor of KN. */

L5:
	    ltri[j + 3 + kn * ltri_dim1] = kt;

/* Store KN as a neighbor of KT, and add a new arc KA. */

L6:
	    ltri[i__ + 3 + kt * ltri_dim1] = kn;
	    if (arcs) {
		++ka;
		ltri[i__ + 6 + kt * ltri_dim1] = ka;
		if (kn != 0) {
		    ltri[j + 6 + kn * ltri_dim1] = ka;
		}
	    }
L7:
	    ;
	}

/* Bottom of loop on triangles. */

L8:
	if (lp2 != lpln1) {
	    goto L1;
	}
/* L9: */
    }

/* No errors encountered. */

    *nt = kt;
    *ier = 0;
    return 0;

/* Invalid input parameter. */

L11:
    *nt = 0;
    *ier = 1;
    return 0;

/* Invalid triangulation data structure:  I1 is a neighbor of */
/*   I2, but I2 is not a neighbor of I1. */

L12:
    *nt = 0;
    *ier = 2;
    return 0;
} /* trlist_ */

/* Subroutine */ int trlprt_(n, x, y, z__, iflag, nrow, nt, ltri, lout)
integer *n;
doublereal *x, *y, *z__;
integer *iflag, *nrow, *nt, *ltri, *lout;
{
    /* Initialized data */

    static integer nmax = 9999;
    static integer nlmax = 58;

    /* Format strings */
    static char fmt_100[] = "(///18x,\002STRIPACK (TRLIST) Output,  N = \002\
,i4)";
    static char fmt_110[] = "(//1x,10x,\002*** Invalid Parameter:  N =\002,i\
5,\002, NROW =\002,i5,\002, NT =\002,i5,\002 ***\002)";
    static char fmt_101[] = "(//8x,\002Node\002,10x,\002X(Node)\002,10x,\002\
Y(Node)\002,10x,\002Z(Node)\002//)";
    static char fmt_108[] = "(///)";
    static char fmt_103[] = "(8x,i4,3d17.6)";
    static char fmt_102[] = "(//16x,\002Node\002,8x,\002Longitude\002,9x,\
\002Latitude\002//)";
    static char fmt_104[] = "(16x,i4,2d17.6)";
    static char fmt_105[] = "(//1x,\002Triangle\002,8x,\002Vertices\002,12x\
,\002Neighbors\002/4x,\002KT\002,7x,\002N1\002,5x,\002N2\002,5x,\002N3\002,4\
x,\002KT1\002,4x,\002KT2\002,4x,\002KT3\002/)";
    static char fmt_106[] = "(//1x,\002Triangle\002,8x,\002Vertices\002,12x\
,\002Neighbors\002,14x,\002Arcs\002/4x,\002KT\002,7x,\002N1\002,5x,\002N2\
\002,5x,\002N3\002,4x,\002KT1\002,4x,\002KT2\002,4x,\002KT3\002,4x,\002KA\
1\002,4x,\002KA2\002,4x,\002KA3\002/)";
    static char fmt_107[] = "(2x,i4,2x,6(3x,i4),3(2x,i5))";
    static char fmt_109[] = "(/1x,\002NB = \002,i4,\002 Boundary Nodes\002,5\
x,\002NA = \002,i5,\002 Arcs\002,5x,\002NT = \002,i5,\002 Triangles\002)";

    /* System generated locals */
    integer ltri_dim1, ltri_offset, i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    static integer i__, k, na, nb, nl, lun;

    /* Fortran I/O blocks */
    static cilist io___437 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___439 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___440 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___442 = { 0, 0, 0, fmt_108, 0 };
    static cilist io___443 = { 0, 0, 0, fmt_103, 0 };
    static cilist io___444 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___445 = { 0, 0, 0, fmt_108, 0 };
    static cilist io___446 = { 0, 0, 0, fmt_104, 0 };
    static cilist io___447 = { 0, 0, 0, fmt_108, 0 };
    static cilist io___448 = { 0, 0, 0, fmt_105, 0 };
    static cilist io___449 = { 0, 0, 0, fmt_106, 0 };
    static cilist io___451 = { 0, 0, 0, fmt_108, 0 };
    static cilist io___452 = { 0, 0, 0, fmt_107, 0 };
    static cilist io___455 = { 0, 0, 0, fmt_109, 0 };



/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/02/98 */

/*   This subroutine prints the triangle list created by Sub- */
/* routine TRLIST and, optionally, the nodal coordinates */
/* (either latitude and longitude or Cartesian coordinates) */
/* on logical unit LOUT.  The numbers of boundary nodes, */
/* triangles, and arcs are also printed. */


/* On input: */

/*       N = Number of nodes in the triangulation. */
/*           3 .LE. N .LE. 9999. */

/*       X,Y,Z = Arrays of length N containing the Cartesian */
/*               coordinates of the nodes if IFLAG = 0, or */
/*               (X and Y only) arrays of length N containing */
/*               longitude and latitude, respectively, if */
/*               IFLAG > 0, or unused dummy parameters if */
/*               IFLAG < 0. */

/*       IFLAG = Nodal coordinate option indicator: */
/*               IFLAG = 0 if X, Y, and Z (assumed to contain */
/*                         Cartesian coordinates) are to be */
/*                         printed (to 6 decimal places). */
/*               IFLAG > 0 if only X and Y (assumed to con- */
/*                         tain longitude and latitude) are */
/*                         to be printed (to 6 decimal */
/*                         places). */
/*               IFLAG < 0 if only the adjacency lists are to */
/*                         be printed. */

/*       NROW = Number of rows (entries per triangle) re- */
/*              served for the triangle list LTRI.  The value */
/*              must be 6 if only the vertex indexes and */
/*              neighboring triangle indexes are stored, or 9 */
/*              if arc indexes are also stored. */

/*       NT = Number of triangles in the triangulation. */
/*            1 .LE. NT .LE. 9999. */

/*       LTRI = NROW by NT array whose J-th column contains */
/*              the vertex nodal indexes (first three rows), */
/*              neighboring triangle indexes (second three */
/*              rows), and, if NROW = 9, arc indexes (last */
/*              three rows) associated with triangle J for */
/*              J = 1,...,NT. */

/*       LOUT = Logical unit number for output.  If LOUT is */
/*              not in the range 0 to 99, output is written */
/*              to unit 6. */

/* Input parameters are not altered by this routine. */

/* On output: */

/*   The triangle list and nodal coordinates (as specified by */
/* IFLAG) are written to unit LOUT. */

/* Modules required by TRLPRT:  None */

/* *********************************************************** */

    /* Parameter adjustments */
    --z__;
    --y;
    --x;
    ltri_dim1 = *nrow;
    ltri_offset = 1 + ltri_dim1 * 1;
    ltri -= ltri_offset;

    /* Function Body */

/* Local parameters: */

/* I =     DO-loop, nodal index, and row index for LTRI */
/* K =     DO-loop and triangle index */
/* LUN =   Logical unit number for output */
/* NA =    Number of triangulation arcs */
/* NB =    Number of boundary nodes */
/* NL =    Number of lines printed on the current page */
/* NLMAX = Maximum number of print lines per page (except */
/*           for the last page which may have two addi- */
/*           tional lines) */
/* NMAX =  Maximum value of N and NT (4-digit format) */

    lun = *lout;
    if (lun < 0 || lun > 99) {
	lun = 6;
    }

/* Print a heading and test for invalid input. */

    io___437.ciunit = lun;
    s_wsfe(&io___437);
    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
    e_wsfe();
    nl = 3;
    if (*n < 3 || *n > nmax || *nrow != 6 && *nrow != 9 || *nt < 1 || *nt > 
	    nmax) {

/* Print an error message and exit. */

	io___439.ciunit = lun;
	s_wsfe(&io___439);
	do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*nrow), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*nt), (ftnlen)sizeof(integer));
	e_wsfe();
	return 0;
    }
    if (*iflag == 0) {

/* Print X, Y, and Z. */

	io___440.ciunit = lun;
	s_wsfe(&io___440);
	e_wsfe();
	nl = 6;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (nl >= nlmax) {
		io___442.ciunit = lun;
		s_wsfe(&io___442);
		e_wsfe();
		nl = 0;
	    }
	    io___443.ciunit = lun;
	    s_wsfe(&io___443);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&y[i__], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&z__[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    ++nl;
/* L1: */
	}
    } else if (*iflag > 0) {

/* Print X (longitude) and Y (latitude). */

	io___444.ciunit = lun;
	s_wsfe(&io___444);
	e_wsfe();
	nl = 6;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (nl >= nlmax) {
		io___445.ciunit = lun;
		s_wsfe(&io___445);
		e_wsfe();
		nl = 0;
	    }
	    io___446.ciunit = lun;
	    s_wsfe(&io___446);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&y[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    ++nl;
/* L2: */
	}
    }

/* Print the triangulation LTRI. */

    if (nl > nlmax / 2) {
	io___447.ciunit = lun;
	s_wsfe(&io___447);
	e_wsfe();
	nl = 0;
    }
    if (*nrow == 6) {
	io___448.ciunit = lun;
	s_wsfe(&io___448);
	e_wsfe();
    } else {
	io___449.ciunit = lun;
	s_wsfe(&io___449);
	e_wsfe();
    }
    nl += 5;
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	if (nl >= nlmax) {
	    io___451.ciunit = lun;
	    s_wsfe(&io___451);
	    e_wsfe();
	    nl = 0;
	}
	io___452.ciunit = lun;
	s_wsfe(&io___452);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	i__2 = *nrow;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    do_fio(&c__1, (char *)&ltri[i__ + k * ltri_dim1], (ftnlen)sizeof(
		    integer));
	}
	e_wsfe();
	++nl;
/* L3: */
    }

/* Print NB, NA, and NT (boundary nodes, arcs, and */
/*   triangles). */

    nb = (*n << 1) - *nt - 2;
    if (nb < 3) {
	nb = 0;
	na = *n * 3 - 6;
    } else {
	na = *nt + *n - 1;
    }
    io___455.ciunit = lun;
    s_wsfe(&io___455);
    do_fio(&c__1, (char *)&nb, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&na, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*nt), (ftnlen)sizeof(integer));
    e_wsfe();
    return 0;

/* Print formats: */

} /* trlprt_ */

/* Subroutine */ int trmesh_(n, x, y, z__, list, lptr, lend, lnew, near__, 
	next, dist, ier)
integer *n;
doublereal *x, *y, *z__;
integer *list, *lptr, *lend, *lnew, *near__, *next;
doublereal *dist;
integer *ier;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    extern logical left_();
    static doublereal d__;
    static integer i__, j, k;
    static doublereal d1, d2;
    static integer nexti;
    static doublereal d3;
    static integer i0;
    extern /* Subroutine */ int addnod_();
    static integer lp, nn, lpl;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   03/04/03 */

/*   This subroutine creates a Delaunay triangulation of a */
/* set of N arbitrarily distributed points, referred to as */
/* nodes, on the surface of the unit sphere.  The Delaunay */
/* triangulation is defined as a set of (spherical) triangles */
/* with the following five properties: */

/*  1)  The triangle vertices are nodes. */
/*  2)  No triangle contains a node other than its vertices. */
/*  3)  The interiors of the triangles are pairwise disjoint. */
/*  4)  The union of triangles is the convex hull of the set */
/*        of nodes (the smallest convex set that contains */
/*        the nodes).  If the nodes are not contained in a */
/*        single hemisphere, their convex hull is the en- */
/*        tire sphere and there are no boundary nodes. */
/*        Otherwise, there are at least three boundary nodes. */
/*  5)  The interior of the circumcircle of each triangle */
/*        contains no node. */

/* The first four properties define a triangulation, and the */
/* last property results in a triangulation which is as close */
/* as possible to equiangular in a certain sense and which is */
/* uniquely defined unless four or more nodes lie in a common */
/* plane.  This property makes the triangulation well-suited */
/* for solving closest-point problems and for triangle-based */
/* interpolation. */

/*   The algorithm has expected time complexity O(N*log(N)) */
/* for most nodal distributions. */

/*   Spherical coordinates (latitude and longitude) may be */
/* converted to Cartesian coordinates by Subroutine TRANS. */

/*   The following is a list of the software package modules */
/* which a user may wish to call directly: */

/*  ADDNOD - Updates the triangulation by appending a new */
/*             node. */

/*  AREAS  - Returns the area of a spherical triangle. */

/*  AREAV  - Returns the area of a Voronoi region associated */
/*           with an interior node without requiring that the */
/*           entire Voronoi diagram be computed and stored. */

/*  BNODES - Returns an array containing the indexes of the */
/*             boundary nodes (if any) in counterclockwise */
/*             order.  Counts of boundary nodes, triangles, */
/*             and arcs are also returned. */

/*  CIRCLE - Computes the coordinates of a sequence of uni- */
/*           formly spaced points on the unit circle centered */
/*           at (0,0). */

/*  CIRCUM - Returns the circumcenter of a spherical trian- */
/*             gle. */

/*  CRLIST - Returns the set of triangle circumcenters */
/*             (Voronoi vertices) and circumradii associated */
/*             with a triangulation. */

/*  DELARC - Deletes a boundary arc from a triangulation. */

/*  DELNOD - Updates the triangulation with a nodal deletion. */

/*  EDGE   - Forces an arbitrary pair of nodes to be connec- */
/*             ted by an arc in the triangulation. */

/*  GETNP  - Determines the ordered sequence of L closest */
/*             nodes to a given node, along with the associ- */
/*             ated distances. */

/*  INSIDE - Locates a point relative to a polygon on the */
/*             surface of the sphere. */

/*  INTRSC - Returns the point of intersection between a */
/*             pair of great circle arcs. */

/*  JRAND  - Generates a uniformly distributed pseudo-random */
/*             integer. */

/*  LEFT   - Locates a point relative to a great circle. */

/*  NEARND - Returns the index of the nearest node to an */
/*             arbitrary point, along with its squared */
/*             distance. */

/*  PROJCT - Applies a perspective-depth projection to a */
/*             point in 3-space. */

/*  SCOORD - Converts a point from Cartesian coordinates to */
/*             spherical coordinates. */

/*  STORE  - Forces a value to be stored in main memory so */
/*             that the precision of floating point numbers */
/*             in memory locations rather than registers is */
/*             computed. */

/*  TRANS  - Transforms spherical coordinates into Cartesian */
/*             coordinates on the unit sphere for input to */
/*             Subroutine TRMESH. */

/*  TRLIST - Converts the triangulation data structure to a */
/*             triangle list more suitable for use in a fin- */
/*             ite element code. */

/*  TRLPRT - Prints the triangle list created by Subroutine */
/*             TRLIST. */

/*  TRMESH - Creates a Delaunay triangulation of a set of */
/*             nodes. */

/*  TRPLOT - Creates a level-2 Encapsulated Postscript (EPS) */
/*             file containing a triangulation plot. */

/*  TRPRNT - Prints the triangulation data structure and, */
/*             optionally, the nodal coordinates. */

/*  VRPLOT - Creates a level-2 Encapsulated Postscript (EPS) */
/*             file containing a Voronoi diagram plot. */


/* On input: */

/*       N = Number of nodes in the triangulation.  N .GE. 3. */

/*       X,Y,Z = Arrays of length N containing the Cartesian */
/*               coordinates of distinct nodes.  (X(K),Y(K), */
/*               Z(K)) is referred to as node K, and K is re- */
/*               ferred to as a nodal index.  It is required */
/*               that X(K)**2 + Y(K)**2 + Z(K)**2 = 1 for all */
/*               K.  The first three nodes must not be col- */
/*               linear (lie on a common great circle). */

/* The above parameters are not altered by this routine. */

/*       LIST,LPTR = Arrays of length at least 6N-12. */

/*       LEND = Array of length at least N. */

/*       NEAR,NEXT,DIST = Work space arrays of length at */
/*                        least N.  The space is used to */
/*                        efficiently determine the nearest */
/*                        triangulation node to each un- */
/*                        processed node for use by ADDNOD. */

/* On output: */

/*       LIST = Set of nodal indexes which, along with LPTR, */
/*              LEND, and LNEW, define the triangulation as a */
/*              set of N adjacency lists -- counterclockwise- */
/*              ordered sequences of neighboring nodes such */
/*              that the first and last neighbors of a bound- */
/*              ary node are boundary nodes (the first neigh- */
/*              bor of an interior node is arbitrary).  In */
/*              order to distinguish between interior and */
/*              boundary nodes, the last neighbor of each */
/*              boundary node is represented by the negative */
/*              of its index. */

/*       LPTR = Set of pointers (LIST indexes) in one-to-one */
/*              correspondence with the elements of LIST. */
/*              LIST(LPTR(I)) indexes the node which follows */
/*              LIST(I) in cyclical counterclockwise order */
/*              (the first neighbor follows the last neigh- */
/*              bor). */

/*       LEND = Set of pointers to adjacency lists.  LEND(K) */
/*              points to the last neighbor of node K for */
/*              K = 1,...,N.  Thus, LIST(LEND(K)) < 0 if and */
/*              only if K is a boundary node. */

/*       LNEW = Pointer to the first empty location in LIST */
/*              and LPTR (list length plus one).  LIST, LPTR, */
/*              LEND, and LNEW are not altered if IER < 0, */
/*              and are incomplete if IER > 0. */

/*       NEAR,NEXT,DIST = Garbage. */

/*       IER = Error indicator: */
/*             IER =  0 if no errors were encountered. */
/*             IER = -1 if N < 3 on input. */
/*             IER = -2 if the first three nodes are */
/*                      collinear. */
/*             IER =  L if nodes L and M coincide for some */
/*                      M > L.  The data structure represents */
/*                      a triangulation of nodes 1 to M-1 in */
/*                      this case. */

/* Modules required by TRMESH:  ADDNOD, BDYADD, COVSPH, */
/*                                INSERT, INTADD, JRAND, */
/*                                LEFT, LSTPTR, STORE, SWAP, */
/*                                SWPTST, TRFIND */

/* Intrinsic function called by TRMESH:  ABS */

/* *********************************************************** */


/* Local parameters: */

/* D =        (Negative cosine of) distance from node K to */
/*              node I */
/* D1,D2,D3 = Distances from node K to nodes 1, 2, and 3, */
/*              respectively */
/* I,J =      Nodal indexes */
/* I0 =       Index of the node preceding I in a sequence of */
/*              unprocessed nodes:  I = NEXT(I0) */
/* K =        Index of node to be added and DO-loop index: */
/*              K > 3 */
/* LP =       LIST index (pointer) of a neighbor of K */
/* LPL =      Pointer to the last neighbor of K */
/* NEXTI =    NEXT(I) */
/* NN =       Local copy of N */

    /* Parameter adjustments */
    --dist;
    --next;
    --near__;
    --lend;
    --z__;
    --y;
    --x;
    --list;
    --lptr;

    /* Function Body */
    nn = *n;
    if (nn < 3) {
	*ier = -1;
	return 0;
    }

/* Store the first triangle in the linked list. */

    if (! left_(&x[1], &y[1], &z__[1], &x[2], &y[2], &z__[2], &x[3], &y[3], &
	    z__[3])) {

/*   The first triangle is (3,2,1) = (2,1,3) = (1,3,2). */

	list[1] = 3;
	lptr[1] = 2;
	list[2] = -2;
	lptr[2] = 1;
	lend[1] = 2;

	list[3] = 1;
	lptr[3] = 4;
	list[4] = -3;
	lptr[4] = 3;
	lend[2] = 4;

	list[5] = 2;
	lptr[5] = 6;
	list[6] = -1;
	lptr[6] = 5;
	lend[3] = 6;

    } else if (! left_(&x[2], &y[2], &z__[2], &x[1], &y[1], &z__[1], &x[3], &
	    y[3], &z__[3])) {

/*   The first triangle is (1,2,3):  3 Strictly Left 1->2, */
/*     i.e., node 3 lies in the left hemisphere defined by */
/*     arc 1->2. */

	list[1] = 2;
	lptr[1] = 2;
	list[2] = -3;
	lptr[2] = 1;
	lend[1] = 2;

	list[3] = 3;
	lptr[3] = 4;
	list[4] = -1;
	lptr[4] = 3;
	lend[2] = 4;

	list[5] = 1;
	lptr[5] = 6;
	list[6] = -2;
	lptr[6] = 5;
	lend[3] = 6;

    } else {

/*   The first three nodes are collinear. */

	*ier = -2;
	return 0;
    }

/* Initialize LNEW and test for N = 3. */

    *lnew = 7;
    if (nn == 3) {
	*ier = 0;
	return 0;
    }

/* A nearest-node data structure (NEAR, NEXT, and DIST) is */
/*   used to obtain an expected-time (N*log(N)) incremental */
/*   algorithm by enabling constant search time for locating */
/*   each new node in the triangulation. */

/* For each unprocessed node K, NEAR(K) is the index of the */
/*   triangulation node closest to K (used as the starting */
/*   point for the search in Subroutine TRFIND) and DIST(K) */
/*   is an increasing function of the arc length (angular */
/*   distance) between nodes K and NEAR(K):  -Cos(a) for arc */
/*   length a. */

/* Since it is necessary to efficiently find the subset of */
/*   unprocessed nodes associated with each triangulation */
/*   node J (those that have J as their NEAR entries), the */
/*   subsets are stored in NEAR and NEXT as follows:  for */
/*   each node J in the triangulation, I = NEAR(J) is the */
/*   first unprocessed node in J's set (with I = 0 if the */
/*   set is empty), L = NEXT(I) (if I > 0) is the second, */
/*   NEXT(L) (if L > 0) is the third, etc.  The nodes in each */
/*   set are initially ordered by increasing indexes (which */
/*   maximizes efficiency) but that ordering is not main- */
/*   tained as the data structure is updated. */

/* Initialize the data structure for the single triangle. */

    near__[1] = 0;
    near__[2] = 0;
    near__[3] = 0;
    for (k = nn; k >= 4; --k) {
	d1 = -(x[k] * x[1] + y[k] * y[1] + z__[k] * z__[1]);
	d2 = -(x[k] * x[2] + y[k] * y[2] + z__[k] * z__[2]);
	d3 = -(x[k] * x[3] + y[k] * y[3] + z__[k] * z__[3]);
	if (d1 <= d2 && d1 <= d3) {
	    near__[k] = 1;
	    dist[k] = d1;
	    next[k] = near__[1];
	    near__[1] = k;
	} else if (d2 <= d1 && d2 <= d3) {
	    near__[k] = 2;
	    dist[k] = d2;
	    next[k] = near__[2];
	    near__[2] = k;
	} else {
	    near__[k] = 3;
	    dist[k] = d3;
	    next[k] = near__[3];
	    near__[3] = k;
	}
/* L1: */
    }

/* Add the remaining nodes */

    i__1 = nn;
    for (k = 4; k <= i__1; ++k) {
	addnod_(&near__[k], &k, &x[1], &y[1], &z__[1], &list[1], &lptr[1], &
		lend[1], lnew, ier);
	if (*ier != 0) {
	    return 0;
	}

/* Remove K from the set of unprocessed nodes associated */
/*   with NEAR(K). */

	i__ = near__[k];
	if (near__[i__] == k) {
	    near__[i__] = next[k];
	} else {
	    i__ = near__[i__];
L2:
	    i0 = i__;
	    i__ = next[i0];
	    if (i__ != k) {
		goto L2;
	    }
	    next[i0] = next[k];
	}
	near__[k] = 0;

/* Loop on neighbors J of node K. */

	lpl = lend[k];
	lp = lpl;
L3:
	lp = lptr[lp];
	j = (i__2 = list[lp], abs(i__2));

/* Loop on elements I in the sequence of unprocessed nodes */
/*   associated with J:  K is a candidate for replacing J */
/*   as the nearest triangulation node to I.  The next value */
/*   of I in the sequence, NEXT(I), must be saved before I */
/*   is moved because it is altered by adding I to K's set. */

	i__ = near__[j];
L4:
	if (i__ == 0) {
	    goto L5;
	}
	nexti = next[i__];

/* Test for the distance from I to K less than the distance */
/*   from I to J. */

	d__ = -(x[i__] * x[k] + y[i__] * y[k] + z__[i__] * z__[k]);
	if (d__ < dist[i__]) {

/* Replace J by K as the nearest triangulation node to I: */
/*   update NEAR(I) and DIST(I), and remove I from J's set */
/*   of unprocessed nodes and add it to K's set. */

	    near__[i__] = k;
	    dist[i__] = d__;
	    if (i__ == near__[j]) {
		near__[j] = nexti;
	    } else {
		next[i0] = nexti;
	    }
	    next[i__] = near__[k];
	    near__[k] = i__;
	} else {
	    i0 = i__;
	}

/* Bottom of loop on I. */

	i__ = nexti;
	goto L4;

/* Bottom of loop on neighbors J. */

L5:
	if (lp != lpl) {
	    goto L3;
	}
/* L6: */
    }
    return 0;
} /* trmesh_ */

/* Subroutine */ int trplot_(lun, pltsiz, elat, elon, a, n, x, y, z__, list, 
	lptr, lend, title, numbr, ier, title_len)
integer *lun;
doublereal *pltsiz, *elat, *elon, *a;
integer *n;
doublereal *x, *y, *z__;
integer *list, *lptr, *lend;
char *title;
logical *numbr;
integer *ier;
ftnlen title_len;
{
    /* Initialized data */

    static logical annot = TRUE_;
    static doublereal fsizn = 10.;
    static doublereal fsizt = 16.;
    static doublereal tol = .5;

    /* Format strings */
    static char fmt_100[] = "(\002%!PS-Adobe-3.0 EPSF-3.0\002/\002%%Bounding\
Box:\002,4i4/\002%%Title:  Triangulation\002/\002%%Creator:  STRIPACK\002\
/\002%%EndComments\002)";
    static char fmt_110[] = "(f12.6,\002 setlinewidth\002)";
    static char fmt_120[] = "(\002306 396 \002,i3,\002 0 360 arc\002)";
    static char fmt_130[] = "(\002stroke\002)";
    static char fmt_140[] = "(2f12.6,\002 translate\002/2f12.6,\002 scale\
\002)";
    static char fmt_150[] = "(\002gsave\002)";
    static char fmt_160[] = "(\0020 0 \002,f12.6,\002 0 360 arc\002)";
    static char fmt_170[] = "(\002clip newpath\002)";
    static char fmt_190[] = "(\002grestore\002)";
    static char fmt_200[] = "(\002/Helvetica findfont\002/f12.6,\002 scalefo\
nt setfont\002)";
    static char fmt_210[] = "(2f12.6,\002 moveto\002)";
    static char fmt_220[] = "(\002(\002,i3,\002) show\002)";
    static char fmt_230[] = "(a80/\002  stringwidth pop 2 div neg \002,f12\
.6,\002 moveto\002)";
    static char fmt_240[] = "(a80/\002  show\002)";
    static char fmt_250[] = "(\002(Window center:  ELAT = \002,f7.2,\002,  E\
LON = \002,f8.2,\002) show\002)";
    static char fmt_260[] = "(\002(Angular extent:  A = \002,f5.2,\002) sho\
w\002)";
    static char fmt_270[] = "(\002stroke\002/\002showpage\002/\002%%EOF\002)";
    static char fmt_280[] = "(a1)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double atan(), sin();
    integer i_dnnt(), s_wsfe(), do_fio(), e_wsfe();
    double cos(), sqrt();

    /* Local variables */
    static integer nseg;
    static doublereal t;
    static integer n0, n1;
    static doublereal p0[3], p1[3], cf, r11, r12, r21, ct, r22, r23, sf;
    static integer ir, lp;
    static doublereal ex, ey, ez, wr, tx, ty;
    extern /* Subroutine */ int drwarc_();
    static integer lpl;
    static doublereal wrs;
    static integer ipx1, ipx2, ipy1, ipy2;

    /* Fortran I/O blocks */
    static cilist io___480 = { 1, 0, 0, fmt_100, 0 };
    static cilist io___482 = { 1, 0, 0, fmt_110, 0 };
    static cilist io___483 = { 1, 0, 0, fmt_120, 0 };
    static cilist io___484 = { 1, 0, 0, fmt_130, 0 };
    static cilist io___488 = { 1, 0, 0, fmt_140, 0 };
    static cilist io___489 = { 1, 0, 0, fmt_110, 0 };
    static cilist io___490 = { 1, 0, 0, fmt_150, 0 };
    static cilist io___491 = { 1, 0, 0, fmt_160, 0 };
    static cilist io___492 = { 1, 0, 0, fmt_170, 0 };
    static cilist io___509 = { 1, 0, 0, fmt_130, 0 };
    static cilist io___510 = { 1, 0, 0, fmt_190, 0 };
    static cilist io___511 = { 1, 0, 0, fmt_200, 0 };
    static cilist io___512 = { 1, 0, 0, fmt_210, 0 };
    static cilist io___513 = { 1, 0, 0, fmt_220, 0 };
    static cilist io___514 = { 1, 0, 0, fmt_200, 0 };
    static cilist io___515 = { 1, 0, 0, fmt_230, 0 };
    static cilist io___516 = { 1, 0, 0, fmt_240, 0 };
    static cilist io___517 = { 1, 0, 0, fmt_210, 0 };
    static cilist io___518 = { 1, 0, 0, fmt_250, 0 };
    static cilist io___519 = { 1, 0, 0, fmt_210, 0 };
    static cilist io___520 = { 1, 0, 0, fmt_260, 0 };
    static cilist io___521 = { 1, 0, 0, fmt_270, 0 };
    static cilist io___522 = { 1, 0, 0, fmt_280, 0 };



/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   03/04/03 */

/*   This subroutine creates a level-2 Encapsulated Post- */
/* script (EPS) file containing a graphical display of a */
/* triangulation of a set of nodes on the surface of the unit */
/* sphere.  The visible portion of the triangulation is */
/* projected onto the plane that contains the origin and has */
/* normal defined by a user-specified eye-position. */


/* On input: */

/*       LUN = Logical unit number in the range 0 to 99. */
/*             The unit should be opened with an appropriate */
/*             file name before the call to this routine. */

/*       PLTSIZ = Plot size in inches.  A circular window in */
/*                the projection plane is mapped to a circu- */
/*                lar viewport with diameter equal to .88* */
/*                PLTSIZ (leaving room for labels outside the */
/*                viewport).  The viewport is centered on the */
/*                8.5 by 11 inch page, and its boundary is */
/*                drawn.  1.0 .LE. PLTSIZ .LE. 8.5. */

/*       ELAT,ELON = Latitude and longitude (in degrees) of */
/*                   the center of projection E (the center */
/*                   of the plot).  The projection plane is */
/*                   the plane that contains the origin and */
/*                   has E as unit normal.  In a rotated */
/*                   coordinate system for which E is the */
/*                   north pole, the projection plane con- */
/*                   tains the equator, and only northern */
/*                   hemisphere nodes are visible (from the */
/*                   point at infinity in the direction E). */
/*                   These are projected orthogonally onto */
/*                   the projection plane (by zeroing the z- */
/*                   component in the rotated coordinate */
/*                   system).  ELAT and ELON must be in the */
/*                   range -90 to 90 and -180 to 180, respec- */
/*                   tively. */

/*       A = Angular distance in degrees from E to the boun- */
/*           dary of a circular window against which the */
/*           triangulation is clipped.  The projected window */
/*           is a disk of radius r = Sin(A) centered at the */
/*           origin, and only visible nodes whose projections */
/*           are within distance r of the origin are included */
/*           in the plot.  Thus, if A = 90, the plot includes */
/*           the entire hemisphere centered at E.  0 .LT. A */
/*           .LE. 90. */

/*       N = Number of nodes in the triangulation.  N .GE. 3. */

/*       X,Y,Z = Arrays of length N containing the Cartesian */
/*               coordinates of the nodes (unit vectors). */

/*       LIST,LPTR,LEND = Data structure defining the trian- */
/*                        gulation.  Refer to Subroutine */
/*                        TRMESH. */

/*       TITLE = Type CHARACTER variable or constant contain- */
/*               ing a string to be centered above the plot. */
/*               The string must be enclosed in parentheses; */
/*               i.e., the first and last characters must be */
/*               '(' and ')', respectively, but these are not */
/*               displayed.  TITLE may have at most 80 char- */
/*               acters including the parentheses. */

/*       NUMBR = Option indicator:  If NUMBR = TRUE, the */
/*               nodal indexes are plotted next to the nodes. */

/* Input parameters are not altered by this routine. */

/* On output: */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered. */
/*             IER = 1 if LUN, PLTSIZ, or N is outside its */
/*                     valid range. */
/*             IER = 2 if ELAT, ELON, or A is outside its */
/*                     valid range. */
/*             IER = 3 if an error was encountered in writing */
/*                     to unit LUN. */

/*   The values in the data statement below may be altered */
/* in order to modify various plotting options. */

/* Module required by TRPLOT:  DRWARC */

/* Intrinsic functions called by TRPLOT:  ABS, ATAN, COS, */
/*                                          DBLE, NINT, SIN, */
/*                                          SQRT */

/* *********************************************************** */


    /* Parameter adjustments */
    --lend;
    --z__;
    --y;
    --x;
    --list;
    --lptr;

    /* Function Body */

/* Local parameters: */

/* ANNOT =     Logical variable with value TRUE iff the plot */
/*               is to be annotated with the values of ELAT, */
/*               ELON, and A */
/* CF =        Conversion factor for degrees to radians */
/* CT =        Cos(ELAT) */
/* EX,EY,EZ =  Cartesian coordinates of the eye-position E */
/* FSIZN =     Font size in points for labeling nodes with */
/*               their indexes if NUMBR = TRUE */
/* FSIZT =     Font size in points for the title (and */
/*               annotation if ANNOT = TRUE) */
/* IPX1,IPY1 = X and y coordinates (in points) of the lower */
/*               left corner of the bounding box or viewport */
/*               box */
/* IPX2,IPY2 = X and y coordinates (in points) of the upper */
/*               right corner of the bounding box or viewport */
/*               box */
/* IR =        Half the width (height) of the bounding box or */
/*               viewport box in points -- viewport radius */
/* LP =        LIST index (pointer) */
/* LPL =       Pointer to the last neighbor of N0 */
/* N0 =        Index of a node whose incident arcs are to be */
/*               drawn */
/* N1 =        Neighbor of N0 */
/* NSEG =      Number of line segments used by DRWARC in a */
/*               polygonal approximation to a projected edge */
/* P0 =        Coordinates of N0 in the rotated coordinate */
/*               system or label location (first two */
/*               components) */
/* P1 =        Coordinates of N1 in the rotated coordinate */
/*               system or intersection of edge N0-N1 with */
/*               the equator (in the rotated coordinate */
/*               system) */
/* R11...R23 = Components of the first two rows of a rotation */
/*               that maps E to the north pole (0,0,1) */
/* SF =        Scale factor for mapping world coordinates */
/*               (window coordinates in [-WR,WR] X [-WR,WR]) */
/*               to viewport coordinates in [IPX1,IPX2] X */
/*               [IPY1,IPY2] */
/* T =         Temporary variable */
/* TOL =       Maximum distance in points between a projected */
/*               triangulation edge and its approximation by */
/*               a polygonal curve */
/* TX,TY =     Translation vector for mapping world coordi- */
/*               nates to viewport coordinates */
/* WR =        Window radius r = Sin(A) */
/* WRS =       WR**2 */


/* Test for invalid parameters. */

    if (*lun < 0 || *lun > 99 || *pltsiz < 1. || *pltsiz > 8.5 || *n < 3) {
	goto L11;
    }
    if (abs(*elat) > 90. || abs(*elon) > 180. || *a > 90.) {
	goto L12;
    }

/* Compute a conversion factor CF for degrees to radians */
/*   and compute the window radius WR. */

    cf = atan(1.) / 45.;
    wr = sin(cf * *a);
    wrs = wr * wr;

/* Compute the lower left (IPX1,IPY1) and upper right */
/*   (IPX2,IPY2) corner coordinates of the bounding box. */
/*   The coordinates, specified in default user space units */
/*   (points, at 72 points/inch with origin at the lower */
/*   left corner of the page), are chosen to preserve the */
/*   square aspect ratio, and to center the plot on the 8.5 */
/*   by 11 inch page.  The center of the page is (306,396), */
/*   and IR = PLTSIZ/2 in points. */

    d__1 = *pltsiz * 36.;
    ir = i_dnnt(&d__1);
    ipx1 = 306 - ir;
    ipx2 = ir + 306;
    ipy1 = 396 - ir;
    ipy2 = ir + 396;

/* Output header comments. */

    io___480.ciunit = *lun;
    i__1 = s_wsfe(&io___480);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&ipx1, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&ipy1, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&ipx2, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&ipy2, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }

/* Set (IPX1,IPY1) and (IPX2,IPY2) to the corner coordinates */
/*   of a viewport box obtained by shrinking the bounding box */
/*   by 12% in each dimension. */

    d__1 = (doublereal) ir * .88;
    ir = i_dnnt(&d__1);
    ipx1 = 306 - ir;
    ipx2 = ir + 306;
    ipy1 = 396 - ir;
    ipy2 = ir + 396;

/* Set the line thickness to 2 points, and draw the */
/*   viewport boundary. */

    t = 2.;
    io___482.ciunit = *lun;
    i__1 = s_wsfe(&io___482);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&t, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }
    io___483.ciunit = *lun;
    i__1 = s_wsfe(&io___483);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&ir, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }
    io___484.ciunit = *lun;
    i__1 = s_wsfe(&io___484);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }

/* Set up an affine mapping from the window box [-WR,WR] X */
/*   [-WR,WR] to the viewport box. */

    sf = (doublereal) ir / wr;
    tx = ipx1 + sf * wr;
    ty = ipy1 + sf * wr;
    io___488.ciunit = *lun;
    i__1 = s_wsfe(&io___488);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&tx, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&ty, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&sf, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&sf, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }

/* The line thickness must be changed to reflect the new */
/*   scaling which is applied to all subsequent output. */
/*   Set it to 1.0 point. */

    t = 1. / sf;
    io___489.ciunit = *lun;
    i__1 = s_wsfe(&io___489);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&t, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }

/* Save the current graphics state, and set the clip path to */
/*   the boundary of the window. */

    io___490.ciunit = *lun;
    i__1 = s_wsfe(&io___490);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }
    io___491.ciunit = *lun;
    i__1 = s_wsfe(&io___491);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&wr, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }
    io___492.ciunit = *lun;
    i__1 = s_wsfe(&io___492);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }

/* Compute the Cartesian coordinates of E and the components */
/*   of a rotation R which maps E to the north pole (0,0,1). */
/*   R is taken to be a rotation about the z-axis (into the */
/*   yz-plane) followed by a rotation about the x-axis chosen */
/*   so that the view-up direction is (0,0,1), or (-1,0,0) if */
/*   E is the north or south pole. */

/*           ( R11  R12  0   ) */
/*       R = ( R21  R22  R23 ) */
/*           ( EX   EY   EZ  ) */

    t = cf * *elon;
    ct = cos(cf * *elat);
    ex = ct * cos(t);
    ey = ct * sin(t);
    ez = sin(cf * *elat);
    if (ct != 0.) {
	r11 = -ey / ct;
	r12 = ex / ct;
    } else {
	r11 = 0.;
	r12 = 1.;
    }
    r21 = -ez * r12;
    r22 = ez * r11;
    r23 = ct;

/* Loop on visible nodes N0 that project to points */
/*   (P0(1),P0(2)) in the window. */

    i__1 = *n;
    for (n0 = 1; n0 <= i__1; ++n0) {
	p0[2] = ex * x[n0] + ey * y[n0] + ez * z__[n0];
	if (p0[2] < 0.) {
	    goto L3;
	}
	p0[0] = r11 * x[n0] + r12 * y[n0];
	p0[1] = r21 * x[n0] + r22 * y[n0] + r23 * z__[n0];
	if (p0[0] * p0[0] + p0[1] * p0[1] > wrs) {
	    goto L3;
	}
	lpl = lend[n0];
	lp = lpl;

/* Loop on neighbors N1 of N0.  LPL points to the last */
/*   neighbor of N0.  Copy the components of N1 into P. */

L1:
	lp = lptr[lp];
	n1 = (i__2 = list[lp], abs(i__2));
	p1[0] = r11 * x[n1] + r12 * y[n1];
	p1[1] = r21 * x[n1] + r22 * y[n1] + r23 * z__[n1];
	p1[2] = ex * x[n1] + ey * y[n1] + ez * z__[n1];
	if (p1[2] < 0.) {

/*   N1 is a 'southern hemisphere' point.  Move it to the */
/*     intersection of edge N0-N1 with the equator so that */
/*     the edge is clipped properly.  P1(3) is set to 0. */

	    p1[0] = p0[2] * p1[0] - p1[2] * p0[0];
	    p1[1] = p0[2] * p1[1] - p1[2] * p0[1];
	    t = sqrt(p1[0] * p1[0] + p1[1] * p1[1]);
	    p1[0] /= t;
	    p1[1] /= t;
	}

/*   If node N1 is in the window and N1 < N0, bypass edge */
/*     N0->N1 (since edge N1->N0 has already been drawn). */

	if (p1[2] >= 0. && p1[0] * p1[0] + p1[1] * p1[1] <= wrs && n1 < n0) {
	    goto L2;
	}

/*   Add the edge to the path.  (TOL is converted to world */
/*     coordinates.) */

	if (p1[2] < 0.) {
	    p1[2] = 0.;
	}
	d__1 = tol / sf;
	drwarc_(lun, p0, p1, &d__1, &nseg);

/* Bottom of loops. */

L2:
	if (lp != lpl) {
	    goto L1;
	}
L3:
	;
    }

/* Paint the path and restore the saved graphics state (with */
/*   no clip path). */

    io___509.ciunit = *lun;
    i__1 = s_wsfe(&io___509);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }
    io___510.ciunit = *lun;
    i__1 = s_wsfe(&io___510);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }
    if (*numbr) {

/* Nodes in the window are to be labeled with their indexes. */
/*   Convert FSIZN from points to world coordinates, and */
/*   output the commands to select a font and scale it. */

	t = fsizn / sf;
	io___511.ciunit = *lun;
	i__1 = s_wsfe(&io___511);
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = do_fio(&c__1, (char *)&t, (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L13;
	}

/* Loop on visible nodes N0 that project to points */
/*   P0 = (P0(1),P0(2)) in the window. */

	i__1 = *n;
	for (n0 = 1; n0 <= i__1; ++n0) {
	    if (ex * x[n0] + ey * y[n0] + ez * z__[n0] < 0.) {
		goto L4;
	    }
	    p0[0] = r11 * x[n0] + r12 * y[n0];
	    p0[1] = r21 * x[n0] + r22 * y[n0] + r23 * z__[n0];
	    if (p0[0] * p0[0] + p0[1] * p0[1] > wrs) {
		goto L4;
	    }

/*   Move to P0 and draw the label N0.  The first character */
/*     will will have its lower left corner about one */
/*     character width to the right of the nodal position. */

	    io___512.ciunit = *lun;
	    i__2 = s_wsfe(&io___512);
	    if (i__2 != 0) {
		goto L13;
	    }
	    i__2 = do_fio(&c__1, (char *)&p0[0], (ftnlen)sizeof(doublereal));
	    if (i__2 != 0) {
		goto L13;
	    }
	    i__2 = do_fio(&c__1, (char *)&p0[1], (ftnlen)sizeof(doublereal));
	    if (i__2 != 0) {
		goto L13;
	    }
	    i__2 = e_wsfe();
	    if (i__2 != 0) {
		goto L13;
	    }
	    io___513.ciunit = *lun;
	    i__2 = s_wsfe(&io___513);
	    if (i__2 != 0) {
		goto L13;
	    }
	    i__2 = do_fio(&c__1, (char *)&n0, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L13;
	    }
	    i__2 = e_wsfe();
	    if (i__2 != 0) {
		goto L13;
	    }
L4:
	    ;
	}
    }

/* Convert FSIZT from points to world coordinates, and output */
/*   the commands to select a font and scale it. */

    t = fsizt / sf;
    io___514.ciunit = *lun;
    i__1 = s_wsfe(&io___514);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&t, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }

/* Display TITLE centered above the plot: */

    p0[1] = wr + t * 3.;
    io___515.ciunit = *lun;
    i__1 = s_wsfe(&io___515);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, title, title_len);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&p0[1], (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }
    io___516.ciunit = *lun;
    i__1 = s_wsfe(&io___516);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, title, title_len);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }
    if (annot) {

/* Display the window center and radius below the plot. */

	p0[0] = -wr;
	p0[1] = -wr - 50. / sf;
	io___517.ciunit = *lun;
	i__1 = s_wsfe(&io___517);
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = do_fio(&c__1, (char *)&p0[0], (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = do_fio(&c__1, (char *)&p0[1], (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L13;
	}
	io___518.ciunit = *lun;
	i__1 = s_wsfe(&io___518);
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = do_fio(&c__1, (char *)&(*elat), (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = do_fio(&c__1, (char *)&(*elon), (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L13;
	}
	p0[1] -= t * 2.;
	io___519.ciunit = *lun;
	i__1 = s_wsfe(&io___519);
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = do_fio(&c__1, (char *)&p0[0], (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = do_fio(&c__1, (char *)&p0[1], (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L13;
	}
	io___520.ciunit = *lun;
	i__1 = s_wsfe(&io___520);
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = do_fio(&c__1, (char *)&(*a), (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L13;
	}
    }

/* Paint the path and output the showpage command and */
/*   end-of-file indicator. */

    io___521.ciunit = *lun;
    i__1 = s_wsfe(&io___521);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }

/* HP's interpreters require a one-byte End-of-PostScript-Job */
/*   indicator (to eliminate a timeout error message): */
/*   ASCII 4. */

    io___522.ciunit = *lun;
    i__1 = s_wsfe(&io___522);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, "\004", (ftnlen)1);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }

/* No error encountered. */

    *ier = 0;
    return 0;

/* Invalid input parameter LUN, PLTSIZ, or N. */

L11:
    *ier = 1;
    return 0;

/* Invalid input parameter ELAT, ELON, or A. */

L12:
    *ier = 2;
    return 0;

/* Error writing to unit LUN. */

L13:
    *ier = 3;
    return 0;
} /* trplot_ */

/* Subroutine */ int trprnt_(n, x, y, z__, iflag, list, lptr, lend, lout)
integer *n;
doublereal *x, *y, *z__;
integer *iflag, *list, *lptr, *lend, *lout;
{
    /* Initialized data */

    static integer nmax = 9999;
    static integer nlmax = 58;

    /* Format strings */
    static char fmt_100[] = "(///15x,\002STRIPACK Triangulation Data \002\
,\002Structure,  N = \002,i5//)";
    static char fmt_110[] = "(1x,10x,\002*** N is outside its valid\002,\002\
 range ***\002)";
    static char fmt_101[] = "(1x,\002Node\002,31x,\002Neighbors of Node\002/\
/)";
    static char fmt_108[] = "(///)";
    static char fmt_104[] = "(1x,i4,4x,14i5/(1x,8x,14i5))";
    static char fmt_107[] = "(1x)";
    static char fmt_102[] = "(1x,\002Node\002,5x,\002Longitude\002,6x,\002La\
titude\002,18x,\002Neighbors of Node\002//)";
    static char fmt_105[] = "(1x,i4,2d15.6,4x,8i5/(1x,38x,8i5))";
    static char fmt_103[] = "(1x,\002Node\002,5x,\002X(Node)\002,8x,\002Y(No\
de)\002,8x,\002Z(Node)\002,11x,\002Neighbors of Node\002//)";
    static char fmt_106[] = "(1x,i4,3d15.6,4x,5i5/(1x,53x,5i5))";
    static char fmt_109[] = "(/1x,\002NB = \002,i4,\002 Boundary Nodes\002,5\
x,\002NA = \002,i5,\002 Arcs\002,5x,\002NT = \002,i5,\002 Triangles\002)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    static integer node, i__, k, nabor[400], na, nb, nd, nl, lp, nn, nt, inc, 
	    lpl, lun;

    /* Fortran I/O blocks */
    static cilist io___527 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___528 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___531 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___539 = { 0, 0, 0, fmt_108, 0 };
    static cilist io___540 = { 0, 0, 0, fmt_104, 0 };
    static cilist io___542 = { 0, 0, 0, fmt_107, 0 };
    static cilist io___543 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___544 = { 0, 0, 0, fmt_108, 0 };
    static cilist io___545 = { 0, 0, 0, fmt_105, 0 };
    static cilist io___546 = { 0, 0, 0, fmt_107, 0 };
    static cilist io___547 = { 0, 0, 0, fmt_103, 0 };
    static cilist io___548 = { 0, 0, 0, fmt_108, 0 };
    static cilist io___549 = { 0, 0, 0, fmt_106, 0 };
    static cilist io___550 = { 0, 0, 0, fmt_107, 0 };
    static cilist io___553 = { 0, 0, 0, fmt_109, 0 };



/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/25/98 */

/*   This subroutine prints the triangulation adjacency lists */
/* created by Subroutine TRMESH and, optionally, the nodal */
/* coordinates (either latitude and longitude or Cartesian */
/* coordinates) on logical unit LOUT.  The list of neighbors */
/* of a boundary node is followed by index 0.  The numbers of */
/* boundary nodes, triangles, and arcs are also printed. */


/* On input: */

/*       N = Number of nodes in the triangulation.  N .GE. 3 */
/*           and N .LE. 9999. */

/*       X,Y,Z = Arrays of length N containing the Cartesian */
/*               coordinates of the nodes if IFLAG = 0, or */
/*               (X and Y only) arrays of length N containing */
/*               longitude and latitude, respectively, if */
/*               IFLAG > 0, or unused dummy parameters if */
/*               IFLAG < 0. */

/*       IFLAG = Nodal coordinate option indicator: */
/*               IFLAG = 0 if X, Y, and Z (assumed to contain */
/*                         Cartesian coordinates) are to be */
/*                         printed (to 6 decimal places). */
/*               IFLAG > 0 if only X and Y (assumed to con- */
/*                         tain longitude and latitude) are */
/*                         to be printed (to 6 decimal */
/*                         places). */
/*               IFLAG < 0 if only the adjacency lists are to */
/*                         be printed. */

/*       LIST,LPTR,LEND = Data structure defining the trian- */
/*                        gulation.  Refer to Subroutine */
/*                        TRMESH. */

/*       LOUT = Logical unit for output.  If LOUT is not in */
/*              the range 0 to 99, output is written to */
/*              logical unit 6. */

/* Input parameters are not altered by this routine. */

/* On output: */

/*   The adjacency lists and nodal coordinates (as specified */
/* by IFLAG) are written to unit LOUT. */

/* Modules required by TRPRNT:  None */

/* *********************************************************** */

    /* Parameter adjustments */
    --lend;
    --z__;
    --y;
    --x;
    --list;
    --lptr;

    /* Function Body */

/* Local parameters: */

/* I =     NABOR index (1 to K) */
/* INC =   Increment for NL associated with an adjacency list */
/* K =     Counter and number of neighbors of NODE */
/* LP =    LIST pointer of a neighbor of NODE */
/* LPL =   Pointer to the last neighbor of NODE */
/* LUN =   Logical unit for output (copy of LOUT) */
/* NA =    Number of arcs in the triangulation */
/* NABOR = Array containing the adjacency list associated */
/*           with NODE, with zero appended if NODE is a */
/*           boundary node */
/* NB =    Number of boundary nodes encountered */
/* ND =    Index of a neighbor of NODE (or negative index) */
/* NL =    Number of lines that have been printed on the */
/*           current page */
/* NLMAX = Maximum number of print lines per page (except */
/*           for the last page which may have two addi- */
/*           tional lines) */
/* NMAX =  Upper bound on N (allows 4-digit indexes) */
/* NODE =  Index of a node and DO-loop index (1 to N) */
/* NN =    Local copy of N */
/* NT =    Number of triangles in the triangulation */

    nn = *n;
    lun = *lout;
    if (lun < 0 || lun > 99) {
	lun = 6;
    }

/* Print a heading and test the range of N. */

    io___527.ciunit = lun;
    s_wsfe(&io___527);
    do_fio(&c__1, (char *)&nn, (ftnlen)sizeof(integer));
    e_wsfe();
    if (nn < 3 || nn > nmax) {

/* N is outside its valid range. */

	io___528.ciunit = lun;
	s_wsfe(&io___528);
	e_wsfe();
	return 0;
    }

/* Initialize NL (the number of lines printed on the current */
/*   page) and NB (the number of boundary nodes encountered). */

    nl = 6;
    nb = 0;
    if (*iflag < 0) {

/* Print LIST only.  K is the number of neighbors of NODE */
/*   that have been stored in NABOR. */

	io___531.ciunit = lun;
	s_wsfe(&io___531);
	e_wsfe();
	i__1 = nn;
	for (node = 1; node <= i__1; ++node) {
	    lpl = lend[node];
	    lp = lpl;
	    k = 0;

L1:
	    ++k;
	    lp = lptr[lp];
	    nd = list[lp];
	    nabor[k - 1] = nd;
	    if (lp != lpl) {
		goto L1;
	    }
	    if (nd <= 0) {

/*   NODE is a boundary node.  Correct the sign of the last */
/*     neighbor, add 0 to the end of the list, and increment */
/*     NB. */

		nabor[k - 1] = -nd;
		++k;
		nabor[k - 1] = 0;
		++nb;
	    }

/*   Increment NL and print the list of neighbors. */

	    inc = (k - 1) / 14 + 2;
	    nl += inc;
	    if (nl > nlmax) {
		io___539.ciunit = lun;
		s_wsfe(&io___539);
		e_wsfe();
		nl = inc;
	    }
	    io___540.ciunit = lun;
	    s_wsfe(&io___540);
	    do_fio(&c__1, (char *)&node, (ftnlen)sizeof(integer));
	    i__2 = k;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		do_fio(&c__1, (char *)&nabor[i__ - 1], (ftnlen)sizeof(integer)
			);
	    }
	    e_wsfe();
	    if (k != 14) {
		io___542.ciunit = lun;
		s_wsfe(&io___542);
		e_wsfe();
	    }
/* L2: */
	}
    } else if (*iflag > 0) {

/* Print X (longitude), Y (latitude), and LIST. */

	io___543.ciunit = lun;
	s_wsfe(&io___543);
	e_wsfe();
	i__1 = nn;
	for (node = 1; node <= i__1; ++node) {
	    lpl = lend[node];
	    lp = lpl;
	    k = 0;

L3:
	    ++k;
	    lp = lptr[lp];
	    nd = list[lp];
	    nabor[k - 1] = nd;
	    if (lp != lpl) {
		goto L3;
	    }
	    if (nd <= 0) {

/*   NODE is a boundary node. */

		nabor[k - 1] = -nd;
		++k;
		nabor[k - 1] = 0;
		++nb;
	    }

/*   Increment NL and print X, Y, and NABOR. */

	    inc = (k - 1) / 8 + 2;
	    nl += inc;
	    if (nl > nlmax) {
		io___544.ciunit = lun;
		s_wsfe(&io___544);
		e_wsfe();
		nl = inc;
	    }
	    io___545.ciunit = lun;
	    s_wsfe(&io___545);
	    do_fio(&c__1, (char *)&node, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&x[node], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&y[node], (ftnlen)sizeof(doublereal));
	    i__2 = k;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		do_fio(&c__1, (char *)&nabor[i__ - 1], (ftnlen)sizeof(integer)
			);
	    }
	    e_wsfe();
	    if (k != 8) {
		io___546.ciunit = lun;
		s_wsfe(&io___546);
		e_wsfe();
	    }
/* L4: */
	}
    } else {

/* Print X, Y, Z, and LIST. */

	io___547.ciunit = lun;
	s_wsfe(&io___547);
	e_wsfe();
	i__1 = nn;
	for (node = 1; node <= i__1; ++node) {
	    lpl = lend[node];
	    lp = lpl;
	    k = 0;

L5:
	    ++k;
	    lp = lptr[lp];
	    nd = list[lp];
	    nabor[k - 1] = nd;
	    if (lp != lpl) {
		goto L5;
	    }
	    if (nd <= 0) {

/*   NODE is a boundary node. */

		nabor[k - 1] = -nd;
		++k;
		nabor[k - 1] = 0;
		++nb;
	    }

/*   Increment NL and print X, Y, Z, and NABOR. */

	    inc = (k - 1) / 5 + 2;
	    nl += inc;
	    if (nl > nlmax) {
		io___548.ciunit = lun;
		s_wsfe(&io___548);
		e_wsfe();
		nl = inc;
	    }
	    io___549.ciunit = lun;
	    s_wsfe(&io___549);
	    do_fio(&c__1, (char *)&node, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&x[node], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&y[node], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&z__[node], (ftnlen)sizeof(doublereal));
	    i__2 = k;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		do_fio(&c__1, (char *)&nabor[i__ - 1], (ftnlen)sizeof(integer)
			);
	    }
	    e_wsfe();
	    if (k != 5) {
		io___550.ciunit = lun;
		s_wsfe(&io___550);
		e_wsfe();
	    }
/* L6: */
	}
    }

/* Print NB, NA, and NT (boundary nodes, arcs, and */
/*   triangles). */

    if (nb != 0) {
	na = nn * 3 - nb - 3;
	nt = (nn << 1) - nb - 2;
    } else {
	na = nn * 3 - 6;
	nt = (nn << 1) - 4;
    }
    io___553.ciunit = lun;
    s_wsfe(&io___553);
    do_fio(&c__1, (char *)&nb, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&na, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nt, (ftnlen)sizeof(integer));
    e_wsfe();
    return 0;

/* Print formats: */

} /* trprnt_ */

/* Subroutine */ int vrplot_(lun, pltsiz, elat, elon, a, n, x, y, z__, nt, 
	listc, lptr, lend, xc, yc, zc, title, numbr, ier, title_len)
integer *lun;
doublereal *pltsiz, *elat, *elon, *a;
integer *n;
doublereal *x, *y, *z__;
integer *nt, *listc, *lptr, *lend;
doublereal *xc, *yc, *zc;
char *title;
logical *numbr;
integer *ier;
ftnlen title_len;
{
    /* Initialized data */

    static logical annot = TRUE_;
    static doublereal fsizn = 10.;
    static doublereal fsizt = 16.;
    static doublereal tol = .5;

    /* Format strings */
    static char fmt_100[] = "(\002%!PS-Adobe-3.0 EPSF-3.0\002/\002%%Bounding\
Box:\002,4i4/\002%%Title:  Voronoi diagram\002/\002%%Creator:  STRIPACK\002\
/\002%%EndComments\002)";
    static char fmt_110[] = "(f12.6,\002 setlinewidth\002)";
    static char fmt_120[] = "(\002306 396 \002,i3,\002 0 360 arc\002)";
    static char fmt_130[] = "(\002stroke\002)";
    static char fmt_140[] = "(2f12.6,\002 translate\002/2f12.6,\002 scale\
\002)";
    static char fmt_150[] = "(\002gsave\002)";
    static char fmt_160[] = "(\0020 0 \002,f12.6,\002 0 360 arc\002)";
    static char fmt_170[] = "(\002clip newpath\002)";
    static char fmt_190[] = "(\002grestore\002)";
    static char fmt_200[] = "(\002/Helvetica findfont\002/f12.6,\002 scalefo\
nt setfont\002)";
    static char fmt_210[] = "(2f12.6,\002 moveto\002)";
    static char fmt_220[] = "(\002(\002,i3,\002) show\002)";
    static char fmt_230[] = "(a80/\002  stringwidth pop 2 div neg \002,f12\
.6,\002 moveto\002)";
    static char fmt_240[] = "(a80/\002  show\002)";
    static char fmt_250[] = "(\002(Window center:  ELAT = \002,f7.2,\002,  E\
LON = \002,f8.2,\002) show\002)";
    static char fmt_260[] = "(\002(Angular extent:  A = \002,f5.2,\002) sho\
w\002)";
    static char fmt_270[] = "(\002stroke\002/\002showpage\002/\002%%EOF\002)";
    static char fmt_280[] = "(a1)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double atan(), sin();
    integer i_dnnt(), s_wsfe(), do_fio(), e_wsfe();
    double cos(), sqrt();

    /* Local variables */
    static integer nseg;
    static doublereal t;
    static integer n0;
    static doublereal p1[3], p2[3], x0, y0, cf, r11, r12, r21, ct, r22, r23, 
	    sf;
    static integer ir, lp;
    static doublereal ex, ey, ez, wr, tx, ty;
    extern /* Subroutine */ int drwarc_();
    static logical in1, in2;
    static integer kv1, kv2, lpl;
    static doublereal wrs;
    static integer ipx1, ipx2, ipy1, ipy2;

    /* Fortran I/O blocks */
    static cilist io___566 = { 1, 0, 0, fmt_100, 0 };
    static cilist io___568 = { 1, 0, 0, fmt_110, 0 };
    static cilist io___569 = { 1, 0, 0, fmt_120, 0 };
    static cilist io___570 = { 1, 0, 0, fmt_130, 0 };
    static cilist io___574 = { 1, 0, 0, fmt_140, 0 };
    static cilist io___575 = { 1, 0, 0, fmt_110, 0 };
    static cilist io___576 = { 1, 0, 0, fmt_150, 0 };
    static cilist io___577 = { 1, 0, 0, fmt_160, 0 };
    static cilist io___578 = { 1, 0, 0, fmt_170, 0 };
    static cilist io___598 = { 1, 0, 0, fmt_130, 0 };
    static cilist io___599 = { 1, 0, 0, fmt_190, 0 };
    static cilist io___600 = { 1, 0, 0, fmt_200, 0 };
    static cilist io___603 = { 1, 0, 0, fmt_210, 0 };
    static cilist io___604 = { 1, 0, 0, fmt_220, 0 };
    static cilist io___605 = { 1, 0, 0, fmt_200, 0 };
    static cilist io___606 = { 1, 0, 0, fmt_230, 0 };
    static cilist io___607 = { 1, 0, 0, fmt_240, 0 };
    static cilist io___608 = { 1, 0, 0, fmt_210, 0 };
    static cilist io___609 = { 1, 0, 0, fmt_250, 0 };
    static cilist io___610 = { 1, 0, 0, fmt_210, 0 };
    static cilist io___611 = { 1, 0, 0, fmt_260, 0 };
    static cilist io___612 = { 1, 0, 0, fmt_270, 0 };
    static cilist io___613 = { 1, 0, 0, fmt_280, 0 };



/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   03/04/03 */

/*   This subroutine creates a level-2 Encapsulated Post- */
/* script (EPS) file containing a graphical depiction of a */
/* Voronoi diagram of a set of nodes on the unit sphere. */
/* The visible portion of the diagram is projected orthog- */
/* onally onto the plane that contains the origin and has */
/* normal defined by a user-specified eye-position. */

/*   The parameters defining the Voronoi diagram may be com- */
/* puted by Subroutine CRLIST. */


/* On input: */

/*       LUN = Logical unit number in the range 0 to 99. */
/*             The unit should be opened with an appropriate */
/*             file name before the call to this routine. */

/*       PLTSIZ = Plot size in inches.  A circular window in */
/*                the projection plane is mapped to a circu- */
/*                lar viewport with diameter equal to .88* */
/*                PLTSIZ (leaving room for labels outside the */
/*                viewport).  The viewport is centered on the */
/*                8.5 by 11 inch page, and its boundary is */
/*                drawn.  1.0 .LE. PLTSIZ .LE. 8.5. */

/*       ELAT,ELON = Latitude and longitude (in degrees) of */
/*                   the center of projection E (the center */
/*                   of the plot).  The projection plane is */
/*                   the plane that contains the origin and */
/*                   has E as unit normal.  In a rotated */
/*                   coordinate system for which E is the */
/*                   north pole, the projection plane con- */
/*                   tains the equator, and only northern */
/*                   hemisphere points are visible (from the */
/*                   point at infinity in the direction E). */
/*                   These are projected orthogonally onto */
/*                   the projection plane (by zeroing the z- */
/*                   component in the rotated coordinate */
/*                   system).  ELAT and ELON must be in the */
/*                   range -90 to 90 and -180 to 180, respec- */
/*                   tively. */

/*       A = Angular distance in degrees from E to the boun- */
/*           dary of a circular window against which the */
/*           Voronoi diagram is clipped.  The projected win- */
/*           dow is a disk of radius r = Sin(A) centered at */
/*           the origin, and only visible vertices whose */
/*           projections are within distance r of the origin */
/*           are included in the plot.  Thus, if A = 90, the */
/*           plot includes the entire hemisphere centered at */
/*           E.  0 .LT. A .LE. 90. */

/*       N = Number of nodes (Voronoi centers) and Voronoi */
/*           regions.  N .GE. 3. */

/*       X,Y,Z = Arrays of length N containing the Cartesian */
/*               coordinates of the nodes (unit vectors). */

/*       NT = Number of Voronoi region vertices (triangles, */
/*            including those in the extended triangulation */
/*            if the number of boundary nodes NB is nonzero): */
/*            NT = 2*N-4. */

/*       LISTC = Array of length 3*NT containing triangle */
/*               indexes (indexes to XC, YC, and ZC) stored */
/*               in 1-1 correspondence with LIST/LPTR entries */
/*               (or entries that would be stored in LIST for */
/*               the extended triangulation):  the index of */
/*               triangle (N1,N2,N3) is stored in LISTC(K), */
/*               LISTC(L), and LISTC(M), where LIST(K), */
/*               LIST(L), and LIST(M) are the indexes of N2 */
/*               as a neighbor of N1, N3 as a neighbor of N2, */
/*               and N1 as a neighbor of N3.  The Voronoi */
/*               region associated with a node is defined by */
/*               the CCW-ordered sequence of circumcenters in */
/*               one-to-one correspondence with its adjacency */
/*               list (in the extended triangulation). */

/*       LPTR = Array of length 3*NT = 6*N-12 containing a */
/*              set of pointers (LISTC indexes) in one-to-one */
/*              correspondence with the elements of LISTC. */
/*              LISTC(LPTR(I)) indexes the triangle which */
/*              follows LISTC(I) in cyclical counterclockwise */
/*              order (the first neighbor follows the last */
/*              neighbor). */

/*       LEND = Array of length N containing a set of */
/*              pointers to triangle lists.  LP = LEND(K) */
/*              points to a triangle (indexed by LISTC(LP)) */
/*              containing node K for K = 1 to N. */

/*       XC,YC,ZC = Arrays of length NT containing the */
/*                  Cartesian coordinates of the triangle */
/*                  circumcenters (Voronoi vertices). */
/*                  XC(I)**2 + YC(I)**2 + ZC(I)**2 = 1. */

/*       TITLE = Type CHARACTER variable or constant contain- */
/*               ing a string to be centered above the plot. */
/*               The string must be enclosed in parentheses; */
/*               i.e., the first and last characters must be */
/*               '(' and ')', respectively, but these are not */
/*               displayed.  TITLE may have at most 80 char- */
/*               acters including the parentheses. */

/*       NUMBR = Option indicator:  If NUMBR = TRUE, the */
/*               nodal indexes are plotted at the Voronoi */
/*               region centers. */

/* Input parameters are not altered by this routine. */

/* On output: */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered. */
/*             IER = 1 if LUN, PLTSIZ, N, or NT is outside */
/*                     its valid range. */
/*             IER = 2 if ELAT, ELON, or A is outside its */
/*                     valid range. */
/*             IER = 3 if an error was encountered in writing */
/*                     to unit LUN. */

/* Module required by VRPLOT:  DRWARC */

/* Intrinsic functions called by VRPLOT:  ABS, ATAN, COS, */
/*                                          DBLE, NINT, SIN, */
/*                                          SQRT */

/* *********************************************************** */


    /* Parameter adjustments */
    --lend;
    --z__;
    --y;
    --x;
    --zc;
    --yc;
    --xc;
    --listc;
    --lptr;

    /* Function Body */

/* Local parameters: */

/* ANNOT =     Logical variable with value TRUE iff the plot */
/*               is to be annotated with the values of ELAT, */
/*               ELON, and A */
/* CF =        Conversion factor for degrees to radians */
/* CT =        Cos(ELAT) */
/* EX,EY,EZ =  Cartesian coordinates of the eye-position E */
/* FSIZN =     Font size in points for labeling nodes with */
/*               their indexes if NUMBR = TRUE */
/* FSIZT =     Font size in points for the title (and */
/*               annotation if ANNOT = TRUE) */
/* IN1,IN2 =   Logical variables with value TRUE iff the */
/*               projections of vertices KV1 and KV2, respec- */
/*               tively, are inside the window */
/* IPX1,IPY1 = X and y coordinates (in points) of the lower */
/*               left corner of the bounding box or viewport */
/*               box */
/* IPX2,IPY2 = X and y coordinates (in points) of the upper */
/*               right corner of the bounding box or viewport */
/*               box */
/* IR =        Half the width (height) of the bounding box or */
/*               viewport box in points -- viewport radius */
/* KV1,KV2 =   Endpoint indexes of a Voronoi edge */
/* LP =        LIST index (pointer) */
/* LPL =       Pointer to the last neighbor of N0 */
/* N0 =        Index of a node */
/* NSEG =      Number of line segments used by DRWARC in a */
/*               polygonal approximation to a projected edge */
/* P1 =        Coordinates of vertex KV1 in the rotated */
/*               coordinate system */
/* P2 =        Coordinates of vertex KV2 in the rotated */
/*               coordinate system or intersection of edge */
/*               KV1-KV2 with the equator (in the rotated */
/*               coordinate system) */
/* R11...R23 = Components of the first two rows of a rotation */
/*               that maps E to the north pole (0,0,1) */
/* SF =        Scale factor for mapping world coordinates */
/*               (window coordinates in [-WR,WR] X [-WR,WR]) */
/*               to viewport coordinates in [IPX1,IPX2] X */
/*               [IPY1,IPY2] */
/* T =         Temporary variable */
/* TOL =       Maximum distance in points between a projected */
/*               Voronoi edge and its approximation by a */
/*               polygonal curve */
/* TX,TY =     Translation vector for mapping world coordi- */
/*               nates to viewport coordinates */
/* WR =        Window radius r = Sin(A) */
/* WRS =       WR**2 */
/* X0,Y0 =     Projection plane coordinates of node N0 or */
/*               label location */


/* Test for invalid parameters. */

    if (*lun < 0 || *lun > 99 || *pltsiz < 1. || *pltsiz > 8.5 || *n < 3 || *
	    nt != 2 * *n - 4) {
	goto L11;
    }
    if (abs(*elat) > 90. || abs(*elon) > 180. || *a > 90.) {
	goto L12;
    }

/* Compute a conversion factor CF for degrees to radians */
/*   and compute the window radius WR. */

    cf = atan(1.) / 45.;
    wr = sin(cf * *a);
    wrs = wr * wr;

/* Compute the lower left (IPX1,IPY1) and upper right */
/*   (IPX2,IPY2) corner coordinates of the bounding box. */
/*   The coordinates, specified in default user space units */
/*   (points, at 72 points/inch with origin at the lower */
/*   left corner of the page), are chosen to preserve the */
/*   square aspect ratio, and to center the plot on the 8.5 */
/*   by 11 inch page.  The center of the page is (306,396), */
/*   and IR = PLTSIZ/2 in points. */

    d__1 = *pltsiz * 36.;
    ir = i_dnnt(&d__1);
    ipx1 = 306 - ir;
    ipx2 = ir + 306;
    ipy1 = 396 - ir;
    ipy2 = ir + 396;

/* Output header comments. */

    io___566.ciunit = *lun;
    i__1 = s_wsfe(&io___566);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&ipx1, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&ipy1, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&ipx2, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&ipy2, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }

/* Set (IPX1,IPY1) and (IPX2,IPY2) to the corner coordinates */
/*   of a viewport box obtained by shrinking the bounding box */
/*   by 12% in each dimension. */

    d__1 = (doublereal) ir * .88;
    ir = i_dnnt(&d__1);
    ipx1 = 306 - ir;
    ipx2 = ir + 306;
    ipy1 = 396 - ir;
    ipy2 = ir + 396;

/* Set the line thickness to 2 points, and draw the */
/*   viewport boundary. */

    t = 2.;
    io___568.ciunit = *lun;
    i__1 = s_wsfe(&io___568);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&t, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }
    io___569.ciunit = *lun;
    i__1 = s_wsfe(&io___569);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&ir, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }
    io___570.ciunit = *lun;
    i__1 = s_wsfe(&io___570);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }

/* Set up an affine mapping from the window box [-WR,WR] X */
/*   [-WR,WR] to the viewport box. */

    sf = (doublereal) ir / wr;
    tx = ipx1 + sf * wr;
    ty = ipy1 + sf * wr;
    io___574.ciunit = *lun;
    i__1 = s_wsfe(&io___574);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&tx, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&ty, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&sf, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&sf, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }

/* The line thickness must be changed to reflect the new */
/*   scaling which is applied to all subsequent output. */
/*   Set it to 1.0 point. */

    t = 1. / sf;
    io___575.ciunit = *lun;
    i__1 = s_wsfe(&io___575);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&t, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }

/* Save the current graphics state, and set the clip path to */
/*   the boundary of the window. */

    io___576.ciunit = *lun;
    i__1 = s_wsfe(&io___576);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }
    io___577.ciunit = *lun;
    i__1 = s_wsfe(&io___577);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&wr, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }
    io___578.ciunit = *lun;
    i__1 = s_wsfe(&io___578);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }

/* Compute the Cartesian coordinates of E and the components */
/*   of a rotation R which maps E to the north pole (0,0,1). */
/*   R is taken to be a rotation about the z-axis (into the */
/*   yz-plane) followed by a rotation about the x-axis chosen */
/*   so that the view-up direction is (0,0,1), or (-1,0,0) if */
/*   E is the north or south pole. */

/*           ( R11  R12  0   ) */
/*       R = ( R21  R22  R23 ) */
/*           ( EX   EY   EZ  ) */

    t = cf * *elon;
    ct = cos(cf * *elat);
    ex = ct * cos(t);
    ey = ct * sin(t);
    ez = sin(cf * *elat);
    if (ct != 0.) {
	r11 = -ey / ct;
	r12 = ex / ct;
    } else {
	r11 = 0.;
	r12 = 1.;
    }
    r21 = -ez * r12;
    r22 = ez * r11;
    r23 = ct;

/* Loop on nodes (Voronoi centers) N0. */
/*   LPL indexes the last neighbor of N0. */

    i__1 = *n;
    for (n0 = 1; n0 <= i__1; ++n0) {
	lpl = lend[n0];

/* Set KV2 to the first (and last) vertex index and compute */
/*   its coordinates P2 in the rotated coordinate system. */

	kv2 = listc[lpl];
	p2[0] = r11 * xc[kv2] + r12 * yc[kv2];
	p2[1] = r21 * xc[kv2] + r22 * yc[kv2] + r23 * zc[kv2];
	p2[2] = ex * xc[kv2] + ey * yc[kv2] + ez * zc[kv2];

/*   IN2 = TRUE iff KV2 is in the window. */

	in2 = p2[2] >= 0. && p2[0] * p2[0] + p2[1] * p2[1] <= wrs;

/* Loop on neighbors N1 of N0.  For each triangulation edge */
/*   N0-N1, KV1-KV2 is the corresponding Voronoi edge. */

	lp = lpl;
L1:
	lp = lptr[lp];
	kv1 = kv2;
	p1[0] = p2[0];
	p1[1] = p2[1];
	p1[2] = p2[2];
	in1 = in2;
	kv2 = listc[lp];

/*   Compute the new values of P2 and IN2. */

	p2[0] = r11 * xc[kv2] + r12 * yc[kv2];
	p2[1] = r21 * xc[kv2] + r22 * yc[kv2] + r23 * zc[kv2];
	p2[2] = ex * xc[kv2] + ey * yc[kv2] + ez * zc[kv2];
	in2 = p2[2] >= 0. && p2[0] * p2[0] + p2[1] * p2[1] <= wrs;

/* Add edge KV1-KV2 to the path iff both endpoints are inside */
/*   the window and KV2 > KV1, or KV1 is inside and KV2 is */
/*   outside (so that the edge is drawn only once). */

	if (! in1 || in2 && kv2 <= kv1) {
	    goto L2;
	}
	if (p2[2] < 0.) {

/*   KV2 is a 'southern hemisphere' point.  Move it to the */
/*     intersection of edge KV1-KV2 with the equator so that */
/*     the edge is clipped properly.  P2(3) is set to 0. */

	    p2[0] = p1[2] * p2[0] - p2[2] * p1[0];
	    p2[1] = p1[2] * p2[1] - p2[2] * p1[1];
	    t = sqrt(p2[0] * p2[0] + p2[1] * p2[1]);
	    p2[0] /= t;
	    p2[1] /= t;
	}

/*   Add the edge to the path.  (TOL is converted to world */
/*     coordinates.) */

	if (p2[2] < 0.) {
	    p2[2] = (float)0.;
	}
	d__1 = tol / sf;
	drwarc_(lun, p1, p2, &d__1, &nseg);

/* Bottom of loops. */

L2:
	if (lp != lpl) {
	    goto L1;
	}
/* L3: */
    }

/* Paint the path and restore the saved graphics state (with */
/*   no clip path). */

    io___598.ciunit = *lun;
    i__1 = s_wsfe(&io___598);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }
    io___599.ciunit = *lun;
    i__1 = s_wsfe(&io___599);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }
    if (*numbr) {

/* Nodes in the window are to be labeled with their indexes. */
/*   Convert FSIZN from points to world coordinates, and */
/*   output the commands to select a font and scale it. */

	t = fsizn / sf;
	io___600.ciunit = *lun;
	i__1 = s_wsfe(&io___600);
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = do_fio(&c__1, (char *)&t, (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L13;
	}

/* Loop on visible nodes N0 that project to points (X0,Y0) in */
/*   the window. */

	i__1 = *n;
	for (n0 = 1; n0 <= i__1; ++n0) {
	    if (ex * x[n0] + ey * y[n0] + ez * z__[n0] < 0.) {
		goto L4;
	    }
	    x0 = r11 * x[n0] + r12 * y[n0];
	    y0 = r21 * x[n0] + r22 * y[n0] + r23 * z__[n0];
	    if (x0 * x0 + y0 * y0 > wrs) {
		goto L4;
	    }

/*   Move to (X0,Y0), and draw the label N0 with the origin */
/*     of the first character at (X0,Y0). */

	    io___603.ciunit = *lun;
	    i__2 = s_wsfe(&io___603);
	    if (i__2 != 0) {
		goto L13;
	    }
	    i__2 = do_fio(&c__1, (char *)&x0, (ftnlen)sizeof(doublereal));
	    if (i__2 != 0) {
		goto L13;
	    }
	    i__2 = do_fio(&c__1, (char *)&y0, (ftnlen)sizeof(doublereal));
	    if (i__2 != 0) {
		goto L13;
	    }
	    i__2 = e_wsfe();
	    if (i__2 != 0) {
		goto L13;
	    }
	    io___604.ciunit = *lun;
	    i__2 = s_wsfe(&io___604);
	    if (i__2 != 0) {
		goto L13;
	    }
	    i__2 = do_fio(&c__1, (char *)&n0, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L13;
	    }
	    i__2 = e_wsfe();
	    if (i__2 != 0) {
		goto L13;
	    }
L4:
	    ;
	}
    }

/* Convert FSIZT from points to world coordinates, and output */
/*   the commands to select a font and scale it. */

    t = fsizt / sf;
    io___605.ciunit = *lun;
    i__1 = s_wsfe(&io___605);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&t, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }

/* Display TITLE centered above the plot: */

    y0 = wr + t * 3.;
    io___606.ciunit = *lun;
    i__1 = s_wsfe(&io___606);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, title, title_len);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, (char *)&y0, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }
    io___607.ciunit = *lun;
    i__1 = s_wsfe(&io___607);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, title, title_len);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }
    if (annot) {

/* Display the window center and radius below the plot. */

	x0 = -wr;
	y0 = -wr - 50. / sf;
	io___608.ciunit = *lun;
	i__1 = s_wsfe(&io___608);
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = do_fio(&c__1, (char *)&x0, (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = do_fio(&c__1, (char *)&y0, (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L13;
	}
	io___609.ciunit = *lun;
	i__1 = s_wsfe(&io___609);
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = do_fio(&c__1, (char *)&(*elat), (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = do_fio(&c__1, (char *)&(*elon), (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L13;
	}
	y0 -= t * 2.;
	io___610.ciunit = *lun;
	i__1 = s_wsfe(&io___610);
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = do_fio(&c__1, (char *)&x0, (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = do_fio(&c__1, (char *)&y0, (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L13;
	}
	io___611.ciunit = *lun;
	i__1 = s_wsfe(&io___611);
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = do_fio(&c__1, (char *)&(*a), (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L13;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L13;
	}
    }

/* Paint the path and output the showpage command and */
/*   end-of-file indicator. */

    io___612.ciunit = *lun;
    i__1 = s_wsfe(&io___612);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }

/* HP's interpreters require a one-byte End-of-PostScript-Job */
/*   indicator (to eliminate a timeout error message): */
/*   ASCII 4. */

    io___613.ciunit = *lun;
    i__1 = s_wsfe(&io___613);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = do_fio(&c__1, "\004", (ftnlen)1);
    if (i__1 != 0) {
	goto L13;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L13;
    }

/* No error encountered. */

    *ier = 0;
    return 0;

/* Invalid input parameter LUN, PLTSIZ, N, or NT. */

L11:
    *ier = 1;
    return 0;

/* Invalid input parameter ELAT, ELON, or A. */

L12:
    *ier = 2;
    return 0;

/* Error writing to unit LUN. */

L13:
    *ier = 3;
    return 0;
} /* vrplot_ */

/* Subroutine */ int random_(ix, iy, iz, rannum)
integer *ix, *iy, *iz;
doublereal *rannum;
{
    static doublereal x;


/*   This routine returns pseudo-random numbers uniformly */
/* distributed in the interval (0,1).  Integer seeds IX, IY, */
/* and IZ should be initialized to values in the range 1 to */
/* 30,000 before the first call to RANDOM, and should not */
/* be altered between subsequent calls (unless a sequence */
/* of random numbers is to be repeated by reinitializing the */
/* seeds). */

/* Reference:  B. A. Wichmann and I. D. Hill, An Efficient */
/*             and Portable Pseudo-random Number Generator, */
/*             Applied Statistics, Vol. 31, No. 2, 1982, */
/*             pp. 188-190. */

    *ix = *ix * 171 % 30269;
    *iy = *iy * 172 % 30307;
    *iz = *iz * 170 % 30323;
    x = (doublereal) (*ix) / 30269. + (doublereal) (*iy) / 30307. + (
	    doublereal) (*iz) / 30323.;
    *rannum = x - (integer) x;
    return 0;
} /* random_ */


#endif
