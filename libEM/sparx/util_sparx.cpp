/**
 * $Id$
 */

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
#ifdef _WIN32
#pragma warning(disable:4819)
#include <malloc.h>
#endif	//_WIN32

#include <cstring>
#include <ctime>
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <sstream>
#include "emdata.h"
#include "util.h"
#include "fundamentals.h"
#include "lapackblas.h"
#include "lbfgsb.h"
using namespace EMAN;
#include "steepest.h"
#include "emassert.h"
#include "randnum.h"
#include "mcqd.h"

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_bessel.h>
#include <cmath>
//#include <omp.h>
using namespace std;
using std::complex;

/* Subroutine */ 
int circum_(double *, double *, double *, double *, int *);
long int left_(double *, double *, double *, double *, double *, double *, double *, double *, double *);
int addnod_(int *, int *, double *, double *, double *, int *, int *, int *, int *, int *);

vector<float> Util::infomask(EMData* Vol, EMData* mask, bool flip = false)
//  flip true:  find statistics under the mask (mask >0.5)
//  flip false: find statistics ourside the mask (mask <0.5)
{
	ENTERFUNC;
	vector<float> stats;
	float *Volptr, *maskptr,MAX,MIN;
	long double Sum1,Sum2;
	long count;

	MAX = -FLT_MAX;
	MIN =  FLT_MAX;
	count = 0L;
	Sum1  = 0.0L;
	Sum2  = 0.0L;

	if (mask == NULL) {
		//Vol->update_stat();
		stats.push_back(Vol->get_attr("mean"));
		stats.push_back(Vol->get_attr("sigma"));
		stats.push_back(Vol->get_attr("minimum"));
		stats.push_back(Vol->get_attr("maximum"));
		return stats;
	}

	/* Check if the sizes of the mask and image are same */

	size_t nx = Vol->get_xsize();
	size_t ny = Vol->get_ysize();
	size_t nz = Vol->get_zsize();

	size_t mask_nx = mask->get_xsize();
	size_t mask_ny = mask->get_ysize();
	size_t mask_nz = mask->get_zsize();

	if  (nx != mask_nx || ny != mask_ny || nz != mask_nz )
		throw ImageDimensionException("The dimension of the image does not match the dimension of the mask!");

 /*       if (nx != mask_nx ||
            ny != mask_ny ||
            nz != mask_nz  ) {
           // should throw an exception here!!! (will clean it up later CY)
           fprintf(stderr, "The dimension of the image does not match the dimension of the mask!\n");
           fprintf(stderr, " nx = %d, mask_nx = %d\n", nx, mask_nx);
           fprintf(stderr, " ny = %d, mask_ny = %d\n", ny, mask_ny);
           fprintf(stderr, " nz = %d, mask_nz = %d\n", nz, mask_nz);
           exit(1);
        }
 */
	Volptr = Vol->get_data();
	maskptr = mask->get_data();

	for (size_t i = 0; i < (size_t)nx*ny*nz; ++i) {
		if ((maskptr[i]>0.5f) == flip) {
			Sum1 += Volptr[i];
			Sum2 += Volptr[i]*double(Volptr[i]);
			MAX = (MAX < Volptr[i])?Volptr[i]:MAX;
			MIN = (MIN > Volptr[i])?Volptr[i]:MIN;
			count++;
		}
	}

	if (count == 0) {
		LOGERR("Invalid mask");
		throw ImageFormatException( "Invalid mask");
	}

	float avg = static_cast<float>(Sum1/count);
	float sig = static_cast<float>(sqrt((Sum2 - Sum1*Sum1/count)/(count-1)));

	stats.push_back(avg);
	stats.push_back(sig);
	stats.push_back(MIN);
	stats.push_back(MAX);

	return stats;
}


vector<double> Util::helixshiftali(vector<EMData*> ctx, vector<vector<float> > pcoords, int nsegms, float maxincline, int kang, int search_rng, int nxc)
{
	int cents, six, incline, kim, ixl, sib;
	double q0, qt, qu, tang, dst, xl, dxl, qm, bang;
	vector<double> result;
	cents = nsegms/2;
	qm = -1.0e23;  
	for ( six =-search_rng; six <= search_rng; six++ ) {
		q0 = ctx[cents]->get_value_at(six+nxc);
		for ( incline = 0; incline <= kang; incline++ ) {
			qt = q0;
			qu = q0;
			if ( kang > 0 ) tang = tan(maxincline/kang*incline);
			else tang = 0.0;
			for ( kim = cents+1; kim < nsegms; kim++ ) {
				dst = sqrt(((double)pcoords[cents][0] - (double)pcoords[kim][0])*((double)pcoords[cents][0] - (double)pcoords[kim][0])+
				((double)pcoords[cents][1] - (double)pcoords[kim][1])*((double)pcoords[cents][1] - (double)pcoords[kim][1]));
				//printf("\ndst=%5.6f", dst); fflush(stdout);//cout <<"dst="<< setprecision(4) << dst<< endl;
				xl = dst * tang + six + nxc;
				ixl = floor(xl);
				dxl = xl - ixl;
				qt += (1.0-dxl)*ctx[kim]->get_value_at(ixl) + dxl*ctx[kim]->get_value_at(ixl+1);
				xl = -dst*tang+six+nxc;
				ixl = floor(xl);
				dxl = xl - ixl;
				qu += (1.0-dxl)*ctx[kim]->get_value_at(ixl) + dxl*ctx[kim]->get_value_at(ixl+1);
			}
			for ( kim = 0; kim < cents; kim++ ) {
			 	dst = sqrt(((double)pcoords[cents][0] - (double)pcoords[kim][0])*((double)pcoords[cents][0] - (double)pcoords[kim][0])+
				((double)pcoords[cents][1] - (double)pcoords[kim][1])*((double)pcoords[cents][1] - (double)pcoords[kim][1]));
				xl = -dst*tang+six+nxc;
				ixl = floor(xl);
				dxl = xl - ixl;
				qt += (1.0-dxl)*ctx[kim]->get_value_at(ixl) + dxl*ctx[kim]->get_value_at(ixl+1);
				xl = dst * tang + six + nxc;
				ixl = floor(xl);
				dxl = xl - ixl;
				qu += (1.0-dxl)*ctx[kim]->get_value_at(ixl) + dxl*ctx[kim]->get_value_at(ixl+1);
			}
			if ( qt > qm ) {
				qm = qt;
				sib = six;
				bang = tang;}
			if ( qu > qm ) {
				qm = qu;
				sib = six;
				bang = -tang;}	
		}
		
	
	
	}
	//printf("\nqm sib bang=%f %f %f", qm, sib, bang);    
    result.push_back((double)sib); 
    result.push_back(bang);
    result.push_back(qm);
    return result;    
}

vector<double> Util::snakeshiftali(vector<EMData*> sccf, vector<vector<float> > pcoords, int nsegms, float resamp_dang, int kang, int search_rng, int nxc, int maxrin)
{
	int cents, six, inplane, kim, ixl, ial, sib, nx, ny;
	double q0, qu0, qt, qu, tang, dst, xl, dxl, al, dal, qm, bang;
	char fln[200];
	vector<double> result;
	nx = sccf[0]->get_xsize();
	ny = sccf[0]->get_ysize();
	printf("nx=%d ny=%d search_rng=%d nxc=%d", nx, ny, search_rng, nxc);fflush(stdout);
	cents = nsegms/2;
	qm = -1.0e23;  
	// for ( int i = 0; i < nsegms; i++ ) {
// 		sprintf(fln,"sccfimage%d.hdf",i);
// 		sccf[i]->write_image(fln);}
	for ( six =-search_rng; six <= search_rng; six++ ) {
		for ( inplane = 0; inplane <= kang; inplane++ ) {
			if ( six+nxc > nx-1 || six+nxc < 0) continue;
			q0 = sccf[cents]->get_value_at(six+nxc,(maxrin-inplane)%maxrin);
			qu0 = sccf[cents]->get_value_at(six+nxc,inplane);
			qt = q0;
			qu = qu0;
			if ( kang > 0 ) tang = tan(resamp_dang*inplane);  //sin(resamp_dang*inplane);  //;
			else tang = 0.0;
			for ( kim = cents+1; kim < nsegms; kim++ ) {
				dst = sqrt(((double)pcoords[cents][0] - (double)pcoords[kim][0])*((double)pcoords[cents][0] - (double)pcoords[kim][0])+
				((double)pcoords[cents][1] - (double)pcoords[kim][1])*((double)pcoords[cents][1] - (double)pcoords[kim][1]));
				//printf("\ndst=%5.6f", dst); fflush(stdout);//cout <<"dst="<< setprecision(4) << dst<< endl;
				xl = -dst * tang + six + nxc;
				ixl = floor(xl);
				dxl = xl - ixl;
				al = (maxrin-inplane)%maxrin ;
				ial = floor(al);
				if (ixl +1 <= nx-1 && ixl >= 0 ) 
					qt += (1.0-dxl)*sccf[kim]->get_value_at(ixl,ial) + dxl*sccf[kim]->get_value_at(ixl+1,ial);
				else 
					qt = -1.0e23;
				//if ( six == -10 && inplane == 0 ) printf("\nkim=%d xl=%f dst=%f coefs=%f", kim, xl, dst, (1.0-dxl)*sccf[kim]->get_value_at(ixl,ial) + dxl*sccf[kim]->get_value_at(ixl+1,ial));
				
				xl =  dst * tang + six + nxc;
				ixl = floor(xl);
				dxl = xl - ixl;
				al = inplane ;
				ial = floor(al);
				if (ixl +1 <= nx-1 && ixl >= 0 ) 
					qu += (1.0-dxl)*sccf[kim]->get_value_at(ixl,ial) + dxl*sccf[kim]->get_value_at(ixl+1,ial);
				else qu = -1.0e23;	
			}
			for ( kim = 0; kim < cents; kim++ ) {
			 	dst = sqrt(((double)pcoords[cents][0] - (double)pcoords[kim][0])*((double)pcoords[cents][0] - (double)pcoords[kim][0])+
				((double)pcoords[cents][1] - (double)pcoords[kim][1])*((double)pcoords[cents][1] - (double)pcoords[kim][1]));
				//printf("\ndst=%5.6f", dst); fflush(stdout);
				xl = dst*tang+six+nxc;
				ixl = floor(xl);
				dxl = xl - ixl;
				al = (maxrin-inplane)%maxrin;
				ial = floor(al);
				if (ixl +1 <= nx-1 && ixl >= 0 )
					qt += (1.0-dxl)*sccf[kim]->get_value_at(ixl,ial) + dxl*sccf[kim]->get_value_at(ixl+1,ial);
				else 
					qt = -1.0e23;
				if ( six == -20 && inplane == 0 ) printf("\nkim=%d xl=%f dst=%f  coefs=%f", kim, xl, dst, (1.0-dxl)*sccf[kim]->get_value_at(ixl,ial) + dxl*sccf[kim]->get_value_at(ixl+1,ial));
				
				xl = -dst * tang + six + nxc;
				ixl = floor(xl);
				dxl = xl - ixl;
				al = inplane;
				ial = floor(al);
				if (ixl +1 <= nx-1 && ixl >= 0 )
					qu += (1.0-dxl)*sccf[kim]->get_value_at(ixl,ial) + dxl*sccf[kim]->get_value_at(ixl+1,ial);
				else 
					qu = -1.0e23;	
			}
			cout << "six= " << six << "inplane= " << (maxrin-inplane)%maxrin << "qt=  " << qt   << endl;
			if ( qt > qm ) {
				cout << " six " << six << "  qm=  " <<qm <<"  qt=  "<< qt << "  inplane=  " << (maxrin-inplane)%maxrin << endl;
				qm = qt;
				sib = six;
				bang = (maxrin-inplane)%maxrin;
			}
			if ( qu > qm ) {
				cout << " six " << six << "qm=  " <<qm <<"  qu=  "<< qu << "  inplane=  " << inplane << endl;
				qm = qu;
				sib = six;
				bang = inplane;
			}
			//cout << "bang" << bang <<endl;		
		}
		
	
	
	}
	//printf("\nqm sib bang=%f %f %f", qm, sib, bang);    
    result.push_back((double)sib); 
    result.push_back(bang);
    result.push_back(qm);
    return result;    
}



vector<float> Util::curhelixshiftali(vector<EMData*> ctx, vector<vector<float> > pcoords, int nsegms, int search_rng, int nx, int ny)
{
	int cents, six, incline, kim, ixl, sib, circnum, circ, nxc, kang;
	float q0, qt, qu, tang, dst, xl, dxl, qm, bang, circd, circR, CircR, x1, x2;
	float dst0, dst1, dst2, xmax, xmin, maxincline, Rsinang, Rcosang, RSinang, RCosang;
	vector<float> result;
	
	nxc = nx/2;
	cents = nsegms/2;
	circnum = nxc;
	qm = -1.0e23; 
	dst1 = sqrt((pcoords[cents][0] - pcoords[0][0])*(pcoords[cents][0] - pcoords[0][0]) + (pcoords[cents][1] - pcoords[0][1])*(pcoords[cents][1] - pcoords[0][1]));
	dst2 = sqrt((pcoords[cents][0] - pcoords[nsegms-1][0])*(pcoords[cents][0] - pcoords[nsegms-1][0]) + (pcoords[cents][1] - pcoords[nsegms-1][1])*(pcoords[cents][1] - pcoords[nsegms-1][1]));
	
	for ( circ= 0; circ < circnum; circ++ ) {	
		circd = 4*(circ+1.0);
		circR = circd/2.0 + ((ny*nsegms)*(ny*nsegms))/(8.0*circd);
		x1 = circR - sqrt(circR*circR-dst1*dst1);
		x2 = circR - sqrt(circR*circR-dst2*dst2);
		xmax = (x1>x2)?x1:x2;
		if ( xmax == x1 ) dst0 = dst1;
		else dst0 = dst2;
				 
		for ( six =-search_rng; six <= search_rng; six++ ) {
			if (nx-(six+nxc+xmax)-2 <= 0) continue;
			xmin = (nx-(six+nxc+xmax)-2<six+nxc-2)?(nx-(six+nxc+xmax)-2):(six+nxc-2);
			maxincline = atan2(xmin,dst0);
			kang = (int)(dst0*tan(maxincline)+0.5);
			q0 = ctx[cents]->get_value_at(six+nxc);
			for ( incline = -kang; incline <= kang; incline++ ) {
				qt = q0;
				//qu = q0;
				if ( kang > 0 ) {
					tang = tan(maxincline/kang*incline);
					Rsinang = circR*sin(maxincline/kang*incline);
					Rcosang = circR*cos(maxincline/kang*incline);}
				else {
					tang = 0.0;
					Rsinang = 0.0;
					Rcosang = 1.0;}
				for ( kim = cents+1; kim < nsegms; kim++ ) {
					dst = sqrt((pcoords[cents][0] - pcoords[kim][0])*(pcoords[cents][0] - pcoords[kim][0])+
					(pcoords[cents][1] - pcoords[kim][1])*(pcoords[cents][1] - pcoords[kim][1]));
					if ( tang >= 0 ) {                     //case 1 and 4.
						xl = Rcosang - sqrt(Rcosang*Rcosang-2*dst*Rsinang-dst*dst)+six+nxc;
						ixl = (int) xl;
						dxl = xl - ixl;
						qt += (1.0-dxl)*ctx[kim]->get_value_at(ixl) + dxl*ctx[kim]->get_value_at(ixl+1);}
						
					if ( tang < 0 ) {
						if ( -dst*tang > circR - sqrt(circR*circR-dst*dst)) {	             //## case 3
							if ( Rcosang*Rcosang+2*dst*Rsinang-dst*dst <0 )
							printf("Rcosang**2+2*dst*Rsinang-dst**2 less than 0 in curhelicalshiftali_MPI"); 
							xl = -Rcosang + sqrt(Rcosang*Rcosang+2*dst*Rsinang-dst*dst)+six+nxc;
							ixl = (int) xl;
							dxl = xl - ixl;
							if ( ixl < 0 ) printf("kim=%d, cents=%d nsegms=%d incline=%d ixl=%d", kim, cents,nsegms, incline,ixl); 
							qt += (1.0-dxl)*ctx[kim]->get_value_at(ixl) + dxl*ctx[kim]->get_value_at(ixl+1);}
						else {                                           //## case 2
							xl = Rcosang - sqrt(Rcosang*Rcosang+2*dst*Rsinang-dst*dst)+six+nxc;
							ixl = (int) xl;
							dxl = xl - ixl;
							qt += (1.0-dxl)*ctx[kim]->get_value_at(ixl) + dxl*ctx[kim]->get_value_at(ixl+1);}}
				}
				for ( kim = 0; kim < cents; kim++ ) {
					dst = sqrt((pcoords[cents][0] - pcoords[kim][0])*(pcoords[cents][0] - pcoords[kim][0])+
					(pcoords[cents][1] - pcoords[kim][1])*(pcoords[cents][1] - pcoords[kim][1]));
					if ( tang >= 0 ) {
						if ( dst*tang > circR - sqrt(circR*circR-dst*dst)) {             
							xl = -Rcosang + sqrt(Rcosang*Rcosang+2*dst*Rsinang-dst*dst)+six+nxc;
							ixl = (int) xl;
							dxl = xl - ixl;
							qt += (1.0-dxl)*ctx[kim]->get_value_at(ixl) + dxl*ctx[kim]->get_value_at(ixl+1);}
						else {                      
							xl = Rcosang - sqrt(Rcosang*Rcosang+2*dst*Rsinang-dst*dst)+six+nxc;
							ixl = (int) xl;
							dxl = xl - ixl;
							qt += (1.0-dxl)*ctx[kim]->get_value_at(ixl) + dxl*ctx[kim]->get_value_at(ixl+1);}}
					if ( tang < 0 ) {	
						xl = Rcosang - sqrt(Rcosang*Rcosang-2*dst*Rsinang-dst*dst)+six+nxc;
						ixl = (int) xl;
						dxl = xl - ixl;
						qt += (1.0-dxl)*ctx[kim]->get_value_at(ixl) + dxl*ctx[kim]->get_value_at(ixl+1);}
				}
				if ( qt > qm ) {
					qm = qt;
					sib = six;
					bang = tang;
					CircR = circR;
					RSinang = Rsinang;
					RCosang = Rcosang;}
			}
		
		}
	}	
	//printf("\nqm sib bang=%f %f %f", qm, sib, bang);    
    result.push_back((float)sib); 
    result.push_back(bang);
    result.push_back(CircR);
    result.push_back(RSinang);
    result.push_back(RCosang);
    result.push_back(qm);
    return result;    
}



//----------------------------------------------------------------------------------------------------------

Dict Util::im_diff(EMData* V1, EMData* V2, EMData* mask)
{
	ENTERFUNC;

	if (!EMUtil::is_same_size(V1, V2)) {
		LOGERR("images not same size");
		throw ImageFormatException( "images not same size");
	}

	size_t nx = V1->get_xsize();
	size_t ny = V1->get_ysize();
	size_t nz = V1->get_zsize();
	size_t size = (size_t)nx*ny*nz;

	EMData *BD = new EMData();
 	BD->set_size(nx, ny, nz);

	float *params = new float[2];

	float *V1ptr, *V2ptr, *MASKptr, *BDptr, A, B;
	long double S1=0.L,S2=0.L,S3=0.L,S4=0.L;
	int nvox = 0L;

   	V1ptr = V1->get_data();
	V2ptr = V2->get_data();
	BDptr = BD->get_data();


	if(!mask){
		EMData * Mask = new EMData();
		Mask->set_size(nx,ny,nz);
		Mask->to_one();
		MASKptr = Mask->get_data();
	} else {
		if (!EMUtil::is_same_size(V1, mask)) {
			LOGERR("mask not same size");
			throw ImageFormatException( "mask not same size");
		}

		MASKptr = mask->get_data();
	}



//	 calculation of S1,S2,S3,S3,nvox

	for (size_t i = 0L;i < size; i++) {
	      if (MASKptr[i]>0.5f) {
	       S1 += V1ptr[i]*V2ptr[i];
	       S2 += V1ptr[i]*V1ptr[i];
	       S3 += V2ptr[i];
	       S4 += V1ptr[i];
	       nvox ++;
	      }
	}

	if ((nvox*S1 - S3*S4) == 0. || (nvox*S2 - S4*S4) == 0) {
		A =1.0f ;
	} else {
		A = static_cast<float>( (nvox*S1 - S3*S4)/(nvox*S2 - S4*S4) );
	}
	B = static_cast<float> (A*S4  -  S3)/nvox;

	// calculation of the difference image

	for (size_t i = 0L;i < size; i++) {
	     if (MASKptr[i]>0.5f) {
	       BDptr[i] = A*V1ptr[i] -  B  - V2ptr[i];
	     }  else  {
               BDptr[i] = 0.f;
	     }
	}

	BD->update();

	params[0] = A;
	params[1] = B;

	Dict BDnParams;
	BDnParams["imdiff"] = BD;
	BDnParams["A"] = params[0];
	BDnParams["B"] = params[1];

	EXITFUNC;
	return BDnParams;
 }

//----------------------------------------------------------------------------------------------------------



EMData *Util::TwoDTestFunc(int Size, float p, float q,  float a, float b, int flag, float alphaDeg) //PRB
{
	ENTERFUNC;
    	int Mid= (Size+1)/2;

    	if (flag==0) { // This is the real function
		EMData* ImBW = new EMData();
	   	ImBW->set_size(Size,Size,1);
	   	ImBW->to_zero();

	   	float tempIm;
	   	float x,y;

	   	for (int ix=(1-Mid);  ix<Mid; ix++){
	        	for (int iy=(1-Mid);  iy<Mid; iy++){
		  		x = (float)ix;
		  		y = (float)iy;
	       	  	tempIm= static_cast<float>( (1/(2*M_PI)) * cos(p*x)* cos(q*y) * exp(-.5*x*x/(a*a))* exp(-.5*y*y/(b*b)) );
		  		(*ImBW)(ix+Mid-1,iy+Mid-1) = tempIm * exp(.5f*p*p*a*a)* exp(.5f*q*q*b*b);
	   		}
	   	}
	   	ImBW->update();
	   	ImBW->set_complex(false);
	   	ImBW->set_ri(true);

	   	return ImBW;
   	}
	else if (flag==1) {  // This is the Fourier Transform
		EMData* ImBWFFT = new EMData();
	   	ImBWFFT ->set_size(2*Size,Size,1);
	   	ImBWFFT ->to_zero();

	   	float r,s;

	   	for (int ir=(1-Mid);  ir<Mid; ir++){
	        	for (int is=(1-Mid);  is<Mid; is++){
		   		r = (float)ir;
		   		s = (float)is;
	       	   	(*ImBWFFT)(2*(ir+Mid-1),is+Mid-1)= cosh(p*r*a*a) * cosh(q*s*b*b) *
		            	exp(-.5f*r*r*a*a)* exp(-.5f*s*s*b*b);
	   		}
	   	}
	   	ImBWFFT->update();
	   	ImBWFFT->set_complex(true);
	   	ImBWFFT->set_ri(true);
	   	ImBWFFT->set_shuffled(true);
	   	ImBWFFT->set_fftodd(true);

	   	return ImBWFFT;
   	}
   	else if (flag==2 || flag==3) { //   This is the projection in Real Space
		float alpha = static_cast<float>( alphaDeg*M_PI/180.0 );
		float C=cos(alpha);
		float S=sin(alpha);
		float D= sqrt(S*S*b*b + C*C*a*a);
		//float D2 = D*D;   PAP - to get rid of warning

		float P = p * C *a*a/D ;
		float Q = q * S *b*b/D ;

		if (flag==2) {
			EMData* pofalpha = new EMData();
			pofalpha ->set_size(Size,1,1);
			pofalpha ->to_zero();

			float Norm0 =  D*(float)sqrt(2*pi);
			float Norm1 =  exp( .5f*(P+Q)*(P+Q)) / Norm0 ;
			float Norm2 =  exp( .5f*(P-Q)*(P-Q)) / Norm0 ;
			float sD;

			for (int is=(1-Mid);  is<Mid; is++){
				sD = is/D ;
				(*pofalpha)(is+Mid-1) =  Norm1 * exp(-.5f*sD*sD)*cos(sD*(P+Q))
                         + Norm2 * exp(-.5f*sD*sD)*cos(sD*(P-Q));
			}
			pofalpha-> update();
			pofalpha-> set_complex(false);
			pofalpha-> set_ri(true);

			return pofalpha;
		}
		if (flag==3) { // This is the projection in Fourier Space
			float vD;

			EMData* pofalphak = new EMData();
			pofalphak ->set_size(2*Size,1,1);
			pofalphak ->to_zero();

			for (int iv=(1-Mid);  iv<Mid; iv++){
				vD = iv*D ;
		 		(*pofalphak)(2*(iv+Mid-1)) =  exp(-.5f*vD*vD)*(cosh(vD*(P+Q)) + cosh(vD*(P-Q)) );
			}
			pofalphak-> update();
			pofalphak-> set_complex(false);
			pofalphak-> set_ri(true);

			return pofalphak;
   		}
   	}
    	else if (flag==4) {
		cout <<" FH under construction";
   		EMData* OutFT= TwoDTestFunc(Size, p, q, a, b, 1);
   		EMData* TryFH= OutFT -> real2FH(4.0);
   		return TryFH;
   	} else {
    		cout <<" flag must be 0,1,2,3, or 4";
    	}

	EXITFUNC;
	return 0;
}


void Util::spline_mat(float *x, float *y, int n,  float *xq, float *yq, int m) //PRB
{

	float x0= x[0];
	float x1= x[1];
	float x2= x[2];
	float y0= y[0];
	float y1= y[1];
	float y2= y[2];
	float yp1 =  (y1-y0)/(x1-x0) +  (y2-y0)/(x2-x0) - (y2-y1)/(x2-x1)  ;
	float xn  = x[n];
	float xnm1= x[n-1];
	float xnm2= x[n-2];
	float yn  = y[n];
	float ynm1= y[n-1];
	float ynm2= y[n-2];
	float ypn=  (yn-ynm1)/(xn-xnm1) +  (yn-ynm2)/(xn-xnm2) - (ynm1-ynm2)/(xnm1-xnm2) ;
	float *y2d = new float[n];
	Util::spline(x,y,n,yp1,ypn,y2d);
	Util::splint(x,y,y2d,n,xq,yq,m); //PRB
	delete [] y2d;
	return;
}


void Util::spline(float *x, float *y, int n, float yp1, float ypn, float *y2) //PRB
{
	int i,k;
	float p, qn, sig, un, *u;
	u = new float[n-1];

	if (yp1 > .99e30){
		y2[0]=u[0]=0.0;
	} else {
		y2[0]=-.5f;
		u[0] =(3.0f/ (x[1] -x[0]))*( (y[1]-y[0])/(x[1]-x[0]) -yp1);
	}

	for (i=1; i < n-1; i++) {
		sig= (x[i] - x[i-1])/(x[i+1] - x[i-1]);
		p = sig*y2[i-1] + 2.0f;
		y2[i]  = (sig-1.0f)/p;
		u[i] = (y[i+1] - y[i] )/(x[i+1]-x[i] ) -  (y[i] - y[i-1] )/(x[i] -x[i-1]);
		u[i] = (6.0f*u[i]/ (x[i+1]-x[i-1]) - sig*u[i-1])/p;
	}

	if (ypn>.99e30){
		qn=0; un=0;
	} else {
		qn= .5f;
		un= (3.0f/(x[n-1] -x[n-2])) * (ypn -  (y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1]= (un - qn*u[n-2])/(qn*y2[n-2]+1.0f);
	for (k=n-2; k>=0; k--){
		y2[k]=y2[k]*y2[k+1]+u[k];
	}
	delete [] u;
}


/***************************************************************************
Description:
  Compute the value at u of ith B_spline Base Function 
  Spline  Knots: [u0,u1,..., u_(m+2k)]. 

Arguments:
  MM = m+k. the total number of B_spline Base Functions at the above interval.   
  i = the ith B_spline Base Function.
  k = The degree of B_spline Base Function.
****************************************************************************/
float Util::bsplineBase(float u)
{
  float ww1, ww2, ww3, sign;
  float value = 0.0;
  
    ww1 = fabs(u);
   if (ww1 >= 2.0) {
        return value;
   }
   ww2 = ww1*ww1;
   ww3 = ww2*ww1;

 
   if (ww1 >=  0.0 && ww1 < 1.0) {
        ww2 = ww1*ww1;
        ww3 = ww2*ww1;
    value = 2.0/3.0 - ww2 + 0.5*ww3;
    return value;
   }

   if (ww1 >=  1.0 && ww1 <= 2.0) {
        ww1 = 2.0 - ww1;
        ww2 = ww1*ww1;
        ww3 = ww2*ww1;
        value = 1.0/6.0*ww3;

   }
	return value;

}

float Util::bsplineBasedu(float u)
{  
float ww1, ww2, ww3, sign;
float values = 0;  
/*
  ww1 = fabs(u);
  if (ww1 >= 2.0) {
    return values;
  }
  ww2 = ww1*ww1;
  
  if (u >= 0.0) {
     sign = 1.0;
  }  else {
     sign = -1.0;
  }

  if (ww1 >=  0.0 && ww1 < 1.0) {
     ww2 = ww1*ww1;
     values = -2.0*u + 1.5*sign*ww2;
     return values;
  }
   
  if (ww1 >=  1.0 && ww1 <= 2.0) {
     ww1 = 2.0 - ww1;
     ww2 = ww1*ww1;
     values = -0.5*sign*ww2;
  }
  return values;
*/
  
  ww1 = fabs(u);
   if (ww1 >= 1.0) {
     return values;
   }

   if ( u >=  -1.0 && u < 0.0) {
     values = u+1.0;
   }

   if ( u >= 0.0 && u <= 1.0 ) {
     values = 1.0 - u;
   }
   
return values;
}


void Util::convertTocubicbsplineCoeffs(vector<float> s, int DataLength, float EPSILON)
//              float   *s,             /* input samples --> output coefficients */
//              int     DataLength,     /* number of samples or coefficients     */
//              float   EPSILON         /* admissible relative error             */

{
int   i, n, ni, ni1, K;
float sum, z1, w1, w2;

n = DataLength + 1;
z1 = sqrt(3.0) - 2.0;
K = log(EPSILON)/log(fabs(z1));
//printf("K = %i\n", K);

// compute initial value s(0)
sum = 0.0;
w2 = pow(z1, 2*n);
if (n < K) {
   for (i = 1; i < n; i++){
      w1 = pow(z1, i);
      sum = sum + s[i-1]*(w1 - w2/w1);
   }
} else {
   for (i = 1; i < n; i++){
      sum = (sum + s[n- i-1])*z1;
   }
}
sum = -sum/(1.0 - w2);

// compute c^{+}
n = DataLength;
s[0]  = s[0] + z1*sum;
for (i = 1; i < n; i++) {
   s[i]  = s[i] + z1*s[i-1];
   //printf("cp[%i] = %e, %f \n", i, cp[i], z1);
}

// compute c^- 
s[n-1] = -z1*s[n-1];
for (i = 1; i < n; i++) {
   ni = n - i;
   ni1 = ni - 1;
   s[ni1]  = z1*(s[ni] - s[ni1]);
}

for (i = 0; i < n; i++) {
   s[i]  = 6.0*s[i];
}

}


void Util::splint( float *xa, float *ya, float *y2a, int n,  float *xq, float *yq, int m) //PRB
{
	int klo, khi, k;
	float h, b, a;

//	klo=0; // can try to put here
	for (int j=0; j<m;j++){
		klo=0;
		khi=n-1;
		while (khi-klo >1) {
			k=(khi+klo) >>1;
			if  (xa[k]>xq[j]){ khi=k;}
			else { klo=k;}
		}
		h=xa[khi]- xa[klo];
		if (h==0.0) printf("Bad XA input to routine SPLINT \n");
		a =(xa[khi]-xq[j])/h;
		b=(xq[j]-xa[klo])/h;
		yq[j]=a*ya[klo] + b*ya[khi]
			+ ((a*a*a-a)*y2a[klo]
			     +(b*b*b-b)*y2a[khi]) *(h*h)/6.0f;
	}
//	printf("h=%f, a = %f, b=%f, ya[klo]=%f, ya[khi]=%f , yq=%f\n",h, a, b, ya[klo], ya[khi],yq[0]);
}


void Util::Radialize(int *PermMatTr, float *kValsSorted,   // PRB
               float *weightofkValsSorted, int Size, int *SizeReturned)
{
	int iMax = (int) floor( (Size-1.0)/2 +.01);
	int CountMax = (iMax+2)*(iMax+1)/2;
	int Count=-1;
	float *kVals     = new float[CountMax];
	float *weightMat = new float[CountMax];
	int *PermMat     = new   int[CountMax];
	SizeReturned[0] = CountMax;

//	printf("Aa \n");	fflush(stdout);
	for (int jkx=0; jkx< iMax+1; jkx++) {
		for (int jky=0; jky< jkx+1; jky++) {
			Count++;
			kVals[Count] = sqrtf((float) (jkx*jkx +jky*jky));
			weightMat[Count]=  1.0;
			if (jkx!=0)  { weightMat[Count] *=2;}
			if (jky!=0)  { weightMat[Count] *=2;}
			if (jkx!=jky){ weightMat[Count] *=2;}
			PermMat[Count]=Count+1;
		}
	}

	int lkVals = Count+1;
//	printf("Cc \n");fflush(stdout);

	sort_mat(&kVals[0],&kVals[Count],
	     &PermMat[0],  &PermMat[Count]);  //PermMat is
				//also returned as well as kValsSorted
	fflush(stdout);

	int newInd;

        for (int iP=0; iP < lkVals ; iP++ ) {
		newInd =  PermMat[iP];
		PermMatTr[newInd-1] = iP+1;
	}

//	printf("Ee \n"); fflush(stdout);

	int CountA=-1;
	int CountB=-1;

	while (CountB< (CountMax-1)) {
		CountA++;
		CountB++;
//		printf("CountA=%d ; CountB=%d \n", CountA,CountB);fflush(stdout);
		kValsSorted[CountA] = kVals[CountB] ;
		if (CountB<(CountMax-1) ) {
			while (fabs(kVals[CountB] -kVals[CountB+1])<.0000001  ) {
				SizeReturned[0]--;
				for (int iP=0; iP < lkVals; iP++){
//					printf("iP=%d \n", iP);fflush(stdout);
					if  (PermMatTr[iP]>CountA+1) {
						PermMatTr[iP]--;
		    			}
		 		}
				CountB++;
	    		}
		}
	}


	for (int CountD=0; CountD < CountMax; CountD++) {
	    newInd = PermMatTr[CountD];
	    weightofkValsSorted[newInd-1] += weightMat[CountD];
        }

}


vector<float>
Util::even_angles(float delta, float t1, float t2, float p1, float p2)
{
	vector<float> angles;
	float psi = 0.0;
	if ((0.0 == t1 && 0.0 == t2)||(t1 >= t2)) {
		t1 =  0.0f;
		t2 = 90.0f;
	}
	if ((0.0 == p1 && 0.0 == p2)||(p1 >= p2)) {
		p1 =   0.0f;
		p2 = 359.9f;
	}
	bool skip = ((t1 < 90.0) && (90.0 == t2) && (0.0 == p1) && (p2 > 180.0));
	for (float theta = t1; theta <= t2; theta += delta) {
		float detphi;
		int lt;
		if ((0.0 == theta)||(180.0 == theta)) {
			detphi = 360.0f;
			lt = 1;
		} else {
			detphi = delta/sin(theta*static_cast<float>(dgr_to_rad));
			lt = int((p2 - p1)/detphi)-1;
			if (lt < 1) lt = 1;
			detphi = (p2 - p1)/lt;
		}
		for (int i = 0; i < lt; i++) {
			float phi = p1 + i*detphi;
			if (skip&&(90.0 == theta)&&(phi > 180.0)) continue;
			angles.push_back(phi);
			angles.push_back(theta);
			angles.push_back(psi);
		}
	}
	return angles;
}


#define  fdata(i,j)      fdata[ i-1 + (j-1)*nxdata ]
/*float Util::quadri(float xx, float yy, int nxdata, int nydata, float* fdata)
{

//  purpose: quadratic interpolation
//
//  parameters:       xx,yy treated as circularly closed.
//                    fdata - image 1..nxdata, 1..nydata
//
//                    f3    fc       f0, f1, f2, f3 are the values
//                     +             at the grid points.  x is the
//                     + x           point at which the function
//              f2++++f0++++f1       is to be estimated. (it need
//                     +             not be in the first quadrant).
//                     +             fc - the outer corner point
//                    f4             nearest x.
c
//                                   f0 is the value of the fdata at
//                                   fdata(i,j), it is the interior mesh
//                                   point nearest  x.
//                                   the coordinates of f0 are (x0,y0),
//                                   the coordinates of f1 are (xb,y0),
//                                   the coordinates of f2 are (xa,y0),
//                                   the coordinates of f3 are (x0,yb),
//                                   the coordinates of f4 are (x0,ya),
//                                   the coordinates of fc are (xc,yc),
c
//                   o               hxa, hxb are the mesh spacings
//                   +               in the x-direction to the left
//                  hyb              and right of the center point.
//                   +
//            ++hxa++o++hxb++o       hyb, hya are the mesh spacings
//                   +               in the y-direction.
//                  hya
//                   +               hxc equals either  hxb  or  hxa
//                   o               depending on where the corner
//                                   point is located.
c
//                                   construct the interpolant
//                                   f = f0 + c1*(x-x0) +
//                                       c2*(x-x0)*(x-x1) +
//                                       c3*(y-y0) + c4*(y-y0)*(y-y1)
//                                       + c5*(x-x0)*(y-y0)
//
//

    float x, y, dx0, dy0, f0, c1, c2, c3, c4, c5, dxb, dyb;
    float quadri;
    int   i, j, ip1, im1, jp1, jm1, ic, jc, hxc, hyc;

    x = xx;
    y = yy;

    // circular closure
	while ( x < 1.0 ) x += nxdata;
	while ( x >= (float)(nxdata+1) )  x -= nxdata;
	while ( y < 1.0 ) y += nydata;
	while ( y >= (float)(nydata+1) )  y -= nydata;


    i   = (int) x;
    j   = (int) y;

    dx0 = x - i;
    dy0 = y - j;

    ip1 = i + 1;
    im1 = i - 1;
    jp1 = j + 1;
    jm1 = j - 1;

    if (ip1 > nxdata) ip1 = ip1 - nxdata;
    if (im1 < 1)      im1 = im1 + nxdata;
    if (jp1 > nydata) jp1 = jp1 - nydata;
    if (jm1 < 1)      jm1 = jm1 + nydata;

    f0  = fdata(i,j);
    c1  = fdata(ip1,j) - f0;
    c2  = (c1 - f0 + fdata(im1,j)) * 0.5;
    c3  = fdata(i,jp1) - f0;
    c4  = (c3 - f0 + fdata(i,jm1)) * 0.5;

    dxb = dx0 - 1;
    dyb = dy0 - 1;

    // hxc & hyc are either 1 or -1
    if (dx0 >= 0) { hxc = 1; } else { hxc = -1; }
    if (dy0 >= 0) { hyc = 1; } else { hyc = -1; }

    ic  = i + hxc;
    jc  = j + hyc;

    if (ic > nxdata) { ic = ic - nxdata; }  else if (ic < 1) { ic = ic + nxdata; }
    if (jc > nydata) { jc = jc - nydata; } else if (jc < 1) { jc = jc + nydata; }

    c5  =  ( (fdata(ic,jc) - f0 - hxc * c1 - (hxc * (hxc - 1.0)) * c2
            - hyc * c3 - (hyc * (hyc - 1.0)) * c4) * (hxc * hyc));

    quadri = f0 + dx0 * (c1 + dxb * c2 + dy0 * c5) + dy0 * (c3 + dyb * c4);

    return quadri;
}*/
float Util::quadri(float xx, float yy, int nxdata, int nydata, float* fdata)
{
//  purpose: quadratic interpolation
//  Optimized for speed, circular closer removed, checking of ranges removed
	float  x, y, dx0, dy0, f0, c1, c2, c3, c4, c5, dxb, dyb;
	float  quadri;
	int    i, j, ip1, im1, jp1, jm1, ic, jc, hxc, hyc;

	x = xx;
	y = yy;

	//     any xx and yy
	while ( x < 1.0 )                 x += nxdata;
	while ( x >= (float)(nxdata+1) )  x -= nxdata;
	while ( y < 1.0 )                 y += nydata;
	while ( y >= (float)(nydata+1) )  y -= nydata;

	i   = (int) x;
	j   = (int) y;

	dx0 = x - i;
	dy0 = y - j;

	ip1 = i + 1;
	im1 = i - 1;
	jp1 = j + 1;
	jm1 = j - 1;

	if (ip1 > nxdata) ip1 -= nxdata;
	if (im1 < 1)	  im1 += nxdata;
	if (jp1 > nydata) jp1 -= nydata;
	if (jm1 < 1)	  jm1 += nydata;

	f0  = fdata(i,j);
	c1  = fdata(ip1,j) - f0;
	c2  = (c1 - f0 + fdata(im1,j)) * 0.5f;
	c3  = fdata(i,jp1) - f0;
	c4  = (c3 - f0 + fdata(i,jm1)) * 0.5f;

	dxb = dx0 - 1;
	dyb = dy0 - 1;

	// hxc & hyc are either 1 or -1
	if (dx0 >= 0) hxc = 1; else hxc = -1;
	if (dy0 >= 0) hyc = 1; else hyc = -1;

	ic  = i + hxc;
	jc  = j + hyc;

	if (ic > nxdata) ic -= nxdata;  else if (ic < 1) ic += nxdata;
	if (jc > nydata) jc -= nydata;  else if (jc < 1) jc += nydata;

	c5  =  ( (fdata(ic,jc) - f0 - hxc * c1 - (hxc * (hxc - 1.0f)) * c2
		- hyc * c3 - (hyc * (hyc - 1.0f)) * c4) * (hxc * hyc));


	quadri = f0 + dx0 * (c1 + dxb * c2 + dy0 * c5) + dy0 * (c3 + dyb * c4);

	return quadri;
}

#undef fdata

#define  fdata(i,j)      fdata[ i-1 + (j-1)*nxdata ]
float Util::quadri_background(float xx, float yy, int nxdata, int nydata, float* fdata, int xnew, int ynew)
{
//  purpose: quadratic interpolation
//  Optimized for speed, circular closer removed, checking of ranges removed
	float  x, y, dx0, dy0, f0, c1, c2, c3, c4, c5, dxb, dyb;
	float  quadri;
	int    i, j, ip1, im1, jp1, jm1, ic, jc, hxc, hyc;

	x = xx;
	y = yy;

	// wrap around is not done circulantly; if (x,y) is not in the image, then x = xnew and y = ynew
	if ( (x < 1.0) || ( x >= (float)(nxdata+1) ) || ( y < 1.0 ) || ( y >= (float)(nydata+1) )){
	      x = (float)xnew;
		  y = (float)ynew;
     }


	i   = (int) x;
	j   = (int) y;

	dx0 = x - i;
	dy0 = y - j;

	ip1 = i + 1;
	im1 = i - 1;
	jp1 = j + 1;
	jm1 = j - 1;

	if (ip1 > nxdata) ip1 -= nxdata;
	if (im1 < 1)	  im1 += nxdata;
	if (jp1 > nydata) jp1 -= nydata;
	if (jm1 < 1)	  jm1 += nydata;

	f0  = fdata(i,j);
	c1  = fdata(ip1,j) - f0;
	c2  = (c1 - f0 + fdata(im1,j)) * 0.5f;
	c3  = fdata(i,jp1) - f0;
	c4  = (c3 - f0 + fdata(i,jm1)) * 0.5f;

	dxb = dx0 - 1;
	dyb = dy0 - 1;

	// hxc & hyc are either 1 or -1
	if (dx0 >= 0) hxc = 1; else hxc = -1;
	if (dy0 >= 0) hyc = 1; else hyc = -1;

	ic  = i + hxc;
	jc  = j + hyc;

	if (ic > nxdata) ic -= nxdata;  else if (ic < 1) ic += nxdata;
	if (jc > nydata) jc -= nydata;  else if (jc < 1) jc += nydata;

	c5  =  ( (fdata(ic,jc) - f0 - hxc * c1 - (hxc * (hxc - 1.0f)) * c2
		- hyc * c3 - (hyc * (hyc - 1.0f)) * c4) * (hxc * hyc));


	quadri = f0 + dx0 * (c1 + dxb * c2 + dy0 * c5) + dy0 * (c3 + dyb * c4);

	return quadri;
}

#undef fdata


float  Util::get_pixel_conv_new(int nx, int ny, int nz, float delx, float dely, float delz, float* data, Util::KaiserBessel& kb) {
	int K = kb.get_window_size();
	int kbmin = -K/2;
	int kbmax = -kbmin;
	int kbc = kbmax+1;

	float pixel = 0.0f;
	float w = 0.0f;

	delx = restrict1(delx, nx);
	int inxold = int(round(delx));
	if ( ny < 2 ) {  //1D
		float tablex1 = kb.i0win_tab(delx-inxold+3);
		float tablex2 = kb.i0win_tab(delx-inxold+2);
		float tablex3 = kb.i0win_tab(delx-inxold+1);
		float tablex4 = kb.i0win_tab(delx-inxold);
		float tablex5 = kb.i0win_tab(delx-inxold-1);
		float tablex6 = kb.i0win_tab(delx-inxold-2);
		float tablex7 = kb.i0win_tab(delx-inxold-3);

		int x1, x2, x3, x4, x5, x6, x7;

	 	if ( inxold <= kbc || inxold >=nx-kbc-2 )  {
			x1 = (inxold-3+nx)%nx;
			x2 = (inxold-2+nx)%nx;
			x3 = (inxold-1+nx)%nx;
			x4 = (inxold  +nx)%nx;
			x5 = (inxold+1+nx)%nx;
			x6 = (inxold+2+nx)%nx;
			x7 = (inxold+3+nx)%nx;
	 	} else {
			x1 = inxold-3;
			x2 = inxold-2;
			x3 = inxold-1;
			x4 = inxold;
			x5 = inxold+1;
			x6 = inxold+2;
			x7 = inxold+3;
	 	}

		pixel = data[x1]*tablex1 + data[x2]*tablex2 + data[x3]*tablex3 +
				data[x4]*tablex4 + data[x5]*tablex5 + data[x6]*tablex6 +
				data[x7]*tablex7 ;

		w = tablex1+tablex2+tablex3+tablex4+tablex5+tablex6+tablex7;
	} else if ( nz < 2 ) {  // 2D
		dely = restrict1(dely, ny);
		int inyold = int(round(dely));
		float tablex1 = kb.i0win_tab(delx-inxold+3);
		float tablex2 = kb.i0win_tab(delx-inxold+2);
		float tablex3 = kb.i0win_tab(delx-inxold+1);
		float tablex4 = kb.i0win_tab(delx-inxold);
		float tablex5 = kb.i0win_tab(delx-inxold-1);
		float tablex6 = kb.i0win_tab(delx-inxold-2);
		float tablex7 = kb.i0win_tab(delx-inxold-3);

		float tabley1 = kb.i0win_tab(dely-inyold+3);
		float tabley2 = kb.i0win_tab(dely-inyold+2);
		float tabley3 = kb.i0win_tab(dely-inyold+1);
		float tabley4 = kb.i0win_tab(dely-inyold);
		float tabley5 = kb.i0win_tab(dely-inyold-1);
		float tabley6 = kb.i0win_tab(dely-inyold-2);
		float tabley7 = kb.i0win_tab(dely-inyold-3);

		int x1, x2, x3, x4, x5, x6, x7, y1, y2, y3, y4, y5, y6, y7;

	 	if ( inxold <= kbc || inxold >=nx-kbc-2 || inyold <= kbc || inyold >=ny-kbc-2 )  {
			x1 = (inxold-3+nx)%nx;
			x2 = (inxold-2+nx)%nx;
			x3 = (inxold-1+nx)%nx;
			x4 = (inxold  +nx)%nx;
			x5 = (inxold+1+nx)%nx;
			x6 = (inxold+2+nx)%nx;
			x7 = (inxold+3+nx)%nx;

			y1 = ((inyold-3+ny)%ny)*nx;
			y2 = ((inyold-2+ny)%ny)*nx;
			y3 = ((inyold-1+ny)%ny)*nx;
			y4 = ((inyold  +ny)%ny)*nx;
			y5 = ((inyold+1+ny)%ny)*nx;
			y6 = ((inyold+2+ny)%ny)*nx;
			y7 = ((inyold+3+ny)%ny)*nx;
		} else {
			x1 = inxold-3;
			x2 = inxold-2;
			x3 = inxold-1;
			x4 = inxold;
			x5 = inxold+1;
			x6 = inxold+2;
			x7 = inxold+3;

			y1 = (inyold-3)*nx;
			y2 = (inyold-2)*nx;
			y3 = (inyold-1)*nx;
			y4 = inyold*nx;
			y5 = (inyold+1)*nx;
			y6 = (inyold+2)*nx;
			y7 = (inyold+3)*nx;
		}

		pixel    = ( data[x1+y1]*tablex1 + data[x2+y1]*tablex2 + data[x3+y1]*tablex3 +
		             data[x4+y1]*tablex4 + data[x5+y1]*tablex5 + data[x6+y1]*tablex6 +
					 data[x7+y1]*tablex7 ) * tabley1 +
				   ( data[x1+y2]*tablex1 + data[x2+y2]*tablex2 + data[x3+y2]*tablex3 +
					 data[x4+y2]*tablex4 + data[x5+y2]*tablex5 + data[x6+y2]*tablex6 +
					 data[x7+y2]*tablex7 ) * tabley2 +
				   ( data[x1+y3]*tablex1 + data[x2+y3]*tablex2 + data[x3+y3]*tablex3 +
					 data[x4+y3]*tablex4 + data[x5+y3]*tablex5 + data[x6+y3]*tablex6 +
					 data[x7+y3]*tablex7 ) * tabley3 +
				   ( data[x1+y4]*tablex1 + data[x2+y4]*tablex2 + data[x3+y4]*tablex3 +
					 data[x4+y4]*tablex4 + data[x5+y4]*tablex5 + data[x6+y4]*tablex6 +
					 data[x7+y4]*tablex7 ) * tabley4 +
				   ( data[x1+y5]*tablex1 + data[x2+y5]*tablex2 + data[x3+y5]*tablex3 +
					 data[x4+y5]*tablex4 + data[x5+y5]*tablex5 + data[x6+y5]*tablex6 +
					 data[x7+y5]*tablex7 ) * tabley5 +
				   ( data[x1+y6]*tablex1 + data[x2+y6]*tablex2 + data[x3+y6]*tablex3 +
					 data[x4+y6]*tablex4 + data[x5+y6]*tablex5 + data[x6+y6]*tablex6 +
					 data[x7+y6]*tablex7 ) * tabley6 +
				   ( data[x1+y7]*tablex1 + data[x2+y7]*tablex2 + data[x3+y7]*tablex3 +
					 data[x4+y7]*tablex4 + data[x5+y7]*tablex5 + data[x6+y7]*tablex6 +
					 data[x7+y7]*tablex7 ) * tabley7;

		w = (tablex1+tablex2+tablex3+tablex4+tablex5+tablex6+tablex7) *
		    (tabley1+tabley2+tabley3+tabley4+tabley5+tabley6+tabley7);
	} else {  //  3D
		dely = restrict1(dely, ny);
		int inyold = int(Util::round(dely));
		delz = restrict1(delz, nz);
		int inzold = int(Util::round(delz));

		float tablex1 = kb.i0win_tab(delx-inxold+3);
		float tablex2 = kb.i0win_tab(delx-inxold+2);
		float tablex3 = kb.i0win_tab(delx-inxold+1);
		float tablex4 = kb.i0win_tab(delx-inxold);
		float tablex5 = kb.i0win_tab(delx-inxold-1);
		float tablex6 = kb.i0win_tab(delx-inxold-2);
		float tablex7 = kb.i0win_tab(delx-inxold-3);

		float tabley1 = kb.i0win_tab(dely-inyold+3);
		float tabley2 = kb.i0win_tab(dely-inyold+2);
		float tabley3 = kb.i0win_tab(dely-inyold+1);
		float tabley4 = kb.i0win_tab(dely-inyold);
		float tabley5 = kb.i0win_tab(dely-inyold-1);
		float tabley6 = kb.i0win_tab(dely-inyold-2);
		float tabley7 = kb.i0win_tab(dely-inyold-3);

		float tablez1 = kb.i0win_tab(delz-inzold+3);
		float tablez2 = kb.i0win_tab(delz-inzold+2);
		float tablez3 = kb.i0win_tab(delz-inzold+1);
		float tablez4 = kb.i0win_tab(delz-inzold);
		float tablez5 = kb.i0win_tab(delz-inzold-1);
		float tablez6 = kb.i0win_tab(delz-inzold-2);
		float tablez7 = kb.i0win_tab(delz-inzold-3);

		int x1, x2, x3, x4, x5, x6, x7, y1, y2, y3, y4, y5, y6, y7, z1, z2, z3, z4, z5, z6, z7;

		if ( inxold <= kbc || inxold >=nx-kbc-2 || inyold <= kbc || inyold >=ny-kbc-2 || inzold <= kbc || inzold >= nz-kbc-2 )  {
			x1 = (inxold-3+nx)%nx;
			x2 = (inxold-2+nx)%nx;
			x3 = (inxold-1+nx)%nx;
			x4 = (inxold  +nx)%nx;
			x5 = (inxold+1+nx)%nx;
			x6 = (inxold+2+nx)%nx;
			x7 = (inxold+3+nx)%nx;

			y1 = ((inyold-3+ny)%ny)*nx;
			y2 = ((inyold-2+ny)%ny)*nx;
			y3 = ((inyold-1+ny)%ny)*nx;
			y4 = ((inyold  +ny)%ny)*nx;
			y5 = ((inyold+1+ny)%ny)*nx;
			y6 = ((inyold+2+ny)%ny)*nx;
			y7 = ((inyold+3+ny)%ny)*nx;

			z1 = ((inzold-3+nz)%nz)*nx*ny;
			z2 = ((inzold-2+nz)%nz)*nx*ny;
			z3 = ((inzold-1+nz)%nz)*nx*ny;
			z4 = ((inzold  +nz)%nz)*nx*ny;
			z5 = ((inzold+1+nz)%nz)*nx*ny;
			z6 = ((inzold+2+nz)%nz)*nx*ny;
			z7 = ((inzold+3+nz)%nz)*nx*ny;
		} else {
			x1 = inxold-3;
			x2 = inxold-2;
			x3 = inxold-1;
			x4 = inxold;
			x5 = inxold+1;
			x6 = inxold+2;
			x7 = inxold+3;

			y1 = (inyold-3)*nx;
			y2 = (inyold-2)*nx;
			y3 = (inyold-1)*nx;
			y4 = inyold*nx;
			y5 = (inyold+1)*nx;
			y6 = (inyold+2)*nx;
			y7 = (inyold+3)*nx;

			z1 = (inzold-3)*nx*ny;
			z2 = (inzold-2)*nx*ny;
			z3 = (inzold-1)*nx*ny;
			z4 = inzold*nx*ny;
			z5 = (inzold+1)*nx*ny;
			z6 = (inzold+2)*nx*ny;
			z7 = (inzold+3)*nx*ny;
	 	}

		pixel  = ( ( data[x1+y1+z1]*tablex1 + data[x2+y1+z1]*tablex2 + data[x3+y1+z1]*tablex3 +
		             data[x4+y1+z1]*tablex4 + data[x5+y1+z1]*tablex5 + data[x6+y1+z1]*tablex6 +
					 data[x7+y1+z1]*tablex7 ) * tabley1 +
				   ( data[x1+y2+z1]*tablex1 + data[x2+y2+z1]*tablex2 + data[x3+y2+z1]*tablex3 +
					 data[x4+y2+z1]*tablex4 + data[x5+y2+z1]*tablex5 + data[x6+y2+z1]*tablex6 +
					 data[x7+y2+z1]*tablex7 ) * tabley2 +
				   ( data[x1+y3+z1]*tablex1 + data[x2+y3+z1]*tablex2 + data[x3+y3+z1]*tablex3 +
					 data[x4+y3+z1]*tablex4 + data[x5+y3+z1]*tablex5 + data[x6+y3+z1]*tablex6 +
					 data[x7+y3+z1]*tablex7 ) * tabley3 +
				   ( data[x1+y4+z1]*tablex1 + data[x2+y4+z1]*tablex2 + data[x3+y4+z1]*tablex3 +
					 data[x4+y4+z1]*tablex4 + data[x5+y4+z1]*tablex5 + data[x6+y4+z1]*tablex6 +
					 data[x7+y4+z1]*tablex7 ) * tabley4 +
				   ( data[x1+y5+z1]*tablex1 + data[x2+y5+z1]*tablex2 + data[x3+y5+z1]*tablex3 +
					 data[x4+y5+z1]*tablex4 + data[x5+y5+z1]*tablex5 + data[x6+y5+z1]*tablex6 +
					 data[x7+y5+z1]*tablex7 ) * tabley5 +
				   ( data[x1+y6+z1]*tablex1 + data[x2+y6+z1]*tablex2 + data[x3+y6+z1]*tablex3 +
					 data[x4+y6+z1]*tablex4 + data[x5+y6+z1]*tablex5 + data[x6+y6+z1]*tablex6 +
					 data[x7+y6+z1]*tablex7 ) * tabley6 +
				   ( data[x1+y7+z1]*tablex1 + data[x2+y7+z1]*tablex2 + data[x3+y7+z1]*tablex3 +
					 data[x4+y7+z1]*tablex4 + data[x5+y7+z1]*tablex5 + data[x6+y7+z1]*tablex6 +
					 data[x7+y7+z1]*tablex7 ) * tabley7 ) *tablez1 +
		         ( ( data[x1+y1+z2]*tablex1 + data[x2+y1+z2]*tablex2 + data[x3+y1+z2]*tablex3 +
		             data[x4+y1+z2]*tablex4 + data[x5+y1+z2]*tablex5 + data[x6+y1+z2]*tablex6 +
					 data[x7+y1+z2]*tablex7 ) * tabley1 +
				   ( data[x1+y2+z2]*tablex1 + data[x2+y2+z2]*tablex2 + data[x3+y2+z2]*tablex3 +
					 data[x4+y2+z2]*tablex4 + data[x5+y2+z2]*tablex5 + data[x6+y2+z2]*tablex6 +
					 data[x7+y2+z2]*tablex7 ) * tabley2 +
				   ( data[x1+y3+z2]*tablex1 + data[x2+y3+z2]*tablex2 + data[x3+y3+z2]*tablex3 +
					 data[x4+y3+z2]*tablex4 + data[x5+y3+z2]*tablex5 + data[x6+y3+z2]*tablex6 +
					 data[x7+y3+z2]*tablex7 ) * tabley3 +
				   ( data[x1+y4+z2]*tablex1 + data[x2+y4+z2]*tablex2 + data[x3+y4+z2]*tablex3 +
					 data[x4+y4+z2]*tablex4 + data[x5+y4+z2]*tablex5 + data[x6+y4+z2]*tablex6 +
					 data[x7+y4+z2]*tablex7 ) * tabley4 +
				   ( data[x1+y5+z2]*tablex1 + data[x2+y5+z2]*tablex2 + data[x3+y5+z2]*tablex3 +
					 data[x4+y5+z2]*tablex4 + data[x5+y5+z2]*tablex5 + data[x6+y5+z2]*tablex6 +
					 data[x7+y5+z2]*tablex7 ) * tabley5 +
				   ( data[x1+y6+z2]*tablex1 + data[x2+y6+z2]*tablex2 + data[x3+y6+z2]*tablex3 +
					 data[x4+y6+z2]*tablex4 + data[x5+y6+z2]*tablex5 + data[x6+y6+z2]*tablex6 +
					 data[x7+y6+z2]*tablex7 ) * tabley6 +
				   ( data[x1+y7+z2]*tablex1 + data[x2+y7+z2]*tablex2 + data[x3+y7+z2]*tablex3 +
					 data[x4+y7+z2]*tablex4 + data[x5+y7+z2]*tablex5 + data[x6+y7+z2]*tablex6 +
					 data[x7+y7+z2]*tablex7 ) * tabley7 ) *tablez2 +
		         ( ( data[x1+y1+z3]*tablex1 + data[x2+y1+z3]*tablex2 + data[x3+y1+z3]*tablex3 +
		             data[x4+y1+z3]*tablex4 + data[x5+y1+z3]*tablex5 + data[x6+y1+z3]*tablex6 +
					 data[x7+y1+z3]*tablex7 ) * tabley1 +
				   ( data[x1+y2+z3]*tablex1 + data[x2+y2+z3]*tablex2 + data[x3+y2+z3]*tablex3 +
					 data[x4+y2+z3]*tablex4 + data[x5+y2+z3]*tablex5 + data[x6+y2+z3]*tablex6 +
					 data[x7+y2+z3]*tablex7 ) * tabley2 +
				   ( data[x1+y3+z3]*tablex1 + data[x2+y3+z3]*tablex2 + data[x3+y3+z3]*tablex3 +
					 data[x4+y3+z3]*tablex4 + data[x5+y3+z3]*tablex5 + data[x6+y3+z3]*tablex6 +
					 data[x7+y3+z3]*tablex7 ) * tabley3 +
				   ( data[x1+y4+z3]*tablex1 + data[x2+y4+z3]*tablex2 + data[x3+y4+z3]*tablex3 +
					 data[x4+y4+z3]*tablex4 + data[x5+y4+z3]*tablex5 + data[x6+y4+z3]*tablex6 +
					 data[x7+y4+z3]*tablex7 ) * tabley4 +
				   ( data[x1+y5+z3]*tablex1 + data[x2+y5+z3]*tablex2 + data[x3+y5+z3]*tablex3 +
					 data[x4+y5+z3]*tablex4 + data[x5+y5+z3]*tablex5 + data[x6+y5+z3]*tablex6 +
					 data[x7+y5+z3]*tablex7 ) * tabley5 +
				   ( data[x1+y6+z3]*tablex1 + data[x2+y6+z3]*tablex2 + data[x3+y6+z3]*tablex3 +
					 data[x4+y6+z3]*tablex4 + data[x5+y6+z3]*tablex5 + data[x6+y6+z3]*tablex6 +
					 data[x7+y6+z3]*tablex7 ) * tabley6 +
				   ( data[x1+y7+z3]*tablex1 + data[x2+y7+z3]*tablex2 + data[x3+y7+z3]*tablex3 +
					 data[x4+y7+z3]*tablex4 + data[x5+y7+z3]*tablex5 + data[x6+y7+z3]*tablex6 +
					 data[x7+y7+z3]*tablex7 ) * tabley7 ) *tablez3 +
		         ( ( data[x1+y1+z4]*tablex1 + data[x2+y1+z4]*tablex2 + data[x3+y1+z4]*tablex3 +
		             data[x4+y1+z4]*tablex4 + data[x5+y1+z4]*tablex5 + data[x6+y1+z4]*tablex6 +
					 data[x7+y1+z4]*tablex7 ) * tabley1 +
				   ( data[x1+y2+z4]*tablex1 + data[x2+y2+z4]*tablex2 + data[x3+y2+z4]*tablex3 +
					 data[x4+y2+z4]*tablex4 + data[x5+y2+z4]*tablex5 + data[x6+y2+z4]*tablex6 +
					 data[x7+y2+z4]*tablex7 ) * tabley2 +
				   ( data[x1+y3+z4]*tablex1 + data[x2+y3+z4]*tablex2 + data[x3+y3+z4]*tablex3 +
					 data[x4+y3+z4]*tablex4 + data[x5+y3+z4]*tablex5 + data[x6+y3+z4]*tablex6 +
					 data[x7+y3+z4]*tablex7 ) * tabley3 +
				   ( data[x1+y4+z4]*tablex1 + data[x2+y4+z4]*tablex2 + data[x3+y4+z4]*tablex3 +
					 data[x4+y4+z4]*tablex4 + data[x5+y4+z4]*tablex5 + data[x6+y4+z4]*tablex6 +
					 data[x7+y4+z4]*tablex7 ) * tabley4 +
				   ( data[x1+y5+z4]*tablex1 + data[x2+y5+z4]*tablex2 + data[x3+y5+z4]*tablex3 +
					 data[x4+y5+z4]*tablex4 + data[x5+y5+z4]*tablex5 + data[x6+y5+z4]*tablex6 +
					 data[x7+y5+z4]*tablex7 ) * tabley5 +
				   ( data[x1+y6+z4]*tablex1 + data[x2+y6+z4]*tablex2 + data[x3+y6+z4]*tablex3 +
					 data[x4+y6+z4]*tablex4 + data[x5+y6+z4]*tablex5 + data[x6+y6+z4]*tablex6 +
					 data[x7+y6+z4]*tablex7 ) * tabley6 +
				   ( data[x1+y7+z4]*tablex1 + data[x2+y7+z4]*tablex2 + data[x3+y7+z4]*tablex3 +
					 data[x4+y7+z4]*tablex4 + data[x5+y7+z4]*tablex5 + data[x6+y7+z4]*tablex6 +
					 data[x7+y7+z4]*tablex7 ) * tabley7 ) *tablez4 +
		         ( ( data[x1+y1+z5]*tablex1 + data[x2+y1+z5]*tablex2 + data[x3+y1+z5]*tablex3 +
		             data[x4+y1+z5]*tablex4 + data[x5+y1+z5]*tablex5 + data[x6+y1+z5]*tablex6 +
					 data[x7+y1+z5]*tablex7 ) * tabley1 +
				   ( data[x1+y2+z5]*tablex1 + data[x2+y2+z5]*tablex2 + data[x3+y2+z5]*tablex3 +
					 data[x4+y2+z5]*tablex4 + data[x5+y2+z5]*tablex5 + data[x6+y2+z5]*tablex6 +
					 data[x7+y2+z5]*tablex7 ) * tabley2 +
				   ( data[x1+y3+z5]*tablex1 + data[x2+y3+z5]*tablex2 + data[x3+y3+z5]*tablex3 +
					 data[x4+y3+z5]*tablex4 + data[x5+y3+z5]*tablex5 + data[x6+y3+z5]*tablex6 +
					 data[x7+y3+z5]*tablex7 ) * tabley3 +
				   ( data[x1+y4+z5]*tablex1 + data[x2+y4+z5]*tablex2 + data[x3+y4+z5]*tablex3 +
					 data[x4+y4+z5]*tablex4 + data[x5+y4+z5]*tablex5 + data[x6+y4+z5]*tablex6 +
					 data[x7+y4+z5]*tablex7 ) * tabley4 +
				   ( data[x1+y5+z5]*tablex1 + data[x2+y5+z5]*tablex2 + data[x3+y5+z5]*tablex3 +
					 data[x4+y5+z5]*tablex4 + data[x5+y5+z5]*tablex5 + data[x6+y5+z5]*tablex6 +
					 data[x7+y5+z5]*tablex7 ) * tabley5 +
				   ( data[x1+y6+z5]*tablex1 + data[x2+y6+z5]*tablex2 + data[x3+y6+z5]*tablex3 +
					 data[x4+y6+z5]*tablex4 + data[x5+y6+z5]*tablex5 + data[x6+y6+z5]*tablex6 +
					 data[x7+y6+z5]*tablex7 ) * tabley6 +
				   ( data[x1+y7+z5]*tablex1 + data[x2+y7+z5]*tablex2 + data[x3+y7+z5]*tablex3 +
					 data[x4+y7+z5]*tablex4 + data[x5+y7+z5]*tablex5 + data[x6+y7+z5]*tablex6 +
					 data[x7+y7+z5]*tablex7 ) * tabley7 ) *tablez5 +
		         ( ( data[x1+y1+z6]*tablex1 + data[x2+y1+z6]*tablex2 + data[x3+y1+z6]*tablex3 +
		             data[x4+y1+z6]*tablex4 + data[x5+y1+z6]*tablex5 + data[x6+y1+z6]*tablex6 +
					 data[x7+y1+z6]*tablex7 ) * tabley1 +
				   ( data[x1+y2+z6]*tablex1 + data[x2+y2+z6]*tablex2 + data[x3+y2+z6]*tablex3 +
					 data[x4+y2+z6]*tablex4 + data[x5+y2+z6]*tablex5 + data[x6+y2+z6]*tablex6 +
					 data[x7+y2+z6]*tablex7 ) * tabley2 +
				   ( data[x1+y3+z6]*tablex1 + data[x2+y3+z6]*tablex2 + data[x3+y3+z6]*tablex3 +
					 data[x4+y3+z6]*tablex4 + data[x5+y3+z6]*tablex5 + data[x6+y3+z6]*tablex6 +
					 data[x7+y3+z6]*tablex7 ) * tabley3 +
				   ( data[x1+y4+z6]*tablex1 + data[x2+y4+z6]*tablex2 + data[x3+y4+z6]*tablex3 +
					 data[x4+y4+z6]*tablex4 + data[x5+y4+z6]*tablex5 + data[x6+y4+z6]*tablex6 +
					 data[x7+y4+z6]*tablex7 ) * tabley4 +
				   ( data[x1+y5+z6]*tablex1 + data[x2+y5+z6]*tablex2 + data[x3+y5+z6]*tablex3 +
					 data[x4+y5+z6]*tablex4 + data[x5+y5+z6]*tablex5 + data[x6+y5+z6]*tablex6 +
					 data[x7+y5+z6]*tablex7 ) * tabley5 +
				   ( data[x1+y6+z6]*tablex1 + data[x2+y6+z6]*tablex2 + data[x3+y6+z6]*tablex3 +
					 data[x4+y6+z6]*tablex4 + data[x5+y6+z6]*tablex5 + data[x6+y6+z6]*tablex6 +
					 data[x7+y6+z6]*tablex7 ) * tabley6 +
				   ( data[x1+y7+z6]*tablex1 + data[x2+y7+z6]*tablex2 + data[x3+y7+z6]*tablex3 +
					 data[x4+y7+z6]*tablex4 + data[x5+y7+z6]*tablex5 + data[x6+y7+z6]*tablex6 +
					 data[x7+y7+z6]*tablex7 ) * tabley7 ) *tablez6 +
		         ( ( data[x1+y1+z7]*tablex1 + data[x2+y1+z7]*tablex2 + data[x3+y1+z7]*tablex3 +
		             data[x4+y1+z7]*tablex4 + data[x5+y1+z7]*tablex5 + data[x6+y1+z7]*tablex6 +
					 data[x7+y1+z7]*tablex7 ) * tabley1 +
				   ( data[x1+y2+z7]*tablex1 + data[x2+y2+z7]*tablex2 + data[x3+y2+z7]*tablex3 +
					 data[x4+y2+z7]*tablex4 + data[x5+y2+z7]*tablex5 + data[x6+y2+z7]*tablex6 +
					 data[x7+y2+z7]*tablex7 ) * tabley2 +
				   ( data[x1+y3+z7]*tablex1 + data[x2+y3+z7]*tablex2 + data[x3+y3+z7]*tablex3 +
					 data[x4+y3+z7]*tablex4 + data[x5+y3+z7]*tablex5 + data[x6+y3+z7]*tablex6 +
					 data[x7+y3+z7]*tablex7 ) * tabley3 +
				   ( data[x1+y4+z7]*tablex1 + data[x2+y4+z7]*tablex2 + data[x3+y4+z7]*tablex3 +
					 data[x4+y4+z7]*tablex4 + data[x5+y4+z7]*tablex5 + data[x6+y4+z7]*tablex6 +
					 data[x7+y4+z7]*tablex7 ) * tabley4 +
				   ( data[x1+y5+z7]*tablex1 + data[x2+y5+z7]*tablex2 + data[x3+y5+z7]*tablex3 +
					 data[x4+y5+z7]*tablex4 + data[x5+y5+z7]*tablex5 + data[x6+y5+z7]*tablex6 +
					 data[x7+y5+z7]*tablex7 ) * tabley5 +
				   ( data[x1+y6+z7]*tablex1 + data[x2+y6+z7]*tablex2 + data[x3+y6+z7]*tablex3 +
					 data[x4+y6+z7]*tablex4 + data[x5+y6+z7]*tablex5 + data[x6+y6+z7]*tablex6 +
					 data[x7+y6+z7]*tablex7 ) * tabley6 +
				   ( data[x1+y7+z7]*tablex1 + data[x2+y7+z7]*tablex2 + data[x3+y7+z7]*tablex3 +
					 data[x4+y7+z7]*tablex4 + data[x5+y7+z7]*tablex5 + data[x6+y7+z7]*tablex6 +
					 data[x7+y7+z7]*tablex7 ) * tabley7 ) *tablez7;

		w = (tablex1+tablex2+tablex3+tablex4+tablex5+tablex6+tablex7) *
		    (tabley1+tabley2+tabley3+tabley4+tabley5+tabley6+tabley7) *
		    (tablez1+tablez2+tablez3+tablez4+tablez5+tablez6+tablez7);
	}
        return pixel/w;
}

float  Util::get_pixel_conv_new_background(int nx, int ny, int nz, float delx, float dely, float delz, float* data, Util::KaiserBessel& kb, int xnew, int ynew) {
	int K = kb.get_window_size();
	int kbmin = -K/2;
	int kbmax = -kbmin;
	int kbc = kbmax+1;

	float pixel =0.0f;
	float w=0.0f;

    float argdelx = delx; // adding this for 2D case where the wrap around is not done circulantly using restrict1.
	delx = restrict1(delx, nx);
	int inxold = int(round(delx));
	if ( ny < 2 ) {  //1D
		float tablex1 = kb.i0win_tab(delx-inxold+3);
		float tablex2 = kb.i0win_tab(delx-inxold+2);
		float tablex3 = kb.i0win_tab(delx-inxold+1);
		float tablex4 = kb.i0win_tab(delx-inxold);
		float tablex5 = kb.i0win_tab(delx-inxold-1);
		float tablex6 = kb.i0win_tab(delx-inxold-2);
		float tablex7 = kb.i0win_tab(delx-inxold-3);

		int x1, x2, x3, x4, x5, x6, x7;

	 	if ( inxold <= kbc || inxold >=nx-kbc-2 )  {
			x1 = (inxold-3+nx)%nx;
			x2 = (inxold-2+nx)%nx;
			x3 = (inxold-1+nx)%nx;
			x4 = (inxold  +nx)%nx;
			x5 = (inxold+1+nx)%nx;
			x6 = (inxold+2+nx)%nx;
			x7 = (inxold+3+nx)%nx;
	 	} else {
			x1 = inxold-3;
			x2 = inxold-2;
			x3 = inxold-1;
			x4 = inxold;
			x5 = inxold+1;
			x6 = inxold+2;
			x7 = inxold+3;
	 	}

		pixel = data[x1]*tablex1 + data[x2]*tablex2 + data[x3]*tablex3 +
				data[x4]*tablex4 + data[x5]*tablex5 + data[x6]*tablex6 +
				data[x7]*tablex7 ;

		w = tablex1+tablex2+tablex3+tablex4+tablex5+tablex6+tablex7;
	} else if ( nz < 2 ) {  // 2D

		delx = argdelx;
		// the wrap around is not done circulantly for 2D case; if (argdelx, argdely) is not in the image, then make them (xnew, ynew) which is definitely in the image
		if ((delx < 0.0f) || (delx >= (float) (nx)) || (dely < 0.0f) || (dely >= (float) (ny)) ){
	        delx = (float)xnew*2.0f;
	        dely = (float)ynew*2.0f;
		}

		int inxold = int(round(delx));
		int inyold = int(round(dely));

		float tablex1 = kb.i0win_tab(delx-inxold+3);
		float tablex2 = kb.i0win_tab(delx-inxold+2);
		float tablex3 = kb.i0win_tab(delx-inxold+1);
		float tablex4 = kb.i0win_tab(delx-inxold);
		float tablex5 = kb.i0win_tab(delx-inxold-1);
		float tablex6 = kb.i0win_tab(delx-inxold-2);
		float tablex7 = kb.i0win_tab(delx-inxold-3);

		float tabley1 = kb.i0win_tab(dely-inyold+3);
		float tabley2 = kb.i0win_tab(dely-inyold+2);
		float tabley3 = kb.i0win_tab(dely-inyold+1);
		float tabley4 = kb.i0win_tab(dely-inyold);
		float tabley5 = kb.i0win_tab(dely-inyold-1);
		float tabley6 = kb.i0win_tab(dely-inyold-2);
		float tabley7 = kb.i0win_tab(dely-inyold-3);

		int x1, x2, x3, x4, x5, x6, x7, y1, y2, y3, y4, y5, y6, y7;

	 	if ( inxold <= kbc || inxold >=nx-kbc-2 || inyold <= kbc || inyold >=ny-kbc-2 )  {
			x1 = (inxold-3+nx)%nx;
			x2 = (inxold-2+nx)%nx;
			x3 = (inxold-1+nx)%nx;
			x4 = (inxold  +nx)%nx;
			x5 = (inxold+1+nx)%nx;
			x6 = (inxold+2+nx)%nx;
			x7 = (inxold+3+nx)%nx;

			y1 = ((inyold-3+ny)%ny)*nx;
			y2 = ((inyold-2+ny)%ny)*nx;
			y3 = ((inyold-1+ny)%ny)*nx;
			y4 = ((inyold  +ny)%ny)*nx;
			y5 = ((inyold+1+ny)%ny)*nx;
			y6 = ((inyold+2+ny)%ny)*nx;
			y7 = ((inyold+3+ny)%ny)*nx;
		} else {
			x1 = inxold-3;
			x2 = inxold-2;
			x3 = inxold-1;
			x4 = inxold;
			x5 = inxold+1;
			x6 = inxold+2;
			x7 = inxold+3;

			y1 = (inyold-3)*nx;
			y2 = (inyold-2)*nx;
			y3 = (inyold-1)*nx;
			y4 = inyold*nx;
			y5 = (inyold+1)*nx;
			y6 = (inyold+2)*nx;
			y7 = (inyold+3)*nx;
		}

		pixel    = ( data[x1+y1]*tablex1 + data[x2+y1]*tablex2 + data[x3+y1]*tablex3 +
		             data[x4+y1]*tablex4 + data[x5+y1]*tablex5 + data[x6+y1]*tablex6 +
					 data[x7+y1]*tablex7 ) * tabley1 +
				   ( data[x1+y2]*tablex1 + data[x2+y2]*tablex2 + data[x3+y2]*tablex3 +
					 data[x4+y2]*tablex4 + data[x5+y2]*tablex5 + data[x6+y2]*tablex6 +
					 data[x7+y2]*tablex7 ) * tabley2 +
				   ( data[x1+y3]*tablex1 + data[x2+y3]*tablex2 + data[x3+y3]*tablex3 +
					 data[x4+y3]*tablex4 + data[x5+y3]*tablex5 + data[x6+y3]*tablex6 +
					 data[x7+y3]*tablex7 ) * tabley3 +
				   ( data[x1+y4]*tablex1 + data[x2+y4]*tablex2 + data[x3+y4]*tablex3 +
					 data[x4+y4]*tablex4 + data[x5+y4]*tablex5 + data[x6+y4]*tablex6 +
					 data[x7+y4]*tablex7 ) * tabley4 +
				   ( data[x1+y5]*tablex1 + data[x2+y5]*tablex2 + data[x3+y5]*tablex3 +
					 data[x4+y5]*tablex4 + data[x5+y5]*tablex5 + data[x6+y5]*tablex6 +
					 data[x7+y5]*tablex7 ) * tabley5 +
				   ( data[x1+y6]*tablex1 + data[x2+y6]*tablex2 + data[x3+y6]*tablex3 +
					 data[x4+y6]*tablex4 + data[x5+y6]*tablex5 + data[x6+y6]*tablex6 +
					 data[x7+y6]*tablex7 ) * tabley6 +
				   ( data[x1+y7]*tablex1 + data[x2+y7]*tablex2 + data[x3+y7]*tablex3 +
					 data[x4+y7]*tablex4 + data[x5+y7]*tablex5 + data[x6+y7]*tablex6 +
					 data[x7+y7]*tablex7 ) * tabley7;

		w = (tablex1+tablex2+tablex3+tablex4+tablex5+tablex6+tablex7) *
		    (tabley1+tabley2+tabley3+tabley4+tabley5+tabley6+tabley7);
	} else {  //  3D
		dely = restrict1(dely, ny);
		int inyold = int(Util::round(dely));
		delz = restrict1(delz, nz);
		int inzold = int(Util::round(delz));

		float tablex1 = kb.i0win_tab(delx-inxold+3);
		float tablex2 = kb.i0win_tab(delx-inxold+2);
		float tablex3 = kb.i0win_tab(delx-inxold+1);
		float tablex4 = kb.i0win_tab(delx-inxold);
		float tablex5 = kb.i0win_tab(delx-inxold-1);
		float tablex6 = kb.i0win_tab(delx-inxold-2);
		float tablex7 = kb.i0win_tab(delx-inxold-3);

		float tabley1 = kb.i0win_tab(dely-inyold+3);
		float tabley2 = kb.i0win_tab(dely-inyold+2);
		float tabley3 = kb.i0win_tab(dely-inyold+1);
		float tabley4 = kb.i0win_tab(dely-inyold);
		float tabley5 = kb.i0win_tab(dely-inyold-1);
		float tabley6 = kb.i0win_tab(dely-inyold-2);
		float tabley7 = kb.i0win_tab(dely-inyold-3);

		float tablez1 = kb.i0win_tab(delz-inzold+3);
		float tablez2 = kb.i0win_tab(delz-inzold+2);
		float tablez3 = kb.i0win_tab(delz-inzold+1);
		float tablez4 = kb.i0win_tab(delz-inzold);
		float tablez5 = kb.i0win_tab(delz-inzold-1);
		float tablez6 = kb.i0win_tab(delz-inzold-2);
		float tablez7 = kb.i0win_tab(delz-inzold-3);

		int x1, x2, x3, x4, x5, x6, x7, y1, y2, y3, y4, y5, y6, y7, z1, z2, z3, z4, z5, z6, z7;

		if ( inxold <= kbc || inxold >=nx-kbc-2 || inyold <= kbc || inyold >=ny-kbc-2 || inzold <= kbc || inzold >= nz-kbc-2 )  {
			x1 = (inxold-3+nx)%nx;
			x2 = (inxold-2+nx)%nx;
			x3 = (inxold-1+nx)%nx;
			x4 = (inxold  +nx)%nx;
			x5 = (inxold+1+nx)%nx;
			x6 = (inxold+2+nx)%nx;
			x7 = (inxold+3+nx)%nx;

			y1 = ((inyold-3+ny)%ny)*nx;
			y2 = ((inyold-2+ny)%ny)*nx;
			y3 = ((inyold-1+ny)%ny)*nx;
			y4 = ((inyold  +ny)%ny)*nx;
			y5 = ((inyold+1+ny)%ny)*nx;
			y6 = ((inyold+2+ny)%ny)*nx;
			y7 = ((inyold+3+ny)%ny)*nx;

			z1 = ((inzold-3+nz)%nz)*nx*ny;
			z2 = ((inzold-2+nz)%nz)*nx*ny;
			z3 = ((inzold-1+nz)%nz)*nx*ny;
			z4 = ((inzold  +nz)%nz)*nx*ny;
			z5 = ((inzold+1+nz)%nz)*nx*ny;
			z6 = ((inzold+2+nz)%nz)*nx*ny;
			z7 = ((inzold+3+nz)%nz)*nx*ny;
		} else {
			x1 = inxold-3;
			x2 = inxold-2;
			x3 = inxold-1;
			x4 = inxold;
			x5 = inxold+1;
			x6 = inxold+2;
			x7 = inxold+3;

			y1 = (inyold-3)*nx;
			y2 = (inyold-2)*nx;
			y3 = (inyold-1)*nx;
			y4 = inyold*nx;
			y5 = (inyold+1)*nx;
			y6 = (inyold+2)*nx;
			y7 = (inyold+3)*nx;

			z1 = (inzold-3)*nx*ny;
			z2 = (inzold-2)*nx*ny;
			z3 = (inzold-1)*nx*ny;
			z4 = inzold*nx*ny;
			z5 = (inzold+1)*nx*ny;
			z6 = (inzold+2)*nx*ny;
			z7 = (inzold+3)*nx*ny;
	 	}

		pixel  = ( ( data[x1+y1+z1]*tablex1 + data[x2+y1+z1]*tablex2 + data[x3+y1+z1]*tablex3 +
		             data[x4+y1+z1]*tablex4 + data[x5+y1+z1]*tablex5 + data[x6+y1+z1]*tablex6 +
					 data[x7+y1+z1]*tablex7 ) * tabley1 +
				   ( data[x1+y2+z1]*tablex1 + data[x2+y2+z1]*tablex2 + data[x3+y2+z1]*tablex3 +
					 data[x4+y2+z1]*tablex4 + data[x5+y2+z1]*tablex5 + data[x6+y2+z1]*tablex6 +
					 data[x7+y2+z1]*tablex7 ) * tabley2 +
				   ( data[x1+y3+z1]*tablex1 + data[x2+y3+z1]*tablex2 + data[x3+y3+z1]*tablex3 +
					 data[x4+y3+z1]*tablex4 + data[x5+y3+z1]*tablex5 + data[x6+y3+z1]*tablex6 +
					 data[x7+y3+z1]*tablex7 ) * tabley3 +
				   ( data[x1+y4+z1]*tablex1 + data[x2+y4+z1]*tablex2 + data[x3+y4+z1]*tablex3 +
					 data[x4+y4+z1]*tablex4 + data[x5+y4+z1]*tablex5 + data[x6+y4+z1]*tablex6 +
					 data[x7+y4+z1]*tablex7 ) * tabley4 +
				   ( data[x1+y5+z1]*tablex1 + data[x2+y5+z1]*tablex2 + data[x3+y5+z1]*tablex3 +
					 data[x4+y5+z1]*tablex4 + data[x5+y5+z1]*tablex5 + data[x6+y5+z1]*tablex6 +
					 data[x7+y5+z1]*tablex7 ) * tabley5 +
				   ( data[x1+y6+z1]*tablex1 + data[x2+y6+z1]*tablex2 + data[x3+y6+z1]*tablex3 +
					 data[x4+y6+z1]*tablex4 + data[x5+y6+z1]*tablex5 + data[x6+y6+z1]*tablex6 +
					 data[x7+y6+z1]*tablex7 ) * tabley6 +
				   ( data[x1+y7+z1]*tablex1 + data[x2+y7+z1]*tablex2 + data[x3+y7+z1]*tablex3 +
					 data[x4+y7+z1]*tablex4 + data[x5+y7+z1]*tablex5 + data[x6+y7+z1]*tablex6 +
					 data[x7+y7+z1]*tablex7 ) * tabley7 ) *tablez1 +
		         ( ( data[x1+y1+z2]*tablex1 + data[x2+y1+z2]*tablex2 + data[x3+y1+z2]*tablex3 +
		             data[x4+y1+z2]*tablex4 + data[x5+y1+z2]*tablex5 + data[x6+y1+z2]*tablex6 +
					 data[x7+y1+z2]*tablex7 ) * tabley1 +
				   ( data[x1+y2+z2]*tablex1 + data[x2+y2+z2]*tablex2 + data[x3+y2+z2]*tablex3 +
					 data[x4+y2+z2]*tablex4 + data[x5+y2+z2]*tablex5 + data[x6+y2+z2]*tablex6 +
					 data[x7+y2+z2]*tablex7 ) * tabley2 +
				   ( data[x1+y3+z2]*tablex1 + data[x2+y3+z2]*tablex2 + data[x3+y3+z2]*tablex3 +
					 data[x4+y3+z2]*tablex4 + data[x5+y3+z2]*tablex5 + data[x6+y3+z2]*tablex6 +
					 data[x7+y3+z2]*tablex7 ) * tabley3 +
				   ( data[x1+y4+z2]*tablex1 + data[x2+y4+z2]*tablex2 + data[x3+y4+z2]*tablex3 +
					 data[x4+y4+z2]*tablex4 + data[x5+y4+z2]*tablex5 + data[x6+y4+z2]*tablex6 +
					 data[x7+y4+z2]*tablex7 ) * tabley4 +
				   ( data[x1+y5+z2]*tablex1 + data[x2+y5+z2]*tablex2 + data[x3+y5+z2]*tablex3 +
					 data[x4+y5+z2]*tablex4 + data[x5+y5+z2]*tablex5 + data[x6+y5+z2]*tablex6 +
					 data[x7+y5+z2]*tablex7 ) * tabley5 +
				   ( data[x1+y6+z2]*tablex1 + data[x2+y6+z2]*tablex2 + data[x3+y6+z2]*tablex3 +
					 data[x4+y6+z2]*tablex4 + data[x5+y6+z2]*tablex5 + data[x6+y6+z2]*tablex6 +
					 data[x7+y6+z2]*tablex7 ) * tabley6 +
				   ( data[x1+y7+z2]*tablex1 + data[x2+y7+z2]*tablex2 + data[x3+y7+z2]*tablex3 +
					 data[x4+y7+z2]*tablex4 + data[x5+y7+z2]*tablex5 + data[x6+y7+z2]*tablex6 +
					 data[x7+y7+z2]*tablex7 ) * tabley7 ) *tablez2 +
		         ( ( data[x1+y1+z3]*tablex1 + data[x2+y1+z3]*tablex2 + data[x3+y1+z3]*tablex3 +
		             data[x4+y1+z3]*tablex4 + data[x5+y1+z3]*tablex5 + data[x6+y1+z3]*tablex6 +
					 data[x7+y1+z3]*tablex7 ) * tabley1 +
				   ( data[x1+y2+z3]*tablex1 + data[x2+y2+z3]*tablex2 + data[x3+y2+z3]*tablex3 +
					 data[x4+y2+z3]*tablex4 + data[x5+y2+z3]*tablex5 + data[x6+y2+z3]*tablex6 +
					 data[x7+y2+z3]*tablex7 ) * tabley2 +
				   ( data[x1+y3+z3]*tablex1 + data[x2+y3+z3]*tablex2 + data[x3+y3+z3]*tablex3 +
					 data[x4+y3+z3]*tablex4 + data[x5+y3+z3]*tablex5 + data[x6+y3+z3]*tablex6 +
					 data[x7+y3+z3]*tablex7 ) * tabley3 +
				   ( data[x1+y4+z3]*tablex1 + data[x2+y4+z3]*tablex2 + data[x3+y4+z3]*tablex3 +
					 data[x4+y4+z3]*tablex4 + data[x5+y4+z3]*tablex5 + data[x6+y4+z3]*tablex6 +
					 data[x7+y4+z3]*tablex7 ) * tabley4 +
				   ( data[x1+y5+z3]*tablex1 + data[x2+y5+z3]*tablex2 + data[x3+y5+z3]*tablex3 +
					 data[x4+y5+z3]*tablex4 + data[x5+y5+z3]*tablex5 + data[x6+y5+z3]*tablex6 +
					 data[x7+y5+z3]*tablex7 ) * tabley5 +
				   ( data[x1+y6+z3]*tablex1 + data[x2+y6+z3]*tablex2 + data[x3+y6+z3]*tablex3 +
					 data[x4+y6+z3]*tablex4 + data[x5+y6+z3]*tablex5 + data[x6+y6+z3]*tablex6 +
					 data[x7+y6+z3]*tablex7 ) * tabley6 +
				   ( data[x1+y7+z3]*tablex1 + data[x2+y7+z3]*tablex2 + data[x3+y7+z3]*tablex3 +
					 data[x4+y7+z3]*tablex4 + data[x5+y7+z3]*tablex5 + data[x6+y7+z3]*tablex6 +
					 data[x7+y7+z3]*tablex7 ) * tabley7 ) *tablez3 +
		         ( ( data[x1+y1+z4]*tablex1 + data[x2+y1+z4]*tablex2 + data[x3+y1+z4]*tablex3 +
		             data[x4+y1+z4]*tablex4 + data[x5+y1+z4]*tablex5 + data[x6+y1+z4]*tablex6 +
					 data[x7+y1+z4]*tablex7 ) * tabley1 +
				   ( data[x1+y2+z4]*tablex1 + data[x2+y2+z4]*tablex2 + data[x3+y2+z4]*tablex3 +
					 data[x4+y2+z4]*tablex4 + data[x5+y2+z4]*tablex5 + data[x6+y2+z4]*tablex6 +
					 data[x7+y2+z4]*tablex7 ) * tabley2 +
				   ( data[x1+y3+z4]*tablex1 + data[x2+y3+z4]*tablex2 + data[x3+y3+z4]*tablex3 +
					 data[x4+y3+z4]*tablex4 + data[x5+y3+z4]*tablex5 + data[x6+y3+z4]*tablex6 +
					 data[x7+y3+z4]*tablex7 ) * tabley3 +
				   ( data[x1+y4+z4]*tablex1 + data[x2+y4+z4]*tablex2 + data[x3+y4+z4]*tablex3 +
					 data[x4+y4+z4]*tablex4 + data[x5+y4+z4]*tablex5 + data[x6+y4+z4]*tablex6 +
					 data[x7+y4+z4]*tablex7 ) * tabley4 +
				   ( data[x1+y5+z4]*tablex1 + data[x2+y5+z4]*tablex2 + data[x3+y5+z4]*tablex3 +
					 data[x4+y5+z4]*tablex4 + data[x5+y5+z4]*tablex5 + data[x6+y5+z4]*tablex6 +
					 data[x7+y5+z4]*tablex7 ) * tabley5 +
				   ( data[x1+y6+z4]*tablex1 + data[x2+y6+z4]*tablex2 + data[x3+y6+z4]*tablex3 +
					 data[x4+y6+z4]*tablex4 + data[x5+y6+z4]*tablex5 + data[x6+y6+z4]*tablex6 +
					 data[x7+y6+z4]*tablex7 ) * tabley6 +
				   ( data[x1+y7+z4]*tablex1 + data[x2+y7+z4]*tablex2 + data[x3+y7+z4]*tablex3 +
					 data[x4+y7+z4]*tablex4 + data[x5+y7+z4]*tablex5 + data[x6+y7+z4]*tablex6 +
					 data[x7+y7+z4]*tablex7 ) * tabley7 ) *tablez4 +
					 ( ( data[x1+y1+z5]*tablex1 + data[x2+y1+z5]*tablex2 + data[x3+y1+z5]*tablex3 +
		             data[x4+y1+z5]*tablex4 + data[x5+y1+z5]*tablex5 + data[x6+y1+z5]*tablex6 +
					 data[x7+y1+z5]*tablex7 ) * tabley1 +
				   ( data[x1+y2+z5]*tablex1 + data[x2+y2+z5]*tablex2 + data[x3+y2+z5]*tablex3 +
					 data[x4+y2+z5]*tablex4 + data[x5+y2+z5]*tablex5 + data[x6+y2+z5]*tablex6 +
					 data[x7+y2+z5]*tablex7 ) * tabley2 +
				   ( data[x1+y3+z5]*tablex1 + data[x2+y3+z5]*tablex2 + data[x3+y3+z5]*tablex3 +
					 data[x4+y3+z5]*tablex4 + data[x5+y3+z5]*tablex5 + data[x6+y3+z5]*tablex6 +
					 data[x7+y3+z5]*tablex7 ) * tabley3 +
				   ( data[x1+y4+z5]*tablex1 + data[x2+y4+z5]*tablex2 + data[x3+y4+z5]*tablex3 +
					 data[x4+y4+z5]*tablex4 + data[x5+y4+z5]*tablex5 + data[x6+y4+z5]*tablex6 +
					 data[x7+y4+z5]*tablex7 ) * tabley4 +
				   ( data[x1+y5+z5]*tablex1 + data[x2+y5+z5]*tablex2 + data[x3+y5+z5]*tablex3 +
					 data[x4+y5+z5]*tablex4 + data[x5+y5+z5]*tablex5 + data[x6+y5+z5]*tablex6 +
					 data[x7+y5+z5]*tablex7 ) * tabley5 +
				   ( data[x1+y6+z5]*tablex1 + data[x2+y6+z5]*tablex2 + data[x3+y6+z5]*tablex3 +
					 data[x4+y6+z5]*tablex4 + data[x5+y6+z5]*tablex5 + data[x6+y6+z5]*tablex6 +
					 data[x7+y6+z5]*tablex7 ) * tabley6 +
				   ( data[x1+y7+z5]*tablex1 + data[x2+y7+z5]*tablex2 + data[x3+y7+z5]*tablex3 +
					 data[x4+y7+z5]*tablex4 + data[x5+y7+z5]*tablex5 + data[x6+y7+z5]*tablex6 +
					 data[x7+y7+z5]*tablex7 ) * tabley7 ) *tablez5 +
		         ( ( data[x1+y1+z6]*tablex1 + data[x2+y1+z6]*tablex2 + data[x3+y1+z6]*tablex3 +
		             data[x4+y1+z6]*tablex4 + data[x5+y1+z6]*tablex5 + data[x6+y1+z6]*tablex6 +
					 data[x7+y1+z6]*tablex7 ) * tabley1 +
				   ( data[x1+y2+z6]*tablex1 + data[x2+y2+z6]*tablex2 + data[x3+y2+z6]*tablex3 +
					 data[x4+y2+z6]*tablex4 + data[x5+y2+z6]*tablex5 + data[x6+y2+z6]*tablex6 +
					 data[x7+y2+z6]*tablex7 ) * tabley2 +
				   ( data[x1+y3+z6]*tablex1 + data[x2+y3+z6]*tablex2 + data[x3+y3+z6]*tablex3 +
					 data[x4+y3+z6]*tablex4 + data[x5+y3+z6]*tablex5 + data[x6+y3+z6]*tablex6 +
					 data[x7+y3+z6]*tablex7 ) * tabley3 +
				   ( data[x1+y4+z6]*tablex1 + data[x2+y4+z6]*tablex2 + data[x3+y4+z6]*tablex3 +
					 data[x4+y4+z6]*tablex4 + data[x5+y4+z6]*tablex5 + data[x6+y4+z6]*tablex6 +
					 data[x7+y4+z6]*tablex7 ) * tabley4 +
				   ( data[x1+y5+z6]*tablex1 + data[x2+y5+z6]*tablex2 + data[x3+y5+z6]*tablex3 +
					 data[x4+y5+z6]*tablex4 + data[x5+y5+z6]*tablex5 + data[x6+y5+z6]*tablex6 +
					 data[x7+y5+z6]*tablex7 ) * tabley5 +
				   ( data[x1+y6+z6]*tablex1 + data[x2+y6+z6]*tablex2 + data[x3+y6+z6]*tablex3 +
					 data[x4+y6+z6]*tablex4 + data[x5+y6+z6]*tablex5 + data[x6+y6+z6]*tablex6 +
					 data[x7+y6+z6]*tablex7 ) * tabley6 +
				   ( data[x1+y7+z6]*tablex1 + data[x2+y7+z6]*tablex2 + data[x3+y7+z6]*tablex3 +
					 data[x4+y7+z6]*tablex4 + data[x5+y7+z6]*tablex5 + data[x6+y7+z6]*tablex6 +
					 data[x7+y7+z6]*tablex7 ) * tabley7 ) *tablez6 +
		         ( ( data[x1+y1+z7]*tablex1 + data[x2+y1+z7]*tablex2 + data[x3+y1+z7]*tablex3 +
		             data[x4+y1+z7]*tablex4 + data[x5+y1+z7]*tablex5 + data[x6+y1+z7]*tablex6 +
					 data[x7+y1+z7]*tablex7 ) * tabley1 +
				   ( data[x1+y2+z7]*tablex1 + data[x2+y2+z7]*tablex2 + data[x3+y2+z7]*tablex3 +
					 data[x4+y2+z7]*tablex4 + data[x5+y2+z7]*tablex5 + data[x6+y2+z7]*tablex6 +
					 data[x7+y2+z7]*tablex7 ) * tabley2 +
				   ( data[x1+y3+z7]*tablex1 + data[x2+y3+z7]*tablex2 + data[x3+y3+z7]*tablex3 +
					 data[x4+y3+z7]*tablex4 + data[x5+y3+z7]*tablex5 + data[x6+y3+z7]*tablex6 +
					 data[x7+y3+z7]*tablex7 ) * tabley3 +
				   ( data[x1+y4+z7]*tablex1 + data[x2+y4+z7]*tablex2 + data[x3+y4+z7]*tablex3 +
					 data[x4+y4+z7]*tablex4 + data[x5+y4+z7]*tablex5 + data[x6+y4+z7]*tablex6 +
					 data[x7+y4+z7]*tablex7 ) * tabley4 +
				   ( data[x1+y5+z7]*tablex1 + data[x2+y5+z7]*tablex2 + data[x3+y5+z7]*tablex3 +
					 data[x4+y5+z7]*tablex4 + data[x5+y5+z7]*tablex5 + data[x6+y5+z7]*tablex6 +
					 data[x7+y5+z7]*tablex7 ) * tabley5 +
				   ( data[x1+y6+z7]*tablex1 + data[x2+y6+z7]*tablex2 + data[x3+y6+z7]*tablex3 +
					 data[x4+y6+z7]*tablex4 + data[x5+y6+z7]*tablex5 + data[x6+y6+z7]*tablex6 +
					 data[x7+y6+z7]*tablex7 ) * tabley6 +
				   ( data[x1+y7+z7]*tablex1 + data[x2+y7+z7]*tablex2 + data[x3+y7+z7]*tablex3 +
					 data[x4+y7+z7]*tablex4 + data[x5+y7+z7]*tablex5 + data[x6+y7+z7]*tablex6 +
					 data[x7+y7+z7]*tablex7 ) * tabley7 ) *tablez7;

		w = (tablex1+tablex2+tablex3+tablex4+tablex5+tablex6+tablex7) *
		    (tabley1+tabley2+tabley3+tabley4+tabley5+tabley6+tabley7) *
		    (tablez1+tablez2+tablez3+tablez4+tablez5+tablez6+tablez7);
	}
        return pixel/w;
}

/*
complex<float> Util::extractpoint2(int nx, int ny, float nuxnew, float nuynew, EMData *fimage, Util::KaiserBessel& kb) {

	int nxreal = nx - 2;
	if (nxreal != ny)
		throw ImageDimensionException("extractpoint requires ny == nx");
	int nhalf = nxreal/2;
	int kbsize = kb.get_window_size();
	int kbmin = -kbsize/2;
	int kbmax = -kbmin;
	bool flip = (nuxnew < 0.f);
	if (flip) {
		nuxnew *= -1;
		nuynew *= -1;
	}
	// put (xnew,ynew) on a grid.  The indices will be wrong for
	// the Fourier elements in the image, but the grid sizing will
	// be correct.
	int ixn = int(Util::round(nuxnew));
	int iyn = int(Util::round(nuynew));
	// set up some temporary weighting arrays
	float* wy0 = new float[kbmax - kbmin + 1];
	float* wy = wy0 - kbmin; // wy[kbmin:kbmax]
	float* wx0 = new float[kbmax - kbmin + 1];
	float* wx = wx0 - kbmin;
	for (int i = kbmin; i <= kbmax; i++) {
			int iyp = iyn + i;
			wy[i] = kb.i0win_tab(nuynew - iyp);
			int ixp = ixn + i;
			wx[i] = kb.i0win_tab(nuxnew - ixp);
	}
	// restrict loops to non-zero elements
	int iymin = 0;
	for (int iy = kbmin; iy <= -1; iy++) {
		if (wy[iy] != 0.f) {
			iymin = iy;
			break;
		}
	}
	int iymax = 0;
	for (int iy = kbmax; iy >= 1; iy--) {
		if (wy[iy] != 0.f) {
			iymax = iy;
			break;
		}
	}
	int ixmin = 0;
	for (int ix = kbmin; ix <= -1; ix++) {
		if (wx[ix] != 0.f) {
			ixmin = ix;
			break;
		}
	}
	int ixmax = 0;
	for (int ix = kbmax; ix >= 1; ix--) {
		if (wx[ix] != 0.f) {
			ixmax = ix;
			break;
		}
	}
	float wsum = 0.0f;
	for (int iy = iymin; iy <= iymax; iy++)
		for (int ix = ixmin; ix <= ixmax; ix++)
			wsum += wx[ix]*wy[iy];

	complex<float> result(0.f,0.f);
	if ((ixn >= -kbmin) && (ixn <= nhalf-1-kbmax) && (iyn >= -nhalf-kbmin) && (iyn <= nhalf-1-kbmax)) {
		// (xin,yin) not within window border from the edge
		for (int iy = iymin; iy <= iymax; iy++) {
			int iyp = iyn + iy;
			for (int ix = ixmin; ix <= ixmax; ix++) {
				int ixp = ixn + ix;
				float w = wx[ix]*wy[iy];
				complex<float> val = fimage->cmplx(ixp,iyp);
				result += val*w;
			}
		}
	} else {
		// points that "stick out"
		for (int iy = iymin; iy <= iymax; iy++) {
			int iyp = iyn + iy;
			for (int ix = ixmin; ix <= ixmax; ix++) {
				int ixp = ixn + ix;
				bool mirror = false;
				int ixt= ixp, iyt= iyp;
				if (ixt < 0) {
					ixt = -ixt;
					iyt = -iyt;
					mirror = !mirror;
				}
				if (ixt > nhalf) {
					ixt = nxreal - ixt;
					iyt = -iyt;
					mirror = !mirror;
				}
				if (iyt > nhalf-1)  iyt -= nxreal;
				if (iyt < -nhalf)   iyt += nxreal;
				float w = wx[ix]*wy[iy];
				complex<float> val = fimage->cmplx(ixt,iyt);
				if (mirror)  result += conj(val)*w;
				else         result += val*w;
			}
		}
	}
	if (flip)  result = conj(result)/wsum;
	else result /= wsum;
	delete [] wx0;
	delete [] wy0;
	return result;
}*/

/*
complex<float> Util::extractpoint2(int nx, int ny, float nuxnew, float nuynew, EMData *fimage, Util::KaiserBessel& kb) {

	int nxreal = nx - 2;
	if (nxreal != ny)
		throw ImageDimensionException("extractpoint requires ny == nx");
	int nhalf = nxreal/2;
	bool flip = false;
	if (nuxnew < 0.f) {
		nuxnew *= -1;
		nuynew *= -1;
		flip = true;
	}
	if (nuynew >= nhalf-0.5)  {
		nuynew -= nxreal;
	} else if (nuynew < -nhalf-0.5) {
		nuynew += nxreal;
	}

	// put (xnew,ynew) on a grid.  The indices will be wrong for
	// the Fourier elements in the image, but the grid sizing will
	// be correct.
	int ixn = int(Util::round(nuxnew));
	int iyn = int(Util::round(nuynew));

	// set up some temporary weighting arrays
	static float wy[7];
	static float wx[7];

	float iynn = nuynew - iyn;
	wy[0] = kb.i0win_tab(iynn+3);
	wy[1] = kb.i0win_tab(iynn+2);
	wy[2] = kb.i0win_tab(iynn+1);
	wy[3] = kb.i0win_tab(iynn);
	wy[4] = kb.i0win_tab(iynn-1);
	wy[5] = kb.i0win_tab(iynn-2);
	wy[6] = kb.i0win_tab(iynn-3);

	float ixnn = nuxnew - ixn;
	wx[0] = kb.i0win_tab(ixnn+3);
	wx[1] = kb.i0win_tab(ixnn+2);
	wx[2] = kb.i0win_tab(ixnn+1);
	wx[3] = kb.i0win_tab(ixnn);
	wx[4] = kb.i0win_tab(ixnn-1);
	wx[5] = kb.i0win_tab(ixnn-2);
	wx[6] = kb.i0win_tab(ixnn-3);

	float wsum = (wx[0]+wx[1]+wx[2]+wx[3]+wx[4]+wx[5]+wx[6])*(wy[0]+wy[1]+wy[2]+wy[3]+wy[4]+wy[5]+wy[6]);

	complex<float> result(0.f,0.f);
	for (int iy = 0; iy < 7; iy++) {
		int iyp = iyn + iy - 3 ;
		for (int ix = 0; ix < 7; ix++) {
			int ixp = ixn + ix - 3;
			float w = wx[ix]*wy[iy];
			complex<float> val = fimage->cmplx(ixp,iyp);
			result += val*w;
		}
	}

	if (flip)  result = conj(result)/wsum;
	else result /= wsum;

	return result;
}*/


complex<float> Util::extractpoint2(int nx, int ny, float nuxnew, float nuynew, EMData *fimage, Util::KaiserBessel& kb) {

	int nxreal = nx - 2;
	if (nxreal != ny)
		throw ImageDimensionException("extractpoint requires ny == nx");
	int nhalf = nxreal/2;
	bool flip = (nuxnew < 0.f);
	if (flip) {
		nuxnew *= -1;
		nuynew *= -1;
	}
	if (nuynew >= nhalf-0.5)  {
		nuynew -= nxreal;
	} else if (nuynew < -nhalf-0.5) {
		nuynew += nxreal;
	}

	// put (xnew,ynew) on a grid.  The indices will be wrong for
	// the Fourier elements in the image, but the grid sizing will
	// be correct.
	int ixn = int(Util::round(nuxnew));
	int iyn = int(Util::round(nuynew));

	// set up some temporary weighting arrays
	static float wy[7];
	static float wx[7];

	float iynn = nuynew - iyn;
	wy[0] = kb.i0win_tab(iynn+3);
	wy[1] = kb.i0win_tab(iynn+2);
	wy[2] = kb.i0win_tab(iynn+1);
	wy[3] = kb.i0win_tab(iynn);
	wy[4] = kb.i0win_tab(iynn-1);
	wy[5] = kb.i0win_tab(iynn-2);
	wy[6] = kb.i0win_tab(iynn-3);

	float ixnn = nuxnew - ixn;
	wx[0] = kb.i0win_tab(ixnn+3);
	wx[1] = kb.i0win_tab(ixnn+2);
	wx[2] = kb.i0win_tab(ixnn+1);
	wx[3] = kb.i0win_tab(ixnn);
	wx[4] = kb.i0win_tab(ixnn-1);
	wx[5] = kb.i0win_tab(ixnn-2);
	wx[6] = kb.i0win_tab(ixnn-3);

	float wsum = (wx[0]+wx[1]+wx[2]+wx[3]+wx[4]+wx[5]+wx[6])*(wy[0]+wy[1]+wy[2]+wy[3]+wy[4]+wy[5]+wy[6]);

	complex<float> result(0.f,0.f);
        if ((ixn >= 3) && (ixn <= nhalf-3) && (iyn >= -nhalf+3) && (iyn <= nhalf-4)) {
		// (xin,yin) not within window border from the edge
		for (int iy = 0; iy < 7; iy++) {
			int iyp = iyn + iy - 3 ;
			for (int ix = 0; ix < 7; ix++) {
				int ixp = ixn + ix - 3;
				float w = wx[ix]*wy[iy];
				complex<float> val = fimage->cmplx(ixp,iyp);
				result += val*w;
			}
		}
	} else {
		// points that "stick out"
		for (int iy = 0; iy < 7; iy++) {
			int iyp = iyn + iy - 3;
			for (int ix = 0; ix < 7; ix++) {
				int ixp = ixn + ix - 3;
				bool mirror = false;
				int ixt = ixp, iyt = iyp;
				if (ixt < 0) {
					ixt = -ixt;
					iyt = -iyt;
					mirror = !mirror;
				}
				if (ixt > nhalf) {
					ixt = nxreal - ixt;
					iyt = -iyt;
					mirror = !mirror;
				}
				if (iyt > nhalf-1)  iyt -= nxreal;
				if (iyt < -nhalf)   iyt += nxreal;
				float w = wx[ix]*wy[iy];
				complex<float> val = fimage->cmplx(ixt,iyt);
				if (mirror)  result += conj(val)*w;
				else         result += val*w;
			}
		}
	}
	if (flip)  result = conj(result)/wsum;
	else result /= wsum;
	return result;
}

/*
complex<float> Util::extractpoint2(int nx, int ny, float nuxnew, float nuynew, EMData *fimage, Util::KaiserBessel& kb) {

	int nxreal = nx - 2;
	if (nxreal != ny)
		throw ImageDimensionException("extractpoint requires ny == nx");
	int nhalf = nxreal/2;
	bool flip = (nuxnew < 0.f);
	if (flip) {
		nuxnew *= -1;
		nuynew *= -1;
	}
	// put (xnew,ynew) on a grid.  The indices will be wrong for
	// the Fourier elements in the image, but the grid sizing will
	// be correct.
	int ixn = int(Util::round(nuxnew));
	int iyn = int(Util::round(nuynew));
	// set up some temporary weighting arrays
	static float wy[7];
	static float wx[7];

	float iynn = nuynew - iyn;
	wy[0] = kb.i0win_tab(iynn+3);
	wy[1] = kb.i0win_tab(iynn+2);
	wy[2] = kb.i0win_tab(iynn+1);
	wy[3] = kb.i0win_tab(iynn);
	wy[4] = kb.i0win_tab(iynn-1);
	wy[5] = kb.i0win_tab(iynn-2);
	wy[6] = kb.i0win_tab(iynn-3);

	float ixnn = nuxnew - ixn;
	wx[0] = kb.i0win_tab(ixnn+3);
	wx[1] = kb.i0win_tab(ixnn+2);
	wx[2] = kb.i0win_tab(ixnn+1);
	wx[3] = kb.i0win_tab(ixnn);
	wx[4] = kb.i0win_tab(ixnn-1);
	wx[5] = kb.i0win_tab(ixnn-2);
	wx[6] = kb.i0win_tab(ixnn-3);

	float wsum = (wx[0]+wx[1]+wx[2]+wx[3]+wx[4]+wx[5]+wx[6])*(wy[0]+wy[1]+wy[2]+wy[3]+wy[4]+wy[5]+wy[6]);

	complex<float> result(0.f,0.f);

	if ((ixn >= 3) && (ixn <= nhalf-3) && (iyn >= -nhalf+3) && (iyn <= nhalf-4)) {
		// (xin,yin) not within window border from the edge
		result = ( fimage->cmplx(ixn-3,iyn-3)*wx[0] +
			   fimage->cmplx(ixn-2,iyn-3)*wx[1] +
			   fimage->cmplx(ixn-1,iyn-3)*wx[2] +
			   fimage->cmplx(ixn+0,iyn-3)*wx[3] +
			   fimage->cmplx(ixn+1,iyn-3)*wx[4] +
			   fimage->cmplx(ixn+2,iyn-3)*wx[5] +
			   fimage->cmplx(ixn+3,iyn-3)*wx[6] )*wy[0] +
			   ( fimage->cmplx(ixn-3,iyn-2)*wx[0] +
			   fimage->cmplx(ixn-2,iyn-2)*wx[1] +
			   fimage->cmplx(ixn-1,iyn-2)*wx[2] +
			   fimage->cmplx(ixn+0,iyn-2)*wx[3] +
			   fimage->cmplx(ixn+1,iyn-2)*wx[4] +
			   fimage->cmplx(ixn+2,iyn-2)*wx[5] +
			   fimage->cmplx(ixn+3,iyn-2)*wx[6] )*wy[1] +
			   ( fimage->cmplx(ixn-3,iyn-1)*wx[0] +
			   fimage->cmplx(ixn-2,iyn-1)*wx[1] +
			   fimage->cmplx(ixn-1,iyn-1)*wx[2] +
			   fimage->cmplx(ixn+0,iyn-1)*wx[3] +
			   fimage->cmplx(ixn+1,iyn-1)*wx[4] +
			   fimage->cmplx(ixn+2,iyn-1)*wx[5] +
			   fimage->cmplx(ixn+3,iyn-1)*wx[6] )*wy[2] +
			   ( fimage->cmplx(ixn-3,iyn+0)*wx[0] +
			   fimage->cmplx(ixn-2,iyn+0)*wx[1] +
			   fimage->cmplx(ixn-1,iyn+0)*wx[2] +
			   fimage->cmplx(ixn+0,iyn+0)*wx[3] +
			   fimage->cmplx(ixn+1,iyn+0)*wx[4] +
			   fimage->cmplx(ixn+2,iyn+0)*wx[5] +
			   fimage->cmplx(ixn+3,iyn+0)*wx[6] )*wy[3] +
			   ( fimage->cmplx(ixn-3,iyn+1)*wx[0] +
			   fimage->cmplx(ixn-2,iyn+1)*wx[1] +
			   fimage->cmplx(ixn-1,iyn+1)*wx[2] +
			   fimage->cmplx(ixn+0,iyn+1)*wx[3] +
			   fimage->cmplx(ixn+1,iyn+1)*wx[4] +
			   fimage->cmplx(ixn+2,iyn+1)*wx[5] +
			   fimage->cmplx(ixn+3,iyn+1)*wx[6] )*wy[4] +
			   ( fimage->cmplx(ixn-3,iyn+2)*wx[0] +
			   fimage->cmplx(ixn-2,iyn+2)*wx[1] +
			   fimage->cmplx(ixn-1,iyn+2)*wx[2] +
			   fimage->cmplx(ixn+0,iyn+2)*wx[3] +
			   fimage->cmplx(ixn+1,iyn+2)*wx[4] +
			   fimage->cmplx(ixn+2,iyn+2)*wx[5] +
			   fimage->cmplx(ixn+3,iyn+2)*wx[6] )*wy[5] +
			   ( fimage->cmplx(ixn-3,iyn+3)*wx[0] +
			   fimage->cmplx(ixn-2,iyn+3)*wx[1] +
			   fimage->cmplx(ixn-1,iyn+3)*wx[2] +
			   fimage->cmplx(ixn+0,iyn+3)*wx[3] +
			   fimage->cmplx(ixn+1,iyn+3)*wx[4] +
			   fimage->cmplx(ixn+2,iyn+3)*wx[5] +
			   fimage->cmplx(ixn+3,iyn+3)*wx[6] )*wy[6];

	} else {
		// points that "stick out"
		for (int iy = 0; iy < 7; iy++) {
			int iyp = iyn + iy - 3;
			for (int ix = 0; ix < 7; ix++) {
				int ixp = ixn + ix - 3;
				bool mirror = false;
				int ixt= ixp, iyt= iyp;
				if (ixt < 0) {
					ixt = -ixt;
					iyt = -iyt;
					mirror = !mirror;
				}
				if (ixt > nhalf) {
					ixt = nxreal - ixt;
					iyt = -iyt;
					mirror = !mirror;
				}
				if (iyt > nhalf-1)  iyt -= nxreal;
				if (iyt < -nhalf)   iyt += nxreal;
				float w = wx[ix]*wy[iy];
				complex<float> val = fimage->cmplx(ixt,iyt);
				if (mirror)  result += conj(val)*w;
				else         result += val*w;
			}
		}
	}
	if (flip)  result = conj(result)/wsum;
	else result /= wsum;
	return result;
}*/


float Util::triquad(float R, float S, float T, float* fdata)
{

    const float C2 = 0.5f;    //1.0 / 2.0;
    const float C4 = 0.25f;   //1.0 / 4.0;
    const float C8 = 0.125f;  //1.0 / 8.0;

    float  RS   = R * S;
    float  ST   = S * T;
    float  RT   = R * T;
    float  RST  = R * ST;

    float  RSQ  = 1-R*R;
    float  SSQ  = 1-S*S;
    float  TSQ  = 1-T*T;

    float  RM1  = (1-R);
    float  SM1  = (1-S);
    float  TM1  = (1-T);

    float  RP1  = (1+R);
    float  SP1  = (1+S);
    float  TP1  = (1+T);

    float triquad =
    (-C8) * RST * RM1  * SM1  * TM1 * fdata[0] +
	( C4) * ST  * RSQ  * SM1  * TM1 * fdata[1] +
	( C8) * RST * RP1  * SM1  * TM1 * fdata[2] +
	( C4) * RT  * RM1  * SSQ  * TM1 * fdata[3] +
	(-C2) * T   * RSQ  * SSQ  * TM1 * fdata[4] +
	(-C4) * RT  * RP1  * SSQ  * TM1 * fdata[5] +
	( C8) * RST * RM1  * SP1  * TM1 * fdata[6] +
	(-C4) * ST  * RSQ  * SP1  * TM1 * fdata[7] +
	(-C8) * RST * RP1  * SP1  * TM1 * fdata[8] +
//
	( C4) * RS  * RM1  * SM1  * TSQ * fdata[9]  +
	(-C2) * S   * RSQ  * SM1  * TSQ * fdata[10] +
	(-C4) * RS  * RP1  * SM1  * TSQ * fdata[11] +
	(-C2) * R   * RM1  * SSQ  * TSQ * fdata[12] +
	              RSQ  * SSQ  * TSQ * fdata[13] +
	( C2) * R   * RP1  * SSQ  * TSQ * fdata[14] +
	(-C4) * RS  * RM1  * SP1  * TSQ * fdata[15] +
	( C2) * S   * RSQ  * SP1  * TSQ * fdata[16] +
	( C4) * RS  * RP1  * SP1  * TSQ * fdata[17] +
 //
	( C8) * RST * RM1  * SM1  * TP1 * fdata[18] +
	(-C4) * ST  * RSQ  * SM1  * TP1 * fdata[19] +
	(-C8) * RST * RP1  * SM1  * TP1 * fdata[20] +
	(-C4) * RT  * RM1  * SSQ  * TP1 * fdata[21] +
	( C2) * T   * RSQ  * SSQ  * TP1 * fdata[22] +
	( C4) * RT  * RP1  * SSQ  * TP1 * fdata[23] +
	(-C8) * RST * RM1  * SP1  * TP1 * fdata[24] +
	( C4) * ST  * RSQ  * SP1  * TP1 * fdata[25] +
	( C8) * RST * RP1  * SP1  * TP1 * fdata[26]   ;
     return triquad;
}

Util::sincBlackman::sincBlackman(int M_, float fc_, int ntable_)
		: M(M_), fc(fc_), ntable(ntable_) {
	// Sinc-Blackman kernel
	build_sBtable();
}

void Util::sincBlackman::build_sBtable() {
	sBtable.resize(ntable+1);
	int ltab = int(round(float(ntable)/1.25f));
	int M2 = M/2;
	fltb = float(ltab)/M2;
	for (int i=ltab+1; i <= ntable; i++) sBtable[i] = 0.0f;
	float x = 1.0e-7f;
	sBtable[0] = (float)(sin(twopi*fc*x)/x*(0.52-0.5*cos(twopi*(x-M2)/M)+0.08*cos(2*twopi*(x-M2)/M)));
	for (int i=1; i <= ltab; i++) {
		x = float(i)/fltb;
		sBtable[i] = (float)(sin(twopi*fc*x)/x*(0.52-0.5*cos(twopi*(x-M2)/M)+0.08*cos(2*twopi*(x-M2)/M)));
		//cout << "  "<<x<<"  "<<sBtable[i] <<endl;
	}
}

Util::KaiserBessel::KaiserBessel(float alpha_, int K_, float r_, float v_,
		                         int N_, float vtable_, int ntable_)
		: alpha(alpha_), v(v_), r(r_), N(N_), K(K_), vtable(vtable_),
		  ntable(ntable_) {
	// Default values are alpha=1.25, K=6, r=0.5, v = K/2
	if (0.f == v) v = float(K)/2;
	if (0.f == vtable) vtable = v;
	alphar = alpha*r;
	fac = static_cast<float>(twopi)*alphar*v;
	vadjust = 1.0f*v;
	facadj = static_cast<float>(twopi)*alphar*vadjust;
	build_I0table();
}

float Util::KaiserBessel::i0win(float x) const {
	float val0 = float(gsl_sf_bessel_I0(facadj));
	float absx = fabs(x);
	if (absx > vadjust) return 0.f;
	float rt = sqrt(1.f - pow(absx/vadjust, 2));
	float res = static_cast<float>(gsl_sf_bessel_I0(facadj*rt))/val0;
	return res;
}

void Util::KaiserBessel::build_I0table() {
	i0table.resize(ntable+1); // i0table[0:ntable]
	int ltab = int(round(float(ntable)/1.25f));
	fltb = float(ltab)/(K/2);
	float val0 = static_cast<float>(gsl_sf_bessel_I0(facadj));
	for (int i=ltab+1; i <= ntable; i++) i0table[i] = 0.f;
	for (int i=0; i <= ltab; i++) {
		float s = float(i)/fltb/N;
		if (s < vadjust) {
			float rt = sqrt(1.f - pow(s/vadjust, 2));
			i0table[i] = static_cast<float>(gsl_sf_bessel_I0(facadj*rt))/val0;
		} else {
			i0table[i] = 0.f;
		}
//		cout << "  "<<s*N<<"  "<<i0table[i] <<endl;
	}
}

float Util::KaiserBessel::I0table_maxerror() {
	float maxdiff = 0.f;
	for (int i = 1; i <= ntable; i++) {
		float diff = fabs(i0table[i] - i0table[i-1]);
		if (diff > maxdiff) maxdiff = diff;
	}
	return maxdiff;
}

float Util::KaiserBessel::sinhwin(float x) const {
	float val0 = sinh(fac)/fac;
	float absx = fabs(x);
	if (0.0 == x) {
		float res = 1.0f;
		return res;
	} else if (absx == alphar) {
		return 1.0f/val0;
	} else if (absx < alphar) {
		float rt = sqrt(1.0f - pow((x/alphar), 2));
		float facrt = fac*rt;
		float res = (sinh(facrt)/facrt)/val0;
		return res;
	} else {
		float rt = sqrt(pow((x/alphar),2) - 1.f);
		float facrt = fac*rt;
		float res = (sin(facrt)/facrt)/val0;
		return res;
	}
}

float Util::FakeKaiserBessel::i0win(float x) const {
	float val0 = sqrt(facadj)*float(gsl_sf_bessel_I1(facadj));
	float absx = fabs(x);
	if (absx > vadjust) return 0.f;
	float rt = sqrt(1.f - pow(absx/vadjust, 2));
	float res = sqrt(facadj*rt)*float(gsl_sf_bessel_I1(facadj*rt))/val0;
	return res;
}

void Util::FakeKaiserBessel::build_I0table() {
	i0table.resize(ntable+1); // i0table[0:ntable]
	int ltab = int(round(float(ntable)/1.1f));
	fltb = float(ltab)/(K/2);
	float val0 = sqrt(facadj)*static_cast<float>(gsl_sf_bessel_I1(facadj));
	for (int i=ltab+1; i <= ntable; i++) i0table[i] = 0.f;
	for (int i=0; i <= ltab; i++) {
		float s = float(i)/fltb/N;
		if (s < vadjust) {
			float rt = sqrt(1.f - pow(s/vadjust, 2));
			i0table[i] = sqrt(facadj*rt)*static_cast<float>(gsl_sf_bessel_I1(facadj*rt))/val0;
		} else {
			i0table[i] = 0.f;
		}
	}
}

float Util::FakeKaiserBessel::sinhwin(float x) const {
	float val0 = sinh(fac)/fac;
	float absx = fabs(x);
	if (0.0 == x) {
		float res = 1.0f;
		return res;
	} else if (absx == alphar) {
		return 1.0f/val0;
	} else if (absx < alphar) {
		float rt = sqrt(1.0f - pow((x/alphar), 2));
		float facrt = fac*rt;
		float res = (sinh(facrt)/facrt)/val0;
		return res;
	} else {
		float rt = sqrt(pow((x/alphar),2) - 1.f);
		float facrt = fac*rt;
		float res = (sin(facrt)/facrt)/val0;
		return res;
	}
}

#if 0 // 1-st order KB window
float Util::FakeKaiserBessel::sinhwin(float x) const {
	//float val0 = sinh(fac)/fac;
	float prefix = 2*facadj*vadjust/float(gsl_sf_bessel_I1(facadj));
	float val0 = prefix*(cosh(facadj) - sinh(facadj)/facadj);
	float absx = fabs(x);
	if (0.0 == x) {
		//float res = 1.0f;
		float res = val0;
		return res;
	} else if (absx == alphar) {
		//return 1.0f/val0;
		return prefix;
	} else if (absx < alphar) {
		float rt = sqrt(1.0f - pow((x/alphar), 2));
		//float facrt = fac*rt;
		float facrt = facadj*rt;
		//float res = (sinh(facrt)/facrt)/val0;
		float res = prefix*(cosh(facrt) - sinh(facrt)/facrt);
		return res;
	} else {
		float rt = sqrt(pow((x/alphar),2) - 1.f);
		//float facrt = fac*rt;
		float facrt = facadj*rt;
		//float res = (sin(facrt)/facrt)/val0;
		float res = prefix*(sin(facrt)/facrt - cos(facrt));
		return res;
	}
}
#endif // 0



#define  circ(i)         circ[i-1]
#define  numr(i,j)       numr[(j-1)*3 + i-1]
#define  xim(i,j)        xim[(j-1)*nsam + i-1]

EMData* Util::Polar2D(EMData* image, vector<int> numr, string cmode){
	int nsam = image->get_xsize();
	int nrow = image->get_ysize();
	int nring = numr.size()/3;
	int lcirc = numr[3*nring-2]+numr[3*nring-1]-1;
	EMData* out = new EMData();
	out->set_size(lcirc,1,1);
	char mode = (cmode == "F" || cmode == "f") ? 'f' : 'h';
	float *xim  = image->get_data();
	float *circ = out->get_data();
/*   alrq(image->get_data(), nsam, nrow, &numr[0], out->get_data(), lcirc, nring, cmode);
   return out;
}
void Util::alrq(float *xim,  int nsam , int nrow , int *numr,
          float *circ, int lcirc, int nring, char mode)
{*/
/*
c
c  purpose:
c
c  resmaple to polar coordinates
c
*/
	//  dimension	      xim(nsam,nrow),circ(lcirc)
	//  integer	      numr(3,nring)

	double dfi, dpi;
	int    ns2, nr2, i, inr, l, nsim, kcirc, lt, j;
	float  yq, xold, yold, fi, x, y;

	ns2 = nsam/2+1;
	nr2 = nrow/2+1;
	dpi = 2.0*atan(1.0);

	for (i=1;i<=nring;i++) {
		// radius of the ring
		inr = numr(1,i);
		yq  = static_cast<float>(inr);
		l   = numr(3,i);
		if (mode == 'h' || mode == 'H')  lt = l/2;
		else                             lt = l/4;

		nsim	       = lt-1;
		dfi	       = dpi/(nsim+1);
		kcirc	       = numr(2,i);
		xold	       = 0.0f;
		yold	       = static_cast<float>(inr);
		circ(kcirc)    = quadri(xold+(float)ns2,yold+(float)nr2,nsam,nrow,xim);
		xold	       = static_cast<float>(inr);
		yold	       = 0.0f;
		circ(lt+kcirc) = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);

		if (mode == 'f' || mode == 'F') {
			xold		  = 0.0f;
			yold		  = static_cast<float>(-inr);
			circ(lt+lt+kcirc) = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);
			xold		  = static_cast<float>(-inr);
			yold		  = 0.0f;
			circ(lt+lt+lt+kcirc) = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);
		}

		for (j=1;j<=nsim;j++) {
			fi		 = static_cast<float>(dfi*j);
			x		 = sin(fi)*yq;
			y		 = cos(fi)*yq;
			xold		 = x;
			yold		 = y;
			circ(j+kcirc)	 = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);
			xold		 =  y;
			yold		 = -x;
			circ(j+lt+kcirc) = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);

			if (mode == 'f' || mode == 'F')  {
				xold		    = -x;
				yold		    = -y;
				circ(j+lt+lt+kcirc) = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);
				xold		    = -y;
				yold		    =  x;
				circ(j+lt+lt+lt+kcirc) = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);
			}
		}
	}
	return  out;
}
/*
 * 10/22/2014
 * Replacement for the previous version which is commented out after this function
 * Optimized:
 * 	*Sin and Cos functions are tabulated for the largest ring
 * 	*Bilinear interpolation
*/
EMData* Util::Polar2Dm(EMData* image, float cnx2, float cny2, vector<int> numr, string cmode){
	int nring = numr.size()/3;
	int r1 = numr(1,1);
	int r2 = numr(1,nring);
	int maxPoints = numr(3,nring);
	
	int nx = image->get_xsize();
	int ny = image->get_ysize();

	int lcirc = numr[3*nring-2]+numr[3*nring-1]-1;
	float xold, yold, xnew, ynew;
	
	EMData* out = new EMData();
	out->set_size(lcirc,1, 1);

	char mode = (cmode == "F" || cmode == "f") ? 'F' : 'H';
	int div = (mode == 'F' ? 4 : 2);

	maxPoints = maxPoints /div - 1;

	double dpi = 2*atan(1.0), dfi;
	dfi = dpi / (maxPoints +1);
//	Table for sin & cos
	vector<float> vsin(maxPoints);
	vector<float> vcos(maxPoints);
	for (int x = 0; x < maxPoints; x++) {
		float ang = static_cast<float>((x+1) * dfi);
		vsin[x] = sin(ang);
		vcos[x] = cos(ang);
	}

	

	float *xim  = image->get_data();
	float* circ = out->get_data();
	for (int r = 0; r < r2 - r1+1; r++) {

		int it = r+1;
		int inr = numr(1,it);

		int iRef = numr(3,it)/div;

		int kcirc = numr(2,it);

		xold  = 0.0f+cnx2;
		yold  = inr+cny2;
		circ(kcirc) = bilinear(xold,yold,nx,ny,xim);    // Sampling on 90 degree

		xold  = inr+cnx2;
		yold  = 0.0f+cny2;
		circ(iRef+kcirc) = bilinear(xold,yold,nx,ny,xim);  // Sampling on 0 degree

		if (mode == 'F') {
			xold = 0.0f+cnx2;
			yold = -inr+cny2;
			circ(2*iRef+kcirc) = bilinear(xold,yold,nx,ny,xim);  // Sampling on 270 degree

			xold = -inr+cnx2;
			yold = 0.0f+cny2;
			circ(3*iRef+kcirc) = bilinear(xold,yold,nx,ny,xim); // Sampling on 180 degree
		}

		int nPoints = iRef-1;
		int mult = (maxPoints+1)/(nPoints+1);
		for (int x = 0; x < nPoints; x++) {
			int jt = x+1;
			int ind = (x+1)*mult - 1;
			xnew    = vsin[ind] * inr;
			ynew    = vcos[ind] * inr;

			xold = xnew+cnx2;
			yold = ynew+cny2;
			circ(jt+kcirc) = bilinear(xold,yold,nx,ny,xim);      // Sampling on the first quadrant

			xold = ynew+cnx2;
			yold = -xnew+cny2;
			circ(jt+iRef+kcirc) = bilinear(xold,yold,nx,ny,xim);	// Sampling on the fourth quadrant

			if (mode == 'F') {
				xold = -xnew+cnx2;
				yold = -ynew+cny2;

				circ(jt+2*iRef+kcirc) = bilinear(xold,yold,nx,ny,xim); // Sampling on the third quadrant

				xold = -ynew+cnx2;
				yold = xnew+cny2;

				circ(jt+3*iRef+kcirc) = bilinear(xold,yold,nx,ny,xim);  // Sampling on the second quadrant
			}
		}
	}
	return out;
}

/*
 * 10/22/2014
 * Previous version
 * Quadratic interpolation*
 EMData* Util::Polar2Dm(EMData* image, float cns2, float cnr2, vector<int> numr, string cmode){
	int nsam = image->get_xsize();
	int nrow = image->get_ysize();
	int nring = numr.size()/3;
	int lcirc = numr[3*nring-2]+numr[3*nring-1]-1;
	EMData* out = new EMData();
	out->set_size(lcirc,1,1);
	char mode = (cmode == "F" || cmode == "f") ? 'f' : 'h';
	float *xim  = image->get_data();
	float *circ = out->get_data();
	double dpi, dfi;
	int    it, jt, inr, l, nsim, kcirc, lt;
	float  xold, yold, fi, x, y;

	//     cns2 and cnr2 are predefined centers
	//     no need to set to zero, all elements are defined
	dpi = 2*atan(1.0);
	for (it=1; it<=nring; it++) {
		// radius of the ring
		inr = numr(1,it);

		// "F" means a full circle interpolation
		// "H" means a half circle interpolation

		l = numr(3,it);
		if ( mode == 'h' || mode == 'H' ) lt = l / 2;
		else                              lt = l / 4;

		nsim  = lt - 1;
		dfi   = dpi / (nsim+1);
		kcirc = numr(2,it);
		xold  = 0.0f+cns2;
		yold  = inr+cnr2;

//		Assert( kcirc <= lcirc );
		circ(kcirc) = quadri(xold,yold,nsam,nrow,xim);    // Sampling on 90 degree

		xold  = inr+cns2;
		yold  = 0.0f+cnr2;
//		Assert( lt+kcirc <= lcirc );
		circ(lt+kcirc) = quadri(xold,yold,nsam,nrow,xim);  // Sampling on 0 degree

		if ( mode == 'f' || mode == 'F' ) {
			xold = 0.0f+cns2;
			yold = -inr+cnr2;
//			Assert( lt+lt+kcirc <= lcirc );
			circ(lt+lt+kcirc) = quadri(xold,yold,nsam,nrow,xim);  // Sampling on 270 degree

			xold = -inr+cns2;
			yold = 0.0f+cnr2;
//			Assert(lt+lt+lt+kcirc <= lcirc );
			circ(lt+lt+lt+kcirc) = quadri(xold,yold,nsam,nrow,xim); // Sampling on 180 degree
		}

		for (jt=1; jt<=nsim; jt++) {
			fi   = static_cast<float>(dfi * jt);
			x    = sin(fi) * inr;
			y    = cos(fi) * inr;

			xold = x+cns2;
			yold = y+cnr2;

//			Assert( jt+kcirc <= lcirc );
			circ(jt+kcirc) = quadri(xold,yold,nsam,nrow,xim);      // Sampling on the first quadrant

			xold = y+cns2;
			yold = -x+cnr2;

//			Assert( jt+lt+kcirc <= lcirc );
			circ(jt+lt+kcirc) = quadri(xold,yold,nsam,nrow,xim);	// Sampling on the fourth quadrant

			if ( mode == 'f' || mode == 'F' ) {
				xold = -x+cns2;
				yold = -y+cnr2;

//				Assert( jt+lt+lt+kcirc <= lcirc );
				circ(jt+lt+lt+kcirc) = quadri(xold,yold,nsam,nrow,xim); // Sampling on the third quadrant

				xold = -y+cns2;
				yold = x+cnr2;

//				Assert( jt+lt+lt+lt+kcirc <= lcirc );
				circ(jt+lt+lt+lt+kcirc) = quadri(xold,yold,nsam,nrow,xim);  // Sampling on the second quadrant
			}
		} // end for jt
	} //end for it
	return out;
}

*/
float Util::bilinear(float xold, float yold, int nsam, int nrow, float* xim)
{
/*
c  purpose: linear interpolation
  Optimized for speed, circular closer removed, checking of ranges removed
*/
    float bilinear;
    int   ixold, iyold;

/*
	float xdif, ydif, xrem, yrem;
	ixold   = (int) floor(xold);
	iyold   = (int) floor(yold);
	ydif = yold - iyold;
	yrem = 1.0f - ydif;
	xdif = xold - ixold;
	xrem = 1.0f- xdif;
//                 RBUF(K) = YDIF*(BUF(NADDR+NSAM)*XREM
//     &                    +BUF(NADDR+NSAM+1)*XDIF)
//     &                    +YREM*(BUF(NADDR)*XREM + BUF(NADDR+1)*XDIF)
	bilinear = ydif*(xim(ixold,iyold+1)*xrem + xim(ixold+1,iyold+1)*xdif) +
	  				yrem*(xim(ixold,iyold)*xrem+xim(ixold+1,iyold)*xdif);

    return bilinear;
}
*/
	float xdif, ydif;
	ixold   = (int) xold;
	iyold   = (int) yold;
	ydif = yold - iyold;
	xdif = xold - ixold;
	bilinear = xim(ixold, iyold) + ydif* (xim(ixold, iyold+1) - xim(ixold, iyold)) +
	           xdif* (xim(ixold+1, iyold) - xim(ixold, iyold) +
			   ydif* (xim(ixold+1, iyold+1) - xim(ixold+1, iyold) - xim(ixold, iyold+1) + xim(ixold, iyold)) );

	return bilinear;
}

void Util::alrl_ms(float *xim, int    nsam, int  nrow, float cns2, float cnr2,
             int  *numr, float *circ, int , int  nring, char  mode) {
	double dpi, dfi;
	int    it, jt, inr, l, nsim, kcirc, lt;
	float  xold, yold, fi, x, y;

	//     cns2 and cnr2 are predefined centers
	//     no need to set to zero, all elements are defined

	dpi = 2*atan(1.0);
	for (it=1; it<=nring; it++) {
		// radius of the ring
		inr = numr(1,it);

		l = numr(3,it);
		if ( mode == 'h' || mode == 'H' ) lt = l / 2;
		else				  lt = l / 4;

		nsim  = lt - 1;
		dfi   = dpi / (nsim+1);
		kcirc = numr(2,it);


		xold  = 0.0f+cns2;
		yold  = inr+cnr2;

		circ(kcirc) = quadri(xold,yold,nsam,nrow,xim);

		xold  = inr+cns2;
		yold  = 0.0f+cnr2;
		circ(lt+kcirc) = quadri(xold,yold,nsam,nrow,xim);

		if ( mode == 'f' || mode == 'F' ) {
			xold = 0.0f+cns2;
			yold = -inr+cnr2;
			circ(lt+lt+kcirc) = quadri(xold,yold,nsam,nrow,xim);

			xold = -inr+cns2;
			yold = 0.0f+cnr2;
			circ(lt+lt+lt+kcirc) = quadri(xold,yold,nsam,nrow,xim);
		}

		for (jt=1; jt<=nsim; jt++) {
			fi   = static_cast<float>(dfi * jt);
			x    = sin(fi) * inr;
			y    = cos(fi) * inr;

			xold = x+cns2;
			yold = y+cnr2;
			circ(jt+kcirc) = quadri(xold,yold,nsam,nrow,xim);

			xold = y+cns2;
			yold = -x+cnr2;
			circ(jt+lt+kcirc) = quadri(xold,yold,nsam,nrow,xim);

			if ( mode == 'f' || mode == 'F' ) {
				xold = -x+cns2;
				yold = -y+cnr2;
				circ(jt+lt+lt+kcirc) = quadri(xold,yold,nsam,nrow,xim);

				xold = -y+cns2;
				yold = x+cnr2;
				circ(jt+lt+lt+lt+kcirc) = quadri(xold,yold,nsam,nrow,xim);
			}
		} // end for jt
	} //end for it
}
/*
void Util::alrl_ms(float *xim, int    nsam, int  nrow, float cns2, float cnr2,
             int  *numr, float *circ, int lcirc, int  nring, char  mode)
{
   double dpi, dfi;
   int    it, jt, inr, l, nsim, kcirc, lt, xold, yold;
   float  yq, fi, x, y;

   //     cns2 and cnr2 are predefined centers
   //     no need to set to zero, all elements are defined

   dpi = 2*atan(1.0);
   for (it=1; it<=nring; it++) {
      // radius of the ring
      inr = numr(1,it);
      yq  = inr;

      l = numr(3,it);
      if ( mode == 'h' || mode == 'H' ) {
         lt = l / 2;
      }
      else { // if ( mode == 'f' || mode == 'F' )
         lt = l / 4;
      }

      nsim  = lt - 1;
      dfi   = dpi / (nsim+1);
      kcirc = numr(2,it);


	xold = (int) (0.0+cns2);
	yold = (int) (inr+cnr2);

	circ(kcirc) = xim(xold, yold);

      xold = (int) (inr+cns2);
      yold = (int) (0.0+cnr2);
      circ(lt+kcirc) = xim(xold, yold);

      if ( mode == 'f' || mode == 'F' ) {
         xold  = (int) (0.0+cns2);
         yold = (int) (-inr+cnr2);
         circ(lt+lt+kcirc) = xim(xold, yold);

         xold  = (int) (-inr+cns2);
         yold = (int) (0.0+cnr2);
         circ(lt+lt+lt+kcirc) = xim(xold, yold);
      }

      for (jt=1; jt<=nsim; jt++) {
         fi   = dfi * jt;
         x    = sin(fi) * yq;
         y    = cos(fi) * yq;

         xold  = (int) (x+cns2);
         yold = (int) (y+cnr2);
         circ(jt+kcirc) = xim(xold, yold);

         xold  = (int) (y+cns2);
         yold = (int) (-x+cnr2);
         circ(jt+lt+kcirc) = xim(xold, yold);

         if ( mode == 'f' || mode == 'F' ) {
            xold  = (int) (-x+cns2);
            yold = (int) (-y+cnr2);
            circ(jt+lt+lt+kcirc) = xim(xold, yold);

            xold  = (int) (-y+cns2);
            yold = (int) (x+cnr2);
            circ(jt+lt+lt+lt+kcirc) = xim(xold, yold);
         }
      } // end for jt
   } //end for it
}
*/
//xim((int) floor(xold), (int) floor(yold))
#undef  xim

EMData* Util::Polar2Dmi(EMData* image, float cns2, float cnr2, vector<int> numr, string cmode, Util::KaiserBessel& kb){
// input image is twice the size of the original image
	int nring = numr.size()/3;
	int lcirc = numr[3*nring-2]+numr[3*nring-1]-1;
	EMData* out = new EMData();
	out->set_size(lcirc,1,1);
	char mode = (cmode == "F" || cmode == "f") ? 'f' : 'h';
	float *circ = out->get_data();
	float *fimage = image->get_data();
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	double dpi, dfi;
	int    it, jt, inr, l, nsim, kcirc, lt;
	float  yq, xold, yold, fi, x, y;

	//     cns2 and cnr2 are predefined centers
	//     no need to set to zero, all elements are defined

	dpi = 2*atan(1.0);
	for (it=1;it<=nring;it++) {
		// radius of the ring
		inr = numr(1,it);
		yq  = static_cast<float>(inr);

		l = numr(3,it);
		if ( mode == 'h' || mode == 'H' )  lt = l / 2;
		else                               lt = l / 4;

		nsim  = lt - 1;
		dfi   = dpi / (nsim+1);
		kcirc = numr(2,it);
		xold  = 0.0f;
		yold  = static_cast<float>(inr);
		circ(kcirc) = get_pixel_conv_new(nx,ny,nz,2*(xold+cns2-1.0f),2*(yold+cnr2-1.0f),0,fimage,kb);
//      circ(kcirc) = image->get_pixel_conv(2*(xold+cns2-1.0f),2*(yold+cnr2-1.0f),0,kb);

		xold  = static_cast<float>(inr);
		yold  = 0.0f;
		circ(lt+kcirc) = get_pixel_conv_new(nx,ny,nz,2*(xold+cns2-1.0f),2*(yold+cnr2-1.0f),0,fimage,kb);
//      circ(lt+kcirc) = image->get_pixel_conv(2*(xold+cns2-1.0f),2*(yold+cnr2-1.0f),0,kb);

	if ( mode == 'f' || mode == 'F' ) {
		xold = 0.0f;
		yold = static_cast<float>(-inr);
		circ(lt+lt+kcirc) = get_pixel_conv_new(nx,ny,nz,2*(xold+cns2-1.0f),2*(yold+cnr2-1.0f),0,fimage,kb);
//         circ(lt+lt+kcirc) = image->get_pixel_conv(2*(xold+cns2-1.0f),2*(yold+cnr2-1.0f),0,kb);

		xold = static_cast<float>(-inr);
		yold = 0.0f;
		circ(lt+lt+lt+kcirc) = get_pixel_conv_new(nx,ny,nz,2*(xold+cns2-1.0f),2*(yold+cnr2-1.0f),0,fimage,kb);
//         circ(lt+lt+lt+kcirc) = image->get_pixel_conv(2*(xold+cns2-1.0f),2*(yold+cnr2-1.0f),0,kb);
	}

	for (jt=1;jt<=nsim;jt++) {
		fi   = static_cast<float>(dfi * jt);
		x    = sin(fi) * yq;
		y    = cos(fi) * yq;

		xold = x;
		yold = y;
		circ(jt+kcirc) = get_pixel_conv_new(nx,ny,nz,2*(xold+cns2-1.0f),2*(yold+cnr2-1.0f),0,fimage,kb);
//         circ(jt+kcirc) = image->get_pixel_conv(2*(xold+cns2-1.0f),2*(yold+cnr2-1.0f),0,kb);

		xold = y;
		yold = -x;
		circ(jt+lt+kcirc) = get_pixel_conv_new(nx,ny,nz,2*(xold+cns2-1.0f),2*(yold+cnr2-1.0f),0,fimage,kb);
//         circ(jt+lt+kcirc) = image->get_pixel_conv(2*(xold+cns2-1.0f),2*(yold+cnr2-1.0f),0,kb);

	if ( mode == 'f' || mode == 'F' ) {
		xold = -x;
		yold = -y;
		circ(jt+lt+lt+kcirc) = get_pixel_conv_new(nx,ny,nz,2*(xold+cns2-1.0f),2*(yold+cnr2-1.0f),0,fimage,kb);
//            circ(jt+lt+lt+kcirc) = image->get_pixel_conv(2*(xold+cns2-1.0f),2*(yold+cnr2-1.0f),0,kb);

		xold = -y;
		yold = x;
		circ(jt+lt+lt+lt+kcirc) = get_pixel_conv_new(nx,ny,nz,2*(xold+cns2-1.0f),2*(yold+cnr2-1.0f),0,fimage,kb);
//            circ(jt+lt+lt+lt+kcirc) = image->get_pixel_conv(2*(xold+cns2-1.0f),2*(yold+cnr2-1.0f),0,kb);
	}
	} // end for jt
	} //end for it
	return  out;
}

/*

        A set of 1-D power-of-two FFTs
	Pawel & Chao 01/20/06

fftr_q(xcmplx,nv)
  single precision

 dimension xcmplx(2,iabs(nv)/2);
 xcmplx(1,1) --- R(0), xcmplx(2,1) --- R(NV/2)
 xcmplx(1,i) --- real, xcmplx(2,i) --- imaginary


fftr_d(xcmplx,nv)
  double precision

 dimension xcmplx(2,iabs(nv)/2);
 xcmplx(1,1) --- R(0), xcmplx(2,1) --- R(NV/2)
 xcmplx(1,i) --- real, xcmplx(2,i) --- imaginary



*/
#define  tab1(i)      tab1[i-1]
#define  xcmplx(i,j)  xcmplx [(j-1)*2 + i-1]
#define  br(i)        br[i-1]
#define  bi(i)        bi[i-1]
//-----------------------------------------
void Util::fftc_d(double *br, double *bi, int ln, int ks)
{
	double rni,sgn,tr1,tr2,ti1,ti2;
	double cc,c,ss,s,t,x2,x3,x4,x5;
	int    b3,b4,b5,b6,b7,b56;
	int    n, k, l, j, i, ix0, ix1, status=0;

	const double tab1[] = {
		9.58737990959775e-5,
		1.91747597310703e-4,
		3.83495187571395e-4,
		7.66990318742704e-4,
		1.53398018628476e-3,
		3.06795676296598e-3,
		6.13588464915449e-3,
		1.22715382857199e-2,
		2.45412285229123e-2,
		4.90676743274181e-2,
		9.80171403295604e-2,
		1.95090322016128e-1,
		3.82683432365090e-1,
		7.07106781186546e-1,
		1.00000000000000,
	};

	n=(int)pow(2.0f,ln);

	k=abs(ks);
	l=16-ln;
	b3=n*k;
	b6=b3;
	b7=k;
	if (ks > 0) {
		sgn=1.0f;
	} else {
		sgn=-1.0f;
		rni=1.0f/(float)(n);
		j=1;
		for (i=1; i<=n; i++) {
			br(j)=br(j)*rni;
			bi(j)=bi(j)*rni;
			j=j+k;
		}
	}

L12:
   b6=b6/2;
   b5=b6;
   b4=2*b6;
   b56=b5-b6;

L14:
   tr1=br(b5+1);
   ti1=bi(b5+1);
   tr2=br(b56+1);
   ti2=bi(b56+1);

   br(b5+1)=tr2-tr1;
   bi(b5+1)=ti2-ti1;
   br(b56+1)=tr1+tr2;
   bi(b56+1)=ti1+ti2;

   b5=b5+b4;
   b56=b5-b6;
   if ( b5 <= b3 )  goto  L14;
   if ( b6 == b7 )  goto  L20;

   b4=b7;
   cc=2.0f*pow(tab1(l),2);
   c=1.0f-cc;
   l++;
   ss=sgn*tab1(l);
   s=ss;

L16:
   b5=b6+b4;
   b4=2*b6;
   b56=b5-b6;

L18:
   tr1=br(b5+1);
   ti1=bi(b5+1);
   tr2=br(b56+1);
   ti2=bi(b56+1);
   br(b5+1)=c*(tr2-tr1)-s*(ti2-ti1);
   bi(b5+1)=s*(tr2-tr1)+c*(ti2-ti1);
   br(b56+1)=tr1+tr2;
   bi(b56+1)=ti1+ti2;

   b5=b5+b4;
   b56=b5-b6;
   if ( b5 <= b3 )  goto  L18;
   b4=b5-b6;
   b5=b4-b3;
   c=-c;
   b4=b6-b5;
   if ( b5 < b4 )  goto  L16;
   b4=b4+b7;
   if ( b4 >= b5 ) goto  L12;

   t=c-cc*c-ss*s;
   s=s+ss*c-cc*s;
   c=t;
   goto  L16;

L20:
   ix0=b3/2;
   b3=b3-b7;
   b4=0;
   b5=0;
   b6=ix0;
   ix1=0;
   if (b6 == b7) goto EXIT;

L22:
   b4=b3-b4;
   b5=b3-b5;
   x2=br(b4+1);
   x3=br(b5+1);
   x4=bi(b4+1);
   x5=bi(b5+1);
   br(b4+1)=x3;
   br(b5+1)=x2;
   bi(b4+1)=x5;
   bi(b5+1)=x4;
   if(b6 < b4)  goto  L22;

L24:
   b4=b4+b7;
   b5=b6+b5;
   x2=br(b4+1);
   x3=br(b5+1);
   x4=bi(b4+1);
   x5=bi(b5+1);
   br(b4+1)=x3;
   br(b5+1)=x2;
   bi(b4+1)=x5;
   bi(b5+1)=x4;
   ix0=b6;

L26:
   ix0=ix0/2;
   ix1=ix1-ix0;
   if( ix1 >= 0)  goto L26;

   ix0=2*ix0;
   b4=b4+b7;
   ix1=ix1+ix0;
   b5=ix1;
   if ( b5 >= b4)  goto  L22;
   if ( b4 < b6)   goto  L24;

EXIT:
   status = 0;
}

// -----------------------------------------------------------------
void Util::fftc_q(float *br, float *bi, int ln, int ks)
{
	//  dimension  br(1),bi(1)

	int b3,b4,b5,b6,b7,b56;
	int n, k, l, j, i, ix0, ix1;
	float rni, tr1, ti1, tr2, ti2, cc, c, ss, s, t, x2, x3, x4, x5, sgn;
	int status=0;

	const float tab1[] = {
		9.58737990959775e-5f,
		1.91747597310703e-4f,
		3.83495187571395e-4f,
		7.66990318742704e-4f,
		1.53398018628476e-3f,
		3.06795676296598e-3f,
		6.13588464915449e-3f,
		1.22715382857199e-2f,
		2.45412285229123e-2f,
		4.90676743274181e-2f,
		9.80171403295604e-2f,
		1.95090322016128e-1f,
		3.82683432365090e-1f,
		7.07106781186546e-1f,
		1.00000000000000f,
	};

	n=(int)pow(2.0f,ln);

	k=abs(ks);
	l=16-ln;
	b3=n*k;
	b6=b3;
	b7=k;
	if( ks > 0 ) {
		sgn=1.0f;
	} else {
		sgn=-1.0f;
		rni=1.0f/(float)n;
		j=1;
		for (i=1; i<=n; i++) {
			br(j)=br(j)*rni;
			bi(j)=bi(j)*rni;
			j=j+k;
		}
	}
L12:
   b6=b6/2;
   b5=b6;
   b4=2*b6;
   b56=b5-b6;
L14:
   tr1=br(b5+1);
   ti1=bi(b5+1);

   tr2=br(b56+1);
   ti2=bi(b56+1);

   br(b5+1)=tr2-tr1;
   bi(b5+1)=ti2-ti1;
   br(b56+1)=tr1+tr2;
   bi(b56+1)=ti1+ti2;

   b5=b5+b4;
   b56=b5-b6;
   if ( b5 <= b3 )  goto  L14;
   if ( b6 == b7 )  goto  L20;

   b4=b7;
   cc=2.0f*pow(tab1(l),2);
   c=1.0f-cc;
   l++;
   ss=sgn*tab1(l);
   s=ss;
L16:
   b5=b6+b4;
   b4=2*b6;
   b56=b5-b6;
L18:
   tr1=br(b5+1);
   ti1=bi(b5+1);
   tr2=br(b56+1);
   ti2=bi(b56+1);
   br(b5+1)=c*(tr2-tr1)-s*(ti2-ti1);
   bi(b5+1)=s*(tr2-tr1)+c*(ti2-ti1);
   br(b56+1)=tr1+tr2;
   bi(b56+1)=ti1+ti2;

   b5=b5+b4;
   b56=b5-b6;
   if(b5 <= b3)  goto L18;
   b4=b5-b6;
   b5=b4-b3;
   c=-c;
   b4=b6-b5;
   if(b5 < b4)  goto  L16;
   b4=b4+b7;
   if(b4 >= b5) goto  L12;

   t=c-cc*c-ss*s;
   s=s+ss*c-cc*s;
   c=t;
   goto  L16;
L20:
   ix0=b3/2;
   b3=b3-b7;
   b4=0;
   b5=0;
   b6=ix0;
   ix1=0;
   if ( b6 == b7) goto EXIT;
L22:
   b4=b3-b4;
   b5=b3-b5;
   x2=br(b4+1);
   x3=br(b5+1);
   x4=bi(b4+1);
   x5=bi(b5+1);
   br(b4+1)=x3;
   br(b5+1)=x2;
   bi(b4+1)=x5;
   bi(b5+1)=x4;
   if (b6 < b4) goto  L22;
L24:
   b4=b4+b7;
   b5=b6+b5;
   x2=br(b4+1);
   x3=br(b5+1);
   x4=bi(b4+1);
   x5=bi(b5+1);
   br(b4+1)=x3;
   br(b5+1)=x2;
   bi(b4+1)=x5;
   bi(b5+1)=x4;
   ix0=b6;
L26:
   ix0=ix0/2;
   ix1=ix1-ix0;
   if(ix1 >= 0)  goto  L26;

   ix0=2*ix0;
   b4=b4+b7;
   ix1=ix1+ix0;
   b5=ix1;
   if (b5 >= b4)  goto  L22;
   if (b4 < b6)   goto  L24;
EXIT:
   status = 0;
}

void  Util::fftr_q(float *xcmplx, int nv)
{
   // dimension xcmplx(2,1); xcmplx(1,i) --- real, xcmplx(2,i) --- imaginary

	int nu, inv, nu1, n, isub, n2, i1, i2, i;
	float ss, cc, c, s, tr, ti, tr1, tr2, ti1, ti2, t;

	const float tab1[] = {
		9.58737990959775e-5f,
		1.91747597310703e-4f,
		3.83495187571395e-4f,
		7.66990318742704e-4f,
		1.53398018628476e-3f,
		3.06795676296598e-3f,
		6.13588464915449e-3f,
		1.22715382857199e-2f,
		2.45412285229123e-2f,
		4.90676743274181e-2f,
		9.80171403295604e-2f,
		1.95090322016128e-1f,
		3.82683432365090e-1f,
		7.07106781186546e-1f,
		1.00000000000000f,
	};

	nu=abs(nv);
	inv=nv/nu;
	nu1=nu-1;
	n=(int)pow(2.f,nu1);
	isub=16-nu1;

	ss=-tab1(isub);
	cc=-2.0f*pow(tab1(isub-1),2.f);
	c=1.0f;
	s=0.0f;
	n2=n/2;
	if ( inv > 0) {
		fftc_q(&xcmplx(1,1),&xcmplx(2,1),nu1,2);
		tr=xcmplx(1,1);
		ti=xcmplx(2,1);
		xcmplx(1,1)=tr+ti;
		xcmplx(2,1)=tr-ti;
		for (i=1;i<=n2;i++) {
			i1=i+1;
			i2=n-i+1;
			tr1=xcmplx(1,i1);
			tr2=xcmplx(1,i2);
			ti1=xcmplx(2,i1);
			ti2=xcmplx(2,i2);
			t=(cc*c-ss*s)+c;
			s=(cc*s+ss*c)+s;
			c=t;
			xcmplx(1,i1)=0.5f*((tr1+tr2)+(ti1+ti2)*c-(tr1-tr2)*s);
			xcmplx(1,i2)=0.5f*((tr1+tr2)-(ti1+ti2)*c+(tr1-tr2)*s);
			xcmplx(2,i1)=0.5f*((ti1-ti2)-(ti1+ti2)*s-(tr1-tr2)*c);
			xcmplx(2,i2)=0.5f*(-(ti1-ti2)-(ti1+ti2)*s-(tr1-tr2)*c);
		}
	} else {
		tr=xcmplx(1,1);
		ti=xcmplx(2,1);
		xcmplx(1,1)=0.5f*(tr+ti);
		xcmplx(2,1)=0.5f*(tr-ti);
		for (i=1; i<=n2; i++) {
			i1=i+1;
			i2=n-i+1;
			tr1=xcmplx(1,i1);
			tr2=xcmplx(1,i2);
			ti1=xcmplx(2,i1);
			ti2=xcmplx(2,i2);
			t=(cc*c-ss*s)+c;
			s=(cc*s+ss*c)+s;
			c=t;
			xcmplx(1,i1)=0.5f*((tr1+tr2)-(tr1-tr2)*s-(ti1+ti2)*c);
			xcmplx(1,i2)=0.5f*((tr1+tr2)+(tr1-tr2)*s+(ti1+ti2)*c);
			xcmplx(2,i1)=0.5f*((ti1-ti2)+(tr1-tr2)*c-(ti1+ti2)*s);
			xcmplx(2,i2)=0.5f*(-(ti1-ti2)+(tr1-tr2)*c-(ti1+ti2)*s);
		}
		fftc_q(&xcmplx(1,1),&xcmplx(2,1),nu1,-2);
	}
}

// -------------------------------------------
void  Util::fftr_d(double *xcmplx, int nv)
{
	// double precision  x(2,1)
	int    i1, i2,  nu, inv, nu1, n, isub, n2, i;
	double tr1,tr2,ti1,ti2,tr,ti;
	double cc,c,ss,s,t;
	const double tab1[] = {
		9.58737990959775e-5,
		1.91747597310703e-4,
		3.83495187571395e-4,
		7.66990318742704e-4,
		1.53398018628476e-3,
		3.06795676296598e-3,
		6.13588464915449e-3,
		1.22715382857199e-2,
		2.45412285229123e-2,
		4.90676743274181e-2,
		9.80171403295604e-2,
		1.95090322016128e-1,
		3.82683432365090e-1,
		7.07106781186546e-1,
		1.00000000000000,
	};

	nu=abs(nv);
	inv=nv/nu;
	nu1=nu-1;
	n=(int)pow(2.0f,nu1);
	isub=16-nu1;
	ss=-tab1(isub);
	cc=-2.0*pow(tab1(isub-1),2);
	c=1.0f;
	s=0.0f;
	n2=n/2;

	if ( inv > 0 ) {
		fftc_d(&xcmplx(1,1),&xcmplx(2,1),nu1,2);
		tr=xcmplx(1,1);
		ti=xcmplx(2,1);
		xcmplx(1,1)=tr+ti;
		xcmplx(2,1)=tr-ti;
		for (i=1;i<=n2;i++) {
			i1=i+1;
			i2=n-i+1;
			tr1=xcmplx(1,i1);
			tr2=xcmplx(1,i2);
			ti1=xcmplx(2,i1);
			ti2=xcmplx(2,i2);
			t=(cc*c-ss*s)+c;
			s=(cc*s+ss*c)+s;
			c=t;
			xcmplx(1,i1)=0.5*((tr1+tr2)+(ti1+ti2)*c-(tr1-tr2)*s);
			xcmplx(1,i2)=0.5*((tr1+tr2)-(ti1+ti2)*c+(tr1-tr2)*s);
			xcmplx(2,i1)=0.5*((ti1-ti2)-(ti1+ti2)*s-(tr1-tr2)*c);
			xcmplx(2,i2)=0.5*(-(ti1-ti2)-(ti1+ti2)*s-(tr1-tr2)*c);
		}
	} else {
		tr=xcmplx(1,1);
		ti=xcmplx(2,1);
		xcmplx(1,1)=0.5*(tr+ti);
		xcmplx(2,1)=0.5*(tr-ti);
		for (i=1; i<=n2; i++) {
			i1=i+1;
			i2=n-i+1;
			tr1=xcmplx(1,i1);
			tr2=xcmplx(1,i2);
			ti1=xcmplx(2,i1);
			ti2=xcmplx(2,i2);
			t=(cc*c-ss*s)+c;
			s=(cc*s+ss*c)+s;
			c=t;
			xcmplx(1,i1)=0.5*((tr1+tr2)-(tr1-tr2)*s-(ti1+ti2)*c);
			xcmplx(1,i2)=0.5*((tr1+tr2)+(tr1-tr2)*s+(ti1+ti2)*c);
			xcmplx(2,i1)=0.5*((ti1-ti2)+(tr1-tr2)*c-(ti1+ti2)*s);
			xcmplx(2,i2)=0.5*(-(ti1-ti2)+(tr1-tr2)*c-(ti1+ti2)*s);
		}
		fftc_d(&xcmplx(1,1),&xcmplx(2,1),nu1,-2);
	}
}
#undef  tab1
#undef  xcmplx
#undef  br
#undef  bi


void Util::Frngs(EMData* circp, vector<int> numr){
	int nring = numr.size()/3;
	float *circ = circp->get_data();
	int i, l;
	for (i=1; i<=nring;i++) {

#ifdef _WIN32
		l = (int)( log((float)numr(3,i))/log(2.0f) );
#else
		l=(int)(log2(numr(3,i)));
#endif	//_WIN32

		fftr_q(&circ(numr(2,i)),l);
	}
}

void Util::Frngs_inv(EMData* circp, vector<int> numr){
	int nring = numr.size()/3;
	float *circ = circp->get_data();
	int i, l;
	for (i=1; i<=nring;i++) {

#ifdef _WIN32
		l = (int)( log((float)numr(3,i))/log(2.0f) );
#else
		l=(int)(log2(numr(3,i)));
#endif	//_WIN32

		fftr_q(&circ(numr(2,i)),-l);
	}
}
#undef  circ

void Util::Applyws(EMData* circp, vector<int> numr, vector<float> wr)
{	/*
	  Apply weights to FTs of rings
	*/
	const int nring = numr.size() / 3;
	const int maxrin = numr.back();
	float *circ = circp->get_data();
	for (int i = 0; i < nring; ++i) {
		const int numr3i = numr[2+i*3];
		const int numr2i = numr[1+i*3]-1;
		const float w = wr[i];
		circ[numr2i] *= w;
		if (numr3i == maxrin)  circ[numr2i+1] *= w;
		else                   circ[numr2i+1] *= 0.5*w;
		for (int j = 2+numr2i; j < numr3i+numr2i; ++j)  circ[j] *= w;
	}
}

void Util::prb3p(double *b, float *pos) {
	//  three points 1D parabola fit,  call this function for random searchea
	double  c2,c3;
	c3 = b[0] - 2.0f*b[1] +b[1];
	if( fabs(c3) > 1.0e-6 ) *pos =  max( min ( static_cast<float>(2*(b[0]-b[2])/c3), 1.0f), -1.0f);
	else  *pos = 0.0f;
	//printf("\n  prb1d    %f    %f    %f    %f\n",b[0] ,b[1] ,b[2], *pos);
}



#define  b(i)            b[i-1]
void Util::prb1d(double *b, int npoint, float *pos) {
	double  c2,c3;
	int 	nhalf;

	nhalf = npoint/2 + 1;
	*pos  = 0.0;

	if (npoint == 7) {
		c2 = 49.*b(1) + 6.*b(2) - 21.*b(3) - 32.*b(4) - 27.*b(5)
		     - 6.*b(6) + 31.*b(7);
		c3 = 5.*b(1) - 3.*b(3) - 4.*b(4) - 3.*b(5) + 5.*b(7);
	}
	else if (npoint == 5) {
		c2 = (74.*b(1) - 23.*b(2) - 60.*b(3) - 37.*b(4)
		   + 46.*b(5) ) / (-70.);
		c3 = (2.*b(1) - b(2) - 2.*b(3) - b(4) + 2.*b(5) ) / 14.0;
	}
	else if (npoint == 3) {
		c2 = (5.*b(1) - 8.*b(2) + 3.*b(3) ) / (-2.0);
		c3 = (b(1) - 2.*b(2) + b(3) ) / 2.0;
	}
	//else if (npoint == 9) {
	else  { // at least one has to be true!!
		c2 = (1708.*b(1) + 581.*b(2) - 246.*b(3) - 773.*b(4)
		     - 1000.*b(5) - 927.*b(6) - 554.*b(7) + 119.*b(8)
		     + 1092.*b(9) ) / (-4620.);
		c3 = (28.*b(1) + 7.*b(2) - 8.*b(3) - 17.*b(4) - 20.*b(5)
		     - 17.*b(6) - 8.*b(7) + 7.*b(8) + 28.*b(9) ) / 924.0;
	}
	if (c3 != 0.0)  *pos = static_cast<float>(c2/(2.0*c3) - nhalf);
}
#undef  b


#define  circ1(i)        circ1[i-1]
#define  circ2(i)        circ2[i-1]
#define  t(i)            t[i-1]
#define  q(i)            q[i-1]
#define  b(i)            b[i-1]
#define  t7(i)           t7[i-1]
Dict Util::Crosrng_e(EMData*  circ1p, EMData* circ2p, vector<int> numr, int neg) {
	//  neg = 0 straight,  neg = 1 mirrored
	int nring = numr.size()/3;
	//int lcirc = numr[3*nring-2]+numr[3*nring-1]-1;
	int maxrin = numr[numr.size()-1];
	double qn;   float  tot;
	float *circ1 = circ1p->get_data();
	float *circ2 = circ2p->get_data();
/*
c checks single position, neg is flag for checking mirrored position
c
c  input - fourier transforms of rings!
c  first set is conjugated (mirrored) if neg
c  circ1 already multiplied by weights!
c       automatic arrays
	dimension         t(maxrin)  removed +2 as it is only needed for other ffts
	double precision  q(maxrin)
	double precision  t7(-3:3)
*/
	float *t;
	double t7[7], *q;
	int    i, j, k, ip, jc, numr3i, numr2i, jtot = 0;
	float  pos;

#ifdef _WIN32
	ip = -(int)(log((float)maxrin)/log(2.0f));
#else
	ip = -(int) (log2(maxrin));
#endif	//_WIN32

	q = (double*)calloc(maxrin, sizeof(double));

//   cout << *qn <<"  " <<*tot<<"  "<<ip<<endl;
	for (i=1; i<=nring; i++) {
		numr3i = numr(3,i);
		numr2i = numr(2,i);

		q(1) += (circ1(numr2i)) * circ2(numr2i);

		if (numr3i != maxrin) {
			 // test .ne. first for speed on some compilers
			q(numr3i+1) += circ1(numr2i+1) * circ2(numr2i+1);

			if (neg) {
				// first set is conjugated (mirrored)
				for (j=3;j<=numr3i;j=j+2) {
					jc = j+numr2i-1;
					q(j)   +=  (circ1(jc))*circ2(jc)   - (circ1(jc+1))*circ2(jc+1);
					q(j+1) += -(circ1(jc))*circ2(jc+1) - (circ1(jc+1))*circ2(jc);
				}
			} else {
				for (j=3;j<=numr3i;j=j+2) {
					jc = j+numr2i-1;
					q(j)   +=  (circ1(jc))*circ2(jc)   + (circ1(jc+1))*circ2(jc+1);
					q(j+1) += -(circ1(jc))*circ2(jc+1) + (circ1(jc+1))*circ2(jc);
				}
			}
		} else {
			q(2) += circ1(numr2i+1) * circ2(numr2i+1);
			if (neg) {
				// first set is conjugated (mirrored)
				for (j=3;j<=maxrin;j=j+2) {
					jc = j+numr2i-1;
					q(j)   +=  (circ1(jc))*circ2(jc)   - (circ1(jc+1))*circ2(jc+1);
					q(j+1) += -(circ1(jc))*circ2(jc+1) - (circ1(jc+1))*circ2(jc);
				}
			} else {
				for (j=3;j<=maxrin;j=j+2) {
					jc = j+numr2i-1;
					q(j)   +=  (circ1(jc))*circ2(jc)   + (circ1(jc+1))*circ2(jc+1);
					q(j+1) += -(circ1(jc))*circ2(jc+1) + (circ1(jc+1))*circ2(jc);
				}
			}
		}
	}

	fftr_d(q,ip);

	qn = -1.0e20;
	for (j=1;j<=maxrin;j++) {
	   if (q(j) >= qn) {
		  qn = q(j); jtot = j;
	   }
	}
	for (k=-3; k<=3; k++) {
		j = (jtot+k+maxrin-1)%maxrin + 1;
		t7(k+4) = q(j);
	}

	prb1d(t7,7,&pos);

	tot = (float)jtot + pos;

	if (q) free(q);

	Dict retvals;
	retvals["qn"]  = qn;
	retvals["tot"] = tot;
	return  retvals;
}

static std::vector<int> shuffled_range(int first, int last)
{
	std::vector<int> v(last - first + 1);
	for (unsigned i = 0; i < v.size(); ++i) v[i] = first + i;
	for (unsigned i = 0; i < v.size(); ++i) {
		unsigned r = Util::get_irand(0,v.size()-1);
		std::swap( v[r], v[i] );
	}
	return v;
}

Dict Util::Crosrng_rand_e(EMData*  circ1p, EMData* circ2p, vector<int> numr, int neg, float previous_max, float an, int psi_pos) {
/*
c checks single position, neg is flag for checking mirrored position
c
c  input - fourier transforms of rings!
c  first set is conjugated (mirrored) if neg
c  circ1 already multiplied by weights!
c       automatic arrays
c   neg = 0 straight,  neg = 1 mirrored
	double precision  q(maxrin)
	double precision  t7(-3:3)
*/
	int nring = numr.size()/3;
	//int lcirc = numr[3*nring-2]+numr[3*nring-1]-1;
	int maxrin = numr[numr.size()-1];
	float *circ1 = circ1p->get_data();
	float *circ2 = circ2p->get_data();
	double p3[3], *q;
	int    i, j, k, ip, jc, numr3i, numr2i;
	float  pos;

#ifdef _WIN32
	ip = -(int)(log((float)maxrin)/log(2.0f));
#else
	ip = -(int) (log2(maxrin));
#endif	//_WIN32

	q = (double*)calloc(maxrin, sizeof(double));

	for (i=1; i<=nring; i++) {
		numr3i = numr(3,i);
		numr2i = numr(2,i);

		q(1) += (circ1(numr2i)) * circ2(numr2i);

		if (numr3i != maxrin) {
			 // test .ne. first for speed on some compilers
			q(numr3i+1) += circ1(numr2i+1) * circ2(numr2i+1);

			if (neg) {
				// first set is conjugated (mirrored)
				for (j=3;j<=numr3i;j=j+2) {
					jc = j+numr2i-1;
					q(j)   +=  (circ1(jc))*circ2(jc)   - (circ1(jc+1))*circ2(jc+1);
					q(j+1) += -(circ1(jc))*circ2(jc+1) - (circ1(jc+1))*circ2(jc);
				}
			} else {
				for (j=3;j<=numr3i;j=j+2) {
					jc = j+numr2i-1;
					q(j)   +=  (circ1(jc))*circ2(jc)   + (circ1(jc+1))*circ2(jc+1);
					q(j+1) += -(circ1(jc))*circ2(jc+1) + (circ1(jc+1))*circ2(jc);
				}
			}
		} else {
			q(2) += circ1(numr2i+1) * circ2(numr2i+1);
			if (neg) {
				// first set is conjugated (mirrored)
				for (j=3;j<=maxrin;j=j+2) {
					jc = j+numr2i-1;
					q(j)   +=  (circ1(jc))*circ2(jc)   - (circ1(jc+1))*circ2(jc+1);
					q(j+1) += -(circ1(jc))*circ2(jc+1) - (circ1(jc+1))*circ2(jc);
				}
			} else {
				for (j=3;j<=maxrin;j=j+2) {
					jc = j+numr2i-1;
					q(j)   +=  (circ1(jc))*circ2(jc)   + (circ1(jc+1))*circ2(jc+1);
					q(j+1) += -(circ1(jc))*circ2(jc+1) + (circ1(jc+1))*circ2(jc);
				}
			}
		}
	}

	fftr_d(q,ip);

	// psi is randomized within +- 2 * an of the anchor position
	const int psi_range = int(2*an/360.0*maxrin + 0.5);
    std::vector<int> rperm(2*psi_range+1);
	for (k=-psi_range; k<=psi_range; k++) {
		j = ( k + psi_pos + maxrin - 1)%maxrin+1;
		rperm[k + psi_range] = j;
		//printf(" JJJ  %d      %d      %d      %d \n",k,rperm[k + psi_range],psi_pos,psi_range);
	}

	for (unsigned i = 0; i < rperm.size(); ++i) {
		unsigned r = Util::get_irand(0,rperm.size()-1);
		std::swap( rperm[r], rperm[i] );
	}

	double qn = -1.0e20;
	int jtot  = -1;
	for (int j_j=0;j_j<rperm.size();j_j++) {
		j = rperm[j_j] + 1;
		if (q(j) > previous_max) {
			qn = q(j);
			jtot = j;
			break;
		}
	}
			//printf(" jtot qn  %d      %f \n",jtot, qn);
	float  tot = 0.0;
	if( jtot > -1 ) {
		for (k=-1; k<=1; k++) {
			j = (jtot+k+maxrin-1)%maxrin + 1;
			p3[k+1] = q(j);
		}
		prb3p(p3,&pos);
		tot = (float)jtot + pos;
	}

	if (q) free(q);
			//printf(" tot pos  %f      %f \n",tot, pos);

	Dict retvals;
	retvals["qn"]  = qn;
	retvals["tot"] = tot;
	return  retvals;
}

Dict Util::Crosrng_ew(EMData*  circ1p, EMData* circ2p, vector<int> numr, vector<float> w, int neg) {
   //  neg = 0 straight,  neg = 1 mirrored
	int nring = numr.size()/3;
	//int lcirc = numr[3*nring-2]+numr[3*nring-1]-1;
	int maxrin = numr[numr.size()-1];
	double qn;   float  tot;
	float *circ1 = circ1p->get_data();
	float *circ2 = circ2p->get_data();
/*
c checks single position, neg is flag for checking mirrored position
c
c  input - fourier transforms of rings!
c  first set is conjugated (mirrored) if neg
c  multiplication by weights!
c       automatic arrays
	dimension         t(maxrin)  removed +2 as it is only needed for other ffts
	double precision  q(maxrin)
	double precision  t7(-3:3)
*/
	float *t;
	double t7[7], *q;
	int    i, j, k, ip, jc, numr3i, numr2i, jtot = 0;
	float  pos;

#ifdef _WIN32
	ip = -(int)(log((float)maxrin)/log(2.0f));
#else
	ip = -(int) (log2(maxrin));
#endif	//_WIN32

	q = (double*)calloc(maxrin, sizeof(double));
	t = (float*)calloc(maxrin, sizeof(float));

//   cout << *qn <<"  " <<*tot<<"  "<<ip<<endl;
	for (i=1;i<=nring;i++) {
		numr3i = numr(3,i);
		numr2i = numr(2,i);

		t(1) = circ1(numr2i) * circ2(numr2i);

		if (numr3i != maxrin) {
			// test .ne. first for speed on some compilers
			t(numr3i+1) = circ1(numr2i+1) * circ2(numr2i+1);
			t(2)  	  = 0.0;

			if (neg) {
				// first set is conjugated (mirrored)
				for (j=3; j<=numr3i; j=j+2) {
					jc = j+numr2i-1;
					t(j)   =  (circ1(jc))*circ2(jc)   - (circ1(jc+1))*circ2(jc+1);
					t(j+1) = -(circ1(jc))*circ2(jc+1) - (circ1(jc+1))*circ2(jc);
				}
			} else {
				for (j=3; j<=numr3i; j=j+2) {
					jc     = j+numr2i-1;
					t(j)   =  (circ1(jc))*circ2(jc)   + (circ1(jc+1))*circ2(jc+1);
					t(j+1) = -(circ1(jc))*circ2(jc+1) + (circ1(jc+1))*circ2(jc);
				}
			}
			for (j=1;j<=numr3i+1;j++) q(j) += t(j)*w[i-1];
		} else {
			t(2) = circ1(numr2i+1) * circ2(numr2i+1);
			if (neg) {
				// first set is conjugated (mirrored)
				for (j=3; j<=maxrin; j=j+2) {
					jc = j+numr2i-1;
					t(j)   =  (circ1(jc))*circ2(jc)   - (circ1(jc+1))*circ2(jc+1);
					t(j+1) = -(circ1(jc))*circ2(jc+1) - (circ1(jc+1))*circ2(jc);
				}
			} else {
				for (j=3; j<=maxrin; j=j+2) {
			    	jc = j+numr2i-1;
			    	t(j)   =  (circ1(jc))*circ2(jc)   + (circ1(jc+1))*circ2(jc+1);
			    	t(j+1) = -(circ1(jc))*circ2(jc+1) + (circ1(jc+1))*circ2(jc);
				}
			}
			for (j = 1; j <= maxrin; j++) q(j) += t(j)*w[i-1];
		}
	}

	fftr_d(q,ip);

	qn = -1.0e20;
	for (j=1;j<=maxrin;j++) {
		//cout << j << "  " << q(j) << endl;
		if (q(j) >= qn) {
			qn = q(j);
			jtot = j;
		}
	}

	for (k=-3; k<=3; k++) {
		j = (jtot+k+maxrin-1)%maxrin + 1;
		t7(k+4) = q(j);
	}

	prb1d(t7,7,&pos);

	tot = (float)jtot + pos;

	//if (q) free(q);
	if (t) free(t);

	Dict retvals;
	//tot = 1;
	//qn = q(1);
	retvals["qn"] = qn;
	retvals["tot"] = tot;

	if (q) free(q);

	return  retvals;
}

Dict Util::Crosrng_ms(EMData* circ1p, EMData* circ2p, vector<int> numr) {
	int nring = numr.size()/3;
	//int lcirc = numr[3*nring-2]+numr[3*nring-1]-1;
	int maxrin = numr[numr.size()-1];
	double qn; float tot; double qm; float tmt;
	float *circ1 = circ1p->get_data();
	float *circ2 = circ2p->get_data();
/*
c
c  checks both straight & mirrored positions
c
c  input - fourier transforms of rings!!
c  circ1 already multiplied by weights!
c
*/

	// dimension		 circ1(lcirc),circ2(lcirc)

	// t(maxrin), q(maxrin), t7(-3:3)  //maxrin+2 removed
	double *t, *q, t7[7];

	int   ip, jc, numr3i, numr2i, i, j, k, jtot = 0;
	float t1, t2, t3, t4, c1, c2, d1, d2, pos;

	qn  = 0.0f;
	qm  = 0.0f;
	tot = 0.0f;
	tmt = 0.0f;
#ifdef _WIN32
	ip = -(int)(log((float)maxrin)/log(2.0f));
#else
 	ip = -(int)(log2(maxrin));
#endif	//_WIN32
  //for (j=1; j<=maxrin;j++) cout <<"  "<<j<<"   "<<circ1(j)<<"   "<<circ2(j) <<endl;

	//  c - straight  = circ1 * conjg(circ2)
	//  zero q array

	q = (double*)calloc(maxrin,sizeof(double));

	//   t - mirrored  = conjg(circ1) * conjg(circ2)
	//   zero t array
	t = (double*)calloc(maxrin,sizeof(double));

   //   premultiply  arrays ie( circ12 = circ1 * circ2) much slower
	for (i=1; i<=nring; i++) {

		numr3i = numr(3,i);   // Number of samples of this ring
		numr2i = numr(2,i);   // The beginning point of this ring

		t1   = circ1(numr2i) * circ2(numr2i);
		q(1) += t1;
		t(1) += t1;

		t1   = circ1(numr2i+1) * circ2(numr2i+1);
		if (numr3i == maxrin)  {
			q(2) += t1;
			t(2) += t1;
		} else {
			q(numr3i+1) += t1;
			t(numr3i+1) += t1;
		}

		for (j=3; j<=numr3i; j += 2) {
			jc     = j+numr2i-1;

// Here, (c1+c2i)*conj(d1+d2i) = (c1*d1+c2*d2)+(-c1*d2+c2*d1)i
//   			          ----- -----    ----- -----
//      			   t1     t2      t3    t4
// Here, conj(c1+c2i)*conj(d1+d2i) = (c1*d1-c2*d2)+(-c1*d2-c2*d1)i
//     		                      ----- -----    ----- -----
//     			               t1    t2       t3    t4

			c1     = circ1(jc);
			c2     = circ1(jc+1);
			d1     = circ2(jc);
			d2     = circ2(jc+1);

			t1     = c1 * d1;
			t2     = c2 * d2;
			t3     = c1 * d2;
			t4     = c2 * d1;

			q(j)   += t1 + t2;
			q(j+1) += -t3 + t4;
			t(j)   += t1 - t2;
			t(j+1) += -t3 - t4;
		}
	}
	//for (j=1; j<=maxrin; j++) cout <<"  "<<j<<"   "<<q(j) <<"   "<<t(j) <<endl;
	fftr_d(q,ip);

	qn  = -1.0e20;
	for (j=1; j<=maxrin; j++) {//cout <<"  "<<j<<"   "<<q(j) <<endl;
		if (q(j) >= qn) {
			qn  = q(j);
			jtot = j;
		}
	}

	for (k=-3; k<=3; k++) {
		j = ((jtot+k+maxrin-1)%maxrin)+1;
		t7(k+4) = q(j);
	}

	// interpolate
	prb1d(t7,7,&pos);
	tot = (float)(jtot)+pos;

	// mirrored
	fftr_d(t,ip);

	// find angle
	qm = -1.0e20;
	for (j=1; j<=maxrin;j++) {//cout <<"  "<<j<<"	"<<t(j) <<endl;
		if ( t(j) >= qm ) {
			qm   = t(j);
			jtot = j;
		}
	}

	for (k=-3; k<=3; k++)  {
		j = ((jtot+k+maxrin-1)%maxrin) + 1;
		t7(k+4) = t(j);
	}

	// interpolate

	prb1d(t7,7,&pos);
	tmt = float(jtot) + pos;

	free(t);
	free(q);

	Dict retvals;
	retvals["qn"] = qn;
	retvals["tot"] = tot;
	retvals["qm"] = qm;
	retvals["tmt"] = tmt;
	return retvals;
}

static Dict Crosrng_rand_ms(EMData* circ1p, EMData* circ2p, vector<int> numr, float previous_max)
{
	int nring = numr.size()/3;
	//int lcirc = numr[3*nring-2]+numr[3*nring-1]-1;
	int maxrin = numr[numr.size()-1];
	float *circ1 = circ1p->get_data();
	float *circ2 = circ2p->get_data();
/*
c
c  checks both straight & mirrored positions
c
c  input - fourier transforms of rings!!
c  circ1 already multiplied by weights!
c
*/

	// dimension		 circ1(lcirc),circ2(lcirc)
	double *t, *q, p3[3];

	int   jc, numr3i, numr2i, i, j, k;
	float t1, t2, t3, t4, c1, c2, d1, d2, pos;

#ifdef _WIN32
	const int ip = -(int)(log((float)maxrin)/log(2.0f));
#else
	const int ip = -(int)(log2(maxrin));
#endif	//_WIN32
  //for (j=1; j<=maxrin;j++) cout <<"  "<<j<<"   "<<circ1(j)<<"   "<<circ2(j) <<endl;

	//  c - straight  = circ1 * conjg(circ2)
	//  zero q array

	q = (double*)calloc(maxrin,sizeof(double));

	//   t - mirrored  = conjg(circ1) * conjg(circ2)
	//   zero t array
	t = (double*)calloc(maxrin,sizeof(double));

   //   premultiply  arrays ie( circ12 = circ1 * circ2) much slower
	for (i=1; i<=nring; i++) {

		numr3i = numr(3,i);   // Number of samples of this ring
		numr2i = numr(2,i);   // The beginning point of this ring

		t1   = circ1(numr2i) * circ2(numr2i);
		q(1) += t1;
		t(1) += t1;

		t1   = circ1(numr2i+1) * circ2(numr2i+1);
		if (numr3i == maxrin)  {
			q(2) += t1;
			t(2) += t1;
		} else {
			q(numr3i+1) += t1;
			t(numr3i+1) += t1;
		}

		for (j=3; j<=numr3i; j += 2) {
			jc     = j+numr2i-1;

// Here, (c1+c2i)*conj(d1+d2i) = (c1*d1+c2*d2)+(-c1*d2+c2*d1)i
//   			                  ----- -----    ----- -----
//      			               t1     t2      t3    t4
// Here, conj(c1+c2i)*conj(d1+d2i) = (c1*d1-c2*d2)+(-c1*d2-c2*d1)i
//     		                          ----- -----    ----- -----
//     			                       t1    t2       t3    t4

			c1     = circ1(jc);
			c2     = circ1(jc+1);
			d1     = circ2(jc);
			d2     = circ2(jc+1);

			t1     = c1 * d1;
			t2     = c2 * d2;
			t3     = c1 * d2;
			t4     = c2 * d1;

			q(j)   += t1 + t2;
			q(j+1) += -t3 + t4;
			t(j)   += t1 - t2;
			t(j+1) += -t3 - t4;
		}
	}
	//for (j=1; j<=maxrin; j++) cout <<"  "<<j<<"   "<<q(j) <<"   "<<t(j) <<endl;
	Util::fftr_d(q,ip);
	// mirrored
	Util::fftr_d(t,ip);

	std::vector<int> angles = shuffled_range(0, 2*maxrin-1);  // [0..maxrin) - not-mirrored, [maxrin..2*maxrin) - mirrored
	bool mirrored = false;
	double qn = -1.0e20;
	int    jtot = -1;
	for (int pos_j=0; pos_j < 2*maxrin; pos_j++) {   //cout <<"  "<<j<<"   "<<q(j-1) <<endl;
		j = angles[pos_j] + 1;
		if (j <= maxrin) {
			// not mirrored
			if (q(j) > previous_max) {
				mirrored = false;
				qn  = q(j);
				jtot = j;
				break;
			}
		} else {
			// mirrored
			j -= maxrin;
			if ( t(j) > previous_max ) {
				mirrored = true;
				qn   = t(j);
				jtot = j;
				break;
			}
		}
	}

	float tot = 0.0;
	if( jtot > -1) {
		// interpolation
		if (! mirrored) {
			for (k=-1; k<=1; k++) {
				j = ((jtot+k+maxrin-1)%maxrin)+1;
				p3[k+1] = q(j);
			}
		} else {
			for (k=-1; k<=1; k++)  {
				j = ((jtot+k+maxrin-1)%maxrin) + 1;
				p3[k+1] = t(j);
			}
		}

		// interpolate
		Util::prb3p(p3,&pos);
		tot = (float)(jtot)+pos;
	}

	free(t);
	free(q);

	Dict retvals;
	retvals["qn"] = qn;
	retvals["tot"] = tot;
	retvals["mirror"] = (mirrored) ? (1) : (0);
	return retvals;
}

Dict Util::Crosrng_ms_delta(EMData* circ1p, EMData* circ2p, vector<int> numr, float delta_start, float delta) {
	int nring = numr.size()/3;
	//int lcirc = numr[3*nring-2]+numr[3*nring-1]-1;
	int maxrin = numr[numr.size()-1];
	double qn; float tot; double qm; float tmt;
	float *circ1 = circ1p->get_data();
	float *circ2 = circ2p->get_data();
/*
c
c  checks both straight & mirrored positions
c
c  input - fourier transforms of rings!!
c  circ1 already multiplied by weights!
c
*/

	// dimension		 circ1(lcirc),circ2(lcirc)

	// t(maxrin), q(maxrin), t7(-3:3)  //maxrin+2 removed
	double *t, *q;

	int   ip, jc, numr3i, numr2i, i, j, jtot = 0;
	float t1, t2, t3, t4, c1, c2, d1, d2;

	qn  = 0.0f;
	qm  = 0.0f;
	tot = 0.0f;
	tmt = 0.0f;
#ifdef _WIN32
	ip = -(int)(log((float)maxrin)/log(2.0f));
#else
 	ip = -(int)(log2(maxrin));
#endif	//_WIN32
  //for (j=1; j<=maxrin;j++) cout <<"  "<<j<<"   "<<circ1(j)<<"   "<<circ2(j) <<endl;

	//  c - straight  = circ1 * conjg(circ2)
	//  zero q array

	q = (double*)calloc(maxrin,sizeof(double));

	//   t - mirrored  = conjg(circ1) * conjg(circ2)
	//   zero t array
	t = (double*)calloc(maxrin,sizeof(double));

   //   premultiply  arrays ie( circ12 = circ1 * circ2) much slower
	for (i=1; i<=nring; i++) {

		numr3i = numr(3,i);   // Number of samples of this ring
		numr2i = numr(2,i);   // The beginning point of this ring

		t1   = circ1(numr2i) * circ2(numr2i);
		q(1) += t1;
		t(1) += t1;

		t1   = circ1(numr2i+1) * circ2(numr2i+1);
		if (numr3i == maxrin)  {
			q(2) += t1;
			t(2) += t1;
		} else {
			q(numr3i+1) += t1;
			t(numr3i+1) += t1;
		}

		for (j=3; j<=numr3i; j += 2) {
			jc     = j+numr2i-1;

// Here, (c1+c2i)*conj(d1+d2i) = (c1*d1+c2*d2)+(-c1*d2+c2*d1)i
//   			          ----- -----    ----- -----
//      			   t1     t2      t3    t4
// Here, conj(c1+c2i)*conj(d1+d2i) = (c1*d1-c2*d2)+(-c1*d2-c2*d1)i
//     		                      ----- -----    ----- -----
//     			               t1    t2       t3    t4

			c1     = circ1(jc);
			c2     = circ1(jc+1);
			d1     = circ2(jc);
			d2     = circ2(jc+1);

			t1     = c1 * d1;
			t2     = c2 * d2;
			t3     = c1 * d2;
			t4     = c2 * d1;

			q(j)   += t1 + t2;
			q(j+1) += -t3 + t4;
			t(j)   += t1 - t2;
			t(j+1) += -t3 - t4;
		}
	}
	//for (j=1; j<=maxrin; j++) cout <<"  "<<j<<"   "<<q(j) <<"   "<<t(j) <<endl;
	fftr_d(q,ip);

	qn  = -1.0e20;

	int jstart = 1+static_cast<int>(delta_start/360.0*maxrin);
	int jstep = static_cast<int>(delta/360.0*maxrin);
	if (jstep < 1) { jstep = 1; }

	for (j=jstart; j<=maxrin; j+=jstep) {//cout <<"  "<<j<<"   "<<q(j) <<endl;
		if (q(j) >= qn) {
			qn  = q(j);
			jtot = j;
		}
	}

	//for (k=-3; k<=3; k++) {
	//	j = ((jtot+k+maxrin-1)%maxrin)+1;
	//	t7(k+4) = q(j);
	//}

	// interpolate
	//prb1d(t7,7,&pos);
	//tot = (float)(jtot)+pos;
	// Do not interpolate
	tot = (float)(jtot);

	// mirrored
	fftr_d(t,ip);

	// find angle
	qm = -1.0e20;
	for (j=jstart; j<=maxrin;j+=jstep) {//cout <<"  "<<j<<"	"<<t(j) <<endl;
		if ( t(j) >= qm ) {
			qm   = t(j);
			jtot = j;
		}
	}

	//for (k=-3; k<=3; k++)  {
	//	j = ((jtot+k+maxrin-1)%maxrin) + 1;
	//	t7(k+4) = t(j);
	//}

	// interpolate

	//prb1d(t7,7,&pos);
	//tmt = float(jtot) + pos;
	// Do not interpolate
	tmt = float(jtot);

	free(t);
	free(q);

	Dict retvals;
	retvals["qn"] = qn;
	retvals["tot"] = tot;
	retvals["qm"] = qm;
	retvals["tmt"] = tmt;
	return retvals;
}


Dict Util::Crosrng_sm_psi(EMData* circ1p, EMData* circ2p, vector<int> numr, float psi, int flag, float psi_max) {
// flag 0 - straight, 1 - mirror

	int nring = numr.size()/3;
	int maxrin = numr[numr.size()-1];
	double qn; float tot;
	float *circ1 = circ1p->get_data();
	float *circ2 = circ2p->get_data();

	double *q;

	int   ip, jc, numr3i, numr2i, i, j, k, jtot = 0;
	float c1, c2, d1, d2;

	qn  = 0.0f;
	tot = 0.0f;
#ifdef _WIN32
	ip = -(int)(log((float)maxrin)/log(2.0f));
#else
 	ip = -(int)(log2(maxrin));
#endif	//_WIN32

	//  c - straight  = circ1 * conjg(circ2)
	//  zero q array

	q = (double*)calloc(maxrin,sizeof(double));
	int neg = 1-2*flag;
   //   premultiply  arrays ie( circ12 = circ1 * circ2) much slower

	for (i=1; i<=nring; i++) {

		numr3i = numr(3,i);   // Number of samples of this ring
		numr2i = numr(2,i);   // The beginning point of this ring

		q(1) += circ1(numr2i) * circ2(numr2i);

		float t1 = circ1(numr2i+1) * circ2(numr2i+1);
		if (numr3i == maxrin)   q(2) += t1;
		else  q(numr3i+1) += t1;

		for (j=3; j<=numr3i; j += 2) {
			jc     = j+numr2i-1;

	// Here, (c1+c2i)*conj(d1+d2i) = (c1*d1+c2*d2)+(-c1*d2+c2*d1)i
	//   			          ----- -----    ----- -----
	//      			   t1     t2      t3    t4

			c1     = circ1(jc);
			c2     = circ1(jc+1);
			d1     = circ2(jc);
			d2     = circ2(jc+1);

			q(j)   +=  c1 * d1 + c2 * d2*neg;
			q(j+1) += -c1 * d2 + c2 * d1*neg;
		}
	}

	fftr_d(q,ip);

	qn  = -1.0e20;
	// psi = 0 should correspond to psi_pos = 1 (meaning no change in in-plane rotation)
	// this has to be verified for mirror
	int psi_pos = (int)fmod((float)Util::round(-psi/360.0*maxrin+10*maxrin),(float)maxrin) + 1;

	const int psi_range = int(psi_max/360.0*maxrin + 0.5);

	for (k=-psi_range; k<=psi_range; k++) {
		j = ( k + psi_pos + maxrin - 1)%maxrin+1;
		if (q(j) >= qn) {
			qn  = q(j);
			jtot = j;
		}
	}

	tot = (float)(jtot);
	free(q);

	Dict retvals;
	retvals["qn"]  = qn;
	retvals["tot"] = tot;
	return retvals;
}

Dict Util::Crosrng_psi(EMData* circ1p, EMData* circ2p, vector<int> numr, float psi, float psi_max) {
// Computes both straight and mirrored

	int nring = numr.size()/3;
	int maxrin = numr[numr.size()-1];
	double qn; float tot; double qm; float tmt;
	float *circ1 = circ1p->get_data();
	float *circ2 = circ2p->get_data();

	double *t, *q;

	int   ip, jc, numr3i, numr2i, i, j, k, jtot = 0;
	float t1, t2, t3, t4, c1, c2, d1, d2;

	qn  = 0.0f;
	qm  = 0.0f;
	tot = 0.0f;
	tmt = 0.0f;
#ifdef _WIN32
	ip = -(int)(log((float)maxrin)/log(2.0f));
#else
 	ip = -(int)(log2(maxrin));
#endif	//_WIN32

	//  c - straight  = circ1 * conjg(circ2)
	//  zero q array

	q = (double*)calloc(maxrin,sizeof(double));
	
	//   t - mirrored  = conjg(circ1) * conjg(circ2)
	//   zero t array
	t = (double*)calloc(maxrin,sizeof(double));
	
   //   premultiply  arrays ie( circ12 = circ1 * circ2) much slower
   	
	for (i=1; i<=nring; i++) {

		numr3i = numr(3,i);   // Number of samples of this ring
		numr2i = numr(2,i);   // The beginning point of this ring

		t1   = circ1(numr2i) * circ2(numr2i);
		q(1) += t1;
		t(1) += t1;
		
		t1   = circ1(numr2i+1) * circ2(numr2i+1);
		if (numr3i == maxrin)  {
			q(2) += t1;
			t(2) += t1;
		} else {
			q(numr3i+1) += t1;
			t(numr3i+1) += t1;
		}

		for (j=3; j<=numr3i; j += 2) {
			jc     = j+numr2i-1;

	// Here, (c1+c2i)*conj(d1+d2i) = (c1*d1+c2*d2)+(-c1*d2+c2*d1)i
	//   			          ----- -----    ----- -----
	//      			   t1     t2      t3    t4

			c1     = circ1(jc);
			c2     = circ1(jc+1);
			d1     = circ2(jc);
			d2     = circ2(jc+1);

			t1     = c1 * d1;
			t3     = c1 * d2;
			t2     = c2 * d2;
			t4     = c2 * d1;

			q(j)   +=  t1 + t2;
			q(j+1) += -t3 + t4;
			t(j)   +=  t1 - t2;
			t(j+1) += -t3 - t4;
		}
	}

	fftr_d(q,ip);

	qn  = -1.0e20;
	// psi = 0 should correspond to psi_pos = 1
	// this has to be verified for mirror
	int psi_pos = (int)fmod((float)Util::round(-psi/360.0*maxrin+10*maxrin),(float)maxrin) + 1;
	const int psi_range = int(psi_max/360.0*maxrin + 0.5);

	for (k=-psi_range; k<=psi_range; k++) {
    j = (k+psi_pos+maxrin-1)%maxrin+1;
		if (q(j) >= qn) {
			qn  = q(j);
			jtot = j;
		}
	}

	tot = (float)(jtot);
	free(q);

    // mirrored
	fftr_d(t,ip);

	qm  = -1.0e20;

	for (k=-psi_range; k<=psi_range; k++) {
		j = (k+psi_pos+maxrin-1)%maxrin+1;
		if (t(j) >= qm) {
			qm  = t(j);
			jtot = j;
		}
	}

	tmt = (float)(jtot);
	free(t);

	Dict retvals;
	retvals["qn"]  = qn;
	retvals["tot"] = tot;
	retvals["qm"]  = qm;
	retvals["tmt"] = tmt;
	return retvals;
}

Dict Util::Crosrng_ns(EMData* circ1p, EMData* circ2p, vector<int> numr) {
	int nring = numr.size()/3;
	//int lcirc = numr[3*nring-2]+numr[3*nring-1]-1;
	int maxrin = numr[numr.size()-1];
	double qn; float tot;
	float *circ1 = circ1p->get_data();
	float *circ2 = circ2p->get_data();
/*
c
c  checks only straight position
c
c  input - fourier transforms of rings!!
c  circ1 already multiplied by weights!
c
*/

	// dimension		 circ1(lcirc),circ2(lcirc)

	// q(maxrin), t7(-3:3)  //maxrin+2 removed
	double *q, t7[7];

	int   ip, jc, numr3i, numr2i, i, j, k, jtot = 0;
	float c1, c2, d1, d2, pos;

	qn  = 0.0;
	tot = 0.0;
#ifdef _WIN32
	ip = -(int)(log((float)maxrin)/log(2.0f));
#else
   ip = -(int)(log2(maxrin));
#endif	//_WIN32
	//for (j=1; j<=maxrin;j++) cout <<"  "<<j<<"   "<<circ1(j)<<"	"<<circ2(j) <<endl;

	//  c - straight  = circ1 * conjg(circ2)
	//  zero q array

	q = (double*)calloc(maxrin,sizeof(double));

			//   premultiply  arrays ie( circ12 = circ1 * circ2) much slower
	for (i=1; i<=nring; i++) {

		numr3i = numr(3,i);   // Number of samples of this ring
		numr2i = numr(2,i);   // The beginning point of this ring

		q(1) += circ1(numr2i) * circ2(numr2i);

		if (numr3i == maxrin)   q(2) += circ1(numr2i+1) * circ2(numr2i+1);
		else                 q(numr3i+1) += circ1(numr2i+1) * circ2(numr2i+1);

		for (j=3; j<=numr3i; j += 2) {
			jc     = j+numr2i-1;

// Here, (c1+c2i)*conj(d1+d2i) = (c1*d1+c2*d2)+(-c1*d2+c2*d1)i
//   			          ----- -----    ----- -----
//      			   t1     t2      t3    t4

			c1     = circ1(jc);
			c2     = circ1(jc+1);
			d1     = circ2(jc);
			d2     = circ2(jc+1);

			q(j)   +=  c1 * d1 + c2 * d2;
			q(j+1) += -c1 * d2 + c2 * d1;
		}
	}
//for (j=1; j<=maxrin; j++) cout <<"  "<<j<<"   "<<q(j) <<endl;
	fftr_d(q,ip);

	qn  = -1.0e20;
	for (j=1; j<=maxrin; j++) {//cout <<"  "<<j<<"   "<<q(j) <<endl;
		if (q(j) >= qn) {
			qn  = q(j);
			jtot = j;
		}
	}

	for (k=-3; k<=3; k++)  {
		j = ((jtot+k+maxrin-1)%maxrin)+1;
		t7(k+4) = q(j);
	}

	// interpolate
	prb1d(t7,7,&pos);
	tot = (float)(jtot)+pos;

	free(q);

	Dict retvals;
	retvals["qn"]  = qn;
	retvals["tot"] = tot;
	return retvals;
}

#define  dout(i,j)        dout[i+maxrin*j]
#define  circ1b(i)        circ1b[i-1]
#define  circ2b(i)        circ2b[i-1]

EMData* Util::Crosrng_msg(EMData* circ1, EMData* circ2, vector<int> numr) {

   // dimension         circ1(lcirc),circ2(lcirc)

	int   ip, jc, numr3i, numr2i, i, j;
	float t1, t2, t3, t4, c1, c2, d1, d2;

	int nring = numr.size()/3;
	//int lcirc = numr[3*nring-2]+numr[3*nring-1]-1;
	int maxrin = numr[numr.size()-1];

	float* circ1b = circ1->get_data();
	float* circ2b = circ2->get_data();

	// t(maxrin), q(maxrin)  // removed +2
	double *t, *q;

	q = (double*)calloc(maxrin,sizeof(double));
	t = (double*)calloc(maxrin,sizeof(double));

#ifdef _WIN32
	ip = -(int)(log((float)maxrin)/log(2.0f));
#else
	ip = -(int)(log2(maxrin));
#endif	//_WIN32

	//  q - straight  = circ1 * conjg(circ2)

	//   t - mirrored  = conjg(circ1) * conjg(circ2)

	//   premultiply  arrays ie( circ12 = circ1 * circ2) much slower

	for (i=1; i<=nring; i++) {

		numr3i = numr(3,i);
		numr2i = numr(2,i);

		t1   = circ1b(numr2i) * circ2b(numr2i);
		q(1) += t1;
		t(1) += t1;

		t1   = circ1b(numr2i+1) * circ2b(numr2i+1);
		if (numr3i == maxrin)  {
			q(2) += t1;
			t(2) += t1;
		} else {
			q(numr3i+1) += t1;
			t(numr3i+1) += t1;
		}

		for (j=3; j<=numr3i; j=j+2) {
			jc     = j+numr2i-1;

			c1     = circ1b(jc);
			c2     = circ1b(jc+1);
			d1     = circ2b(jc);
			d2     = circ2b(jc+1);

			t1     = c1 * d1;
			t3     = c1 * d2;
			t2     = c2 * d2;
			t4     = c2 * d1;

			q(j)   +=  t1 + t2;
			q(j+1) += -t3 + t4;
			t(j)   +=  t1 - t2;
			t(j+1) += -t3 - t4;
		}
	}

	// straight
	fftr_d(q,ip);

	// mirrored
	fftr_d(t,ip);

	EMData* out = new EMData();
	out->set_size(maxrin,2,1);
	float *dout = out->get_data();
	for (int i=0; i<maxrin; i++) {dout(i,0)=static_cast<float>(q[i]); dout(i,1)=static_cast<float>(t[i]);}
	//out->set_size(maxrin,1,1);
	//float *dout = out->get_data();
	//for (int i=0; i<maxrin; i++) {dout(i,0)=q[i];}
	free(t);
	free(q);
	return out;
}


vector<float> Util::Crosrng_msg_vec_p(EMData* circ1, EMData* circ2, vector<int> numr ) {

	int maxrin = numr[numr.size()-1];

        vector<float> r(2*maxrin);

        Crosrng_msg_vec( circ1, circ2, numr, &r[0], &r[maxrin] );

        return r;
}

#define  dout(i,j)        dout[i+maxrin*j]
#define  circ1b(i)        circ1b[i-1]
#define  circ2b(i)        circ2b[i-1]

void Util::Crosrng_msg_vec(EMData* circ1, EMData* circ2, vector<int> numr, float *q, float *t) {

   // dimension         circ1(lcirc),circ2(lcirc)

	int   ip, jc, numr3i, numr2i, i, j;
	float t1, t2, t3, t4, c1, c2, d1, d2;

	int nring = numr.size()/3;
	//int lcirc = numr[3*nring-2]+numr[3*nring-1]-1;
	int maxrin = numr[numr.size()-1];

	float* circ1b = circ1->get_data();
	float* circ2b = circ2->get_data();

#ifdef _WIN32
	ip = -(int)(log((float)maxrin)/log(2.0f));
#else
	ip = -(int)(log2(maxrin));
#endif	//_WIN32
	for (int i=1; i<=maxrin; i++)  {q(i) = 0.0f; t(i) = 0.0f;}

	//  q - straight  = circ1 * conjg(circ2)

	//   t - mirrored  = conjg(circ1) * conjg(circ2)

	for (i=1; i<=nring; i++) {

		numr3i = numr(3,i);
		numr2i = numr(2,i);

		t1   = circ1b(numr2i) * circ2b(numr2i);
		q(1) += t1;
		t(1) += t1;

		t1   = circ1b(numr2i+1) * circ2b(numr2i+1);
		if (numr3i == maxrin)  {
			q(2) += t1;
			t(2) += t1;
		} else {
			q(numr3i+1) += t1;
			t(numr3i+1) += t1;
		}

		for (j=3; j<=numr3i; j=j+2) {
			jc     = j+numr2i-1;

			c1     = circ1b(jc);
			c2     = circ1b(jc+1);
			d1     = circ2b(jc);
			d2     = circ2b(jc+1);

			t1     = c1 * d1;
			t3     = c1 * d2;
			t2     = c2 * d2;
			t4     = c2 * d1;

			q(j)   +=  t1 + t2;
			q(j+1) += -t3 + t4;
			t(j)   +=  t1 - t2;
			t(j+1) += -t3 - t4;
		}
	}
	// straight
	fftr_q(q,ip);
	//for (int i=0; i<maxrin; i++) cout<<i<<"  B    "<<q[i]<<"       "<<t[i]<<endl;

	// mirrored
	fftr_q(t,ip);
}

void Util::Crosrng_msg_vec_snake(EMData* circ1, EMData* circ2, vector<int> numr, float *q) {

   // dimension         circ1(lcirc),circ2(lcirc)

	int   ip, jc, numr3i, numr2i, i, j;
	float t1, t2, t3, t4, c1, c2, d1, d2;

	int nring = numr.size()/3;
	//int lcirc = numr[3*nring-2]+numr[3*nring-1]-1;
	int maxrin = numr[numr.size()-1];

	float* circ1b = circ1->get_data();
	float* circ2b = circ2->get_data();

#ifdef _WIN32
	ip = -(int)(log((float)maxrin)/log(2.0f));
#else
	ip = -(int)(log2(maxrin));
#endif	//_WIN32
	for (int i=1; i<=maxrin; i++)  {q(i) = 0.0f; }

	//  q - straight  = circ1 * conjg(circ2)

	//   t - mirrored  = conjg(circ1) * conjg(circ2)

	for (i=1; i<=nring; i++) {

		numr3i = numr(3,i);
		numr2i = numr(2,i);
		
		t1   = circ1b(numr2i) * circ2b(numr2i);
		q(1) += t1;
		
		t1   = circ1b(numr2i+1) * circ2b(numr2i+1);
		if (numr3i == maxrin)  {
			q(2) += t1;
		} else {
			q(numr3i+1) += t1;
		}

		for (j=3; j<=numr3i; j=j+2) {
			jc     = j+numr2i-1;

			c1     = circ1b(jc);
			c2     = circ1b(jc+1);
			d1     = circ2b(jc);
			d2     = circ2b(jc+1);

			t1     = c1 * d1;
			t3     = c1 * d2;
			t2     = c2 * d2;
			t4     = c2 * d1;

			q(j)   +=  t1 + t2;
			q(j+1) += -t3 + t4;
			
		}
	}
	// straight
	fftr_q(q,ip);
	//for (int i=0; i<maxrin; i++) cout<<i<<"  B    "<<q[i]<<"       "<<t[i]<<endl;


}




EMData* Util::Crosrng_msg_s(EMData* circ1, EMData* circ2, vector<int> numr)
{

	int   ip, jc, numr3i, numr2i, i, j;
	float c1, c2, d1, d2;

	int nring = numr.size()/3;
	int maxrin = numr[numr.size()-1];

	float* circ1b = circ1->get_data();
	float* circ2b = circ2->get_data();

	double *q;

	q = (double*)calloc(maxrin,sizeof(double));

#ifdef _WIN32
	ip = -(int)(log((float)maxrin)/log(2.0f));
#else
	ip = -(int)(log2(maxrin));
#endif	//_WIN32

	 //  q - straight  = circ1 * conjg(circ2)

	for (i=1;i<=nring;i++) {

		numr3i = numr(3,i);
		numr2i = numr(2,i);

		q(1) += circ1b(numr2i) * circ2b(numr2i);

		if (numr3i == maxrin)   q(2) += circ1b(numr2i+1) * circ2b(numr2i+1);
		else             q(numr3i+1) += circ1b(numr2i+1) * circ2b(numr2i+1);

		for (j=3;j<=numr3i;j=j+2) {
			jc     = j+numr2i-1;

			c1     = circ1b(jc);
			c2     = circ1b(jc+1);
			d1     = circ2b(jc);
			d2     = circ2b(jc+1);

			q(j)   +=  c1 * d1 + c2 * d2;
			q(j+1) += -c1 * d2 + c2 * d1;
		}
	}

	// straight
	fftr_d(q,ip);

	EMData* out = new EMData();
	out->set_size(maxrin,1,1);
	float *dout = out->get_data();
	for (int i=0; i<maxrin; i++) dout[i]=static_cast<float>(q[i]);
	free(q);
	return out;

}


EMData* Util::Crosrng_msg_m(EMData* circ1, EMData* circ2, vector<int> numr)
{

	int   ip, jc, numr3i, numr2i, i, j;
	float c1, c2, d1, d2;

	int nring = numr.size()/3;
	int maxrin = numr[numr.size()-1];

	float* circ1b = circ1->get_data();
	float* circ2b = circ2->get_data();

	double *t;

	t = (double*)calloc(maxrin,sizeof(double));

#ifdef _WIN32
	ip = -(int)(log((float)maxrin)/log(2.0f));
#else
	ip = -(int)(log2(maxrin));
#endif	//_WIN32

	 //   t - mirrored  = conjg(circ1) * conjg(circ2)

	for (i=1;i<=nring;i++) {

		numr3i = numr(3,i);
		numr2i = numr(2,i);
		t(1) += circ1b(numr2i) * circ2b(numr2i);

		if (numr3i == maxrin)  t(2) += circ1b(numr2i+1) * circ2b(numr2i+1);
		else          t(numr3i+1) += circ1b(numr2i+1) * circ2b(numr2i+1);

		for (j=3;j<=numr3i;j=j+2) {
			jc     = j+numr2i-1;

			c1     = circ1b(jc);
			c2     = circ1b(jc+1);
			d1     = circ2b(jc);
			d2     = circ2b(jc+1);

			t(j)   +=  c1 * d1 - c2 * d2;
			t(j+1) += -c1 * d2 - c2 * d1;
		}
	}

	// mirrored
	fftr_d(t,ip);

	EMData* out = new EMData();
	out->set_size(maxrin,1,1);
	float *dout = out->get_data();
	for (int i=0; i<maxrin; i++) dout[i]=static_cast<float>(t[i]);
	free(t);
	return out;

}

#undef circ1b
#undef circ2b
#undef dout

#undef  circ1
#undef  circ2
#undef  t
#undef  q
#undef  b
#undef  t7


#define    QUADPI      		    3.141592653589793238462643383279502884197
#define    PI2                      2*QUADPI

float Util::ener(EMData* ave, vector<int> numr) {
	ENTERFUNC;
	long double ener,en;

	int nring = numr.size()/3;
	float *aveptr = ave->get_data();

	ener = 0.0;
	for (int i=1; i<=nring; i++) {
		int numr3i = numr(3,i);
		int np     = numr(2,i)-1;
		float tq = static_cast<float>(PI2*numr(1,i)/numr3i);
		en = tq*(aveptr[np]*aveptr[np]+aveptr[np+1]*aveptr[np+1])*0.5;
		for (int j=np+2; j<np+numr3i-1; j++) en += tq*aveptr[j]*aveptr[j];
		ener += en/numr3i;
	}
	EXITFUNC;
	return static_cast<float>(ener);
}

float Util::ener_tot(const vector<EMData*>& data, vector<int> numr, vector<float> tot) {
	ENTERFUNC;
	long double ener, en;
	float arg, cs, si;

	int nima = data.size();
	int nring = numr.size()/3;
	int maxrin = numr(3,nring);

	ener = 0.0;
	for (int i=1; i<=nring; i++) {
		int numr3i = numr(3,i);
		int np     = numr(2,i)-1;
		float tq = static_cast<float>(PI2*numr(1,i)/numr3i);
		float temp1 = 0.0, temp2 = 0.0;
		for (int kk=0; kk<nima; kk++) {
			float *ptr = data[kk]->get_data();
			temp1 += ptr[np];
			temp2 += static_cast<float>(ptr[np+1]*cos(PI2*(tot[kk]-1.0f)/2.0f*numr3i/maxrin));
		}
		en = tq*(temp1*temp1+temp2*temp2)*0.5;
		for (int j=2; j<numr3i; j+=2) {
			float tempr = 0.0, tempi = 0.0;
			for (int kk=0; kk<nima; kk++) {
				float *ptr = data[kk]->get_data();
				arg = static_cast<float>( PI2*(tot[kk]-1.0)*(j/2)/maxrin );
				cs = cos(arg);
				si = sin(arg);
				tempr += ptr[np + j]*cs - ptr[np + j +1]*si;
				tempi += ptr[np + j]*si + ptr[np + j +1]*cs;
			}
			en += tq*(tempr*tempr+tempi*tempi);
		}
		ener += en/numr3i;
	}
	EXITFUNC;
	return static_cast<float>(ener);
}

void Util::update_fav(EMData* avep, EMData* datp, float tot, int mirror, vector<int> numr) {
	int nring = numr.size()/3;
	float *ave = avep->get_data();
	float *dat = datp->get_data();
	int i, j, numr3i, np;
	float  arg, cs, si;
	int maxrin = numr(3,nring);
	if(mirror == 1) { //for mirrored data has to be conjugated
		for (i=1; i<=nring; i++) {
			numr3i = numr(3,i);
			np     = numr(2,i)-1;
			ave[np]   += dat[np];
			ave[np+1] += static_cast<float>( dat[np+1]*cos(PI2*(tot-1.0f)/2.0f*numr3i/maxrin) );
			for (j=2; j<numr3i; j=j+2) {
				arg = static_cast<float>( PI2*(tot-1.)*(j/2)/maxrin );
				cs = cos(arg);
				si = sin(arg);
				//complex(data[np + j],data[np + j +1])*complex(cos(arg),sin(arg))
				ave[np + j]    += dat[np + j]*cs - dat[np + j +1]*si;
				ave[np + j +1] -= dat[np + j]*si + dat[np + j +1]*cs;
			}
		}
	} else {
		for (i=1; i<=nring; i++) {
			numr3i = numr(3,i);
			np     = numr(2,i)-1;
			ave[np]   += dat[np];
			ave[np+1] += static_cast<float>( dat[np+1]*cos(PI2*(tot-1.0f)/2.0f*numr3i/maxrin) );
			for (j=2; j<numr3i; j=j+2) {
				arg = static_cast<float>( PI2*(tot-1.)*(j/2)/maxrin );
				cs = cos(arg);
				si = sin(arg);
				//complex(data[np + j],data[np + j +1])*complex(cos(arg),sin(arg))
				ave[np + j]    += dat[np + j]*cs - dat[np + j +1]*si;
				ave[np + j +1] += dat[np + j]*si + dat[np + j +1]*cs;
			}
		}
	}
	avep->update();
	EXITFUNC;
}

void Util::sub_fav(EMData* avep, EMData* datp, float tot, int mirror, vector<int> numr) {
	int nring = numr.size()/3;
	float *ave = avep->get_data();
	float *dat = datp->get_data();
	int i, j, numr3i, np;
	float  arg, cs, si;
	int maxrin = numr(3,nring);
	if(mirror == 1) { //for mirrored data has to be conjugated
		for (i=1; i<=nring; i++) {
			numr3i = numr(3,i);
			np     = numr(2,i)-1;
			ave[np]   -= dat[np];
			ave[np+1] -= static_cast<float>( dat[np+1]*cos(PI2*(tot-1.0f)/2.0f*numr3i/maxrin) );
			for (j=2; j<numr3i; j=j+2) {
				arg = static_cast<float>( PI2*(tot-1.)*(j/2)/maxrin );
				cs = cos(arg);
				si = sin(arg);
				//complex(data[np + j],data[np + j +1])*complex(cos(arg),sin(arg))
				ave[np + j]    -= dat[np + j]*cs - dat[np + j +1]*si;
				ave[np + j +1] += dat[np + j]*si + dat[np + j +1]*cs;
			}
		}
	} else {
		for (i=1; i<=nring; i++) {
			numr3i = numr(3,i);
			np     = numr(2,i)-1;
			ave[np]   -= dat[np];
			ave[np+1] -= static_cast<float>( dat[np+1]*cos(PI2*(tot-1.0f)/2.0f*numr3i/maxrin) );
			for (j=2; j<numr3i; j=j+2) {
				arg = static_cast<float>( PI2*(tot-1.)*(j/2)/maxrin );
				cs = cos(arg);
				si = sin(arg);
				//complex(data[np + j],data[np + j +1])*complex(cos(arg),sin(arg))
				ave[np + j]    -= dat[np + j]*cs - dat[np + j +1]*si;
				ave[np + j +1] -= dat[np + j]*si + dat[np + j +1]*cs;
			}
		}
	}
	avep->update();
	EXITFUNC;
}


#undef    QUADPI
#undef    PI2

#undef  numr
#undef  circ


#define	QUADPI	 3.141592653589793238462643383279502884197
#define	PI2      QUADPI*2
#define deg_rad  QUADPI/180.0
#define rad_deg  180.0/QUADPI

struct ori_t
{
    int iphi;
    int itht;
    int id;
};


struct cmpang
{
    bool operator()( const ori_t& a, const ori_t& b )
    {
        if( a.itht != b.itht )
        {
            return a.itht < b.itht;
        }

        return a.iphi < b.iphi;
    }
};


vector<double> Util::cml_weights(const vector<float>& cml){
	static const int NBIN = 100;
	int nline=cml.size()/2;
	vector<double> weights(nline);

	vector<ori_t> angs(nline);
	for( int i=0; i < nline; ++i ) {
		angs[i].iphi = int( NBIN*cml[2*i] );
		angs[i].itht = int( NBIN*cml[2*i+1] );
		if( angs[i].itht == 180*NBIN ) angs[i].itht = 0;
		angs[i].id = i;
	}

	//std::cout << "# of angs: " << angs.size() << std::endl;

	std::sort( angs.begin(), angs.end(), cmpang() );

	vector<float> newphi;
	vector<float> newtht;
	vector< vector<int> > indices;

	int curt_iphi = -1;
	int curt_itht = -1;
	for(unsigned int i=0 ;i < angs.size(); ++i ) {
		if( angs[i].iphi==curt_iphi && angs[i].itht==curt_itht ) {
			Assert( indices.size() > 0 );
			indices.back().push_back(angs[i].id);
		} else {
			curt_iphi = angs[i].iphi;
			curt_itht = angs[i].itht;

			newphi.push_back( float(curt_iphi)/NBIN );
			newtht.push_back( float(curt_itht)/NBIN );
			indices.push_back( vector<int>(1,angs[i].id) );
		}
	}

	//std::cout << "# of indpendent ang: " << newphi.size() << std::endl;


	int num_agl = newphi.size();

	if(num_agl>2) {
		vector<double> w=Util::vrdg(newphi, newtht);

		Assert( w.size()==newphi.size() );
		Assert( indices.size()==newphi.size() );

		for(unsigned int i=0; i < newphi.size(); ++i ) {
		    /*
		    std::cout << "phi,tht,w,n: ";
		    std::cout << boost::format( "%10.3f" ) % newphi[i] << " ";
		    std::cout << boost::format( "%10.3f" ) % newtht[i] << " ";
		    std::cout << boost::format( "%8.6f"  ) % w[i] << " ";
		    std::cout << indices[i].size() << "(";
		    */

		    for(unsigned int j=0; j < indices[i].size(); ++j ) {
		    	    int id = indices[i][j];
		    	    weights[id] = w[i]/indices[i].size();
		    	    //std::cout << id << " ";
		    }

		    //std::cout << ")" << std::endl;

		}
	} else {
		cout<<"warning in Util.cml_weights"<<endl;
		double val = PI2/float(nline);
		for(int i=0; i<nline; i++)  weights[i]=val;
	}

	return weights;

}

/****************************************************
 * New code for common-lines
 ****************************************************/

void Util::set_line(EMData* img, int posline, EMData* line, int offset, int length)
{
	int i;
	int nx=img->get_xsize();
	float *img_ptr  = img->get_data();
	float *line_ptr = line->get_data();
	for (i=0;i<length;i++) img_ptr[nx*posline + i] = line_ptr[offset + i];
	img->update();
}

void Util::cml_prepare_line(EMData* sino, EMData* line, int ilf, int ihf, int pos_line, int nblines){
	int j;
	int nx = sino->get_xsize();
	int i = nx * pos_line;
	float r1, r2;
	float *line_ptr = line->get_data();
	float *sino_ptr = sino->get_data();
	for (j=ilf;j<=ihf; j += 2) {
		r1 = line_ptr[j];
		r2 = line_ptr[j + 1];
		sino_ptr[i + j - ilf] = r1;
		sino_ptr[i + j - ilf + 1] = r2;
		sino_ptr[i + nx * nblines + j - ilf] = r1;
		sino_ptr[i + nx * nblines + j - ilf + 1] = -r2;
	}
	sino->update();
}

vector<double> Util::cml_init_rot(vector<float> Ori){
	int nb_ori = Ori.size() / 4;
	int i, ind;
	float ph, th, ps;
	double cph, cth, cps, sph, sth, sps;
	vector<double> Rot(nb_ori*9);
	for (i=0; i<nb_ori; ++i){
		ind = i*4;
		// spider convention phi=psi-90, psi=phi+90
		ph = Ori[ind+2]-90;
		th = Ori[ind+1];
		ps = Ori[ind]+90;
		ph *= deg_rad;
		th *= deg_rad;
		ps *= deg_rad;
		// pre-calculate some trigo stuffs
		cph = cos(ph);
		cth = cos(th);
		cps = cos(ps);
		sph = sin(ph);
		sth = sin(th);
		sps = sin(ps);
		// fill rotation matrix
		ind = i*9;
		Rot[ind] = cph*cps-cth*sps*sph;
		Rot[ind+1] = cph*sps+cth*cps*sph;
		Rot[ind+2] = sth*sph;
		Rot[ind+3] = -sph*cps-cth*sps*cph;
		Rot[ind+4] = -sph*sps+cth*cps*cph;
		Rot[ind+5] = sth*cph;
		Rot[ind+6] = sth*sps;
		Rot[ind+7] = -sth*cps;
		Rot[ind+8] = cth;
	}

	return Rot;
}

vector<float> Util::cml_update_rot(vector<float> Rot, int iprj, float nph, float th, float nps){
	float ph, ps;
	double cph, cth, cps, sph, sth, sps;
	int ind = iprj*9;
	// spider convention phi=psi-90, psi=phi+90
	ph = nps-90;
	ps = nph+90;
	ph *= deg_rad;
	th *= deg_rad;
	ps *= deg_rad;
	// pre-calculate some trigo stuffs
	cph = cos(ph);
	cth = cos(th);
	cps = cos(ps);
	sph = sin(ph);
	sth = sin(th);
	sps = sin(ps);
	// fill rotation matrix
	Rot[ind] = (float)(cph*cps-cth*sps*sph);
	Rot[ind+1] = (float)(cph*sps+cth*cps*sph);
	Rot[ind+2] = (float)(sth*sph);
	Rot[ind+3] = (float)(-sph*cps-cth*sps*cph);
	Rot[ind+4] = (float)(-sph*sps+cth*cps*cph);
	Rot[ind+5] = (float)(sth*cph);
	Rot[ind+6] = (float)(sth*sps);
	Rot[ind+7] = (float)(-sth*cps);
	Rot[ind+8] = (float)(cth);

	return Rot;
}

vector<int> Util::cml_line_insino(vector<float> Rot, int i_prj, int n_prj){
	vector<int> com(2*(n_prj - 1));
	int a = i_prj*9;
	int i, b, c;
	int n1=0, n2=0;
	float vmax = 1 - 1.0e-6f;
	double r11, r12, r13, r23, r31, r32, r33;

	c = 0;
	for (i=0; i<n_prj; ++i){
		if (i!=i_prj){
			b = i*9;
			// this is equivalent to R = A*B'
			r11 = Rot[a]*Rot[b]+Rot[a+1]*Rot[b+1]+Rot[a+2]*Rot[b+2];
			r12 = Rot[a]*Rot[b+3]+Rot[a+1]*Rot[b+4]+Rot[a+2]*Rot[b+5];
			r13 = Rot[a]*Rot[b+6]+Rot[a+1]*Rot[b+7]+Rot[a+2]*Rot[b+8];
			r23 = Rot[a+3]*Rot[b+6]+Rot[a+4]*Rot[b+7]+Rot[a+5]*Rot[b+8];
			r31 = Rot[a+6]*Rot[b]+Rot[a+7]*Rot[b+1]+Rot[a+8]*Rot[b+2];
			r32 = Rot[a+6]*Rot[b+3]+Rot[a+7]*Rot[b+4]+Rot[a+8]*Rot[b+5];
			r33 = Rot[a+6]*Rot[b+6]+Rot[a+7]*Rot[b+7]+Rot[a+8]*Rot[b+8];
			if (r33 > vmax) {
			    n2 = 270;
			    n1 = 270 + nint180((float)(rad_deg*atan2(r12, r11)));
			}
			else if (r33 < -vmax) {
			    n2 = 270;
			    n1 = 270 - nint180((float)(rad_deg*atan2(r12, r11)));
			} else {
			    n2 = nint180((float)(rad_deg*atan2(r31, -r32)));
			    n1 = nint180((float)(rad_deg*atan2(r13, r23)));
			    if (n1 < 0) {n1 += 360;}
			    if (n2 <= 0) {n2 = abs(n2);}
			    else {n2 = 360 - n2;}
			}

			if (n1 >= 360){n1 = n1 % 360;}
			if (n2 >= 360){n2 = n2 % 360;}

			// store common-lines
			b = c*2;
			com[b] = n1;
			com[b+1] = n2;
			++c;
		}
	}

    return com;

}

vector<int> Util::cml_line_insino_all(vector<float> Rot, vector<int> seq, int, int n_lines) {
	vector<int> com(2*n_lines);
	int a=0, b, c, l;
	int n1=0, n2=0, mem=-1;
	float vmax = 1 - 1.0e-6f;
	double r11, r12, r13, r23, r31, r32, r33;
	c = 0;
	for (l=0; l<n_lines; ++l){
		c = 2*l;
		if (seq[c]!=mem){
		    mem = seq[c];
		    a = seq[c]*9;
		}
		b = seq[c+1]*9;

		// this is equivalent to R = A*B'
		r11 = Rot[a]*Rot[b]+Rot[a+1]*Rot[b+1]+Rot[a+2]*Rot[b+2];
		r12 = Rot[a]*Rot[b+3]+Rot[a+1]*Rot[b+4]+Rot[a+2]*Rot[b+5];
		r13 = Rot[a]*Rot[b+6]+Rot[a+1]*Rot[b+7]+Rot[a+2]*Rot[b+8];
		r23 = Rot[a+3]*Rot[b+6]+Rot[a+4]*Rot[b+7]+Rot[a+5]*Rot[b+8];
		r31 = Rot[a+6]*Rot[b]+Rot[a+7]*Rot[b+1]+Rot[a+8]*Rot[b+2];
		r32 = Rot[a+6]*Rot[b+3]+Rot[a+7]*Rot[b+4]+Rot[a+8]*Rot[b+5];
		r33 = Rot[a+6]*Rot[b+6]+Rot[a+7]*Rot[b+7]+Rot[a+8]*Rot[b+8];
		if (r33 > vmax) {
		    n2 = 270;
		    n1 = 270 + nint180((float)(rad_deg*atan2(r12, r11)));
		}
		else if (r33 < -vmax) {
		    n2 = 270;
		    n1 = 270 - nint180((float)(rad_deg*atan2(r12, r11)));
		} else {
		    n2 = nint180((float)(rad_deg*atan2(r31, -r32)));
		    n1 = nint180((float)(rad_deg*atan2(r13, r23)));
		    if (n1 < 0) {n1 += 360;}
		    if (n2 <= 0) {n2 = abs(n2);}
		    else {n2 = 360 - n2;}
		}
		if (n1 >= 360){n1 = n1 % 360;}
		if (n2 >= 360){n2 = n2 % 360;}

		// store common-lines
		com[c] = n1;
		com[c+1] = n2;
	}

	return com;

}

vector<double> Util::cml_line_in3d(vector<float> Ori, vector<int> seq, int, int nlines){
	// seq is the pairwise index ij: 0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3
	vector<double> cml(2*nlines); // [phi, theta] / line
	float ph1, th1;
	float ph2, th2;
	double nx, ny, nz;
	double norm;
	double sth1=0, sph1=0, cth1=0, cph1=0;
	double sth2, sph2, cth2, cph2;
	int l, ind, c;
	int mem = -1;
	for (l=0; l<nlines; ++l){
		c = 2*l;
		if (seq[c]!=mem){
			mem = seq[c];
			ind = 4*seq[c];
			ph1 = Ori[ind]*deg_rad;
			th1 = Ori[ind+1]*deg_rad;
			sth1 = sin(th1);
			sph1 = sin(ph1);
			cth1 = cos(th1);
			cph1 = cos(ph1);
		}
		ind = 4*seq[c+1];
		ph2 = Ori[ind]*deg_rad;
		th2 = Ori[ind+1]*deg_rad;
		sth2 = sin(th2);
		cth2 = cos(th2);
		sph2 = sin(ph2);
		cph2 = cos(ph2);
		// cross product
		nx = sth1*cph1*cth2 - cth1*sth2*cph2;
		ny = cth1*sth2*sph2 - cth2*sth1*sph1;
		nz = sth1*sph1*sth2*cph2 - sth1*cph1*sth2*sph2;
		norm = sqrt(nx*nx+ny*ny+nz*nz);
		nx /= norm;
		ny /= norm;
		nz /= norm;
		// apply mirror if need
		if (nz<0) {nx=-nx; ny=-ny; nz=-nz;}
		// compute theta and phi
		cml[c+1] = acos(nz);
		if (cml[c+1] == 0) {cml[c] = 0;}
		else {
			cml[c+1] *= rad_deg;
			if (cml[c+1] > 89.99) {cml[c+1] = 89.99;} // this fix some pb in Voronoi
			cml[c] = rad_deg * atan2(nx, ny);
			cml[c] = fmod(360 + cml[c], 360);

		}
	}

	return cml;
}

double Util::cml_disc(const vector<EMData*>& data, vector<int> com, vector<int> seq, vector<float> weights, int n_lines) {
	double res = 0;
	double buf = 0;
	float* line_1;
	float* line_2;
	int i, n, ind;
	int lnlen = data[0]->get_xsize();
	for (n=0; n<n_lines; ++n) {
		ind = n*2;
		line_1 = data[seq[ind]]->get_data() + com[ind] * lnlen;
		line_2 = data[seq[ind+1]]->get_data() + com[ind+1] *lnlen;
		buf = 0;
		for (i=0; i<lnlen; ++i) {
		    buf += (line_1[i]-line_2[i])*(line_1[i]-line_2[i]);
		}
		res += buf * weights[n];
	}

	return res;

}

vector<double> Util::cml_spin_psi(const vector<EMData*>& data, vector<int> com, vector<float> weights, \
				 int iprj, vector<int> iw, int n_psi, int d_psi, int n_prj){
	// res: [best_disc, best_ipsi]
	// seq: pairwise indexes ij, 0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3
	// iw : index to know where is the weight for the common-lines on the current projection in the all weights, [12, 4, 2, 7]
	vector<double> res(2);
	int lnlen = data[0]->get_xsize();
	int end = 2*(n_prj-1);
	double disc, buf, bdisc, tmp;
	int n, i, ipsi, ind, bipsi, c;
	float* line_1;
	float* line_2;
	bdisc = 1.0e6;
	bipsi = -1;
	// loop psi
	for(ipsi=0; ipsi<n_psi; ipsi += d_psi) {
		// discrepancy
		disc = 0;
		c = 0;
		for (n=0; n<n_prj; ++n) {
			if(n!=iprj) {
				ind = 2*c;
				line_1 = data[iprj]->get_data() + com[ind] * lnlen;
				line_2 = data[n]->get_data() + com[ind+1] * lnlen;
				buf = 0;
				for (i=0; i<lnlen; ++i) {
					tmp = line_1[i]-line_2[i];
					buf += tmp*tmp;
				}
				disc += buf * weights[iw[c]];
				++c;
			}
		}
		// select the best value
		if (disc <= bdisc) {
			bdisc = disc;
			bipsi = ipsi;
		}
		// update common-lines
		for (i=0; i<end; i+=2){
			com[i] += d_psi;
			if (com[i] >= n_psi) com[i] = com[i] - n_psi;
		}
	}
	res[0] = bdisc;
	res[1] = float(bipsi);

	return res;
}

vector<double> Util::cml_spin_psi_now(const vector<EMData*>& data, vector<int> com, \
				 int iprj, vector<int> iw, int n_psi, int d_psi, int n_prj){
	// res: [best_disc, best_ipsi]
	// seq: pairwise indexes ij, 0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3
	// iw : index to know where is the weight for the common-lines on the current projection in the all weights, [12, 4, 2, 7]
	vector<double> res(2);
	int lnlen = data[0]->get_xsize();
	int end = 2*(n_prj-1);
	double disc, buf, bdisc, tmp;
	int n, i, ipsi, ind, bipsi, c;
	float* line_1;
	float* line_2;
	bdisc = 1.0e6;
	bipsi = -1;
	// loop psi
	for(ipsi=0; ipsi<n_psi; ipsi += d_psi) {
		// discrepancy
		disc = 0;
		c = 0;
		for (n=0; n<n_prj; ++n) {
			if(n!=iprj) {
				ind = 2*c;
				line_1 = data[iprj]->get_data() + com[ind] * lnlen;
				line_2 = data[n]->get_data() + com[ind+1] * lnlen;
				buf = 0;
				for (i=0; i<lnlen; ++i) {
					tmp = line_1[i]-line_2[i];
					buf += tmp*tmp;
				}
				disc += buf;
				++c;
			}
		}
		// select the best value
		if (disc <= bdisc) {
			bdisc = disc;
			bipsi = ipsi;
		}
		// update common-lines
		for (i=0; i<end; i+=2){
			com[i] += d_psi;
			if (com[i] >= n_psi) com[i] = com[i] - n_psi;
		}
	}
	res[0] = bdisc;
	res[1] = float(bipsi);

	return res;
}

#undef	QUADPI
#undef	PI2
#undef  deg_rad
#undef  rad_deg

/****************************************************
 * END OF NEW CODE FOR COMMON-LINES
 ****************************************************/

// helper function for k-means
Dict Util::min_dist_real(EMData* image, const vector<EMData*>& data) {
	ENTERFUNC;

	int nima = data.size();
	vector<float> res(nima);
	double result = 0.;
	double valmin = 1.0e20;
	int valpos = -1;

	for (int kk=0; kk<nima; kk++){
	result = 0;

	float *y_data = data[kk]->get_data();
	float *x_data = image->get_data();
	long totsize = image->get_xsize()*image->get_ysize();
	for (long i = 0; i < totsize; i++) {
	    double temp = x_data[i]- y_data[i];
	    result += temp*temp;
	}
	result /= totsize;
	res[kk] = (float)result;

	if(result<valmin) {valmin = result; valpos = kk;}

	}

	Dict retvals;
	retvals["dist"] = res;
	retvals["pos"]  = valpos;

	EXITFUNC;
	return retvals;

}

Dict Util::min_dist_four(EMData* image, const vector<EMData*>& data) {
	ENTERFUNC;

	int nima = data.size();
	vector<float> res(nima);
	double result = 0.;
	double valmin = 1.0e20;
	int valpos = -1;

	for (int kk=0; kk<nima; kk++){
	result = 0;
	//validate_input_args(image, data[kk]);

	float *y_data = data[kk]->get_data();
	float *x_data = image->get_data();

	// Implemented by PAP  01/09/06 - please do not change.  If in doubts, write/call me.
	int nx  = data[kk]->get_xsize();
	int ny  = data[kk]->get_ysize();
	nx = (nx - 2 + data[kk]->is_fftodd()); // nx is the real-space size of the input image
	int lsd2 = (nx + 2 - nx%2) ; // Extended x-dimension of the complex image

	int ixb = 2*((nx+1)%2);
	int iyb = ny%2;
	int iz = 0;

	for ( int iy = 0; iy <= ny-1; iy++) {
	    for ( int ix = 2; ix <= lsd2 - 1 - ixb; ix++) {
		int ii = ix + (iy  + iz * ny)* lsd2;
		result += (x_data[ii] - y_data[ii])*double(x_data[ii] - y_data[ii]);
	    }
	}
	for ( int iy = 1; iy <= ny/2-1 + iyb; iy++) {
	    int ii = (iy  + iz * ny)* lsd2;
	    result += (x_data[ii] - y_data[ii])*double(x_data[ii] - y_data[ii]);
	    result += (x_data[ii+1] - y_data[ii+1])*double(x_data[ii+1] - y_data[ii+1]);
	}
	if(nx%2 == 0) {
	    for ( int iy = 1; iy <= ny/2-1 + iyb; iy++) {
		int ii = lsd2 - 2 + (iy  + iz * ny)* lsd2;
		result += (x_data[ii] - y_data[ii])*double(x_data[ii] - y_data[ii]);
		result += (x_data[ii+1] - y_data[ii+1])*double(x_data[ii+1] - y_data[ii+1]);
	    }

	}
	result *= 2;
	result += (x_data[0] - y_data[0])*double(x_data[0] - y_data[0]);
	if(ny%2 == 0) {
	    int ii = (ny/2  + iz * ny)* lsd2;
	    result += (x_data[ii] - y_data[ii])*double(x_data[ii] - y_data[ii]);
	}
	if(nx%2 == 0) {
	    int ii = lsd2 - 2 + (0  + iz * ny)* lsd2;
	    result += (x_data[ii] - y_data[ii])*double(x_data[ii] - y_data[ii]);
	    if(ny%2 == 0) {
		int ii = lsd2 - 2 +(ny/2  + iz * ny)* lsd2;
		result += (x_data[ii] - y_data[ii])*double(x_data[ii] - y_data[ii]);
	    }
	}

	result /= (long int)nx*(long int)ny*(long int)nx*(long int)ny;
	res[kk] = (float)result;

	if(result<valmin) {valmin = result; valpos = kk;}

	}

	Dict retvals;
	retvals["dist"] = res;
	retvals["pos"]  = valpos;

	EXITFUNC;
	return retvals;
}

int Util::k_means_cont_table_(int* group1, int* group2, int* stb, long int s1, long int s2, int flag) {
    long int d2 = group2[s2 - 1] - group2[0];
    long int p2 = 0;
    long int i1 = 0;
    long int i2 = 0;
    long int max = 0;
    long int cont = 0;
    long int i = 0;
    int stop1 = 0;
    int stop2 = 0;

    for (i=0; i<s1; i++) {
	p2 = (long int)(s2 * (double)group1[i] / (double)d2);
	if (p2 >= s2) {p2 = s2 - 1;}
	i1 = p2;
	i2 = p2;
	max = s2;
	if (group1[i] < group2[0] || group1[i] > group2[s2 - 1]) {continue;}

	stop1 = 0;
	stop2 = 0;
	while (max--) {
	    if (group1[i] == group2[i1]) {
		if (flag) {stb[cont] = group1[i];}
		cont++;
		break;
	    }
	    if (group2[i1] < group1[i]) {stop1=1;}
	    if (group1[i] == group2[i2]) {
		if (flag) {stb[cont] = group1[i];}
		cont++;
		break;
	    }
	    if (group2[i2] > group1[i]) {stop2=1;}
	    //printf("i1 %li i2 %li    v2 %i v2 %i   stop1 %i stop2 %i\n", i1, i2, group2[i1], group2[i2], stop1, stop2);

	    if (stop1 & stop2) {break;}
	    i1--;
	    i2++;
	    if (i1 < 0) {i1 = 0;}
	    if (i2 >= s2) {i2 = s2 - 1;}
	}
	//printf("v1: %i    ite: %li   cont: %li\n", group1[i], s2-max, cont);
    }

    return cont;
}



#define old_ptr(i,j,k)          old_ptr[i+(j+(k*ny))*(size_t)nx]
#define new_ptr(iptr,jptr,kptr) new_ptr[iptr+(jptr+(kptr*new_ny))*(size_t)new_nx]
EMData* Util::decimate(EMData* img, int x_step, int y_step, int z_step)
{
	/* Exception Handle */
	if (!img) {
		throw NullPointerException("NULL input image");
	}
	/* ============================== */

	// Get the size of the input image
	int nx=img->get_xsize(),ny=img->get_ysize(),nz=img->get_zsize();
	/* ============================== */


	/* Exception Handle */
	if ((x_step-1 > nx/2 || y_step-1 > ny/2 || z_step-1 > nz/2) || (x_step-1)<0 || (y_step-1)<0 || (z_step-1)<0)
	{
		LOGERR("Parameters for decimation cannot exceed the center of the image.");
		throw ImageDimensionException("Parameters for decimation cannot exceed the center of the image.");
	}
	/* ============================== */


	/*    Calculation of the start point */
	int new_st_x=(nx/2)%x_step, new_st_y=(ny/2)%y_step, new_st_z=(nz/2)%z_step;
	/* ============================*/


	/* Calculation of the size of the decimated image */
	int rx=2*(nx/(2*x_step)), ry=2*(ny/(2*y_step)), rz=2*(nz/(2*z_step));
	int r1=int(ceil((nx-(x_step*rx))/(1.f*x_step))), r2=int(ceil((ny-(y_step*ry))/(1.f*y_step)));
	int r3=int(ceil((nz-(z_step*rz))/(1.f*z_step)));
	if(r1>1){r1=1;}
	if(r2>1){r2=1;}
	if(r3>1){r3=1;}
	int new_nx=rx+r1, new_ny=ry+r2, new_nz=rz+r3;
	/* ===========================================*/


	EMData* img2 = new EMData();
	img2->set_size(new_nx,new_ny,new_nz);
	float *new_ptr = img2->get_data();
	float *old_ptr = img->get_data();
	int iptr, jptr, kptr = 0;
	for (int k=new_st_z; k<nz; k+=z_step) {jptr=0;
		for (int j=new_st_y; j<ny; j+=y_step) {iptr=0;
			for (int i=new_st_x; i<nx; i+=x_step) {
				new_ptr(iptr,jptr,kptr) = old_ptr(i,j,k);
			iptr++;}
		jptr++;}
	kptr++;}
	img2->update();
	return img2;
}
#undef old_ptr
#undef new_ptr

#define inp(i,j,k)  inp[(i+new_st_x)+((j+new_st_y)+((k+new_st_z)*ny))*(size_t)nx]
#define outp(i,j,k) outp[i+(j+(k*new_ny))*(size_t)new_nx]
EMData* Util::window(EMData* img,int new_nx,int new_ny, int new_nz, int x_offset, int y_offset, int z_offset)
{
	/* Exception Handle */
	if (!img) throw NullPointerException("NULL input image");
	/* ============================== */

	// Get the size of the input image
	int nx=img->get_xsize(), ny=img->get_ysize(), nz=img->get_zsize();
	/* ============================== */

	/* Exception Handle */
	if(new_nx>nx || new_ny>ny || new_nz>nz)
		throw ImageDimensionException("The size of the windowed image cannot exceed the input image size. 1");
	if((nx/2)-(new_nx/2)+x_offset<0 || (ny/2)-(new_ny/2)+y_offset<0 || (nz/2)-(new_nz/2)+z_offset<0)
		throw ImageDimensionException("The offset inconsistent with the input image size. 2");
	if(x_offset>((nx-(nx/2))-(new_nx-(new_nx/2))) || y_offset>((ny-(ny/2))-(new_ny-(new_ny/2))) || z_offset>((nz-(nz/2))-(new_nz-(new_nz/2))))
		throw ImageDimensionException("The offset inconsistent with the input image size. 3");
	/* ============================== */

	/*    Calculation of the start point */
	int  new_st_x = nx/2-new_nx/2 + x_offset,
	     new_st_y = ny/2-new_ny/2 + y_offset,
	     new_st_z = nz/2-new_nz/2 + z_offset;
	/* ============================== */

	/* Exception Handle */
	if (new_st_x<0 || new_st_y<0 || new_st_z<0)   //  WHAT HAPPENS WITH THE END POINT CHECK??  PAP
		throw ImageDimensionException("The offset inconsistent with the input image size. 4");
	/* ============================== */

	EMData* wind = img->copy_empty_head();
	wind->set_size(new_nx, new_ny, new_nz);
	float *outp=wind->get_data();
	float *inp=img->get_data();

	for (int k=0; k<new_nz; k++)
		for(int j=0; j<new_ny; j++)
			for(int i=0; i<new_nx; i++)
				outp(i,j,k) = inp(i,j,k);
	wind->update();
	return wind;
}
#undef inp
#undef outp

#define inp(i,j,k) inp[i+(j+(k*ny))*(size_t)nx]
#define outp(i,j,k) outp[(i+new_st_x)+((j+new_st_y)+((k+new_st_z)*new_ny))*(size_t)new_nx]
EMData *Util::pad(EMData* img,int new_nx, int new_ny, int new_nz, int x_offset, int y_offset, int z_offset, const char *params)
{
	/* Exception Handle */
	if (!img)  throw NullPointerException("NULL input image");
	/* ============================== */

	// Get the size of the input image
	int nx=img->get_xsize(),ny=img->get_ysize(),nz=img->get_zsize();
	/* ============================== */

	/* Exception Handle */
	if(new_nx<nx || new_ny<ny || new_nz<nz)
		throw ImageDimensionException("The size of the padded image cannot be lower than the input image size.");
	if((new_nx/2)-(nx/2)+x_offset<0 || (new_ny/2)-(ny/2)+y_offset<0 || (new_nz/2)-(nz/2)+z_offset<0)
		throw ImageDimensionException("The offset inconsistent with the input image size. Solution: Change the offset parameters");
	if(x_offset>((new_nx-(new_nx/2))-(nx-(nx/2))) || y_offset>((new_ny-(new_ny/2))-(ny-(ny/2))) || z_offset>((new_nz-(new_nz/2))-(nz-(nz/2))))
		throw ImageDimensionException("The offset inconsistent with the input image size. Solution: Change the offset parameters");
	/* ============================== */

	EMData* pading = img->copy_head();
	pading->set_size(new_nx, new_ny, new_nz);
	float *inp  = img->get_data();
	float *outp = pading->get_data();


	/* Calculation of the average and the circumference values for background substitution
	=======================================================================================*/
	float background;

	if (strcmp(params,"average")==0) background = img->get_attr("mean");
	else if (strcmp(params,"circumference")==0) {
		float sum1=0.0f;
		size_t cnt=0;
		for(int i=0;i<nx;i++) {
			sum1 += inp(i,0,0) + inp(i,ny-1,nz-1);
			cnt+=2;
		}
		if(nz-1 == 0) {
			for (int j=1;j<ny-1;j++) {
				sum1 += inp(1,j,0) + inp(nx-1,j,0);
				cnt+=2;
			}
		} else {
			for (int k=1;k<nz-1;k++) {
				for (int j=1;j<ny-1;j++) {
					sum1 += inp(1,j,0) + inp(nx-1,j,0);
					cnt+=2;
				}
			}
		}
		background = sum1/cnt;
	} else {
		background = static_cast<float>( atof( params ) );
	}
	/*=====================================================================================*/

	 /*Initial Padding */
	int new_st_x=0,new_st_y=0,new_st_z=0;
	for (int k=0;k<new_nz;k++)
		for(int j=0;j<new_ny;j++)
			for (int i=0;i<new_nx;i++)
				outp(i,j,k)=background;
	/*============================== */

	/*    Calculation of the start point */
	new_st_x=int((new_nx/2-nx/2)  + x_offset);
	new_st_y=int((new_ny/2-ny/2)  + y_offset);
	new_st_z=int((new_nz/2-nz/2)  + z_offset);
	/* ============================== */

	for (int k=0;k<nz;k++)
		for(int j=0;j<ny;j++)
			for(int i=0;i<nx;i++)
				outp(i,j,k)=inp(i,j,k);
	pading->update();
	return pading;
}
#undef inp
#undef outp
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

void Util::colreverse(float* beg, float* end, int nx) {
	float* tmp = new float[nx];
	int n = (end - beg)/nx;
	int nhalf = n/2;
	for (int i = 0; i < nhalf; i++) {
		// swap col i and col n-1-i
		memcpy(tmp, beg+i*nx, nx*sizeof(float));
		memcpy(beg+i*nx, beg+(n-1-i)*nx, nx*sizeof(float));
		memcpy(beg+(n-1-i)*nx, tmp, nx*sizeof(float));
	}
	delete[] tmp;
}

void Util::slicereverse(float *beg, float *end, int nx,int ny)
{
        int nxy = nx*ny;
	colreverse(beg, end, nxy);
}


void Util::cyclicshift(EMData *image, Dict params) {

	if (image->is_complex()) throw ImageFormatException("Real image required for IntegerCyclicShift2DProcessor");

	int dx = params["dx"];
	int dy = params["dy"];
	int dz = params["dz"];

	// The reverse trick we're using shifts to the left (a negative shift)
	int nx = image->get_xsize();
	dx %= nx;
	if (dx < 0) dx += nx;
	int ny = image->get_ysize();
	dy %= ny;
	if (dy < 0) dy += ny;
	int nz = image->get_zsize();
	dz %= nz;
	if (dz < 0) dz += nz;

	int mx = -(dx - nx);
	int my = -(dy - ny);
	int mz = -(dz - nz);

	float* data = image->get_data();
	// x-reverses
	if (mx != 0) {
		for (int iz = 0; iz < nz; iz++)
	               for (int iy = 0; iy < ny; iy++) {
				// reverses for column iy
	        		size_t offset = nx*iy + (size_t)nx*ny*iz; // starting location for column iy in slice iz
				reverse(&data[offset],&data[offset+mx]);
				reverse(&data[offset+mx],&data[offset+nx]);
				reverse(&data[offset],&data[offset+nx]);
	        	}
	}
	// y-reverses
	if (my != 0) {
		for (int iz = 0; iz < nz; iz++) {
	        	size_t offset = (size_t)nx*ny*iz;
			colreverse(&data[offset], &data[offset + my*nx], nx);
			colreverse(&data[offset + my*nx], &data[offset + ny*nx], nx);
			colreverse(&data[offset], &data[offset + ny*nx], nx);
		}
	}
	if (mz != 0) {
		slicereverse(&data[0], &data[(size_t)mz*ny*nx], nx, ny);
		slicereverse(&data[mz*ny*nx], &data[(size_t)nz*ny*nx], nx, ny);
		slicereverse(&data[0], &data[(size_t)nz*ny*nx], nx ,ny);
	}
	image->update();
}

//-----------------------------------------------------------------------------------------------------------------------


vector<float> Util::histogram(EMData* image, EMData* mask, int nbins, float hmin, float hmax)
{
	if (image->is_complex())
                throw ImageFormatException("Cannot do histogram on Fourier image");
	//float hmax, hmin;
	float *imageptr=0, *maskptr=0;
	int nx=image->get_xsize();
	int ny=image->get_ysize();
	int nz=image->get_zsize();

	if(mask != NULL){
		if(nx != mask->get_xsize() || ny != mask->get_ysize() || nz != mask->get_zsize())
			throw ImageDimensionException("The size of mask image should be of same size as the input image");
		maskptr =mask->get_data();
	}
	if( nbins == 0) nbins = nx;
	vector <float> freq(2*nbins, 0.0);

	imageptr=image->get_data();
	if( hmin == hmax ) {
		if(mask == NULL) {
			hmax = image->get_attr("maximum");
			hmin = image->get_attr("minimum");
		} else {
			bool  First = true;
			for (size_t i = 0;i < (size_t)nx*ny*nz; i++) {
	      		if (maskptr[i]>=0.5f) {
					if(First) {
						hmax = imageptr[i];
						hmin = imageptr[i];
						First = false;
					} else {
						hmax = (hmax < imageptr[i])?imageptr[i]:hmax;
						hmin = (hmin > imageptr[i])?imageptr[i]:hmin;
					}
				}
			}
		}
	}
	float hdiff = hmax - hmin;
	float ff = (nbins-1)/hdiff;
	for (int i = 0; i < nbins; i++) freq[nbins+i] = hmin + (float(i)+0.5f)/ff;
	if(mask == NULL) {
		for(size_t i = 0; i < (size_t)nx*ny*nz; i++) {
			int jbin = static_cast<int>((imageptr[i]-hmin)*ff + 1.5);
			if(jbin >= 1 && jbin <= nbins)  freq[jbin-1] += 1.0;
		}
	} else {
		for(size_t i = 0; i < (size_t)nx*ny*nz; i++) {
			if(maskptr[i] >= 0.5) {
				int jbin = static_cast<int>((imageptr[i]-hmin)*ff + 1.5);
				if(jbin >= 1 && jbin <= nbins)  freq[jbin-1] += 1.0;
			}
		}
	}
	return freq;
}

Dict Util::histc(EMData *ref,EMData *img, EMData *mask)
{
	/* Exception Handle */
	if (img->is_complex() || ref->is_complex())
                throw ImageFormatException("Cannot do Histogram on Fourier Image");

	if(mask != NULL){
		if(img->get_xsize() != mask->get_xsize() || img->get_ysize() != mask->get_ysize() || img->get_zsize() != mask->get_zsize())
			throw ImageDimensionException("The size of mask image should be of same size as the input image"); }
	/* ===================================================== */

	/* Image size calculation */
	size_t size_ref = ((size_t)(ref->get_xsize())*(ref->get_ysize())*(ref->get_zsize()));
	size_t size_img = ((size_t)(img->get_xsize())*(img->get_ysize())*(img->get_zsize()));
	/* ===================================================== */

	/* The reference image attributes */
	float *ref_ptr = ref->get_data();
	float ref_h_min = ref->get_attr("minimum");
	float ref_h_max = ref->get_attr("maximum");
	float ref_h_avg = ref->get_attr("mean");
	float ref_h_sig = ref->get_attr("sigma");
	/* ===================================================== */

	/* Input image under mask attributes */
	float *mask_ptr = (mask == NULL)?img->get_data():mask->get_data();

	vector<float> img_data = Util::infomask(img, mask);
	float img_avg = img_data[0];
	float img_sig = img_data[1];

	/* The image under mask -- size calculation */
	int cnt=0;
	for(size_t i=0;i<size_img;++i)
		if (mask_ptr[i]>0.5f)
				cnt++;
	/* ===================================================== */

	/* Histogram of reference image calculation */
	float ref_h_diff = ref_h_max - ref_h_min;

	#ifdef _WIN32
		int hist_len = _cpp_min((unsigned long)size_ref/16,_cpp_min((unsigned long)size_img/16,256lu));
	#else
		int hist_len = std::min((unsigned long)size_ref/16,std::min((unsigned long)size_img/16,256lu));
	#endif	//_WIN32

	float *ref_freq_bin = new float[3*hist_len];

	//initialize value in each bin to zero
	for (int i = 0;i < (3*hist_len);i++) ref_freq_bin[i] = 0.f;

	for (size_t i = 0;i < size_ref;++i) {
		int L = static_cast<int>(((ref_ptr[i] - ref_h_min)/ref_h_diff) * (hist_len-1) + hist_len+1);
		ref_freq_bin[L]++;
	}
	for (int i = 0;i < (3*hist_len);i++) ref_freq_bin[i] *= static_cast<float>(cnt)/static_cast<float>(size_ref);

	//Parameters Calculation (i.e) 'A' x + 'B'
	float A = ref_h_sig/img_sig;
	float B = ref_h_avg - (A*img_avg);

	vector<float> args;
	args.push_back(A);
	args.push_back(B);

	vector<float> scale;
	scale.push_back(1.e-7f*A);
	scale.push_back(-1.e-7f*B);

	vector<float> ref_freq_hist;
	for(int i = 0;i < (3*hist_len);i++) ref_freq_hist.push_back((int)ref_freq_bin[i]);

	vector<float> data;
	data.push_back(ref_h_diff);
	data.push_back(ref_h_min);

	Dict parameter;

	/* Parameters displaying the arguments A & B, and the scaling function and the data's */
	parameter["args"] = args;
	parameter["scale"]= scale;
	parameter["data"] = data;
	parameter["ref_freq_bin"] = ref_freq_hist;
	parameter["size_img"]=(double)size_img;
	parameter["hist_len"]=hist_len;
	/* ===================================================== */

	return parameter;
}


float Util::hist_comp_freq(float PA,float PB,size_t size_img, int hist_len, EMData *img, vector<float> ref_freq_hist, EMData *mask, float ref_h_diff, float ref_h_min)
{
	float *img_ptr = img->get_data();
	float *mask_ptr = (mask == NULL)?img->get_data():mask->get_data();

	int *img_freq_bin = new int[3*hist_len];
	for(int i = 0;i < (3*hist_len);i++) img_freq_bin[i] = 0;
	for(size_t i = 0;i < size_img;++i) {
		if(mask_ptr[i] > 0.5f) {
			float img_xn = img_ptr[i]*PA + PB;
			int L = static_cast<int>(((img_xn - ref_h_min)/ref_h_diff) * (hist_len-1) + hist_len+1);
			if(L >= 0 && L < (3*hist_len)) img_freq_bin[L]++;
		}
	}
	int freq_hist = 0;

	for(int i = 0;i < (3*hist_len);i++) freq_hist += (int)pow((float)((int)ref_freq_hist[i] - (int)img_freq_bin[i]),2.f);
	freq_hist = (-freq_hist);
	return static_cast<float>(freq_hist);
}
//------------------------------------------------------------------------------------------------------------------------------------------------------------------
#define    QUADPI      		        3.141592653589793238462643383279502884197
#define    DGR_TO_RAD    		QUADPI/180
#define    DM(I)         		DM	    [I-1]
#define    SS(I)         		SS	    [I-1]
Dict Util::CANG(float PHI,float THETA,float PSI)
{
	double CPHI,SPHI,CTHE,STHE,CPSI,SPSI;
	vector<float>	DM,SS;

	for(int i =0;i<9;i++) DM.push_back(0);

	for(int i =0;i<6;i++) SS.push_back(0);

	CPHI = cos(double(PHI)*DGR_TO_RAD);
	SPHI = sin(double(PHI)*DGR_TO_RAD);
	CTHE = cos(double(THETA)*DGR_TO_RAD);
	STHE = sin(double(THETA)*DGR_TO_RAD);
	CPSI = cos(double(PSI)*DGR_TO_RAD);
	SPSI = sin(double(PSI)*DGR_TO_RAD);

	SS(1) = float(CPHI);
	SS(2) = float(SPHI);
	SS(3) = float(CTHE);
	SS(4) = float(STHE);
	SS(5) = float(CPSI);
	SS(6) = float(SPSI);

	DM(1) = float(CPHI*CTHE*CPSI-SPHI*SPSI);
	DM(2) = float(SPHI*CTHE*CPSI+CPHI*SPSI);
	DM(3) = float(-STHE*CPSI);
	DM(4) = float(-CPHI*CTHE*SPSI-SPHI*CPSI);
	DM(5) = float(-SPHI*CTHE*SPSI+CPHI*CPSI);
	DM(6) = float(STHE*SPSI);
	DM(7) = float(STHE*CPHI);
	DM(8) = float(STHE*SPHI);
	DM(9) = float(CTHE);

	Dict DMnSS;
	DMnSS["DM"] = DM;
	DMnSS["SS"] = SS;

	return(DMnSS);
}
#undef SS
#undef DM
#undef QUADPI
#undef DGR_TO_RAD
//-----------------------------------------------------------------------------------------------------------------------
struct t_BPCQ_line{
	int rX;     // radius along X axe
	int offset; // offset of voxel in volume - beginning of the line
	float xbb;  // XBB coefficient
	float ybb;  // YBB coefficient
};

void Util::BPCQ( EMData *B, EMData *CUBE, const int radius )
{
	if (B->is_complex()) {
		B->do_ift_inplace();
		B->depad();
	}

	const Transform * transform = B->get_attr("xform.projection");
	Dict transform_params = transform->get_params("spider");

	// ---- build DM matrix (transform matrix) - convert from 3x4 matrix to 2x3 matrix (only 2 first rows are nedeed)
	std::vector<float> DM = transform->get_matrix();
	DM[3+0] = DM[4+0];
	DM[3+1] = DM[4+1];
	DM[3+2] = DM[4+2];

	delete transform;

	const int NSAM = B->get_xsize();
	const int NROW = B->get_ysize();

	// buffer "lines_to_process" should be aligned to size of cache line (usually 64 or 128 bytes)
	t_BPCQ_line * lines_to_process;
#ifdef _WIN32
	if ( (lines_to_process = (t_BPCQ_line *)_aligned_malloc( 4*radius*radius*sizeof(t_BPCQ_line), 256 )) == NULL )
#else
	if ( posix_memalign( reinterpret_cast<void**>(&lines_to_process), 256, 4*radius*radius*sizeof(t_BPCQ_line) ) != 0 )
#endif	//_WIN32
	{
		throw std::bad_alloc();
	}
	t_BPCQ_line * first_free_line = lines_to_process;

	// calculate lines parameters
	{
		//  Unsure about sign of shifts, check later PAP 06/28/09
		const float x_shift_plus_center = float(NSAM/2 +1) + float(transform_params[ "tx" ]);
		const float y_shift_plus_center = float(NROW/2 +1) + float(transform_params[ "ty" ]);

		const int sizeX = CUBE->get_xsize();
		const int sizeY = CUBE->get_ysize();
		const int sizeZ = CUBE->get_zsize();

		const int centerX = sizeX / 2;
		const int centerY = sizeY / 2;
		const int centerZ = sizeZ / 2;

		const int minRZ = (centerZ >= radius) ? (-radius) : (-centerZ);
		const int maxRZ = (sizeZ > centerZ+radius) ? (radius) : (sizeZ-centerZ-1);

		for ( int rZ=minRZ; rZ<=maxRZ; ++rZ ) {
			for ( int rY=-radius; rY<=radius; ++rY ) {
				const int sqRX = radius*radius - rZ*rZ - rY*rY;
				if (sqRX >= 0) {
#ifdef	_WIN32
					first_free_line->rX     = static_cast<int>( floor(sqrtf(sqRX)+0.5) );
#else
					first_free_line->rX     = static_cast<int>( roundf(sqrtf(sqRX)) );
#endif	//_WIN32
					first_free_line->offset = sizeX*( centerY+rY + sizeY*(centerZ+rZ) ) + centerX - first_free_line->rX;
					first_free_line->xbb    = rZ*DM[2] + rY*DM[1] + x_shift_plus_center;
					first_free_line->ybb    = rZ*DM[5] + rY*DM[4] + y_shift_plus_center;
					++first_free_line;
				}
			}
		}
	}

	const float * const Bptr = B->get_data();
	float * const CUBE_begin = CUBE->get_data();

	// update voxels in volume
	// this loop takes more than 95% of calculations time spent in Util::BPCQ function
	for ( t_BPCQ_line * iLine = lines_to_process; iLine < first_free_line; ++iLine ) {
		const int rX_first = -(iLine->rX);
		const int rX_last  =   iLine->rX;
		float  *CUBE_ptr = CUBE_begin + iLine->offset;
		for (int rX=rX_first; rX<=rX_last; ++rX, ++CUBE_ptr) {
			const float XB  = rX * DM[0] + iLine->xbb;
			const float YB  = rX * DM[3] + iLine->ybb;
			const int IQX = int(XB);
			const int IQY = int(YB);
			if ( IQX < 1 || IQX >= NSAM || IQY < 1 || IQY >= NROW )  continue;
			const float DIPX = XB-IQX;
			const float DIPY = YB-IQY;
			const float b00 = Bptr[IQX-1+((IQY-1)*NSAM)];
			const float b01 = Bptr[IQX-1+((IQY-0)*NSAM)];
			const float b10 = Bptr[IQX-0+((IQY-1)*NSAM)];
			const float b11 = Bptr[IQX-0+((IQY-0)*NSAM)];
			*(CUBE_ptr) = *(CUBE_ptr) + b00 + DIPY*(b01-b00) + DIPX*(b10-b00+DIPY*(b11-b10-b01+b00));
		}
	}

	free(lines_to_process);
}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#define    W(i,j) 			Wptr    [i+(j)*Wnx]
#define    PROJ(i,j) 		PROJptr [i+(j)*NNNN]
#define    SS(i,j)         	SS	    [i+(j)*6]

void Util::WTF(EMData* PROJ,vector<float> SS,float SNR,int K)
{

	--K; // now indexes are started from 0

	int NSAM = PROJ->get_xsize();
	int NROW = PROJ->get_ysize();

	if (PROJ->is_fftpadded()) {
		NSAM -= (PROJ->is_fftodd()) ? (1) : (2);  // correction for DFT image
	}

    const int ntotal = NSAM*NROW;
	const float q = 2.0f;
	const float qt = 8.0f/q;
	//  Fix for padding 2x
	const int ipad = 1;
	NSAM *= ipad;
	NROW *= ipad;
	const int NNNN = NSAM+2-(NSAM%2);
	const int NX2 = NSAM/2;
	const int NR2  = NROW/2;

	const int NANG = int(SS.size())/6;

	EMData* W = new EMData();
	const int Wnx = NNNN/2;
	W->set_size(Wnx,NROW,1);
	W->to_zero();
	float *Wptr = W->get_data();

	for (int L=0; L<NANG; L++) {
		const float tmp1 = SS(2,K)*SS(3,L)*(SS(0,K)*SS(0,L) + SS(1,K)*SS(1,L)) - SS(2,L)*SS(3,K);
		const float tmp2 = SS(3,L)*( SS(0,K)*SS(1,L) - SS(0,L)*SS(1,K) );
		float OX = SS(5,K)*tmp2 + SS(4,K)*tmp1;
		float OY = SS(4,K)*tmp2 - SS(5,K)*tmp1;
		if(OX < 0.0f) {
			OX = -OX;
			OY = -OY;
		}

		if( OX > 1.0e-6f || fabs(OY) > 1.0e-6f ) {
			for (int J=0; J<NROW; ++J) {
				const float JY_OY = (J > NR2) ? ((J-NROW)*OY) : (J*OY);
				int xma = NX2;
				int xmi = 0;
				const float fxma = ( q-JY_OY) / OX;
				const float fxmi = (-q-JY_OY) / OX;
				if (fxma < xmi || fxmi > xma ) {
					continue;
				}
				if (fxma < xma) {
					xma = static_cast<int>(fxma+0.5f);
				}
				if (fxmi > xmi) {
					xmi = static_cast<int>(fxmi+0.5f);
				}
				for( int I=xmi; I<=xma; ++I ) {
					const float Y = I*OX + JY_OY;
					W(I,J) += exp(-qt*Y*Y);
				}
			}
		} else {
			for (int J=0; J<NROW; ++J) {
				for (int I=0; I<NNNN/2; ++I) {
					W(I,J) += 1.0f;
				}
			}
		}
	}

    EMData* proj_in = PROJ;

	const bool realOnInput = PROJ->is_real();
    if (realOnInput) {
		// copy input image and run DFT on it
		PROJ = PROJ->norm_pad( false, ipad);
		PROJ->do_fft_inplace();
		PROJ->update();
    }
    float * PROJptr = PROJ->get_data();

	const float osnr = 1.0f/SNR;
	const float WNRMinv = 1.0f/W(0,0);
	for (int J=0; J<NROW; ++J)  {
		float sy = (J > NR2) ? (J - NROW) : (J);
		sy /= NROW;
		sy *= sy;
		for (int I=0; I<NNNN; I+=2) {
			const int KX = I/2;
			const float temp = W(KX,J)*WNRMinv;
			float WW = temp/(temp*temp + osnr);
			// This is supposed to fix fall-off due to Gaussian function in the weighting function
			const float sx = float(KX) / NSAM;
			WW *= exp(qt*(sy + sx*sx));
			PROJ(I,J)   *= WW;
			PROJ(I+1,J) *= WW;
		}
	}
	delete W; W = 0;

	PROJ->do_ift_inplace();
	PROJ->depad();

	if (realOnInput) {
		// copy data back to input image
		float* data_src = PROJ->get_data();
		float* data_dst = proj_in->get_data();
		memcpy( data_dst, data_src, ntotal * sizeof(float) );
		delete PROJ;
	}

	proj_in->update();
}
/*
void Util::WTF(EMData* PROJ,vector<float> SS,float SNR,int K)
{
	int NSAM,NROW,NNNN,NR2,L,JY,KX,NANG;
	float WW,OX,OY,Y;

	NSAM = PROJ->get_xsize();
	NROW = PROJ->get_ysize();
	//  Fix for padding 2x
        int ntotal = NSAM*NROW;
	int ipad = 1;
	NSAM *= ipad;
	NROW *= ipad;
	NNNN = NSAM+2-(NSAM%2);
	NR2  = NROW/2;

	NANG = int(SS.size())/6;

	EMData* W = new EMData();
	int Wnx = NNNN/2;
	W->set_size(Wnx,NROW,1);
	W->to_zero();
	float *Wptr = W->get_data();
	float *PROJptr = PROJ->get_data();
	for (L=1; L<=NANG; L++) {
		float  tmp1 = SS(3,K)*SS(4,L)*(SS(1,K)*SS(1,L) + SS(2,K)*SS(2,L)) - SS(3,L)*SS(4,K);
		float  tmp2 = SS(4,L)*( SS(1,K)*SS(2,L) - SS(1,L)*SS(2,K) ); 
		OX = SS(6,K)*tmp2 + SS(5,K)*tmp1;
		OY = SS(5,K)*tmp2 - SS(6,K)*tmp1;
	//cout << " OX   "<<OX << " OY   "<<OY <<endl;

		if( fabs(OX) > 1.0e-6f || fabs(OY) > 1.0e6f) {
			for(int J=1;J<=NROW;J++) {
				JY = (J-1);
				if(JY > NR2) JY=JY-NROW;
				for(int I=1;I<=NNNN/2;I++) {
					Y =  fabs(OX * (I-1) + OY * JY);
					if(Y < 2.0f) {
					W(I,J) += exp(-4*Y*Y);
	cout << " L   "<<L << " I   "<<I-1 << " JY   "<<JY << " ARG   "<<4*Y*Y<<endl;}
				}
			}
		} else {
			for(int J=1;J<=NROW;J++) for(int I=1;I<=NNNN/2;I++)  W(I,J) += 1.0f;
		}
	}
        EMData* proj_in = PROJ;

	PROJ = PROJ->norm_pad( false, ipad);
	PROJ->do_fft_inplace();
	PROJ->update();
	//cout << " x   "<<PROJ->get_xsize() << " y   "<<PROJ->get_ysize() <<endl;
	PROJptr = PROJ->get_data();

	float WNRMinv,temp;
	float osnr = 1.0f/SNR;
	WNRMinv = 1.0f/W(1,1);
	for(int J=1;J<=NROW;J++)
		for(int I=1;I<=NNNN;I+=2) {
			KX	     = (I+1)/2;
			temp	     = W(KX,J)*WNRMinv;
			WW           = temp/(temp*temp + osnr);
			PROJ(I,J)   *= WW;
			PROJ(I+1,J) *= WW;
		}
	delete W; W = 0;
	PROJ->do_ift_inplace();
	PROJ->depad();

        float* data_src = PROJ->get_data();
        float* data_dst = proj_in->get_data();

        for( int i=0; i < ntotal; ++i )  data_dst[i] = data_src[i];

        proj_in->update();

	delete PROJ;
}
*/
#undef PROJ
#undef W
#undef SS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#define    W(i,j) 			Wptr        [i-1+((j-1)*Wnx)]
#define    PROJ(i,j) 			PROJptr     [i-1+((j-1)*NNNN)]
#define    SS(I,J)         		SS	    [I-1 + (J-1)*6]
#define    RI(i,j)                      RI          [(i-1) + ((j-1)*3)]
#define    CC(i)                        CC          [i-1]
#define    CP(i)                        CP          [i-1]
#define    VP(i)                        VP          [i-1]
#define    VV(i)                        VV          [i-1]
#define    AMAX1(i,j)                   i>j?i:j
#define    AMIN1(i,j)                   i<j?i:j
void Util::WTM(EMData *PROJ,vector<float>SS, int DIAMETER,int NUMP)
{
	float rad2deg =(180.0f/3.1415926f);
	float deg2rad = (3.1415926f/180.0f);

	int NSAM,NROW,NNNN,NR2,NANG,L,JY;

	NSAM = PROJ->get_xsize();
	NROW = PROJ->get_ysize();

	if (PROJ->is_fftpadded()) {
		NSAM -= (PROJ->is_fftodd()) ? (1) : (2);  // correction for DFT image
	}

	NNNN = NSAM+2-(NSAM%2);
	NR2  = NROW/2;
	NANG = int(SS.size())/6;

	float RI[9];
	RI(1,1)=SS(1,NUMP)*SS(3,NUMP)*SS(5,NUMP)-SS(2,NUMP)*SS(6,NUMP);
	RI(2,1)=-SS(1,NUMP)*SS(3,NUMP)*SS(6,NUMP)-SS(2,NUMP)*SS(5,NUMP);
	RI(3,1)=SS(1,NUMP)*SS(4,NUMP);
	RI(1,2)=SS(2,NUMP)*SS(3,NUMP)*SS(5,NUMP)+SS(1,NUMP)*SS(6,NUMP);
	RI(2,2)=-SS(2,NUMP)*SS(3,NUMP)*SS(6,NUMP)+SS(1,NUMP)*SS(5,NUMP);
	RI(3,2)=SS(2,NUMP)*SS(4,NUMP);
	RI(1,3)=-SS(4,NUMP)*SS(5,NUMP);
	RI(2,3)=SS(4,NUMP)*SS(6,NUMP);
	RI(3,3)=SS(3,NUMP);

	float THICK=static_cast<float>( NSAM)/DIAMETER/2.0f ;

	EMData* W = new EMData();
	int Wnx = NNNN/2;
	W->set_size(NNNN/2,NROW,1);
	W->to_one();
	float *Wptr = W->get_data();

	float ALPHA,TMP,FV,RT,FM,CCN,CC[3],CP[2],VP[2],VV[3];

	for (L=1; L<=NANG; L++) {
		if (L != NUMP) {
			CC(1)=SS(2,L)*SS(4,L)*SS(3,NUMP)-SS(3,L)*SS(2,NUMP)*SS(4,NUMP);
			CC(2)=SS(3,L)*SS(1,NUMP)*SS(4,NUMP)-SS(1,L)*SS(4,L)*SS(3,NUMP);
			CC(3)=SS(1,L)*SS(4,L)*SS(2,NUMP)*SS(4,NUMP)-SS(2,L)*SS(4,L)*SS(1,NUMP)*SS(4,NUMP);

			TMP = sqrt(CC(1)*CC(1) +  CC(2)*CC(2) + CC(3)*CC(3));
			CCN=static_cast<float>( AMAX1( AMIN1(TMP,1.0) ,-1.0) );
			ALPHA=rad2deg*float(asin(CCN));
			if (ALPHA>180.0f) ALPHA=ALPHA-180.0f;
			if (ALPHA>90.0f)  ALPHA=180.0f-ALPHA;
			if(ALPHA<1.0E-6) {
				for(int J=1;J<=NROW;J++) for(int I=1;I<=NNNN/2;I++) W(I,J)+=1.0;
			} else {
				FM=THICK/(fabs(sin(ALPHA*deg2rad)));
				CC(1)	= CC(1)/CCN;CC(2)   = CC(2)/CCN;CC(3)	= CC(3)/CCN;
				VV(1)= SS(2,L)*SS(4,L)*CC(3)-SS(3,L)*CC(2);
				VV(2)= SS(3,L)*CC(1)-SS(1,L)*SS(4,L)*CC(3);
				VV(3)= SS(1,L)*SS(4,L)*CC(2)-SS(2,L)*SS(4,L)*CC(1);
				CP(1)	= 0.0;CP(2) = 0.0;
				VP(1)	= 0.0;VP(2) = 0.0;

				CP(1) = CP(1) + RI(1,1)*CC(1) + RI(1,2)*CC(2) + RI(1,3)*CC(3);
				CP(2) = CP(2) + RI(2,1)*CC(1) + RI(2,2)*CC(2) + RI(2,3)*CC(3);
				VP(1) = VP(1) + RI(1,1)*VV(1) + RI(1,2)*VV(2) + RI(1,3)*VV(3);
				VP(2) = VP(2) + RI(2,1)*VV(1) + RI(2,2)*VV(2) + RI(2,3)*VV(3);

				TMP = CP(1)*VP(2)-CP(2)*VP(1);

				//     PREVENT TMP TO BE TOO SMALL, SIGN IS IRRELEVANT
				TMP = AMAX1(1.0E-4f,fabs(TMP));
				float tmpinv = 1.0f/TMP;
				for(int J=1;J<=NROW;J++) {
					JY = (J-1);
					if (JY>NR2)  JY=JY-NROW;
					for(int I=1;I<=NNNN/2;I++) {
						FV     = fabs((JY*CP(1)-(I-1)*CP(2))*tmpinv);
						RT     = 1.0f-FV/FM;
						W(I,J) += ((RT>0.0f)*RT);
					}
				}
			}
		}
	}

    EMData* proj_in = PROJ;
    const bool realOnInput = PROJ->is_real();

	if (realOnInput) {
		// copy input image and run DFT on it
		PROJ = PROJ->norm_pad( false, 1 );
		PROJ->do_fft_inplace();
		PROJ->update();
	}
	float *PROJptr = PROJ->get_data();

	int KX;
	float WW;
	for(int J=1; J<=NROW; J++)
		for(int I=1; I<=NNNN; I+=2) {
			KX	    =  (I+1)/2;
			WW	    =  1.0f/W(KX,J);
			PROJ(I,J)   = PROJ(I,J)*WW;
			PROJ(I+1,J) = PROJ(I+1,J)*WW;
		}
	delete W; W = 0;
	PROJ->do_ift_inplace();
	PROJ->depad();

	if (realOnInput) {
		// copy data back to input image
		float* data_src = PROJ->get_data();
		float* data_dst = proj_in->get_data();
		int ntotal = NSAM*NROW;
		memcpy( data_dst, data_src, ntotal * sizeof(float) );
		delete PROJ;
	}

	proj_in->update();
}
#undef   AMAX1
#undef   AMIN1
#undef   RI
#undef   CC
#undef   CP
#undef   VV
#undef   VP
#undef   W
#undef   SS
#undef   PROJ
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
float Util::tf(float dzz, float ak, float voltage, float cs, float wgh, float b_factor, float sign)
{
	float cst  = cs*1.0e7f;

	wgh /= 100.0;
	float phase = atan(wgh/sqrt(1.0f-wgh*wgh));
	float lambda=12.398f/sqrt(voltage*(1022.0f+voltage));
	float ak2 = ak*ak;
	float g1 = dzz*1.0e4f*lambda*ak2;
	float g2 = cst*lambda*lambda*lambda*ak2*ak2/2.0f;

	float ctfv = static_cast<float>( sin(M_PI*(g1-g2)+phase)*sign );
	if(b_factor != 0.0f)  ctfv *= exp(-b_factor*ak2/4.0f);
	return ctfv;
}

EMData* Util::compress_image_mask(EMData* image, EMData* mask)
{
	/***********
	***get the size of the image for validation purpose
	**************/
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	/********
	***Exception Handle
	*************/
	if(nx != mask->get_xsize() || ny != mask->get_ysize() || nz != mask->get_zsize())
		throw ImageDimensionException("The dimension of the image does not match the dimension of the mask!");

	size_t i, size = (size_t)nx*ny*nz;

	float* img_ptr = image->get_data();
	float* mask_ptr = mask->get_data();

	int ln=0;  //length of the output image = number of points under the mask.
	for(i = 0;i < size;i++) if(mask_ptr[i] > 0.5f) ln++;

	EMData* new_image = new EMData();
	new_image->set_size(ln,1,1); /* set size of the new image */
	float *new_ptr    = new_image->get_data();

	ln=-1;
	for(i = 0;i < size;i++){
		if(mask_ptr[i] > 0.5f) {
			ln++;
			new_ptr[ln]=img_ptr[i];
		}
	}

	return new_image;
}

EMData *Util::reconstitute_image_mask(EMData* image, EMData *mask )
{
	/********
	***Exception Handle
	*************/
	if(mask == NULL)
		throw ImageDimensionException("The mask cannot be an null image");

	/***********
	***get the size of the mask
	**************/
	int nx = mask->get_xsize(),ny = mask->get_ysize(),nz = mask->get_zsize();

	size_t i,size = (size_t)nx*ny*nz;			 /* loop counters */
	/* new image declaration */
	EMData *new_image = new EMData();
	new_image->set_size(nx,ny,nz); 		 /* set the size of new image */
	float *new_ptr  = new_image->get_data(); /* set size of the new image */
	float *mask_ptr = mask->get_data();	 /* assign a pointer to the mask image */
	float *img_ptr  = image->get_data();	 /* assign a pointer to the 1D image */
	int count = 0;
	float sum_under_mask = 0.0 ;
	for(i = 0;i < size;i++){
			if(mask_ptr[i] > 0.5f){
				new_ptr[i] = img_ptr[count];
				sum_under_mask += img_ptr[count];
				count++;
                                if( count > image->get_xsize() ) {
                                    throw ImageDimensionException("Error: in reconstitute_image_mask, the mask doesn't match the image, it is too large");
                                }
			}
	}

        if( count > image->get_xsize() ) {
            throw ImageDimensionException("Error: in reconstitute_image_mask, the mask doesn't match the image, it is too small");
        }

	float avg_under_mask = sum_under_mask / count;
	for(i = 0;i < size;i++) {
		if(mask_ptr[i] <= 0.5f)  new_ptr[i] = avg_under_mask;
	}
	new_image->update();
	return new_image;
}



vector<float> Util::merge_peaks(vector<float> peak1, vector<float> peak2,float p_size)
{
	vector<float>new_peak;
	int n1=peak1.size()/3;
	float p_size2=p_size*p_size;
	for (int i=0;i<n1;++i) {
		vector<float>::iterator it2= peak1.begin()+3*i;
		bool push_back1=true;
		int n2=peak2.size()/3;
		/*cout<<"peak2 size==="<<n2<<"i====="<<i<<endl;
		       cout<<"new peak size==="<<new_peak.size()/3<<endl;*/
		if(n2 ==0) {
			new_peak.push_back(*it2);
			new_peak.push_back(*(it2+1));
			new_peak.push_back(*(it2+2));
		} else  {
			int j=0;
			while (j< n2-1 ) {
				vector<float>::iterator it3= peak2.begin()+3*j;
				float d2=((*(it2+1))-(*(it3+1)))*((*(it2+1))-(*(it3+1)))+((*(it2+2))-(*(it3+2)))*((*(it2+2))-(*(it3+2)));
				if(d2< p_size2 ) {
					if( (*it2)<(*it3) ) {
						new_peak.push_back(*it3);
						new_peak.push_back(*(it3+1));
						new_peak.push_back(*(it3+2));
						peak2.erase(it3);
						peak2.erase(it3);
						peak2.erase(it3);
						push_back1=false;
					} else {
						peak2.erase(it3);
						peak2.erase(it3);
						peak2.erase(it3);
					}
				} else  j=j+1;
				n2=peak2.size()/3;
			}
			if(push_back1) {
				new_peak.push_back(*it2);
				new_peak.push_back(*(it2+1));
				new_peak.push_back(*(it2+2));
			}
		}
	}
	return new_peak;
}

int Util::coveig(int n, float *covmat, float *eigval, float *eigvec)
{
	// n size of the covariance/correlation matrix
	// covmat --- covariance/correlation matrix (n by n)
	// eigval --- returns eigenvalues
	// eigvec --- returns eigenvectors

	ENTERFUNC;

	int i;

	// make a copy of covmat so that it will not be overwritten
	for ( i = 0 ; i < n * n ; i++ )   eigvec[i] = covmat[i];

	char NEEDV = 'V';
	char UPLO = 'U';
	int lwork = -1;
	int info = 0;
	float *work, wsize;

	//  query to get optimal workspace
	ssyev_(&NEEDV, &UPLO, &n, eigvec, &n, eigval, &wsize, &lwork, &info);
	lwork = (int)wsize;

	work = (float *)calloc(lwork, sizeof(float));
	//  calculate eigs
	ssyev_(&NEEDV, &UPLO, &n, eigvec, &n, eigval, work, &lwork, &info);
	free(work);
	EXITFUNC;
	return info;
}

Dict Util::coveig_for_py(int ncov, const vector<float>& covmatpy)
{

	ENTERFUNC;
	int len = covmatpy.size();
	float *eigvec;
	float *eigval;
	float *covmat;
	int status = 0;
        eigval = (float*)calloc(ncov,sizeof(float));
        eigvec = (float*)calloc(ncov*ncov,sizeof(float));
	covmat = (float*)calloc(ncov*ncov, sizeof(float));

	const float *covmat_ptr;
	covmat_ptr = &covmatpy[0];
	for(int i=0;i<len;i++){
	    covmat[i] = covmat_ptr[i];
	}

	status = Util::coveig(ncov, covmat, eigval, eigvec);

	vector<float> eigval_py(ncov);
	const float *eigval_ptr;
	eigval_ptr = &eigval[0];
	for(int i=0;i<ncov;i++){
	    eigval_py[i] = eigval_ptr[i];
	}

	vector<float> eigvec_py(ncov*ncov);
	const float *eigvec_ptr;
	eigvec_ptr = &eigvec[0];
	for(int i=0;i<ncov*ncov;i++){
	    eigvec_py[i] = eigvec_ptr[i];
	}

        Dict res;
        res["eigval"] = eigval_py;
        res["eigvec"] = eigvec_py;

	EXITFUNC;
	return res;
}

vector<float> Util::pw_extract(vector<float>pw, int n, int iswi, float ps)
{
	int k,m,n1,klmd,klm2d,nklmd,n2d,n_larg,l, n2;

	k=(int)pw.size();
	l=0;
	m=k;
	n2=n+2;
	n1=n+1;
	klmd=k+l+m;
	klm2d= k+l+m+2;
	nklmd=k+l+m+n;
	n2d=n+2;
	/*size has to be increased when N is large*/
	n_larg=klmd*2;
	klm2d=n_larg+klm2d;
	klmd=n_larg+klmd;
	nklmd=n_larg+nklmd;
	int size_q=klm2d*n2d;
	int size_cu=nklmd*2;
	static int i__;

	double *q ;
	double *x ;
	double *res;
	double *cu;
	float  *q2;
	float  *pw_;
	long int *iu;
	double *s;
	q   = (double*)calloc(size_q,sizeof(double));
	x   = (double*)calloc(n2d,sizeof(double));
	res = (double*)calloc(klmd,sizeof(double));
	cu  = (double*)calloc(size_cu,sizeof(double));
	s   = (double*)calloc(klmd,sizeof(double));
	q2  = (float*)calloc(size_q,sizeof(float));
	iu  = (long int*)calloc(size_cu,sizeof(long int));
	pw_ = (float*)calloc(k,sizeof(float));

	for( i__ =0; i__<k; ++i__) pw_[i__] = log(pw[i__]);
	long int l_k=k;
	long int l_n=n;
	long int l_iswi=iswi;
	vector<float> cl1_res;
	cl1_res = Util::call_cl1(&l_k, &l_n, &ps, &l_iswi, pw_, q2, q, x, res, cu, s, iu);
	free(q);
	free(x);
	free(res);
	free(s);
	free(cu);
	free(q2);
	free(iu);
	free(pw_);
	return cl1_res;
}
vector<float> Util::call_cl1(long int *k, long int *n, float *ps, long int *iswi, float *pw, float *q2,double *q, double *x, double *res, double *cu, double *s, long int *iu)
{
    long int q2_dim1, q2_offset, q_dim1, q_offset, i__1, i__2;
    float r__1;
    int tmp__i;
    long int i__, j;
    --s;
    --res;
    iu -= 3;
    cu -= 3;
    --x;
    long int klm2d;
    klm2d= *k+*k+2;
    klm2d=klm2d+klm2d;
    q_dim1 = klm2d;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    q2_dim1 = klm2d;
    q2_offset = 1 + q2_dim1;
    q2 -= q2_offset;
    i__2=0;
    i__1 = *n - 1;
    tmp__i=0;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *k;
	tmp__i+=1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    r__1 = float(i__ - 1) /(float) *k / (*ps * 2);
	    q2[i__ + j * q2_dim1] = pow(r__1, tmp__i);
	    }
    }
    for  (i__ = 1; i__ <= i__2; ++i__)
      { q2[i__ + *n * q2_dim1] = 1.f;
	    q2[i__ + (*n + 1) * q2_dim1] = pw[i__-1];
	}
   vector<float> fit_res;
   fit_res=Util::lsfit(k, n, &klm2d, iswi, &q2[q2_offset], &q[q_offset], &x[1], &res[1], &cu[3], &s[1], &iu[3]);
   return fit_res;
}
vector<float> Util::lsfit(long int *ks,long int *n, long int *klm2d, long int *iswi,float *q1,double *q, double *x, double *res, double *cu, double *s, long int *iu)
{
    /* System generated locals */
    long int q_dim1, q_offset, q1_dim1, q1_offset, i__1, i__2;

    /* Local variables */
    long int i__, j, m, n1, ii, jj;
    double tmp;
    vector<float> p;
    --x;
    q_dim1 = *klm2d;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    q1_dim1 = *klm2d;
    q1_offset = 1 + q1_dim1;
    q1 -= q1_offset;
    --s;
    --res;
    iu -= 3;
    cu -= 3;

    /* Function Body */
    long int l = 0;

/* C==ZHONG HUANG,JULY,12,02;L=0,1,2,3,4,5,6 correspond to different equality constraints */
    m = *ks;
    n1 = *n + 1;
    if (*iswi == 1) {
	i__1 = n1;
	for (jj = 1; jj <= i__1; ++jj) {
	    i__2 = *ks;
	    for (ii = 1; ii <= i__2; ++ii) {
	/*	q[ii + jj * q_dim1] = (double) q1[ii + jj * q1_dim1];*/

		q[*ks + ii + jj * q_dim1] = (double) q1[ii + jj * q1_dim1]
			;
	    }
	}
    } else if (*iswi == 2) {
	i__1 = *ks;
	for (ii = 1; ii <= i__1; ++ii) {
	    i__2 = n1;
	    for (jj = 1; jj <= i__2; ++jj) {
		q[ii + jj * q_dim1] = (double) q1[ii + jj * q1_dim1];
		q[*ks + ii + jj * q_dim1] = -((double) q1[ii + jj * q1_dim1]);
	    }
	}
    } else if (*iswi == 3) {
	l = 2;
	i__1 = n1;
	for (jj = 1; jj <= i__1; ++jj) {
	    i__2 = *ks + 2;
	    for (ii = 1; ii <= i__2; ++ii) {
		q[ii + jj * q_dim1] = (double) q1[ii + jj * q1_dim1];
	    }
	    i__2 = *ks;
	    for (ii = 1; ii <= i__2; ++ii) {
		q[*ks + 2 + ii + jj * q_dim1] = (double) q1[ii + jj * q1_dim1];
	    }
	}
    } else if (*iswi == 4) {
	l = 2;
	i__1 = n1;
	for (jj = 1; jj <= i__1; ++jj) {
	    i__2 = *ks + 2;
	    for (ii = 1; ii <= i__2; ++ii) {
		q[ii + jj * q_dim1] = (double) q1[ii + jj * q1_dim1];
	    }
	    i__2 = *ks;
	    for (ii = 1; ii <= i__2; ++ii) {
		q[*ks + 2 + ii + jj * q_dim1] = -((double) q1[ii + jj * q1_dim1]);
	    }
	}
    } else if (*iswi == 5) {
	l = 1;
	i__1 = n1;
	for (jj = 1; jj <= i__1; ++jj) {
	    i__2 = *ks + 1;
	    for (ii = 1; ii <= i__2; ++ii) {
		q[ii + jj * q_dim1] = (double) q1[ii + jj * q1_dim1];
	    }
	    i__2 = *ks;
	    for (ii = 1; ii <= i__2; ++ii) {
		q[*ks + 1 + ii + jj * q_dim1] = -((double) q1[ii + jj * q1_dim1]);
	    }
	}
    } else if (*iswi == 6) {
	l = 1;
	i__1 = n1;
	for (jj = 1; jj <= i__1; ++jj) {
	    i__2 = *ks + 1;
	    for (ii = 1; ii <= i__2; ++ii) {
		q[ii + jj * q_dim1] = (double) q1[ii + jj * q1_dim1];
	    }
	    i__2 = *ks;
	    for (ii = 1; ii <= i__2; ++ii) {
		q[*ks + 1 + ii + jj * q_dim1] = (double) q1[ii + jj * q1_dim1];
	    }
	}
    } else if (*iswi == 7) {
	l = 3;
	i__1 = n1;
	for (jj = 1; jj <= i__1; ++jj) {
	    i__2 = *ks + 3;
	    for (ii = 1; ii <= i__2; ++ii) {
		q[ii + jj * q_dim1] = (double) q1[ii + jj * q1_dim1];
	    }
	    i__2 = *ks;
	    for (ii = 1; ii <= i__2; ++ii) {
		q[*ks + 3 + ii + jj * q_dim1] = -((double) q1[ii + jj * q1_dim1]);
	    }
	}
    } else if (*iswi == 8) {
	l = 4;
	i__1 = n1;
	for (jj = 1; jj <= i__1; ++jj) {
	    i__2 = *ks + 4;
	    for (ii = 1; ii <= i__2; ++ii) {
		q[ii + jj * q_dim1] = (double) q1[ii + jj * q1_dim1];
	    }
	    i__2 = *ks;
	    for (ii = 1; ii <= i__2; ++ii) {
		q[*ks + 4 + ii + jj * q_dim1] = -((double) q1[ii + jj * q1_dim1]);
	    }
	}
    }

    Util::cl1(ks, &l, &m, n, klm2d, &q[q_offset], &x[1], &res[1], &cu[3], &iu[3], &s[1]);
    i__1 = *ks;
    int tmp__j=0;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tmp = 0.f;
	i__2 = *n - 1;
	for (j = 1; j <= i__2; ++j) {
	    tmp__j=j;
	    tmp += pow(q1[i__ + q1_dim1], tmp__j) * x[j];
	    }
	tmp += x[*n];
	p.push_back(static_cast<float>(exp(tmp)));
	p.push_back(q1[i__ + q1_dim1]);
    }
    i__2=*n;
    for (i__=1;i__<=i__2;++i__)
    	{ p.push_back(static_cast<float>(x[i__]));}
    return p;
}
void Util::cl1(long int *k, long int *l, long int *m, long int *n, long int *klm2d,
	double *q, double *x, double *res, double *cu, long int *iu, double *s)
{

    long int q_dim1, q_offset, i__1, i__2;
    double d__1;

    static long int i__, j;
    static double z__;
    static long int n1, n2, ia, ii, kk, in, nk, js;
    static double sn, zu, zv;
    static long int nk1, klm, nkl, jmn, jpn;
    static double cuv;
    static long int klm1, nkl1, klm2, kode, iimn, nklm, iter;
    static float xmin;
    static double xmax;
    static long int iout;
    static double xsum;
    static long int iineg, maxit;
    static double toler;
    static float error;
    static double pivot;
    static long int kforce, iphase;
    static double tpivot;

    --s;
    --res;
    iu -= 3;
    cu -= 3;
    --x;
    q_dim1 = *klm2d;
    q_offset = 1 + q_dim1;
    q -= q_offset;

    /* Function Body */
    maxit = 500;
    kode = 0;
    toler = 1e-4f;
    iter = 0;
    n1 = *n + 1;
    n2 = *n + 2;
    nk = *n + *k;
    nk1 = nk + 1;
    nkl = nk + *l;
    nkl1 = nkl + 1;
    klm = *k + *l + *m;
    klm1 = klm + 1;
    klm2 = klm + 2;
    nklm = *n + klm;
    kforce = 1;
    iter = 0;
    js = 1;
    ia = 0;
/* SET UP LABELS IN Q. */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	q[klm2 + j * q_dim1] = (double) j;
/* L10: */
    }
    i__1 = klm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__ + n2 * q_dim1] = (double) (*n + i__);
	if (q[i__ + n1 * q_dim1] >= 0.f) {
	    goto L30;
	}
	i__2 = n2;
	for (j = 1; j <= i__2; ++j) {
	    q[i__ + j * q_dim1] = -q[i__ + j * q_dim1];
/* L20: */
	}
L30:
	;
    }
/* SET UP PHASE 1 COSTS. */
    iphase = 2;
    i__1 = nklm;
    for (j = 1; j <= i__1; ++j) {
	cu[(j << 1) + 1] = 0.f;
	cu[(j << 1) + 2] = 0.f;
	iu[(j << 1) + 1] = 0;
	iu[(j << 1) + 2] = 0;
/* L40: */
    }
    if (*l == 0) {
	goto L60;
    }
    i__1 = nkl;
    for (j = nk1; j <= i__1; ++j) {
	cu[(j << 1) + 1] = 1.f;
	cu[(j << 1) + 2] = 1.f;
	iu[(j << 1) + 1] = 1;
	iu[(j << 1) + 2] = 1;
/* L50: */
    }
    iphase = 1;
L60:
    if (*m == 0) {
	goto L80;
    }
    i__1 = nklm;
    for (j = nkl1; j <= i__1; ++j) {
	cu[(j << 1) + 2] = 1.f;
	iu[(j << 1) + 2] = 1;
	jmn = j - *n;
	if (q[jmn + n2 * q_dim1] < 0.f) {
	    iphase = 1;
	}
/* L70: */
    }
L80:
    if (kode == 0) {
	goto L150;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if ((d__1 = x[j]) < 0.) {
	    goto L90;
	} else if (d__1 == 0) {
	    goto L110;
	} else {
	    goto L100;
	}
L90:
	cu[(j << 1) + 1] = 1.f;
	iu[(j << 1) + 1] = 1;
	goto L110;
L100:
	cu[(j << 1) + 2] = 1.f;
	iu[(j << 1) + 2] = 1;
L110:
	;
    }
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	jpn = j + *n;
	if ((d__1 = res[j]) < 0.) {
	    goto L120;
	} else if (d__1 == 0) {
	    goto L140;
	} else {
	    goto L130;
	}
L120:
	cu[(jpn << 1) + 1] = 1.f;
	iu[(jpn << 1) + 1] = 1;
	if (q[j + n2 * q_dim1] > 0.f) {
	    iphase = 1;
	}
	goto L140;
L130:
	cu[(jpn << 1) + 2] = 1.f;
	iu[(jpn << 1) + 2] = 1;
	if (q[j + n2 * q_dim1] < 0.f) {
	    iphase = 1;
	}
L140:
	;
    }
L150:
    if (iphase == 2) {
	goto L500;
    }
/* COMPUTE THE MARGINAL COSTS. */
L160:
    i__1 = n1;
    for (j = js; j <= i__1; ++j) {
	xsum = 0.;
	i__2 = klm;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ii = (long int) q[i__ + n2 * q_dim1];
	    if (ii < 0) {
		goto L170;
	    }
	    z__ = cu[(ii << 1) + 1];
	    goto L180;
L170:
	    iineg = -ii;
	    z__ = cu[(iineg << 1) + 2];
L180:
	    xsum += q[i__ + j * q_dim1] * z__;
/*  180       XSUM = XSUM + Q(I,J)*Z */
/* L190: */
	}
	q[klm1 + j * q_dim1] = xsum;
/* L200: */
    }
    i__1 = *n;
    for (j = js; j <= i__1; ++j) {
	ii = (long int) q[klm2 + j * q_dim1];
	if (ii < 0) {
	    goto L210;
	}
	z__ = cu[(ii << 1) + 1];
	goto L220;
L210:
	iineg = -ii;
	z__ = cu[(iineg << 1) + 2];
L220:
	q[klm1 + j * q_dim1] -= z__;
/* L230: */
    }
/* DETERMINE THE VECTOR TO ENTER THE BASIS. */
L240:
    xmax = 0.f;
    if (js > *n) {
	goto L490;
    }
    i__1 = *n;
    for (j = js; j <= i__1; ++j) {
	zu = q[klm1 + j * q_dim1];
	ii = (long int) q[klm2 + j * q_dim1];
	if (ii > 0) {
	    goto L250;
	}
	ii = -ii;
	zv = zu;
	zu = -zu - cu[(ii << 1) + 1] - cu[(ii << 1) + 2];
	goto L260;
L250:
	zv = -zu - cu[(ii << 1) + 1] - cu[(ii << 1) + 2];
L260:
	if (kforce == 1 && ii > *n) {
	    goto L280;
	}
	if (iu[(ii << 1) + 1] == 1) {
	    goto L270;
	}
	if (zu <= xmax) {
	    goto L270;
	}
	xmax = zu;
	in = j;
L270:
	if (iu[(ii << 1) + 2] == 1) {
	    goto L280;
	}
	if (zv <= xmax) {
	    goto L280;
	}
	xmax = zv;
	in = j;
L280:
	;
    }
    if (xmax <= toler) {
	goto L490;
    }
    if (q[klm1 + in * q_dim1] == xmax) {
	goto L300;
    }
    i__1 = klm2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__ + in * q_dim1] = -q[i__ + in * q_dim1];
/* L290: */
    }
    q[klm1 + in * q_dim1] = xmax;
/* DETERMINE THE VECTOR TO LEAVE THE BASIS. */
L300:
    if (iphase == 1 || ia == 0) {
	goto L330;
    }
    xmax = 0.f;
    i__1 = ia;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__ = (d__1 = q[i__ + in * q_dim1], abs(d__1));
	if (z__ <= xmax) {
	    goto L310;
	}
	xmax = z__;
	iout = i__;
L310:
	;
    }
    if (xmax <= toler) {
	goto L330;
    }
    i__1 = n2;
    for (j = 1; j <= i__1; ++j) {
	z__ = q[ia + j * q_dim1];
	q[ia + j * q_dim1] = q[iout + j * q_dim1];
	q[iout + j * q_dim1] = z__;
/* L320: */
    }
    iout = ia;
    --ia;
    pivot = q[iout + in * q_dim1];
    goto L420;
L330:
    kk = 0;
    i__1 = klm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__ = q[i__ + in * q_dim1];
	if (z__ <= toler) {
	    goto L340;
	}
	++kk;
	res[kk] = q[i__ + n1 * q_dim1] / z__;
	s[kk] = (double) i__;
L340:
	;
    }
L350:
    if (kk > 0) {
	goto L360;
    }
    kode = 2;
    goto L590;
L360:
    xmin = static_cast<float>( res[1] );
    iout = (long int) s[1];
    j = 1;
    if (kk == 1) {
	goto L380;
    }
    i__1 = kk;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (res[i__] >= xmin) {
	    goto L370;
	}
	j = i__;
	xmin = static_cast<float>( res[i__] );
	iout = (long int) s[i__];
L370:
	;
    }
    res[j] = res[kk];
    s[j] = s[kk];
L380:
    --kk;
    pivot = q[iout + in * q_dim1];
    ii = (long int) q[iout + n2 * q_dim1];
    if (iphase == 1) {
	goto L400;
    }
    if (ii < 0) {
	goto L390;
    }
    if (iu[(ii << 1) + 2] == 1) {
	goto L420;
    }
    goto L400;
L390:
    iineg = -ii;
    if (iu[(iineg << 1) + 1] == 1) {
	goto L420;
    }
/* 400 II = IABS(II) */
L400:
    ii = abs(ii);
    cuv = cu[(ii << 1) + 1] + cu[(ii << 1) + 2];
    if (q[klm1 + in * q_dim1] - pivot * cuv <= toler) {
	goto L420;
    }
/* BYPASS INTERMEDIATE VERTICES. */
    i__1 = n1;
    for (j = js; j <= i__1; ++j) {
	z__ = q[iout + j * q_dim1];
	q[klm1 + j * q_dim1] -= z__ * cuv;
	q[iout + j * q_dim1] = -z__;
/* L410: */
    }
    q[iout + n2 * q_dim1] = -q[iout + n2 * q_dim1];
    goto L350;
/* GAUSS-JORDAN ELIMINATION. */
L420:
    if (iter < maxit) {
	goto L430;
    }
    kode = 3;
    goto L590;
L430:
    ++iter;
    i__1 = n1;
    for (j = js; j <= i__1; ++j) {
	if (j != in) {
	    q[iout + j * q_dim1] /= pivot;
	}
/* L440: */
    }
/* IF PERMITTED, USE SUBROUTINE COL OF THE DESCRIPTION */
/* SECTION AND REPLACE THE FOLLOWING SEVEN STATEMENTS DOWN */
/* TO AND INCLUDING STATEMENT NUMBER 460 BY.. */
/*     DO 460 J=JS,N1 */
/*        IF(J .EQ. IN) GO TO 460 */
/*        Z = -Q(IOUT,J) */
/*        CALL COL(Q(1,J), Q(1,IN), Z, IOUT, KLM1) */
/* 460 CONTINUE */
    i__1 = n1;
    for (j = js; j <= i__1; ++j) {
	if (j == in) {
	    goto L460;
	}
	z__ = -q[iout + j * q_dim1];
	i__2 = klm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (i__ != iout) {
		q[i__ + j * q_dim1] += z__ * q[i__ + in * q_dim1];
	    }
/* L450: */
	}
L460:
	;
    }
    tpivot = -pivot;
    i__1 = klm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ != iout) {
	    q[i__ + in * q_dim1] /= tpivot;
	}
/* L470: */
    }
    q[iout + in * q_dim1] = 1.f / pivot;
    z__ = q[iout + n2 * q_dim1];
    q[iout + n2 * q_dim1] = q[klm2 + in * q_dim1];
    q[klm2 + in * q_dim1] = z__;
    ii = (long int) abs(z__);
    if (iu[(ii << 1) + 1] == 0 || iu[(ii << 1) + 2] == 0) {
	goto L240;
    }
    i__1 = klm2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__ = q[i__ + in * q_dim1];
	q[i__ + in * q_dim1] = q[i__ + js * q_dim1];
	q[i__ + js * q_dim1] = z__;
/* L480: */
    }
    ++js;
    goto L240;
/* TEST FOR OPTIMALITY. */
L490:
    if (kforce == 0) {
	goto L580;
    }
    if (iphase == 1 && q[klm1 + n1 * q_dim1] <= toler) {
	goto L500;
    }
    kforce = 0;
    goto L240;
/* SET UP PHASE 2 COSTS. */
L500:
    iphase = 2;
    i__1 = nklm;
    for (j = 1; j <= i__1; ++j) {
	cu[(j << 1) + 1] = 0.f;
	cu[(j << 1) + 2] = 0.f;
/* L510: */
    }
    i__1 = nk;
    for (j = n1; j <= i__1; ++j) {
	cu[(j << 1) + 1] = 1.f;
	cu[(j << 1) + 2] = 1.f;
/* L520: */
    }
    i__1 = klm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = (long int) q[i__ + n2 * q_dim1];
	if (ii > 0) {
	    goto L530;
	}
	ii = -ii;
	if (iu[(ii << 1) + 2] == 0) {
	    goto L560;
	}
	cu[(ii << 1) + 2] = 0.f;
	goto L540;
L530:
	if (iu[(ii << 1) + 1] == 0) {
	    goto L560;
	}
	cu[(ii << 1) + 1] = 0.f;
L540:
	++ia;
	i__2 = n2;
	for (j = 1; j <= i__2; ++j) {
	    z__ = q[ia + j * q_dim1];
	    q[ia + j * q_dim1] = q[i__ + j * q_dim1];
	    q[i__ + j * q_dim1] = z__;
/* L550: */
	}
L560:
	;
    }
    goto L160;
L570:
    if (q[klm1 + n1 * q_dim1] <= toler) {
	goto L500;
    }
    kode = 1;
    goto L590;
L580:
    if (iphase == 1) {
	goto L570;
    }
/* PREPARE OUTPUT. */
    kode = 0;
L590:
    xsum = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	x[j] = 0.f;
/* L600: */
    }
    i__1 = klm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	res[i__] = 0.f;
/* L610: */
    }
    i__1 = klm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = (long int) q[i__ + n2 * q_dim1];
	sn = 1.f;
	if (ii > 0) {
	    goto L620;
	}
	ii = -ii;
	sn = -1.f;
L620:
	if (ii > *n) {
	    goto L630;
	}
	x[ii] = sn * q[i__ + n1 * q_dim1];
	goto L640;
L630:
	iimn = ii - *n;
	res[iimn] = sn * q[i__ + n1 * q_dim1];
	if (ii >= n1 && ii <= nk) {
	    xsum += q[i__ + n1 * q_dim1];
	}
L640:
	;
    }
    error = (float)xsum;
    return;
}

float Util::eval(char * images,EMData * img, vector<int> S,int N, int ,int size)
{
	int j,d;
	EMData * e = new EMData();
	float *eptr, *imgptr;
	imgptr = img->get_data();
	float SSE = 0.f;
	for (j = 0 ; j < N ; j++) {
		e->read_image(images,S[j]);
		eptr = e->get_data();
		for (d = 0; d < size; d++) {
			SSE += ((eptr[d] - imgptr[d])*(eptr[d] - imgptr[d]));}
		}
	delete e;
	return SSE;
}


#define		mymax(x,y)		(((x)>(y))?(x):(y))
#define 	mymin(x,y)		(((x)<(y))?(x):(y))
#define		sign(x,y)		(((((y)>0)?(1):(-1))*(y!=0))*(x))


#define		quadpi	 	 	3.141592653589793238462643383279502884197
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
#define ts(i)			ts	[i-1]
#define	thetast(i)		thetast	[i-1]
#define key(i)			key     [i-1]


vector<double> Util::vrdg(const vector<float>& ph, const vector<float>& th)
{

	ENTERFUNC;

	if ( th.size() != ph.size() ) {
		LOGERR("images not same size");
		throw ImageFormatException( "images not same size");
	}

	// rand_seed
	srand(10);

	int i,*key;
	int len = th.size();
	double *theta,*phi,*weight;
	theta   = 	(double*) calloc(len,sizeof(double));
	phi     = 	(double*) calloc(len,sizeof(double));
	weight  = 	(double*) calloc(len,sizeof(double));
	key     =       (int*) calloc(len,sizeof(int));
	const float *thptr, *phptr;

	thptr = &th[0];
	phptr = &ph[0];
	for(i=1;i<=len;i++){
		key(i) = i;
		weight(i) = 0.0;
	}

	for(i = 0;i<len;i++){
		theta[i] = thptr[i];
		phi[i]   = phptr[i];
	}

	//  sort by theta
	Util::hsortd(theta, phi, key, len, 1);

	//Util::voronoidiag(theta,phi, weight, len);
 	Util::voronoi(phi, theta, weight, len);

	//sort by key
	Util::hsortd(weight, weight, key, len, 2);

	free(theta);
	free(phi);
	free(key);
	vector<double> wt;
	double count = 0;
	for(i=1; i<= len; i++)
	{
		wt.push_back(weight(i));
		count += weight(i);
	}

	//if( abs(count-6.28) > 0.1 )
	//{
	//    printf("Warning: SUM OF VORONOI CELLS AREAS IS %lf, should 2*PI\n", count);
	//}

	free(weight);

	EXITFUNC;
	return wt;

}

struct	tmpstruct{
	double theta1,phi1;
	int key1;
	};

void Util::hsortd(double *theta,double *phi,int *key,int len,int option)
{
	ENTERFUNC;
	vector<tmpstruct> tmp(len);
	int i;
	for(i = 1;i<=len;i++)
	{
		tmp[i-1].theta1 = theta(i);
		tmp[i-1].phi1 = phi(i);
		tmp[i-1].key1 = key(i);
	}

	if (option == 1) sort(tmp.begin(),tmp.end(),Util::cmp1);
	if (option == 2) sort(tmp.begin(),tmp.end(),Util::cmp2);

	for(i = 1;i<=len;i++)
	{
		theta(i) = tmp[i-1].theta1;
		phi(i) 	 = tmp[i-1].phi1;
		key(i)   = tmp[i-1].key1;
	}
	EXITFUNC;
}

bool Util::cmp1(tmpstruct tmp1,tmpstruct tmp2)
{
 	return(tmp1.theta1 < tmp2.theta1);
}

bool Util::cmp2(tmpstruct tmp1,tmpstruct tmp2)
{
	return(tmp1.key1 < tmp2.key1);
}

/******************  VORONOI DIAGRAM **********************************/
/*
void Util::voronoidiag(double *theta,double *phi,double* weight,int n)
{
	ENTERFUNC;

	int	*lband;
	double	aat=0.0f,*ts;
	double	aa,acum,area;
	int	last;
	int numth 	= 	1;
	int nbt		=	1;//mymax((int)(sqrt((n/500.0))) , 3);

	int i,it,l,k;
	int nband,lb,low,medium,lhigh,lbw,lenw;


	lband 	=	(int*)calloc(nbt,sizeof(int));
	ts 	=	(double*)calloc(nbt,sizeof(double));

	if(lband == NULL || ts == NULL ){
       		fprintf(stderr,"memory allocation failure!\n");
		exit(1);
	}

	nband=nbt;

	while(nband>0){
		Util::angstep(ts,nband);

		l=1;
		for(i=1;i<=n;i++){
			if(theta(i)>ts(l)){
				lband(l)=i;
				l=l+1;
				if(l>nband)  exit(1);
			}
		}

		l=1;
		for(i=1;i<=n;i++){
			if(theta(i)>ts(l)){
				lband(l)=i;
				l=l+1;
				if(l>nband)  exit(1);
			}
		}

		lband(l)=n+1;
 		acum=0.0;
		for(it=l;it>=1;it-=numth){
			for(i=it;i>=mymax(1,it-numth+1);i--){
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


 			Util::voronoi(&phi(lb),&theta(lb),&weight(lbw),lenw,low,medium,lhigh,last);


			if(nband>1){
				if(i==1)	area=quadpi*2.0*(1.0-cos(ts(1)*dgr_to_rad));
				else		area=quadpi*2.0*(cos(ts(i-1)*dgr_to_rad)-cos(ts(i)*dgr_to_rad));

				aa = 0.0;
				for(k = lbw;k<=lbw+lenw-1;k++)
					aa = aa+weight(k);

				acum=acum+aa;
		   		aat=aa/area;
				}

			}
			for(i=it;mymax(1,it-numth+1);i--){
			if(fabs(aat-1.0)>0.02){
				nband=mymax(0,mymin( (int)(((float)nband) * 0.75) ,nband-1) );
				goto  label2;
				}
			}
		acum=acum/quadpi/2.0;
		exit(1);
label2:

		continue;
		}

	free(ts);
	free(lband);

	}

	EXITFUNC;
}


void Util::angstep(double* thetast,int len){

	ENTERFUNC;

	double t1,t2,tmp;
	int i;
 	if(len>1){
		t1=0;
		for(i=1;i<=len-1;i++){
			tmp=cos(t1)-1.0/((float)len);
			t2=acos(sign(mymin(1.0,fabs(tmp)),tmp));
			thetast(i)=t2 * rad_to_deg;
			t1=t2;
		}
	}
	thetast(len)=90.0;

	EXITFUNC;
}
*/
/*
void Util::voronoi(double *phi, double *theta, double *weight, int lenw, int low, int medium, int nt, int last)
{

	ENTERFUNC;
	int *list, *lptr, *lend, *iwk, *key,*lcnt,*indx,*good;
	int nt6, n, ier,nout,lnew,mdup,nd;
	int i,k,mt,status;


	double *ds, *x, *y, *z;
	double tol=1.0e-8;
	double a;

	if(last){
		if(medium>nt)  n = nt+nt;
		else           n = nt+nt-medium+1;
	}
	else{
		n=nt;
	}

	nt6 = n*6;

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
       		printf("memory allocation failure!\n");
       		exit(1);
	}



	for(i = 1;i<=nt;i++){
		x[i-1] = theta(i);
		y[i-1] = phi(i);
	}



	if (last) {
		for(i=nt+1;i<=n;i++){
			x[i-1]=180.0-x[2*nt-i];
			y[i-1]=180.0+y[2*nt-i];
		}
	}


	Util::disorder2(x,y,key,n);

	Util::ang_to_xyz(x,y,z,n);


	//  Make sure that first three are no colinear
	label1:
	for(k=0; k<2; k++){
		for(i=k+1; i<3; i++){
			if(  x[i]*x[k]+y[i]*y[k]+z[i]*z[k] > 1.0-tol){
				Util::flip23(x, y, z, key, k, n);
				goto label1;
			}
		}
	}


	status = Util::trmsh3_(&n,&tol,x,y,z,&nout,list,lptr,lend,&lnew,indx,lcnt, iwk, good, ds, &ier);


	if (status != 0) {
		printf(" error in trmsh3 \n");
		exit(1);
	}


	mdup=n-nout;
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
			good[nd-1]=k;
		}
	}


	for(i = 1;i<=nout;i++) {
		k=good[i-1];
 		if (key[k-1] >= low && key[k-1]<medium){
			a = Util::areav_(&i,&nout,x,y,z,list,lptr,lend,&ier);
    			if (ier != 0){
				weight[key[k-1]-low] =-1.0;
			}
			else {
				weight[key[k-1]-low]=a/lcnt[i-1];
			}
		}
	}

// Fill out the duplicated weights
	for(i = 1;i<=n;i++){
		mt=-indx[i-1];
		if (mt>0){
			k=good[mt-1];
//  This is a duplicated entry, get the already calculated
//   weight and assign it.
			if (key[i-1]>=low && key[i-1]<medium){
//  Is it already calculated weight??
				if(key[k-1]>=low && key[k-1]<medium){
					weight[key[i-1]-low]=weight[key[k-1]-low];
				}
				else{
//  No, the weight is from the outside of valid region, calculate it anyway
					a = Util::areav_(&mt, &nout, x, y, z, list, lptr, lend, &ier);
					if (ier != 0){
						printf("    *** error in areav:  ier = %d ***\n", ier);
						weight[key[i-1]-low] =-1.0;
					}
					else {
						weight[key[i-1]-low] = a/lcnt[mt-1];
					}
				}
			}
		}
	}


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
	EXITFUNC;
}
*/
void Util::voronoi(double *phi, double *theta, double *weight, int nt)
{

	ENTERFUNC;

	int *list, *lptr, *lend, *iwk, *key,*lcnt,*indx,*good;
	int nt6, n, ier, nout, lnew, mdup, nd;
	int i,k,mt,status;


	double *ds, *x, *y, *z;
	double tol  = 1.0e-8;
	double dtol = 15;
	double a;

	/*if(last){
		if(medium>nt)  n = nt+nt;
		else           n = nt+nt-medium+1;
	}
	else{
		n=nt;
	}*/

	n = nt + nt;

	nt6 = n*6;

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
       		printf("memory allocation failure!\n");
       		exit(1);
	}

        bool colinear=true;
	while(colinear)
	{

	L1:
 	    for(i = 0; i<nt; i++){
		x[i] = theta[i];
		y[i] = phi[i];
		x[nt+i] = 180.0 - x[i];
		y[nt+i] = 180.0 + y[i];
	    }

	    Util::disorder2(x, y, key, n);

	    // check if the first three angles are not close, else shuffle
	    double val;
            for(k=0; k<2; k++){
		for(i=k+1; i<3; i++){
		    val = (x[i]-x[k])*(x[i]-x[k]) + (y[i]-y[k])*(y[i]-y[k]);
		    if( val  < dtol) {
			goto L1;
		    }
		}
	    }

	    Util::ang_to_xyz(x, y, z, n);

	    //  Make sure that first three has no duplication
	    bool dupnode=true;
	    dupnode=true;
  	    while(dupnode)
	    {
	        for(k=0; k<2; k++){
		    for(i=k+1; i<3; i++){
			if(  x[i]*x[k]+y[i]*y[k]+z[i]*z[k] > 1.0-tol) {
				Util::flip23(x, y, z, key, k, n);
				continue;
			}
		    }
	        }
	        dupnode = false;
	    }


       	    ier = 0;

	    status = Util::trmsh3_(&n,&tol,x,y,z,&nout,list,lptr,lend,&lnew, indx, lcnt, iwk, good, ds, &ier);

	    if (status != 0) {
	 	printf(" error in trmsh3 \n");
		exit(1);
	    }

            if (ier > 0) {
		printf("*** Error in TRMESH:  duplicate nodes encountered ***\n");
		exit(1);
	    }

	    mdup=n-nout;
	    if (ier == -2) {
		//printf("in TRMESH:the first three nodes are colinear*** disorder again\n");
	    }
	    else
	    {
	        colinear=false;
	    }
        }


	Assert( ier != -2 );
//  Create a list of unique nodes GOOD, the numbers refer to locations on the full list
//  INDX contains node numbers from the squeezed list
	nd=0;
	for (k=1; k<=n; k++){
		if (indx[k-1]>0) {
			nd++;
			good[nd-1]=k;
		}
	}

//
// *** Compute the Voronoi region areas.
//
	for(i = 1; i<=nout; i++) {
		k=good[i-1];
		//  We only need n weights from hemisphere
		if (key[k-1] <= nt) {
//  CALCULATE THE AREA
 			a = Util::areav_(&i, &nout, x, y, z, list, lptr, lend, &ier);
    			if (ier != 0){
//  We set the weight to -1, this will signal the error in the calling
//   program, as the area will turn out incorrect
				printf("    *** error in areav:  ier = %d ***\n", ier);
				weight[key[k-1]-1] =-1.0;
			} else {
//  Assign the weight
				weight[key[k-1]-1]=a/lcnt[i-1];
			}
		}
	}


// Fill out the duplicated weights
	for(i = 1; i<=n; i++){
		mt -= indx[i-1];
		if (mt>0){
			k = good[mt-1];
//  This is a duplicated entry, get the already calculated
//   weight and assign it.
		//  We only need n weights from hemisphere
			if (key[i-1] <= nt && key[k-1] <= nt) { weight[key[i-1]-1] = weight[key[k-1]-1];}
			}
	}

	free(list);
	free(lend);
	free(iwk);
	free(good);
	free(key);
	free(lptr);
	free(indx);
	free(lcnt);
	free(ds);
	free(x);
	free(y);
	free(z);


	EXITFUNC;
}

void Util::disorder2(double *x,double *y, int *key, int len)
{
	ENTERFUNC;
	int k, i;
	for(i=0; i<len; i++) key[i]=i+1;

	for(i = 0; i<len;i++){
		k = rand()%len;
		std::swap(key[k], key[i]);
		std::swap(x[k], x[i]);
		std::swap(y[k], y[i]);
	}
	EXITFUNC;
}

void Util::ang_to_xyz(double *x,double *y,double *z,int len)
{
	ENTERFUNC;
	double costheta,sintheta,cosphi,sinphi;
	for(int i = 0;  i<len;  i++) {
		cosphi = cos(y[i]*dgr_to_rad);
		sinphi = sin(y[i]*dgr_to_rad);
		if(fabs(x[i]-90.0)< 1.0e-5){
			x[i] = cosphi;
			y[i] = sinphi;
			z[i] = 0.0;
		} else {
			costheta = cos(x[i]*dgr_to_rad);
			sintheta = sin(x[i]*dgr_to_rad);
			x[i] = cosphi*sintheta;
			y[i] = sinphi*sintheta;
			z[i] = costheta;
		}
	}
	EXITFUNC;
}

void Util::flip23(double *x,double *y,double *z,int *key, int k, int len)
{
	ENTERFUNC;
	int i = k;
	while( i == k )  i = rand()%len;
	std::swap(key[i], key[k]);
	std::swap(x[i], x[k]);
	std::swap(y[i], y[k]);
	std::swap(z[i], z[k]);
	EXITFUNC;
}


#undef	mymax
#undef	mymin
#undef	sign
#undef	quadpi
#undef	dgr_to_rad
#undef	deg_to_rad
#undef	rad_to_deg
#undef	rad_to_dgr
#undef	TRUE
#undef	FALSE
#undef 	theta
#undef 	phi
#undef 	weight
#undef 	lband
#undef 	ts
#undef  thetast
#undef 	key


/*################################################################################################
##########  strid.f -- translated by f2c (version 20030320). ###################################
######   You must link the resulting object file with the libraries: #############################
####################	-lf2c -lm   (in that order)   ############################################
################################################################################################*/

/* Common Block Declarations */


#define TRUE_ (1)
#define FALSE_ (0)
#define abs(x) ((x) >= 0 ? (x) : -(x))

struct stcom_{
    double y;
};
stcom_ stcom_1;
#ifdef KR_headers
double floor();
int i_dnnt(x) double *x;
#else
int i_dnnt(double *x)
#endif
{
	return (int)(*x >= 0. ? floor(*x + .5) : -floor(.5 - *x));
}




/* ____________________STRID______________________________________ */
/* Subroutine */ int Util::trmsh3_(int *n0, double *tol, double *x,
	double *y, double *z__, int *n, int *list, int *
	lptr, int *lend, int *lnew, int *indx, int *lcnt,
	int *near__, int *next, double *dist, int *ier)
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static double d__;
    static int i__, j;
    static double d1, d2, d3;
    static int i0, lp, kt, ku, lpl, nku;
    static int nexti;


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

/*       LCNT = Array of int weights (counts) associated */
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

    //printf("pass check colinear\n");

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
/* Subroutine */ int addnod_(int *nst, int *k, double *x,
	double *y, double *z__, int *list, int *lptr, int
	*lend, int *lnew, int *ier)
{
    /* Initialized data */

    static double tol = 0.;

    /* System generated locals */
    int i__1;

    /* Local variables */
    static int l;
    static double p[3], b1, b2, b3;
    static int i1, i2, i3, kk, lp, in1, io1, io2, km1, lpf, ist, lpo1;
    /* Subroutine */ int swap_(int *, int *, int *,
	    int *, int *, int *, int *, int *);
    static int lpo1s;
    /* Subroutine */ int bdyadd_(int *, int *, int *,
	    int *, int *, int *, int *), intadd_(int *,
	    int *, int *, int *, int *, int *, int *,
	    int *), trfind_(int *, double *, int *,
	    double *, double *, double *, int *, int *,
	    int *, double *, double *, double *, int *,
	    int *, int *), covsph_(int *, int *, int *,
	    int *, int *, int *);
    int lstptr_(int *, int *, int *, int *);
    long int swptst_(int *, int *, int *, int *,
	    double *, double *, double *);


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

double angle_(double *v1, double *v2, double *v3)
{
    /* System generated locals */
    double ret_val;

    /* Builtin functions */
    //double sqrt(double), acos(double);

    /* Local variables */
    static double a;
    static int i__;
    static double ca, s21, s23, u21[3], u23[3];


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

double areas_(double *v1, double *v2, double *v3)
{
    /* System generated locals */
    double ret_val;

    /* Builtin functions */
    //double sqrt(double), acos(double);

    /* Local variables */
    static int i__;
    static double a1, a2, a3, s12, s31, s23, u12[3], u23[3], u31[3], ca1,
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

//double areas_(double *, double *, double *);

double Util::areav_(int *k, int *n, double *x, double *y,
	double *z__, int *list, int *lptr, int *lend, int
	*ier)
{
    /* Initialized data */

    static double amax = 6.28;

    /* System generated locals */
    double ret_val;

    /* Local variables */
    static double a, c0[3], c2[3], c3[3];
    static int n1, n2, n3;
    static double v1[3], v2[3], v3[3];
    static int lp, lpl, ierr;
    static double asum;
    static long int first;


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

double areav_new__(int *k, int *n, double *x, double *y,
	double *z__, int *list, int *lptr, int *lend, int
	*ier)
{
    /* System generated locals */
    double ret_val = 0;

    /* Builtin functions */
    //double acos(double);

    /* Local variables */
    static int m;
    static double c1[3], c2[3], c3[3];
    static int n1, n2, n3;
    static double v1[3], v2[3], v3[3];
    static int lp;
    static double c1s[3], c2s[3];
    static int lpl, ierr;
    static double asum;
    double angle_(double *, double *, double *);
    static float areav;


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
    ret_val = asum - (double) (m - 2) * acos(-1.);
    return ret_val;

/* Invalid input. */

L11:
    *ier = 1;
    areav = 0.f;
    return ret_val;

/* K indexes a boundary node. */

L12:
    *ier = 2;
    areav = 0.f;
    return ret_val;

/* Error in CIRCUM. */

L13:
    *ier = 3;
    areav = 0.f;
    return ret_val;
} /* areav_new__ */

/* Subroutine */ int bdyadd_(int *kk, int *i1, int *i2, int *
	list, int *lptr, int *lend, int *lnew)
{
    static int k, n1, n2, lp, lsav, nsav, next;
    /* Subroutine */ int insert_(int *, int *, int *,
	    int *, int *);


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

/* Subroutine */ int bnodes_(int *n, int *list, int *lptr,
	int *lend, int *nodes, int *nb, int *na, int *nt)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int k, n0, lp, nn, nst;


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

/*       NODES = int array of length at least NB */
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
    *nt = nn - (2<<1);
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

/* Subroutine */ int circle_(int *k, double *xc, double *yc,
	int *ier)
{
    /* System generated locals */
    int i__1;

    /* Builtin functions */
    //double atan(double), cos(double), sin(double);

    /* Local variables */
    static double a, c__;
    static int i__;
    static double s;
    static int k2, k3;
    static double x0, y0;
    static int kk, np1;


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
    a = atan(1.) * 2. / (double) kk;
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

/* Subroutine */ int circum_(double *v1, double *v2, double *v3,
	double *c__, int *ier)
{
    /* Builtin functions */
    //double sqrt(double);

    /* Local variables */
    static int i__;
    static double e1[3], e2[3], cu[3], cnorm;


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

/* Subroutine */ int covsph_(int *kk, int *n0, int *list, int
	*lptr, int *lend, int *lnew)
{
    static int k, lp, nst, lsav, next;
    /* Subroutine */ int insert_(int *, int *, int *,
	    int *, int *);


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


/* Subroutine */ int crlist_(int *n, int *ncol, double *x,
	double *y, double *z__, int *list, int *lend, int
	*lptr, int *lnew, int *ltri, int *listc, int *nb,
	double *xc, double *yc, double *zc, double *rc,
	int *ier)
{
    /* System generated locals */
    int i__1, i__2;

    /* Builtin functions */
    //double acos(double);

    /* Local variables */
    static double c__[3], t;
    static int i1, i2, i3, i4, n0, n1, n2, n3, n4;
    static double v1[3], v2[3], v3[3];
    static int lp, kt, nn, nt, nm2, kt1, kt2, kt11, kt12, kt21, kt22, lpl,
	     lpn;
    static long int swp;
    static int ierr;
    int lstptr_(int *, int *, int *, int *);
    long int swptst_(int *, int *, int *, int *,
	    double *, double *, double *);


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

/*       LIST = int array containing the set of adjacency */
/*              lists.  Refer to Subroutine TRMESH. */

/*       LEND = Set of pointers to ends of adjacency lists. */
/*              Refer to Subroutine TRMESH. */

/* The above parameters are not altered by this routine. */

/*       LPTR = Array of pointers associated with LIST.  Re- */
/*              fer to Subroutine TRMESH. */

/*       LNEW = Pointer to the first empty location in LIST */
/*              and LPTR (list length plus one). */

/*       LTRI = int work space array dimensioned 6 by */
/*              NCOL, or unused dummy parameter if NB = 0. */

/*       LISTC = int array of length at least 3*NT, where */
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
/* SWP =       long int variable set to TRUE in each optimiza- */
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

/* Subroutine */ int delarc_(int *n, int *io1, int *io2, int *
	list, int *lptr, int *lend, int *lnew, int *ier)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int n1, n2, n3, lp, lph, lpl;
    /* Subroutine */ int delnb_(int *, int *, int *,
	    int *, int *, int *, int *, int *);
    int lstptr_(int *, int *, int *, int *);


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

/* Subroutine */ int delnb_(int *n0, int *nb, int *n, int *
	list, int *lptr, int *lend, int *lnew, int *lph)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__, lp, nn, lpb, lpl, lpp, lnw;


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

/* Subroutine */ int delnod_(int *k, int *n, double *x,
	double *y, double *z__, int *list, int *lptr, int
	*lend, int *lnew, int *lwk, int *iwk, int *ier)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__, j, n1, n2;
    static double x1, x2, y1, y2, z1, z2;
    static int nl, lp, nn, nr;
    static double xl, yl, zl, xr, yr, zr;
    static int nnb, lp21, lpf, lph, lpl, lpn, iwl, nit, lnw, lpl2;
    static long int bdry;
    static int ierr, lwkl;
    /* Subroutine */ int swap_(int *, int *, int *,
	    int *, int *, int *, int *, int *), delnb_(
	    int *, int *, int *, int *, int *, int *,
	    int *, int *);
    int nbcnt_(int *, int *);
    /* Subroutine */ int optim_(double *, double *, double
	    *, int *, int *, int *, int *, int *, int
	    *, int *);
    static int nfrst;
    int lstptr_(int *, int *, int *, int *);


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

/*       IWK = int work array dimensioned 2 by LWK (or */
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

/* BDRY =    long int variable with value TRUE iff N1 is a */
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
/*      WRITE (*,100) NIT, IERR */
/*  100 FORMAT (//5X,'*** Error in OPTIM (called from ', */
/*     .        'DELNOD):  NIT = ',I4,', IER = ',I1,' ***'/) */
    return 0;

/* Error flag 1 returned by OPTIM. */

L26:
    *ier = 6;
    return 0;
} /* delnod_ */

/* Subroutine */ int drwarc_(int *, double *p, double *q,
	double *tol, int *nseg)
{
    /* System generated locals */
    int i__1;
    double d__1;

    /* Builtin functions */
    //double sqrt(double);

    /* Local variables */
    static int i__, k;
    static double s, p1[3], p2[3], u1, u2, v1, v2;
    static int na;
    static double dp[3], du, dv, pm[3], um, vm, err, enrm;


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

/*       LUN = long int unit number in the range 0 to 99. */

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

    na = (int) (err / *tol + 1.);

/* Initialize for loop on arcs P1-P2, where the intermediate */
/*   points are obtained by normalizing PM = P + k*DP for */
/*   DP = (Q-P)/NA */

    s = 1. / (double) na;
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
/*        WRITE (LUN,100,ERR=5) P1(1), P1(2), P2(1), P2(2) */
/*  100   FORMAT (2F12.6,' moveto',2F12.6,' lineto') */
	p1[0] = p2[0];
	p1[1] = p2[1];
	p1[2] = p2[2];
/* L4: */
    }
/*      WRITE (LUN,100,ERR=5) P1(1), P1(2), Q(1), Q(2) */

/* No error encountered. */

    *nseg = na;
    return 0;

/* Invalid input value of P or Q. */

L5:
    *nseg = 0;
    return 0;
} /* drwarc_ */

/* Subroutine */ int edge_(int *in1, int *in2, double *x,
	double *y, double *z__, int *lwk, int *iwk, int *
	list, int *lptr, int *lend, int *ier)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__, n0, n1, n2;
    static double x0, x1, x2, y0, y1, y2, z0, z1, z2;
    static int nl, lp, nr;
    static double dp12;
    static int lp21, iwc, iwf, lft, lpl, iwl, nit;
    static double dp1l, dp2l, dp1r, dp2r;
    static int ierr;
    /* Subroutine */ int swap_(int *, int *, int *,
	    int *, int *, int *, int *, int *);
    static int next, iwcp1, n1lst, iwend;
    /* Subroutine */ int optim_(double *, double *, double
	    *, int *, int *, int *, int *, int *, int
	    *, int *);
    static int n1frst;


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

/*       IWK = int work array of length at least 2*LWK. */

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
	iwk[(i__ - (1<<1)) + 1] = iwk[(i__ << 1) + 1];
	iwk[(i__ - (1<<1)) + 2] = iwk[(i__ << 1) + 2];
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
	iwk[(i__ + (1<<1)) + 1] = iwk[(i__ << 1) + 1];
	iwk[(i__ + (1<<1)) + 2] = iwk[(i__ << 1) + 2];
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
    iwk[(i__ << 1) + 1] = iwk[(i__ - (1<<1)) + 1];
    iwk[(i__ << 1) + 2] = iwk[(i__ - (1<<1)) + 2];
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

	nit = iwc - (1<<2);
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

	nit = iwend - (iwc<<2);
	i__1 = iwend - iwc;
	optim_(&x[1], &y[1], &z__[1], &i__1, &list[1], &lptr[1], &lend[1], &
		nit, &iwk[(iwc + (1<<1)) + 1], &ierr);
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
/*      WRITE (*,130) IN1, IN2 */
/*  130 FORMAT (//5X,'*** Error in EDGE:  Invalid triangula', */
/*     .        'tion or null triangles on boundary'/ */
/*     .        9X,'IN1 =',I4,', IN2=',I4/) */
    return 0;

/* Error flag (other than 1) returned by OPTIM. */

L34:
    *ier = 4;
/*      WRITE (*,140) NIT, IERR */
/*  140 FORMAT (//5X,'*** Error in OPTIM (called from EDGE):', */
/*     .        '  NIT = ',I4,', IER = ',I1,' ***'/) */
    return 0;

/* Error flag 1 returned by OPTIM. */

L35:
    *ier = 5;
    return 0;
} /* edge_ */

/* Subroutine */ int getnp_(double *x, double *y, double *z__,
	int *list, int *lptr, int *lend, int *l, int *
	npts, double *df, int *ier)
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static int i__, n1;
    static double x1, y1, z1;
    static int nb, ni, lp, np, lm1;
    static double dnb, dnp;
    static int lpl;


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

/* Subroutine */ int insert_(int *k, int *lp, int *list, int *
	lptr, int *lnew)
{
    static int lsav;


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

long int inside_(double *p, int *lv, double *xv, double *yv,
	double *zv, int *nv, int *listv, int *ier)
{
    /* Initialized data */

    static double eps = .001;

    /* System generated locals */
    int i__1;
    long int ret_val = 0;

    /* Builtin functions */
    //double sqrt(double);

    /* Local variables */
    static double b[3], d__;
    static int k, n;
    static double q[3];
    static int i1, i2, k0;
    static double v1[3], v2[3], cn[3], bp, bq;
    static int ni;
    static double pn[3], qn[3], vn[3];
    static int imx;
    static long int lft1, lft2, even;
    static int ierr;
    static long int pinr, qinr;
    static double qnrm, vnrm;
    /* Subroutine */ int intrsc_(double *, double *,
	    double *, double *, int *);


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
/* LFT1,LFT2 = long int variables associated with I1 and I2 in */
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

/* Subroutine */ int intadd_(int *kk, int *i1, int *i2, int *
	i3, int *list, int *lptr, int *lend, int *lnew)
{
    static int k, n1, n2, n3, lp;
    /* Subroutine */ int insert_(int *, int *, int *,
	    int *, int *);
    int lstptr_(int *, int *, int *, int *);


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

/* Subroutine */ int intrsc_(double *p1, double *p2, double *cn,
	double *p, int *ier)
{
    /* Builtin functions */
    //double sqrt(double);

    /* Local variables */
    static int i__;
    static double t, d1, d2, pp[3], ppn;


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

int jrand_(int *n, int *ix, int *iy, int *iz)
{
    /* System generated locals */
    int ret_val;

    /* Local variables */
    static float u, x;


/* *********************************************************** */

/*                                              From STRIPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/28/98 */

/*   This function returns a uniformly distributed pseudo- */
/* random int in the range 1 to N. */


/* On input: */

/*       N = Maximum value to be returned. */

/* N is not altered by this function. */

/*       IX,IY,IZ = int seeds initialized to values in */
/*                  the range 1 to 30,000 before the first */
/*                  call to JRAND, and not altered between */
/*                  subsequent calls (unless a sequence of */
/*                  random numbers is to be repeated by */
/*                  reinitializing the seeds). */

/* On output: */

/*       IX,IY,IZ = Updated int seeds. */

/*       JRAND = Random int in the range 1 to N. */

/* Reference:  B. A. Wichmann and I. D. Hill, "An Efficient */
/*             and Portable Pseudo-random Number Generator", */
/*             Applied Statistics, Vol. 31, No. 2, 1982, */
/*             pp. 188-190. */

/* Modules required by JRAND:  None */

/* Intrinsic functions called by JRAND:  INT, MOD, float */

/* *********************************************************** */


/* Local parameters: */

/* U = Pseudo-random number uniformly distributed in the */
/*     interval (0,1). */
/* X = Pseudo-random number in the range 0 to 3 whose frac- */
/*       tional part is U. */

    *ix = *ix * 171 % 30269;
    *iy = *iy * 172 % 30307;
    *iz = *iz * 170 % 30323;
    x = (float) (*ix) / 30269.f + (float) (*iy) / 30307.f + (float) (*iz) /
	    30323.f;
    u = x - (int) x;
    ret_val = (int) ((float) (*n) * u + 1.f);
    return ret_val;
} /* jrand_ */

long int left_(double *x1, double *y1, double *z1, double *x2,
	double *y2, double *z2, double *x0, double *y0,
	double *z0)
{
    /* System generated locals */
    long int ret_val;


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
    	    *z0 * (*x1 * *y2 - *x2 * *y1) >= -0.000001;


    return ret_val;
} /* left_ */

int lstptr_(int *lpl, int *nb, int *list, int *lptr)
{
    /* System generated locals */
    int ret_val;

    /* Local variables */
    static int nd, lp;


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

int nbcnt_(int *lpl, int *lptr)
{
    /* System generated locals */
    int ret_val;

    /* Local variables */
    static int k, lp;


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

int nearnd_(double *p, int *ist, int *n, double *x,
	double *y, double *z__, int *list, int *lptr, int
	*lend, double *al)
{
    /* System generated locals */
    int ret_val, i__1;

    /* Builtin functions */
    //double acos(double);

    /* Local variables */
    static int l;
    static double b1, b2, b3;
    static int i1, i2, i3, n1, n2, n3, lp, nn, nr;
    static double ds1;
    static int lp1, lp2;
    static double dx1, dx2, dx3, dy1, dy2, dy3, dz1, dz2, dz3;
    static int lpl;
    static double dsr;
    static int nst, listp[25], lptrp[25];
    /* Subroutine */ int trfind_(int *, double *, int *,
	    double *, double *, double *, int *, int *,
	    int *, double *, double *, double *, int *,
	    int *, int *);
    int lstptr_(int *, int *, int *, int *);


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

/* Subroutine */ int optim_(double *x, double *y, double *z__,
	int *na, int *list, int *lptr, int *lend, int *
	nit, int *iwk, int *ier)
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static int i__, n1, n2, lp, io1, io2, nna, lp21, lpl, lpp;
    static long int swp;
    static int iter;
    /* Subroutine */ int swap_(int *, int *, int *,
	    int *, int *, int *, int *, int *);
    static int maxit;
    long int swptst_(int *, int *, int *, int *,
	    double *, double *, double *);


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

/*       IWK = int array dimensioned 2 by NA containing */
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

/* Subroutine */ int projct_(double *px, double *py, double *pz,
	double *ox, double *oy, double *oz, double *ex,
	double *ey, double *ez, double *vx, double *vy,
	double *vz, long int *init, double *x, double *y,
	double *z__, int *ier)
{
    /* Builtin functions */
    //double sqrt(double);

    /* Local variables */
    static double s, sc, xe, ye, ze, xh, yh, zh, xv, yv, zv, xw, yw, zw,
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

/*       INIT = long int switch which must be set to TRUE on */
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

/* Subroutine */ int scoord_(double *px, double *py, double *pz,
	double *plat, double *plon, double *pnrm)
{
    /* Builtin functions */
    //double sqrt(double), atan2(double, double), asin(double);


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

double store_(double *x)
{
    /* System generated locals */
    double ret_val;


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

/* Subroutine */ int swap_(int *in1, int *in2, int *io1, int *
	io2, int *list, int *lptr, int *lend, int *lp21)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int lp, lph, lpsav;
    int lstptr_(int *, int *, int *, int *);


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

long int swptst_(int *n1, int *n2, int *n3, int *n4,
	double *x, double *y, double *z__)
{
    /* System generated locals */
    long int ret_val;

    /* Local variables */
    static double x4, y4, z4, dx1, dx2, dx3, dy1, dy2, dy3, dz1, dz2, dz3;


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

/* Subroutine */ int trans_(int *n, double *rlat, double *rlon,
	double *x, double *y, double *z__)
{
    /* System generated locals */
    int i__1;

    /* Builtin functions */
    //double cos(double), sin(double);

    /* Local variables */
    static int i__, nn;
    static double phi, theta, cosphi;


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

/* Subroutine */ int trfind_(int *nst, double *p, int *n,
	double *x, double *y, double *z__, int *list, int
	*lptr, int *lend, double *b1, double *b2, double *b3,
	int *i1, int *i2, int *i3)
{
    /* Initialized data */

    static int ix = 1;
    static int iy = 2;
    static int iz = 3;

    /* System generated locals */
    int i__1;
    double d__1, d__2;

    /* Local variables */
    static double q[3];
    static int n0, n1, n2, n3, n4, nf;
    static double s12;
    static int nl, lp;
    static double xp, yp, zp;
    static int n1s, n2s;
    static double eps, tol, ptn1, ptn2;
    static int next;
    int jrand_(int *, int *, int *, int *);
    double store_(double *);
    int lstptr_(int *, int *, int *, int *);


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
/* IX,IY,IZ = int seeds for JRAND */
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
		-1e-10) {
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
		-1e-10) {

/*   P is to the right of the boundary edge N0->NF. */

	    n1 = n0;
	    n2 = nf;
	    goto L9;
	}
	if (xp * (y[nl] * z__[n0] - y[n0] * z__[nl]) - yp * (x[nl] * z__[n0]
		- x[n0] * z__[nl]) + zp * (x[nl] * y[n0] - x[n0] * y[nl]) <
		-1e-10) {

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
	    n2] * z__[n0]) + zp * (x[n0] * y[n2] - x[n2] * y[n0]) < -1e-10) {
	goto L7;
    }
    n1 = n2;
    if (n1 != nl) {
	goto L4;
    }
    if (xp * (y[n0] * z__[nf] - y[nf] * z__[n0]) - yp * (x[n0] * z__[nf] - x[
	    nf] * z__[n0]) + zp * (x[n0] * y[nf] - x[nf] * y[n0]) < -1e-10) {
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
		- x[n0] * z__[n1]) + zp * (x[n1] * y[n0] - x[n0] * y[n1]) >
		-1e-10) {
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
     if (*b3 < -1e-10) {

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
		-1e-10) {
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
    if (*b1 < 0.f) {
	*b1 = 0.f;
    }
    if (*b2 < 0.f) {
	*b2 = 0.f;
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
	    >= -1e-10) {

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
		 y[next]) >= -1e-10) {

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

/* Subroutine */ int trlist_(int *n, int *list, int *lptr,
	int *lend, int *nrow, int *nt, int *ltri, int *
	ier)
{
    /* System generated locals */
    int ltri_dim1, ltri_offset, i__1, i__2;

    /* Local variables */
    static int i__, j, i1, i2, i3, n1, n2, n3, ka, kn, lp, kt, nm2, lp2,
	    lpl, isv;
    static long int arcs;
    static int lpln1;


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

/*       LTRI = int array of length at least NROW*NT, */
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

/* ARCS =     long int variable with value TRUE iff are */
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
    ltri_offset = 1 + ltri_dim1;
    ltri -= ltri_offset;

    /* Function Body */
    if (*n < 3 || (*nrow != 6 && *nrow != 9)) {
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

/* Subroutine */ int trlprt_(int *n, double *x, double *y,
	double *z__, int *iflag, int *nrow, int *nt, int *
	ltri, int *lout)
{
    /* Initialized data */

    static int nmax = 9999;
    static int nlmax = 58;

    /* System generated locals */
    int ltri_dim1, ltri_offset, i__1;

    /* Local variables */
    static int i__, k, na, nb, nl, lun;


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
/* on long int unit LOUT.  The numbers of boundary nodes, */
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

/*       LOUT = long int unit number for output.  If LOUT is */
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
    ltri_offset = 1 + ltri_dim1;
    ltri -= ltri_offset;

    /* Function Body */

/* Local parameters: */

/* I =     DO-loop, nodal index, and row index for LTRI */
/* K =     DO-loop and triangle index */
/* LUN =   long int unit number for output */
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

/*      WRITE (LUN,100) N */
    nl = 3;
    if (*n < 3 || *n > nmax || (*nrow != 6 && *nrow != 9) || *nt < 1 || *nt >
	    nmax) {

/* Print an error message and exit. */

/*        WRITE (LUN,110) N, NROW, NT */
	return 0;
    }
    if (*iflag == 0) {

/* Print X, Y, and Z. */

/*        WRITE (LUN,101) */
	nl = 6;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (nl >= nlmax) {
/*            WRITE (LUN,108) */
		nl = 0;
	    }
/*          WRITE (LUN,103) I, X(I), Y(I), Z(I) */
	    ++nl;
/* L1: */
	}
    } else if (*iflag > 0) {

/* Print X (longitude) and Y (latitude). */

/*        WRITE (LUN,102) */
	nl = 6;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (nl >= nlmax) {
/*            WRITE (LUN,108) */
		nl = 0;
	    }
/*          WRITE (LUN,104) I, X(I), Y(I) */
	    ++nl;
/* L2: */
	}
    }

/* Print the triangulation LTRI. */

    if (nl > nlmax / 2) {
/*        WRITE (LUN,108) */
	nl = 0;
    }
    if (*nrow == 6) {
/*        WRITE (LUN,105) */
    } else {
/*        WRITE (LUN,106) */
    }
    nl += 5;
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	if (nl >= nlmax) {
/*          WRITE (LUN,108) */
	    nl = 0;
	}
/*        WRITE (LUN,107) K, (LTRI(I,K), I = 1,NROW) */
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
/*      WRITE (LUN,109) NB, NA, NT */
    return 0;

/* Print formats: */

/*  100 FORMAT (///18X,'STRIPACK (TRLIST) Output,  N = ',I4) */
/*  101 FORMAT (//8X,'Node',10X,'X(Node)',10X,'Y(Node)',10X, */
/*     .        'Z(Node)'//) */
/*  102 FORMAT (//16X,'Node',8X,'Longitude',9X,'Latitude'//) */
/*  103 FORMAT (8X,I4,3D17.6) */
/*  104 FORMAT (16X,I4,2D17.6) */
/*  105 FORMAT (//1X,'Triangle',8X,'Vertices',12X,'Neighbors'/ */
/*     .        4X,'KT',7X,'N1',5X,'N2',5X,'N3',4X,'KT1',4X, */
/*     .        'KT2',4X,'KT3'/) */
/*  106 FORMAT (//1X,'Triangle',8X,'Vertices',12X,'Neighbors', */
/*     .        14X,'Arcs'/ */
/*     .        4X,'KT',7X,'N1',5X,'N2',5X,'N3',4X,'KT1',4X, */
/*     .        'KT2',4X,'KT3',4X,'KA1',4X,'KA2',4X,'KA3'/) */
/*  107 FORMAT (2X,I4,2X,6(3X,I4),3(2X,I5)) */
/*  108 FORMAT (///) */
/*  109 FORMAT (/1X,'NB = ',I4,' Boundary Nodes',5X, */
/*     .        'NA = ',I5,' Arcs',5X,'NT = ',I5, */
/*     .        ' Triangles') */
/*  110 FORMAT (//1X,10X,'*** Invalid Parameter:  N =',I5, */
/*     .        ', NROW =',I5,', NT =',I5,' ***') */
} /* trlprt_ */

/* Subroutine */ int trmesh_(int *n, double *x, double *y,
	double *z__, int *list, int *lptr, int *lend, int
	*lnew, int *near__, int *next, double *dist, int *ier)
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static double d__;
    static int i__, j, k;
    static double d1, d2, d3;
    static int i0, lp, nn, lpl;
    long int left_(double *, double *, double *, double
	    *, double *, double *, double *, double *,
	    double *);
    static int nexti;


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
/*             int. */

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

/* Subroutine */ int trplot_(int *lun, double *pltsiz, double *
	elat, double *elon, double *a, int *n, double *x,
	double *y, double *z__, int *list, int *lptr, int
	*lend, char *, long int *numbr, int *ier, short )
{
    /* Initialized data */

    static long int annot = TRUE_;
    static double fsizn = 10.;
    static double fsizt = 16.;
    static double tol = .5;

    /* System generated locals */
    int i__1, i__2;
    double d__1;

    /* Builtin functions */
    //double atan(double), sin(double);
    //int i_dnnt(double *);
    //double cos(double), sqrt(double);

    /* Local variables */
    static double t;
    static int n0, n1;
    static double p0[3], p1[3], cf, r11, r12, r21, ct, r22, r23, sf;
    static int ir, lp;
    static double ex, ey, ez, wr, tx, ty;
    static int lpl;
    static double wrs;
    static int ipx1, ipx2, ipy1, ipy2, nseg;
    /* Subroutine */ int drwarc_(int *, double *, double *,
	     double *, int *);


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

/*       LUN = long int unit number in the range 0 to 99. */
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

/* ANNOT =     long int variable with value TRUE iff the plot */
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

/*      WRITE (LUN,100,ERR=13) IPX1, IPY1, IPX2, IPY2 */
/*  100 FORMAT ('%!PS-Adobe-3.0 EPSF-3.0'/ */
/*     .        '%%BoundingBox:',4I4/ */
/*     .        '%%Title:  Triangulation'/ */
/*     .        '%%Creator:  STRIPACK'/ */
/*     .        '%%EndComments') */

/* Set (IPX1,IPY1) and (IPX2,IPY2) to the corner coordinates */
/*   of a viewport box obtained by shrinking the bounding box */
/*   by 12% in each dimension. */

    d__1 = (double) ir * .88;
    ir = i_dnnt(&d__1);
    ipx1 = 306 - ir;
    ipx2 = ir + 306;
    ipy1 = 396 - ir;
    ipy2 = ir + 396;

/* Set the line thickness to 2 points, and draw the */
/*   viewport boundary. */

    t = 2.;
/*      WRITE (LUN,110,ERR=13) T */
/*      WRITE (LUN,120,ERR=13) IR */
/*      WRITE (LUN,130,ERR=13) */
/*  110 FORMAT (F12.6,' setlinewidth') */
/*  120 FORMAT ('306 396 ',I3,' 0 360 arc') */
/*  130 FORMAT ('stroke') */

/* Set up an affine mapping from the window box [-WR,WR] X */
/*   [-WR,WR] to the viewport box. */

    sf = (double) ir / wr;
    tx = ipx1 + sf * wr;
    ty = ipy1 + sf * wr;
/*      WRITE (LUN,140,ERR=13) TX, TY, SF, SF */
/*  140 FORMAT (2F12.6,' translate'/ */
/*    .        2F12.6,' scale') */

/* The line thickness must be changed to reflect the new */
/*   scaling which is applied to all subsequent output. */
/*   Set it to 1.0 point. */

    t = 1. / sf;
/*      WRITE (LUN,110,ERR=13) T */

/* Save the current graphics state, and set the clip path to */
/*   the boundary of the window. */

/*      WRITE (LUN,150,ERR=13) */
/*      WRITE (LUN,160,ERR=13) WR */
/*      WRITE (LUN,170,ERR=13) */
/*  150 FORMAT ('gsave') */
/*  160 FORMAT ('0 0 ',F12.6,' 0 360 arc') */
/*  170 FORMAT ('clip newpath') */

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

/*      WRITE (LUN,130,ERR=13) */
/*      WRITE (LUN,190,ERR=13) */
/*  190 FORMAT ('grestore') */
    if (*numbr) {

/* Nodes in the window are to be labeled with their indexes. */
/*   Convert FSIZN from points to world coordinates, and */
/*   output the commands to select a font and scale it. */

	t = fsizn / sf;
/*        WRITE (LUN,200,ERR=13) T */
/*  200   FORMAT ('/Helvetica findfont'/ */
/*     .          F12.6,' scalefont setfont') */

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

/*          WRITE (LUN,210,ERR=13) P0(1), P0(2) */
/*          WRITE (LUN,220,ERR=13) N0 */
/*  210     FORMAT (2F12.6,' moveto') */
/*  220     FORMAT ('(',I3,') show') */
L4:
	    ;
	}
    }

/* Convert FSIZT from points to world coordinates, and output */
/*   the commands to select a font and scale it. */

    t = fsizt / sf;
/*      WRITE (LUN,200,ERR=13) T */

/* Display TITLE centered above the plot: */

    p0[1] = wr + t * 3.;
/*      WRITE (LUN,230,ERR=13) TITLE, P0(2) */
/*  230 FORMAT (A80/'  stringwidth pop 2 div neg ',F12.6, */
/*     .        ' moveto') */
/*      WRITE (LUN,240,ERR=13) TITLE */
/*  240 FORMAT (A80/'  show') */
    if (annot) {

/* Display the window center and radius below the plot. */

	p0[0] = -wr;
	p0[1] = -wr - 50. / sf;
/*        WRITE (LUN,210,ERR=13) P0(1), P0(2) */
/*        WRITE (LUN,250,ERR=13) ELAT, ELON */
	p0[1] -= t * 2.;
/*        WRITE (LUN,210,ERR=13) P0(1), P0(2) */
/*        WRITE (LUN,260,ERR=13) A */
/*  250   FORMAT ('(Window center:  ELAT = ',F7.2, */
/*     .          ',  ELON = ',F8.2,') show') */
/*  260   FORMAT ('(Angular extent:  A = ',F5.2,') show') */
    }

/* Paint the path and output the showpage command and */
/*   end-of-file indicator. */

/*      WRITE (LUN,270,ERR=13) */
/*  270 FORMAT ('stroke'/ */
/*     .        'showpage'/ */
/*     .        '%%EOF') */

/* HP's interpreters require a one-byte End-of-PostScript-Job */
/*   indicator (to eliminate a timeout error message): */
/*   ASCII 4. */

/*      WRITE (LUN,280,ERR=13) CHAR(4) */
/*  280 FORMAT (A1) */

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

/* L13: */
    *ier = 3;
    return 0;
} /* trplot_ */

/* Subroutine */ int trprnt_(int *n, double *x, double *y,
	double *z__, int *iflag, int *list, int *lptr,
	int *lend, int *lout)
{
    /* Initialized data */

    static int nmax = 9999;
    static int nlmax = 58;

    /* System generated locals */
    int i__1;

    /* Local variables */
    static int k, na, nb, nd, nl, lp, nn, nt, inc, lpl, lun, node, nabor[
	    400];


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
/* coordinates) on long int unit LOUT.  The list of neighbors */
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

/*       LOUT = long int unit for output.  If LOUT is not in */
/*              the range 0 to 99, output is written to */
/*              long int unit 6. */

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
/* LUN =   long int unit for output (copy of LOUT) */
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

/*      WRITE (LUN,100) NN */
    if (nn < 3 || nn > nmax) {

/* N is outside its valid range. */

/*        WRITE (LUN,110) */
	return 0;
    }

/* Initialize NL (the number of lines printed on the current */
/*   page) and NB (the number of boundary nodes encountered). */

    nl = 6;
    nb = 0;
    if (*iflag < 0) {

/* Print LIST only.  K is the number of neighbors of NODE */
/*   that have been stored in NABOR. */

/*        WRITE (LUN,101) */
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
/*            WRITE (LUN,108) */
		nl = inc;
	    }
/*          WRITE (LUN,104) NODE, (NABOR(I), I = 1,K) */
/*          IF (K .NE. 14) */
/* 	     WRITE (LUN,107) */
/* L2: */
	}
    } else if (*iflag > 0) {

/* Print X (longitude), Y (latitude), and LIST. */

/*        WRITE (LUN,102) */
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
/*            WRITE (LUN,108) */
		nl = inc;
	    }
/*          WRITE (LUN,105) NODE, X(NODE), Y(NODE), (NABOR(I), I = 1,K) */
/*          IF (K .NE. 8) */
/* 	     PRINT *,K */
/* 	     WRITE (LUN,107) */
/* L4: */
	}
    } else {

/* Print X, Y, Z, and LIST. */

/*        WRITE (LUN,103) */
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
/*            WRITE (LUN,108) */
		nl = inc;
	    }
/*          WRITE (LUN,106) NODE, X(NODE), Y(NODE),Z(NODE), (NABOR(I), I = 1,K) */
/*          IF (K .NE. 5) */
/* 	     print *,K */
/* 	     WRITE (LUN,107) */
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
/*      WRITE (LUN,109) NB, NA, NT */
    return 0;

/* Print formats: */

/*  100 FORMAT (///15X,'STRIPACK Triangulation Data ', */
/*     .        'Structure,  N = ',I5//) */
/*  101 FORMAT (1X,'Node',31X,'Neighbors of Node'//) */
/*  102 FORMAT (1X,'Node',5X,'Longitude',6X,'Latitude', */
/*     .        18X,'Neighbors of Node'//) */
/*  103 FORMAT (1X,'Node',5X,'X(Node)',8X,'Y(Node)',8X, */
/*     .        'Z(Node)',11X,'Neighbors of Node'//) */
/*  104 FORMAT (1X,I4,4X,14I5/(1X,8X,14I5)) */
/*  105 FORMAT (1X,I4,2D15.6,4X,8I5/(1X,38X,8I5)) */
/*  106 FORMAT (1X,I4,3D15.6,4X,5I5/(1X,53X,5I5)) */
/*  107 FORMAT (1X) */
/*  108 FORMAT (///) */
/*  109 FORMAT (/1X,'NB = ',I4,' Boundary Nodes',5X, */
/*     .        'NA = ',I5,' Arcs',5X,'NT = ',I5, */
/*     .        ' Triangles') */
/*  110 FORMAT (1X,10X,'*** N is outside its valid', */
/*     .        ' range ***') */
} /* trprnt_ */

/* Subroutine */ int vrplot_(int *lun, double *pltsiz, double *
	elat, double *elon, double *a, int *n, double *x,
	double *y, double *z__, int *nt, int *listc, int *
	lptr, int *lend, double *xc, double *yc, double *zc,
	char *, long int *numbr, int *ier, short)
{
    /* Initialized data */

    static long int annot = TRUE_;
    static double fsizn = 10.;
    static double fsizt = 16.;
    static double tol = .5;

    /* System generated locals */
    int i__1;
    double d__1;

    /* Builtin functions */
    //double atan(double), sin(double);
    //int i_dnnt(double *);
    //double cos(double), sqrt(double);

    /* Local variables */
    static double t;
    static int n0;
    static double p1[3], p2[3], x0, y0, cf, r11, r12, r21, ct, r22, r23,
	    sf;
    static int ir, lp;
    static double ex, ey, ez, wr, tx, ty;
    static long int in1, in2;
    static int kv1, kv2, lpl;
    static double wrs;
    static int ipx1, ipx2, ipy1, ipy2, nseg;
    /* Subroutine */ int drwarc_(int *, double *, double *,
	     double *, int *);


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

/*       LUN = long int unit number in the range 0 to 99. */
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

/* ANNOT =     long int variable with value TRUE iff the plot */
/*               is to be annotated with the values of ELAT, */
/*               ELON, and A */
/* CF =        Conversion factor for degrees to radians */
/* CT =        Cos(ELAT) */
/* EX,EY,EZ =  Cartesian coordinates of the eye-position E */
/* FSIZN =     Font size in points for labeling nodes with */
/*               their indexes if NUMBR = TRUE */
/* FSIZT =     Font size in points for the title (and */
/*               annotation if ANNOT = TRUE) */
/* IN1,IN2 =   long int variables with value TRUE iff the */
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

/*      WRITE (LUN,100,ERR=13) IPX1, IPY1, IPX2, IPY2 */
/*  100 FORMAT ('%!PS-Adobe-3.0 EPSF-3.0'/ */
/*     .        '%%BoundingBox:',4I4/ */
/*     .        '%%Title:  Voronoi diagram'/ */
/*     .        '%%Creator:  STRIPACK'/ */
/*     .        '%%EndComments') */
/* Set (IPX1,IPY1) and (IPX2,IPY2) to the corner coordinates */
/*   of a viewport box obtained by shrinking the bounding box */
/*   by 12% in each dimension. */

    d__1 = (double) ir * .88;
    ir = i_dnnt(&d__1);
    ipx1 = 306 - ir;
    ipx2 = ir + 306;
    ipy1 = 396 - ir;
    ipy2 = ir + 396;

/* Set the line thickness to 2 points, and draw the */
/*   viewport boundary. */

    t = 2.;
/*      WRITE (LUN,110,ERR=13) T */
/*      WRITE (LUN,120,ERR=13) IR */
/*      WRITE (LUN,130,ERR=13) */
/*  110 FORMAT (F12.6,' setlinewidth') */
/*  120 FORMAT ('306 396 ',I3,' 0 360 arc') */
/*  130 FORMAT ('stroke') */

/* Set up an affine mapping from the window box [-WR,WR] X */
/*   [-WR,WR] to the viewport box. */

    sf = (double) ir / wr;
    tx = ipx1 + sf * wr;
    ty = ipy1 + sf * wr;
/*      WRITE (LUN,140,ERR=13) TX, TY, SF, SF */
/*  140 FORMAT (2F12.6,' translate'/ */
/*     .        2F12.6,' scale') */

/* The line thickness must be changed to reflect the new */
/*   scaling which is applied to all subsequent output. */
/*   Set it to 1.0 point. */

    t = 1. / sf;
/*      WRITE (LUN,110,ERR=13) T */

/* Save the current graphics state, and set the clip path to */
/*   the boundary of the window. */

/*      WRITE (LUN,150,ERR=13) */
/*      WRITE (LUN,160,ERR=13) WR */
/*      WRITE (LUN,170,ERR=13) */
/*  150 FORMAT ('gsave') */
/*  160 FORMAT ('0 0 ',F12.6,' 0 360 arc') */
/*  170 FORMAT ('clip newpath') */

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

	if (! in1 || (in2 && kv2 <= kv1)) {
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
	    p2[2] = 0.f;
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

/*      WRITE (LUN,130,ERR=13) */
/*      WRITE (LUN,190,ERR=13) */
/*  190 FORMAT ('grestore') */
    if (*numbr) {

/* Nodes in the window are to be labeled with their indexes. */
/*   Convert FSIZN from points to world coordinates, and */
/*   output the commands to select a font and scale it. */

	t = fsizn / sf;
/*        WRITE (LUN,200,ERR=13) T */
/*  200   FORMAT ('/Helvetica findfont'/ */
/*     .          F12.6,' scalefont setfont') */

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

/*          WRITE (LUN,210,ERR=13) X0, Y0 */
/*          WRITE (LUN,220,ERR=13) N0 */
/*  210     FORMAT (2F12.6,' moveto') */
/*  220     FORMAT ('(',I3,') show') */
L4:
	    ;
	}
    }

/* Convert FSIZT from points to world coordinates, and output */
/*   the commands to select a font and scale it. */

    t = fsizt / sf;
/*      WRITE (LUN,200,ERR=13) T */

/* Display TITLE centered above the plot: */

    y0 = wr + t * 3.;
/*      WRITE (LUN,230,ERR=13) TITLE, Y0 */
/*  230 FORMAT (A80/'  stringwidth pop 2 div neg ',F12.6, */
/*     .        ' moveto') */
/*      WRITE (LUN,240,ERR=13) TITLE */
/*  240 FORMAT (A80/'  show') */
    if (annot) {

/* Display the window center and radius below the plot. */

	x0 = -wr;
	y0 = -wr - 50. / sf;
/*        WRITE (LUN,210,ERR=13) X0, Y0 */
/*        WRITE (LUN,250,ERR=13) ELAT, ELON */
	y0 -= t * 2.;
/*        WRITE (LUN,210,ERR=13) X0, Y0 */
/*        WRITE (LUN,260,ERR=13) A */
/*  250   FORMAT ('(Window center:  ELAT = ',F7.2, */
/*     .          ',  ELON = ',F8.2,') show') */
/*  260   FORMAT ('(Angular extent:  A = ',F5.2,') show') */
    }

/* Paint the path and output the showpage command and */
/*   end-of-file indicator. */

/*      WRITE (LUN,270,ERR=13) */
/*  270 FORMAT ('stroke'/ */
/*     .        'showpage'/ */
/*     .        '%%EOF') */

/* HP's interpreters require a one-byte End-of-PostScript-Job */
/*   indicator (to eliminate a timeout error message): */
/*   ASCII 4. */

/*      WRITE (LUN,280,ERR=13) CHAR(4) */
/*  280 FORMAT (A1) */

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

/* L13: */
    *ier = 3;
    return 0;
} /* vrplot_ */

/* Subroutine */ int random_(int *ix, int *iy, int *iz,
	double *rannum)
{
    static double x;


/*   This routine returns pseudo-random numbers uniformly */
/* distributed in the interval (0,1).  int seeds IX, IY, */
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
    x = (double) (*ix) / 30269. + (double) (*iy) / 30307. + (
	    double) (*iz) / 30323.;
    *rannum = x - (int) x;
    return 0;
} /* random_ */

#undef TRUE_
#undef FALSE_
#undef abs

/*################################################################################################
##########  strid.f -- translated by f2c (version 20030320). ###################################
######   You must link the resulting object file with the libraries: #############################
####################	-lf2c -lm   (in that order)   ############################################
################################################################################################*/


EMData* Util::mult_scalar(EMData* img, float scalar)
{
	ENTERFUNC;
	/* Exception Handle */
	if (!img) {
		throw NullPointerException("NULL input image");
	}
	/* ============  output = scalar*input  ================== */

	int nx=img->get_xsize(),ny=img->get_ysize(),nz=img->get_zsize();
	size_t size = (size_t)nx*ny*nz;
	EMData * img2   = img->copy_head();
	float *img_ptr  = img->get_data();
	float *img2_ptr = img2->get_data();
	for (size_t i=0;i<size;++i)img2_ptr[i] = img_ptr[i]*scalar;
	img2->update();

	if(img->is_complex()) {
		img2->set_complex(true);
		if(img->is_fftodd()) img2->set_fftodd(true); else img2->set_fftodd(false);
	}
	EXITFUNC;
	return img2;
}

EMData* Util::madn_scalar(EMData* img, EMData* img1, float scalar)
{
	ENTERFUNC;
	/* Exception Handle */
	if (!img) {
		throw NullPointerException("NULL input image");
	}
	/* ==============   output = img + scalar*img1   ================ */

	int nx=img->get_xsize(),ny=img->get_ysize(),nz=img->get_zsize();
	size_t size = (size_t)nx*ny*nz;
	EMData * img2   = img->copy_head();
	float *img_ptr  = img->get_data();
	float *img2_ptr = img2->get_data();
	float *img1_ptr = img1->get_data();
	for (size_t i=0;i<size;++i) img2_ptr[i] = img_ptr[i] + img1_ptr[i]*scalar;
	img2->update();
	if(img->is_complex()) {
		img2->set_complex(true);
		if(img->is_fftodd()) img2->set_fftodd(true); else img2->set_fftodd(false);
	}

	EXITFUNC;
	return img2;
}

EMData* Util::addn_img(EMData* img, EMData* img1)
{
	ENTERFUNC;
	/* Exception Handle */
	if (!img) {
		throw NullPointerException("NULL input image");
	}
	/* ==============   output = img + img1   ================ */

	int nx=img->get_xsize(),ny=img->get_ysize(),nz=img->get_zsize();
	size_t size = (size_t)nx*ny*nz;
	EMData * img2   = img->copy_head();
	float *img_ptr  = img->get_data();
	float *img2_ptr = img2->get_data();
	float *img1_ptr = img1->get_data();
	for (size_t i=0;i<size;++i) img2_ptr[i] = img_ptr[i] + img1_ptr[i];
	img2->update();
	if(img->is_complex()) {
		img2->set_complex(true);
		if(img->is_fftodd()) img2->set_fftodd(true); else img2->set_fftodd(false);
	}

	EXITFUNC;
	return img2;
}

EMData* Util::subn_img(EMData* img, EMData* img1)
{
	ENTERFUNC;
	/* Exception Handle */
	if (!img) {
		throw NullPointerException("NULL input image");
	}
	/* ==============   output = img - img1   ================ */

	int nx=img->get_xsize(),ny=img->get_ysize(),nz=img->get_zsize();
	size_t size = (size_t)nx*ny*nz;
	EMData * img2   = img->copy_head();
	float *img_ptr  = img->get_data();
	float *img2_ptr = img2->get_data();
	float *img1_ptr = img1->get_data();
	for (size_t i=0;i<size;++i) img2_ptr[i] = img_ptr[i] - img1_ptr[i];
	img2->update();
	if(img->is_complex()) {
		img2->set_complex(true);
		if(img->is_fftodd()) img2->set_fftodd(true); else img2->set_fftodd(false);
	}

	EXITFUNC;
	return img2;
}

EMData* Util::muln_img(EMData* img, EMData* img1)
{
	ENTERFUNC;
	/* Exception Handle */
	if (!img) {
		throw NullPointerException("NULL input image");
	}
	/* ==============   output = img * img1   ================ */
	/* =====complex =   output = img * conjg(img1)   ================ */
	int nx=img->get_xsize(),ny=img->get_ysize(),nz=img->get_zsize();
	size_t size = (size_t)nx*ny*nz;
	EMData * img2   = img->copy_head();
	float *img_ptr  = img->get_data();
	float *img2_ptr = img2->get_data();
	float *img1_ptr = img1->get_data();
	if(img->is_complex()) {
		for (size_t i=0; i<size; i+=2) {
			img2_ptr[i]   = img_ptr[i] * img1_ptr[i]   - img_ptr[i+1] * img1_ptr[i+1] ;
			img2_ptr[i+1] = img_ptr[i] * img1_ptr[i+1] + img_ptr[i+1] * img1_ptr[i] ;
		}
		img2->set_complex(true);
		if(img->is_fftodd()) img2->set_fftodd(true); else img2->set_fftodd(false);
	} else {
		for (size_t i=0; i<size; ++i) img2_ptr[i] = img_ptr[i] * img1_ptr[i];
		img2->update();
	}

	EXITFUNC;
	return img2;
}

EMData* Util::divn_img(EMData* img, EMData* img1)
{
	ENTERFUNC;
	/* Exception Handle */
	if (!img) {
		throw NullPointerException("NULL input image");
	}
	/* ==============   output = img / img1   ================ */

	int nx=img->get_xsize(),ny=img->get_ysize(),nz=img->get_zsize();
	size_t size = (size_t)nx*ny*nz;
	EMData * img2   = img->copy_head();
	float *img_ptr  = img->get_data();
	float *img2_ptr = img2->get_data();
	float *img1_ptr = img1->get_data();
	if(img->is_complex()) {
		float  sq2;
		for (size_t i=0; i<size; i+=2) {
			sq2 = 1.0f/(img1_ptr[i] * img1_ptr[i]   + img1_ptr[i+1] * img1_ptr[i+1]);
			img2_ptr[i]   = sq2*(img_ptr[i] * img1_ptr[i]   + img_ptr[i+1] * img1_ptr[i+1]) ;
			img2_ptr[i+1] = sq2*(img_ptr[i+1] * img1_ptr[i] - img_ptr[i] * img1_ptr[i+1]) ;
		}
		img2->set_complex(true);
		if(img->is_fftodd()) img2->set_fftodd(true); else img2->set_fftodd(false);
	} else {
		for (size_t i=0; i<size; ++i) img2_ptr[i] = img_ptr[i] / img1_ptr[i];
		img2->update();
	}

	EXITFUNC;
	return img2;
}

EMData* Util::squaren_img(EMData* img)
{
	ENTERFUNC;
	/* Exception Handle */
	if (!img) {
		throw NullPointerException("NULL input image");
	}
	/* ========= img = |img|^2 ===================== */

	int nx=img->get_xsize(),ny=img->get_ysize(),nz=img->get_zsize();
	size_t size = (size_t)nx*ny*nz;
	EMData * img2   = img->copy_head();
	float *img_ptr  = img->get_data();
	float *img2_ptr = img2->get_data();
	if(img->is_complex()) {
		for (size_t i=0; i<size; i+=2) {
			img2_ptr[i]    = img_ptr[i] * img_ptr[i] + img_ptr[i+1] * img_ptr[i+1];;
			img2_ptr[i+1]  = 0.0f;
		}
		img2->set_complex(true);
		if(img->is_fftodd()) img2->set_fftodd(true); else img2->set_fftodd(false);
	} else {
		for (size_t i=0;i<size;++i) img_ptr[i] *= img_ptr[i];
	}
	img2->update();

	EXITFUNC;
	return img2;
}

EMData* Util::divn_filter(EMData* img, EMData* img1)
{
	ENTERFUNC;
	/* Exception Handle */
	if (!img) {
		throw NullPointerException("NULL input image");
	}
	/* ========= img /= img1 ===================== */

	int nx=img->get_xsize(),ny=img->get_ysize(),nz=img->get_zsize();
	size_t size = (size_t)nx*ny*nz;
	EMData * img2   = img->copy_head();
	float *img_ptr  = img->get_data();
	float *img1_ptr = img1->get_data();
	float *img2_ptr = img2->get_data();
	if(img->is_complex()) {
		for (size_t i=0; i<size; i+=2) {
			if(img1_ptr[i] > 1.e-10f) {
				img2_ptr[i]   = img_ptr[i]  /img1_ptr[i];
				img2_ptr[i+1] = img_ptr[i+1]/img1_ptr[i];
			} else img2_ptr[i] = img2_ptr[i+1] = 0.0f;
		}
	} else  throw ImageFormatException("Only Fourier image allowed");

	img2->update();

	EXITFUNC;
	return img2;
}

void Util::mul_scalar(EMData* img, float scalar)
{
	ENTERFUNC;
	/* Exception Handle */
	if (!img) {
		throw NullPointerException("NULL input image");
	}
	/* ============  output = scalar*input  ================== */

	int nx=img->get_xsize(),ny=img->get_ysize(),nz=img->get_zsize();
	size_t size = (size_t)nx*ny*nz;
	float *img_ptr  =img->get_data();
	for (size_t i=0;i<size;++i) img_ptr[i] *= scalar;
	img->update();

	EXITFUNC;
}

void Util::mad_scalar(EMData* img, EMData* img1, float scalar)
{
	ENTERFUNC;
	/* Exception Handle */
	if (!img) {
		throw NullPointerException("NULL input image");
	}
	/* ==============   img += scalar*img1   ================ */

	int nx=img->get_xsize(),ny=img->get_ysize(),nz=img->get_zsize();
	size_t size = (size_t)nx*ny*nz;
	float *img_ptr  = img->get_data();
	float *img1_ptr = img1->get_data();
	for (size_t i=0;i<size;++i)img_ptr[i] += img1_ptr[i]*scalar;
	img1->update();

	EXITFUNC;
}

void Util::add_img(EMData* img, EMData* img1)
{
	ENTERFUNC;
	/* Exception Handle */
	if (!img || !img1) {
		throw NullPointerException("NULL input image");
	}
	/* ========= img += img1 ===================== */

	int nx=img->get_xsize(),ny=img->get_ysize(),nz=img->get_zsize();
	size_t size = (size_t)nx*ny*nz;
	float *img_ptr  = img->get_data();
	float *img1_ptr = img1->get_data();
	for (size_t i=0;i<size;++i) img_ptr[i] += img1_ptr[i];
	img->update();

	EXITFUNC;
}

void Util::add_img_abs(EMData* img, EMData* img1)
{
	ENTERFUNC;
	/* Exception Handle */
	if (!img) {
		throw NullPointerException("NULL input image");
	}
	/* ========= img += img1 ===================== */

	int nx=img->get_xsize(),ny=img->get_ysize(),nz=img->get_zsize();
	size_t size = (size_t)nx*ny*nz;
	float *img_ptr  = img->get_data();
	float *img1_ptr = img1->get_data();
	for (size_t i=0;i<size;++i) img_ptr[i] += abs(img1_ptr[i]);
	img->update();

	EXITFUNC;
}

void Util::add_img2(EMData* img, EMData* img1)
{
	ENTERFUNC;
	/* Exception Handle */
	if (!img) {
		throw NullPointerException("NULL input image");
	}
	/* ========= img += img1**2 ===================== */

	int nx=img->get_xsize(),ny=img->get_ysize(),nz=img->get_zsize();
	size_t size = (size_t)nx*ny*nz;
	float *img_ptr  = img->get_data();
	float *img1_ptr = img1->get_data();
	if(img->is_complex()) {
		for (size_t i=0; i<size; i+=2) img_ptr[i] += img1_ptr[i] * img1_ptr[i] + img1_ptr[i+1] * img1_ptr[i+1];
	} else {
		for (size_t i=0;i<size;++i) img_ptr[i] += img1_ptr[i]*img1_ptr[i];
	}
	img->update();

	EXITFUNC;
}

void Util::sub_img(EMData* img, EMData* img1)
{
	ENTERFUNC;
	/* Exception Handle */
	if (!img) {
		throw NullPointerException("NULL input image");
	}
	/* ========= img -= img1 ===================== */

	int nx=img->get_xsize(),ny=img->get_ysize(),nz=img->get_zsize();
	size_t size = (size_t)nx*ny*nz;
	float *img_ptr  = img->get_data();
	float *img1_ptr = img1->get_data();
	for (size_t i=0;i<size;++i) img_ptr[i] -= img1_ptr[i];
	img->update();

	EXITFUNC;
}

void Util::mul_img(EMData* img, EMData* img1)
{
	ENTERFUNC;
	/* Exception Handle */
	if (!img) {
		throw NullPointerException("NULL input image");
	}
	/* ========= img *= img1 ===================== */

	int nx=img->get_xsize(),ny=img->get_ysize(),nz=img->get_zsize();
	size_t size = (size_t)nx*ny*nz;
	float *img_ptr  = img->get_data();
	float *img1_ptr = img1->get_data();
	if(img->is_complex()) {
		for (size_t i=0; i<size; i+=2) {
			float tmp     = img_ptr[i] * img1_ptr[i]   - img_ptr[i+1] * img1_ptr[i+1];
			img_ptr[i+1]  = img_ptr[i] * img1_ptr[i+1] + img_ptr[i+1] * img1_ptr[i];
			img_ptr[i]    = tmp;

		}
	} else {
		for (size_t i=0;i<size;++i) img_ptr[i] *= img1_ptr[i];
	}
	img->update();

	EXITFUNC;
}

void Util::div_img(EMData* img, EMData* img1)
{
	ENTERFUNC;
	/* Exception Handle */
	if (!img) {
		throw NullPointerException("NULL input image");
	}
	/* ========= img /= img1 ===================== */

	int nx=img->get_xsize(),ny=img->get_ysize(),nz=img->get_zsize();
	size_t size = (size_t)nx*ny*nz;
	float *img_ptr  = img->get_data();
	float *img1_ptr = img1->get_data();
	if(img->is_complex()) {
		float  sq2;
		for (size_t i=0; i<size; i+=2) {
			sq2 = 1.0f/(img1_ptr[i] * img1_ptr[i]   + img1_ptr[i+1] * img1_ptr[i+1]);
			float tmp    = sq2*(img_ptr[i] * img1_ptr[i]   + img_ptr[i+1] * img1_ptr[i+1]);
			img_ptr[i+1] = sq2*(img_ptr[i+1] * img1_ptr[i] - img_ptr[i] * img1_ptr[i+1]);
			img_ptr[i]   = tmp;
		}
	} else {
		for (size_t i=0; i<size; ++i) img_ptr[i] /= img1_ptr[i];
	}
	img->update();

	EXITFUNC;
}


void Util::square_img(EMData* img)
{
	ENTERFUNC;
	/* Exception Handle */
	if (!img) {
		throw NullPointerException("NULL input image");
	}
	/* ========= img = |img|^2 ===================== */

	int nx=img->get_xsize(),ny=img->get_ysize(),nz=img->get_zsize();
	size_t size = (size_t)nx*ny*nz;
	float *img_ptr  = img->get_data();
	if(img->is_complex()) {
		for (size_t i=0; i<size; i+=2) {
			img_ptr[i]    = img_ptr[i] * img_ptr[i] + img_ptr[i+1] * img_ptr[i+1];;
			img_ptr[i+1]  = 0.0f;
		}
	} else {
		for (size_t i=0;i<size;++i) img_ptr[i] *= img_ptr[i];
	}
	img->update();

	EXITFUNC;
}


void Util::div_filter(EMData* img, EMData* img1)
{
	ENTERFUNC;
	/* Exception Handle */
	if (!img) {
		throw NullPointerException("NULL input image");
	}
	/* ========= img /= img1 ===================== */

	int nx=img->get_xsize(),ny=img->get_ysize(),nz=img->get_zsize();
	size_t size = (size_t)nx*ny*nz;
	float *img_ptr  = img->get_data();
	float *img1_ptr = img1->get_data();
	if(img->is_complex()) {
		for (size_t i=0; i<size; i+=2) {
			if(img1_ptr[i] > 1.e-10f) {
				img_ptr[i]   /= img1_ptr[i];
				img_ptr[i+1] /= img1_ptr[i];
			} else img_ptr[i] = img_ptr[i+1] = 0.0f;
		}
	} else throw ImageFormatException("Only Fourier image allowed");

	img->update();

	EXITFUNC;
}


void Util::set_freq(EMData* freqvol, EMData* temp, EMData* mask, float cutoff, float freq)
{
	ENTERFUNC;
	/* Exception Handle */
	if (!freqvol) {
		throw NullPointerException("NULL input image");
	}

	int nx=freqvol->get_xsize(),ny=freqvol->get_ysize(),nz=freqvol->get_zsize();
	size_t size = (size_t)nx*ny*nz;
	float *freqvol_ptr  = freqvol->get_data();
	float *temp_ptr = temp->get_data();
	float *mask_ptr = mask->get_data();

	for (size_t i=0;i<size;++i) {
		if(mask_ptr[i] >0.5f) {
			if(freqvol_ptr[i]  == 0.0f) {
				if(temp_ptr[i] < cutoff) freqvol_ptr[i] = freq;
			}
		}
	}

	freqvol->update();
	EXITFUNC;
}



#define img_ptr(i,j,k)  img_ptr[2*(i-1)+((j-1)+((k-1)*ny))*(size_t)nxo]

EMData* Util::pack_complex_to_real(EMData* img)
{
	ENTERFUNC;
	/* Exception Handle */
	if (!img) {
		throw NullPointerException("NULL input image");
	}
	/* ==============   img is modulus of a complex image in FFT format (so its imaginary parts are zero),
	                      output is img packed into real image with Friedel part added,   ================ */

	int nxo=img->get_xsize(), ny=img->get_ysize(), nz=img->get_zsize();
	int nx = nxo - 2 + img->is_fftodd();
	int lsd2 = (nx + 2 - nx%2) / 2; // Extended x-dimension of the complex image
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
	//img->set_array_offsets(1,1,1);
	float *img_ptr  = img->get_data();
	for (int iz = 1; iz <= nz; iz++) {
		int jz=iz-1;
		if(jz>=nz2p) jz=jz-nzt;
		for (int iy = 1; iy <= ny; iy++) {
			int jy=iy-1;
			if(jy>=ny2p) jy=jy-nyt;
			for (int ix = 1; ix <= lsd2; ix++) {
				int jx=ix-1;
				if(jx>=nx2p) jx=jx-nx;
				power(jx,jy,jz) = img_ptr(ix,iy,iz); //real(img->cmplx(ix,iy,iz));
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
	power.update();
	power.set_array_offsets(0,0,0);
	return &power;
}
#undef  img_ptr

float Util::ang_n(float peakp, string mode, int maxrin)
{
    if (mode == "f" || mode == "F")
        return fmodf(((peakp-1.0f) / maxrin+1.0f)*360.0f,360.0f);
    else
        return fmodf(((peakp-1.0f) / maxrin+1.0f)*180.0f,180.0f);
}


void Util::Normalize_ring( EMData* ring, const vector<int>& numr )
{
    float* data = ring->get_data();
    float av=0.0;
    float sq=0.0;
    float nn=0.0;
    int nring = numr.size()/3;
    for( int i=0; i < nring; ++i )
    {
        int numr3i = numr[3*i+2];
        int numr2i = numr[3*i+1]-1;
        float w = numr[3*i]*2*M_PI/float(numr[3*i+2]);
        for( int j=0; j < numr3i; ++j )
        {
            int jc = numr2i+j;
            av += data[jc] * w;
            sq += data[jc] * data[jc] * w;
            nn += w;
        }
    }

    float avg = av/nn;
    float sgm = sqrt( (sq-av*av/nn)/nn );
    size_t n = (size_t)ring->get_xsize() * ring->get_ysize() * ring->get_zsize();
    for( size_t i=0; i < n; ++i )
    {
        data[i] -= avg;
        data[i] /= sgm;
    }

    ring->update();
}

vector<float> Util::multiref_polar_ali_2d(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, string mode,
                vector<int>numr, float cnx, float cny) {

    // Manually extract.
/*    vector< EMAN::EMData* > crefim;
    std::size_t crefim_len = PyObject_Length(crefim_list.ptr());
    crefim.reserve(crefim_len);

    for(std::size_t i=0;i<crefim_len;i++) {
        boost::python::extract<EMAN::EMData*> proxy(crefim_list[i]);
        crefim.push_back(proxy());
    }
*/
	const float qv = static_cast<float>( pi/180.0 );
	size_t crefim_len = crefim.size();

// 	int   ky = int(2*yrng/step+0.5)/2;
// 	int   kx = int(2*xrng/step+0.5)/2;
	int lkx = int(xrng[0]/step);
	int rkx = int(xrng[1]/step);
	int lky = int(yrng[0]/step);
	int rky = int(yrng[1]/step);
	
	int   iref, nref=0, mirror=0;
	float iy, ix, sx=0, sy=0;
	float peak = -1.0E23f;
	float ang=0.0f;
	for (int i = -lky; i <= rky; i++) {
		iy = i * step ;
		for (int j = -lkx; j <= rkx; j++) {
			ix = j*step ;
			EMData* cimage = Polar2Dm(image, cnx+ix, cny+iy, numr, mode);

			Normalize_ring( cimage, numr );

			Frngs(cimage, numr);
			//  compare with all reference images
			// for iref in xrange(len(crefim)):
			for ( iref = 0; iref < (int)crefim_len; iref++) {
				Dict retvals = Crosrng_ms(crefim[iref], cimage, numr);
				double qn = retvals["qn"];
				double qm = retvals["qm"];
				if(qn >= peak || qm >= peak) {
					sx = -ix;
					sy = -iy;
					nref = iref;
					if (qn >= qm) {
						ang = ang_n(retvals["tot"], mode, numr[numr.size()-1]);
						peak = static_cast<float>(qn);
						mirror = 0;
					} else {
						ang = ang_n(retvals["tmt"], mode, numr[numr.size()-1]);
						peak = static_cast<float>(qm);
						mirror = 1;
					}
				}
			}  delete cimage; cimage = 0;
		}	
	}

	float co, so, sxs, sys;
	co = static_cast<float>(  cos(ang*qv) );
	so = static_cast<float>( -sin(ang*qv) );
	sxs = sx*co - sy*so;
	sys = sx*so + sy*co;
	vector<float> res;
	res.push_back(ang);
	res.push_back(sxs);
	res.push_back(sys);
	res.push_back(static_cast<float>(mirror));
	res.push_back(static_cast<float>(nref));
	res.push_back(peak);
	return res;
}

vector<float> Util::multiref_polar_ali_2d_peaklist(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, string mode,
                vector<int>numr, float cnx, float cny) {

	const float qv = static_cast<float>( pi/180.0 );
	size_t crefim_len = crefim.size();

//	int   ky = int(2*yrng/step+0.5)/2;
//	int   kx = int(2*xrng/step+0.5)/2;
	int lkx = int(xrng[0]/step);
	int rkx = int(xrng[1]/step);
	int lky = int(yrng[0]/step);
	int rky = int(yrng[1]/step);
		
	float iy, ix;
	vector<float> peak(crefim_len*5);  //, -1.0e23f);
	for (int iref = 0; iref < (int)crefim_len; iref++) {
		peak[iref*5] = -1.0e23f;
		peak[iref*5+1] = 0.0f;
		peak[iref*5+2] = 0.0;
		peak[iref*5+3] = 0.0;
		peak[iref*5+4] = 0;
	}
	
	for (int i = -lky; i <= rky; i++) {
		iy = i * step ;
		for (int j = -lkx; j <= rkx; j++) {
			ix = j*step ;
			EMData* cimage = Polar2Dm(image, cnx+ix, cny+iy, numr, mode);
			Normalize_ring( cimage, numr );
			Frngs(cimage, numr);
			for (int iref = 0; iref < (int)crefim_len; iref++) {
				Dict retvals = Crosrng_ms(crefim[iref], cimage, numr);
				double qn = retvals["qn"];
				double qm = retvals["qm"];
				if(qn >= peak[iref*5] || qm >= peak[iref*5]) {
					if (qn >= qm) {
						peak[iref*5] = static_cast<float>(qn);
						peak[iref*5+1] = ang_n(retvals["tot"], mode, numr[numr.size()-1]);
						peak[iref*5+2] = -ix;
						peak[iref*5+3] = -iy;
						peak[iref*5+4] = 0;
					} else {
						peak[iref*5] = static_cast<float>(qm);
						peak[iref*5+1] = ang_n(retvals["tmt"], mode, numr[numr.size()-1]);
						peak[iref*5+2] = -ix;
						peak[iref*5+3] = -iy;
						peak[iref*5+4] = 1;
					}
				}
			}  delete cimage; cimage = 0;
		}
	}

	for (int iref = 0; iref < (int)crefim_len; iref++) {
		float ang = peak[iref*5+1];
		float sx  = peak[iref*5+2];
		float sy  = peak[iref*5+3];
		float co  =  cos(ang*qv);
		float so  = -sin(ang*qv);
		float sxs = sx*co - sy*so;
		float sys = sx*so + sy*co;
		peak[iref*5+2] = sxs;
		peak[iref*5+3] = sys;
	}
	return peak;
}

vector<float> Util::multiref_polar_ali_2d_peaklist_local(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, float ant, string mode,
                vector<int>numr, float cnx, float cny) {

	size_t crefim_len = crefim.size();
	const float qv = static_cast<float>( pi/180.0 );

	Transform * t = image->get_attr("xform.projection");
	Dict d = t->get_params("spider");
	if(t) {delete t; t=0;}
	float phi   = d["phi"];
	float theta = d["theta"];
//	int   ky    = int(2*yrng/step+0.5)/2;
//	int   kx    = int(2*xrng/step+0.5)/2;
	int lkx = int(xrng[0]/step);
	int rkx = int(xrng[1]/step);
	int lky = int(yrng[0]/step);
	int rky = int(yrng[1]/step);
		
	int   iref, nref=0, mirror=0;
	float imn1 = sin(theta*qv)*cos(phi*qv);
	float imn2 = sin(theta*qv)*sin(phi*qv);
	float imn3 = cos(theta*qv);
	vector<float> n1(crefim_len);
	vector<float> n2(crefim_len);
	vector<float> n3(crefim_len);
	for ( iref = 0; iref < (int)crefim_len; iref++) {
		n1[iref] = crefim[iref]->get_attr("n1");
		n2[iref] = crefim[iref]->get_attr("n2");
		n3[iref] = crefim[iref]->get_attr("n3");
	}

	vector<float> peak;
	for (int i = -lky; i <= rky; i++) {
		float iy = i * step ;
		for (int j = -lkx; j <= rkx; j++) {
			float ix = j*step ;
			EMData* cimage = Polar2Dm(image, cnx+ix, cny+iy, numr, mode);
			Normalize_ring( cimage, numr );
			Frngs(cimage, numr);
			for (int iref = 0; iref < (int)crefim_len; iref++) {
// UUU
//cout<<"  "<<n1[iref]*imn1 + n2[iref]*imn2 + n3[iref]*imn3<<"  "<<ant<<endl;
				if( n1[iref]*imn1 + n2[iref]*imn2 + n3[iref]*imn3 >=ant ) {
					Dict retvals = Crosrng_e(crefim[iref], cimage, numr, 0);
					const float ang = ang_n(retvals["tot"], mode, numr[numr.size()-1]);
					const float sx  = -ix;
					const float sy  = -iy;
					const float co  =  cos(ang*qv);
					const float so  = -sin(ang*qv);
					const float sxs = sx*co - sy*so;
					const float sys = sx*so + sy*co;
					double qn = retvals["qn"];
					peak.push_back(static_cast<float>(qn));
					peak.push_back( ang );
					peak.push_back( sxs );
					peak.push_back( sys );
					peak.push_back( static_cast<float>(iref) );
				}	
			}  delete cimage; cimage = 0;
		}
	}
	return peak;
}

struct assign_groups_comparator {
	const float * values;
	bool operator() (int i,int j) { return (values[i] > values[j]); }
	assign_groups_comparator(const float * v) : values(v) {}
};

vector<int> Util::assign_groups(std::string matrix_address, int nref, int nima)
{
	const float * matrix;
	{ // convert memory address sent as string to pointer to float
		size_t addr = 0;
		for ( std::string::const_iterator i = matrix_address.begin();  i != matrix_address.end();  ++i ) {
			int digit = *i - '0';
			addr *= 10;
			addr += digit;
		}
		matrix = reinterpret_cast<float*>(addr);
	}

	int kt = nref;
	unsigned int maxasi = nima/nref;
	vector< vector<int> > id_list(nref);
	int group, ima;

	{
		int begin = 0;
		// allocate and sort vector of indexes
		std::vector<int> dd(nref*nima);
		for (int i=0; i<nref*nima; i++)  {
			dd[i] = i;
		}
		assign_groups_comparator comparator(matrix);
		sort(dd.begin(), dd.end(), comparator);
		// main loop
		std::vector<bool> del_row(nref, false);
		std::vector<bool> del_column(nima, false);
		while (kt > 0) {
			bool flag = true;
			while (flag) {
				int l = dd[begin];
				group = l/nima;
				ima = l%nima;
				if (del_column[ima] || del_row[group]) begin++;
				else flag = false;
			}

			id_list[group].push_back(ima);
			if (kt > 1) {
				if (id_list[group].size() < maxasi) group = -1;
				else kt -= 1;
			} else {
				if (id_list[group].size() < maxasi+nima%nref) group = -1;
				else kt -= 1;
			}
			del_column[ima] = true;
			if (group != -1) {
				del_row[group] = true;
			}
		}
	}

	vector<int> id_list_1; 
	for (int iref=0; iref<nref; iref++)
		for (unsigned int im=0; im<maxasi; im++)
			id_list_1.push_back(id_list[iref][im]);
	for (unsigned int im=maxasi; im<maxasi+nima%nref; im++)
			id_list_1.push_back(id_list[group][im]);
	id_list_1.push_back(group);

	return id_list_1;
}

int Util::nearest_ang(const vector<float>& vecref, float x, float y, float z) {
	float best_v = -1.0f;
	int best_i = -1;

	for (unsigned int i=0; i<vecref.size()/3; i++) {
		float v = abs(vecref[i*3]*x+vecref[i*3+1]*y+vecref[i*3+2]*z);
		if (v > best_v) {
			best_v = v;
			best_i = i;
		}
	}
	return best_i;
}

struct d_ang {
	float d;
	int i;
	int mirror;
	d_ang(float _d, int _i, int _m):d(_d), i(_i), mirror(_m) {}
	bool operator<(const d_ang& a) const { return d < a.d || (d == a.d && i < a.i); }
};

vector<int> Util::assign_projangles(const vector<float>& projangles, const vector<float>& refangles) {
	int nref = refangles.size()/2;
	int nproj = projangles.size()/2;
	vector<int> asg(nproj);
	vector<float> vecref(nref*3);
	for (int i=0; i<nref; i++)
		getvec(refangles[i*2], refangles[i*2+1], vecref[i*3], vecref[i*3+1], vecref[i*3+2]);
	for (int i=0; i<nproj; i++) {
		float x, y, z;
		getvec(projangles[i*2], projangles[i*2+1], x, y, z);
		asg[i] = nearest_ang(vecref, x, y, z);
	}
	return asg;
}


vector<int> Util::nearestk_to_refdir(const vector<float>& projangles, const vector<float>& refangles, const int howmany) {
	int nref = refangles.size()/2;
	int nproj = projangles.size()/2;
	vector<int> asg(howmany*nref);
	vector<float> vecproj(nproj*3);
	for (int i=0; i<nproj; i++)
		getvec(projangles[i*2], projangles[i*2+1], vecproj[i*3], vecproj[i*3+1], vecproj[i*3+2]);


	vector<bool> taken(nproj);
	for (int k=0; k<nref; k++) {
		for (int i=0; i<nproj; i++)  taken[i] = true;
		float x, y, z;
		getvec(refangles[k*2], refangles[k*2+1], x, y, z);
		for (int h=0; h<howmany; h++) {
			float best_v = -1.0f;
			int best_i = -1;
			for (int i=0; i<nproj; i++) {
				if( taken[i] ) {
					float v = abs(vecproj[i*3]*x+vecproj[i*3+1]*y+vecproj[i*3+2]*z);
					if (v > best_v) {
						best_v = v;
						best_i = i;
					}
				}
			}
			asg[k*howmany + h] = best_i;
			taken[best_i] = false;
		}
	}
	return asg;
}


vector<int> Util::group_proj_by_phitheta(const vector<float>& projangles, const vector<float>& ref_ang, const int img_per_grp) {
	float c = 100.0;
	int L = max(100, img_per_grp);
	int N = projangles.size()/2;

	int sz = ref_ang.size();
	int nref1 = ref_ang[sz-4];
	int nref2 = ref_ang[sz-3];
	int nref3 = ref_ang[sz-2];
	int nref4 = ref_ang[sz-1];
	int nref;

	set<int> pt;
	for (int i=0; i<N; i++) pt.insert(i);
	vector<float> v(N*3, 0.0f);
	for (int i=0; i<N; i++) 
		getvec(projangles[i*2], projangles[i*2+1], v[i*3], v[i*3+1], v[i*3+2], 1);

	int previous_group = -1;
	int previous_zone = 5;
	int max_group = 0;
	vector<float> ref_ang_list;
	vector<float> diff_table;
	map<int, int> diff_table_index;
	vector<int> proj_list;
	vector<int> sg;
	vector<int> remain_index;
	vector<int> asg;
	int mirror;
	for (int grp=0; grp<N/img_per_grp; grp++) {
		int N_remain = N-grp*img_per_grp;
		assert(N_remain == static_cast<int>(pt.size()));
		if (N_remain >= nref4*L) {
			if (previous_zone > 4) {
				ref_ang_list.resize(nref4*2);
				for (int i=0; i<nref4*2; i++)  ref_ang_list[i] = ref_ang[(nref1+nref2+nref3)*2+i];
				nref = nref4;
				previous_group = -1;
				previous_zone = 4;
			}
		} else if (N_remain >= nref3*L) {
			if (previous_zone > 3) {
				ref_ang_list.resize(nref3*2);
				for (int i=0; i<nref3*2; i++)  ref_ang_list[i] = ref_ang[(nref1+nref2)*2+i];
				nref = nref3;
				previous_group = -1;
				previous_zone = 3;
			}
		} else if (N_remain >= nref2*L) {
			if (previous_zone > 2) {
				ref_ang_list.resize(nref2*2);
				for (int i=0; i<nref2*2; i++)  ref_ang_list[i] = ref_ang[nref1*2+i];
				nref = nref2;
				previous_group = -1;
				previous_zone = 2;
			}
		} else if (N_remain >= nref1*L) {
			if (previous_zone > 1) {
				ref_ang_list.resize(nref1*2);
				for (int i=0; i<nref1*2; i++)  ref_ang_list[i] = ref_ang[i];
				nref = nref1;
				previous_group = -1;
				previous_zone = 1;
			}
		} else if (previous_zone > 0) {
			previous_group = -1;
			previous_zone = 0;
		}

		vector<int> index;
		if (N_remain >=  nref1*L) {
			if (previous_group == -1) { // which means it just changed zone
				vector<float> proj_ang_list(N_remain*2, 0.0f);
				remain_index.resize(N_remain);
				int l = 0;
				for (set<int>::const_iterator si = pt.begin(); si != pt.end(); ++si) {
					remain_index[l] = (*si);
					proj_ang_list[l*2] = projangles[(*si)*2];
					proj_ang_list[l*2+1] = projangles[(*si)*2+1];
					l++; 
				}
				assert(N_remain == l);
				asg = assign_projangles(proj_ang_list, ref_ang_list);
				sg.resize(nref);
				for (int i=0; i<nref; i++) sg[i] = 0;
				for (int i=0; i<N_remain; i++)  sg[asg[i]]++;
			}
			int max_group_size = 0;
			for (int i=0; i<nref; i++)
				if (sg[i] > max_group_size)	{ max_group_size = sg[i]; max_group = i; }
			for (unsigned int i=0; i<remain_index.size(); i++)
				if (asg[i] == max_group)  index.push_back(remain_index[i]);
		} else {
			for (set<int>::const_iterator si = pt.begin(); si != pt.end(); ++si) 
				index.push_back(*si);
			max_group = 0;
		}
			
		int Nn = index.size();
		if (max_group != previous_group) {
			diff_table.resize(Nn*Nn);
			diff_table_index.clear();
			for (int i=0; i<Nn-1; i++)
				for (int j=i+1; j<Nn; j++) {
					float diff = ang_diff(v[index[i]*3], v[index[i]*3+1], v[index[i]*3+2], v[index[j]*3], v[index[j]*3+1], v[index[j]*3+2], mirror);
					float q = exp(-c*pow(diff/180.0f*static_cast<float>(M_PI), 2.0f));
					diff_table[i*Nn+j] = q;
					diff_table[j*Nn+i] = q;
				}
			for (int i=0; i<Nn; i++)  {
				diff_table[i*Nn+i] = 0.0f;      // diagonal values
				diff_table_index[index[i]] = i;	
			}
			previous_group = max_group;
		} 

		int diff_table_size = static_cast<int>(sqrt((float)diff_table.size())+0.5f);
		float max_density = -1;
		int max_density_i = -1;
		for (int i=0; i<Nn; i++) {
			float s = 0.0f;
			int z = diff_table_index[index[i]];
			for (int j=0; j<diff_table_size; j++)  s += diff_table[z*diff_table_size+j];
			if (s > max_density) {
				max_density = s;
				max_density_i = i;
			}
		}

		vector<d_ang> dang(Nn, d_ang(0.0, 0, 0));
		for (int i=0; i<Nn; i++) {
			dang[i].d = ang_diff(v[index[i]*3], v[index[i]*3+1], v[index[i]*3+2], 
			  v[index[max_density_i]*3], v[index[max_density_i]*3+1], v[index[max_density_i]*3+2], mirror);
			dang[i].mirror = mirror;
			dang[i].i = i;
		}
		dang[max_density_i].d = -1;
		sort(dang.begin(), dang.end());		

		for (int i=0; i<img_per_grp; i++) {
			int idd = index[dang[i].i];
			mirror = dang[i].mirror;
			for (unsigned int j=0; j<remain_index.size(); j++)
				if (idd == remain_index[j]) asg[j] = -1;
			for (int j=0; j<diff_table_size; j++) {
				diff_table[diff_table_index[idd]*diff_table_size+j] = 0.0f;
				diff_table[diff_table_index[idd]+diff_table_size*j] = 0.0f;
			}
			proj_list.push_back(mirror*idd);
			pt.erase(idd);
		}
		sg[max_group] -= img_per_grp;
	}
	for (set<int>::const_iterator si = pt.begin(); si != pt.end(); ++si) {
		proj_list.push_back(*si);
	}
	return proj_list;
}

vector<float> Util::multiref_polar_ali_2d_delta(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, string mode,
                vector<int>numr, float cnx, float cny, float delta_start, float delta) {

    // Manually extract.
/*    vector< EMAN::EMData* > crefim;
    std::size_t crefim_len = PyObject_Length(crefim_list.ptr());
    crefim.reserve(crefim_len);

    for(std::size_t i=0;i<crefim_len;i++) {
        boost::python::extract<EMAN::EMData*> proxy(crefim_list[i]);
        crefim.push_back(proxy());
    }
*/

	size_t crefim_len = crefim.size();
	const float qv = static_cast<float>( pi/180.0 );

// 	int   ky = int(2*yrng/step+0.5)/2;
// 	int   kx = int(2*xrng/step+0.5)/2;
	int lkx = int(xrng[0]/step);
	int rkx = int(xrng[1]/step);
	int lky = int(yrng[0]/step);
	int rky = int(yrng[1]/step);
	int   iref, nref=0, mirror=0;
	float iy, ix, sx=0, sy=0;
	float peak = -1.0E23f;
	float ang=0.0f;
	
	for (int i = -lky; i <= rky; i++) {
		iy = i * step ;
		for (int j = -lkx; j <= rkx; j++) {
			ix = j*step ;
			EMData* cimage = Polar2Dm(image, cnx+ix, cny+iy, numr, mode);

			Normalize_ring( cimage, numr );

			Frngs(cimage, numr);
			//  compare with all reference images
			// for iref in xrange(len(crefim)):
			for ( iref = 0; iref < (int)crefim_len; iref++) {
				Dict retvals = Crosrng_ms_delta(crefim[iref], cimage, numr, delta_start, delta);
				double qn = retvals["qn"];
				double qm = retvals["qm"];
				if(qn >= peak || qm >= peak) {
					sx = -ix;
					sy = -iy;
					nref = iref;
					if (qn >= qm) {
						ang = ang_n(retvals["tot"], mode, numr[numr.size()-1]);
						peak = static_cast<float>(qn);
						mirror = 0;
					} else {
						ang = ang_n(retvals["tmt"], mode, numr[numr.size()-1]);
						peak = static_cast<float>(qm);
						mirror = 1;
					}
				}
			}  delete cimage; cimage = 0;
		}
	}
	float co, so, sxs, sys;
	co = static_cast<float>(  cos(ang*qv) );
	so = static_cast<float>( -sin(ang*qv) );
	sxs = sx*co - sy*so;
	sys = sx*so + sy*co;
	vector<float> res;
	res.push_back(ang);
	res.push_back(sxs);
	res.push_back(sys);
	res.push_back(static_cast<float>(mirror));
	res.push_back(static_cast<float>(nref));
	res.push_back(peak);
	return res;
}

vector<float> Util::multiref_polar_ali_2d_nom(EMData* image, const vector< EMData* >& crefim,
                float xrng, float yrng, float step, string mode,
                vector< int >numr, float cnx, float cny) {

    // Manually extract.
/*    vector< EMAN::EMData* > crefim;
    std::size_t crefim_len = PyObject_Length(crefim_list.ptr());
    crefim.reserve(crefim_len);

    for(std::size_t i=0;i<crefim_len;i++) {
        boost::python::extract<EMAN::EMData*> proxy(crefim_list[i]);
        crefim.push_back(proxy());
    }
*/
	size_t crefim_len = crefim.size();
	const float qv = static_cast<float>( pi/180.0 );

	int   ky = int(2*yrng/step+0.5)/2;
	int   kx = int(2*xrng/step+0.5)/2;
	int   iref, nref=0;
	float iy, ix, sx=0, sy=0;
	float peak = -1.0E23f;
	float ang=0.0f;
	int nring = numr.size()/3;
	int inr = numr[(nring-1)*3];
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	for (int i = -ky; i <= ky; i++) {
		iy = i * step ;
		//here the centers of the image cny and cnx use Fortran convention which means the index of the image array starts from 1 to nx (ny). 02-23-2015.
		if ( inr + (int)(cny+iy) <= ny - 1  &&  -inr + (int)(cny+iy) >= 1) {
			for (int j = -kx; j <= kx; j++) {
				ix = j*step ;
				if ( inr + (int)(cnx+ix) <= nx - 1  &&  -inr + (int)(cnx+ix) >= 1) {
					EMData* cimage = Polar2Dm(image, cnx+ix, cny+iy, numr, mode);
					Frngs(cimage, numr);
					//  compare with all reference images
					// for iref in xrange(len(crefim)):
					for ( iref = 0; iref < (int)crefim_len; iref++) {
						Dict retvals = Crosrng_ns(crefim[iref], cimage, numr);
						double qn = retvals["qn"];
						if(qn >= peak) {
							sx = -ix;
							sy = -iy;
							nref = iref;
							ang = ang_n(retvals["tot"], mode, numr[numr.size()-1]);
							peak = static_cast<float>(qn);
						}
					}  delete cimage; cimage = 0;
				}
			}		
		}
	}
	float co, so, sxs, sys;
	co = static_cast<float>(  cos(ang*qv) );
	so = static_cast<float>( -sin(ang*qv) );
	sxs = sx*co - sy*so;
	sys = sx*so + sy*co;
	vector<float> res;
	res.push_back(ang);
	res.push_back(sxs);
	res.push_back(sys);
	res.push_back(static_cast<float>(nref));
	res.push_back(peak);
	return res;
}

vector<float> Util::multiref_polar_ali_2d_local(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, float ant, string mode,
                vector<int>numr, float cnx, float cny, string sym) {

	int lkx = int(xrng[0]/step);
	int rkx = int(xrng[1]/step);
	int lky = int(yrng[0]/step);
	int rky = int(yrng[1]/step);
	int   iref, nref=0, mirror=0, report_mirror=0;
	float iy, ix, sx=0, sy=0;
	float peak = -1.0E23f;
	float ang  = 0.0f;

	size_t crefim_len = crefim.size();
	const float qv = static_cast<float>( pi/180.0 );

	Transform * t = image->get_attr("xform.projection");
//	Dict d = t->get_params("spider");
//	float theta1 = d["theta"];
//	mirror = (int)(theta1 > 90.0f);

    // sym is a symmetry string, t is the projection transform of the input image, tsym is a vector of transforms that are input transofrmation multiplied by all symmetries.
    // its length is number of symmetries.
    vector<Transform> tsym = t->get_sym_proj(sym);

    int isym = 0;
    int nsym = tsym.size();


    // Ims structure of 3 floats defined here  NXCt3rIDWnWrG2bj
    vector<Ims> vIms(nsym);

    for (isym = 0; isym < nsym; ++isym) {
        Dict u = tsym[isym].get_params("spider");
        float phi   = u["phi"];
        float theta = u["theta"];
        vIms[isym].ims1 = sin(theta*qv)*cos(phi*qv);
        vIms[isym].ims2 = sin(theta*qv)*sin(phi*qv);
        vIms[isym].ims3 = cos(theta*qv);
    }

	for (iref = 0; iref < crefim_len; iref++) {

	    float ref_theta = crefim[iref]->get_attr("theta");
	    float ref_phi = crefim[iref]->get_attr("phi");
	    float ref_psi = crefim[iref]->get_attr("psi");

		float n1 = crefim[iref]->get_attr("n1");
		float n2 = crefim[iref]->get_attr("n2");
		float n3 = crefim[iref]->get_attr("n3");



		for (isym = 0; isym < nsym; ++isym) {
			float dot_product = n1*vIms[isym].ims1 + n2*vIms[isym].ims2 + n3*vIms[isym].ims3;
			if(abs(dot_product)>=ant) {
				mirror = dot_product < 0;

				for (int i = -lky; i <= rky; i++) {
					iy = i * step ;
					for (int j = -lkx; j <= rkx; j++) {
						ix = j*step;
						EMData* cimage = Polar2Dm(image, cnx+ix, cny+iy, numr, mode);
						Normalize_ring( cimage, numr );
						Frngs(cimage, numr);
						//  compare with all reference images that are on a new list
							Dict retvals = Crosrng_e(crefim[iref], cimage, numr, mirror);
							double qn = retvals["qn"];

							if(qn >= peak) {
								report_mirror = mirror;
								sx = -ix;
								sy = -iy;
								nref = iref;
								ang = ang_n(retvals["tot"], mode, numr[numr.size()-1]);
								peak = static_cast<float>( qn );
							}
						  delete cimage; cimage = 0;
					}
				}
				break;
			}
		}
	}

	float co, so, sxs, sys;
	if(peak == -1.0E23) {
		ang=0.0; sxs=0.0; sys=0.0; mirror=0;
		nref = -1;
	} else {
		co  =  cos(ang*qv);
		so  = -sin(ang*qv);
		sxs = sx*co - sy*so;
		sys = sx*so + sy*co;
		mirror = report_mirror;
	}

	vector<float> res;
	res.push_back(ang);
	res.push_back(sxs);
	res.push_back(sys);
	res.push_back(static_cast<float>(mirror));
	res.push_back(static_cast<float>(nref));
	res.push_back(peak);
	return res;
}


// Needed for debugging
//      printf("\n HHHHHHHHHHHHHHHH \n");
//      const char * myfilename01 = "/users/hvoicu/EMAN2/src/test_scripts/ref_angles02.txt";
//      FILE *fp = fopen(myfilename01, "w");
//
//        for(int i = 0; i < index_crefim.size(); ++i) {
//          float phi = crefim[index_crefim[i]]->get_attr("phi");
//          float theta = crefim[index_crefim[i]]->get_attr("theta");
//          float psi = crefim[index_crefim[i]]->get_attr("psi");
//          cout << phi << " "  << theta << " " << psi << " 0 0"   << endl;
//          fprintf(fp, "%f %f %f 0 0\n", phi, theta, psi);
//        }
//
//      fflush(fp);
//      fclose(fp);
//      exit(0);

vector<float> Util::shc(EMData* image, const vector< EMData* >& crefim,
				vector<float> xrng, vector<float> yrng, float step, float ant, string mode,
				vector<int>numr, float cnx, float cny, string sym) {

	size_t crefim_len = crefim.size();
	const float qv = static_cast<float>( pi/180.0 );
    Transform * t = 0;

    float ang = 0.0f;
    float sx=0.0f, sy=0.0f;
    float peak = -1.0e23f;
    int   nref = -1, mirror = 0;
    bool  found_better = false;
    size_t tiref = 0;
    float an;

	// cout << ant <<endl;
    if( ant > 0.0f) {  //LOCAL SEARCHES
    /*
    Sequence of operations:
    	1. image should have reference projection orientation xform.anchor in addition to current projection orientation
    	2. find a subset of reference projections within ant of the reference projection orientation of image,
    			the similarity should include the symmetry
    	3. shuffle the indexes of the subset
    	4. Do the random matching
    With the above, local searches could be also used for randomization tests as it would be enough to call
    shc with previousmax=-1.0e23 and reference projection direction set to the current angle.
    */
		t = image->get_attr("xform.anchor");
		Dict d = t->get_params("spider");
		//float phi   = d["phi"];
		float theta   = d["theta"];
		float psi   = d["psi"];
        int maxrin = numr[numr.size()-1];

		//printf("\nINPUT   %f   %f   %f   %f   %d \n",phi, theta, psi, 360.0f-psi,mirror);
		mirror = (int)(theta > 90.0f);
        // mirror = 0 then use -psi;
		int psi_pos = (int)fmod((float)Util::round((-psi - mirror*180)/360.0*maxrin+10*maxrin),(float)maxrin) + 1;

		//  Multiply anchor direction object t by all symmetry group rotations.
		vector<Transform> tsym = t->get_sym_proj(sym);
		int isym = 0;
		int nsym = tsym.size();
		vector<Ims> vIms(nsym);

		for (isym = 0; isym < nsym; ++isym) {
			Dict u = tsym[isym].get_params("spider");
			float phi   = u["phi"];
			float theta = u["theta"];
//			float psi = u["psi"];
//			printf("\n%f %f %f\n",phi, theta, psi);
			vIms[isym].ims1 = sin(theta*qv)*cos(phi*qv);
			vIms[isym].ims2 = sin(theta*qv)*sin(phi*qv);
			vIms[isym].ims3 = cos(theta*qv);
		}

		//  extract indexes of reference images that are within predefined angular distance from the anchor direction.
		vector<int> index_crefim;
		vector<int> mirror_crefim;

		for (unsigned i = 0; i < crefim_len; i++) {
			float n1 = crefim[i]->get_attr("n1");
			float n2 = crefim[i]->get_attr("n2");
			float n3 = crefim[i]->get_attr("n3");
			//  for point-group symmetry get any close symmetry-related reference image
			for (isym = 0; isym < nsym; ++isym) {

                float dot_product = n1*vIms[isym].ims1 + n2*vIms[isym].ims2 + n3*vIms[isym].ims3;
                if(abs(dot_product)>=ant) {
                    mirror_crefim.push_back(int(dot_product < 0));
					index_crefim.push_back(i);
					break;
				}
			}
		}

		 an = (float)acos(ant) / qv;

		const float previousmax = image->get_attr("previousmax");
		//printf("\n  previousmax   %f  \n",previousmax);
        crefim_len = index_crefim.size();

        if (crefim_len > 0) {
			const int lkx = int(xrng[0]/step);
			const int rkx = int(xrng[1]/step);
			const int lky = int(yrng[0]/step);
			const int rky = int(yrng[1]/step);

			vector< vector<EMData*> > cimages( lky+rky+1, vector<EMData*>(lkx+rkx+1) );
			for (int i = -lky; i <= rky; i++) {
				const int iy = i * step ;
				for (int j = -lkx; j <= rkx; j++) {
					const int ix = j*step ;
					EMData* cimage = Polar2Dm(image, cnx+ix, cny+iy, numr, mode);
					Normalize_ring( cimage, numr );
					Frngs(cimage, numr);
					cimages[i+lky][j+lkx] = cimage;
				}
			}

			for (unsigned i = 0; i < crefim_len; ++i) {
				unsigned r = Util::get_irand(0,crefim_len-1);
				swap( index_crefim[r], index_crefim[i] );
			}
			bool found_better = false;
			for ( ;  (tiref < crefim_len) && (! found_better); tiref++) {
				const int iref = index_crefim[tiref];
				std::vector<int> shifts = shuffled_range( 0, (lkx+rkx+1) * (lky+rky+1) - 1 );
				for ( unsigned nodeId = 0;  nodeId < shifts.size();  ++nodeId ) {
					const int i = ( shifts[nodeId] % (lky+rky+1) ) - lky;
					const int j = ( shifts[nodeId] / (lky+rky+1) ) - lkx;
					const float iy = i * step;
					const float ix = j * step;
					EMData* cimage = cimages[i+lky][j+lkx];
					Dict retvals = Crosrng_rand_e(crefim[iref], cimage, numr, mirror_crefim[iref], previousmax, an, psi_pos);
					const float new_peak = static_cast<float>( retvals["qn"] );
					//cout << new_peak <<endl;
					if (new_peak > peak) {
						sx = -ix;
						sy = -iy;
						nref = iref;
						ang = ang_n(retvals["tot"], mode, maxrin);
						peak = new_peak;
						mirror = mirror_crefim[iref];
						found_better = (peak > previousmax);
						//cout << found_better <<endl;
						// jump out from search
						if (found_better) break;
					}
				}
			}

			for (unsigned i = 0; i < cimages.size(); ++i) {
				for (unsigned j = 0; j < cimages[i].size(); ++j) {
					delete cimages[i][j];
					cimages[i][j] = NULL;
				}
			}
		}
		//cout << "  JUMPED OUT " <<endl;
		const float co =  cos(ang*qv);
		const float so = -sin(ang*qv);
		const float sxs = sx*co - sy*so;
		const float sys = sx*so + sy*co;
		//printf("\n OUTPUT   %f    %d    %f   %f   %d   %d   %f\n",360.0f - psi, psi_pos, ang, peak, mirror, maxrin,an);
		vector<float> res;
		res.push_back(ang);
		res.push_back(sxs);
		res.push_back(sys);
		res.push_back(static_cast<float>(mirror));
		res.push_back(static_cast<float>(nref));
		res.push_back(peak);
		res.push_back(static_cast<float>(tiref));
		return res;

    } else {  // GLOBAL SEARCHES
		const float previousmax = image->get_attr("previousmax");

		vector<unsigned> listr(crefim_len);
		for (unsigned i = 0; i < crefim_len; ++i) listr[i] = i;
		for (unsigned i = 0; i < crefim_len; ++i) {
			unsigned r = Util::get_irand(0,crefim_len-1);
			swap( listr[r], listr[i] );
		}

		const int lkx = int(xrng[0]/step);
		const int rkx = int(xrng[1]/step);
		const int lky = int(yrng[0]/step);
		const int rky = int(yrng[1]/step);

		vector< vector<EMData*> > cimages( lky+rky+1, vector<EMData*>(lkx+rkx+1) );

		for (int i = -lky; i <= rky; i++) {
			const int iy = i * step ;
			for (int j = -lkx; j <= rkx; j++) {
				const int ix = j*step;
				EMData* cimage = Polar2Dm(image, cnx+ix, cny+iy, numr, mode);
				Normalize_ring( cimage, numr );
				Frngs(cimage, numr);
				cimages[i+lky][j+lkx] = cimage;
			}
		}

		for ( ;  (tiref < crefim_len) && (! found_better); tiref++) {
			const int iref = listr[tiref];
			//  Here check for neighbors
			std::vector<int> shifts = shuffled_range( 0, (lkx+rkx+1) * (lky+rky+1) - 1 );
			for ( unsigned nodeId = 0;  nodeId < shifts.size();  ++nodeId ) {
				const int i = ( shifts[nodeId] % (lky+rky+1) ) - lky;
				const int j = ( shifts[nodeId] / (lky+rky+1) ) - lkx;
				const float iy = i * step;
				const float ix = j * step;
				EMData* cimage = cimages[i+lky][j+lkx];

				Dict retvals = Crosrng_rand_ms(crefim[iref], cimage, numr, previousmax);
				const float new_peak = static_cast<float>( retvals["qn"] );
				//cout << new_peak <<endl;
				if (new_peak > peak) {
					sx = -ix;
					sy = -iy;
					nref = iref;
					ang = ang_n(retvals["tot"], mode, numr[numr.size()-1]);
					peak = new_peak;
					mirror = static_cast<int>( retvals["mirror"] );
					found_better = (peak > previousmax);
					//cout << found_better <<endl;
					// jump out from search
					if (found_better) break;
				}
			}
		}
		//cout << "  JUMPED OUT " <<endl;

		for (unsigned i = 0; i < cimages.size(); ++i) {
			for (unsigned j = 0; j < cimages[i].size(); ++j) {
				delete cimages[i][j];
				cimages[i][j] = NULL;
			}
		}

		const float co =  cos(ang*qv);
		const float so = -sin(ang*qv);
		const float sxs = sx*co - sy*so;
		const float sys = sx*so + sy*co;

		vector<float> res;
		res.push_back(ang);
		res.push_back(sxs);
		res.push_back(sys);
		res.push_back(static_cast<float>(mirror));
		res.push_back(static_cast<float>(nref));
		res.push_back(peak);
		res.push_back(static_cast<float>(tiref));
		return res;

	}
}
	
static std::string toString(int i)
{
	std::ostringstream s;
	s << i;
	return s.str();
}

vector<float> Util::shc_multipeaks(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, float ant, string mode,
                vector<int>numr, float cnx, float cny, int max_peaks_count)
{
	size_t crefim_len = crefim.size();
	const float qv = static_cast<float>( pi/180.0 );

	const float previousmax = image->get_attr("previousmax");
//	const int   ky = int(2*yrng/step+0.5)/2;
//	const int   kx = int(2*xrng/step+0.5)/2;
	const int lkx = int(xrng[0]/step);
	const int rkx = int(xrng[1]/step);
	const int lky = int(yrng[0]/step);
	const int rky = int(yrng[1]/step);
	
	int   iref, nref=0, mirror=0;
	float iy, ix, sx=0, sy=0;
	float ang=0.0f;

	vector< vector<EMData*> > cimages( lky+rky+1, vector<EMData*>(lkx+rkx+1) );
	for (int i = -lky; i <= rky; i++) {
		iy = i * step ;
		for (int j = -lkx; j <= rkx; j++) {
			ix = j*step ;
			EMData* cimage = Polar2Dm(image, cnx+ix, cny+iy, numr, mode);
			Normalize_ring( cimage, numr );
			Frngs(cimage, numr);
			cimages[i+lky][j+lkx] = cimage;
		}
	}

	vector<unsigned> listr(crefim_len);
	for (unsigned i = 0; i < crefim_len; ++i) listr[i] = i;
	for (unsigned i = 0; i < crefim_len; ++i) {
		unsigned r = Util::get_irand(0,crefim_len-1);
		swap( listr[r], listr[i] );
	}

	std::vector<float> results;
	size_t tiref = 0;
	for ( ;  (tiref < crefim_len) && (results.size() / 7 < static_cast<unsigned>(max_peaks_count)); tiref++) {
		iref = listr[tiref];
		float peak = previousmax;

		std::vector<int> shifts = shuffled_range( 0, (lkx+rkx+1) * (lky+rky+1) - 1 );
		//for ( unsigned nodeId = 0;  nodeId < shifts.size();  ++nodeId ) {
		unsigned nodeId = 0;
		bool found_better = false;
		while(nodeId < shifts.size()  &&  !found_better  && (results.size() / 7 < static_cast<unsigned>(max_peaks_count)) ) {
			const int i = ( shifts[nodeId] % (lky+rky+1) ) - lky;
			const int j = ( shifts[nodeId] / (lky+rky+1) ) - lkx;
			const float iy = i * step;
			const float ix = j * step;
			EMData* cimage = cimages[i+lky][j+lkx];
			Dict retvals = Crosrng_rand_ms(crefim[iref], cimage, numr, previousmax);
			const float new_peak = static_cast<float>( retvals["qn"] );
			if (new_peak > peak) {
				found_better = true;//(peak > previousmax);

				sx = -ix;
				sy = -iy;
				ang = ang_n(retvals["tot"], mode, numr[numr.size()-1]);
				peak = new_peak;
				mirror = static_cast<int>( retvals["mirror"] );

				const float co =  cos(ang*qv);
				const float so = -sin(ang*qv);
				const float sxs = sx*co - sy*so;
				const float sys = sx*so + sy*co;

				results.push_back(ang);
				results.push_back(sxs);
				results.push_back(sys);
				results.push_back(static_cast<float>(mirror));
				results.push_back(static_cast<float>(iref));
				results.push_back(peak);
				results.push_back(static_cast<float>(tiref));
			}
			++nodeId;
		}
    }
	for (unsigned i = 0; i < cimages.size(); ++i) {
		for (unsigned j = 0; j < cimages[i].size(); ++j) {
			delete cimages[i][j];
			cimages[i][j] = NULL;
		}
	}
/*
	// sorting
	unsigned no_of_solution = results.size() / 7;
	for (unsigned i = 0; i < no_of_solution; ++i) {
		unsigned max_peak_ind = i;
		float max_peak_val = results[7*i+5];
		for (unsigned j = i+1; j < no_of_solution; ++j) {
			const float peak_j = results[7*j+5];
			if (peak_j > max_peak_val) {
				max_peak_val = peak_j;
				max_peak_ind = j;
			}
		}
		if (max_peak_ind != i) {
			for (unsigned j = 0; j < 7; ++j) {
				std::swap(results[7*i+j], results[7*max_peak_ind+j]);
			}
		}
	}

	// set new previousmax to the middle value
	if ( no_of_solution > 0 ) {
		image->set_attr("previousmax",results[7*(no_of_solution/2)+5]);
	}

	// set weights
	if ( no_of_solution < 2 ) {
		image->set_attr("weight",1.0);
	} else {
		// the idea is like this:
		// - the minimal weight should equal 1/no_of_solution of maximal weight
		// - the sum of weights should equal 1.0
		vector<float> w(no_of_solution);
		for (unsigned i = 0; i < no_of_solution; ++i) {
			w[i] = results[7*i+5];
		}
		//float r = (w.front() - no_of_solution * w.back()) / (no_of_solution - 1);
		float sum = 0.0f;
		for (unsigned i = 0; i < no_of_solution; ++i) {
			//w[i] += r;
			sum += w[i];
		}
		for (unsigned i = 0; i < no_of_solution; ++i) {
			w[i] /= sum;
			if (i == 0) {
				image->set_attr("weight",w[i]);
			} else {
				image->set_attr("weight" + toString(i),w[i]);
			}
		}
	}

	// remove unused parameters from the header
	unsigned i0 = (no_of_solution < 1) ? 1 : no_of_solution;
	for (unsigned i = i0; image->has_attr("weight" + toString(i));  ++i ) {
		image->del_attr("weight" + toString(i));
	}
*/
	return results;
}

vector<float> Util::multiref_polar_ali_2d_local_psi(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, float ant, float psi_max, string mode,
                vector<int>numr, float cnx, float cny) {

    // Manually extract.
/*    vector< EMAN::EMData* > crefim;
    std::size_t crefim_len = PyObject_Length(crefim_list.ptr());
    crefim.reserve(crefim_len);

    for(std::size_t i=0;i<crefim_len;i++) {
        boost::python::extract<EMAN::EMData*> proxy(crefim_list[i]);
        crefim.push_back(proxy());
    }
*/
	size_t crefim_len = crefim.size();
	const float qv = static_cast<float>(pi/180.0);

	Transform* t = image->get_attr("xform.projection");
	Dict d = t->get_params("spider");
	if(t) {delete t; t=0;}
	float phi = d["phi"];
	float theta = d["theta"];
	float psi = d["psi"];
// 	int ky = int(2*yrng/step+0.5)/2;
// 	int kx = int(2*xrng/step+0.5)/2;
	int lkx = int(xrng[0]/step);
	int rkx = int(xrng[1]/step);
	int lky = int(yrng[0]/step);
	int rky = int(yrng[1]/step);
	int iref, nref = 0, mirror = 0;
	float iy, ix, sx = 0, sy = 0;
	float peak = -1.0E23f;
	float ang = 0.0f;
	float imn1 = sin(theta*qv)*cos(phi*qv);
	float imn2 = sin(theta*qv)*sin(phi*qv);
	float imn3 = cos(theta*qv);
	vector<float> n1(crefim_len);
	vector<float> n2(crefim_len);
	vector<float> n3(crefim_len);
	for (iref = 0; iref < (int)crefim_len; iref++) {
			n1[iref] = crefim[iref]->get_attr("n1");
			n2[iref] = crefim[iref]->get_attr("n2");
			n3[iref] = crefim[iref]->get_attr("n3");
	}
	bool nomirror = (theta<90.0) || (theta==90.0);
	if (!nomirror) {
		phi   = fmod(phi+540.0f, 360.0f);
		theta = 180-theta;
		psi   = fmod(540.0f-psi, 360.0f);
	} else { psi = fmod(360.0f-psi, 360.0f); }
	
	for (int i = -lky; i <= rky; i++) {
		iy = i * step ;
		for (int j = -lkx; j <= rkx; j++) {
			ix = j*step ;
			EMData* cimage = Polar2Dm(image, cnx+ix, cny+iy, numr, mode);

			Normalize_ring(cimage, numr);

			Frngs(cimage, numr);
			//  compare with all reference images
			// for iref in xrange(len(crefim)):
			for (iref = 0; iref < (int)crefim_len; iref++) {
				if (abs(n1[iref]*imn1 + n2[iref]*imn2 + n3[iref]*imn3)>=ant) {
					float refpsi = crefim[iref]->get_attr("psi");
					if (nomirror) {
						Dict retvals = Crosrng_sm_psi(crefim[iref], cimage, numr, fmod(360.0+psi+refpsi, 360.0), 0, psi_max);
						double qn = retvals["qn"];
						if (qn >= peak) {
							sx = -ix;
							sy = -iy;
							nref = iref;
							ang = ang_n(retvals["tot"], mode, numr[numr.size()-1]);
							peak = static_cast<float>(qn);
							mirror = 0;
						}
					} else {
						Dict retvals = Crosrng_sm_psi(crefim[iref], cimage, numr, fmod(360.0+psi-refpsi, 360.0), 1, psi_max);
						double qn = retvals["qn"];
						if (qn >= peak) {
							sx = -ix;
							sy = -iy;
							nref = iref;
							ang = ang_n(retvals["tot"], mode, numr[numr.size()-1]);
							peak = static_cast<float>(qn);
							mirror = 1;
						}
					}
				}
			}  delete cimage; cimage = 0;
		}
	}
	float co, so, sxs, sys;
	if(peak == -1.0E23) {
		ang=0.0; sxs=0.0; sys=0.0; mirror=0;
		nref = -1;
	} else {
		co =  cos(ang*qv);
		so = -sin(ang*qv);
		sxs = sx*co - sy*so;
		sys = sx*so + sy*co;
	}
	vector<float> res;
	res.push_back(ang);
	res.push_back(sxs);
	res.push_back(sys);
	res.push_back(static_cast<float>(mirror));
	res.push_back(static_cast<float>(nref));
	res.push_back(peak);
	return res;
}


vector<float> Util::multiref_polar_ali_helical(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, float psi_max, string mode,
                vector<int>numr, float cnx, float cny, int ynumber) {
	
	size_t crefim_len = crefim.size();
	const float qv = static_cast<float>( pi/180.0 );

	int   iref, nref=0, mirror=0;
	float iy, ix, sx=0, sy=0;
	float peak = -1.0E23f;
	float ang=0.0f;
// 	int   kx = int(2*xrng/step+0.5)/2;
	int lkx = int(xrng[0]/step);
	int rkx = int(xrng[1]/step);
	int lky, rky;
	
		
	//if ynumber==-1, use the old code which process x and y direction equally.
	//if ynumber is given, it should be even. We need to check whether it is zero

//	int ky;
	float stepy;
	int kystart;
	
	if (ynumber == -1) {
//	    ky = int(2*yrng/step+0.5)/2;
		lky = int(yrng[0]/step);
		rky = int(yrng[1]/step);
	    stepy = step;
	    kystart = -lky;
	} else if(ynumber == 0) {
	     rky = 0;
	   	 stepy = 0.0f;
	   	 kystart = rky;
	} else {
//	    ky = int(ynumber/2);		
		stepy=(yrng[0]+yrng[1])/ynumber;		//2-25-2015@ming 
		lky = int(yrng[0]/stepy);
		rky = int(yrng[1]/stepy);
//		kystart = -lky + 1;  
		kystart = -lky;   
	}
	//std::cout<<"yrng="<<yrng<<"ynumber="<<ynumber<<"stepy=="<<stepy<<"stepx=="<<step<<std::endl;
	for (int i = kystart; i <= rky; i++) {
		iy = i * stepy ;
		for (int j = -lkx; j <= rkx; j++) {
			ix = j*step ;
			EMData* cimage = Polar2Dm(image, cnx+ix, cny+iy, numr, mode);

			Normalize_ring( cimage, numr );

			Frngs(cimage, numr);
			//  compare with all reference images
			// for iref in xrange(len(crefim)):
			for ( iref = 0; iref < (int)crefim_len; iref++) {
				Dict retvals_0   = Crosrng_psi(crefim[iref], cimage, numr,   0, psi_max);
				Dict retvals_180 = Crosrng_psi(crefim[iref], cimage, numr, 180, psi_max);
				double qn_0   = retvals_0["qn"];
				double qn_180 = retvals_180["qn"];
				double qm_0   = retvals_0["qm"];
				double qm_180 = retvals_180["qm"];
				double qn;
				double qm;
				bool   qn_is_zero = false;
				bool   qm_is_zero = false;

				if (qn_0 >= qn_180) {
					qn = qn_0;
					qn_is_zero = true;
				} else {
					qn = qn_180;
					qn_is_zero = false; 
				}

				if (qm_0 >= qm_180) {
					qm = qm_0;
					qm_is_zero = true;
				} else {
					qm = qm_180;
					qm_is_zero = false; 
				}

				if(qn >= peak || qm >= peak) {
					sx = -ix;
					sy = -iy;
					nref = iref;
					if (qn >= qm) {
						if (qn_is_zero) ang = ang_n(retvals_0["tot"], mode, numr[numr.size()-1]);
						else            ang = ang_n(retvals_180["tot"], mode, numr[numr.size()-1]);
						peak = static_cast<float>(qn);
						mirror = 0;
					} else {
						if (qm_is_zero) ang = ang_n(retvals_0["tmt"], mode, numr[numr.size()-1]);
						else            ang = ang_n(retvals_180["tmt"], mode, numr[numr.size()-1]);
						peak = static_cast<float>(qm);
						mirror = 1;
					}
				}
			}
			delete cimage; cimage = 0;
		}
	}
	float co, so, sxs, sys;
	co = static_cast<float>(  cos(ang*qv) );
	so = static_cast<float>( -sin(ang*qv) );
	sxs = sx*co - sy*so;
	sys = sx*so + sy*co;
	vector<float> res;
	res.push_back(ang);
	res.push_back(sxs);
	res.push_back(sys);
	res.push_back(static_cast<float>(mirror));
	res.push_back(static_cast<float>(nref));
	res.push_back(peak);
	return res;
}

vector<float> Util::multiref_polar_ali_helical_local(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, float ant, float psi_max, string mode,
                vector<int>numr, float cnx, float cny, int ynumber, float yrnglocal) {
	
	size_t crefim_len = crefim.size();
	const float qv = static_cast<float>( pi/180.0 );

	int   iref, nref=-1, mirror=0;
	float iy, ix, sx=0, sy=0;
	float peak = -1.0E23f;
	float jneginf = -2.0E23f;
	float ang=0.0f;
	Transform * t = image->get_attr("xform.projection");
	Dict d = t->get_params("spider");
	if(t) {delete t; t=0;}
	float phi       = d["phi"];
	float theta     = d["theta"];
	float psi       = d["psi"];
	float imn1 = sin(theta*qv)*cos(phi*qv);
	float imn2 = sin(theta*qv)*sin(phi*qv);
	float imn3 = cos(theta*qv);
	vector<float> n1(crefim_len);
	vector<float> n2(crefim_len);
	vector<float> n3(crefim_len);
	for ( iref = 0; iref < (int)crefim_len; iref++) {
		n1[iref] = crefim[iref]->get_attr("n1");
		n2[iref] = crefim[iref]->get_attr("n2");
		n3[iref] = crefim[iref]->get_attr("n3");		
	}
	float nbrinp;
	float nbrinp_mirror;
	bool  use_ref;
	bool  use_ref_mirror;
//	int   kx = int(2*xrng/step+0.5)/2;
//	int   ky;
	int lkx = int(xrng[0]/step);
	int rkx = int(xrng[1]/step);
	int lky, rky;
	
	float stepy;
	
	if (ynumber == 0) {
		lky = rky = 0;
	} else { 

		if (ynumber > 0) stepy=(yrng[0]+yrng[1])/ynumber;	//2-25-2015@ming
		else if (ynumber == -1) stepy = step;
		
		if (yrnglocal >= 0.0) {
//		   	ky = int(yrnglocal/stepy);
		   	lky = min(int(yrng[0]/stepy), int(yrnglocal/stepy));    //2-25-2015@ming
		   	rky = min(int(yrng[1]/stepy), int(yrnglocal/stepy));
		} else { // search range is not restricted
				lky = int(yrng[0]/stepy);                   //2-25-2015@ming
				rky = int(yrng[1]/stepy);
//			if (ynumber > 0) {
//				ky = int(ynumber/2);	
//			} // else {
// 				ky = int(2*yrng/stepy+0.5)/2;	
// 			}	
					
		}
	}

	for (int i = -lky; i <= rky; i++) {
		iy = i * stepy ;
		for (int j = -lkx; j <= rkx; j++) {
			ix = j*step ;
			EMData* cimage = Polar2Dm(image, cnx+ix, cny+iy, numr, mode);
			Normalize_ring( cimage, numr );
	
			Frngs(cimage, numr);
			//  Compare with All reference images within neighborhood ant
			// for iref in xrange(len(crefim)):
			for ( iref = 0; iref < (int)crefim_len; iref++) {

				use_ref = false;
				use_ref_mirror = false;

				// inner product of iref's Eulerian angles with that of the data
				nbrinp = n1[iref]*imn1 + n2[iref]*imn2 + n3[iref]*imn3;
				if (nbrinp >= ant)   use_ref = true;

				// inner product of the mirror of iref's Eulerian angles with that of the data
				nbrinp_mirror = n3[iref]*imn3 - n1[iref]*imn1 - n2[iref]*imn2;
				if (nbrinp_mirror >= ant)  use_ref_mirror = true;
				if(use_ref || use_ref_mirror){

					Dict retvals;
					Dict retvals_mirror;
					if (use_ref_mirror == true) {
						if ((psi-90.0f) < 90.0f) retvals_mirror = Crosrng_sm_psi(crefim[iref], cimage, numr,   0.0f, 1, psi_max);
						else                     retvals_mirror = Crosrng_sm_psi(crefim[iref], cimage, numr, 180.0f, 1, psi_max);
					}
					if (use_ref == true) {
						if ((psi-90.0f) < 90.0f) retvals = Crosrng_sm_psi(crefim[iref], cimage, numr,   0.0f, 0, psi_max);
						else                     retvals = Crosrng_sm_psi(crefim[iref], cimage, numr, 180.0f, 0, psi_max);
					}
					double qn = jneginf;
					if (use_ref) qn = retvals["qn"];
					double qm = jneginf;
					if (use_ref_mirror) qm = retvals_mirror["qn"];

					if(qn >= peak || qm >= peak) {
						sx = -ix;
						sy = -iy;
						nref = iref;
						if (qn >= qm){
							ang = ang_n(retvals["tot"], mode, numr[numr.size()-1]);
							mirror = 0;
							peak = static_cast<float>(qn);
						} else {
							ang = ang_n(retvals_mirror["tot"], mode, numr[numr.size()-1]);
							mirror = 1;
							peak = static_cast<float>(qm);
						}
					}
				}
			}
			delete cimage; cimage = 0;
		}
	}

	float co, so, sxs, sys;
	co = static_cast<float>(  cos(ang*qv) );
	so = static_cast<float>( -sin(ang*qv) );
	sxs = sx*co - sy*so;
	sys = sx*so + sy*co;
	vector<float> res;
	res.push_back(ang);
	res.push_back(sxs);
	res.push_back(sys);
	res.push_back(static_cast<float>(mirror));
	res.push_back(static_cast<float>(nref));
	res.push_back(peak);
	return res;
}


vector<float> Util::multiref_polar_ali_helical_90(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, float psi_max, string mode,
                vector<int>numr, float cnx, float cny, int ynumber) {

	size_t crefim_len = crefim.size();
	const float qv = static_cast<float>( pi/180.0 );

	int   iref, nref=0, mirror=0;
	float iy, ix, sx=0, sy=0;
	float peak = -1.0E23f;
	float ang=0.0f;
//	int   kx = int(2*xrng/step+0.5)/2;
	int lkx = int(xrng[0]/step);
	int rkx = int(xrng[1]/step);
	//if ynumber==-1, use the old code which process x and y direction equally.

	int lky, rky;
	float stepy;
	int kystart;

	if (ynumber == -1) {
//	    ky = int(2*yrng/step+0.5)/2;
		lky = int(yrng[0]/step);
		rky = int(yrng[1]/step);
	    stepy = step;
	    kystart = -lky;
	} else if(ynumber == 0) {
	     rky = 0;
	   	 stepy = 0.0f;
	   	 kystart = rky;
	} else {
//	    ky = int(ynumber/2);	
	    stepy=(yrng[0]+yrng[1])/ynumber;     //2-25-2015@ming
   		lky = int(yrng[0]/stepy);
		rky = int(yrng[1]/stepy); 	
		kystart = -lky;   // + 1;    
	}

	//std::cout<<"yrng="<<yrng<<"ynumber="<<ynumber<<"stepy=="<<stepy<<"stepx=="<<step<<std::endl;
	for (int i = kystart; i <= rky; i++) {
		iy = i * stepy ;
		for (int j = -lkx; j <= rkx; j++) {
			ix = j*step ;
			EMData* cimage = Polar2Dm(image, cnx+ix, cny+iy, numr, mode);

			Normalize_ring( cimage, numr );

			Frngs(cimage, numr);
			//  compare with all reference images
			// for iref in xrange(len(crefim)):
			for ( iref = 0; iref < (int)crefim_len; iref++) {
				Dict retvals_0   = Crosrng_sm_psi(crefim[iref], cimage, numr,   0.0f, 0, psi_max);
				Dict retvals_180 = Crosrng_sm_psi(crefim[iref], cimage, numr, 180.0f, 0, psi_max);
				double qn_0   = retvals_0["qn"];
				double qn_180 = retvals_180["qn"];
				double qn;
				bool qn_is_zero = false;
		
				if (qn_0 >= qn_180){
					qn = qn_0;
					qn_is_zero = true;
				} else {
					qn = qn_180;
					qn_is_zero = false; 
				}

				if(qn >= peak) {
					sx = -ix;
					sy = -iy;
					nref = iref;
			
					if (qn_is_zero){
						ang = ang_n(retvals_0["tot"], mode, numr[numr.size()-1]);
					} else {
						ang = ang_n(retvals_180["tot"], mode, numr[numr.size()-1]);
					}
					peak = static_cast<float>(qn);
					mirror = 0;
				}
			}
			delete cimage; cimage = 0;
		}	
	}	
	float co, so, sxs, sys;
	co = static_cast<float>(  cos(ang*qv) );
	so = static_cast<float>( -sin(ang*qv) );
	sxs = sx*co - sy*so;
	sys = sx*so + sy*co;
	vector<float> res;
	res.push_back(ang);
	res.push_back(sxs);
	res.push_back(sys);
	res.push_back(static_cast<float>(mirror));
	res.push_back(static_cast<float>(nref));
	res.push_back(peak);
	return res;
}

vector<float> Util::multiref_polar_ali_helical_90_local(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, float ant, float psi_max, string mode,
                vector<int>numr, float cnx, float cny, int ynumber, float yrnglocal) {

	size_t crefim_len = crefim.size();
	const float qv = static_cast<float>( pi/180.0 );
	Transform * t = image->get_attr("xform.projection");
	Dict d = t->get_params("spider");
	if(t) {delete t; t=0;}
	float phi   = d["phi"];
	float theta = d["theta"];
	float psi   = d["psi"];
	float imn1 = sin(theta*qv)*cos(phi*qv);
	float imn2 = sin(theta*qv)*sin(phi*qv);
	float imn3 = cos(theta*qv);
	vector<float> n1(crefim_len);
	vector<float> n2(crefim_len);
	vector<float> n3(crefim_len);
	int   iref, nref=-1, mirror=0;
	float iy, ix, sx=0, sy=0;
	float peak = -1.0E23f;
	float ang  = 0.0f;
//	int   kx   = int(2*xrng/step+0.5)/2;
	int lkx = int(xrng[0]/step);
	int rkx = int(xrng[1]/step);

	for ( iref = 0; iref < (int)crefim_len; iref++) {
		n1[iref] = crefim[iref]->get_attr("n1");
		n2[iref] = crefim[iref]->get_attr("n2");
		n3[iref] = crefim[iref]->get_attr("n3");
	}
	
	float stepy;
	int lky, rky;
	
	if (ynumber == 0) {
		lky = rky = 0;
	} else {
		if (ynumber > 0) stepy=(yrng[0]+yrng[1])/ynumber;
		else if (ynumber == -1) stepy = step;
		
		if (yrnglocal >= 0.0) {
//		   	ky = int(yrnglocal/stepy);
			lky = min(int(yrng[0]/stepy),int(yrnglocal/stepy));
			rky = min(int(yrng[1]/stepy),int(yrnglocal/stepy));
		}
		else { // search range is not restricted
			lky = int(yrng[0]/stepy);
			rky = int(yrng[1]/stepy);
//			if (ynumber > 0) {
//				ky = int(ynumber/2);	
//			} // else {
// 				ky = int(2*yrng/stepy+0.5)/2;	
// 			}			
		}
	}

	for (int i = -lky; i <= rky; i++) {
		iy = i * stepy ;
		for (int j = -lkx; j <= rkx; j++) {
			ix = j*step ;
			EMData* cimage = Polar2Dm(image, cnx+ix, cny+iy, numr, mode);

			Normalize_ring( cimage, numr );

			Frngs(cimage, numr);
			//  compare with all reference images
			// for iref in xrange(len(crefim)):
			for ( iref = 0; iref < (int)crefim_len; iref++) {
				if(fabs(n1[iref]*imn1 + n2[iref]*imn2 + n3[iref]*imn3)>=ant) {
					Dict retvals; 
					if ((psi-90.0f) < 90.0f) retvals = Crosrng_sm_psi(crefim[iref], cimage, numr,   0.0f, 0, psi_max);
					else                     retvals = Crosrng_sm_psi(crefim[iref], cimage, numr, 180.0f, 0, psi_max);

					double qn = retvals["qn"];
					if( qn >= peak) {
						sx = -ix;
						sy = -iy;
						nref = iref;
						ang = ang_n(retvals["tot"], mode, numr[numr.size()-1]);
						peak = static_cast<float>(qn);
						mirror = 0;
					}
				}
			}  
			delete cimage; cimage = 0;
		}
	}

	float co, so, sxs, sys;
	co = static_cast<float>(  cos(ang*qv) );
	so = static_cast<float>( -sin(ang*qv) );
	sxs = sx*co - sy*so;
	sys = sx*so + sy*co;
	vector<float> res;
	res.push_back(ang);
	res.push_back(sxs);
	res.push_back(sys);
	res.push_back(static_cast<float>(mirror));
	res.push_back(static_cast<float>(nref));
	res.push_back(peak);
	return res;
}

// HELICON
vector<float> Util::multiref_polar_ali_helicon_local(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, float ant, float psi_max, string mode,
                vector<int>numr, float cnx, float cny, int ynumber, float yrnglocal) {
	
	size_t crefim_len = crefim.size();
	const float qv = static_cast<float>( pi/180.0 );

	int   iref, nref=-1, mirror=0;
	float iy, ix, sx=0, sy=0;
	float peak = -1.0E23f;
	float jneginf = -2.0E23f;
	float ang=0.0f;
	Transform * t = image->get_attr("xform.projection");
	Dict d = t->get_params("spider");
	if(t) {delete t; t=0;}
	float phi       = d["phi"];
	float theta     = d["theta"];
	float psi       = d["psi"];
	float imn1 = sin(theta*qv)*cos(phi*qv);
	float imn2 = sin(theta*qv)*sin(phi*qv);
	float imn3 = cos(theta*qv);
	vector<float> n1(crefim_len);
	vector<float> n2(crefim_len);
	vector<float> n3(crefim_len);
	for ( iref = 0; iref < (int)crefim_len; iref++) {
		n1[iref] = crefim[iref]->get_attr("n1");
		n2[iref] = crefim[iref]->get_attr("n2");
		n3[iref] = crefim[iref]->get_attr("n3");		
	}
	float nbrinp;
	float nbrinp_mirror;
	bool  use_ref;
	bool  use_ref_mirror;
//	int   kx = int(2*xrng/step+0.5)/2;
	int lkx = int(xrng[0]/step);
	int rkx = int(xrng[1]/step);
	int   lky, rky;
	float stepy;
	
	if (ynumber == 0) {
		lky = rky = 0;
	} else {
		if (ynumber > 0) stepy=(yrng[0]+yrng[1])/ynumber;
		else if (ynumber == -1) stepy = step;

		if (yrnglocal >= 0.0) {
//		   	ky = int(yrnglocal/stepy);
			lky = min(int(yrng[0]/stepy),int(yrnglocal/stepy));
			rky = min(int(yrng[1]/stepy),int(yrnglocal/stepy));
			
		} else { // search range is not restricted
			lky = int(yrng[0]/stepy);
			rky = int(yrng[1]/stepy);
						
//			if (ynumber > 0) {
//				ky = int(ynumber/2);	
//			} // else {
// 				ky = int(2*yrng/stepy+0.5)/2;	
// 			}	
		}
	}

	for (int i = -lky; i <= rky; i++) {
		iy = i * stepy ;
		for (int j = -lkx; j <= rkx; j++) {
			ix = j*step ;
			EMData* cimage = Polar2Dm(image, cnx+ix, cny+iy, numr, mode);
			Normalize_ring( cimage, numr );
	
			Frngs(cimage, numr);
			//  Compare with All reference images within neighborhood ant
			// for iref in xrange(len(crefim)):
			for ( iref = 0; iref < (int)crefim_len; iref++) {

				use_ref = false;
				use_ref_mirror = false;

				// inner product of iref's Eulerian angles with that of the data
				nbrinp = n1[iref]*imn1 + n2[iref]*imn2 + n3[iref]*imn3;
				if (nbrinp >= ant)   use_ref = true;

				// inner product of the mirror of iref's Eulerian angles with that of the data
				nbrinp_mirror = n3[iref]*imn3 - n1[iref]*imn1 - n2[iref]*imn2;
				if (nbrinp_mirror >= ant)  use_ref_mirror = true;
				if(use_ref || use_ref_mirror){

					Dict retvals;
					Dict retvals_mirror;
					if (use_ref_mirror == true)  retvals_mirror = Crosrng_sm_psi(crefim[iref], cimage, numr, psi - 90.0f, 1, psi_max);
					if (use_ref == true)         retvals        = Crosrng_sm_psi(crefim[iref], cimage, numr, psi - 90.0f, 0, psi_max);
					double qn = jneginf;
					if (use_ref) qn = retvals["qn"];
					double qm = jneginf;
					if (use_ref_mirror) qm = retvals_mirror["qn"];

					if(qn >= peak || qm >= peak) {
						sx = -ix;
						sy = -iy;
						nref = iref;
						if (qn >= qm){
							ang = ang_n(retvals["tot"], mode, numr[numr.size()-1]);
							mirror = 0;
							peak = static_cast<float>(qn);
						} else {
							ang = ang_n(retvals_mirror["tot"], mode, numr[numr.size()-1]);
							mirror = 1;
							peak = static_cast<float>(qm);
						}
					}
				}
			}
			delete cimage; cimage = 0;
		}
	}

	float co, so, sxs, sys;
	co = static_cast<float>(  cos(ang*qv) );
	so = static_cast<float>( -sin(ang*qv) );
	sxs = sx*co - sy*so;
	sys = sx*so + sy*co;
	vector<float> res;
	res.push_back(ang);
	res.push_back(sxs);
	res.push_back(sys);
	res.push_back(static_cast<float>(mirror));
	res.push_back(static_cast<float>(nref));
	res.push_back(peak);
	return res;
}

//  HELICON

vector<float> Util::multiref_polar_ali_helicon_90_local(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, float ant, float psi_max, string mode,
                vector<int>numr, float cnx, float cny, int ynumber, float yrnglocal) {

	size_t crefim_len = crefim.size();
	const float qv = static_cast<float>( pi/180.0 );

	Transform * t = image->get_attr("xform.projection");
	Dict d = t->get_params("spider");
	if(t) {delete t; t=0;}
	float phi   = d["phi"];
	float theta = d["theta"];
	float psi   = d["psi"];
	float imn1 = sin(theta*qv)*cos(phi*qv);
	float imn2 = sin(theta*qv)*sin(phi*qv);
	float imn3 = cos(theta*qv);
	vector<float> n1(crefim_len);
	vector<float> n2(crefim_len);
	vector<float> n3(crefim_len);
	int   iref, nref=-1, mirror=0;
	float iy, ix, sx=0, sy=0;
	float peak = -1.0E23f;
	float ang  = 0.0f;
//	int   kx   = int(2*xrng/step+0.5)/2;
	int lkx = int(xrng[0]/step);
	int rkx = int(xrng[1]/step);

	for ( iref = 0; iref < (int)crefim_len; iref++) {
		n1[iref] = crefim[iref]->get_attr("n1");
		n2[iref] = crefim[iref]->get_attr("n2");
		n3[iref] = crefim[iref]->get_attr("n3");
	}

	float stepy;
	int lky, rky;

	if (ynumber == 0) {
		lky = rky = 0;
	} else {
		if (ynumber > 0) stepy = (yrng[0]+yrng[1])/ynumber;
		else if (ynumber == -1) stepy = step;

		if (yrnglocal >= 0.0) {
//		   	ky = int(yrnglocal/stepy);
			lky = min(int(yrng[0]/stepy), int(yrnglocal/stepy));
			rky = min(int(yrng[1]/stepy), int(yrnglocal/stepy));
		}
		else { // search range is not restricted
			lky = int(yrng[0]/stepy);
			rky = int(yrng[1]/stepy);
			// if (ynumber > 0) {
// 				ky = int(ynumber/2);
// 			} else {
// 				ky = int(2*yrng/stepy+0.5)/2;
// 			}
		}
	}
	int nring = numr.size()/3;
	for (int i = -lky; i <= rky; i++) {
		iy = i * stepy ;
		for (int j = -lkx; j <= rkx; j++) {
			ix = j*step ;
			EMData* cimage = Polar2Dm(image, cnx+ix, cny+iy, numr, mode);

			Normalize_ring( cimage, numr );

			Frngs(cimage, numr);
			//  compare with all reference images
			// for iref in xrange(len(crefim)):
			for ( iref = 0; iref < (int)crefim_len; iref++) {
				//if(fabs(n1[iref]*imn1 + n2[iref]*imn2 + n3[iref]*imn3)>=ant) {
				if( (n1[iref]*imn1 + n2[iref]*imn2 + n3[iref]*imn3) >= ant ) {
					Dict retvals; 
					retvals = Crosrng_sm_psi(crefim[iref], cimage, numr,  psi - 90.0f, 0, psi_max);

					double qn = retvals["qn"];
					if( qn >= peak) {
						sx = -ix;
						sy = -iy;
						nref = iref;
						ang = ang_n(retvals["tot"], mode, numr[numr.size()-1]);
						peak = static_cast<float>(qn);
						mirror = 0;
					}
				}
			}  
			delete cimage; cimage = 0;
		}
	}

	float co, so, sxs, sys;
	co = static_cast<float>(  cos(ang*qv) );
	so = static_cast<float>( -sin(ang*qv) );
	sxs = sx*co - sy*so;
	sys = sx*so + sy*co;
	vector<float> res;
	res.push_back(ang);
	res.push_back(sxs);
	res.push_back(sys);
	res.push_back(static_cast<float>(mirror));
	res.push_back(static_cast<float>(nref));
	res.push_back(peak);
	return res;
}



void  Util::multiref_peaks_ali2d(EMData* image, EMData* crefim,
			vector<float> xrng, vector<float> yrng, float step, string mode,
			vector< int >numr, float cnx, float cny,
			EMData *peaks, EMData *peakm) {

	int   maxrin = numr[numr.size()-1];

// 	int   ky = int(2*yrng/step+0.5)/2;
// 	int   kx = int(2*xrng/step+0.5)/2;
	int lkx = int(xrng[0]/step);
	int rkx = int(xrng[1]/step);
	int lky = int(yrng[0]/step);
	int rky = int(yrng[1]/step);
	
	peaks->set_size(maxrin, lkx+rkx+3, lky+rky+3);
	float *p_ccf1ds = peaks->get_data();

	peakm->set_size(maxrin, lkx+rkx+3, lky+rky+3);
	float *p_ccf1dm = peakm->get_data();

	for ( int i = 0; i<maxrin*(lkx+rkx+3)*(lky+rky+3); i++) {
		p_ccf1ds[i] = -1.e20f;
		p_ccf1dm[i] = -1.e20f;
	}

	for (int i = -lky; i <= rky; i++) {
		float iy = i * step ;
		for (int j = -lkx; j <= rkx; j++) {
			float ix = j*step ;
			EMData* cimage = Polar2Dm(image, cnx+ix, cny+iy, numr, mode);
			Frngs(cimage, numr);
			Crosrng_msg_vec(crefim, cimage, numr,
			  p_ccf1ds+(j+lkx+1+((i+lky+1)*(lkx+rkx+3)))*maxrin,
			  p_ccf1dm+(j+lkx+1+((i+lky+1)*(lkx+rkx+3)))*maxrin);
			delete cimage; cimage = 0;
		}
	}
	return;
}

void  Util::multiref_peaks_compress_ali2d(EMData* image, EMData* crefim, vector<float> xrng, vector<float> yrng,
     float step, string mode, vector<int>numr, float cnx, float cny, EMData *peaks, EMData *peakm,
     EMData *peaks_compress, EMData *peakm_compress) {

	int   maxrin = numr[numr.size()-1];

// 	int   ky = int(2*yrng/step+0.5)/2;
// 	int   kx = int(2*xrng/step+0.5)/2;
	int lkx = int(xrng[0]/step);
	int rkx = int(xrng[1]/step);
	int lky = int(yrng[0]/step);
	int rky = int(yrng[1]/step);
	
	peaks->set_size(maxrin, lkx+rkx+3, lky+rky+3);
	float *p_ccf1ds = peaks->get_data();

	peakm->set_size(maxrin, lkx+rkx+3, lky+rky+3);
	float *p_ccf1dm = peakm->get_data();

	peaks_compress->set_size(maxrin, 1, 1);
	float *p_ccf1ds_compress = peaks_compress->get_data();

	peakm_compress->set_size(maxrin, 1, 1);
	float *p_ccf1dm_compress = peakm_compress->get_data();

	for ( int i = 0; i<maxrin*(lkx+rkx+3)*(lky+rky+3); i++) {
		p_ccf1ds[i] = -1.e20f;
		p_ccf1dm[i] = -1.e20f;
	}
	for (int i = -lky; i <= rky; i++) {
		float iy = i * step ;
		for (int j = -lkx; j <= rkx; j++) {
			float ix = j*step ;
			EMData* cimage = Polar2Dm(image, cnx+ix, cny+iy, numr, mode);
			Frngs(cimage, numr);
			Crosrng_msg_vec(crefim, cimage, numr,
			  p_ccf1ds+(j+lkx+1+((i+lky+1)*(lkx+rkx+3)))*maxrin,
			  p_ccf1dm+(j+lkx+1+((i+lky+1)*(lkx+rkx+3)))*maxrin);
			delete cimage; cimage = 0;
		}
	}
	for (int x=0; x<maxrin; x++) {
		float maxs = -1.0e22f;
		float maxm = -1.0e22f;
		for (int i=1; i<=lky+rky+1; i++) {
			for (int j=1; j<=lkx+rkx+1; j++) {
				if (p_ccf1ds[(i*(lkx+rkx+3)+j)*maxrin+x] > maxs) maxs = p_ccf1ds[(i*(lkx+rkx+3)+j)*maxrin+x];
				if (p_ccf1dm[(i*(lkx+rkx+3)+j)*maxrin+x] > maxm) maxm = p_ccf1dm[(i*(lkx+rkx+3)+j)*maxrin+x];
			}
		}
		p_ccf1ds_compress[x] = maxs;
		p_ccf1dm_compress[x] = maxm;
	}
	return;
}

struct ccf_point
{
    float value;
    int i;
    int j;
    int k;
    int mirror;
};


struct ccf_value
{
    bool operator()( const ccf_point& a, const ccf_point& b )
    {
        return a.value > b.value;
    }
};


vector<float>  Util::ali2d_ccf_list(EMData* image, EMData* crefim,
			vector<float> xrng, vector<float> yrng, float step, string mode,
			vector< int >numr, float cnx, float cny, double T) {

	int   maxrin = numr[numr.size()-1];

// 	int   ky = int(2*yrng/step+0.5)/2;
// 	int   kx = int(2*xrng/step+0.5)/2;
	int lkx = int(xrng[0]/step);
	int rkx = int(xrng[1]/step);
	int lky = int(yrng[0]/step);
	int rky = int(yrng[1]/step);
	
	float *p_ccf1ds = (float *)malloc(maxrin*sizeof(float));
	float *p_ccf1dm = (float *)malloc(maxrin*sizeof(float));
	int vol = maxrin*(lkx+rkx+1)*(lky+rky+1);
	vector<ccf_point> ccf(2*vol);
	ccf_point temp;

	int index = 0;
	for (int i = -lky; i <= rky; i++) {
		float iy = i * step ;
		for (int j = -lkx; j <= rkx; j++) {
			float ix = j*step ;
			EMData* cimage = Polar2Dm(image, cnx+ix, cny+iy, numr, mode);
			Frngs(cimage, numr);
			Crosrng_msg_vec(crefim, cimage, numr, p_ccf1ds, p_ccf1dm);
			for (int k=0; k<maxrin; k++) {
				temp.value = p_ccf1ds[k];
				temp.i = k;
				temp.j = j;
				temp.k = i;
				temp.mirror = 0;
				ccf[index] = temp;
				index++;
				temp.value = p_ccf1dm[k];
				temp.mirror = 1;
				ccf[index] = temp;
				index++;
			}
			delete cimage; cimage = 0;
		}
	}

	delete p_ccf1ds;
	delete p_ccf1dm;
	std::sort(ccf.begin(), ccf.end(), ccf_value());

	double qt = (double)ccf[0].value;
	vector <double> p(2*vol), cp(2*vol);

	double sump = 0.0;
	for (int i=0; i<2*vol; i++) {
		p[i] = pow(double(ccf[i].value)/qt, 1.0/T);
		sump += p[i];
	}
	for (int i=0; i<2*vol; i++) {
		p[i] /= sump;
	}
	for (int i=1; i<2*vol; i++) {
		p[i] += p[i-1];
	}
	p[2*vol-1] = 2.0;

	float t = get_frand(0.0f, 1.0f);
	int select = 0;
	while (p[select] < t)	select += 1;

	vector<float> a(6);
	a[0] = ccf[select].value;
	a[1] = (float)ccf[select].i;
	a[2] = (float)ccf[select].j;
	a[3] = (float)ccf[select].k;
	a[4] = (float)ccf[select].mirror;
	a[5] = (float)select;
	return a;
}

vector<float>  Util::ali2d_ccf_list_snake(EMData* image, EMData* crefim, vector<float> wr, 
			float xrng, float yrng, float step, string mode,
			vector< int >numr, float cnx, float cny, double T) {

	int   maxrin = numr[numr.size()-1];
	
	int   kx = int(2*xrng/step+0.5)/2;
	
	float *p_ccf1ds = (float *)malloc(maxrin*sizeof(float));
	int vol = maxrin*(2*kx+1);
	vector<float> ccf(vol);
	int index;

	int nring = numr.size()/3;
	int inr = numr[(nring-1)*3];
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	//here the centers of the image cny and cnx use Fortran convention which means the index of the image array starts from 1 to nx (ny). 02-23-2015
	for (int j = -kx; j <= kx; j++) {
		float ix = j*step ;
		if ( inr + (int)(cnx+ix) <= nx - 1  &&  -inr + (int)(cnx+ix) >= 1) {
			EMData* cimage = Polar2Dm(image, cnx+ix, cny, numr, mode);
			Frngs(cimage, numr);
			//Applyws(cimage, numr, wr);
			Crosrng_msg_vec_snake(crefim, cimage, numr, p_ccf1ds);
			for (int k=0; k<maxrin; k++) {
				index = (j+kx)*maxrin+k;
				ccf[index] = p_ccf1ds[k];
			}
			delete cimage; cimage = 0;
		}
				
	}
	

	delete p_ccf1ds;
	return ccf;
	// std::sort(ccf.begin(), ccf.end(), ccf_value());
// 
// 	double qt = (double)ccf[0].value;
// 	vector <double> p(vol), cp(vol);
// 
// 	double sump = 0.0;
// 	for (int i=0; i<vol; i++) {
// 		p[i] = pow(double(ccf[i].value)/qt, 1.0/T);
// 		sump += p[i];
// 	}
// 	for (int i=0; i<vol; i++) {
// 		p[i] /= sump;
// 	}
// 	for (int i=1; i<vol; i++) {
// 		p[i] += p[i-1];
// 	}
// 	p[vol-1] = 2.0;
// 
// 	float t = get_frand(0.0f, 1.0f);
// 	int select = 0;
// 	while (p[select] < t)	select += 1;
// 
// 	vector<float> a(6);
// 	a[0] = ccf[select].value;
// 	a[1] = (float)ccf[select].i;
// 	a[2] = (float)ccf[select].j;
// 	a[3] = (float)ccf[select].k;
// 	a[4] = (float)ccf[select].mirror;
// 	a[5] = (float)select;
//	return a;
}




/*
void Util::multiref_peaks_ali(EMData* image, const vector< EMData* >& crefim,
			float xrng, float yrng, float step, string mode,
			vector< int >numr, float cnx, float cny,
			EMData *peaks, EMData *peakm, int nphi, int ntheta) {

// formerly known as apmq
    // Determine shift and rotation between image and many reference
    // images (crefim, weights have to be applied) quadratic
    // interpolation


    // Manually extract.
*//*    vector< EMAN::EMData* > crefim;
    std::size_t crefim_len = PyObject_Length(crefim_list.ptr());
    crefim.reserve(crefim_len);

    for(std::size_t i=0;i<crefim_len;i++) {
        boost::python::extract<EMAN::EMData*> proxy(crefim_list[i]);
        crefim.push_back(proxy());
    }
*/
/*
	int   maxrin = numr[numr.size()-1];

	size_t crefim_len = crefim.size();

	int   iref;
	int   ky = int(2*yrng/step+0.5)/2;
	int   kx = int(2*xrng/step+0.5)/2;
	int   tkx = 2*kx+3;
	int   tky = 2*ky+3;

	peaks->set_size(maxrin, nphi, ntheta, tkx, tky);
	float *p_ccf1ds = peaks->get_data();

//#define img_ptr(i,j,k,l)  img_ptr[i+(j+(k+(l*nz))*ny))*nx]
//#define img_ptr(i,j,k,l)  img_ptr[i+(j+(k*ny))*nx]
	peakm->set_size(maxrin, nphi, ntheta, tkx, tky);
	float *p_ccf1dm = peakm->get_data();

	for ( int i = 0; i<maxrin*(int)crefim_len*(2*kx+3)*(2*ky+3); i++) {
		p_ccf1ds[i] = -1.e20f;
		p_ccf1dm[i] = -1.e20f;
	}

	float  iy, ix;
	int nring = numr.size()/3;
	int inr = numr[(nring-1)*3];
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	for (int i = -ky; i <= ky; i++) {
		iy = i * step ;
		//here the centers of the image cny and cnx use Fortran convention which means the index of the image array starts from 1 to nx (ny). 02-23-2015
		if ( inr + (int)(cny+iy) <= ny - 1  &&  -inr + (int)(cny+iy) >= 1) {
			for (int j = -kx; j <= kx; j++) {
				ix = j*step ;
				if ( inr + (int)(cnx+ix) <= nx - 1  &&  -inr + (int)(cnx+ix) >= 1) {
					EMData* cimage = Polar2Dm(image, cnx+ix, cny+iy, numr, mode);
					Frngs(cimage, numr);
					//  compare with all reference images
					// for iref in xrange(len(crefim)):
					for ( iref = 0; iref < (int)crefim_len; iref++) {
						Crosrng_msg_vec(crefim[iref], cimage, numr,
							p_ccf1ds+(iref + (j+kx+1+((i+ky+1)*tkx))*(int)crefim_len )*maxrin,
							p_ccf1dm+(iref + (j+kx+1+((i+ky+1)*tkx))*(int)crefim_len )*maxrin);
					}
					delete cimage; cimage = 0;
				}
			}		
		}
	}
	return;
}
*/
vector<float> Util::twoD_fine_ali(EMData* image, EMData *refim, EMData* mask, float ang, float sxs, float sys) {

	EMData *rot;

	const int nmax=3, mmax=3;
	char task[60], csave[60];
	long int lsave[4];
	long int n, m, iprint, nbd[nmax], iwa[3*nmax], isave[44];
	double f, f1, f2, f3, factr, pgtol, x[nmax], l[nmax], u[nmax], g[nmax], dsave[29], wa[2*mmax*nmax+4*nmax+12*mmax*mmax+12*mmax];
	long int SIXTY=60;

	//     We wish to have no output.
	iprint = -1;

	//c     We specify the tolerances in the stopping criteria.
	factr=1.0e1;
	pgtol=1.0e-5;

	//     We specify the dimension n of the sample problem and the number
	//        m of limited memory corrections stored.  (n and m should not
	//        exceed the limits nmax and mmax respectively.)
	n=3;
	m=3;

	//     We now provide nbd which defines the bounds on the variables:
	//                    l   specifies the lower bounds,
	//                    u   specifies the upper bounds.
 	//		      x   specifies the initial guess
	x[0] = ang; nbd[0] = 2; l[0] = ang-2.0; u[0] = ang+2.0;
	x[1] = sxs; nbd[1] = 2; l[1] = sxs-1.5; u[1] = sxs+1.5;
	x[2] = sys; nbd[2] = 2; l[2] = sys-1.5; u[2] = sys+1.5;


	//     We start the iteration by initializing task.
	// (**MUST clear remaining chars in task with spaces (else crash)!**)
	strcpy(task,"START");
	for (int i=5;i<60;i++)	task[i]=' ';

	//     This is the call to the L-BFGS-B code.
	// (* call the L-BFGS-B routine with task='START' once before loop *)
	setulb_(&n,&m,x,l,u,nbd,&f,g,&factr,&pgtol,wa,iwa,task,&iprint,csave,lsave,isave,dsave,SIXTY,SIXTY);
	//int step = 1;

 	// (* while routine returns "FG" or "NEW_X" in task, keep calling it *)
	while (strncmp(task,"FG",2)==0 || strncmp(task,"NEW_X",5)==0) {

		if (strncmp(task,"FG",2)==0) {
	      	//   the minimization routine has returned to request the
		//   function f and gradient g values at the current x

		//        Compute function value f for the sample problem.
		rot = new EMData();
		rot = image->rot_scale_trans2D((float)x[0], (float)x[1], (float)x[2], 1.0f);
		f = rot->cmp("sqeuclidean", refim, Dict("mask", mask));
		//f = -f;
		delete rot;

	      	//        Compute gradient g for the sample problem.
		float dt = 1.0e-3f;
		rot = new EMData();
		rot = image->rot_scale_trans2D((float)x[0]+dt, (float)x[1], (float)x[2], 1.0f);
		f1 = rot->cmp("sqeuclidean", refim, Dict("mask", mask));
		//f1 = -f1;
		g[0] = (f1-f)/dt;
		delete rot;

		dt = 1.0e-2f;
		rot = new EMData();
		rot = image->rot_scale_trans2D((float)x[0], (float)x[1]+dt, (float)x[2], 1.0f);
		f2 = rot->cmp("sqeuclidean", refim, Dict("mask", mask));
		//f2 = -f2;
		g[1] = (f2-f)/dt;
		delete rot;

		rot = new EMData();
		rot = image->rot_scale_trans2D((float)x[0], (float)x[1], (float)x[2]+dt, 1.0f);
		f3 = rot->cmp("sqeuclidean", refim, Dict("mask", mask));
		//f3 = -f3;
		g[2] = (f3-f)/dt;
		delete rot;
   		}

		//c          go back to the minimization routine.
		setulb_(&n,&m,x,l,u,nbd,&f,g,&factr,&pgtol,wa,iwa,task,&iprint,csave,lsave,isave,dsave,SIXTY,SIXTY);
		//step++;
  	}

	//printf("Total step is %d\n", step);
	vector<float> res;
	res.push_back(static_cast<float>(x[0]));
	res.push_back(static_cast<float>(x[1]));
	res.push_back(static_cast<float>(x[2]));
	//res.push_back(step);
	return res;
}

vector<float> Util::twoD_fine_ali_G(EMData* image, EMData *refim, EMData* mask, Util::KaiserBessel& kb, float ang, float sxs, float sys) {

	EMData *rot;

	const int nmax=3, mmax=3;
	char task[60], csave[60];
	long int lsave[4];
	long int n, m, iprint, nbd[nmax], iwa[3*nmax], isave[44];
	double f, f1, f2, f3, factr, pgtol, x[nmax], l[nmax], u[nmax], g[nmax], dsave[29], wa[2*mmax*nmax+4*nmax+12*mmax*mmax+12*mmax];
	long int SIXTY=60;

	//     We wish to have no output.
	iprint = -1;

	//c     We specify the tolerances in the stopping criteria.
	factr=1.0e1;
	pgtol=1.0e-5;

	//     We specify the dimension n of the sample problem and the number
	//        m of limited memory corrections stored.  (n and m should not
	//        exceed the limits nmax and mmax respectively.)
	n=3;
	m=3;

	//     We now provide nbd which defines the bounds on the variables:
	//                    l   specifies the lower bounds,
	//                    u   specifies the upper bounds.
 	//		      x   specifies the initial guess
	x[0] = ang; nbd[0] = 2; l[0] = ang-2.0; u[0] = ang+2.0;
	x[1] = sxs; nbd[1] = 2; l[1] = sxs-1.5; u[1] = sxs+1.5;
	x[2] = sys; nbd[2] = 2; l[2] = sys-1.5; u[2] = sys+1.5;


	//     We start the iteration by initializing task.
	// (**MUST clear remaining chars in task with spaces (else crash)!**)
	strcpy(task,"START");
	for (int i=5;i<60;i++)	task[i]=' ';

	//     This is the call to the L-BFGS-B code.
	// (* call the L-BFGS-B routine with task='START' once before loop *)
	setulb_(&n,&m,x,l,u,nbd,&f,g,&factr,&pgtol,wa,iwa,task,&iprint,csave,lsave,isave,dsave,SIXTY,SIXTY);
	//int step = 1;

 	// (* while routine returns "FG" or "NEW_X" in task, keep calling it *)
	while (strncmp(task,"FG",2)==0 || strncmp(task,"NEW_X",5)==0) {

		if (strncmp(task,"FG",2)==0) {
	      	//   the minimization routine has returned to request the
		//   function f and gradient g values at the current x

		//        Compute function value f for the sample problem.
		rot = new EMData();
		rot = image->rot_scale_conv7((float)(x[0]*pi/180), (float)x[1], (float)x[2], kb, 1.0f);
		f = rot->cmp("sqeuclidean", refim, Dict("mask", mask));
		//f = -f;
		delete rot;

	      	//        Compute gradient g for the sample problem.
		float dt = 1.0e-3f;
		rot = new EMData();
		rot = image->rot_scale_conv7((float)((x[0]+dt)*pi/180), (float)x[1], (float)x[2], kb, 1.0f);
		f1 = rot->cmp("sqeuclidean", refim, Dict("mask", mask));
		//f1 = -f1;
		g[0] = (f1-f)/dt;
		delete rot;

		rot = new EMData();
		rot = image->rot_scale_conv7((float)(x[0]*pi/180), (float)x[1]+dt, (float)x[2], kb, 1.0);
		f2 = rot->cmp("sqeuclidean", refim, Dict("mask", mask));
		//f2 = -f2;
		g[1] = (f2-f)/dt;
		delete rot;

		rot = new EMData();
		rot = image->rot_scale_conv7((float)(x[0]*pi/180), (float)x[1], (float)x[2]+dt, kb, 1.0f);
		f3 = rot->cmp("sqeuclidean", refim, Dict("mask", mask));
		//f3 = -f3;
		g[2] = (f3-f)/dt;
		delete rot;
   		}

		//c          go back to the minimization routine.
		setulb_(&n,&m,x,l,u,nbd,&f,g,&factr,&pgtol,wa,iwa,task,&iprint,csave,lsave,isave,dsave,SIXTY,SIXTY);
		//step++;
  	}

	//printf("Total step is %d\n", step);
	vector<float> res;
	res.push_back(static_cast<float>(x[0]));
	res.push_back(static_cast<float>(x[1]));
	res.push_back(static_cast<float>(x[2]));
	//res.push_back(step);
	return res;
}

vector<float> Util::twoD_to_3D_ali(EMData* volft, Util::KaiserBessel& kb, EMData *refim, EMData* mask, float phi, float theta, float psi, float sxs, float sys) {

	EMData *proj, *proj2;

	const int nmax=5, mmax=5;
	char task[60], csave[60];
	long int lsave[4];
	long int n, m, iprint, nbd[nmax], iwa[3*nmax], isave[44];
	double f, ft, factr, pgtol, x[nmax], l[nmax], u[nmax], g[nmax], dsave[29], wa[2*mmax*nmax+4*nmax+12*mmax*mmax+12*mmax];
	long int SIXTY=60;

	//     We wish to have no output.
	iprint = -1;

	//c     We specify the tolerances in the stopping criteria.
	factr=1.0e1;
	pgtol=1.0e-5;

	//     We specify the dimension n of the sample problem and the number
	//        m of limited memory corrections stored.  (n and m should not
	//        exceed the limits nmax and mmax respectively.)
	n=5;
	m=5;

	//     We now provide nbd which defines the bounds on the variables:
	//                    l   specifies the lower bounds,
	//                    u   specifies the upper bounds.
 	//		      x   specifies the initial guess
	x[0] = phi; 	nbd[0] = 2; 	l[0] = phi-2.0; 	u[0] = phi+2.0;
	x[1] = theta; 	nbd[1] = 2; 	l[1] = theta-2.0;	u[1] = theta+2.0;
	x[2] = psi; 	nbd[2] = 2; 	l[2] = psi-2.0; 	u[2] = psi+2.0;
	x[3] = sxs; 	nbd[3] = 2; 	l[3] = sxs-2.0; 	u[3] = sxs+2.0;
	x[4] = sys; 	nbd[4] = 2; 	l[4] = sys-2.0; 	u[4] = sys+2.0;


	//     We start the iteration by initializing task.
	// (**MUST clear remaining chars in task with spaces (else crash)!**)
	strcpy(task,"START");
	for (int i=5;i<60;i++)	task[i]=' ';

	//     This is the call to the L-BFGS-B code.
	// (* call the L-BFGS-B routine with task='START' once before loop *)
	setulb_(&n,&m,x,l,u,nbd,&f,g,&factr,&pgtol,wa,iwa,task,&iprint,csave,lsave,isave,dsave,SIXTY,SIXTY);
	int step = 1;

 	// (* while routine returns "FG" or "NEW_X" in task, keep calling it *)
	while (strncmp(task,"FG",2)==0 || strncmp(task,"NEW_X",5)==0) {

		if (strncmp(task,"FG",2)==0) {
	      	//   the minimization routine has returned to request the
		//   function f and gradient g values at the current x

		//        Compute function value f for the sample problem.
		proj = new EMData();
		proj2 = new EMData();
		proj = volft->extract_plane(Transform(Dict("type", "SPIDER", "phi", (float)x[0], "theta", (float)x[1], "psi", (float)x[2])), kb);
		proj->fft_shuffle();
		proj->center_origin_fft();
		proj->process_inplace("filter.shift", Dict("x_shift", x[3], "y_shift", x[4], "z_shift", 0.0f));
		proj->do_ift_inplace();
		int M = proj->get_ysize()/2;
		proj2 = proj->window_center(M);
		f = proj2->cmp("sqeuclidean", refim, Dict("mask", mask));
		//f = -f;
		delete proj;
		delete proj2;

	      	//        Compute gradient g for the sample problem.
		float dt = 1.0e-3f;
		proj = new EMData();
		proj2 = new EMData();
		proj = volft->extract_plane(Transform(Dict("type", "SPIDER", "phi", (float)x[0]+dt, "theta", (float)x[1], "psi", (float)x[2])), kb);
		proj->fft_shuffle();
		proj->center_origin_fft();
		proj->process_inplace("filter.shift", Dict("x_shift", x[3], "y_shift", x[4], "z_shift", 0.0f));
		proj->do_ift_inplace();
		proj2 = proj->window_center(M);
		ft = proj2->cmp("sqeuclidean", refim, Dict("mask", mask));
		//ft = -ft;
		delete proj;
		delete proj2;
		g[0] = (ft-f)/dt;

		proj = new EMData();
		proj2 = new EMData();
		proj = volft->extract_plane(Transform(Dict("type", "SPIDER", "phi", (float)x[0], "theta", (float)x[1]+dt, "psi", (float)x[2])), kb);
		proj->fft_shuffle();
		proj->center_origin_fft();
		proj->process_inplace("filter.shift", Dict("x_shift", x[3], "y_shift", x[4], "z_shift", 0.0f));
		proj->do_ift_inplace();
		proj2 = proj->window_center(M);
		ft = proj2->cmp("sqeuclidean", refim, Dict("mask", mask));
		//ft = -ft;
		delete proj;
		delete proj2;
		g[1] = (ft-f)/dt;

		proj = new EMData();
		proj2 = new EMData();
		proj = volft->extract_plane(Transform(Dict("type", "SPIDER", "phi", (float)x[0], "theta", (float)x[1], "psi", (float)x[2]+dt)), kb);
		proj->fft_shuffle();
		proj->center_origin_fft();
		proj->process_inplace("filter.shift", Dict("x_shift", x[3], "y_shift", x[4], "z_shift", 0.0f));
		proj->do_ift_inplace();
		proj2 = proj->window_center(M);
		ft = proj2->cmp("sqeuclidean", refim, Dict("mask", mask));
		//ft = -ft;
		delete proj;
		delete proj2;
		g[2] = (ft-f)/dt;

		proj = new EMData();
		proj2 = new EMData();
		proj = volft->extract_plane(Transform(Dict("type", "SPIDER", "phi", (float)x[0], "theta", (float)x[1], "psi", (float)x[2])), kb);
		proj->fft_shuffle();
		proj->center_origin_fft();
		proj->process_inplace("filter.shift", Dict("x_shift", x[3]+dt, "y_shift", x[4], "z_shift", 0.0f));
		proj->do_ift_inplace();
		proj2 = proj->window_center(M);
		ft = proj2->cmp("sqeuclidean", refim, Dict("mask", mask));
		//ft = -ft;
		delete proj;
		delete proj2;
		g[3] = (ft-f)/dt;

		proj = new EMData();
		proj2 = new EMData();
		proj = volft->extract_plane(Transform(Dict("type", "SPIDER", "phi", (float)x[0], "theta", (float)x[1], "psi", (float)x[2])), kb);
		proj->fft_shuffle();
		proj->center_origin_fft();
		proj->process_inplace("filter.shift", Dict("x_shift", x[3], "y_shift", x[4]+dt, "z_shift", 0.0f));
		proj->do_ift_inplace();
		proj2 = proj->window_center(M);
		ft = proj2->cmp("sqeuclidean", refim, Dict("mask", mask));
		//ft = -ft;
		delete proj;
		delete proj2;
		g[4] = (ft-f)/dt;
   		}

		//c          go back to the minimization routine.
		setulb_(&n,&m,x,l,u,nbd,&f,g,&factr,&pgtol,wa,iwa,task,&iprint,csave,lsave,isave,dsave,SIXTY,SIXTY);
		step++;
  	}

	//printf("Total step is %d\n", step);
	vector<float> res;
	res.push_back(static_cast<float>(x[0]));
	res.push_back(static_cast<float>(x[1]));
	res.push_back(static_cast<float>(x[2]));
	res.push_back(static_cast<float>(x[3]));
	res.push_back(static_cast<float>(x[4]));
	//res.push_back(step);
	return res;
}


vector<float> Util::twoD_fine_ali_SD(EMData* image, EMData *refim, EMData* mask, float ang, float sxs, float sys) {

	double  x[4];
	int n;
	int l = 3;
	int m = 200;
	double e = 1e-9;
	double step = 0.01;
	float (*my_func)(EMData* , EMData* , EMData* , float , float , float) = ccc_images;

	x[1] = ang;
	x[2] = sxs;
	x[3] = sys;

	Steepda(x, step, e, l, m, &n, my_func, image, refim, mask);   // Call steepest descent optimization subroutine
	//printf("Took %d steps\n", n);

	vector<float> res;
	res.push_back(static_cast<float>(x[1]));
	res.push_back(static_cast<float>(x[2]));
	res.push_back(static_cast<float>(x[3]));
	res.push_back(static_cast<float>(n));
	return res;
}

/* Parameters:
 * args - parameters of (L-1) G transformation (G_1, G_2, ..., G_(L-1)), saved as alpha_1, sx_1, sy_1, alpha_2, sx_2, sy_2, ... (transformation G_L is always set to I)
 * all_ali_params - parameters of (L*N) T transformation, saved as flat list: <parameters of N images for G_1 alignment>, <parameters of N images for G_2 alignment>, ...
 *                  where image's parameters consists of four values: alpha, sx, sy, mirror
 * d - particle diameter (2R)
 */
vector<float> Util::multi_align_error(vector<float> args, vector<float> all_ali_params, int d) {
	
	const int nmax=args.size(), mmax=nmax;
	char task[60], csave[60];
	long int lsave[4];
	long int n, m, iprint, isave[44];
	long int* nbd = new long int[nmax];
	long int* iwa = new long int[3*nmax];
	double f, factr, pgtol;
	double* x = new double[nmax];
	double* l = new double[nmax];
	double* u = new double[nmax];
	double* g = new double[nmax];
	double dsave[29];
	double* wa = new double[2*mmax*nmax+4*nmax+12*mmax*mmax+12*mmax];
	long int SIXTY=60;

	int num_ali = nmax/3+1;
	int nima = all_ali_params.size()/(num_ali*4);
	
	//     We wish to have no output.
	iprint = -1;

	//c     We specify the tolerances in the stopping criteria.
	factr=1.0e1;
	pgtol=1.0e-9;

	//     We specify the dimension n of the sample problem and the number
	//        m of limited memory corrections stored.  (n and m should not
	//        exceed the limits nmax and mmax respectively.)
	n=nmax;
	m=mmax;

	//     We now provide nbd which defines the bounds on the variables:
	//                    l   specifies the lower bounds,
	//                    u   specifies the upper bounds.
 	//		      x   specifies the initial guess
	for (int i=0; i<nmax; i++) {
		x[i] = args[i]; 
		nbd[i] = 0;
	}

	//     We start the iteration by initializing task.
	// (**MUST clear remaining chars in task with spaces (else crash)!**)
	strcpy(task,"START");
	for (int i=5;i<60;i++)	task[i]=' ';

	//     This is the call to the L-BFGS-B code.
	// (* call the L-BFGS-B routine with task='START' once before loop *)
	setulb_(&n,&m,x,l,u,nbd,&f,g,&factr,&pgtol,wa,iwa,task,&iprint,csave,lsave,isave,dsave,SIXTY,SIXTY);
	int step = 1;

 	// (* while routine returns "FG" or "NEW_X" in task, keep calling it *)
	while (strncmp(task,"FG",2)==0 || strncmp(task,"NEW_X",5)==0) {

		if (strncmp(task,"FG",2)==0) {
	      	//   the minimization routine has returned to request the
		//   function f and gradient g values at the current x

		//        Compute function value f for the sample problem.
		f = multi_align_error_func(x, all_ali_params, nima, num_ali, d);

	      	//        Compute gradient g for the sample problem.
		multi_align_error_dfunc(x, all_ali_params, nima, num_ali, g, d);

		}
		//c          go back to the minimization routine.
		setulb_(&n,&m,x,l,u,nbd,&f,g,&factr,&pgtol,wa,iwa,task,&iprint,csave,lsave,isave,dsave,SIXTY,SIXTY);
		step++;
  	}

	//printf("Total step is %d\n", step);
	vector<float> res;
	for (int i=0; i<nmax; i++) res.push_back(static_cast<float>(x[i]));
	res.push_back(static_cast<float>(f));

	delete[] nbd;
	delete[] iwa;
	delete[] x;
	delete[] l;
	delete[] u;
	delete[] g;
	delete[] wa;

	return res;

}

double Util::multi_align_error_func(double* x, vector<float> all_ali_params, int nima, int num_ali, int d)
{
	vector<double> sqr_pixel_error = multi_align_error_func2(x, all_ali_params, nima, num_ali, d);
	double sum_sqr_pixel_error = 0.0;
	for (int i=0; i<nima; i++) {
		sum_sqr_pixel_error += sqr_pixel_error[i];
	}
	return ( sum_sqr_pixel_error / nima );
}


vector<double> Util::multi_align_error_func2(double* x, vector<float> ali_params, int nima, int num_ali, int d)
{
	double* args = new double[num_ali*3];
	for (int i=0; i<3*num_ali-3; i++)   args[i] = x[i];
	args[3*num_ali-3] = 0.0;
	args[3*num_ali-2] = 0.0;
	args[3*num_ali-1] = 0.0;
	double* cosa = new double[num_ali];
	double* sina = new double[num_ali];
	for (int i=0; i<num_ali; i++) {
		cosa[i] = cos(args[i*3]*M_PI/180.0);
		sina[i] = sin(args[i*3]*M_PI/180.0);
	}
	double* sx = new double[num_ali];
	double* sy = new double[num_ali];
	
	vector<double> sqr_pixel_error(nima);

	for (int i=0; i<nima; i++) {
		double sum_cosa = 0.0;
		double sum_sina = 0.0;
		for (int j=0; j<num_ali; j++) {
			if (static_cast<int>(ali_params[j*nima*4+i*4+3]) == 0) {
				sum_cosa += cos((args[j*3]+ali_params[j*nima*4+i*4])*M_PI/180.0);
				sum_sina += sin((args[j*3]+ali_params[j*nima*4+i*4])*M_PI/180.0);
				sx[j] = args[j*3+1] + ali_params[j*nima*4+i*4+1]*cosa[j] + ali_params[j*nima*4+i*4+2]*sina[j];
				sy[j] = args[j*3+2] - ali_params[j*nima*4+i*4+1]*sina[j] + ali_params[j*nima*4+i*4+2]*cosa[j];
			} else {
				sum_cosa += cos((-args[j*3]+ali_params[j*nima*4+i*4])*M_PI/180.0);
				sum_sina += sin((-args[j*3]+ali_params[j*nima*4+i*4])*M_PI/180.0);
				sx[j] = -args[j*3+1] + ali_params[j*nima*4+i*4+1]*cosa[j] - ali_params[j*nima*4+i*4+2]*sina[j];
				sy[j] =  args[j*3+2] + ali_params[j*nima*4+i*4+1]*sina[j] + ali_params[j*nima*4+i*4+2]*cosa[j];
			}
		}
		double sqrtP = sqrt(sum_cosa*sum_cosa+sum_sina*sum_sina);
		sqr_pixel_error[i] = d*d/4.0*(1.0-sqrtP/num_ali)+var(sx, num_ali)+var(sy, num_ali);
	}
	
	delete[] args;
	delete[] cosa;
	delete[] sina;
	delete[] sx;
	delete[] sy;
	
	return sqr_pixel_error;
}

void Util::multi_align_error_dfunc(double* x, vector<float> ali_params, int nima, int num_ali, double* g, int d) {

	for (int i=0; i<num_ali*3-3; i++)    g[i] = 0.0;

	double* args = new double[num_ali*3];
	for (int i=0; i<3*num_ali-3; i++)   args[i] = x[i];
	args[3*num_ali-3] = 0.0;
	args[3*num_ali-2] = 0.0;
	args[3*num_ali-1] = 0.0;
	double* cosa = new double[num_ali];
	double* sina = new double[num_ali];
	for (int i=0; i<num_ali; i++) {
		cosa[i] = cos(args[i*3]*M_PI/180.0);
		sina[i] = sin(args[i*3]*M_PI/180.0);
	}
	double* sx = new double[num_ali];
	double* sy = new double[num_ali];
	
	vector<float> sqr_pixel_error(nima);

	for (int i=0; i<nima; i++) {
		double sum_cosa = 0.0;
		double sum_sina = 0.0;
		for (int j=0; j<num_ali; j++) {
			if (static_cast<int>(ali_params[j*nima*4+i*4+3]) == 0) {
				sum_cosa += cos((args[j*3]+ali_params[j*nima*4+i*4])*M_PI/180.0);
				sum_sina += sin((args[j*3]+ali_params[j*nima*4+i*4])*M_PI/180.0);
				sx[j] = args[j*3+1] + ali_params[j*nima*4+i*4+1]*cosa[j] + ali_params[j*nima*4+i*4+2]*sina[j];
				sy[j] = args[j*3+2] - ali_params[j*nima*4+i*4+1]*sina[j] + ali_params[j*nima*4+i*4+2]*cosa[j];
			} else {
				sum_cosa += cos((-args[j*3]+ali_params[j*nima*4+i*4])*M_PI/180.0);
				sum_sina += sin((-args[j*3]+ali_params[j*nima*4+i*4])*M_PI/180.0);
				sx[j] = -args[j*3+1] + ali_params[j*nima*4+i*4+1]*cosa[j] - ali_params[j*nima*4+i*4+2]*sina[j];
				sy[j] =  args[j*3+2] + ali_params[j*nima*4+i*4+1]*sina[j] + ali_params[j*nima*4+i*4+2]*cosa[j];
			}
		}
		double P = sqrt(sum_cosa*sum_cosa+sum_sina*sum_sina);
		sum_cosa /= P;
		sum_sina /= P;
		for (int j=0; j<num_ali-1; j++) {
			double dx = 2.0*(sx[j]-mean(sx, num_ali));
			double dy = 2.0*(sy[j]-mean(sy, num_ali));
			if (static_cast<int>(ali_params[j*nima*4+i*4+3]) == 0) {
				g[j*3] += (d*d/4.0*(sum_cosa*sin((args[j*3]+ali_params[j*nima*4+i*4])*M_PI/180.0) -
					 	    sum_sina*cos((args[j*3]+ali_params[j*nima*4+i*4])*M_PI/180.0)) +
						    dx*(-ali_params[j*nima*4+i*4+1]*sina[j]+ali_params[j*nima*4+i*4+2]*cosa[j])+
						    dy*(-ali_params[j*nima*4+i*4+1]*cosa[j]-ali_params[j*nima*4+i*4+2]*sina[j]))*M_PI/180.0;
				g[j*3+1] += dx;
				g[j*3+2] += dy;
			} else {
				g[j*3] += (d*d/4.0*(-sum_cosa*sin((-args[j*3]+ali_params[j*nima*4+i*4])*M_PI/180.0) +
					 	     sum_sina*cos((-args[j*3]+ali_params[j*nima*4+i*4])*M_PI/180.0)) +
						    dx*( ali_params[j*nima*4+i*4+1]*sina[j]+ali_params[j*nima*4+i*4+2]*cosa[j])+
						    dy*(-ali_params[j*nima*4+i*4+1]*cosa[j]+ali_params[j*nima*4+i*4+2]*sina[j]))*M_PI/180.0;
				g[j*3+1] += -dx;
				g[j*3+2] += dy;
			}
		}
	}
	
	for (int i=0; i<3*num_ali-3; i++)  g[i] /= (num_ali*nima);
	
	delete[] args;
	delete[] cosa;
	delete[] sina;
	delete[] sx;
	delete[] sy;
}

float Util::ccc_images(EMData* image, EMData* refim, EMData* mask, float ang, float sx, float sy) {

	EMData *rot= new EMData();
	float ccc;

	rot = image->rot_scale_trans2D(ang, sx, sy, 1.0);
	ccc = -rot->cmp("sqeuclidean", refim, Dict("mask", mask));
	delete rot;
	return ccc;
}

vector<float> Util::twoD_fine_ali_SD_G(EMData* image, EMData *refim, EMData* mask, Util::KaiserBessel& kb, float ang, float sxs, float sys) {

	double  x[4];
	int n;
	int l = 3;
	int m = 200;
	double e = 1e-9;
	double step = 0.001;
	float (*my_func)(EMData* , EMData* , EMData* , Util::KaiserBessel&, float , float , float) = ccc_images_G;

	x[1] = ang;
	x[2] = sxs;
	x[3] = sys;

	Steepda_G(x, step, e, l, m, &n, my_func, image, refim, mask, kb);   // Call steepest descent optimization subroutine
	//printf("Took %d steps\n", n);

	vector<float> res;
	res.push_back(static_cast<float>(x[1]));
	res.push_back(static_cast<float>(x[2]));
	res.push_back(static_cast<float>(x[3]));
	res.push_back(static_cast<float>(n));
	return res;
}


float Util::ccc_images_G(EMData* image, EMData* refim, EMData* mask, Util::KaiserBessel& kb, float ang, float sx, float sy) {

	EMData *rot= new EMData();
	float ccc;

	rot = image->rot_scale_conv7(static_cast<float>(ang*pi/180.0), sx, sy, kb, 1.0f);
	ccc = -rot->cmp("sqeuclidean", refim, Dict("mask", mask));
	delete rot;
	return ccc;
}

#define img_ptr(i,j,k)  img_ptr[i+(j+(k*ny))*(size_t)nx]
#define img2_ptr(i,j,k) img2_ptr[i+(j+(k*ny))*(size_t)nx]
EMData* Util::move_points(EMData* img, float qprob, int ri, int ro)
{
	ENTERFUNC;
	/* Exception Handle */
	if (!img) {
		throw NullPointerException("NULL input image");
	}

	cout <<"  VERSION  03/19/2015  10:00 AM"<<endl;

	exit(0);

	float dummy;
	dummy = ri;
	cout <<  "   fmod   "<<fmod(qprob,dummy)<<endl;
	cout <<  "   round   "<<round(qprob)<<endl;
	int newx, newy, newz;
	bool  keep_going;
	int nx=img->get_xsize(),ny=img->get_ysize(),nz=img->get_zsize();
	//int size = nx*ny*nz;
	EMData * img2 = new EMData();
	img2->set_size(nx,ny,nz);
	img2->to_zero();
	float *img_ptr  = img->get_data();
	float *img2_ptr = img2->get_data();
	int r2 = ro*ro;
	int r3 = r2*ro;
	int ri2 = ri*ri;
	int ri3 = ri2*ri;

	int n2 = nx/2;

	for (int k=-n2; k<=n2; k++) {		//cout << " k   "<<k <<endl;
		float z2 = static_cast<float>(k*k);
		for (int j=-n2; j<=n2; j++) {
			float y2 = z2 + j*j;
			if(y2 <= r2) {
											//cout << "  j  "<<j <<endl;

				for (int i=-n2; i<=n2; i++) {
					float x2 = y2 + i*i;
					if(x2 <= r3) {
											//cout << "  i  "<<i <<endl;
						int ib = i+n2; int jb = j+n2; int kb = k+n2;
						if(x2 >= ri3) {
							//  this is the outer shell, here points can only vanish
							if(img_ptr(ib,jb,kb) == 1.0f) {
								//cout << "  1  "<<ib <<endl;
								if(Util::get_frand(0.0f, 1.0f) > qprob){
									img2_ptr(ib,jb,kb) = 0.0f;
									keep_going = true;
								//cout << "  try  "<<ib <<endl;
									while(keep_going) {
										newx = Util::get_irand(-ro,ro);
										newy = Util::get_irand(-ro,ro);
										newz = Util::get_irand(-ro,ro);
										if(newx*newx+newy*newy+newz*newz <= r3) {
											newx += n2; newy += n2; newz += n2;
											if( img_ptr(newx,newy,newz) == 0.0f) {
												img2_ptr(newx,newy,newz) = 1.0f;
												keep_going = false; }
										}
									}
								}   else  img2_ptr(ib,jb,kb) = 1.0f;
							}
						}  else  {
							// this is inner shell, the point can only move (or vanish, if all neighbors exist)
							if(img_ptr(ib,jb,kb) == 1.0) {
								if(Util::get_frand(0.0f,1.0f) > qprob) {
									//  find out the number of neighbors
									float  numn = -1.0f;  // we already know the central one is 1
									for (newz = -1; newz <= 1; newz++)
										for (newy = -1; newy <= 1; newy++)
											for (newx = -1; newx <= 1; newx++)
												numn += img_ptr(ib+newx,jb+newy,kb+newz);
									img2_ptr(ib,jb,kb) = 0.0;
									if(numn == 26.0f) {
										//  all neighbors exist, it has to vanish
										keep_going = true;
										while(keep_going) {
											newx = Util::get_irand(-ro,ro);
											newy = Util::get_irand(-ro,ro);
											newz = Util::get_irand(-ro,ro);
											if(newx*newx+newy*newy+newz*newz < r3) {
												newx += n2; newy += n2; newz += n2;
												if( img_ptr(newx,newy,newz) == 0.0f) {
													if(newx*newx+newy*newy+newz*newz < r3) {
														if(newx*newx+newy*newy+newz*newz < r3) {
															newx += n2; newy += n2; newz += n2;
															if( img_ptr(newx,newy,newz) == 0.0f) {
																img2_ptr(newx,newy,newz) = 1.0f;
																keep_going = false; }
														}
													}
												}
											}
										}
									}  else if(numn == 25.0f) {
										// there is only one empty neighbor, move there
										for (newz = -1; newz <= 1; newz++) {
											for (newy = -1; newy <= 1; newy++) {
												for (newx = -1; newx <= 1; newx++) {
													if( newx != 0 && newy != 0 && newz != 0) {
														if(img_ptr(newx+ib,newy+jb,newz+kb) == 0.0f) {
															img2_ptr(newx+ib,newy+jb,newz+kb) = 1.0f;
															}
													}
												}
											}
										}
									}  else {
										//  more than one neighbor is zero, select randomly one and move there
										keep_going = true;
										while(keep_going) {
											newx = Util::get_irand(-1,1);
											newy = Util::get_irand(-1,1);
											newz = Util::get_irand(-1,1);
											if(newx != 0 && newy != 0 && newz != 0)  {
												if(img_ptr(ib+newx,jb+newy,kb+newz) == 0.0f) {
													img2_ptr(ib+newx,jb+newy,kb+newz) = 1.0f;//?????
													keep_going = false;
												}
											}
										}
									}
								}  else img2_ptr(ib,jb,kb) = 1.0f;
							}
						}
					}
				}
			}
		}
	}
	//for (int i=0;i<size;i++) img2_ptr[i] = img_ptr[i] + Util::get_irand(1,7);
	img2->update();

	EXITFUNC;
	return img2;
}
#undef img_ptr
#undef img2_ptr

struct point3d_t
{
	point3d_t( int ix, int iy, int iz ): x(ix), y(iy), z(iz) {}

	int x;
	int y;
	int z;
};


int find_group( int ix, int iy, int iz, int grpid, EMData* mg, EMData* visited )
{
	int offs[][3] = { {-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1} };
	int noff = 6;

	int nx = visited->get_xsize();
	int ny = visited->get_ysize();
	int nz = visited->get_zsize();

	vector< point3d_t > pts;
	pts.push_back( point3d_t(ix, iy, iz) );
	visited->set_value_at( ix, iy, iz, (float)grpid );

	int start = 0;
	int end = pts.size();

	while( end > start ) {
		for(int i=start; i < end; ++i ) {
			int ix = pts[i].x;
			int iy = pts[i].y;
			int iz = pts[i].z;

			for( int j=0; j < noff; ++j ) {
				int jx = ix + offs[j][0];
				int jy = iy + offs[j][1];
				int jz = iz + offs[j][2];

				if( jx < 0 || jx >= nx ) continue;
				if( jy < 0 || jy >= ny ) continue;
				if( jz < 0 || jz >= nz ) continue;


				if( (*mg)(jx, jy, jz)>0 && (*visited)(jx, jy, jz)==0.0 ) {
				    pts.push_back( point3d_t(jx, jy, jz) );
				    visited->set_value_at( jx, jy, jz, (float)grpid );
				}

			}
		}

		start = end;
		end = pts.size();
	}
	return pts.size();
}


EMData* Util::get_biggest_cluster( EMData* mg )
{
	int nx = mg->get_xsize();
	int ny = mg->get_ysize();
	int nz = mg->get_zsize();

	EMData* visited = new EMData();
	visited->set_size( nx, ny, nz );
	visited->to_zero();
	int grpid = 0;
	int maxgrp = 0;
	int maxsize = 0;
	for( int iz=0; iz < nz; ++iz ) {
		for( int iy=0; iy < ny; ++iy ) {
			for( int ix=0; ix < nx; ++ix ) {
				if( (*mg)(ix, iy, iz)==0.0 ) continue;

				if( (*visited)(ix, iy, iz) > 0.0 ) {
					// visited before, must be in other group.
					continue;
				}

				grpid++;
				int grpsize = find_group( ix, iy, iz, grpid, mg, visited );
				if( grpsize > maxsize ) {
					maxsize = grpsize;
					maxgrp = grpid;
				}
			}
		}
	}

	Assert( maxgrp > 0 );

	int npoint = 0;
	EMData* result = new EMData();
	result->set_size( nx, ny, nz );
	result->to_zero();

	for( int iz=0; iz < nz; ++iz ) {
		for( int iy=0; iy < ny; ++iy ) {
			for( int ix=0; ix < nx; ++ix ) {
				if( (*visited)(ix, iy, iz)==maxgrp ) {
					(*result)(ix,iy,iz) = 1.0;
					npoint++;
				}
			}
		}
	}

	Assert( npoint==maxsize );
	delete visited;
	return result;

}

EMData* Util::ctf_img(int nx, int ny, int nz, float dz,float ps,float voltage,float cs, float wgh, float b_factor,float dza, float azz, float sign)
{
	int    ix, iy, iz;
	int    i,  j, k;
	int    nr2, nl2;
	float  ak;
	float  scx, scy, scz;
	int    offset = 2 - nx%2;
	int    lsm = nx + offset;
	EMData* ctf_img1 = new EMData();
	ctf_img1->set_size(lsm, ny, nz);
	float freq = 1.0f/(2.0f*ps);
	scx = 2.0f/float(nx);
	if(ny>=1) scy = 2.0f/float(ny); else scy=0.0f;
	if(nz>=1) scz = 2.0f/float(nz); else scz=0.0f;
	nr2 = ny/2 ;
	nl2 = nz/2 ;
	for ( k=0; k<nz;k++) {
		iz = k;  if(k>nl2) iz=k-nz;
		float oz2 = iz*scz*iz*scz;
		for ( j=0; j<ny;j++) {
			iy = j;  if(j>nr2) iy=j - ny;
			float oy = iy*scy;
			float oy2 = oy*oy;
			for ( i=0; i<lsm/2; i++) {
				ix=i;
				if( dza == 0.0f) {
					ak=pow(ix*ix*scx*scx + oy2 + oz2, 0.5f)*freq;
					(*ctf_img1) (i*2,j,k)   = Util::tf(dz, ak, voltage, cs, wgh, b_factor, sign);
				} else {
					float ox = ix*scx;
					ak=pow(ox*ox + oy2 + oz2, 0.5f)*freq;
					float dzz = dz - dza/2.0f*sin(2*(atan2(oy, ox)+azz*M_PI/180.0f));
					(*ctf_img1) (i*2,j,k)   = Util::tf(dzz, ak, voltage, cs, wgh, b_factor, sign);
				}
				//(*ctf_img1) (i*2+1,j,k) = 0.0f;  PAP  I assumed new EMData sets to zero
			}
		}
	}
	ctf_img1->update();
	ctf_img1->set_complex(true);
	ctf_img1->set_ri(true);
	//ctf_img1->attr_dict["is_complex"] = 1;
	//ctf_img1->attr_dict["is_ri"] = 1;
	if(nx%2==0) ctf_img1->set_fftodd(false); else ctf_img1->set_fftodd(true);
	return ctf_img1;
}



EMData* Util::ctf_rimg(int nx, int ny, int nz, float dz, float ps, float voltage, float cs, float wgh, float b_factor, float dza, float azz, float sign)
{
	int    ix, iy, iz;
	int    i,  j, k;
	float  ak;
	float  scx, scy, scz;
	EMData* ctf_img1 = new EMData();
	ctf_img1->set_size(nx, ny, nz);
	float freq = 1.0f/(2.0f*ps);
	scx = 2.0f/float(nx);
	if(ny>=1) scy = 2.0f/float(ny); else scy=0.0f;
	if(nz>=1) scz = 2.0f/float(nz); else scz=0.0f;
	int ns2 = nx/2 ;
	int nr2 = ny/2 ;
	int nl2 = nz/2 ;
	int nod = nx%2 ;
	int nok = ny%2 ;
	int noz = nz%2 ;
	for ( k=0; k<nz;k++) {
		iz = k - nl2;
		int kz = (nz - k - noz)%nz;
		float oz2 = iz*scz*iz*scz;
		for ( j=0; j<ny;j++) {
			iy = j - nr2;
			int jy = (ny - j - nok)%ny;
			float oy = iy*scy;
			float oy2 = oy*oy;
			for ( i=0; i<=ns2; i++) {
				ix = i - ns2;
				if( dza == 0.0f) {
					ak=pow(ix*ix*scx*scx + oy2 + oz2, 0.5f)*freq;
					(*ctf_img1) (i,j,k)   = Util::tf(dz, ak, voltage, cs, wgh, b_factor, sign);
				} else {
					float ox = ix*scx;
					ak=pow(ox*ox + oy2 + oz2, 0.5f)*freq;
					float dzz = dz - dza/2.0f*sin(2*(atan2(oy, ox)+azz*M_PI/180.0f));
					(*ctf_img1) (i,j,k)   = Util::tf(dzz, ak, voltage, cs, wgh, b_factor, sign);
				}
				ix = nx - i - nod;
				if(ix<nx)  (*ctf_img1) (ix,jy,kz) = (*ctf_img1) (i,j,k);
			}
		}
	}
	ctf_img1->update();
	return ctf_img1;
}




EMData* Util::ctf2_rimg(int nx, int ny, int nz, float dz, float ps, float voltage, float cs, float wgh, float b_factor, float dza, float azz, float sign)
{
	int    ix, iy, iz;
	int    i,  j, k;
	float  ak;
	float  scx, scy, scz;
	EMData* ctf_img1 = new EMData();
	ctf_img1->set_size(nx, ny, nz);
	float freq = 1.0f/(2.0f*ps);
	scx = 2.0f/float(nx);
	if(ny>=1) scy = 2.0f/float(ny); else scy=0.0f;
	if(nz>=1) scz = 2.0f/float(nz); else scz=0.0f;
	int ns2 = nx/2 ;
	int nr2 = ny/2 ;
	int nl2 = nz/2 ;
	int nod = nx%2 ;
	int nok = ny%2 ;
	int noz = nz%2 ;
	for ( k=0; k<nz;k++) {
		iz = k - nl2;
		int kz = (nz - k - noz)%nz;
		float oz2 = iz*scz*iz*scz;
		for ( j=0; j<ny;j++) {
			iy = j - nr2;
			int jy = (ny - j - nok)%ny;
			float oy = iy*scy;
			float oy2 = oy*oy;
			for ( i=0; i<=ns2; i++) {
				ix = i - ns2;
				if( dza == 0.0f) {
					ak=pow(ix*ix*scx*scx + oy2 + oz2, 0.5f)*freq;
					(*ctf_img1) (i,j,k)   = pow(Util::tf(dz, ak, voltage, cs, wgh, b_factor, sign),2);
				} else {
					float ox = ix*scx;
					ak=pow(ox*ox + oy2 + oz2, 0.5f)*freq;
					//az = atan2(oy, ox);
					//float dzz = dz + dza/2.0f*sin(2*(az-azz*M_PI/180.0f-pihalf));
					float dzz = dz - dza/2.0f*sin(2*(atan2(oy, ox)+azz*M_PI/180.0f));
					(*ctf_img1) (i,j,k)   = pow(Util::tf(dzz, ak, voltage, cs, wgh, b_factor, sign),2);
				}
				ix = nx - i - nod;
				if(ix<nx)  (*ctf_img1) (ix,jy,kz) = (*ctf_img1) (i,j,k);
			}
		}
	}
	ctf_img1->update();
	return ctf_img1;
}



/*
#define  cent(i)     out[i+N]
#define  assign(i)   out[i]
vector<float> Util::cluster_pairwise(EMData* d, int K) {

	int nx = d->get_xsize();
	int N = 1 + int((sqrt(1.0 + 8.0*nx)-1.0)/2.0);
	vector<float> out(N+K+2);
	if(N*(N-1)/2 != nx) {
		//print  "  incorrect dimension"
		return out;}
	//  assign random objects as centers
	for(int i=0; i<N; i++) assign(i) = float(i);
	// shuffle
	for(int i=0; i<N; i++) {
		int j = Util::get_irand(0,N-1);
		float temp = assign(i);
		assign(i) = assign(j);
		assign(j) = temp;
	}
	for(int k=0; k<K; k++) cent(k) = float(assign(k));
	//for(int k=0; k<K; k++) cout<<cent(k)<<"    ";cout<<endl;
	//
	for(int i=0; i<N; i++) assign(i) = 0.0f;
	float qm, dispold = 1.1e22f, disp = 1.0e22f, na=0.0f;
	bool change = true;
	int it = -1;
	while(change && disp < dispold) {
		change = false;
		dispold = disp;
		it++;
		//cout<<"Iteration:  "<<it<<endl;
		// dispersion is a sum of distance from objects to object center
		disp = 0.0f;
		for(int i=0; i<N; i++) {
			qm = 1.0e23f;
			for(int k=0; k<K; k++) {
				if(float(i) == cent(k)) {
					qm = 0.0f;
					na = (float)k;
				} else {
					float dt = (*d)(mono(i,int(cent(k))));
					if(dt < qm) {
						qm = dt;
						na = (float)k;
					}
				}
			}
			disp += qm;
			if(na != assign(i)) {
				assign(i) = na;
				change = true;
			}
		}
	//for(int k=0; k<N; k++) cout<<assign(k)<<"    ";cout<<endl;
		//print disp
		//print  assign
		// find centers
		for(int k=0; k<K; k++) {
			qm = 1.0e23f;
			for(int i=0; i<N; i++) {
				if(assign(i) == float(k)) {
					float q = 0.0;
					for(int j=0; j<N; j++) {
						if(assign(j) == float(k)) {
								//it cannot be the same object
							if(i != j)  q += (*d)(mono(i,j));
							//cout<<q<<"   "<<i<<"   "<<j<<"   "<<k<<endl;}
						}
					}
					if(q < qm) {
						//cout<<qm<<"   "<<q<<"   "<<i<<"   "<<k<<endl;
						qm = q;
						cent(k) = float(i);
					}
				}
			}
		}
	//for(int k=0; k<K; k++) cout<<cent(k)<<"    ";cout<<endl;cout<<disp<<endl;
	}
	out[N+K] = disp;
	out[N+K+1] = float(it);
	return  out;
}
#undef  cent
#undef  assign
*/
#define  cent(i)     out[i+N]
#define  assign(i)   out[i]
vector<float> Util::cluster_pairwise(EMData* d, int K, float T, float F) {
	int nx = d->get_xsize();
	int N = 1 + int((sqrt(1.0 + 8.0*nx)-1.0)/2.0);
	vector<float> out(N+K+2);
	if(N*(N-1)/2 != nx) {
		//print  "  incorrect dimension"
		return out;}
	//  assign random objects as centers
	for(int i=0; i<N; i++) assign(i) = float(i);
	// shuffle
	for(int i=0; i<N; i++) {
		int j = Util::get_irand(0,N-1);
		float temp = assign(i);
		assign(i) = assign(j);
		assign(j) = temp;
	}
	for(int k=0; k<K; k++) cent(k) = float(assign(k));
	//for(int k=0; k<K; k++) cout<<cent(k)<<"    ";cout<<endl;
	//
	for(int i=0; i<N; i++) assign(i) = 0.0f;
	float qm, dispold = 1.1e22f, disp = 1.0e22f, na=0.0f;
	bool change = true;
	int it = -1;
	int ct = -1;
	while ((change && disp < dispold) || ct > 0) {

		change = false;
		dispold = disp;
		it++;

		// dispersion is a sum of distance from objects to object center
		disp = 0.0f;
		ct = 0;
		for(int i=0; i<N; i++) {
			qm = 1.0e23f;
			for(int k=0; k<K; k++) {
				if(float(i) == cent(k)) {
					qm = 0.0f;
					na = (float)k;
				} else {
					float dt = (*d)(mono(i,int(cent(k))));
					if(dt < qm) {
						qm = dt;
						na = (float)k;
					}
				}
			}


			// Simulated annealing
			if(exp(-1.0/float(T)) > Util::get_irand(1,1000)/1000.0) {
			    na = (float)(Util::get_irand(0, K));
			    qm = (*d)(mono(i,int(na)));
			    ct++;
			}

			disp += qm;

			if(na != assign(i)) {
				assign(i) = na;
				change = true;
			}
		}

		//cout<<"Iteration:  "<<it<< " disp " <<disp << " T " << T << " disturb " << ct << endl;
		T = T*F;

	//for(int k=0; k<N; k++) cout<<assign(k)<<"    ";cout<<endl;
		//print disp
		//print  assign
		// find centers
		for(int k=0; k<K; k++) {
			qm = 1.0e23f;
			for(int i=0; i<N; i++) {
				if(assign(i) == float(k)) {
					float q = 0.0;
					for(int j=0; j<N; j++) {
						if(assign(j) == float(k)) {
								//it cannot be the same object
							if(i != j)  q += (*d)(mono(i,j));
							//cout<<q<<"   "<<i<<"   "<<j<<"   "<<k<<endl;}
						}
					}
					if(q < qm) {
						//cout<<qm<<"   "<<q<<"   "<<i<<"   "<<k<<endl;
						qm = q;
						cent(k) = float(i);
					}
				}
			}
		}
	//for(int k=0; k<K; k++) cout<<cent(k)<<"    ";cout<<endl;cout<<disp<<endl;
	}
	out[N+K] = disp;
	out[N+K+1] = float(it);
	return  out;
}
#undef  cent
#undef  assign
/*
#define  groupping(i,k)   group[i + k*m]
vector<float> Util::cluster_equalsize(EMData* d, int m) {
	int nx = d->get_xsize();
	int N = 1 + int((sqrt(1.0 + 8.0*nx)-1.0)/2.0);
	int K = N/m;
	//cout<<"  K  "<<K<<endl;
	vector<float> group(N+1);
	if(N*(N-1)/2 != nx) {
		//print  "  incorrect dimension"
		return group;}
	bool active[N];
	for(int i=0; i<N; i++) active[i] = true;

	float dm, qd;
	int   ppi, ppj;
	for(int k=0; k<K; k++) {
		// find two most similiar objects among active
		cout<<"  k  "<<k<<endl;
		dm = 1.0e23;
		for(int i=1; i<N; i++) {
			if(active[i]) {
				for(int j=0; j<i; j++) {
					if(active[j]) {
						qd = (*d)(mono(i,j));
						if(qd < dm) {
							dm = qd;
							ppi = i;
							ppj = j;
						}
					}
				}
			}
		}
		groupping(0,k) = float(ppi);
		groupping(1,k) = float(ppj);
		active[ppi] = false;
		active[ppj] = false;

		// find progressively objects most similar to those in the current list
		for(int l=2; l<m; l++) {
			//cout<<"  l  "<<l<<endl;
			dm = 1.0e23;
			for(int i=0; i<N; i++) {
				if(active[i]) {
					qd = 0.0;
					for(int j=0; j<l; j++) { //j in groupping[k]:
			//cout<<"  groupping(j,k)  "<<groupping(j,k)<<"   "<<j<<endl;
						int jj = int(groupping(j,k));
			//cout<<"   "<<jj<<endl;
						qd += (*d)(mono(i,jj));
					}
					if(qd < dm) {
						dm = qd;
						ppi = i;
					}
				}
			}
			groupping(l,k) = float(ppi);
			active[ppi] = false;
			//cout<<"  k  "<<k<<"  l  "<<l<<"  ppi  "<<ppi<<"  dm  "<<dm<<endl;
		}
		//for(int j=k*m; j<(k+1)*m; j++) cout<<group[j]<<endl;
	}
	// there might be remaining objects when N is not divisible by m, simply put them in one group
	if(N%m != 0) {
		int j = K*m;
		K++;
		//cout<<"  N%m  "<<N%m<<"   "<<K<<endl;
		for(int i=0; i<N; i++) {
			if(active[i]) {
				group[j] = float(i);
		//cout<<"  j,i  "<<j<<"   "<<i<<"   "<<K-1 <<"   "<< group[j]<<endl;
				j++;
			}
		}
	}

	int  cent[K];
	 // find centers
	for(int k=0; k<K; k++) {
		float qm = 1.0e23f;
		for(int i=0; i<N; i++) {
			if(group[i] == float(k)) {
				qd = 0.0;
				for(int j=0; j<N; j++) {
					if(group[j] == float(k)) {
						//it cannot be the same object
						if(i != j)  qd += (*d)(mono(i,j));
					}
				}
				if(qd < qm) {
					qm = qd;
					cent[k] = i;
				}
			}
		}
	}
	// dispersion is a sum of distances from objects to object center
	float disp = 0.0f;
	for(int i=0; i<N; i++) {
		for(int k=0; k<K; k++) if(i == cent[k]) disp += (*d)(mono(i,cent[k]));
	}
	group[N] = disp;
	return  group;
}
#undef  groupping
*/

vector<float> Util::cluster_equalsize(EMData* d) {
	//  WORKS ONLY FOR NUMBER OF OBJECTS N=l^2   !!
	int nx = d->get_xsize();
	int N = 1 + int((sqrt(1.0 + 8.0*nx)-1.0)/2.0);
	int K = N/2;
	vector<float> group(N);
	if(N*(N-1)/2 != nx) {
		//print  "  incorrect dimension"
		return group;}
	//bool active[N];	//this does not compile in VS2005. --Grant Tang
	bool * active = new bool[N];
	for(int i=0; i<N; i++) active[i] = true;

	float dm, qd;
	int   ppi = 0, ppj = 0;
	for(int k=0; k<K; k++) {
		// find pairs of most similiar objects among active
		//cout<<"  k  "<<k<<endl;
		dm = 1.0e23f;
		for(int i=1; i<N; i++) {
			if(active[i]) {
				for(int j=0; j<i; j++) {
					if(active[j]) {
						qd = (*d)(i*(i - 1)/2 + j);
						if(qd < dm) {
							dm = qd;
							ppi = i;
							ppj = j;
						}
					}
				}
			}
		}
		group[2*k] = float(ppi);
		group[1+2*k] = float(ppj);
		active[ppi] = false;
		active[ppj] = false;
	}

	delete [] active;
	active = NULL;
	return  group;
}
/*
#define son(i,j)=i*(i-1)/2+j
vector<float> Util::cluster_equalsize(EMData* d) {
	//  WORKS ONLY FOR NUMBER OF OBJECTS N=l^2   !!
	int nx = d->get_xsize();
	int N = 1 + int((sqrt(1.0 + 8.0*nx)-1.0)/2.0);
	int K = N/2;
	vector<float> group(N);
	if(N*(N-1)/2 != nx) {
		//print  "  incorrect dimension"
		return group;}
	//bool active[N];
	int  active[N];
	for(int i=0; i<N; i++) active[i] = i;

	float dm, qd;
	int   ppi = 0, ppj = 0, ln = N;
	for(int k=0; k<K; k++) {
		// find pairs of most similiar objects among active
		//cout<<"  k:  "<<k<<endl;
		dm = 1.0e23;
		for(int i=1; i<ln; i++) {
			for(int j=0; j<i; j++) {
				//qd = (*d)(mono(active[i],active[j]));
				qd = (*d)(active[i]*(active[i] - 1)/2 + active[j]);
				if(qd < dm) {
					dm = qd;
					ppi = i;
					ppj = j;
				}
			}
		}
		group[2*k]   = float(active[ppi]);
		group[1+2*k] = float(active[ppj]);
		//  Shorten the list
		if(ppi > ln-3 || ppj > ln - 3) {
			if(ppi > ln-3 && ppj > ln - 3) {
			} else if(ppi > ln-3) {
				if(ppi == ln -1) active[ppj] = active[ln-2];
				else             active[ppj] = active[ln-1];
			} else { // ppj>ln-3
				if(ppj == ln -1) active[ppi] = active[ln-2];
				else             active[ppi] = active[ln-1];
			}
		} else {
			active[ppi] = active[ln-1];
			active[ppj] = active[ln-2];
		}
		ln = ln - 2;
	}
	return  group;
}

*/
#define data(i,j) group[i*ny+j]
vector<float> Util::vareas(EMData* d) {
	const float step=0.001f;
	int ny = d->get_ysize();
	//  input emdata should have size 2xN, where N is number of points
	//  output vector should be 2xN, first element is the number of elements
	//  associated with this point, second is 0 is the element is touching the border, 1 if it is interior
	vector<float> group(2*ny);
	for(int i=0; i<2*ny; i++) group[i] = 0.0f;
	int K = int(1.0f/step) +1;
	int hit = 0;
	for(int kx=0; kx<=K; kx++) {
		float tx = kx*step;
		for(int ky=0; ky<=K; ky++) {
			float ty = ky*step;
			float dm = 1.0e23f;
			for(int i=0; i<ny; i++) {
				float qd = pow(tx-(*d)(0,i),2) + pow(ty-(*d)(1,i),2);
				if( qd < dm) {
					dm = qd;
					hit = i;
				}
			}
			data(0,hit) += 1.0f;
			if(kx == 0 || ky == 0 || kx == K || ky == K)  data(1,hit) = 1.0f;
		}
	}
	return  group;
}
#undef data

EMData* Util::get_slice(EMData *vol, int dim, int index) {

	int nx = vol->get_xsize();
	int ny = vol->get_ysize();
	int nz = vol->get_zsize();
	float *vol_data = vol->get_data();
	int new_nx, new_ny;

	if (nz == 1)
		throw ImageDimensionException("Error: Input must be a 3-D object");
	if ((dim < 1) || (dim > 3))
		throw ImageDimensionException("Error: dim must be 1 (x-dimension), 2 (y-dimension) or 3 (z-dimension)");
	if (((dim == 1) && (index < 0 || index > nx-1)) ||
	  ((dim == 1) && (index < 0 || index > nx-1)) ||
	  ((dim == 1) && (index < 0 || index > nx-1)))
		throw ImageDimensionException("Error: index exceeds the size of the 3-D object");

	if (dim == 1) {
		new_nx = ny;
		new_ny = nz;
	} else if (dim == 2) {
		new_nx = nx;
		new_ny = nz;
	} else {
		new_nx = nx;
		new_ny = ny;
	}

	EMData *slice = new EMData();
	slice->set_size(new_nx, new_ny, 1);
	float *slice_data = slice->get_data();

	if (dim == 1) {
		for (int x=0; x<new_nx; x++)
			for (int y=0; y<new_ny; y++)
				slice_data[y*new_nx+x] = vol_data[(y*ny+x)*nx+index];
	} else if (dim == 2) {
		for (int x=0; x<new_nx; x++)
			for (int y=0; y<new_ny; y++)
				slice_data[y*new_nx+x] = vol_data[(y*ny+index)*nx+x];
	} else {
		for (int x=0; x<new_nx; x++)
			for (int y=0; y<new_ny; y++)
				slice_data[y*new_nx+x] = vol_data[((size_t)index*ny+y)*nx+x];
	}

	return slice;
}

void Util::image_mutation(EMData *img, float mutation_rate) {
	int nx = img->get_xsize();
	float min = img->get_attr("minimum");
	float max = img->get_attr("maximum");
	float* img_data = img->get_data();
	array_mutation(img_data, nx*nx, mutation_rate, min, max, 8, 0);
	return;
}


void Util::array_mutation(float *list, int len_list, float mutation_rate, float min_val, float max_val, int L, int is_mirror) {

	if (is_mirror != 0) {
		for (int i=0; i<len_list; i++) {
			int r = rand()%10000;
			float f = r/10000.0f;
			if (f < mutation_rate) list[i] = 1-list[i];
		}
	} else {
		map<int, vector<int> >  graycode;
		map<vector<int>, int> rev_graycode;
		vector <int> gray;

		int K=1;
		for (int i=0; i<L; i++) K*=2;

		for (int k=0; k<K; k++) {
			int shift = 0;
			vector <int> gray;
			for (int i=L-1; i>-1; i--) {
				int t = ((k>>i)%2-shift)%2;
				gray.push_back(t);
				shift += t-2;
			}
			graycode[k] = gray;
			rev_graycode[gray] = k;
		}

		float gap = (K-1)/(max_val-min_val);
		for (int i=0; i<len_list; i++) {
			float val = list[i];
			if (val < min_val) { val = min_val; }
			else if  (val > max_val) { val = max_val; }
			int k = int((val-min_val)*gap+0.5);
			vector<int> gray = graycode[k];
			bool changed = false;
			for (vector<int>::iterator p=gray.begin(); p!=gray.end(); p++) {
				int r = rand()%10000;
				float f = r/10000.0f;
				if (f < mutation_rate) {
					*p = 1-*p;
					changed = true;
				}
			}
			if (changed) {
				k = rev_graycode[gray];
				list[i] = k/gap+min_val;
			}
		}
	}

}

vector<float> Util::list_mutation(vector<float> list, float mutation_rate, float min_val, float max_val, int L, int is_mirror) {

	if (is_mirror != 0) {
		for (vector<float>::iterator q=list.begin(); q!=list.end(); q++) {
			int r = rand()%10000;
			float f = r/10000.0f;
			if (f < mutation_rate) *q = 1-*q;
		}
	} else {
		map<int, vector<int> >  graycode;
		map<vector<int>, int> rev_graycode;
		vector <int> gray;

		int K=1;
		for (int i=0; i<L; i++) K*=2;

		for (int k=0; k<K; k++) {
			int shift = 0;
			vector <int> gray;
			for (int i=L-1; i>-1; i--) {
				int t = ((k>>i)%2-shift)%2;
				gray.push_back(t);
				shift += t-2;
			}
			graycode[k] = gray;
			rev_graycode[gray] = k;
		}

		float gap = (K-1)/(max_val-min_val);
		for (vector<float>::iterator q=list.begin(); q!=list.end(); q++) {
			float val = *q;
			if (val < min_val) { val = min_val; }
			else if  (val > max_val) { val = max_val; }
			int k = int((val-min_val)*gap+0.5);
			vector<int> gray = graycode[k];
			bool changed = false;
			for (vector<int>::iterator p=gray.begin(); p!=gray.end(); p++) {
				int r = rand()%10000;
				float f = r/10000.0f;
				if (f < mutation_rate) {
					*p = 1-*p;
					changed = true;
				}
			}
			if (changed) {
				k = rev_graycode[gray];
				*q = k/gap+min_val;
			}
		}
	}
	return list;
}


bool Util::sanitycheck(int* argParts, int* Indices, int* dimClasses, int nParts, int K, int T, int* output){
	//cout<<"sanitycheck called\n";
	int total_cost = *output;
	int num_matches = *(output+1);

	int cost=0;
	int* intx;
	int intx_size;
	int* intx_next(0);
	int intx_next_size = 0;
	int curclass;
	int curclass_size;
	//cout<<"cost by match: [";
	for(int i = 0; i < num_matches; i++){
		curclass = *(output+2+ i*nParts);
		// check feasibility
		if (*(argParts + Indices[curclass]+1) == -5) {cout<<"infeasible match!\n"; return 0;}
		*(argParts + Indices[curclass]+1) = -5;
		// initialize intx to the number of elements in class curclass of partition 0 of the i-th match
		curclass_size = *(dimClasses+curclass)-2;
		intx = new int[curclass_size];
		for (int ic = 0; ic < curclass_size; ic++) *(intx+ic) = *(argParts + Indices[curclass]+2+ic);
		intx_size = curclass_size;

		for (int j=1; j < nParts; j++){
		      curclass = *(output+2+ i*nParts+j);
		      if (*(argParts + Indices[j*K+curclass]+1)==-5){cout<<"infeasible match!\n"; return 0;}
		      *(argParts + Indices[j*K+curclass]+1)=-5;
		      // compute the intersection of intx and class curclass of partition j of the i-th match
		      intx_next_size = Util::k_means_cont_table_(intx,argParts + Indices[j*K+curclass]+2, intx_next, intx_size, *(dimClasses + j*K+curclass)-2,0);
		      intx_next = new int[intx_next_size];
		      Util::k_means_cont_table_(intx,argParts + Indices[j*K+curclass]+2, intx_next, intx_size, *(dimClasses + j*K+curclass)-2,1);
		      delete[] intx;
		      intx=intx_next;
		      intx_size= intx_next_size;
		}
		delete[] intx_next;

		if (intx_next_size <= T) {cout << "something wrong with solution!\n"; return 0;}
		//cout <<intx_next_size<<",";
		cost = cost + intx_next_size;
	}
	//cout<<"]\n";
	if (cost != total_cost) {cout << "something wrong with solution!\n"; return 0;}

	return 1;

}


// Given J, returns the J matches with the largest weight
// matchlist has room for J matches
// costlist has J elements to record cost of the J largest matches

void Util::search2(int* argParts, int* Indices, int* dimClasses, int nParts, int K, int T, int* matchlist, int* costlist, int J){
	
	// some temp variables
	bool flag = 0;
	int nintx;
	int* dummy(0);
	//int* ret;
	int* curbranch = new int[nParts];
	
	//initialize costlist to all 0
	for(int jit= 0; jit< J; jit++) *(costlist+jit) = 0;
	
	
	for(int a=0; a<K; a++)
	{
	
		// check that class a of partition 0 is active and has greater than T elements. If not the case, then skip to the next class
		if (*(argParts + Indices[a] + 1) < 1) continue;
		if (*(dimClasses + a)-2 <= T) continue;

		// initial pruning: for each partition j>0, set the partition to inactive if its intersection with class a of partition 0 is less than new T

		for( int i=1; i < nParts; i++){
			flag = 0; // if flag stays 0 then no class in this partition has more than T objects in common with a, which implies no feasible match (> T) with class a of part 0 is possible.
			for(int j=0; j < K; j++){
				if (*(argParts + Indices[i*K+j] + 1) < 1) continue;
				nintx = Util::k_means_cont_table_(argParts + Indices[a]+2,argParts + Indices[i*K+j]+2, dummy, *(dimClasses + a)-2, *(dimClasses + i*K+j)-2,0);
				if (nintx > T) flag=1;
				else *(argParts + Indices[i*K+j] + 1) =-4;
			}
			if (flag==0) {break;}
		}

		// explore determines J matchs with the largest weight greater than T where class in partition 0 is class a
		*curbranch = a;

		if (flag > 0) // Each partition has one or more active class
			Util::explore2(argParts, Indices, dimClasses, nParts, K, T, argParts+Indices[a]+2, *(dimClasses+a)-2, argParts+Indices[a]+2,
			*(dimClasses+a)-2,0, J, matchlist, costlist, curbranch);
			
		// take all the classes marked as -4 and remark it as 1 in preparation for next round
		for( int i=1; i < nParts; i++){
			for(int j=0; j < K; j++){
				if (*(argParts + Indices[i*K+j] + 1) == -4) *(argParts + Indices[i*K+j] + 1) =1;

			}
		}
	}
	
	delete[] curbranch;
}

// returns J largest matches
void Util::explore2(int* argParts, int* Indices, int* dimClasses, int nParts, int K, int T, int* curintx, int size_curintx, int* next, int size_next, int depth, int J, int* matchlist, int*costlist, int* curbranch){

// depth is the level which is going to be explored in the current iteration
	int* curintx2(0);
	int nintx = size_curintx;
	
	
	// 2. take the intx of next and cur. Prune if <= T
	if (depth >0){
		nintx = Util::k_means_cont_table_(curintx,next, curintx2, size_curintx, size_next,0);
		if (nintx <= T) return; //prune!
	}

	// 1. we're at a leaf with weight > T, so determine if there is any empty space. If so, put it in. If not, determine if current cost is larger than any of the cost in matchlist, if so, replace the  smallest one in matchlist
	if (depth == (nParts-1)) {
		
		int replace = 0;
		int ind_smallest = -1;
		int smallest_cost = -1;
		
		for (int jit = 0; jit < J; jit++){
			if (*(costlist+jit) < nintx){
				replace = 1;
				if (ind_smallest == -1) {ind_smallest = jit; smallest_cost = *(costlist+jit);}
				if (*(costlist+jit) < smallest_cost) {ind_smallest = jit; smallest_cost = *(costlist+jit);}
			}	
		}
		
		if (replace > 0){
			// replace the smallest cost in matchlist with the current stuff
			*(costlist + ind_smallest) = nintx;
			for (int xit = 0; xit < nParts; xit++)
				*(matchlist + ind_smallest*nParts + xit) = *(curbranch+xit);
				
		}
		
		return;	
	}
	

	// 3. have not yet reached a leaf, and current weight is still greather than T, so keep on going.

	if (depth > 0){
		curintx2 = new int[nintx]; // put the intersection set in here
		Util::k_means_cont_table_(curintx,next,curintx2, size_curintx, size_next,1);
	}

	if (depth == 0){
		// set curintx2 to curintx
		curintx2 = new int[size_curintx];
		for (int cp = 0; cp < size_curintx; cp++) *(curintx2+cp) = *(curintx+cp);
	}


	// recursion (non-leaf case)
	depth=depth+1;
	// we now consider each of the classes in partition depth and recurse upon each of them
	for (int i=0; i < K; i++){

		if (*(argParts + Indices[depth*K+i] + 1) < 1) continue; // class is not active so move on
		size_next = (*(dimClasses + depth*K+i ))-2;
		if (size_next <= T) continue;
		*(curbranch+depth) = i;
		Util::explore2(argParts,Indices, dimClasses, nParts, K, T, curintx2, nintx, argParts + Indices[depth*K+i] + 2, size_next, depth,J, matchlist,
			costlist, curbranch);
		
	}

	delete[] curintx2;
}

void Util::initial_prune(vector <vector <int*> > & Parts, int* dimClasses, int nParts, int K, int T) {
	//cout<<"initial_prune\n";
	// simple initial pruning. For class indClass of partition indPart:
	// For each class of partition which is not indPart, see if there is a class in the partition with which indClass has intersection greater than T
	// If for some partition (not equal to indPart) for which there is no such class, then return 0 immediately.

	// 1. For each class of partition which is not indPart, remove the class from Parts if its intx with indClass of indPart is not gt T

	// remember when calling k_means_cont_table, the first element of each class is an index, and actual number of elements in the class
	// (as stored in dimClasses) should be decremented by 1 accordingly, and second is dummy variable

	int* dummy(0);
	int* cref;
	int cref_size;
	int* ccomp;
	int ccomp_size;
	int nintx;
	for (int i=0; i < nParts; i++){
		for (int j =0; j < K; j++){

			// consider class Parts[i][j]
			cref = Parts[i][j];//incr by 1 since first element is index and second is dummy
			cref_size = dimClasses[i*K+cref[0]]-2;


			if (cref_size <= T){
				cref[0] = -1;
				continue;
			}
			bool done = 0;
			for (int a = 0; a < nParts; a++){
				if (a == i) continue; //consider all classes not in partition i and set to inactive all those classes whose intx with cref is not gt T
				bool hasActive=0;
				for (unsigned int b=0; b < Parts[a].size(); b++){
					// get the card of the intx between Parts[i][j] and Parts[a][b] using k_means_cont_table
					// remember first element of each class is the index of the class
					ccomp = Parts[a][b];
					ccomp_size= dimClasses[a*K+ccomp[0]]-2;
					nintx = Util::k_means_cont_table_(cref+2,ccomp+2, dummy, cref_size, ccomp_size,0);


					if (nintx <= T)
						ccomp[1] = 0; // class Parts[a][b] is 'inactive' for cref
					else{
						ccomp[1] = 1; // class Parts[a][b] is 'active' for cref
						hasActive=1;
					}
				}
				// see if partition a has at least one active class.if not then we're done with cref
				if (hasActive < 1){
				   done=1;
				   break;
				}

			}

			if (done > 0){
				// remove class j from partition i

				cref[0] = -1; // mark for deletion later
				continue; // move on to class Parts[i][j+1]
			}

			// now we adopt more aggressive measures. we find the match with the largest weight which has class Parts[i][j] for partition i.
			// We get rid of Parts[i][j] if this weight is not gt T as no other feasible match containing class Parts[i][j] can be gt T.

			// (To implement later:) To reduce complexity, determine the order the partitions are to be explored based on the cardinality of the active classes of each partition.
			// Note that if there are ALOT of partitions, then sorting the partititions may actually take longer than the time saved by ordering hte
			// partitions. However, we do not have a lot of partitions, (usually around 10), so sorting time will be trivial compared to exploring time.

			// explore returns one if there is a feasible matching containing class cref which has weight gt T, and otherwise returns 0
			//bool found = 1;
			bool found = explore(Parts, dimClasses, nParts, K, T, i, cref+2, cref_size, cref, cref_size,0);

			if (found<1){ // There is NO feasible matching with class j (cref)  with weight greater than T, so delete this class from Parts
				// Parts[i].erase(Parts[i].begin()+j);
				cref[0] = -1;
			}
		}

		// Erase from Parts[i] all the classes that's being designated for erasure

		for (int d = K-1; d > -1; d--){
			if (Parts[i][d][0] < 0) Parts[i].erase(Parts[i].begin()+d);
		}

	}
	//cout <<"number of classes left in each partition after initial prune\n";	
	// Print out how many classes are left in each partition
	//for (int i =0; i < nParts; i++)
	//	cout << Parts[i].size()<<", ";
	//cout << "\n";
}


bool Util::explore(vector <vector <int*> > & Parts, int* dimClasses, int nParts, int K, int T, int partref, int* curintx, int size_curintx, int* next,  int size_next, int depth) {


	if (size_next <= T) return 0;

	// take the intx of next and cur
	int* curintx2(0);
	int nintx = Util::k_means_cont_table_(curintx, next+2, curintx2, size_curintx, size_next,0);
	if (nintx <= T) return 0;

	int old_depth=depth;
	if (depth == partref) depth = depth + 1; // we skip classes in partref
	if (depth == nParts &&  old_depth>0) return 1;

	// have not yet reached a leaf, and current weight is still greather than T, so keep on going.

	curintx2 = new int[nintx]; // put the intersection set in here
	Util::k_means_cont_table_(curintx,next+2,curintx2, size_curintx, size_next,1);

	// if (old_depth > 0) {delete[] curintx;} // don't delete when depth == 0 because curintx in that case is a pointer into Parts

	// we now consider each of the classes in partition (depth+1) in turn
	bool gt_thresh;
	int num_classes = Parts[depth].size(); // (TO DO) have to figure out how many classes partition (depth) has since some may have being removed from before iterations

	for (int i=0; i < num_classes; i++){
		if (Parts[depth][i][1] < 1) continue; // class is not active so move on
		size_next = dimClasses[depth*K + Parts[depth][i][0] ]-2;
		gt_thresh = explore(Parts,dimClasses, nParts, K, T, partref, curintx2,nintx, Parts[depth][i], size_next, depth+1);
		if (gt_thresh) { delete[] curintx2; return 1; }
	}
	delete[] curintx2;
	return 0;
}




vector<int> Util::bb_enumerateMPI_(int* argParts, int* dimClasses, int nParts, int K, int T, int n_guesses, int LARGEST_CLASS, int J,
int max_branching, float stmult, int branchfunc, int LIM) {
	
	
	// Indices is an nParts*K int array storing the index (into argparts) of the first element of the i-th class of the j-th partition
	// So Indices[j*K + i] is the offset from argparts of the first element of the first element  of the i-th class of the j-th partition
	// Make a vector of nParts vectors of K int* each
	 int* Indices = new int[nParts*K];
	 int ind_c = 0;
	 for (int i=0; i < nParts; i++){
		 for(int j = 0; j < K; j++){
			 Indices[i*K + j] = ind_c;
			 ind_c = ind_c + dimClasses[i*K + j];
		 }
	 }

	// do initial pruning on argParts and return the pruned partitions

	// Make a vector of nParts vectors of K int* each
	vector <vector <int*> > Parts(nParts,vector<int*>(K));
	ind_c = 0;
	int argParts_size=0;
	for (int i=0; i < nParts; i++){
		for(int j = 0; j < K; j++){
			Parts[i][j] = argParts + ind_c;
			ind_c = ind_c + dimClasses[i*K + j];
			argParts_size = argParts_size + dimClasses[i*K + j];
		}
	}

	// in the following we call initial_prune with Parts which is a vector. This is not the most
	// efficient since vector does not allow for direct addressing. But since initial_prune doesn't have very high complexity, and
	// the running time for 7 partitions with 288 classes per partition is a couple of minutes at most, i'll just leave it for now.....

	// comment out for testing
	Util::initial_prune(Parts, dimClasses, nParts, K, T);
	for(int i = 0; i < nParts; i++){
		for(int j=0; j < K; j++){
			argParts[Indices[i*K + j]+1] = -1;
		}
	}

	int num_classes;
	int old_index;
	for(int i=0; i<nParts; i++){
		num_classes = Parts[i].size();// number of classes in partition i after pruning
		for (int j=0; j < num_classes; j++){
			old_index = Parts[i][j][0];
			//cout << "old_index: " << old_index<<"\n";
			argParts[Indices[i*K + old_index]+1] = 1;
		}
	}


	// if we're not doing mpi then keep going and call branchMPI and return the output
	//cout <<"begin partition matching\n";
	//int* dummy(0);
	int* output = Util::branchMPI(argParts, Indices,dimClasses, nParts, K, T, 0, n_guesses, LARGEST_CLASS, J, max_branching, stmult, branchfunc, LIM);
	
	//cout<<"total cost: "<<*output<<"\n";
	//cout<<"number of matches: "<<*(output+1)<<"\n";
	// now go check if the matches are sensical! i.e, if the matches are feasible, if the sum of the match weights in output is equal to *output, and if each match in output has weight at least T
	bool correct = Util::sanitycheck(argParts, Indices,dimClasses, nParts, K, T,output);

	delete[] Indices;

	// something is wrong with output of branchMPI!
	if (correct < 1){
		cout << "something is wrong with output of branchMPI!\n";
		vector<int> ret(1);
		ret[0] = -1;
		if (output != 0)  { delete[] output; output = 0; }
		return ret;
	}

	// output is not nonsense, so now put it into a single dimension vector and return
	// output is an int array, the first element is the cost of the output solution, the second element is the total number of matches in the solution
	// and the rest is the list of matches. output is one dimensional

	int output_size = 2 + output[1] * nParts;
	vector<int> ret(output_size);
	for (int i = 0; i < output_size; i++) {
		ret[i]= output[i];
	}
	if (output != 0) { delete[] output; output = 0; }
	return ret;

}


int branch_all=0;
int* Util::branchMPI(int* argParts, int* Indices, int* dimClasses, int nParts, int K, int T, int curlevel,int n_guesses, int
LARGEST_CLASS, int J, int max_branching, float stmult, int branchfunc, int LIM) {

//*************************************
//testing search2
if (1 == 0){
cout <<"begin test search2\n";
int* matchlist = new int[J*nParts];
int* costlist = new int[J];
for (int jit = 0; jit < nParts; jit++) *(costlist+jit) = 0;
Util::search2(argParts,Indices, dimClasses, nParts, K,  T, matchlist,costlist,J);

for (int jit = 0; jit < J; jit++) {
  cout << *(costlist +jit)<<": ";
  for (int yit = 0; yit < nParts; yit++)
  	cout << *(matchlist + jit*nParts + yit)<<",";
  cout <<"\n";	

}
cout <<"end test search2\n";
int* output = new int[1];
output[0] = 1;
delete [] matchlist;
delete [] costlist;
return output;
}
//**************************************

	// Base Case: we're at a leaf, no more feasible matches possible
	if (curlevel > K -1){
		int* output = new int[2];
		output[0] = 0;
		output[1] = 0;
		return output;
	}

	// branch dynamically depending on results of search 2!
	
	int* matchlist = new int[J*nParts];
	int* costlist = new int[J];
	Util::search2(argParts, Indices, dimClasses, nParts, K,  T, matchlist, costlist, J);
	
	
	// each class in the matches found by findTopLargest is encoded by the original index of the first element of the class in argPart
	// each match contains nParts classes, with the i-th class belonging to the i-th partition.

	// if there are no feasible matches with cost gt T, then return 0
	for (int jit = 0; jit < J ; jit++){
	
		if (costlist[jit] > T) break;
		if (jit == J-1){
			int* output = new int[2];
			output[0] = 0;
			output[1] = 0;
			delete[] matchlist;
			delete[] costlist;
			return output;
		}
	}
	

	
	// note that costlist and matchlist are NOT sorted by weight, and branch factor takes care of that...
	if (curlevel==0) branch_all = 0;
	
	int nBranches = -1;

	if (branchfunc == 0)
		nBranches = branch_factor_0(costlist,matchlist,J, T, nParts, curlevel, max_branching, LIM); // branch based on distribution of top J (weighted) matches  with cost > T

	if (branchfunc == 2)
		nBranches = branch_factor_2(costlist,matchlist,J, T, nParts, curlevel, max_branching, LIM); // branch based on distribution of top J (weighted) matches  with cost > T

	if (branchfunc == 3)
		nBranches = branch_factor_3(costlist,matchlist,J, T, nParts, curlevel, max_branching, K, LIM); // branch based on distribution of top J (weighted) matches  with cost > T

	if (branchfunc == 4)
		nBranches = branch_factor_4(costlist,matchlist,J, T, nParts, curlevel, max_branching, stmult); // branch based on distribution of top J (weighted) matches  with cost > T

	int* newcostlist= new int[nBranches];
	int* newmatchlist = new int[nBranches*nParts];
	for (int i=0; i<nBranches; i++){
		newcostlist[i] = costlist[i];
		for (int j=0; j< nParts; j++)
			newmatchlist[i*nParts + j] = matchlist[i*nParts + j];
	}

	delete[] costlist;
	delete[] matchlist;
	
	//int* output = new int[2];//initialize to placeholder
	int* output = new int[2+K*nParts];//initialize to placeholder
	output[0] = 0;
	output[1] = 0;
	// some temporary variables
	int old_index;
	int totalcost;
	int nmatches;
	//int offset;

	for(int i=0; i < nBranches ; i++){

		// consider the i-th match returned by findTopLargest
		//if (newcostlist[i] <= T) continue;

		// mark the classes in the i-th match of matchlist as taken (using the dummy variable and -2), and then call branch again on argParts.
		// After branch returns, compute overall cost, unmark  the classes just marked as 1 again in preparation for next loop.

		for(int j=0; j < nParts; j++){
			// matchlist[i*nParts + j] is the original index of the class belonging to the j-th partition in the i-th match.
			old_index = newmatchlist[i*nParts + j];
			argParts[Indices[j*K+old_index] + 1] = -2;
		}

		
		int* ret = Util::branchMPI(argParts, Indices, dimClasses, nParts, K, T, curlevel+1, n_guesses, LARGEST_CLASS,
		J, max_branching, stmult,branchfunc, LIM);
		
		// first element of ret is the total cost of all the matches in ret, and second element is the number of matches in ret
		totalcost = newcostlist[i] + ret[0];

		//if (curlevel == 0) {
		//	cout <<"totalcost*****************************************************************: "<<totalcost<<", costlist["<<i<<"]="<<newcostlist[i]<<", *ret="<<*ret<<", level: "<<curlevel<<"\n";
			
		//}
		if (totalcost > output[0]) // option 1
		{
			nmatches = 1 + ret[1];
			//delete[] output; // get rid of the old maxreturn
			//output = new int[2+nmatches*nParts];
			output[0] = totalcost;
			output[1] = nmatches;
			int nret = 2+(nmatches-1)*nParts;
			for(int iret=2; iret < nret; iret++) output[iret] = ret[iret];
			for(int imax=0; imax < nParts; imax++) output[nret+imax] = newmatchlist[i*nParts + imax];
		}


		delete[] ret;

		// unmark the marked classes in preparation for the next iteration

		for(int j=0; j < nParts; j++){
			old_index = newmatchlist[i*nParts + j];
			argParts[Indices[j*K+old_index] + 1] = 1;
		}

	}

	delete[] newmatchlist;
	delete[] newcostlist;
	
	return output;
}

int* costlist_global;
// make global costlist
bool jiafunc(int i, int j){
	return (costlist_global[j] < costlist_global[i]) ;

}
// Given J matches, branch always on the first one (i.e., the one with the largest weight, so the worst case we just end up doing greedy).
// Branch on the second one only if it is INFEASIBLE with the first, so you know it will never occur in any branching beginning with the first.
// Branch on subsequent ones only if its infeasible with ALL the ones which we have previously decided to branch on.
// The other option is to use LIM - so we branch on a match if its infeasible with at least LIM matches which we have previously decoded to branch on.
// For now, LIM is defaulted to -1, which means we branch on a match only if it is infeasible to ALL matches we have previously decided to branch on.
int Util::branch_factor_2(int* costlist, int* matchlist, int J, int T, int nParts, int, int max_branching, int LIM){
	
	int ntot=0;
	for (int jit=0; jit < J; jit++){
		if (*(costlist+jit) > T) ntot++;
	}

	int cur;
	// sort matchlist by cost
	int* indx = new int[J];
	for (int jit=0; jit < J; jit++) indx[jit]=jit;
	vector<int> myindx (indx, indx+J);
	vector<int>::iterator it;
	costlist_global=costlist;
	sort(myindx.begin(), myindx.end(), jiafunc);

	// put matchlist in the order of mycost
	int* templist = new int[J];
	int* temp2list = new int[J*nParts];
	int next = 0;
	
	for (it=myindx.begin(); it!=myindx.end();++it){
		cur = *(costlist + *it);
		if (cur > T){
			
			templist[next] = cur;
			for (int vt = 0; vt < nParts; vt++) temp2list[next*nParts + vt] = matchlist[ (*it)*nParts + vt];
			next = next + 1;
		}
	}
	
	for (int jit=0; jit < ntot; jit++){
		*(costlist+jit)=*(templist + jit);
		//cout <<*(costlist+jit)<<", ";
		for (int vit=0; vit < nParts; vit++) matchlist[jit*nParts + vit]= temp2list[jit*nParts + vit];
	}
	//cout <<"\n";
	
	delete [] indx;
	//compute the average 
	
	// partition matches into groups where matches in each group is mutually feasible, and matches from different groups are not feasible
	// The number of branches to explore are the mutually exclusive matches. Take the largest of these
	
	
	int B=1;
	int B_init=B;
	int infeasible=0;
	
	for (int i=B_init; i<ntot; i++){
		if (i==ntot) continue;
		// look at the i-th match. branch on it only if it is infeasible with ALL previous matches we have decided to
		// branch on
		infeasible = 0;
		if (LIM < 0) LIM = B;
		for (int j=0; j<B; j++){
			
			for (int vit=0; vit<nParts; vit++){
				if (temp2list[i*nParts+vit] == matchlist[j*nParts+vit]) {infeasible++; break;}
			}
			if (infeasible >= LIM) break;
		}
		
		if (infeasible >= LIM){
			*(costlist+B)=*(templist+i);
			for (int vit=0; vit < nParts; vit++)
				*(matchlist+B*nParts + vit)=*(temp2list+i*nParts + vit);
			B=B+1;	
		}
	}
	
	delete [] templist;
	delete [] temp2list;
	//cout<<"**************************************** "<<B<<" ***************************\n";
	
	if (branch_all < max_branching){
		if (B>1)
			{branch_all = branch_all + B -1 ; }
	}
	else B=1;
	
	return B;
	

}


// similar to branch_factor_2 except we branch on a match if it is infeasible with all other matches in matchlist (not just the ones we branch on). LIM plays similar role here.
int Util::branch_factor_3(int* costlist, int* matchlist, int J, int T, int nParts, int, int max_branching, int, int LIM){
	
	int ntot=0;
	for (int jit=0; jit < J; jit++){
		if (*(costlist+jit) > T) ntot++;
	}

	int cur;
	// sort matchlist by cost
	int* indx = new int[J];
	for (int jit=0; jit < J; jit++) indx[jit]=jit;
	vector<int> myindx (indx, indx+J);
	vector<int>::iterator it;
	costlist_global=costlist;
	sort(myindx.begin(), myindx.begin()+J, jiafunc);

	// put matchlist in the order of mycost
	int* templist = new int[J];
	int* temp2list = new int[J*nParts];
	int next = 0;
	
	for (it=myindx.begin(); it!=myindx.end();++it){
		cur = *(costlist + *it);
		if (cur > T){
			
			templist[next] = cur;
			for (int vt = 0; vt < nParts; vt++) temp2list[next*nParts + vt] = matchlist[ (*it)*nParts + vt];
			next = next + 1;
		}
	}
	
	for (int jit=0; jit < ntot; jit++){
		*(costlist+jit)=*(templist + jit);
		//cout <<*(costlist+jit)<<", ";
		for (int vit=0; vit < nParts; vit++) matchlist[jit*nParts + vit]= temp2list[jit*nParts + vit];
	}
	//cout <<"\n";
	
	delete [] indx;
	//compute the average 
	
	// partition matches into groups where matches in each group is mutually feasible, and matches from different groups are not feasible
	// The number of branches to explore are the mutually exclusive matches. Take the largest of these
	
	
	int B=1;
	int B_init=B;
	int infeasible=0;
	// if we're near the bottom of the tree then explore more... this is because the larger weights are not likely to change much,
	// whereas the smaller ones can have many permutations
	if (LIM < 0) LIM = ntot-1;
	for (int i=B_init; i<ntot; i++){
		if (i==ntot) continue;
		// look at the i-th match. branch on it only if it is infeasible with ALL previous matches we have decided to
		// branch on
		infeasible = 0;
		
		for (int j=0; j<ntot; j++){
			if (j == i) continue;
			for (int vit=0; vit<nParts; vit++){
				if (temp2list[i*nParts+vit] == temp2list[j*nParts+vit]) {infeasible++; break;}
			}
			if (infeasible >= LIM) break;
		}
		
		if (infeasible >= LIM){
			*(costlist+B)=*(templist+i);
			for (int vit=0; vit < nParts; vit++)
				*(matchlist+B*nParts + vit)=*(temp2list+i*nParts + vit);
			B=B+1;	
		}
	}
	
	delete [] templist;
	delete [] temp2list;
	//cout<<"**************************************** "<<B<<" ***************************\n";
	
	
	if (branch_all < max_branching){
		if (B>1)
			{branch_all = branch_all + B-1;}
	}
	else B=1;
	
	return B;
	

}

// We branch based on distribution of the cost of the J largest matches. Roughly speaking, if there is a match which has significantly larger weight than others, then we branch just on that
// match. Otherwise, we branch on similar weighted matches.
// As before we always branch on the match with the largest cost so worst case we'll get greedy.
// We compute standard dev of the J costs, and then if the difference between the cost of a match and the largest cost is within stmult*standard dev, then we branch on it.
int Util::branch_factor_4(int* costlist, int* matchlist, int J, int T, int nParts, int, int max_branching, float stmult){
	int sum=0;
	float average =0;
	int ntot=0;
	for (int jit=0; jit < J; jit++){
		if (*(costlist+jit) > T) {ntot++; sum = sum +*(costlist+jit);}
	}
	average = ((float)sum)/((float)ntot);
	int cur;
	// sort matchlist by cost
	int* indx = new int[J];
	for (int jit=0; jit < J; jit++) indx[jit]=jit;
	vector<int> myindx (indx, indx+J);
	vector<int>::iterator it;
	costlist_global=costlist;
	sort(myindx.begin(), myindx.begin()+J, jiafunc);

	// put matchlist in the order of mycost
	int* templist = new int[J];
	int* temp2list = new int[J*nParts];
	int next = 0;
	
	for (it=myindx.begin(); it!=myindx.end();++it){
		cur = *(costlist + *it);
		if (cur > T){
			
			templist[next] = cur;
			for (int vt = 0; vt < nParts; vt++) temp2list[next*nParts + vt] = matchlist[ (*it)*nParts + vt];
			next = next + 1;
		}
	}
	
	for (int jit=0; jit < ntot; jit++){
		*(costlist+jit)=*(templist + jit);
		//cout <<*(costlist+jit)<<", ";
		for (int vit=0; vit < nParts; vit++) matchlist[jit*nParts + vit]= temp2list[jit*nParts + vit];
	}
	//cout <<"\n";
	
	delete [] indx;
	delete [] templist;
	delete [] temp2list;
	
	if (ntot == 1) return 1;
	
	// look at the average, standard dev etc. If standard dev very small, i.e., costs very similar, then branch on the similar
	// costs
	float sq_sum=0.0;
	//cout <<"costlist:";
	for (int i=0; i< ntot; i++){
		sq_sum = sq_sum + (float) pow((float) *(costlist+i) - average, (float)2.0);
		//cout <<*(costlist+i)<<", ";
	}	
	//cout <<"\n";
	
	float variance = sq_sum/ntot;
	float stdev = (float)pow((float)variance,(float)0.5);
	
	//cout <<"stdev: "<<int(stdev)<<"\n";
	
	int B=1;
	int largest = *costlist;
	//cout <<"largest: "<<largest<<"\n";
	for (int i=1; i<ntot; i++){
		int cur = *(costlist+i);
		if (largest-cur < (float)(stdev*stmult)) B++;
		else break;
	
	}
	//cout <<"B: "<<B<<"\n";
	if (branch_all < max_branching){
		if (B>1)
			{branch_all = branch_all + B-1;}
	}
	else B=1;
	
	return B;
	

}

int Util::branch_factor_0(int* costlist, int* matchlist, int J, int T, int nParts, int, int max_branching, int LIM){
	
	int ntot=0;
	for (int jit=0; jit < J; jit++){
		if (*(costlist+jit) > T) ntot++;
	}

	int cur;
	// sort matchlist by cost
	int* indx = new int[J];
	for (int jit=0; jit < J; jit++) indx[jit]=jit;
	vector<int> myindx (indx, indx+J);
	vector<int>::iterator it;
	costlist_global=costlist;
	sort(myindx.begin(), myindx.begin()+J, jiafunc);

	// put matchlist in the order of mycost
	int* templist = new int[J];
	int* temp2list = new int[J*nParts];
	int next = 0;
	
	for (it=myindx.begin(); it!=myindx.end();++it){
		cur = *(costlist + *it);
		if (cur > T){
			
			templist[next] = cur;
			for (int vt = 0; vt < nParts; vt++) temp2list[next*nParts + vt] = matchlist[ (*it)*nParts + vt];
			next = next + 1;
		}
	}
	
	for (int jit=0; jit < ntot; jit++){
		*(costlist+jit)=*(templist + jit);
		//cout <<*(costlist+jit)<<", ";
		for (int vit=0; vit < nParts; vit++) matchlist[jit*nParts + vit]= temp2list[jit*nParts + vit];
	}
	//cout <<"\n";
	
	for (int jit=1; jit < ntot; jit++){
	
	     if ((costlist[jit] == costlist[0]) && costlist[jit] > T){
	     
	             for (int vit=0; vit < nParts; vit++){
		             if ( matchlist[jit*nParts + vit] >  matchlist[vit])
			     	 break;
		             if ( matchlist[jit*nParts + vit] ==  matchlist[vit])
			     	 continue;
			     if ( matchlist[jit*nParts + vit] <  matchlist[vit])
			     {
			         // swap
				 for (int swp=0; swp < nParts; swp++){
				       int tmp  = matchlist[swp];
				       matchlist[swp]= matchlist[jit*nParts + swp];
				       matchlist[jit*nParts + swp] = tmp;
				 }
				 break;
			     
			     }	 
		     }
	     }
	
	}
	
	
	delete [] indx;
	//compute the average 
	
	// partition matches into groups where matches in each group is mutually feasible, and matches from different groups are not feasible
	// The number of branches to explore are the mutually exclusive matches. Take the largest of these
	
	
	int B=1;
	int B_init=B;
	int infeasible=0;
	
	for (int i=B_init; i<ntot; i++){
		if (i==ntot) continue;
		// look at the i-th match. branch on it only if it is infeasible with ALL previous matches we have decided to
		// branch on
		infeasible = 0;
		if (LIM < 0) LIM = B;
		for (int j=0; j<B; j++){
			
			for (int vit=0; vit<nParts; vit++){
				if (temp2list[i*nParts+vit] == matchlist[j*nParts+vit]) {infeasible++; break;}
			}
			if (infeasible >= LIM) break;
		}
		
		if (infeasible >= LIM){
			*(costlist+B)=*(templist+i);
			for (int vit=0; vit < nParts; vit++)
				*(matchlist+B*nParts + vit)=*(temp2list+i*nParts + vit);
			B=B+1;	
		}
	}
	
	delete [] templist;
	delete [] temp2list;
	//cout<<"**************************************** "<<B<<" ***************************\n";
	
	if (branch_all < max_branching){
		if (B>1)
			{branch_all = branch_all + B -1 ; }
	}
	else B=1;
	
	return B;

}


Dict Util::get_transform_params(EMData* image, string xform, string convention)
{
	const Transform * transform = image->get_attr(xform);
	Dict transform_params = transform->get_params(convention);
	delete transform;
	return transform_params;
}


static void compose_transform2(float psi1, float sx1, float sy1, float psi2, float sx2, float sy2, float & out_psi, float & out_sx, float & out_sy)
{
	Dict t_params;
	t_params["type"]   = "2D";
	t_params["alpha"]  = psi1;
    t_params["tx"]     = sx1;
    t_params["ty"]     = sy1;
    t_params["mirror"] = 0;
    t_params["scale"]  = 1.0;
	Transform t1(t_params);
	t_params["alpha"]  = psi2;
    t_params["tx"]     = sx2;
    t_params["ty"]     = sy2;
	Transform t2(t_params);
	Transform tt = t2*t1;
	Dict d = tt.get_params("2D");
	out_psi = d["alpha"];
	out_sx  = d["tx"];
	out_sy  = d["ty"];
}


void Util::constrained_helix_exhaustive( vector<EMData*> data, vector<EMData*> fdata, vector<EMData*> refproj, vector<EMData*> rotproj,
			vector<float> dp_dphi_rise_delta, vector<int> nphi_phiwobble_range_ywobble_Dsym_nwx_nwy_nwxc_nwyc,
			bool FindPsi, float psi_max, vector<EMData*> crefim, vector<int> numr, int maxrin, string mode, int cnx, int cny)
{
	if (dp_dphi_rise_delta.size() < 4) {
		printf("Not enough parameters (dp_dphi_rise_delta)");
		return;
	}
	if (nphi_phiwobble_range_ywobble_Dsym_nwx_nwy_nwxc_nwyc.size() < 9) {
		printf("Not enough parameters (nphi_phiwobble_range_ywobble_Dsym_nwx_nwy_nwxc_nwyc)");
		return;
	}
	float dpsi;
	//float dp       = dp_dphi_rise_delta[0];
	float dphi        = dp_dphi_rise_delta[1];
	float rise        = dp_dphi_rise_delta[2];  // Float in pixels!
	float delta       = dp_dphi_rise_delta[3];
	float ywobble     = dp_dphi_rise_delta[4];
	float ystep       = dp_dphi_rise_delta[5];
	int   nphi        = nphi_phiwobble_range_ywobble_Dsym_nwx_nwy_nwxc_nwyc[0];
	int   symrestrict = nphi_phiwobble_range_ywobble_Dsym_nwx_nwy_nwxc_nwyc[1];
	int   phiwobble   = nphi_phiwobble_range_ywobble_Dsym_nwx_nwy_nwxc_nwyc[2];
	int   range       = nphi_phiwobble_range_ywobble_Dsym_nwx_nwy_nwxc_nwyc[3];
	bool  Dsym        = nphi_phiwobble_range_ywobble_Dsym_nwx_nwy_nwxc_nwyc[4];
	int   nwx         = nphi_phiwobble_range_ywobble_Dsym_nwx_nwy_nwxc_nwyc[5];
	int   nwy         = nphi_phiwobble_range_ywobble_Dsym_nwx_nwy_nwxc_nwyc[6];
	int   nwxc        = nphi_phiwobble_range_ywobble_Dsym_nwx_nwy_nwxc_nwyc[7];
	int   nwyc        = nphi_phiwobble_range_ywobble_Dsym_nwx_nwy_nwxc_nwyc[8];

	const int ndata = data.size();

	std::vector< boost::shared_ptr<EMData> > objectsToDelete; // objects added to this vector are automatically deleted at the end of this function

	vector<float> c0 = data[0]->get_attr("ptcl_source_coord");
	vector< vector<EMData*> > ccfs(ndata, vector<EMData*>(nphi));
	vector< vector<EMData*> > ccfr(0);
	if (! Dsym)  ccfr.resize(ndata, vector<EMData*>(nphi));

	for (int im = 0; im < ndata; ++im) {
		for (int iphi = 0; iphi < nphi; ++iphi) {
			std::auto_ptr<EMData> corr( correlation( refproj[iphi], fdata[im], CIRCULANT, true) );
			ccfs[im][iphi] = Util::window( corr.get(), nwx, nwy);
			objectsToDelete.push_back(boost::shared_ptr<EMData>(ccfs[im][iphi]));
			if (! Dsym) {
				std::auto_ptr<EMData> corr2( correlation( rotproj[iphi], fdata[im], CIRCULANT, true) );
				ccfr[im][iphi] = Util::window( corr2.get(), nwx, nwy);
				objectsToDelete.push_back(boost::shared_ptr<EMData>(ccfr[im][iphi]));
			}
		}
	}

	vector<float> dxshiftlocal(ndata, 0);
	vector<float> dyshiftlocal(ndata, 0);
	vector<float> dphilocal(ndata, 0);

	vector<float> xshiftlocal(ndata, 0);
	vector<float> xrshiftlocal(ndata, 0);
	vector<float> mxshiftlocal(ndata, 0);
	vector<float> yshiftlocal(ndata, 0);
	vector<float> yrshiftlocal(ndata, 0);
	vector<float> myshiftlocal(ndata, 0);
	vector<float> philocal(ndata, 0);
	vector<float> phirlocal(ndata, 0);
	vector<float> mphilocal(ndata, 0);
	float dirma = -1.0e23f;
	for (int idir = -1; idir < 2; idir += 2) {
		float tmax = -1.0e23f;
		float mpsi;
		for (int ix = 1; ix < nwx-1; ++ix) {                                         //#  X shift
			//#cout << "im: ", len(ccfs), ix,time()-start_time
			int six = ix - nwxc;
			float siy = -rise/2.0f - ystep + ywobble;
			//for (int iy = 1+ywobble; iy < nwy-ywobble-1; ++iy)
			while( siy < rise/2.0 - ystep - ywobble ) {                     //#  Y shift
				siy += ystep;
				yshiftlocal[0]  = siy;
				yrshiftlocal[0] = siy;

				float yiy = siy + nwyc;
				int   fiy = int(yiy);
				float ydif = yiy - fiy;
				float yrem = 1.0 - ydif;
				int   iuy = fiy;

				float xrem = 1.0f;
				float xdif = 0.0f;


				for (int iphi = 0; iphi < nphi/symrestrict; ++iphi) {               //#  phi search, first segment has phi range restricted by point-group symmetry
					float qphi = iphi*delta;
					philocal[0]  = qphi;
					qphi = 180.0f - qphi;
					if(  qphi > 0.0f )     phirlocal[0] = fmod( qphi , 360.0f );
					else if( qphi < 0.0f ) phirlocal[0] = 360.0f + fmod( qphi , 360.0f );
					//# we use the first segment as a reference, so there is no interpolation, just copy the correlation
					//#  Select largest correlation within +/- range pixels from the location we explore
					float mxm = -1.0e23f;
					float mxr;
					for (int iux = max(1, ix - range); iux < min(nwx - 1, ix+range+1); ++iux) {    //#  X wobble
						//float qcf = ccfs[0][iphi]->get_value_at(iux,iy);
						float qcf =   xrem*yrem*ccfs[0][iphi]->get_value_at(iux,iuy)
									+ xdif*yrem*ccfs[0][iphi]->get_value_at(iux+1,iuy)
									+ xrem*ydif*ccfs[0][iphi]->get_value_at(iux,iuy+1)
									+ xdif*ydif*ccfs[0][iphi]->get_value_at(iux+1,iuy+1);
						if (qcf > mxm) {
							mxm = qcf;
							xshiftlocal[0] = float(iux-nwxc);
						}
					}
					if (! Dsym) {
						mxr = -1.0e23f;
						for (int iux = max(1, ix - range); iux < min(nwx - 1, ix+range+1); ++iux) {     //# Xr wobble
							//float qcf = ccfr[0][iphi]->get_value_at(iux,iy);
							float qcf =   xrem*yrem*ccfr[0][iphi]->get_value_at(iux,iuy)
										+ xdif*yrem*ccfr[0][iphi]->get_value_at(iux+1,iuy)
										+ xrem*ydif*ccfr[0][iphi]->get_value_at(iux,iuy+1)
										+ xdif*ydif*ccfr[0][iphi]->get_value_at(iux+1,iuy+1);
							if (qcf > mxr) {
								mxr = qcf;
								xrshiftlocal[0] = float(iux-nwxc);
							}
						}
					}
					//cout << "  QST "<<iphi<< "  "<<mxm<< "  "<<mxr<<endl;
					for ( int im = 1; im < ndata; ++im) {                                       //#  predicted locations
						//# dst is distance between segment 0 and current segment in pixels
						vector<float> cim = data[im]->get_attr("ptcl_source_coord");
						float dst = sqrt( (c0[0] - cim[0])*(c0[0] - cim[0]) + (c0[1] - cim[1])*(c0[1] - cim[1]));
						//# predict for all remaining segments assuming number 0
						float qd  = round((siy + dst)/rise);
						float piy = siy + dst - rise*qd;                                    //#  predicted y shift
						//cout<<siy<<"  "<<dst<<"  "<<rise<<"  "<<qd<<"  "<<piy<<"  "<<siy<<"  "<<siy<<"  "<<endl;
						float pphi = (philocal[0] + idir*dphi*qd);
						if(pphi >= 360.0f)     pphi = fmod( pphi, 360.0f);
						else if( pphi < 0.0f)  pphi = 360.0f + fmod( pphi, 360.0f);         //#  predicted phi with full angular accuracy, not an integer
						int pix = six;                                                      //# predicted x shift

						int xix = pix + nwxc;
						//#  Local x search
						int fix = int(xix);
						float xdif = xix - fix;
						float xrem = 1.0f - xdif;
						float ciq = -1.0e23f;
						//# interpolate correlation at pphi
						int ttphi = (int( pphi/delta + 0.5))%nphi;
						for (int lphi = -phiwobble; lphi < phiwobble+1; ++lphi) {                                               //#  phi wobble
							int tphi = (ttphi+lphi+nphi)%nphi;
							for (int iux = max(1, fix - range); iux < min(nwx - 1, fix+range+1); ++iux) {                       //#  X wobble
								//for (int iuy = max(1, fiy - ywobble); iuy < min(nwy - 1, fiy+ywobble+1); ++iuy) {               //#  Y wobble
								float tiy = -ywobble - ystep + piy + nwyc;
								while( tiy <= ywobble - ystep + piy + nwyc )  {
									tiy += ystep;
									int iuy = int(tiy);
									if(iuy > -1 && iuy < nwy-1)  {
									
										float ydif = yiy - fiy;
										float yrem = 1.0f - ydif;


										float qcf =   xrem*yrem*ccfs[im][tphi]->get_value_at(iux,iuy)
													+ xdif*yrem*ccfs[im][tphi]->get_value_at(iux+1,iuy)
													+ xrem*ydif*ccfs[im][tphi]->get_value_at(iux,iuy+1)
													+ xdif*ydif*ccfs[im][tphi]->get_value_at(iux+1,iuy+1);
				//cout << lphi<< "  "<< tphi<<  "  "<<iux<< "  "<<iuy<< "  "<<fiy<< "  "<<ywobble<< "  "<<nwy<< "  "<<qcf<<endl;
										if (qcf > ciq) {
											ciq = qcf;
											xshiftlocal[im] = iux + xdif - nwxc;
											yshiftlocal[im] = tiy - nwyc;
											philocal[im]    = tphi * delta;
										}
									}
								}
							}
						}
						//cout << "  St "<<iphi<< "  "<<mxm<< "  "<<ciq<< "  "<<mxm+ciq<< "  "<<endl;
						//cout << im<< "  "<< ttphi<<  "  "<<phiwobble<< "  "<<ywobble<< "  "<<fiy<<endl;
						mxm += ciq;
						//# now for rotated
						if (! Dsym) {
							pphi = (phirlocal[0] + idir*dphi*qd);                               //#  predicted phi for rotated 180 defs with full angular accuracy, not an integer
							if(pphi >= 360.0f)  pphi = fmod( pphi, 360.0f);
							else if( pphi < 0.0f)  pphi = 360.0f + fmod( pphi, 360.0f);         //#  predicted phi with full angular accuracy, not an integer
							ciq = -1.0e23f;
							//# interpolate correlation at pphi
							for (int lphi = -phiwobble; lphi < phiwobble+1; ++lphi) {                                           //#  phi wobble
								qphi = lphi*delta;
								float qtphi = fmod( pphi + qphi + (int(fabs(qphi/360.0f))+1)*360.0f , 360.0f);
								qphi = fmod(540.0f-qtphi, 360.0f) / delta;
								int tphi = (int( qphi + (int(fabs(qphi/nphi))+1)*nphi  + 0.5))%nphi;
								for (int iux = max(1, fix - range); iux < min(nwx - 1, fix+range+1); ++iux) {                   //#  X wobble
									//for (int iuy = max(1, fiy - ywobble); iuy < min(nwy - 1, fiy+ywobble+1); ++iuy) {           //#  Y wobble

									float tiy = -ywobble - ystep + piy + nwyc;
									while( tiy <= ywobble - ystep + piy + nwyc )  {
										tiy += ystep;
										int iuy = int(tiy);
										if(iuy > -1 && iuy < nwy-1)  {
								
											float ydif = yiy - fiy;
											float yrem = 1.0f - ydif;



											float qcf = xrem*yrem*ccfr[im][tphi]->get_value_at(iux,iuy)
														+ xdif*yrem*ccfr[im][tphi]->get_value_at(iux+1,iuy)
														+ xrem*ydif*ccfr[im][tphi]->get_value_at(iux,iuy+1)
														+ xdif*ydif*ccfr[im][tphi]->get_value_at(iux+1,iuy+1);
											if (qcf > ciq) {
												ciq = qcf;
												xrshiftlocal[im] = iux + xdif - nwxc;
												yrshiftlocal[im] = tiy - nwyc;
												phirlocal[im]    = int(qtphi/delta+0.5f)*delta;
											}
										}
									}
								}
							}
							mxr += ciq;
						} else {
							mxr = mxm-1.e5;
						}
					}
					if ( mxr > mxm ) {
						if (mxr > tmax) {
							tmax = mxr;
							mpsi = 270.0f;
							for (int im = 0; im < ndata; ++im) mxshiftlocal[im] = xrshiftlocal[im];
							for (int im = 0; im < ndata; ++im) myshiftlocal[im] = yrshiftlocal[im];
							for (int im = 0; im < ndata; ++im) mphilocal[im]    = fmod(540.0f-phirlocal[im], 360.0f);
						}
					} else {
						if (mxm > tmax) {
							tmax = mxm;
							mpsi = 90.0f;
							for (int im = 0; im < ndata; ++im) mxshiftlocal[im] = xshiftlocal[im];
							for (int im = 0; im < ndata; ++im) myshiftlocal[im] = yshiftlocal[im];
							for (int im = 0; im < ndata; ++im) mphilocal[im]    = philocal[im];
						}
					}
				}
			}
		}
		//cout << idir<< "  "<< tmax<< endl;
		if (tmax > dirma) {
			dirma = tmax;
			dpsi = mpsi;
			for (int im = 0; im < ndata; ++im) dxshiftlocal[im] = mxshiftlocal[im];
			for (int im = 0; im < ndata; ++im) dyshiftlocal[im] = myshiftlocal[im];
			for (int im = 0; im < ndata; ++im) dphilocal[im]    = mphilocal[im];
		}
	}

	for (int im = 0; im < ndata; ++im) {
		float psx  = dxshiftlocal[im];
		float psy  = dyshiftlocal[im];
		float pphi = dphilocal[im];
		float epsi;
		float bestang;
		if (FindPsi) {
			float qphi = pphi/delta;
			int iphi = ( int(qphi + ((int)(abs(qphi/nphi))+1)*nphi + 0.5f))%nphi ;
			//#cout <<  " ref number and current parameters reduced to 2D  ",iphi,0.0, psx, psy
			Dict params = Util::get_transform_params(data[im], "xform.projection", "spider");
			float opsi3 = params["psi"];
			float opx3 = -static_cast<float>(params["tx"]);
			float opy3 = -static_cast<float>(params["ty"]);
			//#cout << " old 3D params in data ",ophi, otheta, opsi3, opx3, opy3
			float gamma = (abs(opsi3 - 90) < abs(opsi3 - 270)) ? (90) : (270);
			float oalpha, osx, osy;
			compose_transform2(0, opx3, opy3, gamma-opsi3, 0, 0, oalpha, osx, osy); //# reduce 3D to 2D
			//#cout << " old 3D params, -> 2D ",oalpha, osx, osy
			//# combine previous with the current in plane
			//#cout << " current 2D combined with old 2D rotation",oalpha, csx, csy
			//#  Find what the shift is without the angle
			float junk, nnsx, nnsy;
			compose_transform2(0, psx, psy, -oalpha, 0, 0, junk, nnsx, nnsy);
			//#cout << " 2D shift without angle ",nnsx, nnsy
			
			//this part below will be changed on a higher level later for Polar2Dm.     2-26-2015@ming
			std::auto_ptr<EMData> cimage( Util::Polar2Dm(data[im], cnx+nnsx, cny+nnsy, numr, mode) );
			Util::Frngs(cimage.get(), numr);
			std::auto_ptr<EMData> temp( Util::Crosrng_msg_s( cimage.get(), crefim[iphi], numr) );

			int ipr = int(psi_max*maxrin/360.0f + 0.5f);
			int incpsi = (dpsi == 270.0f) ? (maxrin/2) : (0);
			float qn = -1.0e23f;
			for (int ips = -ipr; ips < ipr+1; ++ips) {
				int tot = (ips + incpsi + maxrin)%maxrin;
				float tval = temp->get_value_at(tot);
				if (tval > qn) {
					qn = tval;
					bestang = ang_n(tot+1.0f, mode, maxrin);
				}
			}
			//#cout << " best angle ",bestang
			bestang = fmod(bestang - (dpsi-90.0f) + 720.0f, 360.0f);
			//#cout << " angle applied ",bestang
			//#rot_shift2D(data[im],-bestang).write_image("rotated.hdf",im)
			//std::auto_ptr<EMData> rot_data_im( data[im]->rot_scale_trans2D_background(-bestang, 0, 0, 1) );
			//fdata[im] = (rot_data_im->is_complex()) ? (rot_data_im->do_ift()) : (rot_data_im->do_fft());

			//#cout <<  " New composed 3D  ",dpsi,bestang, nnsx, nnsy

			epsi = fmod(bestang+dpsi, 360.0f);
			psx = nnsx;
			psy = nnsy;
			//#cout <<  " New composed 3D  ",pphi,90.0,epsi, psx, psy
		} else {
			epsi = dpsi;
			bestang = 0;
		}
		data[im]->set_attr("bestang",fmod(720.0f-bestang, 360.0f));
		//printf("  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f\n", pphi, 90.0, epsi, psx, psy);
		Dict t_params;
		t_params["type"]  = "spider";
		t_params["phi"]   = pphi;
        t_params["theta"] = 90.0f;
        t_params["psi"]   = epsi;
        t_params["tx"]    = -psx;
        t_params["ty"]    = -psy;
        Transform t(t_params);
        data[im]->set_attr("xform.projection", &t);
	}
}

// fast function to create transform (data in degrees)
inline Transform create_transform(float phi, float theta, float psi, float tx = 0, float ty = 0)
{
	std::vector<float> v(12);
	const float deg2rad = 3.141592653592 / 180.0;

	phi   *= deg2rad;
	theta *= deg2rad;
	psi   *= deg2rad;

	v[0] = cos(psi) * cos(theta) * cos(phi) - sin(psi) * sin(phi);
	v[1] = cos(psi) * cos(theta) * sin(phi) + sin(psi) * cos(phi);
	v[2] = - cos(psi) * sin(theta);
	v[3] = tx;

	v[4] = - sin(psi) * cos(theta) * cos(phi) - cos(psi) * sin(phi);
	v[5] = - sin(psi) * cos(theta) * sin(phi) + cos(psi) * cos(phi);
	v[6] = sin(psi) * sin(theta);
	v[7] = ty;

	v[8] = sin(theta) * cos(phi);
	v[9] = sin(theta) * sin(phi);
	v[10] = cos(theta);
	v[11] = 0;

	Transform t(v);
	return t;
}




// Calculate angles between param1 and param2 rotated by phi, theta, psi
// param1/2: (phi0, theta0, psi0, x0, y0, phi1, theta1, psi1, x1, y1, phi2, ... )
// returns n vales for n given projections
inline void diff_between_3D_parameters_angles( std::vector<float> & param1, std::vector<float> & param2, float phi, float theta, float psi, std::vector<float> & result)
{
	const unsigned n = param1.size() / 5;
	const float deg2rad = 3.141592653592 / 180.0;
	Transform rot = create_transform(phi, theta, psi);
	for ( unsigned i = 0;  i < n;  ++i ) {
		Transform temp = create_transform(param2[5*i+0], param2[5*i+1], param2[5*i+2], param2[5*i+3], param2[5*i+4]) * rot;
		Dict out = temp.get_params("spider");
		// calculate angle between positions on the sphere
		const float phi2   = static_cast<float>(out["phi"  ]) * deg2rad;
		const float theta2 = static_cast<float>(out["theta"]) * deg2rad;
		const float phi1   = param1[5*i+0] * deg2rad;
		const float theta1 = param1[5*i+1] * deg2rad;
		const float st1 = sin(theta1);
		const float st2 = sin(theta2);
		const float ct1 = cos(theta1);
		const float ct2 = cos(theta2);
		const float cp1cp2_sp1sp2 = cos(phi1 - phi2);
		float val = st1 * st2 * cp1cp2_sp1sp2 + ct1 * ct2;
		if (val < -1.0) {
			val = -1.0;
		} else if (val > 1.0) {
			val = 1.0;
		}
		result[i] += acos( val ) / deg2rad;
	}
}



// Calculate average angle error per projection for given sets of configurations and theirs pair-wise rotations
// all_params: [ phi0, theta0, psi0, sx0, sy0, phi1, theta1, psi1, sx1, sy1, phi2, .... ] - contains (5*number_of_projections) elements
// rotations: (phi_1_0, theta_1_0, psi_1_0, phi_2_0, theta_2_0, psi_2_0, phi_2_1, theta_2_1, psi_2_1, ... ) - contains (3*number_of_projections*(number_of_projections-1)/2) elements
// returns: (number_of_projections) elements - average angles' errors (per projection)
std::vector<float> Util::diff_between_matrix_of_3D_parameters_angles( std::vector<float> all_params, std::vector<float> rotations )
{
	// number of configurations (rotations is a bottom half of a square matrix)
	unsigned conf_count = 2;
	for ( const unsigned t = rotations.size() / 3;  conf_count * (conf_count-1) / 2 < t; ++conf_count );
	// number of projections
	const unsigned projs_count = all_params.size() / 5 / conf_count;

	// decode sets of parameters
	std::vector< std::vector<float> > params(conf_count);
	for (unsigned i = 0; i < conf_count; ++i) {
		params[i].assign( all_params.begin() + i*projs_count*5, all_params.begin() + (i+1)*projs_count*5 );
	}

	std::vector<float> avg(projs_count, 0);

	unsigned ir = 0;
	for (unsigned i = 0; i < conf_count; ++i) {
		for (unsigned j = 0; j < i; ++j) {
			diff_between_3D_parameters_angles(params[j], params[i], rotations[ir*3], rotations[ir*3+1], rotations[ir*3+2], avg);
			++ir;
		}
	}

	for (unsigned k = 0; k < projs_count; ++k) {
		avg[k] /= ( conf_count * (conf_count-1) / 2 );
	}

	return avg;
}


std::vector<int> Util::max_clique(std::vector<int> edges)
{
	const int size = *(max_element(edges.begin(), edges.end())) + 1;
	bool **conn = new bool*[size];
	for (int i=0; i < size; i++) {
		conn[i] = new bool[size];
		for (int j = 0; j < size; ++j) {
			conn[i][j] = false;
		}
	}
	for (unsigned i = 1; i < edges.size(); ++i) {
		int v1 = edges[i-1];
		int v2 = edges[i];
		conn[v1][v2] = conn[v2][v1] = true;
	}
	Maxclique m(conn, size);
	int *qmax;
	int qsize;
	m.mcq(qmax, qsize);
	std::vector<int> result(qsize);
	for (int i = 0;  i < qsize; ++i ) {
		result[i] = qmax[i];
	}
	delete [] qmax;
	for (int i=0;i<size;i++) {
		delete [] conn[i];
	}
	delete [] conn;
	return result;
}


#define img_ptr(i,j,k)  img_ptr[i+(j+(k*ny))*(size_t)nx]
#define img2_ptr(i,j,k) img2_ptr[i+(j+(k*ny))*(size_t)nx]
float Util::local_inner_product(EMData* image1, EMData* image2, int lx, int ly, int lz, int w) {
	int nx=image1->get_xsize(),ny=image1->get_ysize(),nz=image1->get_zsize();

	float *img_ptr  = image1->get_data();
	float *img2_ptr = image2->get_data();

	int lop = w/2;

	float lip  = 0.0f;
	float nrm1 = 0.0f;
	float nrm2 = 0.0f;

	for (int k=-lop; k<=lop; k++) {
		int kz = lz + k;
		for (int j=-lop; j<=lop; j++) {
			int ky = ly +j;
			for (int i=-lop; i<=lop; i++) {
				int kx = lx + i;
				nrm1 += img_ptr(kx, ky, kz)  * img_ptr(kx, ky, kz);
				nrm2 += img2_ptr(kx, ky, kz) * img2_ptr(kx, ky, kz);
				lip  += img_ptr(kx, ky, kz)  * img2_ptr(kx, ky, kz);
			}
		}
	}
	return lip/sqrt(nrm1*nrm2);
}
#undef img_ptr
#undef img2_ptr

