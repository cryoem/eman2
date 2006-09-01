/**
 * $Id$
 */
#include <iostream>
#include <stdio.h>

#include "emdata.h"
#include "util.h"

#include "lapackblas.h"

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_bessel.h>

using namespace EMAN;
using namespace std;

vector<float> Util::infomask(EMData* Vol, EMData* mask)
{
	ENTERFUNC;
	vector<float> stats;
	float *Volptr, *maskptr,MAX,MIN;
	long double Sum1,Sum2;	
	long count;
	
	MAX = FLT_MIN;
	MIN = FLT_MAX;
	count = 0L;
	Sum1 = 0.L;
	Sum2 = 0.L;
		
	
	if (mask == NULL)	
          {
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
	
	
	/* Calculation of the Statistics */
			       
	for (size_t i = 0;i < nx*ny*nz; i++)
	    {
	      if (maskptr[i]>0.5f)
	      {
	       Sum1 += Volptr[i];	       
	       Sum2 += Volptr[i]*Volptr[i];	       
	       MAX = (MAX < Volptr[i])?Volptr[i]:MAX;
	       MIN = (MIN > Volptr[i])?Volptr[i]:MIN;
	       count++;	       
	      }
	    }
	
       if (count==0) count++;
    
       float avg = static_cast<float>(Sum1/count);
       float sig2 = static_cast<float>(Sum2/count - avg*avg);
       float sig = sqrt(sig2);
                            
       stats.push_back(avg);
       stats.push_back(sig);
       stats.push_back(MIN);
       stats.push_back(MAX);

       return stats;
}
 
 
//---------------------------------------------------------------------------------------------------------- 
 
Dict Util::im_diff(EMData* V1, EMData* V2, EMData* mask)
{
	ENTERFUNC;
	
	size_t nx = V1->get_xsize();
	size_t ny = V1->get_ysize();
	size_t nz = V1->get_zsize();
	size_t size = nx*ny*nz;
	
	EMData *BD = new EMData();
 	BD->set_size(nx, ny, nz);
	
	float *params = new float[2];	 	 
  	
	float *V1ptr, *V2ptr, *MASKptr, *BDptr, A, B; 
	long double S1=0.L,S2=0.L,S3=0.L,S4=0.L;
	int nvox = 0L;
	
        V1ptr = V1->get_data();
	V2ptr = V2->get_data();
	MASKptr = mask->get_data();
	BDptr = BD->get_data();
	
	
//	 calculation of S1,S2,S3,S3,nvox
			       
	for (size_t i = 0L;i < size; i++) {
	      if (MASKptr[i]>0.5f) {
	       S1 += V1ptr[i]*V2ptr[i];
	       S2 += V2ptr[i]*V2ptr[i];
	       S3 += V2ptr[i]; 
	       S4 += V1ptr[i];
	       nvox ++;
	      }
	}       
	 
			
	A = static_cast<float> (nvox*S1 - S3*S4)/(nvox*S2 - S3*S3);
	B = static_cast<float> (A*S3  -  S4)/nvox;
        
	// calculation of the difference image
	
	for (size_t i = 0L;i < size; i++) {
	     if (MASKptr[i]>0.5f) {
	       BDptr[i] = A*V2ptr[i] -  B  - V1ptr[i];
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
		
	return BDnParams;
 }

//----------------------------------------------------------------------------------------------------------



EMData *Util::TwoDTestFunc(int Size, float p, float q,  float a, float b, int flag, float alphaDeg) //PRB
{
    int Mid= (Size+1)/2;

    if (flag<0 || flag>4) {
    	cout <<" flat must be 0,1,2,3, or 4";
    }
    if (flag==0) { // This is the real function
	   EMData* ImBW = new EMData();
	   ImBW->set_size(Size,Size,1);
	   ImBW->to_zero();
	   
	   float tempIm;
	   float x,y;
	
	   for (int ix=(1-Mid);  ix<Mid; ix++){
	        for (int iy=(1-Mid);  iy<Mid; iy++){
		  x = ix;
		  y = iy;
	       	  tempIm= (1/(2*M_PI)) * cos(p*x)* cos(q*y) * exp(-.5*x*x/(a*a))* exp(-.5*y*y/(b*b)) ;
		  (*ImBW)(ix+Mid-1,iy+Mid-1) = tempIm * exp(.5*p*p*a*a)* exp(.5*q*q*b*b);
	   	}
	   }
	   ImBW->done_data();
	   ImBW->set_complex(false);
	   ImBW->set_ri(true);
	
	
	   return ImBW;
   	}
   	if (flag==1) {  // This is the Fourier Transform
	   EMData* ImBWFFT = new EMData();
	   ImBWFFT ->set_size(2*Size,Size,1);
	   ImBWFFT ->to_zero();
	   
	   float r,s;
	
	   for (int ir=(1-Mid);  ir<Mid; ir++){
	        for (int is=(1-Mid);  is<Mid; is++){
		   r = ir;
		   s = is;
	       	   (*ImBWFFT)(2*(ir+Mid-1),is+Mid-1)= cosh(p*r*a*a) * cosh(q*s*b*b) *
		            exp(-.5*r*r*a*a)* exp(-.5*s*s*b*b);
	   	}
	   }
	   ImBWFFT->done_data();
	   ImBWFFT->set_complex(true);
	   ImBWFFT->set_ri(true);
	   ImBWFFT->set_shuffled(true);
	   ImBWFFT->set_fftodd(true);
	
	   return ImBWFFT;
   	}   		
   	if (flag==2 || flag==3) { //   This is the projection in Real Space
		float alpha =alphaDeg*M_PI/180.0;
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

			float Norm0 =  D*sqrt(2*pi);
			float Norm1 =  exp( .5*(P+Q)*(P+Q)) / Norm0 ;
			float Norm2 =  exp( .5*(P-Q)*(P-Q)) / Norm0 ;
			float sD;

			for (int is=(1-Mid);  is<Mid; is++){
				sD = is/D ;
				(*pofalpha)(is+Mid-1) =  Norm1 * exp(-.5*sD*sD)*cos(sD*(P+Q))
                         + Norm2 * exp(-.5*sD*sD)*cos(sD*(P-Q));
			}
			pofalpha-> done_data();
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
		 		(*pofalphak)(2*(iv+Mid-1)) =  exp(-.5*vD*vD)*(cosh(vD*(P+Q)) + cosh(vD*(P-Q)) );
			}
			pofalphak-> done_data();
			pofalphak-> set_complex(false);
			pofalphak-> set_ri(true);
		
			return pofalphak;
   		}
   	}   		
    if (flag==4) {
		cout <<" FH under construction";
   		EMData* OutFT= TwoDTestFunc(Size, p, q, a, b, 1);
   		EMData* TryFH= OutFT -> real2FH(4.0);
   		return TryFH;
   	}   			
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
	u=new float[n-1];

	if (yp1 > .99e30){
		y2[0]=u[0]=0.0;
	}else{
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
	}}

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
	if ((0.0 == t1)&&(0.0 == t2)||(t1 >= t2)) {
		t1 = 0.0f;
		t2 = 90.0f;
	}
	if ((0.0 == p1)&&(0.0 == p2)||(p1 >= p2)) {
		p1 = 0.0f;
		p2 = 359.9f;
	}
	bool skip = ((t1 < 90.0)&&(90.0 == t2)&&(0.0 == p1)&&(p2 > 180.0));
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

/*	static float Util::trilinear_interpolate(float t, float u, float v, float f[]){
	 return Util:trilinear_interpolate(f[0],f[1], f[2],
				  f[3], f[4], float p6, 
							 float p7, float p8, float t,
								  float u, float v)
	}
*/
float Util::triquad(double r, double s, double t, float f[]) {
	const float c2 = 1.0f / 2.0f;
	const float c4 = 1.0f / 4.0f;
	const float c8 = 1.0f / 8.0f;
	float rs = (float)(r*s);
	float st = (float)(s*t);
	float rt = (float)(r*t);
	float rst = (float)(r*st);
	float rsq = (float)(1 - r*r);
	float ssq = (float)(1 - s*s);
	float tsq = (float)(1 - t*t);
	float rm1 = (float)(1 - r);
	float sm1 = (float)(1 - s);
	float tm1 = (float)(1 - t);
	float rp1 = (float)(1 + r);
	float sp1 = (float)(1 + s);
	float tp1 = (float)(1 + t);

	return (float)(
		(-c8) * rst * rm1  * sm1  * tm1 * f[ 0] +
		( c4) * st	* rsq  * sm1  * tm1 * f[ 1] +
		( c8) * rst * rp1  * sm1  * tm1 * f[ 2] +
		( c4) * rt	* rm1  * ssq  * tm1 * f[ 3] +
		(-c2) * t	* rsq  * ssq  * tm1 * f[ 4] +
		(-c4) * rt	* rp1  * ssq  * tm1 * f[ 5] +
		( c8) * rst * rm1  * sp1  * tm1 * f[ 6] +
		(-c4) * st	* rsq  * sp1  * tm1 * f[ 7] +
		(-c8) * rst * rp1  * sp1  * tm1 * f[ 8] +

		( c4) * rs	* rm1  * sm1  * tsq * f[ 9] +
		(-c2) * s	* rsq  * sm1  * tsq * f[10] +
		(-c4) * rs	* rp1  * sm1  * tsq * f[11] +
		(-c2) * r	* rm1  * ssq  * tsq * f[12] +
					  rsq  * ssq  * tsq * f[13] +
		( c2) * r	* rp1  * ssq  * tsq * f[14] +
		(-c4) * rs	* rm1  * sp1  * tsq * f[15] +
		( c2) * s	* rsq  * sp1  * tsq * f[16] +
		( c4) * rs	* rp1  * sp1  * tsq * f[17] +

		( c8) * rst * rm1  * sm1  * tp1 * f[18] +
		(-c4) * st	* rsq  * sm1  * tp1 * f[19] +
		(-c8) * rst * rp1  * sm1  * tp1 * f[20] +
		(-c4) * rt	* rm1  * ssq  * tp1 * f[21] +
		( c2) * t	* rsq  * ssq  * tp1 * f[22] +
		( c4) * rt	* rp1  * ssq  * tp1 * f[23] +
		(-c8) * rst * rm1  * sp1  * tp1 * f[24] +
		( c4) * st	* rsq  * sp1  * tp1 * f[25] +
		( c8) * rst * rp1  * sp1  * tp1 * f[26]);
}

#define  fdata(i,j)      fdata  [ i-1 + (j-1)*nxdata ]
float Util::quadri(float xx, float yy, int nxdata, int nydata, float* fdata)
{
/*
c  purpose: quadratic interpolation 
c 
c  parameters:       xx,yy treated as circularly closed.
c                    fdata - image 1..nxdata, 1..nydata
c
c                    f3    fc       f0, f1, f2, f3 are the values
c                     +             at the grid points.  x is the
c                     + x           point at which the function
c              f2++++f0++++f1       is to be estimated. (it need
c                     +             not be in the first quadrant).
c                     +             fc - the outer corner point
c                    f4             nearest x.
c
c                                   f0 is the value of the fdata at
c                                   fdata(i,j), it is the interior mesh
c                                   point nearest  x.
c                                   the coordinates of f0 are (x0,y0),
c                                   the coordinates of f1 are (xb,y0),
c                                   the coordinates of f2 are (xa,y0),
c                                   the coordinates of f3 are (x0,yb),
c                                   the coordinates of f4 are (x0,ya),
c                                   the coordinates of fc are (xc,yc),
c
c                   o               hxa, hxb are the mesh spacings
c                   +               in the x-direction to the left
c                  hyb              and right of the center point.
c                   +
c            ++hxa++o++hxb++o       hyb, hya are the mesh spacings
c                   +               in the y-direction.
c                  hya
c                   +               hxc equals either  hxb  or  hxa
c                   o               depending on where the corner
c                                   point is located.
c
c                                   construct the interpolant
c                                   f = f0 + c1*(x-x0) +
c                                       c2*(x-x0)*(x-x1) +
c                                       c3*(y-y0) + c4*(y-y0)*(y-y1)
c                                       + c5*(x-x0)*(y-y0)
c
c
*/
    float x, y, dx0, dy0, f0, c1, c2, c3, c4, c5, dxb, dyb;
    float quadri;
    int   i, j, ip1, im1, jp1, jm1, ic, jc, hxc, hyc;
    
    x = xx;
    y = yy;

    // circular closure
    if (x < 1.0)               x = x+(1 - floor(x) / nxdata) * nxdata;
    if (x > (float)nxdata+0.5) x = fmod(x-1.0f,(float)nxdata) + 1.0f;
    if (y < 1.0)               y = y+(1 - floor(y) / nydata) * nydata;
    if (y > (float)nydata+0.5) y = fmod(y-1.0f,(float)nydata) + 1.0f;


    i   = (int) floor(x);
    j   = (int) floor(y);

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
    if (dx0 >= 0) {
       hxc = 1;
    }
    else {
       hxc = -1;
    }
    if (dy0 >= 0) {
       hyc = 1;
    }
    else {
       hyc = -1;
    }
 
    ic  = i + hxc;
    jc  = j + hyc;

    if (ic > nxdata) {
       ic = ic - nxdata;
    }
    else if (ic < 1) {
       ic = ic + nxdata;
    }

    if (jc > nydata) {
       jc = jc - nydata;
    }
    else if (jc < 1) {
       jc = jc + nydata;
    }

    c5  =  ( (fdata(ic,jc) - f0 - hxc * c1 - (hxc * (hxc - 1.0)) * c2 
            - hyc * c3 - (hyc * (hyc - 1.0)) * c4) * (hxc * hyc));

    quadri = f0 + dx0 * (c1 + dxb * c2 + dy0 * c5) + dy0 * (c3 + dyb * c4);

    return quadri; 
}
#undef fdata

float Util::triquad(float R, float S, float T, float* fdata)
{

    float C2 = 1.0 / 2.0;
    float C4 = 1.0 / 4.0;
    float C8 = 1.0 / 8.0;

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
	float res = gsl_sf_bessel_I0(facadj*rt)/val0;
	return res;
}

void Util::KaiserBessel::build_I0table() {
	i0table.resize(ntable+1); // i0table[0:ntable]
	int ltab = int(round(float(ntable)/1.25f));
	fltb = float(ltab)/(K/2);
	float val0 = gsl_sf_bessel_I0(facadj);
	for (int i=ltab+1; i <= ntable; i++) i0table[i] = 0.f;
	for (int i=0; i <= ltab; i++) {
		float s = float(i)/fltb/N;
		if (s < vadjust) {
			float rt = sqrt(1.f - pow(s/vadjust, 2));
			i0table[i] = gsl_sf_bessel_I0(facadj*rt)/val0;
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
	float val0 = sqrt(facadj)*gsl_sf_bessel_I1(facadj);
	for (int i=ltab+1; i <= ntable; i++) i0table[i] = 0.f;
	for (int i=0; i <= ltab; i++) {
		float s = float(i)/fltb/N;
		if (s < vadjust) {
			float rt = sqrt(1.f - pow(s/vadjust, 2));
			i0table[i] = sqrt(facadj*rt)*gsl_sf_bessel_I1(facadj*rt)/val0;
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



EMData* Util::Polar2D(EMData* image, vector<int> numr, string mode){
   int nsam = image->get_xsize();
   int nrow = image->get_ysize();
   int nring = numr.size()/3;
   int lcirc = numr[3*nring-2]+numr[3*nring-1]-1;
   EMData* out = new EMData();
   char cmode = (mode == "F" || mode == "f") ? 'f' : 'h';
   out->set_size(lcirc,1,1);
   alrq(image->get_data(), nsam, nrow, &numr[0], out->get_data(), lcirc, nring, cmode);
   return out;
}

#define  circ(i)         circ   [(i)-1]
#define  numr(i,j)       numr   [((j)-1)*3 + (i)-1]
#define  xim(i,j)        xim    [((j)-1)*nsam + (i)-1]
void Util::alrq(float *xim,  int nsam , int nrow , int *numr,
          float *circ, int lcirc, int nring, char mode)
{
/* 
c                                                                     
c  purpose:                                                          
c                                                                   
c  resmaple to polar coordinates
c                                                                  
*/
   //  dimension         xim(nsam,nrow),circ(lcirc)
   //  integer           numr(3,nring)

   double dfi, dpi;
   int    ns2, nr2, i, inr, l, nsim, kcirc, lt, j;
   float  yq, xold, yold, fi, x, y;

   ns2 = nsam/2+1;
   nr2 = nrow/2+1;
   dpi = 2.0*atan(1.0);

   for (i=1;i<=nring;i++) {
     // radius of the ring
     inr = numr(1,i);
     yq  = inr;
     l   = numr(3,i);
     if (mode == 'h' || mode == 'H') {
        lt = l/2;
     }
     else  {    //if (mode == 'f' || mode == 'F' )
        lt = l/4;
     }

     nsim           = lt-1;
     dfi            = dpi/(nsim+1);
     kcirc          = numr(2,i);
     xold           = 0.0;
     yold           = inr;
     circ(kcirc)    = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);
     xold           = inr;
     yold           = 0.0;
     circ(lt+kcirc) = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);

     if (mode == 'f' || mode == 'F') {
        xold              = 0.0;
        yold              = -inr;
        circ(lt+lt+kcirc) = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);
        xold              = -inr;
        yold              = 0.0;
        circ(lt+lt+lt+kcirc) = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);
     }

     for (j=1;j<=nsim;j++) {
        fi               = dfi*j;
        x                = sin(fi)*yq;
        y                = cos(fi)*yq;
        xold             = x;
        yold             = y;
        circ(j+kcirc)    = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);
        xold             =  y;
        yold             = -x;
        circ(j+lt+kcirc) = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);

        if (mode == 'f' || mode == 'F')  {
           xold                = -x;
           yold                = -y;
           circ(j+lt+lt+kcirc) = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);
           xold                = -y;
           yold                =  x;
           circ(j+lt+lt+lt+kcirc) = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);
        };
     }
   }
 
}




EMData* Util::Polar2Dm(EMData* image, float cns2, float cnr2, vector<int> numr, string mode){
   int nsam = image->get_xsize();
   int nrow = image->get_ysize();
   int nring = numr.size()/3;
   int lcirc = numr[3*nring-2]+numr[3*nring-1]-1;
   EMData* out = new EMData();
   char cmode = (mode == "F" || mode == "f") ? 'f' : 'h';
   out->set_size(lcirc,1,1);
   for (int i=0; i<60000; i++) alrq_ms(image->get_data(), nsam, nrow, cns2, cnr2, &numr[0], out->get_data(), lcirc, nring, cmode);
   return out;
}
void Util::alrq_ms(float *xim, int    nsam, int  nrow, float cns2, float cnr2,
             int  *numr, float *circ, int lcirc, int  nring, char  mode)
{
   double dpi, dfi;
   int    it, jt, inr, l, nsim, kcirc, lt;
   float  yq, xold, yold, fi, x, y;

   //     cns2 and cnr2 are predefined centers
   //     no need to set to zero, all elements are defined

   dpi = 2*atan(1.0);
   for (it=1;it<=nring;it++) {
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
      xold  = 0.0;
      yold  = inr;

      circ(kcirc) = quadri(xold+cns2,yold+cnr2,nsam,nrow,xim);

      xold  = inr;
      yold  = 0.0;
      circ(lt+kcirc) = quadri(xold+cns2,yold+cnr2,nsam,nrow,xim);

      if ( mode == 'f' || mode == 'F' ) {
         xold = 0.0;
         yold = -inr;
         circ(lt+lt+kcirc) = quadri(xold+cns2,yold+cnr2,nsam,nrow,xim);

         xold = -inr;
         yold = 0.0;
         circ(lt+lt+lt+kcirc) = quadri(xold+cns2,yold+cnr2,nsam,nrow,xim);
      }
      
      for (jt=1;jt<=nsim;jt++) {
         fi   = dfi * jt;
         x    = sin(fi) * yq;
         y    = cos(fi) * yq;

         xold = x;
         yold = y;
         circ(jt+kcirc) = quadri(xold+cns2,yold+cnr2,nsam,nrow,xim);

         xold = y;
         yold = -x;
         circ(jt+lt+kcirc) = quadri(xold+cns2,yold+cnr2,nsam,nrow,xim);

         if ( mode == 'f' || mode == 'F' ) {
            xold = -x;
            yold = -y;
            circ(jt+lt+lt+kcirc) = quadri(xold+cns2,yold+cnr2,nsam,nrow,xim);

            xold = -y;
            yold = x;
            circ(jt+lt+lt+lt+kcirc) = quadri(xold+cns2,yold+cnr2,nsam,nrow,xim);
         }
      } // end for jt
   } //end for it
}
#undef  xim

EMData* Util::Polar2Dmi(EMData* image, float cns2, float cnr2, vector<int> numr, string mode, Util::KaiserBessel& kb){
// input image is twice the size of the original image
   int nring = numr.size()/3;
   int lcirc = numr[3*nring-2]+numr[3*nring-1]-1;
   EMData* out = new EMData();
   char cmode = (mode == "F" || mode == "f") ? 'f' : 'h';
   out->set_size(lcirc,1,1);
   Util::alrq_msi(image, cns2, cnr2, &numr[0], out->get_data(), lcirc, nring, cmode, kb);
   return out;
}

void Util::alrq_msi(EMData* image, float cns2, float cnr2,
             int  *numr, float *circ, int lcirc, int  nring, char  mode, Util::KaiserBessel& kb)
{
   double dpi, dfi;
   int    it, jt, inr, l, nsim, kcirc, lt;
   float  yq, xold, yold, fi, x, y;

   //     cns2 and cnr2 are predefined centers
   //     no need to set to zero, all elements are defined

   dpi = 2*atan(1.0);
   for (it=1;it<=nring;it++) {
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
      xold  = 0.0;
      yold  = inr;
      circ(kcirc) = image->get_pixel_conv(xold+cns2-1.0f,yold+cnr2-1.0f,0,kb);
      
      xold  = inr;
      yold  = 0.0;
      circ(lt+kcirc) = image->get_pixel_conv(xold+cns2-1.0f,yold+cnr2-1.0f,0,kb);

      if ( mode == 'f' || mode == 'F' ) {
         xold = 0.0;
         yold = -inr;
         circ(lt+lt+kcirc) = image->get_pixel_conv(xold+cns2-1.0f,yold+cnr2-1.0f,0,kb);

         xold = -inr;
         yold = 0.0;
         circ(lt+lt+lt+kcirc) = image->get_pixel_conv(xold+cns2-1.0f,yold+cnr2-1.0f,0,kb);
      }
      
      for (jt=1;jt<=nsim;jt++) {
         fi   = dfi * jt;
         x    = sin(fi) * yq;
         y    = cos(fi) * yq;

         xold = x;
         yold = y;
         circ(jt+kcirc) = image->get_pixel_conv(xold+cns2-1.0f,yold+cnr2-1.0f,0,kb);

         xold = y;
         yold = -x;
         circ(jt+lt+kcirc) = image->get_pixel_conv(xold+cns2-1.0f,yold+cnr2-1.0f,0,kb);

         if ( mode == 'f' || mode == 'F' ) {
            xold = -x;
            yold = -y;
            circ(jt+lt+lt+kcirc) = image->get_pixel_conv(xold+cns2-1.0f,yold+cnr2-1.0f,0,kb);

            xold = -y;
            yold = x;
            circ(jt+lt+lt+lt+kcirc) = image->get_pixel_conv(xold+cns2-1.0f,yold+cnr2-1.0f,0,kb);
         }
      } // end for jt
   } //end for it
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
#define  xcmplx(i,j)  xcmplx [((j)-1)*2 + (i)-1]
#define  br(i)        br     [(i)-1]
#define  bi(i)        bi     [(i)-1]
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

   n=(int)pow(2.f,ln);

   k=abs(ks);
   l=16-ln;
   b3=n*k;
   b6=b3;
   b7=k;
   if (ks > 0) {
      sgn=1.0;
   }
   else {
      sgn=-1.0;
      rni=1.0/(float)(n);
      j=1;
      for (i=1;i<=n;i++) {
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
   cc=2.0*pow(tab1(l),2);
   c=1.0-cc;
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

   n=(int)pow(2.f,ln);

   k=abs(ks);
   l=16-ln;
   b3=n*k;
   b6=b3;
   b7=k;
   if( ks > 0 ) {
      sgn=1.0;
   } 
   else {
      sgn=-1.0;
      rni=1.0/(float)n;
      j=1;
      for (i=1; i<=n;i++) {
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
   if (b5 <= b3)  goto  L14;
   if (b6 == b7)  goto  L20;

   b4=b7;
   cc=2.0*pow(tab1(l),2);
   c=1.0-cc;
   l=l+1;
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
   n=(int)pow(2.f,nu1);
   isub=16-nu1;

   ss=-tab1(isub);
   cc=-2.0*pow(tab1(isub-1),2.f);
   c=1.0;
   s=0.0;
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
         xcmplx(1,i1)=0.5*((tr1+tr2)+(ti1+ti2)*c-(tr1-tr2)*s);
         xcmplx(1,i2)=0.5*((tr1+tr2)-(ti1+ti2)*c+(tr1-tr2)*s);
         xcmplx(2,i1)=0.5*((ti1-ti2)-(ti1+ti2)*s-(tr1-tr2)*c);
         xcmplx(2,i2)=0.5*(-(ti1-ti2)-(ti1+ti2)*s-(tr1-tr2)*c);
     }
   }
   else {
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
/*
   tab1(1)=9.58737990959775e-5;
   tab1(2)=1.91747597310703e-4;
   tab1(3)=3.83495187571395e-4;
   tab1(4)=7.66990318742704e-4;
   tab1(5)=1.53398018628476e-3;
   tab1(6)=3.06795676296598e-3;
   tab1(7)=6.13588464915449e-3;
   tab1(8)=1.22715382857199e-2;
   tab1(9)=2.45412285229123e-2;
   tab1(10)=4.90676743274181e-2;
   tab1(11)=9.80171403295604e-2;
   tab1(12)=1.95090322016128e-1;
   tab1(13)=3.82683432365090e-1;
   tab1(14)=7.07106781186546e-1;
   tab1(15)=1.00000000000000;
*/
   nu=abs(nv);
   inv=nv/nu;
   nu1=nu-1;
   n=(int)pow(2.f,nu1);
   isub=16-nu1;
   ss=-tab1(isub);
   cc=-2.0*pow(tab1(isub-1),2);
   c=1.0;
   s=0.0;
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
   }
   else {
      tr=xcmplx(1,1);
      ti=xcmplx(2,1);
      xcmplx(1,1)=0.5*(tr+ti);
      xcmplx(2,1)=0.5*(tr-ti);
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

void Util::Frngs(EMData* circ, vector<int> numr){
   int nring = numr.size()/3;
   for (int i=0; i<60000; i++) frngs(circ->get_data(), &numr[0],  nring);
}
void Util::frngs(float *circ, int *numr, int nring){
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
#undef  circ
//---------------------------------------------------
#define  b(i)            b      [(i)-1]
void Util::prb1d(double *b, int npoint, float *pos)
{
   double  c2,c3;
   int     nhalf;

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
   else if (npoint == 9) {
      c2 = (1708.*b(1) + 581.*b(2) - 246.*b(3) - 773.*b(4)
         - 1000.*b(5) - 927.*b(6) - 554.*b(7) + 119.*b(8)
         + 1092.*b(9) ) / (-4620.);
      c3 = (28.*b(1) + 7.*b(2) - 8.*b(3) - 17.*b(4) - 20.*b(5)
         - 17.*b(6) - 8.*b(7) + 7.*b(8) + 28.*b(9) ) / 924.0;
   }
   if (c3 != 0.0)  *pos = c2/(2.0*c3) - nhalf;
}
#undef  b
boost::tuple<double, float, int> Util::Crosrng_e(EMData*  circ1, EMData* circ2, vector<int> numr, int neg) {
   int nring = numr.size()/3;
   int lcirc = numr[3*nring-2]+numr[3*nring-1]-1;
   int maxrin = numr[numr.size()-1];
   double qn;   float  tot;
   crosrng_e(circ1->get_data(), circ2->get_data(), lcirc, nring, maxrin, &numr[0], 
		      &qn, &tot, neg);
   return boost::make_tuple(qn, tot, neg);
}
#define  circ1(i)        circ1  [(i)-1]
#define  circ2(i)        circ2  [(i)-1]
#define  t(i)            t      [(i)-1]
#define  q(i)            q      [(i)-1]
#define  b(i)            b      [(i)-1]
#define  t7(i)           t7     [(i)-1]


//-----------------------------------------------
void Util::crosrng_e(float *circ1, float *circ2, int lcirc,
                     int    nring, int   maxrin, int *numr,
                     double *qn, float *tot, int neg)
{
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
   int    i, j, k, ip, jc, numr3i, numr2i, jtot;
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

      t(1) = (circ1(numr2i)) * circ2(numr2i);

      if (numr3i != maxrin) {
         // test .ne. first for speed on some compilers
	 t(numr3i+1) = circ1(numr2i+1) * circ2(numr2i+1);
	 t(2)        = 0.0;

         if (neg) {
            // first set is conjugated (mirrored)
	    for (j=3;j<=numr3i;j=j+2) {
	       jc = j+numr2i-1;
	       t(j) =(circ1(jc))*circ2(jc)-(circ1(jc+1))*circ2(jc+1);
	       t(j+1) = -(circ1(jc))*circ2(jc+1)-(circ1(jc+1))*circ2(jc);
	    } 
         } 
         else {
	    for (j=3;j<=numr3i;j=j+2) {
	       jc = j+numr2i-1;
	       t(j) = (circ1(jc))*circ2(jc) + (circ1(jc+1))*circ2(jc+1);
	       t(j+1) = -(circ1(jc))*circ2(jc+1) + (circ1(jc+1))*circ2(jc);
	    }
         } 
         for (j=1;j<=numr3i+1;j++) q(j) = q(j) + t(j);
      }
      else {
	 t(2) = circ1(numr2i+1) * circ2(numr2i+1);
         if (neg) {
            // first set is conjugated (mirrored)
	    for (j=3;j<=maxrin;j=j+2) {
	       jc = j+numr2i-1;
	       t(j) = (circ1(jc))*circ2(jc) - (circ1(jc+1))*circ2(jc+1);
	       t(j+1) = -(circ1(jc))*circ2(jc+1) - (circ1(jc+1))*circ2(jc);
	    }
         }
         else {
	    for (j=3;j<=maxrin;j=j+2) {
	       jc = j+numr2i-1;
	       t(j) = (circ1(jc))*circ2(jc) + (circ1(jc+1))*circ2(jc+1);
	       t(j+1) = -(circ1(jc))*circ2(jc+1) + (circ1(jc+1))*circ2(jc);
	    } 
         }
         for (j = 1; j <= maxrin; j++) q(j) = q(j) + t(j);
      }
   }

   fftr_d(q,ip);

   *qn = -1.0e20;
   for (j=1;j<=maxrin;j++) {
      if (q(j) >= *qn) {
         *qn = q(j);
	 jtot = j;
      }
   } 

   for (k=-3;k<=3;k++) {
      j = (jtot+k+maxrin-1)%maxrin + 1;
      t7(k+4) = q(j);
   }

   prb1d(t7,7,&pos);

   *tot = (float)jtot + pos;

   if (q) free(q);
   if (t) free(t);
}

Dict Util::Crosrng_ms(EMData* circ1, EMData* circ2, vector<int> numr) {
   int nring = numr.size()/3;
   int lcirc = numr[3*nring-2]+numr[3*nring-1]-1;
   int maxrin = numr[numr.size()-1];
   double qn; float tot; double qm; float tmt;
   for (int i=0; i<60000; i++) crosrng_ms(circ1->get_data(), circ2->get_data(), lcirc, nring, maxrin, 
              &numr[0], &qn, &tot, &qm, &tmt);
   Dict retvals;
   retvals["qn"] = qn;
   retvals["tot"] = tot;
   retvals["qm"] = qm;
   retvals["tmt"] = tmt;
   return retvals;
}

//---------------------------------------------------
void Util::crosrng_ms(float *circ1, float *circ2, int  lcirc, int  nring,
                      int   maxrin, int   *numr , double *qn, float *tot,
                      double   *qm, float *tmt)
{
/*
c
c  checks both straight & mirrored positions
c
c  input - fourier transforms of rings!!
c  circ1 already multiplied by weights!
c
c  notes: aug 04 attempted speedup using 
c       premultiply  arrays ie( circ12 = circ1 * circ2) much slower
c       various  other attempts  failed to yield improvement
c       this is a very important compute demand in alignmen & refine.
c       optional limit on angular search should be added.
*/

   // dimension         circ1(lcirc),circ2(lcirc)

   // t(maxrin), q(maxrin), t7(-3:3)  //maxrin+2 removed
   double *t, *q, t7[7];

   int   ip, jc, numr3i, numr2i, i, j, k, jtot;
   float t1, t2, t3, t4, c1, c2, d1, d2, pos;

   *qn  = 0.0;
   *qm  = 0.0;
   *tot = 0.0;
   *tmt = 0.0; 

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

   for (i=1;i<=nring;i++) {

      numr3i = numr(3,i);
      numr2i = numr(2,i);

      t1   = circ1(numr2i) * circ2(numr2i);
      q(1) = q(1)+t1;
      t(1) = t(1)+t1;

      if (numr3i == maxrin)  {
         t1   = circ1(numr2i+1) * circ2(numr2i+1);
         q(2) = q(2)+t1;
         t(2) = t(2)+t1;
      }
      else {
	 t1          = circ1(numr2i+1) * circ2(numr2i+1);
	 q(numr3i+1) = q(numr3i+1)+t1;
      }

      for (j=3;j<=numr3i;j=j+2) {
	 jc     = j+numr2i-1;

 	 c1     = circ1(jc);
 	 c2     = circ1(jc+1);
         d1     = circ2(jc);
         d2     = circ2(jc+1);

  	 t1     = c1 * d1;
 	 t3     = c1 * d2;
 	 t2     = c2 * d2;
 	 t4     = c2 * d1;

	 q(j)   = q(j)   + t1 + t2;
	 q(j+1) = q(j+1) - t3 + t4;
	 t(j)   = t(j)   + t1 - t2;
	 t(j+1) = t(j+1) - t3 - t4;
      } 
  }

  fftr_d(q,ip);

  jtot = 0;
  *qn  = -1.0e20;
  for (j=1; j<=maxrin; j++) {
     if (q(j) >= *qn) {
        *qn  = q(j);
        jtot = j;
     }
  }

 
  for (k=-3;k<=3;k++) {
    j = ((jtot+k+maxrin-1)%maxrin)+1;
    t7(k+4) = q(j);
  }

  // interpolate
  prb1d(t7,7,&pos);
  *tot = (float)(jtot)+pos;

  // Do not interpolate
  //*tot = (float)(jtot);

  // mirrored
  fftr_d(t,ip);

  // find angle
  *qm = -1.0e20;
  for (j=1; j<=maxrin;j++) {
     if ( t(j) >= *qm ) {
        *qm   = t(j);
        jtot = j;
     }
  }


  // find angle
  for (k=-3;k<=3;k++) {
    j       = ((jtot+k+maxrin-1)%maxrin) + 1;
    t7(k+4) = t(j);
  }

  // interpolate

  prb1d(t7,7,&pos);
  *tmt = float(jtot) + pos;

  // Do not interpolate
  //*tmt = float(jtot);
  
  free(t);
  free(q);
}
//  Try rotational gridding

Dict Util::Crosrng_msr(EMData* circ1, EMData* circ2, vector<int> numr) {
   int nring = numr.size()/3;
   int lcirc = numr[3*nring-2]+numr[3*nring-1]-1;
   int maxrin = numr[numr.size()-1];
   float qn; float tot; float qm; float tmt;
   crosrng_msr(circ1->get_data(), circ2->get_data(), lcirc, nring, maxrin, 
              &numr[0], &qn, &tot, &qm, &tmt);
   Dict retvals;
   retvals["qn"] = qn;
   retvals["tot"] = tot;
   retvals["qm"] = qm;
   retvals["tmt"] = tmt;
   return retvals;
}
#define  temp(i)            temp      [(i)-1]

//---------------------------------------------------
void Util::crosrng_msr(float *circ1, float *circ2, int  lcirc, int  nring,
                      int   maxrin, int   *numr , float *qn, float *tot,
                      float   *qm, float *tmt)
{
/*
c
c  checks both straight & mirrored positions
c
c  input - fourier transforms of rings!!
c  circ1 already multiplied by weights!
c
c  notes: aug 04 attempted speedup using 
c       premultiply  arrays ie( circ12 = circ1 * circ2) much slower
c       various  other attempts  failed to yield improvement
c       this is a very important compute demand in alignmen & refine.
c       optional limit on angular search should be added.
*/

   // dimension         circ1(lcirc),circ2(lcirc)

   // t(maxrin), q(maxrin), t7(-3:3)  //maxrin+2 removed
   double *t, *q, t7[7];
   float *temp;

   int   ip, jc, numr3i, numr2i, i, j, k, jtot;
   float t1, t2, t3, t4, c1, c2, d1, d2, pos;

   *qn  = 0.0;
   *qm  = 0.0;
   *tot = 0.0;
   *tmt = 0.0; 

#ifdef WIN32
	ip = -(int)(log((float)maxrin)/log(2.0f));
#else
   ip = -(int)(log2(maxrin));
#endif	//WIN32

   //  c - straight  = circ1 * conjg(circ2)
   //  zero q array
  
   q = (double*)calloc(maxrin,sizeof(double));  

   //   t - mirrored  = conjg(circ1) * conjg(circ2)
   //   zero t array
   t = (double*)calloc(maxrin,sizeof(double));


   temp = (float*)calloc(maxrin,sizeof(float));

   //   premultiply  arrays ie( circ12 = circ1 * circ2) much slower

   for (i=1;i<=nring;i++) {

      numr3i = numr(3,i);
      numr2i = numr(2,i);

      t1   = circ1(numr2i) * circ2(numr2i);
      q(1) = q(1)+t1;
      t(1) = t(1)+t1;

      if (numr3i == maxrin)  {
         t1   = circ1(numr2i+1) * circ2(numr2i+1);
         q(2) = q(2)+t1;
         t(2) = t(2)+t1;
      }
      else {
	 t1          = circ1(numr2i+1) * circ2(numr2i+1);
	 q(numr3i+1) = q(numr3i+1)+t1;
      }

      for (j=3;j<=numr3i;j=j+2) {
	 jc     = j+numr2i-1;

 	 c1     = circ1(jc);
 	 c2     = circ1(jc+1);
         d1     = circ2(jc);
         d2     = circ2(jc+1);

  	 t1     = c1 * d1;
 	 t3     = c1 * d2;
 	 t2     = c2 * d2;
 	 t4     = c2 * d1;

	 q(j)   = q(j)   + t1 + t2;
	 q(j+1) = q(j+1) - t3 + t4;
	 t(j)   = t(j)   + t1 - t2;
	 t(j+1) = t(j+1) - t3 - t4;
      } 
  }

  // straight
  for (i=1; i<=maxrin; i++) {temp(i)=q(i);}
  fftr_q(temp,ip);

  jtot = 0;
  *qn  = -1.0e20;
  for (j=1; j<=maxrin; j++) {
     if (temp(j) >= *qn) {
        *qn  = temp(j);
        jtot = j;
     }
  }
  
 
  for (k=-3;k<=3;k++) {
    j = ((jtot+k+maxrin-1)%maxrin)+1;
    t7(k+4) = temp(j);
  }

  // interpolate
  prb1d(t7,7,&pos);
  *tot = (float)(jtot)+pos;

  // mirrored
  for (i=1; i<=maxrin; i++) {temp(i)=t(i);}
  fftr_q(temp,ip);

  // find angle
  *qm = -1.0e20;
  for (j=1; j<=maxrin;j++) {
     if ( temp(j) >= *qm ) {
        *qm   = temp(j);
        jtot = j;
     }
  }

  // find angle
  for (k=-3;k<=3;k++) {
    j       = ((jtot+k+maxrin-1)%maxrin) + 1;
    t7(k+4) = t(j);
  }

  // interpolate

  prb1d(t7,7,&pos);
  *tmt = float(jtot) + pos;

  free(t);
  free(q);
  free(temp);
}
#undef temp


#define  dout(i,j)   dout[i+maxrin*j]
EMData* Util::Crosrng_msg(EMData* circ1, EMData* circ2, vector<int> numr) {
   int nring = numr.size()/3;
   int lcirc = numr[3*nring-2]+numr[3*nring-1]-1;
   int maxrin = numr[numr.size()-1];

   // t(maxrin), q(maxrin)  // removed +2
   double *t, *q;

   //  q - straight  = circ1 * conjg(circ2)
   //  zero q array
   q = (double*)calloc(maxrin,sizeof(double));

   //   t - mirrored  = conjg(circ1) * conjg(circ2)
   //   zero t array
   t = (double*)calloc(maxrin,sizeof(double));

   crosrng_msg(circ1->get_data(), circ2->get_data(), &q[0], &t[0], lcirc, nring, maxrin, &numr[0]);
   EMData* out = new EMData();
   out->set_size(maxrin,2,1);
   float *dout = out->get_data();
   for (int i=0; i<maxrin; i++) {dout(i,0)=q[i]; dout(i,1)=t[i];}
   /*out->set_size(maxrin,1,1);
   float *dout = out->get_data();
   for (int i=0; i<maxrin; i++) {dout(i,0)=q[i];}*/
   free(t);
   free(q);
   return out;
}
#undef out
//---------------------------------------------------
void Util::crosrng_msg(float *circ1, float *circ2, double *q, double *t, int  lcirc, int  nring,
                      int  maxrin, int   *numr )
{
/*
c
c  checks both straight & mirrored positions
c
c  input - fourier transforms of rings!!
c  circ1 already multiplied by weights!
c
c  notes: aug 04 attempted speedup using 
c       premultiply  arrays ie( circ12 = circ1 * circ2) much slower
c       various  other attempts  failed to yield improvement
c       this is a very important compute demand in alignmen & refine.
c       optional limit on angular search should be added.
*/

   // dimension         circ1(lcirc),circ2(lcirc)

   int   ip, jc, numr3i, numr2i, i, j;
   float t1, t2, t3, t4, c1, c2, d1, d2;

#ifdef _WIN32
	ip = -(int)(log((float)maxrin)/log(2.0f));
#else
	ip = -(int)(log2(maxrin));
#endif	//_WIN32

   //  q - straight  = circ1 * conjg(circ2)

   //   t - mirrored  = conjg(circ1) * conjg(circ2)

   //   premultiply  arrays ie( circ12 = circ1 * circ2) much slower

   for (i=1;i<=nring;i++) {

      numr3i = numr(3,i);
      numr2i = numr(2,i);

      t1   = circ1(numr2i) * circ2(numr2i);
      q(1) = q(1)+t1;
      t(1) = t(1)+t1;

      if (numr3i == maxrin)  {
         t1   = circ1(numr2i+1) * circ2(numr2i+1);
         q(2) = q(2)+t1;
         t(2) = t(2)+t1;
      }
      else {
	 t1          = circ1(numr2i+1) * circ2(numr2i+1);
	 q(numr3i+1) = q(numr3i+1)+t1;
      }

      for (j=3;j<=numr3i;j=j+2) {
	 jc     = j+numr2i-1;

 	 c1     = circ1(jc);
 	 c2     = circ1(jc+1);
         d1     = circ2(jc);
         d2     = circ2(jc+1);

  	 t1     = c1 * d1;
 	 t3     = c1 * d2;
 	 t2     = c2 * d2;
 	 t4     = c2 * d1;

	 q(j)   = q(j)   + t1 + t2;
	 q(j+1) = q(j+1) - t3 + t4;
	 t(j)   = t(j)   + t1 - t2;
	 t(j+1) = t(j+1) - t3 - t4;
      } 
  }
  
  // straight
  fftr_d(q,ip);

  // mirrored
  fftr_d(t,ip);
}
#undef  circ1
#undef  circ2
#undef  t
#undef  q
#undef  b
#undef  t7


#undef  numr



#define old_ptr(i,j,k) old_ptr[(i+(j+(k*ny))*nx)]
#define new_ptr(iptr,jptr,kptr) new_ptr[iptr+(jptr+(kptr*new_ny))*new_nx]
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
		LOGERR("The Parameters for decimation cannot exceed the center of the image.");
		throw ImageDimensionException("The Parameters for decimation cannot exceed the center of the image.");	 
	}
	/* ============================== */
	
	
	/*    Calculation of the start point */
	int new_st_x=(nx/2)%x_step,new_st_y=(ny/2)%y_step,new_st_z=(nz/2)%z_step;
	/* ============================*/
	
	
	/* Calculation of the size of the decimated image */
	int rx=2*(nx/(2*x_step)),ry=2*(ny/(2*y_step)),rz=2*(nz/(2*z_step));
	int r1=int(ceil((nx-(x_step*rx))/(1.f*x_step))),r2=int(ceil((ny-(y_step*ry))/(1.f*y_step)));
	int r3=int(ceil((nz-(z_step*rz))/(1.f*z_step)));
	if(r1>1){r1=1;}
	if(r2>1){r2=1;}
	if(r3>1){r3=1;}
	int new_nx=rx+r1,new_ny=ry+r2,new_nz=rz+r3;
	/* ===========================================*/
	
	
	EMData* img2 = new EMData();
	img2->set_size(new_nx,new_ny,new_nz);
	float *new_ptr=img2->get_data();
	float *old_ptr=img->get_data();
	int iptr,jptr,kptr=0;
	for (int k=new_st_z;k<nz;k+=z_step){jptr=0;
		for (int j=new_st_y;j<ny;j+=y_step){iptr=0;
			for (int i=new_st_x;i<nx;i+=x_step){				
				new_ptr(iptr,jptr,kptr)=old_ptr(i,j,k);
			iptr++;}
		jptr++;}
	kptr++;}
	img2->update();
	return img2;
}
#undef old_ptr
#undef new_ptr

#define inp(i,j,k) inp[(i+new_st_x)+((j+new_st_y)+((k+new_st_z)*ny))*nx]
#define outp(i,j,k) outp[i+(j+(k*new_ny))*new_nx]
EMData* Util::window(EMData* img,int new_nx,int new_ny, int new_nz, int x_offset, int y_offset, int z_offset)
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
	if(new_nx>nx || new_ny>ny || new_nz>nz)
		throw ImageDimensionException("The size of the windowed image cannot exceed the input image size.");
	if((nx/2)-(new_nx/2)+x_offset<0 || (ny/2)-(new_ny/2)+y_offset<0 || (nz/2)-(new_nz/2)+z_offset<0)
		throw ImageDimensionException("The offset imconsistent with the input image size. Solution: Change the offset parameters");
	if(x_offset>((nx-(nx/2))-(new_nx-(new_nx/2))) || y_offset>((ny-(ny/2))-(new_ny-(new_ny/2))) || z_offset>((nz-(nz/2))-(new_nz-(new_nz/2))))
		throw ImageDimensionException("The offset imconsistent with the input image size. Solution: Change the offset parameters");
	/* ============================== */
	
	EMData* wind= new EMData();
	wind->set_size(new_nx,new_ny,new_nz);
	float *outp=wind->get_data();
	float *inp=img->get_data();

	
	/*    Calculation of the start point */
	int new_st_x=int((nx/2-new_nx/2) + x_offset),
	    new_st_y=int((ny/2-new_ny/2) + y_offset),  
	    new_st_z=int((nz/2-new_nz/2) + z_offset);
	/* ============================== */
	    
	/* Exception Handle */
	if (new_st_x<0 || new_st_y<0 || new_st_z<0)   //  WHAT HAPPENS WITH THE END POINT CHECK??  PAP
		throw ImageDimensionException("The offset inconsistent with the input image size. Solution: Change the offset parameters");
	/* ============================== */
	
	
	for (int k=0;k<new_nz;k++)
	    for(int j=0;j<new_ny;j++)
	        for(int i=0;i<new_nx;i++)
		     outp(i,j,k)=inp(i,j,k);		    
	wind->update();
	return wind;
}
#undef inp
#undef outp

#define inp(i,j,k) inp[i+(j+(k*ny))*nx]
#define outp(i,j,k) outp[(i+new_st_x)+((j+new_st_y)+((k+new_st_z)*new_ny))*new_nx]
EMData *Util::pad(EMData* img,Dict params,int new_nx, int new_ny, int new_nz, int x_offset, int y_offset, int z_offset)
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
	if(new_nx<nx || new_ny<ny || new_nz<nz)
		throw ImageDimensionException("The size of the padding image cannot be below the input image size.");
	if((new_nx/2)-(nx/2)+x_offset<0 || (new_ny/2)-(ny/2)+y_offset<0 || (new_nz/2)-(nz/2)+z_offset<0)
		throw ImageDimensionException("The offset imconsistent with the input image size. Solution: Change the offset parameters");
	if(x_offset>((new_nx-(new_nx/2))-(nx-(nx/2))) || y_offset>((new_ny-(new_ny/2))-(ny-(ny/2))) || z_offset>((new_nz-(new_nz/2))-(nz-(nz/2))))
		throw ImageDimensionException("The offset imconsistent with the input image size. Solution: Change the offset parameters");
	/* ============================== */

	
	EMData* pading=new EMData();
	pading->set_size(new_nx,new_ny,new_nz);
	float *inp=img->get_data();
	float *outp=pading->get_data();
		
	
	/* Calculation of the average and the circumference values for background substitution 
	=======================================================================================*/
	float background;
	if (strcmp("average",params["average"])==0){
		background = img->get_attr("mean");
		}
	else if (strcmp("circumference",params["circumference"])==0)
	{
		float sum1=0.f;
		int cnt=0;
		for(int i=0;i<nx;i++){
			sum1 += inp(i,0,0) + inp(i,ny-1,nz-1);
			cnt+=2;}
		if(nz-1 == 0)
		{
			for (int j=1;j<ny-1;j++){
				sum1 += inp(1,j,0) + inp(nx-1,j,0);
				cnt+=2;}
		}
		else
		{
		for (int k=1;k<nz-1;k++){
			for (int j=1;j<ny-1;j++){
				sum1 += inp(1,j,0) + inp(nx-1,j,0);
				cnt+=2;}
		}
		}
		background = sum1/cnt;
	}		
	else{
		background = static_cast<int>(params["background"]);
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
	        for(int i=0;i<nx;i++){
			outp(i,j,k)=inp(i,j,k); 
			}
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
/*
 Performs inplace integer cyclic shift as specified by the "dx","dy","dz" parameters on a 3d volume.
 Implements the inplace swapping using reversals as descibed in  also:
    http://www.csse.monash.edu.au/~lloyd/tildeAlgDS/Intro/Eg01/
    
 
* @author  Phani Ivatury
* @date 18-2006
* @see http://www.csse.monash.edu.au/~lloyd/tildeAlgDS/Intro/Eg01/
*
*
* A[0] A[1] A[2] A[3] A[4] A[5] A[6] A[7] A[8] A[9]
*
* 10   20   30   40   50   60   70   80   90   100
* ------------
*    m = 3 (shift left three places)
*
* Reverse the items from 0..m-1 and m..N-1:
* 
* 30   20   10   100  90   80   70   60   50   40
*
* Now reverse the entire sequence:
*
* 40   50   60   70   80   90   100  10   20   30

    
    cycl_shift() in libpy/fundementals.py calls this function
    
    Usage: 
    EMData *im1 = new EMData();
    im1->set_size(70,80,85);
    im1->to_one();
    Dict params; params["dx"] = 10;params["dy"] = 10000;params["zx"] = -10;
    Utils::cyclicshift(im1,params);
    im1.peak_search(1,1)
*/

if (image->is_complex())
                throw ImageFormatException("Real image required for "
                                                   "IntegerCyclicShift2DProcessor");
 
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
	 
	 
 #ifdef DEBUG
         std::cout << dx << std::endl;
         std::cout << dy << std::endl;
	     std::cout << dz << std::endl;
 #endif
         int mx = -(dx - nx);
         int my = -(dy - ny);
	 	 int mz = -(dz - nz);
	 
         float* data = image->get_data();
         // x-reverses
         if (mx != 0) {
	         for (int iz = 0; iz < nz; iz++)
		 	for (int iy = 0; iy < ny; iy++) {
	                         // reverses for column iy
        	                 int offset = nx*iy + nx*ny*iz; // starting location for column iy in slice iz
                	         reverse(&data[offset],&data[offset+mx]);
                        	 reverse(&data[offset+mx],&data[offset+nx]);
                         	 reverse(&data[offset],&data[offset+nx]);
		         }
         }
         // y-reverses
         if (my != 0) {	 
	         for (int iz = 0; iz < nz; iz++) {
		 	int offset = nx*ny*iz;
            	     	colreverse(&data[offset], &data[offset + my*nx], nx);
                     	colreverse(&data[offset + my*nx], &data[offset + ny*nx], nx);
                     	colreverse(&data[offset], &data[offset + ny*nx], nx);
		}
         }
	 if (mz != 0) {
                 slicereverse(&data[0], &data[mz*ny*nx], nx, ny);
                 slicereverse(&data[my*ny*nx], &data[nz*ny*nx], nx, ny);
                 slicereverse(&data[0], &data[nz*ny*nx], nx ,ny);
         }
	image->done_data();	 
}

//-----------------------------------------------------------------------------------------------------------------------


/*void Util::histogram(EMData* image, EMData* mask)
{
	if (image->is_complex())
                throw ImageFormatException("Cannot do Histogram on Fourier Image");
	float hmax = image->get_attr("maximum");
	float hmin = image->get_attr("minimum");
	float *imageptr,*maskptr;
	int nx=image->get_xsize();
	int ny=image->get_ysize();
	int nz=image->get_zsize();
	
	if(mask != NULL){
		if(nx != mask->get_xsize() || ny != mask->get_ysize() || nz != mask->get_zsize())
			throw ImageDimensionException("The size of mask image should be of same size as the input image");
	 }
	int nbins = 128;
	float *freq = new float[nbins];
		
	for(int i=0;i<nbins;i++)
		freq[i]=0;
	imageptr=image->get_data();
	maskptr =mask->get_data();
	if(mask!=NULL)
	{
		for (int i = 0;i < nx*ny*nz; i++)
	    	{
	      		if (maskptr[i]>=0.5f)
	      		{			              
			       hmax = (hmax < imageptr[i])?imageptr[i]:hmax;
			       hmin = (hmin > imageptr[i])?imageptr[i]:hmin;
			}
		}
	}
	float hdiff = hmax - hmin;
	float ff = (nbins-1)/hdiff;
	float fnumel=0.f,hav=0.f,hav2=0.f;
	
	if(mask!=NULL)
	{
		for(int i = 0;i < nx*ny*nz;i++)
		{
			int jbin = static_cast<int>((imageptr[i]-hmin)*ff + 1.5);
			if(jbin >= 1 && jbin <= nbins)
			{
				freq[jbin-1] += 1.0;
				fnumel += 1;
				hav += imageptr[i];
				hav2 += (double)pow(imageptr[i],2);
			}
		}
	}
	else
	{
		for(int i = 0;i < nx*ny*nz;i++)
		{
			float bin_mode;
			float hist_max = freq[1];
			int max_bin = 0;
			for(int j=1;j<nbins;j++)
			{
				if(freq[j] >= hist_max)
				{
					hist_max = freq[j];
					max_bin = j;
				}
			}
			if(max_bin == 0)
				bin_mode = 0.5;
			else if(max_bin == (nbins-1))
				bin_mode = static_cast<float>(nbins) - 0.5;
			else
			{
				float YM1 = freq[max_bin - 1];
				float YP1 = freq[max_bin + 1];
				bin_mode = static_cast<float>(max_bin-1) + ((YM1 - YP1)*0.5/(YM1 + YP1 - (2.0*hist_max)));
			}
			//float hist_mode = hmin + (bin_mode*bin_size);
			
			double dtop = hav2 - ((hav*hav)/(fnumel=(fnumel==0.f)?1:fnumel));
			
			if(dtop < 0.0)
				throw ImageFormatException("Cannot be negative");
			
			hav = hav/(fnumel=(fnumel==0.f)?1:fnumel);
			//float hsig = sqrt(dtop/(fnumel-1));
			
			if(maskptr[i] >= 0.5)
			{
				int jbin = static_cast<int>((imageptr[i]-hmin)*ff + 1.5);
				if(jbin >= 1 && jbin <= nbins)
				{
					freq[jbin-1] += 1.0;
					fnumel += 1;
					hav += imageptr[i];
					hav2 += (double)pow(imageptr[i],2);
				}
			}
		}
	}
	delete[] freq;
}
*/
			
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
	int size_ref = ((ref->get_xsize())*(ref->get_ysize())*(ref->get_zsize()));
	int size_img = ((img->get_xsize())*(img->get_ysize())*(img->get_zsize()));
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
	
	vector<float> img_data = Util::infomask(img,mask);
	float img_avg = img_data[0];
	float img_sig = img_data[1];
	
	/* The image under mask -- size calculation */
	int cnt=0;
	for(int i=0;i<size_img;i++)
		if (mask_ptr[i]>0.5f)
				cnt++;
	/* ===================================================== */
	
	/* Histogram of reference image calculation */
	float ref_h_diff = ref_h_max - ref_h_min;
		
	#ifdef _WIN32
		int hist_len = _MIN((int)size_ref/16,_MIN((int)size_img/16,256));
	#else
		int hist_len = std::min((int)size_ref/16,std::min((int)size_img/16,256));
	#endif	//_WIN32
	
	float *ref_freq_bin = new float[3*hist_len];

	//initialize value in each bin to zero 
	for (int i = 0;i < (3*hist_len);i++)
		ref_freq_bin[i] = 0.f;
		
	for (int i = 0;i < size_ref;i++)
	{
		int L = static_cast<int>(((ref_ptr[i] - ref_h_min)/ref_h_diff) * (hist_len-1) + hist_len+1);
		ref_freq_bin[L]++;
	}
	for (int i = 0;i < (3*hist_len);i++)
		ref_freq_bin[i] *= static_cast<float>(cnt)/static_cast<float>(size_ref);
		
	//Parameters Calculation (i.e) 'A' x + 'B' 
	float A = ref_h_sig/img_sig;
	float B = ref_h_avg - (A*img_avg);
	
	vector<float> args;
	args.push_back(A);
	args.push_back(B);
	
	vector<float> scale;
	scale.push_back(1.e-7*A);
	scale.push_back(-1.e-7*B);
	
	vector<float> ref_freq_hist;
	for(int i = 0;i < (3*hist_len);i++)
		ref_freq_hist.push_back((int)ref_freq_bin[i]);
		
	vector<float> data;
	data.push_back(ref_h_diff);
	data.push_back(ref_h_min);
	
	Dict parameter;
	
	/* Parameters displaying the arguments A & B, and the scaling function and the data's */
	parameter["args"] = args;
	parameter["scale"]= scale;
	parameter["data"] = data;
	parameter["ref_freq_bin"] = ref_freq_hist;
	parameter["size_img"]=size_img;
	parameter["hist_len"]=hist_len;
	/* ===================================================== */
	
	return parameter;	
}
	
	
float Util::hist_comp_freq(float PA,float PB,int size_img, int hist_len, EMData *img, vector<float> ref_freq_hist, EMData *mask, float ref_h_diff, float ref_h_min)
{
	float *img_ptr = img->get_data();
	float *mask_ptr = (mask == NULL)?img->get_data():mask->get_data();
		
	int *img_freq_bin = new int[3*hist_len];
	for(int i = 0;i < (3*hist_len);i++)
		img_freq_bin[i] = 0;
	for(int i = 0;i < size_img;i++)
	{
		if(mask_ptr[i] > 0.5f)
		{
			float img_xn = img_ptr[i]*PA + PB;
			int L = static_cast<int>(((img_xn - ref_h_min)/ref_h_diff) * (hist_len-1) + hist_len+1);
			if(L >= 0 && L < (3*hist_len))
				img_freq_bin[L]++;
			
		}
	}; 
	int freq_hist = 0;

	for(int i = 0;i < (3*hist_len);i++)
		freq_hist += (int)pow((float)((int)ref_freq_hist[i] - (int)img_freq_bin[i]),2.f);
	freq_hist = (-freq_hist);
	return freq_hist;
}
//------------------------------------------------------------------------------------------------------------------------------------------------------------------
#define    QUADPI      		        3.141592653589793238462643383279502884197
#define    DGR_TO_RAD    		QUADPI/180
#define    DM(I)         		DM	    [I-1]   
#define    SS(I)         		SS	    [I-1]
Dict Util::CANG(float PHI,float THETA,float PSI)
{
 double CPHI,SPHI,CTHE,STHE,CPSI,SPSI;
 vector<float>   DM,SS;
 
 for(int i =0;i<9;i++)
     DM.push_back(0);
     
 for(int i =0;i<6;i++)
     SS.push_back(0);   
  
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
#define    DM(I)         		DM	    [I-1]  
#define    B(i,j) 			Bptr        [i-1+((j-1)*NSAM)] 
#define    CUBE(i,j,k)                  CUBEptr     [(i-1)+((j-1)+((k-1)*NY3D))*NX3D]

void Util::BPCQ(EMData *B,EMData *CUBE, vector<float> DM)
{

 float  *Bptr = B->get_data(); 
 float  *CUBEptr = CUBE->get_data();
 
 int NSAM,NROW,NX3D,NY3D,NZC,KZ,IQX,IQY,LDPX,LDPY,LDPZ,LDPNMX,LDPNMY,NZ1;
 float DIPX,DIPY,XB,YB,XBB,YBB;
 
 NSAM = B->get_xsize();
 NROW = B->get_ysize();
 NX3D = CUBE->get_xsize();
 NY3D = CUBE->get_ysize();
 NZC = CUBE->get_zsize();


 LDPX   = NX3D/2 +1;
 LDPY   = NY3D/2 +1;
 LDPZ   = NZC/2 +1;
 LDPNMX = NSAM/2 +1;
 LDPNMY = NROW/2 +1;
 NZ1    = 1; 
  
 for(int K=1;K<=NZC;K++)
     {
       KZ=K-1+NZ1;
       for(int J=1;J<=NY3D;J++)
           {
	     XBB = (1-LDPX)*DM(1)+(J-LDPY)*DM(2)+(KZ-LDPZ)*DM(3);
             YBB = (1-LDPX)*DM(4)+(J-LDPY)*DM(5)+(KZ-LDPZ)*DM(6);
              for(int I=1;I<=NX3D;I++)
	          {
		     XB  = (I-1)*DM(1)+XBB;
		     IQX = int(XB+float(LDPNMX));
                     if (IQX <1 || IQX >= NSAM) continue;
		     YB   = (I-1)*DM(4)+YBB;
                     IQY  = int(YB+float(LDPNMY));
                     if (IQY<1 || IQY>=NROW)  continue;
                     DIPX = XB+LDPNMX-IQX;
		     DIPY = YB+LDPNMY-IQY;
 
                    CUBE(I,J,K) = CUBE(I,J,K)+B(IQX,IQY)+DIPY*(B(IQX,IQY+1)-B(IQX,IQY))+DIPX*(B(IQX+1,IQY)-B(IQX,IQY)+DIPY*(B(IQX+1,IQY+1)-B(IQX+1,IQY)-B(IQX,IQY+1)+B(IQX,IQY)));
 	          }
           } 
 
    } 
    
   
} 

#undef DM
#undef B
#undef CUBE

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#define    W(i,j) 			Wptr        [i-1+((j-1)*(Wnx))] 
#define    PROJ(i,j) 		        PROJptr     [i-1+((j-1)*(NNNN))] 
#define    SS(I,J)         		SS	    [I-1 + (J-1)*6]

void Util::WTF(EMData* PROJ,vector<float> SS,float SNR,int K,vector<float> exptable)
{
 int NSAM,NROW,NNNN,NR2,L,JY,KX,NANG;
 float WW,OX,OY,Y;
 
 NSAM = PROJ->get_xsize();
 NROW = PROJ->get_ysize(); 
 NNNN   = NSAM+2-(NSAM%2);
 NR2 = NROW/2;
 
 NANG = int(SS.size())/6; 
  
 EMData* W = new EMData();
 int Wnx = NNNN/2;
 W->set_size(Wnx,NROW,1);
 W->to_zero();
 float *Wptr = W->get_data();
 float *PROJptr = PROJ->get_data(); 
 float indcnst = 1000/2.0;
 // we create look-up table for 1001 uniformly distributed samples [0,2];
 
 for (L=1; L<=NANG; L++) {
      OX = SS(6,K)*SS(4,L)*(-SS(1,L)*SS(2,K)+ SS(1,K)*SS(2,L)) + SS(5,K)*(-SS(3,L)*SS(4,K)+SS(3,K)*SS(4,L)*(SS(1,K)*SS(1,L) + SS(2,K)*SS(2,L)));
      OY = SS(5,K)*SS(4,L)*(-SS(1,L)*SS(2,K)+ SS(1,K)*SS(2,L)) - SS(6,K)*(-SS(3,L)*SS(4,K)+SS(3,K)*SS(4,L)*(SS(1,K)*SS(1,L) + SS(2,K)*SS(2,L)));

      if(OX != 0.0f || OY!=0.0f) { 
	 //int count = 0;
        for(int J=1;J<=NROW;J++) {
	      JY = (J-1);
	      if(JY > NR2) JY=JY-NROW;	       
	      for(int I=1;I<=NNNN/2;I++) {
		          Y =  fabs(OX * (I-1) + OY * JY);
                  if(Y < 2.0f) W(I,J) += exptable[int(Y*indcnst)];//exp(-4*Y*Y);//
		  //if(Y < 2.0f) Wptr[count++] += exp(-4*Y*Y);//exptable[int(Y*indcnst)];//
		  }   
	    }
	  } else { 
	    for(int J=1;J<=NROW;J++) for(int I=1;I<=NNNN/2;I++)  W(I,J) += 1.0f;
 	  }
 }

 PROJ->pad_fft();
 PROJ->do_fft_inplace();
 PROJ->update();
 PROJ->done_data();
 PROJptr = PROJ->get_data();
 
 
 float WNRMinv,temp;
 float osnr = 1.0f/SNR;
 WNRMinv = 1/W(1,1);
 for(int J=1;J<=NROW;J++)
    for(int I=1;I<=NNNN;I+=2) {
         KX          = (I+1)/2;
	     temp        = W(KX,J)*WNRMinv;
	     WW          = temp/(temp*temp + osnr);
	     PROJ(I,J)   *= WW;
         PROJ(I+1,J) *= WW;
       }  

PROJ->do_ift_inplace();
PROJ->postift_depad_corner_inplace();
}

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
 float rad2deg =(180.0/3.1415926);
 float deg2rad = (3.1415926/180.0);
 
 int NSAM,NROW,NNNN,NR2,NANG,L,JY;
  
 NSAM = PROJ->get_xsize();
 NROW = PROJ->get_ysize(); 
 NNNN   = NSAM+2-(NSAM%2);
 NR2 = NROW/2;
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
 
 float THICK=NSAM/DIAMETER/2.0;

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
	  CCN=AMAX1( AMIN1(TMP,1.0) ,-1.0);
	  ALPHA=rad2deg*float(asin(CCN));
	  if (ALPHA>180.0) ALPHA=ALPHA-180.0;
	  if (ALPHA>90.0) ALPHA=180.0-ALPHA;
	  if(ALPHA<1.0E-6) {
          for(int J=1;J<=NROW;J++) for(int I=1;I<=NNNN/2;I++) W(I,J)+=1.0;
    } else {
      FM=THICK/(fabs(sin(ALPHA*deg2rad)));
      CC(1)   = CC(1)/CCN;CC(2)   = CC(2)/CCN;CC(3)   = CC(3)/CCN;
      VV(1)= SS(2,L)*SS(4,L)*CC(3)-SS(3,L)*CC(2);
      VV(2)= SS(3,L)*CC(1)-SS(1,L)*SS(4,L)*CC(3);
      VV(3)= SS(1,L)*SS(4,L)*CC(2)-SS(2,L)*SS(4,L)*CC(1);
      CP(1)   = 0.0;CP(2) = 0.0;
      VP(1)   = 0.0;VP(2) = 0.0;
      
	  CP(1) = CP(1) + RI(1,1)*CC(1) + RI(1,2)*CC(2) + RI(1,3)*CC(3);
	  CP(2) = CP(2) + RI(2,1)*CC(1) + RI(2,2)*CC(2) + RI(2,3)*CC(3);
	  VP(1) = VP(1) + RI(1,1)*VV(1) + RI(1,2)*VV(2) + RI(1,3)*VV(3);
	  VP(2) = VP(2) + RI(2,1)*VV(1) + RI(2,2)*VV(2) + RI(2,3)*VV(3);						
      
      TMP = CP(1)*VP(2)-CP(2)*VP(1);

       //     PREVENT TMP TO BE TOO SMALL, SIGN IS IRRELEVANT
       TMP = AMAX1(1.0E-4,fabs(TMP));
	   float tmpinv = 1/TMP;   
       for(int J=1;J<=NROW;J++) {
	     JY = (J-1);
         if (JY>NR2)  JY=JY-NROW;
         for(int I=1;I<=NNNN/2;I++) {
        		FV     = fabs((JY*CP(1)-(I-1)*CP(2))*tmpinv);
        		RT     = 1.0-FV/FM;
        		W(I,J) += ((RT>0.0)*RT);		 
         }
       } 
      }  
	}

 }
 
 PROJ->pad_fft();
 PROJ->do_fft_inplace();
 PROJ->update();
 PROJ->done_data();
 float *PROJptr = PROJ->get_data();
  
 int KX;
 float WW;
 for(int J=1; J<=NROW; J++)
    for(int I=1; I<=NNNN; I+=2) {
         KX          =  (I+1)/2;
         WW          =  1.0f/W(KX,J);
	     PROJ(I,J)   = PROJ(I,J)*WW;
         PROJ(I+1,J) = PROJ(I+1,J)*WW;
    }  

 PROJ->do_ift_inplace();
 PROJ->postift_depad_corner_inplace();  
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
//-----------------------------------------------------------------------------------------------------------------------
Dict Util::ExpMinus4YSqr(float ymax,int nsamples)
{
  //exp(-16) is 1.0E-7 approximately)
  vector<float> expvect;
  
  double inc = double(ymax)/nsamples;
  double temp;
  for(int i =0;i<nsamples;i++)
     {
      temp = exp((-4*(i*inc)*(i*inc)));
      expvect.push_back(float(temp));  
     }
 expvect.push_back(0.0);
 Dict lookupdict;
 lookupdict["table"] = expvect;
 lookupdict["ymax"] = ymax;
 lookupdict["nsamples"] = nsamples;

  return lookupdict;  
}
//------------------------------------------------------------------------------------------------------------------------- 

float Util::tf(float dzz,float ak,float lambda,float cs,float wgh,float b_factor,float sign)  {
return sin(-M_PI*(dzz*lambda*ak*ak-cs*lambda*lambda*lambda*ak*ak*ak*ak/2.)-wgh)*exp(-b_factor*ak*ak)*sign;
}

EMData *Util::ctf_img(int nx, int ny, int nz,float ps,float dz,float cs,float voltage,float dza, float azz,float wgh,float b_factor, float sign)
{               
	int  lsm;
	double ix,iy,iz;
	int i,j,k;    
	int nr2 ,nl2;
	float dzz,az,ak;
	float scx, scy,scz;	  
	if (nx%2==0) lsm=nx+2; else lsm=nx+1;		     
	float lambda=12.398/pow(voltage *(1022.+voltage),.5);	
	cs=cs*1.0e-7f;    
	wgh = atan(wgh/(1.0-wgh));   
	EMData* ctf_img1 = new EMData();
	ctf_img1->set_size(lsm,ny,nz);
	float freq=1./(2.*ps);		    
	scx=2./nx;
	if(ny<=1) scy=2./ny; else scy=0.0;
	if(nz<=1) scz=2./nz; else scz=0.0;
	nr2=ny/2 ;
	nl2=nz/2 ;
	for ( k=0; k<nz;k++) {
	       if(k>nl2) iz=k-float(nz);
	       for ( j=0; j<ny;j++) { 
	     	     if(j>nr2) iy=j-float(ny);
	     	     for ( i=0;i<lsm/2;i++) {
	     		   ix=i;
	     		   ak=pow(ix*ix*scx*scx+iy*scy*iy*scy+iz*scz*iz*scz,.5)*freq;
	     		   if(ak!=0) az=0.0; else az=M_PI;
	     		   dzz=dz+dza/2.*sin(2*(az-azz*M_PI/180.));
			   (*ctf_img1) (i*2,j,k)=tf(dzz,ak,lambda,cs,wgh,b_factor,sign);
	     		   (*ctf_img1) (i*2+1,j,k)=0.0f;
	     	     }
	     	     
	       }

	}
		if(nx%2==0) ctf_img1->set_fftodd(false); else ctf_img1->set_fftodd(true); 
		ctf_img1->set_complex(true);
	    	ctf_img1->set_ri(1);  
	        if(nx%2==0) ctf_img1->set_attr("npad",2); else  ctf_img1->set_attr("npad",1);
		return ctf_img1;
			 			 
} 		
//Return the 1-D image that contains only pixels from a n-d image selected by a mask

EMData* Util::compress_image_mask(EMData* image, EMData* mask)
{
	/***********
	***get the size of the image for validation purpose
	**************/
	int nx = image->get_xsize(),ny = image->get_ysize(),nz = image->get_zsize();  //Aren't  these  implied?  Please check and let me know, PAP.
	/********
	***Exception Handle 
	*************/
	if(nx != mask->get_xsize() || ny != mask->get_ysize() || nz != mask->get_zsize())
		throw ImageDimensionException("The dimension of the image does not match the dimension of the mask!");

	int i, size = nx*ny*nz;
		
	float* img_ptr = image->get_data();
	float* mask_ptr = mask->get_data();

	int ln=0;  //length of the output image = number of points under the mask.
	for(i = 0;i < size;i++){
		if(mask_ptr[i] > 0.5f) ln++;
	}

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

/* Recreates a n-d image using its compressed 1-D form and the mask */
EMData *Util::reconstitute_image_mask(EMData* image, EMData *mask)
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

	int i,size = nx*ny*nz;			 /* loop counters */
	/* new image declaration */
	EMData *new_image = new EMData();
	new_image->set_size(nx,ny,nz); 		 /* set the size of new image */
	float *new_ptr  = new_image->get_data(); /* set size of the new image */
	float *mask_ptr = mask->get_data();	 /* assign a pointer to the mask image */
	float *img_ptr  = image->get_data();	 /* assign a pointer to the 1D image */
	int count = 0;
	for(i = 0;i < size;i++){
		if(mask_ptr[i] > 0.5f){
			new_ptr[i] = img_ptr[count];
			count++;
		}
		else{
			new_ptr[i] = 0.0f;
		}
	}
	new_image->update();
	
	return new_image;
}	
vector<float> Util::merge_peaks(vector<float> peak1, vector<float> peak2,float p_size)
{	
	vector<float>new_peak;
	int n1=peak1.size()/3;
	float p_size2=p_size*p_size;
	for (int i=0;i<n1;++i)
		{
			vector<float>::iterator it2= peak1.begin()+3*i;
			bool push_back1=true;
			int n2=peak2.size()/3;
			/*cout<<"peak2 size==="<<n2<<"i====="<<i<<endl;
			cout<<"new peak size==="<<new_peak.size()/3<<endl;*/
			
			if(n2 ==0) 
				{
					new_peak.push_back(*it2);
					new_peak.push_back(*(it2+1));
					new_peak.push_back(*(it2+2));
				
				}
			else 
				{						
					int j=0;					
					while (j< n2-1 )
						{								
							vector<float>::iterator it3= peak2.begin()+3*j;
							float d2=((*(it2+1))-(*(it3+1)))*((*(it2+1))-(*(it3+1)))+((*(it2+2))-(*(it3+2)))*((*(it2+2))-(*(it3+2)));							
							if(d2< p_size2 )
								{ 	
									if( (*it2)<(*it3))
										{	
											new_peak.push_back(*it3);
											new_peak.push_back(*(it3+1));
											new_peak.push_back(*(it3+2));
											peak2.erase(it3);
											peak2.erase(it3);
											peak2.erase(it3);
											push_back1=false;
										}
									else
										{
											peak2.erase(it3);
											peak2.erase(it3);
											peak2.erase(it3);
									
										}	
								}
							else
								{
									j=j+1;
								}
								n2=peak2.size()/3;						
						}
					if(push_back1)
						{
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
    for ( i = 0 ; i < n * n ; i++ ) {
       eigvec[i] = covmat[i];
    }		
	
    char NEEDV = 'V';
    char UPLO = 'U';
    int lwork = -1;
    int info = 0;
    float *work, wsize;
  
    //	query to get optimal workspace
    ssyev_(&NEEDV, &UPLO, &n, eigvec, &n, eigval, &wsize, 
           &lwork, &info);
    lwork = (int)wsize;

    work = (float *)calloc(lwork, sizeof(float));
    // 	calculate eigs
    ssyev_(&NEEDV, &UPLO, &n, eigvec, &n, eigval, work, 
           &lwork, &info);
    free(work);
    return info;
    EXITFUNC;
}
