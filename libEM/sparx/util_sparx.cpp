/**
 * $Id$
 */
#include <iostream>

#include "emdata.h"
#include "util.h"

using namespace EMAN;
using namespace std;

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
		float D2 = D*D;
			
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


