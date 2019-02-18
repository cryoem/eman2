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


EMData *EMData::real2FH(float OverSamplekB) // PRB
{
	int nx        = get_xsize();
	int ny        = get_ysize();
	int nz        = get_zsize();
	int Center  = (int) floor( (nx+1.0)/2.0 +.01);
#ifdef DEBUG
	printf("nx=%d, ny=%d, nz=%d Center=%d\n", nx,ny,nz, Center);
#endif	//DEBUG
	float ScalFactor=4.1f;
	gsl_set_error_handler_off();

	if ( (nz==1) && (nx==ny) && (!is_complex())  && (Center*2)==(nx+1)){
#ifdef DEBUG
		printf("entered if \n");fflush(stdout);
#endif	//DEBUG
//		MArray2D ImBW = this ->get_2dview();
		EMData*  ImBW = this ;
		int Size=nx;
		int iMax = (int) floor( (Size-1.0)/2 +.01);
		int CountMax = (iMax+2)*(iMax+1)/2;
		int *PermMatTr  = new int[CountMax];
		float *RValsSorted  = new float[CountMax];
		float *weightofkValsSorted = new float[CountMax];
		int *SizeReturned = new int[1];
		Util::Radialize(PermMatTr, RValsSorted,weightofkValsSorted,Size, SizeReturned);
	  	int RIntMax= SizeReturned[0];

		int mMax = (int) floor( ScalFactor*RValsSorted[RIntMax-1]+10.0);

		int kIntMax=2+ (int) floor( RValsSorted[RIntMax-1]*OverSamplekB);
		float *kVec2Use= new float[kIntMax];
		for (int kk=0; kk<kIntMax; kk++){
			kVec2Use[kk]= ((float) kk)/OverSamplekB;}

		float *krVec= new float[kIntMax*RIntMax];
		int Count=0;
		for (int jk=0; jk<kIntMax; jk++ ){
			for (int jR=0; jR<RIntMax; jR++ ){
				krVec[Count]=2.0f*M_PI*RValsSorted[jR]
					*kVec2Use[jk]/( (float) Size);
				Count++;
//				printf("krVec[%d]=%f \n",Count,krVec[Count-1]);fflush(stdout);
		}} // end building up krVec
		float krVecMin= kVec2Use[1]*RValsSorted[1];
		float krVecMax = krVec[kIntMax*RIntMax-1]+krVecMin;
		int Number2Use = (int) floor(OverSamplekB*krVecMax+1.0);
		float *krVec2Use      = new float[Number2Use+1];
		float *sampledBesselJ = new float[Number2Use+1];
#ifdef DEBUG
		printf("Size=%d, iMax=%d, SizeReturned=%d, RIntMax=%d, \n"
		      "mMax=%d, kIntMax=%d, krVecMin=%f, krVecMax=%f,  Number2Use=%d  \n\n",
			Size, iMax, SizeReturned[0], RIntMax, mMax, kIntMax,
			       krVecMin,krVecMax,Number2Use);fflush(stdout);
#endif	//DEBUG
		for (int jkr=0; jkr<= Number2Use; jkr++) {
			krVec2Use[jkr] =((float)jkr)*krVecMax/
			            ((float)Number2Use);
//			printf("krVec2Use[%d]=%f \n",jkr+1,krVec2Use[jkr]);fflush(stdout);
		}


		EMData* rhoOfkmB = copy(); // glibc detected ** malloc(); memory corruption
//		printf("finished O \n");fflush(stdout);
		rhoOfkmB->set_size(2*(mMax+1),kIntMax);
		rhoOfkmB->to_zero();
//		MArray2D rhoOfkmB = FH->get_2dview();

		int CenterM= Center-1; // to convert from Matlab to C++
		std::complex <float> *rhoOfRandmTemp = new std::complex <float>[RIntMax];
		std::complex <float> rhoTemp;

		int PCount=0;

		for (int m=0; m <=mMax; m++){
		//    if m==mMax, tic, end
			std::complex <float> tempF(0.0f,-1.0f);
			std::complex <float> overallFactor = (std::complex <float>) pow(tempF,m);  //(-i)^m ;  % I dropped off the 2 pi
			std::complex <float> mI(0.0f,static_cast<float>(m));
			for (int ii=0; ii< RIntMax; ii++){ rhoOfRandmTemp[ii]=0;}
			for (int jx=0; jx <Center ; jx++) {
				for (int jy=0; jy <=jx; jy++){
					float fjx=float(jx);
					float fjy= float(jy);
          				Count = (jx*jx+jx)/2 +1 +jy;
					PCount = PermMatTr[Count-1];
//					printf("PCount=%d, Count=%d \n", PCount, Count);
  				        rhoTemp =  std::complex <float> ((*ImBW)(CenterM+jx,CenterM+jy)) *exp(mI* std::complex <float> (atan2(+fjy,+fjx)))
				         +   std::complex <float> ((*ImBW)(CenterM+jx,CenterM-jy)) * exp(mI*std::complex <float>(atan2(-fjy,+fjx)))
				         +   std::complex <float> ((*ImBW)(CenterM-jx,CenterM+jy)) * exp(mI*std::complex <float>(atan2(+fjy,-fjx)))
				         +   std::complex <float> ((*ImBW)(CenterM-jx,CenterM-jy)) * exp(mI*std::complex <float>(atan2(-fjy,-fjx)))
			               	 +   std::complex <float> ((*ImBW)(CenterM+jy,CenterM+jx)) * exp(mI*std::complex <float>(atan2(+fjx,+fjy)))
					 +   std::complex <float> ((*ImBW)(CenterM+jy,CenterM-jx)) * exp(mI*std::complex <float>(atan2(-fjx,+fjy)))
					 +   std::complex <float> ((*ImBW)(CenterM-jy,CenterM+jx)) * exp(mI*std::complex <float>(atan2(+fjx,-fjy)))
					 +   std::complex <float> ((*ImBW)(CenterM-jy,CenterM-jx)) * exp(mI*std::complex <float>(atan2(-fjx,-fjy)));
            				if (((jx+jy)==0)&&(m>0) ){
						rhoTemp=0;}
//			printf("m=%d, jx=%d, jy=%d, rhoTemp= %f+ %f i\n", m,jx,jy,(rhoTemp.real()), (rhoTemp.imag()) );fflush(stdout);
//			{" %f,%f %f,%f %f,%f %f,%f \n",
//			       ImBW[CenterM+jx][CenterM+jy] ,ImBW[CenterM+jx][CenterM-jy]  , ImBW[CenterM-jx][CenterM+jy] ,ImBW[CenterM-jx][CenterM-jy],
//			       ImBW[CenterM+jy][CenterM+jx] ,ImBW[CenterM+jy][CenterM-jx]  , ImBW[CenterM-jy][CenterM+jx] ,ImBW[CenterM-jy][CenterM-jx]);
            				rhoOfRandmTemp[PCount-1] +=
				            rhoTemp/((float)pow(2.,(int)( (jx==0)  +(jy==0)+ (jy==jx))));

			}} // end walk through lattice
//			printf("\n m=%d rhoOfRandmTemp" ,m  );fflush(stdout);
//			for (int ss=0; ss< RIntMax; ss++){
//				printf(" %3.1f+ %3.1fi \t",(rhoOfRandmTemp[ss].real()), (rhoOfRandmTemp[ss].imag())   );fflush(stdout);}

// calculate product

			float tempp;
//			printf("\n m=%d sampledBesselJ" ,m  );fflush(stdout);
			for (int st=0; st<= Number2Use; st++){
				tempp=krVec2Use[st];
				sampledBesselJ[st] = static_cast<float>(gsl_sf_bessel_Jn(m,tempp));
//				printf(" %3.2f  \t",sampledBesselJ[st]   );fflush(stdout);
			} // good so far

//			sampledBesselJ  = BesselJ(m,krVec2Use);
			float *tempMB = new float [kIntMax*RIntMax];
			Util::spline_mat(krVec2Use, sampledBesselJ, Number2Use+1,krVec,tempMB,kIntMax*RIntMax );
//			printf("\n tempMB m=%d y2" ,m  );fflush(stdout);
			std::complex <float> *rowV = new std::complex <float> [kIntMax];

//			for (int st=0; st< kIntMax*RIntMax; st++){printf(" %3.2f  \t",tempMB[st]   );fflush(stdout);} // good so far

//   tempMB,krVec is in blocks of RIntMax
//			printf("\n rowV m=%d \t" ,m  );fflush(stdout);
			for (int st=0; st < kIntMax; st++) {
					rowV[st]=0;
					for (int sv=0; sv < RIntMax; sv++) {
						rowV[st]+=  rhoOfRandmTemp[sv] *tempMB[sv+st*RIntMax];
					}
					 rowV[st] *= overallFactor;
//					printf(" %1.3f +%1.3fi \t" , rowV[st].real(), rowV[st].imag() );fflush(stdout);
			}
			for (int st=0; st < kIntMax; st++) {
				(*rhoOfkmB)(2*m  ,st) = rowV[st].real();
				(*rhoOfkmB)(2*m+1,st) = rowV[st].imag();
			}
// 			rowV = overallFactor*rhoOfRandmTemp*tempMBB;
//			rhoOfkmB(m+1,1:kIntMax) = rowV ;

//			if m==mMax, toc, end

// %'final interpolation'
// %     rhoOfkm(m+1,:) = spline(kVec2Use,rowV,RValsSorted); ;


		} // ends m loop

		update();
		rhoOfkmB-> update();
		rhoOfkmB->set_complex(true);
		if(rhoOfkmB->get_ysize()==1 && rhoOfkmB->get_zsize()==1) {
			rhoOfkmB->set_complex_x(true);
		}
    	rhoOfkmB->set_ri(true);
    	rhoOfkmB->set_FH(true);
    	rhoOfkmB->set_fftodd(true);
		return rhoOfkmB;
	} else {
		LOGERR("2D real square odd image expected.");
		throw ImageFormatException("2D real square odd image expected.");
	}
}

EMData *EMData::copy_empty_head() const
{
	ENTERFUNC;
	EMData *ret = new EMData();
	ret->attr_dict = attr_dict;
	ret->flags = flags;
	ret->all_translation = all_translation;
	ret->path = path;
	ret->pathnum = pathnum;

// should these be here? d.woolford I did not comment them out, merely place them here (commented out) to draw attention
// 	ret->xoff = xoff;
// 	ret->yoff = yoff;
// 	ret->zoff = zoff;
// 	ret->changecount = changecount;

	ret->update();

	EXITFUNC;
	return ret;
}


EMData *EMData::FH2F(int Size, float OverSamplekB, int IntensityFlag)  // PRB
{
	int nx=get_xsize();
	int ny=get_ysize();
	int nz=get_zsize();
	float ScalFactor=4.1f;
	int Center = (int) floor((Size+1.0)/2.0 +.1);
	int CenterM= Center-1;
	int CountMax = (Center+1)*Center/2;

	int     *PermMatTr           = new int[CountMax];
	float  *RValsSorted         = new float[CountMax];
	float  *weightofkValsSorted = new float[CountMax];
	int      *SizeReturned        = new int[1];
	Util::Radialize(PermMatTr, RValsSorted,weightofkValsSorted,Size, SizeReturned);
	int RIntMax= SizeReturned[0];  // replaces CountMax; the latter should now never be used.
//	kVec2Use = (0:1/OverSamplek:RValsSorted(RIntMax)+1/OverSamplek); %   in pixels  (otherwise need *2*pi/Size)

	int   mMax = (int) floor( ScalFactor*RValsSorted[RIntMax-1]+10.0);

	int    kIntMax  = 2+ (int) floor( RValsSorted[RIntMax-1]*OverSamplekB);
	float *kVec2Use = new float[kIntMax];
	for (int kk=0; kk<kIntMax; kk++){
		kVec2Use[kk]= ((float) kk)/OverSamplekB;}



#ifdef DEBUG
	printf("nx=%d, ny=%d, nz=%d Center=%d mMax=%d CountMax=%d kIntMax=%d Centerm1=%d  Size=%d\n\n",
	    nx,ny,nz, Center, mMax, CountMax, kIntMax,  CenterM, Size);
#endif

	EMData * rhoOfkmB = this;

//     check mMax's are equal
//     check kIntMax's are equal

	if ( (nx==2*(mMax+1)) && (ny==kIntMax) &&(nz==1) ) {

	EMData *rhoOfkandm = copy();
	rhoOfkandm ->set_size(2*(mMax+1),RIntMax);
	rhoOfkandm ->to_zero();
//	MArray2D rhoOfkandm = tempCopy->get_2dview();  % Just changed Nov 20 2005
//	printf("rhoOfkandm \n");
	for (int mr=0; mr <2*(mMax+1); mr++){
		float *Row= new float[kIntMax];
		float *RowOut= new float[RIntMax];
		for (int ii=0; ii<kIntMax; ii++){ Row[ii]=(*rhoOfkmB)(mr,ii);}
		Util::spline_mat(kVec2Use, Row, kIntMax,  RValsSorted, RowOut, RIntMax );
		for (int ii=0; ii<RIntMax; ii++){
			(*rhoOfkandm)(mr,ii) = RowOut[ii];
//			printf("%3.3f  ",RowOut[ii]);
		}
//		printf(" \n");
//		rhoOfkandm(m+1,:) = spline(kVec2Use,rhoOfkmBReIm(m+1,1:kIntMax),kIntMax,RValsSorted);
	}
	rhoOfkandm ->update();

//          So far so good PRB ....

	EMData* outCopy = rhoOfkandm ->copy();
	outCopy->set_size(2*Size,Size,1);
	outCopy->to_zero();
//	MArray2D ImBWfftRm = outCopy->get_2dview();

	int Count =0, kInt, kIntm1;
	std::complex <float> ImfTemp;
	float kValue, thetak;

	for (int jkx=0; jkx <Center; jkx++) { // These index the outputted picture
		for (int jky=0; jky<=jkx; jky++){
			kInt = PermMatTr[Count];
			kIntm1= kInt-1;
			Count++;
			float fjkx = float(jkx);
			float fjky = float(jky);

			kValue = std::sqrt(fjkx*fjkx +  fjky*fjky )  ;
//        		mMaxR= floor(ScalFactor*kValue +10);

 //                   How many copies

			thetak = atan2(fjky,fjkx);
			ImfTemp = (*rhoOfkandm)(0, kIntm1) ;
        		for (int mm= 1; mm <mMax;mm++) {  // The index for m
				std::complex <float> fact(0,-mm*thetak);
				std::complex <float> expfact= exp(fact);
				std::complex <float> tempRho((*rhoOfkandm)(2*mm,kIntm1),(*rhoOfkandm)(2*mm+1,kIntm1));
				float mmFac = float(1-2*(mm%2));
				if (IntensityFlag==1){ mmFac=1;}
				ImfTemp +=   expfact * tempRho + mmFac  *conj(expfact*tempRho);//pow(float(-1),mm)
        		}
 			(*outCopy)(2*(CenterM+jkx),CenterM+jky)   = ImfTemp.real();
			(*outCopy)(2*(CenterM+jkx)+1,CenterM+jky) = ImfTemp.imag();
//			printf("jkx=%d, jky=%d; %f + %f i \n",jkx,jky,ImfTemp.real(), ImfTemp.imag());

			if (jky>0) {
				thetak = atan2(-fjky,fjkx);
				ImfTemp = (*rhoOfkandm)(0,kIntm1);
				for (int mm= 1; mm<mMax; mm++) { // The index for m
					std::complex <float> fact(0,-mm*thetak);
					std::complex <float> expfact= exp(fact);
					std::complex <float> tempRho((*rhoOfkandm)(2*mm,kIntm1), (*rhoOfkandm)(2*mm+1,kIntm1));
					float mmFac = float(1-2*(mm%2));
					if (IntensityFlag==1){ mmFac=1;}
					ImfTemp +=   expfact * tempRho +  mmFac  *conj(expfact*tempRho);
				}
				(*outCopy)(2*(CenterM+jkx),CenterM-jky)  = ImfTemp.real();

				(*outCopy)(2*(CenterM+jkx)+1,CenterM-jky) = ImfTemp.imag();
			}

			if (jkx>0) {
            			thetak = atan2(fjky,-fjkx);
				ImfTemp = (*rhoOfkandm)(0,kIntm1);
				for (int mm= 1; mm<mMax; mm++) { // The index for m
					std::complex <float> fact(0,-mm*thetak);
					std::complex <float> expfact= exp(fact);
					std::complex <float> tempRho((*rhoOfkandm)(2*mm,kIntm1), (*rhoOfkandm)(2*mm+1,kIntm1));
					float mmFac = float(1-2*(mm%2));
					if (IntensityFlag==1){ mmFac=1;}
					ImfTemp +=   expfact * tempRho +  mmFac *conj(expfact*tempRho);
				}
				(*outCopy)(2*(CenterM-jkx)  ,CenterM+jky) = ImfTemp.real();
				(*outCopy)(2*(CenterM-jkx)+1,CenterM+jky) = ImfTemp.imag();
			}

 			if (jkx>0 && jky>0) {
				thetak = atan2(-fjky,-fjkx);
				ImfTemp = (*rhoOfkandm)(0 , kIntm1);
				for (int mm= 1; mm<mMax; mm++) {  // The index for m
					std::complex <float> fact(0,-mm*thetak);
					std::complex <float> expfact= exp(fact);
					std::complex <float> tempRho( (*rhoOfkandm)(2*mm,kIntm1),(*rhoOfkandm)(2*mm+1,kIntm1) );
					float mmFac = float(1-2*(mm%2));
					if (IntensityFlag==1){ mmFac=1;}
					ImfTemp +=   expfact * tempRho +  mmFac *conj(expfact*tempRho);
				}
				(*outCopy)(2*(CenterM-jkx)  ,CenterM-jky) = ImfTemp.real();
				(*outCopy)(2*(CenterM-jkx)+1,CenterM-jky) = ImfTemp.imag();
			}

			if (jky< jkx) {
				thetak = atan2(fjkx,fjky);
				ImfTemp = (*rhoOfkandm)(0,kIntm1);
				for (int mm= 1; mm<mMax; mm++){ // The index for m
					std::complex <float> fact(0,-mm*thetak);
					std::complex <float> expfact= exp(fact);
					std::complex <float> tempRho((*rhoOfkandm)(2*mm,kIntm1),(*rhoOfkandm)(2*mm+1,kIntm1));
					float mmFac = float(1-2*(mm%2));
					if (IntensityFlag==1){ mmFac=1;}
					ImfTemp +=   expfact * tempRho +  mmFac *conj(expfact*tempRho);
				}
				(*outCopy)(2*(CenterM+jky)  ,CenterM+jkx) = ImfTemp.real();
				(*outCopy)(2*(CenterM+jky)+1,CenterM+jkx) = ImfTemp.imag();

				if (jky>0){
					thetak = atan2(fjkx,-fjky);
					ImfTemp = (*rhoOfkandm)(0, kIntm1);
					for (int mm= 1; mm <mMax; mm++) { // The index for m
						std::complex <float> fact(0,-mm*thetak);
						std::complex <float> expfact= exp(fact);
						std::complex <float> tempRho((*rhoOfkandm)(2*mm,kIntm1),(*rhoOfkandm)(2*mm+1,kIntm1));
					float mmFac = float(1-2*(mm%2));
					if (IntensityFlag==1){ mmFac=1;}
						ImfTemp +=  expfact * tempRho +  mmFac *conj(expfact*tempRho);
					}
					(*outCopy)(2*(CenterM-jky)  ,CenterM+jkx) = ImfTemp.real();
					(*outCopy)(2*(CenterM-jky)+1,CenterM+jkx) = ImfTemp.imag();
				}

				 if (jkx>0) {
					 thetak = atan2(-fjkx,fjky);
					 ImfTemp = (*rhoOfkandm)(0,kIntm1);
					for (int mm= 1; mm <mMax; mm++) { // The index for m
						std::complex <float> fact(0,-mm*thetak);
						std::complex <float> expfact= exp(fact);
						std::complex <float> tempRho((*rhoOfkandm)(2*mm,kIntm1),(*rhoOfkandm)(2*mm+1,kIntm1));
						float mmFac = float(1-2*(mm%2));
						if (IntensityFlag==1){ mmFac=1;}
						ImfTemp +=  expfact * tempRho +  mmFac *conj(expfact*tempRho);
 					}
					(*outCopy)(2*(CenterM+jky)  ,CenterM-jkx) = ImfTemp.real();
					(*outCopy)(2*(CenterM+jky)+1,CenterM-jkx) = ImfTemp.imag();
 				}

	 			if (jkx>0 && jky>0) {
					thetak = atan2(-fjkx,-fjky);
					ImfTemp = (*rhoOfkandm)(0,kIntm1) ;
					for (int mm= 1; mm <mMax; mm++) { // The index for m
						std::complex <float> fact(0,-mm*thetak);
						std::complex <float> expfact= exp(fact);
						std::complex <float> tempRho((*rhoOfkandm)(2*mm,kIntm1) ,(*rhoOfkandm)(2*mm+1,kIntm1) );
						float mmFac = float(1-2*(mm%2));
						if (IntensityFlag==1){ mmFac=1;}
						ImfTemp +=  expfact * tempRho +  mmFac *conj(expfact*tempRho);
					}
					(*outCopy)(2*(CenterM-jky)  ,CenterM-jkx) = ImfTemp.real();
					(*outCopy)(2*(CenterM-jky)+1,CenterM-jkx) = ImfTemp.imag();
 				}
 			} // ends jky <jkx


		} // ends jky
	} // ends jkx
	outCopy->update();
	outCopy->set_complex(true);
	if(outCopy->get_ysize()==1 && outCopy->get_zsize()==1) outCopy->set_complex_x(true);
	outCopy->set_ri(true);
	outCopy->set_FH(false);
	outCopy->set_fftodd(true);
	outCopy->set_shuffled(true);
	return outCopy;
	} else {
		LOGERR("can't be an FH image not this size");
		throw ImageFormatException("something strange about this image: not a FH");

	}
}  // ends FH2F


EMData *EMData::FH2Real(int Size, float OverSamplekB, int)  // PRB
{
	EMData* FFT= FH2F(Size,OverSamplekB,0);
	FFT->process_inplace("xform.fourierorigin.tocorner");
	EMData* eguess= FFT ->do_ift();
	return eguess;
}  // ends FH2F

float dist(int lnlen, const float* line_1, const float* line_2)
{
	float dis2=0.0;
	for( int i=0; i < lnlen; ++i) {
		float tmp = line_1[i] - line_2[i];
		dis2 += tmp*tmp;
	}
	//return static_cast<float>(std::sqrt(dis2));
	return dis2;
}

float dist_r(int lnlen, const float* line_1, const float* line_2)
{
	double dis2 = 0.0;
	for( int i=0; i < lnlen; ++i ) {
		float tmp = line_1[lnlen-1-i] - line_2[i];
		dis2 += tmp*tmp;
	}
	return static_cast<float>(std::sqrt(dis2));
}

/*
float EMData::cm_euc(EMData* sinoj, int n1, int n2, float alpha1, float alpha2)
{
    int lnlen = get_xsize();

	Assert( n1 >=0 && n1 < get_ysize() );
	Assert( n2 >=0 && n2 < get_ysize() );
	Assert( alpha1>=0.0 && alpha1 < 360.0 );
	Assert( alpha2>=0.0 && alpha2 < 360.0 );

	float* line_1 = get_data() + n1*lnlen;
	float* line_2 = sinoj->get_data() + n2*lnlen;
	float just = (alpha1-180.0f)*(alpha2-180.0f);
	if( just > 0.0 ) return dist(lnlen, line_1, line_2);

	if( just == 0.0 ) {
		float dist_1 = dist(lnlen, line_1, line_2);
		float dist_2 = dist_r(lnlen, line_1, line_2);
#ifdef	_WIN32
		return _cpp_min(dist_1, dist_2);
#else
		return std::min(dist_1, dist_2);
#endif	//_WIN32
	}

	Assert( (alpha1-180.0)*(alpha2-180.0) < 0.0 );
	return dist_r(lnlen, line_1, line_2);
}
*/

float EMData::cm_euc(EMData* sinoj, int n1, int n2)
{
	int lnlen = get_xsize();
	float* line_1 = get_data() + n1 * lnlen;
	float* line_2 = sinoj->get_data() + n2 * lnlen;
	return dist(lnlen, line_1, line_2);
}

EMData* EMData::rotavg() {

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

EMData* EMData::rotavg_i() {

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

	avg1D = rotavg();
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


EMData* EMData::mult_radial(EMData* radial) {

	ENTERFUNC;
	if ( ny == 1 && nz == 1 ) {
		LOGERR("Input image must be 2-D or 3-D!");
		throw ImageDimensionException("Input image must be 2-D or 3-D!");
	}

	EMData* result = this->copy_head();

	result->to_zero();
	result->set_array_offsets(-nx/2, -ny/2, -nz/2);
	this->set_array_offsets(-nx/2, -ny/2, -nz/2);
	int rmax = radial->get_xsize();
	int i, j, k, ir;
	float r;
	for ( k = -nz/2; k < nz/2+nz%2; k++) {
		for ( j = -ny/2; j < ny/2+ny%2; j++) {
			for ( i = -nx/2; i < nx/2+nx%2; i++)  {
				r = std::sqrt(float(k*k) + float(j*j) + float(i*i));
				ir = int(r);
				if(ir < rmax-1)  (*result)(i,j,k) = (*this)(i,j,k) * ((*radial)(ir)+((*radial)(ir+1)-(*radial)(ir))*(r - float(ir)));
			}
		}
	}
	result->update();
	result->set_array_offsets(0,0,0);
	this->set_array_offsets(0,0,0);
	EXITFUNC;
	return result;
}

#define rdata(i,j,k) rdata[(i-1)+((j-1)+(k-1)*ny)*(size_t)nx]
#define square(x) ((x)*(x))
vector<float> EMData::cog() {

	vector<float> cntog;
	int ndim = get_ndim();
	int i=1,j=1,k=1;
	float val,sum1=0.f,MX=0.f,RG=0.f,MY=0.f,MZ=0.f,r=0.f;

	if (ndim == 1) {
		for ( i = 1;i <= nx; i++) {
			val   = rdata(i,j,k);
			sum1 += val;
			MX   += ((i-1)*val);
		}
		MX=(MX/sum1);
		for ( i = 1;i <= nx; i++) {
			val   = rdata(i,j,k);
			sum1 += val;
			RG   += val*(square(MX - (i-1)));
		}
		RG=std::sqrt(RG/sum1);
		MX=MX-(nx/2);
		cntog.push_back(MX);
		cntog.push_back(RG);
#ifdef _WIN32
		cntog.push_back((float)Util::round(MX));
#else
		cntog.push_back(round(MX));
#endif	//_WIN32
	} else if (ndim == 2) {
		for (j=1;j<=ny;j++) {
			for (i=1;i<=nx;i++) {
				val = rdata(i,j,k);
				sum1 += val;
				MX   += ((i-1)*val);
				MY   += ((j-1)*val);
			}
		}
		MX=(MX/sum1);
		MY=(MY/sum1);
		sum1=0.f;
		RG=0.f;
		for (j=1;j<=ny;j++) {
			r = (square(MY-(j-1)));
			for (i=1;i<=nx;i++) {
				val = rdata(i,j,k);
				sum1 += val;
				RG   += val*(square(MX - (i-1)) + r);
			}
		}
		RG = std::sqrt(RG/sum1);
		MX = MX - nx/2;
		MY = MY - ny/2;
		cntog.push_back(MX);
		cntog.push_back(MY);
		cntog.push_back(RG);
#ifdef _WIN32
		cntog.push_back((float)Util::round(MX));cntog.push_back((float)Util::round(MY));
#else
		cntog.push_back(round(MX));cntog.push_back(round(MY));
#endif	//_WIN32
	} else {
		for (k = 1;k <= nz;k++) {
			for (j=1;j<=ny;j++) {
				for (i=1;i<=nx;i++) {
					val = rdata(i,j,k);
					sum1 += val;
					MX += ((i-1)*val);
					MY += ((j-1)*val);
					MZ += ((k-1)*val);
				}
			}
		}
		MX = MX/sum1;
		MY = MY/sum1;
		MZ = MZ/sum1;
		sum1=0.f;
		RG=0.f;
		for (k = 1;k <= nz;k++) {
			for (j=1;j<=ny;j++) {
				float r = (square(MY-(j-1)) + square(MZ - (k-1)));
				for (i=1;i<=nx;i++) {
					val = rdata(i,j,k);
					sum1 += val;
					RG   += val*(square(MX - (i-1)) + r);
				}
			}
		}
		RG = std::sqrt(RG/sum1);
		MX = MX - nx/2;
		MY = MY - ny/2;
		MZ = MZ - nz/2;
		cntog.push_back(MX);
		cntog.push_back(MY);
		cntog.push_back(MZ);
		cntog.push_back(RG);
#ifdef _WIN32
		cntog.push_back((float)Util::round(MX));cntog.push_back((float)Util::round(MY));cntog.push_back((float)Util::round(MZ));
#else
		cntog.push_back(round(MX));cntog.push_back(round(MY));cntog.push_back(round(MZ));
#endif	//_WIN32
	}
	return cntog;
}
#undef square
#undef rdata





vector < float >EMData::calc_fourier_shell_correlation(EMData * with, float w)
{
	ENTERFUNC;

/*
 ******************************************************
 *DISCLAIMER
 * 08/16/05 P.A.Penczek
 * The University of Texas
 * Pawel.A.Penczek@uth.tmc.edu
 * Please do not modify the content of calc_fourier_shell_correlation
 ******************************************************/
/*
Fourier Ring/Shell Correlation
Purpose: Calculate CCF in Fourier space as a function of spatial frequency
         between a pair of 2-3D images.
Method: Calculate FFT (if needed), calculate FSC.
Input:  f - real or complex 2-3D image
	g - real or complex 2-3D image
        w - float ring width
Output: 2D 3xk real image.
        k - length of FSC curve, depends on dimensions of the image and ring width
	1 column - FSC,
	2 column - normalized frequency [0,0.5]
	3 column - currently n /error of the FSC = 1/sqrt(n),
                     where n is the number of Fourier coefficients within given shell
*/
	int needfree=0, kz, ky, ii;
	float  argx, argy, argz;

	if (!with) {
		throw NullPointerException("NULL input image");
	}


	EMData *f = this;
	EMData *g = with;

	int nx  = f->get_xsize();
	int ny  = f->get_ysize();
	int nz  = f->get_zsize();

	if (ny==0 && nz==0) {
		throw ImageFormatException( "Cannot calculate FSC for 1D images");
	}

	if (f->is_complex()) nx = (nx - 2 + f->is_fftodd()); // nx is the real-space size of the input image
	int lsd2 = (nx + 2 - nx%2) ; // Extended x-dimension of the complex image

//  Process f if real
	EMData* fpimage = NULL;
	if (f->is_complex()) fpimage = f;
	else {
		fpimage= f->norm_pad(false, 1);
		fpimage->do_fft_inplace();
		needfree|=1; // Extend and do the FFT if f is real
	} 

//  Process g if real
	EMData* gpimage = NULL;
	if (g->is_complex()) gpimage = g;
	else {
		gpimage= g->norm_pad(false, 1);
		gpimage->do_fft_inplace();
		needfree|=2;  // Extend and do the FFT if g is real
	}

	float *d1 = fpimage->get_data();
	float *d2 = gpimage->get_data();

	int nx2 = nx/2;
	int ny2 = ny/2;
	int nz2 = nz/2;

	float dx2 = 1.0f/float(nx2)/float(nx2);
	float dy2 = 1.0f/float(ny2)/float(ny2);

	float dz2 = 1.0f/std::max(float(nz2),1.0f)/std::max(float(nz2),1.0f);
	int inc = Util::round(nx2/w);

	double* ret = new double[inc+1];
	double* n1 = new double[inc+1];
	double* n2 = new double[inc+1];
	float*  lr = new float[inc+1];
	for (int i = 0; i <= inc; i++) {
		ret[i] = 0.0f; n1[i] = 0.0f; n2[i] = 0.0f; lr[i]=0.0f;
	}

	for (int iz = 0; iz <= nz-1; iz++) {
		if(iz>nz2) kz=iz-nz; else kz=iz; argz = float(kz*kz)*dz2;
		for (int iy = 0; iy <= ny-1; iy++) {
			if(iy>ny2) ky=iy-ny; else ky=iy; argy = argz + float(ky*ky)*dy2;
			for (int ix = 0; ix <= lsd2-1; ix+=2) {
			// Skip Friedel related values
				if (ix>0 || (kz>=0 && (ky>=0 || kz!=0))) {
					argx = 0.5f*std::sqrt(argy + float(ix*ix)*0.25f*dx2);
					int r = Util::round(inc*2*argx);
					if(r <= inc) {
						ii = ix + (iy  + iz * ny)* lsd2;
						ret[r] += d1[ii] * double(d2[ii]) + d1[ii + 1] * double(d2[ii + 1]);
						n1[r]  += d1[ii] * double(d1[ii]) + d1[ii + 1] * double(d1[ii + 1]);
						n2[r]  += d2[ii] * double(d2[ii]) + d2[ii + 1] * double(d2[ii + 1]);
						lr[r]  += 2.0f;
					}
				}
			}
		}
	}

	int  linc = 0;
	for (int i = 0; i <= inc; i++) if(lr[i]>0) linc++;

	vector<float> result(linc*3);

	ii = -1;
	for (int i = 0; i <= inc; i++) {
		if(lr[i]>0 && n1[i]>0.0f && n2[i]>0.0f) {
			ii++;
			result[ii]        = float(i)/float(2*inc);
			result[ii+linc]   = float(ret[i] / (std::sqrt(n1[i] * n2[i])));
			result[ii+2*linc] = lr[i]  /*1.0f/sqrt(float(lr[i]))*/;
		}
		else {
			result[i]           = 0.0f;
			result[i+inc+1]     = 0.0f;
			result[i+2*(inc+1)] = 0.0f;}
	}

	if (needfree&1) {
		if (fpimage) {
			delete fpimage;
			fpimage = 0;
		}
	}
	if (needfree&2) {
		if (gpimage) {
			delete gpimage;
			gpimage = 0;
		}
	}
	delete[] ret; delete[]  n1; delete[]  n2; delete[]  lr;

	EXITFUNC;
	return result;
}


vector < float >EMData::scale_factors(EMData * with, int beg, int end)
{
	ENTERFUNC;

/*
 ******************************************************
 *DISCLAIMER
 * 04/20/14 P.A.Penczek
 * The University of Texas
 * Pawel.A.Penczek@uth.tmc.edu
 * Please do not modify
 ******************************************************/
/*
*/
	int needfree=0, kz, ky, ii;
	float  argx, argy, argz;

	if (!with) {
		throw NullPointerException("NULL input image");
	}


	EMData *f = this;
	EMData *g = with;

	int nx  = f->get_xsize();
	int ny  = f->get_ysize();
	int nz  = f->get_zsize();

	if (ny==0 && nz==0) {
		throw ImageFormatException( "Cannot calculate for 1D images");
	}

	if (f->is_complex()) nx = (nx - 2 + f->is_fftodd()); // nx is the real-space size of the input image
	int lsd2 = (nx + 2 - nx%2) ; // Extended x-dimension of the complex image

//  Process f if real
	EMData* fpimage = NULL;
	if (f->is_complex()) fpimage = f;
	else {
		fpimage= f->norm_pad(false, 1); 
		fpimage->do_fft_inplace();
		needfree|=1; // Extend and do the FFT if f is real
	} 

//  Process g if real
	EMData* gpimage = NULL;
	if (g->is_complex()) gpimage = g;
	else {
		gpimage= g->norm_pad(false, 1);
		gpimage->do_fft_inplace();
		needfree|=2;  // Extend and do the FFT if g is real
	}

	float *d1 = fpimage->get_data();
	float *d2 = gpimage->get_data();

	int nx2 = nx/2;
	int ny2 = ny/2;
	int nz2 = nz/2;

	float dx2 = 1.0f/float(nx2)/float(nx2);
	float dy2 = 1.0f/float(ny2)/float(ny2);


	float dz2 = 1.0f/std::max(float(nz2),1.0f)/std::max(float(nz2),1.0f);
	int inc = Util::round(float(std::max(std::max(nx2,ny2),nz2) ));
	int len = end - beg + 1;

	double* ret1 = new double[len];
	double* ret2 = new double[len];
	float*  lr   = new float[len];
	for (int i = 0; i < len; i++) {
		ret1[i] = 0.0f; ret2[i] = 0.0f; lr[i]=0;
	}

	for (int iz = 0; iz <= nz-1; iz++) {
		if(iz>nz2) kz=iz-nz; else kz=iz; argz = float(kz*kz)*dz2;
		for (int iy = 0; iy <= ny-1; iy++) {
			if(iy>ny2) ky=iy-ny; else ky=iy; argy = argz + float(ky*ky)*dy2;
			for (int ix = 0; ix <= lsd2-1; ix+=2) {
			// Skip Friedel related values
				if (ix>0 || (kz>=0 && (ky>=0 || kz!=0))) {
					argx = 0.5f*std::sqrt(argy + float(ix*ix)*0.25f*dx2);
					int r = Util::round(inc*2*argx);
					if((r >= beg) && (r <= end)) {
						ii = ix + (iy  + iz * ny)* lsd2;
						ret1[r-beg] += d1[ii] * double(d2[ii]) + d1[ii + 1] * double(d2[ii + 1]);
						ret2[r-beg] += d2[ii] * double(d2[ii]) + d2[ii + 1] * double(d2[ii + 1]);
						lr[r-beg]   += 2.0f;
					}
				}
			}
		}
	}

	vector<float> result(len*2);

	for (int i = 0; i < len; i++) {
			result[i]       = ret1[i]/lr[i];
			result[i+len]   = ret2[i]/lr[i];
	}

	delete[] ret1; delete[] ret2; delete[]  lr;

	if (needfree&1) {
		if (fpimage) {
			delete fpimage;
			fpimage = 0;
		}
	}
	if (needfree&2) {
		if (gpimage) {
			delete gpimage;
			gpimage = 0;
		}
	}

	EXITFUNC;
	return result;
}

EMData* EMData::symvol(string symString) {
	ENTERFUNC;
	int nsym = Transform::get_nsym(symString); // number of symmetries
	Transform sym;
	// set up output volume
	EMData *svol = new EMData;
	svol->set_size(nx, ny, nz);
	svol->to_zero();
	// actual work -- loop over symmetries and symmetrize
	vector<Transform> rm = sym.get_sym_proj(symString);
	for (int isym = 0; isym < nsym; isym++)  this -> rot_scale_trans_background(rm[isym], svol);
	*svol /=  ((float) nsym);
	svol->update();
	EXITFUNC;
	return svol;
}

EMData* EMData::symfvol(string symString, int radius) {
	ENTERFUNC;
	int nsym = Transform::get_nsym(symString); // number of symmetries
	Transform sym;
	// set up output volume
	//cout<<"  "<<nx<<"  "<<ny<<"  "<<nz<<"  "<<is_complex()<<endl;
	EMData *svol = copy_head();
	svol->set_size(nx, ny, nz);
	svol->to_zero();
	svol->set_array_offsets(0,0,0);
	this->set_array_offsets(0,0,0);
	// actual work -- loop over symmetries and symmetrize
	vector<Transform> rm = sym.get_sym_proj(symString);
	for (int isym = 0; isym < nsym; isym++)  this -> rot_fvol(rm[isym], svol, radius);
	*svol /=  ((float) nsym);
	svol->update();
	EXITFUNC;
	return svol;
}

#define proj(ix,iy,iz)  proj[ix-1+(iy-1+(iz-1)*ny)*(size_t)nx]
#define pnewimg(ix,iy,iz)  pnewimg[ix-1+(iy-1+(iz-1)*ny)*(size_t)nx]
EMData* EMData::average_circ_sub() const
{
//  this is written as though dimensions could be different, but in fact they should be all equal nx=ny=nz,
//                                                           no check of this
	ENTERFUNC;
	const EMData* const image = this;
	EMData* newimg = copy_head();
	float *proj = image->get_data();
	float *pnewimg = newimg->get_data();
	//  Calculate average outside of a circle
	float r2 = static_cast<float>( (nx/2)*(nx/2) );
	float qs=0.0f;
	int m=0;
	int ncz = nz/2 + 1;
	int ncy = ny/2 + 1;
	int ncx = nx/2 + 1;
	for (int iz = 1; iz <= nz; iz++) {
		float yy = static_cast<float>( (iz-ncz)*(iz-ncz) );
		for (int iy = 1; iy <=ny; iy++) { float xx = yy + (iy-ncy)*(iy-ncy);
			for (int ix = 1; ix <= nx; ix++) {
				if ( xx+float((ix-ncx)*(ix-ncx)) > r2 ) {
					qs += proj(ix,iy,iz);
					m++;
				}
			}
		}
	}


	if( m > 0 ) qs /= m;
	for (int iz = 1; iz <= nz; iz++)
		for (int iy = 1; iy <= ny; iy++)
			for (int ix = 1; ix <= nx; ix++)
					pnewimg(ix,iy,iz) = proj(ix,iy,iz) - qs;
	newimg->update();
	return newimg;
	EXITFUNC;
}


//  Helper functions for method nn


void EMData::onelinenn(int j, int n, int n2, EMData* wptr, EMData* bi, const Transform& tf)
{
        //std::cout<<"   onelinenn  "<<j<<"  "<<n<<"  "<<n2<<"  "<<std::endl;
	int nnd4 = n*n/4;
	int jp = (j >= 0) ? j+1 : n+j+1;
	//for(int i = 0; i <= 2; i++){{for(int l = 0; l <= 2; l++) std::cout<<"  "<<tf[l][i];}std::cout<<std::endl;};std::cout<<std::endl;
	// loop over x
	for (int i = 0; i <= n2; i++) {
        	if (((i*i+j*j) < nnd4) && !((0 == i) && (j < 0))) {
//        if ( !((0 == i) && (j < 0))) {
			float xnew = i*tf[0][0] + j*tf[1][0];
			float ynew = i*tf[0][1] + j*tf[1][1];
			float znew = i*tf[0][2] + j*tf[1][2];
			std::complex<float> btq;
			if (xnew < 0.) {
				xnew = -xnew;
				ynew = -ynew;
				znew = -znew;
				btq = conj(bi->cmplx(i,jp));
			} else {
				btq = bi->cmplx(i,jp);
			}
			int ixn = int(xnew + 0.5 + n) - n;
			int iyn = int(ynew + 0.5 + n) - n;
			int izn = int(znew + 0.5 + n) - n;
			
			int iza, iya;
			if (izn >= 0)  iza = izn + 1;
			else	       iza = n + izn + 1;

			if (iyn >= 0) iya = iyn + 1;
			else	      iya = n + iyn + 1;

            //cout <<"  "<<jp<<"  "<<i<<"  "<<j<<"  "<< btq<<endl;

			cmplx(ixn,iya,iza) += btq;
			//std::cout<<"    "<<j<<"  "<<ixn<<"  "<<iya<<"  "<<iza<<"  "<<btq<<std::endl;
			(*wptr)(ixn,iya,iza)++;
			
			/*if ((ixn <= n2) && (iyn >= -n2) && (iyn <= n2)  && (izn >= -n2) && (izn <= n2)) {
				if (ixn >= 0) {
					int iza, iya;
					if (izn >= 0)  iza = izn + 1;
					else	       iza = n + izn + 1;

					if (iyn >= 0) iya = iyn + 1;
					else	      iya = n + iyn + 1;

					cmplx(ixn,iya,iza) += btq;
					//std::cout<<"    "<<j<<"  "<<ixn<<"  "<<iya<<"  "<<iza<<"  "<<btq<<std::endl;
					(*wptr)(ixn,iya,iza)++;
				} else {
					int izt, iyt;
					if (izn > 0) izt = n - izn + 1;
					else	     izt = -izn + 1;

					if (iyn > 0) iyt = n - iyn + 1;
					else	     iyt = -iyn + 1;

					cmplx(-ixn,iyt,izt) += conj(btq);
					//std::cout<<" *  "<<j<<"  "<<ixn<<"  "<<iyt<<"  "<<izt<<"  "<<btq<<std::endl;
					(*wptr)(-ixn,iyt,izt)++;
				}
			}*/
		}
	}
}


void EMData::onelinenn_mult(int j, int n, int n2, EMData* wptr, EMData* bi, const Transform& tf, float mult)
{
        //std::cout<<"   onelinenn  "<<j<<"  "<<n<<"  "<<n2<<"  "<<std::endl;
	int jp = (j >= 0) ? j+1 : n+j+1;
	//for(int i = 0; i <= 1; i++){for(int l = 0; l <= 2; l++){std::cout<<"  "<<tf[i][l]<<"  "<<std::endl;}}
	// loop over x
	for (int i = 0; i <= n2; i++) {
        if (((i*i+j*j) < n*n/4) && !((0 == i) && (j < 0))) {
//        if ( !((0 == i) && (j < 0))) {
			float xnew = i*tf[0][0] + j*tf[1][0];
			float ynew = i*tf[0][1] + j*tf[1][1];
			float znew = i*tf[0][2] + j*tf[1][2];
			std::complex<float> btq;
			if (xnew < 0.) {
				xnew = -xnew;
				ynew = -ynew;
				znew = -znew;
				btq = conj(bi->cmplx(i,jp));
			} else {
				btq = bi->cmplx(i,jp);
			}
			int ixn = int(xnew + 0.5 + n) - n;
			int iyn = int(ynew + 0.5 + n) - n;
			int izn = int(znew + 0.5 + n) - n;
			
			
			int iza, iya;
			if (izn >= 0)  iza = izn + 1;
			else	       iza = n + izn + 1;

			if (iyn >= 0) iya = iyn + 1;
			else	      iya = n + iyn + 1;

			cmplx(ixn,iya,iza) += btq * mult;
			(*wptr)(ixn,iya,iza)+= mult;
			
			/*if ((ixn <= n2) && (iyn >= -n2) && (iyn <= n2)  && (izn >= -n2) && (izn <= n2)) {
				if (ixn >= 0) {
					int iza, iya;
					if (izn >= 0)  iza = izn + 1;
					else	       iza = n + izn + 1;

					if (iyn >= 0) iya = iyn + 1;
					else	      iya = n + iyn + 1;

					cmplx(ixn,iya,iza) += btq*float(mult);
					(*wptr)(ixn,iya,iza)+=float(mult);
				} else {
					int izt, iyt;
					if (izn > 0) izt = n - izn + 1;
					else	     izt = -izn + 1;

					if (iyn > 0) iyt = n - iyn + 1;
					else	     iyt = -iyn + 1;

					cmplx(-ixn,iyt,izt) += conj(btq)*float(mult);
					(*wptr)(-ixn,iyt,izt)+=float(mult);
				}
			}*/
		}
	}
}

void EMData::nn(EMData* wptr, EMData* myfft, const Transform& tf, float mult)
{
	ENTERFUNC;
	int nxc = attr_dict["nxc"]; // # of complex elements along x
	// let's treat nr, bi, and local data as matrices
	vector<int> saved_offsets = get_array_offsets();
	vector<int> myfft_saved_offsets = myfft->get_array_offsets();
	set_array_offsets(0,1,1);
	myfft->set_array_offsets(0,1);
	// loop over frequencies in y
	//for(int i = 0; i <= 2; i++){{for(int l = 0; l <= 2; l++) std::cout<<"  "<<tf[l][i];}std::cout<<std::endl;};std::cout<<std::endl;
	//Dict tt = tf.get_rotation("spider");
	//std::cout << static_cast<float>(tt["phi"]) << " " << static_cast<float>(tt["theta"]) << " " << static_cast<float>(tt["psi"]) << std::endl;
	if( mult == 1 ) {
		for (int iy = -ny/2 + 1; iy <= ny/2; iy++) onelinenn(iy, ny, nxc, wptr, myfft, tf);
	} else {
		for (int iy = -ny/2 + 1; iy <= ny/2; iy++) onelinenn_mult(iy, ny, nxc, wptr, myfft, tf, mult);
	}

	set_array_offsets(saved_offsets);
	myfft->set_array_offsets(myfft_saved_offsets);
	EXITFUNC;
}


void EMData::insert_rect_slice(EMData* w, EMData* myfft, const Transform& trans, int sizeofprojection, float xratio, float yratio, float zratio, int npad, float mult)
{
	ENTERFUNC;
	vector<int> saved_offsets = get_array_offsets();
	vector<int> myfft_saved_offsets = myfft->get_array_offsets();
	set_array_offsets(0,1,1);
	myfft->set_array_offsets(0,1);
	
	// insert rectangular fft from my nn4_rect code

	Vec2f coordinate_2d_square;
	Vec3f coordinate_3dnew;
	Vec3f axis_newx;
	Vec3f axis_newy;
	Vec3f tempv;
	
       	//begin of scaling factor calculation
	//unit vector x,y of 2D fft transformed to new positon after rotation and scaling
	axis_newx[0] = xratio*0.5f*(sizeofprojection*npad)*trans[0][0];
	axis_newx[1] = yratio*0.5f*(sizeofprojection*npad)*trans[0][1];
	axis_newx[2] = zratio*0.5f*(sizeofprojection*npad)*trans[0][2];

	float ellipse_length_x = std::sqrt(axis_newx[0]*axis_newx[0]+axis_newx[1]*axis_newx[1]+axis_newx[2]*axis_newx[2]);
	
	int ellipse_length_x_int = int(ellipse_length_x);
	float ellipse_step_x = 0.5f*(sizeofprojection*npad)/float(ellipse_length_x_int);
	float xscale = ellipse_step_x;//scal increased

  	axis_newy[0] = xratio*0.5f*(sizeofprojection*npad)*trans[1][0];
  	axis_newy[1] = yratio*0.5f*(sizeofprojection*npad)*trans[1][1];
  	axis_newy[2] = zratio*0.5f*(sizeofprojection*npad)*trans[1][2];



	float ellipse_length_y = std::sqrt(axis_newy[0]*axis_newy[0]+axis_newy[1]*axis_newy[1]+axis_newy[2]*axis_newy[2]);
	int ellipse_length_y_int = int(ellipse_length_y);
	float ellipse_step_y = 0.5f*(sizeofprojection*npad)/float(ellipse_length_y_int);
	float yscale = ellipse_step_y;
	//end of scaling factor calculation
	std::complex<float> c1;
	int nxyz = sizeofprojection*npad;

	float r2=0.25f*sizeofprojection*npad*sizeofprojection*npad;
	float r2_at_point;
	
	for(int i=0;i<ellipse_length_x_int;i++) {
		for(int j=-1*ellipse_length_y_int+1; j<=ellipse_length_y_int; j++) {
        	    
			r2_at_point=i*xscale*i*xscale+j*yscale*j*yscale;
			if(r2_at_point<=r2 ) {
				
				
				coordinate_2d_square[0] = xscale*float(i);
				coordinate_2d_square[1] = yscale*float(j);
				float xnew = coordinate_2d_square[0]*trans[0][0] + coordinate_2d_square[1]*trans[1][0];
				float ynew = coordinate_2d_square[0]*trans[0][1] + coordinate_2d_square[1]*trans[1][1];
				float znew = coordinate_2d_square[0]*trans[0][2] + coordinate_2d_square[1]*trans[1][2];
				coordinate_3dnew[0] = xnew*xratio;
				coordinate_3dnew[1] = ynew*yratio;
				coordinate_3dnew[2] = znew*zratio;
				
				//bilinear interpolation
				float xp = coordinate_2d_square[0];
				float yp = ( coordinate_2d_square[1] >= 0) ? coordinate_2d_square[1]+1 : nxyz+coordinate_2d_square[1]+1;
				std::complex<float> lin_interpolated(0,0);
				int xlow=int(xp),xhigh=int(xp)+1;
				int ylow=int(yp),yhigh=int(yp)+1;
				float tx=xp-xlow,ty=yp-ylow;

				
				if(j == -1) {
					
					if(ylow<yp)
						lin_interpolated=myfft->cmplx(xlow,ylow)*(1-tx)*(1-ty) + myfft->cmplx(xlow,yhigh)*(1-tx)*ty
						+ myfft->cmplx(xhigh,ylow)*tx*(1-ty) + myfft->cmplx(xhigh,yhigh)*tx*ty;
					else 
						lin_interpolated=myfft->cmplx(xlow,ylow)*(1-tx)
 		  				+ myfft->cmplx(xhigh,ylow)*tx;
									
				} else {
						lin_interpolated=myfft->cmplx(xlow,ylow)*(1-tx)*(1-ty) + myfft->cmplx(xlow,yhigh)*(1-tx)*ty
						+ myfft->cmplx(xhigh,ylow)*tx*(1-ty)+ myfft->cmplx(xhigh,yhigh)*tx*ty;
					
				}
					
				c1 = lin_interpolated;
				
				//now nearest neighborhood interpolation
				
				std::complex<float> btq;
				if ( coordinate_3dnew[0] < 0.) {
					coordinate_3dnew[0] = -coordinate_3dnew[0];
					coordinate_3dnew[1] = -coordinate_3dnew[1];
					coordinate_3dnew[2] = -coordinate_3dnew[2];
					btq = conj(c1);
					} else {
					btq = c1;
					}
				int ixn = int(coordinate_3dnew[0] + 0.5 + nx) - nx;
				int iyn = int(coordinate_3dnew[1] + 0.5 + ny) - ny;
				int izn = int(coordinate_3dnew[2] + 0.5 + nz) - nz;

				int iza, iya;
				if (izn >= 0)  iza = izn + 1;
				else	       iza = nz + izn + 1;

				if (iyn >= 0) iya = iyn + 1;
				else	      iya = ny + iyn + 1;

				cmplx(ixn,iya,iza) += btq * mult;
				(*w)(ixn,iya,iza) += mult;
					
				}
			}
        		    
		}


	//end insert rectanular fft
		
	set_array_offsets(saved_offsets);
	myfft->set_array_offsets(myfft_saved_offsets);
	EXITFUNC;

}



void EMData::nn_SSNR(EMData* wptr, EMData* wptr2, EMData* myfft, const Transform& tf, float)
{
	ENTERFUNC;
	int nxc = attr_dict["nxc"];

	vector<int> saved_offsets = get_array_offsets();
	vector<int> myfft_saved_offsets = myfft->get_array_offsets();

	set_array_offsets(0,1,1);
       	myfft->set_array_offsets(0,1);

	int iymin = is_fftodd() ? -ny/2 : -ny/2 + 1 ;
	int iymax = ny/2;
	int izmin = is_fftodd() ? -nz/2 : -nz/2 + 1 ;
	int izmax = nz/2;

	for (int iy = iymin; iy <= iymax; iy++) {
		int jp = iy >= 0 ? iy+1 : ny+iy+1; //checked, works for both odd and even
		for (int ix = 0; ix <= nxc; ix++) {
        		if (( 4*(ix*ix+iy*iy) < ny*ny ) && !( ix == 0 && iy < 0 ) ) {
				float xnew = ix*tf[0][0] + iy*tf[1][0];
				float ynew = ix*tf[0][1] + iy*tf[1][1];
				float znew = ix*tf[0][2] + iy*tf[1][2];
				std::complex<float> btq;
				if (xnew < 0.0) {
					xnew = -xnew; // ensures xnew>=0.0
					ynew = -ynew;
					znew = -znew;
					btq = conj(myfft->cmplx(ix,jp));
				} else {
					btq = myfft->cmplx(ix,jp);
				}
				int ixn = int(xnew + 0.5 + nx) - nx; // ensures ixn >= 0
				int iyn = int(ynew + 0.5 + ny) - ny;
				int izn = int(znew + 0.5 + nz) - nz;
				if ((ixn <= nxc) && (iyn >= iymin) && (iyn <= iymax) && (izn >= izmin) && (izn <= izmax)) {
					if (ixn >= 0) {
						int iza, iya;
						if (izn >= 0)  iza = izn + 1;
						else	       iza = nz + izn + 1;

						if (iyn >= 0) iya = iyn + 1;
						else	      iya = ny + iyn + 1;

						cmplx(ixn,iya,iza) += btq;
						(*wptr)(ixn,iya,iza)++;
						(*wptr2)(ixn,iya,iza) += norm(btq);
					} else {
						int izt, iyt;
						if (izn > 0) izt = nz - izn + 1;
						else	     izt = -izn + 1;

						if (iyn > 0) iyt = ny - iyn + 1;
						else	     iyt = -iyn + 1;

						cmplx(-ixn,iyt,izt) += std::conj(btq);
						(*wptr)(-ixn,iyt,izt)++;
						(*wptr2)(-ixn,iyt,izt) += std::norm(btq);
					}
				}
			}
		}
	}
	set_array_offsets(saved_offsets);
	myfft->set_array_offsets(myfft_saved_offsets);
	EXITFUNC;
}



void EMData::symplane0(EMData* wptr) {
	ENTERFUNC;
	int nxc = attr_dict["nxc"];
	int n = nxc*2;

	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,1,1);
	for (int iza = 2; iza <= nxc; iza++) {
		for (int iya = 2; iya <= nxc; iya++) {
			cmplx(0,iya,iza) += conj(cmplx(0,n-iya+2,n-iza+2));
			(*wptr)(0,iya,iza) += (*wptr)(0,n-iya+2,n-iza+2);
			cmplx(0,n-iya+2,n-iza+2) = conj(cmplx(0,iya,iza));
			(*wptr)(0,n-iya+2,n-iza+2) = (*wptr)(0,iya,iza);
			cmplx(0,n-iya+2,iza) += conj(cmplx(0,iya,n-iza+2));
			(*wptr)(0,n-iya+2,iza) += (*wptr)(0,iya,n-iza+2);
			cmplx(0,iya,n-iza+2) = conj(cmplx(0,n-iya+2,iza));
			(*wptr)(0,iya,n-iza+2) = (*wptr)(0,n-iya+2,iza);
		}
	}
	for (int iya = 2; iya <= nxc; iya++) {
		cmplx(0,iya,1) += conj(cmplx(0,n-iya+2,1));
		(*wptr)(0,iya,1) += (*wptr)(0,n-iya+2,1);
		cmplx(0,n-iya+2,1) = conj(cmplx(0,iya,1));
		(*wptr)(0,n-iya+2,1) = (*wptr)(0,iya,1);
	}
	for (int iza = 2; iza <= nxc; iza++) {
		cmplx(0,1,iza) += conj(cmplx(0,1,n-iza+2));
		(*wptr)(0,1,iza) += (*wptr)(0,1,n-iza+2);
		cmplx(0,1,n-iza+2) = conj(cmplx(0,1,iza));
		(*wptr)(0,1,n-iza+2) = (*wptr)(0,1,iza);
	}
	EXITFUNC;
}

void EMData::symplane1(EMData* wptr, EMData* wptr2) {
	ENTERFUNC;
	int nxc = attr_dict["nxc"];
	int n = nxc*2;
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,1,1);
	for (int iza = 2; iza <= nxc; iza++) {
		for (int iya = 2; iya <= nxc; iya++) {
			cmplx(0,iya,iza) += conj(cmplx(0,n-iya+2,n-iza+2));
			(*wptr)(0,iya,iza) += (*wptr)(0,n-iya+2,n-iza+2);
			(*wptr2)(0,iya,iza) += (*wptr2)(0,n-iya+2,n-iza+2);
			cmplx(0,n-iya+2,n-iza+2) = conj(cmplx(0,iya,iza));
			(*wptr)(0,n-iya+2,n-iza+2) = (*wptr)(0,iya,iza);
			(*wptr2)(0,n-iya+2,n-iza+2) = (*wptr2)(0,iya,iza);
			cmplx(0,n-iya+2,iza) += conj(cmplx(0,iya,n-iza+2));
			(*wptr)(0,n-iya+2,iza) += (*wptr)(0,iya,n-iza+2);
			(*wptr2)(0,n-iya+2,iza) += (*wptr2)(0,iya,n-iza+2);
			cmplx(0,iya,n-iza+2) = conj(cmplx(0,n-iya+2,iza));
			(*wptr)(0,iya,n-iza+2) = (*wptr)(0,n-iya+2,iza);
			(*wptr2)(0,iya,n-iza+2) = (*wptr2)(0,n-iya+2,iza);
		}
	}
	for (int iya = 2; iya <= nxc; iya++) {
		cmplx(0,iya,1) += conj(cmplx(0,n-iya+2,1));
		(*wptr)(0,iya,1) += (*wptr)(0,n-iya+2,1);
		(*wptr2)(0,iya,1) += (*wptr2)(0,n-iya+2,1);
		cmplx(0,n-iya+2,1) = conj(cmplx(0,iya,1));
		(*wptr)(0,n-iya+2,1) = (*wptr)(0,iya,1);
		(*wptr2)(0,n-iya+2,1) = (*wptr2)(0,iya,1);
	}
	for (int iza = 2; iza <= nxc; iza++) {
		cmplx(0,1,iza) += conj(cmplx(0,1,n-iza+2));
		(*wptr)(0,1,iza) += (*wptr)(0,1,n-iza+2);
		(*wptr2)(0,1,iza) += (*wptr2)(0,1,n-iza+2);
		cmplx(0,1,n-iza+2) = conj(cmplx(0,1,iza));
		(*wptr)(0,1,n-iza+2) = (*wptr)(0,1,iza);
		(*wptr2)(0,1,n-iza+2) = (*wptr2)(0,1,iza);
	}
	EXITFUNC;
}

void EMData::symplane2(EMData* wptr, EMData* wptr2, EMData* wptr3) {
	ENTERFUNC;
	int nxc = attr_dict["nxc"];
	int n = nxc*2;
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,1,1);
	for (int iza = 2; iza <= nxc; iza++) {
		for (int iya = 2; iya <= nxc; iya++) {
			cmplx(0,iya,iza) += conj(cmplx(0,n-iya+2,n-iza+2));
			(*wptr)(0,iya,iza) += (*wptr)(0,n-iya+2,n-iza+2);
			(*wptr2)(0,iya,iza) += (*wptr2)(0,n-iya+2,n-iza+2);
			(*wptr3)(0,iya,iza) += (*wptr3)(0,n-iya+2,n-iza+2);

			cmplx(0,n-iya+2,n-iza+2) = conj(cmplx(0,iya,iza));
			(*wptr)(0,n-iya+2,n-iza+2) = (*wptr)(0,iya,iza);
			(*wptr2)(0,n-iya+2,n-iza+2) = (*wptr2)(0,iya,iza);
			(*wptr3)(0,n-iya+2,n-iza+2) = (*wptr3)(0,iya,iza);

			cmplx(0,n-iya+2,iza) += conj(cmplx(0,iya,n-iza+2));
			(*wptr)(0,n-iya+2,iza) += (*wptr)(0,iya,n-iza+2);
			(*wptr2)(0,n-iya+2,iza) += (*wptr2)(0,iya,n-iza+2);
			(*wptr3)(0,n-iya+2,iza) += (*wptr3)(0,iya,n-iza+2);

			cmplx(0,iya,n-iza+2) = conj(cmplx(0,n-iya+2,iza));
			(*wptr)(0,iya,n-iza+2) = (*wptr)(0,n-iya+2,iza);
			(*wptr2)(0,iya,n-iza+2) = (*wptr2)(0,n-iya+2,iza);
			(*wptr3)(0,iya,n-iza+2) = (*wptr3)(0,n-iya+2,iza);
		}
	}
	for (int iya = 2; iya <= nxc; iya++) {
		cmplx(0,iya,1) += conj(cmplx(0,n-iya+2,1));
		(*wptr)(0,iya,1) += (*wptr)(0,n-iya+2,1);
		(*wptr2)(0,iya,1) += (*wptr2)(0,n-iya+2,1);
		(*wptr3)(0,iya,1) += (*wptr3)(0,n-iya+2,1);

		cmplx(0,n-iya+2,1) = conj(cmplx(0,iya,1));
		(*wptr)(0,n-iya+2,1) = (*wptr)(0,iya,1);
		(*wptr2)(0,n-iya+2,1) = (*wptr2)(0,iya,1);
		(*wptr3)(0,n-iya+2,1) = (*wptr3)(0,iya,1);
	}
	for (int iza = 2; iza <= nxc; iza++) {
		cmplx(0,1,iza) += conj(cmplx(0,1,n-iza+2));
		(*wptr)(0,1,iza) += (*wptr)(0,1,n-iza+2);
		(*wptr2)(0,1,iza) += (*wptr2)(0,1,n-iza+2);
		(*wptr3)(0,1,iza) += (*wptr3)(0,1,n-iza+2);

		cmplx(0,1,n-iza+2) = conj(cmplx(0,1,iza));
		(*wptr)(0,1,n-iza+2) = (*wptr)(0,1,iza);
		(*wptr2)(0,1,n-iza+2) = (*wptr2)(0,1,iza);
		(*wptr3)(0,1,n-iza+2) = (*wptr3)(0,1,iza);
	}
	EXITFUNC;
}


void EMData::symplane0_ctf(EMData* w) {
	ENTERFUNC;
	int nxc = attr_dict["nxc"];
	int n = nxc*2;
	int nr = nxc-1;
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,1,1);

	for (int iza = 2; iza <= nxc; iza++) {
		for (int iya = 2; iya <= nxc; iya++) {
			cmplx(0,iya,iza) += conj(cmplx(0,n-iya+2,n-iza+2));
			(*w)(0,iya,iza) += (*w)(0,n-iya+2,n-iza+2);
			cmplx(0,n-iya+2,n-iza+2) = conj(cmplx(0,iya,iza));
			(*w)(0,n-iya+2,n-iza+2) = (*w)(0,iya,iza);
			cmplx(0,n-iya+2,iza) += conj(cmplx(0,iya,n-iza+2));
			(*w)(0,n-iya+2,iza) += (*w)(0,iya,n-iza+2);
			cmplx(0,iya,n-iza+2) = conj(cmplx(0,n-iya+2,iza));
			(*w)(0,iya,n-iza+2) = (*w)(0,n-iya+2,iza);
		}
	}
	for (int iya = 2; iya <= nxc; iya++) {
		cmplx(0,iya,1) += conj(cmplx(0,n-iya+2,1));
		(*w)(0,iya,1) += (*w)(0,n-iya+2,1);
		cmplx(0,n-iya+2,1) = conj(cmplx(0,iya,1));
		(*w)(0,n-iya+2,1) = (*w)(0,iya,1);
	}
	for (int iza = 2; iza <= nxc; iza++) {
		cmplx(0,1,iza) += conj(cmplx(0,1,n-iza+2));
		(*w)(0,1,iza) += (*w)(0,1,n-iza+2);
		cmplx(0,1,n-iza+2) = conj(cmplx(0,1,iza));
		(*w)(0,1,n-iza+2) = (*w)(0,1,iza);
	}
	set_array_offsets(saved_offsets);
	EXITFUNC;
}


void EMData::symplane0_odd(EMData* w) {
	ENTERFUNC;
	ny=get_ysize();
	nz=get_zsize();
	if( ny != nz ) {
		LOGERR("symplane0_odd requires ny = nz.");
		throw InvalidValueException(ny, "symplane0_odd requires ny = nz.");
	}
	int nyh = ny/2;
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
	vector<int> saved_offsets_w = w->get_array_offsets();
	w->set_array_offsets(0,0,0);

	int iyb,iza,izc,iyc;
	for (int iz = -nyh; iz <= nyh; iz++) {
		if( iz < 0 )  {
			iyb = 0;
			iza = iz + nz;
			izc = -iz;
		} else if( iz > 0 )  {
			iyb = 1;
			iza = iz;
			izc = nz - iz;
		} else  {
			iyb = 1;
			iza = iz;
			izc = iz;
		}
		for (int iy = iyb; iy <= nyh; iy++) {
			if( iy == 0 ) iyc = iy;
			else  iyc = ny - iy;
			cmplx(0,iy,iza) += conj(cmplx(0,iyc,izc));
			(*w)(0,iy,iza) += (*w)(0,iyc,izc);
			cmplx(0,iyc,izc) = conj(cmplx(0,iy,iza));
			(*w)(0,iyc,izc) = (*w)(0,iy,iza);
		}
	}
	set_array_offsets(saved_offsets);
	w->set_array_offsets(saved_offsets_w);
	EXITFUNC;
}

/*
void EMData::symplane0_odd(EMData* w) {
	ENTERFUNC;
	ny=get_ysize();
	nz=get_zsize();
	if( ny != nz ) {
		LOGERR("symplane0_odd requires ny = nz.");
		throw InvalidValueException(ny, "symplane0_odd requires ny = nz.");
	}
	int nyh = ny/2;
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
	vector<int> saved_offsets_w = w->get_array_offsets();
	w->set_array_offsets(0,0,0);
cout <<"  symplane0_odd   "<<ny<<"   "<<nyh<<endl;

	for (int iza = 1; iza <= nyh; iza++) {
		for (int iya = 1; iya < ny; iya++) {
			cmplx(0,iya,iza) += conj(cmplx(0,ny-iya,nz-iza));
			(*w)(0,iya,iza) += (*w)(0,ny-iya,nz-iza);
			cmplx(0,ny-iya,nz-iza) = conj(cmplx(0,iya,iza));
			(*w)(0,ny-iya,nz-iza) = (*w)(0,iya,iza);
		}
	}
	for (int iya = 1; iya <= nyh; iya++) {
		cmplx(0,iya,0) += conj(cmplx(0,ny-iya,0));
		(*w)(0,iya,0) += (*w)(0,ny-iya,0);
		cmplx(0,ny-iya,1) = conj(cmplx(0,iya,0));
		(*w)(0,ny-iya,1) = (*w)(0,iya,0);
		cmplx(0,0,iya) += conj(cmplx(0,0,ny-iya));
		(*w)(0,0,iya) += (*w)(0,0,ny-iya);
		cmplx(0,0,nz-iya) = conj(cmplx(0,0,iya));
		(*w)(0,0,nz-iya) = (*w)(0,0,iya);
	}
	set_array_offsets(saved_offsets);
	w->set_array_offsets(saved_offsets_w);
	EXITFUNC;
}
*/

void EMData::symplane0_rect(EMData* w) {
	ENTERFUNC;
	nx=get_xsize();
	ny=get_ysize();
	nz=get_zsize();
	int nzc=nz/2;
	int nyc=ny/2;
	
	// let's treat the local data as a matrix
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,1,1);
	for (int iza = 2; iza <= nzc; iza++) {
		for (int iya = 2; iya <= nyc; iya++) {
			cmplx(0,iya,iza) += conj(cmplx(0,ny-iya+2,nz-iza+2));
			(*w)(0,iya,iza) += (*w)(0,ny-iya+2,nz-iza+2);
			cmplx(0,ny-iya+2,nz-iza+2) = conj(cmplx(0,iya,iza));
			(*w)(0,ny-iya+2,nz-iza+2) = (*w)(0,iya,iza);
			cmplx(0,ny-iya+2,iza) += conj(cmplx(0,iya,nz-iza+2));
			(*w)(0,ny-iya+2,iza) += (*w)(0,iya,nz-iza+2);
			cmplx(0,iya,nz-iza+2) = conj(cmplx(0,ny-iya+2,iza));
			(*w)(0,iya,nz-iza+2) = (*w)(0,ny-iya+2,iza);
		}
	}
	for (int iya = 2; iya <= nyc; iya++) {
		cmplx(0,iya,1) += conj(cmplx(0,ny-iya+2,1));
		(*w)(0,iya,1) += (*w)(0,ny-iya+2,1);
		cmplx(0,ny-iya+2,1) = conj(cmplx(0,iya,1));
		(*w)(0,ny-iya+2,1) = (*w)(0,iya,1);
	}
	for (int iza = 2; iza <= nzc; iza++) {
		cmplx(0,1,iza) += conj(cmplx(0,1,nz-iza+2));
		(*w)(0,1,iza) += (*w)(0,1,nz-iza+2);
		cmplx(0,1,nz-iza+2) = conj(cmplx(0,1,iza));
		(*w)(0,1,nz-iza+2) = (*w)(0,1,iza);
	}
	EXITFUNC;
}


class ctf_store
{
public:

    static void init( int winsize, const Ctf* ctf ) {
		Dict params = ctf->to_dict();

		m_winsize = winsize;

		m_voltage = params["voltage"];
		m_pixel   = params["apix"];
		m_cs      = params["cs"];
		m_ampcont = params["ampcont"];
		m_bfactor = params["bfactor"];
		m_defocus = params["defocus"];
		m_dza     = params["dfdiff"];
		m_azz     = params["dfang"];
		m_winsize2= m_winsize*m_winsize;
		m_vecsize = m_winsize2/4;
    }

    static float get_ctf( int r2 ,int i, int j) {
		float  ak = std::sqrt( r2/float(m_winsize2) )/m_pixel;
		if(m_dza == 0.0f)  return Util::tf( m_defocus, ak, m_voltage, m_cs, m_ampcont, m_bfactor, 1);
		else {
			float az = atan2(float(j), float(i));
			float dzz = m_defocus - m_dza/2.0f*sin(2*(az+m_azz*M_PI/180.0f));
			return Util::tf( dzz, ak, m_voltage, m_cs, m_ampcont, m_bfactor, 1);
		}
	}
	
private:

	static int m_winsize, m_winsize2, m_vecsize;
	static float m_cs;
	static float m_voltage;
	static float m_pixel;
	static float m_ampcont;
	static float m_bfactor;
	static float m_defocus;
	static float m_dza;
	static float m_azz;
};

int ctf_store::m_winsize, ctf_store::m_winsize2, ctf_store::m_vecsize;

float ctf_store::m_cs, ctf_store::m_voltage, ctf_store::m_pixel;
float ctf_store::m_ampcont, ctf_store::m_bfactor;
float ctf_store::m_defocus, ctf_store::m_dza, ctf_store::m_azz;

//  This class was for some reason introduced for rectangular reconstruction, which in any case is not used anymore PAP 06/23/2015
class ctf_store_new
{
public:

	static void init( int winsize, const Ctf* ctf ) {
		Dict params = ctf->to_dict();

		m_winsize = winsize;

		m_voltage = params["voltage"];
		m_pixel   = params["apix"];
		m_cs      = params["cs"];
		m_ampcont = params["ampcont"];
		m_bfactor = params["bfactor"];
		m_defocus = params["defocus"];
		m_dza     = params["dfdiff"];
		m_azz     = params["dfang"];
		m_winsize2= m_winsize*m_winsize;
		m_vecsize = m_winsize2/4;
    }

    static float get_ctf( float r2 ) {  //  HAS TO BE CORRECTED AS astigmatism m_dza and m_azz is not used!!  PAP 04/27/2013
		float ak = std::sqrt( r2/float(m_winsize2) )/m_pixel;
		return Util::tf( m_defocus, ak, m_voltage, m_cs, m_ampcont, m_bfactor, 1);
    }

private:

	static int m_winsize, m_winsize2, m_vecsize;
	static float m_cs;
	static float m_voltage;
	static float m_pixel;
	static float m_ampcont;
	static float m_bfactor;
	static float m_defocus;
	static float m_dza;
	static float m_azz;
};


int ctf_store_new::m_winsize, ctf_store_new::m_winsize2, ctf_store_new::m_vecsize;

float ctf_store_new::m_cs, ctf_store_new::m_voltage, ctf_store_new::m_pixel;
float ctf_store_new::m_ampcont, ctf_store_new::m_bfactor;
float ctf_store_new::m_defocus, ctf_store_new::m_dza, ctf_store_new::m_azz;



//  Helper functions for method nn4_ctf
void EMData::onelinenn_ctf(int j, int n, int n2, EMData* w, EMData* bi, const Transform& tf, float mult) {
//std::cout<<"   onelinenn_ctf  "<<j<<"  "<<n<<"  "<<n2<<"  "<<std::endl;

	//int remove = bi->get_attr_default( "remove", 0 );

	int jp = (j >= 0) ? j+1 : n+j+1;
	// loop over x
	for (int i = 0; i <= n2; i++) {
		int r2 = i*i+j*j;
		if ( (r2<n*n/4) && !((0==i) && (j<0)) ) {
			float ctf = ctf_store::get_ctf( r2, i, j ); //This is in 2D projection plane
			float xnew = i*tf[0][0] + j*tf[1][0];
			float ynew = i*tf[0][1] + j*tf[1][1];
			float znew = i*tf[0][2] + j*tf[1][2];
			std::complex<float> btq;
			if (xnew < 0.) {
				xnew = -xnew;
				ynew = -ynew;
				znew = -znew;
				btq = conj(bi->cmplx(i,jp));
			} else  btq = bi->cmplx(i,jp);
			int ixn = int(xnew + 0.5 + n) - n;
			int iyn = int(ynew + 0.5 + n) - n;
			int izn = int(znew + 0.5 + n) - n;

			int iza, iya;
			if (izn >= 0)  iza = izn + 1;
			else           iza = n + izn + 1;

			if (iyn >= 0)  iya = iyn + 1;
			else           iya = n + iyn + 1;

			//if(remove > 0 ) {
			//	cmplx(ixn,iya,iza) -= btq*ctf*mult;
			//	(*w)(ixn,iya,iza)  -= ctf*ctf*mult;
			//} else {
				cmplx(ixn,iya,iza) += btq*ctf*mult;
				(*w)(ixn,iya,iza)  += ctf*ctf*mult;
 			//}

		}
	}
}

//  Helper functions for method nn4_ctfw
void EMData::onelinenn_ctfw(int j, int n, int n2,
		          EMData* w, EMData* bi, EMData* c2, vector<float> bckgnoise, const Transform& tf, float weight) {
//std::cout<<"   onelinenn_ctf  "<<j<<"  "<<n<<"  "<<n2<<"  "<<std::endl;

	int nnd4 = n*n/4;
	int jp = (j >= 0) ? j+1 : n+j+1;
	// loop over x
	for (int i = 0; i <= n2; i++) {
		int r2 = i*i + j*j;
		if ( (r2 < nnd4) && !((0 == i) && (j < 0)) ) {
			float xnew = i*tf[0][0] + j*tf[1][0];
			float ynew = i*tf[0][1] + j*tf[1][1];
			float znew = i*tf[0][2] + j*tf[1][2];
			std::complex<float> btq;
			if (xnew < 0.) {
				xnew = -xnew;
				ynew = -ynew;
				znew = -znew;
				btq = conj(bi->cmplx(i,jp));
			} else  btq = bi->cmplx(i,jp);
			
			float c2val = (*c2)(i,jp);
			
			int ixn = int(xnew + 0.5 + n) - n;
			int iyn = int(ynew + 0.5 + n) - n;
			int izn = int(znew + 0.5 + n) - n;

			int iza, iya;
			if (izn >= 0)  iza = izn + 1;
			else           iza = n + izn + 1;

			if (iyn >= 0) iya = iyn + 1;
			else          iya = n + iyn + 1;

			// linear interpolation of 1D bckgnoise
			float rr = std::sqrt(float(r2));
			int   ir = int(rr);
			float df = rr - float(ir);
			float mult = (1.0f - df)*bckgnoise[ir] + df*bckgnoise[ir+1];
			//cout <<"  "<<jp<<"  "<<i<<"  "<<j<<"  "<<rr<<"  "<<ir<<"  "<<mult<<"  "<<1.0f/mult<<"  "<<btq<<"  "<<weight<<endl;
			// cmplx(ixn, iya, iza) += btq*ctf*mult*weight;
			// (*w)(ixn, iya, iza)  += ctf*ctf*mult*weight;
			cmplx(ixn, iya, iza) +=  btq * mult * weight;
			(*w)(ixn, iya, iza)  += c2val * mult * weight;

		}
	}
}

//  Helper functions for method nn4_ctfw with tri-linear interpolation
void EMData::onelinetr_ctfw(int j, int bign, int n, int n2, int npad,
		          EMData* w, EMData* bi, EMData* c2, vector<float> bckgnoise, const Transform& tf, float weight) {
//std::cout<<"   onelinetr_ctfw  "<<j<<"  "<<n<<"   "<<bign<<"  "<<n<<"  "<<n2<<"  "<<npad<<std::endl;

	int nnd4 = n*n/4;
	int jp = (j >= 0) ? j+1 : n+j+1;
	// loop over x
	for (int i = 0; i <= n2; i++) {
		int r2 = i*i + j*j;
		if ( (r2 < nnd4) && !((0 == i) && (j < 0)) ) {
			float xnew = (i*tf[0][0] + j*tf[1][0])*npad;
			float ynew = (i*tf[0][1] + j*tf[1][1])*npad;
			float znew = (i*tf[0][2] + j*tf[1][2])*npad;
			std::complex<float> btq;
			if (xnew < 0.) {
				xnew = -xnew;
				ynew = -ynew;
				znew = -znew;
				btq = conj(bi->cmplx(i,jp));
			} else  btq = bi->cmplx(i,jp);

			// linear interpolation of 1D bckgnoise
			float rr = std::sqrt(float(r2));
			int   ir = int(rr);
			float df = rr - float(ir);
			float mult = (1.0f - df)*bckgnoise[ir] + df*bckgnoise[ir+1];

			float c2val = (*c2)(i,jp);
			//std::complex<float> numerator = btq * mult * weight;
			std::complex<float> numerator = btq * weight;//  Data is always premultiplied
			float denominator = c2val * mult * weight;


			int ixn = int(xnew + bign); // Here bign has no particular meaning, it only matters it is much larger than -xnew
			int iyn = int(ynew + bign);
			int izn = int(znew + bign);

			float dx = xnew + bign - ixn;
			float dy = ynew + bign - iyn;
			float dz = znew + bign - izn;
			float qdx = 1.0f - dx;
			float qdy = 1.0f - dy;
			float qdz = 1.0f - dz;

			float qq000 = qdx * qdy * qdz;
			float qq010 = qdx *  dy * qdz;
			float qq100 =  dx * qdy * qdz;
			float qq110 =  dx *  dy * qdz;
			float qq001 = qdx * qdy *  dz;
			float qq011 = qdx *  dy *  dz;
			float qq101 =  dx * qdy *  dz;
			float qq111 =  dx *  dy *  dz;

			ixn -= bign;
			iyn -= bign;
			izn -= bign;



/*
			int iza, iya;
			if (izn >= 0)  iza = izn + 1;
			else           iza = n + izn + 1;

			if (iyn >= 0) iya = iyn + 1;
			else          iya = n + iyn + 1;

            cout <<"  "<<jp<<"  "<<i<<"  "<<j<<"  "<<rr<<"  "<<ir<<"  "<<c2val<<"  "<<mult<<"  "<<1.0f/mult<<"  "<<btq<<"  "<<weight<<endl;
			// cmplx(ixn, iya, iza) += btq*ctf*mult*weight;
			// (*w)(ixn, iya, iza)  += ctf*ctf*mult*weight;
			cmplx(ixn, iya, iza) +=  btq * mult * weight;
			(*w)(ixn, iya, iza)  += c2val * mult * weight;
*/
//cout <<"  "<<jp<<"  "<<i<<"  "<<j<<"  "<<ixn<<"  "<<iyn<<"  "<<izn<<"  "<<n<<endl;

			int iza, iya;
			if (iyn >= 0) iya = iyn + 1;
			else          iya = bign + iyn + 1;  // Here bign has correct meaning

			if (izn >= 0)  iza = izn + 1;
			else           iza = bign + izn + 1;

			int ix1 = ixn + 1;

			int iy1 = iya + 1;
			if (iy1 > bign) iy1 -= bign;

			int iz1 = iza + 1;
			if (iz1 > bign) iz1 -= bign;


/*
			int iy1 = iyn +1;
			if (iy1 >= 0) iy1 = iy1 + 1;
			else          iy1 = bign + iy1 + 1;

			int iz1 = izn +1;
			if (iz1 >= 0) iz1 = iz1 + 1;
			else          iz1 = bign + iz1 + 1;
*/
            //cout <<" onetrl "<<jp<<"  "<<i<<"  "<<j<<"  "<<ixn<<"  "<<iya<<"  "<<iza<<"  "<<endl;
			// cmplx(ixn, iya, iza) += btq*ctf*mult*weight;
			// (*w)(ixn, iya, iza)  += ctf*ctf*mult*weight;
/*
	if(ixn<0 or ixn >=n2*npad)  cout<<"  error   ixn  "<<ixn<<endl;
	if(iya<1 or iya >bign)  cout<<"  error   iya  "<<iya<<endl;
	if(iza<1 or iza >bign)  cout<<"  error   iza  "<<iza<<endl;
	if(iy1<1 or iy1 >bign)  cout<<"  error   iy1  "<<iy1<<endl;
	if(iz1<1 or iz1 >bign)  cout<<"  error   iz1  "<<iz1<<endl;
*/
			// numerator
			cmplx(ixn, iya, iza) += qq000 * numerator;
			cmplx(ixn, iy1, iza) += qq010 * numerator;
			cmplx(ix1, iya, iza) += qq100 * numerator;
			cmplx(ix1, iy1, iza) += qq110 * numerator;
			cmplx(ixn, iya, iz1) += qq001 * numerator;
			cmplx(ixn, iy1, iz1) += qq011 * numerator;
			cmplx(ix1, iya, iz1) += qq101 * numerator;
			cmplx(ix1, iy1, iz1) += qq111 * numerator;
			// denominator
			(*w)(ixn, iya, iza) += qq000 * denominator;
			(*w)(ixn, iy1, iza) += qq010 * denominator;
			(*w)(ix1, iy1, iza) += qq110 * denominator;
			(*w)(ix1, iya, iza) += qq100 * denominator;
			(*w)(ixn, iya, iz1) += qq001 * denominator;
			(*w)(ixn, iy1, iz1) += qq011 * denominator;
			(*w)(ix1, iya, iz1) += qq101 * denominator;
			(*w)(ix1, iy1, iz1) += qq111 * denominator;


		}
	}
}

void EMData::onelinenn_ctf_applied(int j, int n, int n2,
		          EMData* w, EMData* bi, const Transform& tf, float mult) {//std::cout<<"   onelinenn_ctf  "<<j<<"  "<<n<<"  "<<n2<<"  "<<std::endl;

	int remove = bi->get_attr_default( "remove", 0 );

	int jp = (j >= 0) ? j+1 : n+j+1;
	// loop over x
	for (int i = 0; i <= n2; i++) {
	        int r2 = i*i + j*j;
		if ( (r2< n*n/4) && !((0==i) && (j< 0)) ) {
			float  ctf = ctf_store::get_ctf(r2, i, j);

			 //	   if ( !((0 == i) && (j < 0))) {
			float xnew = i*tf[0][0] + j*tf[1][0];
			float ynew = i*tf[0][1] + j*tf[1][1];
			float znew = i*tf[0][2] + j*tf[1][2];
			std::complex<float> btq;
			if (xnew < 0.) {
				xnew = -xnew;
				ynew = -ynew;
				znew = -znew;
				btq = conj(bi->cmplx(i,jp));
			} else  btq = bi->cmplx(i,jp);
			int ixn = int(xnew + 0.5 + n) - n;
			int iyn = int(ynew + 0.5 + n) - n;
			int izn = int(znew + 0.5 + n) - n;

			int iza, iya;
			if (izn >= 0)  iza = izn + 1;
			else           iza = n + izn + 1;

			if (iyn >= 0) iya = iyn + 1;
			else          iya = n + iyn + 1;

			if( remove > 0 ) {
				cmplx(ixn,iya,iza) -= btq*mult;
				(*w)(ixn,iya,iza) -= mult*ctf*ctf;
			} else {
				cmplx(ixn,iya,iza) += btq*mult;
				(*w)(ixn,iya,iza) += mult*ctf*ctf;
			}

		}
	}
}

void EMData::onelinenn_ctf_exists(int j, int n, int n2,
		          EMData* w, EMData* bi, EMData* c2, const Transform& tf, float weight) {
	//std::cout<<"   onelinenn_ctf  "<<j<<"  "<<n<<"  "<<n2<<"  "<<std::endl;
    //int remove = bi->get_attr_default( "remove", 0 );

	int nnd4 = n*n/4;
	int jp = (j >= 0) ? j+1 : n+j+1;
	// loop over x
	for (int i = 0; i <= n2; i++) {
	    int r2 = i*i + j*j;
		if ( (r2 < n*n/4) && !((0 == i) && (j < 0)) ) {
			float ctf = ctf_store::get_ctf( r2, i, j ); //This is in 2D projection plane
			//	   if ( !((0 == i) && (j < 0))) {
			float xnew = i*tf[0][0] + j*tf[1][0];
			float ynew = i*tf[0][1] + j*tf[1][1];
			float znew = i*tf[0][2] + j*tf[1][2];
			std::complex<float> btq;
			if (xnew < 0.) {
				xnew = -xnew;
				ynew = -ynew;
				znew = -znew;
				btq = conj(bi->cmplx(i,jp));
			} else  btq = bi->cmplx(i,jp);
			
			float c2val = (*c2)(i,jp);
			
			int ixn = int(xnew + 0.5 + n) - n;
			int iyn = int(ynew + 0.5 + n) - n;
			int izn = int(znew + 0.5 + n) - n;
			
			int iza, iya;
			if (izn >= 0)  iza = izn + 1;
			else           iza = n + izn + 1;

			if (iyn >= 0) iya = iyn + 1;
			else          iya = n + iyn + 1;

			//if( remove > 0 ) {
			//	cmplx(ixn,iya,iza) -= btq*weight;
			//	(*w)(ixn,iya,iza) -= c2val*weight;
			//} else {
				cmplx(ixn, iya, iza) += btq * weight;
				(*w)(ixn, iya, iza) += c2val * weight;
			//}

		}
	}
}

void EMData::nn_ctf(EMData* w, EMData* myfft, const Transform& tf, float mult) {
	ENTERFUNC;
	int nxc = attr_dict["nxc"]; // # of complex elements along x
	// let's treat nr, bi, and local data as matrices
	vector<int> saved_offsets = get_array_offsets();
	vector<int> myfft_saved_offsets = myfft->get_array_offsets();
	set_array_offsets(0,1,1);
	myfft->set_array_offsets(0,1);

    Ctf* ctf = myfft->get_attr("ctf");
    ctf_store::init( ny, ctf );
    if(ctf) {delete ctf; ctf=0;}

	// loop over frequencies in y
	for (int iy = -ny/2 + 1; iy <= ny/2; iy++) onelinenn_ctf(iy, ny, nxc, w, myfft, tf, mult);
	set_array_offsets(saved_offsets);
	myfft->set_array_offsets(myfft_saved_offsets);
	EXITFUNC;
}

void EMData::nn_ctfw(EMData* w, EMData* myfft, EMData* ctf2d2, int npad, vector<float> bckgnoise, const Transform& tf, float weight ) {
	ENTERFUNC;
	//int nxc = attr_dict["nxc"]; // # of complex elements along x
	// let's treat nr, bi, and local data as matrices
	vector<int> saved_offsets = get_array_offsets();
	vector<int> myfft_saved_offsets = myfft->get_array_offsets();
	set_array_offsets(0,1,1);
	myfft->set_array_offsets(0,1);
    
	vector<int> ctf2d2_saved_offsets = ctf2d2->get_array_offsets();
	ctf2d2->set_array_offsets(0,1);
	int mynx = myfft->get_xsize();
	mynx /= 2;
	int myny = myfft->get_ysize();
	//cout<<"  dimensions in nn_ctfw  "<<nx<<"   "<<ny<<"   "<<mynx<<"   "<<myny<<endl;
	// loop over frequencies in y
	for (int iy = -myny/2 + 1; iy <= myny/2; iy++) onelinetr_ctfw(iy, ny, myny, mynx, npad, w, myfft, ctf2d2, bckgnoise, tf, weight);
	set_array_offsets(saved_offsets);
	myfft->set_array_offsets(myfft_saved_offsets);
	ctf2d2->set_array_offsets(ctf2d2_saved_offsets);
	EXITFUNC;
}

void EMData::nn_ctf_applied(EMData* w, EMData* myfft, const Transform& tf, float mult) {
	ENTERFUNC;
	int nxc = attr_dict["nxc"]; // # of complex elements along x
	// let's treat nr, bi, and local data as matrices
	vector<int> saved_offsets = get_array_offsets();
	vector<int> myfft_saved_offsets = myfft->get_array_offsets();
	set_array_offsets(0,1,1);
	myfft->set_array_offsets(0,1);

	Ctf* ctf = myfft->get_attr( "ctf" );
	ctf_store::init( ny, ctf );
	if(ctf) {delete ctf; ctf=0;}

	// loop over frequencies in y
	for (int iy = -ny/2 + 1; iy <= ny/2; iy++) onelinenn_ctf_applied(iy, ny, nxc, w, myfft, tf, mult);
	set_array_offsets(saved_offsets);
	myfft->set_array_offsets(myfft_saved_offsets);
	EXITFUNC;
}

void EMData::nn_ctf_exists(EMData* w, EMData* myfft, EMData* ctf2d2, const Transform& tf, float weight) {
	ENTERFUNC;
	int nxc = attr_dict["nxc"]; // # of complex elements along x
	// let's treat nr, bi, and local data as matrices
	vector<int> saved_offsets = get_array_offsets();
	vector<int> myfft_saved_offsets = myfft->get_array_offsets();
	set_array_offsets(0,1,1);
	myfft->set_array_offsets(0,1);

	vector<int> ctf2d2_saved_offsets = ctf2d2->get_array_offsets();
	ctf2d2->set_array_offsets(0,1);

	// loop over frequencies in y
	for (int iy = -ny/2 + 1; iy <= ny/2; iy++) onelinenn_ctf_exists(iy, ny, nxc, w, myfft, ctf2d2, tf, weight);
	set_array_offsets(saved_offsets);
	myfft->set_array_offsets(myfft_saved_offsets);
	ctf2d2->set_array_offsets(ctf2d2_saved_offsets);
	EXITFUNC;
}

void EMData::insert_rect_slice_ctf(EMData* w, EMData* myfft, const Transform& trans, int sizeofprojection, float xratio, float yratio, float zratio, int npad, float mult)
{
	ENTERFUNC;
	vector<int> saved_offsets = get_array_offsets();
	vector<int> myfft_saved_offsets = myfft->get_array_offsets();
	set_array_offsets(0,1,1);
	myfft->set_array_offsets(0,1);
	
	// insert rectangular fft from my nn4_rect code

	Vec2f coordinate_2d_square;
	Vec3f coordinate_3dnew;
	Vec3f axis_newx;
	Vec3f axis_newy;
	Vec3f tempv;
	
       	//begin of scaling factor calculation
	//unit vector x,y of 2D fft transformed to new positon after rotation and scaling
	axis_newx[0] = xratio*0.5f*(sizeofprojection*npad)*trans[0][0];
	axis_newx[1] = yratio*0.5f*(sizeofprojection*npad)*trans[0][1];
	axis_newx[2] = zratio*0.5f*(sizeofprojection*npad)*trans[0][2];

	float ellipse_length_x = std::sqrt(axis_newx[0]*axis_newx[0]+axis_newx[1]*axis_newx[1]+axis_newx[2]*axis_newx[2]);
	
	int ellipse_length_x_int = int(ellipse_length_x);
	float ellipse_step_x = 0.5f*(sizeofprojection*npad)/float(ellipse_length_x_int);
	float xscale = ellipse_step_x;//scal increased

  	axis_newy[0] = xratio*0.5f*(sizeofprojection*npad)*trans[1][0];
  	axis_newy[1] = yratio*0.5f*(sizeofprojection*npad)*trans[1][1];
  	axis_newy[2] = zratio*0.5f*(sizeofprojection*npad)*trans[1][2];



	float ellipse_length_y = std::sqrt(axis_newy[0]*axis_newy[0]+axis_newy[1]*axis_newy[1]+axis_newy[2]*axis_newy[2]);
	int ellipse_length_y_int = int(ellipse_length_y);
	float ellipse_step_y = 0.5f*(sizeofprojection*npad)/float(ellipse_length_y_int);
	float yscale = ellipse_step_y;
	//end of scaling factor calculation
	std::complex<float> c1;
	int nxyz = sizeofprojection*npad;
	Ctf* ctf = myfft->get_attr( "ctf" );
	ctf_store_new::init( nxyz, ctf );
	if(ctf) {delete ctf; ctf=0;}
	int remove = myfft->get_attr_default( "remove", 0 );

	float r2=0.25f*sizeofprojection*npad*sizeofprojection*npad;
	float r2_at_point;
	
	for(int i=0;i<ellipse_length_x_int;i++) {
		for(int j=-1*ellipse_length_y_int+1; j<=ellipse_length_y_int; j++) {
        	    
			r2_at_point=i*xscale*i*xscale+j*yscale*j*yscale;
			if(r2_at_point<=r2 && ! ((0==i) && (j<0))) {
				
				float ctf_value = ctf_store_new::get_ctf( r2_at_point );
				coordinate_2d_square[0] = xscale*float(i);
				coordinate_2d_square[1] = yscale*float(j);
				float xnew = coordinate_2d_square[0]*trans[0][0] + coordinate_2d_square[1]*trans[1][0];
				float ynew = coordinate_2d_square[0]*trans[0][1] + coordinate_2d_square[1]*trans[1][1];
				float znew = coordinate_2d_square[0]*trans[0][2] + coordinate_2d_square[1]*trans[1][2];
				coordinate_3dnew[0] = xnew*xratio;
				coordinate_3dnew[1] = ynew*yratio;
				coordinate_3dnew[2] = znew*zratio;
				
				//bilinear interpolation
				float xp = coordinate_2d_square[0];
				float yp = ( coordinate_2d_square[1] >= 0) ? coordinate_2d_square[1]+1 : nxyz+coordinate_2d_square[1]+1;
				std::complex<float> lin_interpolated(0,0);
				int xlow=int(xp),xhigh=int(xp)+1;
				int ylow=int(yp),yhigh=int(yp)+1;
				float tx=xp-xlow,ty=yp-ylow;

				
				if(j == -1) {
					
					if(ylow<yp)
						lin_interpolated=myfft->cmplx(xlow,ylow)*(1-tx)*(1-ty) + myfft->cmplx(xlow,yhigh)*(1-tx)*ty
						+ myfft->cmplx(xhigh,ylow)*tx*(1-ty) + myfft->cmplx(xhigh,yhigh)*tx*ty;
					else 
						lin_interpolated=myfft->cmplx(xlow,ylow)*(1-tx)
 		  				+ myfft->cmplx(xhigh,ylow)*tx;
									
				} else {
						lin_interpolated=myfft->cmplx(xlow,ylow)*(1-tx)*(1-ty) + myfft->cmplx(xlow,yhigh)*(1-tx)*ty
						+ myfft->cmplx(xhigh,ylow)*tx*(1-ty)+ myfft->cmplx(xhigh,yhigh)*tx*ty;
					
				}
					
				c1 = lin_interpolated;
				
				//now nearest neighborhood interpolation
				
				std::complex<float> btq;
				if ( coordinate_3dnew[0] < 0.) {
					coordinate_3dnew[0] = -coordinate_3dnew[0];
					coordinate_3dnew[1] = -coordinate_3dnew[1];
					coordinate_3dnew[2] = -coordinate_3dnew[2];
					btq = conj(c1);
					} else {
					btq = c1;
					}
				int ixn = int(coordinate_3dnew[0] + 0.5 + nx) - nx;
				int iyn = int(coordinate_3dnew[1] + 0.5 + ny) - ny;
				int izn = int(coordinate_3dnew[2] + 0.5 + nz) - nz;

				int iza, iya;
				if (izn >= 0)  iza = izn + 1;
				else	       iza = nz + izn + 1;

				if (iyn >= 0) iya = iyn + 1;
				else	      iya = ny + iyn + 1;

				if(remove > 0 ) {
					cmplx(ixn,iya,iza) -= btq*ctf_value * mult;
					(*w)(ixn,iya,iza) -= ctf_value*ctf_value*mult;
					} else {
					cmplx(ixn,iya,iza) += btq*ctf_value * mult;
					(*w)(ixn,iya,iza) += ctf_value*ctf_value*mult;
					}
					
				}
			}
        		    
		}


	//end insert rectanular fft
		
	set_array_offsets(saved_offsets);
	myfft->set_array_offsets(myfft_saved_offsets);
	EXITFUNC;

}


void EMData::insert_rect_slice_ctf_applied(EMData* w, EMData* myfft,const Transform& trans,int sizeofprojection,float xratio,float yratio, float zratio, int npad,float mult)
{
	ENTERFUNC;
	vector<int> saved_offsets = get_array_offsets();
	vector<int> myfft_saved_offsets = myfft->get_array_offsets();
	set_array_offsets(0,1,1);
	myfft->set_array_offsets(0,1);
	
	// insert rectangular fft from my nn4_rect code

	Vec2f coordinate_2d_square;
	Vec3f coordinate_3dnew;
	Vec3f axis_newx;
	Vec3f axis_newy;
	Vec3f tempv;
	
       	//begin of scaling factor calculation
	//unit vector x,y of 2D fft transformed to new positon after rotation and scaling
	axis_newx[0] = xratio*0.5f*(sizeofprojection*npad)*trans[0][0];
	axis_newx[1] = yratio*0.5f*(sizeofprojection*npad)*trans[0][1];
	axis_newx[2] = zratio*0.5f*(sizeofprojection*npad)*trans[0][2];

	float ellipse_length_x = std::sqrt(axis_newx[0]*axis_newx[0]+axis_newx[1]*axis_newx[1]+axis_newx[2]*axis_newx[2]);
	
	int ellipse_length_x_int = int(ellipse_length_x);
	float ellipse_step_x = 0.5f*(sizeofprojection*npad)/float(ellipse_length_x_int);
	float xscale = ellipse_step_x;//scal increased

  	axis_newy[0] = xratio*0.5f*(sizeofprojection*npad)*trans[1][0];
  	axis_newy[1] = yratio*0.5f*(sizeofprojection*npad)*trans[1][1];
  	axis_newy[2] = zratio*0.5f*(sizeofprojection*npad)*trans[1][2];



	float ellipse_length_y = std::sqrt(axis_newy[0]*axis_newy[0]+axis_newy[1]*axis_newy[1]+axis_newy[2]*axis_newy[2]);
	int ellipse_length_y_int = int(ellipse_length_y);
	float ellipse_step_y = 0.5f*(sizeofprojection*npad)/float(ellipse_length_y_int);
	float yscale = ellipse_step_y;
	//end of scaling factor calculation
	std::complex<float> c1;
	int nxyz = sizeofprojection*npad;
	Ctf* ctf = myfft->get_attr( "ctf" );
	ctf_store_new::init( nxyz, ctf );
	if(ctf) {delete ctf; ctf=0;}
	int remove = myfft->get_attr_default( "remove", 0 );

	float r2=0.25f*sizeofprojection*npad*sizeofprojection*npad;
	float r2_at_point;
	
	for(int i=0;i<ellipse_length_x_int;i++) {
		for(int j=-1*ellipse_length_y_int+1; j<=ellipse_length_y_int; j++) {
        	    
			r2_at_point=i*xscale*i*xscale+j*yscale*j*yscale;
			if(r2_at_point<=r2 && ! ((0==i) && (j<0))) {
				
				float ctf_value = ctf_store_new::get_ctf( r2_at_point );
				coordinate_2d_square[0] = xscale*float(i);
				coordinate_2d_square[1] = yscale*float(j);
				float xnew = coordinate_2d_square[0]*trans[0][0] + coordinate_2d_square[1]*trans[1][0];
				float ynew = coordinate_2d_square[0]*trans[0][1] + coordinate_2d_square[1]*trans[1][1];
				float znew = coordinate_2d_square[0]*trans[0][2] + coordinate_2d_square[1]*trans[1][2];
				coordinate_3dnew[0] = xnew*xratio;
				coordinate_3dnew[1] = ynew*yratio;
				coordinate_3dnew[2] = znew*zratio;
				
				//bilinear interpolation
				float xp = coordinate_2d_square[0];
				float yp = ( coordinate_2d_square[1] >= 0) ? coordinate_2d_square[1]+1 : nxyz+coordinate_2d_square[1]+1;
				std::complex<float> lin_interpolated(0,0);
				int xlow=int(xp),xhigh=int(xp)+1;
				int ylow=int(yp),yhigh=int(yp)+1;
				float tx=xp-xlow,ty=yp-ylow;

				
				if(j == -1) {
					
					if(ylow<yp)
						lin_interpolated=myfft->cmplx(xlow,ylow)*(1-tx)*(1-ty) + myfft->cmplx(xlow,yhigh)*(1-tx)*ty
						+ myfft->cmplx(xhigh,ylow)*tx*(1-ty) + myfft->cmplx(xhigh,yhigh)*tx*ty;
					else
						lin_interpolated=myfft->cmplx(xlow,ylow)*(1-tx)
 		  				+ myfft->cmplx(xhigh,ylow)*tx;
									
				} else {
						lin_interpolated=myfft->cmplx(xlow,ylow)*(1-tx)*(1-ty) + myfft->cmplx(xlow,yhigh)*(1-tx)*ty
						+ myfft->cmplx(xhigh,ylow)*tx*(1-ty)+ myfft->cmplx(xhigh,yhigh)*tx*ty;
				}
					
				c1 = lin_interpolated;
				
				//now nearest neighborhood interpolation
				
				std::complex<float> btq;
				if ( coordinate_3dnew[0] < 0.) {
					coordinate_3dnew[0] = -coordinate_3dnew[0];
					coordinate_3dnew[1] = -coordinate_3dnew[1];
					coordinate_3dnew[2] = -coordinate_3dnew[2];
					btq = conj(c1);
					} else {
						btq = c1;
					}
				int ixn = int(coordinate_3dnew[0] + 0.5 + nx) - nx;
				int iyn = int(coordinate_3dnew[1] + 0.5 + ny) - ny;
				int izn = int(coordinate_3dnew[2] + 0.5 + nz) - nz;

				int iza, iya;
				if (izn >= 0)  iza = izn + 1;
				else	       iza = nz + izn + 1;

				if (iyn >= 0) iya = iyn + 1;
				else	      iya = ny + iyn + 1;

				if(remove > 0 ) {
					cmplx(ixn,iya,iza) -= btq * mult;
					(*w)(ixn,iya,iza) -= ctf_value*ctf_value*mult;
					} else {
					cmplx(ixn,iya,iza) += btq * mult;
					(*w)(ixn,iya,iza) += ctf_value*ctf_value*mult;
					}
					
				}
			}
        		    
		}


	//end insert rectanular fft
		
	set_array_offsets(saved_offsets);
	myfft->set_array_offsets(myfft_saved_offsets);
	EXITFUNC;

}


/*
Data::onelinenn_ctf(int j, int n, int n2, EMData* w, EMData* bi, const Transform& tf, int mult) {
//std::cout<<"   onelinenn_ctf  "<<j<<"  "<<n<<"  "<<n2<<"  "<<std::endl;

	int remove = bi->get_attr_default( "remove", 0 );

	int jp = (j >= 0) ? j+1 : n+j+1;
	// loop over x
	for (int i = 0; i <= n2; i++) {
	        int r2 = i*i+j*j;
		if ( (r2<n*n/4) && !( (0==i) && (j<0) ) ) {
			float  ctf = ctf_store::get_ctf( r2 );
			float xnew = i*tf[0][0] + j*tf[1][0];
			float ynew = i*tf[0][1] + j*tf[1][1];
			float znew = i*tf[0][2] + j*tf[1][2];
			std::complex<float> btq;
			if (xnew < 0.) {
				xnew = -xnew;
				ynew = -ynew;
				znew = -znew;
				btq = conj(bi->cmplx(i,jp));
			} else  btq = bi->cmplx(i,jp);
			int ixn = int(xnew + 0.5 + n) - n;
			int iyn = int(ynew + 0.5 + n) - n;
			int izn = int(znew + 0.5 + n) - n;
			if ((ixn <= n2) && (iyn >= -n2) && (iyn <= n2) && (izn >= -n2) && (izn <= n2)) {
				if (ixn >= 0) {
					int iza, iya;
					if (izn >= 0)  iza = izn + 1;
					else           iza = n + izn + 1;

					if (iyn >= 0) iya = iyn + 1;
					else          iya = n + iyn + 1;

					if(remove > 0 ) {
						cmplx(ixn,iya,iza) -= btq*ctf*float(mult);
						(*w)(ixn,iya,iza)  -= ctf*ctf*mult;
					} else {
						cmplx(ixn,iya,iza) += btq*ctf*float(mult);
						(*w)(ixn,iya,iza)  += ctf*ctf*mult;
					}

				       //	std::cout<<"    "<<j<<"  "<<ixn<<"  "<<iya<<"  "<<iza<<"  "<<ctf<<std::endl;
				} else {
					int izt, iyt;
					if (izn > 0) izt = n - izn + 1;
					else         izt = -izn + 1;

					if (iyn > 0) iyt = n - iyn + 1;
					else         iyt = -iyn + 1;

                    if( remove > 0 ) {
					    cmplx(-ixn,iyt,izt) -= conj(btq)*ctf*float(mult);
					    (*w)(-ixn,iyt,izt)  -= ctf*ctf*float(mult);
                    } else {
					    cmplx(-ixn,iyt,izt) += conj(btq)*ctf*float(mult);
					    (*w)(-ixn,iyt,izt)  += ctf*ctf*float(mult);
					}

				        //	std::cout<<" *  " << j << "  " <<-ixn << "  " << iyt << "  " << izt << "  " << ctf <<std::endl;
				}
			}
		}
	}
}
*/


void EMData::nn_SSNR_ctf(EMData* wptr, EMData* wptr2, EMData* wptr3, EMData* myfft, const Transform& tf, float)
{
	/***   Preparing terms for SSNR
	      m_wvolume F^3D Wiener volume
	     wptr   ctf^2
	    wptr5  ctf^2*|P^2D->3D(F^3D)|^2
	   wptr4  2*Real(conj(F_k^2D)*ctf*P^2D->3D(F^3D))
	  wptr2  F_k^2D*conj(F_k^2D) or |F_k^2D|^2
	  Kn is counted in the previous routine, and won't be
	 calculated any more.
	                                            ***/
	ENTERFUNC;
	int nxc = attr_dict["nxc"];
	vector<int> saved_offsets = get_array_offsets();
	vector<int> myfft_saved_offsets = myfft->get_array_offsets();
	set_array_offsets(0,1,1);
	myfft->set_array_offsets(0,1);

	Ctf* ctf = myfft->get_attr("ctf");
	ctf_store::init( ny, ctf );
	int iymin = is_fftodd() ? -ny/2 : -ny/2 + 1;
	int iymax = ny/2;
	int izmin = is_fftodd() ? -nz/2 : -nz/2 + 1;
	int izmax = nz/2;
//	std::complex<float> tmpq, tmp2;
	for (int iy = iymin; iy <= iymax; iy++) {
		int jp = iy >= 0 ? iy+1 : ny+iy+1; //checked, works for both odd and even
		for (int ix = 0; ix <= nxc; ix++) {
			int r2 = ix*ix+iy*iy;
        		if (( 4*r2 < ny*ny ) && !( ix == 0 && iy < 0 ) ) {
				float  ctf = ctf_store::get_ctf( r2, ix, iy )*10.f;// ???PAP
				float xnew = ix*tf[0][0] + iy*tf[1][0];
				float ynew = ix*tf[0][1] + iy*tf[1][1];
				float znew = ix*tf[0][2] + iy*tf[1][2];
				std::complex<float> btq;
				if (xnew < 0.0) {
					xnew = -xnew; // ensures xnew>=0.0
					ynew = -ynew;
					znew = -znew;
					btq = conj(myfft->cmplx(ix,jp));
				} else  {
					btq = myfft->cmplx(ix,jp);
				}
				int ixn = int(xnew + 0.5 + nx) - nx; // ensures ixn >= 0
				int iyn = int(ynew + 0.5 + ny) - ny;
				int izn = int(znew + 0.5 + nz) - nz;
				if ((ixn <= nxc) && (iyn >= iymin) && (iyn <= iymax) && (izn >= izmin) && (izn <= izmax)) {
					if (ixn >= 0) {
						int iza, iya;
						if (izn >= 0) iza = izn + 1;
						else          iza = nz + izn + 1;

						if (iyn >= 0) iya = iyn + 1;
						else          iya = ny + iyn + 1;

						cmplx(ixn,iya,iza)    += btq*ctf;
						(*wptr)(ixn,iya,iza)  += ctf*ctf;
						(*wptr2)(ixn,iya,iza) += std::norm(btq);
						(*wptr3)(ixn,iya,iza) += 1;
					} else {
						int izt, iyt;
						if (izn > 0)  izt = nz - izn + 1;
						else          izt = -izn + 1;

						if (iyn > 0) iyt = ny - iyn + 1;
						else         iyt = -iyn + 1;

						cmplx(-ixn,iyt,izt)    += std::conj(btq)*ctf;
						(*wptr) (-ixn,iyt,izt) += ctf*ctf;
						(*wptr2)(-ixn,iyt,izt) += std::norm(btq);
						(*wptr3)(-ixn,iyt,izt) += 1;
					}
				}
			}
		}
	}
	set_array_offsets(saved_offsets);
	myfft->set_array_offsets(myfft_saved_offsets);
	if(ctf) {delete ctf; ctf=0;}
	EXITFUNC;
}

/*void EMData::nn_wiener(EMData* wptr, EMData* wptr3, EMData* myfft, const Transform& tf, int)
{
     // Wiener volume calculating routine Counting Kn

	ENTERFUNC;
	int nxc = attr_dict["nxc"];
	vector<int> saved_offsets = get_array_offsets();
	vector<int> myfft_saved_offsets = myfft->get_array_offsets();
	set_array_offsets(0,1,1);
       	myfft->set_array_offsets(0,1);
	// if( ! ctf_store::inited() )
        float Cs           = myfft->get_attr( "Cs" );
        float pixel        = myfft->get_attr( "Pixel_size" );
        float voltage      = myfft->get_attr( "voltage");
        float amp_contrast = myfft->get_attr( "amp_contrast" );
        float b_factor     = 0.0;
        ctf_store::init( ny, voltage, pixel, Cs, amp_contrast, b_factor );
        float defocus = myfft->get_attr( "defocus" );
	int iymin = is_fftodd() ? -ny/2 : -ny/2 + 1 ;
	int iymax = ny/2;
	int izmin = is_fftodd() ? -nz/2 : -nz/2 + 1 ;
	int izmax = nz/2;
	for (int iy = iymin; iy <= iymax; iy++) {
		int jp = iy >= 0 ? iy+1 : ny+iy+1; //checked, works for both odd and even
		for (int ix = 0; ix <= nxc; ix++) {
			int r2 = ix*ix+iy*iy;
        		if (( 4*r2 < ny*ny ) && !( ix == 0 && iy < 0 ) )
			{
				float  ctf = ctf_store::get_ctf( defocus, r2 );
				float xnew = ix*tf[0][0] + iy*tf[1][0];
				float ynew = ix*tf[0][1] + iy*tf[1][1];
				float znew = ix*tf[0][2] + iy*tf[1][2];
				std::complex<float> btq;
				if (xnew < 0.0) {
					xnew = -xnew; // ensures xnew>=0.0
					ynew = -ynew;
					znew = -znew;
					btq = conj(myfft->cmplx(ix,jp));
				} else {
					btq = myfft->cmplx(ix,jp);
				}
				int ixn = int(xnew + 0.5 + nx) - nx; // ensures ixn >= 0
				int iyn = int(ynew + 0.5 + ny) - ny;
				int izn = int(znew + 0.5 + nz) - nz;
				if ((ixn <= nxc) && (iyn >= iymin) && (iyn <= iymax) && (izn >= izmin) && (izn <= izmax)) {
					if (ixn >= 0) {
						int iza, iya;
						if (izn >= 0)	iza = izn + 1;
						else			iza = nz + izn + 1;

						if (iyn >= 0)	iya = iyn + 1;
						else			iya = ny + iyn + 1;

						cmplx(ixn,iya,iza)    += btq*ctf;
						(*wptr)(ixn,iya,iza)  += ctf*ctf;
						(*wptr3)(ixn,iya,iza) += 1.0;
					} else {
						int izt, iyt;
						if (izn > 0)	izt = nz - izn + 1;
						else			izt = -izn + 1;

						if (iyn > 0)	iyt = ny - iyn + 1;
						else			iyt = -iyn + 1;

\						cmplx(-ixn,iyt,izt)    += conj(btq)*ctf;
						(*wptr)(-ixn,iyt,izt)  += ctf*ctf;
						(*wptr3)(-ixn,iyt,izt) += 1.0;
					}
				}
			}
		}
	}
	set_array_offsets(saved_offsets);
	myfft->set_array_offsets(myfft_saved_offsets);
	EXITFUNC;
}*/

EMData* EMData::rot_scale_trans2D(float angDeg, float delx, float dely, float scale) { // quadratic, no background, 2D
	float ang=angDeg*M_PI/180.0f;
	if (1 >= ny)
		throw ImageDimensionException("Can't rotate 1D image");
	if (nz<2) {
		vector<int> saved_offsets = get_array_offsets();
		set_array_offsets(0,0,0);
		if (0.0f == scale) scale = 1.0f; // silently fix common user error
		EMData* ret = copy_head();
		delx = restrict2(delx, nx);
		dely = restrict2(dely, ny);
		// center of image
		int xc = nx/2;
		int yc = ny/2;
		// shifted center for rotation
		float shiftxc = xc + delx;
		float shiftyc = yc + dely;
		// trig
		float cang = cos(ang);
		float sang = sin(ang);
			for (int iy = 0; iy < ny; iy++) {
				float y = float(iy) - shiftyc;
				float ycang = y*cang/scale + yc;
				float ysang = -y*sang/scale + xc;
				for (int ix = 0; ix < nx; ix++) {
					float x = float(ix) - shiftxc;
					float xold = x*cang/scale + ysang ;
					float yold = x*sang/scale + ycang ;
					//  quadri is taking care of cyclic count
					(*ret)(ix,iy) = Util::quadri(xold+1.0f, yold+1.0f, nx, ny, get_data());
					   //have to add one as quadri uses Fortran counting
				}
			}
		set_array_offsets(saved_offsets);
		return ret;
	} else {
		throw ImageDimensionException("Volume not currently supported");
	}
}

EMData* EMData::rot_scale_trans2D_background(float angDeg, float delx, float dely, float scale) { // quadratic, no background, 2D
    float ang=angDeg*M_PI/180.0f;
	if (1 >= ny)
		throw ImageDimensionException("Can't rotate 1D image");
	if (nz<2) {
		vector<int> saved_offsets = get_array_offsets();
		set_array_offsets(0,0,0);
		if (0.0f == scale) scale = 1.0f; // silently fix common user error
		EMData* ret = copy_head();
		delx = restrict2(delx, nx);
		dely = restrict2(dely, ny);
		// center of image
		int xc = nx/2;
		int yc = ny/2;
		// shifted center for rotation
		float shiftxc = xc + delx;
		float shiftyc = yc + dely;
		// trig
		float cang = cos(ang);
		float sang = sin(ang);
			for (int iy = 0; iy < ny; iy++) {
				float y = float(iy) - shiftyc;
				float ycang = y*cang/scale + yc;
				float ysang = -y*sang/scale + xc;
				for (int ix = 0; ix < nx; ix++) {
					float x = float(ix) - shiftxc;
					float xold = x*cang/scale + ysang ;
					float yold = x*sang/scale + ycang ;
					//  in quadri_background, wrap around is not done circulantly; if (xold,yold) is not in the image, then it's replaced by (ix,iy)
					(*ret)(ix,iy) = Util::quadri_background(xold+1.0f, yold+1.0f, nx, ny, get_data(),ix+1,iy+1);
					   //have to add one as quadri uses Fortran counting
				}
			}
		set_array_offsets(saved_offsets);
		return ret;
	} else {
		throw ImageDimensionException("Volume not currently supported");
	}
}

#define in(i,j,k)          in[i+(j+(k*ny))*(size_t)nx]
EMData*
EMData::rot_scale_trans(const Transform &RA, EMData* ret) {
//EMData::rot_scale_trans(const Transform &RA) {

//	EMData* ret = copy_head();

    int ret_is_initially_null = ret == NULL; 

	if (ret == NULL) ret = copy_head();
	
	float *in = this->get_data();
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
	Vec3f translations = RA.get_trans();
	Transform RAinv = RA.inverse();

	if (1 >= ny)  throw ImageDimensionException("Can't rotate 1D image");
	if (nz < 2) {
		float  p1, p2, p3, p4;
		float delx = translations.at(0);
		float dely = translations.at(1);
		delx = restrict2(delx, nx);
		dely = restrict2(dely, ny);
		int xc = nx/2;
		int yc = ny/2;
	//         shifted center for rotation
		float shiftxc = xc + delx;
		float shiftyc = yc + dely;
			for (int iy = 0; iy < ny; iy++) {
				float y = float(iy) - shiftyc;
				float ysang = y*RAinv[0][1]+xc;
				float ycang = y*RAinv[1][1]+yc;
				for (int ix = 0; ix < nx; ix++) {
					float x = float(ix) - shiftxc;
					float xold = x*RAinv[0][0] + ysang;
					float yold = x*RAinv[1][0] + ycang;

					xold = restrict1(xold, nx);
					yold = restrict1(yold, ny);

					int xfloor = int(xold);
					int yfloor = int(yold);
					float t = xold-xfloor;
					float u = yold-yfloor;
					if(xfloor == nx -1 && yfloor == ny -1) {

						p1 =in[xfloor   + yfloor*ny];
						p2 =in[ yfloor*ny];
						p3 =in[0];
						p4 =in[xfloor];
					} else if(xfloor == nx - 1) {

						p1 =in[xfloor   + yfloor*ny];
						p2 =in[           yfloor*ny];
						p3 =in[          (yfloor+1)*ny];
						p4 =in[xfloor   + (yfloor+1)*ny];
					} else if(yfloor == ny - 1) {

						p1 =in[xfloor   + yfloor*ny];
						p2 =in[xfloor+1 + yfloor*ny];
						p3 =in[xfloor+1 ];
						p4 =in[xfloor   ];
					} else {
						p1 =in[xfloor   + yfloor*ny];
						p2 =in[xfloor+1 + yfloor*ny];
						p3 =in[xfloor+1 + (yfloor+1)*ny];
						p4 =in[xfloor   + (yfloor+1)*ny];
					}
					if (ret_is_initially_null) (*ret)(ix,iy) = p1 + u * ( p4 - p1) + t * ( p2 - p1 + u *(p3-p2-p4+p1));
					else (*ret)(ix,iy) += p1 + u * ( p4 - p1) + t * ( p2 - p1 + u *(p3-p2-p4+p1));
				} //ends x loop
			} // ends y loop
			set_array_offsets(saved_offsets);
			return ret;
	} else {
	//		 This begins the 3D version tri-linear interpolation.

		float delx = translations.at(0);
		float dely = translations.at(1);
		float delz = translations.at(2);
		delx = restrict2(delx, nx);
		dely = restrict2(dely, ny);
		delz = restrict2(delz, nz);
		int xc = nx/2;
		int yc = ny/2;
		int zc = nz/2;
	//         shifted center for rotation
		float shiftxc = xc + delx;
		float shiftyc = yc + dely;
		float shiftzc = zc + delz;

		for (int iz = 0; iz < nz; iz++) {
			float z = float(iz) - shiftzc;
			float xoldz = z*RAinv[0][2]+xc;
			float yoldz = z*RAinv[1][2]+yc;
			float zoldz = z*RAinv[2][2]+zc;
			for (int iy = 0; iy < ny; iy++) {
				float y = float(iy) - shiftyc;
				float xoldzy = xoldz + y*RAinv[0][1] ;
				float yoldzy = yoldz + y*RAinv[1][1] ;
				float zoldzy = zoldz + y*RAinv[2][1] ;
				for (int ix = 0; ix < nx; ix++) {
					float x = float(ix) - shiftxc;
					float xold = xoldzy + x*RAinv[0][0] ;
					float yold = yoldzy + x*RAinv[1][0] ;
					float zold = zoldzy + x*RAinv[2][0] ;

					xold = restrict1(xold, nx);
					yold = restrict1(yold, ny);
					zold = restrict1(zold, nz);


					int IOX = int(xold);
					int IOY = int(yold);
					int IOZ = int(zold);

					int IOXp1 = std::min( nx-1 ,IOX+1);

					int IOYp1 = std::min( ny-1 ,IOY+1);

					int IOZp1 = std::min( nz-1 ,IOZ+1);

					float dx = xold-IOX;
					float dy = yold-IOY;
					float dz = zold-IOZ;

					float a1 = in(IOX,IOY,IOZ);
					float a2 = in(IOXp1,IOY,IOZ) - in(IOX,IOY,IOZ);
					float a3 = in(IOX,IOYp1,IOZ) - in(IOX,IOY,IOZ);
					float a4 = in(IOX,IOY,IOZp1) - in(IOX,IOY,IOZ);
					float a5 = in(IOX,IOY,IOZ) - in(IOXp1,IOY,IOZ) - in(IOX,IOYp1,IOZ) + in(IOXp1,IOYp1,IOZ);
					float a6 = in(IOX,IOY,IOZ) - in(IOXp1,IOY,IOZ) - in(IOX,IOY,IOZp1) + in(IOXp1,IOY,IOZp1);
					float a7 = in(IOX,IOY,IOZ) - in(IOX,IOYp1,IOZ) - in(IOX,IOY,IOZp1) + in(IOX,IOYp1,IOZp1);
					float a8 = in(IOXp1,IOY,IOZ) + in(IOX,IOYp1,IOZ)+ in(IOX,IOY,IOZp1)
							- in(IOX,IOY,IOZ)- in(IOXp1,IOYp1,IOZ) - in(IOXp1,IOY,IOZp1)
							- in(IOX,IOYp1,IOZp1) + in(IOXp1,IOYp1,IOZp1);
					if (ret_is_initially_null) (*ret)(ix,iy,iz) = a1 + dz*(a4 + a6*dx + (a7 + a8*dx)*dy) + a3*dy + dx*(a2 + a5*dy);
					else (*ret)(ix,iy,iz) += a1 + dz*(a4 + a6*dx + (a7 + a8*dx)*dy) + a3*dy + dx*(a2 + a5*dy);
				} //ends x loop
			} // ends y loop
		} // ends z loop

		set_array_offsets(saved_offsets);
		return ret;
	}
}

// new function added for background option
#define in(i,j,k)          in[i+(j+(k*ny))*(size_t)nx]
EMData*
EMData::rot_scale_trans_background(const Transform &RA, EMData* ret) {

    int ret_is_initially_null = ret == NULL; 

	if (ret == NULL) ret = copy_head();

	float *in = this->get_data();
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
	Vec3f translations = RA.get_trans();
	Transform RAinv = RA.inverse();

	if (1 >= ny)  throw ImageDimensionException("Can't rotate 1D image");
	if (nz < 2) {
		float  p1, p2, p3, p4;
		float delx = translations.at(0);
		float dely = translations.at(1);
		delx = restrict2(delx, nx);
		dely = restrict2(dely, ny);
		int xc = nx/2;
		int yc = ny/2;
	//         shifted center for rotation
		float shiftxc = xc + delx;
		float shiftyc = yc + dely;
			for (int iy = 0; iy < ny; iy++) {
				float y = float(iy) - shiftyc;
				float ysang = y*RAinv[0][1]+xc;
				float ycang = y*RAinv[1][1]+yc;
				for (int ix = 0; ix < nx; ix++) {
					float x = float(ix) - shiftxc;
					float xold = x*RAinv[0][0] + ysang;
					float yold = x*RAinv[1][0] + ycang;

					// if (xold,yold) is outside the image, then let xold = ix and yold = iy

					if ( (xold < 0.0f) || (xold >= (float)(nx)) || (yold < 0.0f) || (yold >= (float)(ny)) ){
						xold = (float)ix;
						yold = (float)iy;
					}

					int xfloor = int(xold);
					int yfloor = int(yold);
					float t = xold-xfloor;
					float u = yold-yfloor;
					if(xfloor == nx -1 && yfloor == ny -1) {

						p1 =in[xfloor   + yfloor*ny];
						p2 =in[ yfloor*ny];
						p3 =in[0];
						p4 =in[xfloor];
					} else if(xfloor == nx - 1) {

						p1 =in[xfloor   + yfloor*ny];
						p2 =in[           yfloor*ny];
						p3 =in[          (yfloor+1)*ny];
						p4 =in[xfloor   + (yfloor+1)*ny];
					} else if(yfloor == ny - 1) {

						p1 =in[xfloor   + yfloor*ny];
						p2 =in[xfloor+1 + yfloor*ny];
						p3 =in[xfloor+1 ];
						p4 =in[xfloor   ];
					} else {

						p1 =in[xfloor   + yfloor*ny];
						p2 =in[xfloor+1 + yfloor*ny];
						p3 =in[xfloor+1 + (yfloor+1)*ny];
						p4 =in[xfloor   + (yfloor+1)*ny];
					}
					if (ret_is_initially_null) (*ret)(ix,iy) = p1 + u * ( p4 - p1) + t * ( p2 - p1 + u *(p3-p2-p4+p1));
					else (*ret)(ix,iy) += p1 + u * ( p4 - p1) + t * ( p2 - p1 + u *(p3-p2-p4+p1));
				} //ends x loop
			} // ends y loop
			set_array_offsets(saved_offsets);
			return ret;
	} else {
//		 This begins the 3D version tri-linear interpolation.

		float delx = translations.at(0);
		float dely = translations.at(1);
		float delz = translations.at(2);
		delx = restrict2(delx, nx);
		dely = restrict2(dely, ny);
		delz = restrict2(delz, nz);
		int xc = nx/2;
		int yc = ny/2;
		int zc = nz/2;
	//         shifted center for rotation
		float shiftxc = xc + delx;
		float shiftyc = yc + dely;
		float shiftzc = zc + delz;

		for (int iz = 0; iz < nz; iz++) {
			float z = float(iz) - shiftzc;
			float xoldz = z*RAinv[0][2]+xc;
			float yoldz = z*RAinv[1][2]+yc;
			float zoldz = z*RAinv[2][2]+zc;
			for (int iy = 0; iy < ny; iy++) {
				float y = float(iy) - shiftyc;
				float xoldzy = xoldz + y*RAinv[0][1] ;
				float yoldzy = yoldz + y*RAinv[1][1] ;
				float zoldzy = zoldz + y*RAinv[2][1] ;
				for (int ix = 0; ix < nx; ix++) {
					float x = float(ix) - shiftxc;
					float xold = xoldzy + x*RAinv[0][0] ;
					float yold = yoldzy + x*RAinv[1][0] ;
					float zold = zoldzy + x*RAinv[2][0] ;

					// if (xold,yold,zold) is outside the image, then let xold = ix, yold = iy and zold=iz

					if ( (xold < 0.0f) || (xold >= (float)(nx)) || (yold < 0.0f) || (yold >= (float)(ny))  || (zold < 0.0f) || (zold >= (float)(nz)) ){
						 xold = (float)ix;
						 yold = (float)iy;
						 zold = (float)iz;
					}

					int IOX = int(xold);
					int IOY = int(yold);
					int IOZ = int(zold);

					int IOXp1 = std::min( nx-1 ,IOX+1);

					int IOYp1 = std::min( ny-1 ,IOY+1);

					int IOZp1 = std::min( nz-1 ,IOZ+1);

					float dx = xold-IOX;
					float dy = yold-IOY;
					float dz = zold-IOZ;

					float a1 = in(IOX,IOY,IOZ);
					float a2 = in(IOXp1,IOY,IOZ) - in(IOX,IOY,IOZ);
					float a3 = in(IOX,IOYp1,IOZ) - in(IOX,IOY,IOZ);
					float a4 = in(IOX,IOY,IOZp1) - in(IOX,IOY,IOZ);
					float a5 = in(IOX,IOY,IOZ) - in(IOXp1,IOY,IOZ) - in(IOX,IOYp1,IOZ) + in(IOXp1,IOYp1,IOZ);
					float a6 = in(IOX,IOY,IOZ) - in(IOXp1,IOY,IOZ) - in(IOX,IOY,IOZp1) + in(IOXp1,IOY,IOZp1);
					float a7 = in(IOX,IOY,IOZ) - in(IOX,IOYp1,IOZ) - in(IOX,IOY,IOZp1) + in(IOX,IOYp1,IOZp1);
					float a8 = in(IOXp1,IOY,IOZ) + in(IOX,IOYp1,IOZ)+ in(IOX,IOY,IOZp1)
							- in(IOX,IOY,IOZ)- in(IOXp1,IOYp1,IOZ) - in(IOXp1,IOY,IOZp1)
							- in(IOX,IOYp1,IOZp1) + in(IOXp1,IOYp1,IOZp1);
					if (ret_is_initially_null) (*ret)(ix,iy,iz) = a1 + dz*(a4 + a6*dx + (a7 + a8*dx)*dy) + a3*dy + dx*(a2 + a5*dy);
					else (*ret)(ix,iy,iz) += a1 + dz*(a4 + a6*dx + (a7 + a8*dx)*dy) + a3*dy + dx*(a2 + a5*dy);
				} //ends x loop
			} // ends y loop
		} // ends z loop

		set_array_offsets(saved_offsets);
		return ret;

	}
}
#undef  in


/*
EMData*
EMData::rot_scale_conv(float ang, float delx, float dely, Util::KaiserBessel& kb) {
	int nxn, nyn, nzn;
	const float scale=0.5;
	float  sum, w;
	if (1 >= ny)
		throw ImageDimensionException("Can't rotate 1D image");
	if (1 < nz)
		throw ImageDimensionException("Volume not currently supported");
	nxn=nx/2;nyn=ny/2;nzn=nz/2;

	int K = kb.get_window_size();
	int kbmin = -K/2;
	int kbmax = -kbmin;
	int kbc = kbmax+1;
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
	EMData* ret = new EMData();
#ifdef _WIN32
	ret->set_size(nxn, _cpp_max(nyn,1), _cpp_max(nzn,1));
#else
	ret->set_size(nxn, std::max(nyn,1), std::max(nzn,1));
#endif	//_WIN32
	ret->to_zero();  //we will leave margins zeroed.
	delx = fmod(delx, float(nxn));
	dely = fmod(dely, float(nyn));
	// center of big image,
	int xc = nxn;
	int ixs = nxn%2;  // extra shift on account of odd-sized images
	int yc = nyn;
	int iys = nyn%2;
	// center of small image
	int xcn = nxn/2;
	int ycn = nyn/2;
	// shifted center for rotation
	float shiftxc = xcn + delx;
	float shiftyc = ycn + dely;
	// bounds if origin at center
	float ymin = -ny/2.0f;
	float xmin = -nx/2.0f;
	float ymax = -ymin;
	float xmax = -xmin;
	if (0 == nx%2) xmax--;
	if (0 == ny%2) ymax--;
	// trig
	float cang = cos(ang);
	float sang = sin(ang);
		for (int iy = 0; iy < nyn; iy++) {
			float y = float(iy) - shiftyc;
			float ycang = y*cang/scale + yc;
			float ysang = -y*sang/scale + xc;
			for (int ix = 0; ix < nxn; ix++) {
				float x = float(ix) - shiftxc;
				float xold = x*cang/scale + ysang-ixs;// have to add the fraction on account on odd-sized images for which Fourier zero-padding changes the center location
				float yold = x*sang/scale + ycang-iys;
				int inxold = int(Util::round(xold)); int inyold = int(Util::round(yold));
				     sum=0.0f;    w=0.0f;
				if(inxold <= kbc || inxold >=nx-kbc-2 || inyold <= kbc || inyold >=ny-kbc-2 )  {
					for (int m2 =kbmin; m2 <=kbmax; m2++){ for (int m1 =kbmin; m1 <=kbmax; m1++) {
						float q = kb.i0win_tab(xold - inxold-m1)*kb.i0win_tab(yold - inyold-m2);
						sum += (*this)((inxold+m1+nx)%nx,(inyold+m2+ny)%ny)*q; w+=q;}
					}
		    	} else {
		    		for (int m2 =kbmin; m2 <=kbmax; m2++){ for (int m1 =kbmin; m1 <=kbmax; m1++) {
		    			float q =kb.i0win_tab(xold - inxold-m1)*kb.i0win_tab(yold - inyold-m2);
		    			sum += (*this)(inxold+m1,inyold+m2)*q;w+=q;}
		    		}
		        }
				(*ret)(ix,iy)=sum/w;
			}
		}
	set_array_offsets(saved_offsets);
	return ret;
}
*/

EMData* EMData::rot_scale_conv(float ang, float delx, float dely, Util::KaiserBessel& kb, float scale_input) {
	int nxn, nyn, nzn;
	if(scale_input == 0.0f) scale_input = 1.0f;
	//const float scale=0.5;
	float  scale = 0.5f*scale_input;
	float  sum, w;
	if (1 >= ny)
		throw ImageDimensionException("Cannot rotate 1D image");
	if (1 < nz)
		throw ImageDimensionException("Volume currently not supported");
	nxn=nx/2; nyn=ny/2; nzn=nz/2;

	int K = kb.get_window_size();
	int kbmin = -K/2;
	int kbmax = -kbmin;
	int kbc = kbmax+1;
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
	EMData* ret = this->copy_head();
	ret->set_size(nxn, std::max(nyn,1), std::max(nzn,1));
	//ret->to_zero();  //we will leave margins zeroed.
	delx = restrict2(delx, nx);
	dely = restrict2(dely, ny);
	// center of big image,
	int xc  = nxn;
	int ixs = nxn%2;  // extra shift on account of odd-sized images
	int yc  = nyn;
	int iys = nyn%2;
	// center of small image
	int xcn = nxn/2;
	int ycn = nyn/2;
	// shifted center for rotation
	float shiftxc = xcn + delx;
	float shiftyc = ycn + dely;
	// bounds if origin at center
	float ymin = -ny/2.0f;
	float xmin = -nx/2.0f;
	float ymax = -ymin;
	float xmax = -xmin;
	if (0 == nx%2) xmax--;
	if (0 == ny%2) ymax--;

	float   *t = (float*)calloc(kbmax-kbmin+1, sizeof(float));

	// trig
	float cang = cos(ang);
	float sang = sin(ang);
	for (int iy = 0; iy < nyn; iy++) {
		float y = float(iy) - shiftyc;
		float ycang = y*cang/scale + yc;
		float ysang = -y*sang/scale + xc;
		for (int ix = 0; ix < nxn; ix++) {
			float x = float(ix) - shiftxc;
			float xold = x*cang/scale + ysang-ixs;// have to add the fraction on account on odd-sized images for which Fourier zero-padding changes the center location
			float yold = x*sang/scale + ycang-iys;

			xold = restrict1(xold, nx);
			yold = restrict1(yold, ny);

			int inxold = int(Util::round(xold)); int inyold = int(Util::round(yold));
			sum=0.0f;    w=0.0f;
			for (int m1 =kbmin; m1 <=kbmax; m1++) t[m1-kbmin] = kb.i0win_tab(xold - inxold-m1);
			if(inxold <= kbc || inxold >=nx-kbc-2 || inyold <= kbc || inyold >=ny-kbc-2 )  {
				for (int m2 =kbmin; m2 <=kbmax; m2++) {
					float qt = kb.i0win_tab(yold - inyold-m2);
					for (int m1 =kbmin; m1 <=kbmax; m1++) {
						float q = t[m1-kbmin]*qt;
						sum += (*this)((inxold+m1+nx)%nx,(inyold+m2+ny)%ny)*q; w+=q;
					}
				}
		    } else {
				for (int m2 =kbmin; m2 <=kbmax; m2++) {
					float qt = kb.i0win_tab(yold - inyold-m2);
			  		for (int m1 =kbmin; m1 <=kbmax; m1++) {
						float q = t[m1-kbmin]*qt;
						sum += (*this)(inxold+m1,inyold+m2)*q; w+=q;}
					}
		    	}
			(*ret)(ix,iy)=sum/w;
		}
	}
	if (t) free(t);
	set_array_offsets(saved_offsets);
	return ret;
}

// Notes by Yang on 10/02/07
// This function is at first just a test, but I found it is slightly faster (about 10%) than rot_scale_conv_new(), so I decided to retain it.
EMData* EMData::rot_scale_conv7(float ang, float delx, float dely, Util::KaiserBessel& kb, float scale_input) {
	int nxn, nyn, nzn;
	float  scale = 0.5f*scale_input;
	float  sum, w;
	if (1 >= ny)
		throw ImageDimensionException("Can't rotate 1D image");
	if (1 < nz)
		throw ImageDimensionException("Volume not currently supported");
	nxn = nx/2; nyn=ny/2; nzn=nz/2;

	int K = kb.get_window_size();
	int kbmin = -K/2;
	int kbmax = -kbmin;
	int kbc = kbmax+1;
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
	EMData* ret = this->copy_head();
	ret->set_size(nxn, std::max(nyn,1), std::max(nzn,1));
	//ret->to_zero();  //we will leave margins zeroed.
	delx = restrict2(delx, nx);
	dely = restrict2(dely, ny);
	// center of big image,
	int xc = nxn;
	int ixs = nxn%2;  // extra shift on account of odd-sized images
	int yc = nyn;
	int iys = nyn%2;
	// center of small image
	int xcn = nxn/2;
	int ycn = nyn/2;
	// shifted center for rotation
	float shiftxc = xcn + delx;
	float shiftyc = ycn + dely;
	// bounds if origin at center
	float ymin = -ny/2.0f;
	float xmin = -nx/2.0f;
	float ymax = -ymin;
	float xmax = -xmin;
	if (0 == nx%2) xmax--;
	if (0 == ny%2) ymax--;

	float   *t = (float*)calloc(kbmax-kbmin+1, sizeof(float));

	// trig
	float cang = cos(ang);
	float sang = sin(ang);
	for (int iy = 0; iy < nyn; iy++) {
		float y = float(iy) - shiftyc;
		float ycang = y*cang/scale + yc;
		float ysang = -y*sang/scale + xc;
		for (int ix = 0; ix < nxn; ix++) {
			float x = float(ix) - shiftxc;
			float xold = x*cang/scale + ysang-ixs;// have to add the fraction on account on odd-sized images for which Fourier zero-padding changes the center location
			float yold = x*sang/scale + ycang-iys;

			xold = restrict1(xold, nx);
			yold = restrict1(yold, ny);

			int inxold = int(Util::round(xold)); int inyold = int(Util::round(yold));
			sum=0.0f;    w=0.0f;

			float tablex1 = kb.i0win_tab(xold-inxold+3);
			float tablex2 = kb.i0win_tab(xold-inxold+2);
			float tablex3 = kb.i0win_tab(xold-inxold+1);
			float tablex4 = kb.i0win_tab(xold-inxold);
			float tablex5 = kb.i0win_tab(xold-inxold-1);
			float tablex6 = kb.i0win_tab(xold-inxold-2);
			float tablex7 = kb.i0win_tab(xold-inxold-3);

			float tabley1 = kb.i0win_tab(yold-inyold+3);
			float tabley2 = kb.i0win_tab(yold-inyold+2);
			float tabley3 = kb.i0win_tab(yold-inyold+1);
			float tabley4 = kb.i0win_tab(yold-inyold);
			float tabley5 = kb.i0win_tab(yold-inyold-1);
			float tabley6 = kb.i0win_tab(yold-inyold-2);
			float tabley7 = kb.i0win_tab(yold-inyold-3);

			int x1, x2, x3, x4, x5, x6, x7, y1, y2, y3, y4, y5, y6, y7;

			if(inxold <= kbc || inxold >=nx-kbc-2 || inyold <= kbc || inyold >=ny-kbc-2 )  {
				x1 = (inxold-3+nx)%nx;
				x2 = (inxold-2+nx)%nx;
				x3 = (inxold-1+nx)%nx;
				x4 = (inxold  +nx)%nx;
				x5 = (inxold+1+nx)%nx;
				x6 = (inxold+2+nx)%nx;
				x7 = (inxold+3+nx)%nx;

				y1 = (inyold-3+ny)%ny;
				y2 = (inyold-2+ny)%ny;
				y3 = (inyold-1+ny)%ny;
				y4 = (inyold  +ny)%ny;
				y5 = (inyold+1+ny)%ny;
				y6 = (inyold+2+ny)%ny;
				y7 = (inyold+3+ny)%ny;
		    	} else {
				x1 = inxold-3;
				x2 = inxold-2;
				x3 = inxold-1;
				x4 = inxold;
				x5 = inxold+1;
				x6 = inxold+2;
				x7 = inxold+3;

				y1 = inyold-3;
				y2 = inyold-2;
				y3 = inyold-1;
				y4 = inyold;
				y5 = inyold+1;
				y6 = inyold+2;
				y7 = inyold+3;
		    	}
			sum    =   ( (*this)(x1,y1)*tablex1 + (*this)(x2,y1)*tablex2 + (*this)(x3,y1)*tablex3 +
		        	     (*this)(x4,y1)*tablex4 + (*this)(x5,y1)*tablex5 + (*this)(x6,y1)*tablex6 +
			   	     (*this)(x7,y1)*tablex7 ) * tabley1 +
			   	   ( (*this)(x1,y2)*tablex1 + (*this)(x2,y2)*tablex2 + (*this)(x3,y2)*tablex3 +
			   	     (*this)(x4,y2)*tablex4 + (*this)(x5,y2)*tablex5 + (*this)(x6,y2)*tablex6 +
		   		     (*this)(x7,y2)*tablex7 ) * tabley2 +
			   	   ( (*this)(x1,y3)*tablex1 + (*this)(x2,y3)*tablex2 + (*this)(x3,y3)*tablex3 +
			   	     (*this)(x4,y3)*tablex4 + (*this)(x5,y3)*tablex5 + (*this)(x6,y3)*tablex6 +
			   	     (*this)(x7,y3)*tablex7 ) * tabley3 +
		   		   ( (*this)(x1,y4)*tablex1 + (*this)(x2,y4)*tablex2 + (*this)(x3,y4)*tablex3 +
			   	     (*this)(x4,y4)*tablex4 + (*this)(x5,y4)*tablex5 + (*this)(x6,y4)*tablex6 +
			   	     (*this)(x7,y4)*tablex7 ) * tabley4 +
			   	   ( (*this)(x1,y5)*tablex1 + (*this)(x2,y5)*tablex2 + (*this)(x3,y5)*tablex3 +
		   		     (*this)(x4,y5)*tablex4 + (*this)(x5,y5)*tablex5 + (*this)(x6,y5)*tablex6 +
			   	     (*this)(x7,y5)*tablex7 ) * tabley5 +
			   	   ( (*this)(x1,y6)*tablex1 + (*this)(x2,y6)*tablex2 + (*this)(x3,y6)*tablex3 +
			   	     (*this)(x4,y6)*tablex4 + (*this)(x5,y6)*tablex5 + (*this)(x6,y6)*tablex6 +
		   		     (*this)(x7,y6)*tablex7 ) * tabley6 +
			   	   ( (*this)(x1,y7)*tablex1 + (*this)(x2,y7)*tablex2 + (*this)(x3,y7)*tablex3 +
			   	     (*this)(x4,y7)*tablex4 + (*this)(x5,y7)*tablex5 + (*this)(x6,y7)*tablex6 +
			   	     (*this)(x7,y7)*tablex7 ) * tabley7;

			w = (tablex1+tablex2+tablex3+tablex4+tablex5+tablex6+tablex7) *
			    (tabley1+tabley2+tabley3+tabley4+tabley5+tabley6+tabley7);

			(*ret)(ix,iy)=sum/w;
		}
	}
	if (t) free(t);
	set_array_offsets(saved_offsets);
	return ret;
}

EMData*
EMData::rot_fvol(const Transform &RA, EMData* ret, int radius) {
//EMData::rot_fvol(const Transform &RA) {
// Note data has to be shifted to corners by n/2
//	EMData* ret = copy_head();

    if( ret == NULL ) {
        ret = copy_head();
	    ret->to_zero();
    }
	ret->set_array_offsets(0,0,0);
    
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
	Transform RAinv = RA.inverse();

	if (1 >= ny)  throw ImageDimensionException("Can't frotate 1D image");
	if (nz < 2) {
		throw ImageDimensionException("Can't frotate 2D image");
	} else {
//		 This begins the 3D version tri-linear interpolation.
		if (is_complex())  {
			if (radius < 1 || radius > nx/2-1) radius = nx/2-1;
			float rm2 = (radius-1)*(radius-1);
			bool iodd = get_attr("is_fftodd");
			int xc = nx/2;
			int yc = ny/2;
			int zc = nz/2;
			int bign = 2*nz;
//cout<<"  "<<nx<<"  "<<ny<<"  "<<nz<<"  "<<xc<<"  "<<yc<<"  "<<zc<<"  "<<bign<<"  "<<rm2<<endl;
			for (int iz = -zc + 1-nz%2; iz < zc; iz++) {
				float xnewz = iz*RAinv[0][2];
				float ynewz = iz*RAinv[1][2];
				float znewz = iz*RAinv[2][2];
				for (int iy = -yc + 1-ny%2; iy < yc; iy++) {
					float xnewzy = xnewz + iy*RAinv[0][1];
					float ynewzy = ynewz + iy*RAinv[1][1];
					float znewzy = znewz + iy*RAinv[2][1];
					for (int ix = 0; ix < xc; ix++) {
						float xnew = xnewzy + ix*RAinv[0][0];
						float ynew = ynewzy + ix*RAinv[1][0];
						float znew = znewzy + ix*RAinv[2][0];
//cout<<"   Before if "<<xnew*xnew + ynew+ynew + znew*znew<<"  "<<rm2<<"  "<<ix<<"  "<<iy<<"  "<<iz<<"  "<<xnew<<"  "<<ynew<<"  "<<znew<<endl;
						if(xnew*xnew + ynew*ynew + znew*znew <= rm2 ) {

							int itz = iz;
							if( itz < 0 ) itz += nz;
							int ity = iy;
							if( ity < 0 ) ity += ny;

							bool flipin;
							if (xnew < 0.) {
								xnew = -xnew;
								ynew = -ynew;
								znew = -znew;
								flipin = true;
							} else  flipin = false;


							int ixn = int(xnew);
							int iyn = int(ynew + bign);
							int izn = int(znew + bign);
							float dx = xnew - ixn;
							float dy = ynew + bign - iyn;
							float dz = znew + bign - izn;
							float qdx = 1.0f - dx;
							float qdy = 1.0f - dy;
							float qdz = 1.0f - dz;

							iyn -= bign;
							izn -= bign;

							int ix1 = ixn + 1;
							if( ix1 < xc) {

								float qq000 = qdx * qdy * qdz;
								float qq010 = qdx *  dy * qdz;
								float qq100 =  dx * qdy * qdz;
								float qq110 =  dx *  dy * qdz;
								float qq001 = qdx * qdy *  dz;
								float qq011 = qdx *  dy *  dz;
								float qq101 =  dx * qdy *  dz;
								float qq111 =  dx *  dy *  dz;

								int iza = izn;
								if (izn < 0)  iza += nz;

								int iya = iyn;
								if (iyn < 0)  iya += ny;

								int iz1 = iza + 1;
								if (iz1 >= nz) iz1 -= nz;

								int iy1 = iya + 1;
								if (iy1 >= ny) iy1 -= ny;
	//cout<<" XXX "<<ix<<"  "<<iy<<"  "<<iz<<"  "<<xnew<<"  "<<ynew<<"  "<<znew<<"  "<<ix<<"  "<<ity<<"  "<<itz<<"  "<<ixn<<"  "<<iya<<"  "<<iza<<endl;
	//if(ixn<0 or ixn >=xc)  cout<<"  error   ixn  "<<ixn<<endl;
	//if(ix1<0 or ix1 >=xc){  cout<<"  error   ix1  "<<ix1<<endl;cout<<" XXX "<<ix<<"  "<<iy<<"  "<<iz<<"  "<<xnew<<"  "<<ynew<<"  "<<znew<<"  "<<ix<<"  "<<ity<<"  "<<itz<<"  "<<ixn<<"  "<<iya<<"  "<<iza<<endl;}
	//if(iya<0 or iya >=ny)  cout<<"  error   iya  "<<iya<<endl;
	//if(iy1<0 or iy1 >=ny)  cout<<"  error   iy1  "<<iy1<<endl;
	//if(iza<0 or iza >=nz)  cout<<"  error   iza  "<<iza<<endl;
	//if(iz1<0 or iz1 >=nz)  cout<<"  error   iz1  "<<iz1<<endl;
								std::complex<float> btq = qq000 * cmplx(ixn, iya, iza) + qq010 * cmplx(ixn, iy1, iza) + qq100 * cmplx(ix1, iya, iza)
										+ qq110 * cmplx(ix1, iy1, iza) + qq001 * cmplx(ixn, iya, iz1) + qq011 * cmplx(ixn, iy1, iz1)
										+ qq101 * cmplx(ix1, iya, iz1) + qq111 * cmplx(ix1, iy1, iz1);

	//if(ix<0 or ix >=xc)  cout<<"  error   ix  "<<ix<<endl;
	//if(ity<0 or ity >=ny)  cout<<"  error   ity  "<<ity<<endl;
	//if(itz<0 or itz >=nz)  cout<<"  error   itz  "<<itz<<endl;
	
								if( flipin )  ret->cmplx(ix,ity,itz) += conj(btq);
                                else  ret->cmplx(ix,ity,itz) += btq;
							}
						}
					} //ends x loop
				} // ends y loop
			} // ends z loop
		} else {
			if (radius < 1 || radius > nx-1) radius = nx-1;
			float rm2 = (radius-1)*(radius-1);
			//cout<< "  radius  "<<radius<<"   "<<rm2<<endl;
			int xc = nx/2;
			int yc = ny/2;
			int zc = nz/2;
			int bign = 2*zc;

			for (int iz = -zc + 1-nz%2; iz < zc; iz++) {
				float xnewz = iz*RAinv[0][2];
				float ynewz = iz*RAinv[1][2];
				float znewz = iz*RAinv[2][2];
				for (int iy = -yc + 1-ny%2; iy < yc; iy++) {
					float xnewzy = xnewz + iy*RAinv[0][1];
					float ynewzy = ynewz + iy*RAinv[1][1];
					float znewzy = znewz + iy*RAinv[2][1];
					for (int ix = 0; ix < nx; ix++) {
						float xnew = xnewzy + ix*RAinv[0][0];
						float ynew = ynewzy + ix*RAinv[1][0];
						float znew = znewzy + ix*RAinv[2][0];
//cout<<"   Before if "<<xnew*xnew + ynew*ynew + znew*znew<<"  "<<rm2<<"  "<<ix<<"  "<<iy<<"  "<<iz<<"  "<<xnew<<"  "<<ynew<<"  "<<znew<<endl;
						if(xnew*xnew + ynew*ynew + znew*znew <= rm2 ) {

							int itz = iz;
							if( itz < 0 ) itz += nz;
							int ity = iy;
							if( ity < 0 ) ity += ny;

							if (xnew < 0.) {
								xnew = -xnew;
								ynew = -ynew;
								znew = -znew;
							}

							int ixn = int(xnew);
							int iyn = int(ynew + bign);
							int izn = int(znew + bign);
							float dx = xnew - ixn;
							float dy = ynew + bign - iyn;
							float dz = znew + bign - izn;
							float qdx = 1.0f - dx;
							float qdy = 1.0f - dy;
							float qdz = 1.0f - dz;

							iyn -= bign;
							izn -= bign;

							int ix1 = ixn + 1;
							if( ix1 < nx)  {

								float qq000 = qdx * qdy * qdz;
								float qq010 = qdx *  dy * qdz;
								float qq100 =  dx * qdy * qdz;
								float qq110 =  dx *  dy * qdz;
								float qq001 = qdx * qdy *  dz;
								float qq011 = qdx *  dy *  dz;
								float qq101 =  dx * qdy *  dz;
								float qq111 =  dx *  dy *  dz;

								int iza = izn;
								if (izn < 0)  iza += nz;

								int iya = iyn;
								if (iyn < 0)  iya += ny;

								int iz1 = iza + 1;
								if (iz1 >= nz) iz1 -= nz;

								int iy1 = iya + 1;
								if (iy1 >= ny) iy1 -= ny;
	//cout<<" XXX "<<ix<<"  "<<iy<<"  "<<iz<<" new: "<<xnew<<"  "<<ynew<<"  "<<znew<<"  "<<ix<<"  "<<ity<<"  "<<itz<<"  "<<ixn<<"  "<<iya<<"  "<<iza<<endl;
	//btq = cmplx(ixn, iya, iza);

                                (*ret)(ix,ity,itz) += qq000 * (*this)(ixn, iya, iza) + qq010 * (*this)(ixn, iy1, iza) + qq100 * (*this)(ix1, iya, iza)
                                            + qq110 * (*this)(ix1, iy1, iza) + qq001 * (*this)(ixn, iya, iz1) + qq011 * (*this)(ixn, iy1, iz1)
                                            + qq101 * (*this)(ix1, iya, iz1) + qq111 * (*this)(ix1, iy1, iz1);
                                
							}

						}
					} //ends x loop
				} // ends y loop
			} // ends z loop
		}
		set_array_offsets(saved_offsets);
		return ret;
	}
}

EMData* EMData::downsample(Util::sincBlackman& kb, float scale) {

	/*int M = kb.get_sB_size();
	int kbmin = -M/2;
	int kbmax = -kbmin;*/

	int nxn, nyn, nzn;
	nxn = (int)(nx*scale+0.5f); nyn = (int)(ny*scale+0.5f); nzn = (int)(nz*scale+0.5f);

	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
	EMData* ret = this->copy_head();
	ret->set_size(nxn, std::max(nyn,1), std::max(nzn,1));
	ret->to_zero();  //we will leave margins zeroed.

	// scan new, find pixels in old
	if(nz == 1)  {
		for (int iy =0; iy < nyn; iy++) {
			float y = float(iy)/scale;
			for (int ix = 0; ix < nxn; ix++) {
				float x = float(ix)/scale;
				(*ret)(ix,iy) = this->get_pixel_filtered(x, y, 1.0f, kb);
			}
		}
	}  else  {
		
		for (int iz =0; iz < nzn; iz++) {
			float z = float(iz)/scale;
			for (int iy =0; iy < nyn; iy++) {
				float y = float(iy)/scale;
				for (int ix = 0; ix < nxn; ix++) {
					float x = float(ix)/scale;
					(*ret)(ix,iy,iz) = this->get_pixel_filtered(x, y, z, kb);
				}
			}
		}
	}
	set_array_offsets(saved_offsets);
	return ret;
}


EMData* EMData::rot_scale_conv_new(float ang, float delx, float dely, Util::KaiserBessel& kb, float scale_input) {

	if (scale_input == 0.0f) scale_input = 1.0f;
	float  scale = 0.5f*scale_input;

	if (1 >= ny)
		throw ImageDimensionException("Can't rotate 1D image");
	if (1 < nz)
		throw ImageDimensionException("Use rot_scale_conv_new_3D for volumes");
	int nxn = nx/2; int nyn = ny/2; int nzn = nz/2;

	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
	EMData* ret = this->copy_head();
	ret->set_size(nxn, std::max(nyn,1), std::max(nzn,1));
	//ret->to_zero();  //we will leave margins zeroed.
	delx = restrict2(delx, nx);
	dely = restrict2(dely, ny);
	// center of big image,
	int xc = nxn;
	int ixs = nxn%2;  // extra shift on account of odd-sized images
	int yc = nyn;
	int iys = nyn%2;
	// center of small image
	int xcn = nxn/2;
	int ycn = nyn/2;
	// shifted center for rotation
	float shiftxc = xcn + delx;
	float shiftyc = ycn + dely;
	// bounds if origin at center
	float ymin = -ny/2.0f;
	float xmin = -nx/2.0f;
	float ymax = -ymin;
	float xmax = -xmin;
	if (0 == nx%2) xmax--;
	if (0 == ny%2) ymax--;

	float* data = this->get_data();

	float cang = cos(ang);
	float sang = sin(ang);
	for (int iy = 0; iy < nyn; iy++) {
		float y = float(iy) - shiftyc;
		float ycang = y*cang/scale + yc;
		float ysang = -y*sang/scale + xc;
		for (int ix = 0; ix < nxn; ix++) {
			float x = float(ix) - shiftxc;
			float xold = x*cang/scale + ysang-ixs;// have to add the fraction on account on odd-sized images for which Fourier zero-padding changes the center location
			float yold = x*sang/scale + ycang-iys;

			(*ret)(ix,iy) = Util::get_pixel_conv_new(nx, ny, 1, xold, yold, 1, data, kb);
		}
	}
	set_array_offsets(saved_offsets);
	return ret;
}

EMData* EMData::rot_scale_conv_new_3D(float phi, float theta, float psi, float delx, float dely, float delz, Util::KaiserBessel& kb, float scale_input, bool wrap) {

	if (scale_input == 0.0f) scale_input = 1.0f;
	float  scale = 0.5f*scale_input;

	if (1 >= ny)
		throw ImageDimensionException("Can't rotate 1D image");
	int nxn = nx/2; int nyn = ny/2; int nzn = nz/2;

	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
	EMData* ret = this->copy_head();
	ret->set_size(nxn, std::max(nyn,1), std::max(nzn,1));
	//ret->to_zero();  //we will leave margins zeroed.
	if(wrap){
		delx = restrict2(delx, nx);
		dely = restrict2(dely, ny);
		delz = restrict2(delz, nz);
	}
	// center of big image,
	int xc = nxn;
	int ixs = nxn%2;  // extra shift on account of odd-sized images
	int yc = nyn;
	int iys = nyn%2;
	int zc = nzn;
	int izs = nzn%2;
	// center of small image
	int xcn = nxn/2;
	int ycn = nyn/2;
	int zcn = nzn/2;
	// shifted center for rotation
	float shiftxc = xcn + delx;
	float shiftyc = ycn + dely;
	float shiftzc = zcn + delz;
	// bounds if origin at center
	float zmin = -nz/2.0f;
	float ymin = -ny/2.0f;
	float xmin = -nx/2.0f;
	float zmax = -zmin;
	float ymax = -ymin;
	float xmax = -xmin;
	if (0 == nx%2) xmax--;
	if (0 == ny%2) ymax--;
	if (0 == nz%2) zmax--;

	float* data = this->get_data();

	float cf = cos(phi);   float sf = sin(phi);
	float ct = cos(theta); float st = sin(theta);
	float cp = cos(psi);   float sp = sin(psi);
	// rotation matrix (the transpose is used in the loop to get (xold,yold,zold)):
	float a11 =  cp*ct*cf-sp*sf; float a12 =  cp*ct*sf+sp*cf; float a13 = -cp*st;
	float a21 = -sp*ct*cf-cp*sf; float a22 = -sp*ct*sf+cp*cf; float a23 =  sp*st;
	float a31 =  st*cf;          float a32 =  st*sf;          float a33 =  ct;
	for (int iz = 0; iz < nzn; iz++) {
		float z = (float(iz) - shiftzc)/scale;
		float zco1 = a31*z+xc;
		float zco2 = a32*z+yc;
		float zco3 = a33*z+zc;
		for (int iy = 0; iy < nyn; iy++) {
			float y = (float(iy) - shiftyc)/scale;
			float yco1 = zco1+a21*y;
			float yco2 = zco2+a22*y;
			float yco3 = zco3+a23*y;
			for (int ix = 0; ix < nxn; ix++) {
				float x = (float(ix) - shiftxc)/scale;
				float xold = yco1+a11*x-ixs; //have to add the fraction on account of odd-sized images for which Fourier zero-padding changes the center location
				float yold = yco2+a12*x-iys;
				float zold = yco3+a13*x-izs;
				if(!wrap && (xold<0.0 || xold>nx-1 || yold<0.0 || yold>ny-1 || zold<0.0 || zold>nz-1))
					(*ret)(ix,iy,iz) = 0.0;
				else
					(*ret)(ix,iy,iz) = Util::get_pixel_conv_new(nx, ny, nz, xold, yold, zold, data, kb);
			}
		}
	}
	set_array_offsets(saved_offsets);
	return ret;
}

EMData* EMData::rot_scale_conv_new_background(float ang, float delx, float dely, Util::KaiserBessel& kb, float scale_input) {

	int nxn, nyn, nzn;

	if (scale_input == 0.0f) scale_input = 1.0f;
	float  scale = 0.5f*scale_input;

	if (1 >= ny)
		throw ImageDimensionException("Cannot rotate 1D image");
	if (1 < nz)
		throw ImageDimensionException("Use rot_scale_conv_new_background_3D for volumes");
	nxn = nx/2; nyn = ny/2; nzn = nz/2;

	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
	EMData* ret = this->copy_head();
	ret->set_size(nxn, std::max(nyn,1), std::max(nzn,1));
	//ret->to_zero();  //we will leave margins zeroed.
	delx = restrict2(delx, nx);
	dely = restrict2(dely, ny);
	// center of big image,
	int xc = nxn;
	int ixs = nxn%2;  // extra shift on account of odd-sized images
	int yc = nyn;
	int iys = nyn%2;
	// center of small image
	int xcn = nxn/2;
	int ycn = nyn/2;
	// shifted center for rotation
	float shiftxc = xcn + delx;
	float shiftyc = ycn + dely;
	// bounds if origin at center
	float ymin = -ny/2.0f;
	float xmin = -nx/2.0f;
	float ymax = -ymin;
	float xmax = -xmin;
	if (0 == nx%2) xmax--;
	if (0 == ny%2) ymax--;

	float* data = this->get_data();

	// trig
	float cang = cos(ang);
	float sang = sin(ang);
	for (int iy = 0; iy < nyn; iy++) {
		float y = float(iy) - shiftyc;
		float ycang = y*cang/scale + yc;
		float ysang = -y*sang/scale + xc;
		for (int ix = 0; ix < nxn; ix++) {
			float x = float(ix) - shiftxc;
			float xold = x*cang/scale + ysang-ixs;// have to add the fraction on account on odd-sized images for which Fourier zero-padding changes the center location
			float yold = x*sang/scale + ycang-iys;

			(*ret)(ix,iy) = Util::get_pixel_conv_new_background(nx, ny, 1, xold, yold, 1, data, kb, ix, iy);
		}
	}
	set_array_offsets(saved_offsets);
	return ret;
}

//  This function returns twice resampled image.  It is needed for gridding CCF applications
//  It is slightly modified version of rot_scale_conv_new_background
EMData* EMData::rot_scale_conv_new_background_twice(float ang, float delx, float dely, Util::KaiserBessel& kb, float scale_input) {

	int nxn, nyn, nzn;

	if (scale_input == 0.0f) scale_input = 1.0f;
	float  scale = scale_input;

	if (1 >= ny)
		throw ImageDimensionException("Cannot rotate 1D image");
	if (1 < nz)
		throw ImageDimensionException("Use rot_scale_conv_new_background_3D for volumes");
	nxn = nx; nyn = ny; nzn = nz;

	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
	EMData* ret = this->copy_head();
	ret->set_size(nxn, std::max(nyn,1), std::max(nzn,1));
	//ret->to_zero();  //we will leave margins zeroed.
	delx = restrict2(delx, nx);
	dely = restrict2(dely, ny);
	// center of big image,
	int xc = nxn/2;
	int ixs = nxn%2;  // extra shift on account of odd-sized images
	int yc = nyn/2;
	int iys = nyn%2;
	// center of small image
	int xcn = nxn/2;
	int ycn = nyn/2;
	// shifted center for rotation
	float shiftxc = xcn + delx;
	float shiftyc = ycn + dely;
	// bounds if origin at center
	float ymin = -ny/2.0f;
	float xmin = -nx/2.0f;
	float ymax = -ymin;
	float xmax = -xmin;
	if (0 == nx%2) xmax--;
	if (0 == ny%2) ymax--;

	float* data = this->get_data();

	// trig
	float cang = cos(ang);
	float sang = sin(ang);
	for (int iy = 0; iy < nyn; iy++) {
		float y = float(iy) - shiftyc;
		float ycang = y*cang/scale + yc;
		float ysang = -y*sang/scale + xc;
		for (int ix = 0; ix < nxn; ix++) {
			float x = float(ix) - shiftxc;
			float xold = x*cang/scale + ysang-ixs;// have to add the fraction on account on odd-sized images for which Fourier zero-padding changes the center location
			float yold = x*sang/scale + ycang-iys;

			//(*ret)(ix,iy) = Util::get_pixel_conv_new_background(nx, ny, 1, xold, yold, 1, data, kb, ix, iy);
//cout << "   "<<ix<< "   "<<iy<< "   "<<xold<< "   "<<yold<< "   "<<nx<< "   "<<ny<<endl;
//cout << "   "<<xc<< "   "<<yc<< "   "<<xcn<< "   "<<ycn<< "   "<<shiftxc<< "   "<<shiftyc<<endl;
			(*ret)(ix,iy) = Util::get_pixel_conv_new(nx, ny, 1, xold, yold, 1, data, kb);
		}
	}
	set_array_offsets(saved_offsets);
	return ret;
}


EMData* EMData::rot_scale_conv_new_background_3D(float phi, float theta, float psi, float delx, float dely, float delz, Util::KaiserBessel& kb, float scale_input, bool wrap) {

	if (scale_input == 0.0f) scale_input = 1.0f;
	float  scale = 0.5f*scale_input;

	if (1 >= ny)
		throw ImageDimensionException("Can't rotate 1D image");
	int nxn = nx/2; int nyn = ny/2; int nzn = nz/2;

	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
	EMData* ret = this->copy_head();
	ret->set_size(nxn, std::max(nyn,1), std::max(nzn,1));
	//ret->to_zero();  //we will leave margins zeroed.
	if (wrap){
		delx = restrict2(delx, nx);
		dely = restrict2(dely, ny);
		delz = restrict2(delz, nz);
	}
	// center of big image,
	int xc = nxn;
	int ixs = nxn%2;  // extra shift on account of odd-sized images
	int yc = nyn;
	int iys = nyn%2;
	int zc = nzn;
	int izs = nzn%2;
	// center of small image
	int xcn = nxn/2;
	int ycn = nyn/2;
	int zcn = nzn/2;
	// shifted center for rotation
	float shiftxc = xcn + delx;
	float shiftyc = ycn + dely;
	float shiftzc = zcn + delz;
	// bounds if origin at center
	float zmin = -nz/2.0f;
	float ymin = -ny/2.0f;
	float xmin = -nx/2.0f;
	float zmax = -zmin;
	float ymax = -ymin;
	float xmax = -xmin;
	if (0 == nx%2) xmax--;
	if (0 == ny%2) ymax--;
	if (0 == nz%2) zmax--;

	float* data = this->get_data();

	float cf = cos(phi);   float sf = sin(phi);
	float ct = cos(theta); float st = sin(theta);
	float cp = cos(psi);   float sp = sin(psi);
	// rotation matrix (the transpose is used in the loop to get (xold,yold,zold)):
	float a11 =  cp*ct*cf-sp*sf; float a12 =  cp*ct*sf+sp*cf; float a13 = -cp*st;
	float a21 = -sp*ct*cf-cp*sf; float a22 = -sp*ct*sf+cp*cf; float a23 =  sp*st;
	float a31 =  st*cf;          float a32 =  st*sf;          float a33 =  ct;
	for (int iz = 0; iz < nzn; iz++) {
		float z = (float(iz) - shiftzc)/scale;
		float zco1 = a31*z+xc;
		float zco2 = a32*z+yc;
		float zco3 = a33*z+zc;
		for (int iy = 0; iy < nyn; iy++) {
			float y = (float(iy) - shiftyc)/scale;
			float yco1 = zco1+a21*y;
			float yco2 = zco2+a22*y;
			float yco3 = zco3+a23*y;
			for (int ix = 0; ix < nxn; ix++) {
				float x = (float(ix) - shiftxc)/scale;
				float xold = yco1+a11*x-ixs; //have to add the fraction on account on odd-sized images for which Fourier zero-padding changes the center location
				float yold = yco2+a12*x-iys;
				float zold = yco3+a13*x-izs;
				if(!wrap && (xold<0.0 || xold>nx-1 || yold<0.0 || yold>ny-1 || zold<0.0 || zold>nz-1))
					(*ret)(ix,iy,iz) = 0.0;
				else
					(*ret)(ix,iy,iz) = Util::get_pixel_conv_new_background(nx, ny, nz, xold, yold, zold, data, kb, ix, iy);
			}
		}
	}
	set_array_offsets(saved_offsets);
	return ret;
}


float  EMData::get_pixel_conv(float delx, float dely, float delz, Util::KaiserBessel& kb) {
//  here counting is in C style, so coordinates of the pixel delx should be [0-nx-1]

	int K     = kb.get_window_size();
	int kbmin = -K/2;
	int kbmax = -kbmin;
	int kbc   = kbmax+1;

	float pixel =0.0f;
	float w=0.0f;

	delx = restrict2(delx, nx);
	int inxold = int(Util::round(delx));
	if(ny<2) {  //1D
	 	if(inxold <= kbc || inxold >=nx-kbc-2 )  {
	 		//  loop for ends
         		for (int m1 =kbmin; m1 <=kbmax; m1++) {
	 			float q = kb.i0win_tab(delx - inxold-m1);
	 			pixel += (*this)((inxold+m1+nx)%nx)*q; w+=q;
			}
	 	} else {
         		for (int m1 =kbmin; m1 <=kbmax; m1++) {
	 			float q = kb.i0win_tab(delx - inxold-m1);
	 			pixel += (*this)(inxold+m1)*q; w+=q;
			}
	 	}

	} else if(nz<2) {  // 2D
		dely = restrict2(dely, ny);
		int inyold = int(Util::round(dely));
	 	if(inxold <= kbc || inxold >=nx-kbc-2 || inyold <= kbc || inyold >=ny-kbc-2 )  {
	 		//  loop for strips
         		for (int m2 =kbmin; m2 <=kbmax; m2++){ for (int m1 =kbmin; m1 <=kbmax; m1++) {
	 			float q = kb.i0win_tab(delx - inxold-m1)*kb.i0win_tab(dely - inyold-m2);
	 			pixel += (*this)((inxold+m1+nx)%nx,(inyold+m2+ny)%ny)*q; w+=q;}
			}
	 	} else {
         		for (int m2 =kbmin; m2 <=kbmax; m2++){ for (int m1 =kbmin; m1 <=kbmax; m1++) {
	 			float q = kb.i0win_tab(delx - inxold-m1)*kb.i0win_tab(dely - inyold-m2);
	 			pixel += (*this)(inxold+m1,inyold+m2)*q; w+=q;}
			}
	 	}
	} else {  //  3D
		dely = restrict2(dely, ny);
		int inyold = int(Util::round(dely));
		delz = restrict2(delz, nz);
		int inzold = int(Util::round(delz));
		    //cout << inxold<<"  "<< kbc<<"  "<< nx-kbc-2<<"  "<< endl;
	 	if(inxold <= kbc || inxold >=nx-kbc-2 || inyold <= kbc || inyold >=ny-kbc-2  || inzold <= kbc || inzold >=nz-kbc-2 )  {
	 		//  loop for strips
         		for (int m3 =kbmin; m3 <=kbmax; m3++){ for (int m2 =kbmin; m2 <=kbmax; m2++){ for (int m1 =kbmin; m1 <=kbmax; m1++) {
	 			float q = kb.i0win_tab(delx - inxold-m1)*kb.i0win_tab(dely - inyold-m2)*kb.i0win_tab(delz - inzold-m3);
				//cout << "BB  "<<m1<<"  "<< m2<<"  "<< m3<<"  "<< q<<"  "<< q<<"  "<<(*this)((inxold+m1+nx)%nx,(inyold+m2+ny)%ny,(inzold+m3+nz)%nz)<< endl;
	 			pixel += (*this)((inxold+m1+nx)%nx,(inyold+m2+ny)%ny,(inzold+m3+nz)%nz)*q ;w+=q;}}
			}
	 	} else {
         		for (int m3 =kbmin; m3 <=kbmax; m3++){ for (int m2 =kbmin; m2 <=kbmax; m2++){ for (int m1 =kbmin; m1 <=kbmax; m1++) {
	 			float q = kb.i0win_tab(delx - inxold-m1)*kb.i0win_tab(dely - inyold-m2)*kb.i0win_tab(delz - inzold-m3);
				//cout << "OO  "<<m1<<"  "<< m2<<"  "<< m3<<"  "<< q<<"  "<< q<<"  "<<(*this)(inxold+m1,inyold+m2,inzold+m3)<< endl;
	 			pixel += (*this)(inxold+m1,inyold+m2,inzold+m3)*q; w+=q;}}
			}
	 	}
	}
        return pixel/w;
}


float  EMData::get_pixel_filtered(float delx, float dely, float delz, Util::sincBlackman& kb) {
//  here counting is in C style, so coordinates of the pixel delx should be [0-nx-1]

	int K     = kb.get_sB_size();
	int kbmin = -K/2;
	int kbmax = -kbmin;
	int kbc   = kbmax+1;

	float pixel = 0.0f;
	float w     = 0.0f;

	//delx = restrict2(delx, nx);	//  In this function the old location is always within the	image
	int inxold = int(Util::round(delx));
	/*if(ny<2) {  //1D
	 	if(inxold <= kbc || inxold >=nx-kbc-2 )  {
	 		//  loop for ends
         		for (int m1 =kbmin; m1 <=kbmax; m1++) {
	 			float q = kb.sBwin_tab(delx - inxold-m1);
	 			pixel += (*this)((inxold+m1+nx)%nx)*q; w+=q;
			}
	 	} else {
         		for (int m1 =kbmin; m1 <=kbmax; m1++) {
	 			float q = kb.sBwin_tab(delx - inxold-m1);
	 			pixel += (*this)(inxold+m1)*q; w+=q;
			}
	 	}

	} else */
	if(nz<2) {  
		//dely = restrict2(dely, ny);
		int inyold = int(Util::round(dely));
	 	if(inxold <= kbc || inxold >=nx-kbc-2 || inyold <= kbc || inyold >=ny-kbc-2 )  {
	 		//  loop for strips
	 		//cout<<"   FIRST  "<<inxold<<"  "<< inyold<<"  "<< kbc<<"  "<< nx-kbc-2<<"  "<< ny-kbc-2<< endl;
         	for (int m2 =kbmin; m2 <=kbmax; m2++){
				float t = kb.sBwin_tab(dely - inyold-m2);
				for (int m1 =kbmin; m1 <=kbmax; m1++) {
	 				float q = kb.sBwin_tab(delx - inxold-m1)*t;
	 				pixel += (*this)((inxold+m1+nx)%nx,(inyold+m2+ny)%ny)*q;
					w += q;
					float zzz = (*this)((inxold+m1+nx)%nx,(inyold+m2+ny)%ny);
					//if(inxold == 0 && inyold == 5)  cout<<"   OO  "<<m1<<"  "<< m2<<"  "<< (inxold+m1+nx)%nx<<"  "<< (inyold+m2+ny)%ny<<"  "<< zzz<<"  "<< q<<"  "<< pixel<< endl;
					//if(inxold == 1 && inyold == 5)  cout<<"   JJ  "<<m1<<"  "<< m2<<"  "<< (inxold+m1+nx)%nx<<"  "<< (inyold+m2+ny)%ny<<"  "<< zzz<<"  "<< q<<"  "<< pixel<< endl;
				}
			}
	 	} else {
	 		//cout<<"   SECOND  "<<inxold<<"  "<< inyold<<"  "<< kbc<<"  "<< nx-kbc-2<<"  "<< ny-kbc-2<< endl;
			for (int m2 =kbmin; m2 <=kbmax; m2++){
				float t = kb.sBwin_tab(dely - inyold-m2);
				for (int m1 =kbmin; m1 <=kbmax; m1++) {
					float q = kb.sBwin_tab(delx - inxold-m1)*t;
					pixel += (*this)(inxold+m1,inyold+m2)*q;
					w += q;
				}
			}
	 	}
	} else {  //  3D
		//std::cout<<"pixel_filtered 3D"<<std::endl;
		dely = restrict2(dely, ny);
		int inyold = int(Util::round(dely));
		delz = restrict2(delz, nz);
		int inzold = int(Util::round(delz));
		//cout << inxold<<"  "<< kbc<<"  "<< nx-kbc-2<<"  "<< endl;
	 	if(inxold <= kbc || inxold >=nx-kbc-2 || inyold <= kbc || inyold >=ny-kbc-2  || inzold <= kbc || inzold >=nz-kbc-2 )  {
	 		//  loop for strips
         		for (int m3 =kbmin; m3 <=kbmax; m3++){ for (int m2 =kbmin; m2 <=kbmax; m2++){ for (int m1 =kbmin; m1 <=kbmax; m1++) {
	 			float q = kb.sBwin_tab(delx - inxold-m1)*kb.sBwin_tab(dely - inyold-m2)*kb.sBwin_tab(delz - inzold-m3);
				//cout << "BB  "<<m1<<"  "<< m2<<"  "<< m3<<"  "<< q<<"  "<< q<<"  "<<(*this)((inxold+m1+nx)%nx,(inyold+m2+ny)%ny,(inzold+m3+nz)%nz)<< endl;
	 			pixel += (*this)((inxold+m1+nx)%nx,(inyold+m2+ny)%ny,(inzold+m3+nz)%nz)*q ;w+=q;}}
			}
	 	} else {
         		for (int m3 =kbmin; m3 <=kbmax; m3++){ for (int m2 =kbmin; m2 <=kbmax; m2++){ for (int m1 =kbmin; m1 <=kbmax; m1++) {
	 			float q = kb.sBwin_tab(delx - inxold-m1)*kb.sBwin_tab(dely - inyold-m2)*kb.sBwin_tab(delz - inzold-m3);
				//cout << "OO  "<<m1<<"  "<< m2<<"  "<< m3<<"  "<< q<<"  "<< q<<"  "<<(*this)(inxold+m1,inyold+m2,inzold+m3)<< endl;
	 			pixel += (*this)(inxold+m1,inyold+m2,inzold+m3)*q; w+=q;}}
			}
	 	}
	}
	return pixel/w;
}

// Note by Yang on 10/02/07
// get_pixel_conv7() is equivalent to get_pixel_conv_new(), however, it is written in this way such that it can be used in python directly
// By the way, get_pixel_conv_new() is a faster version of get_pixel_conv(), I have done a lot of testing and show that their results are the same.
float  EMData::get_pixel_conv7(float delx, float dely, float delz, Util::KaiserBessel& kb) {
//  here counting is in C style, so coordinates of the pixel delx should be [0-nx-1]

	float *image=(this->get_data());
	int nx = this->get_xsize();
	int ny = this->get_ysize();
	int nz = this->get_zsize();

	float result;

	result = Util::get_pixel_conv_new(nx,ny,nz,delx,dely,delz,image,kb);
	return result;
}

float EMData::getconvpt2d_kbi0(float x, float y, Util::KaiserBessel::kbi0_win win, int size) {
	const int nxhalf = nx/2;
	const int nyhalf = ny/2;
	const int bd = size/2;
	float* wxarr = new float[size];
	float* wyarr = new float[size];
	float* wx = wxarr + bd; // wx[-bd] = wxarr[0]
	float* wy = wyarr + bd;
	int ixc = int(x + 0.5f*Util::sgn(x));
	int iyc = int(y + 0.5f*Util::sgn(y));
	if (abs(ixc) > nxhalf)
		throw InvalidValueException(ixc, "getconv: X value out of range");
	if (abs(iyc) > nyhalf)
		throw InvalidValueException(ixc, "getconv: Y value out of range");
	for (int i = -bd; i <= bd; i++) {
		int iyp = iyc + i;
		wy[i] = win(y - iyp);
		int ixp = ixc + i;
		wx[i] = win(x - ixp);
	}
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(-nxhalf, -nyhalf);
	float conv = 0.f, wsum = 0.f;
	for (int iy = -bd; iy <= bd; iy++) {
		int iyp = iyc + iy;
		for (int ix = -bd; ix <= bd; ix++) {
			int ixp = ixc + ix;
			float wg = wx[ix]*wy[iy];
			conv += (*this)(ixp,iyp)*wg;
			wsum += wg;
		}
	}
	set_array_offsets(saved_offsets);
	delete [] wxarr;
	delete [] wyarr;
	//return conv/wsum;
	return conv;
}

std::complex<float> EMData::extractpoint(float nuxnew, float nuynew, Util::KaiserBessel& kb) {
	if (2 != get_ndim())
		throw ImageDimensionException("extractpoint needs a 2-D image.");
	if (!is_complex())
		throw ImageFormatException("extractpoint requires a fourier image");
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
	std::complex<float> result(0.f,0.f);
	if ((ixn >= -kbmin) && (ixn <= nhalf-1-kbmax) && (iyn >= -nhalf-kbmin) && (iyn <= nhalf-1-kbmax)) {
		// (xin,yin) not within window border from the edge
		for (int iy = iymin; iy <= iymax; iy++) {
			int iyp = iyn + iy;
			for (int ix = ixmin; ix <= ixmax; ix++) {
				int ixp = ixn + ix;
				float w = wx[ix]*wy[iy];
				std::complex<float> val = cmplx(ixp,iyp);
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
				std::complex<float> val = this->cmplx(ixt,iyt);
				if (mirror)  result += conj(val)*w;
				else         result += val*w;
			}
		}
	}
	if (flip)  result = conj(result)/wsum;
	else       result /= wsum;
	delete [] wx0;
	delete [] wy0;
	return result;
}

EMData* EMData::extractline(Util::KaiserBessel& kb, float nuxnew, float nuynew)
{
	if (!is_complex())
		throw ImageFormatException("extractline requires a fourier image");
	if (nx%2 != 0)
		throw ImageDimensionException("extractline requires nx to be even");
	int nxreal = nx - 2;
	if (nxreal != ny)
		throw ImageDimensionException("extractline requires ny == nx");
	// build complex result image
	EMData* res = new EMData();
	res->set_size(nx,1,1);
	res->to_zero();
	res->set_complex(true);
	res->set_fftodd(false);
	res->set_fftpad(true);
	res->set_ri(true);
	// Array offsets: (0..nhalf,-nhalf..nhalf-1)
	int n = nxreal;
	int nhalf = n/2;
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,-nhalf,-nhalf);

	// set up some temporary weighting arrays
	int kbsize = kb.get_window_size();
	int kbmin = -kbsize/2;
	int kbmax = -kbmin;
	float* wy0 = new float[kbmax - kbmin + 1];
	float* wy = wy0 - kbmin; // wy[kbmin:kbmax]
	float* wx0 = new float[kbmax - kbmin + 1];
	float* wx = wx0 - kbmin;

	int   count = 0;
	float wsum = 0.f;
	bool  flip = (nuxnew < 0.f);

	for (int jx = 0; jx <= nhalf; jx++) {
		float xnew = jx*nuxnew, ynew = jx*nuynew;
		count++;
		std::complex<float> btq(0.f,0.f);
		if (flip) {
			xnew = -xnew;
			ynew = -ynew;
		}
		int ixn = int(Util::round(xnew));
		int iyn = int(Util::round(ynew));
		// populate weight arrays
		for (int i=kbmin; i <= kbmax; i++) {
			int iyp = iyn + i;
			wy[i] = kb.i0win_tab(ynew - iyp);
			int ixp = ixn + i;
			wx[i] = kb.i0win_tab(xnew - ixp);
		}
		// restrict weight arrays to non-zero elements

		int lnby = 0;
		for (int iy = kbmin; iy <= -1; iy++) {
			if (wy[iy] != 0.f) {
				lnby = iy;
				break;
			}
		}
		int lney = 0;
		for (int iy = kbmax; iy >= 1; iy--) {
			if (wy[iy] != 0.f) {
				lney = iy;
				break;
			}
		}
		int lnbx = 0;
		for (int ix = kbmin; ix <= -1; ix++) {
			if (wx[ix] != 0.f) {
				lnbx = ix;
				break;
			}
		}
		int lnex = 0;
		for (int ix = kbmax; ix >= 1; ix--) {
			if (wx[ix] != 0.f) {
				lnex = ix;
				break;
			}
		}
		if (ixn >= -kbmin && ixn <= nhalf-1-kbmax
				&& iyn >= -nhalf-kbmin && iyn <= nhalf-1-kbmax) {
			// interior points
			for (int ly=lnby; ly<=lney; ly++) {
				int iyp = iyn + ly;
				for (int lx=lnbx; lx<=lnex; lx++) {
					int ixp = ixn + lx;
					float wg = wx[lx]*wy[ly];
					btq += cmplx(ixp,iyp)*wg;
					wsum += wg;
				}
			}
		} else {
			// points "sticking out"
			for (int ly=lnby; ly<=lney; ly++) {
				int iyp = iyn + ly;
				for (int lx=lnbx; lx<=lnex; lx++) {
					int ixp = ixn + lx;
					float wg = wx[lx]*wy[ly];
					bool mirror = false;
					int ixt(ixp), iyt(iyp);
					if (ixt > nhalf || ixt < -nhalf) {
						ixt = Util::sgn(ixt)*(n - abs(ixt));
						iyt = -iyt;
						mirror = !mirror;
					}
					if (iyt >= nhalf || iyt < -nhalf) {
						if (ixt != 0) {
							ixt = -ixt;
							iyt = Util::sgn(iyt)*(n - abs(iyt));
							mirror = !mirror;
						} else {
							iyt -= n*Util::sgn(iyt);
						}
					}
					if (ixt < 0) {
						ixt = -ixt;
						iyt = -iyt;
						mirror = !mirror;
					}
					if (iyt == nhalf) iyt = -nhalf;
					if (mirror) btq += conj(cmplx(ixt,iyt))*wg;
					else        btq += cmplx(ixt,iyt)*wg;
					wsum += wg;
				}
			}
		}
		if (flip) res->cmplx(jx) = conj(btq);
		else      res->cmplx(jx) = btq;
	}
	for (int jx = 0; jx <= nhalf; jx++)  res->cmplx(jx) *= count/wsum;

	delete[] wx0; delete[] wy0;
	set_array_offsets(saved_offsets);
	res->set_array_offsets(0,0,0);
	return res;
}


/** Helper function for EMData::fft_shuffle, below */
inline void swapx(float* a, float* b, float* temp, size_t nbytes) {
	memcpy(temp, a, nbytes);
	memcpy(a, b, nbytes);
	memcpy(b, temp, nbytes);
}

void EMData::fft_shuffle() {
	//  swap Fourier image columns
	//  such that center is moved from (0,0,0) to (0, ny/2, nz/2)
	if (!is_complex())
		throw ImageFormatException("fft_shuffle requires a fourier image");
	vector<int> offsets = get_array_offsets();
	set_array_offsets(); // clear offsets before shuffling
	EMData& self = *this;
	int nyhalf = ny/2;
	int nzhalf = nz/2;
	int nbytes = nx*sizeof(float);
	float* temp = new float[nx];
	for (int iz=0; iz < nz; iz++)
		for (int iy=0; iy < nyhalf; iy++)
			swapx(&self(0,iy,iz),&self(0,iy+nyhalf,iz),temp,nbytes);
	if (nz > 1) {
		for (int iy=0; iy < ny; iy++)
			for (int iz=0; iz < nzhalf; iz++)
				swapx(&self(0,iy,iz),&self(0,iy,iz+nzhalf),temp,nbytes);
	}
	set_shuffled(!is_shuffled()); // toggle
	set_array_offsets(offsets); // reset offsets
	update();
	delete[] temp;
}

void EMData::pad_corner(float *pad_image) {
	size_t nbytes = nx*sizeof(float);
	for (int iy=0; iy<ny; iy++)
		memcpy(&(*this)(0,iy), pad_image+3+(iy+3)*nx, nbytes);
}

void EMData::shuffle_pad_corner(float *pad_image) {
	int nyhalf = ny/2;
	size_t nbytes = nx*sizeof(float);
	for (int iy = 0; iy < nyhalf; iy++)
		memcpy(&(*this)(0,iy), pad_image+6+(iy+nyhalf+3)*nx, nbytes);
	for (int iy = nyhalf; iy < ny; iy++)
		memcpy(&(*this)(0,iy), pad_image+6+(iy-nyhalf+3)*nx, nbytes);
}

EMData* EMData::ft2polar(int ring_length, int nb, int ne) {
	if (2 != get_ndim())
		throw ImageDimensionException("ft2polar requires a 2-D image.");
	EMData *cimage = NULL;
	if ( is_complex() )  ImageFormatException("ft2polar requires a real image.");
	else {
		cimage = this->copy();
		cimage->set_attr("npad",1);
		cimage->div_sinc(1);
		cimage = cimage->norm_pad(false, 1);
		cimage->do_fft_inplace();
		cimage->center_origin_fft();
	}

	int nx = cimage->get_xsize();
	int ny = cimage->get_ysize();
	int nc = ny/2;
	int lcirc = ne-nb+1;
	if( ne > nc-2) throw InvalidValueException(ne, "Maximum radius too large.");
	//printf("  dims   nx ny   %d   %d   %d   %d\n",nx,ny,ring_length,lcirc);

	EMData* rings = new EMData(2*ring_length-2, lcirc, 1, false);
	rings->set_fftpad(0);
	//int lx = rings->get_xsize();
	//int ly = rings->get_ysize();
	//printf("  dirns  nx ny   %d   %d   %d   %d\n",lx,ly,ring_length,lcirc);
	//float* xim  = cimage->get_data();
	//for(unsigned int inr = 0; inr <= 2*nx; inr++) printf("ft2polar   %d   %f\n",inr,xim[inr]);

	float dfi;
	dfi = TWOPI / ring_length;
//	Table for sin & cos
	vector<float> vsin(ring_length/2);
	vector<float> vcos(ring_length/2);
	for (int x = 0; x < ring_length/2; x++) {
		float ang = static_cast<float>(x * dfi);
		vsin[x] = sin(ang);
		vcos[x] = cos(ang);
		//printf("trigtab   %d      %f      %f  %f\n",x,ang,vsin[x],vcos[x]);
	}
	for (unsigned int inr = nb; inr <= ne; inr++) {
		for (unsigned int it = 0; it < ring_length/2; it++) {
			float xnew = vsin[it] * inr;
			float ynew = vcos[it] * inr;
			complex<float> v1 = (complex<float>)cimage->get_complex_at_interp(xnew,ynew);
			//printf("   %d   %d       %f  %f      (%f , %f)\n",it,inr,xnew,ynew, std::real(v1), std::imag(v1));
			rings->set_value_at(2*it,inr-nb,std::real(v1));
			rings->set_value_at(2*it+1,inr-nb,std::imag(v1));
			rings->set_value_at(2*it+ring_length,inr-nb,std::real(v1));
			rings->set_value_at(2*it+1+ring_length,inr-nb,-std::imag(v1));
		}
	}

	delete cimage;
	cimage = 0;

	EXITFUNC;
	return rings;
}


EMData* EMData::ft2polargrid(int ring_length, int nb, int ne, Util::KaiserBessel& kb) {
	if(get_ndim() != 2) throw ImageDimensionException("ft2polargrid requires a 2-D image.");
	if(!is_complex()) throw ImageFormatException("ft2polargrid requires a fourier image");

	int nxreal = nx - 2 + int(is_fftodd());
	if(nxreal != ny) throw ImageDimensionException("ft2polargrid requires ny == nx(real)");
	if(nxreal%2 != 0) throw ImageDimensionException("ft2polargrid needs an even image.");

	if(!is_shuffled()) fft_shuffle();
	set_array_offsets(0,-ny/2);

	int nc = ny/2;
	int lcirc = ne-nb+1;
	if( ne > nc-2) throw InvalidValueException(ne, "Maximum radius too large.");
	//printf("  dims   nx ny   %d   %d   %d   %d\n",nx,ny,ring_length,lcirc);

	EMData* rings = new EMData(2*ring_length-2, lcirc, 1, false);
	rings->set_fftpad(0);
	float dfi;
	dfi = TWOPI / ring_length;
	//	Table for sin & cos
	vector<float> vsin(ring_length/2);
	vector<float> vcos(ring_length/2);
	for (int x = 0; x < ring_length/2; x++) {
		float ang = static_cast<float>(x * dfi);
		vsin[x] = sin(ang);
		vcos[x] = cos(ang);
		//printf("trigtab   %d      %f      %f  %f\n",x,ang,vsin[x],vcos[x]);
	}

	for (unsigned int inr = nb; inr <= ne; inr++) {
		for (unsigned int it = 0; it < ring_length/2; it++) {
			float nuxold = vsin[it] * 2*inr;
			float nuyold = vcos[it] * 2*inr;
		//if(it ==0 && inr==nb) {nuxold=0.0f; nuyold=0.0f;}
			complex<float> v1 = Util::extractpoint2(nx, ny, nuxold, nuyold, this, kb);
			rings->cmplx(it,inr-nb) = v1;
			rings->cmplx(it+ring_length/2,inr-nb) = std::conj(v1);
			//printf("   %d   %d       %f  %f      (%f , %f)\n",2*it,inr,nuxold,nuyold, std::real(v1), std::imag(v1));
			
		}
	}

	rings->update();
	set_array_offsets();
	fft_shuffle(); // reset input to an unshuffled complex image
	return rings;
}

EMData* EMData::fourier_rotate_shift2d(float ang, float sx, float sy, int npad) {
	if(get_ndim() != 2) throw ImageDimensionException("fourier_rotate_shift2d requires a 2-D image.");
	EMData *cimage = NULL;
	int nxreal, nxhalf, nyhalf, nyf;
	int ix,iy,np;
	if( is_complex() )  {
		cimage = this;
		nxreal = ny/npad;
		nxhalf = nxreal/2;
		nyhalf = ny/2/npad;
		nyf = ny/npad;
	} else {
		nxreal = nx;
		nxhalf = nxreal/2;
		nyhalf = ny/2;
		nyf = ny;
		EMData *padimage = Util::pad(this, 2*nx, 2*ny, 1, 0,0,0, "0.0");
		padimage->set_attr("npad",npad);
		padimage->div_sinc(1);
		cimage = padimage->norm_pad(false, 1);
		delete padimage; padimage=0;
		cimage->do_fft_inplace();
		cimage->center_origin_fft();
		cimage->fft_shuffle();
		cimage->set_attr("npad",npad);
	}

	float cir = (float)((nxhalf-1)*(nxhalf-1));

	EMData* result = new EMData(nxreal, nxreal, 1, false);
	result->to_zero();

	float fang = ang*(float)deg_rad;
	float cang = cos(fang);
	float sang = sin(fang);
	float temp = -2.0f*M_PI/nxreal;
	for (int iy = -nyhalf; iy < nyhalf; iy++) {
		int it = iy + nyhalf;
		float ycang = iy*cang;
		float ysang = iy*sang;
		for (int ix = 0; ix <= nxhalf; ix++) {
			float nuxold = ix*cang - ysang;
			float nuyold = ix*sang + ycang;
			if(nuxold*nuxold+nuyold*nuyold<cir) {
				nuyold += nyhalf;
				if(nuyold>nyhalf)  nuyold -= nyf;
				nuxold *= npad;
				nuyold *= npad;
				complex<float> v1 = (complex<float>)cimage->get_complex_at_interp(nuxold, nuyold);
				float phase_ang = temp*(sx*ix+sy*iy);
				result->cmplx(ix,it) = v1*complex<float>(cos(phase_ang), sin(phase_ang));
				//complex<float> v0 = cimage->cmplx(ix,it);
				//complex<float> v2 = result->cmplx(ix,it);
				//printf("indexes   %d      %d     %d         %f   %f           %f   %f        %f   %f       %f   %f\n",ix,iy,it,std::real(v0),std::imag(v0),std::real(v1),std::imag(v1),nuxold,nuyold,std::real(v2),std::imag(v2));
			}
		}
	}

	result->fft_shuffle(); // reset to an un shuffled result
	result->set_shuffled(false);
	result->center_origin_fft();
	result->set_fftpad(true);
	result->set_attr("npad", 1);
	result->set_fftodd(false);
	if( !is_complex() )  {
		result->do_ift_inplace();
		result->depad();
		delete cimage; cimage = 0;
	}
	result->update();
	return result;
}

// We tried to pad the Fourier image to reduce the stick out points, howover it is not very efficient.
/*
EMData* EMData::fouriergridrot2d(float ang, float scale, Util::KaiserBessel& kb) {
	if (2 != get_ndim())
		throw ImageDimensionException("fouriergridrot2d needs a 2-D image.");
	if (!is_complex())
		throw ImageFormatException("fouriergridrot2d requires a fourier image");
	int nxreal = nx - 2 + int(is_fftodd());
	if (nxreal != ny)
		throw ImageDimensionException("fouriergridrot2d requires ny == nx(real)");
	if (0 != nxreal%2)
		throw ImageDimensionException("fouriergridrot2d needs an even image.");
	if (scale == 0.0f) scale = 1.0f;
	int nxhalf = nxreal/2;
	int nyhalf = ny/2;

	EMData *pad_this = new EMData();
	pad_this->set_size(nx+12, ny+6);
	//pad_this->to_zero();
	float* pad_image = pad_this-> get_data();

	if (!is_shuffled()) {
		shuffle_pad_corner(pad_image);
	} else {
		pad_corner(pad_image);
	}
	pad_this -> set_array_offsets(-6, -nyhalf-3);

	EMData* result = copy_head();
	set_array_offsets(0,-nyhalf);
	result->set_array_offsets(0,-nyhalf);

	ang = ang*deg_rad;
	float cang = cos(ang);
	float sang = sin(ang);
	for (int iy = -nyhalf; iy < nyhalf; iy++) {
		float ycang = iy*cang;
		float ysang = iy*sang;
		for (int ix = 0; ix <= nxhalf; ix++) {
			float nuxold = (ix*cang - ysang)*scale;
			float nuyold = (ix*sang + ycang)*scale;
			result->cmplx(ix,iy) = Util::extractpoint2(nx, ny, nuxold, nuyold, pad_this, kb);
		}
	}
	result->set_array_offsets();
	result->fft_shuffle(); // reset to an unshuffled result
	result->update();
	set_array_offsets();
	fft_shuffle(); // reset to an unshuffled complex image
	return result;
}*/


EMData* EMData::fouriergridrot2d(float ang, float scale, Util::KaiserBessel& kb) {
	if (2 != get_ndim())
		throw ImageDimensionException("fouriergridrot2d needs a 2-D image.");
	if (!is_complex())
		throw ImageFormatException("fouriergridrot2d requires a fourier image");
	int nxreal = nx - 2 + int(is_fftodd());
	if (nxreal != ny)
		throw ImageDimensionException("fouriergridrot2d requires ny == nx(real)");
	if (0 != nxreal%2)
		throw ImageDimensionException("fouriergridrot2d needs an even image.");
	if (scale == 0.0f) scale = 1.0f;
	int nxhalf = nxreal/2;
	int nyhalf = ny/2;
	float cir = (float)((nxhalf-1)*(nxhalf-1));

	if (!is_shuffled()) fft_shuffle();

	EMData* result = copy_head();
	set_array_offsets(0,-nyhalf);
	result->set_array_offsets(0,-nyhalf);

	ang = ang*(float)deg_rad;
	float cang = cos(ang);
	float sang = sin(ang);
	for (int iy = -nyhalf; iy < nyhalf; iy++) {
		float ycang = iy*cang;
		float ysang = iy*sang;
		for (int ix = 0; ix <= nxhalf; ix++) {
			float nuxold = (ix*cang - ysang)*scale;
			float nuyold = (ix*sang + ycang)*scale;
			if(nuxold*nuxold+nuyold*nuyold<cir) result->cmplx(ix,iy) = Util::extractpoint2(nx, ny, nuxold, nuyold, this, kb);
			//result->cmplx(ix,iy) = extractpoint(nuxold, nuyold, kb);
		}
	}
	result->set_array_offsets();
	result->fft_shuffle(); // reset to an unshuffled result
	result->update();
	set_array_offsets();
	fft_shuffle(); // reset input to an unshuffled complex image
	return result;
}

EMData* EMData::fouriergridrot_shift2d(float ang, float sx, float sy, Util::KaiserBessel& kb) {
	if (2 != get_ndim())
		throw ImageDimensionException("fouriergridrot_shift2d needs a 2-D image.");
	if (!is_complex())
		throw ImageFormatException("fouriergridrot_shift2d requires a fourier image");
	int nxreal = nx - 2 + int(is_fftodd());
	if (nxreal != ny)
		throw ImageDimensionException("fouriergridrot_shift2d requires ny == nx(real)");
	if (0 != nxreal%2)
		throw ImageDimensionException("fouriergridrot_shift2d needs an even image.");
	int nxhalf = nxreal/2;
	int nyhalf = ny/2;

	if (!is_shuffled()) fft_shuffle();

	EMData* result = copy_head();
	set_array_offsets(0, -nyhalf);
	result->set_array_offsets(0, -nyhalf);

	ang = ang*(float)deg_rad;
	float cang = cos(ang);
	float sang = sin(ang);
	float temp = -2.0f*M_PI/nxreal;
	for (int iy = -nyhalf; iy < nyhalf; iy++) {
		float ycang = iy*cang;
		float ysang = iy*sang;
		for (int ix = 0; ix <= nxhalf; ix++) {
			float nuxold = ix*cang - ysang;
			float nuyold = ix*sang + ycang;
			result->cmplx(ix,iy) = Util::extractpoint2(nx, ny, nuxold, nuyold, this, kb);
			//result->cmplx(ix,iy) = extractpoint(nuxold, nuyold, kb);
			float phase_ang = temp*(sx*ix+sy*iy);
			result->cmplx(ix,iy) *= complex<float>(cos(phase_ang), sin(phase_ang));
		}
	}
	result->set_array_offsets();
	result->fft_shuffle(); // reset to an unshuffled result
	result->update();
	set_array_offsets();
	fft_shuffle(); // reset to an unshuffled complex image
	return result;
}

void EMData::divkbsinh(const Util::KaiserBessel& kb) {
	
	if (is_complex())
		throw ImageFormatException("divkbsinh requires a real image.");
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
	// Note that the following loops will work for 1-, 2-, and 3-D
	// images, since the "extra" weights will be 1.0.  (For example,
	// for a 2-d image iz=0, nz=1, so iz-nz/2 = 0 - 1/2 = 0, since
	// the division is an integer division.)
	for (int iz=0; iz < nz; iz++) {
		float wz = kb.sinhwin(static_cast<float>(iz-nz/2));
		for (int iy=0; iy < ny; iy++) {
			float wy = kb.sinhwin(static_cast<float>(iy-ny/2));
			for (int ix=0; ix < nx; ix++) {
				float wx = kb.sinhwin(static_cast<float>(ix-nx/2));
				float w = wx*wy*wz;
				(*this)(ix,iy,iz) /= w;
			}
		}
	}
	set_array_offsets(saved_offsets);
}

void EMData::divkbsinh_rect(const Util::KaiserBessel& kbx, const Util::KaiserBessel& kby, const Util::KaiserBessel& kbz) {

	if (is_complex())
		throw ImageFormatException("divkbsinh requires a real image.");
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
	// Note that the following loops will work for 1-, 2-, and 3-D
	// images, since the "extra" weights will be 1.0.  (For example,
	// for a 2-d image iz=0, nz=1, so iz-nz/2 = 0 - 1/2 = 0, since
	// the division is an integer division.)
	for (int iz=0; iz < nz; iz++) {
		float wz = kbz.sinhwin(static_cast<float>(iz-nz/2));
		for (int iy=0; iy < ny; iy++) {
			float wy = kby.sinhwin(static_cast<float>(iy-ny/2));
			for (int ix=0; ix < nx; ix++) {
				float wx = kbx.sinhwin(static_cast<float>(ix-nx/2));
				float w = wx*wy*wz;
				(*this)(ix,iy,iz) /= w;
			}
		}
	}
	
	set_array_offsets(saved_offsets);
}

/* OBSOLETED  PAP
Dict EMData::masked_stats(const EMData* mask) {
	if (is_complex())
		throw ImageFormatException(
				"Complex images not supported by EMData::masked_stats");
	float* ptr = get_data();
	float* mptr = mask->get_data();
	long double sum1 = 0.L;
	long double sum2 = 0.L;
	long nmask = 0L;
	for (long i = 0; i < nx*ny*nz; i++,ptr++,mptr++) {
		if (*mptr > 0.5f) {
			nmask++;
			sum1 += *ptr;
			sum2 += (*ptr)*(*ptr);
		}
	}
	float avg = static_cast<float>(sum1/nmask);
	float sig2 = static_cast<float>(sum2/nmask - avg*avg);
	float sig = sqrt(sig2);
	Dict mydict;
	mydict["avg"] = avg; mydict["sigma"] = sig; mydict["nmask"] = int(nmask);
	return mydict;
}
*/

EMData* EMData::extract_plane(const Transform& tf, Util::KaiserBessel& kb) {
	if (!is_complex())
		throw ImageFormatException("extractplane requires a complex image");
	if (nx%2 != 0)
		throw ImageDimensionException("extractplane requires nx to be even");
	int nxreal = nx - 2;
	if (nxreal != ny || nxreal != nz)
		throw ImageDimensionException("extractplane requires ny == nx == nz");
	// build complex result image
	EMData* res = new EMData();
	res->set_size(nx,ny,1);
	res->to_zero();
	res->set_complex(true);
	res->set_fftodd(false);
	res->set_fftpad(true);
	res->set_ri(true);
	// Array offsets: (0..nhalf,-nhalf..nhalf-1,-nhalf..nhalf-1)
	int n = nxreal;
	int nhalf = n/2;
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,-nhalf,-nhalf);
	res->set_array_offsets(0,-nhalf,0);
	// set up some temporary weighting arrays
	int kbsize =  kb.get_window_size();
	int kbmin  = -kbsize/2;
	int kbmax  = -kbmin;
	float* wy0 = new float[kbmax - kbmin + 1];
	float* wy  = wy0 - kbmin; // wy[kbmin:kbmax]
	float* wx0 = new float[kbmax - kbmin + 1];
	float* wx  = wx0 - kbmin;
	float* wz0 = new float[kbmax - kbmin + 1];
	float* wz  = wz0 - kbmin;
	float rim = nhalf*float(nhalf);
	int count = 0;
	float wsum = 0.f;
	Transform tftrans = tf; // need transpose of tf here for consistency
	tftrans.invert();      // with spider
	for (int jy = -nhalf; jy < nhalf; jy++) 
	{
		for (int jx = 0; jx <= nhalf; jx++) 
		{
			Vec3f nucur((float)jx, (float)jy, 0.f);
			Vec3f nunew = tftrans*nucur;
			float xnew = nunew[0], ynew = nunew[1], znew = nunew[2];
			if (xnew*xnew+ynew*ynew+znew*znew <= rim) 
			{
				count++;
				std::complex<float> btq(0.f,0.f);
				bool flip = false;
				if (xnew < 0.f) {
					flip = true;
					xnew = -xnew;
					ynew = -ynew;
					znew = -znew;
				}
				int ixn = int(Util::round(xnew));
				int iyn = int(Util::round(ynew));
				int izn = int(Util::round(znew));
				// populate weight arrays
				for (int i=kbmin; i <= kbmax; i++) {
					int izp = izn + i;
					wz[i] = kb.i0win_tab(znew - izp);
					int iyp = iyn + i;
					wy[i] = kb.i0win_tab(ynew - iyp);
					int ixp = ixn + i;
					wx[i] = kb.i0win_tab(xnew - ixp);

				}
				// restrict weight arrays to non-zero elements
				int lnbz = 0;
				for (int iz = kbmin; iz <= -1; iz++) {
					if (wz[iz] != 0.f) {
						lnbz = iz;
						break;
					}
				}
				int lnez = 0;
				for (int iz = kbmax; iz >= 1; iz--) {
					if (wz[iz] != 0.f) {
						lnez = iz;
						break;
					}
				}
				int lnby = 0;
				for (int iy = kbmin; iy <= -1; iy++) {
					if (wy[iy] != 0.f) {
						lnby = iy;
						break;
					}
				}
				int lney = 0;
				for (int iy = kbmax; iy >= 1; iy--) {
					if (wy[iy] != 0.f) {
						lney = iy;
						break;
					}
				}
				int lnbx = 0;
				for (int ix = kbmin; ix <= -1; ix++) {
					if (wx[ix] != 0.f) {
						lnbx = ix;
						break;
					}
				}
				int lnex = 0;
				for (int ix = kbmax; ix >= 1; ix--) {
					if (wx[ix] != 0.f) {
						lnex = ix;
						break;
					}
				}
				if    (ixn >= -kbmin      && ixn <= nhalf-1-kbmax
				   && iyn >= -nhalf-kbmin && iyn <= nhalf-1-kbmax
				   && izn >= -nhalf-kbmin && izn <= nhalf-1-kbmax) {
					// interior points
					for (int lz = lnbz; lz <= lnez; lz++) {
						int izp = izn + lz;
						for (int ly=lnby; ly<=lney; ly++) {
							int iyp = iyn + ly;
							float ty = wz[lz]*wy[ly];
							for (int lx=lnbx; lx<=lnex; lx++) {
								int ixp = ixn + lx;
								float wg = wx[lx]*ty;
								btq += cmplx(ixp,iyp,izp)*wg;
								wsum += wg;
							}
						}
					}
				} else {
					// points "sticking out"
					for (int lz = lnbz; lz <= lnez; lz++) {
						int izp = izn + lz;
						for (int ly=lnby; ly<=lney; ly++) {
							int iyp = iyn + ly;
							float ty = wz[lz]*wy[ly];
							for (int lx=lnbx; lx<=lnex; lx++) {
								int ixp = ixn + lx;
								float wg = wx[lx]*ty;
								bool mirror = false;
								int ixt(ixp), iyt(iyp), izt(izp);
								if (ixt > nhalf || ixt < -nhalf) {
									ixt = Util::sgn(ixt)
										  *(n - abs(ixt));
									iyt = -iyt;
									izt = -izt;
									mirror = !mirror;
								}
								if (iyt >= nhalf || iyt < -nhalf) {
									if (ixt != 0) {
										ixt = -ixt;
										iyt = Util::sgn(iyt)
											  *(n - abs(iyt));
										izt = -izt;
										mirror = !mirror;
									} else {
										iyt -= n*Util::sgn(iyt);
									}
								}
								if (izt >= nhalf || izt < -nhalf) {
									if (ixt != 0) {
										ixt = -ixt;
										iyt = -iyt;
										izt = Util::sgn(izt)
											  *(n - abs(izt));
										mirror = !mirror;
									} else {
										izt -= Util::sgn(izt)*n;
									}
								}
								if (ixt < 0) {
									ixt = -ixt;
									iyt = -iyt;
									izt = -izt;
									mirror = !mirror;
								}
								if (iyt == nhalf) iyt = -nhalf;
								if (izt == nhalf) izt = -nhalf;
								if (mirror)   btq += conj(cmplx(ixt,iyt,izt))*wg;
								else          btq += cmplx(ixt,iyt,izt)*wg;
								wsum += wg;
							}
						}
					}
				}
				if (flip)  res->cmplx(jx,jy) = conj(btq);
				else       res->cmplx(jx,jy) = btq;
			}
		}
	}
	for (int jy = -nhalf; jy < nhalf; jy++)
		for (int jx = 0; jx <= nhalf; jx++)
			res->cmplx(jx,jy) *= count/wsum;
	delete[] wx0; delete[] wy0; delete[] wz0;
	set_array_offsets(saved_offsets);
	res->set_array_offsets(0,0,0);
	res->set_shuffled(true);
	return res;
}



EMData* EMData::extract_plane_rect(const Transform& tf, Util::KaiserBessel& kbx,Util::KaiserBessel& kby, Util::KaiserBessel& kbz) {
	

	if (!is_complex())
		throw ImageFormatException("extractplane requires a complex image");
	if (nx%2 != 0)
		throw ImageDimensionException("extractplane requires nx to be even");

	int nxfromxyz = max( max(nx-2,ny), nz) + 2;
	//int nxfromz = nz+2;
	//int nxcircal = nxfromz - 2;
	int nxcircal = nxfromxyz - 2;
	EMData* res = new EMData();
	//res->set_size(nxfromz,nz,1);
	res->set_size(nxfromxyz,nxcircal,1);
	res->to_zero();
	res->set_complex(true);
	res->set_fftodd(false);
	res->set_fftpad(true);
	res->set_ri(true);
	// Array offsets: (0..nhalf,-nhalf..nhalf-1,-nhalf..nhalf-1)
	int n = nxcircal;
	int nhalf = n/2;
	int nxhalf = (nx-2)/2;
	int nyhalf = ny/2;
	int nzhalf = nz/2;
	
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0, -nyhalf, -nzhalf);
	res->set_array_offsets(0,-nhalf,0);
	// set up some temporary weighting arrays
	int kbxsize =  kbx.get_window_size();
	int kbxmin  = -kbxsize/2;
	int kbxmax  = -kbxmin;

	int kbysize =  kby.get_window_size();
	int kbymin  = -kbysize/2;
	int kbymax  = -kbymin;

	int kbzsize =  kbz.get_window_size();
	int kbzmin  = -kbzsize/2;
	int kbzmax  = -kbzmin;

	//std::cout<<"kb size x,y,z=="<<kbxsize<<" "<<kbysize<<" "<<kbzsize<<std::endl;
	float* wy0 = new float[kbymax - kbymin + 1];
	float* wy  = wy0 - kbymin; // wy[kbmin:kbmax]
	float* wx0 = new float[kbxmax - kbxmin + 1];
	float* wx  = wx0 - kbxmin;
	float* wz0 = new float[kbzmax - kbzmin + 1];
	float* wz  = wz0 - kbzmin;
	float rim = nhalf*float(nhalf);
	int count = 0;
	float wsum = 0.f;
	Transform tftrans = tf; // need transpose of tf here for consistency
	tftrans.invert();      // with spider
	float xratio=float(nx-2)/float(nxcircal);
	float yratio=float(ny)/float(nxcircal);
	float zratio=float(nz)/float(nxcircal);
	//std::cout<<"xratio,yratio=="<<xratio<<" "<<yratio<<std::endl;
	for (int jy = -nhalf; jy < nhalf; jy++) 
	{
		for (int jx = 0; jx <= nhalf; jx++) 
		{
			Vec3f nucur((float)jx, (float)jy, 0.f);
			Vec3f nunew = tftrans*nucur;
			float xnew = nunew[0]*xratio, ynew = nunew[1]*yratio, znew = nunew[2]*zratio;
			
			if (nunew[0]*nunew[0]+nunew[1]*nunew[1]+nunew[2]*nunew[2] <= rim)
			{
				count++;
				std::complex<float> btq(0.f,0.f);
				bool flip = false;
				if (xnew < 0.f) {
					flip = true;
					xnew = -xnew;
					ynew = -ynew;
					znew = -znew;
				}
				int ixn = int(Util::round(xnew));
				int iyn = int(Util::round(ynew));
				int izn = int(Util::round(znew));
				// populate weight arrays
				for (int i=kbzmin; i <= kbzmax; i++) {
					int izp = izn + i;
					wz[i] = kbz.i0win_tab(znew - izp);
					}
				for (int i=kbymin; i <= kbymax; i++) {
					int iyp = iyn + i;
					wy[i] = kby.i0win_tab(ynew - iyp);
					}
				for (int i=kbxmin; i <= kbxmax; i++) {
					int ixp = ixn + i;
					wx[i] = kbx.i0win_tab(xnew - ixp);
					}
		

				
				// restrict weight arrays to non-zero elements
				int lnbz = 0;
				for (int iz = kbzmin; iz <= -1; iz++) {
					if (wz[iz] != 0.f) {
						lnbz = iz;
						break;
					}
				}
				int lnez = 0;
				for (int iz = kbzmax; iz >= 1; iz--) {
					if (wz[iz] != 0.f) {
						lnez = iz;
						break;
					}
				}
				int lnby = 0;
				for (int iy = kbymin; iy <= -1; iy++) {
					if (wy[iy] != 0.f) {
						lnby = iy;
						break;
					}
				}
				int lney = 0;
				for (int iy = kbymax; iy >= 1; iy--) {
					if (wy[iy] != 0.f) {
						lney = iy;
						break;
					}
				}
				int lnbx = 0;
				for (int ix = kbxmin; ix <= -1; ix++) {
					if (wx[ix] != 0.f) {
						lnbx = ix;
						break;
					}
				}
				int lnex = 0;
				for (int ix = kbxmax; ix >= 1; ix--) {
					if (wx[ix] != 0.f) {
						lnex = ix;
						break;
					}
				}
				if    (ixn >= -kbxmin      && ixn <= nxhalf-1-kbxmax
				   && iyn >= -nyhalf-kbymin && iyn <= nyhalf-1-kbymax
				   && izn >= -nzhalf-kbzmin && izn <= nzhalf-1-kbzmax) {
					// interior points
					for (int lz = lnbz; lz <= lnez; lz++) {
						int izp = izn + lz;
						for (int ly=lnby; ly<=lney; ly++) {
							int iyp = iyn + ly;
							float ty = wz[lz]*wy[ly];
							for (int lx=lnbx; lx<=lnex; lx++) {
								int ixp = ixn + lx;
								float wg = wx[lx]*ty;
								btq += cmplx(ixp,iyp,izp)*wg;
								wsum += wg;
							}
						}
					}
				} 
				else {
					// points "sticking out"
					for (int lz = lnbz; lz <= lnez; lz++) {
						int izp = izn + lz;
						for (int ly=lnby; ly<=lney; ly++) {
							int iyp = iyn + ly;
							float ty = wz[lz]*wy[ly];
							for (int lx=lnbx; lx<=lnex; lx++) {
								int ixp = ixn + lx;
								float wg = wx[lx]*ty;
								bool mirror = false;
								int ixt(ixp), iyt(iyp), izt(izp);
								if (ixt > nxhalf || ixt < -nxhalf) {
									ixt = Util::sgn(ixt)
										  *(nx-2-abs(ixt));
									iyt = -iyt;
									izt = -izt;
									mirror = !mirror;
								}
								if (iyt >= nyhalf || iyt < -nyhalf) {
									if (ixt != 0) {
										ixt = -ixt;
										iyt = Util::sgn(iyt)
											  *(ny - abs(iyt));
										izt = -izt;
										mirror = !mirror;
									} else {
										iyt -= ny*Util::sgn(iyt);
									}
								}
								if (izt >= nzhalf || izt < -nzhalf) {
									if (ixt != 0) {
										ixt = -ixt;
										iyt = -iyt;
										izt = Util::sgn(izt)
											  *(nz - abs(izt));
										mirror = !mirror;
									} else {
										izt -= Util::sgn(izt)*nz;
									}
								}
								if (ixt < 0) {
									ixt = -ixt;
									iyt = -iyt;
									izt = -izt;
									mirror = !mirror;
								}
								if (iyt == nyhalf) iyt = -nyhalf;
								if (izt == nzhalf) izt = -nzhalf;
								if (mirror)   btq += conj(cmplx(ixt,iyt,izt))*wg;
								else          btq += cmplx(ixt,iyt,izt)*wg;
								wsum += wg;
							}
						}
					}
				}
				if (flip)  res->cmplx(jx,jy) = conj(btq);
				else       res->cmplx(jx,jy) = btq;
			}
		}
	}
	for (int jy = -nhalf; jy < nhalf; jy++)
		for (int jx = 0; jx <= nhalf; jx++)
			res->cmplx(jx,jy) *= count/wsum;
	delete[] wx0; delete[] wy0; delete[] wz0;
	set_array_offsets(saved_offsets);
	res->set_array_offsets(0,0,0);
	res->set_shuffled(true);
	return res;
}



EMData* EMData::extract_plane_rect_fast(const Transform& tf, Util::KaiserBessel& kbx,Util::KaiserBessel& kby, Util::KaiserBessel& kbz) {
	
 	

	if (!is_complex())
		throw ImageFormatException("extractplane requires a complex image");
	if (nx%2 != 0)
		throw ImageDimensionException("extractplane requires nx to be even");

	int nxfromz=nz+2;
	int nxcircal = nxfromz - 2;

	// build complex result image
	float xratio=float(nx-2)/float(nz);
	float yratio=float(ny)/float(nz);
	Vec3f axis_newx,axis_newy;
	axis_newx[0] = xratio*0.5f*nz*tf[0][0];
	axis_newx[1] = yratio*0.5f*nz*tf[0][1];
	axis_newx[2] = 0.5f*nz*tf[0][2];


	float ellipse_length_x=std::sqrt(axis_newx[0]*axis_newx[0]+axis_newx[1]*axis_newx[1]+axis_newx[2]*axis_newx[2]);
	
	int ellipse_length_x_int=int(ellipse_length_x);
	float ellipse_step_x=0.5f*nz/float(ellipse_length_x_int);
	float xscale=ellipse_step_x;//scal increased

  	axis_newy[0] = xratio*0.5f*nz*tf[1][0];
  	axis_newy[1] = yratio*0.5f*nz*tf[1][1];
  	axis_newy[2] = 0.5f*nz*tf[1][2];


	float ellipse_length_y=std::sqrt(axis_newy[0]*axis_newy[0]+axis_newy[1]*axis_newy[1]+axis_newy[2]*axis_newy[2]);
	int ellipse_length_y_int=int(ellipse_length_y);
	float ellipse_step_y=0.5f*nz/float(ellipse_length_y_int);
	float yscale=ellipse_step_y;
	//end of scaling factor calculation
	int nx_e=ellipse_length_x_int*2;
	int ny_e=ellipse_length_y_int*2;
	int nx_ec=nx_e+2;	

	EMData* res = new EMData();
	res->set_size(nx_ec,ny_e,1);
	res->to_zero();
	res->set_complex(true);
	res->set_fftodd(false);
	res->set_fftpad(true);
	res->set_ri(true);
	//std::cout<<"cpp fast extract_plane is called"<<std::endl;
	//std::cout<<"nx_e,ny_e===="<<nx_e<<"  "<<ny_e<<std::endl;
	// Array offsets: (0..nhalf,-nhalf..nhalf-1,-nhalf..nhalf-1)
	int n = nxcircal;
	int nhalf = n/2;
	int nhalfx_e = nx_e/2;
	int nhalfy_e = ny_e/2;
	int nxhalf=(nx-2)/2;
	int nyhalf=ny/2;
	int nzhalf=nz/2;
	//std::cout<<"nhalf,nxhalf,nyhalf,nzhalf=="<<nhalf<<" "<<nxhalf<<" "<<nyhalf<<" "<<nzhalf<<std::endl;
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,-nyhalf,-nzhalf);
	res->set_array_offsets(0,-nhalfy_e,0);
	// set up some temporary weighting arrays
	int kbxsize =  kbx.get_window_size();
	int kbxmin  = -kbxsize/2;
	int kbxmax  = -kbxmin;

	int kbysize =  kby.get_window_size();
	int kbymin  = -kbysize/2;
	int kbymax  = -kbymin;

	int kbzsize =  kbz.get_window_size();
	int kbzmin  = -kbzsize/2;
	int kbzmax  = -kbzmin;

	//std::cout<<"kb size x,y,z=="<<kbxsize<<" "<<kbysize<<" "<<kbzsize<<std::endl;
	float* wy0 = new float[kbymax - kbymin + 1];
	float* wy  = wy0 - kbymin; // wy[kbmin:kbmax]
	float* wx0 = new float[kbxmax - kbxmin + 1];
	float* wx  = wx0 - kbxmin;
	float* wz0 = new float[kbzmax - kbzmin + 1];
	float* wz  = wz0 - kbzmin;
	float rim = nhalf*float(nhalf);
	int count = 0;
	float wsum = 0.f;
	Transform tftrans = tf; // need transpose of tf here for consistency
	tftrans.invert();      // with spider

	//std::cout<<"xratio,yratio=="<<xratio<<" "<<yratio<<std::endl;
	for (int jy = -nhalfy_e; jy < nhalfy_e; jy++) 
	{
		for (int jx = 0; jx <= nhalfx_e; jx++) 
		{
			Vec3f nucur((float)jx, (float)jy, 0.f);
			nucur[0]=nucur[0]*xscale;nucur[1]=nucur[1]*yscale;;
			Vec3f nunew = tftrans*nucur;
			float xnew = nunew[0]*xratio, ynew = nunew[1]*yratio, znew = nunew[2];

			if (nunew[0]*nunew[0]+nunew[1]*nunew[1]+nunew[2]*nunew[2] <= rim)
			{
				count++;
				std::complex<float> btq(0.f,0.f);
				bool flip = false;
				if (xnew < 0.f) {
					flip = true;
					xnew = -xnew;
					ynew = -ynew;
					znew = -znew;
				}
				int ixn = int(Util::round(xnew));
				int iyn = int(Util::round(ynew));
				int izn = int(Util::round(znew));
				// populate weight arrays
				for (int i=kbzmin; i <= kbzmax; i++) {
					int izp = izn + i;
					wz[i] = kbz.i0win_tab(znew - izp);
					}
				for (int i=kbymin; i <= kbymax; i++) {
					int iyp = iyn + i;
					wy[i] = kby.i0win_tab(ynew - iyp);
					}
				for (int i=kbxmin; i <= kbxmax; i++) {
					int ixp = ixn + i;
					wx[i] = kbx.i0win_tab(xnew - ixp);
					}
		

				
				// restrict weight arrays to non-zero elements
				int lnbz = 0;
				for (int iz = kbzmin; iz <= -1; iz++) {
					if (wz[iz] != 0.f) {
						lnbz = iz;
						break;
					}
				}
				int lnez = 0;
				for (int iz = kbzmax; iz >= 1; iz--) {
					if (wz[iz] != 0.f) {
						lnez = iz;
						break;
					}
				}
				int lnby = 0;
				for (int iy = kbymin; iy <= -1; iy++) {
					if (wy[iy] != 0.f) {
						lnby = iy;
						break;
					}
				}
				int lney = 0;
				for (int iy = kbymax; iy >= 1; iy--) {
					if (wy[iy] != 0.f) {
						lney = iy;
						break;
					}
				}
				int lnbx = 0;
				for (int ix = kbxmin; ix <= -1; ix++) {
					if (wx[ix] != 0.f) {
						lnbx = ix;
						break;
					}
				}
				int lnex = 0;
				for (int ix = kbxmax; ix >= 1; ix--) {
					if (wx[ix] != 0.f) {
						lnex = ix;
						break;
					}
				}
				if    (ixn >= -kbxmin      && ixn <= nxhalf-1-kbxmax
				   && iyn >= -nyhalf-kbymin && iyn <= nyhalf-1-kbymax
				   && izn >= -nzhalf-kbzmin && izn <= nzhalf-1-kbzmax) {
					// interior points
					for (int lz = lnbz; lz <= lnez; lz++) {
						int izp = izn + lz;
						for (int ly=lnby; ly<=lney; ly++) {
							int iyp = iyn + ly;
							float ty = wz[lz]*wy[ly];
							for (int lx=lnbx; lx<=lnex; lx++) {
								int ixp = ixn + lx;
								float wg = wx[lx]*ty;
								btq += cmplx(ixp,iyp,izp)*wg;
								wsum += wg;
							}
						}
					}
				} 
				else {
					// points "sticking out"
					for (int lz = lnbz; lz <= lnez; lz++) {
						int izp = izn + lz;
						for (int ly=lnby; ly<=lney; ly++) {
							int iyp = iyn + ly;
							float ty = wz[lz]*wy[ly];
							for (int lx=lnbx; lx<=lnex; lx++) {
								int ixp = ixn + lx;
								float wg = wx[lx]*ty;
								bool mirror = false;
								int ixt(ixp), iyt(iyp), izt(izp);
								if (ixt > nxhalf || ixt < -nxhalf) {
									ixt = Util::sgn(ixt)
										  *(nx-2-abs(ixt));
									iyt = -iyt;
									izt = -izt;
									mirror = !mirror;
								}
								if (iyt >= nyhalf || iyt < -nyhalf) {
									if (ixt != 0) {
										ixt = -ixt;
										iyt = Util::sgn(iyt)
											  *(ny - abs(iyt));
										izt = -izt;
										mirror = !mirror;
									} else {
										iyt -= ny*Util::sgn(iyt);
									}
								}
								if (izt >= nzhalf || izt < -nzhalf) {
									if (ixt != 0) {
										ixt = -ixt;
										iyt = -iyt;
										izt = Util::sgn(izt)
											  *(nz - abs(izt));
										mirror = !mirror;
									} else {
										izt -= Util::sgn(izt)*nz;
									}
								}
								if (ixt < 0) {
									ixt = -ixt;
									iyt = -iyt;
									izt = -izt;
									mirror = !mirror;
								}
								if (iyt == nyhalf) iyt = -nyhalf;
								if (izt == nzhalf) izt = -nzhalf;
								if (mirror)   btq += conj(cmplx(ixt,iyt,izt))*wg;
								else          btq += cmplx(ixt,iyt,izt)*wg;
								wsum += wg;
							}
						}
					}
				}
				if (flip)  res->cmplx(jx,jy) = conj(btq);
				else       res->cmplx(jx,jy) = btq;
			}
		}
	}
	for (int jy = -nhalfy_e; jy < nhalfy_e; jy++)
		for (int jx = 0; jx <= nhalfx_e; jx++)
			res->cmplx(jx,jy) *= count/wsum;
	delete[] wx0; delete[] wy0; delete[] wz0;
	set_array_offsets(saved_offsets);
	res->set_array_offsets(0,0,0);
	res->set_shuffled(true);
	return res;
}




bool EMData::peakcmp(const Pixel& p1, const Pixel& p2) {
    return (p1.value > p2.value);
}

ostream& operator<< (ostream& os, const Pixel& peak) {
    os <<  peak.x <<  peak.y << peak.z  << peak.value;
    return os;
}

/*vector<float> EMData::max_search() {

	EMData& buf = *this;

	int nx = buf.get_xsize();
	int ny = buf.get_ysize();
	int nz = buf.get_zsize();

	int dim = buf.get_ndim();

	vector<float> result;

	if (dim == 1) {
		float max = -1e20f;
		int index = -1;
		for (int i=0; i<nx; i++) {
			if (buf(i)>max) {
				max = buf(i);
				index = i;
			}
		}
		result.push_back((float)index);
		result.push_back(max);
	} else if (dim == 2) {
		float max = -1e20f;
		int index1 = -1;
		int index2 = -1;
		for (int i=0; i<nx; i++) {
			for (int j=0; j<ny; j++) {
				if (buf(i, j)>max) {
					max = buf(i, j);
					index1 = i;
					index2 = j;
				}
			}
		}
		result.push_back((float)index1);
		result.push_back((float)index2);
		result.push_back(max);
	} else {
		float max = -1e20f;
		int index1 = -1;
		int index2 = -1;
		int index3 = -1;
		for (int i=0; i<nx; i++) {
			for (int j=0; j<ny; j++) {
				for (int k=0; k<nz; k++) {
					if (buf(i, j, k)>max) {
						max = buf(i, j, k);
						index1 = i;
						index2 = j;
						index3 = k;
					}
				}
			}
		}
		result.push_back((float)index1);
		result.push_back((float)index2);
		result.push_back((float)index3);
		result.push_back(max);
	}
	return result;
}*/

/*
Divide volume by 
 interpolate_method = 0  -> sinc
 interpolate_method = 1  -> sinc^2
 We assume volume has equal sizes!!
*/
void EMData::div_sinc(int interpolate_method) {
	ENTERFUNC;
	if (is_complex()) {
		LOGERR("Real image expected. Input image is complex.");
		throw ImageFormatException("Real image expected. Input image is complex.");
	}
    vector<int> saved_offsets = get_array_offsets();
	int nx = this->get_xsize();
	int ny = this->get_ysize();
	int nz = this->get_zsize();
	if (nz>1 && (nx != ny || ny != nz))
		throw ImageDimensionException("div_sinc requires ny == nx == nz");
	else if(nx != ny)
		throw ImageDimensionException("div_sinc requires ny == nx");
	int npad = this->get_attr("npad");
	float cdf = M_PI/(nx);
/*
	int IP = nx/2+1;

	//  tabulate sinc function
	float* sincx = new float[IP];
	sincx[0] = 1.0f;
	//  It is 1/sinc
	if(interpolate_method == 0) {
		for (int i = 1; i < IP; ++i)  sincx[i] = (i*cdf)/sin(i*cdf);
	} else {
		for (int i = 1; i < IP; ++i)  sincx[i] = pow((i*cdf)/sin(i*cdf),2);	
	}
*/
	set_array_offsets(-nx/2,-ny/2,-nz/2);

/*
	for (int k = -nz/2; k < nz/2 + nz%2; k++) {
		float kkp = sincx[abs(k)];
		for (int j = -ny/2; j < ny/2 + ny%2; j++) {
			cdf = sincx[abs(j)]*kkp;
			for (int i = -nx/2; i < nx/2 + nx%2; i++) (*this)(i,j,k) *= sincx[abs(i)]*cdf;
		}
	}
*/
	int nzb = nz/2/npad;
	int nyb = ny/2/npad;
	int nxb = nx/2/npad;
	int nze = nzb + (nz/npad)%2;
	int nye = nyb + (ny/npad)%2;
	int nxe = nxb + (nx/npad)%2;
	if(nz == 1) nze = 1;
	for (int k = -nzb; k < nze; k++) {
		for (int j = -nyb; j < nye; j++) {
			for (int i = -nxb; i < nxe; i++) {
				float  rrr = std::sqrt(k*k+j*j+float(i*i));
				if(rrr>0.0f)  {
					if( interpolate_method == 1 ) (*this)(i,j,k) *= pow((rrr*cdf)/sin(rrr*cdf),2);
					else  (*this)(i,j,k) *= (rrr*cdf)/sin(rrr*cdf);
				}
			}
		}
	}
	set_array_offsets(saved_offsets);
	update();
	EXITFUNC;
}


EMData* EMData::extract_section(const Transform& tf, int interpolate_method) {
	if (!is_complex())
		throw ImageFormatException("extract_section requires a complex image");
	if (nx%2 != 0)
		throw ImageDimensionException("extract_section requires nx to be even");
	int nxreal = nx - 2;
	if (nxreal != ny || nxreal != nz)
		throw ImageDimensionException("extract_section requires ny == nx == nz");
	// build complex result image
	EMData* res = new EMData();
	res->set_size(nx,ny,1);
	res->to_zero();
	res->set_complex(true);
	res->set_fftodd(false);
	res->set_fftpad(true);
	res->set_ri(true);
	// Array offsets: (0..nhalf,-nhalf..nhalf-1,-nhalf..nhalf-1)
	int n = nxreal;
	int nhalf = n/2;
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,-nhalf,-nhalf);
	res->set_array_offsets(0,-nhalf,0);


	float rim = nhalf*float(nhalf);
	int count = 0;
	float wsum = 0.f;
	Transform tftrans = tf; // need transpose of tf here for consistency
	tftrans.invert();      // with spider
	if( interpolate_method == 0 ) {
		for (int jy = -nhalf; jy < nhalf; jy++)  {
			for (int jx = 0; jx <= nhalf; jx++) {
				Vec3f nucur((float)jx, (float)jy, 0.f);
				Vec3f nunew = tftrans*nucur;
				float xnew = nunew[0], ynew = nunew[1], znew = nunew[2];
				if (xnew*xnew+ynew*ynew+znew*znew <= rim) {
					int ixn = int(Util::round(xnew));
					int iyn = int(Util::round(ynew));
					int izn = int(Util::round(znew));
					bool flip = false;
					if (ixn < 0) {
						flip = true;
						ixn = -ixn;
						iyn = -iyn;
						izn = -izn;
					}
					if(   ixn >= 0      && ixn <= nhalf-1
					   && iyn >= -nhalf && iyn <= nhalf-1
					   && izn >= -nhalf && izn <= nhalf-1) {
						//std::cout<<"    2D: "<<jx<<"  "<<jy<<"        3D: "<<ixn<<"  "<<iyn<<"  "<<izn<<"  "<<flip<<"  "<<std::endl;
						if (flip)  res->cmplx(jx,jy) = conj(cmplx(ixn,iyn,izn));
						else       res->cmplx(jx,jy) = cmplx(ixn,iyn,izn);
					} else {
						// points "sticking out"
						bool mirror = false;
						int ixt(ixn), iyt(iyn), izt(izn);
						if (ixt > nhalf || ixt < -nhalf) {
							ixt = Util::sgn(ixt)*(n - abs(ixt));
							iyt = -iyt;
							izt = -izt;
							mirror = !mirror;
						}
						if (iyt >= nhalf || iyt < -nhalf) {
							if (ixt != 0) {
								ixt = -ixt;
								iyt = Util::sgn(iyt)*(n - abs(iyt));
								izt = -izt;
								mirror = !mirror;
							} else {
								iyt -= n*Util::sgn(iyt);
							}
						}
						if (izt >= nhalf || izt < -nhalf) {
							if (ixt != 0) {
								ixt = -ixt;
								iyt = -iyt;
								izt = Util::sgn(izt)*(n - abs(izt));
								mirror = !mirror;
							} else {
								izt -= Util::sgn(izt)*n;
							}
						}
						if (ixt < 0) {
							ixt = -ixt;
							iyt = -iyt;
							izt = -izt;
							mirror = !mirror;
						}
						if (iyt == nhalf) iyt = -nhalf;
						if (izt == nhalf) izt = -nhalf;
						if (mirror)   res->cmplx(jx,jy) = conj(cmplx(ixt,iyt,izt));
						else          res->cmplx(jx,jy) = cmplx(ixt,iyt,izt);
					}
				}
			}
		}
	} else {
		// tri-linear interpolation
		for (int jy = -nhalf; jy < nhalf; jy++)  {
			for (int jx = 0; jx <= nhalf; jx++) {
				Vec3f nucur((float)jx, (float)jy, 0.f);
				Vec3f nunew = tftrans*nucur;
				float xnew = nunew[0], ynew = nunew[1], znew = nunew[2];
				if (xnew*xnew+ynew*ynew+znew*znew <= rim) {
					bool flip = false;
					if (xnew < 0.f) {
						flip = true;
						xnew = -xnew;
						ynew = -ynew;
						znew = -znew;
					}

					int ixn = int(xnew + n);
					int iyn = int(ynew + n);
					int izn = int(znew + n);

					float dx = xnew + n - ixn;
					float dy = ynew + n - iyn;
					float dz = znew + n - izn;

					ixn -= n;
					iyn -= n;
					izn -= n;

					if(   ixn >= 0      && ixn <= nhalf-2
					   && iyn >= -nhalf && iyn <= nhalf-2
					   && izn >= -nhalf && izn <= nhalf-2) {
						//std::cout<<"    2D: "<<jx<<"  "<<jy<<"        3D: "<<ixn<<"  "<<iyn<<"  "<<izn<<"  "<<flip<<"  "<<std::endl;

						std::complex<float> a1 = cmplx(ixn,iyn,izn);
						std::complex<float> a2 = cmplx(ixn+1,iyn,izn) - cmplx(ixn,iyn,izn);
						std::complex<float> a3 = cmplx(ixn,iyn+1,izn) - cmplx(ixn,iyn,izn);
						std::complex<float> a4 = cmplx(ixn,iyn,izn+1) - cmplx(ixn,iyn,izn);
						std::complex<float> a5 = cmplx(ixn,iyn,izn) - cmplx(ixn+1,iyn,izn) - cmplx(ixn,iyn+1,izn) + cmplx(ixn+1,iyn+1,izn);
						std::complex<float> a6 = cmplx(ixn,iyn,izn) - cmplx(ixn+1,iyn,izn) - cmplx(ixn,iyn,izn+1) + cmplx(ixn+1,iyn,izn+1);
						std::complex<float> a7 = cmplx(ixn,iyn,izn) - cmplx(ixn,iyn+1,izn) - cmplx(ixn,iyn,izn+1) + cmplx(ixn,iyn+1,izn+1);
						std::complex<float> a8 = cmplx(ixn+1,iyn,izn) + cmplx(ixn,iyn+1,izn)+ cmplx(ixn,iyn,izn+1)
								- cmplx(ixn,iyn,izn)- cmplx(ixn+1,iyn+1,izn) - cmplx(ixn+1,iyn,izn+1)
								- cmplx(ixn,iyn+1,izn+1) + cmplx(ixn+1,iyn+1,izn+1);
						std::complex<float> btq = a1 + dz*(a4 + a6*dx + (a7 + a8*dx)*dy) + a3*dy + dx*(a2 + a5*dy);
						if (flip)  res->cmplx(jx,jy) = conj(btq);
						else       res->cmplx(jx,jy) = btq;
					} else {
						// points "sticking out" just do using linear
						bool mirror = false;
						int ixt(ixn), iyt(iyn), izt(izn);
						if (ixt > nhalf || ixt < -nhalf) {
							ixt = Util::sgn(ixt)*(n - abs(ixt));
							iyt = -iyt;
							izt = -izt;
							mirror = !mirror;
						}
						if (iyt >= nhalf || iyt < -nhalf) {
							if (ixt != 0) {
								ixt = -ixt;
								iyt = Util::sgn(iyt)*(n - abs(iyt));
								izt = -izt;
								mirror = !mirror;
							} else {
								iyt -= n*Util::sgn(iyt);
							}
						}
						if (izt >= nhalf || izt < -nhalf) {
							if (ixt != 0) {
								ixt = -ixt;
								iyt = -iyt;
								izt = Util::sgn(izt)*(n - abs(izt));
								mirror = !mirror;
							} else {
								izt -= Util::sgn(izt)*n;
							}
						}
						if (ixt < 0) {
							ixt = -ixt;
							iyt = -iyt;
							izt = -izt;
							mirror = !mirror;
						}
						if (iyt == nhalf) iyt = -nhalf;
						if (izt == nhalf) izt = -nhalf;
						if (mirror)   res->cmplx(jx,jy) = conj(cmplx(ixt,iyt,izt));
						else          res->cmplx(jx,jy) = cmplx(ixt,iyt,izt);
					}
				}
			}
		}
	}
	set_array_offsets(saved_offsets);
	res->set_array_offsets(0,0,0);
	res->set_shuffled(true);
	return res;
}


EMData* EMData::extract_section2(const Transform& tf, int interpolate_method) {
	if (!is_complex())
		throw ImageFormatException("extract_section2 requires a complex image");
	if (nx%2 != 0)
		throw ImageDimensionException("extract_section2 requires nx to be even");
	int nxreal = nx - 2;
	if (nxreal != ny || nxreal != nz)
		throw ImageDimensionException("extract_section2 requires ny == nx == nz");
	// build complex result image, the assumption is that incoming volume is even and padded 2x
	int nxo = nxreal/2 +2 ;
	int nyo = ny/2;
	EMData* res = new EMData();
	res->set_size(nxo,nyo,1);
	res->to_zero();
	res->set_complex(true);
	res->set_fftodd(false);
	res->set_fftpad(false);
	res->set_ri(true);
	// Array offsets: (0..nhalf,-nhalf..nhalf-1,-nhalf..nhalf-1)
	int n = nxreal;
	int nhalf = ny/2;
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,-nhalf,-nhalf);
	int nhalfo = nyo/2;
	res->set_array_offsets(0,-nhalfo,0);

//	int ixg =0, iyg=0, izg=0;
	float rim = nhalfo*float(nhalfo);
	Transform tftrans = tf; // need transpose of tf here for consistency
	tftrans.invert();      // with spider
	if( interpolate_method == 0 ) {
		for (int jy = -nhalf; jy < nhalf; jy++)  {
			for (int jx = 0; jx <= nhalf; jx++) {
				Vec3f nucur((float)jx, (float)jy, 0.f);
				Vec3f nunew = tftrans*nucur;
				float xnew = nunew[0], ynew = nunew[1], znew = nunew[2];
				if (xnew*xnew+ynew*ynew+znew*znew <= rim) {
					int ixn = int(Util::round(xnew));
					int iyn = int(Util::round(ynew));
					int izn = int(Util::round(znew));
					bool flip = false;
					if (ixn < 0) {
						flip = true;
						ixn = -ixn;
						iyn = -iyn;
						izn = -izn;
					}
					if(   ixn >= 0      && ixn <= nhalf-1
					   && iyn >= -nhalf && iyn <= nhalf-1
					   && izn >= -nhalf && izn <= nhalf-1) {
						//std::cout<<"    2D: "<<jx<<"  "<<jy<<"        3D: "<<ixn<<"  "<<iyn<<"  "<<izn<<"  "<<flip<<"  "<<std::endl;
						if (flip)  res->cmplx(jx,jy) = conj(cmplx(ixn,iyn,izn));
						else       res->cmplx(jx,jy) = cmplx(ixn,iyn,izn);
					} else {
						// points "sticking out"
						bool mirror = false;
						int ixt(ixn), iyt(iyn), izt(izn);
						if (ixt > nhalf || ixt < -nhalf) {
							ixt = Util::sgn(ixt)*(n - abs(ixt));
							iyt = -iyt;
							izt = -izt;
							mirror = !mirror;
						}
						if (iyt >= nhalf || iyt < -nhalf) {
							if (ixt != 0) {
								ixt = -ixt;
								iyt = Util::sgn(iyt)*(n - abs(iyt));
								izt = -izt;
								mirror = !mirror;
							} else {
								iyt -= n*Util::sgn(iyt);
							}
						}
						if (izt >= nhalf || izt < -nhalf) {
							if (ixt != 0) {
								ixt = -ixt;
								iyt = -iyt;
								izt = Util::sgn(izt)*(n - abs(izt));
								mirror = !mirror;
							} else {
								izt -= Util::sgn(izt)*n;
							}
						}
						if (ixt < 0) {
							ixt = -ixt;
							iyt = -iyt;
							izt = -izt;
							mirror = !mirror;
						}
						if (iyt == nhalf) iyt = -nhalf;
						if (izt == nhalf) izt = -nhalf;
						if (mirror)   res->cmplx(jx,jy) = conj(cmplx(ixt,iyt,izt));
						else          res->cmplx(jx,jy) = cmplx(ixt,iyt,izt);
					}
				}
			}
		}
	} else {
		// tri-linear interpolation
		for (int jy = -nhalfo; jy < nhalfo; jy++)  {
			for (int jx = 0; jx <= nhalfo; jx++) {
				if (jx*jx+jy*jy <= rim) {
					Vec3f nucur((float)(2*jx), (float)(2*jy), 0.f);
					Vec3f nunew = tftrans*nucur;
					float xnew = nunew[0], ynew = nunew[1], znew = nunew[2];
//	std::cout<<"INDX  "<<jy+nhalf<<"  "<<jy<<"  "<<jx<<"  "<<xnew<<"  "<<ynew<<"  "<<znew<<std::endl;
					bool flip = false;
					bool xmirror;
					if (xnew < 0.f) {
						flip = true;
						xnew = -xnew;
						ynew = -ynew;
						znew = -znew;
					}

					int ixn = int(xnew + n);
					int iyn = int(ynew + n);
					int izn = int(znew + n);

					float dx = xnew + n - ixn;
					float dy = ynew + n - iyn;
					float dz = znew + n - izn;

					ixn -= n;
					iyn -= n;
					izn -= n;

					if(   ixn >= 0      && ixn <= nhalf-2
					   && iyn >= -nhalf && iyn <= nhalf-2
					   && izn >= -nhalf && izn <= nhalf-2) {
				/*		std::cout<<"    2D: "<<jx<<"  "<<jy<<"        3D: "<<ixn<<"  "<<iyn<<"  "<<izn<<"  "<<flip<<"  "<<std::endl;
						ixg = ixn;
						iyg = iyn;
						izg = izn; */
						std::complex<float> a1 = cmplx(ixn,iyn,izn);
						std::complex<float> a2 = cmplx(ixn+1,iyn,izn) - cmplx(ixn,iyn,izn);
						std::complex<float> a3 = cmplx(ixn,iyn+1,izn) - cmplx(ixn,iyn,izn);
						std::complex<float> a4 = cmplx(ixn,iyn,izn+1) - cmplx(ixn,iyn,izn);
						std::complex<float> a5 = cmplx(ixn,iyn,izn) - cmplx(ixn+1,iyn,izn) - cmplx(ixn,iyn+1,izn) + cmplx(ixn+1,iyn+1,izn);
						std::complex<float> a6 = cmplx(ixn,iyn,izn) - cmplx(ixn+1,iyn,izn) - cmplx(ixn,iyn,izn+1) + cmplx(ixn+1,iyn,izn+1);
						std::complex<float> a7 = cmplx(ixn,iyn,izn) - cmplx(ixn,iyn+1,izn) - cmplx(ixn,iyn,izn+1) + cmplx(ixn,iyn+1,izn+1);
						std::complex<float> a8 = cmplx(ixn+1,iyn,izn) + cmplx(ixn,iyn+1,izn)+ cmplx(ixn,iyn,izn+1)
								- cmplx(ixn,iyn,izn)- cmplx(ixn+1,iyn+1,izn) - cmplx(ixn+1,iyn,izn+1)
								- cmplx(ixn,iyn+1,izn+1) + cmplx(ixn+1,iyn+1,izn+1);
						std::complex<float> btq = a1 + dz*(a4 + a6*dx + (a7 + a8*dx)*dy) + a3*dy + dx*(a2 + a5*dy);
						if (flip)  res->cmplx(jx,jy) = conj(btq);
						else       res->cmplx(jx,jy) = btq;
					} else {
						// points "sticking out" just do using linear
						bool mirror = false;
						int ixt(ixn), iyt(iyn), izt(izn);
						if (ixt > nhalf || ixt < -nhalf) {
							ixt = Util::sgn(ixt)*(n - abs(ixt));
							iyt = -iyt;
							izt = -izt;
							mirror = !mirror;
						}
						if (iyt >= nhalf || iyt < -nhalf) {
							if (ixt != 0) {
								ixt = -ixt;
								iyt = Util::sgn(iyt)*(n - abs(iyt));
								izt = -izt;
								mirror = !mirror;
							} else {
								iyt -= n*Util::sgn(iyt);
							}
						}
						if (izt >= nhalf || izt < -nhalf) {
							if (ixt != 0) {
								ixt = -ixt;
								iyt = -iyt;
								izt = Util::sgn(izt)*(n - abs(izt));
								mirror = !mirror;
							} else {
								izt -= Util::sgn(izt)*n;
							}
						}
						if (ixt < 0) {
							ixt = -ixt;
							iyt = -iyt;
							izt = -izt;
							mirror = !mirror;
						}
						if (iyt == nhalf) iyt = -nhalf;
						if (izt == nhalf) izt = -nhalf;
						if (mirror)   res->cmplx(jx,jy) = conj(cmplx(ixt,iyt,izt));
						else          res->cmplx(jx,jy) = cmplx(ixt,iyt,izt);
						
						
					/*	xmirror = mirror;
						ixg = ixt;
						iyg = iyt;
						izg = izt; */


					}
	//std::cout<<"PROJ  "<<jy+nhalfo<<"  "<<jy<<"  "<<jx<<" | "<<xnew<<"  "<<ynew<<"  "<<znew<<" | "<<ixn<<"  "<<iyn<<"  "<<izn<<"  |  "<<std::real(res->cmplx(jx,jy))<<"  " <<std::imag(res->cmplx(jx,jy))<<" | "<<xmirror<<" | "<<std::real(cmplx(ixg,iyg,izg))<<"  "<<std::imag(cmplx(ixg,iyg,izg))<<std::endl;
				}
			}
		}
	}
	set_array_offsets(saved_offsets);
	res->set_array_offsets(0,0,0);
	res->set_shuffled(true);
	return res;
}


EMData* EMData::extract_rings(const Transform& tf, int interpolate_method, int ring_length, int nb, int ne) {
	if (!is_complex())
		throw ImageFormatException("extract_rings requires a complex image");
	if (nx%2 != 0)
		throw ImageDimensionException("extract_rings requires nx to be even");
	int nxreal = nx - 2;
	if (nxreal != ny || nxreal != nz)
		throw ImageDimensionException("extract_rings requires ny == nx == nz");
	// build complex ringsult image, the assumption is that incoming volume is even and padded 2x
	int nxo = nxreal/2 +2 ;
	int nyo = ny/2;


	int lcirc = ne-nb+1;
	int nc = ny/2;// ?????


	float dfi;
	dfi = TWOPI / ring_length;
//	Table for sin & cos
	vector<float> vsin(ring_length/2);
	vector<float> vcos(ring_length/2);
	for (int x = 0; x < ring_length/2; x++) {
		float ang = static_cast<float>(x * dfi);
		vsin[x] = sin(ang);
		vcos[x] = cos(ang);
		//printf("trigtab   %d      %f  %f\n",x,vsin[x],vcos[x]);
	}

	
	EMData* rings = new EMData(2*ring_length-2, lcirc, 1, false);
	rings->set_fftpad(0);
	rings->set_complex(true);
	rings->set_fftodd(false);
	rings->set_fftpad(false);
	rings->set_ri(true);
	// Array offsets: (0..nhalf,-nhalf..nhalf-1,-nhalf..nhalf-1)
	int n = nxreal;
	int nhalf = ny/2;
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,-nhalf,-nhalf);
	int nhalfo = nyo/2;
	//rings->set_array_offsets(0,-nb,0);

//	int ixg =0, iyg=0, izg=0;
	float rim = nhalfo*float(nhalfo);
	Transform tftrans = tf; // need transpose of tf here for consistency
	tftrans.invert();      // with spider
	if( interpolate_method == 0 ) {
		for (unsigned int inr = nb; inr <= ne; inr++) {
			for (unsigned int it = 0; it < ring_length/2; it++) {
				Vec3f nucur(vsin[it] * inr * 2, vcos[it] * inr * 2, 0.f);
				Vec3f nunew = tftrans*nucur;
				float xnew = nunew[0], ynew = nunew[1], znew = nunew[2];
				int ixn = int(Util::round(xnew));
				int iyn = int(Util::round(ynew));
				int izn = int(Util::round(znew));
				bool flip = false;
				if (ixn < 0) {
					flip = true;
					ixn = -ixn;
					iyn = -iyn;
					izn = -izn;
				}
				if(   ixn >= 0      && ixn <= nhalf-1
				   && iyn >= -nhalf && iyn <= nhalf-1
				   && izn >= -nhalf && izn <= nhalf-1) {
					//std::cout<<"    2D: "<<jx<<"  "<<jy<<"        3D: "<<ixn<<"  "<<iyn<<"  "<<izn<<"  "<<flip<<"  "<<std::endl;
					std::complex<float> btq = cmplx(ixn,iyn,izn);
					if (flip)  btq = conj(btq);
					rings->set_value_at(2*it,inr-nb,std::real(btq));
					rings->set_value_at(2*it+1,inr-nb,std::imag(btq));
					rings->set_value_at(2*it+ring_length,inr-nb,std::real(btq));
					rings->set_value_at(2*it+1+ring_length,inr-nb,-std::imag(btq));
				} else {
					printf(" NN  WHY WOULD IT BE HERE??\n");
					// points "sticking out"
					bool mirror = false;
					int ixt(ixn), iyt(iyn), izt(izn);
					if (ixt > nhalf || ixt < -nhalf) {
						ixt = Util::sgn(ixt)*(n - abs(ixt));
						iyt = -iyt;
						izt = -izt;
						mirror = !mirror;
					}
					if (iyt >= nhalf || iyt < -nhalf) {
						if (ixt != 0) {
							ixt = -ixt;
							iyt = Util::sgn(iyt)*(n - abs(iyt));
							izt = -izt;
							mirror = !mirror;
						} else {
							iyt -= n*Util::sgn(iyt);
						}
					}
					if (izt >= nhalf || izt < -nhalf) {
						if (ixt != 0) {
							ixt = -ixt;
							iyt = -iyt;
							izt = Util::sgn(izt)*(n - abs(izt));
							mirror = !mirror;
						} else {
							izt -= Util::sgn(izt)*n;
						}
					}
					if (ixt < 0) {
						ixt = -ixt;
						iyt = -iyt;
						izt = -izt;
						mirror = !mirror;
					}
					if (iyt == nhalf) iyt = -nhalf;
					if (izt == nhalf) izt = -nhalf;
					std::complex<float> btq = cmplx(ixt,iyt,izt);
					if (flip)  btq = conj(btq);
					rings->set_value_at(2*it,inr-nb,std::real(btq));
					rings->set_value_at(2*it+1,inr-nb,std::imag(btq));
					rings->set_value_at(2*it+ring_length,inr-nb,std::real(btq));
					rings->set_value_at(2*it+1+ring_length,inr-nb,-std::imag(btq));
				}
			}
		}
	} else {
		// tri-linear interpolation
		for (unsigned int inr = nb; inr <= ne; inr++) {
			for (unsigned int it = 0; it < ring_length/2; it++) {
				Vec3f nucur(vsin[it] * inr * 2, vcos[it] * inr * 2, 0.f);
				Vec3f nunew = tftrans*nucur;
				float xnew = nunew[0], ynew = nunew[1], znew = nunew[2];
//	std::cout<<"INDX  "<<jy+nhalf<<"  "<<jy<<"  "<<jx<<"  "<<xnew<<"  "<<ynew<<"  "<<znew<<std::endl;
				bool flip = false;
				bool xmirror;
				if (xnew < 0.f) {
					flip = true;
					xnew = -xnew;
					ynew = -ynew;
					znew = -znew;
				}

				int ixn = int(xnew + n);
				int iyn = int(ynew + n);
				int izn = int(znew + n);

				float dx = xnew + n - ixn;
				float dy = ynew + n - iyn;
				float dz = znew + n - izn;

				ixn -= n;
				iyn -= n;
				izn -= n;

				if(   ixn >= 0      && ixn <= nhalf-2
				   && iyn >= -nhalf && iyn <= nhalf-2
				   && izn >= -nhalf && izn <= nhalf-2) {
			/*		std::cout<<"    2D: "<<jx<<"  "<<jy<<"        3D: "<<ixn<<"  "<<iyn<<"  "<<izn<<"  "<<flip<<"  "<<std::endl;
					ixg = ixn;
					iyg = iyn;
					izg = izn; */
					std::complex<float> a1 = cmplx(ixn,iyn,izn);
					std::complex<float> a2 = cmplx(ixn+1,iyn,izn) - cmplx(ixn,iyn,izn);
					std::complex<float> a3 = cmplx(ixn,iyn+1,izn) - cmplx(ixn,iyn,izn);
					std::complex<float> a4 = cmplx(ixn,iyn,izn+1) - cmplx(ixn,iyn,izn);
					std::complex<float> a5 = cmplx(ixn,iyn,izn) - cmplx(ixn+1,iyn,izn) - cmplx(ixn,iyn+1,izn) + cmplx(ixn+1,iyn+1,izn);
					std::complex<float> a6 = cmplx(ixn,iyn,izn) - cmplx(ixn+1,iyn,izn) - cmplx(ixn,iyn,izn+1) + cmplx(ixn+1,iyn,izn+1);
					std::complex<float> a7 = cmplx(ixn,iyn,izn) - cmplx(ixn,iyn+1,izn) - cmplx(ixn,iyn,izn+1) + cmplx(ixn,iyn+1,izn+1);
					std::complex<float> a8 = cmplx(ixn+1,iyn,izn) + cmplx(ixn,iyn+1,izn)+ cmplx(ixn,iyn,izn+1)
							- cmplx(ixn,iyn,izn)- cmplx(ixn+1,iyn+1,izn) - cmplx(ixn+1,iyn,izn+1)
							- cmplx(ixn,iyn+1,izn+1) + cmplx(ixn+1,iyn+1,izn+1);
					std::complex<float> btq = a1 + dz*(a4 + a6*dx + (a7 + a8*dx)*dy) + a3*dy + dx*(a2 + a5*dy);
					/*if (flip)  rings->cmplx(it,inr) = conj(btq);
					else       rings->cmplx(it,inr) = btq;
					printf(" it inr  flip  %d  %d  %d  %f  %f  \n",it,inr,flip,std::real(rings->cmplx(it,inr)),std::imag(rings->cmplx(it,inr)) );
					rings->cmplx(it+ring_length/2,inr) = conj(rings->cmplx(it,inr));
					*/
					if (flip)  btq = conj(btq);
					rings->set_value_at(2*it,inr-nb,std::real(btq));
					rings->set_value_at(2*it+1,inr-nb,std::imag(btq));
					rings->set_value_at(2*it+ring_length,inr-nb,std::real(btq));
					rings->set_value_at(2*it+1+ring_length,inr-nb,-std::imag(btq));

//std::cout<<"PROJ  "<<jy+nhalfo<<"  "<<jy<<"  "<<jx<<" | "<<xnew<<"  "<<ynew<<"  "<<znew<<" | "<<ixn<<"  "<<iyn<<"  "<<izn<<"  |  "<<std::real(rings->cmplx(jx,jy))<<"  " <<std::imag(rings->cmplx(jx,jy))<<" | "<<xmirror<<" | "<<std::real(cmplx(ixg,iyg,izg))<<"  "<<std::imag(cmplx(ixg,iyg,izg))<<std::endl;
				}  else printf(" WHY WOULD IT BE HERE??\n");
			}
		}
	}
	set_array_offsets(saved_offsets);
	//rings->set_array_offsets(0,0,0);
	//rings->set_shuffled(true);
	return rings;
}


vector<float> EMData::peak_search(int ml, float invert) {

	EMData& buf = *this;
	vector<Pixel> peaks;
	int img_dim;
	int i, j, k;
	int i__1, i__2;
	int j__1, j__2;
	//int k__1, k__2;
 	bool peak_check;
 	img_dim=buf.get_ndim();
 	vector<int> ix, jy, kz;
	vector<float>res;
 	int nx = buf.get_xsize();
 	int ny = buf.get_ysize();
 	int nz = buf.get_zsize();
	if(invert <= 0.0f)  invert=-1.0f;
	else                invert=1.0f ;
	int count = 0;
 	switch (img_dim)  {
 	case(1):
		for(i=0;i<=nx-1;++i)  {
 		 	i__1 = (i-1+nx)%nx;
 		 	i__2 = (i+1)%nx;
		 	// Commented by Yang on 05/14/07
		 	// I changed the following line from > to >=, or in some rare cases (the peak happens to be flat), it will fail to find the peak.
		 	//  03/07/08  I undid the change.  If you change the comparison, it changes the meaning of peak definition.
			float qbf = buf(i)*invert;
 			peak_check = qbf > buf(i__1)*invert && qbf > buf(i__2)*invert;
			if(peak_check) {
	    			if(count < ml) {
					count++;
					peaks.push_back( Pixel(i, 0, 0, qbf) );
					if(count == ml-1) sort(peaks.begin(), peaks.end(), peakcmp);
				} else {
					if( qbf > (peaks.back()).value ) {
						//  do the switch and sort again
						peaks.pop_back();
						peaks.push_back( Pixel(i, 0, 0, qbf) );
						if(ml > 1) sort(peaks.begin(), peaks.end(), peakcmp);
					}
				}
			}
		}
 	break;
 	case(2):
	/*  Removed boundary conditions, PAP 03/10/08
		for(j=0;j<=ny-1;++j)  {
	 		j__1 = (j-1+ny)%ny;
 			j__2 = (j+1)%ny;
	 		for(i=0;i<=nx-1;++i) {
	        		i__1 = (i-1+nx)%nx;
	        		i__2 = (i+1)%nx;
	*/
		for(j=1;j<=ny-2;++j)  {
	 		j__1 = j-1;
 			j__2 = j+1;
	 		for(i=1;i<=nx-2;++i) {
	        		i__1 = i-1;
	        		i__2 = i+1;
				float qbf = buf(i,j)*invert;
	        		peak_check = (qbf > buf(i,j__1)*invert) && (qbf > buf(i,j__2)*invert);
				if(peak_check) {
					peak_check = (qbf > buf(i__1,j)*invert) && (qbf > buf(i__2,j)*invert);
	        			if(peak_check) {
						peak_check = (qbf > buf(i__1,j__1)*invert) && (qbf > buf(i__1,j__2)*invert);
	        				if(peak_check) {
							peak_check = (qbf > buf(i__2,j__1)*invert) && (qbf > buf(i__2,j__2)*invert);
							if(peak_check) {
	    							if(count < ml) {
									count++;
									peaks.push_back( Pixel(i, j, 0, qbf) );
									if(count == ml-1) sort(peaks.begin(), peaks.end(), peakcmp);
								} else {
									if( qbf > (peaks.back()).value ) {
										//  do the switch and sort again
										peaks.pop_back();
										peaks.push_back( Pixel(i, j, 0, qbf) );
										if(ml > 1) sort(peaks.begin(), peaks.end(), peakcmp);
									}
								}
							}
						}
					}
				}
	        	}
		}
 	break;
 	case(3):  //looks ugly, but it is the best I can do,  PAP 03/07/08
	/*  Removed boundary conditions, PAP 03/10/08
		for(k=0;k<=nz-1;++k) {
		 	kz.clear();
			k__1 = (k-1+nz)%nz;
			k__2 = (k+1)%nz;
			kz.push_back(k__1);
			kz.push_back(k);
			kz.push_back(k__2);
 			for(j=0;j<=ny-1;++j) {
				jy.clear();
				j__1 = (j-1+ny)%ny;
 				j__2 = (j+1)%ny;
				jy.push_back(j__1);
				jy.push_back(j);
				jy.push_back(j__2);
 				for(i=0;i<=nx-1;++i) {
			        	ix.clear();
			        	i__1 = (i-1+nx)%nx;
			        	i__2 = (i+1)%nx;
	*/
		for(k=1; k<=nz-2; ++k) {
		 	kz.clear();
			kz.push_back(k-1);
			kz.push_back(k);
			kz.push_back(k+1);
 			for(j=1; j<=ny-2; ++j) {
				jy.clear();
				jy.push_back(j-1);
				jy.push_back(j);
				jy.push_back(j+1);
 				for(i=1; i<=nx-2; ++i) {
			        	ix.clear();
			        	ix.push_back(i-1);
			        	ix.push_back(i);
			        	ix.push_back(i+1);
					float qbf = buf(i,j,k)*invert;
			        	peak_check = qbf > buf(ix[0],jy[0],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[1])*invert;
					//if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[2])*invert;
					if(peak_check) {
	    					if(count < ml) {
							count++;
							peaks.push_back( Pixel(i, j, k, qbf) );
							if(count == ml-1) sort(peaks.begin(), peaks.end(), peakcmp);
						} else {
							if( qbf > (peaks.back()).value ) {
								//  do the switch and sort again
								//cout << qbf<<"   "<< (peaks.back()).value <<"   "<< (*peaks.begin()).value <<endl;
								peaks.pop_back();
								peaks.push_back( Pixel(i, j, k, qbf) );
								if(ml > 1) sort(peaks.begin(), peaks.end(), peakcmp);
							}
						}
					}
					}}}}}}}}}}}}}}}}}}}}}}}}}
				}
			}
		}
		//  Add circular closure for x direction: needed for circular ccf,
		//  should not have adverse impact on other code.  PAP -7/22/08
		for(k=1; k<=nz-2; ++k) {
		 	kz.clear();
			kz.push_back(k-1);
			kz.push_back(k);
			kz.push_back(k+1);
 			for(j=1; j<=ny-2; ++j) {
				jy.clear();
				jy.push_back(j-1);
				jy.push_back(j);
				jy.push_back(j+1);
 				for(i=0; i<=0; ++i) {
			        	ix.clear();
			        	ix.push_back(nx-1);
			        	ix.push_back(i);
			        	ix.push_back(i+1);
					float qbf = buf(i,j,k)*invert;
			        	peak_check = qbf > buf(ix[0],jy[0],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[1])*invert;
					//if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[2])*invert;
					if(peak_check) {
	    					if(count < ml) {
							count++;
							peaks.push_back( Pixel(i, j, k, qbf) );
							if(count == ml-1) sort(peaks.begin(), peaks.end(), peakcmp);
						} else {
							if( qbf > (peaks.back()).value ) {
								//  do the switch and sort again
								//cout << qbf<<"   "<< (peaks.back()).value <<"   "<< (*peaks.begin()).value <<endl;
								peaks.pop_back();
								peaks.push_back( Pixel(i, j, k, qbf) );
								if(ml > 1) sort(peaks.begin(), peaks.end(), peakcmp);
							}
						}
					}
					}}}}}}}}}}}}}}}}}}}}}}}}}
				}
 				for(i=nx-1; i<=nx-1; ++i) {
			        	ix.clear();
			        	ix.push_back(i-1);
			        	ix.push_back(i);
			        	ix.push_back(0);
					float qbf = buf(i,j,k)*invert;
			        	peak_check = qbf > buf(ix[0],jy[0],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[1])*invert;
					//if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[2])*invert;
					if(peak_check) {
	    					if(count < ml) {
							count++;
							peaks.push_back( Pixel(i, j, k, qbf) );
							if(count == ml-1) sort(peaks.begin(), peaks.end(), peakcmp);
						} else {
							if( qbf > (peaks.back()).value ) {
								//  do the switch and sort again
								//cout << qbf<<"   "<< (peaks.back()).value <<"   "<< (*peaks.begin()).value <<endl;
								peaks.pop_back();
								peaks.push_back( Pixel(i, j, k, qbf) );
								if(ml > 1) sort(peaks.begin(), peaks.end(), peakcmp);
							}
						}
					}
					}}}}}}}}}}}}}}}}}}}}}}}}}
				}
			}
		}
	break;
/* 	case(5):  //looks ugly, but it is the best I can do,  PAP 03/07/08
 	int nu = buf.get_usize();
 	int nv = buf.get_vsize();
	vector<int> lu, mv;
		for(m=1; m<=nv-2; ++m) {
		 	mv.clear();
			mv.push_back(m-1);
			mv.push_back(m);
			mv.push_back(m+1);
		for(l=1; l<=nu-2; ++l) {
		 	lu.clear();
			lu.push_back(l-1);
			lu.push_back(l);
			lu.push_back(l+1);
		for(k=1; k<=nz-2; ++k) {
		 	kz.clear();
			kz.push_back(k-1);
			kz.push_back(k);
			kz.push_back(k+1);
 			for(j=1; j<=ny-2; ++j) {
				jy.clear();
				jy.push_back(j-1);
				jy.push_back(j);
				jy.push_back(j+1);
 				for(i=1; i<=nx-2; ++i) {
			        	ix.clear();
			        	ix.push_back(i-1);
			        	ix.push_back(i);
			        	ix.push_back(i+1);
					float qbf = buf(i,j,k,l,m)*invert;
			        	peak_check = qbf > buf(ix[0],jy[0],kz[0],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[0],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[0],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[0],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[0],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[0],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[0],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[0],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[0],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[1],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[1],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[1],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[1],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[1],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[1],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[1],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[1],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[1],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[2],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[2],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[2],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[2],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[2],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[2],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[2],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[2],lu[0],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[2],lu[0],mv[0])*invert;

					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[0],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[0],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[0],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[0],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[0],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[0],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[0],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[0],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[0],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[1],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[1],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[1],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[1],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[1],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[1],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[1],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[1],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[1],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[2],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[2],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[2],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[2],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[2],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[2],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[2],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[2],lu[1],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[2],lu[1],mv[0])*invert;

					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[0],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[0],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[0],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[0],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[0],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[0],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[0],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[0],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[0],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[1],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[1],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[1],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[1],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[1],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[1],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[1],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[1],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[1],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[2],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[2],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[2],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[2],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[2],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[2],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[2],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[2],lu[2],mv[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[2],lu[2],mv[0])*invert;


					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[0],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[0],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[0],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[0],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[0],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[0],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[0],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[0],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[0],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[1],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[1],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[1],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[1],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[1],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[1],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[1],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[1],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[1],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[2],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[2],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[2],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[2],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[2],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[2],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[2],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[2],lu[0],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[2],lu[0],mv[1])*invert;

					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[0],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[0],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[0],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[0],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[0],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[0],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[0],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[0],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[0],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[1],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[1],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[1],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[1],lu[1],mv[1])*invert;
					//if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[1],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[1],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[1],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[1],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[1],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[2],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[2],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[2],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[2],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[2],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[2],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[2],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[2],lu[1],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[2],lu[1],mv[1])*invert;

					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[0],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[0],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[0],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[0],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[0],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[0],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[0],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[0],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[0],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[1],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[1],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[1],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[1],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[1],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[1],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[1],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[1],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[1],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[2],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[2],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[2],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[2],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[2],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[2],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[2],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[2],lu[2],mv[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[2],lu[2],mv[1])*invert;


					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[0],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[0],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[0],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[0],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[0],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[0],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[0],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[0],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[0],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[1],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[1],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[1],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[1],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[1],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[1],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[1],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[1],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[1],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[2],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[2],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[2],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[2],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[2],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[2],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[2],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[2],lu[0],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[2],lu[0],mv[2])*invert;

					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[0],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[0],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[0],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[0],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[0],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[0],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[0],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[0],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[0],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[1],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[1],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[1],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[1],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[1],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[1],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[1],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[1],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[1],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[2],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[2],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[2],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[2],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[2],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[2],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[2],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[2],lu[1],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[2],lu[1],mv[2])*invert;

					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[0],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[0],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[0],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[0],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[0],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[0],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[0],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[0],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[0],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[1],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[1],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[1],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[1],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[1],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[1],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[1],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[1],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[1],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[2],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[2],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[2],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[2],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[2],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[2],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[2],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[2],lu[2],mv[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[2],lu[2],mv[2])*invert;
					if(peak_check) {
	    					if(count < ml) {
							count++;
							//peaks.push_back( Pixel(i, j, k, l, m, qbf) );
							if(count == ml-1) sort(peaks.begin(), peaks.end(), peakcmp);
						} else {
							if( qbf > (peaks.back()).value ) {
								//  do the switch and sort again
								//cout << qbf<<"   "<< (peaks.back()).value <<"   "<< (*peaks.begin()).value <<endl;
								peaks.pop_back();
								//peaks.push_back( Pixel(i, j, k, l, m, qbf) );
								if(ml > 1) sort(peaks.begin(), peaks.end(), peakcmp);
							}
						}
					}
					}}}}}}}}}}}}}}}}}}}}}}}}}}
					}}}}}}}}}}}}}}}}}}}}}}}}}}}
					}}}}}}}}}}}}}}}}}}}}}}}}}}}
					}}}}}}}}}}}}}}}}}}}}}}}}}}}
					}}}}}}}}}}}}}}}}}}}}}}}}}}
					}}}}}}}}}}}}}}}}}}}}}}}}}}}
					}}}}}}}}}}}}}}}}}}}}}}}}}}}
					}}}}}}}}}}}}}}}}}}}}}}}}}}}
					}}}}}}}}}}}}}}}}}}}}}}}}}}}
				}
			}
		}
		}
		}
		//  Add circular closure for x, y, and z directions: needed for circular ccf,
		//  should not have adverse impact on other code.  PAP 11/7/08
		for(m=1; m<=nv-2; ++m) {
		 	mv.clear();
			mv.push_back(m-1);
			mv.push_back(m);
			mv.push_back(m+1);
		for(l=1; l<=nu-2; ++l) {
		 	lu.clear();
			lu.push_back(l-1);
			lu.push_back(l);
			lu.push_back(l+1);
		for(k=1; k<=nz-2; ++k) {
		 	kz.clear();
			kz.push_back(k-1);
			kz.push_back(k);
			kz.push_back(k+1);
 			for(j=1; j<=ny-2; ++j) {
				jy.clear();
				jy.push_back(j-1);
				jy.push_back(j);
				jy.push_back(j+1);
 				for(i=0; i<=0; ++i) {
			        	ix.clear();
			        	ix.push_back(nx-1);
			        	ix.push_back(i);
			        	ix.push_back(i+1);
					float qbf = buf(i,j,k)*invert;
			        	peak_check = qbf > buf(ix[0],jy[0],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[1])*invert;
					//if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[2])*invert;
					if(peak_check) {
	    					if(count < ml) {
							count++;
							peaks.push_back( Pixel(i, j, k, qbf) );
							if(count == ml-1) sort(peaks.begin(), peaks.end(), peakcmp);
						} else {
							if( qbf > (peaks.back()).value ) {
								//  do the switch and sort again
								//cout << qbf<<"   "<< (peaks.back()).value <<"   "<< (*peaks.begin()).value <<endl;
								peaks.pop_back();
								peaks.push_back( Pixel(i, j, k, qbf) );
								if(ml > 1) sort(peaks.begin(), peaks.end(), peakcmp);
							}
						}
					}
					}}}}}}}}}}}}}}}}}}}}}}}}}
				}
 				for(i=nx-1; i<=nx-1; ++i) {
			        	ix.clear();
			        	ix.push_back(i-1);
			        	ix.push_back(i);
			        	ix.push_back(0);
					float qbf = buf(i,j,k)*invert;
			        	peak_check = qbf > buf(ix[0],jy[0],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[0])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[1])*invert;
					//if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[1])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[0],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[0],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[0],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[1],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[1],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[1],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[0],jy[2],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[1],jy[2],kz[2])*invert;
					if( peak_check ) {peak_check = qbf > buf(ix[2],jy[2],kz[2],3,3)*invert;
					if(peak_check) {
	    					if(count < ml) {
							count++;
							peaks.push_back( Pixel(i, j, k, qbf) );
							if(count == ml-1) sort(peaks.begin(), peaks.end(), peakcmp);
						} else {
							if( qbf > (peaks.back()).value ) {
								//  do the switch and sort again
								//cout << qbf<<"   "<< (peaks.back()).value <<"   "<< (*peaks.begin()).value <<endl;
								peaks.pop_back();
								peaks.push_back( Pixel(i, j, k, qbf) );
								if(ml > 1) sort(peaks.begin(), peaks.end(), peakcmp);
							}
						}
					}
					}}}}}}}}}}}}}}}}}}}}}}}}}
				}
			}
		}
		}
		}
	break;*/
	}
	// do we have a peak list yet?
	if (peaks.begin() != peaks.end()) {
	  // yes. sort it
	  sort(peaks.begin(), peaks.end(), peakcmp);

	  int count = 0;

	  float xval = (*peaks.begin()).value;
	  // loop over all peaks
	  for (vector<Pixel>::iterator it = peaks.begin(); it != peaks.end(); it++)  {
	    // current peak count
	    count++;
	    // is current peak count below max?
	    if(count <= ml) {
	      // yes, so append it
	      res.push_back((*it).value);
	      res.push_back(static_cast<float>((*it).x));

	      if(img_dim > 1) {
		res.push_back(static_cast<float>((*it).y));
		if(nz > 1) res.push_back(static_cast<float>((*it).z));
	      }

	      if(xval != 0.0) res.push_back((*it).value/xval);
	      else            res.push_back((*it).value);
	      res.push_back((*it).x-float(int(nx/2)));
	      if(img_dim >1) {
		res.push_back((*it).y-float(int(ny/2)));
		if(nz>1)   res.push_back((*it).z-float(nz/2));
	      }
	    }
	  }
	  res.insert(res.begin(),1,img_dim);
	} else {
	  // no peak list. build empty list
	  res.push_back(buf(0,0,0));
	  res.insert(res.begin(),1,0.0);
	}

	// return results list
	return res;
}

#define rdata(i,j,k) rdata[(i-1)+((j-1)+(k-1)*ny)*(size_t)nx]
#define X(i) X[i-1]
#define Y(j) Y[j-1]
#define Z(k) Z[k-1]
vector<float> EMData::phase_cog()
{
	vector<float> ph_cntog;
	int i=1,j=1,k=1;
	float C=0.f,S=0.f,P=0.f,F1=0.f,SNX;
	if (get_ndim()==1) {
		P = 8*atan(1.0f)/nx;
		for (i=1;i<=nx;i++) {
			C += cos(P * (i-1)) * rdata(i,j,k);
			S += sin(P * (i-1)) * rdata(i,j,k);
		}
		F1 = atan2(S,C);
		if (F1 < 0.0)  F1 += 8*atan(1.0f);
		SNX = F1/P +1.0f;
		SNX = SNX - ((nx/2)+1);
		ph_cntog.push_back(SNX);
#ifdef _WIN32
		ph_cntog.push_back((float)Util::round(SNX));
#else
		ph_cntog.push_back(round(SNX));
#endif //_WIN32
	} else if (get_ndim()==2)  {
#ifdef _WIN32
		float SNY;
		float T=0.0f;
		vector<float> X;
		X.resize(nx);
#else
		float SNY,X[nx],T=0.f;
#endif	//_WIN32
		for ( i=1;i<=nx;i++) X(i)=0.0;
                P = 8*atan(1.0f)/ny;
		for(j=1;j<=ny;j++) {
			T=0.f;
		        for(i=1;i<=nx;i++) {
				T += rdata(i,j,k);
		        	X(i)+=rdata(i,j,k);
			}
		        C += cos(P*(j-1))*T;
		        S += sin(P*(j-1))*T;
		}
		F1=atan2(S,C);
		if(F1<0.0)  F1 += 8*atan(1.0f);
		SNY = F1/P +1.0f;
		C=0.f;  S=0.f;
		P = 8*atan(1.0f)/nx;
		for(i=1;i<=nx;i++) {
			C += cos(P*(i-1))*X(i);
		        S += sin(P*(i-1))*X(i);
		}
	        F1=atan2(S,C);
		if(F1<0.0) F1 += 8*atan(1.0f);
		SNX = F1/P +1.0f;
		SNX = SNX - ((nx/2)+1);
		SNY = SNY - ((ny/2)+1);
		ph_cntog.push_back(SNX); ph_cntog.push_back(SNY);
#ifdef _WIN32
		ph_cntog.push_back((float)Util::round(SNX)); ph_cntog.push_back((float)Util::round(SNY));
#else
		ph_cntog.push_back(round(SNX)); ph_cntog.push_back(round(SNY));
#endif	//_WIN32
	} else {
#ifdef _WIN32
		float val=0.f,sum1=0.f, SNY,SNZ;
		vector<float> X;
		X.resize(nx);
		vector<float> Y;
		Y.resize(ny);
		vector<float> Z;
		Z.resize(nz);
#else
		float val=0.f, sum1=0.f, X[nx], Y[ny], Z[nz], SNY, SNZ;
#endif	//_WIN32
		 for (i=1;i<=nx;i++)  X(i)=0.0;
		 for (j=1;j<=ny;j++)  Y(j)=0.0;
		 for (k=1;k<=nz;k++)  Z(k)=0.0;
		 for(k=1;k<=nz;k++)  {
		 	for(j=1;j<=ny;j++) {
				sum1=0.f;
				for(i=1;i<=nx;i++)  {
					val = rdata(i,j,k);
					sum1 += val;
					X(i) += val;
				}
				Y(j) += sum1;
				Z(k) += sum1;
			}
		}
		P = 8*atan(1.0f)/nx;
		for (i=1;i<=nx;i++) {
			C += cos(P*(i-1))*X(i);
		        S += sin(P*(i-1))*X(i);
		}
		F1=atan2(S,C);
		if(F1<0.0) F1 += 8*atan(1.0f);
		SNX = F1/P +1.0f;
		C=0.f;  S=0.f;
		P = 8*atan(1.0f)/ny;
		for(j=1;j<=ny;j++) {
			C += cos(P*(j-1))*Y(j);
		        S += sin(P*(j-1))*Y(j);
		}
		F1=atan2(S,C);
		if(F1<0.0)  F1 += 8*atan(1.0f);
		SNY = F1/P +1.0f;
		C=0.f;  S=0.f;
		P = 8*atan(1.0f)/nz;
		for(k=1;k<=nz;k++) {
			C += cos(P*(k-1))*Z(k);
			S += sin(P*(k-1))*Z(k);
		}
		F1=atan2(S,C);
		if(F1<0.0)  F1 += 8*atan(1.0f);
		SNZ = F1/P +1.0f;
		SNX = SNX - ((nx/2)+1);
		SNY = SNY - ((ny/2)+1);
		SNZ = SNZ - ((nz/2)+1);
		ph_cntog.push_back(SNX); ph_cntog.push_back(SNY); ph_cntog.push_back(SNZ);
#ifdef _WIN32
		ph_cntog.push_back((float)Util::round(SNX)); ph_cntog.push_back((float)Util::round(SNY)); ph_cntog.push_back((float)Util::round(SNZ));
#else
		ph_cntog.push_back(round(SNX)); ph_cntog.push_back(round(SNY));ph_cntog.push_back(round(SNZ));
#endif
	}
	return ph_cntog;
}
#undef rdata
#undef X
#undef Y
#undef Z

#define avagadro (6.023*(double)pow(10.0,23.0))
#define density_protein (1.36)
#define R (0.61803399f)
#define C (1.f-R)
float EMData::find_3d_threshold(float mass, float pixel_size)
{
	/* Exception Handle */
	if(get_ndim()!=3)
		throw ImageDimensionException("The image should be 3D");
	/* ===============================================================*/

	/* Calculation of the volume of the voxels */
	float density_1_mole, vol_1_mole, vol_angstrom;
	int  vol_voxels;
	density_1_mole = static_cast<float>( (mass*1000.0f)/avagadro );
	vol_1_mole =  static_cast<float>( density_1_mole/density_protein );
	vol_angstrom =  static_cast<float>( vol_1_mole*(double)pow((double)pow(10.0,8),3) );
	vol_voxels = static_cast<int> (vol_angstrom/(double)pow(pixel_size,3));
	/* ===============================================================*/


	float thr1 = get_attr("maximum");
	float thr3 = get_attr("minimum");
	float thr2 = (thr1-thr3)/2 + thr3;
	size_t size = (size_t)nx*ny*nz;
	float x0 = thr1,x3 = thr3,x1,x2,THR=0;

		int ILE = std::min(nx*ny*nx,std::max(1,vol_voxels));

	if (abs(thr3-thr2)>abs(thr2-thr1)) {
		x1=thr2;
		x2=thr2+C*(thr3-thr2);
	} else {
		x2=thr2;
		x1=thr2-C*(thr2-thr1);
	}

	int cnt1=0,cnt2=0;
	for (size_t i=0;i<size;++i) {
		if(rdata[i]>=x1)  cnt1++;
		if(rdata[i]>=x2)  cnt2++;
	}
	float LF1 = static_cast<float>( cnt1 - ILE );
	float F1 = LF1*LF1;
	float LF2 = static_cast<float>( cnt2 - ILE );
	float F2 = LF2*LF2;

	while ((LF1 != 0 || LF2 != 0) && (fabs(LF1-LF2) >= 1.f) && (abs(x1-x2) > (double)pow(10.0,-5) && abs(x1-x3) > (double)pow(10.0,-5) && abs(x2-x3) > (double)pow(10.0,-5)))
	{
		if(F2 < F1) {
			x0=x1;
			x1=x2;
			x2 = R*x1 + C*x3;
			F1=F2;
			int cnt=0;
			for(size_t i=0;i<size;++i)
				if(rdata[i]>=x2)
					cnt++;
			LF2 = static_cast<float>( cnt - ILE );
			F2 = LF2*LF2;
		} else {
			x3=x2;
			x2=x1;
			x1=R*x2 + C*x0;
			F2=F1;
			int cnt=0;
			for(size_t i=0;i<size;++i)
				if(rdata[i]>=x1)
					cnt++;
			LF1 = static_cast<float>( cnt - ILE );
			F1 = LF1*LF1;
		}
	}

	if(F1 < F2) {
		ILE = static_cast<int> (LF1 + ILE);
		THR = x1;
	} else {
		ILE = static_cast<int> (LF2 + ILE);
		THR = x2;
	}
	return THR;

}
#undef avagadro
#undef density_protein
#undef R
#undef C


// reworked peak_ccf uses max queue length for peak objects, i.e. lowest
//    peaks are deleted if queue length is exceeded and a new peak is inserted
//    instead.


vector<float> EMData::peak_ccf(float hf_p)
{

  // cout << "peak ccf starting up" << endl;

  EMData & buf = *this;
  vector<Pixel> peaks;
  int half=int(hf_p);
  float hf_p2 = hf_p*hf_p;
  int i,j;
  int i__1,i__2;
  int j__1,j__2;
  vector<float>res;
  int nx = buf.get_xsize()-half;
  int ny = buf.get_ysize()-half;
  // iterate over image
  for(i=half; i<=nx; ++i) {
    // static assignment so we don't have to re-evaluate
    i__1 = i-1;
    i__2 = i+1;
    for (j=half;j<=ny;++j) {
      j__1 = j-1;
      j__2 = j+1;

      if((buf(i,j)>0.0f)&&buf(i,j)>buf(i,j__1)) {
	if(buf(i,j)>buf(i,j__2)) {
	  if(buf(i,j)>buf(i__1,j)) {
	    if(buf(i,j)>buf(i__2,j)) {
	      if(buf(i,j)>buf(i__1,j__1)) {
		if((buf(i,j))> buf(i__1,j__2)) {
		  if(buf(i,j)>buf(i__2,j__1)) {
		    if(buf(i,j)> buf(i__2,j__2)) {

		      // found a peak
		      // empty list?
		      if (peaks.size()==0) {
			// yes, so just push the peak onto the list
			peaks.push_back(Pixel(i,j,0,buf(i,j)));

		      } else {
			// not empty list. check neighbourhood for peaks
			// logical not in the name is awkward. renamed to overlap
			bool overlap = false;
			//int  size = peaks.size();

			// list of peaks to be deleted, if the current peak is the largest (see below).
			//    list contains iterators to the original list, which will have to be processed
			//    back to front (i.e. LIFO: stl::stack)
			std::stack <vector<Pixel>::iterator> delete_stack;

			// loop over all peaks found so far. this would be nicer with iterators
			for (vector<Pixel>::iterator it=peaks.begin();it!=peaks.end();++it) {
			// for ( int kk= 0; kk< size; kk++) {
			//  vector<Pixel>::iterator it = peaks.begin()+kk;

			  // calc L2 distance
			  float radius=((*it).x-float(i))*((*it).x-float(i))+((*it).y-float(j))*((*it).y-float(j));
			  if (radius <= hf_p2 ) {
			    // peaks overlap
			    if( buf(i,j) > (*it).value) {
			      // this peak (indexed by (i,j)) is larger, mark the old for deletion
			      //    however, we have to be careful. if there is a larger peak within the vicinity of
			      //    the new one, this new peak is not marked as such, and the deletion of prior low
			      //    peaks should not continued. to make sure this deletion does not happen, we have
			      //    to make sure we cycle through all peaks within the vicinity, and only delete smaller
			      //    peaks if this new one is the largest in the vicinity.
			      delete_stack.push(it);

			      //(*it).x = -half; // this marks entry to be deleted, since it's smaller than the new one


			    } else {
			      overlap = true;
			      // old peak is larger, ignore this one. since it's enough to know there is some peak larger
			      //    than this one, we can break out of the peak list loop, instead of continuing.
			      break;
			    }
			  }
			}

			// check whether we need to delete anything. this is marked by the flag overlap == false
			// loop over all peaks and clean out redundant ones
			if (false == overlap) {
			  vector<Pixel>::iterator delete_iterator;
			  while (!delete_stack.empty()) {
			    // pop empties the stack from the back. since we are dealing with iterators, we need to delete
			    //    from the back, so as to keep the rest stack intact upon deletion.
			    delete_iterator = delete_stack.top();
			    peaks.erase(delete_iterator);
			    delete_stack.pop();
			  }
			  // before pushing the peak, we need to check whether max queue length is exceeded and delete
			  //     peaks if necessary.
			  // XXX: remove hardcoded value!
			  if (! (peaks.size() < 2000 )) {

			    //cout << ".";
			    // we need to delete a peak first.
			    // - resort list to get lowest peak at the back
			    sort(peaks.begin(), peaks.end(), peakcmp);

			    // - remove lowest peak
			    peaks.pop_back();
			  }

			  // push the new peak onto the list of peaks
			  peaks.push_back(Pixel(i,j,0,buf(i,j)));
			  //cout << "done." << endl;

			} else {
			  // this peak too small and is ignored, so delete_list is ignored as well. make sure delete_list
			  //    is empty. probably redundant because of scope, but better safe than sorry.....
			  while (!delete_stack.empty()) delete_stack.pop();
			}
		      }
		    }
		  }}}}}}}
    }
  }

  // we have peaks, so build a results vector.
  if(peaks.size()>0) {
    // sort peaks by size
    sort(peaks.begin(),peaks.end(), peakcmp);
    // and push all peaks to the results vector
    for (vector<Pixel>::iterator it = peaks.begin(); it != peaks.end(); it++) {
      // XXX: this format is necessary for Boost to work???
      res.push_back((*it).value);
      res.push_back(static_cast<float>((*it).x));
      res.push_back(static_cast<float>((*it).y));
    }
  } else {
    // only one or zero (?) entries
    res.push_back(buf(0,0,0));
    res.insert(res.begin(),1,0.0);
  }
  return res;
}

EMData* EMData::get_pow(float n_pow)
{
	EMData* buf_new = this->copy_head();
	float *in  = this->get_data();
	float *out = buf_new->get_data();
	for(size_t i=0; i<(size_t)nx*ny*nz; ++i) out[i] = pow(in[i],n_pow);
	return buf_new;
}

EMData* EMData::conjg()
{
	if(this->is_complex()) {
		EMData* buf_new = this->copy_head();
 		float *in  = this->get_data();
		float *out = buf_new->get_data();
		for(size_t i=0; i<(size_t)nx*ny*nz; i+=2) {out[i] = in[i]; out[i+1] = -in[i+1];}
		return buf_new;
	} else throw ImageFormatException("image has to be complex");
}

EMData* EMData::delete_disconnected_regions(int ix, int iy, int iz) {
	if (3 != get_ndim())
		throw ImageDimensionException("delete_disconnected_regions needs a 3-D image.");
	if (is_complex())
		throw ImageFormatException("delete_disconnected_regions requires a real image");
	if ((*this)(ix+nx/2,iy+ny/2,iz+nz/2) == 0)
		throw ImageDimensionException("delete_disconnected_regions starting point is zero.");

	EMData* result = this->copy_head();
	result->to_zero();
	(*result)(ix+nx/2,iy+ny/2,iz+nz/2) = (*this)(ix+nx/2,iy+ny/2,iz+nz/2);
	bool kpt = true;
	//cout << "  delete   "<<(*result)(ix+nx/2,iy+ny/2,iz+nz/2)<<endl;
	while(kpt) {
		kpt = false;
		for (int cz = 1; cz < nz-1; cz++) {
			for (int cy = 1; cy < ny-1; cy++) {
				for (int cx = 1; cx < nx-1; cx++) {
					if((*result)(cx,cy,cz) == 1) {
						for (int lz = -1; lz <= 1; lz++) {
							for (int ly = -1; ly <= 1; ly++) {
								for (int lx = -1; lx <= 1; lx++) {
									if(((*this)(cx+lx,cy+ly,cz+lz) == 1) && ((*result)(cx+lx,cy+ly,cz+lz) == 0))  {
										(*result)(cx+lx,cy+ly,cz+lz) = 1;
										kpt = true;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	result->update();
	return result;
}

EMData* EMData::helicise(float pixel_size, float dp, float dphi, float section_use, float radius, float minrad) {
	if(3 != get_ndim())               throw ImageDimensionException("helicise needs a 3-D image.");
	if(is_complex())                  throw ImageFormatException("helicise requires a real image");
	if(int(section_use*nz+0.5)>=nz-2) throw ImageFormatException("Reduce section used for helicise");

	EMData* result = this->copy_head();
	result->to_zero();


	int nxc = nx/2;
	int nyc = ny/2;
	int nzc = nz/2;
	//  calculations are done in Angstroms
	float volen = nz*pixel_size;
	float nzcp  = nzc*pixel_size;
	float sectl = nz*pixel_size*section_use;  //In Angstroms
	float nb = nzcp - sectl/2.0f;
	float ne = nzcp + sectl/2.0f;
	int numst = int( nz*pixel_size/dp );  // A
	int numri = int(sectl/dp);             // pix
	if(numri < 1)   throw ImageFormatException("Increase section used for helicise");

	float r2, ir;
	if(radius < 0.0f) r2 = (float)((nxc-1)*(nxc-1));
	else r2 = radius*radius;
	if(minrad < 0.0f) ir = 0.0f;
	else ir = minrad*minrad;

	for (int k = 0; k<nz; k++) {
		int nq = 0;
		for (int ist = 0; ist<numst; ist++) {
			float z = k*pixel_size + ist*dp;
			float phi = ist*dphi;
			if( z >= volen ) {
				z = k*pixel_size + (ist-numst)*dp;
				phi = (ist-numst)*dphi;
			}
			float ca = cos(phi*(float)deg_rad);
			float sa = sin(phi*(float)deg_rad);
			if((z >= nb) && (z <= ne )) {
				nq++;
				if( nq > numri ) break;
				float zz = z/pixel_size;
				//cout <<" zz  "<<zz<<"  k  "<<k<<"  phi  "<<phi<<endl;
				for (int j=0; j<ny; j++) {
					int jy = j - nyc;
					int jj = jy*jy;
					for (int i=0; i<nx; i++) {
						int ix = i - nxc;
						float d2 = float((ix*ix + jj));
						if(d2 <= r2 && d2>=ir) {
							float xx =  ix*ca + jy*sa + nxc;
							float yy = -ix*sa + jy*ca + nyc;


	//  Do tri-linear interpolation
	int IOX = int(xx);
	int IOY = int(yy);
	int IOZ = int(zz);

	int IOXp1 = std::min( nx-1 ,IOX+1);

	int IOYp1 = std::min( ny-1 ,IOY+1);

	int IOZp1 = std::min( nz-1 ,IOZ+1);

	float dx = xx-IOX;
	float dy = yy-IOY;
	float dz = zz-IOZ;

	float a1 = (*this)(IOX,IOY,IOZ);
	float a2 = (*this)(IOXp1,IOY,IOZ) - (*this)(IOX,IOY,IOZ);
	float a3 = (*this)(IOX,IOYp1,IOZ) - (*this)(IOX,IOY,IOZ);
	float a4 = (*this)(IOX,IOY,IOZp1) - (*this)(IOX,IOY,IOZ);
	float a5 = (*this)(IOX,IOY,IOZ) - (*this)(IOXp1,IOY,IOZ) - (*this)(IOX,IOYp1,IOZ) + (*this)(IOXp1,IOYp1,IOZ);
	float a6 = (*this)(IOX,IOY,IOZ) - (*this)(IOXp1,IOY,IOZ) - (*this)(IOX,IOY,IOZp1) + (*this)(IOXp1,IOY,IOZp1);
	float a7 = (*this)(IOX,IOY,IOZ) - (*this)(IOX,IOYp1,IOZ) - (*this)(IOX,IOY,IOZp1) + (*this)(IOX,IOYp1,IOZp1);
	float a8 = (*this)(IOXp1,IOY,IOZ) + (*this)(IOX,IOYp1,IOZ)+ (*this)(IOX,IOY,IOZp1)
			- (*this)(IOX,IOY,IOZ)- (*this)(IOXp1,IOYp1,IOZ) - (*this)(IOXp1,IOY,IOZp1)
			- (*this)(IOX,IOYp1,IOZp1) + (*this)(IOXp1,IOYp1,IOZp1);



							(*result)(i,j,k) += a1 + dz*(a4 + a6*dx + (a7 + a8*dx)*dy) + a3*dy + dx*(a2 + a5*dy);
						}
					}
				}
			}
		}
		if(nq < numri) throw InvalidValueException(nq, "Helicise: incorrect number of repeats encoutered.");
	}
	for (int k = 0; k<nz; k++) for (int j = 0; j<ny; j++) for (int i = 0; i<nx; i++) (*result)(i,j,k) /= numst ;

	result->update();
	return result;
}



EMData* EMData::helicise_grid(float pixel_size, float dp, float dphi, Util::KaiserBessel& kb, float section_use, float radius, float minrad) {
	std::cout<<"111111"<<std::endl;
	if (3 != get_ndim())
		throw ImageDimensionException("helicise needs a 3-D image.");
	if (is_complex())
		throw ImageFormatException("helicise requires a real image");
	//begin griding
	//if (scale_input == 0.0f) scale_input = 1.0f;
	float  scale = 0.5f;//*scale_input;

	
	int nxn = nx/2; int nyn = ny/2; int nzn = nz/2;

	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
	EMData* ret = this->copy_head();
	ret->set_size(nxn, std::max(nyn,1), std::max(nzn,1));
	ret->to_zero();  //we will leave margins zeroed.

	// center of big image,
	int xc = nxn;
	int ixs = nxn%2;  // extra shift on account of odd-sized images
	int yc = nyn;
	int iys = nyn%2;
	int zc = nzn;
	int izs = nzn%2;
	// center of small image
	int xcn = nxn/2;
	int ycn = nyn/2;
	int zcn = nzn/2;
	// shifted center for rotation
	float shiftxc = xcn; // + delx;
	float shiftyc = ycn; // + dely;
	float shiftzc = zcn; // + delz;
	// bounds if origin at center
	float zmin = -nz/2.0f;
	float ymin = -ny/2.0f;
	float xmin = -nx/2.0f;
	float zmax = -zmin;
	float ymax = -ymin;
	float xmax = -xmin;
	if (0 == nx%2) xmax--;
	if (0 == ny%2) ymax--;
	if (0 == nz%2) zmax--;

	float* data = this->get_data();

	
	// rotation matrix (the transpose is used in the loop to get (xold,yold,zold)):
	 
	//float a13 = -0.0f;	float a23 =  0.0f;
	//float a31 =  0.0f;    float a32 =  0.0f;        float a33 =  1.0f;
		
	//end gridding

	
	int nyc = nyn/2;
	int nxc = nxn/2;
	int nb = int(nzn*(1.0f - section_use)/2.);
	int ne = nzn - nb -1;
	int numst = int(nzn*section_use*pixel_size/dp);
	// how many steps needed total, fewer will be used, only those that fall between nb and ne
	int nst = int(nzn*pixel_size/dp);
	float r2, ir;
	if(radius < 0.0f) r2 = (float)((nxc-1)*(nxc-1));
	else r2 = radius*radius;
	if(minrad < 0.0f) ir = 0.0f;
	else ir = minrad*minrad;
	
	for (int k = 0; k<nzn; k++) {
		for (int j = 0; j<nyn; j++) {
			int jy = j - nyc;
			int jj = jy*jy;
			for (int i = 0; i<nxn; i++) {
				int ix = i - nxc;
				float d2 = (float)(ix*ix + jj);
				if(d2 <= r2 && d2>=ir) {
					int nq = 0;
					for ( int ist = -nst; ist <= nst; ist++) {
						float zold = (k*pixel_size + ist*dp)/pixel_size;
						int IOZ = int(zold);
						if(IOZ >= nb && IOZ <= ne) {
						
							float cphi = ist*dphi*(float)deg_rad;
							float ca = cos(cphi);
							float sa = sin(cphi);
							
							float xold = ix*ca - jy*sa + nxc;
							float yold = ix*sa + jy*ca + nyc;
							
							float xold_big = (xold-shiftxc)/scale - ixs + xc;
							float yold_big = (yold-shiftyc)/scale - iys + yc;
							float zold_big = (zold-shiftzc)/scale - izs + zc;
							
							/*float a11 =  ca; float a12 =  sa;
							float a21 = -sa; float a22 = ca;
							
							float z = (zold - shiftzc)/scale;
							float zco1 = a31*z+xc;
							float zco2 = a32*z+yc;
							float zco3 = a33*z+zc;
														
							float y = (float(j) - shiftyc)/scale;
							float yco1 = zco1+a21*y;
							float yco2 = zco2+a22*y;
							float yco3 = zco3+a23*y;
							
							float x = (float(i) - shiftxc)/scale;
							float xold_big = yco1+a11*x-ixs; //have to add the fraction on account on odd-sized images for which Fourier zero-padding changes the center location
							float yold_big = yco2+a12*x-iys;
							float zold_big = yco3+a13*x-izs;*/
							
												
							nq++;
							
								
							(*ret)(i,j,k) += Util::get_pixel_conv_new(nx, ny, nz, xold_big, yold_big, zold_big, data, kb);
							
							
							if(nq == numst) break;
						}
					}
					if(nq != numst)
						throw InvalidValueException(nq, "Helicise: incorrect number of repeats encoutered.");
				}
			}
		}
	}
	
	for (int k = 0; k<nzn; k++) for (int j = 0; j<nyn; j++) for (int i = 0; i<nxn; i++) (*ret)(i,j,k) /= numst ;
	set_array_offsets(saved_offsets);
	ret->update();
	return ret;
}


/*
Purpose: Depad and remove FT extension from a real image.
Method: Depad and remove FT extension from a real image.
return new image.
Input: f real n-dimensional image
Output: depadded image
 */
void EMData::depad() {
	if (is_complex())
		throw ImageFormatException("Depadding of complex images not supported");
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
	int npad = attr_dict["npad"];
	if (0 == npad) npad = 1;
	int offset = is_fftodd() ? 1 : 2;
	int nxold = (nx - offset)/npad;
	int nyold = std::max<int>(ny/npad, 1);
	int nzold = std::max<int>(nz/npad, 1);
	int xstart = 0, ystart = 0, zstart = 0;
	if( npad > 1) {
		xstart = (nx - offset - nxold)/2 + nxold%2;
		if(ny > 1) {
			ystart = (ny - nyold)/2 + nyold%2;
			if(nz > 1) {
				zstart = (nz - nzold)/2 + nzold%2;
			}
		}
	}
	int bytes = nxold*sizeof(float);
	float* dest = get_data();
	for (int iz=0; iz < nzold; iz++) {
		for (int iy = 0; iy < nyold; iy++) {
			memmove(dest, &(*this)(xstart,iy+ystart,iz+zstart), bytes);
			dest += nxold;
		}
	}
	set_size(nxold, nyold, nzold);
	set_attr("npad", 1);
	set_fftpad(false);
	set_fftodd(false);
	set_complex(false);
	if(ny==1 && nz==1) set_complex_x(false);
	set_array_offsets(saved_offsets);
	update();
	EXITFUNC;
}

/*
Purpose: Depad and remove FT extension from a real image.
Method: Depad and remove FT extension from a real image.
return new image.
Input: f real n-dimensional image
Output: depadded image
 */
void EMData::depad_corner() {
	if(is_complex())
		throw ImageFormatException("Depadding of complex images not allowed");
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
	int npad = attr_dict["npad"];
	if(0 == npad) npad = 1;
	int offset = is_fftodd() ? 1 : 2;
	int nxold = (nx - offset)/npad;
	int nyold = std::max<int>(ny/npad, 1);
	int nzold = std::max<int>(nz/npad, 1);
	size_t bytes = nxold*sizeof(float);
	float* dest = get_data();
	for (int iz=0; iz < nzold; iz++) {
		for (int iy = 0; iy < nyold; iy++) {
			memmove(dest, &(*this)(0,iy,iz), bytes);
			dest += nxold;
		}
	}
	set_size(nxold, nyold, nzold);
	set_attr("npad", 1);
	set_fftpad(false);
	set_fftodd(false);
	set_complex(false);
	if(ny==1 && nz==1) set_complex_x(false);
	set_array_offsets(saved_offsets);
	update();
	EXITFUNC;
}


// calculate circumference of the surrounding 1 pixel.
float circumference( EMData* emdata, int npixel )
{
	int nx = emdata->get_xsize();
	int ny = emdata->get_ysize();
	int nz = emdata->get_zsize();

        float* data = emdata->get_data();
	if( ny==1 && nz==1 ) {
        	// 1d case
        	float sumf=0.0;
        	int   sumn=0;
        	for( int i=0; i < npixel; ++i ) {
        		sumf += data[i];
        		sumf += data[nx-1-i];
        		sumn += 2;
        	}
        	return sumf/sumn;
        }

        if( nz==1 ) {
        	float sumf=0.0;
        	int   sumn=0;
        	int   id=0;
        	for( int iy=0; iy < ny; ++iy ) {
        		for( int ix=0; ix < nx; ++ix ) {
        			if( iy<npixel || iy>ny-1-npixel || ix<npixel || ix>nx-1-npixel ) {
        			    sumf += data[id];
        			    sumn += 1;
        			}
        			id++;
        		}
        	}

        	Assert( id==nx*ny  );
        	Assert( sumn == nx*ny - (nx-2*npixel)*(ny-2*npixel) );
        	return sumf/sumn;
        }

        // 3d cases;

        float sumf = 0.0;
        size_t   sumn = 0;
        size_t   id = 0;
        for( int iz=0; iz < nz; ++iz) {
        	for( int iy=0; iy < ny; ++iy) {
        		for( int ix=0; ix < nx; ++ix ) {
        			if( iz<npixel||iz>nz-1-npixel||iy<npixel||iy>ny-1-npixel||ix<npixel||ix>nx-1-npixel) {
        				sumf += data[id];
        				sumn += 1;
        			}
        			id++;
        		}
        	}
        }


        Assert( id==(size_t)nx*ny*nz);
        Assert( sumn==(size_t)nx*ny*nz-(size_t)(nx-2*npixel)*(ny-2*npixel)*(nz-2*npixel) );
        return sumf/sumn;
}
/*
Purpose: Create a new [normalized] [zero-padded] image.
Method: Normalize, pad with zero or circumference, extend for fft,
return new image.
Input: real n-dimensional image
flag specify normalize, pad, and/or extend
Output: zero-padded, ft-extended, normalized input image
 */
EMData* EMData::norm_pad(bool donorm, int npad, int valtype) {
	if (this->is_complex())
		throw ImageFormatException("Padding of complex images not supported");
	int nx = this->get_xsize();
	int ny = this->get_ysize();
	int nz = this->get_zsize();
	float mean = 0., stddev = 1.;
	if(donorm) { // Normalization requested
		mean = this->get_attr("mean");
		stddev = this->get_attr("sigma");
	}
	// sanity check
	if (npad < 1) npad = 1;
	int nxpad = npad*nx;
	int nypad = npad*ny;
	int nzpad = npad*nz;
	if (1 == ny) {
		// 1-d image, don't want to pad along y or z
		// Also, assuming that we can't have an image sized as nx=5, ny=1, nz=5.
		nypad = ny;
		nzpad = nz;
	} else if (nz == 1) {
		// 2-d image, don't want to pad along z
		nzpad = nz;
	}
	size_t bytes;
	size_t offset;
	// Not currently fft-extended, so we want to extend for ffts
	offset = 2 - nxpad%2;
	bytes = nx*sizeof(float);
	EMData* fpimage = copy_head();
	fpimage->set_size(nxpad+offset, nypad, nzpad);
	int xstart = 0, ystart = 0, zstart = 0;
	if( npad > 1) {
		if( valtype==0 ) {
			fpimage->to_zero();
		} else {
			float val = circumference(this, 1);
			float* data = fpimage->get_data();
			int nxyz = (nxpad+offset)*nypad*nzpad;
			for( int i=0; i < nxyz; ++i )  data[i] = val;
		}

		xstart = (nxpad - nx)/2 + nx%2;
		if(ny > 1) {
			ystart = (nypad - ny)/2 + ny%2;
			if(nz > 1) {
				zstart = (nzpad - nz)/2 + nz%2;
			}
		}
	} else fpimage->to_zero();


	vector<int> saved_offsets = this->get_array_offsets();
	this->set_array_offsets( 0, 0, 0 );
	for (int iz = 0; iz < nz; iz++) {
		for (int iy = 0; iy < ny; iy++) {
			memcpy(&(*fpimage)(xstart,iy+ystart,iz+zstart), &(*this)(0,iy,iz), bytes);
		}
	}
	this->set_array_offsets( saved_offsets );


	//  Perform normalization (only on the non-zero section of the image)
	if (donorm) { // Normalization requested
		for (int iz = zstart; iz < nz+zstart; iz++)
			for (int iy = ystart; iy < ny+ystart; iy++)
				for (int ix = xstart; ix < nx+xstart; ix++)
					(*fpimage)(ix,iy,iz) = ((*fpimage)(ix,iy,iz)-mean)/stddev;
	}

	fpimage->set_fftpad(true);
	fpimage->set_attr("npad", npad);
	if (offset == 1) fpimage->set_fftodd(true);
	else             fpimage->set_fftodd(false);
	return fpimage;
}

void EMData::center_origin()
{
	ENTERFUNC;
	if (is_complex()) {
		LOGERR("Real image expected. Input image is complex.");
		throw ImageFormatException("Real image expected. Input image is complex.");
	}
	for (int iz = 0; iz < nz; iz++) {
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				// next line multiplies by +/- 1
				(*this)(ix,iy,iz) *= -2*((ix+iy+iz)%2) + 1;
			}
		}
	}
	update();
	EXITFUNC;
}

void EMData::center_origin_yz()
{
	ENTERFUNC;
	if (is_complex()) {
		LOGERR("Real image expected. Input image is complex.");
		throw ImageFormatException("Real image expected. Input image is complex.");
	}
	for (int iz = 0; iz < nz; iz++) {
		for (int iy = (iz+1)%2; iy < ny; iy+=2) {
			for (int ix = 0; ix < nx; ix++) {
				(*this)(ix,iy,iz) *= -1;
			}
		}
	}
	update();
	EXITFUNC;
}

void EMData::center_origin_fft()
{
	ENTERFUNC;
	if (!is_complex()) {
		LOGERR("complex image expected. Input image is real image.");
		throw ImageFormatException("complex image expected. Input image is real image.");
	}

	if (!is_ri()) {
		LOGWARN("Only RI should be used. ");
	}
	vector<int> saved_offsets = get_array_offsets();
	// iz in [1,nz], iy in [1,ny], ix in [0,nx/2], but nx comes in as extended and is the same for odd
	//                                                 and even, so we can ignore the difference...
	//                         in short, as nx is extended, it should be  ix in [0,(nx-2)/2],  corrected PAP 05/20
	set_array_offsets(0,1,1);
	int nxc = nx/2;

	if (is_fftodd()) {

//cout<<"  center_origin_fft odd volume dimensions  "<<nxc<<"   "<<ny<<"  "<<nz<<endl;
		for (int iz = 1; iz <= nz; iz++) {
			for (int iy = 1; iy <= ny; iy++) {
				for (int ix = 0; ix < nxc; ix++) {
					cmplx(ix,iy,iz) *= float(-2*((ix+iy+iz)%2) + 1);
					float temp = float(iz-1+iy-1+ix)/float(ny)*M_PI;
					complex<float> temp2 = complex<float>(cos(temp), -sin(temp));
					cmplx(ix,iy,iz) *= temp2;
				}
			}
		}
	} else {

//cout<<"  center_origin_fft dimensions  "<<nxc<<"   "<<ny<<"  "<<nz<<endl;
		for (int iz = 1; iz <= nz; iz++) {
			for (int iy = 1; iy <= ny; iy++) {
				for (int ix = 0; ix < nxc; ix++) {
					// next line multiplies by +/- 1
					cmplx(ix,iy,iz) *= float(-2*((ix+iy+iz)%2) + 1);
				}
			}
		}
	}
	set_array_offsets(saved_offsets);
	update();
	EXITFUNC;
}


#define  fint(i,j,k)  fint[(i-1) + ((j-1) + (k-1)*ny)*(size_t)lsd]
#define  fout(i,j,k)  fout[(i-1) + ((j-1) + (k-1)*nyn)*(size_t)lsdn]
EMData *EMData::FourInterpol(int nxn, int nyn, int nzn, bool RetReal, bool normalize) {

	int lsd, lsdn, inx, iny, inz, nxo;
	int i, j, k;
	EMData *temp_ft = new EMData();
	if (is_complex()) {
		temp_ft = this->copy();
		lsd = nx;
		nxo = lsd - 2 + this->is_fftodd();
	} else {
		//  do out of place ft
		temp_ft = do_fft();
		lsd = nx + 2 - nx%2;
		nxo = nx;;
	}
	lsdn = nxn + 2 - nxn%2;
	if(nxn<nxo || nyn<ny || nzn<nz)	throw ImageDimensionException("Cannot reduce the image size");
	EMData *ret = this->copy_head();
	ret->set_size(lsdn, nyn, nzn);
	ret->to_zero();
	float *fout = ret->get_data();
	float *fint = temp_ft->get_data();
	if( normalize ) {
		float  anorm = (float) nxn* (float) nyn* (float) nzn/(float) nxo/ (float) ny/ (float) nz;
		for (i = 0; i < lsd*ny*nz; i++)  fint[i] *= anorm;
	}
	inx = nxn-nxo; iny = nyn - ny; inz = nzn - nz;
	for (k=1; k<=nz/2+1; k++) for (j=1; j<=ny/2+1; j++) for (i=1; i<=lsd; i++) fout(i,j,k) = fint(i,j,k);
	if(nyn>1) {
	//cout << "  " <<nxn<<"  " <<nyn<<" A " <<nzn<<endl;
		for (k=1; k<=nz/2+1; k++) for (j=ny/2+2+iny; j<=nyn; j++) for (i=1; i<=lsd; i++) fout(i,j,k) = fint(i,j-iny,k);
		if(nzn>1) {
			for (k=nz/2+2+inz; k<=nzn; k++) {
				for (j=1; j<=ny/2+1; j++) {
					for (i=1; i<=lsd; i++) {
						fout(i,j,k)=fint(i,j,k-inz);
					}
				}
				for (j=ny/2+2+iny; j<=nyn; j++) {
					for (i=1; i<=lsd; i++) {
						fout(i,j,k)=fint(i,j-iny,k-inz);
					}
				}
			}
		}
	}
//  TO KEEP EXACT VALUES ON THE ORIGINAL GRID ONE SHOULD USE
//  SQ2     = 2.0. HOWEVER, TOTAL ENERGY WILL NOT BE CONSERVED
	float  sq2 = 1.0f/std::sqrt(2.0f);
//       WEIGHTING FACTOR USED FOR EVEN NSAM. REQUIRED SINCE ADDING ZERO FOR
//       INTERPOLATION WILL INTRODUCE A COMPLEX CONJUGATE FOR NSAM/2'TH
//       ELEMENT.
	if(nxo%2 == 0 && inx !=0) {
		for (k=1; k<=nzn; k++) {
			for (j=1; j<=nyn; j++) {
				fout(nxo+1,j,k) *= sq2;
				fout(nxo+2,j,k) *= sq2;
			}
		}
		if(nyn>1) {
			for (k=1; k<=nzn; k++) {
			  for (i=1; i<=lsd; i++) {
			    fout(i,ny/2+1+iny,k) = sq2*fout(i,ny/2+1,k);
			    fout(i,ny/2+1,k) *= sq2;
			  }
			}
			if(nzn>1) {
				for (j=1; j<=nyn; j++) {
					for (i=1; i<=lsd; i++) {
						fout(i,j,nz/2+1+inz) = sq2*fout(i,j,nz/2+1);
						fout(i,j,nz/2+1) *= sq2;
					}
				}
			}
		}
	}
	ret->set_complex(true);
/*
//  For padding from odd to even dimension additional shift by 1 pixel is necessary.
	float  xshift = 0.f, yshift = 0.f, zshift = 0.f;
	int nyn2, nzn2;
	if(nxn > nx && nx%2 == 1)  xshift = 1.0f;
	if(ny > 1) {
		if(nyn > ny && ny%2 == 1)  yshift = 1.0f;
		nyn2 = nyn/2;
		if(nz > 1) {
			if(nzn > nz && nz%2 == 1)  zshift = 1.0f;
			nzn2 = nzn/2;
		}  else {
			nzn2 = 0;
		}
	} else {
		nyn2 = 0; nzn2 = 0;
	}
	if(xshift == 1.0 || yshift == 1.0 || zshift == 1.0)  {
		ret->set_array_offsets(1,1,1);
		int  lsdn2 = lsd/2;
		for (int iz = 1; iz <= nzn; iz++) {
			int jz=iz-1; if(jz>nzn2) jz=jz-nzn;
			for (int iy = 1; iy <= nyn; iy++) {
				int jy=iy-1; if(jy>nyn2) jy=jy-nyn;
				for (int ix = 1; ix <= lsdn2; ix++) {
					int jx=ix-1;
					ret->cmplx(ix,iy,iz) *=
					exp(-float(twopi)*iimag*(xshift*jx/nxn + yshift*jy/nyn+ zshift*jz/nzn));
				}
			}
		}
		ret->set_array_offsets(0,0,0);
	}*/
	ret->set_ri(1);
	ret->set_fftpad(true);
	ret->set_attr("npad", 1);
	if (nxn%2 == 1) {ret->set_fftodd(true);} else {ret->set_fftodd(false);}
	if(RetReal) {
		ret->do_ift_inplace();
		ret->depad();
	}
	ret->update();

	/*Dict d1 = temp_ft->get_attr_dict();
	Dict d2 = ret->get_attr_dict();
	printf("-----------------Attribute Dict for temp_ft--------------\n");
	EMUtil::dump_dict(d1);
	printf("-----------------Attribute Dict for ret--------------\n");
	EMUtil::dump_dict(d2);*/
	delete temp_ft;
	temp_ft = 0;
	return ret;
}

EMData *EMData::FourTruncate(int nxn, int nyn, int nzn, bool RetReal, bool normalize) {

	int lsd, lsdn, inx, iny, inz;
	int i, j, k;
	float  *fint;
	EMData *temp_ft = NULL;
	//if (is_complex())
	//	throw ImageFormatException("Input image has to be real");

	if (is_complex()) {
		if(nxn>nx - 2 + nx%2 || nyn>ny || nzn>nz)	throw ImageDimensionException("Cannot increase the image size");
		lsd = nx;
		lsdn = nxn + 2 - nxn%2;
		fint = get_data();
	} else {
		//  do out of place ft
		if(nxn>nx || nyn>ny || nzn>nz)	throw ImageDimensionException("Cannot increase the image size");
		lsd = nx + 2 - nx%2;
		lsdn = nxn + 2 - nxn%2;
		temp_ft = do_fft();
		fint = temp_ft->get_data();
	}
	EMData *ret = this->copy_head();
	ret->set_size(lsdn, nyn, nzn);
	float *fout = ret->get_data();
	inx = nx - nxn;  iny = ny - nyn;  inz = nz - nzn;
	for (k=1; k<=nzn/2+1; k++) for (j=1; j<=nyn/2+1; j++) for (i=1; i<=lsdn; i++) fout(i,j,k)=fint(i,j,k);
	if(nyn>1) {
		for (k=1; k<=nzn/2+1; k++) for (j=nyn/2+2; j<=nyn; j++) for (i=1; i<=lsdn; i++) fout(i,j,k)=fint(i,j+iny,k);
		if(nzn>1) {
			for (k=nzn/2+2; k<=nzn; k++) {
				for (j=1; j<=nyn/2+1; j++) {
					for (i=1; i<=lsdn; i++) {
						fout(i,j,k)=fint(i,j,k+inz);
					}
				}
				for (j=nyn/2+2; j<=nyn; j++) {
					for (i=1; i<=lsdn; i++) {
						fout(i,j,k)=fint(i,j+iny,k+inz);
					}
				}
			}
		}
	}
//       WEIGHTING FACTOR USED FOR EVEN NSAM. REQUIRED SINCE ADDING ZERO FOR
//       INTERPOLATION WILL INTRODUCE A COMPLEX CONJUGATE FOR NSAM/2'TH
//       ELEMENT.
	/*
	if(nxn%2 == 0 && inx !=0) {
		for (k=1; k<=nzn; k++) {
			for (j=1; j<=nyn; j++) {
				fout(nxn+1,j,k) *= sq2;
				fout(nxn+2,j,k) *= sq2;
			}
		}
		if(nyn>1) {
			for (k=1; k<=nzn; k++) {
			  for (i=1; i<=lsdn; i++) {
			    fout(i,nyn/2+1+iny,k) = sq2*fout(i,nyn/2+1,k);
			    fout(i,nyn/2+1,k) *= sq2;
			  }
			}
			if(nzn>1) {
				for (j=1; j<=nyn; j++) {
					for (i=1; i<=lsdn; i++) {
						fout(i,j,nzn/2+1+inz) = sq2*fout(i,j,nzn/2+1);
						fout(i,j,nzn/2+1) *= sq2;
					}
				}
			}
		}
	}*/
//  TO KEEP EXACT VALUES ON THE ORIGINAL GRID ONE SHOULD USE
//  SQ2     = 2.0. HOWEVER, TOTAL ENERGY WILL NOT BE CONSERVED
	//float  sq2 = std::sqrt(2.0f);
	if( normalize )  {
		float  anorm = (float) nxn* (float) nyn* (float) nzn/(float) nx/ (float) ny/ (float) nz;
		for (i = 0; i < lsdn*nyn*nzn; i++)  fout[i] *= anorm;
	}
	ret->set_complex(true);
	ret->set_ri(1);
	ret->set_fftpad(true);
	ret->set_attr("npad", 1);
	if (nxn%2 == 1) {ret->set_fftodd(true);} else {ret->set_fftodd(false);}
	if(RetReal) {
		ret->do_ift_inplace();
		ret->depad();
	}
	ret->update();

	/*Dict d1 = temp_ft->get_attr_dict();
	Dict d2 = ret->get_attr_dict();
	printf("-----------------Attribute Dict for temp_ft--------------\n");
	EMUtil::dump_dict(d1);
	printf("-----------------Attribute Dict for ret--------------\n");
	EMUtil::dump_dict(d2);*/
	if (!is_complex()) {
		delete temp_ft;
		temp_ft = 0;
	}
	return ret;
}
/*
EMData *EMData::FourInterpol_i(int nxn, int nyni, int nzni, bool RetReal) {

	int nyn, nzn, lsd, lsdn, inx, iny, inz;
	int i, j, k;

	if(ny > 1) {
		nyn = nyni;
		if(nz > 1) {
			nzn = nzni;
		}  else {
			nzn = 1;
		}
	} else {
		nyn = 1; nzn = 1;
	}
	if(nxn<nx || nyn<ny || nzn<nz)	throw ImageDimensionException("Cannot reduce the image size");
	lsd = nx-2 + 2 - nx%2;
	lsdn = nxn + 2 - nxn%2;
//  do out of place ft
	EMData *temp_ft = this->copy();
	EMData *ret = this->copy();
	ret->set_size(lsdn, nyn, nzn);
	ret->to_zero();
	float *fout = ret->get_data();
	float *fint = temp_ft->get_data();
//  TO KEEP EXACT VALUES ON THE ORIGINAL GRID ONE SHOULD USE
//  SQ2     = 2.0. HOWEVER, TOTAL ENERGY WILL NOT BE CONSERVED
	float  sq2 = 1.0f/std::sqrt(2.0f);
	float  anorm = (float) nxn* (float) nyn* (float) nzn/(float) nx/ (float) ny/ (float) nz;
	for (i = 0; i < lsd*ny*nz; i++)  fint[i] *= anorm;
	inx = nxn-(nx-2); iny = nyn - ny; inz = nzn - nz;
	for (k=1; k<=nz/2+1; k++) for (j=1; j<=ny/2+1; j++) for (i=1; i<=lsd; i++) fout(i,j,k)=fint(i,j,k);
	if(nyn>1) {
	//cout << "  " <<nxn<<"  " <<nyn<<" A " <<nzn<<endl;
		for (k=1; k<=nz/2+1; k++) for (j=ny/2+2+iny; j<=nyn; j++) for (i=1; i<=lsd; i++) fout(i,j,k)=fint(i,j-iny,k);
		if(nzn>1) {
			for (k=nz/2+2+inz; k<=nzn; k++) {
				for (j=1; j<=ny/2+1; j++) {
					for (i=1; i<=lsd; i++) {
						fout(i,j,k)=fint(i,j,k-inz);
					}
				}
				for (j=ny/2+2+iny; j<=nyn; j++) {
					for (i=1; i<=lsd; i++) {
						fout(i,j,k)=fint(i,j-iny,k-inz);
					}
				}
			}
		}
	}
//       WEIGHTING FACTOR USED FOR EVEN NSAM. REQUIRED SINCE ADDING ZERO FOR
//       INTERPOLATION WILL INTRODUCE A COMPLEX CONJUGATE FOR NSAM/2'TH
//       ELEMENT.
	if(nx%2 == 0 && inx !=0) {
		for (k=1; k<=nzn; k++) {
			for (j=1; j<=nyn; j++) {
				fout(nx-2+1,j,k) *= sq2;
				fout(nx-2+2,j,k) *= sq2;
			}
		}
		if(nyn>1) {
			for (k=1; k<=nzn; k++) {
			  for (i=1; i<=lsd; i++) {
			    fout(i,ny/2+1+iny,k) = sq2*fout(i,ny/2+1,k);
			    fout(i,ny/2+1,k) *= sq2;
			  }
			}
			if(nzn>1) {
				for (j=1; j<=nyn; j++) {
					for (i=1; i<=lsd; i++) {
						fout(i,j,nz/2+1+inz) = sq2*fout(i,j,nz/2+1);
						fout(i,j,nz/2+1) *= sq2;
					}
				}
			}
		}
	}
	ret->set_complex(true);
	ret->set_ri(1);
	ret->set_fftpad(true);
	ret->set_attr("npad", 1);
	if (nxn%2 == 1) {ret->set_fftodd(true);} else {ret->set_fftodd(false);}
	if(RetReal) {
		ret->do_ift_inplace();
		ret->depad();
	}
	ret->update();

	delete temp_ft;
	temp_ft = 0;
	return ret;
}
*/

EMData *EMData::Four_ds(int nxn, int nyni, int nzni, bool RetReal) {

	int nyn, nzn, lsd, lsdn, inx, iny, inz;
	int i, j;

	if(ny > 1) {
		nyn = nyni;
		if(nz > 1) {
			nzn = nzni;
		}  else {
			nzn = 1;
		}
	} else {
		nyn = 1; nzn = 1;
	}
	lsd = nx-2 + 2 - nx%2;
	lsdn = nxn + 2 - nxn%2;
//  do out of place ft
	EMData *temp_ft = this->copy();
	EMData *ret = this->copy();
	ret->set_size(lsdn, nyn, nzn);
	ret->to_zero();
	float *fout = ret->get_data();
	float *fint = temp_ft->get_data();
//  TO KEEP EXACT VALUES ON THE ORIGINAL GRID ONE SHOULD USE
//  SQ2     = 2.0. HOWEVER, TOTAL ENERGY WILL NOT BE CONSERVED
//	float  sq2 = 1.0f/std::sqrt(2.0f);
	float  anorm = (float) nxn* (float) nyn* (float) nzn/(float) nx/ (float) ny/ (float) nz;
	for (i = 0; i < lsd*ny*nz; i++)  fint[i] *= anorm;
	inx = nxn-(nx-2); iny = nyn - ny; inz = nzn - nz;
	for (j=1; j<=nyn; j++)
		for (i=1; i<=lsdn; i++)
			fout(i,j,1)=fint((i-1)/2*4+2-i%2,j*2-1,1);
	ret->set_complex(true);
	ret->set_ri(1);
	//ret->set_fftpad(true);
	//ret->set_attr("npad", 1);
	if (nxn%2 == 1) {ret->set_fftodd(true);} else {ret->set_fftodd(false);}
	if(RetReal) {
		ret->do_ift_inplace();
		ret->depad();
	}
	ret->update();

	delete temp_ft;
	temp_ft = 0;
	return ret;
}

EMData *EMData::Four_shuf_ds_cen_us(int nxn, int nyni, int, bool RetReal) {

	int nyn, nzn, lsd, lsdn, inx, iny, inz;
	int i, j;

	nyn = nyni;
	nzn = 1;
	lsd = nx;
	lsdn = nxn + 2 - nxn%2;

	EMData *temp_ft = this->copy();
	EMData *ret = this->copy();
	ret->set_size(lsdn, nyn, nzn);
	ret->to_zero();
	float *fout = ret->get_data();
	float *fint = temp_ft->get_data();
//  TO KEEP EXACT VALUES ON THE ORIGINAL GRID ONE SHOULD USE
//  SQ2     = 2.0. HOWEVER, TOTAL ENERGY WILL NOT BE CONSERVED
	float  sq2 = 1.0f/std::sqrt(2.0f);

	for (size_t i = 0; i < (size_t)lsd*ny*nz; i++)  fint[i] *= 4;

	inx = nxn-(nx-2); iny = nyn - ny; inz = nzn - nz;
	for (j=1; j<=ny/4; j++)
		for (i=1; i<=(nx-2)/2+2; i++) {
			int g = (i-1)/2+1;
			if ((g+j)%2 == 0) {
				fout(i,j,1)=fint(g*4-2-i%2,j*2-1+ny/2,1);
			} else {
				fout(i,j,1)=-fint(g*4-2-i%2,j*2-1+ny/2,1);
			}
		}

	for (j=ny/4+1; j<=ny/4+1; j++)
		for (i=1; i<=(nx-2)/2+2; i++) {
			int g = (i-1)/2+1;
			if ((g+j)%2 == 0) {
				fout(i,j,1)=fint(g*4-2-i%2,j*2-1-ny/2,1);
			} else {
				fout(i,j,1)=-fint(g*4-2-i%2,j*2-1-ny/2,1);
			}
		}

	for (j=ny/4+2; j<=ny/2; j++)
		for (i=1; i<=(nx-2)/2+2; i++) {
			int g = (i-1)/2+1;
			if ((g+j)%2 == 0) {
				fout(i,j+ny/2,1)=fint(g*4-2-i%2,j*2-1-ny/2,1);
			} else {
				fout(i,j+ny/2,1)=-fint(g*4-2-i%2,j*2-1-ny/2,1);
			}
		}

	if (nx%2 == 0) {
		for (j=1; j<=nyn; j++) {
			fout((nx-2)/2+1,j,1) *= sq2;
			fout((nx-2)/2+2,j,1) *= sq2;
		}
		for (i=1; i<=lsd/2+1; i++) {
			fout(i,ny/4+1+ny/2,1) = sq2*fout(i,ny/4+1,1);
			fout(i,ny/4+1,1) *= sq2;
		}
	}

	ret->set_complex(true);
	ret->set_ri(1);

	if (nxn%2 == 1) {ret->set_fftodd(true);} else {ret->set_fftodd(false);}
	if(RetReal) {
		ret->do_ift_inplace();
		ret->depad();
	}
	ret->update();

	delete temp_ft;
	temp_ft = 0;
	return ret;
}

#undef fint
#undef fout

#define  fint(jx,jy,jz)  fint[jx + (jy + jz*ny)*(size_t)nox]
EMData *EMData::filter_by_image(EMData* image, bool RetReal) {


	bool   complex_input = this->is_complex();
	nx  = this->get_xsize();
	ny  = this->get_ysize();
	nz  = this->get_zsize();
	int nox;
	if (complex_input) nox = (nx - 2 + this->is_fftodd()); else nox = nx;

	int lsd2 = (nox + 2 - nox%2) / 2; // Extended x-dimension of the complex image

	EMData* fp = NULL; // output image
	if(complex_input) {
		// fimage must remain pristine
		fp = this->copy();
	} else {
		fp = this->norm_pad( false, 1);
		fp->do_fft_inplace();
	}
	fp->set_array_offsets(1,1,1);
	int nx2 = nox/2;
	int ny2 = ny/2;
	int nz2 = nz/2;
	float *fint = image->get_data();
	for ( int iz = 1; iz <= nz; iz++) {
		int jz=nz2-iz+1; if(jz<0) jz += nz;
		for ( int iy = 1; iy <= ny; iy++) {
			int jy=ny2-iy+1; if(jy<0) jy += ny;
			for ( int ix = 1; ix <= lsd2; ix++) {
				int jx = nx2-ix+1;
				fp->cmplx(ix,iy,iz) *= fint(jx,jy,jz);
			}
		}
	}

	fp->set_ri(1);
	fp->set_fftpad(true);
	fp->set_attr("npad", 1);
	if (nx%2 == 1) fp->set_fftodd(true);
	else fp->set_fftodd(false);
	if(RetReal) {
		fp->do_ift_inplace();
		fp->depad();
	}
	fp->set_array_offsets(0,0,0);
	fp->update();

	return fp;
}
#undef   fint
#define  fint(jx,jy,jz)  fint[jx + (jy + jz*ny)*(size_t)nx]
#define  fout(jx,jy,jz)  fout[jx + (jy + jz*ny)*(size_t)nx]
EMData *EMData::replace_amplitudes(EMData* image, bool RetReal) {


	bool   complex_input = this->is_complex();
	nx  = this->get_xsize();
	ny  = this->get_ysize();
	nz  = this->get_zsize();
	int nox;
	if (complex_input) nox = (nx - 2 + this->is_fftodd()); else nox = nx;

	EMData* fp = NULL; // output image
	if(complex_input) {
		// fimage must remain pristine
		fp = this->copy();
	} else {
		fp = this->norm_pad( false, 1);
		fp->do_fft_inplace();
	}
	float *fout = fp->get_data();
	float *fint = image->get_data();
	for ( int iz = 0; iz < nz; iz++) {
		for ( int iy = 0; iy < ny; iy++) {
			for ( int ix = 0; ix < nx; ix+=2) {
				float qt = fint(ix,iy,iz)*fint(ix,iy,iz)+fint(ix+1,iy,iz)*fint(ix+1,iy,iz);
				float rt = fout(ix,iy,iz)*fout(ix,iy,iz)+fout(ix+1,iy,iz)*fout(ix+1,iy,iz);
				if(rt > 1.0e-20) {
						fout(ix,iy,iz) *= (qt/rt);
						fout(ix+1,iy,iz) *= (qt/rt);
				} else {
						qt = std::sqrt(qt/2.0f);
						fout(ix,iy,iz) = qt;
						fout(ix+1,iy,iz) = qt;
				}
			}
		}
	}

	fp->set_ri(1);
	fp->set_fftpad(true);
	fp->set_attr("npad", 1);
	if (nx%2 == 1) fp->set_fftodd(true);
	else fp->set_fftodd(false);
	if(RetReal) {
		fp->do_ift_inplace();
		fp->depad();
	}
	fp->set_array_offsets(0,0,0);
	fp->update();

	return fp;
}
#undef fint
#undef fout


#undef QUADPI
#undef deg_rad
