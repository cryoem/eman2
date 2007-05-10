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

#include "emdata.h"
#include <iostream>
#include "math.h"
 
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <vector>
#include <boost/shared_ptr.hpp>

using boost::shared_ptr;
using std::vector;
using std::cout;
using namespace EMAN;
using namespace std;

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
			std::complex <float> overallFactor = pow(tempF,m);  //(-i)^m ;  % I dropped off the 2 pi
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
		done_data();
		rhoOfkmB-> done_data();
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
	rhoOfkandm ->done_data();

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
	outCopy->done_data();
	outCopy->set_complex(true);
	if(outCopy->get_ysize()==1 && outCopy->get_zsize()==1) {
		outCopy->set_complex_x(true);
	}
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


EMData *EMData::FH2Real(int Size, float OverSamplekB, int IntensityFlag)  // PRB
{
	EMData* FFT= FH2F(Size,OverSamplekB,0);
	FFT->process_inplace("eman1.xform.fourierorigin");
	EMData* eguess= FFT ->do_ift();
	return eguess;
}  // ends FH2F

float dist(int lnlen, const float* line_1, const float* line_2)
{
    double dis2=0.0;
    for( int i=0; i < lnlen; ++i)
    {
       float tmp = line_1[i] - line_2[i];
       dis2 += tmp*tmp;
    }
    return std::sqrt( dis2 );
}   

float dist_r(int lnlen, const float* line_1, const float* line_2)
{
    double dis2 = 0.0;
    for( int i=0; i < lnlen; ++i )
    {
        float tmp = line_1[lnlen-1-i] - line_2[i];
        dis2 += tmp*tmp;
    }
    return std::sqrt(dis2);
}


float EMData::cm_euc(EMData* sinoj, int n1, int n2, float alpha1, float alpha2)
{
    int lnlen = get_xsize();
    int nline = get_ysize();

    assert( n1 >=0 && n1 < nline );
    assert( n2 >=0 && n2 < nline );
    assert( alpha1>=0.0 && alpha1 < 360.0 );
    assert( alpha2>=0.0 && alpha2 < 360.0 );

    float* line_1 = get_data() + n1*lnlen;
    float* line_2 = sinoj->get_data() + n2*lnlen;
    float just = (alpha1-180.0)*(alpha2-180.0);
    if( just > 0.0 )
    {
        return dist(lnlen, line_1, line_2);
    }

    if( just == 0.0 )
    {
        float dist_1 = dist(lnlen, line_1, line_2);
	float dist_2 = dist_r(lnlen, line_1, line_2);
	return std::min(dist_1, dist_2);
    }

    assert( (alpha1-180.0)*(alpha2-180.0) < 0.0 );
    return dist_r(lnlen, line_1, line_2);
}



EMData* EMData::rotavg() {

	int rmax;

	ENTERFUNC;

	if (ny<2 && nz <2) {
		LOGERR("No 1D images.");
		throw ImageDimensionException("No 1D images!");
	}
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(-nx/2,-ny/2,-nz/2);
#ifdef _WIN32
	//int rmax = _MIN(nx/2 + nx%2, ny/2 + ny%2);
	if ( nz == 1 ) {
		rmax = _MIN(nx/2 + nx%2, ny/2 + ny%2);
	} else {
		rmax = _MIN(nx/2 + nx%2, _MIN(ny/2 + ny%2, nz/2 + nz%2));
	}
#else
	//int rmax = std::min(nx/2 + nx%2, ny/2 + ny%2);
	if ( nz == 1 ) {
		rmax = std::min(nx/2 + nx%2, ny/2 + ny%2);	
	} else {
		rmax = std::min(nx/2 + nx%2, std::min(ny/2 + ny%2, nz/2 + nz%2));
	}
#endif	//_WIN32
	EMData* ret = new EMData();
	ret->set_size(rmax+1, 1, 1);
	ret->to_zero();
	vector<float> count(rmax+1);
	for (int k = -nz/2; k < nz/2 + nz%2; k++) {
	   if (abs(k) > rmax) continue;
	   for (int j = -ny/2; j < ny/2 + ny%2; j++) {
		if (abs(j) > rmax) continue;
		for (int i = -nx/2; i < nx/2 + nx%2; i++) {
			float r = std::sqrt(float(k*k) + float(j*j) + float(i*i));
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
	for (int ir = 0; ir <= rmax; ir++) {
	#ifdef _WIN32
		(*ret)(ir) /= _MAX(count[ir],1.0f);
	#else
		(*ret)(ir) /= std::max(count[ir],1.0f);
	#endif	//_WIN32
	}

	set_array_offsets(saved_offsets);
	ret->update();
	ret->done_data();
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
	float padded_value = 0.0, number_of_pixel = 0.0, r, frac;
	int i, j, k, ir;
	for ( k = -nz/2; k < nz/2 + nz%2; k++) {
		if (abs(k) > rmax) continue;
		for ( j = -ny/2; j < ny/2 + ny%2; j++) {
			if (abs(j) > rmax) continue;
			for (i = -nx/2; i < nx/2 + nx%2; i++) {
				r = std::sqrt(float(k*k) + float(j*j) + float(i*i));
				ir = int(r);
				if (ir > rmax || ir < rmax-2 ) continue ;				 
				else
      					{
	      					padded_value += (*avg1D)(ir) ;
	      					number_of_pixel += 1.0 ;
					}
			}
		}
	}
	padded_value /= number_of_pixel ;
	for ( k = -nz/2; k < nz/2 + nz%2; k++) 
	{
		for ( j = -ny/2; j < ny/2 + ny%2; j++) 
			{
					for ( i = -nx/2; i < nx/2 + nx%2; i++) 
					{
						r = std::sqrt(float(k*k) + float(j*j) + float(i*i));
						ir = int(r);
						if (ir >= rmax) (*result)(i,j,k) = padded_value ;
						else  
						{
							frac = r - float(ir); 
							(*result)(i,j,k) = (*avg1D)(ir)*(1.0f - frac)+(*avg1D)(ir+1)*frac;
						}	
						
					}
			}
	}				
	result->set_array_offsets(0,0,0);
	EXITFUNC;
	return result;
}

#define rdata(i,j,k) rdata[(i-1)+((j-1)+(k-1)*ny)*nx]
#define square(x) ((x)*(x))
vector<float> EMData::cog() {
	
	vector<float> cntog;
	int ndim = get_ndim();
	int i=1,j=1,k=1;
	float val,sum1=0.f,MX=0.f,RG=0.f,MY=0.f,MZ=0.f,r=0.f;
	
	if (ndim == 1)
	{
			for ( i = 1;i <= nx; i++)
			{
				val   = rdata(i,j,k);
				sum1 += val;
				MX   += ((i-1)*val);
			}
			MX=(MX/sum1);
			for ( i = 1;i <= nx; i++)
			{
				val   = rdata(i,j,k);
				sum1 += val;
				RG   += val*(square(MX - (i-1)));
			}
			RG=std::sqrt(RG/sum1);
			MX=MX-(nx/2);
			cntog.push_back(MX);
			cntog.push_back(RG);
#ifdef _WIN32
			cntog.push_back(Util::round(MX));
#else
			cntog.push_back(round(MX));
#endif	//_WIN32
	}	
	else if (ndim == 2)
	{	
			for (j=1;j<=ny;j++)
				{
					for (i=1;i<=nx;i++)
					{
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
			for (j=1;j<=ny;j++)
				{
					r = (square(MY-(j-1)));
					for (i=1;i<=nx;i++)
					{
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
			cntog.push_back(Util::round(MX));cntog.push_back(Util::round(MY));
#else
			cntog.push_back(round(MX));cntog.push_back(round(MY));
#endif	//_WIN32
	}
	else 
	{		
			for (k = 1;k <= nz;k++)
			{
				for (j=1;j<=ny;j++)
				{
					for (i=1;i<=nx;i++)
					{
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
			for (k = 1;k <= nz;k++)
			{
				for (j=1;j<=ny;j++)
				{
					float r = (square(MY-(j-1)) + square(MZ - (k-1)));
					for (i=1;i<=nx;i++)
					{
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
			cntog.push_back(Util::round(MX));cntog.push_back(Util::round(MY));cntog.push_back(Util::round(MZ));
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
	int needfree=0, nx, ny, nz, nx2, ny2, nz2, ix, iy, iz, kz, ky, ii;
	float  dx2, dy2, dz2, argx, argy, argz;

	if (!with) {
		throw NullPointerException("NULL input image");
	}

	
	EMData *f = this;
	EMData *g = with;
	
	nx  = f->get_xsize();
	ny  = f->get_ysize();
	nz  = f->get_zsize();

	if (ny==0 && nz==0) {
		throw ImageFormatException( "Cannot calculate FSC for 1D images");
	}

	if (!equalsize(f, g)) {
		LOGERR("FSC requires congruent images");
		throw ImageDimensionException("FSC requires congruent images");
	}

	if (f->is_complex()) nx = (nx - 2 + f->is_fftodd()); // nx is the real-space size of the input image
	int lsd2 = (nx + 2 - nx%2) ; // Extended x-dimension of the complex image

//  Process f if real
	EMData* fpimage = NULL;
	if(f->is_complex()) fpimage = f;
	else {fpimage= norm_pad_ft(f, false, false); needfree|=1;} // Extend and do the FFT if f is real


//  Process g if real
	EMData* gpimage = NULL;
	if(g->is_complex()) gpimage = g; 
	else {gpimage= norm_pad_ft(g, false, false); needfree|=2;} // Extend and do the FFT if f is real


	float *d1 = fpimage->get_data();
	float *d2 = gpimage->get_data();

	nx2=nx/2; ny2 = ny/2; nz2 = nz/2;
	dx2 = 1.0f/float(nx2)/float(nx2); 
	dy2 = 1.0f/float(ny2)/float(ny2);

#ifdef _WIN32
	dz2 = 1.0f / _MAX(float(nz2),1.0f)/_MAX(float(nz2),1.0f);
	
	int inc = Util::round(float( _MAX( _MAX(nx2,ny2),nz2) )/w );
#else
	dz2 = 1.0f/std::max(float(nz2),1.0f)/std::max(float(nz2),1.0f);
	
	int inc = Util::round(float(std::max(std::max(nx2,ny2),nz2))/w);
#endif	//_WIN32

	double *ret = new double[inc+1];
	double *n1  = new double[inc+1];
	double *n2  = new double[inc+1];
	float  *lr  = new float[inc+1];
	for (int i = 0; i <= inc; i++) {
		ret[i] = 0; n1[i] = 0; n2[i] = 0; lr[i]=0;
	}

	for ( iz = 0; iz <= nz-1; iz++) {
		if(iz>nz2) kz=iz-nz; else kz=iz; argz = float(kz*kz)*dz2;
		for ( iy = 0; iy <= ny-1; iy++) {
			if(iy>ny2) ky=iy-ny; else ky=iy; argy = argz + float(ky*ky)*dy2;
			for ( ix = 0; ix <= lsd2-1; ix+=2) {
			// Skip Friedel related values
			   if(ix>0 || (kz>=0 && (ky>=0 || kz!=0))) {
				argx = 0.5f*std::sqrt(argy + float(ix*ix)*0.25f*dx2);
				int r = Util::round(inc*2*argx);
				if(r <= inc) {
					ii = ix + (iy  + iz * ny)* lsd2;
					ret[r] += d1[ii] * double(d2[ii]) + d1[ii + 1] * double(d2[ii + 1]);
					n1[r]  += d1[ii] * double(d1[ii]) + d1[ii + 1] * double(d1[ii + 1]);
					n2[r]  += d2[ii] * double(d2[ii]) + d2[ii + 1] * double(d2[ii + 1]);
					lr[r]  +=2;
				}
			   }
			}
		}
	}


	int  linc = 0;
	for (int i = 0; i <= inc; i++) if(lr[i]>0) linc++;
	

	vector < float >result(linc*3);

	ii = -1;
	for (int i = 0; i <= inc; i++) {
		if(lr[i]>0) {
			ii++;
			result[ii]        = float(i)/float(2*inc);
			result[ii+linc]   = float(ret[i] / (std::sqrt(n1[i] * n2[i])));
			result[ii+2*linc] = lr[i]  /*1.0f/sqrt(float(lr[i]))*/;}
		/*else {
			result[i]           = 0.0f;
			result[i+inc+1]     = 0.0f;
			result[i+2*(inc+1)] = 0.0f;}*/
	}

	if( ret )
	{
		delete[]ret;
		ret = 0;
	}

	if( n1 )
	{
		delete[]n1;
		n1 = 0;
	}
	if( n2 )
	{
		delete[]n2;
		n2 = 0;
	}

	if (needfree&1)
	{
		if( fpimage )
		{
			delete fpimage;
			fpimage = 0;
		}
	}
	if (needfree&2)
	{
		if( gpimage )
		{
			delete gpimage;
			gpimage = 0;
		}
	}

	EXITFUNC;
	return result;
}



EMData* EMData::symvol(string symString) {
	ENTERFUNC;
	int nsym = Transform3D::get_nsym(symString); // number of symmetries
	Transform3D sym;
	// set up output volume
	EMData *svol = new EMData;
	svol->set_size(nx, ny, nz);
	svol->to_zero();
	// set up new copy
	EMData* symcopy = new EMData;
	symcopy->set_size(nx, ny, nz);
	// set up coord grid
	// actual work -- loop over symmetries and symmetrize
	for (int isym = 0; isym < nsym; isym++) {
		Transform3D rm = sym.get_sym(symString, isym);
		symcopy = this -> rot_scale_trans(rm);
		*svol += (*symcopy);
	}
	*svol /=  ((float) nsym);
	svol->done_data();
	svol->update();
	EXITFUNC;
	return svol;
}

#define proj(ix,iy,iz)  proj[ix-1+(iy-1+(iz-1)*ny)*nx]
#define pnewimg(ix,iy,iz)  pnewimg[ix-1+(iy-1+(iz-1)*ny)*nx]
EMData* EMData::average_circ_sub() {
//  this is written as though dimensions could be different, but in fact they should be all equal nx=ny=nz,
//                                                           no check of this	
	ENTERFUNC;
	EMData* image = this;
	EMData* newimg = copy_head();
	newimg->set_size(nx,ny,nz);
	float *proj = image->get_data();
	float *pnewimg = newimg->get_data();
	//  Calculate average outside of a circle
	float r2 = (nx/2)*(nx/2);
	float qs=0.0f;
	int m=0;
	int ncz = nz/2 + 1;
	int ncy = ny/2 + 1;
	int ncx = nx/2 + 1;
	for (int iz = 1; iz <= nz; iz++) { float yy = (iz-ncz)*(iz-ncz);
		for (int iy = 1; iy <=ny; iy++) { float xx = yy + (iy-ncy)*(iy-ncy);
			for (int ix = 1; ix <= nx; ix++) {
				if ( xx+float((ix-ncx)*(ix-ncx)) > r2 ) {
					qs += proj(ix,iy,iz);
					m++;
				}
			}
		}
	}
	qs /= m;
	for (int iz = 1; iz <= nz; iz++) 
		for (int iy = 1; iy <= ny; iy++) 
			for (int ix = 1; ix <= nx; ix++)
					pnewimg(ix,iy,iz) = proj(ix,iy,iz) - qs;
	newimg->done_data();
	return newimg;
	EXITFUNC;
}


//  Helper functions for method nn


void EMData::onelinenn(int j, int n, int n2, 
		       EMData* wptr, EMData* bi, const Transform3D& tf)
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
			if ((ixn <= n2) && (iyn >= -n2) && (iyn <= n2)
				            && (izn >= -n2) && (izn <= n2)) {
				if (ixn >= 0) {
					int iza, iya;
					if (izn >= 0) {
						iza = izn + 1;
					} else {
						iza = n + izn + 1;
					}
					if (iyn >= 0) {
						iya = iyn + 1;
					} else {
						iya = n + iyn + 1;
					}
					cmplx(ixn,iya,iza) += btq;
					//std::cout<<"    "<<j<<"  "<<ixn<<"  "<<iya<<"  "<<iza<<"  "<<btq<<std::endl;
					(*wptr)(ixn,iya,iza)++;
				} else {
					int izt, iyt;
					if (izn > 0) {
						izt = n - izn + 1;
					} else {
						izt = -izn + 1;
					}
					if (iyn > 0) {
						iyt = n - iyn + 1;
					} else {
						iyt = -iyn + 1;
					}
					cmplx(-ixn,iyt,izt) += conj(btq);
					//std::cout<<" *  "<<j<<"  "<<ixn<<"  "<<iyt<<"  "<<izt<<"  "<<btq<<std::endl;
					(*wptr)(-ixn,iyt,izt)++;
				}
			}
		}
	}
}


void EMData::onelinenn_mult(int j, int n, int n2, 
		       EMData* wptr, EMData* bi, const Transform3D& tf, int mult)
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
			if ((ixn <= n2) && (iyn >= -n2) && (iyn <= n2)
				            && (izn >= -n2) && (izn <= n2)) {
				if (ixn >= 0) {
					int iza, iya;
					if (izn >= 0) {
						iza = izn + 1;
					} else {
						iza = n + izn + 1;
					}
					if (iyn >= 0) {
						iya = iyn + 1;
					} else {
						iya = n + iyn + 1;
					}
					cmplx(ixn,iya,iza) += btq*float(mult);
					(*wptr)(ixn,iya,iza)+=float(mult);
				} else {
					int izt, iyt;
					if (izn > 0) {
						izt = n - izn + 1;
					} else {
						izt = -izn + 1;
					}
					if (iyn > 0) {
						iyt = n - iyn + 1;
					} else {
						iyt = -iyn + 1;
					}
					cmplx(-ixn,iyt,izt) += conj(btq)*float(mult);
					(*wptr)(-ixn,iyt,izt)+=float(mult);
				}
			}
		}
	}
}

void EMData::nn(EMData* wptr, EMData* myfft, const Transform3D& tf, int mult) 
{
	ENTERFUNC;
	int nxc = attr_dict["nxc"]; // # of complex elements along x
	// let's treat nr, bi, and local data as matrices
	vector<int> saved_offsets = get_array_offsets();
	vector<int> myfft_saved_offsets = myfft->get_array_offsets();
	set_array_offsets(0,1,1);
	myfft->set_array_offsets(0,1);
	// loop over frequencies in y
	if( mult == 1 )
	{
	    for (int iy = -ny/2 + 1; iy <= ny/2; iy++) onelinenn(iy, ny, nxc, wptr, myfft, tf);
	}
	else
	{
	    for (int iy = -ny/2 + 1; iy <= ny/2; iy++) onelinenn_mult(iy, ny, nxc, wptr, myfft, tf, mult);
        }

        set_array_offsets(saved_offsets);
	myfft->set_array_offsets(myfft_saved_offsets);
	EXITFUNC;
}

void EMData::nn_SSNR(EMData* wptr, EMData* wptr2, EMData* myfft, const Transform3D& tf, int mult)
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
						if (izn >= 0) {
							iza = izn + 1;
						} else {
							iza = nz + izn + 1;
						}
						if (iyn >= 0) {
							iya = iyn + 1;
						} else {
							iya = ny + iyn + 1;
						}
						cmplx(ixn,iya,iza) += btq;
						(*wptr)(ixn,iya,iza)++;
						(*wptr2)(ixn,iya,iza) += norm(btq);
					} else {
						int izt, iyt;
						if (izn > 0) {
							izt = nz - izn + 1;
						} else {
							izt = -izn + 1;
						}
						if (iyn > 0) {
							iyt = ny - iyn + 1;
						} else {
							iyt = -iyn + 1;
						}
						cmplx(-ixn,iyt,izt) += conj(btq);
						(*wptr)(-ixn,iyt,izt)++;
						(*wptr2)(-ixn,iyt,izt) += norm(btq);
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
	// let's treat the local data as a matrix
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


class ctf_store
{
public:

    static void init( int winsize, float voltage, float pixel, float Cs, float amp_contrast, float b_factor )
    {
        m_inited = true;
        m_winsize = winsize;
	m_voltage = voltage;
	m_pixel   = pixel;
	m_Cs      = Cs * 1.0e7f;
	m_amp_contrast = amp_contrast;
	m_b_factor = b_factor;

        m_winsize2= m_winsize*m_winsize;
        m_vecsize = m_winsize2/4;
        m_wgh=atan( amp_contrast/(1.0-amp_contrast) );
	m_lambda = 12.398f/std::sqrt(voltage*(1022.f+voltage));
    }

    static bool inited( )
    {
        return m_inited;
    }

    static float get_ctf( float defocus, int r2 )
    {
    /*
        assert( m_inited );

        shared_ptr< vector<float> > ptr;

        map< float, shared_ptr< vector<float> > >::iterator i = m_store.find( defocus );
 
        if( i == m_store.end() )
	{
	    ptr = shared_ptr< vector<float> >( new vector<float>(m_vecsize, 0.0) );
	    m_store[defocus] = ptr;
	}    
        else
	{
	    ptr = i->second;
	}

        float ctf = ptr->at(r2);

	if( ctf == 0.0 )
	{
     */
	    float ak = std::sqrt( r2/float(m_winsize2) )/m_pixel;
	    float a = m_lambda*ak*ak;
	    float b = m_lambda*a*a;
	    float ctf = -sin(-M_PI*(defocus*a-m_Cs*b/2.)-m_wgh);
     //	    ptr->at(r2) = ctf;
     //    }

        return ctf;
    }

private:
 
    static bool m_inited;

    static int m_winsize, m_winsize2, m_vecsize;

    static float m_Cs, m_voltage, m_pixel, m_amp_contrast, m_b_factor;

    static float m_lambda, m_wgh;

    static map< float, shared_ptr< vector<float> > > m_store;
};

bool ctf_store::m_inited = false;

int ctf_store::m_winsize, ctf_store::m_winsize2, ctf_store::m_vecsize;

float ctf_store::m_Cs, ctf_store::m_voltage, ctf_store::m_pixel, ctf_store::m_amp_contrast, ctf_store::m_b_factor;

float ctf_store::m_lambda, ctf_store::m_wgh;

std::map< float, shared_ptr< std::vector<float> > > ctf_store::m_store;


//  Helper functions for method nn4_ctf
void EMData::onelinenn_ctf(int j, int n, int n2, 
		          EMData* w, EMData* bi, const Transform3D& tf, float defocus, int mult) {//std::cout<<"   onelinenn_ctf  "<<j<<"  "<<n<<"  "<<n2<<"  "<<std::endl;

	int jp = (j >= 0) ? j+1 : n+j+1;
	// loop over x
	for (int i = 0; i <= n2; i++) {
	        int r2 = i*i+j*j;
		if ( (r2<n*n/4) && !( (0==i) && (j<0) ) ) {
		        float  ctf = ctf_store::get_ctf( defocus, r2 );
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
					if (izn >= 0) {
						iza = izn + 1;
					} else {
						iza = n + izn + 1;
					}
					if (iyn >= 0) {
						iya = iyn + 1;
					} else {
						iya = n + iyn + 1;
					}
					cmplx(ixn,iya,iza) += btq*ctf*float(mult);
				       //	std::cout<<"    "<<j<<"  "<<ixn<<"  "<<iya<<"  "<<iza<<"  "<<ctf<<std::endl;
					(*w)(ixn,iya,iza) += ctf*ctf*mult;
				} else {
					int izt, iyt;
					if (izn > 0) {
						izt = n - izn + 1;
					} else {
						izt = -izn + 1;
					}
					if (iyn > 0) {
						iyt = n - iyn + 1;
					} else {
						iyt = -iyn + 1;
					}
					cmplx(-ixn,iyt,izt) += conj(btq)*ctf*float(mult);
				        //	std::cout<<" *  " << j << "  " <<-ixn << "  " << iyt << "  " << izt << "  " << ctf <<std::endl;
					(*w)(-ixn,iyt,izt) += ctf*ctf*float(mult);
				}
			}
		}
	}
}

void EMData::onelinenn_ctf_applied(int j, int n, int n2, 
		          EMData* w, EMData* bi, const Transform3D& tf, float defocus, int mult) {//std::cout<<"   onelinenn_ctf  "<<j<<"  "<<n<<"  "<<n2<<"  "<<std::endl;

	int jp = (j >= 0) ? j+1 : n+j+1;
	// loop over x
	for (int i = 0; i <= n2; i++) {
	        int r2 = i*i + j*j;
		if ( (r2< n*n/4) && !((0==i) && (j< 0)) ) {
                        float  ctf = ctf_store::get_ctf( defocus, r2);

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

			if ((ixn <= n2) && (iyn >= -n2) && (iyn <= n2) && (izn >= -n2) && (izn <= n2)) {
				if (ixn >= 0) {
					int iza, iya;
					if (izn >= 0) {
						iza = izn + 1;
					} else {
						iza = n + izn + 1;
					}
					if (iyn >= 0) {
						iya = iyn + 1;
					} else {
						iya = n + iyn + 1;
					}
					cmplx(ixn,iya,iza) += btq*float(mult);
					(*w)(ixn,iya,iza) += mult*ctf*ctf;
				} else {
					int izt, iyt;
					if (izn > 0) {
						izt = n - izn + 1;
					} else {
						izt = -izn + 1;
					}
					if (iyn > 0) {
						iyt = n - iyn + 1;
					} else {
						iyt = -iyn + 1;
					}
					cmplx(-ixn,iyt,izt) += conj(btq)*float(mult);
					(*w)(-ixn,iyt,izt) += mult*ctf*ctf;
					//std::cout<<" *  "<<j<<"  "<<ixn<<"  "<<iyt<<"  "<<izt<<"  "<<btq<<std::endl;
				}
			}
		}
	}
}

void
EMData::nn_ctf(EMData* w, EMData* myfft, const Transform3D& tf, int mult) {
	ENTERFUNC;
	int nxc = attr_dict["nxc"]; // # of complex elements along x
	// let's treat nr, bi, and local data as matrices
	vector<int> saved_offsets = get_array_offsets();
	vector<int> myfft_saved_offsets = myfft->get_array_offsets();
	set_array_offsets(0,1,1);
	myfft->set_array_offsets(0,1);

	// if( ! ctf_store::inited() )
	{
            float Cs = myfft->get_attr( "Cs" );
            float pixel = myfft->get_attr( "Pixel_size" );
            float voltage = myfft->get_attr("voltage");
            float amp_contrast = myfft->get_attr( "amp_contrast" );
            float b_factor = 0.0;
            ctf_store::init( ny, voltage, pixel, Cs, amp_contrast, b_factor );
	}

	// loop over frequencies in y
        float defocus = myfft->get_attr( "defocus" );
	for (int iy = -ny/2 + 1; iy <= ny/2; iy++) onelinenn_ctf(iy, ny, nxc, w, myfft, tf, defocus, mult);
	set_array_offsets(saved_offsets);
	myfft->set_array_offsets(myfft_saved_offsets);
	EXITFUNC;
}

void
EMData::nn_ctf_applied(EMData* w, EMData* myfft, const Transform3D& tf, int mult) {
	ENTERFUNC;
	int nxc = attr_dict["nxc"]; // # of complex elements along x
	// let's treat nr, bi, and local data as matrices
	vector<int> saved_offsets = get_array_offsets();
	vector<int> myfft_saved_offsets = myfft->get_array_offsets();
	set_array_offsets(0,1,1);
	myfft->set_array_offsets(0,1);

	// if( ! ctf_store::inited() )
	{
            float Cs= myfft->get_attr( "Cs" );
            float pixel = myfft->get_attr( "Pixel_size" );
            float voltage = myfft->get_attr("voltage");
            float amp_contrast = myfft->get_attr( "amp_contrast" );
            float b_factor=0.0;
            ctf_store::init( ny, voltage, pixel, Cs, amp_contrast, b_factor );
	}

	// loop over frequencies in y
	float defocus = myfft->get_attr( "defocus" );
	for (int iy = -ny/2 + 1; iy <= ny/2; iy++) onelinenn_ctf_applied(iy, ny, nxc, w, myfft, tf, defocus, mult);
	set_array_offsets(saved_offsets);
	myfft->set_array_offsets(myfft_saved_offsets);
	EXITFUNC;
}



void EMData::nn_SSNR_ctf(EMData* wptr, EMData* wptr2, EMData* wptr3, EMData* wptr4, EMData* wptr5, EMData* myfft, EMData* m_wvolume, const Transform3D& tf, int mult)
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
	std::complex<float> tmpq, tmp2;
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
				if (xnew < 0.0) 
				{
					xnew = -xnew; // ensures xnew>=0.0
					ynew = -ynew;
					znew = -znew;
					btq = conj(myfft->cmplx(ix,jp));
				} 
				else 
				{
					btq = myfft->cmplx(ix,jp);
				}
				int ixn = int(xnew + 0.5 + nx) - nx; // ensures ixn >= 0
				int iyn = int(ynew + 0.5 + ny) - ny;
				int izn = int(znew + 0.5 + nz) - nz;
				if ((ixn <= nxc) && (iyn >= iymin) && (iyn <= iymax) && (izn >= izmin) && (izn <= izmax)) 
				{
					if (ixn >= 0) {
						int iza, iya;
						if (izn >= 0) {
							iza = izn + 1;
						} else {
							iza = nz + izn + 1;
						}
						if (iyn >= 0) {
							iya = iyn + 1;
						} else {
							iya = ny + iyn + 1;
						}
						tmpq = (*m_wvolume)(ixn,iya,iza);
						cmplx(ixn,iya,iza)    += btq*ctf;
						(*wptr)(ixn,iya,iza)  += ctf*ctf;						
						(*wptr5)(ixn,iya,iza) += ctf*ctf*std::norm(tmpq);
						(*wptr2)(ixn,iya,iza) += std::norm(btq);
						tmp2 = tmpq*ctf*conj(btq);
						(*wptr4)(ixn,iya,iza) += -2.*std::real(tmp2);
					} else {
						int izt, iyt;
						if (izn > 0) {
							izt = nz - izn + 1;
						} else {
							izt = -izn + 1;
						}
						if (iyn > 0) {
							iyt = ny - iyn + 1;
						} else {
							iyt = -iyn + 1;
						}
						tmpq = (*m_wvolume)(-ixn,iyt,izt);
						cmplx(-ixn,iyt,izt)     += std::conj(btq)*ctf;
						(*wptr) (-ixn,iyt,izt)  += ctf*ctf;
						(*wptr5)(-ixn,iyt,izt)  += ctf*ctf*std::norm(tmpq);
						(*wptr2)(-ixn,iyt,izt)  += std::norm(btq);
						tmp2 = tmpq*ctf*conj(btq);
						(*wptr4)(-ixn,iyt,izt)  += -2.*std::real(tmp2);
					}
				}
			}
		}
	}
	set_array_offsets(saved_offsets);
	myfft->set_array_offsets(myfft_saved_offsets);
	EXITFUNC;
}

void EMData::nn_wiener(EMData* wptr, EMData* wptr3, EMData* myfft, const Transform3D& tf, int mult)
{
     /*** Wiener volume calculating routine
          Counting Kn  
                                        ***/
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
				if (xnew < 0.0) 
				{
					xnew = -xnew; // ensures xnew>=0.0
					ynew = -ynew;
					znew = -znew;
					btq = conj(myfft->cmplx(ix,jp));
				} else 
				{
					btq = myfft->cmplx(ix,jp);
				}
				int ixn = int(xnew + 0.5 + nx) - nx; // ensures ixn >= 0
				int iyn = int(ynew + 0.5 + ny) - ny;
				int izn = int(znew + 0.5 + nz) - nz;
				if ((ixn <= nxc) && (iyn >= iymin) && (iyn <= iymax) && (izn >= izmin) && (izn <= izmax)) {
					if (ixn >= 0) 
					{
						int iza, iya;
						if (izn >= 0) 
						{
							iza = izn + 1;
						} else 
						{
							iza = nz + izn + 1;
						}
						if (iyn >= 0) 
						{
							iya = iyn + 1;
						} else 
						{
							iya = ny + iyn + 1;
						}
						cmplx(ixn,iya,iza)    += btq*ctf;						
						(*wptr)(ixn,iya,iza)  += ctf*ctf;
						(*wptr3)(ixn,iya,iza) += 1.0;
					} 
					else 
					{
						int izt, iyt;
						if (izn > 0) 
						{
							izt = nz - izn + 1;
						} else 
						{
							izt = -izn + 1;
						}
						if (iyn > 0) 
						{
							iyt = ny - iyn + 1;
						} else 
						{
							iyt = -iyn + 1;
						}
						cmplx(-ixn,iyt,izt)    += conj(btq)*ctf;
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
}

void EMData::symplane0_ctf(EMData* w) {
	ENTERFUNC;
	int nxc = attr_dict["nxc"];
	int n = nxc*2;
	// let's treat the local data as a matrix
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
	EXITFUNC;
}


EMData*
EMData::rot_scale_trans2D(float angDeg, float delx,float dely, float scale) { // quadratic, no background, 2D
	float ang=angDeg*M_PI/180.0f;
	if (1 >= ny)
		throw ImageDimensionException("Can't rotate 1D image");
	if (nz==1) { 
		vector<int> saved_offsets = get_array_offsets();
		set_array_offsets(0,0,0);
		if (0.f == scale) scale = 1.f; // silently fix common user error
		EMData* ret = copy_head();
		if(delx >= 0.0f) { delx = fmod(delx, float(nx));} else {delx = -fmod(-delx, float(nx));}
		if(dely >= 0.0f) { dely = fmod(dely, float(ny));} else {dely = -fmod(-dely, float(ny));}
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

					if (xold < 0.0f) xold = fmod(float(nx) - fmod(-xold, float(nx)), float(nx));
					else if (xold > (float) (nx-1) ) xold = fmod(xold, float(nx));
					if (yold < 0.0f) yold = fmod(float(ny) - fmod(-yold, float(ny)), float(ny));
					else if (yold > (float) (ny-1) ) yold = fmod(yold, float(ny));

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


EMData*
EMData::rot_scale_trans(const Transform3D &RA) {
	
	EMData* ret = copy_head();
	float *in = this->get_data();
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
	Vec3f  translations = RA.get_posttrans();
	Transform3D RAinv; // = new Transform3D();
	RAinv= RA.inverse();

	if (1 >= ny)
		throw ImageDimensionException("Can't rotate 1D image");
	if (nz==1) { 
	float  p1, p2, p3, p4;
	float delx = translations.at(0);
	float dely = translations.at(1);
	if(delx >= 0.0f) { delx = fmod(delx, float(nx));} else {delx = -fmod(-delx, float(nx));}
	if(dely >= 0.0f) { dely = fmod(dely, float(ny));} else {dely = -fmod(-dely, float(ny));}
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

				if (xold < 0.0f) xold = fmod(float(nx) - fmod(-xold, float(nx)), float(nx));
				else if (xold > (float) (nx-1) ) xold = fmod(xold, float(nx));
				if (yold < 0.0f) yold = fmod(float(ny) - fmod(-yold, float(ny)), float(ny));
				else if (yold > (float) (ny-1) ) yold = fmod(yold, float(ny));

				int xfloor = int(xold); int yfloor = int(yold);
				float t=xold-xfloor; float u = yold-yfloor;
				if(xfloor == nx -1 && yfloor == ny -1) {

				p1 =in[xfloor   + yfloor*ny];
				p2 =in[ yfloor*ny];
				p3 =in[0];
				p4 =in[xfloor];}

				else if(xfloor == nx - 1) {
				
				p1 =in[xfloor   + yfloor*ny];
				p2 =in[           yfloor*ny];
				p3 =in[          (yfloor+1)*ny];
				p4 =in[xfloor   + (yfloor+1)*ny];}

				else if(yfloor == ny - 1) {
				
				p1 =in[xfloor   + yfloor*ny];
				p2 =in[xfloor+1 + yfloor*ny];
				p3 =in[xfloor+1 ];
				p4 =in[xfloor   ];}
				
				else {
				p1 =in[xfloor   + yfloor*ny];
				p2 =in[xfloor+1 + yfloor*ny];
				p3 =in[xfloor+1 + (yfloor+1)*ny];
				p4 =in[xfloor   + (yfloor+1)*ny];}
				(*ret)(ix,iy) = p1 + u * ( p4 - p1) + t * ( p2 - p1 + u *(p3-p2-p4+p1));
			} //ends x loop
		} // ends y loop
		set_array_offsets(saved_offsets);
		return ret;
	} else {
//		 This begins the 3D version

	float delx = translations.at(0);
	float dely = translations.at(1);
	float delz = translations.at(2);
	if(delx >= 0.0f) { delx = fmod(delx, float(nx));} else {delx = -fmod(-delx, float(nx));}
	if(dely >= 0.0f) { dely = fmod(dely, float(ny));} else {dely = -fmod(-dely, float(ny));}
	if(dely >= 0.0f) { delz = fmod(delz, float(nz));} else {delz = -fmod(-delz, float(nz));}
	int xc = nx/2;
	int yc = ny/2;
	int zc = nz/2;
//         shifted center for rotation
	float shiftxc = xc + delx;
	float shiftyc = yc + dely;
	float shiftzc = zc + delz;
//                  set up array to use later
//
		int xArr[27];
		int yArr[27];
		int zArr[27];
		float fdata[27];
		
		for (int iL=0; iL<27 ; iL++){  // need this indexing array later
			xArr[iL]  =  (int) (fmod(iL,3.0f) - 1);
			yArr[iL]  =  (int)( fmod( ((int) (iL/3) ),3.0f)- 1);
			zArr[iL]  = ((int) (iL/9)  ) -1;
//			printf("iL=%d, \t xArr=%d, \t yArr=%d, \t zArr=%d \n",iL, xArr[iL],yArr[iL],zArr[iL]);
		}
					
//		for (int iz = 0; iz < nz; iz++) {for (int iy = 0; iy < ny; iy++) {for (int ix = 0; ix < nx; ix++) {
//		      (*ret)(ix,iy,iz) = 0;}}}   // initialize returned data
		
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


				if (xold < 0.0f) xold = fmod((int(xold/float(nx))+1)*nx-xold, float(nx));
				else if (xold > (float) (nx-1) ) xold = fmod(xold, float(nx));
				if (yold < 0.0f) yold =fmod((int(yold/float(ny))+1)*ny-yold, float(ny));
				else if (yold > (float) (ny-1) ) yold = fmod(yold, float(ny));
				if (zold < 0.0f) zold =fmod((int(zold/float(nz))+1)*nz-zold, float(nz));
				else if (zold > (float) (nz-1) ) zold = fmod(zold, float(nz));

//         This is currently coded the way  SPIDER coded it,
//            changing floor to round  in the next 3 lines below may be better
//					int IOX = (int) floor(xold); // This is the center of the array
//					int IOY = (int) floor(yold ); // In the next loop we interpolate
//					int IOZ = (int) floor(zold ); //  If floor is used dx is positive
					int IOX = int(xold);
					int IOY = int(yold);
					int IOZ = int(zold);

					float dx = xold-IOX; //remainder(xold,1);  //  now |dx| <= .5 
					float dy = yold-IOY; //remainder(yold,1);
					float dz = zold-IOZ; //remainder(zold,1);
					
//					printf(" IOX=%d \t IOY=%d \t IOZ=%d \n", IOX, IOY, IOZ);
//					if (IOX>=0 && IOX<nx  && IOY>=0 && IOY < ny && IOZ >= 0 && IOZ < nz ) {
//                                      	ROTATED POSITION IS INSIDE OF VOLUME
//						FIND INTENSITIES ON 3x3x3 COORDINATE GRID;
//                                     Solution is wrapped
						for  (int iL=0; iL<27 ; iL++){
							int xCoor = (int) fmod(IOX+xArr[iL] + nx + .0001f, (float) nx);
							int yCoor = (int) fmod(IOY+yArr[iL] + ny + .0001f, (float) ny);
							int zCoor = (int) fmod(IOZ+zArr[iL] + nz + .0001f, (float) nz);
							fdata[iL] = (*this)( xCoor, yCoor ,zCoor );
//							if (iy==iz && iz==0){printf(" fdata=%f \n", fdata[iL]);}
//						}
					}

					(*ret)(ix,iy,iz) = Util::triquad(dx, dy, dz, fdata);
//					(*ret)(ix,iy,iz) = Util:: trilinear_interpolate(fdata[13],fdata[14],fdata[16],
//											fdata[17],fdata[22],fdata[23],
//											fdata[25],fdata[26],dx, dy, dz);
//	p1 iL=13,   xArr= 0,         yArr= 0,         zArr= 0
//	p2 iL=14,   xArr= 1,         yArr= 0,         zArr= 0
//	p3 iL=16,   xArr= 0,         yArr= 1,         zArr= 0
//	p4 iL=17,   xArr= 1,         yArr= 1,         zArr= 0
//	p5 iL=22,   xArr= 0,         yArr= 0,         zArr= 1
//	p6 iL=23,   xArr= 1,         yArr= 0,         zArr= 1
//	p7 iL=25,   xArr= 0,         yArr= 1,         zArr= 1
//	p8 iL=26,   xArr= 1,         yArr= 1,         zArr= 1



				} //ends x loop
			} // ends y loop
		} // ends z loop

		set_array_offsets(saved_offsets);
		return ret;
/*		static inline float trilinear_interpolate(float p1, 
                     float p2, float p3,float p4, float p5, float p6,float p7, float p8, float t, float u, float v)
*/
//		throw ImageDimensionException("Volume not currently supported");
	}
}

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
	ret->set_size(nxn, _MAX(nyn,1), _MAX(nzn,1));
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
		                  sum += (*this)((inxold+m1+nx)%nx,(inyold+m2+ny)%ny)*q; w+=q;}}
		                }else{
                                  for (int m2 =kbmin; m2 <=kbmax; m2++){ for (int m1 =kbmin; m1 <=kbmax; m1++) {
				  float q =kb.i0win_tab(xold - inxold-m1)*kb.i0win_tab(yold - inyold-m2);
		                  sum += (*this)(inxold+m1,inyold+m2)*q;w+=q;}}
		                }
				(*ret)(ix,iy)=sum/w;
			}
		}
	set_array_offsets(saved_offsets);
	return ret;
}
*/

EMData*
EMData::rot_scale_conv(float ang, float delx, float dely, Util::KaiserBessel& kb, float scale_input) {
	int nxn, nyn, nzn;
	if(scale_input == 0.0f) scale_input = 1.0f;
	//const float scale=0.5;
	float  scale = 0.5*scale_input;
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
	EMData* ret = this->copy_head();
#ifdef _WIN32
	ret->set_size(nxn, _MAX(nyn,1), _MAX(nzn,1));
#else
	ret->set_size(nxn, std::max(nyn,1), std::max(nzn,1));
#endif	//_WIN32 
	//ret->to_zero();  //we will leave margins zeroed.
	if(delx >= 0.0f) { delx = fmod(delx, float(nx));} else {delx = -fmod(-delx, float(nx));}
	if(dely >= 0.0f) { dely = fmod(dely, float(ny));} else {dely = -fmod(-dely, float(ny));}
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

		      if (xold < 0.0f) xold = fmod(float(nx) - fmod(-xold, float(nx)), float(nx));
		      else if (xold > (float) (nx-1) ) xold = fmod(xold, float(nx));
		      if (yold < 0.0f) yold = fmod(float(ny) - fmod(-yold, float(ny)), float(ny));
		      else if (yold > (float) (ny-1) ) yold = fmod(yold, float(ny));

			int inxold = int(Util::round(xold)); int inyold = int(Util::round(yold));
			sum=0.0f;    w=0.0f;
			for (int m1 =kbmin; m1 <=kbmax; m1++) t[m1-kbmin] = kb.i0win_tab(xold - inxold-m1);
			if(inxold <= kbc || inxold >=nx-kbc-2 || inyold <= kbc || inyold >=ny-kbc-2 )  {
				for (int m2 =kbmin; m2 <=kbmax; m2++) { float qt = kb.i0win_tab(yold - inyold-m2);
				for (int m1 =kbmin; m1 <=kbmax; m1++) {
					float q = t[m1-kbmin]*qt;
					sum += (*this)((inxold+m1+nx)%nx,(inyold+m2+ny)%ny)*q; w+=q;}}
		    	} else {
				for (int m2 =kbmin; m2 <=kbmax; m2++) { float qt = kb.i0win_tab(yold - inyold-m2);
			  	for (int m1 =kbmin; m1 <=kbmax; m1++) {
					float q = t[m1-kbmin]*qt;
					sum += (*this)(inxold+m1,inyold+m2)*q; w+=q;}}
		    	}
			(*ret)(ix,iy)=sum/w;
		}
	}
	if (t) free(t);
	set_array_offsets(saved_offsets);
	return ret;
}

EMData* EMData::rot_scale_conv_new(float ang, float delx, float dely, Util::KaiserBessel& kb, float scale_input) {

	int nxn, nyn, nzn;
	
	if (scale_input == 0.0f) scale_input = 1.0f;
	float  scale = 0.5*scale_input;

	if (1 >= ny)
		throw ImageDimensionException("Can't rotate 1D image");
	if (1 < nz) 
		throw ImageDimensionException("Volume not currently supported");
	nxn = nx/2; nyn = ny/2; nzn = nz/2;

	int K = kb.get_window_size();
	int kbmin = -K/2;
	int kbmax = -kbmin;

	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
	EMData* ret = this->copy_head();
#ifdef _WIN32
	ret->set_size(nxn, _MAX(nyn,1), _MAX(nzn,1));
#else
	ret->set_size(nxn, std::max(nyn,1), std::max(nzn,1));
#endif	//_WIN32 
	//ret->to_zero();  //we will leave margins zeroed.
	if(delx >= 0.0f) { delx = fmod(delx, float(nx));} else {delx = -fmod(-delx, float(nx));}
	if(dely >= 0.0f) { dely = fmod(dely, float(ny));} else {dely = -fmod(-dely, float(ny));}
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
	
	float *t = (float*)calloc(kbmax-kbmin+1, sizeof(float));
	
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
			
			xold = xold/2.0;
			yold = yold/2.0;
			(*ret)(ix,iy) = Util::get_pixel_conv_new(nx,ny,1,xold,yold,1,data,kb);
			
		}
	}
	if (t) free(t);
	set_array_offsets(saved_offsets);
	return ret;
}


float  EMData::get_pixel_conv(float delx, float dely, float delz, Util::KaiserBessel& kb) {
//  here counting is in C style, so coordinates of the pixel delx should be [0-nx-1] 

	int K = kb.get_window_size();
	int kbmin = -K/2;
	int kbmax = -kbmin;
	int kbc = kbmax+1;

	float pixel =0.0f;
	float w=0.0f;
	
	delx = fmod(2*delx, float(nx));
	int inxold = int(Util::round(delx));
	if(ny<2) {  //1D
	 		 if(inxold <= kbc || inxold >=nx-kbc-2 )  {
	 //  loop for ends
         		   for (int m1 =kbmin; m1 <=kbmax; m1++) {
	 		     float q = kb.i0win_tab(delx - inxold-m1);
	 		     pixel += (*this)((inxold+m1+nx)%nx)*q;w+=q;}
	 		 }else{
         		   for (int m1 =kbmin; m1 <=kbmax; m1++) {
	 		     float q = kb.i0win_tab(delx - inxold-m1);
	 		     pixel += (*this)(inxold+m1)*q;w+=q;}
	 		 }
	
	} else if(nz<2) {  // 2D
	dely = fmod(2*dely, float(ny));
	int inyold = int(Util::round(dely));
	 		 if(inxold <= kbc || inxold >=nx-kbc-2 || inyold <= kbc || inyold >=ny-kbc-2 )  {
	 //  loop for strips
         		   for (int m2 =kbmin; m2 <=kbmax; m2++){ for (int m1 =kbmin; m1 <=kbmax; m1++) {
	 		     float q = kb.i0win_tab(delx - inxold-m1)*kb.i0win_tab(dely - inyold-m2);
	 		     pixel += (*this)((inxold+m1+nx)%nx,(inyold+m2+ny)%ny)*q;w+=q;}}
	 		 }else{
         		   for (int m2 =kbmin; m2 <=kbmax; m2++){ for (int m1 =kbmin; m1 <=kbmax; m1++) {
	 		     float q = kb.i0win_tab(delx - inxold-m1)*kb.i0win_tab(dely - inyold-m2);
	 		     pixel += (*this)(inxold+m1,inyold+m2)*q;w+=q;}}
	 		 }
	} else {  //  3D
	dely = fmod(2*dely, float(ny));
	int inyold = int(Util::round(dely));
	delz = fmod(2*delz, float(nz));
	int inzold = int(Util::round(delz));
			     //cout << inxold<<"  "<< kbc<<"  "<< nx-kbc-2<<"  "<< endl;
	 		 if(inxold <= kbc || inxold >=nx-kbc-2 || inyold <= kbc || inyold >=ny-kbc-2  || inzold <= kbc || inzold >=nz-kbc-2 )  {
	 //  loop for strips
         		   for (int m3 =kbmin; m3 <=kbmax; m3++){ for (int m2 =kbmin; m2 <=kbmax; m2++){ for (int m1 =kbmin; m1 <=kbmax; m1++) {
	 		     float q = kb.i0win_tab(delx - inxold-m1)*kb.i0win_tab(dely - inyold-m2)*kb.i0win_tab(delz - inzold-m3);
			     //cout << "BB  "<<m1<<"  "<< m2<<"  "<< m3<<"  "<< q<<"  "<< q<<"  "<<(*this)((inxold+m1+nx)%nx,(inyold+m2+ny)%ny,(inzold+m3+nz)%nz)<< endl;
	 		     pixel += (*this)((inxold+m1+nx)%nx,(inyold+m2+ny)%ny,(inzold+m3+nz)%nz)*q;w+=q;}}}
	 		 } else {
         		   for (int m3 =kbmin; m3 <=kbmax; m3++){ for (int m2 =kbmin; m2 <=kbmax; m2++){ for (int m1 =kbmin; m1 <=kbmax; m1++) {
	 		     float q = kb.i0win_tab(delx - inxold-m1)*kb.i0win_tab(dely - inyold-m2)*kb.i0win_tab(delz - inzold-m3);
			     //cout << "OO  "<<m1<<"  "<< m2<<"  "<< m3<<"  "<< q<<"  "<< q<<"  "<<(*this)(inxold+m1,inyold+m2,inzold+m3)<< endl;
	 		     pixel += (*this)(inxold+m1,inyold+m2,inzold+m3)*q;w+=q;}}}
	 		 }
	}
        return pixel/w;
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
	// displacements of (xnew,ynew) from the grid
	float nuxdispl = nuxnew - ixn;
	float nuydispl = nuynew - iyn;
	// set up some temporary weighting arrays
	float* wy0 = new float[kbmax - kbmin + 1];
	float* wy = wy0 - kbmin; // wy[kbmin:kbmax]
	float* wx0 = new float[kbmax - kbmin + 1];
	float* wx = wx0 - kbmin;
	for (int i = kbmin; i <= kbmax; i++) {
		wy[i] = kb.i0win_tab(nuydispl - i);
		//wy[i] = (0 == i) ? 1.f : 0.f; // FIXME: remove after debugging
		wx[i] = kb.i0win_tab(nuxdispl - i);
		//wx[i] = (0 == i) ? 1.f : 0.f; // FIXME: remove after debugging
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
	double wsum = 0.f;
	for (int iy = iymin; iy <= iymax; iy++)
		for (int ix = ixmin; ix <= ixmax; ix++)
			wsum += wx[ix]*wy[iy];
	std::complex<float> result(0.f,0.f);
	if ((ixn >= -kbmin) && (ixn <= nhalf-1-kbmax)
			&& (iyn >= -nhalf-kbmin) && (iyn <= nhalf-1-kbmax)) {
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
				if ((ixt > nhalf) || (ixt < -nhalf)) {
					ixt = Util::sgn(ixt)*(nxreal - abs(ixt));
					iyt *= -1;
					mirror = !mirror;
				}
				if ((iyt >= nhalf) || (iyt < -nhalf)) {
					if (ixt != 0) {
						ixt = -ixt;
						iyt = Util::sgn(iyt)*(nxreal-abs(iyt));
						mirror = !mirror;
					} else {
						iyt -= Util::sgn(iyt)*nxreal;
					}
				}
				if (ixt < 0) {
					ixt = -ixt;
					iyt = -iyt;
					mirror = !mirror;
				}
				if (iyt == nhalf) iyt = -nhalf;
				float w = wx[ix]*wy[iy];
				std::complex<float> val = this->cmplx(ixt,iyt);
				if (mirror) 
					result += conj(val)*w;
				else
					result += val*w;
			}
		}
	}
	if (flip) 
		result = conj(result)/static_cast<float>(wsum);
	else
		result /= static_cast<float>(wsum);
	delete [] wx0;
	delete [] wy0;
	return result;
}

void EMData::center_padded() {
	int npad = get_attr("npad");
	if (1 == npad) return;
	int nxreal = nx;
	if (is_fftpadded())
		nxreal = nx - 2 + int(is_fftodd());
	EMData& self = *this;
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets();
	int nxorig = nxreal/npad;
	int nyorig = ny/npad;
	int nzorig = nz/npad;
	int nxcorner = (nxreal - nxorig)/2 + nxorig%2;
	int nycorner = (ny - nyorig)/2 + nyorig%2;
	int nzcorner = (nz - nzorig)/2 + nzorig%2;
	switch (get_ndim()) {
		case 1:
			// 1-d, so no padding along y or z
			nyorig = ny; nzorig = nz;
			nycorner = nzcorner = 0;
			break;
		case 2:
			// 2-d, so no padding along z
			nzorig = nz;
			nzcorner = 0;
			break;
		case 3:
			break;
		default:
			throw ImageDimensionException("center_padded needs a 1-,"
					                      "2-, or 3-d image.");
	}
	for (int iz = nzorig-1; iz >= 0; iz--)
		for (int iy = nyorig-1; iy >= 0; iy--) 
			for (int ix = nxorig-1; ix >= 0; ix--) 
				std::swap(self(nxcorner+ix,nycorner+iy,nzcorner+iz),
						  self(ix,iy,iz));
	set_array_offsets(saved_offsets);
}

/** Helper function for EMData::fft_shuffle, below */
inline void swapx(float* a, float* b, float* temp, size_t nbytes) {
	memcpy(temp, a, nbytes);
	memcpy(a, b, nbytes);
	memcpy(b, temp, nbytes);
}

void EMData::fft_shuffle() {
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
	done_data();
	delete[] temp;
}

EMData* EMData::fouriergridrot2d(float ang, Util::KaiserBessel& kb) {
	if (2 != get_ndim())
		throw ImageDimensionException("fouriergridrot2d needs a 2-D image.");
	if (!is_complex()) 
		throw ImageFormatException("fouriergridrot2d requires a fourier image");
	int nxreal = nx - 2 + int(is_fftodd());
	if (nxreal != ny)
		throw ImageDimensionException("fouriergridrot2d requires ny == nx(real)");
	if (0 != nxreal%2)
		throw ImageDimensionException("fouriergridrot2d needs an even image.");
	int nxhalf = nxreal/2;
	//cmplx(0,0) = 0.;
	if (!is_shuffled()) 
		fft_shuffle();

	//float nxhalf2 = nxhalf*float(nxhalf);
	int nyhalf = ny/2;
	EMData* result = copy();
	set_array_offsets(0,-nyhalf);
	result->set_array_offsets(0,-nyhalf);
	float cang = cos(ang);
	float sang = -sin(ang);
	for (int iy = -nyhalf; iy < nyhalf; iy++) {
		float ycang = iy*cang;
		//float ysang = -iy*sang;
		float ysang = iy*sang;
		//float iy2 = iy*float(iy);
		for (int ix = 0; ix <= nxhalf; ix++) {
			//float ix2 = ix*float(ix);
#if 0 // old version
			if (ix2 + iy2 <= nxhalf2) {
				float nuyold = ix*sang + ycang;
				float nuxold = ix*cang + ysang;
//				result->cmplx(ix,iy) = extractpoint(nuxold,nuyold,kb);
			} else {
//				result->cmplx(ix,iy) = complex<float>(0.f,0.f);
			}
#endif // 0
			float nuyold = -ix*sang + ycang;
			float nuxold = ix*cang + ysang;
			result->cmplx(ix,iy) = extractpoint(nuxold,nuyold,kb);
		}
	}
	result->set_array_offsets();
	result->fft_shuffle(); // reset to an unshuffled result
	result->done_data();
	set_array_offsets();
	fft_shuffle(); // reset to an unshuffled complex image
	//result->cmplx(0,0) = 0.;
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
		float wz = kb.sinhwin(iz-nz/2);
		for (int iy=0; iy < ny; iy++) {
			float wy = kb.sinhwin(iy-ny/2);
			for (int ix=0; ix < nx; ix++) {
				float wx = kb.sinhwin(ix-nx/2);
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

EMData*  
EMData::extractplane(const Transform3D& tf, Util::KaiserBessel& kb) {
	if (!is_complex()) 
		throw ImageFormatException("extractplane requires a fourier image");
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
	int kbsize = kb.get_window_size();
	int kbmin = -kbsize/2;
	int kbmax = -kbmin;
	float* wy0 = new float[kbmax - kbmin + 1];
	float* wy = wy0 - kbmin; // wy[kbmin:kbmax]
	float* wx0 = new float[kbmax - kbmin + 1];
	float* wx = wx0 - kbmin;
	float* wz0 = new float[kbmax - kbmin + 1];
	float* wz = wz0 - kbmin;
	float rim = nhalf*float(nhalf);
	int count = 0;
	float wsum = 0.f;
	Transform3D tftrans = tf; // need transpose of tf here for consistency
	tftrans.transpose();      // with spider
	for (int jy = -nhalf; jy < nhalf; jy++) {
		for (int jx = 0; jx <= nhalf; jx++) {
			Vec3f nucur(jx, jy, 0.f);
			Vec3f nunew = tftrans*nucur;
			float xnew = nunew[0], ynew = nunew[1], znew = nunew[2];
			if (xnew*xnew+ynew*ynew+znew*znew <= rim) {
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
				if (ixn >= -kbmin && ixn <= nhalf-1-kbmax
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
								if (mirror) 
									btq += conj(cmplx(ixt,iyt,izt))*wg;
								else
									btq += cmplx(ixt,iyt,izt)*wg;
								wsum += wg;
							}
						}
					}
				}
				if (flip) 
					res->cmplx(jx,jy) = conj(btq);
				else
					res->cmplx(jx,jy) = btq;
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


bool EMData::peakcmp(const Pixel& p1, const Pixel& p2) {
    return (p1.value > p2.value);
}

ostream& operator<< (ostream& os, const Pixel& peak) {
    os <<  peak.x <<  peak.y << peak.z
       << peak.value;
    return os;
}

vector<float> EMData::peak_search(int ml, float invert)
{
 	 EMData& buf = *this;
	 vector<Pixel> peaks;
 	 int img_dim;
 	 int i,j,k,itx,ity,itz;
 	 int i__1,i__2;
 	 int j__1,j__2;
 	 int k__1,k__2;
 	 bool peak_check;
 	 img_dim=buf.get_ndim();
 	 vector<int>ix,jy,kz;
	 vector<float>res;
 	 int nx = buf.get_xsize();
 	 int ny = buf.get_ysize();
 	 int nz = buf.get_zsize();
	 if(invert <= 0) 
	     { invert=-1.;}
	 else 
	     { invert=1. ;}
 	 switch (img_dim)
     {
 	 case(1):
		for(i=0;i<=nx-1;++i)
 	  	{
 	   		i__1=(i-1+nx)%nx;
 	   		i__2=(i+1)%nx;
 	      	peak_check=buf(i)*invert>buf(i__1)*invert && buf(i)*invert>buf(i__2)*invert;
	 	  	if(peak_check)
		  		{peaks.push_back(Pixel(i, 0, 0, buf(i)*invert));}  
	 	}
 	 break;
 	 case(2):
		for(j=0;j<=ny-1;++j)
 	    {  
 	    	j__1=(j-1+ny)%ny;
 		j__2=(j+1)%ny;
 	        for(i=0;i<=nx-1;++i)
			{ 
				i__1=(i-1+nx)%nx;
			  	i__2=(i+1)%nx;
			  	peak_check=(buf(i,j)*invert>buf(i,j__1)*invert) && (buf(i,j)*invert>buf(i,j__2)*invert);
			  	peak_check=peak_check && (buf(i,j)*invert>buf(i__1,j)*invert) && (buf(i,j)*invert>buf(i__2,j)*invert);
			  	peak_check=peak_check && (buf(i,j)*invert>buf(i__1,j__1)*invert) && ((buf(i,j)*invert)> buf(i__1,j__2)*invert);
			  	peak_check=peak_check && (buf(i,j)*invert>buf(i__2,j__1)*invert) && (buf(i,j)*invert> buf(i__2,j__2)*invert);
			  	if(peak_check)
 			    {peaks.push_back(Pixel(i, j, 0, buf(i,j)*invert));}
			}
 		}
 	 break;
 	 case(3):
		for(k=0;k<=nz-1;++k)
	   	{  
	   		kz.clear();
	      	k__1=(k-1+nz)%nz;
	      	k__2=(k+1)%nz;
	      	kz.push_back(k__1);
	      	kz.push_back(k);
	      	kz.push_back(k__2);
 	     	for(j=0;j<=ny-1;++j)
		    {
		    	jy.clear();
		     	j__1=(j-1+ny)%ny;
 	            j__2=(j+1)%ny;
		     	jy.push_back(j__1);
		     	jy.push_back(j);
		     	jy.push_back(j__2);
 	            for(i=0;i<=nx-1;++i)
			   	{ 
			   		ix.clear();
			     	i__1=(i-1+nx)%nx;
			     	i__2=(i+1)%nx;
			     	ix.push_back(i__1);
			     	ix.push_back(i);
			     	ix.push_back(i__2);
			     	peak_check=true ;
				 	for(itx=0;itx<= 2; ++itx)
				    {  
				    	for(ity=0;ity<= 2; ++ity)
					    {  
					    	for(itz=0;itz<= 2; ++itz) 
						 	{ 
								if((buf(i,j,k)*invert<=buf(ix[itx],jy[ity],kz[itz])*invert) && i !=ix[itx] && j !=jy[ity] && k !=kz[itz])			   	
						        {peak_check=false;}						   
						  	}
					     }		 
					}
					if(peak_check)
				    {
						peaks.push_back(Pixel(i, j, k, buf(i,j,k)*invert));
					} 
			      }
			   }
			}
		break;
     	}
   if(peaks.begin()!=peaks.end())
   	     	    	 
   { sort(peaks.begin(),peaks.end(), peakcmp);   
    int count=0;
    float xval=(*peaks.begin()).value;  
    	   for (vector<Pixel>::iterator it = peaks.begin(); it != peaks.end(); it++)  
    	   {count=count+1;
    		if(count<=ml)
    		   {
    		   res.push_back((*it).value);
    		   res.push_back((*it).x);
    		   if(img_dim!=1) res.push_back((*it).y);
    			 
    		   if(nz!=1) {res.push_back((*it).z);}   
    			 
    		   if(xval != 0.0)
    			 {res.push_back((*it).value/xval);}
    		   else
    			 {res.push_back((*it).value);}
    		   res.push_back((*it).x-float(int(nx/2)));
    		   if(img_dim!=1)
    		       {res.push_back((*it).y-float(int(ny/2)));}
    		   if(nz!=1)
    		      {res.push_back((*it).z-float(nz/2));} 
    		   
    		 }
    	    }
    	    res.insert(res.begin(),1,img_dim);
	  } 
 else 
    {
    res.push_back(buf(0,0,0));
    res.insert(res.begin(),1,0.0);     
	                      }	   
  return res;
 }	
	  
#define rdata(i,j,k) rdata[(i-1)+((j-1)+(k-1)*ny)*nx]
#define X(i) X[i-1]
#define Y(j) Y[j-1]
#define Z(k) Z[k-1]
vector<float> EMData::phase_cog()
{	
	vector<float> ph_cntog;
	int i=1,j=1,k=1;
	float C=0.f,S=0.f,P=0.f,F1=0.f,SNX;
	if (get_ndim()==1)
		{P = 8*atan(1.0f)/nx;
		for (i=1;i<=nx;i++)
			{C += cos(P * (i-1)) * rdata(i,j,k);
		         S += sin(P * (i-1)) * rdata(i,j,k);}
		F1 = atan2(S,C);
		if (F1 < 0.0){ F1 += 8*atan(1.0f); }
		SNX = F1/P +1.0;
		SNX = SNX - ((nx/2)+1);
		ph_cntog.push_back(SNX);
#ifdef _WIN32
		ph_cntog.push_back(Util::round(SNX));
#else
		ph_cntog.push_back(round(SNX));
#endif //_WIN32
		}
	else if (get_ndim()==2)
#ifdef _WIN32
		{
		float SNY;
		float T=0.0f;
		vector<float> X;
		X.resize(nx);
#else
		{float SNY,X[nx],T=0.f;
#endif	//_WIN32
		for ( i=1;i<=nx;i++)
			{X(i)=0.0;}			
                 P = 8*atan(1.0f)/ny;
		 for(j=1;j<=ny;j++)
			{T=0.f;
			 for(i=1;i<=nx;i++)
				{T += rdata(i,j,k);
				 X(i)+=rdata(i,j,k);}
			 C += cos(P*(j-1))*T;
			 S += sin(P*(j-1))*T;}
		 F1=atan2(S,C);
		 if(F1<0.0){ F1 += 8*atan(1.0f); }
		 SNY = F1/P +1.0;
		 C=0.f;S=0.f;
		 P = 8*atan(1.0f)/nx;
		 for(i=1;i<=nx;i++)
			{C += cos(P*(i-1))*X(i);
			 S += sin(P*(i-1))*X(i);}
	         F1=atan2(S,C);
		 if(F1<0.0){ F1 += 8*atan(1.0f); }
		 SNX = F1/P +1.0;
		 SNX = SNX - ((nx/2)+1);
		 SNY = SNY - ((ny/2)+1);
		 ph_cntog.push_back(SNX); ph_cntog.push_back(SNY);	
#ifdef _WIN32
		 ph_cntog.push_back(Util::round(SNX)); ph_cntog.push_back(Util::round(SNY));
#else
		 ph_cntog.push_back(round(SNX)); ph_cntog.push_back(round(SNY));
#endif	//_WIN32
		}
	else
#ifdef _WIN32
		{float val=0.f,sum1=0.f, SNY,SNZ;
		vector<float> X;
		X.resize(nx);
		vector<float> Y;
		Y.resize(ny);
		vector<float> Z;
		Z.resize(nz);
#else
		{float val=0.f,sum1=0.f,X[nx],Y[ny],Z[nz],SNY,SNZ;
#endif	//_WIN32
		 for (i=1;i<=nx;i++)
			{X(i)=0.0;}
		 for (j=1;j<=ny;j++)
			{Y(j)=0.0;}
		 for (k=1;k<=nz;k++)
			{Z(k)=0.0;}
		 for(k=1;k<=nz;k++)
			{for(j=1;j<=ny;j++)
				{sum1=0.f;
				 for(i=1;i<=nx;i++)
					{val = rdata(i,j,k);
					 sum1 += val;
					 X(i) += val;}
				 Y(j) += sum1;
				 Z(k) += sum1;}
			}
		 P = 8*atan(1.0f)/nx;
		 for (i=1;i<=nx;i++)
			{C += cos(P*(i-1))*X(i);
			 S += sin(P*(i-1))*X(i);}
		 F1=atan2(S,C);
		 if(F1<0.0){ F1 += 8*atan(1.0); }
		 SNX = F1/P +1.0;
		 C=0.f;S=0.f;
		 P = 8*atan(1.0f)/ny;
		 for(j=1;j<=ny;j++)
			{C += cos(P*(j-1))*Y(j);
			 S += sin(P*(j-1))*Y(j);}
		 F1=atan2(S,C);
		 if(F1<0.0){ F1 += 8*atan(1.0f); }
		 SNY = F1/P +1.0;
		 C=0.f;S=0.f;
		 P = 8*atan(1.0f)/nz;
		 for(k=1;k<=nz;k++)
			{C += cos(P*(k-1))*Z(k);
		         S += sin(P*(k-1))*Z(k);}
		 F1=atan2(S,C);
		 if(F1<0.0){ F1 += 8*atan(1.0f); }
		 SNZ = F1/P +1.0;	
		 SNX = SNX - ((nx/2)+1);
		 SNY = SNY - ((ny/2)+1);
		 SNZ = SNZ - ((nz/2)+1);		
		 ph_cntog.push_back(SNX); ph_cntog.push_back(SNY); ph_cntog.push_back(SNZ);
#ifdef _WIN32
		 ph_cntog.push_back(Util::round(SNX)); ph_cntog.push_back(Util::round(SNY)); ph_cntog.push_back(Util::round(SNZ));
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
#define R (0.61803399)
#define C (1.f-R)
float EMData::find_3d_threshold(float mass,float pixel_size)
{
	/* Exception Handle */
	if(get_ndim()!=3)
		throw ImageDimensionException("The image should be 3D");
	/* ===============================================================*/
	
	/* Calculation of the volume of the voxels */
	float density_1_mole,vol_1_mole,vol_angstrom;
	int vol_voxels;
	density_1_mole = (float)(mass*1000.0f)/avagadro;
	vol_1_mole = density_1_mole/density_protein;
	vol_angstrom = vol_1_mole*(double)pow((double)pow(10.0,8),3);
	vol_voxels = static_cast<int> (vol_angstrom/(double)pow(pixel_size,3));
	/* ===============================================================*/

	
	float thr1 = get_attr("maximum");
	float thr3 = get_attr("minimum");
	float thr2 = (thr1-thr3)/2 + thr3;
	int size = nx*ny*nz;
	float x0 = thr1,x3 = thr3,x1,x2,THR=0;
	
	#ifdef _WIN32
		int ILE = _MIN(nx*ny*nx,_MAX(1,vol_voxels));
	#else
		int ILE = std::min(nx*ny*nx,std::max(1,vol_voxels));
	#endif	//_WIN32
	
	if (abs(thr3-thr2)>abs(thr2-thr1))
	{	x1=thr2;
		x2=thr2+C*(thr3-thr2);}
	else
	{	x2=thr2;
		x1=thr2-C*(thr2-thr1);	}
		
	int cnt1=0,cnt2=0;
	for (int i=0;i<size;i++)
	{	if(rdata[i]>=x1)
			cnt1++;
		if(rdata[i]>=x2)
			cnt2++;
	}
	float LF1 = cnt1 - ILE;
	float F1 = LF1*LF1;
	float LF2 = cnt2 - ILE;
	float F2 = LF2*LF2;
	
	while ((LF1 != 0 || LF2 != 0) && (fabs(LF1-LF2) >= 1.f) && (abs(x1-x2) > (double)pow(10.0,-5) && abs(x1-x3) > (double)pow(10.0,-5) && abs(x2-x3) > (double)pow(10.0,-5)))
	{
		if(F2 < F1)
		{
			x0=x1;
			x1=x2;
			x2 = R*x1 + C*x3;
			F1=F2;
			int cnt=0;
			for(int i=0;i<size;i++)
				if(rdata[i]>=x2)
					cnt++;
			LF2 = cnt - ILE;
			F2 = LF2*LF2;
		}
		else
		{
			x3=x2;
			x2=x1;
			x1=R*x2 + C*x0;
			F2=F1;
			int cnt=0;
			for(int i=0;i<size;i++)
				if(rdata[i]>=x1)
					cnt++;
			LF1 = cnt - ILE;
			F1 = LF1*LF1;
		}
	}
	
	if(F1 < F2)
	{
		ILE = static_cast<int> (LF1 + ILE);
		THR = x1;
	}
	else
	{
		ILE = static_cast<int> (LF2 + ILE);
		THR = x2;
	}
	return THR;
	
}
#undef avagadro
#undef density_protein
#undef R
#undef C			
vector<float> EMData::peak_ccf(float hf_p)
{   
    	 EMData & buf = *this;
	 vector<Pixel> peaks;	 
	 int half=int(hf_p);
	 float hf_p2=hf_p*hf_p;
 	 int i,j;
 	 int i__1,i__2;
 	 int j__1,j__2;
 	 bool peak_found;
	 bool not_overlap;
	 vector<float>res;
 	 int nx = buf.get_xsize()-half;
 	 int ny = buf.get_ysize()-half; 
	 int n_peak=0;
	 for(i=half;i<=nx;++i)																																																													 
 		{	
			i__1=i-1;		
 		  	i__2=i+1;	 																	 							    
			for (j=half;j<=ny;++j)
				{  
			    	 	j__1=j-1;
				 	j__2=j+1;   																      
					peak_found=(buf(i,j)>buf(i,j__1)) && (buf(i,j)>buf(i,j__2));
					peak_found=peak_found && (buf(i,j)>buf(i__1,j)) && (buf(i,j)>buf(i__2,j));
					peak_found=peak_found && (buf(i,j)>buf(i__1,j__1)) && ((buf(i,j))> buf(i__1,j__2));
					peak_found=peak_found && (buf(i,j)>buf(i__2,j__1)) && (buf(i,j)> buf(i__2,j__2));				
			  		if(peak_found)
 			       			{	
							if (peaks.size()==0) 
								{	 
									peaks.push_back(Pixel(i,j,0,buf(i,j)));
							 		n_peak=n_peak+1;
								}
							else									
								{	
									not_overlap=true;							
									bool higher_peak=false;	
						 			int size=peaks.size();
						 			for ( int kk= 0; kk< size; kk++)								 
										{	
											vector<Pixel>::iterator it= peaks.begin()+kk;
											float radius=((*it).x-float(i))*((*it).x-float(i))+((*it).y-float(j))*((*it).y-float(j));																	
							  				if (radius <= hf_p2 )
												{	
													not_overlap=false;
													if( buf(i,j) > (*it).value)
														{	
															peaks.erase(it);												 													 
															higher_peak=true;																							
														}																																	
									      												
													else
														{	
															higher_peak=false; 
															break;
														}
												}
																
										}
													
									if(not_overlap|higher_peak)
										{										
											n_peak=n_peak+1;
								 			peaks.push_back(Pixel(i,j,0,buf(i,j)));
										}
								}
							                                
					}
			}						
	}	
	if(peaks.size()>=1)
   	     	    	 
  		{	sort(peaks.begin(),peaks.end(), peakcmp); 
			for (vector<Pixel>::iterator it = peaks.begin(); it != peaks.end(); it++) 
				{
					res.push_back((*it).value);
					res.push_back((*it).x);
					res.push_back((*it).y);			
			
		 		}
		}
	else
		{	
			res.push_back(buf(0,0,0));
 			res.insert(res.begin(),1,0.0); 	
		}
	return res;
}

EMData* EMData::get_pow(float n_pow)
{   
	vector<int> saved_offsets = get_array_offsets();
	EMData* buf_new = new EMData();
 	int nx = this->get_xsize();
 	int ny = this->get_ysize(); 
	int nz = this->get_zsize(); 
	buf_new->set_size(nx,ny,nz);
	float *in = this->get_data();
	float *out = buf_new->get_data();
	for(int i=0; i<nx*ny*nz; i++) out[i]=pow(in[i],n_pow);
	return buf_new;
}						

EMData* EMData::extractline(Util::KaiserBessel& kb,float nuxnew,float nuynew) 
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

	int count = 0;
	float wsum = 0.f;
	bool flip = (nuxnew < 0.f);

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
			}
			else {
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
								iyt = Util::sgn(iyt)
									  *(n - abs(iyt));
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
						if (mirror) 
							btq += conj(cmplx(ixt,iyt))*wg;
						else
							btq += cmplx(ixt,iyt)*wg;
							wsum += wg;
					}
				}
			}
		if (flip) 
			res->cmplx(jx) = conj(btq);
		else
			res->cmplx(jx) = btq;
	}
	for (int jx = 0; jx <= nhalf; jx++) 
		res->cmplx(jx) *= count/wsum;
	
	delete[] wx0; delete[] wy0;
	set_array_offsets(saved_offsets);
	res->set_array_offsets(0,0,0);
	//res->set_shuffled(true);
	return res;
}

EMData* EMData::ctf_img(int nx, int ny, int nz,float dz,float ps,float voltage,float cs, float wgh,float b_factor,float dza, float azz, float sign)
{               
	int   lsm;
	int   ix, iy, iz;
	int   i,  j, k;    
	int   nr2, nl2;
	float  dzz, az, ak;
	float  scx, scy, scz;
	int offset = 2 - nx%2;  
	lsm = nx + offset;		     	   
	EMData* ctf_img1 = new EMData();
	ctf_img1->set_size(lsm, ny, nz);
	float freq=1./(2.*ps);		    
	scx=2./float(nx);
	if(ny<=1) scy=2./ny; else scy=0.0;
	if(nz<=1) scz=2./nz; else scz=0.0;
	nr2 = ny/2 ;
	nl2 = nz/2 ;
	for ( k=0; k<nz;k++) {
	       iz = k;  if(k>nl2) iz=k-nz;
	       for ( j=0; j<ny;j++) { 
	     	     iy = j;  if(j>nr2) iy=j - ny;
	     	     for ( i=0; i<lsm/2; i++) {
	     		   ix=i;
	     		   ak=pow(ix*ix*scx*scx+iy*scy*iy*scy+iz*scz*iz*scz, 0.5f)*freq;
	     		   if(ak!=0) az=0.0; else az=M_PI;
	     		   dzz=dz+dza/2.*sin(2*(az-azz*M_PI/180.));
			       (*ctf_img1) (i*2,j,k)   = Util::tf(dzz, ak, voltage, cs, wgh, b_factor, sign);
	     		   (*ctf_img1) (i*2+1,j,k) = 0.0f;
	     	     }
	     	     
	       }

	}
		ctf_img1->done_data();
	//	ctf_img1->set_complex(true);//
		ctf_img1->attr_dict["is_complex"] = 1;
		ctf_img1->attr_dict["is_ri"] = 1;
		if(nx%2==0) ctf_img1->set_fftodd(false); else ctf_img1->set_fftodd(true);		
		return ctf_img1;
			 			 
} 		
