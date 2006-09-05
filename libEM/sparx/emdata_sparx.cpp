/**
 * $Id$
 */
#include "emdata.h"
#include <iostream>
#include "math.h"
 
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <vector>
using std::vector;
using std::cout;
using namespace EMAN;
using namespace std;

EMData *EMData::real2FH(float OverSamplekB) // PRB
{
	int nx=get_xsize();
	int ny=get_ysize();
	int nz=get_zsize();
	int Center = (int) floor( (nx+1.0)/2.0 +.01);
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

	int   *PermMatTr           = new int[CountMax];
	float *RValsSorted         = new float[CountMax];
	float *weightofkValsSorted = new float[CountMax];
	int   *SizeReturned        = new int[1];
	Util::Radialize(PermMatTr, RValsSorted,weightofkValsSorted,Size, SizeReturned);
	int RIntMax= SizeReturned[0];  // replaces CountMax; the latter should now never be used.
//	kVec2Use = (0:1/OverSamplek:RValsSorted(RIntMax)+1/OverSamplek); %   in pixels  (otherwise need *2*pi/Size)

	int mMax = (int) floor( ScalFactor*RValsSorted[RIntMax-1]+10.0);

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
						float mmFac = mmFac;
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


EMData* EMData::rotavg()
{
	ENTERFUNC;

	if (ny<2 && nz <2) {
		LOGERR("No 1D images.");
		throw ImageDimensionException("No 1D images!");
	}
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(-nx/2,-ny/2,-nz/2);
#ifdef _WIN32
	int rmax = _MIN(nx/2 + nx%2, ny/2 + ny%2);
#else
	int rmax = std::min(nx/2 + nx%2, ny/2 + ny%2);
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

#define rdata(i,j,k) rdata[(i-1)+((j-1)+(k-1)*ny)*nx]
#define square(x) ((x)*(x))
vector<float> EMData::cog()
{
	
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



EMData* EMData::symvol(string symmetry) {
	ENTERFUNC;
	int nsym = Transform3D::get_nsym(symmetry); // number of symmetries
	Transform3D sym;
	int llim = -nx/2;
	int ulim = (nx/2) -1 + (nx % 2);
	// set up output volume
	EMData* svol = this->copy();
	// set up coord grid
	const int nsize = 27;
	int x[nsize], y[nsize], z[nsize];
	float f[nsize];
	for (int i = 0; i < nsize; i+=3) {
		x[i] = -1;
		x[i+1] = 0;
		x[i+2] = 1;
		int imod = (i/3) % 3;
		y[i] = imod - 1;
		y[i+1] = imod - 1;
		y[i+2] = imod - 1;
	}
	for (int i = 0; i < nsize; i++) z[i] = (i/9) - 1;
	// calc radius within which the rotation is valid
	int iradrt = (0 == nx % 2) ? nx/2 - 1 : nx / 2;
	int iradi = (iradrt-1)*(iradrt-1);
	// actual work -- loop over symmetries and symmetrize
	for (int isym = 0; isym < nsym; isym++) {
		Transform3D rm = sym.get_sym(symmetry, isym);
		if ((1.0 == rm[0][0]) && (1.0 == rm[1][1]) && (1.0 == rm[2][2])) {
			// symmetry has to include the identity matrix and this was taken care of by the copy statement above
			/*for (int iz = 0; iz < nz; iz++)
				for (int iy = 0; iy < ny; iy++)
					for (int ix = 0; ix < nx; ix++)
						(*svol)(ix,iy,iz) += (*this)(ix,iy,iz);*/
		} else {
			// symmetry is something interesting
			Vec3f qrt, qr;
			for (int iz = llim; iz <= ulim; iz++) {
				for (int iy = llim; iy <= ulim; iy++) {
					qrt[0] = rm[0][1]*iy + rm[0][2]*iz;
					qrt[1] = rm[1][1]*iy + rm[1][2]*iz;
					qrt[2] = rm[2][1]*iy + rm[2][2]*iz;
					for (int ix = llim; ix <= ulim; ix++) {
						int icrd = ix*ix + iy*iy + iz*iz;
						if (icrd <= iradi) {
							qr[0] = qrt[0] + rm[0][0]*ix;
							qr[1] = qrt[1] + rm[1][0]*ix;
							qr[2] = qrt[2] + rm[2][0]*ix;
							// iox -- integer location in -nx/2...nx/2
							int iox = static_cast<int>(floorf(qr[0]));
							int ioy = static_cast<int>(floorf(qr[1]));
							int ioz = static_cast<int>(floorf(qr[2]));
							// dx -- offset from integer array
							float dx = qr[0] - iox;
							float dy = qr[1] - ioy;
							float dz = qr[2] - ioz;
							// find intensities on 3x3x3 grid
							for (int i = 0; i < nsize; i++) {
								int jx = iox + x[i] - llim;
								int jy = ioy + y[i] - llim;
								int jz = ioz + z[i] - llim;
								f[i] = (*this)(jx,jy,jz);
							}
							// eval intensity at px, py, pz
							int jx = ix - llim;
							int jy = iy - llim;
							int jz = iz - llim;
							(*svol)(jx,jy,jz) += Util::triquad(dx, dy, dz, f);
						} else {
							// rotated position is outside volume
							int jx = ix - llim;
							int jy = iy - llim;
							int jz = iz - llim;
							(*svol)(jx,jy,jz) += (*this)(jx,jy,jz);
						}
					}
				}
			}
		}
	}
	// normalize
	for (int iz = 0; iz < nz; iz++)
		for (int iy = 0; iy < ny; iy++)
			for (int ix = 0; ix < nx; ix++)
				(*svol)(ix,iy,iz) /= nsym;
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
	int nc = nx/2 + 1;
	for (int iz = 1; iz <= nz; iz++) { float yy = (iz-nc)*(iz-nc);
		for (int iy = 1; iy <=ny; iy++) { float xx = yy + (iy-nc)*(iy-nc);
			for (int ix = 1; ix <= nx; ix++) {
				if(xx+float((ix-nc)*(ix-nc)) > r2 ) {
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
		          EMArray<int>& nr, EMData* bi, const Transform3D& tf) {//std::cout<<"   onelinenn  "<<j<<"  "<<n<<"  "<<n2<<"  "<<std::endl;
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
					nr(ixn,iya,iza)++;
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
					nr(-ixn,iyt,izt)++;
				}
			}

		}
	}
}

void EMData::nn(EMArray<int>& nr, EMData* myfft, const Transform3D& tf) {
	ENTERFUNC;
	int nxc = attr_dict["nxc"]; // # of complex elements along x
	// let's treat nr, bi, and local data as matrices
	vector<int> saved_offsets = get_array_offsets();
	vector<int> myfft_saved_offsets = myfft->get_array_offsets();
	set_array_offsets(0,1,1);
	myfft->set_array_offsets(0,1);
	// loop over frequencies in y
	for (int iy = -ny/2 + 1; iy <= ny/2; iy++) onelinenn(iy, ny, nxc, nr, myfft, tf);
	set_array_offsets(saved_offsets);
	myfft->set_array_offsets(myfft_saved_offsets);
	EXITFUNC;
}

void EMData::symplane0(EMArray<int>& w) {
	ENTERFUNC;
	int nxc = attr_dict["nxc"];
	int n = nxc*2;
	// let's treat the local data as a matrix
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,1,1);
	for (int iza = 2; iza <= nxc; iza++) {
		for (int iya = 2; iya <= nxc; iya++) {
			cmplx(0,iya,iza) += conj(cmplx(0,n-iya+2,n-iza+2));
			w(0,iya,iza) += w(0,n-iya+2,n-iza+2);
			cmplx(0,n-iya+2,n-iza+2) = conj(cmplx(0,iya,iza));
			w(0,n-iya+2,n-iza+2) = w(0,iya,iza);
			cmplx(0,n-iya+2,iza) += conj(cmplx(0,iya,n-iza+2));
			w(0,n-iya+2,iza) += w(0,iya,n-iza+2);
			cmplx(0,iya,n-iza+2) = conj(cmplx(0,n-iya+2,iza));
			w(0,iya,n-iza+2) = w(0,n-iya+2,iza);
		}
	}
	for (int iya = 2; iya <= nxc; iya++) {
		cmplx(0,iya,1) += conj(cmplx(0,n-iya+2,1));
		w(0,iya,1) += w(0,n-iya+2,1);
		cmplx(0,n-iya+2,1) = conj(cmplx(0,iya,1));
		w(0,n-iya+2,1) = w(0,iya,1);
	}
	for (int iza = 2; iza <= nxc; iza++) {
		cmplx(0,1,iza) += conj(cmplx(0,1,n-iza+2));
		w(0,1,iza) += w(0,1,n-iza+2);
		cmplx(0,1,n-iza+2) = conj(cmplx(0,1,iza));
		w(0,1,n-iza+2) = w(0,1,iza);
	}
	EXITFUNC;
}


//  Helper functions for method nn4_ctf
void EMData::onelinenn_ctf(int j, int n, int n2, 
		          EMArray<float>& w, EMData* bi, const Transform3D& tf, float dz) {//std::cout<<"   onelinenn_ctf  "<<j<<"  "<<n<<"  "<<n2<<"  "<<std::endl;
	//  CTF parameters
	float ps = 4.2f, voltage = 200.0f, cs= 2.0*1.0e-7f, ww=0.1, b_factor=0.0f, sign=-1.0f;
	float wgh=atan(ww/(1.0-ww)), lambda = 12.398f/std::sqrt(voltage *(1022.f+voltage));
	int jp = (j >= 0) ? j+1 : n+j+1;
	//for(int i = 0; i <= 1; i++){for(int l = 0; l <= 2; l++){std::cout<<"  "<<tf[i][l]<<"  "<<std::endl;}}
	float  dy2 = pow(float(j)/float(n),2);
	// loop over x
	for (int i = 0; i <= n2; i++) {
		if (((i*i+j*j) < n*n/4) && !((0 == i) && (j < 0))) {
			 //	   if ( !((0 == i) && (j < 0))) {
			//  Calculate absolute frequency
			float  ak = std::sqrt(pow(float(i)/float(2*n2),2)+dy2)/ps;
			//float  ctf = Util::tf(dz, freq, 12.398f/sqrt(lambda *(1022.f+lambda)), cs*1.0e-7f, atan(wgh/(1.0-wgh)), b_factor, sign);
			float  ctf = -sin(-M_PI*(dz*lambda*ak*ak-cs*lambda*lambda*lambda*ak*ak*ak*ak/2.)-wgh);
			float xnew = i*tf[0][0] + j*tf[1][0];
			float ynew = i*tf[0][1] + j*tf[1][1];
			float znew = i*tf[0][2] + j*tf[1][2];
			std::complex<float> btq;
			if (xnew < 0.) {
				xnew = -xnew;
				ynew = -ynew;
				znew = -znew;
				btq = conj(bi->cmplx(i,jp))*ctf;
			} else  btq = bi->cmplx(i,jp)*ctf;
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
					cmplx(ixn,iya,iza) += btq;
					//std::cout<<"    "<<j<<"  "<<ixn<<"  "<<iya<<"  "<<iza<<"  "<<btq<<std::endl;
					w(ixn,iya,iza) += ctf*ctf;
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
					w(-ixn,iyt,izt) += ctf*ctf;
				}
			}
		}
	}
}

void EMData::onelinenn_ctf_applied(int j, int n, int n2, 
		          EMArray<float>& w, EMData* bi, const Transform3D& tf, float dz) {//std::cout<<"   onelinenn_ctf  "<<j<<"  "<<n<<"  "<<n2<<"  "<<std::endl;
	//  CTF parameters
	float ps = 4.2f, voltage = 200.0f, cs= 2.0*1.0e-7f, ww=0.1, b_factor=0.0f, sign=-1.0f;
	float wgh=atan(ww/(1.0-ww)), lambda = 12.398f/std::sqrt(voltage *(1022.f+voltage));
	int jp = (j >= 0) ? j+1 : n+j+1;
	//for(int i = 0; i <= 1; i++){for(int l = 0; l <= 2; l++){std::cout<<"  "<<tf[i][l]<<"  "<<std::endl;}}
	float  dy2 = pow(float(j)/float(n),2);
	// loop over x
	for (int i = 0; i <= n2; i++) {
		if (((i*i+j*j) < n*n/4) && !((0 == i) && (j < 0))) {
			 //	   if ( !((0 == i) && (j < 0))) {
			//  Calculate absolute frequency
			float  ak = std::sqrt(pow(float(i)/float(2*n2),2)+dy2)/ps;
			//float  ctf = Util::tf(dz, freq, 12.398f/sqrt(lambda *(1022.f+lambda)), cs*1.0e-7f, atan(wgh/(1.0-wgh)), b_factor, sign);
			float  ctf = -sin(-M_PI*(dz*lambda*ak*ak-cs*lambda*lambda*lambda*ak*ak*ak*ak/2.)-wgh);
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
					cmplx(ixn,iya,iza) += btq;
					//std::cout<<"    "<<j<<"  "<<ixn<<"  "<<iya<<"  "<<iza<<"  "<<btq<<std::endl;
					w(ixn,iya,iza) += ctf*ctf;
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
					w(-ixn,iyt,izt) += ctf*ctf;
				}
			}
		}
	}
}

void
EMData::nn_ctf(EMArray<float>& w, EMData* myfft, const Transform3D& tf, float defocus) {
	ENTERFUNC;
	int nxc = attr_dict["nxc"]; // # of complex elements along x
	// let's treat nr, bi, and local data as matrices
	vector<int> saved_offsets = get_array_offsets();
	vector<int> myfft_saved_offsets = myfft->get_array_offsets();
	set_array_offsets(0,1,1);
	myfft->set_array_offsets(0,1);
	// loop over frequencies in y
	for (int iy = -ny/2 + 1; iy <= ny/2; iy++) onelinenn_ctf(iy, ny, nxc, w, myfft, tf, defocus);
	set_array_offsets(saved_offsets);
	myfft->set_array_offsets(myfft_saved_offsets);
	EXITFUNC;
}

void
EMData::nn_ctf_applied(EMArray<float>& w, EMData* myfft, const Transform3D& tf, float defocus) {
	ENTERFUNC;
	int nxc = attr_dict["nxc"]; // # of complex elements along x
	// let's treat nr, bi, and local data as matrices
	vector<int> saved_offsets = get_array_offsets();
	vector<int> myfft_saved_offsets = myfft->get_array_offsets();
	set_array_offsets(0,1,1);
	myfft->set_array_offsets(0,1);
	// loop over frequencies in y
	for (int iy = -ny/2 + 1; iy <= ny/2; iy++) onelinenn_ctf_applied(iy, ny, nxc, w, myfft, tf, defocus);
	set_array_offsets(saved_offsets);
	myfft->set_array_offsets(myfft_saved_offsets);
	EXITFUNC;
}

void
EMData::symplane0_ctf(EMArray<float>& w) {
	ENTERFUNC;
	int nxc = attr_dict["nxc"];
	int n = nxc*2;
	// let's treat the local data as a matrix
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,1,1);
	for (int iza = 2; iza <= nxc; iza++) {
		for (int iya = 2; iya <= nxc; iya++) {
			cmplx(0,iya,iza) += conj(cmplx(0,n-iya+2,n-iza+2));
			w(0,iya,iza) += w(0,n-iya+2,n-iza+2);
			cmplx(0,n-iya+2,n-iza+2) = conj(cmplx(0,iya,iza));
			w(0,n-iya+2,n-iza+2) = w(0,iya,iza);
			cmplx(0,n-iya+2,iza) += conj(cmplx(0,iya,n-iza+2));
			w(0,n-iya+2,iza) += w(0,iya,n-iza+2);
			cmplx(0,iya,n-iza+2) = conj(cmplx(0,n-iya+2,iza));
			w(0,iya,n-iza+2) = w(0,n-iya+2,iza);
		}
	}
	for (int iya = 2; iya <= nxc; iya++) {
		cmplx(0,iya,1) += conj(cmplx(0,n-iya+2,1));
		w(0,iya,1) += w(0,n-iya+2,1);
		cmplx(0,n-iya+2,1) = conj(cmplx(0,iya,1));
		w(0,n-iya+2,1) = w(0,iya,1);
	}
	for (int iza = 2; iza <= nxc; iza++) {
		cmplx(0,1,iza) += conj(cmplx(0,1,n-iza+2));
		w(0,1,iza) += w(0,1,n-iza+2);
		cmplx(0,1,n-iza+2) = conj(cmplx(0,1,iza));
		w(0,1,n-iza+2) = w(0,1,iza);
	}
	EXITFUNC;
}








EMData*
EMData::rot_trans2D(float ang, float delx, float dely) {
	if (1 >= ny) 
		throw ImageDimensionException("Can't rotate 1D image");
	if (1 < nz) 
		throw ImageDimensionException("Volume not currently supported");
	if (0.f == ang) {
		EMData* ret = copy();
		return ret;
	}
	update();
	float background = get_attr("mean");
	if (ang > pi) ang -= static_cast<float>(twopi);
	if (ang < -pi) ang += static_cast<float>(twopi);
	float cang = cos(ang);
	float sang = sin(ang);
	EMData* ret = copy_head();
	// center of the image
	int xc = nx/2;
	int yc = ny/2;
	// shift center for rotation (if desired)
	float shiftxc = xc + delx;
	float shiftyc = yc + dely;
	for (int iy = 0; iy < ny; iy++) {
		float y = float(iy) - shiftyc;
		float ycang = y*cang + shiftyc;
		float ysang = -y*sang + shiftxc;
		for (int ix = 0; ix < nx; ix++) {
			(*ret)(ix,iy) = background;
			float x = float(ix) - shiftxc;
			float xold = x*cang + ysang;
			float yold = x*sang + ycang;
			int iyold = int(yold);
			float q = yold - float(iyold);
			float qcomp = 1.f - q;
			int ixold = int(xold);
			// Note: nx-2 or ny-2 below because need room for
			// (forward) interpolation
			if ((yold>=0 && iyold<=(ny-2)) && (xold>=0 && ixold<=(nx-2))) {
				// inside boundaries of input image
				float p = xold - ixold;
				float pcomp = 1.f - p;
				(*ret)(ix,iy) = q*(pcomp*(*this)(ixold,iyold+1)
						         + p*(*this)(ixold+1,iyold+1))
					        + qcomp*(pcomp*(*this)(ixold,iyold)
									 + p*(*this)(ixold+1,iyold));
			}
		}
	}
	ret->done_data();
	ret->update();
	return ret;
}

EMData*
EMData::rot_scale_trans2D(float angDeg, float delx,float dely, float scale) {
	float ang=angDeg*M_PI/180.0f;
	if (1 >= ny)
		throw ImageDimensionException("Can't rotate 1D image");
	if (nz==1) { 
		vector<int> saved_offsets = get_array_offsets();
		set_array_offsets(0,0,0);
		if (0.f == scale) scale = 1.f; // silently fix common user error
		EMData* ret = copy_head();
		delx = fmod(delx, float(nx));
		dely = fmod(dely, float(ny));
		// center of image
		int xc = nx/2;
		int yc = ny/2;
		// shifted center for rotation
		float shiftxc = xc + delx;
		float shiftyc = yc + dely;
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
			for (int iy = 0; iy < ny; iy++) {
				float y = float(iy) - shiftyc;
			#ifdef _WIN32
				if (y < ymin) y = _MIN(y+ny,ymax);
				if (y > ymax) y = _MAX(y-ny,ymin);
			#else
				if (y < ymin) y = std::min(y+ny,ymax);
				if (y > ymax) y = std::max(y-ny,ymin);
			#endif	//_WIN32
				float ycang = y*cang/scale + yc;
				float ysang = -y*sang/scale + xc;
				for (int ix = 0; ix < nx; ix++) {
					float x = float(ix) - shiftxc;
				#ifdef _WIN32
					if (x < xmin) x = _MIN(x+nx,xmax);
					if (x > xmax) x = _MAX(x-nx,xmin);
				#else
					if (x < xmin) x = std::min(x+nx,xmax);
					if (x > xmax) x = std::max(x-nx,xmin);
				#endif	//_WIN32
					float xold = x*cang/scale + ysang ;
					float yold = x*sang/scale + ycang ;
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
	
	Vec3f  translations = RA.get_posttrans();
	float delx = translations.at(0);
	float dely = translations.at(1);
	float delz = translations.at(2);
	EMData* ret = copy_head();
	delx = fmod(delx, float(nx));
	dely = fmod(dely, float(ny));
	delz = fmod(delz, float(nz));
	vector<int> saved_offsets = get_array_offsets();
	set_array_offsets(0,0,0);
	int xc = nx/2;
	int yc = ny/2;
	int zc = nz/2;
//         shifted center for rotation
	float shiftxc = xc + delx;
	float shiftyc = yc + dely;
	float shiftzc = zc + delz;
// 	bounds if origin at center
	float xmin = -nx/2.0f;
	float ymin = -ny/2.0f;
	float zmin = -nz/2.0f;
	float xmax = -xmin;
	float ymax = -ymin;
	float zmax = -zmin;
	if (0 == nx%2) xmax--;
	if (0 == ny%2) ymax--;
	if (0 == nz%2) zmax--;

	Transform3D RAinv; // = new Transform3D();
	RAinv= RA.inverse();

	if (1 >= ny)
		throw ImageDimensionException("Can't rotate 1D image");
	if (nz==1) { 
		for (int iy = 0; iy < ny; iy++) {
			float y = float(iy) - shiftyc;
		#ifdef _WIN32
			if (y < ymin) y = _MIN(y+ny,ymax);
			if (y > ymax) y = _MAX(y-ny,ymin);
		#else
			if (y < ymin) y = std::min(y+ny,ymax);
			if (y > ymax) y = std::max(y-ny,ymin);
		#endif	//_WIN32
			float xoldy = y*RAinv[0][1]+xc; 
			float yoldy = y*RAinv[1][1]+yc;
			for (int ix = 0; ix < nx; ix++) {
				float x = float(ix) - shiftxc;
			#ifdef _WIN32
				if (x < xmin) x = _MIN(x+nx,xmax);
				if (x > xmax) x = _MAX(x-nx,xmin);
			#else
				if (x < xmin) x = std::min(x+nx,xmax);
				if (x > xmax) x = std::max(x-nx,xmin);
			#endif	//_WIN32
				float xold = xoldy + x*RAinv[0][0]; 
				float yold = yoldy + x*RAinv[1][0];
//                    		printf("\t\t xold = %f, yold=%f \n",xold,yold);
				(*ret)(ix,iy) = Util::quadri(xold+1.0f, yold+1.0f, nx, ny, get_data());
				   //have to add one as quadri uses Fortran counting
			} //ends x loop
		} // ends y loop
		set_array_offsets(saved_offsets);
		return ret;
	} else {
//		 This begins the 3D version
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
					
		for (int iz = 0; iz < nz; iz++) {for (int iy = 0; iy < ny; iy++) {for (int ix = 0; ix < nx; ix++) {
		      (*ret)(ix,iy,iz) = 0;}}}   // initialize returned data
		
		for (int iz = 0; iz < nz; iz++) {
			float z = float(iz) - shiftzc;
		#ifdef _WIN32
			if (z < zmin) z = _MIN(z+nz,zmax);
			if (z > zmax) z = _MAX(z-nz,zmin);
		#else
			if (z < zmin) z = std::min(z+nz,ymax);
			if (z > zmax) z = std::max(z-nz,ymin);
		#endif	//_WIN32
			float xoldz = z*RAinv[0][2]+xc;
			float yoldz = z*RAinv[1][2]+yc;
			float zoldz = z*RAinv[2][2]+zc;
			for (int iy = 0; iy < ny; iy++) {
				float y = float(iy) - shiftyc;
			#ifdef _WIN32
				if (y < ymin) y = _MIN(y+ny,ymax);
				if (y > ymax) y = _MAX(y-ny,ymin);
			#else
				if (y < ymin) y = std::min(y+ny,ymax);
				if (y > ymax) y = std::max(y-ny,ymin);
			#endif	//_WIN32
				float xoldzy = xoldz + y*RAinv[0][1] ;
				float yoldzy = yoldz + y*RAinv[1][1] ;
				float zoldzy = zoldz + y*RAinv[2][1] ;
				for (int ix = 0; ix < nx; ix++) {
					float x = float(ix) - shiftxc;
				#ifdef _WIN32
					if (x < xmin) x = _MIN(x+nx,xmax);
					if (x > xmax) x = _MAX(x-nx,xmin);
				#else
					if (x < xmin) x = std::min(x+nx,xmax);
					if (x > xmax) x = std::max(x-nx,xmin);
				#endif	//_WIN32
					float xold = xoldzy + x*RAinv[0][0] ;
					float yold = yoldzy + x*RAinv[1][0] ;
					float zold = zoldzy + x*RAinv[2][0] ;

//         This is currently coded the way  SPIDER coded it,
//            changing floor to round  in the next 3 lines below may be better
					int IOX = (int) floor(xold); // This is the center of the array
					int IOY = (int) floor(yold ); // In the next loop we interpolate
					int IOZ = (int) floor(zold ); //  If floor is used dx is positive

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
	//ret->to_zero();  //we will leave margins zeroed.
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
