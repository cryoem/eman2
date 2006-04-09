/**
 * $Id$
 */
#include "emdata.h"
 
#include "gsl_sf_bessel.h"
#include "gsl_errno.h"
#include <vector>
using std::vector;
using std::cout;
using namespace EMAN;

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

			kValue =sqrt(fjkx*fjkx +  fjky*fjky )  ;
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
			float r = sqrt(float(k*k) + float(j*j) + float(i*i));
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
vector<float> EMData::cog()
{
	
	vector<float> cntog;
	int ndim = get_ndim();
	int i=1,j=1,k=1;
	float val,sum1=0.f,MX=0.f,RX=0.f,MY=0.f,RY=0.f,MZ=0.f,RZ=0.f;
	
	if (ndim == 1)
	{
			for ( i = 1;i <= nx; i++)
			{
				val=rdata(i,j,k);
				sum1 += val;
				MX   += (i*val);
				RX   += ((i^2)*val);
			}
			cntog.push_back(MX);
			cntog.push_back(RX);
	
	}	
	else if (ndim == 2)
	{	
			for (j=1;j<=ny;j++)
				{
					for (i=1;i<=nx;i++)
					{
						val = rdata(i,j,k);
						sum1 += val;
						MX += (i*val);
						MY += (j*val);
						RX += ((i^2)*val);
						RY += ((j^2)*val);
					}
				}
			cntog.push_back(MX);
			cntog.push_back(RX);
			cntog.push_back(MY);
			cntog.push_back(RY);
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
						MX += (i*val);
						MY += (j*val);
						MZ += (k*val);
						RX += ((i^2)*val);
						RY += ((j^2)*val);
						RZ += ((k^2)*val);
					}
				}
			}
			cntog.push_back(MX);
			cntog.push_back(RX);
			cntog.push_back(MY);
			cntog.push_back(RY);
			cntog.push_back(MZ);
			cntog.push_back(RZ);
	}	 
	return cntog;
}
#undef rdata
