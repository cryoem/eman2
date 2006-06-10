/**
 * $Id$
 */
#ifndef eman_native_fft_h__
#define eman_native_fft_h__

#ifdef NATIVE_FFT

namespace EMAN
{
	class Nativefft
	{
	public:
		/**
		*  D REAL MIXED RADIX FFT.
		*  INPUT:  X(N) - REAL ARRAY
		*  OUTPUT: N EVEN  X(N+2)
		*  ORDER OF ELEMENTS:
		*  R(0),0.0, R(1), I(1), R(2), I(2), ....., R(N/2-1), I(N/2-1), R(N/2),0.0
		*
		*         N ODD  X(N+1)
		*  R(0),0.0, R(1), I(1), R(2), I(2), ....., R(N/2-1), I(N/2-1), R(N/2),I(N/2)
		* 
		*  HERE WE FOLLOW THE CONVENTION THAT INTEGER DIVISION 
		*  IS ROUNDED DOWN, E.G. 5/2 =2)
		*
		*  INV: +1 FORWARD FFT
		*       -1 INVERSE FFT
		*  ON OUTPUT INV=0 INDICATES ERROR
		*/
		
		//----------2D FFT INPLACE---------------------------
		// 1d inplace FFT
		static int fmrs_1rf(float *x, float *work, int nsam);
		
		// 1d inplace IFFT
		static int fmrs_1rb(float *x, float *work, int nsam);
		
		//----------2D FFT INPLACE---------------------------
		// 2d inplace fft
		static int fmrs_2rf(float *x, float *work, int lda, int nsam, int nrow);
		
		// 2d inplace IFFT
		static int fmrs_2rb(float *y, float *work, int lda, int nsam, int nrow);
		
		//---------------3D INPLACE FFT----------------------
		// 3d inplace FFT
		static int fmrs_3rf(float *b, float *work, int lda, int nsam, int nrow, int nslice);
		
		// 3d inplace IFFT
		static int fmrs_3rb(float *b, float *work, int lda, int nsam, int nrow, int nslice);
		
		//  2D FFT OUT OF PLACE
		static int Nativefft::ftp_2rf(float *reald, float *complexd, int lda, int nsam, int nrow);
		
		//  2D IFT This one is part of out of place, but overwrites the input by the real image (not resized!)....
		static int Nativefft::ftp_2rb(float *complexd, int lda, int nsam, int nrow);

		//  3D FFT OUT OF PLACE
		static int Nativefft::ftp_3rf(float *reald, float *complexd, int lda, int nsam, int nrow, int nslice);

		//  3D IFT This one is part of out of place, but overwrites the input by the real image (not resized!)....
		static int Nativefft::ftp_3rb(float *complexd, int lda, int nsam, int nrow, int nslice);

		
	private:
		/* chao modify integer to int */
		typedef int integer;
		
		static int fftmcf_(float *a, float *b, integer *ntot, integer *n, integer *nspan, integer *isn);
	};



//#define __cplusplus
//#ifdef __cplusplus
//extern "C" {
//#endif

//#define DOUBLE

#ifdef DOUBLE
#define Treal double
#else
#define Treal float
#endif

 void cfftf(int N, Treal data[], Treal wrk[]);
 void cfftb(int N, Treal data[], Treal wrk[]);
 void cffti(int N, Treal wrk[]);

 void rfftf(int N, Treal data[], Treal wrk[]);
 void rfftb(int N, Treal data[], Treal wrk[]);
 void rffti(int N, Treal wrk[]);

//#ifdef __cplusplus
//}
#//endif



}

#endif	//NATIVE_FFT

#endif	//eman_native_fft_h__
