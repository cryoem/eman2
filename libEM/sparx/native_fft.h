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
		static int ftp_2rf(float *reald, float *complexd, int lda, int nsam, int nrow);
		
		//  2D IFT This one is part of out of place, but overwrites the input by the real image (not resized!)....
		static int ftp_2rb(float *complexd, int lda, int nsam, int nrow);

		//  3D FFT OUT OF PLACE
		static int ftp_3rf(float *reald, float *complexd, int lda, int nsam, int nrow, int nslice);

		//  3D IFT This one is part of out of place, but overwrites the input by the real image (not resized!)....
		static int ftp_3rb(float *complexd, int lda, int nsam, int nrow, int nslice);

		
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
