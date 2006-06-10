/**
 * $Id$
 */
#include <string>
#include "emfft.h"
#include "log.h"


#include <iostream>
using std::cout;
using std::endl;

#include "util.h"


#ifdef DJBFFT
extern "C" {
	#include <fftr4.h>
}
#endif	//DJBFFT

using namespace EMAN;

namespace {
	int get_rank(int ny, int nz)
	{
		int rank = 3;
		if (ny == 1) {
			rank = 1;
		}
		else if (nz == 1) {
			rank = 2;
		}
		return rank;
	}
}

#ifdef FFTW2

int EMfft::real_to_complex_1d(float *real_data, float *complex_data, int n)
{
#ifdef DJBFFT
	switch(n)
	{
		if ( n==2 || n==4 || n==8 || n==16 || n==32 || n==64 || n==128
			|| n==256 || n==512 || n==1024 || n==2048 || n==4096 || n==8192 )
		{
			memcpy( complex_data, real_data, n * sizeof(float) );		
		}
		
		case 2:
			fftr4_2( (real4 *)complex_data );
			
		case 4: 
			fftr4_4( (real4 *)complex_data );
			
		case 8:
			fftr4_8( (real4 *)complex_data );
			 
		case 16:
			fftr4_16( (real4 *)complex_data ); 
			
		case 32:
			fftr4_32( (real4 *)complex_data ); 
			
		case 64:
			fftr4_64( (real4 *)complex_data ); 
			
		case 128: 
			fftr4_128( (real4 *)complex_data );
			
		case 256:
			fftr4_256( (real4 *)complex_data ); 
			
		case 512: 
			fftr4_512( (real4 *)complex_data );
			
		case 1024: 
			fftr4_1024( (real4 *)complex_data );
			
		case 2048: 
			fftr4_2048( (real4 *)complex_data );
			
		case 4096: 
			fftr4_4096( (real4 *)complex_data );
			
		case 8192:
			fftr4_8192( (real4 *)complex_data );
			
		default:
			const int complex_n = n + 2 - n%2;
			float * fft_data = new float[n];
			rfftw_plan p = rfftw_create_plan(n, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
			rfftw_one(p, (fftw_real *) real_data, (fftw_real *) complex_data);
			rfftw_destroy_plan(p);
			for(int i=0; i<complex_n; ++i) {
				if(i%2==0) {	//copy real part of complex array
					complex_data[i] = fft_data[i/2];
				}
				else {	//copy imaginary part of complex array
					if(i==1) {
						complex_data[i] = 0.0f;
					}
					else {
						if(n%2 == 0 && i == complex_n-1 ) {
							complex_data[i] = 0.0f;
						}
						else {
							complex_data[i] = fft_data[n-i/2];
						}
					}
				}
			}
	
			delete [] fft_data;
	}	
#else
	const int complex_n = n + 2 - n%2;
	float * fft_data = new float[n];
	rfftw_plan p = rfftw_create_plan(n, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);	
	rfftw_one(p, (fftw_real *) real_data, (fftw_real *) fft_data);
	rfftw_destroy_plan(p);
	
	for(int i=0; i<complex_n; ++i) {
		if(i%2==0) {	//copy real part of complex array
			complex_data[i] = fft_data[i/2];
		}
		else {	//copy imaginary part of complex array
			if(i==1) {
				complex_data[i] = 0.0f;
			}
			else {
				if(n%2 == 0 && i == complex_n-1 ) {
					complex_data[i] = 0.0f;
				}
				else {
					complex_data[i] = fft_data[n-i/2];
				}
			}
		}
	}
	
	delete [] fft_data;
#endif	//DJBFFT	

	return 0;
}

int EMfft::complex_to_real_1d(float *complex_data, float *real_data, int n)
{
#ifdef DJBFFT
	switch(n)
	{
		if ( n==2 || n==4 || n==8 || n==16 || n==32 || n==64 || n==128
			|| n==256 || n==512 || n==1024 || n==2048 || n==4096 || n==8192 )
		{
			memcpy( real_data, complex_data, n * sizeof(float) );		
		}
		
		case 2:	
			fftr4_scale2( (real4 *)real_data );
     		fftr4_un2( (real4 *)real_data );
			
		case 4: 
			fftr4_scale4( (real4 *)real_data );
     		fftr4_un4( (real4 *)real_data );
			
		case 8:
			fftr4_scale8( (real4 *)real_data );
     		fftr4_un8( (real4 *)real_data );
			 
		case 16:
			fftr4_scale16( (real4 *)real_data );
     		fftr4_un16( (real4 *)real_data );
			
		case 32:
			fftr4_scale32( (real4 *)real_data );
     		fftr4_un32( (real4 *)real_data );
			
		case 64:
			fftr4_scale64( (real4 *)real_data );
     		fftr4_un64( (real4 *)real_data );
			
		case 128: 
			fftr4_scale128( (real4 *)real_data );
     		fftr4_un128( (real4 *)real_data );
			
		case 256:
			fftr4_scale256( (real4 *)real_data );
     		fftr4_un256( (real4 *)real_data );
			
		case 512: 
			fftr4_scale512( (real4 *)real_data );
     		fftr4_un512( (real4 *)real_data );
			
		case 1024: 
			fftr4_scale1024( (real4 *)real_data );
     		fftr4_un1024( (real4 *)real_data );
			
		case 2048: 
			fftr4_scale2048( (real4 *)real_data );
     		fftr4_un2048( (real4 *)real_data );
			
		case 4096: 
			fftr4_scale4096( (real4 *)real_data );
     		fftr4_un4096( (real4 *)real_data );
			
		case 8192:
			fftr4_scale8192( (real4 *)real_data );
     		fftr4_un8192( (real4 *)real_data );
			
		default:
			const int complex_n = n + 2 - n%2;
			float * fft_data = new float[n];
			
			for(int i=0; i<complex_n; ++i) {
				if(i%2 == 0) {	//copy real part of complex array
					fft_data[i/2] = complex_data[i];
				}
				else {	//copy imaginary part of complex array
					if(i==1) {continue;}
					if(!(n%2 == 0 && i == complex_n-1)) {
						fft_data[n-i/2] = complex_data[i];
					}
				}
			}
			
			rfftw_plan p = rfftw_create_plan(n, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
			rfftw_one(p, (fftw_real *) fft_data, (fftw_real *) real_data);
			rfftw_destroy_plan(p);
	}
#else
	const int complex_n = n + 2 - n%2;
	float * fft_data = new float[n];
	
	for(int i=0; i<complex_n; ++i) {
		if(i%2 == 0) {	//copy real part of complex array
			fft_data[i/2] = complex_data[i];
		}
		else {	//copy imaginary part of complex array
			if(i==1) {continue;}
			if(!(n%2 == 0 && i == complex_n-1)) {
				fft_data[n-i/2] = complex_data[i];
			}
		}
	}
	rfftw_plan p = rfftw_create_plan(n, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
	rfftw_one(p, (fftw_real *) fft_data, (fftw_real *) real_data);	
	rfftw_destroy_plan(p);
	
	delete [] fft_data;
#endif	//DJBFFT
	
	return 0;
	
}

int EMfft::real_to_complex_nd(float *real_data, float *complex_data, int nx, int ny, int nz)
{
//	const int complex_nx = nx + 2 - nx%2;
	const int rank = get_rank(ny, nz);
	int dims[3];
	dims[0] = nz;
	dims[1] = ny;
	dims[2] = nx;
	
	switch(rank) {
		case 1:
			real_to_complex_1d(real_data, complex_data, nx);
			break;
		
		case 2: 
		case 3:
		{
			rfftwnd_plan plan_nd;
			if(real_data == complex_data) {
				plan_nd = rfftwnd_create_plan(rank, dims + (3 - rank), 
						FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
			}
			else {
				plan_nd = rfftwnd_create_plan(rank, dims + (3 - rank), 
						FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
			}
			rfftwnd_one_real_to_complex(plan_nd, (fftw_real *) real_data,
										(fftw_complex *) complex_data);
			rfftwnd_destroy_plan(plan_nd);
		}
			break;
		
		default:
			LOGERR("Should NEVER be here!!!");
			break;
	}
	
	return 0;
}

int EMfft::complex_to_real_nd(float *complex_data, float *real_data, int nx, int ny, int nz)
{
//	const int complex_nx = nx + 2 - nx%2;
	const int rank = get_rank(ny, nz);
	int dims[3];
	dims[0] = nz;
	dims[1] = ny;
	dims[2] = nx;
	
	switch(rank) {
		case 1:
			complex_to_real_1d(complex_data, real_data, nx);
			break;
		
		case 2:
		case 3:
		{
			rfftwnd_plan plan_nd;
			if(real_data == complex_data) {
				plan_nd = rfftwnd_create_plan(rank, dims + (3 - rank), 
						FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);
			}
			else {
				plan_nd = rfftwnd_create_plan(rank, dims + (3 - rank), 
						FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
			}	
			rfftwnd_one_complex_to_real(plan_nd, (fftw_complex *) complex_data,	
				(fftw_real *) real_data);
			rfftwnd_destroy_plan(plan_nd);
		}		
			break;
		
		default:
			LOGERR("Should NEVER be here!!!");
			break;
	}			

	return 0;
}
#endif	//FFTW2


#ifdef FFTW3

int EMfft::real_to_complex_1d(float *real_data, float *complex_data, int n)
{//cout<<"doing fftw3"<<endl;
	fftwf_plan plan = fftwf_plan_dft_r2c_1d(n, real_data, (fftwf_complex *) complex_data,
											FFTW_ESTIMATE);
	fftwf_execute(plan);
	fftwf_destroy_plan(plan);
	return 0;
}

int EMfft::complex_to_real_1d(float *complex_data, float *real_data, int n)
{
	fftwf_plan plan = fftwf_plan_dft_c2r_1d(n, (fftwf_complex *) complex_data, real_data,
											FFTW_ESTIMATE);
	fftwf_execute(plan);
	fftwf_destroy_plan(plan);
	return 0;
}

int EMfft::real_to_complex_nd(float *real_data, float *complex_data, int nx, int ny, int nz)
{
	const int rank = get_rank(ny, nz);
	int dims[3];
	dims[0] = nz;
	dims[1] = ny;
	dims[2] = nx;
	
	switch(rank) {
		case 1:
			real_to_complex_1d(real_data, complex_data, nx);
			break;
		
		case 2:
		case 3:
		{
			fftwf_plan plan = fftwf_plan_dft_r2c(rank, dims + (3 - rank), 
					real_data, (fftwf_complex *) complex_data, FFTW_ESTIMATE);
			fftwf_execute(plan);
			fftwf_destroy_plan(plan);
		}
			break;
		
		default:
			LOGERR("Should NEVER be here!!!");
			break;
	}
	
	return 0;
}

int EMfft::complex_to_real_nd(float *complex_data, float *real_data, int nx, int ny, int nz)
{
	const int rank = get_rank(ny, nz);
	int dims[3];
	dims[0] = nz;
	dims[1] = ny;
	dims[2] = nx;
	
	switch(rank) {
		case 1:
			complex_to_real_1d(complex_data, real_data, nx);
			break;
		
		case 2:
		case 3:
		{
			fftwf_plan plan = fftwf_plan_dft_c2r(rank, dims + (3 - rank), 
					(fftwf_complex *) complex_data, real_data, FFTW_ESTIMATE);
			fftwf_execute(plan);
			fftwf_destroy_plan(plan);
		}
			break;
			
		default:
			LOGERR("Should NEVER be here!!!");
			break;
	}
	
	return 0;
}

#endif	//FFTW3

#ifdef NATIVE_FFT
#include "sparx/native_fft.h"
int EMfft::real_to_complex_1d(float *real_data, float *complex_data, int n)
{
	//int complex_size = n + 2 - n%2;
	float * work = (float*) malloc((2*n+15)*sizeof(float));
	if (!work) {
		fprintf(stderr,"real_to_complex_1d: failed to allocate work\n");
		LOGERR("real_to_complex_1d: failed to allocate work\n");	
	}
	
	rffti(n, work);
	memcpy(&complex_data[1], real_data, n * sizeof(float));
	rfftf(n, &complex_data[1], work);
	complex_data[0] = complex_data[1] ;
	complex_data[1] = 0.0f ;
	if (n%2 == 0)  complex_data[n+1] = 0.0f ;

	free(work);
	return 0;
}
/*{
	//int complex_size = n + 2 - n%2;
	
	//memcpy(complex_data, real_data, n * sizeof(float));
	//float * work = (float*) malloc((2*n+15)*sizeof(float));
	//if (!work) {
	//	fprintf(stderr,"real_to_complex_1d: failed to allocate work\n");
	//	LOGERR("real_to_complex_1d: failed to allocate work\n");	
	//}
	//Nativefft::fmrs_1rf(complex_data, work, n);
	//cout<<"doing rfftf"<<endl;
	//rffti(n, work);
	
	//rfftf(n, real_data, work);
	
	//cout<<"doing fftr_q"<<endl;
	int l=(int)(log2(n));
        Util::fftr_q(real_data,l);

	//free(work);
	return 0;
}//
{
	int complex_size = n + 2 - n%2;
	
	memcpy(complex_data, real_data, n * sizeof(float));
	float * work = (float*) malloc(complex_size*sizeof(float));
	if (!work) {
		fprintf(stderr,"real_to_complex_1d: failed to allocate work\n");
		LOGERR("real_to_complex_1d: failed to allocate work\n");	
	}
	
	Nativefft::fmrs_1rf(complex_data, work, n);
	
	free(work);
	return 0;
}*/

int EMfft::complex_to_real_1d(float *complex_data, float *real_data, int n)
{
	//here, n is the "logical" size of DFT, not the array size
	float * work = (float*) malloc((2*n+15)*sizeof(float));
	if (!work) {
		fprintf(stderr,"real_to_complex_1d: failed to allocate work\n");
		LOGERR("complex_to_real_1d: failed to allocate work\n");	
	}
	
	rffti(n, work);
	
	memcpy(&real_data[1], &complex_data[2], (n-1) * sizeof(float));
	real_data[0] = complex_data[0] ;
	rfftb(n, real_data, work);
	//  Normalize
	float nrm = 1.0f/float(n);
	for (int i = 0; i<n; i++) real_data[i] *= nrm;
	free(work);
	return 0;
}
/*{
	int complex_size = n + 2 - n%2;
	
	//here, n is the "logical" size of DFT, not the array size
	memcpy(real_data, complex_data, complex_size * sizeof(float));
	float * work = (float*)malloc(complex_size*sizeof(float));
	if (!work) {
		fprintf(stderr,"real_to_complex_1d: failed to allocate work\n");
		LOGERR("complex_to_real_1d: failed to allocate work\n");	
	}
	
	Nativefft::fmrs_1rb(real_data, work, n);
	
	free(work);
	return 0;
}*/

int EMfft::real_to_complex_nd(float *real_data, float *complex_data, int nx, int ny, int nz)
{
	const int rank = get_rank(ny, nz);
	const int complex_nx = nx + 2 - nx%2;
	
	switch(rank) {
		case 1:		//for 1D fft
			real_to_complex_1d(real_data, complex_data, nx);
			return 0;
		case 2:		//for 2D fft
		{
			/*if(real_data != complex_data) {
				for (int j = 0; j < ny; j++) {
					memcpy(&complex_data[complex_nx*j], &real_data[nx*j], nx*sizeof(float));	
				}
			}
			float * work = (float*) malloc(complex_nx*sizeof(float));
	   		if (!work) {
	   			fprintf(stderr,"real_to_complex_nd(2df): failed to allocate work\n");
	   			LOGERR("real_to_complex_nd(2df): failed to allocate work\n");
	   		}
	   		
	   		// 2d inplace fft, overwrite y
	   		int status = Nativefft::fmrs_2rf(complex_data, work, complex_nx, nx, ny);
	   		if (status !=0) {
	      		fprintf(stderr, "real_to_complex_nd(2df): status = %d\n", status);
	      		LOGWARN("real_to_complex_nd(2df): status = %d\n", status);
	   		}
	   	 
	   		free(work);*/
			
	   		int status = Nativefft::ftp_2rf(real_data, complex_data, complex_nx, nx, ny);
	   		if (status !=0) {
	      		fprintf(stderr, "real_to_complex_nd(2df): status = %d\n", status);
	      		LOGWARN("real_to_complex_nd(2df): status = %d\n", status);
	   		}
	   		return 0;
		}
	   	case 3:		//for 3D fft
	   	{
			/*if(real_data != complex_data) {
				for (int k = 0; k<nz; k++) {
		      		for (int j = 0; j < ny; j++) {
		         		memcpy(&complex_data[complex_nx*ny*k+j*complex_nx], &real_data[nx*ny*k+j*nx], nx*sizeof(float));
		      		}
				}
			}
			float * work = (float*) malloc(1*sizeof(float));//malloc(complex_nx*sizeof(float));
	   		if (!work) {
	   			fprintf(stderr,"real_to_complex_nd(3df): failed to allocate work\n");
	   			LOGERR("real_to_complex_nd(3df): failed to allocate work\n");
	   		}*/
	   		
	   		// 3d inplace fft, overwrite complex_data
	   		int status = Nativefft::ftp_3rf(real_data, complex_data, complex_nx, nx, ny, nz);
	   		if (status !=0) {
	      		fprintf(stderr, "real_to_complex_nd(3df): status = %d\n", status);
	      		LOGWARN("real_to_complex_nd(3df): status = %d\n", status);
	   		}
	   		
	   		//free(work);
	   		return 0;
	   	}
	   	default:
	   		LOGERR("Never should be here...\n");
	   		return -1;	
	}
}

int EMfft::complex_to_real_nd(float *complex_data, float *real_data, int nx, int ny, int nz)
{
	const int rank = get_rank(ny, nz);
	const int complex_nx = nx + 2 - nx%2;
	
	switch(rank) {
		case 1:		//for 1D ift
			complex_to_real_1d(complex_data, real_data, nx);
			return 0;
		case 2:		//for 2D ift
		/*{	
			if(real_data != complex_data) {
				memcpy(real_data, complex_data, complex_nx*ny*sizeof(float));
			}
			
			float * work = (float*) malloc(complex_nx*sizeof(float));
			if (!work) {
				fprintf(stderr,"complex_to_real_nd(2db): failed to allocate work\n");
				LOGERR("complex_to_real_nd(2db): failed to allocate work\n");
			}
			
			// 2d inplace ift, overwrite real_data
	   		int status = Nativefft::fmrs_2rb(real_data, work, complex_nx, nx, ny);
	   		if (status !=0) {
	      		fprintf(stderr, "complex_to_real_nd(2db): status = %d\n", status);
	      		LOGWARN("complex_to_real_nd(2db): status = %d\n", status);
			}
			
			free(work);
			return 0;
		}*/
		{	//  Only out of place!
			memcpy(real_data, complex_data, complex_nx*ny*sizeof(float));
			
			// 2d inplace ift, overwrite real_data
	   		int status = Nativefft::ftp_2rb(real_data, complex_nx, nx, ny);

	   		if (status !=0) {
	      		fprintf(stderr, "complex_to_real_nd(2db): status = %d\n", status);
	      		LOGWARN("complex_to_real_nd(2db): status = %d\n", status);
			}
			return 0;
		}
		case 3:		//for 3D ift
		{	// Only out of place!
			memcpy(real_data, complex_data, complex_nx*ny*nz*sizeof(float));
	   		
	   		// 3d inplace fft, overwrite real_data
	   		int status = Nativefft::ftp_3rb(real_data, complex_nx, nx, ny, nz);
	   		if (status !=0) {
	      		fprintf(stderr, "complex_to_real_nd(3db): status = %d\n", status);
	      		LOGWARN("complex_to_real_nd(3db): status = %d\n", status);
	   		}
	   		return 0;
		}
	   	default:
	   		LOGERR("Never should be here...\n");
	   		return -1;
	}
}
#endif	//NATIVE_FFT

#ifdef 	ACML
#include <iostream>
using std::cout;
using std::endl;

int EMfft::real_to_complex_1d(float *real_data, float *complex_data, int n)
{
	const int complex_n = n + 2 - n%2;	//the size for 1D complex array
	int info;
	
	/* Allocate communication work array */
	float * comm = (float *)malloc((3*n+100)*sizeof(float));
	/* Allocate work array to store ACML complex array*/
	float * fft_data = (float *)malloc(n * sizeof(float));
	
	//copy real_data to complex_data then apply inplace FFT on complex data
	memcpy(fft_data, real_data, n * sizeof(float));
	
	/* Initialize communication work array */
	scfft(0, n, fft_data, comm, &info);
	if(info != 0) {
		LOGERR("Error happening in Initialize communication work array: %d", info);
	}
	
	/* Compute a real --> Hermitian transform */
	scfft(1, n, fft_data, comm, &info);
	if(info != 0) {
		LOGERR("Error happening in Initialize communication work array: %d", info);
	}
	
	/**ACML fft real to complex 1D result store as:
	 * let X be an array of length N and with first index 0,
	 * - X(i) contains the real part of Z(i) for i = 0, ..., N/2
	 * - X(N-i) contains the imaginary part of Z(i) for i=1, ..., (N-1)/2 
	 * so we need re-organize the data layout and time all data by sqrt(N)
	 * to make the reault consistent with FFTW  */
	transform(fft_data, fft_data+n, fft_data, time_sqrt_n(n));
	
	for(int i=0; i<complex_n; ++i) {
		if(i%2==0) {	//copy real part of complex array
			complex_data[i] = fft_data[i/2];
		}
		else {	//copy imaginary part of complex array
			if(i==1) {
				complex_data[i] = 0.0f;
			}
			else {
				if(n%2 == 0 && i == complex_n-1 ) {
					complex_data[i] = 0.0f;
				}
				else {
					complex_data[i] = fft_data[n-i/2];
				}
			}
		}
	}
	
	free(fft_data);
	free(comm);
	return 0;
}

int EMfft::complex_to_real_1d(float *complex_data, float *real_data, int n)
{
	int complex_n = n + 2 - n%2;	//the size for 1D complex array
	int info;
	
	/* Allocate communication work array */
	float * comm = (float *)malloc((3*n+100)*sizeof(float));
		
	for(int i=0; i<complex_n; ++i) {
		if(i%2 == 0) {	//copy real part of complex array
			real_data[i/2] = complex_data[i];
		}
		else {	//copy imaginary part of complex array
			if(i==1) {continue;}
			if(!(n%2 == 0 && i == complex_n-1)) {
				real_data[n-i/2] = complex_data[i];
			}
		}
	}
	transform(real_data, real_data+n, real_data, divide_sqrt_n(n));
	
	/* Initialize communication work array */
	csfft(0, n, real_data, comm, &info);
	if(info != 0) {
		LOGERR("Error happening in Initialize communication work array: %d", info);
	}
	
	/* Conjugate the Vector X to simulate inverse transform */
	for (int j = n/2+1; j < n; j++) {
    	real_data[j] = -real_data[j];
	}
	
	/* Compute a Hermitian --> real transform */
	csfft(1, n, real_data, comm, &info);
	if(info != 0) {
		LOGERR("Error happening in Initialize communication work array: %d", info);
	}
	
	free(comm);
	return 0;
}

int EMfft::real_to_complex_2d(float *real_data, float *complex_data, int nx, int ny)
{
	const int complex_nx = nx + 2 - nx%2;
	int info;
	/* Allocate communication work array */
	float * comm = (float *)malloc((3*nx+100)*sizeof(float));
	
	/* Allocate work array to store ACML complex array*/
	float * fft_data = (float *)malloc(complex_nx * ny * sizeof(float));
	cout << "fft_data after allocation:" << endl;
	for(int j=0; j<ny; ++j) {
		for(int i=0; i<complex_nx; ++i) {
			cout << "data[" << i << "][" << j << "] = " << fft_data[i+j*complex_nx] << "\t";
		}
		cout << endl;
	}
	
	/*copy real_data to complex_data then apply inplace FFT on complex data*/
	for(int i=0; i<ny; ++i) {
		memcpy(fft_data+i*complex_nx, real_data+i*nx, nx*sizeof(float));
	}
	//memcpy(fft_data, real_data, nx * ny * sizeof(float));
	cout << endl << "real_data array: " << endl;
	for(int j=0; j<ny; ++j) {
		for(int i=0; i<nx; ++i) {
			cout << real_data[i+j*nx] << "\t";
		}
		cout << endl;
	}
	cout << endl << "the fft_data array: " << endl;
	for(int j=0; j<ny; ++j) {
		for(int i=0; i<complex_nx; ++i) {
			cout << fft_data[i+j*complex_nx] << "\t";
		}
		cout << endl;
	}
	
	//do a multiple 1d real-to-complex transform on x direction
	scfftm(ny, nx, fft_data, comm, &info);
	
	cout << endl << "the fft_data array after x dim transform: " << endl;
	for(int j=0; j<ny; ++j) {
		for(int i=0; i<complex_nx; ++i) {
			cout << fft_data[i+j*complex_nx] << "\t";
		}
		cout << endl;
	}
	
/*	cout << "original fft array" << endl;
/	cout << "complex_nx = " << complex_nx << " ny = " << n*(fft_data2+i+j*2*ny+1) << endl;
	for(int i=0; i<ny; ++i) {
		for(int j=0; j<complex_nx; ++j) {
			cout << *(fft_data+j+i*complex_nx) << "\t";
		}
		cout << endl;
	}
*/	
	//do a multiple 1d complex to complex transformation on y direction
	float * fft_data2 = (float *)malloc(complex_nx * ny * sizeof(float));
/*	cout << "fft array rearranged in y dimension" << endl;
	for(int i=0; i<complex_nx; ++i) {
		for(int j=0; j<ny; ++j) {
			*(fft_data2+i+j*2*ny) = *(fft_data+i+j*complex_nx);
			*(fft_data2+i+j*2*ny+1) = *(fft_data+i+complex_nx/2+j*complex_nx);
			cout << *(fft_data2+i+j*2*ny) << "\t" << *(fft_data2+i+j*2*ny+1) << "\t";
		}
		cout << endl;
	}
*/	
	if(info != 0) {
		LOGERR("Error happening in scfftm: %d", info);
	}
	
	return 0;	
}

int EMfft::complex_to_real_2d(float *complex_data, float *real_data, int nx, int ny)
{
	return 0;
}

int EMfft::real_to_complex_3d(float *real_data, float *complex_data, int nx, int ny, int nz)
{
	return 0;
}

int EMfft::complex_to_real_3d(float *complex_data, float *real_data, int nx, int ny, int nz)
{
	return 0;
}

int EMfft::real_to_complex_nd(float *real_data, float *complex_data, int nx, int ny, int nz)
{
	const int rank = get_rank(ny, nz);
	
	switch(rank) {
		case 1:
			return real_to_complex_1d(real_data, complex_data, nx);
		case 2:
			return real_to_complex_2d(real_data, complex_data, nx, ny);
		case 3:
			return real_to_complex_3d(real_data, complex_data, nx, ny, nz);		
		default:
			LOGERR("Never should be here...\n");
	   		return -1;
	}	
}

int EMfft::complex_to_real_nd(float *complex_data, float *real_data, int nx, int ny, int nz)
{
	const int rank = get_rank(ny, nz);
	
	switch(rank) {
		case 1:
			return complex_to_real_1d(complex_data, real_data, nx);
		case 2:
			return complex_to_real_2d(complex_data, real_data, nx, ny);
		case 3:
			return complex_to_real_3d(complex_data, real_data, nx, ny, nz);
		default:
			LOGERR("Never should be here...\n");
	   		return -1;
	}
}


#endif	//ACML
