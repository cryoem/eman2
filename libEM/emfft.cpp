/**
 * $Id$
 */
#include <string>
#include "emfft.h"
#include "log.h"

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

//list < FftwPlan * >EMfft::fwplans = list < FftwPlan * >();
//FftwPlan EMfft::planf_1d = FftwPlan();
//FftwPlan EMfft::planr_1d = FftwPlan();

//rfftw_plan FftwPlan::plan_1d = 0;
//rfftwnd_plan FftwPlan::plan_nd = 0;
/*
FftwPlan::FftwPlan()
:	rank(0), direction(UNKNOWN_DIRECTION)
{
	memset(dims, 0, NDIMS * sizeof(int));
}


FftwPlan::FftwPlan(int n, FftwDirection dir, int f)
	:rank(1), direction(dir), flags(f)
{
	dims[NDIMS - 1] = n;
	for (int i = 0; i < NDIMS - 1; i++) {
		dims[i] = 1;
	}

	if (dir == REAL_TO_COMPLEX) {
		plan_1d = rfftw_create_plan(n, FFTW_REAL_TO_COMPLEX, flags);
	}
	else if (dir == COMPLEX_TO_REAL) {
		plan_1d = rfftw_create_plan(n, FFTW_COMPLEX_TO_REAL, flags);
	}
}


FftwPlan::FftwPlan(int ndim, int nx, int ny, int nz, FftwDirection dir, int f)
	:rank(ndim), direction(dir), flags(f)
{
	dims[0] = nz;
	dims[1] = ny;
	dims[2] = nx;

	if (dir == REAL_TO_COMPLEX) {
		plan_nd = rfftwnd_create_plan(rank, dims + (NDIMS - rank), FFTW_REAL_TO_COMPLEX, flags);
	}
	else if (dir == COMPLEX_TO_REAL) {
		plan_nd = rfftwnd_create_plan(rank, dims + (NDIMS - rank), FFTW_COMPLEX_TO_REAL, flags);
	}
}



FftwPlan::~FftwPlan()
{
#if 0
	if (plan_nd) {
		rfftwnd_destroy_plan(plan_nd);
		plan_nd = 0;
	}

	if (plan_1d) {
		rfftw_destroy_plan(plan_1d);
		plan_1d = 0;
	}
#endif
}


bool EMAN::operator==(const FftwPlan & p1, const FftwPlan & p2)
{
	if (p1.rank == p2.rank && p1.direction == p2.direction && p1.flags == p2.flags) {
		for (int i = 0; i < FftwPlan::NDIMS; i++) {
			if (p1.dims[i] != p2.dims[i]) {
				return false;
			}
		}
		return true;
	}
	return false;
}
*/

int EMfft::real_to_complex_1d(float *real_data, float *complex_data, int n)
{
/*	if (planf_1d.get_dim(0) != n) {
		planf_1d = FftwPlan(n, REAL_TO_COMPLEX, FFTW_ESTIMATE);
	}
*/
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
			rfftw_one(planf_1d.get_plan_1d(), (fftw_real *) real_data, (fftw_real *) complex_data);
	}	
#else
	rfftw_plan p = rfftw_create_plan(n, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
	
	rfftw_one(p, (fftw_real *) real_data, (fftw_real *) complex_data);
	
	rfftw_destroy_plan(p);
#endif	//DJBFFT	

	return 0;
}

int EMfft::complex_to_real_1d(float *complex_data, float *real_data, int n)
{
/*	if (planr_1d.get_dim(0) != n) {
		planr_1d = FftwPlan(n, COMPLEX_TO_REAL, FFTW_ESTIMATE);
	}
*/
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
			rfftw_one(planr_1d.get_plan_1d(), (fftw_real *) complex_data, (fftw_real *) real_data);
	}
#else
	rfftw_plan p = rfftw_create_plan(n, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
	rfftw_one(p, (fftw_real *) complex_data, (fftw_real *) real_data);
	
	rfftw_destroy_plan(p);
#endif	//DJBFFT
	
	return 0;
	
}



int EMfft::real_to_complex_nd(float *real_data, float *complex_data, int nx, int ny, int nz)
{
/*	// in-place or out-of-place?
	FFTPLACE fftplace = (real_data == complex_data) 
		? FFT_IN_PLACE
		: FFT_OUT_OF_PLACE;
	FftwPlan *fwplan = make_plan(nx, ny, nz, REAL_TO_COMPLEX, fftplace);
	rfftwnd_one_real_to_complex(fwplan->get_plan_nd(), (fftw_real *) real_data,
								(fftw_complex *) complex_data);
	
	delete fwplan;
	rfftwnd_destroy_plan(fwplan->get_plan_nd());
*/
/*	FFTPLACE fftplace = (real_data == complex_data) 
		? FFT_IN_PLACE
		: FFT_OUT_OF_PLACE;
	
	int rank = get_rank(ny, nz);
	int flags = FFTW_ESTIMATE | FFTW_USE_WISDOM | FFTW_THREADSAFE;
	if(FFT_OUT_OF_PLACE == fftplace) {
		flags |= FFTW_OUT_OF_PLACE;
	}
	else {
		flags |= FFTW_IN_PLACE;
	}
*/
	int rank = get_rank(ny, nz);
	int dims[3];
	dims[0] = nz;
	dims[1] = ny;
	dims[2] = nx;
	rfftwnd_plan plan_nd = 
		rfftwnd_create_plan(rank, dims + (3 - rank), FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
	rfftwnd_one_real_to_complex(plan_nd, (fftw_real *) real_data,
								(fftw_complex *) complex_data);
	rfftwnd_destroy_plan(plan_nd);
	return 0;
}

int EMfft::complex_to_real_nd(float *complex_data, float *real_data, int nx, int ny, int nz)
{
	int rank = get_rank(ny, nz);
	int dims[3];
	dims[0] = nz;
	dims[1] = ny;
	dims[2] = nx;
	
//	FftwPlan *fwplan = make_plan(nx, ny, nz, COMPLEX_TO_REAL);
	
	rfftwnd_plan plan_nd = 
		rfftwnd_create_plan(rank, dims + (3 - rank), FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);	
								
	rfftwnd_one_complex_to_real(plan_nd, (fftw_complex *) complex_data,
								(fftw_real *) real_data);					

	rfftwnd_destroy_plan(plan_nd);
	return 0;
}
/*
FftwPlan *EMfft::make_plan(int nx, int ny, int nz, 
		                   FftwDirection dir,
						   FFTPLACE fftplace)
{
	int rank = get_rank(ny, nz);
	int flags = FFTW_ESTIMATE | FFTW_USE_WISDOM | FFTW_THREADSAFE;
	if (dir == REAL_TO_COMPLEX) {
		if (FFT_OUT_OF_PLACE == fftplace) {
			flags |= FFTW_OUT_OF_PLACE;
		} else {
			flags |= FFTW_IN_PLACE;
		}
	}
	else if (dir == COMPLEX_TO_REAL) {
		// Always in place because the inverse transform
		// destroys the input array.
		flags |= FFTW_IN_PLACE;
	}


	// plan caching seems to fail if out-of-place transform
	// is followed by an in-place transform, since the in-place
	// transform then uses the out-of-place plan (which
	// has the wrong flags)
	list < FftwPlan * >::iterator qi = fwplans.begin();
	for (; qi != fwplans.end(); qi++) {
		if (**qi == FftwPlan(rank, nx, ny, nz, dir, flags)) {
			break;
		}
	}
	if (qi != fwplans.end()) {
		return *qi;
	}

	if (fwplans.size() == MAX_NUM_PLANS) {
		delete fwplans.front();
		fwplans.pop_front();
	}

	FftwPlan *new_plan = new FftwPlan(rank, nx, ny, nz, dir, flags);
	fwplans.push_back(new_plan);
	return new_plan;
}
*/
#endif	//FFTW2

#ifdef FFTW3

int EMfft::real_to_complex_1d(float *real_data, float *complex_data, int n)
{
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
	int rank = get_rank(ny, nz);
	int dims[3];
	dims[0] = nz;
	dims[1] = ny;
	dims[2] = nx;

	fftwf_plan plan = fftwf_plan_dft_r2c(rank, dims + (3 - rank), real_data,
										 (fftwf_complex *) complex_data, FFTW_ESTIMATE);

	fftwf_execute(plan);
	fftwf_destroy_plan(plan);
	return 0;
}

int EMfft::complex_to_real_nd(float *complex_data, float *real_data, int nx, int ny, int nz)
{
	int rank = get_rank(ny, nz);
	int dims[3];
	dims[0] = nz;
	dims[1] = ny;
	dims[2] = nx;

	fftwf_plan plan = fftwf_plan_dft_c2r(rank, dims + (3 - rank), (fftwf_complex *) complex_data,
										 real_data, FFTW_ESTIMATE);

	fftwf_execute(plan);
	fftwf_destroy_plan(plan);
	return 0;
}

#endif	//FFTW3

#ifdef NATIVE_FFT
#include "sparx/native_fft.h"
int EMfft::real_to_complex_1d(float *real_data, float *complex_data, int n)
{
	int complex_size = n + 2 - n%2;
	
	memcpy(complex_data, real_data, n * sizeof(float));
	float * work = (float*) malloc(complex_size*sizeof(float));
	if (!work) {
		fprintf(stderr,"real_to_complex_1d: failed to allocate work\n");
		LOGERR("real_to_complex_1d: failed to allocate work\n");	
	}
	
	Nativefft::fmrs_1rf(complex_data, work, n);
	
	//the complex part need negate to be consistent with fftw result
	int i;
	for(i=1; i<complex_size; i+=2) {
		complex_data[i] *= -1;
	}
	
	free(work);
	return 0;
}

int EMfft::complex_to_real_1d(float *complex_data, float *real_data, int n)
{
	int complex_size = n + 2 - n%2;
	
	//the complex part need negate to be consistent with fftw result
	int i;
	for(i=1; i<complex_size; i+=2) {
		complex_data[i] *= -1;
	}
	
	//here, n is the "logical" size of DFT, not the array size, n is the real array size
	memcpy(real_data, complex_data, complex_size * sizeof(float));
	float * work = (float*)malloc(complex_size*sizeof(float));
	if (!work) {
		fprintf(stderr,"real_to_complex_1d: failed to allocate work\n");
		LOGERR("complex_to_real_1d: failed to allocate work\n");	
	}
	
	Nativefft::fmrs_1rb(real_data, work, n);
	
	free(work);
	return 0;
}

int EMfft::real_to_complex_nd(float *real_data, float *complex_data, int nx, int ny, int nz)
{
	const int rank = get_rank(ny, nz);
	const int complex_nx = nx + 2 - nx%2;
	
	switch(rank) {
		case 1:		//for 1D fft
			real_to_complex_1d(real_data, complex_data, nx);
			return 0;;
		case 2:		//for 2D fft
		{
			int j;
			for (j = 0; j < ny; j++) {
				memcpy(&complex_data[complex_nx*j], &real_data[nx*j], nx*sizeof(float));	
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
	   		
	   		free(work);
	   		return 0;
		}
	   	case 3:		//for 3D fft
	   	{
			int j, k;
			for (k = 0; k<nz; k++) {
	      		for (j = 0; j < ny; j++) {
	         		memcpy(&complex_data[complex_nx*ny*k+j*complex_nx], &real_data[nx*ny*k+j*nx], nx*sizeof(float));
	      		}
			}
			float * work = (float*) malloc(complex_nx*sizeof(float));
	   		if (!work) {
	   			fprintf(stderr,"real_to_complex_nd(3df): failed to allocate work\n");
	   			LOGERR("real_to_complex_nd(3df): failed to allocate work\n");
	   		}
	   		
	   		// 3d inplace fft, overwrite complex_data
	   		int status = Nativefft::fmrs_3rf(complex_data, work, complex_nx, nx, ny, nz);
	   		if (status !=0) {
	      		fprintf(stderr, "real_to_complex_nd(3df): status = %d\n", status);
	      		LOGWARN("real_to_complex_nd(3df): status = %d\n", status);
	   		}
	   		
	   		free(work);
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
		{
			memcpy(real_data, complex_data, complex_nx*ny*sizeof(float));
			
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
		}
		case 3:		//for 3D ift
		{
			memcpy(real_data, complex_data, complex_nx*ny*nz*sizeof(float));
			
			float * work = (float*) malloc(complex_nx*sizeof(float));
	   		if (!work) {
	   			fprintf(stderr,"complex_to_real_nd(3db): failed to allocate work\n");
	   			LOGERR("complex_to_real_nd(3db): failed to allocate work\n");
	   		}
	   		
	   		// 3d inplace fft, overwrite real_data
	   		int status = Nativefft::fmrs_3rb(real_data, work, complex_nx, nx, ny, nz);
	   		if (status !=0) {
	      		fprintf(stderr, "complex_to_real_nd(3db): status = %d\n", status);
	      		LOGWARN("complex_to_real_nd(3db): status = %d\n", status);
	   		}
	   		
	   		free(work);
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
	int complex_n = n + 2 - n%2;	//the size for 1D complex array
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
