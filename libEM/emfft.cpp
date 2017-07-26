/**
 * $Id$
 */
 
/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
 * Copyright (c) 2000-2006 Baylor College of Medicine
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
 * */

#include <string>
#include <cstring>
#include "emfft.h"
#include "log.h"

#include <iostream>
using std::cout;
using std::endl;

#include "util.h"

#ifdef EMAN2_USING_CUDA
#include "cuda/cuda_emfft.h"
#endif 

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
#ifdef _WIN32
MUTEX fft_mutex;
#else
pthread_mutex_t fft_mutex=PTHREAD_MUTEX_INITIALIZER;
#endif

#ifdef FFTW_PLAN_CACHING
// The only thing important about these constants is that they don't equal each other
const int EMfft::EMAN2_REAL_2_COMPLEX = 1;
const int EMfft::EMAN2_COMPLEX_2_REAL = 2;


#ifdef FFTW3
EMfft::EMfftw3_cache::EMfftw3_cache() :
		num_plans(0)
{
	for(int i = 0; i < EMFFTW3_CACHE_SIZE; ++i)
	{
		rank[i] = 0;
		plan_dims[i][0] = 0; plan_dims[i][1] = 0; plan_dims[i][2] = 0;
		r2c[i] = -1;
		ip[i] = -1;
		fftwplans[i] = NULL;
	}
}

void EMfft::EMfftw3_cache::debug_plans()
{
	for(int i = 0; i < EMFFTW3_CACHE_SIZE; ++i)
	{
		cout << "Plan " << i << " has dims " << plan_dims[i][0] << " " 
				<< plan_dims[i][1] << " " << 
				plan_dims[i][2] << ", rank " <<
				rank[i] << ", rc flag " 
				<< r2c[i] << ", ip flag " << ip[i] << endl;
	}
}

EMfft::EMfftw3_cache::~EMfftw3_cache()
{
	for(int i = 0; i < EMFFTW3_CACHE_SIZE; ++i)
	{
		if (fftwplans[i] != NULL)
		{
			int mrt = Util::MUTEX_LOCK(&fft_mutex);
			fftwf_destroy_plan(fftwplans[i]);
			mrt = Util::MUTEX_UNLOCK(&fft_mutex);
			fftwplans[i] = NULL;
		}
	}
}

fftwf_plan EMfft::EMfftw3_cache::get_plan(const int rank_in, const int x, const int y, const int z, const int r2c_flag, const int ip_flag, fftwf_complex* complex_data, float* real_data )
{

	if ( rank_in > 3 || rank_in < 1 ) throw InvalidValueException(rank_in, "Error, can not get an FFTW plan using rank out of the range [1,3]");
	if ( r2c_flag != EMAN2_REAL_2_COMPLEX && r2c_flag != EMAN2_COMPLEX_2_REAL ) throw InvalidValueException(r2c_flag, "The real two complex flag is not supported");
	
// 	static int num_added = 0;
// 	cout << "Was asked for " << rank_in << " " << x << " " << y << " " << z << " " << r2c_flag << endl;
	
	int dims[3];
	dims[0] = z;
	dims[1] = y;
	dims[2] = x;
	
	// First check to see if we already have the plan
	int i;
	for (i=0; i<num_plans; i++) {
		if (plan_dims[i][0]==x && plan_dims[i][1]==y && plan_dims[i][2]==z 
				  && rank[i]==rank_in && r2c[i]==r2c_flag && ip[i]==ip_flag) return fftwplans[i];
	}
	
	int mrt = Util::MUTEX_LOCK(&fft_mutex);
	
	fftwf_plan plan;
	// Create the plan
	if ( y == 1 && z == 1 )
	{
		if ( r2c_flag == EMAN2_REAL_2_COMPLEX )
			plan = fftwf_plan_dft_r2c_1d(x, real_data, complex_data, FFTW_ESTIMATE);
		else // r2c_flag == EMAN2_COMPLEX_2_REAL, this is guaranteed by the error checking at the beginning of the function
			plan = fftwf_plan_dft_c2r_1d(x, complex_data, real_data, FFTW_ESTIMATE);
	}
	else
	{
		if ( r2c_flag == EMAN2_REAL_2_COMPLEX )
			plan = fftwf_plan_dft_r2c(rank_in, dims + (3 - rank_in), real_data, complex_data, FFTW_ESTIMATE);
		else // r2c_flag == EMAN2_COMPLEX_2_REAL, this is guaranteed by the error checking at the beginning of the function
			plan = fftwf_plan_dft_c2r(rank_in, dims + (3 - rank_in), complex_data, real_data, FFTW_ESTIMATE);
	}

	if (fftwplans[EMFFTW3_CACHE_SIZE-1] != NULL )
	{
		fftwf_destroy_plan(fftwplans[EMFFTW3_CACHE_SIZE-1]);
		fftwplans[EMFFTW3_CACHE_SIZE-1] = NULL;
	}
	
	mrt = Util::MUTEX_UNLOCK(&fft_mutex);
				
	int upper_limit = num_plans;
	if ( upper_limit == EMFFTW3_CACHE_SIZE ) upper_limit -= 1;
	for (int i=upper_limit-1; i>0; i--)
	{
		fftwplans[i]=fftwplans[i-1];
		rank[i]=rank[i-1];
		r2c[i]=r2c[i-1];
		ip[i]=ip[i-1];
		plan_dims[i][0]=plan_dims[i-1][0];
		plan_dims[i][1]=plan_dims[i-1][1];
		plan_dims[i][2]=plan_dims[i-1][2];
	}
		//dimplan[0]=-1;

	plan_dims[0][0]=x;
	plan_dims[0][1]=y;
	plan_dims[0][2]=z;
	r2c[0]=r2c_flag;
	ip[0]=ip_flag;
	fftwplans[0] = plan;
	rank[0]=rank_in;
	if (num_plans<EMFFTW3_CACHE_SIZE) num_plans++;
// 			debug_plans();
// 			cout << "Created plan 0" << endl;
// 	++num_added;
// 	cout << "I have created " << num_added << " plans" << endl;
	return fftwplans[0];

}

// Static init
EMfft::EMfftw3_cache EMfft::plan_cache;

#endif // FFTW3

#endif // FFTW_PLAN_CACHING

//#ifdef CUDA_FFT
//int EMfft::real_to_complex_1d(float *real_data, float *complex_data, int n)
//{
//	return  cuda_fft_real_to_complex_1d(real_data,complex_data,n);
//}

//int EMfft::complex_to_real_1d(float *complex_data, float *real_data, int n)
//{
//	return cuda_fft_complex_to_real_1d(complex_data,real_data,n);
//}

//int EMfft::real_to_complex_nd(float *real_data, float *complex_data, int nx, int ny, int nz)
//{
//	return cuda_fft_real_to_complex_nd(real_data,complex_data,nx,ny,nz);
//}

//int EMfft::complex_to_real_nd(float *complex_data, float *real_data, int nx, int ny, int nz)
//{
//	return cuda_fft_complex_to_real_nd(complex_data,real_data,nx,ny,nz);
//}

//#endif

#ifdef FFTW3

int EMfft::real_to_complex_1d(float *real_data, float *complex_data, int n)
{//cout<<"doing fftw3"<<endl;
#ifdef FFTW_PLAN_CACHING
	bool ip = ( complex_data == real_data );
	fftwf_plan plan = plan_cache.get_plan(1,n,1,1,EMAN2_REAL_2_COMPLEX,ip,(fftwf_complex *) complex_data, real_data);
	// According to FFTW3, this is making use of the "guru" interface - this is necessary if plans are to be reused
	fftwf_execute_dft_r2c(plan, real_data,(fftwf_complex *) complex_data);
#else
	int mrt = Util::MUTEX_LOCK(&fft_mutex);
	fftwf_plan plan = fftwf_plan_dft_r2c_1d(n, real_data, (fftwf_complex *) complex_data,
											FFTW_ESTIMATE);
	mrt = Util::MUTEX_UNLOCK(&fft_mutex);

	fftwf_execute(plan);
	mrt = Util::MUTEX_LOCK(&fft_mutex);
	fftwf_destroy_plan(plan);
	mrt = Util::MUTEX_UNLOCK(&fft_mutex);
#endif // FFTW_PLAN_CACHING
	return 0;
};

int EMfft::complex_to_real_1d(float *complex_data, float *real_data, int n)
{
#ifdef FFTW_PLAN_CACHING
	bool ip = ( complex_data == real_data );
	fftwf_plan plan = plan_cache.get_plan(1,n,1,1,EMAN2_COMPLEX_2_REAL,ip,(fftwf_complex *) complex_data, real_data);
	// According to FFTW3, this is making use of the "guru" interface - this is necessary if plans are to be reused
	fftwf_execute_dft_c2r(plan, (fftwf_complex *) complex_data, real_data);
#else
	int mrt = Util::MUTEX_LOCK(&fft_mutex);
	fftwf_plan plan = fftwf_plan_dft_c2r_1d(n, (fftwf_complex *) complex_data, real_data,
											FFTW_ESTIMATE);
	mrt = Util::MUTEX_UNLOCK(&fft_mutex);
	fftwf_execute(plan);
	mrt = Util::MUTEX_LOCK(&fft_mutex);
	fftwf_destroy_plan(plan);
	mrt = Util::MUTEX_UNLOCK(&fft_mutex);
#endif // FFTW_PLAN_CACHING
	
	return 0;
}


// ming add c->c fft with fftw3 library//
int EMfft::complex_to_complex_1d(float *complex_data_in, float *complex_data_out, int n)
{
	fftwf_plan p;
	fftwf_complex *in=(fftwf_complex *) complex_data_in;
	fftwf_complex *out=(fftwf_complex *) complex_data_out;
	int mrt = Util::MUTEX_LOCK(&fft_mutex);
	p=fftwf_plan_dft_1d(n/2,in,out, FFTW_FORWARD, FFTW_ESTIMATE);
	mrt = Util::MUTEX_UNLOCK(&fft_mutex);
	fftwf_execute(p);
	mrt = Util::MUTEX_LOCK(&fft_mutex);
	fftwf_destroy_plan(p);
	mrt = Util::MUTEX_UNLOCK(&fft_mutex);
	return 0;
}


int EMfft::complex_to_complex_nd(float *in, float *out, int nx,int ny,int nz)
{
	const int rank = get_rank(ny, nz);
	int dims[3];
	dims[0] = nz;
	dims[1] = ny;
	dims[2] = nx;

	switch(rank) {
		case 1:
			complex_to_complex_1d(in, out, nx);
			break;

		case 2:
		case 3:
		{
			fftwf_plan p;

			int mrt = Util::MUTEX_LOCK(&fft_mutex);

			if(out == in) {
				p=fftwf_plan_dft_3d(nx/2,ny,nz,(fftwf_complex *) in,(fftwf_complex *) out, FFTW_FORWARD, FFTW_ESTIMATE);
			}
			else {

				p=fftwf_plan_dft_3d(nx/2,ny,nz,(fftwf_complex *) in,(fftwf_complex *) out, FFTW_FORWARD, FFTW_ESTIMATE);
			}
			mrt = Util::MUTEX_UNLOCK(&fft_mutex);

			fftwf_execute(p);
			
			mrt = Util::MUTEX_LOCK(&fft_mutex);
			fftwf_destroy_plan(p);
			mrt = Util::MUTEX_UNLOCK(&fft_mutex);

		}
	}
	return 0;
}
// end ming


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
#ifdef FFTW_PLAN_CACHING
			bool ip = ( complex_data == real_data );
			fftwf_plan plan = plan_cache.get_plan(rank,nx,ny,nz,EMAN2_REAL_2_COMPLEX,ip,(fftwf_complex *) complex_data, real_data);
			// According to FFTW3, this is making use of the "guru" interface - this is necessary if plans are to be re-used
			fftwf_execute_dft_r2c(plan, real_data,(fftwf_complex *) complex_data );
#else
			int mrt = Util::MUTEX_LOCK(&fft_mutex);
			fftwf_plan plan = fftwf_plan_dft_r2c(rank, dims + (3 - rank), 
					real_data, (fftwf_complex *) complex_data, FFTW_ESTIMATE);
			mrt = Util::MUTEX_UNLOCK(&fft_mutex);
			
			fftwf_execute(plan);
			
			mrt = Util::MUTEX_LOCK(&fft_mutex);
			fftwf_destroy_plan(plan);
			mrt = Util::MUTEX_UNLOCK(&fft_mutex);

#endif // FFTW_PLAN_CACHING
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
#ifdef FFTW_PLAN_CACHING
			bool ip = ( complex_data == real_data );
			fftwf_plan plan = plan_cache.get_plan(rank,nx,ny,nz,EMAN2_COMPLEX_2_REAL,ip,(fftwf_complex *) complex_data, real_data);
			// According to FFTW3, this is making use of the "guru" interface - this is necessary if plans are to be re-used
			fftwf_execute_dft_c2r(plan, (fftwf_complex *) complex_data, real_data);
#else
			int mrt = Util::MUTEX_LOCK(&fft_mutex);
			fftwf_plan plan = fftwf_plan_dft_c2r(rank, dims + (3 - rank), 
					(fftwf_complex *) complex_data, real_data, FFTW_ESTIMATE);
			mrt = Util::MUTEX_UNLOCK(&fft_mutex);

			fftwf_execute(plan);
			
			mrt = Util::MUTEX_LOCK(&fft_mutex);
			fftwf_destroy_plan(plan);
			mrt = Util::MUTEX_UNLOCK(&fft_mutex);

#endif // FFTW_PLAN_CACHING
			
			
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
