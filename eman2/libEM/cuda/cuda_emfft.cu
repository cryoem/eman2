#include <cufft.h>
#include <cuda_runtime_api.h>
#include "cuda_defs.h"
#include "cuda_util.h"

const int real_2_complex = 1;
const int complex_2_real = 2;
const int EMCUDA_FFT_CACHE_SIZE = 8;
int EMCUDA_FFT_DD_CACHE_NUM_PLANS = -1;

// This struct is the device to device cufft plan cache. This is useful when your data exists solely on the GPU
// In an ideal world the cufft_hh_plan_cache would have one of these as a variable (a lot of copying is occuring)
struct cufft_dd_plan_cache {
	int rank; // Store the rank of the plan
	int plan_dims[3]; // Store the dimensions of the plan (always in 3D, if dimensions are "unused" they are taken to be 1)
	int r2c; // store whether the plan was real to complex or vice versa
	cufftHandle handle;  // Store the plans themselves
	int ip; // Store whether or not the plan was inplace
	int batch;
//	int batch; // Batch fft number - applicable only to 1D
};

cufft_dd_plan_cache* CudaDDFftPlanCache[EMCUDA_FFT_CACHE_SIZE];

void print_cuda_fft_fail(cufftResult result) {
	if ( result == CUFFT_SUCCESS) {
		//printf("CUDA FFT sucess!\n");
		return;
	}
	if (result == CUFFT_INVALID_PLAN ) {
		printf("CUFFT invalid plan\n");
	} else if ( result == CUFFT_INVALID_VALUE ) {
		printf("CUFFT invalid value\n");
	}else if ( result == CUFFT_SETUP_FAILED ) {
		printf("CUFFT setup failed\n");
	}else if ( result == CUFFT_EXEC_FAILED ) {
		printf("CUFFT exec failed\n");
	}else if ( result == CUFFT_ALLOC_FAILED ) {
		printf("CUFFT alloc failed\n");
	}else if ( result == CUFFT_INVALID_TYPE ) {
		printf("CUFFT invalid type\n");
	}else if ( result == CUFFT_INTERNAL_ERROR) {
		printf("CUFFT internal error\n");
	}else if ( result == CUFFT_INVALID_SIZE ) {
	   printf("CUFFT invalid size\n");
	}else {
	printf("CUFFT success\n");}
}

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

int cuda_fft_dd_plan_cache_params_match(const cufft_dd_plan_cache* const c, const int x, const int y, const int z, const int rank_in, const int r2c_flag, const int ip_flag, const int batch) 
{
	
	if (c->plan_dims[0]==x && c->plan_dims[1]==y && c->plan_dims[2]==z 
					  && c->rank==rank_in && c->r2c==r2c_flag && c->ip==ip_flag && c->batch==batch) {
		return 1;
	}
	return 0;
}

cufft_dd_plan_cache* get_cuda_fft_dd_plan_cache(const int rank_in, const int x, const int y, const int z, const int r2c_flag, const int ip_flag, const int batch) 
{
	// First check to see if we already have the plan, last one on array is first in line(make managment faster)

	for (int i=EMCUDA_FFT_DD_CACHE_NUM_PLANS; i>=0; i--) {
		//printf("looping %d\n", i);
		cufft_dd_plan_cache* c = CudaDDFftPlanCache[i];
		if (cuda_fft_dd_plan_cache_params_match(c,x,y,z,rank_in,r2c_flag,ip_flag, batch)) {
			//printf("Already have that cache\n");
			return c;
		}
	}	

	cufftHandle plan;
	cufftResult result;

	if ( y == 1 )
	{ /* 1D */
		if ( r2c_flag == real_2_complex ) {
			result = cufftPlan1d(&plan,x,CUFFT_R2C, batch);
		}
		else {
			// r2c_flag == complex_2_real, this is guaranteed by the error checking at the beginning of the function
			result = cufftPlan1d(&plan,x,CUFFT_C2R, batch);
		}
	} 
	else if ( z == 1 )
	{
		if ( r2c_flag == real_2_complex ) {
			result = cufftPlan2d(&plan,x,y,CUFFT_R2C);
		}
		else {
			// r2c_flag == complex_2_real, this is guaranteed by the error checking at the beginning of the function
			result = cufftPlan2d(&plan,x,y,CUFFT_C2R);
		}
	}
	else /* 3D */ {
		if ( r2c_flag == real_2_complex ) {
			result = cufftPlan3d(&plan,x,y,z,CUFFT_R2C);
		}
		else {
			// r2c_flag == complex_2_real, this is guaranteed by the error checking at the beginning of the function 
			result = cufftPlan3d(&plan,x,y,z,CUFFT_C2R);
		}
	}
	print_cuda_fft_fail(result);
	
	if (EMCUDA_FFT_DD_CACHE_NUM_PLANS == EMCUDA_FFT_CACHE_SIZE )
	{
		cufftDestroy(CudaDDFftPlanCache[0]->handle);
		delete CudaDDFftPlanCache[0];
		for(int i=0; i<(EMCUDA_FFT_DD_CACHE_NUM_PLANS-1); i++) {
			CudaDDFftPlanCache[i] = CudaDDFftPlanCache[i+1];
		}	
	}else{
		EMCUDA_FFT_DD_CACHE_NUM_PLANS++;
	}

	cufft_dd_plan_cache* c = new cufft_dd_plan_cache;
	c->plan_dims[0]=x;
	c->plan_dims[1]=y;
	c->plan_dims[2]=z;
	c->r2c=r2c_flag;
	c->ip=ip_flag;
	c->handle = plan;
	c->rank =rank_in;
	c->batch = batch;
	//c->batch =batch;

	CudaDDFftPlanCache[EMCUDA_FFT_DD_CACHE_NUM_PLANS] = c;
	
	return c;
}

void do_cuda_fft_cache_destroy()
{
	for (int i = 0; i <= EMCUDA_FFT_DD_CACHE_NUM_PLANS; i++) {
		cufftDestroy(CudaDDFftPlanCache[i]->handle);
		delete CudaDDFftPlanCache[i];
	}
}

int cuda_dd_fft_real_to_complex_nd(float *real_data, float *complex_data, int nx, int ny, int nz, int batch)
{
	//device_init(); //this breaks CUDA FFT!!!! (caused by setCudaDevice)
	const int rank = get_rank(ny, nz);
	bool ip;
	cufft_dd_plan_cache* cache = 0;
	//printf("Step 1, nx %d, ny %d, nz %d, rank %d\n",nx,ny,nz,rank);
	//cudaError_t error = cudaGetLastError();
	
	ip = ( complex_data == real_data );

	if ( !ip ) {
		cache = get_cuda_fft_dd_plan_cache(rank,nx,ny,nz,real_2_complex,ip,batch);
		cufftResult result = cufftExecR2C(cache->handle, real_data, (cufftComplex*)complex_data );
		print_cuda_fft_fail(result);
		cudaThreadSynchronize();
	}
	else {
		cache = get_cuda_fft_dd_plan_cache(rank,nx,ny,nz,real_2_complex,ip,batch);
		/// CHECK LATER - Not sure if this will work
		cufftResult result = cufftExecR2C(cache->handle, (cufftReal*)complex_data, (cufftComplex*)complex_data );
		print_cuda_fft_fail(result);
		cudaThreadSynchronize();
	}

	
	return 0;
}

int cuda_dd_fft_complex_to_real_nd(float *complex_data, float *real_data, int nx, int ny, int nz, int batch)
{
	//device_init();
	const int rank = get_rank(ny, nz);
	bool ip;
	cufftResult result;
	//printf("using nx2 %d\n",nx2);
	cufft_dd_plan_cache* cache = 0;

	ip = ( complex_data == real_data );
	cache = get_cuda_fft_dd_plan_cache(rank,nx,ny,nz,complex_2_real,ip,batch);
	result = cufftExecC2R(cache->handle, (cufftComplex*)complex_data, real_data);
	print_cuda_fft_fail(result);
	cudaThreadSynchronize();
			
	return 0;
}
