

#include <cufft.h>
#include <cuda_runtime_api.h>
#include "cuda_defs.h"
#include "cuda_util.h"


const int EMCUDA_FFT_CACHE_SIZE = 30;
int EMCUDA_FFT_HH_CACHE_NUM_PLANS = 0;
int EMCUDA_FFT_DD_CACHE_NUM_PLANS = 0;


void print_cuda_fft_fail(cufftResult result) {
	if ( result == CUFFT_SUCCESS) return;
	
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
	}/*else if ( result == CUFFT_SHUTDOWN_FAILED ) {
	printf("CUFFT shutdown failed\n");
}*/else if ( result == CUFFT_INVALID_SIZE ) {
	   printf("CUFFT invalid size\n");
}else {
	printf("CUFFT success\n");
}
}


// This struct is the device to device cufft plan cache. This is useful when your data exists solely on the GPU
// In an ideal world the cufft_hh_plan_cache would have one of these as a variable (a lot of copying is occuring)
struct cufft_dd_plan_cache {
	int rank; // Store the rank of the plan
	int plan_dims[3]; // Store the dimensions of the plan (always in 3D, if dimensions are "unused" they are taken to be 1)
	int r2c; // store whether the plan was real to complex or vice versa
	cufftHandle handle;  // Store the plans themselves
	int ip; // Store whether or not the plan was inplace
	int batch; // Batch fft number - applicable only to 1D
};

// This struct is the host to host cufft plan cache. This means that the front end off the FFT mechanism takes and returns data that is valid on the host.
// Thus, a lot of memory copying results. As a result each cache has two device pointers to store the real and complex data
struct cufft_hh_plan_cache {
	cufft_dd_plan_cache dd_cache; // Generic cache parts
	cufftReal* real; // Device real pointer
	cufftComplex* complex; // Device complex pointer 
};


cufft_hh_plan_cache CudaHHFftPlanCache[EMCUDA_FFT_CACHE_SIZE];
cufft_dd_plan_cache CudaDDFftPlanCache[EMCUDA_FFT_CACHE_SIZE];


void reset_cuda_fft_dd_plan_cache(cufft_dd_plan_cache* c) {
	c->rank = 0;
	c->plan_dims[0] = 0; c->plan_dims[1] = 0; c->plan_dims[2] = 0;
	c->r2c = -1;
	c->ip = -1;
	if (c->handle != 0) {
// 		printf("Destroying plan %d %d\n",c->handle,EMCUDA_FFT_DD_CACHE_NUM_PLANS);
		cufftDestroy(c->handle);
	}
	c->handle = 0;
	c->batch = 0;
}

void print_fft_dd(cufft_dd_plan_cache* c, const int i) {
	printf("%d: %d %d %d %d %d %d %d %d\n",i,c->rank,c->plan_dims[0],c->plan_dims[1],c->plan_dims[2],c->r2c,c->ip, c->handle,c->batch );
}


void reset_cuda_fft_hh_plan_cache(cufft_hh_plan_cache* c) {
	reset_cuda_fft_dd_plan_cache(&(c->dd_cache));
	if (c->real != 0) {
		cudaFree(c->real);
		c->real = 0;
	}
	if (c->complex != 0) {
		cudaFree(c->complex);
		c->complex = 0;
	}
}

void copy_cuda_fft_dd_plan_cache(cufft_dd_plan_cache* to,  cufft_dd_plan_cache* from) {
	to->rank = from->rank;
	to->plan_dims[0] = from->plan_dims[0];
	to->plan_dims[1] = from->plan_dims[1];
	to->plan_dims[2] = from->plan_dims[2];
	to->r2c = from->r2c;
	to->ip = from->ip;
	to->handle = from->handle;
	to->batch = from->batch;
}

void copy_cuda_fft_hh_plan_cache(cufft_hh_plan_cache* to,  cufft_hh_plan_cache* from) {

	copy_cuda_fft_dd_plan_cache(&(to->dd_cache),&(from->dd_cache));
	to->real = from->real;
	to->complex = from->complex;
}

void init_cuda_fft_dd_plan_cache() {
	for(int i = 0; i < EMCUDA_FFT_CACHE_SIZE; ++i)
	{
		reset_cuda_fft_dd_plan_cache(&(CudaDDFftPlanCache[i]));
	}
}

void init_cuda_fft_hh_plan_cache() {
	for(int i = 0; i < EMCUDA_FFT_CACHE_SIZE; ++i)
	{
		reset_cuda_fft_hh_plan_cache(&(CudaHHFftPlanCache[i]));
	}
}

void cleanup_cuda_fft_dd_plan_cache()
{	
	for(int i = 0; i < EMCUDA_FFT_CACHE_SIZE; ++i)
	{
		reset_cuda_fft_dd_plan_cache(&CudaDDFftPlanCache[i]);
	}
}
void debug_fft_dd_plan_cache()
{	
	for(int i = 0; i < EMCUDA_FFT_CACHE_SIZE; ++i)
	{
		print_fft_dd(&CudaDDFftPlanCache[i],i);
	}
}

void cleanup_cuda_fft_hh_plan_cache()
{	
	for(int i = 0; i < EMCUDA_FFT_CACHE_SIZE; ++i)
	{
		reset_cuda_fft_hh_plan_cache(&CudaHHFftPlanCache[i]);
	}
}

int cuda_fft_dd_plan_cache_params_match(const cufft_dd_plan_cache* const c, const int x, const int y, const int z, const int rank_in, const int r2c_flag, const int ip_flag, const int batch ) {
	
	if (c->plan_dims[0]==x && c->plan_dims[1]==y && c->plan_dims[2]==z 
					  && c->rank==rank_in && c->r2c==r2c_flag && c->ip==ip_flag && c->batch == batch) {
		return 1;
	}
	return 0;
}

int real_2_complex = 1;
int complex_2_real = 2;

cufft_hh_plan_cache* get_cuda_fft_hh_plan_cache(const int rank_in, const int x, const int y, const int z, const int r2c_flag, const int ip_flag,const int batch) {

	if ( rank_in > 3 || rank_in < 1 ) throw; //InvalidValueException(rank_in, "Error, can not get an FFTW plan using rank out of the range [1,3]")
	if ( r2c_flag != real_2_complex && r2c_flag != complex_2_real ) throw; //InvalidValueException(r2c_flag, "The real two complex flag is not supported");
	
	// First check to see if we already have the plan
	int i;
	for (i=0; i<EMCUDA_FFT_HH_CACHE_NUM_PLANS; i++) {
		cufft_dd_plan_cache* c = &(CudaHHFftPlanCache[i].dd_cache);
		if (cuda_fft_dd_plan_cache_params_match(c,x,y,z,rank_in,r2c_flag,ip_flag,batch)) {
			//printf("Already have that cache\n");
			return &CudaHHFftPlanCache[i];
		}
	}
	//printf("Returning a new cache\n");
	cufftHandle plan;
	cufftReal* real; // Device real pointer
	cufftComplex* complex; // Device complex pointer 
	
	if (r2c_flag == complex_2_real) {
		int x2 = x + 2;
		int complex_mem_size = sizeof(cufftComplex) * x2 * y * z/2;
		int real_mem_size = sizeof(cufftReal) * x2 * y * z;
		cudaMalloc((void**)&complex, complex_mem_size);
		cudaMalloc((void**)&real, real_mem_size);
	} else {
		int offset = 2 - x%2;
		int x2 = x + offset;
		int complex_mem_size = sizeof(cufftComplex) * x2 * y * z/2;
		int real_mem_size = sizeof(cufftReal) * x * y * z;
		cudaMalloc((void**)&complex, complex_mem_size);
		cudaMalloc((void**)&real, real_mem_size);
	}
	// Create the plan
	if ( y == 1 && z == 1 )
	{
		if ( r2c_flag == real_2_complex ) {
			cufftPlan1d(&plan,x,CUFFT_R2C,batch);
		}
		else { // r2c_flag == complex_2_real, this is guaranteed by the error checking at the beginning of the function
			cufftPlan1d(&plan, x, CUFFT_C2R,batch);
		}
	}
	else if ( z == 1 )
	{
		if ( r2c_flag == real_2_complex ) {
			cufftPlan2d(&plan,x,y,CUFFT_R2C);
		}
		else // r2c_flag == complex_2_real, this is guaranteed by the error checking at the beginning of the function
			cufftPlan2d(&plan,x,y,CUFFT_C2R);
	}
	else /* 3D */ {
		if ( r2c_flag == real_2_complex ) {
			cufftPlan3d(&plan,x,y,z,CUFFT_R2C);
		}
		else // r2c_flag == complex_2_real, this is guaranteed by the error checking at the beginning of the function
			cufftPlan3d(&plan,x,y,z,CUFFT_C2R);
	}

	if (EMCUDA_FFT_HH_CACHE_NUM_PLANS == EMCUDA_FFT_CACHE_SIZE )
	{
		reset_cuda_fft_hh_plan_cache(&(CudaHHFftPlanCache[EMCUDA_FFT_CACHE_SIZE-1]));
	}
				
	int upper_limit = EMCUDA_FFT_HH_CACHE_NUM_PLANS+1;
	if ( upper_limit == EMCUDA_FFT_CACHE_SIZE ) upper_limit -= 1;
	
	for (int i=upper_limit-1; i>0; i--)
	{
		copy_cuda_fft_hh_plan_cache(&(CudaHHFftPlanCache[i]),&(CudaHHFftPlanCache[i-1]));
		
	}
	
	cufft_hh_plan_cache* c = &CudaHHFftPlanCache[0];
	
	c->dd_cache.plan_dims[0]=x;
	c->dd_cache.plan_dims[1]=y;
	c->dd_cache.plan_dims[2]=z;
	c->dd_cache.r2c=r2c_flag;
	c->dd_cache.ip=ip_flag;
	c->dd_cache.handle = plan;
	c->dd_cache.rank =rank_in;
	c->dd_cache.batch = batch;
	c->complex = complex;
	c->real = real;

	if (EMCUDA_FFT_HH_CACHE_NUM_PLANS<EMCUDA_FFT_CACHE_SIZE) EMCUDA_FFT_HH_CACHE_NUM_PLANS++;

	
	return c;
}

cufft_dd_plan_cache* get_cuda_fft_dd_plan_cache(const int rank_in, const int x, const int y, const int z, const int r2c_flag, const int ip_flag, const int batch) {

	if ( rank_in > 3 || rank_in < 1 ) throw; //InvalidValueException(rank_in, "Error, can not get an FFTW plan using rank out of the range [1,3]")
	if ( r2c_flag != real_2_complex && r2c_flag != complex_2_real ) throw; //InvalidValueException(r2c_flag, "The real two complex flag is not supported");
	
	// First check to see if we already have the plan
	int i;
	for (i=0; i<EMCUDA_FFT_DD_CACHE_NUM_PLANS; i++) {
		cufft_dd_plan_cache* c = &CudaDDFftPlanCache[i];
		if (cuda_fft_dd_plan_cache_params_match(c,x,y,z,rank_in,r2c_flag,ip_flag,batch)) {
			//printf("Already have that cache\n");
			return c;
		}
	}
	//printf("Returning a new cache\n");
	cufftHandle plan;
	// Create the plan
	cufftResult result;
	if ( y == 1 && z == 1 )
	{
		if ( r2c_flag == real_2_complex ) {
			result = cufftPlan1d(&plan,x,CUFFT_R2C,batch);
		}
		else {
			// r2c_flag == complex_2_real, this is guaranteed by the error checking at the beginning of the function
			result = cufftPlan1d(&plan, x, CUFFT_C2R,batch);
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
	if ( plan == 0 ) {
		printf("The plan is 0\n");
	}
	
	if (EMCUDA_FFT_DD_CACHE_NUM_PLANS == EMCUDA_FFT_CACHE_SIZE )
	{
		reset_cuda_fft_dd_plan_cache(&(CudaDDFftPlanCache[EMCUDA_FFT_CACHE_SIZE-1]));
	}
				
	int upper_limit = EMCUDA_FFT_DD_CACHE_NUM_PLANS+1;
	if ( upper_limit == EMCUDA_FFT_CACHE_SIZE ) upper_limit -= 1;
	
	for (int i=upper_limit-1; i>0; i--)
	{
		copy_cuda_fft_dd_plan_cache(&(CudaDDFftPlanCache[i]),&(CudaDDFftPlanCache[i-1]));
		
	}
	
	cufft_dd_plan_cache* c = &CudaDDFftPlanCache[0];
	
	c->plan_dims[0]=x;
	c->plan_dims[1]=y;
	c->plan_dims[2]=z;
	c->r2c=r2c_flag;
	c->ip=ip_flag;
	c->handle = plan;
	c->rank =rank_in;
	c->batch =batch;
	if (EMCUDA_FFT_DD_CACHE_NUM_PLANS<EMCUDA_FFT_CACHE_SIZE) EMCUDA_FFT_DD_CACHE_NUM_PLANS++;
// 	debug_fft_dd_plan_cache();
	return c;
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

int cuda_hh_fft_real_to_complex_1d(float *real_data, float *complex_data, int n,int batch)
{
	device_init();
	bool ip = false;
	int offset = 2 - n%2;
	int n2 = n + offset;
	int complex_mem_size = sizeof(cufftComplex) * n2/2;
	int real_mem_size = sizeof(cufftReal) * n;
	cufft_hh_plan_cache* cache =  get_cuda_fft_hh_plan_cache(1,n,1,1,real_2_complex,ip,batch);
	CUDA_SAFE_CALL(cudaMemcpy(cache->real,real_data, real_mem_size, cudaMemcpyHostToDevice));
	cufftExecR2C(cache->dd_cache.handle, cache->real, cache->complex );
	CUDA_SAFE_CALL(cudaMemcpy(complex_data,cache->complex, complex_mem_size, cudaMemcpyDeviceToHost));

	return 0;
};

int cuda_hh_fft_complex_to_real_1d(float *complex_data, float *real_data, int n,int batch)
{
	device_init();
	bool ip = false;
	int offset = 2 - n%2;
	int n2 = n + offset;
	int complex_mem_size = sizeof(cufftComplex) * n2/2;
	int real_mem_size = sizeof(cufftReal) * n2;
	cufft_hh_plan_cache* cache = get_cuda_fft_hh_plan_cache(1,n,1,1,complex_2_real,ip,batch);
	CUDA_SAFE_CALL(cudaMemcpy(cache->complex,complex_data, complex_mem_size, cudaMemcpyHostToDevice));
	cufftExecC2R(cache->dd_cache.handle, cache->complex, cache->real);
	CUDA_SAFE_CALL(cudaMemcpy(real_data,cache->real, real_mem_size, cudaMemcpyDeviceToHost));

	return 0;
}

int cuda_hh_fft_real_to_complex_nd(float *real_data, float *complex_data, int nx, int ny, int nz)
{
	device_init();
	const int rank = get_rank(ny, nz);
	bool ip;
	int offset = 2 - nx%2;
	int nx2 = nx + offset;
	int complex_mem_size = sizeof(cufftComplex) * nx2 * ny * nz/2;
	int real_mem_size = sizeof(cufftReal) * nx * ny * nz;
	cufft_hh_plan_cache* cache = 0;
#ifdef EMAN2_USING_CUDA_MALLOC	
	cudaStream_t stream1,stream2;
#endif
	switch(rank) {
		case 1:
			cuda_hh_fft_real_to_complex_1d(real_data, complex_data, nx,1);
			break;
		
		case 2:
		case 3:
			ip = ( complex_data == real_data );
			//ip = false;
			
			if ( !ip ) {
				cache = get_cuda_fft_hh_plan_cache(rank,nx,ny,nz,real_2_complex,ip,1);
#ifdef EMAN2_USING_CUDA_MALLOC	
				
				cudaStreamCreate(&stream1);
				cudaMemcpyAsync(cache->real,real_data, real_mem_size, cudaMemcpyHostToDevice,stream1);
				cudaThreadSynchronize();
				cudaStreamDestroy(stream1);
#else
				cudaMemcpy(cache->real,real_data, real_mem_size, cudaMemcpyHostToDevice);
#endif
				cufftExecR2C(cache->dd_cache.handle, cache->real, cache->complex );
				cudaThreadSynchronize();
#ifdef EMAN2_USING_CUDA_MALLOC	
				cudaStreamCreate(&stream2);
				cudaMemcpyAsync(complex_data,cache->complex, complex_mem_size, cudaMemcpyDeviceToHost,stream2);
				cudaThreadSynchronize();
				cudaStreamDestroy(stream2);
#else
				cudaMemcpy(complex_data,cache->complex, complex_mem_size,cudaMemcpyDeviceToHost);
#endif
			}
			else {
				cache = get_cuda_fft_hh_plan_cache(rank,nx,ny,nz,real_2_complex,ip,1);
#ifdef EMAN2_USING_CUDA_MALLOC	
				cudaStreamCreate(&stream1);
				cudaMemcpyAsync(cache->complex,real_data, complex_mem_size, cudaMemcpyHostToDevice,stream1);
				cudaThreadSynchronize();
				cudaStreamDestroy(stream1);
#else
				cudaMemcpy(cache->complex,real_data, complex_mem_size, cudaMemcpyHostToDevice);
#endif
				cufftExecR2C(cache->dd_cache.handle, (cufftReal*)cache->complex, cache->complex );
				cudaThreadSynchronize();
#ifdef EMAN2_USING_CUDA_MALLOC	
				cudaStreamCreate(&stream2);
				cudaMemcpyAsync(complex_data,cache->complex, complex_mem_size, cudaMemcpyDeviceToHost,stream2);
				cudaThreadSynchronize();
				cudaStreamDestroy(stream2);
#else
				cudaMemcpy(complex_data,cache->complex, complex_mem_size,cudaMemcpyDeviceToHost);
#endif
			}
		break;

		
		default:throw;
	}
	
	return 0;
}

int cuda_hh_fft_complex_to_real_nd(float *complex_data, float *real_data, int nx, int ny, int nz)
{
	device_init();
	const int rank = get_rank(ny, nz);
	bool ip;
	int offset = 2 - nx%2;
	int nx2 = nx + offset;
	//printf("using nx2 %d\n",nx2);
	int complex_mem_size = sizeof(cufftComplex) * nx2 * ny * nz/2;
	int real_mem_size = sizeof(cufftReal) * nx2 * ny * nz;
	cufft_hh_plan_cache* cache = 0;
#ifdef EMAN2_USING_CUDA_MALLOC
	cudaStream_t stream1,stream2;
#endif
	switch(rank) {
		case 1:
			cuda_hh_fft_complex_to_real_1d(complex_data, real_data, nx,1);
			break;
		
		case 2:
		case 3:
			ip = ( complex_data == real_data );
			cache = get_cuda_fft_hh_plan_cache(rank,nx,ny,nz,complex_2_real,ip,1);
	
#ifdef EMAN2_USING_CUDA_MALLOC			
			cudaStreamCreate(&stream1);
			cudaMemcpyAsync(cache->complex,complex_data, complex_mem_size, cudaMemcpyHostToDevice,stream1);
			cudaThreadSynchronize();
			cudaStreamDestroy(stream1);
#else
			cudaMemcpy(cache->complex,complex_data, complex_mem_size, cudaMemcpyHostToDevice);
#endif
			cufftExecC2R(cache->dd_cache.handle, cache->complex, cache->real);
			cudaThreadSynchronize();
	
#ifdef EMAN2_USING_CUDA_MALLOC
			cudaStreamCreate(&stream2);
			cudaMemcpyAsync(real_data,cache->real, real_mem_size, cudaMemcpyDeviceToHost,stream2);
			cudaThreadSynchronize();
			cudaStreamDestroy(stream2);
#else
			cudaMemcpy(real_data,cache->real, real_mem_size, cudaMemcpyDeviceToHost);
#endif
			break;
			
		default:throw;
	}
	
	return 0;
}


int cuda_dd_fft_real_to_complex_1d(float *real_data, float *complex_data, int n,int batch)
{
	device_init();
	bool ip = false;
	cufft_dd_plan_cache* cache =  get_cuda_fft_dd_plan_cache(1,n,1,1,real_2_complex,ip,batch);
	cufftResult result = cufftExecR2C(cache->handle, real_data, (cufftComplex*)complex_data );
	print_cuda_fft_fail(result);
	cudaThreadSynchronize();
	return 0;
};

int cuda_dd_fft_complex_to_real_1d(float *complex_data, float *real_data, int n, int batch)
{
	device_init();
	bool ip = false;
	cufft_dd_plan_cache* cache = get_cuda_fft_dd_plan_cache(1,n,1,1,complex_2_real,ip,batch);
	cufftResult result = cufftExecC2R(cache->handle, (cufftComplex*)complex_data,real_data);
	print_cuda_fft_fail(result);
	cudaThreadSynchronize();
	return 0;
}

int cuda_dd_fft_real_to_complex_nd(float *real_data, float *complex_data, int nx, int ny, int nz)
{
	device_init();
	const int rank = get_rank(ny, nz);
	bool ip;
	cufft_dd_plan_cache* cache = 0;
	cufftResult result;
	switch(rank) {
		case 1:
			cuda_dd_fft_real_to_complex_1d(real_data, complex_data, nx,1);
			break;
		
		case 2:
		case 3:
			ip = ( complex_data == real_data );
			//ip = false;
			
			if ( !ip ) {
				cache = get_cuda_fft_dd_plan_cache(rank,nx,ny,nz,real_2_complex,ip,1);
				cufftResult result = cufftExecR2C(cache->handle, real_data, (cufftComplex*)complex_data );
				cudaThreadSynchronize();
				print_cuda_fft_fail(result);
			}
			else {
				cache = get_cuda_fft_dd_plan_cache(rank,nx,ny,nz,real_2_complex,ip,1);
				/// CHECK LATER - Not sure if this will work
				result = cufftExecR2C(cache->handle, (cufftReal*)complex_data, (cufftComplex*)complex_data );
				print_cuda_fft_fail(result);
				cudaThreadSynchronize();

			}
		break;

		
		default:throw;
	}
	
	return 0;
}


int cuda_dd_fft_complex_to_real_nd(float *complex_data, float *real_data, int nx, int ny, int nz)
{
	device_init();
	const int rank = get_rank(ny, nz);
	bool ip;
	cufftResult result;
	//printf("using nx2 %d\n",nx2);
	cufft_dd_plan_cache* cache = 0;
	switch(rank) {
		case 1:
			cuda_dd_fft_complex_to_real_1d(complex_data, real_data, nx,1);
			break;
		
		case 2:
		case 3:
			ip = ( complex_data == real_data );
			cache = get_cuda_fft_dd_plan_cache(rank,nx,ny,nz,complex_2_real,ip,1);
			//printf("Executing handle %d\n",cache->handle);
			result = cufftExecC2R(cache->handle, (cufftComplex*)complex_data, real_data);
			print_cuda_fft_fail(result);
			cudaThreadSynchronize();
			break;
			
		default:throw;
	}
	return 0;
}
