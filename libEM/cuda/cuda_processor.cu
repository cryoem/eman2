

#include "cuda_util.h"



typedef unsigned int uint;
__global__ void mult_kernel(float *data,const float scale, const int z, const int xsize, const int xysize)
{

	uint x=threadIdx.x;
	uint y=blockIdx.x;


	data[x+y*xsize+z*xysize] *= scale;
}


void emdata_processor_mult( const EMDataForCuda* cuda_data, const float& mult) {
	
	const dim3 blockSize(cuda_data->ny,1, 1);
	const dim3 gridSize(cuda_data->nx,1,1);
		
	for (int i = 0; i < cuda_data->nz; ++i) {
		mult_kernel<<<blockSize,gridSize>>>(cuda_data->data,mult,i,cuda_data->nx,cuda_data->nx*cuda_data->ny);
	}
	//CUDA_SAFE_CALL(cuCtxSynchronize());
	cudaThreadSynchronize();	
}



__global__ void correlation_kernel(float *ldata, float* rdata, const int z,const int xsize, const int xysize)
{

	uint x=threadIdx.x;
	uint y=blockIdx.x;

	uint idx = 2*x + y*xsize+z*xysize;
	
	float v1 = ldata[idx];
	float v2 = ldata[idx+1];
	float u1 = rdata[idx];
	float u2 = rdata[idx+1];
	
	ldata[idx] = v1*u1 - v2*u2;
	ldata[idx+1] = v1*u2 + v2*u1;
}


void emdata_processor_correlation( const EMDataForCuda* left,const EMDataForCuda* right) {
	const dim3 blockSize(left->ny,1, 1);
	const dim3 gridSize(left->nx/2,1,1);
		
	for (int i = 0; i < left->nz; ++i) {
		correlation_kernel<<<blockSize,gridSize>>>(left->data,right->data,i,left->nx,left->nx*left->ny);
	}
	cudaThreadSynchronize();
}