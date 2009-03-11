

#include "cuda_util.h"
#include <stdio.h>
// Global texture
extern texture<float, 3, cudaReadModeElementType> tex;
extern texture<float, 2, cudaReadModeElementType> tex2d;

typedef unsigned int uint;

typedef unsigned int uint;
__global__ void mult_kernel(float *data,const float scale,const int num_threads)
{

	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	data[x+y*num_threads] *= scale;
}

void emdata_processor_mult( EMDataForCuda* cuda_data, const float& mult) {
	
	int max_threads = 512;

	int num_calcs = cuda_data->nx*cuda_data->ny*cuda_data->nz;
	
	int grid_y = num_calcs/max_threads;
	int res_y = num_calcs - (grid_y*max_threads);
	
	if ( grid_y > 0 ) {
		const dim3 blockSize(max_threads,1, 1);
		const dim3 gridSize(grid_y,1,1);
		mult_kernel<<<gridSize,blockSize>>>(cuda_data->data,mult,max_threads);
	}
	
	if ( res_y > 0 ) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		mult_kernel<<<gridSize,blockSize>>>(cuda_data->data+grid_y*max_threads,mult,0);
	}

	cudaThreadSynchronize();	
}

__global__ void correlation_kernel(float *ldata, float* rdata,const int num_threads)
{

	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	const uint idx = 2*x + y*num_threads;
	const uint idxp1 = idx+1;
	
	const float v1 = ldata[idx];
	const float v2 = ldata[idxp1];
	const float u1 = rdata[idx];
	const float u2 = rdata[idxp1];
	
	ldata[idx] = v1*u1 + v2*u2;
	ldata[idxp1] = v1*u2 - v2*u1;
}

__global__ void auto_correlation_kernel(float *ldata, float* rdata,const int num_threads)
{

	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	const uint idx = 2*x + y*num_threads;
	const uint idxp1 = idx+1;
	
	const float v1 = ldata[idx];
	const float v2 = ldata[idxp1];
	const float u1 = rdata[idx];
	const float u2 = rdata[idxp1];
	
	ldata[idx] = v1*u1 + v2*u2;
	ldata[idxp1] = 0;
}


__global__ void correlation_kernel_texture_2D(float *ldata,const int num_threads,const int xsize,const int offset)
{

	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	const uint idx = 2*x + y*num_threads;
	const uint idxp1 = idx+1;
	
	const uint tex_idx = offset+idx;
	const uint tx = tex_idx % xsize;
	const uint ty = tex_idx / xsize;
	
	const float v1 = ldata[idx];
	const float v2 = ldata[idxp1];
	const float u1 = tex2D(tex2d,tx,ty);
	const float u2 =  tex2D(tex2d,tx+1,ty);
	
	ldata[idx] = v1*u1 + v2*u2;
	ldata[idxp1] = v1*u2 - v2*u1;
}


__global__ void correlation_kernel_texture_3D(float *ldata,const int num_threads, const int xsize, const int xysize, const int offset)
{

	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	const uint idx = 2*x + y*num_threads;
	const uint idxp1 = idx+1;
	
	const uint tex_idx = offset+idx;
	const uint tx = tex_idx % xsize;
	const uint tz = tex_idx / xysize;
	const uint ty = (tex_idx - tz*xysize)/xsize;
	
	const float v1 = ldata[idx];
	const float v2 = ldata[idxp1];
	const float u1 = tex3D(tex,tx,ty,tz);
	const float u2 = tex3D(tex,tx+1,ty,tz);
	
	ldata[idx] = v1*u1 + v2*u2;
	ldata[idxp1] = v1*u2 - v2*u1;
}



void emdata_processor_correlation_texture( const EMDataForCuda* cuda_data) {
	int max_threads = 512; // I halve the threads because each kernel access memory in two locations

	int num_calcs = cuda_data->nx*cuda_data->ny*cuda_data->nz;
	
	int grid_y = num_calcs/(2*max_threads);
	int res_y = (num_calcs - (2*grid_y*max_threads))/2;
	
// 	printf("Grid %d, Res %d, dims %d %d %d\n",grid_y,res_y,cuda_data->nx,cuda_data->ny,cuda_data->nz);
	
	if ( grid_y > 0 ) {
		const dim3 blockSize(max_threads,1, 1);
		const dim3 gridSize(grid_y,1,1);
		if (cuda_data->nz == 1) {
			correlation_kernel_texture_2D<<<gridSize,blockSize>>>(cuda_data->data,2*max_threads,cuda_data->nx,0);
		} else {
			correlation_kernel_texture_3D<<<gridSize,blockSize>>>(cuda_data->data,2*max_threads,cuda_data->nx,cuda_data->nx*cuda_data->ny,0);
		}
	}
// 	res_y = 0;
	if ( res_y > 0 ) {
		const dim3 blockSize(res_y,1,1);
		const dim3 gridSize(1,1,1);
		int inc = 2*grid_y*max_threads;
// 		printf("Res %d, inc %d\n",res_y,inc);
		if (cuda_data->nz == 1) {
			correlation_kernel_texture_2D<<<gridSize,blockSize>>>(cuda_data->data+inc,0,cuda_data->nx,inc);
		} else {
			correlation_kernel_texture_3D<<<gridSize,blockSize>>>(cuda_data->data+inc,0,cuda_data->nx,cuda_data->nx*cuda_data->ny,inc);
		}
	}
	
	cudaThreadSynchronize();
}


void emdata_processor_correlation( const EMDataForCuda* left, const EMDataForCuda* right ) {
	int max_threads = 512;

	int num_calcs = left->nx*left->ny*left->nz;
	
	int grid_y = num_calcs/(2*max_threads);
	int res_y = (num_calcs - (2*grid_y*max_threads))/2;
	
	//printf("Grid y %d, res %d, dims %d %d %d\n", grid_y,res_y,left->nx,left->ny,left->nz);
	
	if ( grid_y > 0 ) {
		const dim3 blockSize(max_threads,1, 1);
		const dim3 gridSize(grid_y,1,1);
		if (left->data != right->data) {
			correlation_kernel<<<gridSize,blockSize>>>(left->data,right->data,2*max_threads);
		} else {
			auto_correlation_kernel<<<gridSize,blockSize>>>(left->data,right->data,2*max_threads);
		}
	}
	
	if ( res_y > 0 ) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		int inc = 2*grid_y*max_threads;
		if (left->data != right->data) {
			correlation_kernel<<<gridSize,blockSize>>>(left->data+inc,right->data+inc,0);
		} else {
			auto_correlation_kernel<<<gridSize,blockSize>>>(left->data+inc,right->data+inc,0);
		}
	}
	cudaThreadSynchronize();
}


