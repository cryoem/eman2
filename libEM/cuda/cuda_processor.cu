

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

	const uint idx = 2*x + y*num_threads+offset;
	const uint idxp1 = idx+1;
	
	const uint tx = idx % xsize;
	const uint ty = idx / xsize;
	
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

	const uint idx = 2*x + y*num_threads + offset;
	const uint idxp1 = idx+1;
	
	const uint tx = idx % xsize;
	const uint tz = idx / xysize;
	const uint ty = (idx - tz*xysize)/xsize;
	
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
			correlation_kernel_texture_2D<<<gridSize,blockSize>>>(cuda_data->data,0,cuda_data->nx,inc);
		} else {
			correlation_kernel_texture_3D<<<gridSize,blockSize>>>(cuda_data->data,0,cuda_data->nx,cuda_data->nx*cuda_data->ny,inc);
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

__global__ void unwrap_kernel(float* dptr, const int num_threads, const int r1, const float p, const int nx, const int ny, const int nxp, const int offset) {
	const uint x=threadIdx.x;
	const uint y=blockIdx.x;
	
	const uint idx = x + y*num_threads+offset;
	
	const uint tx = idx % nxp;
	const uint ty = idx / nxp;
	
	float si = sinf(tx * 3.141592653589793 * p );
	float co = cosf(tx * 3.141592653589793 * p );

	float xx = (ty + r1) * co + nx / 2;
	float yy = (ty + r1) * si + ny / 2;
	dptr[idx] = tex2D(tex2d,xx+0.5,yy+0.5)*(ty+r1);
}

EMDataForCuda* emdata_unwrap(int r1, int r2, int xs, int do360,int nx, int ny) {
	int p = 1;
	if (do360) {
		p = 2;
	}

	if (xs < 1) {
		//xs = (int) Util::fast_floor(p * M_PI * ny / 4);
		xs = (int) (p * 3.1415926535897931 * ny / 4);
		xs -= xs % 8;
		if (xs<=8) xs=16;
	}

	if (r1 < 0) {
		r1 = 4;
	}

	int rr = ny / 2 - 2;
	
	rr-=rr%2;
	if (r2 <= r1 || r2 > rr) {
		r2 = rr;
	}

	if ( (r2-r1) < 0 ) throw;
	float* dptr;
	int n = xs*(r2-r1);
	cudaError_t error = cudaMalloc((void**)&dptr,n*sizeof(float));
	if ( error != cudaSuccess ) {
		const char* s = cudaGetErrorString(error);
		printf("Cuda malloc failed in emdata_unwrap: %s\n",s);
		throw;
	}
	
	int max_threads = 512;
	int num_calcs = n;
	
	int grid_y = num_calcs/(max_threads);
	int res_y = num_calcs - grid_y*max_threads;
	
	//printf("Grid %d, res %d, n %d, p %f \n",grid_y,res_y,n, p/xs);
	
	if ( grid_y > 0 ) {
		const dim3 blockSize(max_threads,1, 1);
		const dim3 gridSize(grid_y,1,1);
		unwrap_kernel<<<gridSize,blockSize>>>(dptr,max_threads,r1,(float) p/ (float)xs, nx,ny,xs,0);	
	}
	
	if ( res_y > 0 ) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		unwrap_kernel<<<gridSize,blockSize>>>(dptr,max_threads,r1, (float) p/ (float)xs, nx,ny,xs,grid_y*max_threads);	
	}
	
	EMDataForCuda* tmp = (EMDataForCuda*) malloc( sizeof(EMDataForCuda) );
	tmp->data = dptr;
	tmp->nx = xs;
	tmp->ny = r2-r1;
	tmp->nz = 1;
	return tmp;
}
