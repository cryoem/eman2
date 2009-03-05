

#include "cuda_util.h"

// Global texture
extern texture<float, 3, cudaReadModeElementType> tex;
extern texture<float, 2, cudaReadModeElementType> tex2d;

typedef unsigned int uint;
__global__ void mult_kernel(float *data,const float scale, const int z, const int xsize, const int xysize)
{

	uint x=threadIdx.x;
	uint y=blockIdx.x;


	data[x+y*xsize+z*xysize] *= scale;
}

__global__ void mult_kernel_offset(float *data,const float scale, const int z, const int xsize, const int xysize, const int offset)
{

	uint x=threadIdx.x;
	uint y=blockIdx.x;


	data[x+offset+y*xsize+z*xysize] *= scale;
}



void emdata_processor_mult( const EMDataForCuda* cuda_data, const float& mult) {
	
	
// // // 	if ( cuda_data->nx <= 256 ) {
// // // 		const dim3 blockSize(2*cuda_data->nx,1, 1);
// // // 		const dim3 gridSize(cuda_data->ny/2,1,1);
// // // 			
// // // 		for (int i = 0; i < cuda_data->nz; ++i) {
// // // 			mult_kernel<<<gridSize,blockSize>>>(cuda_data->data,mult,i,2*cuda_data->nx,cuda_data->nx*cuda_data->ny);
// // // 		}
// // // 	}
// // // 	else {
	if (cuda_data->nx <= 512) {
		const dim3 blockSize(cuda_data->nx,1, 1);
		const dim3 gridSize(cuda_data->ny,1,1);
			
		for (int i = 0; i < cuda_data->nz; ++i) {
			mult_kernel<<<gridSize,blockSize>>>(cuda_data->data,mult,i,cuda_data->nx,cuda_data->nx*cuda_data->ny);
		}
	} else {
		int offset = 0; 
		while(offset < cuda_data->nx) {
			int block_size = 512;
			if ( (block_size + offset) > cuda_data->nx ) {
				block_size = cuda_data->nx - offset;
			}
			const dim3 blockSize(block_size,1, 1);
			const dim3 gridSize(cuda_data->ny,1,1);
			
			for (int i = 0; i < cuda_data->nz; ++i) {
				mult_kernel_offset<<<gridSize,blockSize>>>(cuda_data->data,mult,i,cuda_data->nx,cuda_data->nx*cuda_data->ny,offset);
			}
			offset += 512;
		}
	}
		//CUDA_SAFE_CALL(cuCtxSynchronize());
	cudaThreadSynchronize();	
}



__global__ void correlation_kernel_2D(float *ldata, float* rdata,const int xsize)
{

	uint x=threadIdx.x;
	uint y=blockIdx.x;

	uint idx = 2*x + y*xsize;
	uint idxp1 = idx+1;
	
	float v1 = ldata[idx];
	float v2 = ldata[idxp1];
	float u1 = rdata[idx];
	float u2 = rdata[idxp1];
	
	ldata[idx] = v1*u1 + v2*u2;
	ldata[idxp1] = v1*u2 - v2*u1;
}

__global__ void correlation_kernel_3D(float *ldata, const int z,const int xsize, const int xysize)
{

	uint x=threadIdx.x;
	uint y=blockIdx.x;

	uint twox = 2*x;
	uint idx = twox + y*xsize+z*xysize;
	uint idxp1 = idx+1;
	
	float v1 = ldata[idx];
	float v2 = ldata[idxp1];
	float u1 = tex3D(tex,twox,y,z);
	float u2 = tex3D(tex,twox+1,y,z);
	
	ldata[idx] = v1*u1 + v2*u2;
	ldata[idxp1] = v1*u2 - v2*u1;
}



__global__ void correlation_kernel_2D_2(float *ldata,const int xsize)
{

	uint x=threadIdx.x;
	uint y=blockIdx.x;

	uint twox = 2*x;
	uint idx = twox + y*xsize;
	uint idxp1 = idx+1;
	
	float v1 = ldata[idx];
	float v2 = ldata[idxp1];
	float u1 = tex2D(tex2d,twox,y);
	float u2 =  tex2D(tex2d,twox+1,y);
	
	ldata[idx] = v1*u1 + v2*u2;
	ldata[idxp1] = v1*u2 - v2*u1;
}


__global__ void correlation_kernel_2D_2_offset(float *ldata,const int xsize, const int offset)
{

	uint x=threadIdx.x;
	uint y=blockIdx.x;

	uint twox = 2*x + offset;
	uint idx = twox + y*xsize;
	uint idxp1 = idx+1;
	
	float v1 = ldata[idx];
	float v2 = ldata[idxp1];
	float u1 = tex2D(tex2d,twox,y);
	float u2 =  tex2D(tex2d,twox+1,y);
	
	ldata[idx] = v1*u1 + v2*u2;
	ldata[idxp1] = v1*u2 - v2*u1;
}

void emdata_processor_correlation( const EMDataForCuda* left) {
	const dim3 gridSize(left->ny,1, 1);
	const dim3 blockSize(left->nx/2,1,1);
	int nz = left->nz;
	if (nz > 1) {
		for (int i = 0; i < left->nz; ++i) {
			correlation_kernel_3D<<<gridSize,blockSize>>>(left->data,i,left->nx,left->nx*left->ny);
		}
	}
	else {
		if (left->nx <= 512) {
			correlation_kernel_2D_2<<<gridSize,blockSize>>>(left->data,left->nx);
		} else {
			int offset = 0; 
			while(offset < left->nx) {
				int block_size = 512;
				if ( (block_size + offset) > left->nx ) {
					block_size = left->nx - offset;
				}
				const dim3 gridSize(left->ny,1, 1);
				const dim3 blockSize(block_size/2,1,1);
				
				for (int i = 0; i < left->nz; ++i) {
					correlation_kernel_2D_2_offset<<<gridSize,blockSize>>>(left->data,left->nx,offset);
				}
				offset += 512;
			}
			
		}
		//correlation_kernel_2D<<<blockSize,gridSize>>>(left->data,right->data,left->nx);
	}	
	cudaThreadSynchronize();
}