// Currently an empty file


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda.h>
#include "cuda_defs.h"
#include "cuda_util.h"

// Global texture
extern texture<float, 3, cudaReadModeElementType> tex3d;

typedef unsigned int uint;
__global__ void proj_kernel(float *out,float size, float size_on_two, float3 mxx,float3 mxy, float3 mxz)
{
//	uint x=(threadIdx.x&0xfffe)+(blockIdx.x&1);
//	uint y=(blockIdx.x&0xfffe)+(threadIdx.x&1);
//	uint x=(threadIdx.x>>1)+((blockIdx.x&1)<<7);
//	uint y=(blockIdx.x&0xfffe)+(threadIdx.x&1);
	uint x=threadIdx.x;
	uint y=blockIdx.x;
	
	float fx=x-size_on_two;
	float fy=y-size_on_two;

	float tx,ty,tz;

	float sum=0;
	// The 0.5f offsets for x,y and z are required - Read section D.2 in Appendix D of the CUDA
	// Programming Guide (version 2.0).
	// Thankyou http://sites.google.com/site/cudaiap2009/cookbook-1
	for (float fz=-size_on_two; fz<size_on_two-.1; fz+=1.0) {
		tx=fx*mxx.x+fy*mxx.y+fz*mxx.z+size_on_two+0.5;
		ty=fx*mxy.x+fy*mxy.y+fz*mxy.z+size_on_two+0.5;
		tz=fx*mxz.x+fy*mxz.y+fz*mxz.z+size_on_two+0.5;
		sum += tex3D(tex3d, tx,ty,tz);
	}

	out[x+y*(int)size]=sum;
}

void standard_project(const float* const matrix,const EMDataForCuda* const data) 
{
	//device_init();
	
	//int idx = stored_cuda_array(rdata,nx,ny,nz);
	//bind_cuda_texture(idx);
	
	const dim3 blockSize(data->ny,1, 1);
	const dim3 gridSize(data->nx,1,1);
	
	float3 mxx,mxy,mxz;
	
	mxx.x=matrix[0];
	mxx.y=matrix[4];
	mxx.z=matrix[8];
	mxy.x=matrix[1];
	mxy.y=matrix[5];
	mxy.z=matrix[9];
	mxz.x=matrix[2];
	mxz.y=matrix[6];
	mxz.z=matrix[10];
		
	proj_kernel<<<blockSize,gridSize>>>(data->data,(float)data->nx,(float)data->nx/2,mxx,mxy,mxz);
	//CUDA_SAFE_CALL(cuCtxSynchronize());
	cudaThreadSynchronize();
}

