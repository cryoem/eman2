// Currently an empty file


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda.h>
#include "cuda_defs.h"
#include "cuda_util.h"

// Global texture
extern texture<float, 3, cudaReadModeElementType> tex;

typedef unsigned int uint;
__global__ void proj_kernel(float *out,float size,float3 mxx,float3 mxy, float3 mxz)
{
//	uint x=(threadIdx.x&0xfffe)+(blockIdx.x&1);
//	uint y=(blockIdx.x&0xfffe)+(threadIdx.x&1);
//	uint x=(threadIdx.x>>1)+((blockIdx.x&1)<<7);
//	uint y=(blockIdx.x&0xfffe)+(threadIdx.x&1);
	uint x=threadIdx.x;
	uint y=blockIdx.x;
	float fx=x-size/2;
	float fy=y-size/2;

	float tx,ty,tz;

	float sum=0;
	for (float fz=-size/2.0; fz<size/2.0-.1; fz+=1.0) {
		tx=fx*mxx.x+fy*mxx.y+fz*mxx.z+size/2.0;
		ty=fx*mxy.x+fy*mxy.y+fz*mxy.z+size/2.0;
		tz=fx*mxz.x+fy*mxz.y+fz*mxz.z+size/2.0;
		sum += CUDA_SAFE_CALL(tex3D(tex, tx,ty,tz));
	}

	out[x+y*(int)size]=sum;
}

void standard_project(const float* const matrix,const float* const rdata, const int nx, const int ny, const int nz, float*const d) 
{
	//device_init();
	
	//int idx = stored_cuda_array(rdata,nx,ny,nz);
	//bind_cuda_texture(idx);
	
	const dim3 blockSize(nx,1, 1);
	const dim3 gridSize(nx,1,1);

	float *memout=0;
	CUDA_SAFE_CALL(cudaMalloc((void **)&memout, nx*ny*sizeof(float)));
	
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
		
	proj_kernel<<<blockSize,gridSize>>>(memout,(float)nx,mxx,mxy,mxz);
	//CUDA_SAFE_CALL(cuCtxSynchronize());
	cudaThreadSynchronize();
	cudaMemcpy(d, memout, nx*ny*sizeof(float), cudaMemcpyDeviceToHost);
	cudaFree(memout);
}

