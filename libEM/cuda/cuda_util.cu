
#include <cuda.h>
#include <stdio.h>

#include "cuda_defs.h"
#include "cuda_emfft.h"
#include "cuda_util.h"

texture<float, 3, cudaReadModeElementType> tex;
texture<float, 2, cudaReadModeElementType> tex2d;

void cuda_bind_texture_3d(texture<float, 3, cudaReadModeElementType> &texture,const cudaArray * const array, const bool interp_mode) {
	texture.normalized = 0;
	if (interp_mode) texture.filterMode = cudaFilterModeLinear;
	else texture.filterMode = cudaFilterModePoint;
	texture.addressMode[0] = cudaAddressModeClamp;
	texture.addressMode[1] = cudaAddressModeClamp;
	texture.addressMode[2] = cudaAddressModeClamp;
	
	cudaBindTextureToArray(texture, array);
}

void cuda_bind_texture_2d(texture<float, 2, cudaReadModeElementType> &texture, const cudaArray * const array, const bool interp_mode) {
	texture.normalized = 0;
	if (interp_mode) texture.filterMode = cudaFilterModeLinear;
	else texture.filterMode = cudaFilterModePoint;
	// tex.filterMode = cudaFilterModePoint;
	texture.addressMode[0] = cudaAddressModeClamp;
	texture.addressMode[1] = cudaAddressModeClamp;
// 	tex.addressMode[2] = cudaAddressModeClamp;
	
	cudaBindTextureToArray(texture, array);
}

cudaArray* get_cuda_array(const float * const data,const int nx, const int ny, const int nz, const cudaMemcpyKind mem_cpy_flag)
{
	cudaArray *array = 0;
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	
	if (nz > 1) {
		cudaExtent VS = make_cudaExtent(nx,ny,nz);
// 		printf("It's a 3D one %d %d %d %d\n",VS.width,VS.height,nx,ny);
		cudaMalloc3DArray(&array, &channelDesc, VS);
// 		printf("It's a 3D one %d %d %d %d %d\n",VS.width,data,nx,ny,nz);
		cudaMemcpy3DParms copyParams = {0};
		copyParams.srcPtr   = make_cudaPitchedPtr((void*)data, VS.width*sizeof(float), VS.width, VS.height);
		copyParams.dstArray = array;
		copyParams.extent   = VS;
		copyParams.kind     = mem_cpy_flag;
		cudaError_t error =  cudaMemcpy3D(&copyParams);
		if ( error != cudaSuccess) {
			const char* e = cudaGetErrorString(error);
			printf("CUDA error from cudaMemcpy3D: %s\n",e);
			cudaFreeArray(array);
			return 0;	
		}
	} else if ( ny > 1) {
// 		printf("It's a 2D one\n");d
		cudaMallocArray(&array,&channelDesc,nx,ny);
		cudaExtent VS = make_cudaExtent(nx,ny,nz);
// 		printf("It's a 3D one %d %d %d %d\n",VS.width,VS.height,nx,ny);
		cudaMalloc3DArray(&array, &channelDesc, VS);
		cudaError_t error = cudaMemcpyToArray(array, 0, 0, data, nx*ny*nz*sizeof(float), mem_cpy_flag);
		if ( error != cudaSuccess)
		{
			const char* e = cudaGetErrorString(error);
			printf("CUDA error from cudaMemcpyToArray: %s\n",e);
			cudaFreeArray(array);
			return 0;	
		}
	} else throw;
	
	return array;
}


cudaArray* get_cuda_array_device(const float * const data,const int nx, const int ny, const int nz)
{
	return get_cuda_array(data,nx,ny,nz,cudaMemcpyDeviceToDevice);
}

cudaArray* get_cuda_array_host(const float * const data,const int nx, const int ny, const int nz)
{
	return get_cuda_array(data,nx,ny,nz,cudaMemcpyHostToDevice);
}

void bind_cuda_array_to_texture( const cudaArray* const array, const int ndims,const bool interp_mode) {
	if (ndims == 3) {
		//printf("Binding 3D texture\n");
		cuda_bind_texture_3d(tex,array,interp_mode);
	} else if ( ndims == 2) {
		cuda_bind_texture_2d(tex2d,array,interp_mode);
	} else throw;
	//printf("Done bind\n");
}

void unbind_cuda_texture(const int ndims) {
	if (ndims == 3) {
		cudaUnbindTexture(&tex);
	}else if ( ndims == 2) {
		cudaUnbindTexture(&tex2d);
	} else throw;
}

void device_init() {
	static bool init = true;
	
	if (init) {
		int deviceCount;
		cudaGetDeviceCount(&deviceCount);
		
		if (deviceCount == 0) throw;
		
// 		if (deviceCount > 1) {
// 			printf("%d CUDA devices detected\n",deviceCount);
// 		} else { // must be one
// 			printf("1 CUDA device detected\n");
// 		}
		
		cudaDeviceProp deviceProp;
		cudaGetDeviceProperties(&deviceProp, 0);
		if (deviceProp.major < 1) exit(2);
		
		cudaSetDevice(0);
		init_cuda_fft_hh_plan_cache(); // should only be performed if the host is using Cuda ffts, which is unlikey. Will do after development has progressed.
		init_cuda_fft_dd_plan_cache();
		init = false; //Force init everytikme
	}
}

__global__ void  calc_max_location_wrap(int* const soln, const float* data,const int maxdx, const int maxdy, const int maxdz, const int nx, const int ny, const int nz) {
	int maxshiftx = maxdx, maxshifty = maxdy, maxshiftz = maxdz;
	if (maxdx == -1) maxshiftx = nx/4;
	if (maxdy == -1) maxshifty = ny/4;
	if (maxdz == -1) maxshiftz = nz/4;

	float max_value = -10000000000000;
	int nxy = nx*ny;

	for (int k = -maxshiftz; k <= maxshiftz; k++) {
		for (int j = -maxshifty; j <= maxshifty; j++) {
			for (int i = -maxshiftx; i <= maxshiftx; i++) {
				
				int kk = k;
				if (kk < 0) {
					kk = nz+kk;
				}
				int jj = j;
				if (jj < 0) {
					jj = ny+jj;
				}
				
				int ii = i;
				if (ii < 0) {
					ii = nx+ii;
				}
				float value = data[ii+jj*nx+kk*nxy];

				if (value > max_value) {
					max_value = value;
					soln[0] = i;
					soln[1] = j;
					soln[2] = k;
				}
			}
		}
	}
}

int* calc_max_location_wrap_cuda(const EMDataForCuda* data, const int maxdx, const int maxdy, const int maxdz) {
	
	int * device_soln=0;
	cudaError_t error = cudaMalloc((void **)&device_soln, 3*sizeof(int));
	if ( error != cudaSuccess ){
		printf("Cuda malloc failed in calc_max_location_wrap_cuda");
		return 0;	
	}
		
	int * host_soln = 0;
	host_soln = (int*) malloc(3*sizeof(int));
	
	const dim3 blockSize(1,1, 1);
	const dim3 gridSize(1,1,1);
	
	calc_max_location_wrap<<<blockSize,gridSize>>>(device_soln,data->data,maxdx,maxdy,maxdz,data->nx,data->ny,data->nz);
	cudaThreadSynchronize();
	cudaMemcpy(host_soln,device_soln,3*sizeof(int),cudaMemcpyDeviceToHost);
	cudaFree(device_soln);
	return host_soln;
}


__global__ void cut_slice_kernel(float *out,int size, float3 mxx,float3 mxy, float3 mxz, float3 trans)
{
	unsigned int x=threadIdx.x;
	unsigned int y=blockIdx.x;
	
	float fx=x-size/2.0;
	float fy=y-size/2.0;

	// The 0.5f offsets for x,y and z are required - Read section D.2 in Appendix D of the CUDA
	// Programming Guide (version 2.0).
	// Thankyou http://sites.google.com/site/cudaiap2009/cookbook-1
	float tx=fx*mxx.x+fy*mxx.y+size/2.0+trans.x+0.5;
	float ty=fx*mxy.x+fy*mxy.y+size/2.0+trans.y+0.5;
	float tz=fx*mxz.x+fy*mxz.y+size/2.0+trans.z+0.5;

	out[x+y*(int)size]=tex3D(tex, tx,ty,tz);
}



void cut_slice_cuda_(const EMDataForCuda* to_data,const float* const matrix)
{
	
	const dim3 blockSize(to_data->ny,1, 1);
	const dim3 gridSize(to_data->nx,1,1);
	
	float3 mxx,mxy,mxz,trans;
	
	mxx.x=matrix[0];
	mxx.y=matrix[1];
	mxx.z=matrix[2];
	mxy.x=matrix[4];
	mxy.y=matrix[5];
	mxy.z=matrix[6];
	mxz.x=matrix[8];
	mxz.y=matrix[9];
	mxz.z=matrix[10];
	trans.x =matrix[3];
	trans.y =matrix[7];
	trans.z =matrix[11];
	
		
	cut_slice_kernel<<<blockSize,gridSize>>>(to_data->data,to_data->nx,mxx,mxy,mxz,trans);
	//CUDA_SAFE_CALL(cuCtxSynchronize());
	cudaThreadSynchronize();
// 	cudaMemcpy(d, memout, nx*ny*sizeof(float), cudaMemcpyDeviceToHost);
// 	cudaFree(memout);
	
}
typedef unsigned int uint;
__global__ void column_sum(float* sum, int ny, int num_threads, int offset ) {
	
	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	float s = 0.0;
	const uint idx_x = x + y*num_threads+offset; /* This is always an x index */
	for(int i =0; i < ny; ++i) {
		s += tex2D(tex2d,idx_x,i);
	}
	sum[idx_x] = s; 
}

void emdata_column_sum(const EMDataForCuda* sum_target,const int ny) {
	int max_threads = 192;
	if (max_threads > sum_target->nx) max_threads = sum_target->nx;
	
	int num_calcs = sum_target->nx;
		
	int grid_y = num_calcs/(max_threads);
	int res_y = (num_calcs - (grid_y*max_threads));
	
	if ( grid_y > 0 ) {
		const dim3 blockSize(max_threads,1, 1);
		const dim3 gridSize(grid_y,1,1);
		column_sum<<<gridSize,blockSize>>>(sum_target->data,ny,max_threads,0);
	}
	
	if ( res_y > 0 ) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		column_sum<<<gridSize,blockSize>>>(sum_target->data,ny,max_threads,grid_y*max_threads);
	}
	cudaThreadSynchronize();
}
