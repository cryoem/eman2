
#include <cuda.h>
#include <stdio.h>

#include "cuda_defs.h"
#include "cuda_emfft.h"
#include "cuda_util.h"

#define MAX_THREADS 128

texture<float, 3, cudaReadModeElementType> texA;
texture<float, 2, cudaReadModeElementType> texA2d;
texture<float, 3, cudaReadModeElementType> texB;
texture<float, 2, cudaReadModeElementType> texB2d;

#include "cuda_processor.cu"
#include "cuda_cmp.cu"
#include "cuda_projector.cu"
#include "cuda_reconstructor.cu"

void cuda_bind_texture_3d(texture<float, 3, cudaReadModeElementType> &texture, const cudaArray * const array, const bool interp_mode) {
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
	texture.addressMode[0] = cudaAddressModeClamp;
	texture.addressMode[1] = cudaAddressModeClamp;

//	printf("Bound 2D texture to array %x\n", array);
	cudaBindTextureToArray(texture, array);

}

cudaArray* get_cuda_array(const int nx, const int ny, const int nz)
{
	cudaArray *array = 0;
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	
	if (nz > 1) {
		cudaExtent VS = make_cudaExtent(nx,ny,nz);
		cudaError_t error = cudaMalloc3DArray(&array, &channelDesc, VS);
		if ( error != cudaSuccess) {
			printf("Could not allocate array\n");
			return 0;
		}
	} else if ( ny > 1) {
		cudaError_t error = cudaMallocArray(&array,&channelDesc,nx,ny);
		if ( error != cudaSuccess) {
			printf("Could not allocate array\n");
			return 0;
		}
	}

	return array;
}

bool copy_to_array(const float * data, cudaArray * array, const int nx, const int ny, const int nz, const cudaMemcpyKind memkind)
{

	if (nz > 1) {
		cudaExtent VS = make_cudaExtent(nx,ny,nz);
		cudaMemcpy3DParms copyParams = {0};
		copyParams.srcPtr   = make_cudaPitchedPtr((void*)data, VS.width*sizeof(float), VS.width, VS.height);
		copyParams.dstArray = array;
		copyParams.extent   = VS;
		copyParams.kind     = memkind;
		cudaError_t error =  cudaMemcpy3D(&copyParams);
		if ( error != cudaSuccess) {
			const char* e = cudaGetErrorString(error);
			printf("CUDA error from cudaMemcpy3D: %s\n",e);
			return 0;	
		}
	} else if ( ny > 1) {
		cudaError_t error = cudaMemcpyToArray(array, 0, 0, data, nx*ny*nz*sizeof(float), memkind);
		if ( error != cudaSuccess)
		{
			const char* e = cudaGetErrorString(error);
			printf("CUDA error from cudaMemcpyToArray: %s\n",e);
			return 0;	
		}
	}
	return 1;
}

void bind_cuda_array_to_textureA( const cudaArray* const array, const int ndims,const bool interp_mode) {
	
	if (ndims == 3) {
		cuda_bind_texture_3d(texA,array,interp_mode);
	} else {
		cuda_bind_texture_2d(texA2d,array,interp_mode);
	} 
	
}

void unbind_cuda_textureA(const int ndims) {
	if (ndims == 3) {
		cudaUnbindTexture(&texA);
	}else {
		cudaUnbindTexture(&texA2d);
	}
}

void bind_cuda_array_to_textureB( const cudaArray* const array, const int ndims,const bool interp_mode) {
	
	if (ndims == 3) {
		cuda_bind_texture_3d(texB,array,interp_mode);
	} else {
		cuda_bind_texture_2d(texB2d,array,interp_mode);
	} 
	
}

void unbind_cuda_textureB(const int ndims) {
	if (ndims == 3) {
		cudaUnbindTexture(&texB);
	}else {
		cudaUnbindTexture(&texB2d);
	}
}

bool device_init() {
	static bool init = true;
	
	if (init) {
		int deviceCount;
		cudaGetDeviceCount(&deviceCount);
		
		if (deviceCount == 0){
			printf("WARNING NO CUDA DEVICES FOUND, NOT USING CUDA\n");
			return 0;
		}
		
 		if (deviceCount > 1) {
 			printf("%d CUDA devices detected\n",deviceCount);
 		} else { // must be one
 			printf("1 CUDA device detected\n");
 		}
		
		cudaDeviceProp deviceProp;
		cudaGetDeviceProperties(&deviceProp, 0);
		if (deviceProp.major < 1) exit(2);
		
		int device = (getenv("EMANUSECUDA") == NULL || atoi(getenv("EMANUSECUDA")) == 0) ? 0 : atoi(getenv("EMANUSECUDA")) - 1;	
		cudaError_t cudareturn = cudaSetDevice(device); 
		if(cudareturn != cudaSuccess) {
			printf("\nERROR in cudaSetDevice.... %s\n", cudaGetErrorString(cudareturn));
			exit(2);
		} else {
			int curdev;
			cudaGetDevice(&curdev);
			printf("Using CUDA device %d\n", curdev);
		}

		init = false; //Force init everytime
	}
	return 1;
}

__global__ void get_edgemean_kernal(const float* data, float* edgemean, const int nx, const int ny, const int nz)
{
	int di = 0;
	float edge_sum = 0;
	float edge_mean = 0;
	size_t nxy = nx * ny;
	if (nz == 1) {
		for (int i = 0, j = (ny - 1) * nx; i < nx; ++i, ++j) {
			edge_sum += data[i] + data[j];
		}
		for (size_t i = 0, j = nx - 1; i < nxy; i += nx, j += nx) {
			edge_sum += data[i] + data[j];
		}
		edge_mean = (float)edge_sum / (nx * 2 + ny * 2);
	}
	else {
		if (nx == ny && nx == nz * 2 - 1) {
			for (size_t j = (nxy * (nz - 1)); j < nxy * nz; ++j, ++di) {
				edge_sum += data[j];
			}
		}
		else {
			for (size_t i = 0, j = (nxy * (nz - 1)); i < nxy; ++i, ++j, ++di) {
				edge_sum += data[i] + data[j];
			}
		}

		int nxy2 = nx * (ny - 1);
		for (int k = 1; k < nz - 1; ++k) {
			size_t k2 = k * nxy;
			size_t k3 = k2 + nxy2;
			for (int i = 0; i < nx; ++i, ++di) {
				edge_sum += data[i + k2] + data[i + k3];
			}
		}
		for (int k = 1; k < nz - 1; ++k) {
			size_t k2 = k * nxy;
			size_t k3 = nx - 1 + k2;
			for (int i = 1; i < ny - 1; ++i, ++di) {
				edge_sum += data[i * nx + k2] + data[i * nx + k3];
			}
		}

		edge_mean = (float)edge_sum / (di * 2);
	}
	*edgemean = edge_mean;
} 
float get_edgemean_cuda(const float* data, const int nx, const int ny, const int nz)
{

	const dim3 blockSize(1,1,1);
	const dim3 gridSize(1,1,1);

	float * d_edgemean=0;
	cudaMalloc((void **)&d_edgemean, sizeof(float));
	float * h_edgemean = 0;
	h_edgemean = (float*) malloc(sizeof(float));

	get_edgemean_kernal<<<gridSize,blockSize>>>(data,d_edgemean, nx, ny, nz);
	cudaThreadSynchronize();
	cudaMemcpy(h_edgemean,d_edgemean,sizeof(float),cudaMemcpyDeviceToHost);
	cudaFree(d_edgemean);
	
	float result = *h_edgemean;
	free(h_edgemean);
	
	return result;
}

__global__ void tovalue_kernal(float* data, const float value, const int num_threads, const int offset)
{

	const uint idx = threadIdx.x + blockIdx.x*num_threads + offset;

	data[idx] = value;

}

void to_value_cuda(float* data, const float value, const int nx, const int ny, const int nz)
{

	int num_calcs = ny*nx*nz;

	int grid_y = num_calcs/(MAX_THREADS);
	int res_y = num_calcs - grid_y*MAX_THREADS;

	if (grid_y > 0) {
		const dim3 blockSize(MAX_THREADS,1, 1);
		const dim3 gridSize(grid_y,1,1);
		tovalue_kernal<<<gridSize,blockSize>>>(data, value, MAX_THREADS,0);
	}
	if (res_y) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		int offset = grid_y*MAX_THREADS;
		tovalue_kernal<<<gridSize,blockSize>>>(data, value, MAX_THREADS,offset);
	}
	cudaThreadSynchronize();

	return;
}

void to_zero_cuda(float* data, const int nx, const int ny, const int nz)
{
	to_value_cuda(data, 0.0, nx, ny, nz);
}