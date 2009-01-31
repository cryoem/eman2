
#include <cuda.h>
#include <stdio.h>

#include "cuda_defs.h"
#include "cuda_emfft.h"
texture<float, 3, cudaReadModeElementType> tex;

void cudaBindTexture(texture<float, 3, cudaReadModeElementType> &tex,cudaArray *array) {
	tex.normalized = false;
	tex.filterMode = cudaFilterModeLinear;
	tex.addressMode[0] = cudaAddressModeClamp;
	tex.addressMode[1] = cudaAddressModeClamp;
	tex.addressMode[2] = cudaAddressModeClamp;
	
	CUDA_SAFE_CALL(cudaBindTextureToArray(tex, array));
}


struct CudaEMDataArray {
	cudaArray* array;
	const float* data;
};

const int max_cuda_arrays = 100;
int num_cuda_arrays = 0;
CudaEMDataArray cuda_arrays[max_cuda_arrays];

void init_cuda_emdata_arrays() {
	for(int i = 0; i < max_cuda_arrays; ++i ) {
		CudaEMDataArray c =  { 0, 0 };
		cuda_arrays[i] = c;
	}
}

int stored_cuda_array(const float * data,const int nx, const int ny, const int nz) {
	
	for(int i = 0; i < num_cuda_arrays; ++i ) {
		if (cuda_arrays[i].data == data ) {
			//printf("Found that cuda arrary\n");
			return i;
		}
	}
	
	//printf("Making a new cuda array\n");
	// If we make it here then it doesn't exist
	cuda_arrays[num_cuda_arrays].data = data;
	cudaExtent VS = make_cudaExtent(nx,ny,nz);
	
	cudaArray *array = 0;
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	CUDA_SAFE_CALL(cudaMalloc3DArray(&array, &channelDesc, VS));
	
	cudaMemcpy3DParms copyParams = {0};
	copyParams.srcPtr   = make_cudaPitchedPtr((void*)data, VS.width*sizeof(float), VS.width, VS.height);
	copyParams.dstArray = array;
	copyParams.extent   = VS;
	copyParams.kind     = cudaMemcpyHostToDevice;
	CUDA_SAFE_CALL(cudaMemcpy3D(&copyParams));
	
	//cudaBindTexture(tex,array);
	cuda_arrays[num_cuda_arrays].array = array;
	num_cuda_arrays++;
	return num_cuda_arrays-1;
}

int delete_cuda_array(const int idx) {
	CUDA_SAFE_CALL(cudaFree(cuda_arrays[idx].array));
	cuda_arrays[idx].array = 0;
	/*
	for (int i = idx; i < num_cuda_arrays-1; ++i ) {
		
	}
	*/
	return 0;
}

int delete_cuda_memory(float*p) {
	CUDA_SAFE_CALL(cudaFree(p));
	p = 0;
	return 0;
}

void bind_cuda_texture(const int idx) {
	//printf("Binding texture\n");
	CUDA_SAFE_CALL(cudaBindTexture(tex,cuda_arrays[idx].array));
	//printf("Done bind\n");
}

void device_init() {
	static bool init = true;
	
	if (init) {
		int deviceCount;
		CUDA_SAFE_CALL(cudaGetDeviceCount(&deviceCount));
		printf("%d CUDA devices detected\n",deviceCount);
		if (deviceCount == 0) exit(1);
		
		cudaDeviceProp deviceProp;
		CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, 0));
		if (deviceProp.major < 1) exit(2);
		
		CUDA_SAFE_CALL(cudaSetDevice(0));
		init_cuda_emdata_arrays();
#ifdef FFTW_PLAN_CACHING
		init_cuda_emfft_cache();
#endif //FFTW_PLAN_CACHING
		init = false; //Force init everytikme
	}
}