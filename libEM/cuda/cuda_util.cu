
#include <cuda.h>
#include <stdio.h>

texture<float, 3, cudaReadModeElementType> tex;

void cudaBindTexture(texture<float, 3, cudaReadModeElementType> &tex,cudaArray *array) {
	

	tex.normalized = false;
	tex.filterMode = cudaFilterModeLinear;
	tex.addressMode[0] = cudaAddressModeClamp;
	tex.addressMode[1] = cudaAddressModeClamp;
	tex.addressMode[2] = cudaAddressModeClamp;
	
	cudaBindTextureToArray(tex, array);
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
			printf("Found that cuda arrary\n");
			//cudaBindTexture(tex,cuda_arrays[i].array);
			return i;
		}
	}
	
	printf("Making a new cuda array\n");
	// If we make it here then it doesn't exist
	cuda_arrays[num_cuda_arrays].data = data;
	cudaExtent VS = make_cudaExtent(nx,ny,nz);
	
	cudaArray *array = 0;
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	cudaMalloc3DArray(&array, &channelDesc, VS);
	
	cudaMemcpy3DParms copyParams = {0};
	copyParams.srcPtr   = make_cudaPitchedPtr((void*)data, VS.width*sizeof(float), VS.width, VS.height);
	copyParams.dstArray = array;
	copyParams.extent   = VS;
	copyParams.kind     = cudaMemcpyHostToDevice;
	cudaMemcpy3D(&copyParams);
	
	
	printf("Now bind\n");
	//cudaBindTexture(tex,array);
	cuda_arrays[num_cuda_arrays].array = array;
	num_cuda_arrays++;
	return num_cuda_arrays-1;
}

void bind_cuda_texture(const int idx) {
	cudaBindTexture(tex,cuda_arrays[idx].array);
}

void device_init() {
	static bool init = true;
	
	if (init) {
		int deviceCount;
		cudaGetDeviceCount(&deviceCount);
		printf("%d devices detected\n",deviceCount);
		if (deviceCount == 0) exit(1);
		
		cudaDeviceProp deviceProp;
		cudaGetDeviceProperties(&deviceProp, 0);
		if (deviceProp.major < 1) exit(2);
		
		cudaSetDevice(0);
		init_cuda_emdata_arrays();
		init = false; //Force init everytikme
	}
}