
#include <cuda.h>
#include <stdio.h>

#include "cuda_defs.h"
#include "cuda_emfft.h"

#include "emcudautil.h"

texture<float, 3, cudaReadModeElementType> tex;

void cudaBindTexture(texture<float, 3, cudaReadModeElementType> &tex,cudaArray *array) {
	tex.normalized = 0;
	tex.filterMode = cudaFilterModeLinear;
	tex.addressMode[0] = cudaAddressModeClamp;
	tex.addressMode[1] = cudaAddressModeClamp;
	tex.addressMode[2] = cudaAddressModeClamp;
	
	CUDA_SAFE_CALL(cudaBindTextureToArray(tex, array));
}


struct CudaEMDataArray {
	cudaArray* array;
	const float* data; /*This one may be unecessary*/
	void* emdata_pointer;
};

void copy_array_data(CudaEMDataArray* to, CudaEMDataArray* from) {
	to->array = from->array;
	to->data = from->data;
	to->emdata_pointer = from->emdata_pointer;
}

void set_array_data_null(CudaEMDataArray* p)
{
	p->array = 0;
	p->data = 0;
	p->emdata_pointer = 0;
}

const int max_cuda_arrays = 2;
int num_cuda_arrays = 0;
CudaEMDataArray cuda_arrays[max_cuda_arrays];

void init_cuda_emdata_arrays() {
	for(int i = 0; i < max_cuda_arrays; ++i ) {
		CudaEMDataArray c =  { 0, 0 };
		cuda_arrays[i] = c;
	}
}


void debug_arrays()
{
	printf("debug\n");
	for(int i=0; i < max_cuda_arrays; ++i) 
	{
		CudaEMDataArray t = cuda_arrays[i];
		printf("%d %d %d\n",t.array,t.data,t.emdata_pointer);
	}
	printf("end debug\n");
}

int make_cuda_array_space_0_free() {
	//printf("Freeing space 0\n");
	//debug_arrays();
	int n = num_cuda_arrays-1;
	set_emdata_cuda_array_handle(-1,cuda_arrays[n].emdata_pointer);
	CUDA_SAFE_CALL(cudaFree(cuda_arrays[n].array));
	cuda_arrays[n].array = 0;
	
	for (int i = 0; i < num_cuda_arrays-1; ++i ) {
		CudaEMDataArray* to = &cuda_arrays[i+1];
		CudaEMDataArray* from = &cuda_arrays[i];
		copy_array_data(to,from);
		set_emdata_cuda_array_handle(i+1,to->emdata_pointer);
	}
	set_array_data_null(&cuda_arrays[0]);
	
	//debug_arrays();
	return 0;
}


int get_cuda_array_handle(const float * data,const int nx, const int ny, const int nz, void* emdata_pointer) {
	
	//printf("Get cuda array %d\n", emdata_pointer);
	//debug_arrays();
	for(int i = 0; i < num_cuda_arrays; ++i ) {
		if (cuda_arrays[i].emdata_pointer == emdata_pointer ) {
			//printf("Found that cuda arrary\n");
			return i;
		}
	}
	int idx = num_cuda_arrays;
	if (num_cuda_arrays == max_cuda_arrays) {
		make_cuda_array_space_0_free();
		idx = 0;
	}
	//printf("Making a new cuda array\n");
	// If we make it here then it doesn't exist
	cuda_arrays[idx].data = data;
	cuda_arrays[idx].emdata_pointer = emdata_pointer;
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
	
	cuda_arrays[idx].array = array;
	if (num_cuda_arrays != max_cuda_arrays) num_cuda_arrays++;
	return idx;
}

int delete_cuda_array(const int idx) {
	//printf("Deleting a cuda array\n");
	CUDA_SAFE_CALL(cudaFree(cuda_arrays[idx].array));
	cuda_arrays[idx].array = 0;
	
	for (int i = idx; i < num_cuda_arrays; ++i ) {
		CudaEMDataArray* to = &cuda_arrays[i];
		CudaEMDataArray* from = &cuda_arrays[i+1];
		copy_array_data(to,from);
		set_emdata_cuda_array_handle(i,to->emdata_pointer);
	}
	set_array_data_null(&cuda_arrays[num_cuda_arrays-1]);
	num_cuda_arrays--;
	
	return 0;
}

// highly specialized, use cautiously


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
		init_cuda_emfft_cache();
		init = false; //Force init everytikme
	}
}

void* cuda_malloc(const size_t size)
{
	device_init();
	void *mem=0;
	cudaMallocHost((void **)&mem, size);
	return mem;
}

void cuda_free(void* mem)
{
	device_init();
	cudaFreeHost(mem);
}
void cuda_memset(void* mem,int value, size_t size) {
	device_init();
	cudaMemset(mem,value,size);
}

void cuda_memcpy(void* dst, const void* const src, size_t count) {
	device_init();
	cudaMemcpy(dst,src,count,cudaMemcpyHostToHost);
// 	cudaStream_t stream;
// 	cudaStreamCreate(&stream);
// 	cudaMemcpyAsync(dst,src,count,cudaMemcpyHostToHost,stream);
// 	cudaStreamSynchronize(stream);
// 	cudaStreamDestroy(stream);
}

