
#include <cuda.h>
#include <stdio.h>

#include "cuda_defs.h"
#include "cuda_emfft.h"
#include "cuda_util.h"

#include "emcudautil.h"

texture<float, 3, cudaReadModeElementType> tex;
texture<float, 2, cudaReadModeElementType> tex2d;

void cuda_bind_texture_3d(texture<float, 3, cudaReadModeElementType> &tex,cudaArray *array) {
	tex.normalized = 0;
	tex.filterMode = cudaFilterModeLinear;
	tex.addressMode[0] = cudaAddressModeClamp;
	tex.addressMode[1] = cudaAddressModeClamp;
	tex.addressMode[2] = cudaAddressModeClamp;
	
	cudaBindTextureToArray(tex, array);
}

void cuda_bind_texture_2d(texture<float, 2, cudaReadModeElementType> &tex,cudaArray *array) {
	tex.normalized = 0;
	tex.filterMode = cudaFilterModeLinear;
	tex.addressMode[0] = cudaAddressModeClamp;
	tex.addressMode[1] = cudaAddressModeClamp;
// 	tex.addressMode[2] = cudaAddressModeClamp;
	
	cudaBindTextureToArray(tex, array);
}


struct CudaEMDataArray {
	cudaArray* array;
	const float* data; /*This one may be unecessary*/
	void* emdata_pointer;
	int nx;
	int ny;
	int nz;
};

void copy_array_data(CudaEMDataArray* to, CudaEMDataArray* from) {
	to->array = from->array;
	to->data = from->data;
	to->emdata_pointer = from->emdata_pointer;
	to->nx = from->nx;
	to->ny = from->ny;
	to->nz = from->nz;
}

void set_array_data_null(CudaEMDataArray* p)
{
	p->array = 0;
	p->data = 0;
	p->emdata_pointer = 0;
	p->nx = 0;
	p->ny = 0;
	p->nz = 0;
}

const int max_cuda_arrays = 10;
int num_cuda_arrays = 0;
CudaEMDataArray cuda_arrays[max_cuda_arrays];

void init_cuda_emdata_arrays() {
	for(int i = 0; i < max_cuda_arrays; ++i ) {
		CudaEMDataArray c =  { 0, 0, 0, 0, 0, 0};
		cuda_arrays[i] = c;
	}
}

int make_cuda_array_space_0_free() {
	//printf("Freeing space 0\n");
	//debug_arrays();
	int n = num_cuda_arrays-1;
	if (cuda_arrays[n].emdata_pointer != 0) set_emdata_cuda_array_handle(-1,cuda_arrays[n].emdata_pointer);
	CUDA_SAFE_CALL(cudaFree(cuda_arrays[n].array));
	cuda_arrays[n].array = 0;
	
	for (int i = 0; i < num_cuda_arrays-1; ++i ) {
		CudaEMDataArray* to = &cuda_arrays[i+1];
		CudaEMDataArray* from = &cuda_arrays[i];
		copy_array_data(to,from);
		if (to->emdata_pointer != 0) set_emdata_cuda_array_handle(i+1,to->emdata_pointer);
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
// 	printf("Making a new cuda array\n");
	// If we make it here then it doesn't exist
	CudaEMDataArray* c = &cuda_arrays[idx];
	c->data = data;
	c->emdata_pointer = emdata_pointer;
	c->nx = nx;
	c->ny = ny;
	c->nz = nz;
	
	cudaArray *array = 0;
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	
	if (nz > 1) {
		cudaExtent VS = make_cudaExtent(nx,ny,nz);
// 		printf("It's a 3D one\n");
		cudaMalloc3DArray(&array, &channelDesc, VS);
		
		cudaMemcpy3DParms copyParams = {0};
		copyParams.srcPtr   = make_cudaPitchedPtr((void*)data, VS.width*sizeof(float), VS.width, VS.height);
		copyParams.dstArray = array;
		copyParams.extent   = VS;
		copyParams.kind     = cudaMemcpyDeviceToDevice;
		cudaMemcpy3D(&copyParams);
	} else if ( ny > 1) {
// 		printf("It's a 2D one\n");d
		cudaMallocArray(&array,&channelDesc,nx,ny);
		cudaMemcpy2DToArray(array, 0, 0, data, nx*sizeof(float), nx*sizeof(float), ny, cudaMemcpyDeviceToDevice);
	} else throw;
	
	c->array = array;
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
		if (to->emdata_pointer != 0) set_emdata_cuda_array_handle(i,to->emdata_pointer);
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
	CudaEMDataArray* c = &cuda_arrays[idx];
	if (c->nz > 1) {
		cuda_bind_texture_3d(tex,cuda_arrays[idx].array);
	} else if ( c->ny > 1 ) {
		cuda_bind_texture_2d(tex2d,cuda_arrays[idx].array);
	} else throw;
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
		init_cuda_fft_hh_plan_cache(); // should only be performed if the host is using Cuda ffts, which is unlikey. Will do after development has progressed.
		init_cuda_fft_dd_plan_cache();
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


void cuda_memcpy_host_to_device(const void* const host_rdata, void* device_rdata, const size_t num_bytes )
{
	device_init();
	cudaMemcpy(device_rdata,host_rdata,num_bytes,cudaMemcpyHostToDevice);
}

void cuda_malloc_device(void** device_rdata, const size_t num_bytes)
{
	device_init();
	int val = cudaMalloc(device_rdata, num_bytes);
	printf("Success? %d\n",(val==cudaSuccess));
}

void cuda_free_device(void* device_rdata)
{
	device_init(); // Technically unecessary seeing as the device must have been initialized to have allocate the pointer
	cudaFree(device_rdata);
}

__global__ void  calc_max_location_wrap(int* const soln, const float* data,const int maxdx, const int maxdy, const int maxdz, const int nx, const int ny, const int nz) {
	int maxshiftx = maxdx, maxshifty = maxdy, maxshiftz = maxdz;
	if (maxdx == -1) maxshiftx = nx/4;
	if (maxdy == -1) maxshifty = ny/4;
	if (maxdz == -1) maxshiftz = nz/4;

	float max_value = -10000000000000;

	for (int k = -maxshiftz; k <= maxshiftz; k++) {
		for (int j = -maxshifty; j <= maxshifty; j++) {
			for (int i = -maxshiftx; i <= maxshiftx; i++) {
				
				int kk = k;
				if (kk < 0) {
					kk = nz-kk;
				}
				int jj = j;
				if (jj < 0) {
					jj = nz-jj;
				}
				int ii = i;
				if (ii < 0) {
					ii = nz-ii;
				}
				float value = data[ii+jj*nx+kk*nx*ny];

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
	cudaMalloc((void **)&device_soln, 3*sizeof(int));
		
	int * host_soln = 0;
	host_soln = (int*) malloc(3*sizeof(int));
	
	const dim3 blockSize(1,1, 1);
	const dim3 gridSize(1,1,1);
		
	calc_max_location_wrap<<<blockSize,gridSize>>>(device_soln,data->data,maxdx,maxdy,maxdz,data->nx,data->ny,data->nz);
	
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

