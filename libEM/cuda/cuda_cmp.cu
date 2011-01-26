#include "cuda_util.h"
#include <sm_11_atomic_functions.h>

// these mean, mean square kernals work by sweeping over the first array dim
__global__ void mean_kernal(const float * data, float * device_stats, const int size, const int num_calcs, const int num_threads, const int offset)
{
	float mean = 0.0f;
	
	const uint x=threadIdx.x;
	const uint y=blockIdx.x;
	int idx = x + y*num_threads + offset;
	
	for(int i = 0; i < size; i++){
		int index = i*size + idx % size + ((idx/size)*size*size); //for coalesing
		mean += data[index];
	}
	
	device_stats[idx] = mean/size;
}

__global__ void meansquare_kernal(const float * data, float * device_stats, const int size, const int num_calcs, const int num_threads, const int offset)
{
	float meansq = 0.0f;
	
	const uint x=threadIdx.x;
	const uint y=blockIdx.x;
	int idx = x + y*num_threads + offset;
	
	for(int i = 0; i < size; i++){
		int index = i*size + idx % size + ((idx/size)*size*size); //for coalesing
		meansq += data[index]*data[index];
	}
	
	device_stats[idx] = meansq/size;
}

__global__ void sumup_kernal(const float * data, float * device_stats, const int size, const int dim2size, const int num_threads, const int offset)
{
	float sum = 0.0f;
	
	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	int idx = x + y*num_threads + offset;
	
	for(int i = 0; i < size; i++){
		int index = i*dim2size + idx % dim2size;
		sum += data[index];
	}
	
	device_stats[idx] = sum/size; 
}


__global__ void dot_cmp_kernaldm(const float* data1, const float* data2, const float* dm, float* device_soln, const int size, const int num_threads, const int offset)
{
	float dot = 0.0f;
	float nnn = 0.0f;

	int idx = threadIdx.x + blockIdx.x*num_threads + offset;
	
	for(int i = 0; i < size; i++){
		int index = i*size + idx % size + ((idx/size)*size*size); //for coalesing
		if(dm[index] > 0.5){
			dot += data1[index]*data2[index];
			nnn += 1.0f;
		}
	}
	
	device_soln[idx] = dot/nnn;	
	
}

__global__ void dot_cmp_kernal(const float* data1, const float* data2, float* device_soln, const int size, const int num_threads, const int offset)
{
	float dot = 0.0f;

	int idx = threadIdx.x + blockIdx.x*num_threads + offset;
	
	for(int i = 0; i < size; i++){
		int index = i*size + idx % size + ((idx/size)*size*size); //for coalesing
		dot += data1[index]*data2[index];
	}
	
	device_soln[idx] = dot/size;	
	
}

__global__ void dot_cmp_kernal_atomic(const float* data1, const float* data2, float* device_soln, const int num_threads, const int offset)
{

	float dot = 0.0f;
	
	const uint x=threadIdx.x;
	const uint y=blockIdx.x;
	
	int idx = x + y*num_threads + offset;
	
	dot = data1[idx]*data1[idx];
	
	atomicAdd(device_soln, dot); //VERY SLOW!!!!!

}

__shared__ float sdata[MAX_THREADS];
__global__ void dot_cmp_kernal_reduce(float *g_idata1, float *g_idata2, float *g_odata) 
{
	extern __shared__ float sdata[];
	// each thread loads one element from global to shared mem
	unsigned int tid = threadIdx.x;
	unsigned int i = blockIdx.x*(blockDim.x*2) + threadIdx.x;
	sdata[tid] = g_idata1[i]*g_idata2[i] + g_idata1[i+blockDim.x]*g_idata2[i+blockDim.x];
	__syncthreads();

	// do reduction in shared mem
	for (unsigned int s=blockDim.x/2; s>0; s>>=1) {
		if (tid < s) {
		sdata[tid] += sdata[tid + s];
		}
		__syncthreads();
	}

	// write result for this block to global mem
	if (tid == 0) g_odata[blockIdx.x] = sdata[0];

}

__global__ void kernal_reduce(float *g_idata, float *g_odata) 
{
	extern __shared__ float sdata[];
	// each thread loads one element from global to shared mem
	unsigned int tid = threadIdx.x;
	unsigned int i = blockIdx.x*(blockDim.x*2) + threadIdx.x;
	sdata[tid] = g_idata[i] + g_idata[i+blockDim.x];
	__syncthreads();

	// do reduction in shared mem
	for (unsigned int s=blockDim.x/2; s>0; s>>=1) {
		if (tid < s) {
		sdata[tid] += sdata[tid + s];
		}
		__syncthreads();
	}

	// write result for this block to global mem
	if (tid == 0) g_odata[blockIdx.x] = sdata[0];

}

float * device_results=0;
float dot_cmp_cuda(float* data1, float* data2, const float* dm, const int &nx, const int &ny, const int &nz)
{
	/*
	int num_calcs = nx*ny*nz;
	int rsize = num_calcs/MAX_THREADS;
	//float * device_results=0;
	//cudaMalloc((void **)&device_results, rsize*sizeof(float));
	const dim3 blockSize(128,1,1);
	const dim3 gridSize(rsize,1,1);
	dot_cmp_kernal_reduce<<<gridSize, blockSize>>>(data1, data2, device_results);
	
	//now reduce
	for(int i = rsize/MAX_THREADS; i > 0 ; i /= MAX_THREADS)
	{
		const dim3 blockSizei(128,1,1);
		const dim3 gridSizei(i,1,1);
		kernal_reduce<<<gridSizei, blockSizei>>>(device_results, device_results);
	}
	
	float * host_results = (float*) malloc(sizeof(float));
	cudaMemcpy(host_results,device_results,sizeof(float),cudaMemcpyDeviceToHost);
	//cudaFree(device_results);
	
	float dot = host_results[0];
	dot /= (nx*ny*nz);
	free(host_results);
	
	return dot;
	*/
	/* This method is absoutly terrible!!!!
	float * device_results=0;
	cudaMalloc((void **)&device_results, sizeof(float));
	
	int num_calcs = nx*ny*nz;
	int grid_y = num_calcs/(MAX_THREADS);
	int res_y = num_calcs - grid_y*MAX_THREADS;
	
	const dim3 blockSize(MAX_THREADS,1,1);
	const dim3 gridSize(grid_y,1,1);
	dot_cmp_kernal_atomic<<<gridSize, blockSize>>>(data1, data2, device_results, MAX_THREADS, 0);
	cudaThreadSynchronize();
	
	float * host_results = (float*) malloc(sizeof(float));
	cudaMemcpy(host_results,device_results,sizeof(float),cudaMemcpyDeviceToHost);
	cudaFree(device_results);
	
	float dot = host_results[0]/num_calcs;
	free(host_results);
	
	return dot;
	*/
	/*
	//This method is slower
	//first calc the mean and store it there
	int darraysize = ny*nz + nz + 1;
	float * device_results=0;
	cudaMalloc((void **)&device_results, darraysize*sizeof(float));
	
	//First calculate the dot
	int num_calcs = ny*nz;
	int grid_y = num_calcs/(MAX_THREADS);
	int res_y = num_calcs - grid_y*MAX_THREADS;
	int ardim;
	if(nz==1){
		ardim = 1;
	}else{
		ardim = nz+1;
	}
	
	if (grid_y > 0) {
		const dim3 blockSize(MAX_THREADS,1,1);
		const dim3 gridSize(grid_y,1,1);
		dot_cmp_kernal<<<gridSize, blockSize>>>(data1, data2, device_results+ardim, nx, MAX_THREADS, 0);
	}
	if (res_y) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		int offset = grid_y*MAX_THREADS;
		dot_cmp_kernal<<<gridSize, blockSize>>>(data1, data2, device_results+ardim, nx, MAX_THREADS, offset);
		printf("Z\n");
	}
	cudaThreadSynchronize();
	
	if(nz == 1){
		const dim3 blockSize(1,1,1);
		const dim3 gridSize(1,1,1);
		sumup_kernal<<<gridSize, blockSize>>>(device_results+1, device_results, ny, 1, MAX_THREADS, 0);
	}else{
		grid_y = nz/(MAX_THREADS);
		res_y = nz - grid_y*MAX_THREADS;
		if (grid_y > 0) {
			const dim3 blockSize(MAX_THREADS,1,1);
			const dim3 gridSize(grid_y,1,1);
			sumup_kernal<<<gridSize, blockSize>>>(device_results+ardim, device_results+1, ny, nz, MAX_THREADS, 0);
		}
		if (res_y){
			const dim3 blockSize(res_y,1,1);
			const dim3 gridSize(1,1,1); //obviously only need one block here
			int offset = grid_y*MAX_THREADS;
			sumup_kernal<<<gridSize, blockSize>>>(device_results+ardim, device_results+1, ny, nz, MAX_THREADS, offset);
		}
		cudaThreadSynchronize();
		const dim3 blockSize(1,1,1);
		const dim3 gridSize(1,1,1);
		sumup_kernal<<<gridSize, blockSize>>>(device_results+1, device_results, nz, 1, MAX_THREADS, 0);
	}
	cudaThreadSynchronize();	
	
	float * host_results = (float*) malloc(sizeof(float));
	cudaMemcpy(host_results,device_results,sizeof(float),cudaMemcpyDeviceToHost);
	cudaFree(device_results);
	
	float dot = host_results[0];
	free(host_results);
	
	return dot;
	*/
	
	
	//this version is faster
	const int size = nx;
	int num_calcs = ny*nz;
	int grid_y = num_calcs/(MAX_THREADS);
	int res_y = num_calcs - grid_y*MAX_THREADS;
	
	float * device_soln=0;
	cudaMalloc((void **)&device_soln, num_calcs*sizeof(float));
	float * host_soln = 0;
	host_soln = (float*) malloc(num_calcs*sizeof(float));

	if (grid_y > 0) {
		const dim3 blockSize(MAX_THREADS,1,1);
		const dim3 gridSize(grid_y,1,1);
		if(dm){
			dot_cmp_kernaldm<<<gridSize, blockSize>>>(data1, data2, dm, device_soln, size, MAX_THREADS, 0);
		}else{
			dot_cmp_kernal<<<gridSize, blockSize>>>(data1, data2, device_soln, size, MAX_THREADS, 0);
		}
	}else{
		const dim3 blockSize(num_calcs,1,1);
		const dim3 gridSize(1,1,1); //obviously only need one block here
		if(dm){
			dot_cmp_kernaldm<<<gridSize, blockSize>>>(data1, data2, dm, device_soln, size, MAX_THREADS, 0);
		}else{
			dot_cmp_kernal<<<gridSize, blockSize>>>(data1, data2, device_soln, size, MAX_THREADS, 0);
		}
	}
	if (res_y) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		int offset = grid_y*MAX_THREADS;
		if(dm){
			dot_cmp_kernaldm<<<gridSize,blockSize>>>(data1, data2, dm, device_soln, size, MAX_THREADS, offset);
		}else{
			dot_cmp_kernal<<<gridSize,blockSize>>>(data1, data2, device_soln, size, MAX_THREADS, offset);
		}
	}
	
	cudaThreadSynchronize();
	cudaMemcpy(host_soln,device_soln,num_calcs*sizeof(float),cudaMemcpyDeviceToHost);
	cudaFree(device_soln);
	
	float dot = 0.0f;
	float nnn = 0.0f;
	for(int i = 0; i < num_calcs; i++){
		
		if(host_soln[i] == host_soln[i]){ //this wil return false when host_soln[i] is NaN
			dot += host_soln[i];
			nnn++;
		}
	}
	
	free(host_soln);

	
	dot = dot/nnn;

	return dot;
	
}

__global__ void norm_kernal(float * data, float mean, float var, int num_threads, int offset)
{
	
	int idx = threadIdx.x + blockIdx.x*num_threads + offset;
	
	data[idx] = (data[idx] - mean)/ var;
	
}

void normalize_cuda(float * data, float mean, float var, const int nx, const int ny, const int nz)
{
	
	const int num_calcs = nx*ny*nz;
	const int grid_y = num_calcs/(MAX_THREADS);
	const int res_y = num_calcs - grid_y*MAX_THREADS;
	
	if (grid_y > 0) {
		const dim3 blockSize(MAX_THREADS,1,1);
		const dim3 gridSize(grid_y,1,1);
		norm_kernal<<<gridSize, blockSize>>>(data, mean, var, MAX_THREADS, 0);
	}
	if (res_y) {	
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		int offset = grid_y*MAX_THREADS;
		norm_kernal<<<gridSize,blockSize>>>(data, mean, var, MAX_THREADS, offset);	
	}
}

__global__ void stats_kernal(const float *data, float * device_soln, const int size, const int num_calcs, const int num_threads, const int offset)
{
	
	float mean = 0.0f;
	float var = 0.0f;
	
	int idx = threadIdx.x + blockIdx.x*num_threads + offset;
	
	for(int i = 0; i < size; i++){
		int index = i*size + idx % size + ((idx/size)*size*size); //for coalesing
		float datum = data[index]; //so we dno't need multiple accesses to global mem, I would think that the compiler would optimize this, but the manual said to program like this....
		mean += datum;
		var += datum*datum;
	}
	
	device_soln[idx] = mean/size;
	device_soln[idx + num_calcs] = var/size;
	
}

float2 get_stats_cuda(const float * data, const int nx, const int ny, const int nz)
{
	
	
	const int size = nx;
	const int num_calcs = ny*nz;
	const int grid_y = num_calcs/(MAX_THREADS);
	const int res_y = num_calcs - grid_y*MAX_THREADS;
	
	float * device_soln=0;
	cudaMalloc((void **)&device_soln, 2*num_calcs*sizeof(float));
	float * host_soln = 0;
	host_soln = (float*) malloc(2*num_calcs*sizeof(float));
	
	if (grid_y > 0) {
		const dim3 blockSize(MAX_THREADS,1,1);
		const dim3 gridSize(grid_y,1,1);
		stats_kernal<<<gridSize, blockSize>>>(data, device_soln, size, num_calcs, MAX_THREADS, 0);
	}else{
		const dim3 blockSize(num_calcs,1,1);
		const dim3 gridSize(1,1,1); //obviously only need one block here
		stats_kernal<<<gridSize, blockSize>>>(data, device_soln, size, num_calcs, MAX_THREADS, 0);
	}
	if (res_y) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		int offset = grid_y*MAX_THREADS;
		stats_kernal<<<gridSize,blockSize>>>(data, device_soln, size, num_calcs, MAX_THREADS, offset);
	}
	
	cudaThreadSynchronize();
	cudaMemcpy(host_soln,device_soln,2*num_calcs*sizeof(float),cudaMemcpyDeviceToHost);
	cudaFree(device_soln);
	
	float mean = 0.0f;
	float var = 0.0f;
	float nnn = 0.0f;
	for(int i = 0; i < num_calcs; i++){
		
		if(host_soln[i] == host_soln[i]){ //this will return false when host_soln[i] is NaN, in this case only need to check once, if mean = NaN then so will var
			mean += host_soln[i];
			var += host_soln[i + num_calcs];
			nnn++;
		}
	}
	
	free(host_soln);

	float2 stats;
	stats.x = mean/nnn;
	stats.y = var/nnn - stats.x*stats.x;
	
	
	return stats;
	
}

float getvalueat_cuda(float * data, int tx, int ty, int tz, int nx, int ny, int nz)
{
	float * host_soln = 0;
	host_soln = (float*) malloc(sizeof(float));
	cudaMemcpy(host_soln,data,sizeof(float),cudaMemcpyDeviceToHost);
	float peak = host_soln[0];
	free(host_soln);
	
	return peak;
	
}
__global__ void ccc_cmp_kernaldm(const float* data1, const float* data2, const float* dm, float* device_soln, const int size, const int num_calcs, const int num_threads, const int offset)
{
	float avg1 = 0.0f;
	float avg2 = 0.0f;
	float var1 = 0.0f;
	float var2 = 0.0f;
	float ccc = 0.0f;
	float nnn = 0.0f;
	
	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	int idx = x + y*num_threads + offset;
	
	for(int i = 0; i < size; i++){
		int index = i*size + idx % size + ((idx/size)*size*size); //for coalesing
		if(dm[index] > 0.5){
			avg1 += data1[index];
			avg2 += data2[index];
			var1 += data1[index]*data1[index];
			var2 += data2[index]*data2[index];
			ccc += data1[index]*data2[index];
			nnn += 1.0f;
		}
	}
	
	device_soln[idx] = avg1/nnn;
	device_soln[idx + num_calcs] = avg2/nnn;
	device_soln[idx + 2*num_calcs] = var1/nnn;
	device_soln[idx + 3*num_calcs] = var2/nnn;
	device_soln[idx + 4*num_calcs] = ccc/nnn;
	
	
}

__global__ void ccc_cmp_kernal(const float* data1, const float* data2, float* device_soln, const int size, const int num_calcs, const int num_threads, const int offset)
{
	float avg1 = 0.0f;
	float avg2 = 0.0f;
	float var1 = 0.0f;
	float var2 = 0.0f;
	float ccc = 0.0f;
	
	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	int idx = x + y*num_threads + offset;
	
	for(int i = 0; i < size; i++){
		int index = i*size + idx % size + ((idx/size)*size*size); //for coalesing
		avg1 += data1[index];
		avg2 += data2[index];
		var1 += data1[index]*data1[index];
		var2 += data2[index]*data2[index];
		ccc += data1[index]*data2[index];
	}
	
	device_soln[idx] = avg1/size;
	device_soln[idx + num_calcs] = avg2/size;
	device_soln[idx + 2*num_calcs] = var1/size;
	device_soln[idx + 3*num_calcs] = var2/size;
	device_soln[idx + 4*num_calcs] = ccc/size;
	
	
}

float ccc_cmp_cuda(const float* data1, const float* data2, const float* dm, const int &nx, const int &ny, const int &nz)
{
	
	const int size = nx;
	int num_calcs = ny*nz;
	int grid_y = num_calcs/(MAX_THREADS);
	int res_y = num_calcs - grid_y*MAX_THREADS;	
	
	float * device_soln=0;
	cudaMalloc((void **)&device_soln, 5*num_calcs*sizeof(float));
	float * host_soln = 0;
	host_soln = (float*) malloc(5*num_calcs*sizeof(float));

	if (grid_y > 0) {
		const dim3 blockSize(MAX_THREADS,1,1);
		const dim3 gridSize(grid_y,1,1);
		if(dm){
			ccc_cmp_kernaldm<<<gridSize, blockSize>>>(data1, data2, dm, device_soln, size, num_calcs, MAX_THREADS,0);
		}else{;
			ccc_cmp_kernal<<<gridSize, blockSize>>>(data1, data2, device_soln, size, num_calcs, MAX_THREADS, 0);
		}
	}else{
		const dim3 blockSize(num_calcs,1,1);
		const dim3 gridSize(1,1,1); //obviously only need one block here
		if(dm){
			ccc_cmp_kernaldm<<<gridSize, blockSize>>>(data1, data2, dm, device_soln, size, num_calcs, MAX_THREADS,0);
		}else{
			ccc_cmp_kernal<<<gridSize, blockSize>>>(data1, data2, device_soln, size, num_calcs, MAX_THREADS, 0);
		}
	}
	if (res_y) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		int offset = grid_y*MAX_THREADS;
		if(dm){
			ccc_cmp_kernaldm<<<gridSize,blockSize>>>(data1, data2, dm, device_soln, size, num_calcs, MAX_THREADS, offset);
		}else{
			ccc_cmp_kernal<<<gridSize,blockSize>>>(data1, data2, device_soln, size, num_calcs, MAX_THREADS, offset);
		}
	}
	
	cudaThreadSynchronize();
	cudaMemcpy(host_soln,device_soln,5*num_calcs*sizeof(float),cudaMemcpyDeviceToHost);
	cudaFree(device_soln);
	
	float avg1 = 0.0f;
	float avg2 = 0.0f;
	float var1 = 0.0f;
	float var2 = 0.0f;
	float ccc = 0.0f;
	float nnn = 0.0f;
	for(int i = 0; i < num_calcs; i++){	
		if(host_soln[i] == host_soln[i]){ //this wil return false when host_soln[i] is NaN
			avg1 += host_soln[i];
			avg2 += host_soln[i + num_calcs];
			var1 += host_soln[i + 2*num_calcs];
			var2 += host_soln[i + 3*num_calcs];
			ccc += host_soln[i + 4*num_calcs];
			nnn++;
		}
	}
	
	free(host_soln);

	avg1 = avg1/nnn;
	var1 = var1/nnn - avg1*avg1;
	avg2 = avg2/nnn;
	var2 = var2/nnn - avg2*avg2;
	ccc = ccc/nnn - avg1*avg2;
	ccc /= sqrt(var1*var2);
	//ccc = avg2;

	return ccc;
	
}

/*
float* calc_fourier_shell_correlation_cuda(const int nx, const int ny, const int nz, const int d)
{
	
}
*/

/*
typedef unsigned int uint;

__global__ void phase_cmp_weights(float* out, int num_threads, int nx, int ny, int nz, int nxy, float np, int offset) {

	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	int idx = x + y*num_threads + offset;

	int tx = idx % nx;
	int tz = idx / nxy;
	int ty = (idx - tz*nxy)/nx;

	tx /= 2;

	if (ty > ny/2) ty = ny-ty;
	if (tz > nz/2) tz = nz-tz;

	float rad = sqrtf(tx*tx+ty*ty+tz*tz);
	float x2 = 10.0*rad/np;
	out[idx] = x2*exp(-x2);
}

void calc_phase_weights_cuda(EMDataForCuda* t,float np) {

	int max_threads = 192;
	int num_calcs = t->ny*t->nx*t->nz;

	int grid_y = num_calcs/(max_threads);
	int res_y = num_calcs - grid_y*max_threads;

	if (grid_y > 0) {
		const dim3 blockSize(max_threads,1, 1);
		const dim3 gridSize(grid_y,1,1);
		phase_cmp_weights<<<gridSize,blockSize>>>(t->data,max_threads,t->nx,t->ny,t->nz,t->nx*t->ny,np,0);
	}

	if (res_y) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		phase_cmp_weights<<<gridSize,blockSize>>>(t->data,max_threads,t->nx,t->ny,t->nz,t->nx*t->ny,np,grid_y*max_threads);
	}

	cudaThreadSynchronize();

}

__global__ void mean_phase_error_kernel(float *ldata,float *rdata,int num_threads){
	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	const uint idx = x + y*num_threads;
	const uint idx2 = 2*idx;
	// This could be optimized using shared memory  -
	// Coalesce the reading of ldata memory into shared memory
	// Then operate on shared memory
	// Then coaless the writing of data to rdata
	// See the ri2ap etc for an implementation of this idea
	rdata[idx] = ldata[idx2]+ldata[idx2+1];

}


void histogram_sum(EMDataForCuda* hist, const int num_hist){
	for(int i = 0; i < num_hist-1; ++i) {
		EMDataForCuda* l = hist+i;
		EMDataForCuda* r = hist+(i+1);
		int max_threads = 192;
		int num_calcs = (l->nx*l->ny*l->nz)/2;
		num_calcs += num_calcs%2;
		int grid_y = num_calcs/(max_threads);
		int res_y = num_calcs - grid_y*max_threads;

		if (grid_y > 0) {
			const dim3 blockSize(max_threads,1, 1);
			const dim3 gridSize(grid_y,1,1);
			mean_phase_error_kernel<<<gridSize,blockSize>>>(l->data,r->data,max_threads);
		}

		if (res_y) {
			const dim3 blockSize(res_y,1, 1);
			const dim3 gridSize(1,1,1);
			int inc = grid_y*max_threads;
			mean_phase_error_kernel<<<gridSize,blockSize>>>(l->data+inc,r->data+inc,max_threads);
		}
		cudaThreadSynchronize();
	}

}

__device__ float angle_sub_2pi(float x, float y)
{
	float twopi = 6.2831853071795862;
	float r = fmod(fabs(x - y), twopi );
	if (r >  3.1415926535897931 ) {
		r = (float) (twopi - r);
	}

	return r;
}

__global__ void mean_phase_error_kernel(float *ldata,float *rdata,float *wdata,float* hdata,float* ndata, int num_threads)
{
	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	const uint idx = 2*x + y*num_threads;
	const uint idxp1 = idx+1;

	float l1 = ldata[idx];
	float l2 = ldata[idxp1];
	float amp = sqrtf(l1*l1+l2*l2);
	float phase1 = atan2(l2,l1);
	float r1 = rdata[idx];
	float r2 = rdata[idxp1];
	float phase2 = atan2(r2,r1);

	float a = wdata[idx] * amp;
	float f = angle_sub_2pi(phase1,phase2) * a;

	const uint idx2 = idx/2;
	hdata[idx2] = f;
	ndata[idx2] = a;
}

__global__ void single_normal_kernel(float *ldata,float *rdata)
{
	ldata[0] = (*ldata)/(*rdata);
}

void mean_phase_error_cuda(EMDataForCuda* left, EMDataForCuda* right, EMDataForCuda* wt,EMDataForCuda* hist, EMDataForCuda* norm, const int num_hist) {
	int max_threads = 192; // I halve the threads because each kernel access memory in two locations

	int num_calcs = left->nx*left->ny*left->nz;

	int grid_y = num_calcs/(2*max_threads);
	int res_y = (num_calcs - (2*grid_y*max_threads))/2;

	if ( grid_y > 0 ) {
		const dim3 blockSize(max_threads,1, 1);
		const dim3 gridSize(grid_y,1,1);
		mean_phase_error_kernel<<<gridSize,blockSize>>>(left->data,right->data,wt->data,hist[0].data,norm[0].data,2*max_threads);
	}

	if ( res_y > 0 ) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		int inc = 2*grid_y*max_threads;
		int inc2 = grid_y*max_threads;
		mean_phase_error_kernel<<<gridSize,blockSize>>>(left->data+inc,right->data+inc,wt->data+inc2,hist[0].data+inc2,norm[0].data+inc2,2*max_threads);
	}
	cudaThreadSynchronize();

	histogram_sum(hist,num_hist);
	histogram_sum(norm,num_hist);

	single_normal_kernel<<<1,1>>>(hist[num_hist-1].data,norm[num_hist-1].data);
	cudaThreadSynchronize();
}
*/

