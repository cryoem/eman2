#include "cuda_util.h"

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

__global__ void norm_kernal(float * data, float mean, float var, int totaltc)
{
	
	const uint index = threadIdx.x + (blockIdx.x + gridDim.x*blockIdx.y)*MAX_THREADS;
	
	if(index < totaltc){
		data[index] = (data[index] - mean)/var;
	}
	
}

void normalize_cuda(float * data, float mean, float var, const int nx, const int ny, const int nz)
{
	
	int grid = int(ceil(sqrt(nx*ny*nz/MAX_THREADS)));
	
	const dim3 blockSize(MAX_THREADS,1,1);
	const dim3 gridSize(grid,grid,1);
	norm_kernal<<<gridSize, blockSize>>>(data, mean, var, nx*ny*nz);
	
	cudaThreadSynchronize();
}

__global__ void stats_kernal(const float *data, float * device_soln, const int size, const int num_calcs, const int num_threads, const int offset)
{
	
	float sum = 0.0f;
	float sum_sq = 0.0f;
	
	int idx = threadIdx.x + blockIdx.x*num_threads + offset;
	
	for(int i = 0; i < size; i++){
		int index = i*size + idx % size + ((idx/size)*size*size); //for coalesing
		float datum = data[index]; //so we dno't need multiple accesses to global mem, I would think that the compiler would optimize this, but the manual said to program like this....
		sum += datum;
		sum_sq += datum*datum;
	}
	
	device_soln[idx] = sum;
	device_soln[idx + num_calcs] = sum_sq;
	
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
	
	double sum = 0.0f;
	double sum_sq = 0.0f;
	double nnn = 0.0f;
	for(int i = 0; i < num_calcs; i++){
		if(host_soln[i] == host_soln[i]){ //this will return false when host_soln[i] is NaN, in this case only need to check once, if mean = NaN then so will var
			sum += host_soln[i];
			sum_sq += host_soln[i + num_calcs];
			nnn+=size;
		}
	}
	
	free(host_soln);

	float2 stats;
	stats.x = float(sum/nnn);
	stats.y = float((sum_sq - sum*stats.x)/(nnn-1));
	
	
	return stats;
	
}

float get_value_at_wrap_cuda(float * data, int tx, int ty, int tz, int nx, int ny, int nz)
{
	float * host_soln = 0;
	int lx = tx;
	int ly = ty;
	int lz = tz;

	if (lx < 0) lx = nx + lx;
	if (ly < 0) ly = ny + ly;
	if (lz < 0) lz = nz + lz;

	int index = (lx + ly * nx + lz * nx * ny);
	host_soln = (float*) malloc(sizeof(float));
	cudaMemcpy(host_soln,data+index,sizeof(float),cudaMemcpyDeviceToHost);
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

__global__ void fsc_tomo_cmp_kernal(const float* data1, const float* data2, float* device_soln, const float data1threshold, const float data2threshold, const int nx, const int ny, const int nz, const int offset)
{
	
	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	int idx = x + y*MAX_THREADS + offset;
	
	float sum_data1_amps = 0.0;
	float sum_data2_amps = 0.0;
	float top = 0.0;
	for(int i = 0; i < ny; i++){
		//int index = i*nx + idx % nx + ((idx/nx)*ny*nz); //for coalesing
		int rindex = i*nx + 2*(idx % nx/2) + (2*idx/nx)*ny*nz;
		int iindex = i*nx + 2*(idx % nx/2)+ 1 + (2*idx/nx)*ny*nz;
		float data1_r = data1[rindex];
		float data1_i = data1[iindex];
		float data2_r = data2[rindex];
		float data2_i = data2[iindex];
		if((data1_r* data1_r +  data1_i*data1_i) > data1threshold && (data2_r* data2_r +  data2_i*data2_i) > data2threshold){
			sum_data1_amps += (data1_r* data1_r +  data1_i*data1_i);
			sum_data2_amps += (data2_r* data2_r +  data2_i*data2_i);
			top += (data1_r*data2_r + data1_i*data2_i);
		}
	}
	device_soln[idx*3] = top;
	device_soln[idx*3 +1] = sum_data1_amps;
	device_soln[idx*3 +2] = sum_data2_amps;
	
}

float fsc_tomo_cmp_cuda(const float* data1, const float* data2, const float data1threshold, const float data2threshold, const float minres, const float maxres, const int nx, const int ny, const int nz)
{
	
	int num_calcs = (nx/2)*nz;
	int grid_y = num_calcs/(MAX_THREADS);
	int res_y = num_calcs - grid_y*MAX_THREADS;	
	
	float * device_soln=0;
	cudaMalloc((void **)&device_soln, 3*num_calcs*sizeof(float));
	float * host_soln = 0;
	host_soln = (float*) malloc(3*num_calcs*sizeof(float));

	if (grid_y > 0) {
		const dim3 blockSize(MAX_THREADS,1,1);
		const dim3 gridSize(grid_y,1,1);
		fsc_tomo_cmp_kernal<<<gridSize, blockSize>>>(data1, data2, device_soln, data1threshold, data2threshold, nx, ny, nz,0);

	}else{
		const dim3 blockSize(num_calcs,1,1);
		const dim3 gridSize(1,1,1); //obviously only need one block here
		fsc_tomo_cmp_kernal<<<gridSize, blockSize>>>(data1, data2, device_soln, data1threshold, data2threshold, nx, ny, nz,0);
	}
	if (res_y) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		int offset = grid_y*MAX_THREADS;
		fsc_tomo_cmp_kernal<<<gridSize, blockSize>>>(data1, data2, device_soln, data1threshold, data2threshold, nx, ny, nz, offset);
	}
	
	cudaThreadSynchronize();
	cudaMemcpy(host_soln,device_soln,3*num_calcs*sizeof(float),cudaMemcpyDeviceToHost);
	cudaFree(device_soln);
	
	double top = 0.0;
	double sum_data1_amps = 0.0;
	double sum_data2_amps = 0.0;
	for(int i = 0; i < num_calcs; i++){
		if(!isnan(host_soln[i])){ //this wil return false when host_soln[i] is NaN
			top += host_soln[i*3];
			sum_data1_amps += host_soln[i*3 + 1];
			sum_data1_amps += host_soln[i*3 + 2];
		}
	}
	
	free(host_soln);
	
	return top/sqrt(sum_data1_amps*sum_data2_amps);
	
}