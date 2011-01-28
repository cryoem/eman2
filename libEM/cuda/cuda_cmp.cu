#include "cuda_util.h"

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


