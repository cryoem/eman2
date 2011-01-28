

#include "cuda_util.h"
#include <stdio.h>

// Global texture
extern texture<float, 3, cudaReadModeElementType> tex;
extern texture<float, 2, cudaReadModeElementType> tex2d;

typedef unsigned int uint;

#ifdef WIN32
	#define M_PI 3.14159265358979323846f
#endif	//WIN32

#define MAX_THREADS 192

__global__ void mult_kernel(float *data,const float scale,const int num_threads)
{

	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	data[x+y*num_threads] *= scale;
}

void emdata_processor_mult( EMDataForCuda* cuda_data, const float& mult) {

	int num_calcs = cuda_data->nx*cuda_data->ny*cuda_data->nz;

	int grid_y = num_calcs/MAX_THREADS;
	int res_y = num_calcs - (grid_y*MAX_THREADS);

	if ( grid_y > 0 ) {
		const dim3 blockSize(MAX_THREADS,1, 1);
		const dim3 gridSize(grid_y,1,1);
		mult_kernel<<<gridSize,blockSize>>>(cuda_data->data,mult,MAX_THREADS);
	}

	if ( res_y > 0 ) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		mult_kernel<<<gridSize,blockSize>>>(cuda_data->data+grid_y*MAX_THREADS,mult,0);
	}

	cudaThreadSynchronize();
}

__global__ void add_kernel(float *data,const float add,const int num_threads)
{

	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	data[x+y*num_threads] += add;
}

void emdata_processor_add( EMDataForCuda* cuda_data, const float& add) {


	int num_calcs = cuda_data->nx*cuda_data->ny*cuda_data->nz;

	int grid_y = num_calcs/MAX_THREADS;
	int res_y = num_calcs - (grid_y*MAX_THREADS);

	if ( grid_y > 0 ) {
		const dim3 blockSize(MAX_THREADS,1, 1);
		const dim3 gridSize(grid_y,1,1);
		add_kernel<<<gridSize,blockSize>>>(cuda_data->data,add,MAX_THREADS);
	}

	if ( res_y > 0 ) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		add_kernel<<<gridSize,blockSize>>>(cuda_data->data+grid_y*MAX_THREADS,add,0);
	}

	cudaThreadSynchronize();
}

__global__ void assignment_kernel(float *data,const float value,const int num_threads)
{

	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	data[x+y*num_threads] = value;
}

void emdata_processor_to_value( EMDataForCuda* cuda_data, const float& value) {

	int num_calcs = cuda_data->nx*cuda_data->ny*cuda_data->nz;

	int grid_y = num_calcs/MAX_THREADS;
	int res_y = num_calcs - (grid_y*MAX_THREADS);

	if ( grid_y > 0 ) {
		const dim3 blockSize(MAX_THREADS,1, 1);
		const dim3 gridSize(grid_y,1,1);
		assignment_kernel<<<gridSize,blockSize>>>(cuda_data->data,value,MAX_THREADS);
	}

	if ( res_y > 0 ) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		assignment_kernel<<<gridSize,blockSize>>>(cuda_data->data+grid_y*MAX_THREADS,value,0);
	}

	cudaThreadSynchronize();
}


__global__ void phaseorigin_to_center_fourier(float* data, const int num_threads, const int nx, const int ny, const int nz, const int offset)
{
	const uint x=threadIdx.x;
	const uint y=blockIdx.x;
	// This optimization didn't seem to make much difference - the advantage was that shared memory was being used, essentially as a way force colesced memory accesses
//
// 	__shared__ float shared_data[MAX_THREADS];
//
// 	uint idx = x+y*num_threads+offset;
// 	shared_data[x] = data[idx];
// 	__syncthreads();
//
// 	const uint nxy = nx*ny;
// 	if ( idx % 4 == 0 ) {
// 		uint zz = idx/(nxy);
// 		uint yy = (idx-zz*nxy)/nx;
// 		uint xoff = ((yy+zz)%2==0?2:0);
// 		shared_data[x+xoff] *= -1;
// 		shared_data[x+xoff+1] *= -1;
//
// 	}
// 	__syncthreads();
// 	data[idx] = shared_data[x];

	uint idx = x+y*num_threads+offset;
	uint nxon4 = nx/4;
	const uint nxy = nxon4*ny;
	uint zz = idx/(nxy);
	uint yy = (idx-zz*nxy)/(nxon4);

	const uint xx = 4*(idx%(nxon4));

	const uint rnxy = nx*ny;
	const uint xoff = ((yy+zz)%2==0?2:0);
	const uint didx = zz*rnxy+yy*nx+xx+xoff;
	data[didx] *= -1;
	data[didx+1] *= -1;
}

void emdata_phaseorigin_to_center_fourier(const EMDataForCuda* cuda_data) {
	int nx = cuda_data->nx;
	int ny = cuda_data->ny;
	int nz = cuda_data->nz;
	float* data = cuda_data->data;

	if ( nx%2==0 && (ny%2==0 || ny==1 ) && (nz%2==0 || nz==1 ) ) {

		int num_calcs = nz*ny*(nx/4);

		int grid_y = num_calcs/(MAX_THREADS);
		int res_y = num_calcs - grid_y*MAX_THREADS;

		//int odd_offset=0;
		//if (((ny/2)%2)+((nz/2)%2)==1) odd_offset=1;
		if (grid_y > 0) {
			const dim3 blockSize(MAX_THREADS,1, 1);
			const dim3 gridSize(grid_y,1,1);
			phaseorigin_to_center_fourier<<<gridSize,blockSize>>>(data,MAX_THREADS,nx,ny,nz,0);
		}

		if (res_y > 0) {
			const dim3 blockSize(res_y,1, 1);
			const dim3 gridSize(1,1,1);
			phaseorigin_to_center_fourier<<<gridSize,blockSize>>>(data,MAX_THREADS,nx,ny,nz,grid_y*MAX_THREADS);
		}
		cudaThreadSynchronize();
	} else {
		throw;
	}
}

__global__ void correlation_kernel(float *ldata, float* rdata,const int num_threads)
{

	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	__shared__ float shared_rdata[MAX_THREADS];
	__shared__ float shared_ldata[MAX_THREADS];

	shared_ldata[x] = ldata[x+y*num_threads];
// 	__syncthreads(); // Not sure if this is necessary
	shared_rdata[x] = rdata[x+y*num_threads];
	__syncthreads();


	if (x % 2 == 0) {
		float v1 = shared_ldata[x];
		float v2 = shared_ldata[x+1];

		float u1 = shared_rdata[x];
		float u2 = shared_rdata[x+1];

		shared_ldata[x] = v1*u1 + v2*u2;
		shared_ldata[x+1] = v2*u1 - v1*u2;
	}
	__syncthreads();

	ldata[x+y*num_threads] = shared_ldata[x];
}

__global__ void auto_correlation_kernel(float *ldata, float* rdata,const int num_threads)
{

	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	__shared__ float shared_rdata[MAX_THREADS];
	__shared__ float shared_ldata[MAX_THREADS];

	shared_ldata[x] = ldata[x+y*num_threads];
	shared_rdata[x] = rdata[x+y*num_threads];

	__syncthreads();


	if (x % 2 == 0) {
		float v1 = shared_ldata[x];
		float v2 = shared_ldata[x+1];

		float u1 = shared_rdata[x];
		float u2 = shared_rdata[x+1];

		shared_ldata[x] = v1*u1 + v2*u2;
		shared_ldata[x+1] = 0;
	}
	__syncthreads();

	ldata[x+y*num_threads] = shared_ldata[x];
// 	const uint x=threadIdx.x;
// 	const uint y=blockIdx.x;
//
// 	const uint idx = 2*x + y*num_threads;
// 	const uint idxp1 = idx+1;
//
// 	const float v1 = ldata[idx];
// 	const float v2 = ldata[idxp1];
// 	const float u1 = rdata[idx];
// 	const float u2 = rdata[idxp1];
//
// 	ldata[idx] = v1*u1 + v2*u2;
// 	ldata[idxp1] = 0;
}

__global__ void correlation_kernel_texture_2D(float *ldata,const int num_threads,const int xsize,const int offset)
{
	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	__shared__ float shared_ldata[MAX_THREADS];

	shared_ldata[x] = ldata[x+y*num_threads+offset];

	__syncthreads();


	if (x % 2 == 0) {
		float v1 = shared_ldata[x];
		float v2 = shared_ldata[x+1];

		const uint idx = x + y*num_threads+offset;
// 		const uint idxp1 = idx+1;

		const uint tx = idx % xsize;
		const uint ty = idx / xsize;

		const float u1 = tex2D(tex2d,tx,ty);
		const float u2 =  tex2D(tex2d,tx+1,ty);

		shared_ldata[x] = v1*u1 + v2*u2;
		shared_ldata[x+1] = v2*u1 - v1*u2;
	}
	__syncthreads();

	ldata[x+y*num_threads+offset] = shared_ldata[x];
}


__global__ void correlation_kernel_texture_3D(float *ldata,const int num_threads, const int xsize, const int xysize, const int offset)
{

	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	__shared__ float shared_ldata[MAX_THREADS];

	shared_ldata[x] = ldata[x+y*num_threads+offset];

	__syncthreads();


	if (x % 2 == 0) {
		float v1 = shared_ldata[x];
		float v2 = shared_ldata[x+1];

		const uint idx = x + y*num_threads+offset;
// 		const uint idxp1 = idx+1;

		const uint tx = idx % xsize;
		const uint tz = idx / xysize;
		const uint ty = (idx - tz*xysize)/xsize;

		const float u1 = tex3D(tex,tx,ty,tz);
		const float u2 = tex3D(tex,tx+1,ty,tz);

		shared_ldata[x] = v1*u1 + v2*u2;
		shared_ldata[x+1] = v2*u1 - v1*u2;
	}
	__syncthreads();

	ldata[x+y*num_threads+offset] = shared_ldata[x];

// 	const uint idx = 2*x + y*num_threads + offset;
// 	const uint idxp1 = idx+1;
//
// 	const uint tx = idx % xsize;
// 	const uint tz = idx / xysize;
// 	const uint ty = (idx - tz*xysize)/xsize;
//
// 	const float v1 = ldata[idx];
// 	const float v2 = ldata[idxp1];
// 	const float u1 = tex3D(tex,tx,ty,tz);
// 	const float u2 = tex3D(tex,tx+1,ty,tz);
//
// 	ldata[idx] = v1*u1 + v2*u2;
// 	ldata[idxp1] = v2*u1 - v1*u2;
}

void emdata_processor_correlation_texture( const EMDataForCuda* cuda_data, const int center ) {
	int num_calcs = cuda_data->nx*cuda_data->ny*cuda_data->nz;

	int grid_y = num_calcs/(MAX_THREADS);
	int res_y = num_calcs - grid_y*MAX_THREADS;

// 	printf("Grid %d, Res %d, dims %d %d %d\n",grid_y,res_y,cuda_data->nx,cuda_data->ny,cuda_data->nz);

	if ( grid_y > 0 ) {
		const dim3 blockSize(MAX_THREADS,1, 1);
		const dim3 gridSize(grid_y,1,1);
		if (cuda_data->nz == 1) {
			correlation_kernel_texture_2D<<<gridSize,blockSize>>>(cuda_data->data,MAX_THREADS,cuda_data->nx,0);
		} else {
			correlation_kernel_texture_3D<<<gridSize,blockSize>>>(cuda_data->data,MAX_THREADS,cuda_data->nx,cuda_data->nx*cuda_data->ny,0);
		}
	}
// 	res_y = 0;
	if ( res_y > 0 ) {
		const dim3 blockSize(res_y,1,1);
		const dim3 gridSize(1,1,1);
		int inc = grid_y*MAX_THREADS;
// 		printf("Res %d, inc %d\n",res_y,inc);
		if (cuda_data->nz == 1) {
			correlation_kernel_texture_2D<<<gridSize,blockSize>>>(cuda_data->data,0,cuda_data->nx,inc);
		} else {
			correlation_kernel_texture_3D<<<gridSize,blockSize>>>(cuda_data->data,0,cuda_data->nx,cuda_data->nx*cuda_data->ny,inc);
		}
	}

	cudaThreadSynchronize();
	if (center) {
		emdata_phaseorigin_to_center_fourier(cuda_data);
	}
}


void emdata_processor_correlation( const EMDataForCuda* left, const EMDataForCuda* right, const int center) {

	int num_calcs = left->nx*left->ny*left->nz;

	int grid_y = num_calcs/(MAX_THREADS);
	int res_y = num_calcs - grid_y*MAX_THREADS;

	//printf("Grid y %d, res %d, dims %d %d %d\n", grid_y,res_y,left->nx,left->ny,left->nz);

	if ( grid_y > 0 ) {
		const dim3 blockSize(MAX_THREADS,1, 1);
		const dim3 gridSize(grid_y,1,1);
		if (left->data != right->data) {
			correlation_kernel<<<gridSize,blockSize>>>(left->data,right->data,MAX_THREADS);
		} else {
			auto_correlation_kernel<<<gridSize,blockSize>>>(left->data,right->data,MAX_THREADS);
		}
	}

	if ( res_y > 0 ) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		int inc = grid_y*MAX_THREADS;
		if (left->data != right->data) {
			correlation_kernel<<<gridSize,blockSize>>>(left->data+inc,right->data+inc,0);
		} else {
			auto_correlation_kernel<<<gridSize,blockSize>>>(left->data+inc,right->data+inc,0);
		}
	}
	cudaThreadSynchronize();

	if (center) {
		emdata_phaseorigin_to_center_fourier(left);
	}
}

__global__ void unwrap_kernel(float* dptr, const int num_threads, const int r1, const float p, const int nx, const int ny, const int nxp, const int dx,const int dy,const int weight_radial,const int offset) {
	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	const uint idx = x + y*num_threads+offset;

	const uint tx = idx % nxp;
	const uint ty = idx / nxp;

	float ang = tx * M_PI * p;
	float si = sinf(ang);
	float co = cosf(ang);

	float ypr1 = ty + r1;
	float xx = ypr1 * co + nx / 2 + dx;
	float yy = ypr1 * si + ny / 2 + dy;
	if ( weight_radial ) dptr[idx] = tex2D(tex2d,xx+0.5,yy+0.5)*ypr1;
	else dptr[idx] = tex2D(tex2d,xx+0.5,yy+0.5);
}


void emdata_unwrap(EMDataForCuda* data, int r1, int r2, int xs, int num_pi, int dx, int dy, int weight_radial, int nx, int ny) {

// 	float* dptr;
	int n = xs*(r2-r1);
// 	cudaError_t error = cudaMalloc((void**)&dptr,n*sizeof(float));
// 	if ( error != cudaSuccess ) {
// 		const char* s = cudaGetErrorString(error);
// 		printf("Cuda malloc failed in emdata_unwrap: %s\n",s);
// 		throw;
// 	}
	int num_calcs = n;

	int grid_y = num_calcs/(MAX_THREADS);
	int res_y = num_calcs - grid_y*MAX_THREADS;

	//printf("Grid %d, res %d, n %d, p %f \n",grid_y,res_y,n, p/xs);

	if ( grid_y > 0 ) {
		const dim3 blockSize(MAX_THREADS,1, 1);
		const dim3 gridSize(grid_y,1,1);
		unwrap_kernel<<<gridSize,blockSize>>>(data->data,MAX_THREADS,r1,(float) num_pi/ (float)xs, nx,ny,xs,dx,dy,weight_radial,0);
	}

	if ( res_y > 0 ) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		unwrap_kernel<<<gridSize,blockSize>>>(data->data,MAX_THREADS,r1, (float) num_pi/ (float)xs, nx,ny,xs,dx,dy,weight_radial,grid_y*MAX_THREADS);
	}

// 	EMDataForCuda* tmp = (EMDataForCuda*) malloc( sizeof(EMDataForCuda) );
// 	tmp->data = dptr;
// 	tmp->nx = xs;
// 	tmp->ny = r2-r1;
// 	tmp->nz = 1;
// 	return tmp;
}



__global__ void swap_bot_left_top_right(float* data, const int num_threads, const int nx, const int ny, const int xodd, const int yodd, const int offset) {
	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	const uint gpu_idx = x+y*num_threads+offset;
	const uint c = gpu_idx % (nx/2);
	const uint r = gpu_idx / (nx/2);

	const uint idx1 = r*nx + c;
	const uint idx2 = (r+ny/2+yodd)*nx + c + nx/2+xodd;
	float tmp = data[idx1];
	data[idx1] = data[idx2];
	data[idx2] = tmp;
}

__global__ void swap_top_left_bot_right(float* data, const int num_threads, const int nx, const int ny, const int xodd, const int yodd, const int offset) {
	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	const uint gpu_idx = x+y*num_threads+offset;
	const uint c = gpu_idx % (nx/2);
	const uint r = gpu_idx / (nx/2) + ny/2+yodd;

	const uint idx1 = r*nx + c;
	const uint idx2 = (r-ny/2-yodd)*nx + c + nx/2+xodd;
	float tmp = data[idx1];
	data[idx1] = data[idx2];
	data[idx2] = tmp;
}

__global__ void swap_middle_row(float* data, const int num_threads, const int nx, const int ny, const int xodd, const int yodd, const int offset) {
	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	const uint c = x+y*num_threads+offset;
	int r = ny/2;
	int idx1 = r*nx + c;
	int idx2 = r*nx + c + nx/2+ xodd;
	float tmp = data[idx1];
	data[idx1] = data[idx2];
	data[idx2] = tmp;
}

// Iterate along the central column, swapping values where appropriate
__global__ void swap_middle_column(float* data, const int num_threads, const int nx, const int ny, const int xodd, const int yodd, const int offset) {
	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	const uint r = x+y*num_threads+offset;
	int c = nx/2;
	int idx1 = r*nx + c;
	int idx2 = (r+ny/2+yodd)*nx + c;
	float tmp = data[idx1];
	data[idx1] = data[idx2];
	data[idx2] = tmp;
}

void swap_central_slices_180(EMDataForCuda* cuda_data)
{
	int nx = cuda_data->nx;
	int ny = cuda_data->ny;
	int nz = cuda_data->nz;

	int xodd = (nx % 2) == 1;
	int yodd = (ny % 2) == 1;
	//int zodd = (nz % 2) == 1;

	//int nxy = nx * ny;
	float *data = cuda_data->data;

	if ( ny == 1 && nz == 1 ){
		throw;
	}
	else if ( nz == 1 ) {
		if ( yodd ) {
			// Iterate along middle row, swapping values where appropriate
			int num_calcs = nx/2;

			int grid_y = num_calcs/(MAX_THREADS);
			int res_y = num_calcs - grid_y*MAX_THREADS;

			if (grid_y > 0) {
				const dim3 blockSize(MAX_THREADS,1, 1);
				const dim3 gridSize(grid_y,1,1);
				swap_middle_row<<<gridSize,blockSize>>>(data,MAX_THREADS,nx,ny,xodd,yodd,0);
			}

			if (res_y) {
				const dim3 blockSize(res_y,1, 1);
				const dim3 gridSize(1,1,1);
				swap_middle_row<<<gridSize,blockSize>>>(data,MAX_THREADS,nx,ny,xodd,yodd,grid_y*MAX_THREADS);
			}
		}

		if ( xodd )	{
			// Iterate along the central column, swapping values where appropriate
			int num_calcs = ny/2;

			int grid_y = num_calcs/(MAX_THREADS);
			int res_y = num_calcs - grid_y*MAX_THREADS;

			if (grid_y > 0) {
				const dim3 blockSize(MAX_THREADS,1, 1);
				const dim3 gridSize(grid_y,1,1);
				swap_middle_column<<<gridSize,blockSize>>>(data,MAX_THREADS,nx,ny,xodd,yodd,0);
			}

			if (res_y) {
				const dim3 blockSize(res_y,1, 1);
				const dim3 gridSize(1,1,1);
				swap_middle_column<<<gridSize,blockSize>>>(data,MAX_THREADS,nx,ny,xodd,yodd,grid_y*MAX_THREADS);
			}

		}
	}
	else // nx && ny && nz are greater than 1
	{
		throw;
	}
}

void swap_corners_180(EMDataForCuda* cuda_data)
{
	int nx = cuda_data->nx;
	int ny = cuda_data->ny;
	int nz = cuda_data->nz;

	int xodd = (nx % 2) == 1;
	int yodd = (ny % 2) == 1;
	//int zodd = (nz % 2) == 1;

	//int nxy = nx * ny;

	float *data = cuda_data->data;

	if ( ny == 1 && nz == 1 ){
		throw;
	}
	else if ( nz == 1 ) {
		int num_calcs = ny/2*nx/2;

		int grid_y = num_calcs/(MAX_THREADS);
		int res_y = num_calcs - grid_y*MAX_THREADS;

		//printf("Grid %d, res %d, n %d\n",grid_y,res_y,num_calcs );
		// Swap bottom left and top right
		if (grid_y > 0) {
			const dim3 blockSize(MAX_THREADS,1, 1);
			const dim3 gridSize(grid_y,1,1);
			swap_bot_left_top_right<<<gridSize,blockSize>>>(data,MAX_THREADS,nx,ny,xodd,yodd,0);
		}

		if (res_y) {
			const dim3 blockSize(res_y,1, 1);
			const dim3 gridSize(1,1,1);
			swap_bot_left_top_right<<<gridSize,blockSize>>>(data,MAX_THREADS,nx,ny,xodd,yodd,grid_y*MAX_THREADS);
		}

		num_calcs = (ny-ny/2+yodd)*nx/2;
		//printf("Grid %d, res %d, n %d\n",grid_y,res_y,num_calcs );

		grid_y = num_calcs/(MAX_THREADS);
		res_y = num_calcs - grid_y*MAX_THREADS;
		// Swap the top left and bottom right corners
		if (grid_y > 0) {
			const dim3 blockSize(MAX_THREADS,1, 1);
			const dim3 gridSize(grid_y,1,1);
			swap_top_left_bot_right<<<gridSize,blockSize>>>(data,MAX_THREADS,nx,ny,xodd,yodd,0);
		}

		if (res_y) {
			const dim3 blockSize(res_y,1, 1);
			const dim3 gridSize(1,1,1);
			swap_top_left_bot_right<<<gridSize,blockSize>>>(data,MAX_THREADS,nx,ny,xodd,yodd,grid_y*MAX_THREADS);

		}
	}
	else // nx && ny && nz are greater than 1
	{
		throw;
	}
}

__global__ void middle_to_right(float* data, const int nx, const int ny)
{
	float tmp;
	for ( int r  = 0; r < ny; ++r ) {
		float last_val = data[r*nx+nx/2];
		for ( int c = nx-1; c >=  nx/2; --c ){
			int idx = r*nx+c;
			tmp = data[idx];
			data[idx] = last_val;
			last_val = tmp;
		}
	}
}

__global__ void middle_to_top(float* data, const int nx, const int ny)
{
	float tmp;
	for ( int c = 0; c < nx; ++c ) {
		// Get the value in the top row
		float last_val = data[ny/2*nx + c];
		for ( int r = ny-1; r >= ny/2; --r ){
			int idx = r*nx+c;
			tmp = data[idx];
			data[idx] = last_val;
			last_val = tmp;
		}
	}
}


void emdata_phaseorigin_to_center(EMDataForCuda* cuda_data) {
	int xodd = (cuda_data->nx % 2) == 1;
	int yodd = (cuda_data->ny % 2) == 1;
	//int zodd = (cuda_data->nz % 2) == 1;

	//int nxy = nx * ny;
	if ( cuda_data->nz == 1 && cuda_data->ny > 1 ){
		// The order in which these operations occur literally undoes what the
		// PhaseToCornerProcessor did to the image.
		// First, the corners sections of the image are swapped appropriately
		swap_corners_180(cuda_data);
		// Second, central pixel lines are swapped
		swap_central_slices_180(cuda_data);

		// Third, appropriate sections of the image are cyclically shifted by one pixel
		if (xodd) {
			// Transfer the middle column to the far right
			// Shift all from the far right to (but not including the) middle one to the left
			middle_to_right<<<1,1>>>(cuda_data->data,cuda_data->nx,cuda_data->ny);
		}
		if (yodd) {
			// Tranfer the middle row to the top,
			// shifting all pixels from the top row down one, until  but not including the) middle
			middle_to_top<<<1,1>>>(cuda_data->data,cuda_data->nx,cuda_data->ny);
		}
		cudaThreadSynchronize();
	} else {
		throw;
	}
}

__global__ void transform_kernel_3D(float *out,int nx,int ny,int nz,int num_threads, float3 mxx,float3 mxy, float3 mxz, float3 trans,int offset)
{
	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	const uint idx = x + y*num_threads + offset;

	const uint nxy = nx*ny;
	const float fx = idx % nx - nx/2.0;
	float fz = idx / nxy;
	const float fy = (idx - ((int)fz)*nxy)/nx - ny/2.0;
	fz -= nz/2.0;

	// The 0.5f offsets for x,y and z are required - Read section D.2 in Appendix D of the CUDA
	// Programming Guide (version 2.0).
	// Thankyou http://sites.google.com/site/cudaiap2009/cookbook-1
	float tx=fx*mxx.x+fy*mxx.y+fz*mxx.z+nx/2.0+trans.x+0.5;
	float ty=fx*mxy.x+fy*mxy.y+fz*mxy.z+ny/2.0+trans.y+0.5;
	float tz=fx*mxz.x+fy*mxz.y+fz*mxz.z+nz/2.0+trans.z+0.5;

	out[idx]=tex3D(tex, tx,ty,tz);
}

__global__ void transform_kernel_2D(float *out,int nx,int ny,int num_threads, float3 mxx,float3 mxy, float3 mxz, float3 trans,int offset)
{
	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	const uint idx = x + y*num_threads + offset;

	const float fx = idx % nx - nx/2.0;
	const float fy = idx/nx - ny/2.0;


	// The 0.5f offsets for x,y and z are required - Read section D.2 in Appendix D of the CUDA
	// Programming Guide (version 2.0).
	// Thankyou http://sites.google.com/site/cudaiap2009/cookbook-1
	float tx=fx*mxx.x+fy*mxx.y+nx/2.0+trans.x+0.5;
	float ty=fx*mxy.x+fy*mxy.y+ny/2.0+trans.y+0.5;

	out[idx]=tex2D(tex2d, tx,ty);
}

EMDataForCuda* emdata_transform_cuda(const float* const matrix,const int nx,const int ny,const int nz) {
	EMDataForCuda* t = (EMDataForCuda*) malloc(sizeof(EMDataForCuda));

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

	int num_calcs = ny*nx*nz;

	int grid_y = num_calcs/(MAX_THREADS);
	int res_y = num_calcs - grid_y*MAX_THREADS;

	cudaError_t error = cudaMalloc((void**)&t->data,nx*ny*nz*sizeof(float));
	if ( error != cudaSuccess ) {
		return 0; //Calling function should know something went wrong
	}

        //debugging stuff
        //printf("CUDA ptr = %d\n", t->data);
 
	//int asize = 10*sizeof(float);
	//float * rdatax;
	//float * ddata;
	//rdatax = (float*) malloc(asize);

        //cudaMalloc((void**)&ddata,asize);
        //cudaError_t cerror = cudaMemcpy(rdatax,ddata,asize,cudaMemcpyDeviceToHost);
	//if (cerror != cudaSuccess ){printf("CUDA failed at low level\n");}
	//printf("Success %d, error %d\n", cudaSuccess, cerror);

	if (grid_y > 0) {
		const dim3 blockSize(MAX_THREADS,1, 1);
		const dim3 gridSize(grid_y,1,1);
		if ( nz > 1 ) {
			transform_kernel_3D<<<gridSize,blockSize>>>(t->data,nx,ny,nz,MAX_THREADS,mxx,mxy,mxz,trans,0);
		}
		else if ( ny > 1 ) {
			transform_kernel_2D<<<gridSize,blockSize>>>(t->data,nx,ny,MAX_THREADS,mxx,mxy,mxz,trans,0);
		} else throw;
	}

	if (res_y) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		if ( nz > 1 ) {
			transform_kernel_3D<<<gridSize,blockSize>>>(t->data,nx,ny,nz,MAX_THREADS,mxx,mxy,mxz,trans,grid_y*MAX_THREADS);
		}
		else if ( ny > 1 ) {
			transform_kernel_2D<<<gridSize,blockSize>>>(t->data,nx,ny,MAX_THREADS,mxx,mxy,mxz,trans,grid_y*MAX_THREADS);
		} else throw;
	}
	cudaThreadSynchronize();

	t->nx = nx;
	t->ny = ny;
	t->nz = nz;

	return t;
}

__global__ void ri2inten_kernel(float *data,const int num_threads)
{

	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	__shared__ float shared_data[MAX_THREADS];

	shared_data[x] = data[x+y*num_threads];

	__syncthreads();

	if (x % 2 == 0) {
		float a = shared_data[x];
		float b = shared_data[x+1];
		shared_data[x]= a*a + b*b;
		shared_data[x+1] = 0;
	}

	__syncthreads();

	data[x+y*num_threads] = shared_data[x];
}

void emdata_ri2inten( EMDataForCuda* cuda_data) {

	int num_calcs = cuda_data->nx*cuda_data->ny*cuda_data->nz;

	int grid_y = num_calcs/MAX_THREADS;
	int res_y = num_calcs - (grid_y*MAX_THREADS);

	if ( grid_y > 0 ) {
		const dim3 blockSize(MAX_THREADS,1, 1);
		const dim3 gridSize(grid_y,1,1);
		ri2inten_kernel<<<gridSize,blockSize>>>(cuda_data->data,MAX_THREADS);
	}

	if ( res_y > 0 ) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		ri2inten_kernel<<<gridSize,blockSize>>>(cuda_data->data+grid_y*MAX_THREADS,0);
	}

	cudaThreadSynchronize();
}



__global__ void ap2ri_kernel(float *data,const int num_threads)
{

	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	__shared__ float shared_data[MAX_THREADS];

	shared_data[x] = data[x+y*num_threads];

	__syncthreads();

	if (x % 2 == 0) {
		float a = shared_data[x];
		float b = shared_data[x+1];

		shared_data[x+1]= a*sin(b);
		shared_data[x] =a*cos(b);
	}

	__syncthreads();

	data[x+y*num_threads] = shared_data[x];
}


void emdata_ap2ri( EMDataForCuda* cuda_data) {

	int num_calcs = cuda_data->nx*cuda_data->ny*cuda_data->nz;

	int grid_y = num_calcs/MAX_THREADS;
	int res_y = num_calcs - (grid_y*MAX_THREADS);

	if ( grid_y > 0 ) {
		const dim3 blockSize(MAX_THREADS,1, 1);
		const dim3 gridSize(grid_y,1,1);
		ap2ri_kernel<<<gridSize,blockSize>>>(cuda_data->data,MAX_THREADS);
	}

	if ( res_y > 0 ) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		ap2ri_kernel<<<gridSize,blockSize>>>(cuda_data->data+grid_y*MAX_THREADS,0);
	}

	cudaThreadSynchronize();
}


__global__ void ri2ap_kernel(float *data,const int num_threads)
{

	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	__shared__ float shared_data[MAX_THREADS];

	shared_data[x] = data[x+y*num_threads];

	__syncthreads();

	if (x % 2 == 0) {
		float a = shared_data[x];
		float b = shared_data[x+1];

		float amp = hypotf(a,b);
		float phase = atan2(b, a);
		shared_data[x] = amp;
		shared_data[x+1] = phase;
	}

	__syncthreads();

	data[x+y*num_threads] = shared_data[x];
}

void emdata_ri2ap( EMDataForCuda* cuda_data) {

	int num_calcs = cuda_data->nx*cuda_data->ny*cuda_data->nz;

	int grid_y = num_calcs/MAX_THREADS;
	int res_y = num_calcs - (grid_y*MAX_THREADS);

	if ( grid_y > 0 ) {
		const dim3 blockSize(MAX_THREADS,1, 1);
		const dim3 gridSize(grid_y,1,1);
		ri2ap_kernel<<<gridSize,blockSize>>>(cuda_data->data,MAX_THREADS);
	}

	if ( res_y > 0 ) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		ri2ap_kernel<<<gridSize,blockSize>>>(cuda_data->data+grid_y*MAX_THREADS,0);
	}

	cudaThreadSynchronize();
}

__global__ void binarize_fourier_kernel(float *data,const int num_threads, const float threshold)
{
	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	__shared__ float shared_data[MAX_THREADS];

	shared_data[x] = data[x+y*num_threads];

	__syncthreads();

	if (x % 2 == 0) {
		float a = shared_data[x];
		shared_data[x+1] = 0;
		if ( a >= threshold ) {
			shared_data[x] = 1;
		} else {
			shared_data[x] = 0;
		}
	}

	__syncthreads();

	data[x+y*num_threads] = shared_data[x];
}


void binarize_fourier_amp_processor(EMDataForCuda* cuda_data,const float& threshold)
{
	int num_calcs = cuda_data->nx*cuda_data->ny*cuda_data->nz;

	int grid_y = num_calcs/MAX_THREADS;
	int res_y = num_calcs - (grid_y*MAX_THREADS);

	if ( grid_y > 0 ) {
		const dim3 blockSize(MAX_THREADS,1, 1);
		const dim3 gridSize(grid_y,1,1);
		binarize_fourier_kernel<<<gridSize,blockSize>>>(cuda_data->data,MAX_THREADS,threshold);
	}

	if ( res_y > 0 ) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		binarize_fourier_kernel<<<gridSize,blockSize>>>(cuda_data->data+grid_y*MAX_THREADS,MAX_THREADS,threshold);
	}

	cudaThreadSynchronize();
}


__global__ void rotate_180( float* data,int nx, int nxy, int offset, unsigned int size) {

	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	__shared__ float shared_lower_data[MAX_THREADS];
	__shared__ float shared_upper_data[MAX_THREADS];

	shared_lower_data[x] = data[x+y*MAX_THREADS+offset];
	shared_upper_data[x] = data[nxy + x+(-y-1)*MAX_THREADS-offset];
	__syncthreads();


	if (size == 0) {
		float tmp = shared_lower_data[x];
		shared_lower_data[x] = shared_upper_data[MAX_THREADS-x-1];
		shared_upper_data[MAX_THREADS-x-1] = tmp;
	} else {
		if ( x < size ) {
			float tmp = shared_lower_data[x];
			shared_lower_data[x] = shared_upper_data[MAX_THREADS-x-1];
			shared_upper_data[MAX_THREADS-x-1]= tmp;

		}
	}

	__syncthreads();
	if (size == 0) {
		data[x+y*MAX_THREADS+offset] = shared_lower_data[x];
		data[nxy+x+(-y-1)*MAX_THREADS-offset] = shared_upper_data[x];
	} else {
		if ( x < size ) {
			data[nxy-x-1+(-y)*MAX_THREADS-offset] = shared_upper_data[MAX_THREADS-x-1];
			data[x+y*MAX_THREADS+offset] = shared_lower_data[x];
		}
	}

}

 void emdata_rotate_180( EMDataForCuda* cuda_data) {

	// Only works for 2D images
	int nx = cuda_data->nx;
	int ny = cuda_data->ny;
	int nxy = nx*ny;
	// Half of the pixels in the image are swapped with the pixels of the other half, hence the divide by two
	int num_mem_swaps = 0;
	int offset = 0;
	if (nx%2 == 0 && ny%2 == 0) {
		num_mem_swaps = (nx-2)*(ny+1); // If you figure it out, it's (nx-2)*ny + nx-1 -1
		offset = nx+1;
	} else if ( nx % 2 == 0 ) {
		num_mem_swaps = nx*ny -2;
		offset = 1;
	}else if ( ny % 2 == 0 ) {
		num_mem_swaps = (nx-1)*ny-1;
		offset =  cuda_data->nx;
	}else {
		num_mem_swaps = nx*ny-1;
		offset =  0;
	}

	int num_calcs = num_mem_swaps/2;
	nxy -= offset;

	unsigned int grid_y = num_calcs/MAX_THREADS;
	unsigned int res_y = num_calcs - (grid_y*MAX_THREADS);

	if ( grid_y > 0 ) {
		const dim3 blockSize(MAX_THREADS,1, 1);
		const dim3 gridSize(grid_y,1,1);
		rotate_180<<<gridSize,blockSize>>>(cuda_data->data+offset,nx,nxy,0,0);
	}

	if ( res_y > 0 ) {
		const dim3 blockSize(MAX_THREADS,1, 1);
		const dim3 gridSize(1,1,1);
		rotate_180<<<gridSize,blockSize>>>(cuda_data->data+offset,nx,nxy,grid_y*MAX_THREADS,res_y);
	}

	cudaThreadSynchronize();
 }


