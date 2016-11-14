

typedef unsigned int uint;

#ifdef WIN32
	#define M_PI 3.14159265358979323846f
#endif	//WIN32

__global__ void transform_kernel_3D(float *out,int nx,int ny,int nz, float3 mxx,float3 mxy, float3 mxz, float3 trans, int totaltc)
{
	const uint idx = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*MAX_THREADS;

	if (idx < totaltc){
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

		out[idx]=tex3D(texA, tx,ty,tz);
	}
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

	out[idx]=tex2D(texA2d, tx,ty);
}

float* emdata_transform_cuda(const float* const matrix,const int nx,const int ny,const int nz) {

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

	float* data;
	cudaError_t error = cudaMalloc((void**)&data,nx*ny*nz*sizeof(float));
	if ( error != cudaSuccess ) {
		return 0; //Calling function should know something went wrong
	}

	if ( nz > 1 ) {
		int grid = int(ceil(sqrt(nx*ny*nz/MAX_THREADS)));
		const dim3 blockSize(MAX_THREADS,1, 1);
		const dim3 gridSize(grid,grid,1);
		transform_kernel_3D<<<gridSize,blockSize>>>(data,nx,ny,nz,mxx,mxy,mxz,trans,nx*ny*nz);
		cudaThreadSynchronize();
		
		return data;
	}
		
	if (grid_y > 0) {
		const dim3 blockSize(MAX_THREADS,1, 1);
		const dim3 gridSize(grid_y,1,1);
		transform_kernel_2D<<<gridSize,blockSize>>>(data,nx,ny,MAX_THREADS,mxx,mxy,mxz,trans,0);
	}

	if (res_y) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		transform_kernel_2D<<<gridSize,blockSize>>>(data,nx,ny,MAX_THREADS,mxx,mxy,mxz,trans,grid_y*MAX_THREADS);
	}
	cudaThreadSynchronize();

	return data;
}

__global__ void complexmult_conj_kernal(float *afft, const float *bfft, int totaltc)
{

	const uint ridx = 2*(threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*MAX_THREADS);
	
	//maybe use float2 to improve coalessing....
	if (ridx < totaltc){
		const uint iidx = ridx + 1;
		float afftr = afft[ridx];
		float affti = afft[iidx];
		float bfftr = bfft[ridx];
		float bffti = bfft[iidx];
	
		afft[ridx] = afftr*bfftr + affti*bffti;  //real portion
		afft[iidx] = affti*bfftr - afftr*bffti; //imaginary portion
	}

}

void calc_ccf_cuda(float* afft, const float* bfft, const int nx, const int ny, const int nz)
{

	int grid = int(ceil(sqrt((nx/2)*ny*nz/MAX_THREADS)));
	
	const dim3 blockSize(MAX_THREADS,1, 1);
	const dim3 gridSize(grid,grid,1);
	complexmult_conj_kernal<<<gridSize,blockSize>>>(afft, bfft, nx*ny*nz);

	cudaThreadSynchronize();

}

__global__ void complexmult_kernal(float *afft, const float *bfft, int totaltc)
{

	const uint ridx = 2*(threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*MAX_THREADS);
	
	if(ridx < totaltc){
		const uint iidx = ridx + 1;
		//maybe use float2 to improve coalessing....

		float afftr = afft[ridx];
		float affti = afft[iidx];
		float bfftr = bfft[ridx];
		float bffti = bfft[iidx];

		afft[ridx] = afftr*bfftr - affti*bffti;  //real portion
		afft[iidx] = affti*bfftr + afftr*bffti; //imaginary portion
	}

}

void calc_conv_cuda(float* afft, const float* bfft, const int nx, const int ny, const int nz)
{

	int grid = int(ceil(sqrt((nx/2)*ny*nz/MAX_THREADS)));

	const dim3 blockSize(MAX_THREADS,1, 1);
	const dim3 gridSize(grid,grid,1);
	complexmult_kernal<<<gridSize,blockSize>>>(afft, bfft, nx*ny*nz);

	cudaThreadSynchronize();

}

// This is a bit rubbish way of doing things. This function is not really suited to CUDA the only reason I have it here is to prevent
// copying back to host. See below for an improved version...
__global__ void  calc_max_location_wrap(const float* in, CudaPeakInfo* soln, const int nx, const int ny, const int nz, const int maxdx, const int maxdy, const int maxdz) {
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
				float value = in[ii+jj*nx+kk*nxy];

				if (value > max_value) {
					max_value = value;
					soln->px = i;
					soln->py = j;
					soln->pz = k;
					soln->peak = max_value;
				}
			}
		}
	}
}

//This is basically the above kernal, but for 2D images
__global__ void calc_max_location_wrapv2(const float* in, CudaPeakInfo* soln, const int nx, const int ny, const int nz, const int maxdx, const int maxdy, const int maxdz)
{
	int maxshiftx = maxdx, maxshifty = maxdy;
	if (maxdx == -1) maxshiftx = nx/4;
	if (maxdy == -1) maxshifty = ny/4;

	float max_value = -10000000000000;
	for (int j = -maxshifty; j <= maxshifty; j++) {
		for (int i = -maxshiftx; i <= maxshiftx; i++) {
				
			int jj = j;
			if (jj < 0) {
				jj = ny+jj;
			}
				
			int ii = i;
			if (ii < 0) {
				ii = nx+ii;
			}
			float value = in[ii+jj*nx];

			if (value > max_value) {
				max_value = value;
				soln->px = i;
				soln->py = j;
				soln->pz = 0;
				soln->peak = max_value;
			}
		}
	}

	
}

//so far only for 2D, and this needs some tweaking to get it to work, currently we are doing it the dumb way (see above)
// But then again, the real bottleneck, for CUDA, it the 2D FFTs. To get a significant speedup for 2D images, we need to
// process images in parallel.
__global__ void calc_max_location_wrapv3(const float* in, CudaPeakInfo* solns, const int num_threads, const int nx, const int maxdx, const int maxdy, const int maxdz, const int offset)
{

	int idx = threadIdx.x + blockIdx.x*num_threads + offset;
	
	float bestval = -10000000000.0f;
	int yidx = idx % maxdx + (idx/maxdx)*(nx-maxdx);
	for(int i = -maxdx; i < maxdx; i ++)
	{	
		//int iabs = (i+(i>>31))^(i>>31); //find abs w/o branching
		int index = yidx*nx + i;
		float val = in[index];
		if(val > bestval)
		{	 
			bestval = val;
			solns[idx].px = i;
			solns[idx].py = yidx - nx*(2*yidx/nx) - (nx - 1 - i)/nx; //works unless search range is full image
			solns[idx].pz = index;
			solns[idx].peak = bestval;
		}
	}
	
}

CudaPeakInfo* calc_max_location_wrap_cuda(const float* in, const int nx, const int ny, const int nz, const int maxdx, const int maxdy, const int maxdz)
{
	
	if(nz > 1)
	{
		const dim3 blockSize(1,1,1);
		const dim3 gridSize(1,1,1);
	
		CudaPeakInfo * device_soln=0;
		cudaMalloc((void **)&device_soln, sizeof(CudaPeakInfo));
		CudaPeakInfo * host_soln = 0;
		host_soln = (CudaPeakInfo*) malloc(sizeof(CudaPeakInfo));

		calc_max_location_wrap<<<blockSize,gridSize>>>(in,device_soln, nx, ny, nz, maxdx, maxdy, maxdz);
		cudaThreadSynchronize();
		cudaMemcpy(host_soln,device_soln,sizeof(CudaPeakInfo),cudaMemcpyDeviceToHost);
		cudaFree(device_soln);

		return host_soln;
	}else{
		//printf("2D\n");

		const dim3 blockSize(1,1,1);
		const dim3 gridSize(1,1,1);
	
		CudaPeakInfo * device_soln=0;
		cudaMalloc((void **)&device_soln, sizeof(CudaPeakInfo));
		CudaPeakInfo * host_soln = 0;
		host_soln = (CudaPeakInfo*) malloc(sizeof(CudaPeakInfo));

		calc_max_location_wrapv2<<<blockSize,gridSize>>>(in,device_soln, nx, ny, nz, maxdx, maxdy, maxdz);
		cudaThreadSynchronize();
		cudaMemcpy(host_soln,device_soln,sizeof(CudaPeakInfo),cudaMemcpyDeviceToHost);
		cudaFree(device_soln);

		return host_soln;
		
		/* This code needs to be fixed......(In theory it should be faster)
		int num_calcs = 2*maxdx;
		int grid_y = num_calcs/(MAX_THREADS);
		int res_y = num_calcs - grid_y*MAX_THREADS;

		CudaPeakInfo * device_soln=0;
		cudaMalloc((void **)&device_soln, num_calcs*sizeof(CudaPeakInfo));
		CudaPeakInfo * host_soln = 0;
		host_soln = (CudaPeakInfo*) malloc(num_calcs*sizeof(CudaPeakInfo));
	
		if (grid_y > 0) {
			const dim3 blockSize(MAX_THREADS,1,1);
			const dim3 gridSize(grid_y,1,1);
			calc_max_location_wrapv3<<<gridSize, blockSize>>>(in,device_soln, MAX_THREADS, nx, maxdx, maxdy, maxdz, 0);
		}else{
			const dim3 blockSize(num_calcs,1,1);
			const dim3 gridSize(1,1,1);
			calc_max_location_wrapv3<<<gridSize, blockSize>>>(in,device_soln, MAX_THREADS, nx, maxdx, maxdy, maxdz, 0);
		}
		if(res_y){
			const dim3 blockSize(res_y,1,1);
			const dim3 gridSize(1,1,1);
			calc_max_location_wrapv3<<<gridSize, blockSize>>>(in,device_soln, MAX_THREADS, nx, maxdx, maxdy, maxdz, grid_y*MAX_THREADS);
		}

		cudaThreadSynchronize();

		cudaMemcpy(host_soln,device_soln,num_calcs*sizeof(CudaPeakInfo),cudaMemcpyDeviceToHost);
		cudaFree(device_soln);

		float bestpeak = -1000000000000.0f;
		CudaPeakInfo * soln = (CudaPeakInfo*) malloc(sizeof(CudaPeakInfo));
		for(int i = 0; i < num_calcs; i++)
		{
			if(host_soln[i].peak > bestpeak)
			{
				bestpeak = host_soln[i].peak;
				soln->px = host_soln[i].px;
				soln->py = host_soln[i].py;
				soln->pz = 0;
				soln->peak = bestpeak;
			}	
			//printf("x %d y %d z %d peak %f\n",host_soln[i].px,host_soln[i].py,host_soln[i].pz,host_soln[i].peak);
		}
		free(host_soln);
		return soln;
		*/

	}
}

__device__ float dget_value_at_wrap_cuda(const float* data, const int nx, const int ny, const int nz, int x, int y, int z)
{
	
	if (x < 0) x = nx + x;
	if (y < 0) y = ny + y;
	if (z < 0) z = nz + z;

	return data[x + y*nx + z*nx*ny];

}

__global__ void  calc_max_location_wrap_intp(const float* in, CudaPeakInfoFloat* soln, const int nx, const int ny, const int nz, const int maxdx, const int maxdy, const int maxdz) {
	int maxshiftx = maxdx, maxshifty = maxdy, maxshiftz = maxdz;
	if (maxdx == -1) maxshiftx = nx/4;
	if (maxdy == -1) maxshifty = ny/4;
	if (maxdz == -1) maxshiftz = nz/4;

	float max_value = -10000000000000;
	int nxy = nx*ny;

	int px = 0;
	int py = 0;
	int pz = 0;
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
				float value = in[ii+jj*nx+kk*nxy];

				if (value > max_value) {
					max_value = value;
					px = i;
					py = j;
					pz = k;
					soln->peak = max_value;
				}
			}
		}
	}

	float cmx = 0.0f; float cmy = 0.0f; float cmz = 0.0f;
	float sval = 0.0f;
	for (float x = float(px)-2.0f; x <= float(px)+2.0f; x++) {
		for (float y = float(py)-2.0f; y <= float(py)+2.0f; y++) {
			for (float z = float(pz)-2.0f; z <= float(pz)+2.0f; z++) {
				//Compute center of mass
				float val = dget_value_at_wrap_cuda(in,nx,ny,nz,x,y,z);
				cmx += x*val;
				cmy += y*val;
				cmz += z*val;
				sval += val;
			}
		}
	}
	soln->xintp = cmx/sval;
	soln->yintp = cmy/sval;
	soln->zintp = cmz/sval;

/**
	// Now do the intepolation...
	float x2 = float(px);
	float x1 = x2-1.0f;
	float x3 = x2+1.0f;
	float y2 = float(py);
	float y1 = y2-1.0f;
	float y3 = y2+1.0f;
	float z2 = float(pz);
	float z1 = z2-1.0f;
	float z3 = z2+1.0f;

	float yx1 = dget_value_at_wrap_cuda(in, nx, ny, nz, x1 , y2, z2);
	float yx2 = dget_value_at_wrap_cuda(in, nx, ny, nz, x2 , y2, z2);
	float yx3 = dget_value_at_wrap_cuda(in, nx, ny, nz, x3 , y2, z2);
	float yy1 = dget_value_at_wrap_cuda(in, nx, ny, nz, x2 , y1, z2);
	float yy2 = dget_value_at_wrap_cuda(in, nx, ny, nz, x2 , y2, z2);
	float yy3 = dget_value_at_wrap_cuda(in, nx, ny, nz, x2 , y3, z2);
	float yz1 = dget_value_at_wrap_cuda(in, nx, ny, nz, x2 , y2, z1);
	float yz2 = dget_value_at_wrap_cuda(in, nx, ny, nz, x2 , y2, z2);
	float yz3 = dget_value_at_wrap_cuda(in, nx, ny, nz, x2 , y2, z3);	

	// Fit peak in X to y = ax^2 + bx +c
	float bx = ((yx1 - yx2)*(x2*x2 - x3*x3)/(x1*x1 - x2*x2) - (yx2-yx3))/(-(x2 - x3) + (x1 - x2)*(x2*x2 - x3*x3)/(x1*x1 - x2*x2));
	float ax = ((yx1 - yx2) - bx*(x1 - x2))/(x1*x1 - x2*x2);
	//Find minima
	soln->xintp = -bx/(2*ax);
	
	// Fit peak in X to y = ax^2 + bx +c
	float by = ((yy1 - yy2)*(x2*x2 - x3*x3)/(x1*x1 - x2*x2) - (yy2-yy3))/(-(x2 - x3) + (x1 - x2)*(x2*x2 - x3*x3)/(x1*x1 - x2*x2));
	float ay = ((yy1 - yy2) - by*(x1 - x2))/(x1*x1 - x2*x2);
	//Find minima
	soln->yintp = -by/(2*ay);
	
	// Fit peak in X to y = ax^2 + bx +c
	float bz = ((yz1 - yz2)*(x2*x2 - x3*x3)/(x1*x1 - x2*x2) - (yz2-yz3))/(-(x2 - x3) + (x1 - x2)*(x2*x2 - x3*x3)/(x1*x1 - x2*x2));
	float az = ((yz1 - yz2) - bz*(x1 - x2))/(x1*x1 - x2*x2);
	//Find minima
	soln->zintp = -bz/(2*az);
**/
}

CudaPeakInfoFloat* calc_max_location_wrap_intp_cuda(const float* in, const int nx, const int ny, const int nz, const int maxdx, const int maxdy, const int maxdz)
{
	
	if(nz > 1)
	{
		const dim3 blockSize(1,1,1);
		const dim3 gridSize(1,1,1);
	
		CudaPeakInfoFloat * device_soln=0;
		cudaMalloc((void **)&device_soln, sizeof(CudaPeakInfoFloat));
		CudaPeakInfoFloat * host_soln = 0;
		host_soln = (CudaPeakInfoFloat*) malloc(sizeof(CudaPeakInfoFloat));

		calc_max_location_wrap_intp<<<blockSize,gridSize>>>(in,device_soln, nx, ny, nz, maxdx, maxdy, maxdz);
		cudaThreadSynchronize();
		cudaMemcpy(host_soln,device_soln,sizeof(CudaPeakInfoFloat),cudaMemcpyDeviceToHost);
		cudaFree(device_soln);

		return host_soln;
	} else {
		printf("Cannot to INTP on 2D yet\n");
		exit(1);
	}
}

__global__ void mult_kernel(float* data, const float scale, const int realtc)
{

	const uint index = threadIdx.x + (blockIdx.x + gridDim.x*blockIdx.y)*MAX_THREADS;
	
	if (index < realtc){
		data[index] *= scale;
	}
}

void emdata_processor_mult(float* data, const float& mult, const int nx, const int ny, const int nz) {

	int grid = int(ceil(sqrt(nx*ny*nz/MAX_THREADS)));
	const dim3 blockSize(MAX_THREADS,1, 1);
	const dim3 gridSize(grid,grid,1);
	mult_kernel<<<gridSize,blockSize>>>(data,mult,nx*ny*nz);
	
	cudaThreadSynchronize();
}

__global__ void mult_complex_eff_kernal(float* data, const float* src_data, const int nx, const int nxy, const int size)
{
	int idx = threadIdx.z*nxy + threadIdx.y*nx + threadIdx.x;

	data[idx] *= src_data[idx];
	data[size-idx-1] *= src_data[size-idx-1];
}

//this version, unlike the CPU version assumes the same image size, I don't know why you would ever need to mult two images of diffrent size....
void mult_complex_efficient_cuda(float* data, const float* src_data, const int nx, const int ny, const int nz, const int radius)
{

	int nxy = nx*ny;
	int size = nxy*nz;

	// now only for 2D and 3D, works on 3D till radius is <= 10 and for 2D radius is <= 33 (the radius is, in my experiance always small than 10)
	if(nz == 1)
	{
		const dim3 blockSize(radius,radius,1);
		const dim3 gridSize(1,1,1);
		mult_complex_eff_kernal<<<gridSize,blockSize>>>(data,src_data,nx,nxy,size);
	}else{
		const dim3 blockSize(radius,radius,radius);
		const dim3 gridSize(1,1,1);
		mult_complex_eff_kernal<<<gridSize,blockSize>>>(data,src_data,nx,nxy,size);
	}
	cudaThreadSynchronize();

}

__global__ void mcfauto_kernal(const float* data1, float* data2, const int totaltc)
{
	int idx = 2*(threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*MAX_THREADS);
	
	if(idx < totaltc){
		data2[idx] = sqrt(data1[idx] * data2[idx] + data1[idx + 1] * data2[idx + 1]);
		data2[idx + 1] = 0;
	}
}

void mcf_cuda(const float* data1, float* data2, const int nx, const int ny, const int nz)
{
	
	int grid = int(ceil(sqrt((nx/2)*ny*nz/MAX_THREADS)));

	const dim3 blockSize(MAX_THREADS,1, 1);
	const dim3 gridSize(grid,grid,1);
	mcfauto_kernal<<<gridSize,blockSize>>>(data1,data2,nx*ny*nz);

	cudaThreadSynchronize();
	
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

void swap_central_slices_180(float* data, const int nx, const int ny, const int nz)
{
	int xodd = (nx % 2) == 1;
	int yodd = (ny % 2) == 1;
	//int zodd = (nz % 2) == 1;

	//int nxy = nx * ny;

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
		cudaThreadSynchronize();
	}
	else // nx && ny && nz are greater than 1
	{
		throw;
	}
}

void swap_corners_180(float* data, const int nx, const int ny, const int nz)
{
	int xodd = (nx % 2) == 1;
	int yodd = (ny % 2) == 1;
	//int zodd = (nz % 2) == 1;

	//int nxy = nx * ny;

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
		cudaThreadSynchronize();
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

void emdata_phaseorigin_to_center(float* data, const int nx, const int ny, const int nz) 
{
	int xodd = (nx % 2) == 1;
	int yodd = (ny % 2) == 1;
	//int zodd = (nz % 2) == 1; // not yet supported

	//int nxy = nx * ny;
	if ( nz == 1 && ny > 1 ){
		// The order in which these operations occur literally undoes what the
		// PhaseToCornerProcessor did to the image.
		// First, the corners sections of the image are swapped appropriately
		swap_corners_180(data, nx, ny, nz);
		// Second, central pixel lines are swapped
		swap_central_slices_180(data, nx, ny, nz);

		// Third, appropriate sections of the image are cyclically shifted by one pixel
		if (xodd) {
			// Transfer the middle column to the far right
			// Shift all from the far right to (but not including the) middle one to the left
			middle_to_right<<<1,1>>>(data,nx,ny);
		}
		if (yodd) {
			// Tranfer the middle row to the top,
			// shifting all pixels from the top row down one, until  but not including the) middle
			middle_to_top<<<1,1>>>(data,nx,ny);
		}
		cudaThreadSynchronize();
	} else {
		throw;
	}
}

// This function just does the opposite of what emdata_phaseorigin_to_center does (it undoes its effects
void emdata_phaseorigin_to_corner(float* data, const int nx, const int ny, const int nz) 
{

	int xodd = (nx % 2) == 1;
	int yodd = (ny % 2) == 1;
	//int zodd = (nz % 2) == 1; // not yet supported

	if ( nz == 1 && ny > 1 ){
		// First ppropriate sections of the image are cyclically shifted by one pixel (this may not be right, I'll need to test, but it will work fine for even images...
		if (yodd) {
			// Tranfer the middle row to the top,
			// shifting all pixels from the top row down one, until  but not including the) middle
			middle_to_top<<<1,1>>>(data,nx,ny);
		}
		if (xodd) {
			// Transfer the middle column to the far right
			// Shift all from the far right to (but not including the) middle one to the left
			middle_to_right<<<1,1>>>(data,nx,ny);
		}
		// Second, central pixel lines are swapped
		swap_central_slices_180(data, nx, ny, nz);
		// Third, the corners sections of the image are swapped appropriately
		swap_corners_180(data, nx, ny, nz);

		cudaThreadSynchronize();
	} else {
		throw;
	}
}

__global__ void subtract_kernal(float* data, float f, const int totaltc)
{
	
	int idx = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*MAX_THREADS;
	
	if(idx < totaltc){
		data[idx] = data[idx] - f;
	}
}

void subtract_cuda(float* data, float f, const int nx, const int ny, const int nz)
{

	int grid = int(ceil(sqrt(nx*ny*nz/MAX_THREADS)));

	const dim3 blockSize(MAX_THREADS,1, 1);
	const dim3 gridSize(grid,grid,1);
	subtract_kernal<<<gridSize,blockSize>>>(data,f,nx*ny*nz);

	cudaThreadSynchronize();

}

__global__ void unwrap_kernel(float* dptr, const int num_threads, const int r1, const float p, const int nx, const int ny, const int nxp, const int dx,const int dy,const int weight_radial,const int offset) {

	const uint idx = threadIdx.x + blockIdx.x*num_threads+offset;

	const uint tx = idx % nxp;
	const uint ty = idx / nxp;

	float ang = tx * M_PI * p;
	float si = sinf(ang);
	float co = cosf(ang);

	float ypr1 = ty + r1;
	float xx = ypr1 * co + nx / 2 + dx;
	float yy = ypr1 * si + ny / 2 + dy;
	if ( weight_radial ) dptr[idx] = tex2D(texA2d,xx+0.5,yy+0.5)*ypr1;
	else dptr[idx] = tex2D(texA2d,xx+0.5,yy+0.5);
}


void emdata_unwrap(float* data, int r1, int r2, int xs, int num_pi, int dx, int dy, int weight_radial, int nx, int ny) {

	int n = xs*(r2-r1);
	int num_calcs = n;

	int grid_y = num_calcs/(MAX_THREADS);
	int res_y = num_calcs - grid_y*MAX_THREADS;

	//printf("Grid %d, res %d, n %d, p %f \n",grid_y,res_y,n, p/xs);

	if ( grid_y > 0 ) {
		const dim3 blockSize(MAX_THREADS,1, 1);
		const dim3 gridSize(grid_y,1,1);
		unwrap_kernel<<<gridSize,blockSize>>>(data,MAX_THREADS,r1,(float) num_pi/ (float)xs, nx,ny,xs,dx,dy,weight_radial,0);
	}

	if ( res_y > 0 ) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		unwrap_kernel<<<gridSize,blockSize>>>(data,MAX_THREADS,r1, (float) num_pi/ (float)xs, nx,ny,xs,dx,dy,weight_radial,grid_y*MAX_THREADS);
	}
	cudaThreadSynchronize();

}

__global__ void column_sum(const float* data, float* sum, int nx, int ny, int num_threads, int offset ) {
	
	float s = 0.0;
	const uint idx = threadIdx.x + blockIdx.x*num_threads+offset; 
	for(int i =0; i < ny; i++) {
		s += data[idx + i*nx];
	}
	sum[idx] = s; 
}

float* emdata_column_sum(const float* data, const int nx, const int ny) {

        int max_threads = MAX_THREADS;
	if (max_threads > nx) max_threads = nx;
	
	int num_calcs = nx;

	float* sum = 0;
	cudaMalloc((void **)&sum, num_calcs*sizeof(float));
	
	int grid_y = num_calcs/(max_threads);
	int res_y = (num_calcs - (grid_y*max_threads));
	
	if ( grid_y > 0 ) {
		const dim3 blockSize(max_threads,1, 1);
		const dim3 gridSize(grid_y,1,1);
		column_sum<<<gridSize,blockSize>>>(data,sum,nx,ny,max_threads,0);
	}
	
	if ( res_y > 0 ) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		column_sum<<<gridSize,blockSize>>>(data,sum,nx,ny,max_threads,grid_y*max_threads);
	}
	cudaThreadSynchronize();

	return sum;
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

void emdata_rotate_180(float* data, const int nx, const int ny) {

	// Only works for 2D images

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
		offset =  nx;
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
		rotate_180<<<gridSize,blockSize>>>(data+offset,nx,nxy,0,0);
	}

	if ( res_y > 0 ) {
		const dim3 blockSize(MAX_THREADS,1, 1);
		const dim3 gridSize(1,1,1);
		rotate_180<<<gridSize,blockSize>>>(data+offset,nx,nxy,grid_y*MAX_THREADS,res_y);
	}

	cudaThreadSynchronize();
}
