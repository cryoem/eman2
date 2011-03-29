

typedef unsigned int uint;

#ifdef WIN32
	#define M_PI 3.14159265358979323846f
#endif	//WIN32

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

	out[idx]=tex3D(texA, tx,ty,tz);
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

	if (grid_y > 0) {
		const dim3 blockSize(MAX_THREADS,1, 1);
		const dim3 gridSize(grid_y,1,1);
		if ( nz > 1 ) {
			transform_kernel_3D<<<gridSize,blockSize>>>(data,nx,ny,nz,MAX_THREADS,mxx,mxy,mxz,trans,0);
		}
		else if ( ny > 1 ) {
			transform_kernel_2D<<<gridSize,blockSize>>>(data,nx,ny,MAX_THREADS,mxx,mxy,mxz,trans,0);
		} else throw;
	}

	if (res_y) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		if ( nz > 1 ) {
			transform_kernel_3D<<<gridSize,blockSize>>>(data,nx,ny,nz,MAX_THREADS,mxx,mxy,mxz,trans,grid_y*MAX_THREADS);
		}
		else if ( ny > 1 ) {
			transform_kernel_2D<<<gridSize,blockSize>>>(data,nx,ny,MAX_THREADS,mxx,mxy,mxz,trans,grid_y*MAX_THREADS);
		} else throw;
	}
	cudaThreadSynchronize();

	return data;
}

__global__ void complexmult_conj_kernal(float *afft, const float *bfft, int num_threads, int offset)
{

	const uint ridx = 2*(threadIdx.x + blockIdx.x*num_threads + offset);
	const uint iidx = ridx + 1;
	//maybe use float2 to improve coalessing....

	float afftr = afft[ridx];
	float affti = afft[iidx];
	float bfftr = bfft[ridx];
	float bffti = bfft[iidx];

	afft[ridx] = afftr*bfftr + affti*bffti;  //real portion
	afft[iidx] = affti*bfftr - afftr*bffti; //imaginary portion

}

void calc_ccf_cuda(float* afft, const float* bfft, const int nx, const int ny, const int nz)
{

	int num_calcs = ny*(nx/2)*nz;
	int grid_y = num_calcs/(MAX_THREADS);
	int res_y = num_calcs - grid_y*MAX_THREADS;

	if (grid_y > 0) {
		const dim3 blockSize(MAX_THREADS,1, 1);
		const dim3 gridSize(grid_y,1,1);
		complexmult_conj_kernal<<<gridSize,blockSize>>>(afft, bfft, MAX_THREADS, 0);
	}
	if (res_y) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		complexmult_conj_kernal<<<gridSize,blockSize>>>(afft, bfft, MAX_THREADS, grid_y*MAX_THREADS);
	}
	cudaThreadSynchronize();

}

__global__ void complexmult_kernal(float *afft, const float *bfft, int num_threads, int offset)
{

	const uint ridx = 2*(threadIdx.x + blockIdx.x*num_threads + offset);
	const uint iidx = ridx + 1;
	//maybe use float2 to improve coalessing....

	float afftr = afft[ridx];
	float affti = afft[iidx];
	float bfftr = bfft[ridx];
	float bffti = bfft[iidx];

	afft[ridx] = afftr*bfftr - affti*bffti;  //real portion
	afft[iidx] = affti*bfftr + afftr*bffti; //imaginary portion

}

void calc_conv_cuda(float* afft, const float* bfft, const int nx, const int ny, const int nz)
{

	int num_calcs = ny*(nx/2)*nz;
	int grid_y = num_calcs/(MAX_THREADS);
	int res_y = num_calcs - grid_y*MAX_THREADS;

	if (grid_y > 0) {
		const dim3 blockSize(MAX_THREADS,1, 1);
		const dim3 gridSize(grid_y,1,1);
		complexmult_kernal<<<gridSize,blockSize>>>(afft, bfft, MAX_THREADS, 0);
	}
	if (res_y) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		complexmult_kernal<<<gridSize,blockSize>>>(afft, bfft, MAX_THREADS, grid_y*MAX_THREADS);
	}
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

__global__ void mult_kernel(float* data, const float scale, const int num_threads)
{

	const uint x=threadIdx.x;
	const uint y=blockIdx.x;

	data[x+y*num_threads] *= scale;
}

void emdata_processor_mult(float* data, const float& mult, const int nx, const int ny, const int nz) {

	int num_calcs = nx*ny*nz;

	int grid_y = num_calcs/MAX_THREADS;
	int res_y = num_calcs - (grid_y*MAX_THREADS);

	if ( grid_y > 0 ) {
		const dim3 blockSize(MAX_THREADS,1, 1);
		const dim3 gridSize(grid_y,1,1);
		mult_kernel<<<gridSize,blockSize>>>(data,mult,MAX_THREADS);
	}

	if ( res_y > 0 ) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		mult_kernel<<<gridSize,blockSize>>>(data+grid_y*MAX_THREADS,mult,0);
	}

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

__global__ void mcfauto_kernal(const float* data1, float* data2, const int num_threads, const int offset)
{
	int idx = 2*(threadIdx.x + blockIdx.x*num_threads + offset);
	
	data2[idx] = sqrt(data1[idx] * data2[idx] + data1[idx + 1] * data2[idx + 1]);
	data2[idx + 1] = 0;
}

void mcf_cuda(const float* data1, float* data2, const int nx, const int ny, const int nz)
{
	
	int num_calcs = nx*ny*nz/2;

	int grid_y = num_calcs/MAX_THREADS;
	int res_y = num_calcs - (grid_y*MAX_THREADS);

	if ( grid_y > 0 ) {
		const dim3 blockSize(MAX_THREADS,1, 1);
		const dim3 gridSize(grid_y,1,1);
		mcfauto_kernal<<<gridSize,blockSize>>>(data1,data2,MAX_THREADS,0);
	}
	if ( res_y > 0 ) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		mcfauto_kernal<<<gridSize,blockSize>>>(data1,data2,MAX_THREADS,grid_y*MAX_THREADS);
	}
	
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

__global__ void subtract_kernal(float* data, float f, const int num_threads, const int offset)
{

	int idx = threadIdx.x + blockIdx.x*num_threads + offset;
	
	data[idx] = data[idx] - f;
}

void subtract_cuda(float* data, float f, const int nx, const int ny, const int nz)
{

	int num_calcs = nx*ny*nz;

	int grid_y = num_calcs/MAX_THREADS;
	int res_y = num_calcs - (grid_y*MAX_THREADS);

	if ( grid_y > 0 ) {
		const dim3 blockSize(MAX_THREADS,1, 1);
		const dim3 gridSize(grid_y,1,1);
		subtract_kernal<<<gridSize,blockSize>>>(data,f,MAX_THREADS,0);
	}
	if ( res_y > 0 ) {
		const dim3 blockSize(res_y,1, 1);
		const dim3 gridSize(1,1,1);
		subtract_kernal<<<gridSize,blockSize>>>(data,f,MAX_THREADS,grid_y*MAX_THREADS);
	}

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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//NOTHING BUT RUBBISH BELOW HERE

/*
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

void emdata_phaseorigin_to_center_fourier(float* data, const int nx, const int ny, const int nz) {

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

*/
