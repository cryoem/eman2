
typedef unsigned int uint;
__device__ void add_complex_at_df(float* vol, float* tmp_data, const int x, const int y, const int z,const int nx, const int ny, const int nz, const float dtr, const float dti, const float gg)
{

	// we could use shared memory here, might speed thingup as vol tmp_data are acessed at least twice per thread on avg
	uint idx;
	if (x<0) {
		idx=-x*2+(y<=0?-y:ny-y)*nx+(z<=0?-z:nz-z)*nx*ny;
		vol[idx]+=dtr;
		vol[idx+1]-=dti;
		tmp_data[idx/2]+=gg;
		//atomicAdd(vol + idx, dtr);
		//atomicAdd(vol + idx + 1, -dti);
		//atomicAdd(tmp_data + idx/2, gg);
		return;
	}

	idx=x*2+(y<0?ny+y:y)*nx+(z<0?nz+z:z)*nx*ny;
	vol[idx]+=dtr;
	vol[idx+1]+=dti;
	tmp_data[idx/2]+=gg;
	//atomicAdd(vol + idx, dtr);
	//atomicAdd(vol + idx + 1, dti);
	//atomicAdd(tmp_data + idx/2, gg);
}

__const__ float hdecay = 1.0f/((float) (4.0/(M_PI*M_PI)));
__device__ void insert_pixel_df(float* vol, float* tmp_data, const int nx, const int ny, const int nz, const float xx, const float yy, const float zz, const float dtr, const float dti, const float weight)
{

	int x0 = (int) floor(xx);
	int y0 = (int) floor(yy);
	int z0 = (int) floor(zz);
	

	if (x0<-nx/2-2 || y0<-ny/2-1 || z0<-nz/2-1 || x0>nx/2-1 || y0>ny/2 || z0>nz/2 ) return;

	// no error checking on add_complex_fast, so we need to be careful here
	int x1=x0+1;
	int y1=y0+1;
	int z1=z0+1;
	if (x0<-(nx/2-1)) x0=-(nx/2-1);
	if (x1>(nx/2-1)) x1=(nx/2-1);
	if (y0<-ny/2) y0=-ny/2;
	if (y1>ny/2) y1=ny/2;
	if (z0<-nz/2) z0=-nz/2;
	if (z1>nz/2) z1=nz/2;

	//unroll loops
	float r = (float(x0) - xx)*(float(x0) - xx) + (float(y0) - yy)*(float(y0) - yy) + (float(z0) - zz)*(float(z0) - zz);
	float gg = expf(-r*hdecay)*weight;
	add_complex_at_df(vol, tmp_data, x0, y0, z0, nx, ny, nz, dtr*gg, dti*gg, gg);

	r = (float(x1) - xx)*(float(x1) - xx) + (float(y0) - yy)*(float(y0) - yy) + (float(z0) - zz)*(float(z0) - zz);
	gg = expf(-r*hdecay)*weight;
	add_complex_at_df(vol, tmp_data, x1, y0, z0, nx, ny, nz, dtr*gg, dti*gg, gg);

	r = (float(x0) - xx)*(float(x0) - xx) + (float(y1) - yy)*(float(y1) - yy) + (float(z0) - zz)*(float(z0) - zz);
	gg = expf(-r*hdecay)*weight;
	add_complex_at_df(vol, tmp_data, x0, y1, z0, nx, ny, nz, dtr*gg, dti*gg, gg);

	r = (float(x1) - xx)*(float(x1) - xx) + (float(y1) - yy)*(float(y1) - yy) + (float(z0) - zz)*(float(z0) - zz);
	gg = expf(-r*hdecay)*weight;
	add_complex_at_df(vol, tmp_data, x1, y1, z0, nx, ny, nz, dtr*gg, dti*gg, gg);

	r = (float(x0) - xx)*(float(x0) - xx) + (float(y0) - yy)*(float(y0) - yy) + (float(z1) - zz)*(float(z1) - zz);
	gg = expf(-r*hdecay)*weight;
	add_complex_at_df(vol, tmp_data, x0, y0, z1, nx, ny, nz, dtr*gg, dti*gg, gg);

	r = (float(x1) - xx)*(float(x1) - xx) + (float(y0) - yy)*(float(y0) - yy) + (float(z1) - zz)*(float(z1) - zz);
	gg = expf(-r*hdecay)*weight;
	add_complex_at_df(vol, tmp_data, x1, y0, z1, nx, ny, nz, dtr*gg, dti*gg, gg);

	r = (float(x0) - xx)*(float(x0) - xx) + (float(y1) - yy)*(float(y1) - yy) + (float(z1) - zz)*(float(z1) - zz);
	gg = expf(-r*hdecay)*weight;
	add_complex_at_df(vol, tmp_data, x0, y1, z1, nx, ny, nz, dtr*gg, dti*gg, gg);

	r = (float(x1) - xx)*(float(x1) - xx) + (float(y1) - yy)*(float(y1) - yy) + (float(z1) - zz)*(float(z1) - zz);
	gg = expf(-r*hdecay)*weight;
	add_complex_at_df(vol, tmp_data, x1, y1, z1, nx, ny, nz, dtr*gg, dti*gg, gg);

}

__global__ void insert_slice_kernal(const float* const slice_data, float* vol, float* tmp_data, const int inx, const int iny, const int nx, const int ny, const int nz, float3 mxx, float3 myy, float3 mzz, const float weight)
{

	//Here is where we do each pixel
	int x=threadIdx.x; 
	int y=blockIdx.x - iny/2; //need b/c......

	float rx = (float) x/(inx-2.0f);	// coords relative to Nyquist=.5
	float ry = (float) y/iny;

	float xx = rx*mxx.x + ry*mxx.y;
	float yy = rx*myy.x + ry*myy.y;
	float zz = rx*mzz.x + ry*mzz.y;

	xx=xx*(nx-2);
	yy=yy*ny;
	zz=zz*nz;

	//get_complex_at function
	float dtr, dti;
	if (fabsf(x)>=inx/2 || fabsf(y)>iny/2) {
		dtr = 0.0f; dti = 0.0f;
	} else if (x>=0 && y>=0) {
		dtr = slice_data[ x*2+y*inx]; dti = slice_data[x*2+y*inx+1];
	} else if (x>0 && y<0) {
		dtr = slice_data[ x*2+(iny+y)*inx]; dti = slice_data[x*2+(iny+y)*inx+1];
	} else if (x<0 && y>0) {
		dtr = slice_data[-x*2+(iny-y)*inx]; dti = -slice_data[-x*2+(iny-y)*inx+1];
	} else {
		dtr = slice_data[-x*2-y*inx]; dti = -slice_data[-x*2+-y*inx+1];
	}

	//Here is where we actuall insert the pixel, for each inserter, we'll need a new kernal
	insert_pixel_df(vol, tmp_data, nx, ny, nz, xx, yy, zz, dtr, dti, weight);
}

__global__ void stevens_insert_slice_kernal(const float* const slice_data, float* vol, float* tmp_data, const int inx, const int iny, const int nx, const int ny, const int nz, float3 mxx, float3 myy, float3 mzz, const float weight)
{

	//Here is where we do each pixel
	int x=threadIdx.x; 
	int y=blockIdx.x - iny/4; //need b/c......

	x*=2;
	y*=2;
	for(int i = 0; i < 2; i++) {
		x+=i;
		for(int j = 0; j < 2; j++) {
			y+=j;

			float rx = (float) x/(inx-2.0f);	// coords relative to Nyquist=.5
			float ry = (float) y/iny;

			float xx = rx*mxx.x + ry*mxx.y;
			float yy = rx*myy.x + ry*myy.y;
			float zz = rx*mzz.x + ry*mzz.y;

			xx=xx*(nx-2);
			yy=yy*ny;
			zz=zz*nz;

			//get_complex_at function
			float dtr, dti;
			if (fabsf(x)>=inx/2 || fabsf(y)>iny/2) {
				dtr = 0.0f; dti = 0.0f;
			} else if (x>=0 && y>=0) {
				dtr = slice_data[ x*2+y*inx]; dti = slice_data[x*2+y*inx+1];
			} else if (x>0 && y<0) {
				dtr = slice_data[ x*2+(iny+y)*inx]; dti = slice_data[x*2+(iny+y)*inx+1];
			} else if (x<0 && y>0) {
				dtr = slice_data[-x*2+(iny-y)*inx]; dti = -slice_data[-x*2+(iny-y)*inx+1];
			} else {
				dtr = slice_data[-x*2-y*inx]; dti = -slice_data[-x*2+-y*inx+1];
			}

			//Here is where we actuall insert the pixel, for each inserter, we'll need a new kernal
			insert_pixel_df(vol, tmp_data, nx, ny, nz, xx, yy, zz, dtr, dti, weight);
		}
		y-=1;
	}
}

void insert_slice_cuda(const float* const matrix, const float* const slice_data, float* vol, float* tmp_data, const int inx, const int iny, const int nx, const int ny, const int nz, const float weight)
{

	const dim3 gridSize(iny/2,1, 1);
	const dim3 blockSize((inx/2 + 1)/2,1,1);

	float3 mxx,myy,mzz;
	
	mxx.x=matrix[0];
	mxx.y=matrix[4];
	mxx.z=matrix[8];
	myy.x=matrix[1];
	myy.y=matrix[5];
	myy.z=matrix[9];
	mzz.x=matrix[2];
	mzz.y=matrix[6];
	mzz.z=matrix[10];

	stevens_insert_slice_kernal<<<gridSize,blockSize>>>(slice_data, vol, tmp_data, inx ,iny, nx, ny, nz, mxx, myy, mzz, weight);
	cudaThreadSynchronize();
}

__device__ int get_complex_index_kernal(const int x, const int y, const int z, const int nx, const int ny, const int nz)
{
	if (x<0) {
		return x*-2+(y<=0?-y:ny-y)*nx+(z<=0?-z:nz-z)*nx*ny;
	}
	return x*2+(y<0?ny+y:y)*nx+(z<0?nz+z:z)*nx*ny;
}

__const__ float hagg = (float) (4.0/(M_PI*M_PI));
__device__ bool pixel_at(float* vol, float* tmp_data, const int nx, const int ny, const int nz, const float xx, const float yy, const float zz, float3 &dt)
{
	int x0 = (int) floor(xx);
	int y0 = (int) floor(yy);
	int z0 = (int) floor(zz);

	//intialize vars
	float normsum=0,normsum2=0;
	dt.x = 0.0f; dt.y = 0.0f; dt.z = 0.0f;

	if (x0<-nx/2-2 || y0<-ny/2-1 || z0<-nz/2-1 || x0>nx/2-1 || y0>ny/2 || z0>nz/2 ) return false;

	// no error checking on add_complex_fast, so we need to be careful here
	int x1=x0+1;
	int y1=y0+1;
	int z1=z0+1;
	if (x0<-(nx/2-1)) x0=-(nx/2-1);
	if (x1>(nx/2-1)) x1=(nx/2-1);
	if (y0<-ny/2) y0=-ny/2;
	if (y1>ny/2) y1=ny/2;
	if (z0<-nz/2) z0=-nz/2;
	if (z1>nz/2) z1=nz/2;

	
	// Unroll loops
	float r = powf((float(x0) - xx),2) + powf((float(y0) - yy),2) + powf((float(z0) - zz),2);
	float gg = expf(-r/hagg);
	int idx = get_complex_index_kernal(x0, y0, z0, nx, ny, nz);
	dt.x +=gg*vol[idx];
	dt.y+=(x0<0?-1.0f:1.0f)*gg*vol[idx+1];
	r = tmp_data[idx/2]; // I am recycling the r register to save registers
	dt.z+=r*gg;
	normsum2+=gg;
	normsum+=gg*r;

	r = powf((float(x1) - xx),2) + powf((float(y0) - yy),2) + powf((float(z0) - zz),2);
	gg = expf(-r/hagg);
	idx = get_complex_index_kernal(x1, y0, z0, nx, ny, nz);
	dt.x +=gg*vol[idx];
	dt.y+=(x1<0?-1.0f:1.0f)*gg*vol[idx+1];
	r = tmp_data[idx/2]; // I am recycling the r register to save registers
	dt.z+=r*gg;
	normsum2+=gg;
	normsum+=gg*r;

	r = powf((float(x0) - xx),2) + powf((float(y1) - yy),2) + powf((float(z0) - zz),2);
	gg = expf(-r/hagg);
	idx = get_complex_index_kernal(x0, y1, z0, nx, ny, nz);
	dt.x +=gg*vol[idx];
	dt.y+=(x0<0?-1.0f:1.0f)*gg*vol[idx+1];
	r = tmp_data[idx/2]; // I am recycling the r register to save registers
	dt.z+=r*gg;
	normsum2+=gg;
	normsum+=gg*r;

	r = powf((float(x1) - xx),2) + powf((float(y1) - yy),2) + powf((float(z0) - zz),2);
	gg = expf(-r/hagg);
	idx = get_complex_index_kernal(x1, y1, z0, nx, ny, nz);
	dt.x +=gg*vol[idx];
	dt.y+=(x1<0?-1.0f:1.0f)*gg*vol[idx+1];
	r = tmp_data[idx/2]; // I am recycling the r register to save registers
	dt.z+=r*gg;
	normsum2+=gg;
	normsum+=gg*r;

	r = powf((float(x0) - xx),2) + powf((float(y0) - yy),2) + powf((float(z1) - zz),2);
	gg = expf(-r/hagg);
	idx = get_complex_index_kernal(x0, y0, z1, nx, ny, nz);
	dt.x +=gg*vol[idx];
	dt.y+=(x0<0?-1.0f:1.0f)*gg*vol[idx+1];
	r = tmp_data[idx/2]; // I am recycling the r register to save registers
	dt.z+=r*gg;
	normsum2+=gg;
	normsum+=gg*r;

	r = powf((float(x1) - xx),2) + powf((float(y0) - yy),2) + powf((float(z1) - zz),2);
	gg = expf(-r/hagg);
	idx = get_complex_index_kernal(x1, y0, z1, nx, ny, nz);
	dt.x +=gg*vol[idx];
	dt.y+=(x1<0?-1.0f:1.0f)*gg*vol[idx+1];
	r = tmp_data[idx/2]; // I am recycling the r register to save registers
	dt.z+=r*gg;
	normsum2+=gg;
	normsum+=gg*r;

	r = powf((float(x0) - xx),2) + powf((float(y1) - yy),2) + powf((float(z1) - zz),2);
	gg = expf(-r/hagg);
	idx = get_complex_index_kernal(x0, y1, z1, nx, ny, nz);
	dt.x +=gg*vol[idx];
	dt.y+=(x0<0?-1.0f:1.0f)*gg*vol[idx+1];
	r = tmp_data[idx/2]; // I am recycling the r register to save registers
	dt.z+=r*gg;
	normsum2+=gg;
	normsum+=gg*r;

	r = powf((float(x1) - xx),2) + powf((float(y1) - yy),2) + powf((float(z1) - zz),2);
	gg = expf(-r/hagg);
	idx = get_complex_index_kernal(x1, y1, z1, nx, ny, nz);
	dt.x +=gg*vol[idx];
	dt.y+=(x1<0?-1.0f:1.0f)*gg*vol[idx+1];
	r = tmp_data[idx/2]; // I am recycling the r register to save registers
	dt.z+=r*gg;
	normsum2+=gg;
	normsum+=gg*r;

	if (normsum==0) return false;

	dt.x/=normsum;
	dt.y/=normsum;
	dt.z/=normsum2;

	return true;
}

__global__ void  determine_slice_kernal(const float* const slice_data, float* vol, float* tmp_data, float* results, const int inx, const int iny, const int nx, const int ny, const int nz, float3 mxx, float3 myy, float3 mzz, const float weight)
{

	//Here is where we do each pixel
	int x=threadIdx.x; //need b/c......
	int y=blockIdx.x - iny/2;
	int idx;
	int size = blockDim.x*gridDim.x;

	float rx = (float) x/(inx-2.0f);	// coords relative to Nyquist=.5
	float ry = (float) y/iny;

	float xx = rx*mxx.x + ry*mxx.y;
	float yy = rx*myy.x + ry*myy.y;
	float zz = rx*mzz.x + ry*mzz.y;

	if (fabsf(xx)>0.5 || fabsf(yy)>=0.5 || fabsf(zz)>=0.5) {
		idx = blockIdx.x*blockDim.x + threadIdx.x;
		results[idx] = 0.0;
		results[size + idx] = 0.0;
		results[2*size + idx] = 0.0;
		results[3*size + idx] = 0.0;
		return;
	}

	xx=xx*(nx-2);
	yy=yy*ny;
	zz=zz*nz;

	float2 dt2; float3 dt;
	idx = (int)(x * 2 + inx*(y<0?iny+y:y));
	dt2.x = slice_data[idx];
	dt2.y = slice_data[idx+1];

	//Here is where we actuall insert the pixel, for each inserter, we'll need a new kernal
	if (!pixel_at(vol, tmp_data, nx, ny, nz, xx, yy, zz, dt)){
		idx = blockIdx.x*blockDim.x + threadIdx.x;
		results[idx] = 0.0;
		results[size + idx] = 0.0;
		results[2*size + idx] = 0.0;
		results[3*size + idx] = 0.0;
		return;
	}

	idx = blockIdx.x*blockDim.x + threadIdx.x;
	results[idx]=(dt.x*dt2.x+dt.y*dt2.y)*dt.z;
	results[size + idx]=dt.z;
	results[2*size + idx]=(dt.x*dt.x+dt.y*dt.y)*dt.z;
	results[3*size + idx]=(dt2.x*dt2.x+dt2.y*dt2.y)*dt.z;

	
}

__global__ void dsa_sumup_kernal(float* results, const int inx, const int iny)
{
	int tidx = threadIdx.x;
	int bd = blockDim.x;
	int size = iny*(inx/2 + 1);

	float dot = 0.0f; float vweight = 0.0f; float power = 0.0f; float power2 = 0.0f;
	for (int i = 0; i < (inx/2 + 1); i++) {
		int idx = i*bd + tidx;
		dot += results[idx];
		vweight += results[size + idx];
		power += results[2*size + idx];
		power2 += results[3*size + idx];
	}

	results[tidx] = dot;
	results[size + tidx] = vweight;
	results[2*size + tidx] = power;
	results[3*size + tidx] = power2;

}

float4 determine_slice_agreement_cuda(const float* const matrix, const float* const slice_data, float* vol, float* tmp_data, const int inx, const int iny, const int nx, const int ny, const int nz, const float weight)
{
	const dim3 gridSize(iny,1, 1);
	const dim3 blockSize((inx/2 + 1),1,1);

	float* results=0;
	cudaMalloc((void **)&results, 4*(inx/2+1)*iny*sizeof(float));

	float3 mxx,myy,mzz;
	
	mxx.x=matrix[0];
	mxx.y=matrix[4];
	mxx.z=matrix[8];
	myy.x=matrix[1];
	myy.y=matrix[5];
	myy.z=matrix[9];
	mzz.x=matrix[2];
	mzz.y=matrix[6];
	mzz.z=matrix[10];

	//Here is where we actuall work for each pixel
	determine_slice_kernal<<<gridSize,blockSize>>>(slice_data, vol, tmp_data, results, inx ,iny, nx, ny, nz, mxx, myy, mzz, weight);
	cudaThreadSynchronize();

	//now do the reductions
	float4 stats;
	stats.x = 0.0f; stats.y = 0.0f; stats.z = 0.0f; stats.w = 0.0f;
	dsa_sumup_kernal<<<1,iny>>>(results, inx, iny); // a bit slight of hand here....

	//copy back from device
	float* host_results = (float*) malloc(4*(inx/2+1)*iny*sizeof(float));
	cudaMemcpy(host_results,results,4*(inx/2+1)*iny*sizeof(float),cudaMemcpyDeviceToHost);
	cudaFree(results);

	int size = iny*(inx/2 + 1);
	for (int i = 0; i < iny; i++) {
		//printf("step %d dot %f vweight %f power %f power2 %f\n",i,host_results[i],host_results[size + i],host_results[2*size + i],host_results[3*size + i]); 
		stats.x += host_results[i];
		stats.y += host_results[size + i];
		stats.z += host_results[2*size + i];
		stats.w += host_results[3*size + i];
	}
	free(host_results);

	return stats;
}
