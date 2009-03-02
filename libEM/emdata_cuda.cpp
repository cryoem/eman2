/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
 * Copyright (c) 2000-2006 Baylor College of Medicine
 *
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * */



#ifdef EMAN2_USING_CUDA

#include "emdata.h"
#include "exception.h"
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "cuda/cuda_util.h"
#include "cuda/cuda_processor.h"

using namespace EMAN;
// Static init
EMData::CudaDeviceEMDataCache EMData::cuda_rw_cache(100);

float* EMData::get_cuda_data() const {
	if (cuda_pointer_handle==-1 || EMDATA_GPU_NEEDS_UPDATE & flags) {
		if (cuda_pointer_handle != -1 && !(EMDATA_GPU_RO_NEEDS_UPDATE & flags ) ) {
			throw;
			cuda_rw_cache.copy_ro_to_rw_data(cuda_pointer_handle);
		}
		cuda_pointer_handle = cuda_rw_cache.cache_data(this,rdata,nx,ny,nz);
		flags &= ~EMDATA_GPU_NEEDS_UPDATE;
	}
	return cuda_rw_cache.get_rw_data(cuda_pointer_handle);
}

void EMData::bind_cuda_array() {
	if (cuda_array_handle==-1|| EMDATA_GPU_RO_NEEDS_UPDATE & flags) {
		if (cuda_array_handle!=-1)  {
			delete_cuda_array(cuda_array_handle);
			
		}
		cuda_array_handle = get_cuda_array_handle(get_cuda_data(),nx,ny,nz,this);
		
		flags &= ~EMDATA_GPU_RO_NEEDS_UPDATE;
	}
	bind_cuda_texture(cuda_array_handle);
}

void EMData::cuda_cache_lost_imminently() const {
	get_data(); // This causes cuda memory to be copied to cpu memory
	cuda_pointer_handle = -1;
}

void EMData::gpu_update() {
	flags |= EMDATA_NEEDUPD | EMDATA_CPU_NEEDS_UPDATE | EMDATA_GPU_RO_NEEDS_UPDATE;
}

void EMData::mult_cuda(const float& val) {
	Dict d("scale",(float)val);
	process_inplace("cuda.math.mult",d);
}

EMData* EMData::calc_ccf_cuda( EMData*  image ) {
	
	EMData* tmp;
	if (is_complex()) {
		tmp = new EMData(*this);
	} else {
		tmp = do_fft_cuda();
	}
	
	Dict d;
	EMData* with = 0;
	if (!image->is_complex()) {
		with = image->do_fft_cuda();
		d["with"] = (EMData*) with;
	} else {
		d["with"] = (EMData*)image;
	}
	
	tmp->process_inplace("cuda.correlate",d);
	
	if (with != 0) {
		delete with;
		with = 0;
	}

	EMData* soln = tmp->do_ift_cuda();
	soln->gpu_update();
	delete tmp;
	tmp = 0;
	
	return soln;
}

void EMData::free_cuda_array() const {
	if (cuda_array_handle != -1){
		delete_cuda_array(cuda_array_handle);
		cuda_array_handle = -1;
	}
}

void EMData::free_cuda_memory() const {
	if (cuda_pointer_handle!=-1) {
		cuda_rw_cache.clear_item(cuda_pointer_handle);
		cuda_pointer_handle = -1;
	}
}



EMData::CudaDeviceEMDataCache::CudaDeviceEMDataCache(const unsigned int size) : cache_size(size),current_insert_idx(0)
{
	device_init();
	rw_cache = new float *[cache_size];
	caller_cache = new const EMData*[cache_size];
	ro_cache = new cudaArray *[cache_size];
	
	for(unsigned int i = 0; i < cache_size; ++ i ) {
		rw_cache[i] = 0;
		caller_cache[i] = 0;
		ro_cache[i] = 0;
	}
}

EMData::CudaDeviceEMDataCache::~CudaDeviceEMDataCache()
{
	for (unsigned int i = 0; i < cache_size; i++) {
		clear_item(i);
	}

	if( rw_cache )
	{
		delete[]rw_cache;
		rw_cache = 0;
	}
	
	// No deletion responsibility for the caller_cache
	if( caller_cache )
	{
		delete[]caller_cache;
		caller_cache = 0;
	}
}

unsigned int EMData::CudaDeviceEMDataCache::cache_data(const EMData* const emdata, const float* const data,const int nx, const int ny, const int nz)
{
	const EMData* previous = caller_cache[current_insert_idx];
	if (previous != 0) {
		previous->cuda_cache_lost_imminently();
		clear_item(current_insert_idx);
	}
	
	float* cuda_rw_data = alloc_rw_data(nx,ny,nz);
	
	if (data != 0 ) { // If rdata is zero it means we're working exclusively on the GPU
		size_t num_bytes = nx*ny*nz*sizeof(float);
		cudaMemcpy(cuda_rw_data,data,num_bytes,cudaMemcpyHostToDevice);
	}
	
	rw_cache[current_insert_idx] = cuda_rw_data;
	caller_cache[current_insert_idx] = emdata;
	ro_cache[current_insert_idx] = 0;
	
	unsigned int ret = current_insert_idx;
	current_insert_idx += 1;
	current_insert_idx %= cache_size; // Potentially inefficient to do this every time, the alternative is an if statement. Which is faster?
	return ret;
}

float* EMData::CudaDeviceEMDataCache::alloc_rw_data(const int nx, const int ny, const int nz) {
	float* cuda_rw_data;
	size_t num_bytes = nx*ny*nz*sizeof(float);

	cudaMalloc((void**)&cuda_rw_data,num_bytes);
	
	if (cuda_rw_data == 0) {
		stringstream ss;
		string gigs;
		ss << (float) num_bytes/1000000000.0;
		ss >> gigs;
		string message = "GPU - Can't allocate " + gigs + " GB - not enough memory.";
		throw BadAllocException(message);
		//throw BadAllocException("Cuda malloc failed, memory may be exhaused");
	}
	return cuda_rw_data;
	
}

float* EMData::CudaDeviceEMDataCache::get_rw_data(const unsigned int idx) {
	if (rw_cache[idx] != 0) return rw_cache[idx]; // This is a common case

	// If there is no current copy of the rw data then the ro should definitely exist - that is the only way the EMData object should have the slot index (idx);
	if (ro_cache[idx] == 0) throw UnexpectedBehaviorException("The CUDA read write data was requested, it was 0, so there was an attempt to copy the CUDA RO data, but this was 0.");
	
	const EMData* d = caller_cache[idx];
	int nx = d->get_xsize();
	int ny = d->get_ysize();
	int nz = d->get_zsize();
	
	rw_cache[idx] = alloc_rw_data(nx,ny,nz);
	
	if (nz == 1) {
		cudaMemcpyFromArray(rw_cache[idx], ro_cache[idx], nx,ny, nx*ny*sizeof(float), cudaMemcpyDeviceToDevice);
	} else if (nz > 1) {
		throw;
		/*
		cudaExtent ca_extent;
		ca_extent.width  = nx;
		ca_extent.height = ny;
		ca_extent.depth  = nz;
		
		cudaMemcpy3DParms copyParams = {0};
		copyParams.srcArray   = ro_cache[idx];
		copyParams.dstPtr = rw_cache[idx];
		copyParams.extent   = ca_extent;
		copyParams.kind     = cudaMemcpyDeviceToDevice;
		
		cudaMemcpy3D(&copyParams);
		*/
	}
}

unsigned int EMData::CudaDeviceEMDataCache::cache_ro_data(const EMData* const emdata, const float* const data,const int nx, const int ny, const int nz) {
	
	const EMData* previous = caller_cache[current_insert_idx];
	if (previous != 0) {
		previous->cuda_cache_lost_imminently();
		clear_item(current_insert_idx);
	}
		
	cudaArray *array = get_cuda_arrary_host(data,nx,ny,nz);
	
	rw_cache[current_insert_idx] = 0;
	caller_cache[current_insert_idx] = emdata;
	ro_cache[current_insert_idx] = array;
	
	unsigned int ret = current_insert_idx;
	current_insert_idx += 1;
	current_insert_idx %= cache_size; // Potentially inefficient to do this everytime, the alternative is an if statement. Which is faster?
	return ret;		
}


void EMData::CudaDeviceEMDataCache::clear_item(const unsigned int idx) {
	float** pointer_to_cuda_pointer = &rw_cache[idx];
	if  ( (*pointer_to_cuda_pointer) != 0) {
		cudaFree(*pointer_to_cuda_pointer);
	} 
	*pointer_to_cuda_pointer = 0;
	caller_cache[idx] = 0;
	
}

int EMData::CudaDeviceEMDataCache::copy_ro_to_rw_data(const unsigned int idx)
{
	return 0;
	//float** pointer_to_cuda_pointer = &rw_cache[idx];
	
}


#endif //EMAN2_USING_CUDA