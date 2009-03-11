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
#include <driver_functions.h>
#include <cuda.h>
#include "cuda/cuda_util.h"
#include "cuda/cuda_processor.h"
#include "cuda/cuda_emfft.h"

using namespace EMAN;
// Static init
EMData::CudaDeviceEMDataCache EMData::cuda_cache(100);

float* EMData::get_cuda_data() const {
	if (cuda_cache_handle==-1 || EMDATA_GPU_NEEDS_UPDATE & flags) {
		if (cuda_cache_handle != -1 && gpu_ro_is_current() ) {
			cuda_cache.copy_ro_to_rw(cuda_cache_handle);
		}
		cuda_cache_handle = cuda_cache.cache_rw_data(this,rdata,nx,ny,nz);
		flags &= ~EMDATA_GPU_NEEDS_UPDATE;
	}
	return cuda_cache.get_rw_data(cuda_cache_handle);
}

bool EMData::gpu_rw_is_current() const {
	if (cuda_cache_handle !=-1 || !(EMDATA_GPU_NEEDS_UPDATE & flags)) return cuda_cache.has_rw_data(cuda_cache_handle);
	else return false;
}

bool EMData::cpu_rw_is_current() const {
	if 	(!(EMDATA_CPU_NEEDS_UPDATE & flags) && rdata != 0) return true;
	return false;
}

bool EMData::gpu_ro_is_current() const {
	if (cuda_cache_handle !=-1 || !(EMDATA_GPU_RO_NEEDS_UPDATE & flags)) return cuda_cache.has_ro_data(cuda_cache_handle);
	else return false;
}

void EMData::bind_cuda_texture(const bool interp_mode) {
	check_cuda_array_update();
	bind_cuda_array_to_texture(cuda_cache.get_ro_data(cuda_cache_handle),cuda_cache.get_ndim(cuda_cache_handle),interp_mode);
}

cudaArray* EMData::get_cuda_array() {
	check_cuda_array_update();
	return cuda_cache.get_ro_data(cuda_cache_handle);
}

void EMData::check_cuda_array_update() {
	if (cuda_cache_handle==-1 || EMDATA_GPU_RO_NEEDS_UPDATE & flags) {
		if (cuda_cache_handle !=- 1 && gpu_rw_is_current() )  {
			cuda_cache.copy_rw_to_ro(cuda_cache_handle);
		} else {
			cuda_cache_handle = cuda_cache.cache_ro_data(this,rdata,nx,ny,nz);
		}
		flags &= ~EMDATA_GPU_RO_NEEDS_UPDATE;
	}
}

void EMData::cuda_cache_lost_imminently() const {
	get_data(); // This causes cuda memory to be copied to cpu memory
	cuda_cache_handle = -1;
}



void EMData::mult_cuda(const float& val) {
// 	Dict d("scale",(float)val);
// 	process_inplace("cuda.math.mult",d);
	EMDataForCuda tmp = get_data_struct_for_cuda();
	emdata_processor_mult(&tmp,val);
	gpu_update();
}

EMData* EMData::calc_ccf_cuda( EMData*  image, bool use_texturing ) {
	
	EMData* tmp;
	if (is_complex()) {
// 		cout << "Tmp is a copy of this" << endl;
		tmp = new EMData(*this);
	} else {
// 		cout << "Tmp is this fftd" << endl;
		tmp = do_fft_cuda();
	}
	
	Dict d;
	EMData* with = 0;
	if (image == this) {
		d["with"] = (EMData*) tmp;
	} else {
		if (!image->is_complex()) {
			int wnx = image->get_xsize(); int wny = image->get_ysize(); int wnz = image->get_zsize();
			if ( wnx != nx || wny != ny || wnz != nz ) {
				
				Region r;
				if (nz > 1) {
					r = Region((wnx-nx)/2, (wny-ny)/2, (wnz-nz)/2,nx,ny,nz);
				}
				else if (ny > 1) {
					r = Region((wnx-nx)/2, (wny-ny)/2,nx,ny);
				}
				else throw UnexpectedBehaviorException("Calc_ccf_cuda doesn't work on 1D images");
				EMData* tmp = image->get_clip(r);
				with = tmp->do_fft_cuda();
				delete tmp;
			}else {
				cout << "With is the fft of the input" << endl;
				with = image->do_fft_cuda();
			}
			d["with"] = (EMData*) with;
		} else {
	// 		cout << "With is the input image" << endl;
			d["with"] = (EMData*)image;
		}
	}
	
	
	EMDataForCuda left = tmp->get_data_struct_for_cuda();
	if (use_texturing) {
	 	((EMData*)d["with"])->bind_cuda_texture(false);
	 	emdata_processor_correlation_texture(&left);
	} else {
		EMDataForCuda right = ((EMData*)d["with"])->get_data_struct_for_cuda();
		emdata_processor_correlation(&left,&right);
	}
	tmp->gpu_update();
	
// 	tmp->process_inplace("cuda.correlate",d);
// 	return tmp;
	if (with != 0 && image != this) {
		delete with;
		with = 0;
	}

	EMData* soln = tmp->do_ift_cuda();
	soln->gpu_update();
	delete tmp;
	tmp = 0;
	
	return soln;
}

void EMData::free_cuda_memory() const {
	if (cuda_cache_handle!=-1) {
		cuda_cache.clear_item(cuda_cache_handle);
		cuda_cache_handle = -1;
	}
}

/// THIS functoin was only ever meant for testing. Just use get_data() instead, where it's all automatic
void EMData::copy_gpu_rw_to_cpu() {
	get_data();
}

void EMData::copy_cpu_to_gpu_rw() {
	get_cuda_data();
}

void EMData::copy_cpu_to_gpu_ro() {
	get_cuda_array();
}

void EMData::copy_gpu_rw_to_gpu_ro() {
	cuda_cache.copy_rw_to_ro(cuda_cache_handle);
}

void EMData::copy_gpu_ro_to_gpu_rw() {
	cuda_cache.copy_ro_to_rw(cuda_cache_handle);
}



EMData::CudaDeviceEMDataCache::CudaDeviceEMDataCache(const int size) : cache_size(size), current_insert_idx(0), mem_allocated(0)
{
	device_init();
	rw_cache = new float *[cache_size];
	caller_cache = new const EMData*[cache_size];
	ro_cache = new cudaArray *[cache_size];
	
	for(int i = 0; i < cache_size; ++ i ) {
		rw_cache[i] = 0;
		caller_cache[i] = 0;
		ro_cache[i] = 0;
	}
}

EMData::CudaDeviceEMDataCache::~CudaDeviceEMDataCache()
{
	for (int i = 0; i < cache_size; i++) {
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
	// This might need some thinking
	cleanup_cuda_fft_dd_plan_cache();
}

int EMData::CudaDeviceEMDataCache::cache_rw_data(const EMData* const emdata, const float* const data,const int nx, const int ny, const int nz)
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
	
	int ret = current_insert_idx;
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
	mem_allocated += num_bytes;
// 	cout << "Allocation went up, it is currently " << (float)mem_allocated/1000000.0f << " MB " << endl;
	return cuda_rw_data;
	
}

int EMData::CudaDeviceEMDataCache::cache_ro_data(const EMData* const emdata, const float* const data,const int nx, const int ny, const int nz) {
	const EMData* previous = caller_cache[current_insert_idx];
	if (previous != 0) {
		previous->cuda_cache_lost_imminently();
		clear_item(current_insert_idx);
	}
		
	cudaArray *array = get_cuda_array_host(data,nx,ny,nz);
	if (array != 0) {
		mem_allocated += nx*ny*nz*sizeof(float);
// 		cout << "Allocation went up, it is currently " << (float)mem_allocated/1000000.0f << " MB " << endl;
		rw_cache[current_insert_idx] = 0;
		caller_cache[current_insert_idx] = emdata;
		ro_cache[current_insert_idx] = array;
		
		int ret = current_insert_idx;
		current_insert_idx += 1;
		current_insert_idx %= cache_size; // Potentially inefficient to do this everytime, the alternative is an if statement. Which is faster?
		return ret; 
	}
	else {
		throw BadAllocException("The allocation of the CUDA array failed");
	}
}


void  EMData::CudaDeviceEMDataCache::copy_rw_to_ro(const int idx) {
	if (rw_cache[idx] == 0) throw UnexpectedBehaviorException("Can not update RO CUDA data: RW data is null.");

	if (ro_cache[idx] != 0)  {
		cudaFreeArray(ro_cache[idx]);
		ro_cache[idx] = 0;
	}
	
	const EMData* d = caller_cache[idx];
	int nx = d->get_xsize();
	int ny = d->get_ysize();
	int nz = d->get_zsize();
	
	cudaArray *array = get_cuda_array_device(rw_cache[idx],nx,ny,nz);
	if (array == 0) throw BadAllocException("The allocation of the CUDA array failed");
	ro_cache[idx] = array;
}

void  EMData::CudaDeviceEMDataCache::copy_ro_to_rw(const int idx) {
	if (ro_cache[idx] == 0) throw UnexpectedBehaviorException("Can not update RW CUDA data: RO data is null.");

	if (rw_cache[idx] != 0)  {
		cudaFree(rw_cache[idx]);
		rw_cache[idx] = 0;
	}

	const EMData* d = caller_cache[idx];
	int nx = d->get_xsize();
	int ny = d->get_ysize();
	int nz = d->get_zsize();
	size_t num_bytes = nx*ny*nz*sizeof(float);
	
	float* cuda_rw_data = alloc_rw_data(nx,ny,nz);
	
	if (nz > 1) {
		cudaExtent extent;
		extent.width  = nx;
		extent.height = ny;
		extent.depth  = nz;
		cudaMemcpy3DParms copyParams = {0};
		copyParams.srcArray   = ro_cache[idx];
		copyParams.dstPtr = make_cudaPitchedPtr((void*)cuda_rw_data, extent.width*sizeof(float), extent.width, extent.height);
		copyParams.extent   = extent;
		copyParams.kind     = cudaMemcpyDeviceToDevice;
		cudaMemcpy3D(&copyParams);
		
	} else if ( ny > 1 ) {
		cudaMemcpyFromArray(cuda_rw_data,ro_cache[idx],0,0,num_bytes,cudaMemcpyDeviceToDevice);
		
	} else throw UnexpectedBehaviorException("Cuda infrastructure has not been designed to work on 1D data");
	
	
	rw_cache[idx] = cuda_rw_data;
}

void EMData::CudaDeviceEMDataCache::clear_item(const int idx) {
	float** pointer_to_cuda_pointer = &rw_cache[idx];
	if  ( (*pointer_to_cuda_pointer) != 0) {
		mem_allocated -= get_emdata_bytes(idx);
		cudaFree(*pointer_to_cuda_pointer);
	}
	*pointer_to_cuda_pointer = 0;
	
	cudaArray** pointer_to_cuda_ro = &ro_cache[idx];
	if  ( (*pointer_to_cuda_ro) != 0) {
		mem_allocated -= get_emdata_bytes(idx);
		cudaFreeArray(*pointer_to_cuda_ro);
	}
	*pointer_to_cuda_ro = 0;
	
	caller_cache[idx] = 0;
	
// 	cout << "Allocation went down, it is currently " << (float)mem_allocated/1000000.0f << " MB " << endl;
	
}



#endif //EMAN2_USING_CUDA