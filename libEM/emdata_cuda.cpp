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
	if (cuda_cache_handle !=-1 && !(EMDATA_GPU_NEEDS_UPDATE & flags)) return cuda_cache.has_rw_data(cuda_cache_handle);
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

void EMData::bind_cuda_texture(const bool interp_mode) const {
	check_cuda_array_update();
	bind_cuda_array_to_texture(cuda_cache.get_ro_data(cuda_cache_handle),cuda_cache.get_ndim(cuda_cache_handle),interp_mode);
}

void EMData::unbind_cuda_texture() const {
	::unbind_cuda_texture(cuda_cache.get_ndim(cuda_cache_handle));
}

cudaArray* EMData::get_cuda_array() const {
	check_cuda_array_update();
	return cuda_cache.get_ro_data(cuda_cache_handle);
}

void EMData::check_cuda_array_update() const {
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
	flags |=  EMDATA_GPU_NEEDS_UPDATE| EMDATA_GPU_RO_NEEDS_UPDATE;
	cuda_cache_handle = -1;
}


bool EMData::gpu_operation_preferred() const {
	bool cpu = cpu_rw_is_current();
	bool gpu = gpu_rw_is_current();
	if ( cpu==0 &&  gpu==0 ) {		
		cout << (!(EMDATA_CPU_NEEDS_UPDATE & flags) && rdata != 0) << " " << (cuda_cache_handle !=-1 && !(EMDATA_GPU_NEEDS_UPDATE & flags) && cuda_cache.has_rw_data(cuda_cache_handle)) << endl;
		cout << "GPU flag " << !(EMDATA_GPU_NEEDS_UPDATE & flags) << endl;
		cout << "CPU flag " << !(EMDATA_CPU_NEEDS_UPDATE & flags) << endl;
		cout << "Rdata " << rdata << endl;
		cout << "Cuda handle " << cuda_cache_handle << endl;
		throw UnexpectedBehaviorException("Neither the CPU and GPU data are current");
	}
	if (gpu) return true;
	return false;
}

EMData* EMData::calc_ccf_cuda( EMData*  image, bool use_texturing,bool center ) const {
	
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
	 	emdata_processor_correlation_texture(&left,center);
	} else {
		EMDataForCuda right = ((EMData*)d["with"])->get_data_struct_for_cuda();
		emdata_processor_correlation(&left,&right,center);
	}
	tmp->gpu_update();
	
// 	tmp->process_inplace("cuda.correlate",d);
// 	return tmp;
	if (with != 0 && image != this) {
		delete with;
		with = 0;
	}

	EMData* soln = tmp->do_ift_cuda(false);
	soln->gpu_update();
	delete tmp;
	tmp = 0;
	
	return soln;
}

EMData *EMData::make_rotational_footprint_cuda( bool unwrap)
{
	ENTERFUNC;
//
//	update_stat();
	//float edge_mean = get_edge_mean();
	float edge_mean = 0;
//	//if ( rot_fp != 0 && unwrap == true) {
//	//	return new EMData(*rot_fp);
//	//}
//	
//	//static EMData obj_filt;
//	//EMData* filt = &obj_filt;
//	//filt->set_complex(true);
//// 	Region filt_region;
//
//// 	if (nx & 1) {
//// 		LOGERR("even image xsize only");		throw ImageFormatException("even image xsize only");
//// 	}
//
	int cs = (((nx * 7 / 4) & 0xfffff8) - nx) / 2; // this pads the image to 1 3/4 * size with result divis. by 8

	Region r1;
	if (nz == 1) {
		r1 = Region(-cs, -cs, nx + 2 * cs, ny + 2 * cs);
	}
	else {
		r1 = Region(-cs, -cs, -cs, nx + 2 * cs, ny + 2 * cs, nz + 2 * cs);
	}
//	// It is important to set all newly established pixels around the boundaries to the mean
//	// If this is not done then the associated rotational alignment routine breaks, in fact
//	// everythin just goes foo. 
	EMData *clipped = get_clip(r1,edge_mean);
//// 	EMData *clipped = copy()
//	
//	// The filter object is nothing more than a cached high pass filter
//	// Ultimately it is used an argument to the EMData::mult(EMData,prevent_complex_multiplication (bool)) 
//	// function in calc_mutual_correlation. Note that in the function the prevent_complex_multiplication 
//	// set to true, which is used for speed reasons. 
//	if (filt->get_xsize() != clipped->get_xsize() +2-(clipped->get_xsize()%2) || filt->get_ysize() != clipped->get_ysize() ||
//		   filt->get_zsize() != clipped->get_zsize()) {
//		filt->set_size(clipped->get_xsize() + 2-(clipped->get_xsize()%2), clipped->get_ysize(), clipped->get_zsize());
//		filt->to_one();
//		filt->process_inplace("eman1.filter.highpass.gaussian", Dict("highpass", 1.5f/nx));
//	}
//	
	//cout << "Clip 1" << endl;
	EMData *mc = clipped->calc_ccf_cuda(clipped,false,true);
	//mc->sub(mc->get_edge_mean());
	//mc->process_inplace("xform.phaseorigin.tocenter");
	//mc->write_image("mc.hdf");
	if( clipped ) {
		delete clipped;
		clipped = 0;
	}

	Region r2;
	if (nz == 1) {
		r2 = Region(cs - nx / 4, cs - ny / 4, nx * 3 / 2, ny * 3 / 2);
	}
	else {
		r2 = Region(cs - nx / 4, cs - ny / 4, cs - nz / 4, nx * 3 / 2, ny * 3 / 2, nz * 3 / 2);
	}
	EMData* clipped_mc = mc->get_clip(r2);
	//clipped_mc->write_image("clipped_mc.hdf");
	if( mc ) {
		delete mc;
		mc = 0;
	}
//	
	EMData * result = NULL;
//
//	if (nz == 1) {
	if (!unwrap) {
		//clipped_mc->process_inplace("mask.sharp", Dict("outer_radius", -1, "value", 0));
		result = clipped_mc;
	}
	else {
		result = clipped_mc->unwrap();
		if( clipped_mc ) {
			delete clipped_mc;
			clipped_mc = 0;
		}
	}
	
	result->gpu_update();
	
	return result;
//	}
//	else {
//		// I am not sure why there is any consideration of non 2D images, but it was here
//		// in the first port so I kept when I cleaned this function up (d.woolford)
//		result = clipped_mc;
//	}
//
//	EXITFUNC;
//	if ( unwrap == true)
//	{ // this if statement reflects a strict policy of caching in only one scenario see comments at beginning of function block
//		
//		// Note that the if statement at the beginning of this function ensures that rot_fp is not zero, so there is no need
//		// to throw any exception
//		// if ( rot_fp != 0 ) throw UnexpectedBehaviorException("The rotational foot print is only expected to be cached if it is not NULL");
//		
//		// Here is where the caching occurs - the rot_fp takes ownsherhip of the pointer, and a deep copied EMData object is returned.
//		// The deep copy invokes a cost in terms of CPU cycles and memory, but prevents the need for complicated memory management (reference counting)
//		rot_fp = result;
//		return new EMData(*rot_fp);
//	}
//	else return result;
}

void EMData::set_gpu_rw_data(float* data, const int x, const int y, const int z) {
	if (cuda_cache_handle!=-1) {
		cuda_cache.replace_gpu_rw(cuda_cache_handle,data);
	} else {
		cuda_cache_handle = cuda_cache.store_rw_data(this,data);
	}
	nx = x; ny = y; nz = z;
	nxy = nx*ny;
	gpu_update();
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
	clear_current();
	
	float* cuda_rw_data = alloc_rw_data(nx,ny,nz);
	
	if (data != 0 ) { // If rdata is zero it means we're working exclusively on the GPU
		size_t num_bytes = nx*ny*nz*sizeof(float);
		cudaError_t error = cudaMemcpy(cuda_rw_data,data,num_bytes,cudaMemcpyHostToDevice);
		if ( error != cudaSuccess) throw UnexpectedBehaviorException( "CudaMemcpy (host to device) error:" + string(cudaGetErrorString(error)));	
	}
	
	return force_store_rw_data(emdata,cuda_rw_data);
}

int EMData::CudaDeviceEMDataCache::force_store_rw_data(const EMData* const emdata, float*  cuda_rw_data)
{	
	rw_cache[current_insert_idx] = cuda_rw_data;
	caller_cache[current_insert_idx] = emdata;
	ro_cache[current_insert_idx] = 0;
	
	int ret = current_insert_idx;
	current_insert_idx += 1;
	current_insert_idx %= cache_size; // Potentially inefficient to do this every time, the alternative is an if statement. Which is faster?
	return ret;
}

int EMData::CudaDeviceEMDataCache::store_rw_data(const EMData* const emdata, float* cuda_rw_data)
{	
	clear_current();
	
	int nx = emdata->get_xsize();
	int ny = emdata->get_ysize();
	int nz = emdata->get_zsize();
	size_t num_bytes = nx*ny*nz*sizeof(float);
	mem_allocated += num_bytes;
	
	return force_store_rw_data(emdata, cuda_rw_data);
}

void EMData::CudaDeviceEMDataCache::replace_gpu_rw(const int idx,float* cuda_rw_data)
{
	clear_item(idx); // The ro data goes out of date anyway
	
	const EMData* d = caller_cache[idx];
	int nx = d->get_xsize();
	int ny = d->get_ysize();
	int nz = d->get_zsize();
	size_t num_bytes = nx*ny*nz*sizeof(float);
	mem_allocated += num_bytes;
	
	rw_cache[current_insert_idx] = cuda_rw_data;
}

void EMData::CudaDeviceEMDataCache::clear_current() {
	const EMData* previous = caller_cache[current_insert_idx];
	if (previous != 0) {
		previous->cuda_cache_lost_imminently();
		clear_item(current_insert_idx);
	}
}

float* EMData::CudaDeviceEMDataCache::alloc_rw_data(const int nx, const int ny, const int nz) {
	float* cuda_rw_data;
	size_t num_bytes = nx*ny*nz*sizeof(float);

	cudaError_t error = cudaMalloc((void**)&cuda_rw_data,num_bytes);
	if ( error != cudaSuccess) throw UnexpectedBehaviorException( "cudaMalloc error:" + string(cudaGetErrorString(error)));	

	
//	float* testing;
//	size_t pitch;
//	cudaMallocPitch( (void**)&testing, &pitch, nx*sizeof(float), ny*nz);
//	cout << "The pitch of that malloc as " << pitch << endl;
//	cudaFree(testing);
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
		cudaError_t error = cudaFreeArray(ro_cache[idx]);
		if ( error != cudaSuccess) throw UnexpectedBehaviorException( "CudaFreeArray error " + string(cudaGetErrorString(error)));
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
		cudaError_t error = cudaFree(rw_cache[idx]);
		if ( error != cudaSuccess)
			throw UnexpectedBehaviorException( "CudaFree error " + string(cudaGetErrorString(error)));
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
		cudaError_t error = cudaMemcpy3D(&copyParams);
		if ( error != cudaSuccess)
			throw UnexpectedBehaviorException( "Copying device array to device pointer - CudaMemcpy3D error : " + string(cudaGetErrorString(error)));
		
	} else if ( ny > 1 ) {
		cudaError_t error = cudaMemcpyFromArray(cuda_rw_data,ro_cache[idx],0,0,num_bytes,cudaMemcpyDeviceToDevice);
		if ( error != cudaSuccess)
			throw UnexpectedBehaviorException( "Copying device array to device pointer - cudaMemcpyFromArray error : " + string(cudaGetErrorString(error)));
	} else throw UnexpectedBehaviorException("Cuda infrastructure has not been designed to work on 1D data");
	
	
	rw_cache[idx] = cuda_rw_data;
}

void EMData::CudaDeviceEMDataCache::clear_item(const int idx) {
	
	if  ( rw_cache[idx] != 0) {
		mem_allocated -= get_emdata_bytes(idx);
		cudaError_t error = cudaFree(rw_cache[idx]);
		if ( error != cudaSuccess)
			throw UnexpectedBehaviorException( "CudaFree error " + string(cudaGetErrorString(error)));
	}
	rw_cache[idx] = 0;
	
	if  ( ro_cache[idx] != 0) {
		mem_allocated -= get_emdata_bytes(idx);
		cudaError_t error = cudaFreeArray(ro_cache[idx]);
		if ( error != cudaSuccess) throw UnexpectedBehaviorException( "CudaFreeArray error " + string(cudaGetErrorString(error)));

	}
	ro_cache[idx] = 0;
	
	caller_cache[idx] = 0;
	
// 	cout << "Allocation went down, it is currently " << (float)mem_allocated/1000000.0f << " MB " << endl;
	
}



#endif //EMAN2_USING_CUDA