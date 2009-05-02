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
EMData::CudaCache EMData::cuda_cache(100);

float* EMData::get_cuda_data() const {
	if (get_size() == 0 ) throw UnexpectedBehaviorException("The size of the data is 0?");
	if (cuda_cache_handle==-1 || EMDATA_GPU_NEEDS_UPDATE & flags) {
		if (cuda_cache_handle != -1 && gpu_ro_is_current() ) {
			cuda_cache.copy_ro_to_rw(cuda_cache_handle);
		} else {
			if (cuda_cache_handle !=-1 ) {
				cuda_cache.clear_item(cuda_cache_handle);
			}
			cuda_cache_handle = cuda_cache.cache_rw_data(this,rdata,nx,ny,nz);
			if (cuda_cache_handle == -1) throw;
		}
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
	if (cuda_cache_handle !=-1 && !(EMDATA_GPU_RO_NEEDS_UPDATE & flags)) return cuda_cache.has_ro_data(cuda_cache_handle);
	else return false;
}

void EMData::bind_cuda_texture(const bool interp_mode) const {
	check_cuda_array_update();
	cuda_cache.lock(cuda_cache_handle);
	bind_cuda_array_to_texture(cuda_cache.get_ro_data(cuda_cache_handle),cuda_cache.get_ndim(cuda_cache_handle),interp_mode);
}

void EMData::unbind_cuda_texture() const {
	::unbind_cuda_texture(cuda_cache.get_ndim(cuda_cache_handle));
	cuda_cache.unlock(cuda_cache_handle);
}

cudaArray* EMData::get_cuda_array() const {
	if (get_size() == 0 ) throw UnexpectedBehaviorException("The size of the data is 0?");
	check_cuda_array_update();
	return cuda_cache.get_ro_data(cuda_cache_handle);
}

void EMData::check_cuda_array_update() const {
	if (cuda_cache_handle==-1 || EMDATA_GPU_RO_NEEDS_UPDATE & flags) {
		if (cuda_cache_handle !=- 1 && gpu_rw_is_current() )  {
			cuda_cache.copy_rw_to_ro(cuda_cache_handle);
		} else {
			if (cuda_cache_handle !=-1 ) cuda_cache.clear_item(cuda_cache_handle);
			cuda_cache_handle = cuda_cache.cache_ro_data(this,rdata,nx,ny,nz);
			if (cuda_cache_handle >=50 ) throw InvalidValueException(cuda_cache_handle,"In get cuda data, the handle is strange");
			if (cuda_cache_handle == -1) throw;
		}
		flags &= ~EMDATA_GPU_RO_NEEDS_UPDATE;
	}
}

void EMData::cuda_cache_lost_imminently() const {
	//scout << "In cache lost " << cuda_cache_handle << " " << nx << " " << ny << " " << nz << endl;
	get_data(); // This causes cuda memory to be copied to cpu memory
	flags |=  EMDATA_GPU_NEEDS_UPDATE| EMDATA_GPU_RO_NEEDS_UPDATE;
	cuda_cache_handle = -1;
}
void EMData::cuda_lock() const {
	if (cuda_cache_handle == -1) throw UnexpectedBehaviorException("No cuda handle, can't lock");
	cuda_cache.lock(cuda_cache_handle);
	//cuda_cache.debug_print();
}
void EMData::cuda_unlock() const {
	//cout << " " << cuda_cache_handle << endl;
	//cuda_cache.debug_print();
	if (cuda_cache_handle == -1) throw UnexpectedBehaviorException("No cuda handle, can't lock");
	cuda_cache.unlock(cuda_cache_handle);
}
EMDataForCuda EMData::get_data_struct_for_cuda() const {
	EMDataForCuda tmp = {get_cuda_data(),nx,ny,nz};
	return tmp;
}

bool EMData::gpu_operation_preferred() const {
	bool cpu = cpu_rw_is_current();
	bool gpu = gpu_rw_is_current();
	if ( cpu==0 &&  gpu==0 ) {
		// This is what happens when set_size doesn't allocate
		return false;
//		cout << (!(EMDATA_CPU_NEEDS_UPDATE & flags) && rdata != 0) << " " << (cuda_cache_handle !=-1 && !(EMDATA_GPU_NEEDS_UPDATE & flags) && cuda_cache.has_rw_data(cuda_cache_handle)) << endl;
//		cout << "GPU flag " << !(EMDATA_GPU_NEEDS_UPDATE & flags) << endl;
//		cout << "CPU flag " << !(EMDATA_CPU_NEEDS_UPDATE & flags) << endl;
//		cout << "Rdata " << rdata << endl;
//		cout << "Cuda handle " << cuda_cache_handle << endl;
//		throw UnexpectedBehaviorException("Neither the CPU or GPU data are current");
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
				with = image->do_fft_cuda();
			}
			d["with"] = (EMData*) with;
		} else {
	// 		cout << "With is the input image" << endl;
			d["with"] = (EMData*)image;
		}
	}


	EMDataForCuda left = tmp->get_data_struct_for_cuda();
	CudaDataLock lock(tmp);
	if (use_texturing) {
	 	((EMData*)d["with"])->bind_cuda_texture(false);
	 	emdata_processor_correlation_texture(&left,center);
	 	((EMData*)d["with"])->unbind_cuda_texture();
	} else {
		EMDataForCuda right = ((EMData*)d["with"])->get_data_struct_for_cuda();
		CudaDataLock lock2((EMData*)d["with"]);
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
//	float edge_mean = get_edge_mean();
	float edge_mean = 0;
	CudaDataLock(this);
	if ( rot_fp != 0 && unwrap == true) {
		return new EMData(*rot_fp);
	}
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

	static EMData big_clip;
	int big_x = nx+2*cs;
	int big_y = ny+2*cs;
	int big_z = 1;
	if ( nz != 1 ) {
		big_z = nz+2*cs;
	}


	if ( big_clip.get_xsize() != big_x || big_clip.get_ysize() != big_y || big_clip.get_zsize() != big_z ) {
		big_clip.set_size_cuda(big_x,big_y,big_z);
		big_clip.get_cuda_data();
		big_clip.cuda_lock(); // Just lock for the entire duration of the program, it's static anyway...
	}
	big_clip.to_value(edge_mean);

	if (nz != 1) {
		big_clip.insert_clip(this,IntPoint(cs,cs,cs));
	} else  {
		big_clip.insert_clip(this,IntPoint(cs,cs,0));
	}
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
	EMData *mc = big_clip.calc_ccf_cuda(&big_clip,false,true);
	mc->sub(mc->get_edge_mean());

	static EMData sml_clip;
	int sml_x = nx * 3 / 2;
	int sml_y = ny * 3 / 2;
	int sml_z = 1;
	if ( nz != 1 ) {
		sml_z = nz * 3 / 2;
	}

	if ( sml_clip.get_xsize() != sml_x || sml_clip.get_ysize() != sml_y || sml_clip.get_zsize() != sml_z ) {
		sml_clip.set_size_cuda(sml_x,sml_y,sml_z);
		sml_clip.get_cuda_data();
		sml_clip.cuda_lock(); // Just lock for the entire duration of the program, it's static anyway...
	}
	if (nz != 1) {
		sml_clip.insert_clip(mc,IntPoint(-cs+nx/4,-cs+ny/4,-cs+nz/4));
	} else {
		sml_clip.insert_clip(mc,IntPoint(-cs+nx/4,-cs+ny/4,0));
	}

	EMData * result = NULL;

	if (!unwrap || nz != 1) {
		//clipped_mc->process_inplace("mask.sharp", Dict("outer_radius", -1, "value", 0));
		result = new EMData(sml_clip);
	}
	else {
		result = sml_clip.unwrap();
	}

	result->gpu_update();

	if ( unwrap == true)
	{ // this if statement reflects a strict policy of caching in only one scenario see comments at beginning of function block

		// Note that the if statement at the beginning of this function ensures that rot_fp is not zero, so there is no need
		// to throw any exception
		// if ( rot_fp != 0 ) throw UnexpectedBehaviorException("The rotational foot print is only expected to be cached if it is not NULL");

		// Here is where the caching occurs - the rot_fp takes ownsherhip of the pointer, and a deep copied EMData object is returned.
		// The deep copy invokes a cost in terms of CPU cycles and memory, but prevents the need for complicated memory management (reference counting)
		rot_fp = result;
		return new EMData(*rot_fp);
	}
	else return result;
}

EMData* EMData::calc_ccfx_cuda( EMData * const with, int y0, int y1, bool no_sum)
{
	ENTERFUNC;
// 	cout << "calc_ccfx cuda" << endl;
	if (!with) {
		LOGERR("NULL 'with' image. ");
		throw NullPointerException("NULL input image");
	}

	if (!EMUtil::is_same_size(this, with)) {
		LOGERR("images not same size: (%d,%d,%d) != (%d,%d,%d)",
			   nx, ny, nz,
			   with->get_xsize(), with->get_ysize(), with->get_zsize());
		throw ImageFormatException("images not same size");
	}
	if (get_ndim() > 2) {
		LOGERR("2D images only");
		throw ImageDimensionException("2D images only");
	}

	if (y1 <= y0) {
		y1 = ny;
	}

	if (y0 >= y1) {
		y0 = 0;
	}

	if (y0 < 0) {
		y0 = 0;
	}

	if (y1 > ny) {
		y1 = ny;
	}

	static int nx_device_fft = 0;
	static int ny_defice_fft = 0;
	static EMData f1;
	static EMData f2;
	static EMData rslt;

	int height = y1-y0;
	int width = (nx+2-(nx%2));
	if (width != nx_device_fft || height != ny_defice_fft ) {
		f1.set_size_cuda(width,height);
		f2.set_size_cuda(width,height);
		rslt.set_size_cuda(nx,height);
		nx_device_fft = width;
		ny_defice_fft = height;
	}

	{// Make a local scope so that the locks are destructed
		float * cd = get_cuda_data();
		CudaDataLock lock(this);
		float * f1cd = f1.get_cuda_data();
		CudaDataLock lock2(&f1);
		cuda_dd_fft_real_to_complex_1d(cd,f1cd,nx,height);
	}
	{// Make a local scope so that the locks are destructed
		float * wcd = with->get_cuda_data();
		CudaDataLock lock(this);
		float * f2cd = f2.get_cuda_data();
		CudaDataLock lock2(&f2);
		cuda_dd_fft_real_to_complex_1d(wcd,f2cd,nx,height);
	}

	EMDataForCuda left = f1.get_data_struct_for_cuda();
	CudaDataLock lock(&f1);

	bool use_texturing = false;
	bool center = false;
	if (use_texturing) {
	 	f2.bind_cuda_texture(false);
	 	emdata_processor_correlation_texture(&left,center);
	 	f2.unbind_cuda_texture();
	} else {
		EMDataForCuda right = f2.get_data_struct_for_cuda();
		CudaDataLock lock2(&f2);
		emdata_processor_correlation(&left,&right,center);
	}

	{// Make a local scope so that the locks are destructed
		float* rcd = rslt.get_cuda_data();
		CudaDataLock rlock(&rslt);
		float * f1cd = f1.get_cuda_data();
		CudaDataLock lock2(&f1);
		cuda_dd_fft_complex_to_real_1d(f1cd,rcd,nx,height);
	}

	if (no_sum) {
		rslt.gpu_update();
		EXITFUNC;
		return new EMData(rslt);
	}
	else {
		EXITFUNC;
		return rslt.column_sum_cuda();
	}

}

EMData* EMData::column_sum_cuda() const {
	ENTERFUNC;
	if (get_ndim() != 2) throw ImageDimensionException("Column sum cuda has been prgogrammed work exclusively with 2D data.");
	EMData *cf = new EMData();
	cf->set_size_cuda(nx, 1, 1);
	EMDataForCuda left = cf->get_data_struct_for_cuda();
	CudaDataLock llock(cf);
	bind_cuda_texture(false);
	emdata_column_sum(&left,ny);
	unbind_cuda_texture();
	cf->gpu_update();
	EXITFUNC;
	return cf;
}

void EMData::set_gpu_rw_data(float* data, const int x, const int y, const int z) {
	nx = x; ny = y; nz = z;
	nxy = nx*ny;
	if (cuda_cache_handle!=-1) {
		cuda_cache.replace_gpu_rw(cuda_cache_handle,data);
	} else {
		cuda_cache_handle = cuda_cache.store_rw_data(this,data);
	}
	gpu_update();
}

void EMData::free_cuda_memory() const {
// 	cout << "Death comes to " << this << " " << cuda_cache_handle << endl;
	if (cuda_cache_handle!=-1) {
		cuda_cache.clear_item(cuda_cache_handle);
		cuda_cache_handle = -1;
	}
}

/// THIS function was only ever meant for testing. Just use get_data() instead, where it's all automatic
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

void EMData::copy_gpu_ro_to_gpu_rw() const {
	cuda_cache.copy_ro_to_rw(cuda_cache_handle);
}

void EMData::copy_gpu_ro_to_cpu() const {
	cuda_cache.copy_ro_to_cpu(cuda_cache_handle,rdata);
}


EMData::CudaCache::CudaCache(const int size) : cache_size(size), current_insert_idx(0), mem_allocated(0), locked(size,0)
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

EMData::CudaCache::~CudaCache()
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

void EMData::CudaCache::lock(const int idx) {
	if (idx < 0 || idx >= cache_size) throw InvalidValueException(idx,"The idx is beyond the cache size");
	locked[idx] += 1;
// 	debug_print();
}
void EMData::CudaCache::unlock(const int idx) {
	if (idx < 0 || idx >= cache_size) throw InvalidValueException(idx,"The idx is beyond the cache size");
	if (locked[idx] == 0) {
// // 		cout << "Warning - unlocked something that didn't need it" << endl;
		return;

// 		 throw UnexpectedBehaviorException("Can't unlock, it wasn't locked!");
	}
	locked[idx] -=1;
}

int EMData::CudaCache::cache_rw_data(const EMData* const emdata, const float* const data,const int nx, const int ny, const int nz)
{
	ensure_slot_space();

	float* cuda_rw_data = alloc_rw_data(nx,ny,nz);

	if (data != 0 ) { // If rdata is zero it means we're working exclusively on the GPU
		size_t num_bytes = nx*ny*nz*sizeof(float);
		cudaError_t error = cudaMemcpy(cuda_rw_data,data,num_bytes,cudaMemcpyHostToDevice);
		if ( error != cudaSuccess) throw UnexpectedBehaviorException( "CudaMemcpy (host to device) error:" + string(cudaGetErrorString(error)));
	}

	return blind_store_rw_data(emdata,cuda_rw_data);
}

int EMData::CudaCache::blind_store_rw_data(const EMData* const emdata, float*  cuda_rw_data)
{
// 	debug_print();
	rw_cache[current_insert_idx] = cuda_rw_data;
	caller_cache[current_insert_idx] = emdata;
	ro_cache[current_insert_idx] = 0;

	int ret = current_insert_idx;
	current_insert_idx += 1;
	current_insert_idx %= cache_size; // Potentially inefficient to do this every time, the alternative is an if statement. Which is faster?
	if ( current_insert_idx > cache_size ) throw;// This is just for debug
// 	cout << "Inserted at " << ret  << " inc to " << current_insert_idx << " size " << get_emdata_bytes(ret)/sizeof(float) << endl;
	return ret;
}

int EMData::CudaCache::store_rw_data(const EMData* const emdata, float* cuda_rw_data)
{
	ensure_slot_space();

	int nx = emdata->get_xsize();
	int ny = emdata->get_ysize();
	int nz = emdata->get_zsize();
	size_t num_bytes = nx*ny*nz*sizeof(float);
	mem_allocated += num_bytes;

	return blind_store_rw_data(emdata, cuda_rw_data);
}

void EMData::CudaCache::debug_print() const {
	cout << "Cuda device cache debug. Total mem allocated: " << static_cast<float>(mem_allocated)/1000000.0 << "MB" << endl;
	for(int i = 0; i < cache_size; ++i) {
		int handle = -1;
		int nx = 0;
		int ny = 0;
		int nz = 0;
		if (caller_cache[i] != 0) {
			handle = caller_cache[i]->cuda_cache_handle;
			nx = caller_cache[i]->get_xsize();
			ny = caller_cache[i]->get_ysize();
			nz = caller_cache[i]->get_zsize();
		}
		cout << i << ": " << handle << " " << caller_cache[i] << " dims: " << nx << " " << ny << " " << nz << " locked: " << locked[i] << " rw " << rw_cache[i] << " ro " << ro_cache[i] << endl;
// 		}
	}
}

void EMData::CudaCache::replace_gpu_rw(const int idx,float* cuda_rw_data)
{
	//clear_item(idx); // The ro data goes out of date anyway
	if  ( rw_cache[idx] != 0) {
		mem_allocated -= get_emdata_bytes(idx);
		cudaError_t error = cudaFree(rw_cache[idx]);
		if ( error != cudaSuccess)
			throw UnexpectedBehaviorException( "CudaFree error : " + string(cudaGetErrorString(error)));
	}
	rw_cache[idx] = 0;

	const EMData* d = caller_cache[idx];
	int nx = d->get_xsize();
	int ny = d->get_ysize();
	int nz = d->get_zsize();
	size_t num_bytes = nx*ny*nz*sizeof(float);
	mem_allocated += num_bytes;

	rw_cache[idx] = cuda_rw_data;
}

void EMData::CudaCache::ensure_slot_space() {

	int checked_entries = 0;
	while ( checked_entries < cache_size) {
		const EMData* previous = caller_cache[current_insert_idx];
		if (previous != 0 ) {
			if ( locked[current_insert_idx] == 0 ) {
// 				cout << "Sending imminent lost sig " << current_insert_idx  << endl;
				previous->cuda_cache_lost_imminently();
// 				cout << "Clear..." << endl;
				clear_item(current_insert_idx);
				break;
			} else {
// 				cout <<  "Lucky it was locked! " << current_insert_idx << endl;
				current_insert_idx++;
				current_insert_idx %= cache_size;
// 				cout <<  "Incremented to " << current_insert_idx << endl;
				checked_entries++;
			}
		} else break; // There IS space!
	}

	if (checked_entries == cache_size) {
		throw UnexpectedBehaviorException("All of the data objects in the cuda cache are locked! There is no space.");
	}
}

float* EMData::CudaCache::alloc_rw_data(const int nx, const int ny, const int nz) {
	float* cuda_rw_data;
	size_t num_bytes = nx*ny*nz*sizeof(float);

	cudaError_t error = cudaMalloc((void**)&cuda_rw_data,num_bytes);
	if ( error != cudaSuccess) {
		debug_print();
		throw BadAllocException( "cudaMalloc error :" + string(cudaGetErrorString(error)));
	}


//	float* testing;
//	size_t pitch;
//	cudaMallocPitch( (void**)&testing, &pitch, nx*sizeof(float), ny*nz);
//	cout << "The pitch of that malloc as " << pitch << endl;
//	cudaFree(testing);

	mem_allocated += num_bytes;
// 	cout << "Allocation went up, it is currently " << (float)mem_allocated/1000000.0f << " MB " << endl;
	return cuda_rw_data;

}

int EMData::CudaCache::cache_ro_data(const EMData* const emdata, const float* const data,const int nx, const int ny, const int nz) {
	ensure_slot_space();

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
// 		cout << "Inserted at " << ret  << " inc to " << current_insert_idx << " size " << nx*ny*nz << endl;
		return ret;
	}
	else {
		throw BadAllocException("The allocation of the CUDA array failed");
	}
}


void  EMData::CudaCache::copy_rw_to_ro(const int idx) {
// 	cout << "Copy rw to ro " << idx << endl;
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

void  EMData::CudaCache::copy_ro_to_rw(const int idx) {
// 	cout << "Copy ro to rw " << idx << endl;
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


void  EMData::CudaCache::copy_ro_to_cpu(const int idx,float* data) {
	if (ro_cache[idx] == 0) throw UnexpectedBehaviorException("Can not update RW CUDA data: RO data is null.");
	if (data == 0) throw NullPointerException("The cpu data pointer is NULL in copy_ro_to_cpu");

	const EMData* d = caller_cache[idx];
	int nx = d->get_xsize();
	int ny = d->get_ysize();
	int nz = d->get_zsize();
	size_t num_bytes = nx*ny*nz*sizeof(float);

	if (nz > 1) {
		cudaExtent extent;
		extent.width  = nx;
		extent.height = ny;
		extent.depth  = nz;
		cudaMemcpy3DParms copyParams = {0};
		copyParams.srcArray   = ro_cache[idx];
		copyParams.dstPtr = make_cudaPitchedPtr((void*)data, extent.width*sizeof(float), extent.width, extent.height);
		copyParams.extent   = extent;
		copyParams.kind     = cudaMemcpyDeviceToHost;
		cudaError_t error = cudaMemcpy3D(&copyParams);
		if ( error != cudaSuccess)
			throw UnexpectedBehaviorException( "Copying device array to device pointer - CudaMemcpy3D error : " + string(cudaGetErrorString(error)));

	} else if ( ny > 1 ) {
		cudaError_t error = cudaMemcpyFromArray(data,ro_cache[idx],0,0,num_bytes,cudaMemcpyDeviceToHost);
		if ( error != cudaSuccess)
			throw UnexpectedBehaviorException( "Copying device array to device pointer - cudaMemcpyFromArray error : " + string(cudaGetErrorString(error)));
	} else throw UnexpectedBehaviorException("Cuda infrastructure has not been designed to work on 1D data");

}
void EMData::CudaCache::clear_item(const int idx) {
// 	debug_print();
	if  ( rw_cache[idx] != 0) {
		mem_allocated -= get_emdata_bytes(idx);
		cudaError_t error = cudaFree(rw_cache[idx]);
		if ( error != cudaSuccess)
			throw UnexpectedBehaviorException( "CudaFree error : " + string(cudaGetErrorString(error)));
	}
	rw_cache[idx] = 0;

	if  ( ro_cache[idx] != 0) {
		mem_allocated -= get_emdata_bytes(idx);
		cudaError_t error = cudaFreeArray(ro_cache[idx]);
		if ( error != cudaSuccess) throw UnexpectedBehaviorException( "CudaFreeArray error : " + string(cudaGetErrorString(error)));

	}
	ro_cache[idx] = 0;

	caller_cache[idx] = 0;

	locked[idx] = 0;
}


EMData::CudaDataLock::CudaDataLock(const EMData* const emdata) : data_cuda_handle(-1)
{
	emdata->set_gpu_rw_current();
	data_cuda_handle = emdata->cuda_cache_handle;
	EMData::cuda_cache.lock(data_cuda_handle);
}

EMData::CudaDataLock::~CudaDataLock() {
	EMData::cuda_cache.unlock(data_cuda_handle);
}


#endif //EMAN2_USING_CUDA
