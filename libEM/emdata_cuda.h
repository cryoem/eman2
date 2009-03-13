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

#ifndef eman__emdatacuda_h__
#define eman__emdatacuda_h__ 1

#ifdef EMAN2_USING_CUDA

public:

	/** Bind either a 2D or a 3D global scope CUDA texture which may be accessed using tex3D/tex2D
	 * Automatically binds to the correct dimensionality using the dimenisons of this EMData object
	 * Presently the 3D texture is called "tex" and the 2D texture is called "tex2d", these are currently
	 * defined in cuda_util.cu. The number of available textures and naming is likely to change
	 * to accommodate the need for accessing many textures using a single kernel
	 * @param interp_mode if true the texture will be bound using the cudaFilterModeLinear filtermode,  otherwise cudaFilterModePoint is usedddd
	 */
	void bind_cuda_texture(const bool interp_mode =true) const;
	
	void unbind_cuda_texture() const;
	
	/** Get the cuda device pointer to the raw data
	 * May cause an allocation, reflecting the lazy strategy employed in EMAN2
	 * Will copy the CPU rdata pointer if it is up to date.
	 * @return a cuda device allocated float pointer
	 */
	float* get_cuda_data() const;
	
	/** A convenient way to get the EMData object as an EMDataForCuda struct
	 * @return an EMDataForCuda struct storing vital information
	 */
	inline EMDataForCuda get_data_struct_for_cuda() const { 
		EMDataForCuda tmp = {get_cuda_data(),nx,ny,nz};
		return tmp;
	}
	
	/** Calc the correlation of this EMData with another using the CUDA device
	 * Does this using Fourier methods
	 * @param image the image to perform the correlation operation with
	 * @return a real space correlation image
	 */
	EMData* calc_ccf_cuda(EMData* image, bool use_texturing ) const;
	
	/** Multiply the image by a constant value using the CUDA device
	 * @param val the amount by which to multiply each pixel in the image
	 */
	void mult_cuda(const float& val);
	
	EMData* unwrap_cuda(int r1 = -1, int r2 = -1, int xs = -1, int dx = 0,
							   int dy = 0, bool do360 = false) const;
	
	/** Explicitly register that the raw data on the GPU has changed in some/any way.
	 * An important part of the EMAN2 device/host framework.
	 * This will set a flag telling the EMData object that both the CPU raw data pointer and
	 * read only array on the CUDA device are out of date and need to be updated at a later point.
	 * Also sets a flag dictating that the cached statistics of the image need to be updated.
	 */
	inline void gpu_update() {
		flags |= EMDATA_NEEDUPD | EMDATA_CPU_NEEDS_UPDATE | EMDATA_GPU_RO_NEEDS_UPDATE;
	}
	
	
	/** Explicitly force a copy from the cuda device pointer to the CPU
	 * Originally added to facilitate testing only
	 */
	void copy_gpu_rw_to_cpu();

	void copy_cpu_to_gpu_rw();

	void copy_cpu_to_gpu_ro();

	void copy_gpu_rw_to_gpu_ro();

	void copy_gpu_ro_to_gpu_rw();
	
	
private:
	
	void set_gpu_rw_data(float* data, const int x, const int y, const int z) ;
	
	/** Check whether the CUDA-cached read-write version of the data pointer is current
	 * Used to double check before copying the cuda rw data. It might be the case that the
	 * cuda_cache_handle is non-zero but that the cuda rw is actually not available.
	 * @return true or false
	 */
	bool gpu_rw_is_current() const;
	
	bool cpu_rw_is_current() const;
	
	/** Check whether the CUDA-cached read-only version of the data pointer is current
	 * Used to double check before copying the cuda ro data. It might be the case that the
	 * cuda_cache_handle is non-zero but that the cuda ro is actually not available.
	 * @return true or false
	 */
	bool gpu_ro_is_current() const;
	
	void check_cuda_array_update() const;
	cudaArray* get_cuda_array() const;
		
	/** Free CUDA memory
	 */
	void free_cuda_memory() const;
	
	/// A handle which may used to retrieve the device float pointer from the CudaDeviceEMDataCache using its [] operator
	mutable int cuda_cache_handle;
	
	/** CudaDeviceEMDataCache stores the CUDA device pointers and cudaArray pointers associated with EMData objects
	 * The cudaArray pointers may be bound to CUDA textures and used as read only memory. They are sometimes referred to as read-only or ro data.
	 * The cuda device pointers may be used for general GPU applications. They are sometimes referred to as read-write or rw data.
	 * This is a "snake cache" that eats its own tail once the available slots are all taken. 
	 * In the event of a random slot becoming free it is simply left that way, waiting empty until the snake head devours it.
	 * @author David Woolford
	 * @date February 2009
	*/
	class CudaDeviceEMDataCache {
		/** Class EMData is a friend to reflect the tight coupling of these two classes
		 * The EMData's use of this class' protected functions is somewhat specialized -
		 * It is presently hard to envisage any other class calling these functions directly. 
		 * It could be possible to move some of the code in EMData::get_cuda_data
		 * and EMData::bind_cuda_texture into this object, however it remains as is because
		 * most of that code is based on the EMData.flags member variable.
		 */
		friend class EMData;
	public:
		/** Constructor
		 * @param size the size of the cache
		 */
		CudaDeviceEMDataCache(const int size);
		
		/** Destructor, frees any non zero CUDA memory 
		 */
		~CudaDeviceEMDataCache();
	protected:
		
		/** Cache a read-write (CUDA device pointer) version of the  data and return an index reflecting its position in the cache
		 * @param emdata the EMData object making the call
		 * @param data the raw data pointer of the EMData object which is making the call
		 * @param nx the length of the x dimension of the raw data pointer
		 * @param ny the length of the y dimension of the raw data pointer
		 * @param nz the length of the z dimension of the raw data pointer
		 * @return the index which can be used to retrieve the CUDA device pointer using get_rw_data
		 */
		int cache_rw_data(const EMData* const emdata, const float* const data,const int nx, const int ny, const int nz);
		
		/** Cache a read-only (cudaArray) version of the  data and return an index reflecting its position in the cache
		 * @param emdata the EMData object making the call
		 * @param data the raw data pointer of the EMData object which is making the call
		 * @param nx the length of the x dimension of the raw data pointer
		 * @param ny the length of the y dimension of the raw data pointer
		 * @param nz the length of the z dimension of the raw data pointer
		 * @return the index which can be used to retrieve the CUDA device pointer using get_ro_data
		 */
		int cache_ro_data(const EMData* const emdata, const float* const data,const int nx, const int ny, const int nz);
		
		/** Get the read-write (CUDA device) pointer from a specific location in the cache
		 * @param idx the index used to retrieve the cached data
		 * @return a CUDA device pointer that may be used for reading and writing
		 */
		inline float* get_rw_data(const int idx) const { return rw_cache[idx]; }
		
		/** Get the read-only (cudaArray) pointer from a specific location in the cache
		 * @param idx the index used to retrieve the cached data
		 * @return a cudaArray pointer that may be used for binding a CUDA texture
		 */
		inline cudaArray* get_ro_data(const int idx) const { return ro_cache[idx]; }
		
		/** Check whether the cache has read-write (CUDA device) data at the given index
		 * @param idx the cache index that will be checked
		 * @return true or false
		 */
		inline bool has_rw_data(const int idx) const { return (rw_cache[idx] != 0); }
		
		/** Check whether the cache has read-only (cudaArray) data at the given index
		 * @param idx the cache index that will be checked
		 * @return true or false
		 */
		inline bool has_ro_data(const int idx) const { return (ro_cache[idx] != 0); }
		
		/** Clear a cache slot at a given index
		 * Performs device memory freeing if stored data is non zero
		 * @param idx the index of the stored CUDA device pointer
		 */
		void clear_item(const int idx);
		
		/** Copy the read-write data to the read-only data at the given cache index
		 * This function is used by the EMData object when it knows the read-write data is present and up to date
		 * Deletes the old read-only data
		 *  @param idx the cache index defining where the copying operation will occur
		 */
		void copy_rw_to_ro(const int idx);
		
		/** Copy the read-only data to the read-write data at the given cache index
		 * This function is used by the EMData object when it knows the read-only data is present and up to date
		 * Deletes the old read-write data
		 *  @param idx the cache index defining where the copying operation will occur
		 */
		void copy_ro_to_rw(const int idx);
		
		/** Get the number of dimensions of the EMData object the is associated with a specific cached object
		  * @param idx the cache index
		  * @return the dimensionality of the EMData object that is associated with the cached data 
		  */
		inline int get_ndim(const int idx) {
			return caller_cache[idx]->get_ndim();
		}
	private:
		/** Allocate a CUDA device pointer using cudaMalloc
		 * @param nx the length of the x dimension
		 * @param ny the length of the y dimension
		 * @param nz the length of the z dimension
		 * @return the cuda malloced device pointer
		 */
		float* alloc_rw_data(const int nx, const int ny, const int nz);
		
		/// The size of the cache
		int cache_size;
		/// The current position of the "snake head" in terms of the cache_size
		int current_insert_idx;
		/// Keep track of how much memory has been allocated
		size_t mem_allocated;
		
		inline size_t get_emdata_bytes(const int idx) {
			// This function can be called at program exit, in which case the EMData objects could
			// Already be deallocated.
// 			try {
				const EMData* e = caller_cache[idx];
				return e->get_size()*sizeof(float);
// 			} except (...) {
// 				return 0;	
// 			}
		}
		
		void clear_current();
		int store_rw_data(const EMData* const emdata, float* cuda_rw_data);
		int force_store_rw_data(const EMData* const emdata, float*  cuda_rw_data);
		void replace_gpu_rw(const int handle, float* cuda_rw_data);
		
		/// The CUDA device rw pointer cache
		float** rw_cache;
		/// A cache of the objects that called the cache_data function, so EMData::cuda_cache_lost_imminently can be called if necessary
		const EMData** caller_cache;
		/// The CUDA  read only cache
		cudaArray ** ro_cache;
	};
	
	/// CudaDeviceEMDataCache is a friend because it calls cuda_cache_lost_imminently, a specialized function that should not be made public.
	friend class CudaDeviceEMDataCache;
	
	/** Called by the CudaDeviceEMDataCache just before it frees the associated cuda device pointer
	 * Internally the EMData object will update its CPU version of the raw data if it is out of data as a result of GPU processing.
	 * Currently based loosely on the assumption that the host will always have more memory than the device, this may change, but works for the time being
	 */
	void cuda_cache_lost_imminently() const;
	
	/// Cuda device pointer cache
	static CudaDeviceEMDataCache cuda_cache;
	
#endif // EMAN2_USING_CUDA
	
#endif //eman__emdatacuda_h__ 1

