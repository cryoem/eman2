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

	/** Class CudaDataLock is a very basic object that you can use to
	 * temporarily 'lock' the cached cuda data (that is associated with a particular
	 * EMData), thus preventing it from being deleted inadvertently (see the CudaCache).
	 * When an object of this type is destructed it automatically unlocks the
	 * associated data.
	 * The following code demontrates use of this object
	 *@code
	 *    {//make a temporary scope
     *    EMData e("something.mrc");
     *    e.set_gpu_rw_current(); // Copy data to GPU
     *    CudaDataLock tmp_lock(&e); // GPU data is locked
     *    EMData f("somethingelse.mrc");
     * 	  f.set_gpu_rw_current(); // This could have inadvertently deleted e's GPU data if it was not locked
     *    // Cuda processing
     *    }// tmp_lock is destructed as it goes out of scope ( which unlocks e's CUDA data)
     *@endcode
     * You could achieve the same thing using the EMData member functions cuda_lock and cuda_unlock,
     * demonstrated as follows
     * @code
     *    e.cuda_lock()
     *    // processing
	 *    e.cuda_unlock()
	 * @endcode
	 *
	 * Thoughts: The need for locking is due solely to the CudaCache, which would indiscriminately devour
	 * whatever occupied the next slot if there were no mechanism for stopping it. Is there a better solution?
	 * @date April 6th 2009
	 * @author David Woolford
	 */
	class CudaDataLock {
	public:
		/** Constructor
		 * @param that the EMData object which will have its cuda data locked.
		 */
		CudaDataLock(const EMData* const that);

		/** Destructor
		 * Unlocks the associated CUDA data
		 */
		~CudaDataLock();

	private:
		/** Disallow assignment by declaring the method private
		 * Supporting this would add unnecessary complexity
		 */
		CudaDataLock& operator=(const CudaDataLock&);
		/** Disallow copy construction by declaring the method private
		 * Supporting this would add unnecessary complexity
		 */
		CudaDataLock(const CudaDataLock& );
		/** Disallow default construction
		 *  Supporting this does not make sense
		 */
		CudaDataLock();

		/// The CudaCache index of the associated EMData
		int data_cuda_handle;
	};

public:

	/** Bind either a 2D or a 3D global scope CUDA texture which may be accessed using tex3D/tex2D
	 * Automatically binds to the correct dimensionality using the dimenisons of this EMData object
	 * Presently the 3D texture is called "tex" and the 2D texture is called "tex2d", these are currently
	 * defined in cuda_util.cu. The number of available textures and naming is likely to change
	 * to accommodate the need for accessing many textures using a single kernel
	 * @param interp_mode if true the texture will be bound using the cudaFilterModeLinear filtermode,
	 * otherwise cudaFilterModePoint is used
	 */
	void bind_cuda_texture(const bool interp_mode =true) const;

	/** Unbind the cuda texture that was bound by the call to bind_cuda_texture
	 * Should be called once the texture is no longer needed
	 */
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
	EMDataForCuda get_data_struct_for_cuda() const;

	/** Calc the correlation of this EMData with another using the CUDA device
	 * Does this using Fourier methods
	 * @param image the image to perform the correlation operation with
	 * @param use_texturing
	 * @param center
	 * @return a real space correlation image
	 */
	EMData* calc_ccf_cuda(EMData* image, bool use_texturing,bool center=false ) const;

	EMData* calc_ccfx_cuda( EMData * const with, int y0=0, int y1=-1, bool no_sum=false);

	EMData * make_rotational_footprint_cuda( bool unwrap=true);

	/** Explicitly register that the raw data on the GPU has changed in some/any way.
	 * An important part of the EMAN2 device/host framework.
	 * This will set a flag telling the EMData object that both the CPU raw data pointer and
	 * read only array on the CUDA device are out of date and need to be updated at a later point.
	 * Also sets a flag dictating that the cached statistics of the image need to be updated.
	 */
	inline void gpu_update() {
		flags |= EMDATA_NEEDUPD | EMDATA_CPU_NEEDS_UPDATE | EMDATA_GPU_RO_NEEDS_UPDATE;
	}

	/** Get the column sum as an EMData object using CUDA
	 * Return object exists solely on the GPU
	 * @exception ImageDimensionException if this image is not 2D
	 * @return an EMData object that stores the column sum of this image and exists on the GPU
	 */
	EMData* column_sum_cuda() const;


	/** Explicitly force a copy from the cuda device pointer to the CPU
	 * Originally added to facilitate testing only
	 */
	void copy_gpu_rw_to_cpu();

	void copy_cpu_to_gpu_rw();

	// A long term solution?
	inline void set_gpu_rw_current() const {
		get_cuda_data();
	}

	bool gpu_operation_preferred() const;

	void copy_cpu_to_gpu_ro();

	void copy_gpu_rw_to_gpu_ro();

	void copy_gpu_ro_to_gpu_rw() const;

	void copy_gpu_ro_to_cpu() const;

	void print_this() const { cout << "this " << this << " " << cuda_cache_handle << endl; }




	/** Lock the associated cuda data
	 * A method for prevented the deletion of associated CUDA memory
	 * Note that this does not prevent its alteration (be it rw or ro)
	 */
	void cuda_lock() const;

	/** Unlock the assocated cuda data
	 * Tells the CudaCache that it no longer needs to prevent deletion
	 * of GPU memory associated with this object.
	 * Used exclusively in conjuntion with cuda_unlock
	 */
	void cuda_unlock() const;

	/** Check whether the CUDA-cached read-write version of the data pointer is current
	 * Used to double check before copying the cuda rw data. It might be the case that the
	 * cuda_cache_handle is non-zero but that the cuda rw is actually not available.
	 * @return true or false
	 */
	bool gpu_rw_is_current() const;

	/** Check whether the cpu data is both available and up-to-date. It may be the case
	 * that the CPU data is null (i.e. GPU operations are being used exclusively). Alternatively
	 * it may be that the EMDATA_CPU_NEEDS_UPDATE flag is toggled, which means
	 * the GPU version of the image is the current one and it is different to that which
	 * is stored on the CPU
	 */
	bool cpu_rw_is_current() const;


	void set_gpu_rw_data(float* data, const int x, const int y, const int z);
private:
	/** Get the cuda handle
	 * This is the entry index in the CudaCache that contains
	 * the associated GPU ro/rw data.
	 */
	int get_cuda_handle() const { return cuda_cache_handle; };


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

	/// A handle which may used to retrieve the device float pointer from the CudaCache using its [] operator
	mutable int cuda_cache_handle;


	/** CudaCache stores the CUDA device pointers and cudaArray pointers associated with EMData objects
	 * The cudaArray pointers may be bound to CUDA textures and used as read only memory. They are sometimes referred to as read-only or ro data.
	 * The cuda device pointers may be used for general GPU applications. They are sometimes referred to as read-write or rw data.
	 * This is a "snake cache" that eats its own tail once the available slots are all taken.
	 * In the event of a random slot becoming free it is simply left that way, waiting empty until the snake head devours it.
	 * Slots may be locked to prevent the snake head from consuming them using the lock member function
	 * Slots that are locked must be unlocked using the unlock member function. Note that when an EMData
	 * object is destructed it calls unlock, so calling unlock after calling lock is not strictly necessary. But it IS recommended as
	 * a good practice.
	 * @author David Woolford
	 * @date February 2009
	*/
	class CudaCache {
		/** Class EMData is a friend to reflect the tight coupling of these two classes
		 * The EMData's use of this class' protected functions is somewhat specialized -
		 * It is presently hard to envisage any other class calling these functions directly.
		 * It could be possible to move some of the code in EMData::get_cuda_data
		 * and EMData::bind_cuda_texture into this object, however it remains as is because
		 * most of that code is based on the EMData::flags member variable.
		 */
		friend class EMData;
		friend class CudaDataLock;
	public:
		/** Constructor
		 * @param size the size of the cache
		 */
		CudaCache(const int size);

		/** Destructor, indiscriminately frees any non zero CUDA memory
		 */
		~CudaCache();
	protected:
		/// I separated functions into protected and private groups - protected functions are those
		/// which are called by the (friend) EMData object. Private functions are called internally by
		/// member functions exclusively. I know this is not conventional, but it does document something usefully.

		/** Cache a read-write (CUDA device pointer) version of the  data and return an index reflecting its position in the cache
		 * @param emdata the EMData object making the call
		 * @param data the raw data pointer of the EMData object which is making the call
		 * @param nx the length of the x dimension of the raw data pointer
		 * @param ny the length of the y dimension of the raw data pointer
		 * @param nz the length of the z dimension of the raw data pointer
		 * @return the index which can be used to retrieve the CUDA device pointer using get_rw_data
		 */
		int cache_rw_data(const EMData* const emdata, const float* const data,const int nx, const int ny, const int nz);


		/** A way of storing CUDA allocated rw memory directly
		 * called from EMData::set_gpu_rw_data
		 * @param emdata the associated EMData object
		 * @param cuda_rw_data device allocated rw memory
		 * @return the index where the data was stored
		 */
		int store_rw_data(const EMData* const emdata, float* cuda_rw_data);

		/** A way of replacing CUDA allocated rw memory directly
		 * called from EMData::set_gpu_rw_data
		 * @param emdata the associated EMData object
		 * @param handle the index of the slot that will have its gpu rw data replaced
		 */
		void replace_gpu_rw(const int handle, float* cuda_rw_data);

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
		inline bool has_rw_data(const int idx) const {
			if (idx < 0 || idx >= cache_size) throw InvalidValueException(idx,"The idx is beyond the cache size");
			return (rw_cache[idx] != 0);
		}

		/** Check whether the cache has read-only (cudaArray) data at the given index
		 * @param idx the cache index that will be checked
		 * @return true or false
		 */
		inline bool has_ro_data(const int idx) const {
			if (idx < 0 || idx >= cache_size) throw InvalidValueException(idx,"The idx is beyond the cache size");
			return (ro_cache[idx] != 0);
		}

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

		/** Copy the read-only data to given the CPU read-write data
		 * This is primary for testing purposes
		 * @param idx the cache index defining where the copying operation will occur
		 * @param data the cpu allocated data that will be copied into
		 */
		void copy_ro_to_cpu(const int idx,float* data);

		/** Get the number of dimensions of the EMData object the is associated with a specific cached object
		  * @param idx the cache index
		  * @return the dimensionality of the EMData object that is associated with the cached data
		  */
		inline int get_ndim(const int idx) {
			if (idx < 0 || idx >= cache_size) throw InvalidValueException(idx,"The idx is beyond the cache size");
			return caller_cache[idx]->get_ndim();
		}

		/** Lock the data stored at the given index, preventing it from inadvertent deletion
		 * Locking that same item more than once is allowable. The lock is really an index that is incremented
		 * when this function is called.
		 * @param idx the index corresponding the cache slot that will be locked
		 * @exception InvalidValueException if the idx is outside the limits of the internal cache
		 */
		void lock(const int idx);

		/** Unlock data stored at the given index, signalling that it can be deleted by the traveling snake head.
		 * NOTE that this may not truly unlock the slot, for all that really happens is an integer that tallies the number
		 * of times lock has been called is decremented.
		 * @param the index corresponding to the cache slot that will be unlocked
		 * @exception InvalidValueException if the idx is outside the limits of the internal cache
		 */
		void unlock(const int idx);

		/** Print the contents of the cache to standard out
		 * Useful for debug
		 */
		void debug_print() const;
	private:
		/** Prevent copying it doesn't make sense */
		CudaCache(const CudaCache&);
		/** Prevent assignment it doesn't make sense */
		CudaCache& operator=(const CudaCache&);
		/** Prevent the default constructor it doesn't make sense*/
		CudaCache();

		/** Allocate a CUDA device pointer using cudaMalloc
		 * @param nx the length of the x dimension
		 * @param ny the length of the y dimension
		 * @param nz the length of the z dimension
		 * @return the cuda malloced device pointer
		 * @exception BadAllocException if the cudaMalloc call failed
		 */
		float* alloc_rw_data(const int nx, const int ny, const int nz);


		/** Get the size in bytes of the EMData object that is associated with the slot at a particular index
		 * @param idx the slot index
		 * @return the size in bytes
		 */
		inline size_t get_emdata_bytes(const int idx) {
			if (idx < 0 || idx >= cache_size) throw InvalidValueException(idx,"The idx is beyond the cache size");

			const EMData* e = caller_cache[idx];
			return e->get_size()*sizeof(float);
		}

		/** Called internally to make sure that the slot corresponding to current_insert_idx is empty/available
		 * May cause emptying of a slot.
		 * May changed current_insert_idx if locked slots are encountered
		 */
		void ensure_slot_space();

		/** Store the rw_data directly into the slot at current_insert_idx
		 * This is an encapsulation of a commonly used fragment of code.
		 * Called "blind" because it doesn't check to make sure the slot is empty,
		 * this is okay because the functions that call this have already called ensure_slot_space.
		 * If this function was used naively it could cause a memory leak.
		 * @param emdata the associated EMData object
		 * @param cuda_rw_data device allocated rw memory
		 * @return the index where the data was stored
		 */
		int blind_store_rw_data(const EMData* const emdata, float*  cuda_rw_data);

		/// The size of the cache
		int cache_size;
		/// The current position of the "snake head" in terms of the cache_size
		int current_insert_idx;
		/// Keep track of how much memory has been allocated
		size_t mem_allocated;
		/// The CUDA device rw pointer cache
		float** rw_cache;
		/// A cache of the objects that called the cache_data function, so EMData::cuda_cache_lost_imminently can be called if necessary
		const EMData** caller_cache;
		/// The CUDA  read only cache
		cudaArray ** ro_cache;
		/// An array storing locking flags - everytime a cache entry is locked the value is incremented
		vector<int> locked;
	};

	/// CudaCache is a friend because it calls cuda_cache_lost_imminently, a specialized function that should not be made public.
	friend class CudaCache;

	/// CudaDataLock is a friend because it needs to know the cuda_cache_handle
	friend class CudaDataLock;

	/** Called by the CudaCache just before it frees the associated cuda device pointer
	 * Internally the EMData object will update its CPU version of the raw data if it is out of data as a result of GPU processing.
	 * Currently based loosely on the assumption that the host will always have more memory than the device, this may change, but works for the time being
	 */
	void cuda_cache_lost_imminently() const;

	/// Cuda device pointer cache
	static CudaCache cuda_cache;

public:
	void cuda_cache_debug_print() const { cuda_cache.debug_print(); }

#endif // EMAN2_USING_CUDA

#endif //eman__emdatacuda_h__ 1

