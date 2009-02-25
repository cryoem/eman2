
#ifndef eman__cuda_util_h__
#define eman__cuda_util_h__ 1

// Various utility functions
/** Initialize the cuda device
 * Can be called any number of times but the actual initialization occurs only the first time
 */
void device_init();

/** Get the stored cuda arrary corresponding to the input arguments
 */
int get_cuda_array_handle(const float* ,const int,const int,const int, void*);

int delete_cuda_array(const int idx);

int delete_cuda_memory(float *p);



/** Texture binding
 */
void bind_cuda_texture(const int);

void* cuda_malloc(const size_t size);
void cuda_free(void*);
void cuda_memset(void*,int value, size_t size);
void cuda_memcpy(void* dst,const void* const, size_t count);

/** Wrapping memory allocation functions
 */

void cuda_memcpy_host_to_device(const void* const host_rdata, void* device_rdata, const size_t num_bytes );
void cuda_malloc_device(void** device_rdata, const size_t num_bytes);
void cuda_free_device(void* device_rdata);

/** A struct for passing EMData objects to and from things like processors
*/
struct EMDataForCuda {
	float * data; // Cuda device pointer
	int nx; // Number of pixels in the x dimension
	int ny; // Number of pixels in the y dimension
	int nz; // Number of pixels in the z dimension
};


int* calc_max_location_wrap_cuda(const EMDataForCuda* data, const int maxdx, const int maxdy, const int maxdz);

void cut_slice_cuda_(const EMDataForCuda* data,const float* const);


#endif // eman__cuda_util_h__




