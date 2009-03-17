
#ifndef eman__cuda_util_h__
#define eman__cuda_util_h__ 1

// Various utility functions
/** Initialize the cuda device
 * Can be called any number of times but the actual initialization occurs only the first time
 */
void device_init();

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

cudaArray* get_cuda_array_host(const float * const data,const int nx, const int ny, const int nz);

cudaArray* get_cuda_array_device(const float * const data,const int nx, const int ny, const int nz);

void bind_cuda_array_to_texture( const cudaArray* const array, const int ndims, const bool interp_mode);

void unbind_cuda_texture(const int ndims);

void emdata_column_sum(const EMDataForCuda* sum_target,const int ny);

#endif // eman__cuda_util_h__




