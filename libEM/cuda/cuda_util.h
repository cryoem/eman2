
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

struct CudaPeakInfo {
	int px;
	int py;
	int pz;
	float peak;
};

bool copy_to_array(const float * data, cudaArray * array, const int nx, const int ny, const int n, const cudaMemcpyKind memkindz);

//int* calc_max_location_wrap_cuda(const EMDataForCuda* data, const int maxdx, const int maxdy, const int maxdz);

//void cut_slice_cuda_(const EMDataForCuda* data,const float* const);

cudaArray* get_cuda_array(const int nx, const int ny, const int nz);

void bind_cuda_array_to_textureA( const cudaArray* const array, const int ndims, const bool interp_mode);

void unbind_cuda_textureA(const int ndims);

void bind_cuda_array_to_textureB( const cudaArray* const array, const int ndims, const bool interp_mode);

void unbind_cuda_textureB(const int ndims);

float get_edgemean_cuda(const float* data, const int nx, const int ny, const int nz);

void to_value_cuda(float* data, const float value, const int nx, const int ny, const int nz);

#endif // eman__cuda_util_h__




