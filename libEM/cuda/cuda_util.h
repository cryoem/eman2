
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


/** Texture binding
 */
void bind_cuda_texture(const int);


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




