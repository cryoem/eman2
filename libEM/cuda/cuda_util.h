
#ifndef eman__cuda_util_h__
#define eman__cuda_util_h__ 1

// Various utility functions
/** Initialize the cuda device
 * Can be called any number of times but the actual initialization occurs only the first time
 */
void device_init();

/** Get the stored cuda arrary corresponding to the input arguments
 */
int stored_cuda_array(const float* ,const int,const int,const int);

int delete_cuda_array(const int idx);

int delete_cuda_memory(float *p);

/** Texture binding
 */
void bind_cuda_texture(const int);

float* cuda_malloc_float(const int size);
void cuda_free(float*);
void cuda_memset(float*,int value, int size);

#endif // eman__cuda_util_h__