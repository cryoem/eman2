
#ifndef eman__cuda_emfft_h__
#define eman__cuda_emfft_h__ 1

void init_cuda_emfft_cache();
void cleanup_cuda_emfft_cache();

int cuda_fft_real_to_complex_1d(float *real_data, float *complex_data, int n);
int cuda_fft_complex_to_real_1d(float *complex_data, float *real_data, int n);
int cuda_fft_real_to_complex_nd(float *real_data, float *complex_data, int nx, int ny, int nz);
int cuda_fft_complex_to_real_nd(float *complex_data, float *real_data, int nx, int ny, int nz);

#endif //eman__cuda_emfft_h__

