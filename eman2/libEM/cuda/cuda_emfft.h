
#ifndef eman__cuda_emfft_h__
#define eman__cuda_emfft_h__ 1

//void init_cuda_fft_hh_plan_cache();
//void cleanup_cuda_fft_hh_plan_cache();

//int cuda_hh_fft_real_to_complex_1d(float *real_data, float *complex_data, int n, int batch);
//int cuda_hh_fft_complex_to_real_1d(float *complex_data, float *real_data, int n, int batch);
//int cuda_hh_fft_real_to_complex_nd(float *real_data, float *complex_data, int nx, int ny, int nz);
//int cuda_hh_fft_complex_to_real_nd(float *complex_data, float *real_data, int nx, int ny, int nz);

//void init_cuda_fft_dd_plan_cache();
//void cleanup_cuda_fft_dd_plan_cache();

//int cuda_dd_fft_real_to_complex_1d(float *real_data, float *complex_data, int n, int batch);
//int cuda_dd_fft_complex_to_real_1d(float *complex_data, float *real_data, int n, int batch);
int cuda_dd_fft_real_to_complex_nd(float *real_data, float *complex_data, int nx, int ny, int nz, int batch);
int cuda_dd_fft_complex_to_real_nd(float *complex_data, float *real_data, int nx, int ny, int nz, int batch);
void do_cuda_fft_cache_destroy();

#endif //eman__cuda_emfft_h__

