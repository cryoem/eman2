
#ifndef cuda_cmp_h__
#define cuda_cmp_h__ 1

float ccc_cmp_cuda(const float* data1, const float* data2, const float* dm, const int &nx, const int &ny, const int &nz);

float dot_cmp_cuda(float* data1, float* data2, const float* dm, const int &nx, const int &ny, const int &nz);

float2 get_stats_cuda(const float * data, const int nx, const int ny, const int nz);

void normalize_cuda(float * data, float mean, float var, const int nx, const int ny, const int nz);

float get_value_at_wrap_cuda(float * data, int tx, int ty, int tz, int nx, int ny, int nz);

float* calc_fourier_shell_correlation_cuda(const int nx, const int ny, const int nz, const int d);

float fsc_tomo_cmp_cuda(const float* data1, const float* data2, const float data1threshold, const float data2threshold, const float minres, const float maxres, const int nx, const int ny, const int nz);

#endif //cuda_processor_h__