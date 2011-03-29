
#ifndef cuda_processor_h__
#define cuda_processor_h__ 1

#include "cuda_util.h"
	
void emdata_processor_mult(float* data, const float& mult, const int nx, const int ny, const int nz);

void emdata_processor_add( EMDataForCuda* cuda_data, const float& sub);

void emdata_processor_to_value( EMDataForCuda* cuda_data, const float& value);

void emdata_processor_correlation_texture( const EMDataForCuda* left,const int center);

void emdata_processor_correlation( const EMDataForCuda* left, const EMDataForCuda* right,const int center);

void emdata_unwrap(float* data, int r1, int r2, int xs, int num_pi, int dx, int dy, int weight_radial, int nx, int ny);

float* emdata_transform_cuda(const float* const m, const int nx, const int ny, const int nz);

void calc_ccf_cuda(float* afft, const float* bfft, const int nx, const int ny, const int nz);

void calc_conv_cuda(float* afft, const float* bfft, const int nx, const int ny, const int nz);

CudaPeakInfo* calc_max_location_wrap_cuda(const float* in, const int nx, const int ny, const int nz, const int maxdx, const int maxdy, const int maxdz);

CudaPeakInfoFloat* calc_max_location_wrap_intp_cuda(const float* in, const int nx, const int ny, const int nz, const int maxdx, const int maxdy, const int maxdz);

void emdata_phaseorigin_to_center(float* data, const int nx, const int ny, const int nz);

void emdata_phaseorigin_to_corner(float* data, const int nx, const int ny, const int nz) ;

void emdata_ri2ap( EMDataForCuda* cuda_data);

void emdata_ap2ri( EMDataForCuda* cuda_data);

void emdata_ri2inten( EMDataForCuda* cuda_data);

void binarize_fourier_amp_processor(EMDataForCuda* cuda_data,const float& threshold);

void mult_complex_efficient_cuda(float* data, const float* src_data, const int nx, const int ny, const int nz, const int radius);

void mcf_cuda(const float* data1, float* data2, const int nx, const int ny, const int nz);

void subtract_cuda(float* data, float f, const int nx, const int ny, const int nz);

float* emdata_column_sum(const float* data, const int nx, const int ny); 

/** Rotates by 180 degrees using memory swapping, uses shared memory for efficiency
 * Works on 2D images - they can be odd in any dimension
 * @param cuda_data an EMDataForCuda struct - should have the data from a 2D image - doesn't check this, it's assumed that the calling function knows what it's doing
 * no return, processes the data inplace
 */
void emdata_rotate_180(float* data, const int nx, const int ny);
#endif // cuda_processor_h__
