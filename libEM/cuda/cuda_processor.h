
#ifndef cuda_processor_h__
#define cuda_processor_h__ 1

#include "cuda_util.h"

void emdata_processor_mult( EMDataForCuda* cuda_data, const float& mult);

void emdata_processor_add( EMDataForCuda* cuda_data, const float& sub);

void emdata_processor_to_value( EMDataForCuda* cuda_data, const float& value);

void emdata_processor_correlation_texture( const EMDataForCuda* left,const int center);

void emdata_processor_correlation( const EMDataForCuda* left, const EMDataForCuda* right,const int center);

//EMDataForCuda* emdata_unwrap( int r1, int r2, int xs, int do360, int nx, int ny);
void emdata_unwrap(EMDataForCuda* data, int r1, int r2, int xs, int num_pi, int dx, int dy, int weight_radial, int nx, int ny);

void emdata_phaseorigin_to_center(EMDataForCuda* cuda_data);

EMDataForCuda* emdata_transform_cuda(const float* const m,const int nx, const int ny, const int nz);

void emdata_ri2ap( EMDataForCuda* cuda_data);

void emdata_ap2ri( EMDataForCuda* cuda_data);

void emdata_ri2inten( EMDataForCuda* cuda_data);

void binarize_fourier_amp_processor(EMDataForCuda* cuda_data,const float& threshold);

/** Rotates by 180 degrees using memory swapping, uses shared memory for efficiency
 * Works on 2D images - they can be odd in any dimension
 * @param cuda_data an EMDataForCuda struct - should have the data from a 2D image - doesn't check this, it's assumed that the calling function knows what it's doing
 * no return, processes the data inplace
 */
void emdata_rotate_180( EMDataForCuda* cuda_data);
#endif // cuda_processor_h__
