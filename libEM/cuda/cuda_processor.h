
#ifndef cuda_processor_h__
#define cuda_processor_h__ 1

#include "cuda_util.h"

void emdata_processor_mult( EMDataForCuda* cuda_data, const float& mult);

void emdata_processor_add( EMDataForCuda* cuda_data, const float& sub);

void emdata_processor_to_value( EMDataForCuda* cuda_data, const float& value);

void emdata_processor_correlation_texture( const EMDataForCuda* left);

void emdata_processor_correlation( const EMDataForCuda* left, const EMDataForCuda* right);

EMDataForCuda* emdata_unwrap( int r1, int r2, int xs, int do360, int nx, int ny);

void emdata_phaseorigin_to_center(EMDataForCuda* cuda_data);

#endif // cuda_processor_h__