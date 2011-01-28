
#ifndef cuda_processor_h__
#define cuda_processor_h__ 1


void calc_phase_weights_cuda(EMDataForCuda* data,float np);

void mean_phase_error_cuda(EMDataForCuda* left, EMDataForCuda* right,EMDataForCuda* wt, EMDataForCuda* hist, EMDataForCuda* norm,int num_hist);
		
#endif //cuda_processor_h__