#ifndef cuda_reconstruct_h__
#define cuda_reconstruct_h__ 1

void insert_slice_cuda(const float* const matrix, const float* const slice_data, float* vol, float* tmp_data, const int inx, const int iny, const int nx, const int ny, const int nz, const float weight);

float4 determine_slice_agreement_cuda(const float* const matrix, const float* const slice_data, float* vol, float* tmp_data, const int inx, const int iny, const int nx, const int ny, const int nz, const float weight);

#endif //  cuda_project_h__ 1