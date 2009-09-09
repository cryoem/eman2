#ifndef cuda_mpi_kmeans_h__
#define cuda_mpi_kmeans_h__ 1

int cuda_mpi_init(float* h_IM, float* d_IM, float* d_AVE, float* d_DIST, int size_IM, int size_AVE, int size_DIST, int rnd);

int cuda_mpi_kmeans(float* h_AVE, float* d_AVE, float* h_DIST, float* d_DIST, float* d_IM, float* h_IM2, float* h_AVE2, unsigned short int* h_ASG, unsigned int* h_NC, int* flag_stop_part, int K, int N, int m);

int cuda_mpi_kmeans_SA(float* h_AVE, float* d_AVE, float* h_DIST, float* d_DIST, float* d_IM, float* h_IM2, float* h_AVE2, unsigned short int* h_ASG, unsigned int* h_NC, int* flag_stop_part, int* ct_im_mv, int K, int N, int m, float T0, int BLOCK_SIZE, int NB, int ins_BLOCK);

int cuda_mpi_shutdown(float* d_IM, float* d_AVE, float* d_DIST);

#endif // cuda_mpi_kmeans_h__ 1

