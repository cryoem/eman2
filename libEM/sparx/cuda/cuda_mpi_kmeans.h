#ifndef cuda_mpi_kmeans_h__
#define cuda_mpi_kmeans_h__ 1

int cuda_inittest(int numdev);
int cuda_readinit();
int cuda_mpi_init(float* h_IM, float** hd_IM, float** hd_AVE, float** hd_DIST, int size_IM, int size_AVE, int size_DIST, int rnd, int numdev);
int cuda_mpi_kmeans(float* h_AVE, float* d_AVE, float* h_DIST, float* d_DIST, float* d_IM, float* h_IM2, float* h_AVE2, unsigned short int* h_ASG, unsigned int* h_NC, int* params);
int cuda_mpi_shutdown(float* d_IM, float* d_AVE, float* d_DIST);
int cuda_mpi_kmeans_SA(float* h_AVE, float* d_AVE, float* h_DIST, float* d_DIST, float* d_IM, float* h_IM2, float* h_AVE2, unsigned short int* h_ASG, unsigned int* h_NC, float T0, int* params);

#endif // cuda_mpi_kmeans_h__ 1

