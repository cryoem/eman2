#ifndef sparx_cuda_kmeans_h__
#define sparx_cuda_kmeans_h__ 1

int cuda_kmeans(float* h_IM, float* h_AVE, unsigned short int* h_ASG, float* h_INFO, int N, int m, int K, int maxite, float F, float T0, int rnd);

#endif // sparx_cuda_kmeans_h__ 1

