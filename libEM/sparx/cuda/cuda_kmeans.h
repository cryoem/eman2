#ifndef sparx_cuda_kmeans_h__
#define sparx_cuda_kmeans_h__ 1

int cuda_kmeans(float* h_IM, float* h_AVE, unsigned short int* h_PART, float* h_INFO, int N, int m, int K, int maxite, float F, int rnd);

#endif // sparx_cuda_kmeans_h__ 1

