// Includes, system
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cublas.h>
#include <cuda.h>
#include <cuda_runtime_api.h>

// prototype function
//int rnd_asg(unsigned short int*, unsigned int*, int, int);
//void criterion_part(float*, unsigned short int*, unsigned int*, float*, float*, float*, float*, int, int, int);

// ERROR system
#define EXIT_OK (0)
#define ERROR_HOST_MEM (1) 
#define ERROR_DEVICE_MEM (2)
#define ERROR_DEVICE (3)
#define ERROR_INIT (4)
#define ERROR_EMPTY (5)
#define ERROR_SETDEVICE (6)
#define EXIT_DONE (255)

// kernel to calculate the exp
__global__ void kmeans_exp_kernel(float* DIST, float pw)
{
    register int idx = blockIdx.x * blockDim.x + threadIdx.x;
    register float arg = DIST[idx] * pw;
    if (arg < -70) arg = -70;
    DIST[idx] = exp(arg);
}

int cuda_inittest(int numdev) {
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    printf("NB device %d\n", deviceCount);
    return cudaSetDevice(numdev);

}

int cuda_readinit() {
    int deviceid;
    cudaGetDevice(&deviceid);
    return deviceid;
}

////// Init mem device //////////////////
int cuda_mpi_init(float* h_IM, float** hd_IM, float** hd_AVE, float** hd_DIST, int size_IM, int size_AVE, int size_DIST, int rnd, int numdev) {
    int status;
    float* d_AVE = NULL;
    float* d_IM = NULL;
    float* d_DIST = NULL;

    // Select a device
    if (numdev != -1) {
	if (EXIT_OK != cudaSetDevice(numdev)) return ERROR_SETDEVICE;
    }

    // Initialize CUBLAS
    srand(rnd);
    if (CUBLAS_STATUS_SUCCESS != cublasInit()) return ERROR_DEVICE;
    // Allocate device memory for images
    status = cublasAlloc(size_IM, sizeof(float), (void**)&d_IM);
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE_MEM;
    // Allocate device memory for averages
    status = cublasAlloc(size_AVE, sizeof(float), (void**)&d_AVE);
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE_MEM;
    // Allocate device memory for distances
    status = cublasAlloc(size_DIST, sizeof(float), (void**)&d_DIST);
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE_MEM;
    // Load images to the device
    status = cublasSetVector(size_IM, sizeof(float), h_IM, 1, d_IM, 1);
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE;

    //printf("\tthe value of d_AVE is %p\n", d_AVE);
    //printf("\tthe value of d_IM is %p\n", d_IM);
    //printf("\tthe value of d_DIST is %p\n", d_DIST);
    *hd_AVE = d_AVE;
    *hd_IM = d_IM;
    *hd_DIST = d_DIST;
    
    return EXIT_OK;

}

// k-means one step
int cuda_mpi_kmeans(float* h_AVE, float* d_AVE, float* h_DIST, float* d_DIST, float* d_IM, float* h_IM2, float* h_AVE2, unsigned short int* h_ASG, unsigned int* h_NC, int* params) {
    //           0, 1, 2,         3,        4
    // params = [N, m, K, flag_stop, ct_im_mv]
    int N = params[0];
    int m = params[1];
    int K = params[2];

     // status vars
    cublasStatus status;
    
    // prog vars GPU
    float alpha = 1.0f;
    float beta = 0.0f;
    int size_AVE = m * K;
    int size_DIST = N * K;

    // re-usable vars
    register int i, j, ind, c;
    float vmin = 0.0f;

    ////// GPU CODE (DISTANCES) ///////////////////////////

    // load averages to the device memory
    status = cublasSetVector(size_AVE, sizeof(float), h_AVE, 1, d_AVE, 1);
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE;
    // compute the distances
    cublasGetError();
    cublasSgemm('t', 'n', K, N, m, alpha, d_AVE, m, d_IM, m, beta, d_DIST, K);
    status = cublasGetError();
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE;
    // read the results from device to host memory
    status = cublasGetVector(size_DIST, sizeof(float), d_DIST, 1, h_DIST, 1);
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE;
    
    ////// K-MEANS //////////////////////////////////////
    params[3] = 1;                        // ask to stop k-means
    params[4] = 0;                        // set the re-assignment counter to zero
	
    // compute distances
    for (i = 0; i < N; i++) {
	vmin = 1.0e10;
	c = 0;
	ind = -1;
	for (j = i * K; j < (i + 1) * K; j++) {
	    h_DIST[j] = h_IM2[i] + h_AVE2[c] - 2 * h_DIST[j];  // compute distances
	    if (h_DIST[j] < vmin) {
		vmin = h_DIST[j];                              // select the smaller distance
		ind = c;
	    }
	    c++;
	}
	
	if (ind != h_ASG[i]) {
	    h_NC[h_ASG[i]]--;                                  // remove one object of the current class
	    h_ASG[i] = ind;                                    // assign the new object
	    h_NC[ind]++;                                       // add one object of the new class
	    params[3] = 0;                                     // if one object moves, ask to continue k-means
	    params[4]++; // ct_im_mv
	}
    }

    ////// NB CLUSTER CONTROL ///////////////////////////////////////

    // check the number of objects per class
    for (i = 0; i < K; i++) {
	if (h_NC[i] < 1) {
	    return ERROR_EMPTY;
	}
    }

    // check K-means is done
    if (params[3] == 1) {return EXIT_DONE;}

    return EXIT_OK;
}

// k-means SA one step
int cuda_mpi_kmeans_SA(float* h_AVE, float* d_AVE, float* h_DIST, float* d_DIST, float* d_IM, float* h_IM2, float* h_AVE2, unsigned short int* h_ASG, unsigned int* h_NC, float T, int* params) {
    //           0  1  2          3         4           5   6          7
    // params = [N, m, K, flag_stop, ct_im_mv, BLOCK_SIZE, NB, ins_BLOCK]
    int N = params[0];
    int m = params[1];
    int K = params[2];
    int BLOCK_SIZE = params[5];
    int NB = params[6];
    int ins_BLOCK = params[7];

    // status vars
    cublasStatus status;
    
    // prog vars GPU
    int size_DIST = N * K;
    int size_AVE = m * K;
    float alpha = 1.0f;
    float beta = 0.0f;

    // re-usable vars
    register int i, j, c;
    float buf = 0.0f;
    float vmax = 0.0f;
    float vmin = 0.0f;
    //float Je = 0.0f;

    ////// GPU CODE (DISTANCES) ///////////////////////////

    // load averages to the device memory
    status = cublasSetVector(size_AVE, sizeof(float), h_AVE, 1, d_AVE, 1);
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE;

    // compute the distances
    cublasGetError();
    cublasSgemm('t', 'n', K, N, m, alpha, d_AVE, m, d_IM, m, beta, d_DIST, K);
    status = cublasGetError();
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE;

    // read the results from device to host memory
    status = cublasGetVector(size_DIST, sizeof(float), d_DIST, 1, h_DIST, 1);
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE;

    // compute Je
    //Je = 0.0f;
    //for (i = 0; i < N; i++) Je += ((h_IM2[i] + h_AVE2[h_ASG[i]] - 2 * h_DIST[i * K + h_ASG[i]]) / (float)m);
	   
    ////// SIMULATED ANNEALING /////////////////////////////////

    // compute dJe normalize
    for (i = 0; i < N; i++) {
	// Ji
	buf = (float)h_NC[h_ASG[i]] / (float)(h_NC[h_ASG[i]] - 1) * (h_IM2[i] + h_AVE2[h_ASG[i]] - 2 * h_DIST[i * K + h_ASG[i]]);
	vmax = -1.0e6;
	vmin = 1.0e6;
	c = 0;
	for (j = i * K; j < (i + 1) * K; j++) {
	    if (c != h_ASG[i]) {
		// dJe = Ji - Jj
		h_DIST[j] = buf - (float)h_NC[c] / (float)(h_NC[c] + 1) * (h_IM2[i] + h_AVE2[c] - 2 * h_DIST[j]);
	    } else {
		// Ji = Jj so dJe = 0.0
		h_DIST[j] = 0.0f; 
	    }
	    if (h_DIST[j] > vmax) vmax = h_DIST[j];
	    if (h_DIST[j] < vmin) vmin = h_DIST[j];
	    c++;
	}
	// Normalize dJe [0, 1]
	for (j = i * K; j < (i + 1) * K; j++) h_DIST[j] = 1.0 - (h_DIST[j] - vmin) / (vmax - vmin);     
    }
		
    ////// GPU CODE (PROB) ///////////////////////////////////////

    // load dJe norm (DIST) to the device memory
    status = cublasSetVector(size_DIST, sizeof(float), h_DIST, 1, d_DIST, 1);
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE;

    // prepare kernel
    dim3 threads(BLOCK_SIZE);
    dim3 grid(NB);

    // execute the kernel probability for each distances
    kmeans_exp_kernel<<< grid, threads >>>(d_DIST, (-1.0 / T));

    // read the results from device to host memory
    status = cublasGetVector(size_DIST, sizeof(float), d_DIST, 1, h_DIST, 1);
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE;

    // if not enough block, complete the probability with c-code
    if (ins_BLOCK != size_DIST) {
	for (i = ins_BLOCK; i < size_DIST; i++) {
	    buf = h_DIST[i] * (-1.0 / T);
	    if (buf < -70) buf = -70;
	    h_DIST[i] = exp(buf);
	}
    }

    //////////////////////////////////////////////////////////////

    // select new assignment
    params[3] = 1;                                               // ask to stop k-means
    params[4] = 0;                                               // set the re-assignment counter to zero
    for (i = 0; i < N; i++) {
	c = 0;
	buf = 0;
	for (j = i * K; j < (i + 1) * K; j++) buf += h_DIST[j];  // sum of p along all classes for the given image i
	for (j = i * K; j < (i + 1) * K; j++) {
	    h_DIST[j] /= buf;                                    // norm of p (sum = 1)
	    if (c > 0) h_DIST[j] += h_DIST[j - 1];               // pk = pk + p(k-1)
	    c++;
	}
	buf = rand() / (float)RAND_MAX;
	j = 0;
	while (h_DIST[i * K + j] < buf && j < (K- 1)) j++;       // select new class, add j < K to resolve an issue due to
	if (j != h_ASG[i]) {                                     // the comparaison between approximation of buf and h_DIST,
	    h_NC[h_ASG[i]]--;                                    // decreasse the number of objects to the current class
	    h_ASG[i] = j;                                        // assign object in a new class
	    h_NC[h_ASG[i]]++;                                    // increase the number of objects to the new class
	    params[3] = 0;                                       // if one object moves, ask to continue k-means
	    params[4]++;                                         // count the number of re-assignments
	}
    }

    ////// NB CLUSTER CONTROL ///////////////////////////////////////

    // check the number of objects per class
    for (i = 0; i < K; i++) {
	if (h_NC[i] < 1) {
	    return ERROR_EMPTY;
	}
    }

    return EXIT_OK;
}

// Compute ji
int cuda_mpi_dist(float *h_AVE, float* d_AVE, float* h_DIST, float* d_DIST, float* d_IM, int N, int K, int m) { 
    // status vars
    cublasStatus status;

    // prog vars GPU
    int size_DIST = N * K;
    int size_AVE = m * K;
    float alpha = 1.0f;
    float beta = 0.0f;

    // load averages to the device memory
    status = cublasSetVector(size_AVE, sizeof(float), h_AVE, 1, d_AVE, 1);
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE;

    // compute the distances
    cublasGetError();
    cublasSgemm('t', 'n', K, N, m, alpha, d_AVE, m, d_IM, m, beta, d_DIST, K);
    status = cublasGetError();
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE;

    // read the results from device to host memory
    status = cublasGetVector(size_DIST, sizeof(float), d_DIST, 1, h_DIST, 1);
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE;

    return EXIT_OK;
}





// Shutdown and release memory
int cuda_mpi_shutdown(float* d_IM, float* d_AVE, float* d_DIST) {
    int status;
    // device memory
    status = cublasFree(d_IM);
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE_MEM;
    status = cublasFree(d_AVE);
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE_MEM;
    status = cublasFree(d_DIST);
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE_MEM;

    status = cublasShutdown();
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE;
    
    return EXIT_OK;
}


