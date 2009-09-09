// Includes, system
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cublas.h>

// prototype function
int rnd_asg(unsigned short int*, unsigned int*, int, int);
void criterion_part(float*, unsigned short int*, unsigned int*, float*, float*, float*, float*, int, int, int);

// ERROR system
#define EXIT_OK (0)
#define ERROR_HOST_MEM (1) 
#define ERROR_DEVICE_MEM (2)
#define ERROR_DEVICE (3)
#define ERROR_INIT (4)
#define ERROR_EMPTY (5)

// kernel to calculate the exp
__global__ void exp_kernel(float* DIST, float pw)
{
    register int idx = blockIdx.x * blockDim.x + threadIdx.x;
    register float arg = DIST[idx] * pw;
    if (arg < -70) arg = -70;
    DIST[idx] = exp(arg);
}

void print_int(unsigned int* data, int size) {
    for (int i = 0; i < size; i++) printf(" %u", data[i]);
    printf("\n");
}

void print_float(float* data, int size) {
    for (int i = 0; i < size; i++) printf(" %f", data[i]);
    printf("\n");
}

// Init mem device
int cuda_mpi_init(float* h_IM, float* d_IM, float* d_AVE, float* d_DIST, int size_IM, int size_AVE, int size_DIST, int rnd) {

    ////// Initialize CUBLAS /////////////////////////
    //printf("simpleCUBLAS test running..\n");
    srand(rnd);
    status = cublasInit();
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE;

    ////// ALLOCATION DEVICE MEMORY ///////////////////////////

    // allocate device memory for images
    status = cublasAlloc(size_IM, sizeof(float), (void**)&d_IM);
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE_MEM;
    // allocate device memory for averages
    status = cublasAlloc(size_AVE, sizeof(float), (void**)&d_AVE);
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE_MEM;
    // allocate device memory for distances
    status = cublasAlloc(size_DIST, sizeof(float), (void**)&d_DIST);
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE_MEM;

    ////// LOADING IMAGES TO THE DEVICE MEMORY /////////////////
    
    status = cublasSetVector(size_IM, sizeof(float), h_IM, 1, d_IM, 1);
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE;

    return EXIT_OK;

}

// k-means one step
int cuda_mpi_kmeans(float* h_AVE, float* d_AVE, float* h_DIST, float* d_DIST, float* d_IM, float* h_IM2, float* h_AVE2, unsigned short int* h_ASG, unsigned int* h_NC, int* flag_stop_part, int K, int N, int m) {
     // status vars
    cublasStatus status;
    
    // prog vars GPU
    float alpha = 1.0f;
    float beta = 0.0f;
    int size_AVE = m * K;
    int isze_DIST = n * K;

    // re-usable vars
    register int i, j, ind, c;
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

    ////// K-MEANS //////////////////////////////////////
    flag_stop_part = 1;                        // ask to stop k-means
    count_im_moved = 0;                        // set the re-assignment counter to zero
	
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
	    flag_stop_part = 0;                                // if one object moves, ask to continue k-means
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


// k-means SA one step
int cuda_mpi_kmeans_SA(float* h_AVE, float* d_AVE, float* h_DIST, float* d_DIST, float* d_IM, float* h_IM2, float* h_AVE2, unsigned short int* h_ASG, unsigned int* h_NC, int* flag_stop_part, int* ct_im_mv, int K, int N, int m, float T0, int BLOCK_SIZE, int NB, int ins_BLOCK) {
     // status vars
    cublasStatus status;
    
    // prog vars GPU
    int size_DIST = N * K;
    int size_AVE = m * K;
    float alpha = 1.0f;
    float beta = 0.0f;

    // re-usable vars
    register int i, j, ind, c;
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
    exp_kernel<<< grid, threads >>>(d_DIST, (-1.0 / T));

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
    flag_stop_part = 1;                                           // ask to stop k-means
    ct_im_mv = 0;                                                 // set the re-assignment counter to zero
    for (i = 0; i < N; i++) {
	c = 0;
	buf = 0;
	for (j = i * K; j < (i + 1) * K; j++) buf += h_DIST[j];         // sum of p along all classes for the given image i
	for (j = i * K; j < (i + 1) * K; j++) {
	    h_DIST[j] /= buf;                                           // norm of p (sum = 1)
	    if (c > 0) h_DIST[j] += h_DIST[j - 1];                      // pk = pk + p(k-1)
	    c++;
	}
	buf = rand() / (float)RAND_MAX;
	j = 0;
	while (h_DIST[i * K + j] < buf && j < (K- 1)) j++;                   // select new class, add j < K to resolve an issue due to
	if (j != h_ASG[i]) {                                            // the comparaison between approximation of buf and h_DIST,
	    h_NC[h_ASG[i]]--;                                           // decreasse the number of objects to the current class
	    h_ASG[i] = j;                                               // assign object in a new class
	    h_NC[h_ASG[i]]++;                                           // increase the number of objects to the new class
	    flag_stop_part = 0;                                         // if one object moves, ask to continue k-means
	    count_im_moved++;                                           // count the number of re-assignments
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

int cuda_mpi_shutdown(float* d_IM, float* d_AVE, float* d_DIST) {
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

    /*

    ////// COMPUTE CRITERION  /////////////////////////////////////////
    // compute the distances again with GPU
    cublasGetError();
    cublasSgemm('t', 'n', K, N, m, alpha, d_AVE, m, d_IM, m, beta, d_DIST, K);
    status = cublasGetError();
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE;

    // read the results from device to host memory
    status = cublasGetVector(size_DIST, sizeof(float), d_DIST, 1, h_DIST, 1);
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE;

    // compute some criterion and store to h_INFO
    criterion_part(h_DIST, h_ASG, h_NC, h_AVE, h_IM2, h_AVE2, h_INFO, N, m, K);


}
    */
////// SOME FUNCTIONS /////////////////////////////////////////////

/* DATA STRUCTURE OF INFO
 * SEQ INFO = [Je, Coleman, Davies-Bouldin, Harabasz, number of iterations, running time]
 *
 * compute some criterions (Je, C, DB, H, and noi)
 */
void criterion_part(float* DIST, unsigned short int* ASG, unsigned int* NC, \
		    float* AVE, float* IM2, float* AVE2, float* INFO, \
		    int N, int m, int K) {
    float buf = 0.0f;
    float Je = 0.0f;
    float Tr_AVE = 0.0f;
    float v_max = 0.0f;
    float* Ji = (float*)calloc(K, sizeof(float));
    float* S_AVE2 = (float*)calloc(m, sizeof(float));
    float* S_AVE = (float*)calloc(m, sizeof(float));
    int i, j, k;

    // first compute Je
    for (i = 0; i < N; i++) Ji[ASG[i]] += (IM2[i] + AVE2[ASG[i]] - 2 * DIST[i * K + ASG[i]]);
    for (i = 0; i < K; i++) Je += (Ji[i] / (float)m);
    INFO[JE_INFO] = Je;
    
    // second compute Trace(ave)
    for (i =0; i < K; i++) {
	for (j = 0; j < m; j++) {
	    S_AVE[j] += AVE[i * m + j];
	    S_AVE2[j] += (AVE[i * m + j] * AVE[i * m + j]);
	}
    }
    buf = 1 / (float)K;
    for (i = 0; i < m; ++i) Tr_AVE += (buf * (S_AVE2[i] - buf * S_AVE[i] * S_AVE[i]));

    // compute Coleman's criteria (Tr_AVE * Je)
    INFO[C_INFO] = Tr_AVE * Je;

    // compute Harabasz's criteria (N-K)/(K-1) * Tr_AVE / Je
    //INFO[ipart * SEQ_INFO + H_INFO] = (Tr_AVE / (float)(K - 1)) / (Je / (float)(N - K));
    INFO[H_INFO] = (Tr_AVE * (float)(N - K)) / (Je * (float)(K - 1));
    //INFO[ipart * SEQ_INFO + H_INFO] = ((float)(N - K) / (float)(K - 1)) * (Tr_AVE / INFO[ipart * SEQ_INFO + JE_INFO]);

    // compute Davies-Bouldin's criteria
    for (i = 0; i < K; i++) {
	v_max = 0.0f;
	for (j = 0; j < K; j++) {
	    if (j != i) {
		buf = 0.0f;
		for (k = 0; k < m; k++) buf += ((AVE[j * m + k] - AVE[i * m + k]) * (AVE[j * m + k] - AVE[i * m + k]));
		buf = (Ji[i] / (float)NC[i] + Ji[j] / (float)NC[j]) * ((float)m / buf);
	    }
	    if (buf > v_max) v_max = buf;
	}
	INFO[DB_INFO] += v_max;
    }
    INFO[DB_INFO] /= K;

    // free memory
    free(Ji);
    free(S_AVE);
    free(S_AVE2);
}
