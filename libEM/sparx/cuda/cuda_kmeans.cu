// Includes, system
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
// includes, project
//#include <cutil_inline.h>
// includes, cublas
#include <cublas.h>

// prototype function
int rnd_asg(unsigned short int*, unsigned int*, int, int);
void criterion_part(float*, unsigned short int*, unsigned int*, float*, float*, float*, float*, int, int, int);

/* DATA STRUCTURE OF INFO
 * SEQ INFO = [Je, Coleman, Davies-Bouldin, Harabasz, number of iterations, running time]
 * 
 */
#define JE_INFO (0)
#define C_INFO (1)
#define DB_INFO (2)
#define H_INFO (3)
#define NOI_INFO (4)
#define TIME_INFO (5)

// ERROR system
#define EXIT_OK (0)
#define ERROR_HOST_MEM (1) 
#define ERROR_DEVICE_MEM (2)
#define ERROR_DEVICE (3)
#define ERROR_INIT (4)
#define ERROR_EMPTY (5)

// number of iterations to calcualte the average of the number of re-assignments
#define SAMPLE_SA (5)
// number of elements in the schedule first temperature
#define N_SCHEDULE (37)

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

// Main, k-means core
int cuda_kmeans(float* h_IM, float* h_AVE, unsigned short int* h_ASG, float* h_INFO, int N, int m, int K, int maxite, float F, float T0, int rnd) {    
    // host vars
    float* h_DIST;
    float* h_AVE2;
    float* h_IM2;
    unsigned int* h_NC;

    // device vars
    float* d_IM;
    float* d_AVE;
    float* d_DIST;

    // status vars
    cublasStatus status;

    // timing vars
    time_t mt1, mt2;
    
    // prog vars GPU
    int size_IM = m * N;
    int size_AVE = m * K;
    int size_DIST = N * K;
    const int BLOCK_SIZE = 512;
    int NB = N * K / BLOCK_SIZE;
    int ins_BLOCK = NB * BLOCK_SIZE;
    float alpha = 1.0f;
    float beta = 0.0f;

    // prog vars flow control
    int ite;
    int flag_stop_part;
    int flag_stop_SA;
    int count_im_moved;

    // SA
    float T;
    float ref_SA = (float) N * 0.97;
    float inc_SA [N_SCHEDULE] = {100, 90, 80, 70, 60, 50, 40, 30, 20, 10, \
				 9, 8, 7, 6, 5, 4, 3, 2, 1,		\
				 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, \
				 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01};
    float val_SA;
    int flag_servoing_SA;
    int c_SA;
        
    // re-usable vars
    register int i, j, ind, c, d;
    float buf = 0.0f;
    float vmax = 0.0f;
    float vmin = 0.0f;

    srand(rnd);
    //float Je = 0.0f;
    
    // Initialize CUBLAS
    //printf("simpleCUBLAS test running..\n");

    status = cublasInit();
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE;

    ////// ALLOCATION HOST MEMORY ///////////////////////////

    // allocate the memory for the sum squared of averages
    h_AVE2 = (float*)malloc(K * sizeof(float));
    if (h_AVE2 == 0) return ERROR_HOST_MEM;
    // allocate the memory for the sum squared of images
    h_IM2 = (float*)malloc(N * sizeof(float));
    if (h_IM2 == 0) return ERROR_HOST_MEM;
    // allocate the memory for the distances
    h_DIST = (float*)malloc(size_DIST * sizeof(float));
    if (h_DIST == 0) return ERROR_HOST_MEM;
    // allocate the memory for the number of images per class
    h_NC = (unsigned int*)malloc(K * sizeof(unsigned int));
    if (h_NC == 0) return ERROR_HOST_MEM;

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

    ////// PRECALCULATION SUM SQUARED IMAGES /////////////////////
    
    for (i = 0; i < N; i++) {
	h_IM2[i] = 0.0f;
	for (j = 0; j < m; j++) h_IM2[i] += (h_IM[i * m + j] * h_IM[i * m + j]);
    }

    // set the first assignment randomly
    if (rnd_asg(h_ASG, h_NC, K, N) != EXIT_SUCCESS) return ERROR_INIT;

    // set flag to run the partition
    flag_stop_part = 0;

    // init simulated annealing
    if (F != 0.0) {
	flag_stop_SA = 0;     // set flag to run SA
	if (T0 == -1) {
	    flag_servoing_SA = 1; // ask to auto determine the first temperature 
	    T = inc_SA[0];        // assign the temperature T = 100 (first in the schedule)
	} else {
	    flag_servoing_SA = 0;
	    T = T0;
	}
    } else {
	flag_stop_SA = 1;     // set flag to run only k-means (without SA)
	flag_servoing_SA = 0; // disable
    }
    val_SA = 0;           // reset the counter for random assignment
    c_SA = 0;             // init the temperature schedule index

    // start main timer
    mt1 = time(NULL);

    ////////////////////////////////////////////////////////////
    ////// HERE IS THE LOOP OVER ITERATION /////////////////////
    ////////////////////////////////////////////////////////////
    ite = 0;
    while (ite < maxite && !flag_stop_part) {
	ite++;
	    
	// compute the averages according ASG
	for (i = 0; i < size_AVE; ++i) h_AVE[i] = 0.0f;                          // set averages to zero
	for (i = 0; i < N; ++i) {
	    c = h_ASG[i] * m;
	    d = i * m;
	    for (j = 0; j < m; ++j) h_AVE[c + j] += h_IM[d + j];                 // accumulate images
	}
	for (i = 0; i < K; i++) {
	    buf = 0.0f;
	    for (j= 0 ; j < m; j++) {
		ind = i * m + j;
		h_AVE[ind] /= (float)h_NC[i];                                    // compute average 
		buf += (h_AVE[ind] * h_AVE[ind]);                                // do sum squared AVE 
	    }
	    h_AVE2[i] = buf;
	}

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
	   
	////////////////////////////////////////////////////////////
	////// SIMULATED ANNEALING /////////////////////////////////
	////////////////////////////////////////////////////////////

	if (flag_stop_SA == 0) {

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
	    flag_stop_part = 1;                                                 // ask to stop k-means
	    count_im_moved = 0;                                                 // set the re-assignment counter to zero
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

	    ////// FLOW CONTROL //////////////////////////////////////////

	    // SA determine the first temperature automatically
	    if (flag_servoing_SA) {
		val_SA += count_im_moved;              // accumulate the number of re-assignments for 5 iterations
		if (ite % SAMPLE_SA == 0) {            // for every 5 iterations
		    val_SA /= (float)SAMPLE_SA;        // do the average
		    if (val_SA > ref_SA) {             // if the value is still too high
			c_SA++;                        // inc the index in the schedule
			if (c_SA < N_SCHEDULE) {       // if the index is not ouside the schedule (37)
			    T = inc_SA[c_SA];          // assign new T from the schedule
			    val_SA = 0;                // set to zero the accumulator
			} else {                       // if the index is outside the schedule
			    T = inc_SA[N_SCHEDULE - 1];// set T = 0.01, with the lower value from the schedule (36)
			    flag_servoing_SA = 0;      // start SA exploitation
			}
		    } else {                           // if the value is lower than 97 % of N
			if (c_SA != 0) {               // if the index schedule is not the first temperature
			    T = inc_SA[c_SA - 1];      // select the previous one temperature
			} else {
			    T = inc_SA[0];             // if the index schedule is the first temperature T = 100
			}         
			flag_servoing_SA = 0;          // start SA exploitation
		    }
		}

            // SA exploitation
	    } else {
		T *= F;                                // cooling temperature
		if (T < 0.0001) flag_stop_SA = 1;      // turn off SA, we consider this T value is equivalent to k-means
	    }

        ////////////////////////////////////////////////////////////
        ////// K-MEANS /////////////////////////////////////////////
        ////////////////////////////////////////////////////////////
	} else {
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
		    count_im_moved++;                                  // count the number of re-assignments
		}
	    }

	} // k-means

        ////// FLOW CONTROL //////////////////////////////////////////

	// check the number of objects per class
	for (i = 0; i < K; i++) {
	    if (h_NC[i] < 1) {
		// clean host memory
		free(h_DIST);
		free(h_AVE2);
		free(h_IM2);
		free(h_NC);
		// clean device memory
		status = cublasFree(d_IM);
		if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE_MEM;
		status = cublasFree(d_AVE);
		if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE_MEM;
		status = cublasFree(d_DIST);
		if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE_MEM;
		// shutdown
		status = cublasShutdown();
		if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE;
		
		return ERROR_EMPTY;
	    }
	}
	
	    
	//printf("ite: %i T: %f CRE: %i         Je: %e   ", ite, T, count_im_moved, Je);
	//printf("\n");
	//print_int(h_NC, K);

    //////////////////////////////////////////////////////////////////
    } // end of loop over iteration
    //////////////////////////////////////////////////////////////////
        
	
    // save running time of the current partition							
    mt2 = time(NULL);
    //printf("time: %f\n", (float)(mt2 - mt1));
    h_INFO[TIME_INFO] = (float)(mt2 - mt1);

    // save number of iterations
    h_INFO[NOI_INFO] = (float)ite;

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

    ////// MEMORY CLEAN UP ///////////////////////////////////////////

    // device memory
    status = cublasFree(d_IM);
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE_MEM;
    status = cublasFree(d_AVE);
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE_MEM;
    status = cublasFree(d_DIST);
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE_MEM;

    // host memory
    free(h_DIST);
    free(h_AVE2);
    free(h_IM2);
    free(h_NC);

    ////// SHUTDOWN ///////////////////////////////////////////
    
    status = cublasShutdown();
    if (status != CUBLAS_STATUS_SUCCESS) return ERROR_DEVICE;
    
    return EXIT_OK;

}

////// SOME FUNCTIONS /////////////////////////////////////////////

// selected randomize assignment and count the number of images per class
int rnd_asg(unsigned short int* ASG, unsigned int* NC, int K, int N) {
    int ret = 20;
    int flag = 0;
    int n, k;
    for (k = 0; k < K; k++) NC[k] = 0;
    while (ret > 0) {
	ret--;
	for (n = 0; n < N; n++) {
	    ASG[n] = rand() % K;
	    NC[ASG[n]]++;
	}
	flag = 1;
	k = K;
	while (k > 0 && flag) {
	    k--;
	    if (NC[k] <= 1) {
		flag = 0;
		if (ret == 0) {
		    //printf("Erreur randomize assignment\n");
		    return EXIT_FAILURE;
		}
		for (k = 0; k < K; k++) NC[k] = 0;
	    }
	if (flag == 1) ret = 0;
	}
    }
    return EXIT_SUCCESS;
}

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
