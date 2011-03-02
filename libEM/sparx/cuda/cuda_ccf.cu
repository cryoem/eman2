/*
 * Copyright 1993-2007 NVIDIA Corporation.  All rights reserved.
 *
 * NOTICE TO USER:
 *
 * This source code is subject to NVIDIA ownership rights under U.S. and
 * international Copyright laws.  Users and possessors of this source code
 * are hereby granted a nonexclusive, royalty-free license to use this code
 * in individual and commercial software.
 *
 * NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE
 * CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR
 * IMPLIED WARRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.
 * IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL,
 * OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS,  WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION,  ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOURCE CODE.
 *
 * U.S. Government End Users.   This source code is a "commercial item" as
 * that term is defined at  48 C.F.R. 2.101 (OCT 1995), consisting  of
 * "commercial computer  software"  and "commercial computer software
 * documentation" as such terms are  used in 48 C.F.R. 12.212 (SEPT 1995)
 * and is provided to the U.S. Government only as a commercial end item.
 * Consistent with 48 C.F.R.12.212 and 48 C.F.R. 227.7202-1 through
 * 227.7202-4 (JUNE 1995), all U.S. Government End Users acquire the
 * source code with only those rights set forth herein.
 *
 * Any use of this source code in individual and commercial software must
 * include, in the user documentation and internal comments to the code,
 * the above Disclaimer and U.S. Government End Users Notice.
 */


/* Includes, system */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Includes, cuda, cufft */
#include <cufft.h>

/* Matrix size */
#define PI (3.14159265358979f)

__global__ void complex_mul(float *ccf, int BLOCK_SIZE, int NRING, int NIMAGE, int KX, int KY);
__global__ void resample_to_polar(float* image, float dx, float dy, int NX, int NY, int RING_LENGTH, int NRING, int offset);
__global__ void mul_ctf(float *image, int nx, int ny, float defocus, float cs, float voltage, float apix, float bfactor, float ampcont);
__global__ void add_img(float *image_padded, float *ave1, float *ave2, int nx, int ny, int nima);
__global__ void rotate_shift(float *image, int nx, int ny, int offset);
__global__ void normalize_rings(float *image, int RING_LENGTH, int NRING, int offset);

texture<float, 2, cudaReadModeElementType> tex;
texture<float, 1, cudaReadModeElementType> texim_subject;
texture<float, 1, cudaReadModeElementType> texim_ref;
texture<float, 1, cudaReadModeElementType> texim_points;
texture<float, 1, cudaReadModeElementType> texim_shifts;
texture<float, 1, cudaReadModeElementType> texim_ali_params;
texture<float, 1, cudaReadModeElementType> texim_polar;
 
/* Main */
void cudasetup(int id) {
	cudaSetDevice(id);
	cudaError_t cudaError;
	cudaError = cudaGetLastError();
	if (cudaError != cudaSuccess) {
		fprintf(stderr, "CUDA Runtime API Error reported: %s.   Devide id = %d\n", cudaGetErrorString(cudaError), id);
		exit(0);
	}
}

void calculate_ccf(float *subject_image, float *ref_image, float *ccf, int NIMAGE, int NX, int NY, int RING_LENGTH, int NRING, int OU, float STEP, int KX, int KY, float *sx, float *sy, int silent)
{    
	float *d_subject_image_polar, *d_ref_image_polar;
	float *d_ccf;
	float *points, *d_points;
	float *shifts_ref, *d_shifts_ref;
	float *shifts_subject, *d_shifts_subject;
	int i, j, k, index;
	int ccf_base_addr;
	float x, y, ang;
	size_t* offset = (size_t *)malloc(sizeof(size_t));

	int BLOCK_SIZE = RING_LENGTH/2+1;
	int POINTS_PER_IMAGE = (RING_LENGTH+2)*NRING;

	// For a texture reference bound to a two-dimensional CUDA array,
	// the maximum width is 2^16 (=65536) and the maximum height is 2^15 (=32768)
	int NIMAGE_ROW = 32768/NX;
	int NROW = NIMAGE/NIMAGE_ROW;
	int NIMAGE_LEFT = NIMAGE%NIMAGE_ROW;

	// For a texture reference bound to linear memory, the maximum width is 2^27
	int NIMAGE_IN_TEXTURE = (1<<27)/POINTS_PER_IMAGE*9/10;  
	int NTEXTURE = NIMAGE/NIMAGE_IN_TEXTURE;
	int NIMAGE_LEFT_TEXTURE = NIMAGE%NIMAGE_IN_TEXTURE;

	int IMAGE_PER_BATCH1 = 65535/NRING;            // The maximum FFT per batch is 65535
	int IMAGE_BATCH1 = NIMAGE/IMAGE_PER_BATCH1;
	int IMAGE_LEFT_BATCH1 = NIMAGE%IMAGE_PER_BATCH1;

	int POS_PER_IMAGE = 2*(2*KX+1)*(2*KY+1);
	int IMAGE_PER_BATCH2 = 65535/POS_PER_IMAGE; // The maximum FFT per batch is 65535
	int IMAGE_BATCH2 = NIMAGE/IMAGE_PER_BATCH2;
	int IMAGE_LEFT_BATCH2 = NIMAGE%IMAGE_PER_BATCH2;

	cudaArray *ref_image_array, *subject_image_array[NROW], *subject_image_array_left;
	dim3 GridSize1(NRING, NIMAGE_ROW);
	dim3 GridSize2(NRING, NIMAGE_LEFT);
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	const textureReference *texPtr;

	tex.normalized = false;
	tex.filterMode = cudaFilterModeLinear;
	tex.addressMode[0] = cudaAddressModeWrap;
	tex.addressMode[1] = cudaAddressModeWrap;

	cufftHandle plan_subject_image, plan_subject_image_rest, plan_ref_image, plan_ccf, plan_ccf_rest;
	cufftPlan1d(&plan_subject_image, RING_LENGTH, CUFFT_R2C, NRING*IMAGE_PER_BATCH1);
	if (IMAGE_LEFT_BATCH1 != 0)
		cufftPlan1d(&plan_subject_image_rest, RING_LENGTH, CUFFT_R2C, NRING*IMAGE_LEFT_BATCH1);
	cufftPlan1d(&plan_ref_image, RING_LENGTH, CUFFT_R2C, NRING);
	cufftPlan1d(&plan_ccf, RING_LENGTH, CUFFT_C2R, POS_PER_IMAGE*IMAGE_PER_BATCH2);
	if (IMAGE_LEFT_BATCH2 != 0)
		cufftPlan1d(&plan_ccf_rest, RING_LENGTH, CUFFT_C2R, POS_PER_IMAGE*IMAGE_LEFT_BATCH2);

	/* Allocate host memory for the coordinates of sampling points */
	points = (float *)malloc(RING_LENGTH*NRING*2*sizeof(float));
	if (points == 0) {
		fprintf (stderr, "Host memory allocation error!\n");
		return;
	}

	/* Allocate host memory for the coordinates of all shifts */
	shifts_ref = (float *)malloc(2*sizeof(float));
	if (shifts_ref == 0) {
		fprintf (stderr, "Host memory allocation error!\n");
		return;
	}

	shifts_subject = (float *)malloc(NIMAGE*2*sizeof(float));
	if (shifts_subject == 0) {
		fprintf(stderr, "Host memory allocation error!\n");
		return;
	}

	if (silent == 0) {
		printf("Initialization on the host memory done.\n");

		printf("\nMemory to be allocated on the video card:\n");
		printf("For %5d subject images                      : %10.3f MB\n", NIMAGE, NIMAGE*NX*NY*4.0/1000000.0);
		printf("For reference image			      : %10.3f KB\n", NX*NY*4/1000.0);
		printf("For %5d subject images in polar coordinates : %10.3f MB\n", NIMAGE, NIMAGE*POINTS_PER_IMAGE*4.0/1000000.0);
		printf("For reference image in polar coordinates      : %10.3f KB\n", (RING_LENGTH+2)*NRING*4/1000.0);
		printf("For all cross-correlation functions (CCF)     : %10.3f MB\n", (RING_LENGTH+2)*NIMAGE*POS_PER_IMAGE*4.0/1000000.0);
		printf("Total memory used			      : %10.3f MB\n\n", ((NIMAGE+1)*(NX*NY+POINTS_PER_IMAGE)+(RING_LENGTH+2)*NIMAGE*POS_PER_IMAGE)*4.0/1000000.0);
	}


	/* Allocate the matrix for all NIMAGE subject images on the video card */
	for (k=0; k<NROW; k++)
		cudaMallocArray(&subject_image_array[k], &channelDesc, NX, NY*NIMAGE_ROW);
	if (NIMAGE_LEFT != 0)
		cudaMallocArray(&subject_image_array_left, &channelDesc, NX, NY*NIMAGE_LEFT);

	/* Allocate the matrix for the reference image on the video card */
	cudaMallocArray(&ref_image_array, &channelDesc, NX, NY);

	/* Allocate the matrix for all NIMAGE subject images in polar coordinates on the video card */
	cudaMalloc((void**)&d_subject_image_polar, NIMAGE*POINTS_PER_IMAGE*sizeof(float));

	/* Allocate the matrix for the reference image in polar coordinates on the video card */
	cudaMalloc((void**)&d_ref_image_polar, POINTS_PER_IMAGE*sizeof(float));

	/* Allocate the matrix for the ccf on the video card */
	cudaMalloc((void**)&d_ccf, (RING_LENGTH+2)*NIMAGE*POS_PER_IMAGE*sizeof(float));

	/* Allocate the matrix for the coordinates of the sampling points on the video card */
	cudaMalloc((void**)&d_points, RING_LENGTH*NRING*2*sizeof(float));

	/* Allocate the matrix for the coordinates of the shifts of reference image on the video card */
	cudaMalloc((void**)&d_shifts_ref, 2*sizeof(float));

	/* Allocate the matrix for the coordinates of the shifts of subject images on the video card */
	cudaMalloc((void**)&d_shifts_subject, NIMAGE*2*sizeof(float));

	/* Fill the matrix for the coordinates of sampling points */
	for (i = 0; i < NRING; i++) 
		for (j = 0; j < RING_LENGTH; j++) {
			index = i*RING_LENGTH+j;
			ang = float(j)/RING_LENGTH*PI*2;
			points[index*2] = (i+1)*sinf(ang)*float(OU)/float(NRING);
			points[index*2+1] = (i+1)*cosf(ang)*float(OU)/float(NRING);
		}

	/* Fill the matrix for the coordinates of shifts for reference images (currently hard-wired to 0.0) */
	shifts_ref[0] = 0.0;
	shifts_ref[1] = 0.0;

	/* Copy the matrix for the coordinates of sampling points to the video card */
	cudaMemcpy(d_points, points, RING_LENGTH*NRING*2*sizeof(float), cudaMemcpyHostToDevice);
	cudaGetTextureReference(&texPtr, "texim_points");
	cudaBindTexture(0, texPtr, d_points, &channelDesc, RING_LENGTH*NRING*2*sizeof(float));

	/* Copy the matrix for the coordinates of shifts to the video card */
	cudaMemcpy(d_shifts_ref, shifts_ref, 2*sizeof(float), cudaMemcpyHostToDevice);
	cudaGetTextureReference(&texPtr, "texim_shifts");
	cudaBindTexture(offset, texPtr, d_shifts_ref, &channelDesc, 2*sizeof(float));

	/* Copy the matrix for the reference image to the video card */
	cudaMemcpyToArray(ref_image_array, 0, 0, ref_image,  NX*NY*sizeof(float), cudaMemcpyHostToDevice);
	cudaBindTextureToArray(tex, ref_image_array, channelDesc);
	
	/* Convert the reference image to polar coordinates */
	resample_to_polar<<<NRING, RING_LENGTH>>>(d_ref_image_polar, 0.0, 0.0, NX, NY, RING_LENGTH, NRING, *offset);
	cudaThreadSynchronize();
	cudaUnbindTexture(tex);

	/* Conduct FFT of the reference image in polar coordinates */
	cufftExecR2C(plan_ref_image, (cufftReal *)d_ref_image_polar, (cufftComplex *)d_ref_image_polar);

	/* Fill the matrix for the coordinates of shifts for subject images */
	for (i = 0; i < NIMAGE; i++) {
		shifts_subject[i*2] = sx[i];
		shifts_subject[i*2+1] = sy[i];
	} 

	/* Copy the matrix for the coordinates of shifts to the video card */
	cudaMemcpy(d_shifts_subject, shifts_subject, NIMAGE*2*sizeof(float), cudaMemcpyHostToDevice);

	/* Copy the matrix for NIMAGE_ROW*NROW subject images to the video card */
	for (k=0; k<NROW; k++)  
	    cudaMemcpyToArray(subject_image_array[k], 0, 0, subject_image+k*NIMAGE_ROW*NX*NY, NIMAGE_ROW*NX*NY*sizeof(float), cudaMemcpyHostToDevice);
	/* Copy the matrix for NIMAGE_LEFT subject images to the video card */
	if (NIMAGE_LEFT != 0) 
	    cudaMemcpyToArray(subject_image_array_left, 0, 0, subject_image+NROW*NIMAGE_ROW*NX*NY, NIMAGE_LEFT*NX*NY*sizeof(float), cudaMemcpyHostToDevice);

	for (i=-KX; i<=KX; i++) {
		for (j=-KY; j<=KY; j++) {
			
			x = i*STEP;
			y = j*STEP;

			cudaGetTextureReference(&texPtr, "texim_shifts");
		    	for (k=0; k<NROW; k++) { 
		    		cudaBindTextureToArray(tex, subject_image_array[k], channelDesc);
				cudaBindTexture(offset, texPtr, d_shifts_subject+2*NIMAGE_ROW*k, &channelDesc, NIMAGE_ROW*2*sizeof(float));

				/* Convert NIMAGE_ROW subject images to the polar coordinates */
				resample_to_polar<<<GridSize1, RING_LENGTH>>>(d_subject_image_polar+k*NIMAGE_ROW*POINTS_PER_IMAGE, x, y, NX, NY, RING_LENGTH, NRING, *offset);
		    		cudaThreadSynchronize();
		    		cudaUnbindTexture(tex);
			}
			if (NIMAGE_LEFT != 0) {
		    		cudaBindTextureToArray(tex, subject_image_array_left, channelDesc);
				cudaBindTexture(offset, texPtr, d_shifts_subject+2*NIMAGE_ROW*NROW, &channelDesc, NIMAGE_LEFT*2*sizeof(float));

				/* Convert NIMAGE_LEFT subject images to the polar coordinates */
				resample_to_polar<<<GridSize2, RING_LENGTH>>>(d_subject_image_polar+NROW*NIMAGE_ROW*POINTS_PER_IMAGE, x, y, NX, NY, RING_LENGTH, NRING, *offset);
		    		cudaThreadSynchronize();
		    		cudaUnbindTexture(tex);
			}
			cudaUnbindTexture(texim_shifts);
			
			/* Conduct FFT for all subject images */
			for (k=0; k<IMAGE_BATCH1; k++)
				cufftExecR2C(plan_subject_image, (cufftReal *)(d_subject_image_polar+k*IMAGE_PER_BATCH1*POINTS_PER_IMAGE), (cufftComplex*)(d_subject_image_polar+k*IMAGE_PER_BATCH1*POINTS_PER_IMAGE));
			if (IMAGE_LEFT_BATCH1 != 0)
				cufftExecR2C(plan_subject_image_rest, (cufftReal *)(d_subject_image_polar+IMAGE_BATCH1*IMAGE_PER_BATCH1*POINTS_PER_IMAGE), (cufftComplex*)(d_subject_image_polar+IMAGE_BATCH1*IMAGE_PER_BATCH1*POINTS_PER_IMAGE));
			
			cudaGetTextureReference(&texPtr, "texim_ref");
			cudaBindTexture(0, texPtr, d_ref_image_polar, &channelDesc, POINTS_PER_IMAGE*sizeof(float));

			cudaGetTextureReference(&texPtr, "texim_subject");
		    	ccf_base_addr = ((j+KY)*(2*KX+1)+(i+KX))*NIMAGE*(RING_LENGTH+2);
			for (k=0; k<NTEXTURE; k++) {
				cudaBindTexture(0, texPtr, d_subject_image_polar+k*NIMAGE_IN_TEXTURE*POINTS_PER_IMAGE, &channelDesc, NIMAGE_IN_TEXTURE*POINTS_PER_IMAGE*sizeof(float));
				complex_mul<<<NIMAGE_IN_TEXTURE, BLOCK_SIZE>>>(d_ccf+ccf_base_addr+k*NIMAGE_IN_TEXTURE*(RING_LENGTH+2), BLOCK_SIZE, NRING, NIMAGE, KX, KY);
				cudaThreadSynchronize();
			}
			if (NIMAGE_LEFT_TEXTURE != 0) {
				cudaBindTexture(0, texPtr, d_subject_image_polar+NTEXTURE*NIMAGE_IN_TEXTURE*POINTS_PER_IMAGE, &channelDesc, NIMAGE_LEFT_TEXTURE*POINTS_PER_IMAGE*sizeof(float));
				complex_mul<<<NIMAGE_LEFT_TEXTURE, BLOCK_SIZE>>>(d_ccf+ccf_base_addr+NTEXTURE*NIMAGE_IN_TEXTURE*(RING_LENGTH+2), BLOCK_SIZE, NRING, NIMAGE, KX, KY);
				cudaThreadSynchronize();
			}
			cudaUnbindTexture(texim_subject);
			cudaUnbindTexture(texim_ref);
		}
	}
	cudaUnbindTexture(texim_points);

	for (i=0; i<IMAGE_BATCH2; i++)
	    	cufftExecC2R(plan_ccf, (cufftComplex *)(d_ccf+i*POS_PER_IMAGE*IMAGE_PER_BATCH2*(RING_LENGTH+2)), (cufftReal *)(d_ccf+i*POS_PER_IMAGE*IMAGE_PER_BATCH2*(RING_LENGTH+2)));
	if (IMAGE_LEFT_BATCH2 != 0) 
	    	cufftExecC2R(plan_ccf_rest, (cufftComplex *)(d_ccf+IMAGE_BATCH2*POS_PER_IMAGE*IMAGE_PER_BATCH2*(RING_LENGTH+2)), (cufftReal *)(d_ccf+IMAGE_BATCH2*POS_PER_IMAGE*IMAGE_PER_BATCH2*(RING_LENGTH+2)));


	cudaMemcpy(ccf, d_ccf, (RING_LENGTH+2)*NIMAGE*POS_PER_IMAGE*sizeof(float), cudaMemcpyDeviceToHost);

	/* Memory clean up */
	for (k=0; k<NROW; k++)
		cudaFreeArray(subject_image_array[k]);
	if (NIMAGE_LEFT!=0)
		cudaFreeArray(subject_image_array_left);
	cudaFreeArray(ref_image_array);
	cudaFree(d_subject_image_polar);
	cudaFree(d_ref_image_polar);
	cudaFree(d_ccf);
	cudaFree(d_points);
	cudaFree(d_shifts_ref);
	cudaFree(d_shifts_subject);
	cufftDestroy(plan_subject_image);
	cufftDestroy(plan_subject_image_rest);
	cufftDestroy(plan_ref_image);
	cufftDestroy(plan_ccf);
	cufftDestroy(plan_ccf_rest);

	free(offset);
	free(points);
	free(shifts_ref);
	free(shifts_subject);

	return;
}

void calculate_multiref_ccf(float *subject_image, float *ref_image, float *ccf, int NIMAGE, int NREF, int NX, int NY, int RING_LENGTH, int NRING, int OU, float STEP, int KX, int KY, float *sx, float *sy, int silent)
{    
	float *d_subject_image_polar, *d_ref_image_polar;
	float *d_ccf;
	float *points, *d_points;
	float *shifts_ref, *d_shifts_ref;
	float *shifts_subject, *d_shifts_subject;
	int i, j, k, l, index;
	int ccf_base_addr;
	float x, y, ang;
	size_t* offset = (size_t *)malloc(sizeof(size_t));

	int BLOCK_SIZE = RING_LENGTH/2+1;
	int POINTS_PER_IMAGE = (RING_LENGTH+2)*NRING;
	
	// For a texture reference bound to a two-dimensional CUDA array,
	// the maximum width is 2^16 (=65536) and the maximum height is 2^15 (=32768)
	int NIMAGE_ROW = 32768/NX;
	int NROW = NIMAGE/NIMAGE_ROW;
	int NIMAGE_LEFT = NIMAGE%NIMAGE_ROW;
	int NROW_REF = NREF/NIMAGE_ROW;
	int NREF_LEFT = NREF%NIMAGE_ROW;

	// For a texture reference bound to linear memory, the maximum width is 2^27
	int NIMAGE_IN_TEXTURE = (1<<27)/POINTS_PER_IMAGE*9/10;  
	//int NTEXTURE = NIMAGE/NIMAGE_IN_TEXTURE;
	//int NIMAGE_LEFT_TEXTURE = NIMAGE%NIMAGE_IN_TEXTURE;
	int NTEXTURE_REF = NREF/NIMAGE_IN_TEXTURE;
	int NREF_LEFT_TEXTURE = NREF%NIMAGE_IN_TEXTURE;

	int NIMAGE_PER_BATCH_FFT = 65535/NRING;            // The maximum FFT per batch is 65535
	int NIMAGE_BATCH_FFT = NIMAGE/NIMAGE_PER_BATCH_FFT;
	int NIMAGE_LEFT_BATCH_FFT = NIMAGE%NIMAGE_PER_BATCH_FFT;
	int NREF_BATCH_FFT = NREF/NIMAGE_PER_BATCH_FFT;
	int NREF_LEFT_BATCH_FFT = NREF%NIMAGE_PER_BATCH_FFT;

	int POS_PER_IMAGE = 2*(2*KX+1)*(2*KY+1);
	int NIMAGE_PER_BATCH_IFFT = 65535/POS_PER_IMAGE;    // The maximum FFT per batch is 65535
	int NREF_BATCH_IFFT = NREF/NIMAGE_PER_BATCH_IFFT;
	int NREF_LEFT_BATCH_IFFT = NREF%NIMAGE_PER_BATCH_IFFT;

/*	int device, device_count;
	cudaGetDevice(&device);
	cudaGetDeviceCount(&device_count);
	cudaDeviceProp Prop;
	cudaGetDeviceProperties(&Prop, device);
	printf("In multi-ccf: device = %d    name = %s   device count = %d\n", device, Prop.name, device_count);

	cudaError_t cudaError;
	cudaError = cudaGetLastError();
	if( cudaError != cudaSuccess ) {
		fprintf(stderr, "CUDA Runtime API Error reported in multi-ccf: %s\n", cudaGetErrorString(cudaError));
	}
*/

	cudaArray *subject_image_array[NROW], *subject_image_array_left, *ref_image_array[NROW_REF], *ref_image_array_left;
	dim3 GridSize1(NRING, NIMAGE_ROW);
	dim3 GridSize2(NRING, NIMAGE_LEFT);
	dim3 GridSize3(NRING, NREF_LEFT);
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	const textureReference *texPtr;

	tex.normalized = false;
	tex.filterMode = cudaFilterModeLinear;
	tex.addressMode[0] = cudaAddressModeWrap;
	tex.addressMode[1] = cudaAddressModeWrap;

	cufftHandle plan_image, plan_image_left, plan_ref, plan_ref_left, plan_ccf, plan_ccf_left;
	cufftPlan1d(&plan_image, RING_LENGTH, CUFFT_R2C, NRING*NIMAGE_PER_BATCH_FFT);
	if (NIMAGE_LEFT_BATCH_FFT != 0)
		cufftPlan1d(&plan_image_left, RING_LENGTH, CUFFT_R2C, NRING*NIMAGE_LEFT_BATCH_FFT);
	cufftPlan1d(&plan_ref, RING_LENGTH, CUFFT_R2C, NRING*NIMAGE_PER_BATCH_FFT);
	if (NREF_LEFT_BATCH_FFT != 0)
		cufftPlan1d(&plan_ref_left, RING_LENGTH, CUFFT_R2C, NRING*NREF_LEFT_BATCH_FFT);
	cufftPlan1d(&plan_ccf, RING_LENGTH, CUFFT_C2R, POS_PER_IMAGE*NIMAGE_PER_BATCH_IFFT);
	if (NREF_LEFT_BATCH_IFFT != 0)
		cufftPlan1d(&plan_ccf_left, RING_LENGTH, CUFFT_C2R, POS_PER_IMAGE*NREF_LEFT_BATCH_IFFT);

	/* Allocate host memory for the coordinates of sampling points */
	points = (float *)malloc(RING_LENGTH*NRING*2*sizeof(float));
	if (points == 0) {
		fprintf (stderr, "Host memory allocation error!\n");
		return;
	}

	/* Allocate host memory for the coordinates of all shifts */
	shifts_ref = (float *)malloc(NREF*2*sizeof(float));
	if (shifts_ref == 0) {
		fprintf (stderr, "Host memory allocation error!\n");
		return;
	}

	shifts_subject = (float *)malloc(NIMAGE*2*sizeof(float));
	if (shifts_subject == 0) {
		fprintf(stderr, "Host memory allocation error!\n");
		return;
	}

	if (silent == 0) {
		printf("Initialization on the host memory done.\n");

		printf("\nMemory to be allocated on the video card:\n");
		printf("For %5d subject images                        : %10.3f MB\n", NIMAGE, NIMAGE*NX*NY*4.0/1000000.0);
		printf("For %5d reference images                      : %10.3f MB\n", NREF, NREF*NX*NY*4.0/1000000.0);
		printf("For %5d subject images in polar coordinates   : %10.3f MB\n", NIMAGE, NIMAGE*POINTS_PER_IMAGE*4.0/1000000.0);
		printf("For %5d reference images in polar coordinates : %10.3f MB\n", NREF, NREF*POINTS_PER_IMAGE*4.0/1000000.0);
		printf("For all cross-correlation functions (CCF)       : %10.3f MB\n", (RING_LENGTH+2)*NIMAGE*NREF*POS_PER_IMAGE*4.0/1000000.0);
		printf("Total memory used			        : %10.3f MB\n\n", ((NIMAGE+NREF)*(NX*NY+POINTS_PER_IMAGE)+(RING_LENGTH+2)*NIMAGE*NREF*POS_PER_IMAGE)*4.0/1000000.0);
	}

	/* Allocate the matrix for all NIMAGE subject images on the video card */
	for (k=0; k<NROW; k++)
		cudaMallocArray(&subject_image_array[k], &channelDesc, NX, NY*NIMAGE_ROW);
	if (NIMAGE_LEFT != 0)
		cudaMallocArray(&subject_image_array_left, &channelDesc, NX, NY*NIMAGE_LEFT);

	/* Allocate the matrix for all NREF reference images on the video card */
	for (k=0; k<NROW_REF; k++)
		cudaMallocArray(&ref_image_array[k], &channelDesc, NX, NY*NIMAGE_ROW);
	if (NREF_LEFT != 0)
		cudaMallocArray(&ref_image_array_left, &channelDesc, NX, NY*NREF_LEFT);

	/* Allocate the matrix for all NIMAGE subject images in polar coordinates on the video card */
	cudaMalloc((void**)&d_subject_image_polar, NIMAGE*POINTS_PER_IMAGE*sizeof(float));

	/* Allocate the matrix for all NREF reference images in polar coordinates on the video card */
	cudaMalloc((void**)&d_ref_image_polar, NREF*POINTS_PER_IMAGE*sizeof(float));

	/* Allocate the matrix for the ccf on the video card */
	cudaMalloc((void**)&d_ccf, (RING_LENGTH+2)*NIMAGE*NREF*POS_PER_IMAGE*sizeof(float));

	/* Allocate the matrix for the coordinates of the sampling points on the video card */
	cudaMalloc((void**)&d_points, RING_LENGTH*NRING*2*sizeof(float));

	/* Allocate the matrix for the coordinates of the shifts of reference image on the video card */
	cudaMalloc((void**)&d_shifts_ref, NREF*2*sizeof(float));

	/* Allocate the matrix for the coordinates of the shifts of subject images on the video card */
	cudaMalloc((void**)&d_shifts_subject, NIMAGE*2*sizeof(float));

	/* Fill the matrix for the coordinates of sampling points */
	for (i = 0; i < NRING; i++) 
		for (j = 0; j < RING_LENGTH; j++) {
			index = i*RING_LENGTH+j;
			ang = float(j)/RING_LENGTH*PI*2;
			points[index*2] = (i+1)*sinf(ang)*float(OU)/float(NRING);
			points[index*2+1] = (i+1)*cosf(ang)*float(OU)/float(NRING);
		}

	/* Fill the matrix for the coordinates of shifts for reference images (currently hard-wired to 0.0) */
	for (i = 0; i < NREF; i++) {
		shifts_ref[i*2] = 0.0;
		shifts_ref[i*2+1] = 0.0;
	}

	/* Copy the matrix for the coordinates of sampling points to the video card */
	cudaMemcpy(d_points, points, RING_LENGTH*NRING*2*sizeof(float), cudaMemcpyHostToDevice);
	cudaGetTextureReference(&texPtr, "texim_points");
	cudaBindTexture(0, texPtr, d_points, &channelDesc, RING_LENGTH*NRING*2*sizeof(float));

	/* Copy the matrix for the coordinates of shifts to the video card */
	cudaMemcpy(d_shifts_ref, shifts_ref, NREF*2*sizeof(float), cudaMemcpyHostToDevice);
	//cudaGetTextureReference(&texPtr, "texim_shifts");
	//cudaBindTexture(offset, texPtr, d_shifts_ref, &channelDesc, NREF*2*sizeof(float));

	/* Fill the matrix for the coordinates of shifts for subject images */
	for (i = 0; i < NIMAGE; i++) {
		shifts_subject[i*2] = sx[i];
		shifts_subject[i*2+1] = sy[i];
	} 

	/* Copy the matrix for the coordinates of shifts to the video card */
	cudaMemcpy(d_shifts_subject, shifts_subject, NIMAGE*2*sizeof(float), cudaMemcpyHostToDevice);

	/* Copy the matrix for NIMAGE_ROW*NROW subject images to the video card */
	for (k=0; k<NROW; k++)  
	    cudaMemcpyToArray(subject_image_array[k], 0, 0, subject_image+k*NIMAGE_ROW*NX*NY, NIMAGE_ROW*NX*NY*sizeof(float), cudaMemcpyHostToDevice);
	/* Copy the matrix for NIMAGE_LEFT subject images to the video card */
	if (NIMAGE_LEFT != 0) 
	    cudaMemcpyToArray(subject_image_array_left, 0, 0, subject_image+NROW*NIMAGE_ROW*NX*NY, NIMAGE_LEFT*NX*NY*sizeof(float), cudaMemcpyHostToDevice);

	/* Copy the matrix for NIMAGE_ROW*NROW_REF reference images to the video card */
	for (k=0; k<NROW_REF; k++)  
	    cudaMemcpyToArray(ref_image_array[k], 0, 0, ref_image+k*NIMAGE_ROW*NX*NY, NIMAGE_ROW*NX*NY*sizeof(float), cudaMemcpyHostToDevice);
	/* Copy the matrix for NREF_LEFT reference images to the video card */
	if (NREF_LEFT != 0)
	    cudaMemcpyToArray(ref_image_array_left, 0, 0, ref_image+NROW_REF*NIMAGE_ROW*NX*NY, NREF_LEFT*NX*NY*sizeof(float), cudaMemcpyHostToDevice);

	cudaGetTextureReference(&texPtr, "texim_shifts");
	for (k=0; k<NROW_REF; k++) { 
		cudaBindTextureToArray(tex, ref_image_array[k], channelDesc);
		cudaBindTexture(offset, texPtr, d_shifts_ref+2*NIMAGE_ROW*k, &channelDesc, NIMAGE_ROW*2*sizeof(float));

		/* Convert NIMAGE_ROW reference images to the polar coordinates */
		resample_to_polar<<<GridSize1, RING_LENGTH>>>(d_ref_image_polar+k*NIMAGE_ROW*POINTS_PER_IMAGE, 0.0f, 0.0f, NX, NY, RING_LENGTH, NRING, *offset);
		cudaThreadSynchronize();
		cudaUnbindTexture(tex);
	}
	if (NREF_LEFT != 0) {
		cudaBindTextureToArray(tex, ref_image_array_left, channelDesc);
		cudaBindTexture(offset, texPtr, d_shifts_ref+2*NIMAGE_ROW*NROW_REF, &channelDesc, NREF_LEFT*2*sizeof(float));

		/* Convert NREF_LEFT reference images to the polar coordinates */
		resample_to_polar<<<GridSize3, RING_LENGTH>>>(d_ref_image_polar+NROW_REF*NIMAGE_ROW*POINTS_PER_IMAGE, 0.0f, 0.0f, NX, NY, RING_LENGTH, NRING, *offset);
		cudaThreadSynchronize();
		cudaUnbindTexture(tex);
	}
	cudaUnbindTexture(texim_shifts);

	cudaGetTextureReference(&texPtr, "texim_polar");
	for (k=0; k<NTEXTURE_REF; k++) {
		cudaBindTexture(offset, texPtr, d_ref_image_polar+k*NIMAGE_IN_TEXTURE*POINTS_PER_IMAGE, &channelDesc, NIMAGE_IN_TEXTURE*POINTS_PER_IMAGE*sizeof(float));
		/* Normalize all NIMAGE_IN_TEXTURE images */
		normalize_rings<<<NIMAGE_IN_TEXTURE, 1>>>(d_ref_image_polar+k*NIMAGE_IN_TEXTURE*POINTS_PER_IMAGE, RING_LENGTH, NRING, *offset);
		cudaThreadSynchronize();
	}
	if (NREF_LEFT_TEXTURE!=0) {
		cudaBindTexture(offset, texPtr, d_ref_image_polar+NTEXTURE_REF*NIMAGE_IN_TEXTURE*POINTS_PER_IMAGE, &channelDesc, NREF_LEFT_TEXTURE*POINTS_PER_IMAGE*sizeof(float));
		/* Normalize all NREF_LEFT_TEXTURE images */
		normalize_rings<<<NREF_LEFT_TEXTURE, 1>>>(d_ref_image_polar+NTEXTURE_REF*NIMAGE_IN_TEXTURE*POINTS_PER_IMAGE, RING_LENGTH, NRING, *offset);
		cudaThreadSynchronize();
	}
	cudaUnbindTexture(texim_polar);

	/* Conduct FFT for all reference images */
	for (k=0; k<NREF_BATCH_FFT; k++)
		cufftExecR2C(plan_ref, (cufftReal *)(d_ref_image_polar+k*NIMAGE_PER_BATCH_FFT*POINTS_PER_IMAGE), (cufftComplex *)(d_ref_image_polar+k*NIMAGE_PER_BATCH_FFT*POINTS_PER_IMAGE));
	if (NREF_LEFT_BATCH_FFT != 0)
		cufftExecR2C(plan_ref_left, (cufftReal *)(d_ref_image_polar+NREF_BATCH_FFT*NIMAGE_PER_BATCH_FFT*POINTS_PER_IMAGE), (cufftComplex *)(d_ref_image_polar+NREF_BATCH_FFT*NIMAGE_PER_BATCH_FFT*POINTS_PER_IMAGE));

	for (i=-KX; i<=KX; i++) {
		for (j=-KY; j<=KY; j++) {
			x = i*STEP;
			y = j*STEP;

			cudaGetTextureReference(&texPtr, "texim_shifts");
		    	for (k=0; k<NROW; k++) { 
		    		cudaBindTextureToArray(tex, subject_image_array[k], channelDesc);
				cudaBindTexture(offset, texPtr, d_shifts_subject+2*NIMAGE_ROW*k, &channelDesc, NIMAGE_ROW*2*sizeof(float));

				/* Convert NIMAGE_ROW subject images to the polar coordinates */
				resample_to_polar<<<GridSize1, RING_LENGTH>>>(d_subject_image_polar+k*NIMAGE_ROW*POINTS_PER_IMAGE, x, y, NX, NY, RING_LENGTH, NRING, *offset);
		    		cudaThreadSynchronize();
		    		cudaUnbindTexture(tex);
			}
			if (NIMAGE_LEFT != 0) {
		    		cudaBindTextureToArray(tex, subject_image_array_left, channelDesc);
				cudaBindTexture(offset, texPtr, d_shifts_subject+2*NIMAGE_ROW*NROW, &channelDesc, NIMAGE_LEFT*2*sizeof(float));

				/* Convert NIMAGE_LEFT subject images to the polar coordinates */
				resample_to_polar<<<GridSize2, RING_LENGTH>>>(d_subject_image_polar+NROW*NIMAGE_ROW*POINTS_PER_IMAGE, x, y, NX, NY, RING_LENGTH, NRING, *offset);
		    		cudaThreadSynchronize();
		    		cudaUnbindTexture(tex);
			}
			cudaUnbindTexture(texim_shifts);
			
			/* Conduct FFT for all subject images */
			for (k=0; k<NIMAGE_BATCH_FFT; k++)
				cufftExecR2C(plan_image, (cufftReal *)(d_subject_image_polar+k*NIMAGE_PER_BATCH_FFT*POINTS_PER_IMAGE), (cufftComplex *)(d_subject_image_polar+k*NIMAGE_PER_BATCH_FFT*POINTS_PER_IMAGE));
			if (NIMAGE_LEFT_BATCH_FFT != 0)
				cufftExecR2C(plan_image_left, (cufftReal *)(d_subject_image_polar+NIMAGE_BATCH_FFT*NIMAGE_PER_BATCH_FFT*POINTS_PER_IMAGE), (cufftComplex *)(d_subject_image_polar+NIMAGE_BATCH_FFT*NIMAGE_PER_BATCH_FFT*POINTS_PER_IMAGE));
			
			for (k=0; k<NIMAGE; k++) {
				cudaGetTextureReference(&texPtr, "texim_ref");
				cudaBindTexture(0, texPtr, d_subject_image_polar+k*POINTS_PER_IMAGE, &channelDesc, POINTS_PER_IMAGE*sizeof(float));
			    	ccf_base_addr = ((j+KY)*(2*KX+1)+(i+KX))*NREF*(RING_LENGTH+2)+k*NREF*POS_PER_IMAGE*(RING_LENGTH+2);
				for (l=0; l<NTEXTURE_REF; l++) {
					cudaGetTextureReference(&texPtr, "texim_subject");
					cudaBindTexture(0, texPtr, d_ref_image_polar+l*NIMAGE_IN_TEXTURE*POINTS_PER_IMAGE, &channelDesc, NIMAGE_IN_TEXTURE*POINTS_PER_IMAGE*sizeof(float));
					complex_mul<<<NIMAGE_IN_TEXTURE, BLOCK_SIZE>>>(d_ccf+ccf_base_addr+l*NIMAGE_IN_TEXTURE*(RING_LENGTH+2), BLOCK_SIZE, NRING, NREF, KX, KY);
					cudaThreadSynchronize();
					cudaUnbindTexture(texim_subject);
				}
				if (NREF_LEFT_TEXTURE != 0) {
					cudaGetTextureReference(&texPtr, "texim_subject");
					cudaBindTexture(0, texPtr, d_ref_image_polar+NTEXTURE_REF*NIMAGE_IN_TEXTURE*POINTS_PER_IMAGE, &channelDesc, NIMAGE_IN_TEXTURE*POINTS_PER_IMAGE*sizeof(float));
					complex_mul<<<NREF_LEFT_TEXTURE, BLOCK_SIZE>>>(d_ccf+ccf_base_addr+NTEXTURE_REF*NIMAGE_IN_TEXTURE*(RING_LENGTH+2), BLOCK_SIZE, NRING, NREF, KX, KY);
					cudaThreadSynchronize();
					cudaUnbindTexture(texim_subject);
				}
				cudaUnbindTexture(texim_ref);
			}
		}
	}
	cudaUnbindTexture(texim_points);

	for (k=0; k<NIMAGE; k++) {
		ccf_base_addr = k*NREF*POS_PER_IMAGE*(RING_LENGTH+2);
		for (i=0; i<NREF_BATCH_IFFT; i++)
		    	cufftExecC2R(plan_ccf, (cufftComplex *)(d_ccf+ccf_base_addr+i*POS_PER_IMAGE*NIMAGE_PER_BATCH_IFFT*(RING_LENGTH+2)), (cufftReal *)(d_ccf+ccf_base_addr+i*POS_PER_IMAGE*NIMAGE_PER_BATCH_IFFT*(RING_LENGTH+2)));
		if (NREF_LEFT_BATCH_IFFT != 0) 
		    	cufftExecC2R(plan_ccf_left, (cufftComplex *)(d_ccf+ccf_base_addr+NREF_BATCH_IFFT*POS_PER_IMAGE*NIMAGE_PER_BATCH_IFFT*(RING_LENGTH+2)), (cufftReal *)(d_ccf+ccf_base_addr+NREF_BATCH_IFFT*POS_PER_IMAGE*NIMAGE_PER_BATCH_IFFT*(RING_LENGTH+2)));
	}
	
	cudaMemcpy(ccf, d_ccf, (RING_LENGTH+2)*NIMAGE*NREF*POS_PER_IMAGE*sizeof(float), cudaMemcpyDeviceToHost);

	/* Memory clean up */
	for (k=0; k<NROW; k++)
		cudaFreeArray(subject_image_array[k]);
	if (NIMAGE_LEFT!=0)
		cudaFreeArray(subject_image_array_left);
	for (k=0; k<NROW_REF; k++)
		cudaFreeArray(ref_image_array[k]);
	if (NREF_LEFT!=0)
		cudaFreeArray(ref_image_array_left);
	cudaFree(d_subject_image_polar);
	cudaFree(d_ref_image_polar);
	cudaFree(d_ccf);
	cudaFree(d_points);
	cudaFree(d_shifts_ref);
	cudaFree(d_shifts_subject);
	cufftDestroy(plan_image);
	cufftDestroy(plan_image_left);
	cufftDestroy(plan_ref);
	cufftDestroy(plan_ref_left);
	cufftDestroy(plan_ccf);
	cufftDestroy(plan_ccf_left);

	free(offset);
	free(points);
	free(shifts_ref);
	free(shifts_subject);

	return;
}


void filter_image(float *image_in, float *image_out, int NIMA, int NX, int NY, float *params) {
	
	float *image_padded, *d_image_padded;
	int padded_size = (NX*2+2)*(NY*2);
	int im, iy;
	
	cufftHandle plan_R2C, plan_C2R;
	cufftPlan2d(&plan_R2C, NX*2, NY*2, CUFFT_R2C);
	cufftPlan2d(&plan_C2R, NX*2, NY*2, CUFFT_C2R);
	
	image_padded = (float *)malloc(padded_size*NIMA*sizeof(float));
	memset(image_padded, 0, padded_size*NIMA*sizeof(float));
	
	cudaMalloc((void **)&d_image_padded, padded_size*NIMA*sizeof(float));
	cudaMemset(d_image_padded, 0, padded_size*NIMA*sizeof(float));

	for (im=0; im<NIMA; im++) 
		for (iy=0; iy<NY; iy++)
			memmove(image_padded+im*padded_size+(iy+NY/2)*(NX*2+2)+NX/2, image_in+im*NX*NY+iy*NX, NX*sizeof(float));

	cudaMemcpy(d_image_padded, image_padded, padded_size*NIMA*sizeof(float), cudaMemcpyHostToDevice);

	for (im=0; im<NIMA; im++) {
		int base_address = im*padded_size;
		cufftExecR2C(plan_R2C, (cufftReal *)(d_image_padded+base_address), (cufftComplex *)(d_image_padded+base_address));
		mul_ctf<<<NX+1, NY*2>>>(d_image_padded+base_address, NX*2, NY*2, params[im*6], params[im*6+1], params[im*6+2], params[im*6+3], params[im*6+4], params[im*6+5]);
		cudaThreadSynchronize();
		cufftExecC2R(plan_C2R, (cufftComplex *)(d_image_padded+base_address), (cufftReal *)(d_image_padded+base_address));
	}
	
	cudaMemcpy(image_padded, d_image_padded, padded_size*NIMA*sizeof(float), cudaMemcpyDeviceToHost);

	for (im=0; im<NIMA; im++)
		for (iy=0; iy<NY; iy++)
			memmove(image_out+im*NX*NY+iy*NX, image_padded+im*padded_size+(iy+NY/2)*(NX*2+2)+NX/2, NX*sizeof(float));	
		
	cufftDestroy(plan_R2C);
	cufftDestroy(plan_C2R);
	cudaFree(d_image_padded);
	
	free(image_padded);

	return;
}

void rot_filt_sum(float *image, int NIMA, int NX, int NY, int CTF, float *ctf_params, float *ali_params, float *ave1, float *ave2) {

	float *d_image_padded, *d_ave1, *d_ave2;
	float *d_ali_params;
	size_t* offset = (size_t *)malloc(sizeof(size_t));
	int padded_size = (NX*2+2)*(NY*2);
	int k, im;

	int NIMAGE_ROW = 32768/NX;
	int NROW = NIMA/NIMAGE_ROW;
	int NIMAGE_LEFT = NIMA%NIMAGE_ROW;

	cufftHandle plan_R2C, plan_C2R;
	cufftPlan2d(&plan_R2C, NX*2, NY*2, CUFFT_R2C);
	cufftPlan2d(&plan_C2R, NX*2, NY*2, CUFFT_C2R);

	cudaArray *image_array[NROW], *image_array_left;
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	const textureReference *texPtr;
	dim3 GridSize1(NIMAGE_ROW, NY);
	dim3 GridSize2(NIMAGE_LEFT, NY);

	tex.normalized = false;
	tex.filterMode = cudaFilterModeLinear;
	tex.addressMode[0] = cudaAddressModeWrap;
	tex.addressMode[1] = cudaAddressModeWrap;

	/* Allocate the matrix for all NIMA subject images on the video card */
	for (k=0; k<NROW; k++)
		cudaMallocArray(&image_array[k], &channelDesc, NX, NY*NIMAGE_ROW);
	if (NIMAGE_LEFT != 0)
		cudaMallocArray(&image_array_left, &channelDesc, NX, NY*NIMAGE_LEFT);

	cudaMalloc((void **)&d_image_padded, padded_size*NIMA*sizeof(float));
	cudaMalloc((void **)&d_ali_params, NIMA*4*sizeof(float));
	cudaMalloc((void **)&d_ave1, NX*NY*sizeof(float));
	cudaMalloc((void **)&d_ave2, NX*NY*sizeof(float));

	/* Copy the matrix for NIMAGE_ROW*NROW subject images to the video card */
	for (k=0; k<NROW; k++)  
		cudaMemcpyToArray(image_array[k], 0, 0, image+k*NIMAGE_ROW*NX*NY, NIMAGE_ROW*NX*NY*sizeof(float), cudaMemcpyHostToDevice);
	/* Copy the matrix for NIMAGE_LEFT subject images to the video card */
	if (NIMAGE_LEFT != 0)
		cudaMemcpyToArray(image_array_left, 0, 0, image+NROW*NIMAGE_ROW*NX*NY, NIMAGE_LEFT*NX*NY*sizeof(float), cudaMemcpyHostToDevice);
	/* Copy the matrix for the alignment parameters to the video card */
	cudaMemcpy(d_ali_params, ali_params, NIMA*4*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemset(d_image_padded, 0, padded_size*NIMA*sizeof(float));

	cudaGetTextureReference(&texPtr, "texim_ali_params");

	for (k=0; k<NROW; k++) { 
		cudaBindTextureToArray(tex, image_array[k], channelDesc);
		cudaBindTexture(offset, texPtr, d_ali_params+4*NIMAGE_ROW*k, &channelDesc, NIMAGE_ROW*4*sizeof(float));

		/* Rotate and shift NIMAGE_ROW subject images */
		rotate_shift<<<GridSize1, NX>>>(d_image_padded+k*NIMAGE_ROW*padded_size, NX, NY, *offset);
	        cudaThreadSynchronize();
		cudaUnbindTexture(tex);
	}
	if (NIMAGE_LEFT != 0) {
		cudaBindTextureToArray(tex, image_array_left, channelDesc);
		cudaBindTexture(offset, texPtr, d_ali_params+4*NIMAGE_ROW*NROW, &channelDesc, NIMAGE_LEFT*4*sizeof(float));

		/* Rotate and shift NIMAGE_LEFT subject images */
		rotate_shift<<<GridSize2, NX>>>(d_image_padded+NROW*NIMAGE_ROW*padded_size, NX, NY, *offset);
	       	cudaThreadSynchronize();
		cudaUnbindTexture(tex);
	}

	if (CTF == 1) {
		/* Transform to Fourier space, multiple CTF, and transform back to real space */
		for (im=0; im<NIMA; im++) {
			int base_address = im*padded_size;
			cufftExecR2C(plan_R2C, (cufftReal *)(d_image_padded+base_address), (cufftComplex *)(d_image_padded+base_address));
			mul_ctf<<<NX+1, NY*2>>>(d_image_padded+base_address, NX*2, NY*2, ctf_params[im*6], ctf_params[im*6+1], ctf_params[im*6+2], ctf_params[im*6+3], ctf_params[im*6+4], ctf_params[im*6+5]);
			cudaThreadSynchronize();
			cufftExecC2R(plan_C2R, (cufftComplex *)(d_image_padded+base_address), (cufftReal *)(d_image_padded+base_address));
		}
	}

	add_img<<<NX, NY>>>(d_image_padded, d_ave1, d_ave2, NX, NY, NIMA);
	cudaThreadSynchronize();
	cudaMemcpy(ave1, d_ave1, NX*NY*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(ave2, d_ave2, NX*NY*sizeof(float), cudaMemcpyDeviceToHost);

	/* memory cleanup */
	cufftDestroy(plan_R2C);
	cufftDestroy(plan_C2R);
	cudaFree(d_image_padded);
	cudaFree(d_ave1);
	cudaFree(d_ave2);
	cudaFree(d_ali_params);
	for (k=0; k<NROW; k++)
		cudaFreeArray(image_array[k]);
	if (NIMAGE_LEFT != 0)
		cudaFreeArray(image_array_left);
	free(offset);

	return;
}

__global__ void add_img(float *image_padded, float *ave1, float *ave2, int nx, int ny, int nima) {
	
	// Block index
	int bx = blockIdx.x;

	// Thread index
	int tx = threadIdx.x;

	float sum1 = 0.0;
	float sum2 = 0.0;
	int index = tx+bx*nx;
	int index2 = tx+(nx>>1)+(bx+(ny>>1))*(nx*2+2);    

	for (int i=0; i<nima; i+=2) sum1 += image_padded[index2+i*(nx*2+2)*ny*2];
	for (int i=1; i<nima; i+=2) sum2 += image_padded[index2+i*(nx*2+2)*ny*2];
	ave1[index] = sum1;
	ave2[index] = sum2;

	return;
}


__global__ void mul_ctf(float *image, int nx, int ny, float defocus, float cs, float voltage, float apix, float bfactor, float ampcont) {

	// Block index
	int bx = blockIdx.x;

	// Thread index
	int tx = threadIdx.x;

	float x, y;

	x = float(bx);
	if (tx >= ny>>1) y = float(tx-ny);
	else y = float(tx);
	int index = bx*2+tx*(nx+2);

	float ak = sqrt(x*x+y*y)/nx/apix;
	float cst = cs*1.0e7f;
	float wgh = ampcont/100.0f;
	float phase = atan(wgh/sqrt(1.0f-wgh*wgh));
	float lambda = 12.398f/sqrt(voltage*(1022.f+voltage));
	float ak2 = ak*ak;
	float g1 = defocus*1.0e4f*lambda*ak2;
	float g2 = cst*lambda*lambda*lambda*ak2*ak2/2.0f;
	float ctfv = sin(PI*(g1-g2)+phase);
	if (bfactor != 0.0f)  ctfv *= exp(-bfactor*ak2/4.0f);

	image[index] *= ctfv;
	image[index+1] *= ctfv;
}

__global__ void complex_mul(float *ccf, int BLOCK_SIZE, int NRING, int NIMAGE, int KX, int KY) {
	
	// Block index
	int bx = blockIdx.x;

	// Thread index
	int tx = threadIdx.x;

	float ccf_s_real = 0.0;
	float ccf_s_imag = 0.0;
	float ccf_m_real = 0.0;
	float ccf_m_imag = 0.0;
	float sub_real, sub_imag, ref_real, ref_imag;
	int i_block, b_image, subject_index, ref_index, ccf_index;
	float rr_sr, ri_si, rr_si, ri_sr;

	b_image = bx*BLOCK_SIZE*NRING;
	for (int i = 0; i < NRING; i++) {
		i_block = i*BLOCK_SIZE;
		subject_index = (b_image+i_block+tx)*2;
		ref_index = (i_block+tx)*2;

		sub_real = tex1Dfetch(texim_subject, subject_index);
		sub_imag = tex1Dfetch(texim_subject, subject_index+1);
		ref_real = tex1Dfetch(texim_ref, ref_index);
		ref_imag = tex1Dfetch(texim_ref, ref_index+1);
		rr_sr = ref_real*sub_real;
		ri_si = ref_imag*sub_imag;
		rr_si = ref_real*sub_imag;
		ri_sr = ref_imag*sub_real;
		ccf_s_real += (rr_sr+ri_si)*(i+1);
		ccf_s_imag += (-rr_si+ri_sr)*(i+1);
		ccf_m_real += (rr_sr-ri_si)*(i+1);
		ccf_m_imag += (-rr_si-ri_sr)*(i+1);
	}

	ccf_index = (bx*BLOCK_SIZE+tx)*2;
	int ccf_OFFSET = BLOCK_SIZE*2*NIMAGE*(2*KX+1)*(2*KY+1);
	ccf[ccf_index] = ccf_s_real;
	ccf[ccf_index+1] = ccf_s_imag;
	ccf[ccf_index+ccf_OFFSET] = ccf_m_real;
	ccf[ccf_index+ccf_OFFSET+1] = ccf_m_imag;
}

__global__ void resample_to_polar(float* image, float sx, float sy, int NX, int NY, int RING_LENGTH, int NRING, int offset) {
	
	// Block index
	int bx = blockIdx.x;
	int by = blockIdx.y;

	// Thread index
	int tx = threadIdx.x;

	float cnx =	  (NX>>1)+sx+0.5f;
	float cny = by*NY+(NY>>1)+sy+0.5f;

	int img_index = (by*NRING+bx)*(RING_LENGTH+2)+tx;
	int index = bx*RING_LENGTH+tx;
	float points_x = tex1Dfetch(texim_points, index*2);
	float points_y = tex1Dfetch(texim_points, index*2+1);
	float shifts_x = tex1Dfetch(texim_shifts, by*2+offset/sizeof(float));
	float shifts_y = tex1Dfetch(texim_shifts, by*2+1+offset/sizeof(float));

	image[img_index] = tex2D(tex, cnx+points_x+shifts_x, cny+points_y+shifts_y);
}

__global__ void normalize_rings(float *image, int RING_LENGTH, int NRING, int offset) {

	// Block index
	int bx = blockIdx.x;

	int img_base = bx*(RING_LENGTH+2)*NRING;
	int i, j;

	float sum = 0.0f, var = 0.0f, p = 0.0f;
	float w = 0.0f;
	for (i=0; i<NRING; i++) 
		for (j=0; j<RING_LENGTH; j++) {
			p = tex1Dfetch(texim_polar, img_base+i*(RING_LENGTH+2)+j+offset/sizeof(float));
			sum += p*(i+1);
			var += p*p*(i+1);
			w += (i+1);
			//sum += p;
			//var += p*p;
		}
	//int n = RING_LENGTH*NRING;
	//var = sqrt((var-sum*sum/n)/(n-1));
	//sum /= n;
	var = sqrt((var-sum*sum/w)/w);
	sum /= w; 
	
	for (i=0; i<NRING; i++) 
		for (j=0; j<RING_LENGTH; j++) {
			image[img_base+i*(RING_LENGTH+2)+j] -= sum;
			image[img_base+i*(RING_LENGTH+2)+j] /= var;
		}
}

__global__ void rotate_shift(float* image, int nx, int ny, int offset) {
	
	// Block index
	int bx = blockIdx.x;
	int by = blockIdx.y;
	
	// Thread index
	int tx = threadIdx.x;
	
	float alpha = tex1Dfetch(texim_ali_params, bx*4+offset/sizeof(float));
	float sx = tex1Dfetch(texim_ali_params, bx*4+1+offset/sizeof(float));
	float sy = tex1Dfetch(texim_ali_params, bx*4+2+offset/sizeof(float));
	float mirror = tex1Dfetch(texim_ali_params, bx*4+3+offset/sizeof(float));

	float x, y;
	
	int img_index = bx*(nx*2+2)*(ny*2)+(by+(ny>>1))*(nx*2+2)+tx+(nx>>1);
	
	if  (mirror > 0.5) x = nx-tx-sx-(nx>>1);
	else x = tx-sx-(nx>>1);
	y = by-sy-(ny>>1);
	
	float cosa = __cosf(alpha*PI/180.0f);
	float sina = __sinf(alpha*PI/180.0f);
	
	image[img_index] = tex2D(tex, x*cosa-y*sina+(nx>>1)+0.5f, x*sina+y*cosa+(ny>>1)+bx*ny+0.5f); 
}
