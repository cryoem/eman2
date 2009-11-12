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

/* This example demonstrates how to use the CUBLAS library
 * to realize apmq.
 */

/* Includes, system */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Includes, cuda, cufft */
#include <cufft.h>

/* Matrix size */
#define PI (3.14159265358979)
#define NIMAGE_ROW (512)
#define NIMAGE_IN_TEXTURE (15000)

__global__ void complex_mul(float *ccf, int BLOCK_SIZE, int NRING, int NIMAGE, int KX, int KY);
__global__ void resample_to_polar(float* image, float dx, float dy, int NX, int NY, int RING_LENGTH, int NRING);
__global__ void mul_ctf(float *image, int nx, int ny, float defocus, float cs, float voltage, float apix, float bfactor, float ampcont);

texture<float, 2, cudaReadModeElementType> tex;
texture<float, 1, cudaReadModeElementType> texim_subject;
texture<float, 1, cudaReadModeElementType> texim_ref;
texture<float, 1, cudaReadModeElementType> texim_points;
texture<float, 1, cudaReadModeElementType> texim_shifts;
 
/* Main */
void calculate_ccf(float *subject_image, float *ref_image, float *ccf, int NIMAGE, int NX, int NY, int RING_LENGTH, int NRING, int OU, float STEP, int KX, int KY, float *sx, float *sy, int id, int silent)
{    
    float *d_subject_image_polar, *d_ref_image_polar;
    float *d_ccf;
    float *points, *d_points;
    float *shifts_ref, *d_shifts_ref;
    float *shifts_subject, *d_shifts_subject;
    int i, j, k, index;
    int ccf_base_addr;
    float x, y, ang;

    int BLOCK_SIZE = RING_LENGTH/2+1;
    int NROW = NIMAGE/NIMAGE_ROW;
    int NIMAGE_LEFT = NIMAGE%NIMAGE_ROW;
    int NTEXTURE = NIMAGE/NIMAGE_IN_TEXTURE;
    int NIMAGE_LEFT_TEXTURE = NIMAGE%NIMAGE_IN_TEXTURE;

    int IMAGE_PER_BATCH1 = 65535/NRING;
    int IMAGE_BATCH1 = NIMAGE/IMAGE_PER_BATCH1;
    int IMAGE_LEFT_BATCH1 = NIMAGE%IMAGE_PER_BATCH1;
    int POINTS_PER_IMAGE = 2*(2*KX+1)*(2*KY+1);
    int IMAGE_PER_BATCH2 = 65535/POINTS_PER_IMAGE;
    int IMAGE_BATCH2 = NIMAGE/IMAGE_PER_BATCH2;
    int IMAGE_LEFT_BATCH2 = NIMAGE%IMAGE_PER_BATCH2;

    // Unblock this line you have multiple GPU
    //cudaSetDevice(id); 

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
    cufftPlan1d(&plan_ccf, RING_LENGTH, CUFFT_C2R, POINTS_PER_IMAGE*IMAGE_PER_BATCH2);
    if (IMAGE_LEFT_BATCH2 != 0)
	cufftPlan1d(&plan_ccf_rest, RING_LENGTH, CUFFT_C2R, POINTS_PER_IMAGE*IMAGE_LEFT_BATCH2);
    
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
	printf("For %5d subject images                      : %10.3f MB\n", NIMAGE, NIMAGE*NX*NY*4/1000000.0);
	printf("For reference image                           : %10.3f KB\n", NX*NY*4/1000.0);
	printf("For %5d subject images in polar coordinates : %10.3f MB\n", NIMAGE, NIMAGE*(RING_LENGTH+2)*NRING*4/1000000.0);
	printf("For reference image in polar coordinates      : %10.3f KB\n", (RING_LENGTH+2)*NRING*4/1000.0);
	printf("For all cross-correlation functions (CCF)     : %10.3f MB\n", (RING_LENGTH+2)*NIMAGE*POINTS_PER_IMAGE*4/1000000.0);
	printf("Total memory used                             : %10.3f MB\n\n", ((NIMAGE+1)*(NX*NY+(RING_LENGTH+2)*NRING)+(RING_LENGTH+2)*NIMAGE*POINTS_PER_IMAGE)*4/1000000.0);
    }


    /* Allocate the matrix for all NIMAGE subject images on the video card */
    for (k=0; k<NROW; k++)
       cudaMallocArray(&subject_image_array[k], &channelDesc, NX, NY*NIMAGE_ROW);
    if (NIMAGE_LEFT != 0)
       cudaMallocArray(&subject_image_array_left, &channelDesc, NX, NY*NIMAGE_LEFT);
 
    /* Allocate the matrix for the reference image on the video card */
    cudaMallocArray(&ref_image_array, &channelDesc, NX, NY);
 
    /* Allocate the matrix for all NIMAGE subject images in polar coordinates on the video card */
    cudaMalloc((void**)&d_subject_image_polar, NIMAGE*(RING_LENGTH+2)*NRING*sizeof(float));
   
    /* Allocate the matrix for the reference image in polar coordinates on the video card */
    cudaMalloc((void**)&d_ref_image_polar, (RING_LENGTH+2)*NRING*sizeof(float));
 
    /* Allocate the matrix for the ccf on the video card */
    cudaMalloc((void**)&d_ccf, (RING_LENGTH+2)*NIMAGE*POINTS_PER_IMAGE*sizeof(float));

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
    for (i = 0; i < 1; i++) {
	shifts_ref[i*2] = 0.0;
	shifts_ref[i*2+1] = 0.0;
    }

    /* Copy the matrix for the coordinates of sampling points to the video card */
    cudaMemcpy(d_points, points, RING_LENGTH*NRING*2*sizeof(float), cudaMemcpyHostToDevice);
    cudaGetTextureReference(&texPtr, "texim_points");
    cudaBindTexture(0, texPtr, d_points, &channelDesc, RING_LENGTH*NRING*2*sizeof(float));

    /* Copy the matrix for the coordinates of shifts to the video card */
    cudaMemcpy(d_shifts_ref, shifts_ref, 2*sizeof(float), cudaMemcpyHostToDevice);
    cudaGetTextureReference(&texPtr, "texim_shifts");
    cudaBindTexture(0, texPtr, d_shifts_ref, &channelDesc, 2*sizeof(float));

    /* Copy the matrix for the reference image to the video card */
    cudaMemcpyToArray(ref_image_array, 0, 0, ref_image,  NX*NY*sizeof(float), cudaMemcpyHostToDevice);
    cudaBindTextureToArray(tex, ref_image_array, channelDesc);
 
    /* Convert the reference image to polar coordinates */
    resample_to_polar<<<NRING, RING_LENGTH>>>(d_ref_image_polar, 0.0, 0.0, NX, NY, RING_LENGTH, NRING);
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
			cudaBindTexture(0, texPtr, d_shifts_subject+2*NIMAGE_ROW*k, &channelDesc, NIMAGE_ROW*2*sizeof(float));

			/* Convert NIMAGE_ROW subject images to the polar coordinates */
			resample_to_polar<<<GridSize1, RING_LENGTH>>>(d_subject_image_polar+k*NIMAGE_ROW*(RING_LENGTH+2)*NRING, x, y, NX, NY, RING_LENGTH, NRING);
                	cudaThreadSynchronize();
     			cudaUnbindTexture(tex);
		}
		if (NIMAGE_LEFT != 0) {
    			cudaBindTextureToArray(tex, subject_image_array_left, channelDesc);
 			cudaBindTexture(0, texPtr, d_shifts_subject+2*NIMAGE_ROW*NROW, &channelDesc, NIMAGE_LEFT*2*sizeof(float));
  		
			/* Convert NIMAGE_LEFT subject images to the polar coordinates */
			resample_to_polar<<<GridSize2, RING_LENGTH>>>(d_subject_image_polar+NROW*NIMAGE_ROW*(RING_LENGTH+2)*NRING, x, y, NX, NY, RING_LENGTH, NRING);
                	cudaThreadSynchronize();
    			cudaUnbindTexture(tex);
		}
		cudaUnbindTexture(texim_shifts);
		
		/* Conduct FFT for all subject images */
		for (k=0; k<IMAGE_BATCH1; k++)
			cufftExecR2C(plan_subject_image, (cufftReal *)(d_subject_image_polar+k*NRING*IMAGE_PER_BATCH1*(RING_LENGTH+2)), (cufftComplex *)(d_subject_image_polar+k*NRING*IMAGE_PER_BATCH1*(RING_LENGTH+2)));
		if (IMAGE_LEFT_BATCH1 != 0)
			cufftExecR2C(plan_subject_image_rest, (cufftReal *)(d_subject_image_polar+IMAGE_BATCH1*NRING*IMAGE_PER_BATCH1*(RING_LENGTH+2)), (cufftComplex *)(d_subject_image_polar+IMAGE_BATCH1*NRING*IMAGE_PER_BATCH1*(RING_LENGTH+2)));
		
		cudaGetTextureReference(&texPtr, "texim_ref");
		cudaBindTexture(0, texPtr, d_ref_image_polar, &channelDesc, (RING_LENGTH+2)*NRING*sizeof(float));
   	
		cudaGetTextureReference(&texPtr, "texim_subject");
                ccf_base_addr = ((j+KY)*(2*KX+1)+(i+KX))*NIMAGE*(RING_LENGTH+2);
		for (k=0; k<NTEXTURE; k++) {
			cudaBindTexture(0, texPtr, d_subject_image_polar+k*NIMAGE_IN_TEXTURE*(RING_LENGTH+2)*NRING, &channelDesc, NIMAGE_IN_TEXTURE*(RING_LENGTH+2)*NRING*sizeof(float));
			complex_mul<<<NIMAGE_IN_TEXTURE, BLOCK_SIZE>>>(d_ccf+ccf_base_addr+k*NIMAGE_IN_TEXTURE*(RING_LENGTH+2), BLOCK_SIZE, NRING, NIMAGE, KX, KY);
		}
		if (NIMAGE_LEFT_TEXTURE != 0) {
			cudaBindTexture(0, texPtr, d_subject_image_polar+NTEXTURE*NIMAGE_IN_TEXTURE*(RING_LENGTH+2)*NRING, &channelDesc, NIMAGE_LEFT_TEXTURE*(RING_LENGTH+2)*NRING*sizeof(float));
			complex_mul<<<NIMAGE_LEFT_TEXTURE, BLOCK_SIZE>>>(d_ccf+ccf_base_addr+NTEXTURE*NIMAGE_IN_TEXTURE*(RING_LENGTH+2), BLOCK_SIZE, NRING, NIMAGE, KX, KY);
		}
		cudaUnbindTexture(texim_subject);
		cudaUnbindTexture(texim_ref);
	}
    }
    cudaUnbindTexture(texim_points);

    for (i=0; i<IMAGE_BATCH2; i++)
	    cufftExecC2R(plan_ccf, (cufftComplex *)(d_ccf+i*POINTS_PER_IMAGE*IMAGE_PER_BATCH2*(RING_LENGTH+2)), (cufftReal *)(d_ccf+i*POINTS_PER_IMAGE*IMAGE_PER_BATCH2*(RING_LENGTH+2)));
    if (IMAGE_LEFT_BATCH2 != 0) 
	    cufftExecC2R(plan_ccf_rest, (cufftComplex *)(d_ccf+IMAGE_BATCH2*POINTS_PER_IMAGE*IMAGE_PER_BATCH2*(RING_LENGTH+2)), (cufftReal *)(d_ccf+IMAGE_BATCH2*POINTS_PER_IMAGE*IMAGE_PER_BATCH2*(RING_LENGTH+2)));


    cudaMemcpy(ccf, d_ccf, (RING_LENGTH+2)*NIMAGE*POINTS_PER_IMAGE*sizeof(float), cudaMemcpyDeviceToHost);

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
    
    free(points);
    free(shifts_ref);
    free(shifts_subject);
    
    return;
}

void filter_image(float *image, float *filter_image, int NIMA, int NX, int NY, float *params) {
	
	float *image_padded, *d_image_padded;
	
	cufftHandle plan_R2C, plan_C2R;
	cufftPlan2d(&plan_R2C, NX*2, NY*2, CUFFT_R2C);
	cufftPlan2d(&plan_C2R, NX*2, NY*2, CUFFT_C2R);
	
	int padded_size = (NX*2+2)*(NY*2);
	
	image_padded = (float *)malloc(padded_size*NIMA*sizeof(float));
	for (int i=0; i<NIMA*padded_size; i++)   image_padded[i] = 0.0f;
	
	cudaMalloc((void **)&d_image_padded, padded_size*NIMA*sizeof(float));

	for (int im=0; im<NIMA; im++) 
		for (int iy=0; iy<NY; iy++)
			memcpy(image_padded+im*padded_size+(iy+NY/2)*(NX*2+2)+NX/2, image+im*NX*NY+iy*NX, NX*sizeof(float));

	cudaMemcpy(d_image_padded, image_padded, padded_size*NIMA*sizeof(float), cudaMemcpyHostToDevice);
	
	for (int im=0; im<NIMA; im++) {
		int base_address = im*padded_size;
		cufftExecR2C(plan_R2C, (cufftReal *)(d_image_padded+base_address), (cufftComplex *)(d_image_padded+base_address));
		mul_ctf<<<NX+1, NY*2>>>(d_image_padded+base_address, NX*2, NY*2, params[im*6], params[im*6+1], params[im*6+2], params[im*6+3], params[im*6+4], params[im*6+5]);
		cufftExecC2R(plan_C2R, (cufftComplex *)(d_image_padded+base_address), (cufftReal *)(d_image_padded+base_address));
	}
	
	cudaMemcpy(image_padded, d_image_padded, padded_size*NIMA*sizeof(float), cudaMemcpyDeviceToHost);

	for (int im=0; im<NIMA; im++)
		for (int iy=0; iy<NY; iy++)
			memcpy(filter_image+im*NX*NY+iy*NX, image_padded+im*padded_size+(iy+NY/2)*(NX*2+2)+NX/2, NX*sizeof(float));	
		
	cufftDestroy(plan_R2C);
	cufftDestroy(plan_C2R);
	cudaFree(d_image_padded);
	
	free(image_padded);

	return;
}

__global__ void mul_ctf(float *image, int nx, int ny, float defocus, float cs, float voltage, float apix, float bfactor, float ampcont) {

    // Block index
    int bx = blockIdx.x;

    // Thread index
    int tx = threadIdx.x;

    float x, y;
    
    x = float(bx);
    if (tx >= ny/2) y = float(tx-ny);
    else y = float(tx);

    float ak = sqrt(x*x+y*y)/nx/apix;
    float cst = cs*1.0e7f;
    float wgh = ampcont/100.0;
    float phase = atan(wgh/sqrt(1.0f-wgh*wgh));
    float lambda = 12.398f/sqrt(voltage*(1022.f+voltage));
    float ak2 = ak*ak;
    float g1 = defocus*1.0e4f*lambda*ak2;
    float g2 = cst*lambda*lambda*lambda*ak2*ak2/2.0f;
    float ctfv = static_cast<float>(sin(PI*(g1-g2)+phase));

    if (bfactor != 0.0f)  ctfv *= exp(-bfactor*ak2/4.0f);
    image[(bx+tx*(nx/2+1))*2] *= ctfv;
    image[(bx+tx*(nx/2+1))*2+1] *= ctfv;
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

__global__ void resample_to_polar(float* image, float sx, float sy, int NX, int NY, int RING_LENGTH, int NRING) {

    // Block index
    int bx = blockIdx.x;
    int by = blockIdx.y;
    
    // Thread index
    int tx = threadIdx.x;

    float cnx =       NX/2+sx+0.5;
    float cny = by*NY+NY/2+sy+0.5;

    int img_index = (by*NRING+bx)*(RING_LENGTH+2)+tx;
    int index = bx*RING_LENGTH+tx;
    float points_x = tex1Dfetch(texim_points, index*2);
    float points_y = tex1Dfetch(texim_points, index*2+1);
    float shifts_x = tex1Dfetch(texim_shifts, by*2);
    float shifts_y = tex1Dfetch(texim_shifts, by*2+1);

    image[img_index] = tex2D(tex, cnx+points_x+shifts_x, cny+points_y+shifts_y);
}

