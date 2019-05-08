/*
 * Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
 * Please do not copy or modify this file without written consent of the author.
 * Copyright (c) 2000-2019 The University of Texas - Houston Medical School
 *
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */
/* This header file includes all the function delarations needed to run HyBR in C */

#define integer int
#define frand() ((float) rand() / (RAND_MAX+1.0))

/* from matrix.c */
extern float **matrix(int m, int n);
extern float **submatrix(float **A, int m, int n, int newm, int newn);
extern void free_matrix(float **A);
extern void fill_rand(float **A, int m, int n);
extern void print_matrix(float **A, int m, int n);
extern void transpose_matrix(float **A, int m, int n);
extern void column_orient(float **A, int m, int n, float ** B);
extern void row_sum(float **A, int m, int n, float *e);
extern void copy_matrix(float **A, float **B, int m, int n);

extern float *myvector(int n);
extern float *subvector(float *x, int n, int newn);
extern void free_vector(float * x);
extern void print_vector(float *x, int n);
extern void fill_ones(float *x, int n);
extern float sum_vector(float *x, int n);
extern void copy_vector(float *x, float *y, int n);
extern void fillveczeros(float *x, int n);
extern float norm(float *x, int n);
extern void alphaxpy(float alpha, float *x, float *y, int n);

extern void matvec_mult(float ** A, float *x, int m, int n, float *b);
extern void matvec_multT(float **A, float *x, int m, int n, float *b);
extern void get_col(float ** A, int m, int k, float * ak);
extern void get_row(float ** A, int n, int k, float * ak);
extern void put_col(float ** A, int m, float * v, int k);


/* from direct_methods.c */
extern void svd(float **A, int m, int n, float **U, float *x, float **V);
extern void tsvd(float **U, float *s, float **V, float *b, int m, int n, float omega, float *x, float *tol);
extern void gcv_tsvd(float *s, float *bhat, int m, int n, float omega, float * tol);
extern void tikhonov(float **U, float *s, float **V, float *b, int m, int n, float omega, float *x, float * alpha);
extern float GCVmin_Tik(float *s, float *bhat, int m, int n, int MaxIts, float omega);
extern float GCVfun_Tik(float alpha, float *s, float * bhat, int m, int n, float omega);
extern void generate_noise(double *E, int m, int n);



