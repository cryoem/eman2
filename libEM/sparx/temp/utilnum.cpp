#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <float.h>

#include "lapackblas.h"
#include "utilnum.h"

int asta2(float *img, int nx, int ny, int ri, double *abaloc, int *klploc)
{
    // calculate the sum of the background pixels, and then set the intensity
    // of these pixels to zero. Also count the number of background pixels.
    // A background pixel is a pixel outside of the circle with radius ri
    // This is done for images assigned to a processor
    //
    int xcent = (nx / 2) + 1;
    int ycent = (ny / 2) + 1;
    int r_squared = ri*ri;

    int x_summand, y_summand;

    for ( int i = 0 ; i < nx ; ++i ) {
	x_summand = (i-xcent) * (i-xcent);
	for ( int j = 0 ; j < ny ; ++j ) {
	    y_summand = (j-ycent) * (j-ycent);
	    if ( x_summand + y_summand > r_squared ) {
		*abaloc += (double) img[j*nx + i];
		//chao set the background to zero
		img[j*nx+i]=0.0;
		++*klploc;
	    }
	}
    }

    return 0;
}

int ifix(float a)
{
    int ia = 0;
    if (a >= 0.0) {
       ia = (int)floor(a);
    }
    else {
       ia = (int)ceil(a);
    }
    return ia;
}

// ================= imported from direct_methods.c =================

/* direct_methods.c
 * ------------------------------------------------------------------
 *  This file contains SVD based methods to solve discrete
 *  ill-posed problems.
 *
 *------------------------------------------------------------------
 */


/*-------SVD--------------------------------------------------------*/
void svd(float **A, int m, int n, float **U, float *s, float **V)
/*
 * Computes the SVD of an m-by-n matrix A:
 *
 *        A = USV^T
 *
 *  NOTE: We assume m >= n
 *
 *  Input:  A    : m-by-n matrix
 *          m, n : row and column dimensions of A
 *
 *  Output: U    : Orthogonal matrix containing left singular vectors.
 *          s    : vector containing singular values of A.
 *               : That is, s = diag(S)
 *          V    : Orthogonal matrix containing right singular vectors.
 */
 
{ 
  char jobu = 'A';
  char jobv = 'A';
  
  int i;

// chao something funny about min/max here... hardwire lwork to 100 for now  
//  integer lwork = max(3*min(m,n) + max(m,n), 5*min(m,n));
  integer lwork = 100;
  integer lda=m; 
  integer ldu=m; 
  integer ldv=n;
  integer info;
  integer mm = m;
  integer nn = n;  

  float ** A_column = matrix(n,m);
  float * work;
  
  work = myvector((int) lwork);

  /* We have to column orient the matrix going in because dgesvd expects the 
   * matrix to be column oriented and in C, the matrix is row oriented. */
 column_orient(A, m, n, A_column);
  
 /* The LAPACK routines destroy the input matrix, so use a
  * copy of it.
  *
  * Check to see if we are using single or float precision,
  * and use appropriate LAPACK routine to compute SVD.
  */

 sgesvd_(&jobu, &jobv, &mm, &nn, *A_column, &lda, s, *U, &ldu, *V, &ldv, 
         work,  &lwork, &info);

//sgesvd_(jobu, jobv, mm, nn, *A_column, lda, s, *U, ldu, *V, ldv, &info);

  /* We have to transpose the U matrix to make it row oriented again. 
   * We don't worry about V because we would have had to transpose it anyways. */
	transpose_matrix(U, m, m);
  
  free_vector(work);
  free_matrix(A_column);

  
}

/*-------TSVD--------------------------------------------------------*/
void tsvd(float **U, float *s, float **V, float *b, int m, int n, float omega, float *x, float *tol)
   /*
    *  Computes the TSVD solution
    *
    *  NOTE: We assume m >= n
    *
    *  Input:  U    : Orthogonal matrix containing left singular vectors.
    *          s    : vector containing singular values of A.
    *               : That is, s = diag(S)
    *          V    : Orthogonal matrix containing right singular vectors.
    *          b    : Right hand side vector
    *          m, n : row and column dimensions of A.
    *          tol  : Truncation tolerance.  If tol < 0, then
    *                   GCV is used to choose it.
    *
    *  Output: x    : TSVD solution
    *          tol  : Truncation tolerance.
    */
{    
    int i, k;
    float * sinv = myvector(n);
    float * bhat = myvector(m);

    matvec_multT(U, b, m, m, bhat);

//printf("m=%d, n=%d\n", m,n);
     /*
      * Use GCV to choose tol.
      */
    k = 0;
    if (*tol < (float) 0.0)
      gcv_tsvd(s, bhat, m, n, omega, tol);

    /*
     * Now we need to compute the the filtered inverse of the
     * singular values.  That is, 1/s(i) if s(i) >= tol, and
     * 0 otherwise.  I'm not sure this is the best way to do this
     * in Fortran -- it would be better to not have to check each
     * time through the loop.  But, I do not want to assume the
     * the singular values are ordered -- they may not be in 
     * some situations. Maybe a pre-sort might be helpful???
     */
 
    for (i=0; i<n; i++){
      if (fabs(s[i]) >= *tol){
        sinv[i] = (float) 1.0 / s[i];
        k++;
      }
    }
*tol = k;
//printf("truncindex = %f\n", *tol);
//print_vector(s,n);

    /*
     * Now the TSVD solution is easy to compute.
     */
    for (i=0; i<n; i++)
      bhat[i] = sinv[i]*bhat[i];
    
    matvec_mult(V, bhat, n, n, x);
    
    
    free_vector(sinv);
    free_vector(bhat);
}
    
    
    
/*-----GCV_TSVD------------------------------------------------*/

void gcv_tsvd(float *s, float *bhat, int m, int n, float omega, float * tol)
   /*
    *  This function uses the GCV function to choose a default truncation tolerance for TSVD regularization
    *
    *  Use GCV to compute default regularization parameter.
    *
    * Input:    s : Vector containing singular values.
    *        bhat : Vector containing spectral coefficients of right
    *                 hand side vector.  That is, bhat = U^T*b.
    *
    * Output: tol : Scalar truncation tolerance for TSVD.
    */
{
    float * bhat_copy = myvector(m);
    float * rho = myvector(m-1);
    float * G = myvector(m-1);
    int i, kk=0;
    float ni;
    /*
     *  NEED TO: Make sure singular values are sorted from
     *           largest to smallest, and do the same for the
     *           bhat values.  In MATLAB this looks like:
     *              [s, idx] = sort(abs(s)); 
     *              s = flipud(s); 
     *              idx = flipud(idx);
     *              bhat = abs( bhat(idx) );
     *
     * Initialize residual and GCV function values.
     */
    
    copy_vector(bhat, bhat_copy, m);
    for (i=0; i<m; i++)
      bhat_copy[i] = fabs( bhat[i] );
    
    rho[m-2] = bhat_copy[m-1]*bhat_copy[m-1];

    /*
     * Recursively compute residual and GCV function values.
     */
    for (i = m-3; i>-1; i--){
      rho[i] = rho[i+1] + bhat_copy[i+1]*bhat_copy[i+1];
    }

//    print_vector(bhat_copy, m-1);
    for(i=0; i< m-1; i++){
      ni = (float) m - omega*(i+1);
      G[i] = (float) n* rho[i]/(ni*ni);
    }
 // print_vector(G,m-1);
    /*
     * Now find minimum of the discrete GCV function.
     */
    for (i=0; i < m-2; i++){
      if (G[i+1]<G[i])
        kk = i+1;
    }

    *tol = s[kk];
 // printf("kk = %d, s[kk] = %f\n", kk, *tol);

    free_vector(bhat_copy);
    free_vector(rho);
    free_vector(G);
}


/*------TIKHONOV------------------------------------------------------------*/

void tikhonov(float **U, float *s, float **V, float *b, int m, int n, float omega, float *x, float * alpha)
  /*
   *  Computes the Tikhonov solution
   *
   *  NOTE: We assume m >= n
   *
   *  Input:  U    : Orthogonal matrix containing left singular vectors.
   *          s    : vector containing singular values of A.
   *               : That is, s = diag(S)
   *          V    : Orthogonal matrix containing right singular vectors.
   *          b    : Right hand side vector
   *          m, n : row and column dimensions of A.
   *        alpha  : Regularization parameter.  If alpha < 0, then
   *                   GCV is used to choose it.
   *        omega  : Omega parameter
   *
   *  Output: x    : Tikhonov solution
   *        alpha  : Regularization parameter
   */
{
    int i, MaxIts = 500;
    float * bhat = myvector(m);
    float * xhat = myvector(n);
    float * D = myvector(n);

    matvec_multT(U, b, m, m, bhat);

/*printf("b:\n");
    print_vector(b, m);
    printf("m = %d, n=%d, omega=%f\n", m,n,omega);
printf("bhat:\n");
    print_vector(bhat, m);
printf("s:\n");
    print_vector(s, n);
*/


     /*
      * Use GCV to choose alpha.
      */
    
    if (*alpha < (float)0.0)
       *alpha = GCVmin_Tik(s, bhat, m, n, MaxIts, omega);

   /*
    * Now we need to compute the the filtered Tikhonov solution
    */
    for (i=0; i<n;i++){
     D[i] = fabs(s[i])*fabs(s[i]) + (*alpha)*(*alpha);
     xhat[i] = (s[i] * bhat[i])/D[i];
    }
 
    matvec_mult(V, xhat, n, n, x);

    free_vector(bhat);
    free_vector(xhat);
    free_vector(D);
}

/*------GCVmin_Tik------------------------------------------------------------*/

float GCVmin_Tik(float *s, float *bhat, int m, int n, int MaxIts, float omega)
   /*
    *
    *  Use a combination of Golden Section Search and Parabolic
    *  Interpolation to compute the minimum of the GCV 
    *  function to find a default regularization parameter (see MATLAB's fminbnd)
    *
    * Input:    s : Vector containing singular or spectral values of A.
    *        bhat : Vector containing spectral coefficients of right
    *                 hand side vector.  That is, bhat = U^T*b.
    *       MaxIts : Maximum number of iterations to find the minimum
    *
    * Output: xf : minimum of the GCV function for Tikhonov.
    */
{
    float ax, bx, seps, c, a, b, v, w, sf, d, e, x, fx, fv, fw, xm, tol, tol1, tol2, r, q, p, si, fu, fval, xf;
    int funccount, iter, maxfun, gs;
    
/*float tt;

tt = GCVfun_Tik(0.01, s, bhat, m,n, omega);
printf("alpha = 0.01;  GCV(a) = %f\n", tt);
tt = GCVfun_Tik(0.59, s, bhat, m,n, omega);
printf("alpha = 0.59;  GCV(a) = %f\n", tt);
tt = GCVfun_Tik(.999, s, bhat, m,n, omega);
printf("alpha = .999;  GCV(a) = %f\n", tt);
tt = GCVfun_Tik(1.2, s, bhat, m,n, omega);
printf("alpha = 1.2;  GCV(a) = %f\n", tt);
*/

    ax = (float) 0.0;
 //   bx = (float)1.0;
   bx = s[0];

    /*  Initialize some values */
    seps = sqrt(FLT_EPSILON);
    c = (float)0.5*((float)3.0 - sqrt((float)5.0));
    a = ax;
    b = bx;
    v = a + c*(b-a);
    w = v;
    xf = v;
    d = (float)0.0;
    e = (float)0.0;
    x = xf; 

    fx = GCVfun_Tik(x, s, bhat, m,n, omega);

    fv = fx;
    fw = fx;
    xm = (float)0.5*(a+b);
    tol = (float).0001;
    tol1 = seps*fabs(xf) + tol/(float)3.0;
    tol2 = (float)2.0*tol1;
    funccount = 0;
    iter = 0;
    maxfun = 500;

    for (iter = 1; iter< MaxIts; iter++){
      gs = 1;
      /* Is a parabolic fit possible*/
      if (fabs(e) > tol1){
        gs = 0;
        r = (xf-w)*(fx-fv);
        q = (xf-v)*(fx-fw);
        p = (xf-v)*q-(xf-w)*r;
        q = (float)2.0*(q-r);
        if (q > (float)0.0)
          p = -p;
        q = fabs(q);
        r = e;
        e = d;
        /*Is the parabola acceptable*/
        if ( (fabs(p)<fabs((float)0.5*q*r)) && (p>q*(a-xf)) && (p<q*(b-xf)) ){
          /* Yes, parabolic interpolation step*/
          d = p/q;
          x = xf+d;
          if (((x-a) < tol2) || ((b-x) < tol2)){
            if(xm-xf <0)
              si = -(float)1.0;
            else
              si = (float)1.0;
            
            d = tol1*si;
          }
        }
        else{
          /* Not acceptable, must do a golden section step*/
          gs = 1;
        }
      }
      if (gs){
        /* A golden-section step is required*/
        if (xf >= xm)
          e = a-xf;
        else
          e = b-xf;

        d = c*e;
      }
      /* The function must not be evaluated too close to xf*/
     if(d <0)
       si = -(float)1.0;
     else
       si = (float)1.0;
     
      x = xf + si * max( fabs(d), tol1 );
      fu = GCVfun_Tik( x, s, bhat, m,n, omega);
      funccount = funccount + 1;
      if (fu <= fx){
        if (x >= xf)
          a = xf;
        else
          b = xf;
   
        v = w;
        fv = fw;
        w = xf;
        fw = fx;
        xf = x;
        fx = fu;
    }
      else{
        if (x < xf)
          a = x;
        else 
          b = x;
        
        if ( (fu <= fw) || (w == xf) ){
          v = w;
          fv = fw;
          w = x;
          fw = fu;
        }
        else if ((fu <= fv) || (v == xf) || (v == w)){
          v = x;
          fv = fu;
        }
      }
      xm = (float)0.5*(a+b);
      tol1 = seps*fabs(xf) + tol/(float)3.0;
      tol2 = (float)2.0*tol1;
      if (funccount >= maxfun || fabs(xf-xm) <= (tol2 - (float)0.5*(b-a))){
        fval = fx;
        return xf;
      }
}
    fval = fx;
    return xf;
}


/*---------GCVFUN_TIK---------------------------------------------------------*/

float GCVfun_Tik(float alpha, float * s, float * bhat, int m, int n, float omega)
    /*
     *   This function evaluates the GCV function for Tikhonov regularization.
     *   Assume the problem has the form b = Ax + noise, with [U, S, V] = svd(A).
     *
     * Input: alpha: value at which to evaluate the GCV function
     *           s : Vector containing singular or spectral values of A.
     *        bhat : Vector containing spectral coefficients of right
     *                 hand side vector.  That is, bhat = U^T*b.
     *
     * Output: G : GCV function for Tikhonov at alpha
     */
{
    float t0, G;
    float * hold = myvector(m-n);
    float * tt= myvector(n);
    float * t1= myvector(n);
    float * t2= myvector(n);
    float * t3= myvector(n);
    int i;
    
    for (i=0; i< m-n; i++)
      hold[i] = fabs(bhat[n])*fabs(bhat[n]);
    t0 = sum_vector(hold,m-n);

    for (i = 0; i<n; i++){
      tt[i] = (float)1.0 /(s[i]*s[i] +alpha*alpha);
      t1[i] = alpha*alpha * tt[i];
      t2[i] = fabs(bhat[i] * t1[i]) * fabs(bhat[i] * t1[i]);
      t3[i] = t1[i] + ((float)1.0 - omega)* (s[i]*s[i]) * tt[i];
    }
    G = (float) n*(sum_vector(t2,n) + t0) / ((sum_vector(t3,n) + (float)m - (float)n)* (sum_vector(t3,n) + (float)m - (float)n));

//printf(" svec = %f, to = %f, Num = %f; Dem = %f \n", sum_vector(t2,n), t0, (float) n*(sum_vector(t2,n) + t0), ((sum_vector(t3,n) + (float)m - (float)n)* (sum_vector(t3,n) + (float)m - (float)n)));

    free_vector(hold);
    free_vector(tt);
    free_vector(t1);
    free_vector(t2);
    free_vector(t3);
    
    return G;
}

/* --------GENERATE_NOISE CODE--------------------------------------------*/

void generate_noise(float *E, int m, int n)
   /*
    * Generate random noise.
    *
    * Input:  
    *       m,n : Dimensions of E (needs to be a vector!)
    *
    * Output: E : Array containing random noise, with norm = 1.
    */
{
    float **U = matrix(m,12);
    float nm, a = -6;
    int i, sz = 1, j;
    float *o = myvector(m);

    fill_ones(o,m);
       
    fill_rand(U, m, 12);
    row_sum(U, m, 12, E);
    
    for(j=0; j<m; j++)
      E[j] = E[j] + a*o[j];

    //daxpy_(&m, &a, o, &sz, E, &sz);    
    //norm = dnrm2_(&m, E, &sz);
    nm = norm(E, m);

    for(i = 0; i < m; i++)
      E[i] = E[i] / nm;
    
    free_matrix(U);
    free_vector(o);

}


// ======================= imported from matrix.cpp =================

/* matrix.c contains all the functions that pertain to matrices and vectors
 *
 * The following are MATRIX functions:
 *
 * 	float **matrix(int m, int n);
 *  float **submatrix(float **A, int m, int n, int newm, int newn);
 *  void free_matrix(float **A);
 *  void fill_rand(float **A, int m, int n);
 *  void print_matrix(float **A, int m, int n);
 *  void transpose_matrix(float **A, int m, int n);
 *  void row_sum(float **A, int m, int n, float *e);
 *  void copy_matrix(float **A, float **B, int m, int n);
 *
 *  The following are VECTOR functions:
 *
 *  float *myvector(int n);
 *  float *subvector(float *x, int n, int newn);
 *  void free_vector(float * x);
 *  void print_vector(float *x, int n);
 *  void fill_ones(float *x, int n);
 *  float sum_vector(float *x, int n);
 *  void copy_vector(float *x, float *y, int n);
 *  void fillveczeros(float *x, int n);
 *  float norm(float *x, int n);
 *  void alphaxpy(float alpha, float *x, float *y, int n)
 *
 *  The following are MATRIX-VECTOR functions:
 *
 *  void matvec_mult(float **A, float *x, int m, int n, float *b);
 *  void matvec_multT(float **A, float *x, int m, int n, float *b);
 *  float * get_col(float ** A, int m, int k);
 *  float * get_row(float ** A, int n, int k);
 *  void put_col(float ** A, int m, float * v, int k);
 */


/* -------------MATRIX FUNCTIONS ------------------------------------*/

/* Allocates a float matrix A of size m x n with subscript range A[0..m-1][0..n-1] */
float **matrix(int m, int n)
{
  int i;
  float **A=(float **) calloc(m, sizeof(float*));
  if (A == NULL)
  {
    fprintf(stderr, "Allocation failure in matrix()");
    exit (1);
  }
  A[0]=(float *) calloc(m * n, sizeof(float));
  if (A[0] == NULL)
  {
    fprintf(stderr, "Allocation failure in matrix()");
    exit (1);
  }
  for(i=1; i < m; i++)
    A[i]=A[0] + i * n;
  return A;
}

/* Points a submatrix B[0...(newm-1)][0... (newn-1)] to A[0..(newm-1)][0..(newn-1)] */
float **submatrix(float **A, int m, int n, int newm, int newn)
{
  int i,j;
  float **B = matrix(newm, newn);
  if (newm > m || newn > n){
    printf("error in submatrix\n");

    return B;
  }
  for(i=0;i<newm;i++)
    for(j=0;j<newn;j++)
      B[i][j]=A[i][j];
 
  return B;
}

/* Frees a matrix allocated by matrix() */
void free_matrix(float **A)
{
  free(A[0]);
  free(A);
}

/* Fills a matrix A[0..m-1][0..n-1] with random numbers. */
void fill_rand(float **A, int m, int n)
{
  int i,j;
  for(i = 0; i < m; i++) {
    for (j = 0; j < n; j++)
      A[i][j]=frand();
    
  }
}

/* Prints a matrix A[0..m-1][0..n-1] */
void print_matrix(float **A, int m, int n)
{
  int i,j;
  for (i=0; i< m; i++)
  {
    for (j=0; j< n; j++)
      printf ("%.15f   ", A[i][j]);
      putchar ('\n');
  }
}

/* Does a transpose of matrix A and stores it in matrix A*/
void transpose_matrix(float **A, int m, int n)
{
  float save;
  int i,j;
   for (i = 0; i < m; i++){
     for (j = i+1; j < n; j++)
      {   
        save = A[i][j];
        A[i][j] = A[j][i];
        A[j][i] = save;
      }
   }
}


/* Rearranges matrix A to be in column-major form and stores it in matrix B*/
void column_orient(float **A, int m, int n, float ** B)
{
  int i,j;
   for (j = 0; j < n; j++){
     for (i = 0; i < m; i++)
      {   
        B[j][i] = A[i][j];
      }
   }
}

/* Computes the row sums of matrix a and stores it in vector e*/
void row_sum(float **A, int m, int n, float *e)
{
  int i,j;
   for (i = 0; i < m; i++){
     for (j = 0; j < n; j++)
      {   
        e[i] = e[i]+ A[i][j];
      }
   }
}

/* Copies matrix A into matrix B. */
void copy_matrix(float **A, float **B, int m, int n)
{
  int i,j;
  for(i = 0; i < m; i++) {
    for (j = 0; j < n; j++)
      B[i][j]=A[i][j];
    
  }
}


/* -------VECTOR FUNCTIONS ---------------------------------------*/

/* Allocates a float vector with subscript range x[0..n-1] */
float *myvector(int n)
{
  float *x=(float *) calloc(n, sizeof(float));
  if (x == NULL)
  {
    fprintf(stderr, "Allocation failure in matrix()");
    exit (1);
  }
  return x;
}

/* Points a subvector y[0...(newn-1)] to x[0..(newn-1)] */
float *subvector(float *x, int n, int newn)
{
  int i;
  float *y = myvector(newn);
  if (newn > n){
    printf("error in subvector");
    return y;
  }
  for(i=0;i<newn;i++)
      y[i]=x[i];
  
  return y;
}

/* Frees a vector allocated by vector() */
void free_vector(float *x)
{
  free(x);
}

/* Prints a vector x[0..n-1] */
void print_vector(float *x, int n)
{
  int i;
  for (i=0; i< n; i++)
  {
      printf ("%.15f  \n ", x[i]);
  }
}

/* Fills a vector with all ones*/
void fill_ones(float *x, int m)
{
  int i;
   for (i = 0; i < m; i++){
     x[i] = (float) 1.0;
   }
}

/* sums a vector*/
float sum_vector(float *x, int n)
{
  int i;
  float sum = 0.0;
   for (i = 0; i < n; i++){
     sum = sum + x[i];
   }
  return sum;
}

/* Copies vector x into vector y. */
void copy_vector(float *x, float *y, int n)
{
  int i;
  for(i = 0; i < n; i++) 
      y[i]=x[i];    
}

/* Fills a vector x with all zeros */
void fillveczeros(float *x, int n)
{
 int i;
 for (i=0; i<n; i++)
   x[i] = 0.0;
}

/* Computes the Euclidean norm of vector x */
float norm(float *x, int n)
{
  float x_norm = 0.0;
  int i;
  for (i=0; i<n; i++)
	x_norm = x_norm + x[i]*x[i];

  x_norm = sqrt(x_norm);

//x_norm = cblas_dnrm2(n, x, 1);

  return x_norm;
}

/* Computes a daxpy: y = alpha*x + y  */
void alphaxpy(float alpha, float *x, float *y, int n)
{ 
 int i;
 for (i=0; i<n; i++)
   y[i] = alpha*x[i] + y[i];
}

/* ------MATRIX-VECTOR FUNCTIONS ------------------------------------*/

/* Does a matrix vector multiply: b = A*x */
void matvec_mult(float **A, float *x, int m, int n, float *b)
{
int i,j;
for(i=0; i< m; i++){
	b[i] = 0.0;
	for(j=0; j<n; j++){
		b[i] = b[i] + A[i][j]*x[j];
}
}

//cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, *A, n, x, 1, 0.0, b, 1);
}

/* Does a matrix-TRANSPOSE vector multiply: b = A'*x */
void matvec_multT(float **A, float *x, int m, int n, float *b)
{
	int i,j;
for(i=0; i< n; i++){
	b[i] = 0.0;
	for(j=0; j<m; j++){
		b[i] = b[i] + A[j][i]*x[j];
}
}

//cblas_dgemv(CblasRowMajor, CblasTrans, m, n, 1.0, *A, m, x, 1, 0.0, b, 1);
}

/* The following function gets the kth column of A and 
 * copies it into the vector ak, which is returned.
 */
void get_col(float ** A, int m, int k, float * ak)
{
  int i;
  for (i = 0; i < m; i++)
    ak[i] = A[i][k];

}

/* The following function gets the kth row of A and 
 * copies it into the vector ak, which is returned.
 */
void get_row(float ** A, int n, int k, float * ak)
{
  int i;
  for (i = 0; i < n; i++)
    ak[i] = A[k][i];

}

/* The following function puts vector v into the kth column of A*/
void put_col(float ** A, int m, float * v, int k)
{
  int i;
  for (i = 0; i < m; i++)
    A[i][k] = v[i];

}

