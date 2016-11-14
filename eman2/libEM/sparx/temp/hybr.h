/* This header file includes all the function delarations needed to run HyBR in C */

#define integer int
#define frand() ((float) rand() / (RAND_MAX+1.0))
//#define max(a,b) a > b ? a : b
//#define min(a,b) ((a) > (b) ? (a) : (b))

#define max(a,b) a > b ? a : b
#define min(a,b) a < b ? a : b

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



