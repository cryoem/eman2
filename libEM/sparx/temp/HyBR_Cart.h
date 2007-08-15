#ifndef HYBR_CART_H
#define HYBR_CART_H
#include "mpi.h"

#ifndef PI
#define PI 3.141592653589793
#endif

#define frand() ((float) rand() / (RAND_MAX+1.0))

//#define max(a,b) a > b ? a : b
//#define min(a,b) a < b ? a : b

#include "mpi.h"
#include "emdata.h"

using namespace EMAN;

int recons3d_HyBR_mpi_Cart(MPI_Comm comm_2d, MPI_Comm comm_row, MPI_Comm comm_col, EMData ** images, float * angleshift, EMData *& xvol, int nangloc, int radius = -1, int maxit = 100, std::string symmetry = "c1", int insolve = 1);
void LBD_Cart(MPI_Comm comm_2d, MPI_Comm comm_row, MPI_Comm comm_col, float *angleshift, Vec3i volsize, int nraysloc, int nnzloc, Vec3i origin, int radius, int *ptrs, int myptrstart, int *cord, int nangloc, int nx, int ny, float *Ukloc, float **B, float **Vloc,int  m, int n, int k, int maxiter, std::string symmetry);
float findomega(float * bhat, float *s, int m, int n, int insolve);
float gcvstopfun(float alpha, float * u, float *s, float beta, int k, int m, int n, int insolve);

void recons3d_CGLS_mpi_Cart(MPI_Comm comm_2d, MPI_Comm comm_row, MPI_Comm comm_col, EMData ** images, float * angleshift, EMData *& xvol, int nangloc, int radius, int maxiter, std::string symmetry);

#endif // HYBR_CART_H
