
#ifndef PROJECT3D_CART_H
#define PROJECT3D_CART_H

#include "mpi.h"
#include "emdata.h"

int sphpart(MPI_Comm comm_2d, int nrays, int *ptrs, int *nnzbase, int *ptrstart);
int getcb2sph(Vec3i volsize, int ri, Vec3i origin, int nnz0, int *ptrs, int *cord);

int fwdpj3_Cart(Vec3i volsize, int nraysloc, int nnzloc, float *dm, Vec3i origin, int ri, int *ptrs, int *cord, int myptrstart, float *x, float *y);

int bckpj3_Cart(Vec3i volsize, int nraysloc, int nnzloc, float *dm, Vec3i origin, int ri, int *ptrs, int *cord, int myptrstart, float *x, float *y);

#endif // PROJECT3D_CART_H

