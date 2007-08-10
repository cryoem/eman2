#ifndef UTILCOMM_CART_H
#define UTILCOMM_CART_H

#include "mpi.h"
#include "stdlib.h"
#include "emdata.h"
#include "alignoptions.h"

using namespace EMAN;

int ReadStackandDist_Cart(MPI_Comm comm, EMData ***images2D, char *stackfname, int *nloc);
int CleanStack_Cart(MPI_Comm comm, EMData ** image_stack, int nloc, int ri, Vec3i volsize, Vec3i origin);
int ReadAngTrandDist_Cart(MPI_Comm comm_2d, MPI_Comm comm_row, int *dim, float *angleshift, char *angfname, int nloc);
int setpart_gc1(MPI_Comm comm_2d, int nangs, int *psize, int *nbase);
int setpart_gr1(MPI_Comm comm_2d, int nnz, int *nnzpart, int *nnzbase);

int sphpart(MPI_Comm comm_2d, int nrays, int *ptrs, int *nnzbase, int *ptrstart);
int getcb2sph(Vec3i volsize, int ri, Vec3i origin, int nnz0, int *ptrs, int *cord);

int fwdpj3_Cart(Vec3i volsize, int nraysloc, int nnzloc, float *dm, Vec3i origin, int ri, int *ptrs, int *cord, int myptrstart, float *x, float *y);

int bckpj3_Cart(Vec3i volsize, int nraysloc, int nnzloc, float *dm, Vec3i origin, int ri, int *ptrs, int *cord, int myptrstart, float *x, float *y);


#endif // UTILCOMM_CART_H

