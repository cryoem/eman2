#ifndef UTILCOMM_H
#define UTILCOMM_H

#include "mpi.h"
#include "stdlib.h"
#include "emdata.h"
#include "alignoptions.h"

using namespace EMAN;

int ReadVandBcast(MPI_Comm comm, EMData *volume, char *volfname);
int ReadStackandDist(MPI_Comm comm, EMData ***images2D, char *stackfname, int *nloc);
int ReadAngTrandDist(MPI_Comm comm, float *angleshift, char *paramfname, int nloc);
int CleanStack(MPI_Comm comm, EMData ** image_stack, int nloc, int ri, Vec3i volsize, Vec3i origin);
int setpart(MPI_Comm comm, int nima, int *psize, int *nbase);
int ParseAlignOptions(MPI_Comm comm, AlignOptions& options, char* optionsfname, int nvoxels, EMData*& mask3D);
int asta2(float *img, int nx, int ny, int ri, double *abaloc, int *klploc);
int setpart(MPI_Comm comm, int nang, int *psize, int *nbase);

#endif // UTILCOMM_H
