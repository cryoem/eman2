#ifndef ALI3D_UNIFIED_MPI_H
#define ALI3D_UNIFIED_MPI_H
#include "mpi.h"

#define NUMBER float

#ifndef PI
#define PI 3.141592653589793
#endif

using namespace EMAN;
using namespace std;

int unified(MPI_Comm comm, EMData * volume, EMData **projdata, 
            float *angleshift, int nang, int max_iter, char *fname_base);

int fgcalc(MPI_Comm comm, float *volsph, Vec3i volsize, int nnz, 
           int nrays    , Vec3i  origin, int        ri, int *ptrs, 
           int *cord    , float *angtrs, int      nang, float *rhs, 
           float     aba, NUMBER  *fval, float   *grad, char *fname_base);
int fcalc(float *volsph, Vec3i volsize, 
           int nnz, int nrays, Vec3i origin, int ri, 
           int *ptrs, int *cord, float *angtrs, int nang, 
	  float *rhs, NUMBER *fval);

int ifix(float a);
int setpart(MPI_Comm comm, int nang, int *psize, int *nbase);
int fwdpj3(Vec3i volsize, int nrays, int   nnz, float *dm, 
           Vec3i  origin, int    ri, int *ptrs, int *cord, 
           float      *x, float  *y);
int bckpj3(Vec3i volsize, int nrays, int   nnz, float *dm, 
           Vec3i  origin, int    ri, int *ptrs, int *cord, 
           float      *x, float  *y);

int sph2cb(float *sphere, Vec3i volsize, int  nrays, int    ri, 
           int      nnz0, int     *ptrs, int  *cord, float *cube);
int cb2sph(float *cube, Vec3i volsize, int    ri, Vec3i origin, 
           int    nnz0, int     *ptrs, int *cord, float *sphere);
int getnnz(Vec3i volsize, int ri, Vec3i origin, int *nrays, int *nnz);

int asta2(float *img, int nx, int ny, int ri, double *abaloc, int *klploc);
int make_proj_mat(float phi, float theta, float psi, float * dm);

#endif // ALI3D_UNIFIED_MPI_H
