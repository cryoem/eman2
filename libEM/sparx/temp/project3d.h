#ifndef PROJECT3D
#define PROJECT3D

using namespace EMAN;

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
int ifix(float a);
int make_proj_mat(float phi, float theta, float psi, float * dm);

#endif // PROJECT3D
