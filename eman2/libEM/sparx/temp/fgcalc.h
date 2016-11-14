#define NUMBER float

int fgcalc(MPI_Comm comm, float *volsph, Vec3i volsize, int nnz, 
           int nrays    , Vec3i  origin, int        ri, int *ptrs, 
           int *cord    , float *angtrs, int      nang, float *rhs, 
           float     aba, NUMBER  *fval, float   *grad, char *fname_base);

int fcalc(float *volsph, Vec3i volsize, 
           int nnz, int nrays, Vec3i origin, int ri, 
           int *ptrs, int *cord, float *angtrs, int nang, 
	  float *rhs, NUMBER *fval);
