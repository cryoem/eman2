#ifndef SIRT_H
#define SIRT_H
#include "mpi.h"

#ifndef PI
#define PI 3.141592653589793
#endif

int recons3d_sirt_mpi(MPI_Comm comm, EMData ** images, float * angleshift, EMData *& xvol, int nangloc, int radius = -1, float lam = 1.0e-4, int maxit = 100, std::string symmetry = "c1", float tol = 1.0e-3);

#endif // SIRT_H
