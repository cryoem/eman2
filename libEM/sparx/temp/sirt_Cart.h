#ifndef SIRT_CART_H
#define SIRT_CART_H
#include "mpi.h"
#include "emdata.h"

using namespace EMAN;

#ifndef PI
#define PI 3.141592653589793
#endif

int recons3d_sirt_mpi_Cart(MPI_Comm comm_2d, MPI_Comm comm_row, MPI_Comm comm_col, EMData ** images, float * angleshift, EMData *& xvol, int nangloc, int radius = -1, float lam = 1.0e-4, int maxit = 100, std::string symmetry = "c1", float tol = 1.0e-3);

#endif // SIRT_CART_H
