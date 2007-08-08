#ifndef ALI3D_UNIFIED_MPI_H
#define ALI3D_UNIFIED_MPI_H
#include "mpi.h"

#define NUMBER float

#ifndef PI
#define PI 3.141592653589793
#endif

using namespace EMAN;
//using namespace std;

int unified(MPI_Comm comm, EMData * volume, EMData **projdata, 
            float *angleshift, int nang, int max_iter, char *fname_base);


#endif // ALI3D_UNIFIED_MPI_H
