#ifndef ALI3D_D_MPI_H
#define ALI3D_D_MPI_H

#include "mpi.h"
#include "emdata.h"
#include "projector.h"

#include "alignoptions.h"

#ifndef PI
#define PI 3.141592653589793
#endif

using namespace EMAN;

std::vector<int> Numrinit(int first_ring, int last_ring, int skip = 1, std::string mode = "F");

std::vector<float> ringwe(std::vector<int> numr, std::string mode = "F");

int Applyws(EMData * circ, std::vector<int> numr, std::vector<float> wr);

EMData * recons3d_4nn(std::string stack_name, std::vector<int> list_proj, int npad = 4);

int ali3d_d(MPI_Comm comm, EMData*& volume, EMData** projdata, EMData** cleandata,
            float *angleshift, int nloc, AlignOptions& options, char* fname_base);

#endif // ALI3D_D_MPI_H
