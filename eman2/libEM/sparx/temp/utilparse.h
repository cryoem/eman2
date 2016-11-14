#ifndef UTILPARSE_H
#define UTILPARSE_H

#include "mpi.h"
#include "emdata.h"
#include "alignoptions.h"

using namespace EMAN;

int ParseAlignOptions(MPI_Comm comm, AlignOptions& options, char* optionsfname, int nvoxels, EMData*& mask3D);

#endif // UTILPARSE_H
