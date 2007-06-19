#include "mpi.h"
#include "stdlib.h"

// include EMAN2
#include "emdata.h"
#include "assert.h"

#include "ali3d_d_mpi.h"
#include "ali3d_unified_mpi.h"
#include "alignoptions.h"
#include "utilcomm.h"

int main(int argc, char *argv[])
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int ncpus, mypid, ierr;
    int nloc; 
    double t0;

    MPI_Status mpistatus;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(comm,&ncpus);
    MPI_Comm_rank(comm,&mypid);
    printf("mypid = %d, ncpus = %d\n", mypid, ncpus);
    char  volfname[100], paramfname[100], stackfname[100],voutfname[100];
    EMData **expimages;

    // parse the command line and set filenames	
    if (argc < 4) {
      if (mypid == 0) {
          printf("Not enough arguments to the command...\n");
          printf("Usage: rununified -data=<imagestack> ");
          printf("-model=<initial 3D volume filename> ");
          printf("-param=<initial angles&shifts> ");
          printf("-out=<output volume filename>\n"); 
      }
      ierr = MPI_Finalize();
      exit(1);
    }
    int ia=0;
    while (ia < argc) {
       if ( !strncmp(argv[ia],"-data",5) ) {
          strcpy(stackfname,&argv[ia][6]);
          if (mypid == 0)
             printf("image stack        = %s\n", stackfname);
       }
       else if ( !strncmp(argv[ia],"-model",6) ) {
          strcpy(volfname,&argv[ia][7]);
          if (mypid == 0)
             printf("initial model      = %s\n", volfname);
       }
       else if ( !strncmp(argv[ia],"-param",6) ) {
          strcpy(paramfname,&argv[ia][7]);
          if (mypid == 0)
             printf("initial parameters = %s\n", paramfname);
       }
       else if ( !strncmp(argv[ia],"-out",4) ) {
          strcpy(voutfname,&argv[ia][5]);
          if (mypid == 0)
             printf("output volume file = %s\n", voutfname);
       }
       ia++;
    }

    // read and broadcast the initial model
    t0 = MPI_Wtime();
    EMData *volume = new EMData();
    ierr = ReadVandBcast(comm, volume, volfname);
    if (mypid == 0) {
       printf("Finished reading and replicating volume\n");
       printf("I/O time for reading volume = %11.3e\n",
              MPI_Wtime() - t0);
    }
    
    // read and distribute a stack of experimental images
    t0 = MPI_Wtime();
    nloc = ReadStackandDist(comm, &expimages, stackfname);
    if (mypid == 0) {
       printf("Finished reading and distributing image stack\n");
       printf("I/O time for reading image stack = %11.3e\n",
              MPI_Wtime() - t0);
    }

    Vec3i volsize;
    Vec3i origin;
    volsize[0] = volume->get_xsize();
    volsize[1] = volume->get_ysize();
    volsize[2] = volume->get_zsize();
    origin[0] = volume->get_xsize()/2 + 1;
    origin[1] = volume->get_ysize()/2 + 1;
    origin[2] = volume->get_zsize()/2 + 1;
    int    ri = volume->get_xsize()/2 - 2;
    ierr = CleanStack(comm, expimages, nloc, ri, volsize, origin);

    float * angleshift = new float[5*nloc];
    
    int maxiter = 10;

    try {
	unified(comm, volume, expimages, angleshift, nloc, 
                maxiter, voutfname);
    }
    catch (std::exception const& e) {
	printf("%s\n", e.what());
    }

    EMDeletePtr(volume);
    for ( int i = 0 ; i < nloc; ++i ) {
	EMDeletePtr(expimages[i]);
    }
    EMDeleteArray(expimages);
    EMDeleteArray(angleshift);
    ierr = MPI_Finalize();
    return 0;
}

