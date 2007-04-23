#include "mpi.h"
#include "stdlib.h"

// include EMAN2
#include "emdata.h"
#include "assert.h"

#include "ali3d_d_mpi.h"
#include "ali3d_unified_mpi.h"
#include "alignoptions.h"

int ReadVandBcast(MPI_Comm comm, EMData *volume, char *volfname);
int ReadStackandDist(MPI_Comm comm, EMData ***expimages, char *stackfname);
int setpart(MPI_Comm comm, int nima, int *psize, int *nbase);

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
    if (argc < 3) {
      if (mypid == 0) {
          printf("Not enough arguments to the command...\n");
          printf("Usage: sxhrefine -data=<imagestack> ");
          printf("-model=<initial 3D volume>\n"); 
      }
      ierr = MPI_Finalize();
      exit(1);
    }
    int ia=0;
    while (ia < argc) {
       if ( !strncmp(argv[ia],"-data",5) ) {
          strcpy(stackfname,&argv[ia][6]);
       }
       else if ( !strncmp(argv[ia],"-model",6) ) {
          strcpy(volfname,&argv[ia][7]);
       }
       else if ( !strncmp(argv[ia],"-param",6) ) {
          strcpy(paramfname,&argv[ia][7]);
       }
       else if ( !strncmp(argv[ia],"-out",4) ) {
          strcpy(voutfname,&argv[ia][5]);
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

    float * angleshift = new float[5*nloc];
    
    int max_refine_cycle = 2;
    int max_iter_ali3d   = 1;
    int max_iter_unified = 5;

    char out_fname[64];

    AlignOptions options;
    options.set_mask3D(NULL);
    options.set_first_ring(1);
    options.set_last_ring(26);
    options.set_rstep(1);
    options.set_ri(30);
    options.set_xrng(2.0);
    options.set_yrng(2.0);
    options.set_step(1.0);
    options.set_dtheta(15.0);
    options.set_snr(1.0);
    options.set_symmetry("c1");
    options.set_CTF(false);
    options.set_have_angles(false);

    try {
	for ( int iter = 0; iter < max_refine_cycle ; ++iter ) {
            if (mypid == 0) printf("refinement cycle: %d\n", iter+1);
	    sprintf(out_fname, "major%d", iter);

	    ali3d_d(comm, volume, expimages, angleshift, nloc, options, 
                    max_iter_ali3d, out_fname);

	    unified(comm, volume, expimages, angleshift, nloc, 
                    max_iter_unified, out_fname);
	    options.set_have_angles(true);
	}
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

