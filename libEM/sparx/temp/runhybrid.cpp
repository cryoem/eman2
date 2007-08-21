#include "mpi.h"
#include "stdlib.h"

// include EMAN2
#include "emdata.h"
#include "emassert.h"

#include "ali3d_d_mpi.h"
#include "ali3d_unified_mpi.h"
#include "alignoptions.h"
#include "utilparse.h"
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
    char  volfname[100], stackfname[100],voutfname[100],optionsfname[100];
    voutfname[0] = 0; // default to empty string
    optionsfname[0] = 0;
    EMData **expimages;

    // parse the command line and set filenames	
    if (argc < 3) {
      if (mypid == 0) {
          printf("Not enough arguments to the command...\n");
          printf("Usage: sxhrefine -data=<imagestack> ");
          printf("-model=<initial 3D volume> "); 
	  printf("-out=<output filename base string> ");
	  printf("-options=<options filename>\n");
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
       else if ( !strncmp(argv[ia],"-out",4) ) {
          strcpy(voutfname,&argv[ia][5]);
       }
       else if ( !strncmp(argv[ia],"-options",8) ) {
	   strcpy(optionsfname,&argv[ia][9]);
       }
       ia++;
    }

    // read and broadcast the initial model
    t0 = MPI_Wtime();
    EMData *volume = new EMData();
    ierr = ReadVandBcast(comm, volume, volfname);
    if (ierr == 0) {
	if (mypid == 0) {
	   printf("Finished reading and replicating volume\n");
	   printf("I/O time for reading volume = %11.3e\n",
		  MPI_Wtime() - t0);
	}
    }
    else {
       if (mypid == 0) printf("Failed to read the model volume %s! exit...\n", volfname);
       ierr = MPI_Finalize();
       exit(1);
    }

    // read and distribute a stack of experimental images
    t0 = MPI_Wtime();
    ierr = ReadStackandDist(comm, &expimages, stackfname, &nloc);
    if (ierr == 0) {
	if (mypid == 0) {
	   printf("Finished reading and distributing image stack\n");
	   printf("I/O time for reading image stack = %11.3e\n",
		  MPI_Wtime() - t0);
	}
    }
    else {
       if (mypid == 0) printf("Failed to read the image stack %s! exit...\n",stackfname);
       ierr = MPI_Finalize();
       exit(1);
    }

    Vec3i volsize;
    Vec3i origin;
    volsize[0] = volume->get_xsize();
    volsize[1] = volume->get_ysize();
    volsize[2] = volume->get_zsize();
    origin[0] = volume->get_xsize()/2 + 1;
    origin[1] = volume->get_ysize()/2 + 1;
    origin[2] = volume->get_zsize()/2 + 1;
    int ri = volume->get_xsize()/2 - 2;

    // make a copy of the images for removing the background; this stack will be used for reconstruction
    EMData** cleanimages = new EMData*[nloc];
    for ( int i = 0 ; i < nloc ; ++i) {
	cleanimages[i] = expimages[i]->copy();
    }
    ierr = CleanStack(comm, cleanimages, nloc, ri, volsize, origin);


    float * angleshift = new float[5*nloc];
    

    // these hard-coded numbers need to be removed eventually
    int max_refine_cycle = 0;
    int max_iter_unified = 10;

    char out_fname[128];

    AlignOptions options(volsize);


    std::ifstream option_stream;
    std::string current_option;
    char option_buffer[100];
    int int_option;
    float float_option;

    EMData * mask3D = NULL;
    if ( optionsfname[0] ) {
	ParseAlignOptions(comm, options, optionsfname, volsize[0]*volsize[1]*volsize[2], mask3D);
    } else {
	if ( mypid == 0 ) printf("Using default alignment options\n");
    }
    options.set_have_angles(false);

    // the below is a bit of a kluge; maxit in the options file indicates total number of refinement cycles, 
    // options.maxit is number of projection matching steps within each refinement cycle
    max_refine_cycle = options.get_maxit();
    options.set_maxit(1);

    try {
	for ( int iter = 0; iter < max_refine_cycle ; ++iter ) {
            if (mypid == 0) printf("refinement cycle: %d\n", iter+1);
	    sprintf(out_fname, "%smajor%d", voutfname, iter);

	    ali3d_d(comm, volume, expimages, cleanimages, angleshift, nloc, options, out_fname);

	    unified(comm, volume, cleanimages, angleshift, nloc, 
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
    for ( int i = 0 ; i < nloc; ++i ) {
	EMDeletePtr(cleanimages[i]);
    }
    EMDeleteArray(cleanimages);
    EMDeleteArray(angleshift);
    if (mask3D != NULL) {
	EMDeletePtr(mask3D);
    }
    ierr = MPI_Finalize();
    return 0;
}

