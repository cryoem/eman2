#include "mpi.h"

#include "emdata.h"

#include "sirt_Cart.h"
#include "utilcomm_Cart.h"
#include "project3d.h"
#include "HyBR_Cart.h"

#define PI 3.14159265358979
#define ROW 0
#define COL 1

using namespace EMAN;

// Reconstruct the 3-D density of a macromolecule from
// a stack of 2-D images using the SIRT algorithm
//
// MPI Parallel Implementation using Cartesian virtual topology

int main(int argc, char ** argv)
{
   MPI_Comm comm = MPI_COMM_WORLD;
   int ncpus, mypid, ierr, mpierr=0;
   int nloc, maxit=0, method=0; 
   float lam = 0.0, tol = 0.0;
   char symstr[10]="none";
   std::string symmetry;
   double t0;
   FILE *fp;
   MPI_Comm comm_2d, comm_row, comm_col;

   MPI_Status mpistatus;
   MPI_Init(&argc,&argv);
   MPI_Comm_size(comm,&ncpus);
   MPI_Comm_rank(comm,&mypid);

   int dims[2], periods[2], my2dpid, mycoords[2];
   int srpid, srcoords[2], keep_dims[2];

   char  stackfname[100],voutfname[100], paramfname[100];
   EMData **expimages;

   // parse the command line and set filenames	
   if (argc < 5) {
     if (mypid == 0) {
         printf("Not enough arguments to the command...\n");
         printf("Usage: runsirt -data=<imagestack> ");
         printf("-param=<parameter file that contains angles and shifts> "); 
         printf("-out=<output filename base string> ");
         printf("-maxit=<maximum number of iterations> ");
         printf("-lam=<lambda> ");
         printf("-tol=<convergence tolerance> ");
         printf("-sym=<symmtry> \n");
         printf("-rowdim=<row dimension of Cartesian topology> ");
         printf("-coldim=<column dimension of Cartesian topology> ");
         printf("-method=<method to use: 0 SIRT (default), 1 Hybrid>\n");
     }
     ierr = MPI_Finalize();
     exit(1);
   }
   int ia=1;
   while (ia < argc) {
      if ( !strncmp(argv[ia],"-data",5) ) {
         strcpy(stackfname,&argv[ia][6]);
      }
      else if ( !strncmp(argv[ia],"-param",6) ) {
         strcpy(paramfname,&argv[ia][7]);
      }
      else if ( !strncmp(argv[ia],"-out",4) ) {
         strcpy(voutfname,&argv[ia][5]);
      }
      else if ( !strncmp(argv[ia],"-rowdim",7) ) {
         dims[ROW] = atoi(&argv[ia][8]); // Row dimension of the topology
      }
      else if ( !strncmp(argv[ia],"-coldim",7) ) {
         dims[COL] = atoi(&argv[ia][8]); // Column dimension of the topology
      }
      else if ( !strncmp(argv[ia],"-maxit",6) ) {
         maxit = atoi(&argv[ia][7]);
      }
      else if ( !strncmp(argv[ia],"-lam",4) ) {
         lam = atof(&argv[ia][5]);
      }
      else if ( !strncmp(argv[ia],"-tol",4) ) {
         tol = atof(&argv[ia][5]);
      }
      else if ( !strncmp(argv[ia],"-sym",4) ) {
         strcpy(symstr, &argv[ia][5]);
      }
      else if ( !strncmp(argv[ia],"-method",7) ) {
         method = atoi(&argv[ia][8]);
      }
      else {
         if (mypid ==0) printf("invalid option: %s\n", argv[ia]);
         ierr = MPI_Finalize();
         exit(1);
      }
      ia++;
   }

  if (dims[ROW]*dims[COL] != ncpus){
	printf("ERROR: rowdim*coldim != ncpus\n");
        ierr = MPI_Finalize();
	return -1;
  }

  // Set up the Cartesian virtual topology: comm_2d
  periods[ROW] = periods[COL] = 1; // Set the periods for wrap-around
  MPI_Cart_create(comm, 2, dims, periods, 1, &comm_2d);
  MPI_Comm_rank(comm_2d, &my2dpid); //Get my pid in the new 2D topology
  MPI_Cart_coords(comm_2d, my2dpid, 2, mycoords); // Get my coordinates
  
  // printf("MPI_2d: mypid = %d, my2dpid = %d, mycoords = [%d, %d] \n", 
  //mypid, my2dpid, mycoords[ROW], mycoords[COL]);

  /* Create the row-based sub-topology */ 
  keep_dims[ROW] = 0; 
  keep_dims[COL] = 1; 
  MPI_Cart_sub(comm_2d, keep_dims, &comm_row); 

  /* Create the column-based sub-topology */ 
  keep_dims[ROW] = 1; 
  keep_dims[COL] = 0; 
  MPI_Cart_sub(comm_2d, keep_dims, &comm_col); 

   // read and distribute a stack of experimental images along row processors
   t0 = MPI_Wtime();
   ierr = ReadStackandDist_Cart(comm_2d, &expimages, stackfname, &nloc);
   if (ierr == 0) {
	if (mypid == 0) {
	   printf("Finished reading and distributing image stack onto ");
           printf("Cartesian topology\n");
	   printf("I/O time for reading image stack = %11.3e\n",
		  MPI_Wtime() - t0);
	}
   }
   else {
      if (mypid == 0) 
         printf("Failed to read the image stack %s! exit...\n",stackfname);
      ierr = MPI_Finalize();
      exit(1);
   }

   // make a copy of the images for removing the background; 
   // this stack will be used for reconstruction
   EMData** cleanimages = new EMData*[nloc];
   for ( int i = 0 ; i < nloc ; ++i) {
	cleanimages[i] = expimages[i]->copy();
   }

   int nx = cleanimages[0]->get_xsize();
   int ny = cleanimages[0]->get_ysize();
   int nz = nx;

   Vec3i volsize;
   Vec3i origin;
   volsize[0] = nx;
   volsize[1] = ny;
   volsize[2] = nz;
   origin[0] = nx/2 + 1;
   origin[1] = ny/2 + 1;
   origin[2] = nz/2 + 1;
   int ri = nx/2 - 2;

   ierr = CleanStack_Cart(comm_col, cleanimages, nloc, ri, volsize, origin);

   // read angle and shift data and distribute along first column
   float * angleshift = new float[5*nloc];

   ierr = ReadAngTrandDist_Cart(comm_2d, comm_row, dims, angleshift, 
                                paramfname, nloc);
   if (ierr!=0) { 
      mpierr = MPI_Finalize();
      return 1;
   }

   // Use xvol to hold reconstructed volume
   EMData * xvol = new EMData();

   // set SIRT parameters
  
   if (maxit==0) maxit = 10;
   if (lam == 0.0) lam = 5.0e-6;
   if (tol == 0.0) tol = 1.0e-3;
   if ( !strncmp(symstr,"none",4) ) {
      symmetry = string("c1");
   }
   else {
      symmetry = string(symstr);
   }
   // write out all options 
   if (mypid == 0) {
      printf("SIRT options used:\n");
      printf("maxit  = %d\n", maxit);
      printf("lam    = %11.3e\n", lam);
      printf("tol    = %11.3e\n", tol);
      printf("sym    = %s\n", symmetry.c_str()); 
      printf("rowdim = %d\n", dims[ROW]); 
      printf("coldim = %d\n", dims[COL]); 
   }

   if (method!=1) method = 0; //use SIRT unless method=1


   // call SIRT to reconstruct
   t0 = MPI_Wtime();
   if (method == 0) {
      recons3d_sirt_mpi_Cart(comm_2d, comm_row, comm_col, cleanimages, 
                             angleshift, xvol, nloc, ri, lam, maxit, 
                             symmetry, tol);
   }
   else {
      int insolve = 1;
      recons3d_HyBR_mpi_Cart(comm_2d, comm_row, comm_col, cleanimages, 
                             angleshift, xvol, nloc, ri, maxit, symmetry, insolve);
   }

   if ( my2dpid == 0 ) 
       printf("Done with reconstruction\n");

   // write the reconstructed volume to disk
   EMUtil::ImageType WRITE_SPI = EMUtil::IMAGE_SINGLE_SPIDER;
   if ( my2dpid == 0 ) {
	xvol->write_image(voutfname, 0, WRITE_SPI);
   }

   // cleanup
   for ( int i = 0 ; i < nloc; ++i ) {
       EMDeletePtr(expimages[i]);
   }
   EMDeleteArray(expimages);
   for ( int i = 0 ; i < nloc; ++i ) {
       EMDeletePtr(cleanimages[i]);
   }
   EMDeleteArray(cleanimages);

   EMDeletePtr(xvol);
   EMDeleteArray(angleshift);

   ierr = MPI_Finalize();

   return 0; // main
}
