#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <string.h>

// GLOBAL Variables
int nangs, nangsloc, nx, ny, myptrstart;
int nnz, nrays, nnzloc, nraysloc;
float *imagesloc, *images;
float *angles, *anglesloc;

// Additional functions that are found in backproj_Cart.cpp
int setpart_gc1(MPI_Comm comm_2d, int nangs, int *psize, int *nbase);
int getnnz(int *volsize, int ri, int *origin, int *nrays, int *nnz);
int setpart_gr1(MPI_Comm comm_2d, int nnz, int *nnzpart, int *nnzbase);
int getcb2sph(int *volsize, int ri, int *origin, int nnz0, int *ptrs, int *cord);
int sphpart(MPI_Comm comm_2d, int nrays, int *ptrs, int *nnzbase, int *ptrstart);
int make_proj_mat(float phi, float theta, float psi, float * dm);
int ifix(float a);
int bckpj3_Cart(int *volsize, int nrays, int nnz, float *dm, int *origin, int ri, int *ptrs, int *cord, int myptrstart, float *x, float *y);

// ---------------- MAIN FUNCTION --------------------------------
// This code will read in the 2D images and perform the backprojection operation in parallel

int main(int argc, char *argv[])
{
  int ncpus, mypid, nrem, ierr;
  MPI_Status mpistatus;
  MPI_Comm comm = MPI_COMM_WORLD;
  
  // Variables needed for the Cartesian topology.
  int ROW = 0, COL = 1;
  int dims[2], periods[2], keep_dims[2];
  int my2dpid, mycoords[2], srcoords[2], otherpid;
  MPI_Comm comm_2d, comm_row, comm_col; 
		  
// Initialize MPI.
  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &ncpus);
  MPI_Comm_rank(comm, &mypid);
  //printf("mypid = %d, ncpus = %d\n", mypid, ncpus);
  
  if ( argc < 3  ) {
	  printf ("ERROR: %s requires Cartesian dimensions input\n", argv[0]);
          MPI_Finalize();
	  return -1;
  }

// Set up a Cartesian virtual topology and get the rank and coordinates of the processes in the topology. 
  dims[ROW] = atoi(argv[1]); // Row dimension of the topology
  dims[COL] = atoi(argv[2]); // Column dimension of the topology
  
  if (dims[ROW]*dims[COL] != ncpus){
	printf("ERROR: Row dim and col dim not equal to ncpus\n");
        MPI_Finalize();
	return -1;
  }
  
  periods[ROW] = periods[COL] = 1; // Set the periods for wrap-around
  
  MPI_Cart_create(comm, 2, dims, periods, 1, &comm_2d);
  MPI_Comm_rank(comm_2d, &my2dpid); //Get my pid in the new 2D topology
  MPI_Cart_coords(comm_2d, my2dpid, 2, mycoords); // Get my coordinates
  
  /* Create the row-based sub-topology */ 
  keep_dims[ROW] = 0; 
  keep_dims[COL] = 1; 
  MPI_Cart_sub(comm_2d, keep_dims, &comm_row); 
 
  /* Create the column-based sub-topology */ 
  keep_dims[ROW] = 1; 
  keep_dims[COL] = 0; 
  MPI_Cart_sub(comm_2d, keep_dims, &comm_col); 

// STEP 1: Have processor (0,0) read in the entire set of 2D images, divide up the images, and send corresponding images to processors in processor group: g_c_0  Do the same for the angles.
  
  if (mycoords[ROW] == 0 && mycoords[COL] == 0){ //I'm processor (0,0)
    FILE *fp=NULL, *fpa=NULL;
    char imagefname[80]="tf2d84.raw", anglesfname[80]="angles.dat";
	  
    fp = fopen(imagefname,"r");
    if (!fp ) {
	printf("ERROR: cannot open %s\n", imagefname);
        MPI_Finalize();
	return -1;
    }
    fread(&nangs, sizeof(int), 1, fp);
    fread(&nx, sizeof(int), 1, fp);
    fread(&ny, sizeof(int), 1, fp);
    
    images = new float[nx*ny*nangs];
    fread(images, sizeof(float), nx*ny*nangs, fp);
    fclose(fp);
    
    fpa = fopen(anglesfname,"r");
    if (!fpa ) {
	printf("ERROR: cannot open %s\n", anglesfname);
        MPI_Finalize();
	return -1;
    }
    angles = new float[3*nangs];
    for (int i = 0; i< 3*nangs; i++)
      fscanf(fpa, "%f",&angles[i]);
       
    fclose(fpa);
  }
  
  // Broadcast variables nangs, nx, ny to all processors
  srcoords[ROW] = srcoords[COL] = 0;
  MPI_Cart_rank(comm_2d, srcoords, &otherpid); 
  
  MPI_Bcast (&nangs, 1, MPI_INT, otherpid, comm_2d);
  MPI_Bcast (&nx, 1, MPI_INT, otherpid, comm_2d);
  MPI_Bcast (&ny, 1, MPI_INT, otherpid, comm_2d);
  
  // Send images and angles from Processor (0,0) to processors in group g_c_0
  int *psize = new int[dims[ROW]];
  int *nbase = new int[dims[ROW]];
  nangsloc = setpart_gc1(comm_2d, nangs, psize, nbase);
  imagesloc = new float[psize[mycoords[ROW]]*nx*ny];
  anglesloc = new float[psize[mycoords[ROW]]*3];
  
  if (mycoords[COL] == 0 && mycoords[ROW] == 0) { //I'm Proc. (0,0)
    for(int ip = 0; ip < dims[ROW]; ++ip){
      int begidx = nbase[ip]*nx*ny;
      if (ip !=0){ // Proc (0,0) sends images and angle data to other processors
	 srcoords[COL] = 0;
	 srcoords[ROW] = ip;
	 MPI_Cart_rank(comm_2d, srcoords, &otherpid);
	 MPI_Send(&images[begidx],psize[ip]*nx*ny, MPI_FLOAT, otherpid, otherpid, comm_2d);
	 MPI_Send(&angles[nbase[ip]*3],psize[ip]*3, MPI_FLOAT, otherpid, otherpid, comm_2d);
	 //printf("Sent data to processor %d\n", otherpid);
      }
      else{ // ip = 0: Proc (0,0) needs to copy images and angles into its imagesloc and anglesloc
	for (int i = 0; i < psize[ip]*nx*ny; i++){
	  imagesloc[i] = images[begidx+i];
	}
	for (int i = 0; i < psize[ip]*3; i++){
		anglesloc[i] = images[nbase[ip]*3 + i];
	}
	//printf("Finished copying to Proc (0,0) local");
      }
    } //End for loop
  } //End if
  
  if (mycoords[COL] == 0 && mycoords[ROW] != 0) { //I'm in g_c_0 and I'm not Processor (0,0) so I should receive data.
    MPI_Recv(imagesloc, psize[mycoords[ROW]]*nx*ny, MPI_FLOAT, 0, mypid, comm_2d, &mpistatus);
    MPI_Recv(anglesloc, psize[mycoords[ROW]]*3, MPI_FLOAT, 0, mypid, comm_2d, &mpistatus);
    //printf("received data for processor %d\n", mypid);
  }
  // Now have all the processors in group g_c_0 broadcast the images and angles along the row communicator
  srcoords[ROW] = 0;
  MPI_Cart_rank(comm_row, srcoords, &otherpid);
  MPI_Bcast(imagesloc, nangsloc*nx*ny, MPI_FLOAT, otherpid , comm_row);
  MPI_Bcast(anglesloc, nangsloc*3, MPI_FLOAT, otherpid , comm_row);
 
// Now distribute the volume (in spherical format) among columns of processors and use nnz to determine the splitting.  Note: ptrs and coord are on all processors
   int radius;
   int volsize[3], origin[3];
   volsize[0] = nx;
   volsize[1] = nx;
   volsize[2] = nx;
   origin[0] = nx/2+1;
   origin[1] = nx/2+1;
   origin[2] = nx/2+1;
   radius = nx/2-1;
   
   ierr = getnnz( volsize, radius, origin, &nrays, &nnz);
   if (mycoords[ROW] == 0 && mycoords[COL] == 0)
      printf("nrays is %d, nnz is %d\n", nrays, nnz);
   
  int * ptrs = new int[nrays+1];
  int * cord = new int[3*nrays];
  ierr = getcb2sph(volsize, radius, origin, nnz, ptrs, cord);
  
  int *nnzpart = new int[dims[COL]];
  int *nnzbase = new int[dims[COL]+1]; 
  nnzloc = setpart_gr1(comm_2d, nnz, nnzpart, nnzbase);
  
  int *ptrstart = new int[dims[COL]+1];
  nraysloc = sphpart(comm_2d, nrays, ptrs, nnzbase, ptrstart);
  
  myptrstart = ptrstart[mycoords[COL]];
  nnzloc = ptrs[ptrstart[mycoords[COL]+1]] - ptrs[myptrstart];
   
  // Print some stuff.
  printf("My coords are (%d,%d) and nraysloc = %d, myptrstart = %d , nnzloc = %d\n", mycoords[ROW], mycoords[COL], nraysloc, myptrstart, nnzloc);
  
  float *bvol_loc = new float[nnzloc];
  for (int i=0; i< nnzloc; i++)
    bvol_loc[i] = 0.0;
  
  // STEP 2: Have everyone perform the backprojection operation for their assigned images and portions of the volume.  Then perform an all_reduce along the columns.
  
  /*float phi, theta, psi;
  float dm[8];
  
  for (int i=0; i<nangsloc; i++){
    phi = anglesloc[3*i+0];
    theta = anglesloc[3*i+1];
    psi = anglesloc[3*i+2];
    dm[6] = 0;
    dm[7] = 0;
    
    make_proj_mat(phi, theta, psi, dm);
  
    ierr = bckpj3_Cart(volsize, nraysloc, nnzloc, dm, origin, radius, ptrs, cord, myptrstart, &imagesloc[nx*ny*i], bvol_loc);
  }
  
  */
  
  MPI_Comm_free(&comm_2d);
  MPI_Comm_free(&comm_row);
  MPI_Comm_free(&comm_col);
  
  MPI_Finalize();
}
