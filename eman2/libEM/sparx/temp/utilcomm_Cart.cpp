#include "mpi.h"
#include "stdlib.h"
#include "emdata.h"

#include "ali3d_unified_mpi.h"
#include "utilcomm_Cart.h"
#include "utilcomm.h"

using namespace EMAN;

int ReadStackandDist_Cart(MPI_Comm comm_2d, EMData ***images2D, char *stackfname, int *nloc)
{
    int ncpus, ierr, mpierr=0;
    MPI_Status mpistatus;
    EMUtil *my_util;
    int nima; 
    FILE *fp=NULL;
    int ROW = 0, COL = 1;
    int my2dpid, mycoords[2], dims[2], periods[2];
    int srpid, srcoords[2]; // Send/receive info

    // Get dims and my coordinates
    MPI_Cart_get(comm_2d, 2, dims, periods, mycoords);
    MPI_Comm_rank(comm_2d, &my2dpid); //Get my pid in the new 2D topology
    ncpus = dims[ROW]*dims[COL];

    ierr = 0;
    if (mycoords[ROW] == 0 && mycoords[COL] == 0) {
        fp = fopen(stackfname,"r");
        if (!fp) {
 	    ierr = 1;
            printf("failed to open %s\n", stackfname);
        }
        else {
            fclose(fp);
            nima = my_util->get_image_count(stackfname);
/**********************************************/
//    nima = 84;
/****************************/
        }
    }
    mpierr = MPI_Bcast(&ierr, 1, MPI_INT, 0, comm_2d);
    if (ierr == 0) {
	int *psize = new int[dims[ROW]];
	int *nbase = new int[dims[ROW]];
    
        mpierr = MPI_Bcast(&nima, 1, MPI_INT, 0, comm_2d);
    
	*nloc = setpart_gc1(comm_2d, nima, psize, nbase);
	*images2D = new EMData*[*nloc]; // NB!: whoever calls ReadStackandDist must delete this!

	EMData *img_ptr;
	int img_index;
	float *imgdata;

	// read the first image to get size
	img_ptr = new EMData();
	img_ptr->read_image(stackfname, 0);
	int nx = img_ptr->get_xsize();
	int ny = img_ptr->get_ysize();

	float s2x, s2y;

	if (mycoords[COL] == 0 && mycoords[ROW] == 0) {
	    printf("Master node reading and distributing %d images along the first column of processors\n", nima);
	    for ( int ip = 0 ; ip < dims[ROW] ; ++ip ) {
		for ( int i = 0 ; i < psize[ip] ; ++i ) {
		    img_index = nbase[ip] + i;

		    if (ip != 0) {
			img_ptr->read_image(stackfname, img_index);
			// get a pointer to the image's data
			imgdata = img_ptr->get_data();
			// find the x/y shift values if it has them, otherwise set them to 0.0
			try {
			    s2x = (*images2D)[i]->get_attr("s2x");
			} catch ( std::exception& e ) {
			    s2x = 0.0;
			}
			try {
			    s2y = (*images2D)[i]->get_attr("s2y");
			} catch ( std::exception& e ) {
			    s2y = 0.0;
			}
			// send these to processor (ip,0)
                        srcoords[COL] = 0;
	                srcoords[ROW] = ip;
	                MPI_Cart_rank(comm_2d, srcoords, &srpid);

			MPI_Send(imgdata, nx*ny, MPI_FLOAT, srpid, srpid, comm_2d);
			MPI_Send(&s2x, 1, MPI_FLOAT, srpid, srpid, comm_2d);
			MPI_Send(&s2y, 1, MPI_FLOAT, srpid, srpid, comm_2d);
		    } else { // ip == 0				    
			(*images2D)[i] = new EMData();
			(*images2D)[i]->read_image(stackfname, img_index);
			try {
			    s2x = (*images2D)[i]->get_attr("s2x");
			} catch ( std::exception& e ) {
			    (*images2D)[i]->set_attr("s2x",0.0);
			}
			try {
			    s2y = (*images2D)[i]->get_attr("s2y");
			} catch ( std::exception& e ) {
			    (*images2D)[i]->set_attr("s2y",0.0);
			}
		    }
		}
		//printf("finished reading data for processor %d\n", ip);
	    }

	} else if (mycoords[COL] == 0 && mycoords[ROW] != 0) {
	    for ( int i = 0 ; i < psize[mycoords[ROW]] ; ++i ) {
		(*images2D)[i] = new EMData();
		(*images2D)[i]->set_size(nx, ny, 1);
		imgdata = (*images2D)[i]->get_data();
		
                MPI_Recv(imgdata, nx*ny, MPI_FLOAT, 0, my2dpid, comm_2d, &mpistatus);
		MPI_Recv(&s2x, 1, MPI_FLOAT, 0, my2dpid, comm_2d, &mpistatus);
		MPI_Recv(&s2y, 1, MPI_FLOAT, 0, my2dpid, comm_2d, &mpistatus);
		(*images2D)[i]->set_attr("s2x",s2x);
		(*images2D)[i]->set_attr("s2y",s2y);
	    }
	    //printf("received %d images for processor (%d,%d)\n", *nloc, mycoords[ROW], mycoords[COL]);
        }

// Now have all the processors in group g_c_0 broadcast the images and angles along the row communicator
    if (mycoords[COL] == 0 && mycoords[ROW] == 0)
      printf("First column of processors distributing images to other columns..\n");

    if (mycoords[COL] == 0) {
      for (int img_index = 0; img_index < *nloc; img_index++){
        imgdata = (*images2D)[img_index]->get_data();
	
	// find the x/y shift values if it has them, otherwise set them to 0.0
	try {
	    s2x = (*images2D)[img_index]->get_attr("s2x");
	} catch ( std::exception& e ) {
	    s2x = 0.0;
	}
	try {
	    s2y = (*images2D)[img_index]->get_attr("s2y");
	} catch ( std::exception& e ) {
	    s2y = 0.0;
	}

        for(int ip=1; ip< dims[COL]; ip++){
// printf("Proc (%d, %d) sending image %d of %d to Proc (%d,%d) \n", mycoords[ROW], mycoords[COL], img_index, *nloc, mycoords[ROW], ip); 
          srcoords[ROW] = mycoords[ROW];
          srcoords[COL] = ip;
          MPI_Cart_rank(comm_2d, srcoords, &srpid);

          MPI_Send(imgdata, nx*ny, MPI_FLOAT, srpid, img_index, comm_2d);
	  MPI_Send(&s2x, 1, MPI_FLOAT, srpid, img_index, comm_2d);
	  MPI_Send(&s2y, 1, MPI_FLOAT, srpid, img_index, comm_2d);
        }
      }
    }
    if (mycoords[COL] != 0) {
  //      printf("Proc (%d, %d) receiving...", mycoords[ROW], mycoords[COL]);
        for(int img_index = 0; img_index < *nloc; img_index++){
          (*images2D)[img_index] = new EMData();
	  (*images2D)[img_index]->set_size(nx, ny, 1);
	  imgdata = (*images2D)[img_index]->get_data();
          
          srcoords[ROW] = mycoords[ROW];
          srcoords[COL] = 0;
          MPI_Cart_rank(comm_2d, srcoords, &srpid);

          MPI_Recv(imgdata, nx*ny, MPI_FLOAT, srpid, img_index, comm_2d, &mpistatus);
	  MPI_Recv(&s2x, 1, MPI_FLOAT, srpid, img_index, comm_2d, &mpistatus);
	  MPI_Recv(&s2y, 1, MPI_FLOAT, srpid, img_index, comm_2d, &mpistatus);

          (*images2D)[img_index]->set_attr("s2x",s2x);
	  (*images2D)[img_index]->set_attr("s2y",s2y);
         }
    }

        ierr = mpierr; 

	EMDeletePtr(img_ptr);
	EMDeleteArray(psize);
	EMDeleteArray(nbase);
    }
    return ierr;
}

int ReadAngTrandDist_Cart(MPI_Comm comm_2d, MPI_Comm comm_row, int *dims, float *angleshift, char *angfname, int nloc)
{
   int ierr = 0;
   int ROW = 0, COL = 1;
   int my2dpid, mycoords[2];
   int srpid, srcoords[2], keep_dims[2];
   MPI_Status mpistatus;
   FILE *fp = NULL;

   int nimgs=0;

   MPI_Comm_rank(comm_2d, &my2dpid); //Get my pid in the new 2D topology
   MPI_Cart_coords(comm_2d, my2dpid, 2, mycoords); // Get my coordinates

   float * iobuffer   = new float[5*nloc];
   if (!iobuffer) {
      fprintf(stderr,"failed to allocate buffer to read angles shifts\n");
      ierr = -1;
      goto EXIT;
   }

   if (mycoords[COL] == 0 && mycoords[ROW] == 0) { //I am Proc (0,0)
      fp = fopen(angfname,"r");
      if (!fp)  ierr = 1;
   }
   MPI_Bcast(&ierr, 1, MPI_INT, 0, comm_2d);

   if ( ierr ) {
      if (mycoords[COL] == 0 && mycoords[ROW] == 0) 
          fprintf(stderr,"failed to open %s\n", angfname);
      ierr = MPI_Finalize();
      goto EXIT;
   }
   else {
       if (mycoords[COL] == 0 && mycoords[ROW] == 0) { //I am Proc (0,0)
	  for (int iproc = 0; iproc < dims[ROW]; iproc++) {
	     // figure out the number of images assigned to processor (iproc,0)
	     if (iproc > 0) {
		srcoords[COL] = 0;
		srcoords[ROW] = iproc;
		MPI_Cart_rank(comm_2d, srcoords, &srpid);

		MPI_Recv(&nimgs, 1, MPI_INT, srpid, srpid, comm_2d, &mpistatus);

		// Read the next nimgs set of angles and shifts
		for (int i = 0; i < nimgs; i++) {
		   fscanf(fp,"%f %f %f %f %f", 
			  &iobuffer[5*i+0],
			  &iobuffer[5*i+1],
			  &iobuffer[5*i+2],
			  &iobuffer[5*i+3],
			  &iobuffer[5*i+4]);
		}
		MPI_Send(iobuffer,5*nimgs,MPI_FLOAT,srpid,srpid,comm_2d);
	     }
	     else {
		for (int i = 0; i < nloc; i++) {
		   fscanf(fp,"%f %f %f %f %f", 
			  &angleshift[5*i+0],
			  &angleshift[5*i+1],
			  &angleshift[5*i+2],
			  &angleshift[5*i+3],
			  &angleshift[5*i+4]);
		}
	     }
	  }
	  fclose(fp);
       }
       else if (mycoords[COL] == 0 && mycoords[ROW] != 0) { //I am in the first column
	  // send image count to the master processor (mypid = 0)

	  MPI_Send(&nloc, 1, MPI_INT, 0, my2dpid, comm_2d);
	  // Receive angleshifts
	  MPI_Recv(angleshift, 5*nloc, MPI_FLOAT, 0, my2dpid, comm_2d, &mpistatus);
       }
  }

  // Now have all the processors in group g_c_0 broadcast the angles along the row communicator
  srcoords[ROW] = 0;
  MPI_Cart_rank(comm_row, srcoords, &srpid);
  MPI_Bcast(angleshift, 5*nloc, MPI_FLOAT, srpid, comm_row);

  EMDeleteArray(iobuffer);

EXIT:
   return ierr;
}

int CleanStack_Cart(MPI_Comm comm_col, EMData ** image_stack, int nloc, int ri, Vec3i volsize, Vec3i origin)
{
    int nx = volsize[0];
    int ny = volsize[1];
    	
    float * rhs = new float[nx*ny];
    float * imgdata;

    // Calculate average "background" from all pixels strictly outside of radius
    double aba, abaloc; // average background
    aba = 0.0;
    abaloc = 0.0;
    int klp, klploc; // number of pixels in background
    klp = 0;
    klploc = 0;

    // calculate avg background in parallel
    for ( int i = 0 ; i < nloc ; ++i ) {
 	imgdata = image_stack[i]->get_data();
	asta2(imgdata, nx, ny, ri, &abaloc, &klploc);
    }
	
    // All reduce using the column-based sub-topology

    MPI_Allreduce(&abaloc, &aba, 1, MPI_DOUBLE, MPI_SUM, comm_col);
    MPI_Allreduce(&klploc, &klp, 1, MPI_INT, MPI_SUM, comm_col);

    aba /= klp;

    // subtract off the average background from pixels weakly inside of radius
    int x_summand, y_summand;
    int r_squared = ri * ri;
    for ( int i = 0 ; i < nloc ; ++i ) {
	imgdata = image_stack[i]->get_data();
	for ( int j = 0 ; j < nx ; ++j) {
	    x_summand = (j - origin[0]) *  (j - origin[0]);
	    for ( int k = 0 ; k < ny ; ++k ) {
		y_summand = (k - origin[1]) *  (k - origin[1]);
		if ( x_summand + y_summand <= r_squared) {
		    imgdata[j*ny + k] -= aba;
		}
	    }
	}
    }
    EMDeleteArray(rhs);
    return 0;
}

int setpart_gc1(MPI_Comm comm_2d, int nangs, int *psize, int *nbase)
/* This function provides the partition for the set of 2D images along the first column of processors

	Input:
	    comm_2d - Cartesian communicator
	      nangs - number of 2D images
	Output:
	      psize - vector containing the partition distribution
	      nbase - vector containing the base of the partitions
	   nangsloc - number of local angles
*/
{
	int ROW = 0, COL = 1;
	int dims[2], periods[2], mycoords[2];
	int nangsloc, nrem;

	// Get information associated with comm_2d
	MPI_Cart_get(comm_2d, 2, dims, periods, mycoords);
 
	nangsloc = nangs/dims[ROW];
	nrem = nangs - dims[ROW]*nangsloc;
	if (mycoords[ROW] < nrem) nangsloc++;

	for (int i = 0; i < dims[ROW]; i++) {
		psize[i] = nangs/dims[ROW];
		if (i < nrem) psize[i]++;
	}
 
	for (int i = 0; i < dims[ROW]; i++) {
          if(i==0)
            nbase[i] = 0;
          else
	    nbase[i] = nbase[i-1] + psize[i-1];
	}
   
//printf("I am [%d,%d], nloc = %d, psize = [%d, %d], nbase = [%d, %d]", mycoords[ROW], mycoords[COL],nangsloc,psize[0],psize[1], nbase[0], nbase[1]);

	return nangsloc;
}


int setpart_gr1(MPI_Comm comm_2d, int nnz, int *nnzpart, int *nnzbase)
/* This function provides the ideal partition for nnz in the 3D volume along the first row of processors.

	Input:
	    comm_2d - Cartesian communicator
	        nnz - total number of non-zeros in the volume
	Output:
	    nnzpart - vector containing the partition distribution
	    nnzbase - vector containing the base of the partitions
	     nnzloc - initial local nnz
*/
{
	int ROW = 0, COL = 1;
	int dims[2], periods[2], mycoords[2];
	int nnzloc, nrem;
	
	// Get information associated with comm_2d
	MPI_Cart_get(comm_2d, 2, dims, periods, mycoords);
 
	nnzloc = nnz/dims[COL];
	nrem = nnz - dims[COL]*nnzloc;
	if (mycoords[COL] < nrem) nnzloc++;

	for (int i = 0; i < dims[COL]; i++) {
		nnzpart[i] = nnz/dims[COL];
		if (i < nrem) nnzpart[i]++;
	}
 
	nnzbase[0] = 0; 
	for (int i = 1; i < dims[COL]; i++) {
		nnzbase[i] = nnzbase[i-1] + nnzpart[i-1];
	}
 	nnzbase[dims[COL]] = nnz;
	return nnzloc;
}


