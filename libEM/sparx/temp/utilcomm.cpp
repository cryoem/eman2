#include "mpi.h"
#include "stdlib.h"
#include "emdata.h"

using namespace EMAN;

int ReadVandBcast(MPI_Comm comm, EMData *volume, char *volfname);
int ReadStackandDist(MPI_Comm comm, EMData ***images2D, char *stackfname);
int setpart(MPI_Comm comm, int nima, int *psize, int *nbase);

int ReadVandBcast(MPI_Comm comm, EMData *volume, char *volfname)
{
    int mypid, ierr, ndim, nx, ny, nz;
    MPI_Comm_rank(comm,&mypid);

    if (mypid == 0) {
	volume->read_image(volfname);
	ndim = volume->get_ndim();
	nx   = volume->get_xsize();
	ny   = volume->get_ysize();
	nz   = volume->get_zsize();
//debug
	printf("ndim = %d, nx = %d, ny = %d, nz = %d\n", ndim, nx, ny, nz);
    }
    ierr = MPI_Bcast(&nx, 1, MPI_INT, 0, comm);
    ierr = MPI_Bcast(&ny, 1, MPI_INT, 0, comm);
    ierr = MPI_Bcast(&nz, 1, MPI_INT, 0, comm);
    if (mypid !=0) volume->set_size(nx,ny,nz);
	
    float * voldata = volume->get_data();
    ierr = MPI_Bcast(voldata, nx*ny*nz, MPI_FLOAT, 0, comm);
    return ierr;
}

int ReadStackandDist(MPI_Comm comm, EMData ***images2D, char *stackfname)
{
    int ncpus, mypid, ierr;
    MPI_Status mpistatus;
    EMUtil *my_util;
    int nima, nloc;

    MPI_Comm_size(comm,&ncpus);
    MPI_Comm_rank(comm,&mypid);
    if (mypid == 0) nima = my_util->get_image_count(stackfname);
    ierr = MPI_Bcast(&nima, 1, MPI_INT, 0, comm);
       
    int *psize = new int[ncpus];
    int *nbase = new int[ncpus];
    nloc = setpart(comm, nima, psize, nbase);

    *images2D = new (EMData*)[nloc];
	
    EMData *img_ptr;
    int img_index;
    float *imgdata;
    
    // read the first image to get size
    img_ptr = new EMData();
    img_ptr->read_image(stackfname, 0);
    int nx = img_ptr->get_xsize();
    int ny = img_ptr->get_ysize();

    if (mypid == 0) {
	for ( int ip = 0 ; ip < ncpus ; ++ip ) {
	    for ( int i = 0 ; i < psize[ip] ; ++i ) {
		img_index = nbase[ip] + i;
		if (ip != 0) {
		    img_ptr->read_image(stackfname, img_index);
		    imgdata = img_ptr->get_data();
		    MPI_Send(imgdata, nx*ny, MPI_FLOAT, ip, ip, comm);
		} else { // ip == 0				    
		    (*images2D)[i] = new EMData();
		    (*images2D)[i]->read_image(stackfname, img_index);
		    (*images2D)[i]->set_attr("s2x",0.0);
		    (*images2D)[i]->set_attr("s2y",0.0);
		}
	    }
	    printf("finished reading data for processor %d\n", ip);
	}
    } else { // mypid != 0 : everyone else receives and reads in their data
	for ( int i = 0 ; i < psize[mypid] ; ++i ) {
	    (*images2D)[i] = new EMData();
	    (*images2D)[i]->set_size(nx, ny, 1);
	    imgdata = (*images2D)[i]->get_data();
	    MPI_Recv(imgdata, nx*ny, MPI_FLOAT, 0, mypid, comm, &mpistatus);
	    (*images2D)[i]->set_attr("s2x",0.0);
	    (*images2D)[i]->set_attr("s2y",0.0);
	}
	printf("received data for processor %d\n", mypid);
    }
    if (mypid == 0) printf("finished reading and distributing data\n");
    EMDeletePtr(img_ptr);
    return nloc;
}


int setpart(MPI_Comm comm, int nang, int *psize, int *nbase)
{
   int ncpus, mypid, nangloc, nrem;

   MPI_Comm_size(comm,&ncpus);
   MPI_Comm_rank(comm,&mypid);

   nangloc = nang/ncpus;
   nrem    = nang - ncpus*nangloc;
   if (mypid < nrem) nangloc++;

   for (int i = 0; i < ncpus; i++) {
      psize[i] = nang/ncpus;
      if (i < nrem) psize[i]++;
   }
 
   nbase[0] = 0; 
   for (int i = 1; i < ncpus; i++) {
      nbase[i] = nbase[i-1] + psize[i-1];
   }
   
   return nangloc;
}
