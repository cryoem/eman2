#include "mpi.h"
#include "stdlib.h"
#include "emdata.h"

#include "ali3d_unified_mpi.h"
#include "utilcomm.h"

using namespace EMAN;

int ReadVandBcast(MPI_Comm comm, EMData *volume, char *volfname)
{
    int mypid, ierr, ndim, nx, ny, nz, mpierr=0;
    MPI_Comm_rank(comm,&mypid);
    FILE *fp=NULL;

    ierr = 0;   
    if (mypid == 0) {
        // check to see if the file exists
        fp = fopen(volfname,"r");
        if (!fp) {
 	    ierr = 1;
            printf("failed to open %s\n", volfname);
        }
        else {
            fclose(fp);
	    volume->read_image(volfname);
	    ndim = volume->get_ndim();
	    nx   = volume->get_xsize();
	    ny   = volume->get_ysize();
	    nz   = volume->get_zsize();
	    printf("ndim = %d, nx = %d, ny = %d, nz = %d\n", ndim, nx, ny, nz);
        }
    }
    mpierr = MPI_Bcast(&ierr, 1, MPI_INT, 0, comm);
    if (ierr==0) {
        mpierr = MPI_Bcast(&nx, 1, MPI_INT, 0, comm);
        mpierr = MPI_Bcast(&ny, 1, MPI_INT, 0, comm);
        mpierr = MPI_Bcast(&nz, 1, MPI_INT, 0, comm);
        if (mypid !=0) volume->set_size(nx,ny,nz);
	
        float * voldata = volume->get_data();
        mpierr = MPI_Bcast(voldata, nx*ny*nz, MPI_FLOAT, 0, comm);
        ierr = mpierr;
    }
    return ierr;
}

int ReadStackandDist(MPI_Comm comm, EMData ***images2D, char *stackfname, int *nloc)
{
    int ncpus, mypid, ierr, mpierr=0;
    MPI_Status mpistatus;
    EMUtil *my_util;
    int nima; 
    FILE *fp=NULL;

    MPI_Comm_size(comm,&ncpus);
    MPI_Comm_rank(comm,&mypid);

    ierr = 0;
    if (mypid == 0) {
        fp = fopen(stackfname,"r");
        if (!fp) {
 	    ierr = 1;
            printf("failed to open %s\n", stackfname);
        }
        else {
            fclose(fp);
            nima = my_util->get_image_count(stackfname);
        }
    }
    mpierr = MPI_Bcast(&ierr, 1, MPI_INT, 0, comm);
    if (ierr == 0) {
	int *psize = new int[ncpus];
	int *nbase = new int[ncpus];
    
        mpierr = MPI_Bcast(&nima, 1, MPI_INT, 0, comm);
    
	*nloc = setpart(comm, nima, psize, nbase);
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

	if (mypid == 0) {
	    printf("Master node reading and distributing %d images...\n", nima);
	    for ( int ip = 0 ; ip < ncpus ; ++ip ) {
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
			// send these to processor ip
			MPI_Send(imgdata, nx*ny, MPI_FLOAT, ip, ip, comm);
			MPI_Send(&s2x, 1, MPI_FLOAT, ip, ip, comm);
			MPI_Send(&s2y, 1, MPI_FLOAT, ip, ip, comm);
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
		printf("finished reading data for processor %d\n", ip);
	    }
	} else { // mypid != 0 : everyone else receives and reads in their data
	    for ( int i = 0 ; i < psize[mypid] ; ++i ) {
		(*images2D)[i] = new EMData();
		(*images2D)[i]->set_size(nx, ny, 1);
		imgdata = (*images2D)[i]->get_data();
		MPI_Recv(imgdata, nx*ny, MPI_FLOAT, 0, mypid, comm, &mpistatus);
		MPI_Recv(&s2x, 1, MPI_FLOAT, 0, mypid, comm, &mpistatus);
		MPI_Recv(&s2y, 1, MPI_FLOAT, 0, mypid, comm, &mpistatus);
		(*images2D)[i]->set_attr("s2x",s2x);
		(*images2D)[i]->set_attr("s2y",s2y);
	    }
	    printf("received %d images for processor %d\n", *nloc, mypid);
	}
	if (mypid == 0) printf("finished reading and distributing data\n");
        ierr = mpierr; 

	EMDeletePtr(img_ptr);
	EMDeleteArray(psize);
	EMDeleteArray(nbase);
    }
    return ierr;
}

int ReadAngTrandDist(MPI_Comm comm, float *angleshift, char *paramfname, int nloc)
{
    // read a parameter file that contains a list of angles and shifts and distribute 
    // them evenly among different processors

    int mypid, ierr=0, mpierr=0, nimgs=0, ncpus;
    double t0;
    MPI_Status mpistatus;
    FILE *fp=NULL;

    float *iobuffer = new float[5*nloc];
    if (!iobuffer) {
       fprintf(stderr,"failed to allocate buffer to read angles shifts\n");
       ierr = -1;
       goto EXIT;
    }

    MPI_Comm_rank(comm,&mypid);
    MPI_Comm_size(comm,&ncpus);

    t0 = MPI_Wtime();

    if (mypid ==0) {
       fp = fopen(paramfname,"r");
       if (!fp)  ierr = 1;
    }
    MPI_Bcast(&ierr, 1, MPI_INT, 0, comm);
 
    if ( ierr ) {
       if (mypid ==0) fprintf(stderr,"failed to open %s\n", paramfname);
    }
    else {
       if (mypid == 0) {
          for (int iproc = 0; iproc < ncpus; iproc++) {
             // figure out the number of images assigned to processor iproc
  	     if (iproc > 0) {
  	        MPI_Recv(&nimgs, 1, MPI_INT, iproc, iproc, comm, &mpistatus);
                // Read the next nimgs set of angles and shifts
                for (int i = 0; i < nimgs; i++) {
                   fscanf(fp,"%f %f %f %f %f", 
                          &iobuffer[5*i+0],
                          &iobuffer[5*i+1],
                          &iobuffer[5*i+2],
                          &iobuffer[5*i+3],
                          &iobuffer[5*i+4]);
                }
                MPI_Send(iobuffer,5*nimgs,MPI_FLOAT,iproc,iproc,comm);
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
       else {
          // send image count to the master processor (mypid = 0)
          MPI_Send(&nloc, 1, MPI_INT, 0, mypid, comm);
          // Receive angleshifts
          MPI_Recv(angleshift, 5*nloc, MPI_FLOAT, 0, mypid, comm, &mpistatus);
       }
    }
    EMDeleteArray(iobuffer);
    if (mypid == 0)
       printf("I/O time for reading angles & shifts = %11.3e\n",
              MPI_Wtime() - t0);
EXIT:
    return ierr;
}

int CleanStack(MPI_Comm comm, EMData ** image_stack, int nloc, int ri, Vec3i volsize, Vec3i origin)
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
	
    MPI_Allreduce(&abaloc, &aba, 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(&klploc, &klp, 1, MPI_INT, MPI_SUM, comm);

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

