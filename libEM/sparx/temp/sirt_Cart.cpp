#include "mpi.h"

#include "emdata.h"

#include "sirt_Cart.h"
#include "utilcomm_Cart.h"
#include "project3d_Cart.h"
#include "project3d.h"

using namespace EMAN;

//int CleanStack(MPI_Comm comm, EMData ** image_stack, int nloc, int ri, Vec3i volsize, Vec3i origin);

int recons3d_sirt_mpi_Cart(MPI_Comm comm_2d, MPI_Comm comm_row, MPI_Comm comm_col, EMData ** images, float * angleshift, EMData *& xvol, int nangloc, int radius, float lam, int maxit, std::string symmetry, float tol)
{
    MPI_Status mpistatus;
    int ROW = 0, COL = 1;
    int ncpus, my2dpid, ierr;
    int mycoords[2], dims[2], periods[2];
    int srpid, srcoords[2]; // Send/receive info

    // Get dims and my coordinates
    MPI_Cart_get(comm_2d, 2, dims, periods, mycoords);
    MPI_Comm_rank(comm_2d, &my2dpid); //Get my pid in the new 2D topology
    ncpus = dims[ROW]*dims[COL];
    
    int * psize;
    int * nbase;

    int nangles;
    psize = new int[dims[ROW]];
    nbase = new int[dims[ROW]];
    MPI_Allreduce(&nangloc, &nangles, 1, MPI_INT, MPI_SUM, comm_col);
    

    int nsym = 0;
    // get image size from first image
    int nx = images[0]->get_xsize();
    if ( radius == -1 ) radius = nx/2 - 1; // make radius as large as possible if the user didn't provide one
    
    Vec3i volsize, origin;
    volsize[0] = nx;
    volsize[1] = nx;
    volsize[2] = nx;
    origin[0] = nx/2+1;
    origin[1] = nx/2+1;
    origin[2] = nx/2+1;
//    this is not currently needed, because the stack that gets passed to sirt will have its background subtracted already
//    ierr = CleanStack(comm, images, nangloc, radius, volsize, origin);

    // vector of symmetrized angles
    std::vector<float> symangles(3,0.0); 
	
    float old_rnorm = 1.00001; // kluge, make sure if its 1.0 + epsilon it still works;
	

// Now distribute the volume (in spherical format) among columns of processors and use nnz to determine the splitting.  Note: ptrs and coord are on all processors, but voldata is only on Proc 0

    int nrays, nnz;

    ierr = getnnz(volsize, radius, origin, &nrays, &nnz);

    int nnzloc, nraysloc;
    int * ptrs = new int[nrays+1];
    int * cord = new int[3*nrays];
    int *nnzpart = new int[dims[COL]];
    int *nnzbase = new int[dims[COL]+1]; 
    int *ptrstart = new int[dims[COL]+1];

    ierr = getcb2sph(volsize, radius, origin, nnz, ptrs, cord);
    nnzloc = setpart_gr1(comm_2d, nnz, nnzpart, nnzbase);
    nraysloc = sphpart(comm_2d, nrays, ptrs, nnzbase, ptrstart);
  
    int myptrstart = ptrstart[mycoords[COL]];
    int nnzall[dims[COL]];
    for (int i = 0; i<dims[COL]; i++)
      nnzall[i] = ptrs[ptrstart[i+1]] - ptrs[ptrstart[i]];
  
    nnzloc = nnzall[mycoords[COL]];

    float *bvol_loc = new float[nnzloc];
    float *bvol = new float[nnzloc];   
    float *xvol_sphloc = new float[nnzloc];
    float *pxvol_loc = new float[nnzloc];
    float *pxvol = new float[nnzloc];
    float * grad_loc = new float[nnzloc];
    for (int i=0; i< nnzloc; i++){
      xvol_sphloc[i] = 0.0;
      bvol[i] = 0.0;
      bvol_loc[i] = 0.0;
      pxvol_loc[i] = 0.0;
      pxvol[i] = 0.0;
      grad_loc[i] = 0.0;
}

    // this is just to set ptrs and cord, voldata is all 0.0 at this point
//    ierr = cb2sph(voldata, volsize, radius, origin, nnz, ptrs, cord, xvol_sph); 

    // arrays to hold volume data:
    // initial backprojected volume, local
    //float * bvol_loc = new float[nnz];
    // initial backprojected volume, global
    //float * bvol = new float[nnz];
    // P^T * P * xvol, local
    //float * pxvol = new float[nnz];
    // P^T * P * xvol, global
    //float * pxvol_loc = new float[nnz];
    
   /* for ( int i = 0 ; i < nnz ; ++i ) {
	xvol_sph[i] = 0.0;
	bvol_loc[i] = 0.0;
	bvol[i] = 0.0;
	pxvol[i] = 0.0;
	pxvol_loc[i] = 0.0;
    }*/

    EMData * current_image;
    float phi, theta, psi;
    Transform3D RA;
    Transform3D Tf;
    nsym = Tf.get_nsym(symmetry);
    Transform3D::EulerType EULER_SPIDER = Transform3D::SPIDER;
    Dict angdict;
//printf("nsym = %d\n", nsym);
    int iter = 1;
	
    double rnorm = 0.0, rnorm_loc = 0.0;
    double bnorm = 0.0, bnorm_loc = 0.0;
    //float * grad = new float[nnz];
    float * image_data;
    float * projected_data_loc = new float[nangloc*nx*nx];
    float * projected_data = new float[nangloc*nx*nx];

    float dm[8];
    
    int restarts = 0;
    
    while (iter <= maxit) {
	if ( iter == 1 ) {
	    if ( restarts == 0 ) {
		// only do this if we aren't restarting due to lam being too large
		for ( int i = 0 ; i < nangloc ; ++i ) {
		    current_image = images[i];
		    image_data = current_image->get_data(); 
		    // retrieve the angles and shifts associated with each image from the array
		    // angleshift.
		    phi   = angleshift[5*i + 0];
		    theta = angleshift[5*i + 1];
		    psi   = angleshift[5*i + 2];
		    dm[6] = angleshift[5*i + 3] * -1.0;
		    dm[7] = angleshift[5*i + 4] * -1.0;
		    // make an instance of Transform3D with the angles
		    RA = Transform3D(EULER_SPIDER, phi, theta, psi);
		    for ( int ns = 1 ; ns < nsym + 1 ; ++ns ) {
			// compose it with each symmetry rotation in turn
			// shifts (stored in dm[6] and dm[7] remain fixed
			Tf = Tf.get_sym(symmetry, ns) * RA;
			angdict = Tf.get_rotation(EULER_SPIDER);
			phi   = (float) angdict["phi"]   * PI/180.0;
			theta = (float) angdict["theta"] * PI/180.0;
			psi   = (float) angdict["psi"]   * PI/180.0;
			make_proj_mat(phi, theta, psi, dm); 
			// accumulate the back-projected images in bvol_loc
			ierr = bckpj3_Cart(volsize, nraysloc, nnzloc, dm, origin, radius, ptrs, cord, myptrstart, image_data, bvol_loc);
		    }
		}
		// reduce bvol_loc so each processor has the full volume
		//ierr = MPI_Allreduce(bvol_loc, bvol, nnz, MPI_FLOAT, MPI_SUM, comm);
                // Now an all reduce along the columns
                ierr = MPI_Allreduce (bvol_loc, bvol, nnzloc, MPI_FLOAT, MPI_SUM, comm_col);

	    }

	    // calculate the norm of the backprojected volume
            bnorm_loc = 0.0;
	    for ( int j = 0 ; j < nnzloc ; ++j ) {
		bnorm_loc += bvol[j] * (double) bvol[j];
		grad_loc[j]= bvol[j];
	    }
            ierr = MPI_Allreduce (&bnorm_loc, &bnorm, 1, MPI_DOUBLE, MPI_SUM, comm_row);

	    bnorm /= nnz;
	    bnorm = sqrt(bnorm);
	  //  printf("bnorm = %f\n", bnorm);
	} else {

          for (int i=0; i<nangloc*nx*nx; i++){
	     projected_data_loc[i] = 0.0;
	     projected_data[i] = 0.0;
          }

	    for ( int i = 0 ; i < nangloc ; ++i ) {
		// retrieve the angles and shifts from angleshift
		RA = Transform3D(EULER_SPIDER, angleshift[5*i + 0], angleshift[5*i + 1], angleshift[5*i + 2]);
		dm[6] = angleshift[5*i + 3] * -1.0;
		dm[7] = angleshift[5*i + 4] * -1.0;
		for ( int ns = 1 ; ns < nsym + 1 ; ++ns ) {
		    // iterate over symmetries
		    Tf = Tf.get_sym(symmetry, ns) * RA;
		    angdict = Tf.get_rotation(EULER_SPIDER);

		    phi   = (float) angdict["phi"]   * PI/180.0;
		    theta = (float) angdict["theta"] * PI/180.0;
		    psi   = (float) angdict["psi"]   * PI/180.0;
		    make_proj_mat(phi, theta, psi, dm); 
		    // accumulate P^TPxvol in pxvol_loc
		    ierr = fwdpj3_Cart(volsize, nraysloc, nnzloc, dm, origin, radius, ptrs, cord, myptrstart, xvol_sphloc, &projected_data_loc[nx*nx*i]);
		}
	    }
	    // and reduce the accumulated pxvol_loc's
	      // Now an all reduce along the rows
            ierr = MPI_Allreduce(projected_data_loc, projected_data, nangloc*nx*nx, MPI_FLOAT, MPI_SUM, comm_row);

	    for ( int i = 0 ; i < nangloc ; ++i ) {
		// retrieve the angles and shifts from angleshift
		RA = Transform3D(EULER_SPIDER, angleshift[5*i + 0], angleshift[5*i + 1], angleshift[5*i + 2]);
		dm[6] = angleshift[5*i + 3] * -1.0;
		dm[7] = angleshift[5*i + 4] * -1.0;
		for ( int ns = 1 ; ns < nsym + 1 ; ++ns ) {
		    // iterate over symmetries
		    Tf = Tf.get_sym(symmetry, ns) * RA;
		    angdict = Tf.get_rotation(EULER_SPIDER);
		    // reset the array in which projected data are stored
		    phi   = (float) angdict["phi"]   * PI/180.0;
		    theta = (float) angdict["theta"] * PI/180.0;
		    psi   = (float) angdict["psi"]   * PI/180.0;
		    make_proj_mat(phi, theta, psi, dm); 
		    // accumulate P^TPxvol in pxvol_loc
		    ierr = bckpj3_Cart(volsize, nraysloc, nnzloc, dm, origin, radius, ptrs, cord, myptrstart, &projected_data[nx*nx*i], pxvol_loc);
		}
	    }
    
          ierr = MPI_Allreduce(pxvol_loc, pxvol, nnzloc, MPI_FLOAT, MPI_SUM, comm_col);
	    
	    for ( int j = 0 ; j < nnzloc ; ++j ) {
		grad_loc[j] = bvol[j];
		grad_loc[j] -= pxvol[j];
	    }
	}

        rnorm_loc = 0.0;
	for ( int j = 0 ; j < nnzloc ; ++j ) {
	    rnorm_loc += grad_loc[j]* (double) grad_loc[j];
	}
   //    printf("I am (%d, %d), with rnorm = %f\n", mycoords[ROW], mycoords[COL], rnorm);
        ierr = MPI_Allreduce (&rnorm_loc, &rnorm, 1, MPI_DOUBLE, MPI_SUM, comm_row);
 	rnorm /= nnz;
	rnorm = sqrt(rnorm);
	if ( my2dpid == 0 ) printf("iter = %3d, rnorm / bnorm = %6.3f, rnorm = %6.3f\n", iter, rnorm / bnorm, rnorm);
	// if on the second pass, rnorm is greater than bnorm, lam is probably set too high
	// reduce it by a factor of 2 and start over
	//	if ( iter == 2 && rnorm / bnorm > old_rnorm ) {
	if ( rnorm / bnorm > old_rnorm ) {
	    // but don't do it more than 20 times
	    if ( restarts > 20 ) {
		if ( my2dpid == 0 ) printf("Failure to converge, even with lam = %f\n", lam);
		break;
	    } else {
		++restarts;
		iter = 1;
		lam /= 2.0;
		// reset these 
		old_rnorm = 1.0001; // kluge, make sure if its 1.0 + epsilon it still works
		for ( int j = 0 ; j < nnzloc ; ++j ) {
		    xvol_sphloc[j]  = 0.0;
		    pxvol_loc[j] = 0.0; 
		}
		if ( my2dpid == 0 ) printf("reducing lam to %f, restarting\n", lam);
		continue;
	    }
	}
	// if changes are sufficiently small, or if no further progress is made, terminate
	if ( rnorm / bnorm < tol || rnorm / bnorm > old_rnorm ) {
	    if ( my2dpid == 0 ) printf("Terminating with rnorm/bnorm = %f, tol = %f, old_rnorm = %f\n",rnorm/bnorm, tol, old_rnorm);
	    break;
	}
	// update the termination threshold
	old_rnorm = rnorm / bnorm;
	// update the reconstructed volume

	for ( int j = 0 ; j < nnzloc ; ++j ) {
	    xvol_sphloc[j] += lam * grad_loc[j];
 	    pxvol_loc[j] = 0.0; // reset it so it's ready to accumulate for the next iteration
	}

	++iter;
    }

// Bring all parts of the spherical volume back together and turn it into cube format.
if (mycoords[ROW] == 0 && mycoords[COL] != 0 ){
  MPI_Send(xvol_sphloc, nnzloc, MPI_FLOAT, 0, my2dpid, comm_2d);
}

if (mycoords[ROW] == 0 && mycoords[COL] == 0 ){

  xvol->set_size(nx, nx, nx);
  xvol->to_zero();
  float * voldata = xvol->get_data();
  float * xvol_sph = new float[nnz];

  for(int i=0; i<nnzloc; i++)
    xvol_sph[i] = xvol_sphloc[i];

  for(int i=1; i< dims[COL]; i++){
    srcoords[ROW] = 0;
    srcoords[COL] = i;
    MPI_Cart_rank(comm_2d, srcoords, &srpid);

    MPI_Recv(&xvol_sph[ptrs[ptrstart[i]]], nnzall[i] , MPI_FLOAT, srpid, srpid, comm_2d, &mpistatus);
  }

    // unpack the spherical volume back out into the original EMData object
    ierr = sph2cb(xvol_sph, volsize, nrays, radius, nnz, ptrs, cord, voldata);
    EMDeleteArray(xvol_sph);
}
    EMDeleteArray(grad_loc);
    EMDeleteArray(pxvol_loc);
    EMDeleteArray(bvol_loc);

    EMDeleteArray(ptrs);
    EMDeleteArray(cord);
    //EMDeleteArray(projected_data);
    
    EMDeleteArray(psize);
    EMDeleteArray(nbase);

    delete [] projected_data_loc;

    return 0; // recons3d_sirt_mpi
}

// #define PROJECTOR_TYPE "chao"

// int main(int argc, char ** argv)
// {
//     MPI_Comm comm = MPI_COMM_WORLD;
//     int ncpus, mypid, ierr;
//     int nloc; 

//     MPI_Status mpistatus;
//     MPI_Init(&argc,&argv);
//     MPI_Comm_size(comm,&ncpus);
//     MPI_Comm_rank(comm,&mypid);
//     printf("mypid = %d, ncpus = %d\n", mypid, ncpus);

//     EMData * xvol = new EMData();

//     EMData *volume = new EMData();
//     if ( mypid == 0 ) {
// 	volume->read_image("vol001.tfc");
//     }
//     int nx = volume->get_xsize();
//     ierr = MPI_Bcast(&nx, 1, MPI_INT, 0, comm);
//     float dtheta = 5.0;
//     std::vector<float> ref_angles = Util::even_angles(dtheta);

//     int nangles = ref_angles.size() / 3;
//     if ( mypid == 0 ) printf("Using %d images\n", nangles);
//     int radius = 30;
//     int maxit = 50;

//     int * psize;
//     int * nbase;
//     int nangloc;
//     psize = new int[ncpus];
//     nbase = new int[ncpus];
//     nangloc = setpart(comm, nangles, psize, nbase);

//     EMData ** images = (EMData **) malloc(nangloc * sizeof(EMData *));
//     std::vector<float> anglelist(3,0.0);
//     Dict volparams;
//     volparams["angletype"] = "SPIDER";

//     float * imgdata;
//     EMData * send_projection;
//     volparams["radius"] = radius;

//     float * angleshift = new float[5 * nloc];

//     if ( mypid == 0 ) {
// 	for ( int i = 0 ; i < nangloc ; ++i ) {	
// 	    anglelist[0] = ref_angles[i*3 + 0];
// 	    anglelist[1] = ref_angles[i*3 + 1];
// 	    anglelist[2] = ref_angles[i*3 + 2];
// 	    volparams["anglelist"] = anglelist;
// 	    images[i] = volume->project(PROJECTOR_TYPE, volparams);
// // 	    images[i]->set_attr("phi", anglelist[0]);
// // 	    images[i]->set_attr("theta", anglelist[1]);
// // 	    images[i]->set_attr("psi", anglelist[2]);
// // 	    images[i]->set_attr("s2x", 0.0);
// // 	    images[i]->set_attr("s2y", 0.0);
// 	    angleshift[5*i + 0] = anglelist[0];
// 	    angleshift[5*i + 1] = anglelist[1];
// 	    angleshift[5*i + 2] = anglelist[2];
// 	    angleshift[5*i + 3] = 0.0;
// 	    angleshift[5*i + 4] = 0.0;
// 	}
// 	for ( int np = 1 ; np < ncpus ; ++np ) {
// 	    for ( int i = 0 ; i < psize[np] ; ++i ) {	
// 		anglelist[0] = ref_angles[(nbase[np] + i)*3 + 0];
// 		anglelist[1] = ref_angles[(nbase[np] + i)*3 + 1];
// 		anglelist[2] = ref_angles[(nbase[np] + i)*3 + 2];
// 		volparams["anglelist"] = anglelist;
// 		MPI_Send(&anglelist[0], 3, MPI_FLOAT, np, np, comm);
// 		send_projection = volume->project(PROJECTOR_TYPE, volparams);
// 		imgdata = send_projection->get_data();
// 		MPI_Send(imgdata, nx*nx, MPI_FLOAT, np, np, comm);
// 		EMDeletePtr(send_projection);
// 	    }
		
// 	}
//     } else { // mypid != 0
// 	for ( int i = 0 ; i < nangloc ; ++i ) {
// 	    images[i] = new EMData();
// 	    images[i]->set_size(nx,nx);
// 	    MPI_Recv(&anglelist[0], 3, MPI_FLOAT, 0, mypid, comm, &mpistatus);
// // 	    images[i]->set_attr("phi", anglelist[0]);
// // 	    images[i]->set_attr("theta", anglelist[1]);
// // 	    images[i]->set_attr("psi", anglelist[2]);
// // 	    images[i]->set_attr("s2x", 0.0);
// // 	    images[i]->set_attr("s2y", 0.0);
// 	    angleshift[5*i + 0] = anglelist[0];
// 	    angleshift[5*i + 1] = anglelist[1];
// 	    angleshift[5*i + 2] = anglelist[2];
// 	    angleshift[5*i + 3] = 0.0;
// 	    angleshift[5*i + 4] = 0.0;
// 	    imgdata = images[i]->get_data();
// 	    MPI_Recv(imgdata, nx*nx, MPI_FLOAT, 0, mypid, comm, &mpistatus);
// 	}
//     }
//     if ( mypid == 0 ) printf("Done with data distribution\n");
	
//     float lam = 1.0e-4;
//     float tol = 1.0e-3;
//     std::string symmetry = "c1";

//     recons3d_sirt_mpi(comm, images, angleshift, xvol, nangloc, radius, lam, maxit, symmetry, tol);
//     if ( mypid == 0 ) printf("Done with SIRT\n");

//     EMUtil::ImageType WRITE_SPI = EMUtil::IMAGE_SINGLE_SPIDER;
//     if ( mypid == 0 ) {
// 	xvol->write_image("volume_sirt.spi", 0, WRITE_SPI);
//     }


//     for ( int i = 0 ; i < nangloc ; ++i ) {	
// 	EMDeletePtr(images[i]);
//     }		
//     free(images);
//     EMDeletePtr(xvol);
//     EMDeletePtr(volume);
//     EMDeleteArray(angleshift);

//     ierr = MPI_Finalize();

//     return 0; // main
// }
