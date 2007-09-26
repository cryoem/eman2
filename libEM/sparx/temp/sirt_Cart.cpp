#include "mpi.h"

#include "emdata.h"

#include "sirt_Cart.h"
#include "utilcomm_Cart.h"
#include "project3d_Cart.h"
#include "project3d.h"

using namespace EMAN;
#define ROW  0
#define COL  1

//int CleanStack(MPI_Comm comm, EMData ** image_stack, int nloc, int ri, Vec3i volsize, Vec3i origin);

int recons3d_sirt_mpi_Cart(MPI_Comm comm_2d    , MPI_Comm comm_row, 
                           MPI_Comm comm_col   , EMData ** images , 
                           float * angleshift  , EMData *& xvol   , 
                           int nangloc         , int radius       , 
                           float lam           , int maxit        , 
                           std::string symmetry, float tol)
{
    MPI_Status mpistatus;
    int ncpus, my2dpid, ierr;
    int mycoords[2], dims[2], periods[2];
    int srpid, srcoords[2]; // Send/receive info
    double t0;

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

    // make radius as large as possible if the user didn't provide one
    if ( radius == -1 ) radius = nx/2 - 1; 
    
    Vec3i volsize, origin;
    volsize[0] = nx;
    volsize[1] = nx;
    volsize[2] = nx;
    origin[0] = nx/2+1;
    origin[1] = nx/2+1;
    origin[2] = nx/2+1;

    // this is not currently needed, because the stack that gets passed to 
    // sirt will have its background subtracted already
    // ierr = CleanStack(comm, images, nangloc, radius, volsize, origin);

    // vector of symmetrized angles
    std::vector<float> symangles(3,0.0); 
	
    // kluge, make sure if its 1.0 + epsilon it still works;
    float old_rnorm = 1.00001; 
	

    // Now distribute the volume (in spherical format) among columns of 
    // processors and use nnz to determine the splitting.  
    // Note: ptrs and coord are on all processors, but voldata is only on 
    // Proc 0

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
   
    t0 = MPI_Wtime();
 
    while (iter <= maxit) {
	if ( iter == 1 ) {
	    if ( restarts == 0 ) {
		// only do this if we aren't restarting due to lam being 
                // too large
		for ( int i = 0 ; i < nangloc ; ++i ) {
		    current_image = images[i];
		    image_data = current_image->get_data(); 
		    // retrieve the angles and shifts associated with
                    // each image from the array angleshift.
		    phi   = angleshift[5*i + 0];
		    theta = angleshift[5*i + 1];
		    psi   = angleshift[5*i + 2];

                    // need to change signs here because the input shifts
                    // are shifts associated with 2-D images. Because
                    // the projection operator actually shifts the volume
                    // the signs should be negated here
		    dm[6] = -angleshift[5*i + 3];
		    dm[7] = -angleshift[5*i + 4];
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
			ierr = bckpj3_Cart(volsize, nraysloc, nnzloc, dm, 
                                           origin, radius, ptrs, cord, 
                                           myptrstart, image_data, bvol_loc);
		    }
		}
		// reduce bvol_loc so each processor has the full volume
		//ierr = MPI_Allreduce(bvol_loc, bvol, nnz, MPI_FLOAT, MPI_SUM,
                //                     comm);
                // Now an all reduce along the columns
                ierr = MPI_Allreduce (bvol_loc, bvol, nnzloc, MPI_FLOAT, 
                                      MPI_SUM, comm_col);

	    }

	    // calculate the norm of the backprojected volume
            bnorm_loc = 0.0;
	    for ( int j = 0 ; j < nnzloc ; ++j ) {
		bnorm_loc += bvol[j] * (double) bvol[j];
		grad_loc[j]= bvol[j];
	    }
            ierr = MPI_Allreduce (&bnorm_loc, &bnorm, 1, MPI_DOUBLE, MPI_SUM,
                                  comm_row);

	    bnorm /= nnz;
	    bnorm = sqrt(bnorm);
	    // printf("bnorm = %f\n", bnorm);
	} else {

            // iterate over symmetry
            for ( int ns = 1 ; ns < nsym + 1 ; ++ns ) {
 	       // reset local and global 2-D projections to zeros
               for (int i=0; i<nangloc*nx*nx; i++){
                  projected_data_loc[i] = 0.0;
                  projected_data[i] = 0.0;
               }
   
               // project from 3-D to 2-D
	       for ( int i = 0 ; i < nangloc ; ++i ) {
                  // retrieve the angles and shifts from angleshift
                  RA = Transform3D(EULER_SPIDER, angleshift[5*i + 0], 
                                angleshift[5*i + 1], angleshift[5*i + 2]);

                  // need to change signs here because the input shifts
                  // are shifts associated with 2-D images. Because
                  // the projection operator actually shifts the volume
                  // the signs should be negated here
                  dm[6] = -angleshift[5*i + 3];
                  dm[7] = -angleshift[5*i + 4];
   
                  // apply symmetry transformation
                  Tf = Tf.get_sym(symmetry, ns) * RA;
                  angdict = Tf.get_rotation(EULER_SPIDER);
   
                  phi   = (float) angdict["phi"]   * PI/180.0;
                  theta = (float) angdict["theta"] * PI/180.0;
                  psi   = (float) angdict["psi"]   * PI/180.0;
                  make_proj_mat(phi, theta, psi, dm); 
   
                  ierr = fwdpj3_Cart(volsize, nraysloc, nnzloc, dm, 
                                     origin, radius, ptrs, cord, myptrstart, 
                                     xvol_sphloc, &projected_data_loc[nx*nx*i]);
               }


               // Now perform global sum on 2-D images along the rows
               ierr = MPI_Allreduce(projected_data_loc, projected_data, 
                                    nangloc*nx*nx, MPI_FLOAT, MPI_SUM, comm_row);

               // backproject from 2-D to 3-D
               for ( int i = 0 ; i < nangloc ; ++i ) {
                  // retrieve the angles and shifts from angleshift
                  RA = Transform3D(EULER_SPIDER, angleshift[5*i + 0], 
                                   angleshift[5*i + 1], angleshift[5*i + 2]);

                  // need to change signs here because the input shifts
                  // are shifts associated with 2-D images. Because
                  // the projection operator actually shifts the volume
                  // the signs should be negated here
                  dm[6] = -angleshift[5*i + 3];
                  dm[7] = -angleshift[5*i + 4];

                  // apply symmetry transformation
                  Tf = Tf.get_sym(symmetry, ns) * RA;
                  angdict = Tf.get_rotation(EULER_SPIDER);

                  // reset the array in which projected data are stored
                  phi   = (float) angdict["phi"]   * PI/180.0;
                  theta = (float) angdict["theta"] * PI/180.0;
                  psi   = (float) angdict["psi"]   * PI/180.0;
                  make_proj_mat(phi, theta, psi, dm); 

                  // accumulate P^TPxvol in pxvol_loc
                  ierr = bckpj3_Cart(volsize, nraysloc, nnzloc, dm, origin, 
                                     radius, ptrs, cord, myptrstart, 
                                     &projected_data[nx*nx*i], pxvol_loc);
	       }
	    } // endfor (ns=1,nsym)

            ierr = MPI_Allreduce(pxvol_loc, pxvol, nnzloc, MPI_FLOAT, MPI_SUM, 
                                 comm_col);

	    for ( int j = 0 ; j < nnzloc ; ++j ) {
		grad_loc[j] = bvol[j];
		grad_loc[j] -= pxvol[j];
	    }
	}

        rnorm_loc = 0.0;
	for ( int j = 0 ; j < nnzloc ; ++j ) {
	    rnorm_loc += grad_loc[j]* (double) grad_loc[j];
	}
        // printf("I am (%d, %d), with rnorm = %11.3e\n", 
        // mycoords[ROW], mycoords[COL], rnorm);
        ierr = MPI_Allreduce (&rnorm_loc, &rnorm, 1, MPI_DOUBLE, MPI_SUM, 
                              comm_row);
 	rnorm /= nnz;
	rnorm = sqrt(rnorm);
	if ( my2dpid == 0 ) 
            printf("iter = %3d, rnorm / bnorm = %11.3e, rnorm = %11.3e\n", 
                   iter, rnorm / bnorm, rnorm);
        // if on the second pass, rnorm is greater than bnorm, 
        // lam is probably set too high reduce it by a factor of 2 
        // and start over
        //	if ( iter == 2 && rnorm / bnorm > old_rnorm ) {
        if ( rnorm / bnorm > old_rnorm ) {
           // but don't do it more than 20 times
	   if ( restarts > 20 ) {
               if ( my2dpid == 0 ) 
                  printf("Failure to converge, even with lam = %f\n", lam);
                  break;
           } else {
               ++restarts;
               iter = 1;
               lam /= 2.0;
               // reset these 
               // kluge, make sure if its 1.0 + epsilon it still works
               old_rnorm = 1.0001; 
               for ( int j = 0 ; j < nnzloc ; ++j ) {
                  xvol_sphloc[j]  = 0.0;
                  pxvol_loc[j] = 0.0; 
               }
               if ( my2dpid == 0 ) 
                  printf("reducing lam to %11.3e, restarting\n", lam);
               continue;
           } // endif (restarts)
        } // endif (rnorm/bnorm) 

        // if changes are sufficiently small, 
        // or if no further progress is made, terminate
	if ( rnorm / bnorm < tol || rnorm / bnorm > old_rnorm ) {
	   if ( my2dpid == 0 ) 
              printf("Terminating with rnorm/bnorm = %11.3e, ");
              printf("tol = %11.3e, old_rnorm = %11.3e\n",
                     rnorm/bnorm, tol, old_rnorm);
	   break;
	}
	// update the termination threshold
	old_rnorm = rnorm / bnorm;
	// update the reconstructed volume

	for ( int j = 0 ; j < nnzloc ; ++j ) {
	    xvol_sphloc[j] += lam * grad_loc[j];
            // reset it so it's ready to accumulate for the next iteration
 	    pxvol_loc[j] = 0.0; 
	}
	++iter;
    }

    if (my2dpid == 0) printf("Total time in SIRT = %11.3e\n", MPI_Wtime()-t0);

    // Bring all parts of the spherical volume back together and turn it 
    // into cube format.
    if (mycoords[ROW] == 0 && mycoords[COL] != 0 ){
       srcoords[ROW] = 0;
       srcoords[COL] = 0;
       MPI_Cart_rank(comm_2d, srcoords, &srpid);
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
          MPI_Recv(&xvol_sph[ptrs[ptrstart[i]]-1], nnzall[i], MPI_FLOAT, 
                   srpid, srpid, comm_2d, &mpistatus);
       }

       // unpack the spherical volume back out into the original EMData object
       ierr = sph2cb(xvol_sph, volsize, nrays, radius, nnz, ptrs, cord, voldata);
       EMDeleteArray(xvol_sph);
    } // endif (mycord[ROW]== 0...)
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

