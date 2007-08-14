#include "mpi.h"

#include "emdata.h"

#include "project3d.h"
#include "sirt.h"
#include "utilcomm.h"

using namespace EMAN;

int recons3d_sirt_mpi(MPI_Comm comm , EMData ** images, float * angleshift  , 
                      EMData *& xvol, int nangloc     , int radius          , 
                      float lam     , int maxit       , std::string symmetry, 
                      float tol)
{
    int ncpus, mypid, ierr;
    double t0;

    MPI_Status mpistatus;
    MPI_Comm_size(comm,&ncpus);
    MPI_Comm_rank(comm,&mypid);
    
    int * psize;
    int * nbase;

    int nangles;
    psize = new int[ncpus];
    nbase = new int[ncpus];
    MPI_Allreduce(&nangloc, &nangles, 1, MPI_INT, MPI_SUM, comm);
    

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

    // this is not currently needed, because the stack that gets passed to sirt 
    // will have its background subtracted already
    // ierr = CleanStack(comm, images, nangloc, radius, volsize, origin);

    xvol->set_size(nx, nx, nx);
    xvol->to_zero();
    float * voldata = xvol->get_data();
	
    // vector of symmetrized angles
    std::vector<float> symangles(3,0.0); 
	
    // kluge, make sure if its 1.0 + epsilon it still works;
    float old_rnorm = 1.00001; 
	
    int nrays, nnz;

    ierr = getnnz(volsize, radius, origin, &nrays, &nnz);

    int * ptrs = new int[nrays+1];
    int * cord = new int[3*nrays];

    float * xvol_sph = new float[nnz];

    // this is just to set ptrs and cord, voldata is all 0.0 at this point
    ierr = cb2sph(voldata, volsize, radius, origin, nnz, ptrs, cord, xvol_sph); 

    // arrays to hold volume data:
    // initial backprojected volume, local
    float * bvol_loc = new float[nnz];
    // initial backprojected volume, global
    float * bvol = new float[nnz];
    // P^T * P * xvol, local
    float * pxvol = new float[nnz];
    // P^T * P * xvol, global
    float * pxvol_loc = new float[nnz];
    
    for ( int i = 0 ; i < nnz ; ++i ) {
	xvol_sph[i] = 0.0;
	bvol_loc[i] = 0.0;
	bvol[i] = 0.0;
	pxvol[i] = 0.0;
	pxvol_loc[i] = 0.0;
    }

    EMData * current_image;
    float phi, theta, psi;
    Transform3D RA;
    Transform3D Tf;
    nsym = Tf.get_nsym(symmetry);
    Transform3D::EulerType EULER_SPIDER = Transform3D::SPIDER;
    Dict angdict;

    int iter = 1;
	
    double rnorm = 0.0;
    double bnorm = 0.0;
    float * grad = new float[nnz];
    float * image_data;
    float * projected_data = new float[nx*nx];
    float dm[8];
    
    int restarts = 0;

    t0 = MPI_Wtime();    
    while (iter <= maxit) {
	if ( iter == 1 ) {
	    if ( restarts == 0 ) {
		// only do this if we aren't restarting due to lam being too large
		for ( int i = 0 ; i < nangloc ; ++i ) {
		    current_image = images[i];
		    image_data = current_image->get_data(); 
		    // retrieve the angles and shifts associated with each image 
                    // from the array angleshift.
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
			ierr = bckpj3(volsize, nrays, nnz, dm, origin, radius, ptrs, cord, 
				      image_data, bvol_loc);
		    }
		}
		// reduce bvol_loc so each processor has the full volume
		ierr = MPI_Allreduce(bvol_loc, bvol, nnz, MPI_FLOAT, MPI_SUM, comm);

	    }

	    // calculate the norm of the backprojected volume
	    for ( int j = 0 ; j < nnz ; ++j ) {
		bnorm += bvol[j] * (double) bvol[j];
		grad[j] = bvol[j];
	    }
	    bnorm /= nnz;
	    bnorm = sqrt(bnorm);
	    
	} else {
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
		for ( int ns = 1 ; ns < nsym + 1 ; ++ns ) {
		    // iterate over symmetries
		    Tf = Tf.get_sym(symmetry, ns) * RA;
		    angdict = Tf.get_rotation(EULER_SPIDER);
		    // reset the array in which projected data are stored
		    for ( int j = 0 ; j < nx*nx ; ++j ) {
			projected_data[j] = 0.0;
		    }
		    phi   = (float) angdict["phi"]   * PI/180.0;
		    theta = (float) angdict["theta"] * PI/180.0;
		    psi   = (float) angdict["psi"]   * PI/180.0;
		    make_proj_mat(phi, theta, psi, dm); 
		    // accumulate P^TPxvol in pxvol_loc
		    ierr = fwdpj3(volsize, nrays, nnz, dm, origin, radius, ptrs, 
                                  cord, xvol_sph, projected_data);
		    ierr = bckpj3(volsize, nrays, nnz, dm, origin, radius, ptrs, 
                                  cord, projected_data, pxvol_loc);
		}
	    }
	    // and reduce the accumulated pxvol_loc's
	    ierr = MPI_Allreduce(pxvol_loc, pxvol, nnz, MPI_FLOAT, MPI_SUM, comm);
	    
	    for ( int j = 0 ; j < nnz ; ++j ) {
		grad[j] = bvol[j];
		grad[j] -= pxvol[j];
	    }
	}
	rnorm = 0.0;
	for ( int j = 0 ; j < nnz ; ++j ) {
	    rnorm += grad[j]* (double) grad[j];
	}
	rnorm /= nnz;
	rnorm = sqrt(rnorm);
	if ( mypid == 0 ) printf("iter = %3d, rnorm / bnorm = %11.3e, rnorm = %11.3e\n", 
                                 iter, rnorm / bnorm, rnorm);
	// if on the second pass, rnorm is greater than bnorm, 
        // lam is probably set too high reduce it by a factor of 2 and start over
	//	if ( iter == 2 && rnorm / bnorm > old_rnorm ) {
	if ( rnorm / bnorm > old_rnorm ) {
	    // but don't do it more than 20 times
	    if ( restarts > 20 ) {
		if ( mypid == 0 ) 
                   printf("Failure to converge, even with lam = %11.3e\n", lam);
		break;
	    } else {
		++restarts;
		iter = 1;
		lam /= 2.0;
		// reset these 
                // kluge, make sure if its 1.0 + epsilon it still works
		old_rnorm = 1.0001; 
		for ( int j = 0 ; j < nnz ; ++j ) {
		    xvol_sph[j]  = 0.0;
		    pxvol_loc[j] = 0.0; 
		}
		if ( mypid == 0 ) printf("reducing lam to %11.3e, restarting\n", lam);
		continue;
	    }
	}
	// if changes are sufficiently small, or if no further progress is made, terminate
	if ( rnorm / bnorm < tol || rnorm / bnorm > old_rnorm ) {
	    if ( mypid == 0 ) 
               printf("Terminating with rnorm/bnorm = %11.3e, tol = %11.3e, ");
               printf("old_rnorm = %11.3e\n", rnorm/bnorm, tol, old_rnorm);
	    break;
	}
	// update the termination threshold
	old_rnorm = rnorm / bnorm;
	// update the reconstructed volume
	for ( int j = 0 ; j < nnz ; ++j ) {
	    xvol_sph[j] += lam * grad[j];
            // reset it so it's ready to accumulate for the next iteration
 	    pxvol_loc[j] = 0.0; 
	}

	++iter;
    }
    if (mypid == 0) printf("Total time in SIRT = %11.3e\n", MPI_Wtime()-t0);

    // unpack the spherical volume back out into the original EMData object
    ierr = sph2cb(xvol_sph, volsize, nrays, radius, nnz, ptrs, cord, voldata);

    EMDeleteArray(grad);
    EMDeleteArray(pxvol);
    EMDeleteArray(bvol);
    EMDeleteArray(pxvol_loc);
    EMDeleteArray(bvol_loc);

    EMDeleteArray(ptrs);
    EMDeleteArray(cord);
    EMDeleteArray(xvol_sph);
    EMDeleteArray(projected_data);
    
    EMDeleteArray(psize);
    EMDeleteArray(nbase);

    return 0; // recons3d_sirt_mpi
}

