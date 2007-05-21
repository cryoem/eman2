#include "mpi.h"

#include "emdata.h"

#include "ali3d_unified_mpi.h"
#include "sirt.h"

#define PROJECTOR_TYPE "chao"
#define PI 3.14159265358979
using namespace EMAN;


int CleanStack(MPI_Comm comm, EMData ** image_stack, int nloc, int ri, Vec3i volsize, Vec3i origin);

int recons3d_sirt_mpi(MPI_Comm comm, EMData ** images, EMData *& xvol, int nangloc, int radius, float lam, int maxit, std::string symmetry, float tol)
{
    int ncpus, mypid, ierr;

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
    if ( radius == -1 ) radius = nx/2 - 1; // make radius as large as possible if the user didn't provide one
    
    Vec3i volsize, origin;
    volsize[0] = nx;
    volsize[1] = nx;
    volsize[2] = nx;
    origin[0] = nx/2+1;
    origin[1] = nx/2+1;
    origin[2] = nx/2+1;
    ierr = CleanStack(comm, images, nangloc, radius, volsize, origin);

    xvol->set_size(nx, nx, nx);
    xvol->to_zero();
    float * voldata = xvol->get_data();
	
	
    // vector of symmetrized angles
    std::vector<float> symangles(3,0.0); 
    std::vector<float> angles;
    std::vector<float> shifts;
	
    float old_rnorm = 1.0;
	
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
    
    while (iter <= maxit) {
	if ( iter == 1 ) {
	    if ( restarts == 0 ) {
		// only do this if we aren't restarting due to lam being too large
		for ( int i = 0 ; i < nangloc ; ++i ) {
		    current_image = images[i];
		    image_data = current_image->get_data(); 
		    // store the angles and shifts associated with each image in the vectors
		    // angles and shifts.
		    phi   = current_image->get_attr("phi"); 
		    angles.push_back(phi);
		    theta = current_image->get_attr("theta"); 
		    angles.push_back(theta);
		    psi   = current_image->get_attr("psi"); 
		    angles.push_back(psi);
		    dm[6] = current_image->get_attr("s2x");
		    dm[6] *= -1.0;
		    shifts.push_back(dm[6]);
		    dm[7] = current_image->get_attr("s2y");
		    dm[7] *= -1.0;
		    shifts.push_back(dm[7]);
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
		// retrieve the angles and shifts from the vectors that store them
		RA = Transform3D(EULER_SPIDER, angles[3*i + 0], angles[3*i + 1], angles[3*i + 2]);
		dm[6] = shifts[2*i + 0];
		dm[7] = shifts[2*i + 1];
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
		    ierr = fwdpj3(volsize, nrays, nnz, dm, origin, radius, ptrs, cord, xvol_sph, projected_data);
		    ierr = bckpj3(volsize, nrays, nnz, dm, origin, radius, ptrs, cord, 
				  projected_data, pxvol_loc);
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
	if ( mypid == 0 ) printf("iter = %3d, rnorm / bnorm = %6.3f, rnorm = %6.3f\n", iter, rnorm / bnorm, rnorm);
	// if on the second pass, rnorm is greater than bnorm, lam is probably set too high
	// reduce it by a factor of 10 and start over
	if ( iter == 2 && rnorm / bnorm > old_rnorm ) {
	    // but don't do it more than 10 times
	    if ( restarts > 20 ) {
		if ( mypid == 0 ) printf("Failure to converge, even with lam = %f\n",lam);
		break;
	    } else {
		++restarts;
		iter = 1;
		lam /= 2.0;
		// reset these 
		old_rnorm = 1.0;
		for ( int j = 0 ; j < nnz ; ++j ) {
		    xvol_sph[j]  = 0.0;
		    pxvol_loc[j] = 0.0; 
		    // bvol_loc[j]  = 0.0; // can reuse this instead, catch restarts != 0 in primary loop
		}
		if ( mypid == 0 ) printf("reducing lam to %f, restarting\n", lam);
		continue;
	    }
	}
	// if changes are sufficiently small, or if no further progress is made, terminate
	if ( rnorm / bnorm < tol || rnorm / bnorm > old_rnorm ) {
	    break;
	}
	// update the termination threshold
	old_rnorm = rnorm / bnorm;
	// update the reconstructed volume
	for ( int j = 0 ; j < nnz ; ++j ) {
	    xvol_sph[j] += lam * grad[j];
 	    pxvol_loc[j] = 0.0; // reset it so it's ready to accumulate for the next iteration
	}

	++iter;
    }

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
	
    return 0; // recons3d_sirt_mpi
}
// 
// int main(int argc, char ** argv)
// {
//     MPI_Comm comm = MPI_COMM_WORLD;
//     int ncpus, mypid, ierr;
//     int nloc; 
// 
//     MPI_Status mpistatus;
//     MPI_Init(&argc,&argv);
//     MPI_Comm_size(comm,&ncpus);
//     MPI_Comm_rank(comm,&mypid);
//     printf("mypid = %d, ncpus = %d\n", mypid, ncpus);
// 
//     EMData * xvol = new EMData();
// 
//     EMData *volume = new EMData();
//     if ( mypid == 0 ) {
// 	volume->read_image("vol001.tfc");
//     }
//     int nx = volume->get_xsize();
//     ierr = MPI_Bcast(&nx, 1, MPI_INT, 0, comm);
//     float dtheta = 5.0;
//     std::vector<float> ref_angles = Util::even_angles(dtheta);
// 
//     int nangles = ref_angles.size() / 3;
//     if ( mypid == 0 ) printf("Using %d images\n", nangles);
//     int radius = 30;
//     int maxit = 50;
// 
//     int * psize;
//     int * nbase;
//     int nangloc;
//     psize = new int[ncpus];
//     nbase = new int[ncpus];
//     nangloc = setpart(comm, nangles, psize, nbase);
// 
//     EMData ** images = (EMData **) malloc(nangloc * sizeof(EMData *));
//     std::vector<float> anglelist(3,0.0);
//     Dict volparams;
//     volparams["angletype"] = "SPIDER";
// 
//     float * imgdata;
//     EMData * send_projection;
//     volparams["radius"] = radius;
//     if ( mypid == 0 ) {
// 	for ( int i = 0 ; i < nangloc ; ++i ) {	
// 	    anglelist[0] = ref_angles[i*3 + 0];
// 	    anglelist[1] = ref_angles[i*3 + 1];
// 	    anglelist[2] = ref_angles[i*3 + 2];
// 	    volparams["anglelist"] = anglelist;
// 	    images[i] = volume->project(PROJECTOR_TYPE, volparams);
// 	    images[i]->set_attr("phi", anglelist[0]);
// 	    images[i]->set_attr("theta", anglelist[1]);
// 	    images[i]->set_attr("psi", anglelist[2]);
// 	    images[i]->set_attr("s2x", 0.0);
// 	    images[i]->set_attr("s2y", 0.0);
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
// 		
// 	}
//     } else { // mypid != 0
// 	for ( int i = 0 ; i < nangloc ; ++i ) {
// 	    images[i] = new EMData();
// 	    images[i]->set_size(nx,nx);
// 	    MPI_Recv(&anglelist[0], 3, MPI_FLOAT, 0, mypid, comm, &mpistatus);
// 	    images[i]->set_attr("phi", anglelist[0]);
// 	    images[i]->set_attr("theta", anglelist[1]);
// 	    images[i]->set_attr("psi", anglelist[2]);
// 	    images[i]->set_attr("s2x", 0.0);
// 	    images[i]->set_attr("s2y", 0.0);
// 	    imgdata = images[i]->get_data();
// 	    MPI_Recv(imgdata, nx*nx, MPI_FLOAT, 0, mypid, comm, &mpistatus);
// 	}
//     }
//     if ( mypid == 0 ) printf("Done with data distribution\n");
// 	
//     float lam = 1.0e-4;
//     float tol = 1.0e-3;
//     std::string symmetry = "c1";
// 
//     recons3d_sirt_mpi(comm, images, xvol, nangloc, radius, lam, maxit, symmetry, tol);
//     if ( mypid == 0 ) printf("Done with SIRT\n");
// 
//     EMUtil::ImageType WRITE_SPI = EMUtil::IMAGE_SINGLE_SPIDER;
//     if ( mypid == 0 ) {
// 	xvol->write_image("volume_sirt.spi", 0, WRITE_SPI);
//     }
// 
// 
//     for ( int i = 0 ; i < nangloc ; ++i ) {	
// 	EMDeletePtr(images[i]);
//     }		
//     free(images);
//     EMDeletePtr(xvol);
//     EMDeletePtr(volume);
// 
//     ierr = MPI_Finalize();
// 
//     return 0; // main
// }
