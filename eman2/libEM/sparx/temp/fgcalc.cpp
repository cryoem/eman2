#include "mpi.h"
#include "emdata.h"
#include "project3d.h"
#include "fgcalc.h"
#include "utilcomm.h"

#define gradloc(i) gradloc[(i)]
#define rhs(i,j)   rhs[((j)-1)*nx*ny + (i) - 1]
#define rvec(i)    rvec[(i) - 1]
#define dt 1.0e-4
 
int fgcalc(MPI_Comm comm, float *volsph, Vec3i volsize, 
           int nnz, int nrays, Vec3i origin, int ri, 
           int *ptrs, int *cord, float *angtrs, int nang, 
           float *rhs, float aba, NUMBER *fval, float *grad, char * fname_base)
{
	int mypid, ncpus, nvars;
	int *psize, *nbase;

	int	   nangloc, nx, ny, nz, jglb, ierr, ibeg;
	float  phi, theta, psi, sx, sy;
	NUMBER cphi, sphi, cthe, sthe, cpsi, spsi;
	float  dm[8];
	float  *rvec, *gvec1, *gvec2, *gvec3, *gvec4, *gvec5, *gradloc;
	NUMBER fvalloc= 0.0;
	
	MPI_Comm_rank(comm,&mypid);
	MPI_Comm_size(comm,&ncpus);

	nvars = nnz + nang*5;
	*fval = 0.0;
	for (int i = 0; i < nvars; i++) grad[i]=0.0;

	nx = volsize[0];
	ny = volsize[1];
	nz = volsize[2];
 
	psize = new int[ncpus];
	nbase = new int[ncpus];
	rvec  = new float[nx*ny];
	gvec1  = new float[nx*ny];
	gvec2  = new float[nx*ny];
	gvec3  = new float[nx*ny];
	gvec4  = new float[nx*ny];
	gvec5  = new float[nx*ny];
	gradloc = new float[nvars]; 

	
	// float ref_params[5]; // Phi:: stores matrix parameters so they won't get +dt/-dt jitters
	                     // Could do this with just one float instead of an array
	float ref_param; // Phi:: stores matrix parameters so they won't get +dt/-dt jitters
	// ref_param shouldn't be needed, but it seemed that there were some strange things going on nonetheless...

	// Phi:: gradloc was not initialized
	for ( int i = 0 ; i < nvars ; ++i ) {
		gradloc[i] = 0.0;
	}

	std::ofstream res_out;
	char out_fname[64];
	sprintf(out_fname, "res_%s.dat", fname_base);
	res_out.open(out_fname);

	nangloc = setpart(comm, nang, psize, nbase);
	float * img_res = new float[nangloc]; // here we'll store the residuals for each data image
	
	dm[6] = 0.0;
	dm[7] = 0.0;
	
	for (int j = 0; j < nangloc; j++) {
	    img_res[j] = 0.0; // initialize; in the inner loop we accumulate the residual at each pixel
	    for (int i = 0; i<nx*ny; i++) { 
		rvec[i] = 0.0;
		gvec1[i] = 0.0; 
		gvec2[i] = 0.0; 
		gvec3[i] = 0.0; 
		gvec4[i] = 0.0; 
		gvec5[i] = 0.0; 
	    }
	    // rvec <-- P(theta)*f
	    jglb  = nbase[mypid] + j;
		 
	    phi	  = angtrs[5*jglb+0];
	    theta = angtrs[5*jglb+1];
	    psi	  = angtrs[5*jglb+2];
	    sx	  = angtrs[5*jglb+3];
	    sy	  = angtrs[5*jglb+4];

	    cphi=cos(phi);
	    sphi=sin(phi);
	    cthe=cos(theta);
	    sthe=sin(theta);
	    cpsi=cos(psi);
	    spsi=sin(psi);
		
	    // sincos(phi, &sphi, &cphi);
	    // sincos(theta, &sthe, &cthe);
	    // sincos(psi, &spsi, &cpsi);

	    dm[0]=cphi*cthe*cpsi-sphi*spsi;
	    dm[1]=sphi*cthe*cpsi+cphi*spsi;
	    dm[2]=-sthe*cpsi;
	    dm[3]=-cphi*cthe*spsi-sphi*cpsi;
	    dm[4]=-sphi*cthe*spsi+cphi*cpsi;
	    dm[5]=sthe*spsi;
	    dm[6]=sx;
	    dm[7]=sy;

	    ierr = fwdpj3(volsize, nrays, nnz, dm, origin, ri, ptrs, cord, volsph, rvec);

	    // gvec1 <--- P(phi+dt)*f
	    ref_param = phi;
	    phi = phi + dt;
	    // sincos(phi, &sphi, &cphi);
	    cphi = cos(phi);
	    sphi = sin(phi);
	    phi = ref_param;

	    dm[0]=cphi*cthe*cpsi-sphi*spsi;
	    dm[1]=sphi*cthe*cpsi+cphi*spsi;
	    dm[2]=-sthe*cpsi;
	    dm[3]=-cphi*cthe*spsi-sphi*cpsi;
	    dm[4]=-sphi*cthe*spsi+cphi*cpsi;
	    dm[5]=sthe*spsi;

	    ierr = fwdpj3(volsize, nrays, nnz, dm, origin, ri, ptrs, cord, 
			  volsph, gvec1);

	    // phi = phi - dt;
	    // phi = ref_params[0];
	    cphi = cos(phi);
	    sphi = sin(phi);
	    // sincos(phi, &sphi, &cphi);

	    // gvec1 <--- (gvec1 - rvec)/dt
	    for (int i = 0; i<nx*ny; i++) gvec1[i] = (gvec1[i] - rvec[i])/dt;

	    // gvec2 <--- P(theta+dt)*f
	    ref_param = theta;
	    theta = theta + dt;
	    cthe  = cos(theta);
	    sthe  = sin(theta);
	    // sincos(theta, &sthe, &cthe);
	    theta = ref_param;

	    dm[0]=cphi*cthe*cpsi-sphi*spsi;
	    dm[1]=sphi*cthe*cpsi+cphi*spsi;
	    dm[2]=-sthe*cpsi;
	    dm[3]=-cphi*cthe*spsi-sphi*cpsi;
	    dm[4]=-sphi*cthe*spsi+cphi*cpsi;
	    dm[5]=sthe*spsi;

	    ierr = fwdpj3(volsize, nrays, nnz, dm, origin, ri, ptrs, cord, 
			  volsph, gvec2);

	    // theta = theta - dt;
	    // theta = ref_params[1];
	    cthe  = cos(theta);
	    sthe  = sin(theta);
	    // sincos(theta, &sthe, &cthe);
  
	    // gvec2 <--- (gvec2 - rvec)/dt
	    for (int i = 0; i<nx*ny; i++) gvec2[i] = (gvec2[i]-rvec[i])/dt;

	    // gvec3 <--- P(psi+dt)*f
	    ref_param = psi;
	    psi	 = psi + dt;
	    cpsi = cos(psi);
	    spsi = sin(psi);  
	    // sincos(psi, &spsi, &cpsi);
	    psi = ref_param;

	    dm[0]=cphi*cthe*cpsi-sphi*spsi;
	    dm[1]=sphi*cthe*cpsi+cphi*spsi;
	    dm[2]=-sthe*cpsi;
	    dm[3]=-cphi*cthe*spsi-sphi*cpsi;
	    dm[4]=-sphi*cthe*spsi+cphi*cpsi;
	    dm[5]=sthe*spsi;

	    ierr = fwdpj3(volsize, nrays, nnz, dm, origin, ri, ptrs, cord, 
			  volsph, gvec3);

	    // psi	 = psi - dt;
	    // psi = ref_params[2];
	    cpsi = cos(psi);
	    spsi = sin(psi);  
	    // sincos(psi, &spsi, &cpsi);

	    for (int i = 0; i < nx*ny; i++) gvec3[i] = (gvec3[i] - rvec[i])/dt;
		
	    // Phi:: Needed to put this in too; still working with shifted psi in dm before
	    dm[0]=cphi*cthe*cpsi-sphi*spsi;
	    dm[1]=sphi*cthe*cpsi+cphi*spsi;
	    dm[2]=-sthe*cpsi;
	    dm[3]=-cphi*cthe*spsi-sphi*cpsi;
	    dm[4]=-sphi*cthe*spsi+cphi*cpsi;
	    dm[5]=sthe*spsi;

	    // gvec4 <--- P(phi,theta,psi)*f(sx+dt,sy)
	    //Phi:: pass shift in the last two entries of dm
	    ref_param = dm[6];
	    dm[6] += dt;		
		
	    ierr = fwdpj3(volsize, nrays, nnz, dm, origin, ri, ptrs, cord, 
			  volsph, gvec4);
					
	    // dm[6] = ref_params[3];
	    dm[6] = ref_param;
		
	    // gvec4 <--- (gvec4 - rvec)/dt
	    for (int i = 0; i < nx*ny; i++) gvec4[i] = (gvec4[i]-rvec[i])/dt;
		
	    ref_param = dm[7];
	    dm[7] += dt;
	    // gvec5 <--- P(phi,theta,psi)*f(sx,sy+dt)
	    // Phi:: changed typo: last arg was gvec4, is now gvec5
	    ierr = fwdpj3(volsize, nrays, nnz, dm, origin, ri, ptrs, cord, 
			  volsph, gvec5);
					
	    // dm[7] = ref_params[4];
	    dm[7] = ref_param;

	    // gvec5 <--- (gvec5 - rvec)/dt
	    for (int i = 0; i < nx*ny; i++) gvec5[i] = (gvec5[i]-rvec[i])/dt;
	    // std::cout << "\n" << rhs(1,j+1) << "\n" << std::endl;		
	    // rvec <--- rvec - rhs
	    for (int i = 1; i <= nx*ny; i++) {
		rvec(i) = rvec(i) - rhs(i,j+1);
	    }
	    // std::cout << "\n" << rvec[0] << "\n" << std::endl;		
	    ierr = bckpj3(volsize, nrays, nnz, dm, origin, ri, ptrs, cord, 
			  rvec, gradloc);

	    // fval <--- norm(rvec)^2/2 
	    // ibeg = nnz + 5*(jglb-1);
	    // Phi:: Indexing error, the -1 should be outside the parenths
	    // also, the parenths can be removed
	    ibeg = nnz + 5*(jglb) - 1;
	    for (int i = 0; i < nx*ny; i++) {

		fvalloc += rvec[i]*rvec[i]; 
		img_res[j] += rvec[i] * rvec[i];
		gradloc(ibeg+1) += rvec[i]*gvec1[i];
		gradloc(ibeg+2) += rvec[i]*gvec2[i];
		gradloc(ibeg+3) += rvec[i]*gvec3[i];
		gradloc(ibeg+4) += rvec[i]*gvec4[i];
		gradloc(ibeg+5) += rvec[i]*gvec5[i];
	    }
	    // 

	}
	
	ierr = MPI_Allreduce(gradloc, grad, nvars, MPI_FLOAT, MPI_SUM, comm);
	ierr = MPI_Allreduce(&fvalloc, fval, 1, MPI_FLOAT, MPI_SUM, comm);

// 	// in turn, send each array of image residuals to master, then write to disk
	MPI_Status mpistatus;
	if (mypid == 0) {
	    for ( int j = 0 ; j < psize[0] ; ++j ) {
		res_out << std::scientific << img_res[j] << std::endl;
	    }
	    for ( int j = 1 ; j < ncpus ; ++j ) {
		ierr = MPI_Recv(img_res, psize[j], MPI_FLOAT, j, j, comm, &mpistatus);
		for ( int i = 0 ; i < psize[j] ; ++i ) {
		    res_out << std::scientific << img_res[i] << std::endl;
		}
	    }
	} else { // mypid != 0 , send my data to master to write out
	    ierr = MPI_Send(img_res, nangloc, MPI_FLOAT, 0, mypid, comm);
	}


	res_out.close();

	EMDeleteArray(psize);
	EMDeleteArray(nbase);
	EMDeleteArray(rvec);
	EMDeleteArray(gvec1);
	EMDeleteArray(gvec2);
	EMDeleteArray(gvec3);
	EMDeleteArray(gvec4);
	EMDeleteArray(gvec5);
	EMDeleteArray(gradloc);
	EMDeleteArray(img_res);
	
	return 0;
}

 
int fcalc(float *volsph, Vec3i volsize, 
           int nnz, int nrays, Vec3i origin, int ri, 
           int *ptrs, int *cord, float *angtrs, int nang, 
           float *rhs, NUMBER *fval)
{
    int	   nx, ny, nz, jglb, ierr;
    float  phi, theta, psi, sx, sy;
    NUMBER cphi, sphi, cthe, sthe, cpsi, spsi;
    float  dm[8];
    float * rvec;
    *fval = 0.0;

    nx = volsize[0];
    ny = volsize[1];
    nz = volsize[2];
 
    rvec  = new float[nx*ny];

    for (int j = 0; j < nang; j++) {
	for (int i = 0; i<nx*ny; i++) { 
	    rvec[i] = 0.0;
	}
	// rvec <-- P(theta)*f

	phi	  = angtrs[5*j+0];
	theta = angtrs[5*j+1];
	psi	  = angtrs[5*j+2];
	sx	  = angtrs[5*j+3];
	sy	  = angtrs[5*j+4];

	cphi=cos(phi);
	sphi=sin(phi);
	cthe=cos(theta);
	sthe=sin(theta);
	cpsi=cos(psi);
	spsi=sin(psi);

	dm[0]=cphi*cthe*cpsi-sphi*spsi;
	dm[1]=sphi*cthe*cpsi+cphi*spsi;
	dm[2]=-sthe*cpsi;
	dm[3]=-cphi*cthe*spsi-sphi*cpsi;
	dm[4]=-sphi*cthe*spsi+cphi*cpsi;
	dm[5]=sthe*spsi;
	dm[6]=sx;
	dm[7]=sy;

	ierr = fwdpj3(volsize, nrays, nnz, dm, origin, ri, ptrs, cord, volsph, rvec);

	// rvec <--- rvec - rhs
	for (int i = 1; i <= nx*ny; i++) {
	    rvec(i) = rvec(i) - rhs(i,j+1);
	}

	// fval <--- norm(rvec)^2/2 
	for (int i = 0; i < nx*ny; i++) {
	    *fval += rvec[i]*rvec[i]; 
	}
    }
	
    EMDeleteArray(rvec);
	
    return 0;
}
