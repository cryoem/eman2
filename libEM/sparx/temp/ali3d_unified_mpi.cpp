#include "mpi.h"
#include "stdlib.h"

// include EMAN2
#include "emdata.h"
#include "assert.h"
#include "projector.h"

#include "ali3d_unified_mpi.h"

// include CCTBX
#include <scitbx/lbfgsb/raw.h>
#include <scitbx/array_family/shared.h>

using namespace EMAN;
using namespace std;

namespace {

    using scitbx::af::shared;
    using scitbx::fn::pow2;
    using namespace scitbx::lbfgsb::raw;

    template <typename ElementType>
    ref1<ElementType>
    make_ref1(shared<ElementType>& a)
    {
	return ref1<ElementType>(a.begin(), a.size());
    }
} // namespace <anonymous>


#define grad(i) grad[(i)-1]
int  unified(MPI_Comm comm, EMData *volume, EMData **projdata, 
             float *angleshift, int nloc, int max_iter, char *fname_base)
{
    int mypid, ncpus;
    MPI_Status mpistatus;

    float *voldata, *imgdata;
    int	  *ptrs, *cord;
    int	  ndim, nx, ny, nz, ierr, nang, iglb, nangloc;
    float psi, theta, phi, sx, sy;
    int	  nnz, nrays, ri;
    EMData *imgbuf;
    int	   *psize, *nbase;
    char  filename[100];
    Vec3i volsize;
    Vec3i origin;

    float *rhs;

    // get parallel configuration

    MPI_Comm_rank(comm,&mypid);
    MPI_Comm_size(comm,&ncpus);

    // nang is the total number of images
    MPI_Allreduce(&nloc, &nang, 1, MPI_INT, MPI_SUM, comm);

    // we need psize and nbase to place distributed angleshift into the x array for optimization
    psize = new int[ncpus];
    nbase = new int[ncpus];
    nangloc = setpart(comm, nang, psize, nbase);
    if (nangloc != nloc) {
       printf("   ali3d_unified: nloc does not match with nangloc on Proc %d , exit...\n", mypid);
       ierr = MPI_Finalize();
       exit(1);
    }
    for (int i = 0; i < ncpus; i++) {
       psize[i] = psize[i]*5;
       nbase[i] = nbase[i]*5;
    }

    ndim = volume->get_ndim();
    nx   = volume->get_xsize();
    ny   = volume->get_ysize();
    nz   = volume->get_zsize();

    volsize[0] = nx;
    volsize[1] = ny;
    volsize[2] = nz;
    origin[0] = nx/2 + 1;
    origin[1] = ny/2 + 1;
    origin[2] = nz/2 + 1;

    voldata = volume->get_data();

    // the following will have to be cleaned up later...
    rhs = new float[nx*ny*nloc];
	
    for ( int i = 0 ; i < nloc ; ++i ) {
	imgdata = projdata[i]->get_data();
	for ( int j = 0 ; j < nx * ny ; ++j ) {
	    rhs[nx*ny*i + j] = *(imgdata + j);
	}
    }
	
    double aba;
    ri = nx/2 - 1;

    ierr = getnnz(volsize, ri, origin, &nrays, &nnz);
    if (mypid ==0) printf("    nnz = %d\n", nnz);
	
    Dict  myparams;

    // LBFGS-B setup
    std::string task, csave;
    shared<bool> lsave_(4);
    ref1<bool> lsave = make_ref1(lsave_);
    int n=nnz+5*nang;
    int m=5;
    int iprint;
    shared<int> nbd_(n);
    ref1<int> nbd = make_ref1(nbd_);
    shared<int> iwa_(3*n);
    ref1<int> iwa = make_ref1(iwa_);
    shared<int> isave_(44);
    ref1<int> isave = make_ref1(isave_);
    NUMBER f, factr, pgtol;
    shared<NUMBER> x_(n);
    ref1<NUMBER> x = make_ref1(x_);
    shared<NUMBER> l_(n);
    ref1<NUMBER> l = make_ref1(l_);
    shared<NUMBER> u_(n);
    ref1<NUMBER> u = make_ref1(u_);
    shared<NUMBER> g_(n);
    ref1<NUMBER> g = make_ref1(g_);
    shared<NUMBER> dsave_(29);
    ref1<NUMBER> dsave = make_ref1(dsave_);
    shared<NUMBER> wa_(2*m*n+4*n+12*m*m+12*m);
    ref1<NUMBER> wa = make_ref1(wa_);
    NUMBER t1, t2;
    iprint = 0; // don't print out the whole initial guess
    factr=1.0e+1;
    pgtol=1.0e-5;
    for(int i=1;i<=n;i++) {
	nbd(i)=0;
    }

    // pack spherically masked volume into x
    ptrs = new int[nrays+1];
    cord = new int[3*nrays];
    ierr = cb2sph(voldata, volsize, ri, origin, nnz, ptrs, cord, &(x(1))); 

    // pack orientation parameters angleshift into x
    ierr = MPI_Allgatherv(angleshift, 5*nloc, MPI_FLOAT, &(x(nnz + 1)), psize, nbase, MPI_FLOAT, comm);

    for (int i = 0; i<nang; i++) {
	x((nnz+5*i+0)+1) *= (PI/180.0);
	x((nnz+5*i+1)+1) *= (PI/180.0);
	x((nnz+5*i+2)+1) *= (PI/180.0);
	x((nnz+5*i+3)+1) *= -1.0;
	x((nnz+5*i+4)+1) *= -1.0;		
    }

   
    if (mypid == 0) {
	printf(
            "\n"
            "     Structure refinement by the unified approach...\n"
            "\n");
    }

    int numfg = 0;
    char res_fname_base[64];
    task = "START";
 lbl_111:
    setulb(
	n,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint,
	csave,lsave,isave,dsave);
    iprint = 1; // print out info after every iteration 
    if (task.substr(0,2) == "FG") {

	if (numfg >= max_iter) goto Terminate;
	sprintf(res_fname_base, "%s_fg%d" , fname_base, numfg);
	ierr =  fgcalc(comm, &(x(1)), volsize, nnz, nrays, origin, ri, 
		       ptrs, cord, &(x(nnz + 1)), nang, rhs, aba, &f, &(g(1)), res_fname_base); 

	numfg++;
	goto lbl_111;
    }
    if (task.substr(0,5) == "NEW_X") goto lbl_111;
 Terminate:
    sph2cb(&(x(1)), volsize, nrays, ri, nnz, ptrs, cord, voldata);

    EMUtil::ImageType WRITE_SPI = EMUtil::IMAGE_SINGLE_SPIDER;
    char out_fname[64];
    sprintf(out_fname, "vol_r_%s.spi", fname_base);
    if (mypid == 0) volume->write_image(out_fname, 0, WRITE_SPI);

    // copy new angles and shifts to angleshift from x
    for ( int i = 0 ; i < 5 * nangloc ; ++i ) {
	angleshift[i] = x(nnz + 1 + nbase[mypid] + i);
    }
    // and convert back to ali3d_d's format
    for ( int i = 0 ; i < nangloc ; ++i ) {
	angleshift[5*i + 0] *= (180.0/PI);
	angleshift[5*i + 1] *= (180.0/PI);
	angleshift[5*i + 2] *= (180.0/PI);
	angleshift[5*i + 3] *= -1.0;
	angleshift[5*i + 4] *= -1.0;		
    }

    EMDeleteArray(rhs);
    EMDeleteArray(ptrs);
    EMDeleteArray(cord);
    EMDeleteArray(psize);
    EMDeleteArray(nbase);

    return 0;
	
}

int asta2(float *img, int nx, int ny, int ri, double *abaloc, int *klploc)
{
	int xcent = (nx / 2) + 1;
	int ycent = (ny / 2) + 1;
	int r_squared = ri*ri;

	int x_summand, y_summand;
	
	for ( int i = 0 ; i < nx ; ++i ) {
		x_summand = (i-xcent) * (i-xcent);
		for ( int j = 0 ; j < ny ; ++j ) {
			y_summand = (j-ycent) * (j-ycent);
			if ( x_summand + y_summand > r_squared ) {
				*abaloc += (double) img[j*nx + i];
//chao set the background to zero
                                img[j*nx+i]=0.0;
				++*klploc;
			}
		}
	}
	
	return 0;
}

//#define cube(i,j,k) cube[ ((k-1)*ny + j-1)*nx + i-1 ] 
// Phi:: put parenths around params to cube
// Phi:: switch index order of cube to try to get it to match spider/fortran version
#define cube(i,j,k) cube[ (((k)-1)*ny + (j)-1)*nx + (i)-1 ] 
// #define cube(k,j,i) cube[ (((k)-1)*ny + (j)-1)*nx + (i)-1 ] 
#define sphere(i)   sphere[(i)-1]
#define cord(i,j)   cord[((j)-1)*3 + (i) -1]
#define ptrs(i)     ptrs[(i)-1]
#define dm(i)       dm[(i)-1]

int cb2sph(float *cube, Vec3i volsize, int    ri, Vec3i origin, 
           int    nnz0, int     *ptrs, int *cord, float *sphere) 
{
    int    xs, ys, zs, xx, yy, zz, rs, r2;
    int    ix, iy, iz, jnz, nnz, nrays;
    int    ftm = 0, status = 0;  

    int xcent = (int)origin[0];
    int ycent = (int)origin[1];
    int zcent = (int)origin[2];

    int nx = (int)volsize[0];
    int ny = (int)volsize[1];
    int nz = (int)volsize[2];

    r2      = ri*ri;
    nnz     = 0;
    nrays    = 0;
    ptrs(1) = 1;

    for (ix = 1; ix <= nx; ix++) {
       xs  = ix-xcent;
       xx  = xs*xs;
       for ( iy = 1; iy <= ny; iy++ ) {
           ys = iy-ycent;
           yy = ys*ys;
           jnz = 0;

           ftm = 1;
           // not the most efficient implementation
           for (iz = 1; iz <= nz; iz++) {
               zs = iz-zcent;
               zz = zs*zs;
               rs = xx + yy + zz;
               if (rs <= r2) {
                  jnz++;
                  nnz++;
                  sphere(nnz) = cube(iz, iy, ix); 

                  //  record the coordinates of the first nonzero ===
                  if (ftm) {
  		     nrays++;
                     cord(1,nrays) = iz; 
                     cord(2,nrays) = iy; 
                     cord(3,nrays) = ix;
                     ftm = 0;
                  }
               }
            } // end for (iz..)
            if (jnz > 0) {
		ptrs(nrays+1) = ptrs(nrays) + jnz;
	    }  // endif (jnz)
       } // end for iy
    } // end for ix
    if (nnz != nnz0) status = -1;
    return status;
}

// decompress sphere into a cube
int sph2cb(float *sphere, Vec3i volsize, int  nrays, int    ri, 
           int      nnz0, int     *ptrs, int  *cord, float *cube)
{
    int       status=0;
    int       r2, i, j, ix, iy, iz,  nnz;

    int nx = (int)volsize[0];
    int ny = (int)volsize[1];
    // int nz = (int)volsize[2];

    r2      = ri*ri;
    nnz     = 0;
    ptrs(1) = 1;

    // no need to initialize
    // for (i = 0; i<nx*ny*nz; i++) cube[i]=0.0;

    nnz = 0;
    for (j = 1; j <= nrays; j++) {
       iz = cord(1,j);
       iy = cord(2,j);
       ix = cord(3,j);
       for (i = ptrs(j); i<=ptrs(j+1)-1; i++, iz++) {
           nnz++;
	   cube(iz,iy,ix) = sphere(nnz);
       }
    }
    if (nnz != nnz0) status = -1;
    return status;
}

#define x(i)        x[(i)-1]
#define y(i,j)      y[((j)-1)*nx + (i) - 1]

// project from 3D to 2D (single image)
int fwdpj3(Vec3i volsize, int nrays, int   nnz, float *dm, 
           Vec3i  origin, int    ri, int *ptrs, int *cord, 
           float      *x, float  *y)
{
    /*
        purpose:  y <--- proj(x)
        input  :  volsize  the size (nx,ny,nz) of the volume
                  nrays    number of rays within the compact spherical 
                           representation
                  nnz      number of voxels within the sphere
                  dm       an array of size 9 storing transformation 
                           associated with the projection direction
                  origin   coordinates of the center of the volume
                  ri       radius of the sphere
                  ptrs     the beginning address of each ray
                  cord     the coordinates of the first point in each ray
                  x        3d input volume
                  y        2d output image 
    */

    int    iqx, iqy, i, j, xc, yc, zc;
    float  ct, dipx, dipy, dipx1m, dipy1m, xb, yb, dm1, dm4;
    int    status = 0;
    
    // Phi: adding the shift parameters that get passed in as the last two entries of dm
    float sx, sy;

    sx = dm(7);
    sy = dm(8);

    int xcent = origin[0];
    int ycent = origin[1];
    int zcent = origin[2];

    int nx = volsize[0];
    int ny = volsize[1];

    dm1 = dm(1);
    dm4 = dm(4);
 
    if ( nx > 2*ri ) {
	for (i = 1; i <= nrays; i++) {

            zc = cord(1,i)-zcent;
            yc = cord(2,i)-ycent;
            xc = cord(3,i)-xcent;
            xb = zc* dm(1) +yc* dm(2) +xc* dm(3) + xcent + sx;
            yb = zc* dm(4) +yc* dm(5) +xc* dm(6) + ycent + sy;

            for (j = ptrs(i); j< ptrs(i+1); j++) {
               iqx = ifix(xb);
               iqy = ifix(yb);

  	       ct   = x(j);

               // dipx =  xb - (float)(iqx);
               // dipy = (yb - (float)(iqy)) * ct;
	           dipx =  xb - iqx;
	           dipy = (yb - iqy) * ct;

               dipy1m = ct - dipy;
               dipx1m = 1.0 - dipx;

			if (iqx <= nx && iqy <= ny && iqx >= 1 && iqy >= 1) 
               // y(iqx  ,iqy)   = y(iqx  ,iqy)   + dipx1m*dipy1m;
               y(iqx  ,iqy)   +=  dipx1m*dipy1m;
			if (iqx + 1 <= nx && iqy <= ny && iqx >= 0 && iqy >= 1) 
               // y(iqx+1,iqy)   = y(iqx+1,iqy)   + dipx*dipy1m; 
               y(iqx+1,iqy)   +=  dipx*dipy1m; 
			if (iqx + 1 <= nx && iqy + 1 <= ny && iqx >= 0 && iqy >= 0) 
               // y(iqx+1,iqy+1) = y(iqx+1,iqy+1) + dipx*dipy;         
               y(iqx+1,iqy+1) +=  dipx*dipy;         
			if (iqx <= nx && iqy + 1 <= ny && iqx >= 1 && iqy >= 0) 
               // y(iqx  ,iqy+1) = y(iqx  ,iqy+1) + dipx1m*dipy;
               y(iqx  ,iqy+1) +=  dipx1m*dipy;
               xb += dm1;
               yb += dm4;
	   }
	}
    }
    else {
	fprintf(stderr, " nx must be greater than 2*ri\n");
        exit(1);
    }
    return status;
}
#undef x
#undef y

#define y(i)        y[(i)-1]
#define x(i,j)      x[((j)-1)*nx + (i) - 1]

// backproject from 2D to 3D for a single image
int bckpj3(Vec3i volsize, int nrays, int   nnz, float *dm, 
           Vec3i  origin, int    ri, int *ptrs, int *cord, 
           float      *x, float *y)
{
    int       i, j, iqx,iqy, xc, yc, zc;
    float     xb, yb, dx, dy, dx1m, dy1m, dxdy;
    int       status = 0; 

    int xcent = origin[0];
    int ycent = origin[1];
    int zcent = origin[2];

    int nx = volsize[0];
    int ny = volsize[1];

    // Phi: adding the shift parameters that get passed in as the last two entries of dm
    float sx, sy;

    sx = dm(7);
    sy = dm(8);


    if ( nx > 2*ri) {
	for (i = 1; i <= nrays; i++) {
	    zc = cord(1,i) - zcent;
	    yc = cord(2,i) - ycent;
            xc = cord(3,i) - xcent;

            xb = zc*dm(1)+yc*dm(2)+xc*dm(3) + xcent + sx;
            yb = zc*dm(4)+yc*dm(5)+xc*dm(6) + ycent + sy;

            for (j = ptrs(i); j <ptrs(i+1); j++) {
		iqx = ifix(xb);
		iqy = ifix(yb);

		dx = xb - iqx;
		dy = yb - iqy;
		dx1m = 1.0 - dx;
		dy1m = 1.0 - dy;
		dxdy = dx*dy;
/*
c               y(j) = y(j) + dx1m*dy1m*x(iqx  , iqy)
c     &                     + dx1m*dy  *x(iqx  , iqy+1)
c     &                     + dx  *dy1m*x(iqx+1, iqy)
c     &                     + dx  *dy  *x(iqx+1, iqy+1)  
c
c              --- faster version of the above commented out
c                  code (derived by summing the following table 
c                  of coefficients along  the colunms) ---
c
c                        1         dx        dy      dxdy
c                     ------   --------  --------  -------
c                      x(i,j)   -x(i,j)   -x(i,j)    x(i,j)  
c                                        x(i,j+1) -x(i,j+1)
c                              x(i+1,j)           -x(i+1,j)
c                                                x(i+1,j+1) 
c
*/
		// Phi: add index checking, now that shifts are being used
		if ( iqx <= nx && iqy <= ny && iqx >= 1 && iqy >= 1 ) {
		    y(j) += x(iqx,iqy);
		    if ( iqx + 1 <= nx && iqx + 1 >= 1 ) {
			y(j) += dx*(-x(iqx,iqy)+x(iqx+1,iqy));
		    }
		    if ( iqy + 1 <= ny && iqy + 1 >= 1 ) {
			y(j) += dy*(-x(iqx,iqy)+x(iqx,iqy+1));
		    }
		    if ( iqx + 1 <= nx && iqy + 1 <= ny && iqx + 1 >= 1 && iqy + 1 >= 1 ) {
			y(j) += dxdy*( x(iqx,iqy) - x(iqx,iqy+1) -x(iqx+1,iqy) + x(iqx+1,iqy+1) );
		    }
		}

//                y(j) += x(iqx,iqy)
//                     +  dx*(-x(iqx,iqy)+x(iqx+1,iqy))
//                     +  dy*(-x(iqx,iqy)+x(iqx,iqy+1))
//                     +  dxdy*( x(iqx,iqy) - x(iqx,iqy+1) 
//                              -x(iqx+1,iqy) + x(iqx+1,iqy+1) );

               xb += dm(1);
               yb += dm(4);
	    } // end for j
	} // end for i
     }
    else {
	fprintf(stderr, "bckpj3: nx must be greater than 2*ri\n");
    }

    return status;
}

#undef x
#undef y
#undef dm

int ifix(float a) 
{
    int ia = 0;
    if (a >= 0.0) {
       ia = (int)floor(a);
    }
    else {
       ia = (int)ceil(a);
    }
    return ia;
}

int getnnz(Vec3i volsize, int ri, Vec3i origin, int *nrays, int *nnz) 
/*
   purpose: count the number of voxels within a sphere centered
            at origin and with a radius ri.

     input:
     volsize contains the size information (nx,ny,nz) about the volume
     ri      radius of the object embedded in the cube.
     origin  coordinates for the center of the volume
   
     output:
     nnz    total number of voxels within the sphere (of radius ri)
     nrays  number of rays in z-direction. 
*/
{
    int  ix, iy, iz, rs, r2, xs, ys, zs, xx, yy, zz;
    int  ftm=0, status = 0;

    r2    = ri*ri;
    *nnz  = 0;
    *nrays = 0;
    int nx = volsize[0];
    int ny = volsize[1];
    int nz = volsize[2];
    // int nx = (int)volsize[0];
    // int ny = (int)volsize[1];
    // int nz = (int)volsize[2];

    int xcent = origin[0]; 
    int ycent = origin[1]; 
    int zcent = origin[2]; 
    // int xcent = (int)origin[0]; 
    // int ycent = (int)origin[1]; 
    // int zcent = (int)origin[2]; 

    // need to add some error checking 
    for (ix = 1; ix <=nx; ix++) { 
	xs  = ix-xcent;
	xx  = xs*xs;
        for (iy = 1; iy <= ny; iy++) {
            ys = iy-ycent;
            yy = ys*ys; 
            ftm = 1;
            for (iz = 1; iz <= nz; iz++) {
		zs = iz-zcent;
		zz = zs*zs;
		rs = xx + yy + zz;
		if (rs <= r2) {
		    (*nnz)++;
		    if (ftm) {
                       (*nrays)++;
                       ftm = 0;
                    }
                }
            }
	} // end for iy
    } // end for ix
    return status;
}

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

// This might be a convenient function to use in fgcalc, rather than repeating the same code six times
int make_proj_mat(float phi, float theta, float psi, float * dm)
{
//     float cphi=cos(phi);
//     float sphi=sin(phi);
//     float cthe=cos(theta);
//     float sthe=sin(theta);
//     float cpsi=cos(psi);
//     float spsi=sin(psi);

    double cphi, sphi, cthe, sthe, cpsi, spsi;
    double dphi = phi;
    double dthe = theta;
    double dpsi = psi;
    sincos(dphi, &sphi, &cphi);
    sincos(dthe, &sthe, &cthe);
    sincos(dpsi, &spsi, &cpsi);

    dm[0]=cphi*cthe*cpsi-sphi*spsi;
    dm[1]=sphi*cthe*cpsi+cphi*spsi;
    dm[2]=-sthe*cpsi;
    dm[3]=-cphi*cthe*spsi-sphi*cpsi;
    dm[4]=-sphi*cthe*spsi+cphi*cpsi;
    dm[5]=sthe*spsi;

    return 0;
}

