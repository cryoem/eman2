#include "mpi.h"

#include "emdata.h"
#include "project3d.h"

#define cube(i,j,k) cube[ (((k)-1)*ny + (j)-1)*nx + (i)-1 ] 
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
