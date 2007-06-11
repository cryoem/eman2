#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <string.h>

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
 
	nbase[0] = 0; 
	for (int i = 1; i < dims[ROW]; i++) {
		nbase[i] = nbase[i-1] + psize[i-1];
	}
   
	return nangsloc;
}

int getnnz(int *volsize, int ri, int *origin, int *nrays, int *nnz) 
/*
   purpose: count the number of voxels within a sphere centered
            at origin and with a radius ri.

   input:
     volsize vector containing size of the volume/cube
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

#define cord(i,j)   cord[((j)-1)*3 + (i) -1]
#define ptrs(i)     ptrs[(i)-1]
#define dm(i)       dm[(i)-1]

int getcb2sph(int *volsize, int ri, int *origin, int nnz0, int *ptrs, int *cord) 
/* This function is similar to 'cb2sph' except we don't actually deal with the cube and sphere data. We only use this function to get vectors ptrs and coord.

	Input:
	 volsize - dimensions of 3D cube
	      ri - radius of sphere
	  origin - origin of sphere
   	    nnz0 - number of non-zeros

	Output:
	      ptrs - vector containing pointers to the beginning of each ray
	      cord - vector containing the x,y,z coordinates of the beginning of each ray
*/

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
     //             sphere(nnz) = cube(iz, iy, ix); 

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

int sphpart(MPI_Comm comm_2d, int nrays, int *ptrs, int *nnzbase, int *ptrstart)
/* This function divides the volume among the column processors by computing nraysloc for each processor and the starting pointer for the first ray.

	Input:
	    comm_2d - Cartesian communicator
	      nrays - total number of rays
	       ptrs - vector containing pointers
	    nnzbase - ideal volume partition of nnz
	Output:
	   ptrstart - vector containing all the starting pointers for each column processor group
	   nraysloc - actual number of local rays
*/
{
	int ROW = 0, COL = 1;
	int dims[2], periods[2], mycoords[2];
	int nraysloc = 0;
	
	// Get information associated with comm_2d
	MPI_Cart_get(comm_2d, 2, dims, periods, mycoords);

	int count = 1;
	if (mycoords[COL] == 0){ // First column starts out with the first ray
	  nraysloc++;
	}
	ptrstart[0] = 0;

	for(int i=1; i<nrays ; i++){
	  if (ptrs[i]-1 <= nnzbase[count] && ptrs[i+1]-1 >= nnzbase[count]){ 
		//nnzbase is between or equal to ptrs

	    if (nnzbase[count]-(ptrs[i]-1)>= ptrs[i+1]-1-nnzbase[count]){
	      if(mycoords[COL] == count-1){  // ray belongs to count-1
		nraysloc++;
	      }
	      ptrstart[count] = i+1;
	      count++;

	    } else { //nnzbase[count]-(ptrs[i]-1)>= ptrs[i+1]-1-nnzbase[count]
	        if(mycoords[COL] == count){// ray belongs to count and it's a new starting ray
		  nraysloc++;
		}
		ptrstart[count] = i;
	        count++;
	    }

	}
	else {  //ptrs[i]-1 > nnzbase[count] so ray belongs to count-1
	  if(mycoords[COL] == count-1){
	    nraysloc++;
	  }

	}
	} // end for
	ptrstart[dims[COL]] = nrays;
	return nraysloc;

}

int make_proj_mat(float phi, float theta, float psi, float * dm)
/* This function computes the projection matrix provided the Euler angles.
	Input:
	  phi, theta, psi - Euler angles
	Output:
	  dm - projection matrix
*/
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

int ifix(float a)
/* This function will floor non-negative numbers and ceil negative numbers.
*/
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

#define x(i)        x[(i)-1]
#define y(i,j)      y[((j)-1)*nx + (i) - 1]

int fwdpj3_Cart(int *volsize, int nraysloc, int nnzloc, float *dm, int *origin, int ri, int *ptrs, int *cord, int myptrstart, float *x, float *y)
/* This function will project from 3D to 2D for a single image for one portion of the 3D volume.  This function is to be run on every processor. y = P*x

	Input:
	  volsize - size (nx,ny,nz) of 3D cube volume
	  nraysloc - local number of rays within the compact spherical representation
	  nnzloc - local number of voxels
	  dm - an array of size 9 storing transformation 
	  origin, ri - origin and radius of sphere
	  ptrs, cord - pointers and coordinates for first ray
	  myptrstart - determines starting rays
	  x - portion of the 3D input volume
	Output:
	  y - 2D output image
*/
{
    int    iqx, iqy, i, j, xc, yc, zc, jj;
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
	for (i = 1; i <= nraysloc; i++) {

            zc = cord(1,myptrstart+i)-zcent;
            yc = cord(2,myptrstart+i)-ycent;
            xc = cord(3,myptrstart+i)-xcent;
            xb = zc* dm(1) +yc* dm(2) +xc* dm(3) + xcent + sx;
            yb = zc* dm(4) +yc* dm(5) +xc* dm(6) + ycent + sy;

            for (j = ptrs(myptrstart+i); j< ptrs(myptrstart+i+1); j++) {
	       jj = j-ptrs[myptrstart]+1;
               iqx = ifix(xb);
               iqy = ifix(yb);

		// jj is the index on the local volume
  	       ct   = x(jj);

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

int bckpj3_Cart(int *volsize, int nraysloc, int nnzloc, float *dm, int *origin, int ri, int *ptrs, int *cord, int myptrstart, float *x, float *y)
/* This function will backproject from 2D to 3D for a single image for one portion of the 3D volume.  This function is to be run on every processor. y = P'x

	Input:
	  volsize - size of 3D cube volume
	  nraysloc - local number of rays
	  nnzloc - local number of voxels
	  dm - projection matrix
	  origin, ri - origin and radius of sphere
	  ptrs, cord - pointers and coordinates for first ray
	  myptrstart - determines starting rays
	  x - 2D image
	Output:
	  y - portion of the 3D volume
*/
{
    int       i, j, iqx,iqy, xc, yc, zc, jj;
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
	for (i = 1; i <= nraysloc; i++) {
	    zc = cord(1,myptrstart+i) - zcent;
	    yc = cord(2,myptrstart+i) - ycent;
            xc = cord(3,myptrstart+i) - xcent;

            xb = zc*dm(1)+yc*dm(2)+xc*dm(3) + xcent + sx;
            yb = zc*dm(4)+yc*dm(5)+xc*dm(6) + ycent + sy;

            for (j = ptrs(myptrstart+i); j <ptrs(myptrstart+i+1); j++) {
		jj = j-ptrs[myptrstart]+1;

		// jj is the index on the local volume
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
		    y(jj) += x(iqx,iqy);
		    if ( iqx + 1 <= nx && iqx + 1 >= 1 ) {
			y(jj) += dx*(-x(iqx,iqy)+x(iqx+1,iqy));
		    }
		    if ( iqy + 1 <= ny && iqy + 1 >= 1 ) {
			y(jj) += dy*(-x(iqx,iqy)+x(iqx,iqy+1));
		    }
		    if ( iqx + 1 <= nx && iqy + 1 <= ny && iqx + 1 >= 1 && iqy + 1 >= 1 ) {
			y(jj) += dxdy*( x(iqx,iqy) - x(iqx,iqy+1) -x(iqx+1,iqy) + x(iqx+1,iqy+1) );
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




#define x(i)        x[(i)-1]
#define y(i,j)      y[((j)-1)*nx + (i) - 1]

// project from 3D to 2D (single image)
int fwdpj3(int *volsize, int nrays, int   nnz, float *dm, int *origin, int ri, int *ptrs, int *cord, float *x, float  *y)
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
int bckpj3(int *volsize, int nrays, int   nnz, float *dm, int *origin, int ri, int *ptrs, int *cord, float *x, float *y)
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

