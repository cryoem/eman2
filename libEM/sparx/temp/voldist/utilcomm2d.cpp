#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <string.h>

int setpart_gc1(MPI_Comm comm_2d, int nangs, int *psize, int *nbase)
/* This function provides the 1-D partition for the set of 2D images along the first column of processors

	Input:
	    comm_2d - Cartesian communicator
	      nangs - total number of 2D images
	Output:
	      psize - vector containing the partition distribution
                  pisze[i] gives the number of images assigned to processor (i,0) on
                  the 2-D processor grid

	      nbase - vector containing the base of the partitions
                  nbase[i]+1 gives the global index of the first image assigned 
                  to processor (i,0) on the 2-D processor grid.
 
	   nangsloc - number of local angles (number of angles assigned to the calling (local) processor)
*/
{
    int ROW = 0, COL = 1;
    int dims[2], periods[2], mycoords[2];
    int nangsloc, nrem;

    // Get information associated with comm_2d
    MPI_Cart_get(comm_2d, 2, dims, periods, mycoords);

    int nrows = dims[ROW];

    nangsloc = nangs/nrows;
    nrem = nangs - nrows*nangsloc;
    if (mycoords[ROW] < nrem) nangsloc++;

    for (int i = 0; i < nrows; i++) {
       psize[i] = nangs/nrows;
       if (i < nrem) psize[i]++;
    }

    nbase[0] = 0; 
    for (int i = 1; i < nrows; i++) {
	    nbase[i] = nbase[i-1] + psize[i-1];
    }

    return nangsloc;
}

int getnnz(int *volsize, int ri, int *origin, int *nrays, int *nnz) 
/*
   purpose: count the number of voxels within a spheric mask centered
              at origin and with a radius ri.

   input:
     volsize vector containing size of the volume/cube
     ri      radius of the spherical mask embedded in the cube.
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
/* This function provides the ideal 1-D partition for the 3D volume that contains nnz voxels
   along the first row of processors.

	Input:
	    comm_2d - Cartesian communicator
	        nnz - total number of voxels in the volume
	Output:
	    nnzpart - vector containing the partition distribution
                nnzpart[i] gives the number of voxels assigned to processor (0,i)
     
	    nnzbase - vector containing the base of the partitions
                nnzbase[i]+1 gives the global index of the first voxel assigned to 
                processor (0,i)        
  
	     nnzloc - initial local nnz
*/
{
    int ROW = 0, COL = 1;
    int dims[2], periods[2], mycoords[2];
    int nnzloc, nrem;

    // Get information associated with comm_2d
    MPI_Cart_get(comm_2d, 2, dims, periods, mycoords);

    int ncols = dims[COL];
    nnzloc = nnz/ncols;
    nrem = nnz - ncols*nnzloc;
    if (mycoords[COL] < nrem) nnzloc++;

    for (int i = 0; i < ncols; i++) {
       nnzpart[i] = nnz/ncols;
       if (i < nrem) nnzpart[i]++;
    }

    nnzbase[0] = 0; 
    for (int i = 1; i < ncols; i++) {
	    nnzbase[i] = nnzbase[i-1] + nnzpart[i-1];
    }
    nnzbase[ncols] = nnz;
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
/* This function divides the volume among the column processors by computing nraysloc 
   for each processor and the starting pointer for the first ray.

	Input:
	    comm_2d - Cartesian communicator
	      nrays - total number of rays
	       ptrs - vector containing pointers
                      ptr[i] gives the global index of the first voxel
                      along the i-th ray

	    nnzbase - ideal volume partition of nnz
                      nnzbase[i]+1 gives global index of the first
                      voxel assigned to the i-th column processor group

	Output:
	   ptrstart - vector containing all the starting pointers for each column processor group
               ptrstart[i] gives the global index of the first ray assigned to the i-th column group
	   nraysloc - actual number of rays assigned to the local processor
*/
{
    int ROW = 0, COL = 1;
    int dims[2], periods[2], mycoords[2];
    int nraysloc = 0;

    // Get information associated with comm_2d
    MPI_Cart_get(comm_2d, 2, dims, periods, mycoords);
    int ncols = dims[COL];

    if (mycoords[COL] == 0){ // First column group starts out with the first ray
      nraysloc++;
    }
    ptrstart[0] = 0;

    int cgrpid = 1; // cgrpid is used to represent the column group number
    // go through all other rays 
    for(int i=1; i<nrays ; i++){
        if (ptrs[i]-1 <= nnzbase[cgrpid] && ptrs[i+1]-1 >= nnzbase[cgrpid]){ 
           // the i+1st ray won't fit in the cgrpid-1 column group according
           // to the ideal partition limit. Need to decide whether to
           // assign it to cgrpid-1 or to cgrpid

	   if (nnzbase[cgrpid]-(ptrs[i]-1)>= ptrs[i+1]-1-nnzbase[cgrpid]){
              // assign it to cgrpid-1
              if  (mycoords[COL] == cgrpid-1){  
                 nraysloc++;
              }
              ptrstart[cgrpid] = i+1;
              cgrpid++;
	   } 
           else { 
              //nnzbase[cgrpid]-(ptrs[i]-1)>= ptrs[i+1]-1-nnzbase[cgrpid]
              // assign it to cgrpid
              if(mycoords[COL] == cgrpid){// ray belongs to cgrpid and it's a new starting ray
	         nraysloc++;
              }
              ptrstart[cgrpid] = i;
              cgrpid++;
	   }
       } 
       else {  
           //ptrs[i]-1 > nnzbase[cgrpid], so the i-th belongs to cgrpid-1
           if(mycoords[COL] == cgrpid-1){  // ray belongs to cgrpid
	      nraysloc++;
           }
       }
    } // end for
    ptrstart[ncols] = nrays;
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


#define y(i)        y[(i)-1]
#define x(i,j)      x[((j)-1)*nx + (i) - 1]

int bckpj3_Cart(int *volsize, int nraysloc, int nnzloc, float *dm, 
           int *origin, int ri, int *ptrs, int *cord, int myptrstart,
           float *x, float *y)
/* This function will backproject from 2D to 3D for a single image for one portion of the 3D volume.  This function is to be run on every processor.

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

    int  mystartvol = ptrs[myptrstart];

    if ( nx > 2*ri) {
	for (i = 1; i <= nraysloc; i++) {
	    zc = cord(1,myptrstart+i) - zcent;
	    yc = cord(2,myptrstart+i) - ycent;
            xc = cord(3,myptrstart+i) - xcent;

            xb = zc*dm(1)+yc*dm(2)+xc*dm(3) + xcent + sx;
            yb = zc*dm(4)+yc*dm(5)+xc*dm(6) + ycent + sy;

            for (j = ptrs(myptrstart+i)-mystartvol; j <ptrs(myptrstart+i+1); j++) {
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
