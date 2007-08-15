#include "mpi.h"
#include "emdata.h"

using namespace EMAN;

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <float.h>

#include "HyBR_Cart.h"
#include "hybr.h"
#include "utilcomm_Cart.h"
#include "project3d_Cart.h"
#include "project3d.h"

int recons3d_HyBR_mpi_Cart(MPI_Comm comm_2d, MPI_Comm comm_row   , MPI_Comm comm_col, 
                           EMData ** images, float * angleshift  , EMData *& xvol   , 
                           int nangloc     , int radius          , int maxiter      , 
                           std::string symmetry, int insolve)
 /*
  *   recons3d_HyBR_mpi_Cart(MPI_Comm comm_2d, MPI_Comm comm_row, 
  *     MPI_Comm comm_col, EMData ** images, float * angleshift, 
  *     EMData *& xvol, int nangloc, int radius, int maxiter, 
  *     std::string symmetry, int insolve)
  *
  * HyBR is a Hybrid Bidiagonalization Regularization method used for 
  * solving large-scale, ill-posed inverse problems of the form:
  *             b = A*x + noise
  * The method combines an iterative Lanczos Bidiagonalization (LBD) Method 
  * with a SVD-based regularization method to stabilize the semiconvergence
  * behavior that is characteristic of many ill-posed problems.
  *
  * Note: Here we use a "stripped-down" version where all options are set 
  *   to default.
  *
  * This code was written specifically for the cryo-EM project,to be used 
  *   with MPI parallel implementation on a Cartesian virtual topology 
  *  
  *  Input:  
  *   comm_2d, comm_row, comm_col : MPI communicators
  *                        images : 2D projection images (distributed)
  *                    angleshift : vector of angles (distributed)
  *                       nangloc : local number of 2D images
  *                         radius: used to determine spherical format
  *                       maxiter : maximum number of HyBR iterations
  *                      symmetry : symmetry parameter
  *                       insolve : inner solver (0=TSVD, 1=Tikhonov)
  *
  *  Output:              xvol   : HyBR solution
  */
{
  MPI_Status mpistatus;
  int ROW = 0, COL = 1;
  int my2dpid, ierr;
  int mycoords[2], dims[2], periods[2];
  int srpid, srcoords[2]; // Send/receive info

  // Get dims and my coordinates
  MPI_Cart_get(comm_2d, 2, dims, periods, mycoords);
  MPI_Comm_rank(comm_2d, &my2dpid); //Get my pid in the new 2D topology

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

  // vector of symmetrized angles
  std::vector<float> symangles(3,0.0); 

  // Now distribute the volume (in spherical format) among columns 
  // of processors and use nnz to determine the splitting.
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

  int m = nx*nx*nangles;
  int n = nnz; 

  int iterations_save, warning;
  float w, omega,curmin, g, alpha, beta_loc, beta;

  float * Ukloc = myvector(nx*nx*nangloc);
  float ** Vloc = matrix(nnzloc, maxiter+1);
  float ** BB = matrix(maxiter+2, maxiter+1);

  float ** Ub = matrix(maxiter+2, maxiter+2);
  float ** Vb = matrix(maxiter+1, maxiter+1);
  float * Sb = myvector(maxiter+1);

  float * vec = myvector(maxiter+2);
  float * f = myvector(maxiter+2);
  float * Omegas = myvector(maxiter+1);
  float * GCV = myvector(maxiter+1);

  float * x_save = myvector(nnzloc);
  float * x = myvector(nnzloc);

  float * UbTvector;
  float ** BB2;
  float ** Ub2;
  float * Sb2, *f2, *Utb;
  float ** Vb2;
  float * vec2;
  float * row1Ub2 = myvector(maxiter+2);

  for (int i=0; i< nnzloc; i++){
    x[i] = 0.0;
    x_save[i] = 0.0;
    for(int j=0; j<maxiter+1; j++)
      Vloc[i][j] = 0.0;
  }

  EMData * current_image;
  float phi, theta, psi;
  Transform3D RA;
  Transform3D Tf;
  nsym = Tf.get_nsym(symmetry);
  Transform3D::EulerType EULER_SPIDER = Transform3D::SPIDER;
  Dict angdict;

  float * image_data;
  float dm[8];

  //  beta = norm(b, m);
  beta_loc = 0.0;
  for (int i=0; i<nangloc; i++){
    current_image = images[i];
    image_data = current_image->get_data(); 

    for(int j=0; j<nx*nx; j++)
      beta_loc += image_data[j] * image_data[j];
  }
  ierr = MPI_Allreduce (&beta_loc, &beta, 1, MPI_FLOAT, MPI_SUM, comm_col);
  beta = sqrt(beta);

  for (int i=0; i<nangloc; i++){
    current_image = images[i];
    image_data = current_image->get_data(); 

    for(int j=0; j<nx*nx; j++)
      Ukloc[nx*nx*i+j] = image_data[j] / beta;
  }

  vec[0] = beta;
  warning = 0;
  
  double t0;
  t0 = MPI_Wtime();
  for (int i = 0; i< maxiter+1; i++){
    LBD_Cart(comm_2d, comm_row, comm_col, angleshift, volsize, nraysloc, nnzloc, origin, radius, ptrs, myptrstart, cord,nangloc, nx, nx, Ukloc, BB, Vloc, m, n, i, maxiter, symmetry);

    if (i+1 >= 2){ // Use inner solver and adaptive GCV to solve projected problem
      BB2 = submatrix(BB, maxiter+2, maxiter+1, i+2, i+1);
      Ub2 = submatrix(Ub, maxiter+2, maxiter+2, i+2, i+2);
      Vb2 = submatrix(Vb, maxiter+1, maxiter+1, i+1,i+1);
      Sb2 = subvector(Sb, maxiter+1, i+1);
      f2 = subvector(f, maxiter+2 ,i+1);
      vec2 = subvector(vec, maxiter+2, i+2);

      svd(BB2, i+2, i+1, Ub2, Sb2, Vb2); 

      /*Use the adaptive, modified GCV method*/
      UbTvector = myvector(i+2);
      matvec_multT(Ub2,vec2, i+2, i+2, UbTvector);  
      omega = findomega(UbTvector, Sb2, i+2, i+1, insolve);
      Omegas[i-1] = min((float)1.0, omega);

      w = sum_vector(Omegas, i)/((float)(i+1) -(float)1.0);

      alpha = -(float)1.0;
      //Solve the projected problem
      if(insolve == 0)
        tsvd(Ub2, Sb2, Vb2, vec2, i+2, i+1, w, f2, &alpha);
      else
        tikhonov(Ub2, Sb2, Vb2, vec2, i+2, i+1, w, f2, &alpha);

      /*Compute the GCV value used to find the stopping criteria*/
      get_row(Ub2, i+2, 0, row1Ub2);
      g = gcvstopfun(alpha, row1Ub2, Sb2, beta, i+1, m, n, insolve);

      GCV[i-1] = g;

      /*Determine if GCV wants us to stop*/
      if (i+1 > 2){
        /*-------- GCV curve is flat, we stop ------------------------*/
        if (fabs((GCV[i-1]-GCV[i-2]))/GCV[0] < (pow((float)10.0, - (float)6.0))){
          matvec_mult(Vloc, f2, nnzloc, i+1, x);
          // Return volume in cube format.
          if (mycoords[ROW] == 0 && mycoords[COL] != 0 ){
            srcoords[ROW] = 0;
            srcoords[COL] = 0;
            MPI_Cart_rank(comm_2d, srcoords, &srpid);

            MPI_Send(x, nnzloc, MPI_FLOAT, srpid, my2dpid, comm_2d);
          }

          if (mycoords[ROW] == 0 && mycoords[COL] == 0 ){
            printf("Exit Code 1: Time = %11.3e\n", MPI_Wtime()-t0);
            xvol->set_size(nx, nx, nx);
            xvol->to_zero();
            float * voldata = xvol->get_data();
            float * xvol_sph = new float[nnz];

            for(int i=0; i< dims[COL]; i++){
              if (i==0){
                for(int i=0; i<nnzloc; i++)
                  xvol_sph[i] = x[i];
              }else{
                srcoords[ROW] = 0;
                srcoords[COL] = i;
                MPI_Cart_rank(comm_2d, srcoords, &srpid);

                MPI_Recv(&xvol_sph[ptrs[ptrstart[i]]-1], nnzall[i], MPI_FLOAT, srpid, srpid, comm_2d, &mpistatus);
              }
            }

            // unpack the spherical volume back out into the original EMData object
            ierr = sph2cb(xvol_sph, volsize, nrays, radius, nnz, ptrs, cord, voldata);
            EMDeleteArray(xvol_sph);
          }
            return 0;
        }
        /*--- Have warning: Avoid bumps by using a window of 4 its --*/
        else if (warning && i > iterations_save + 3){
          curmin = GCV[iterations_save];
          for(int j = 1; j< 3; j++){
            if (GCV[iterations_save + j]< curmin){
              curmin = GCV[iterations_save+ j];
            }
          }
          if (GCV[iterations_save] <= curmin){
            /* We should have stopped at iterations_save.*/
            copy_vector(x_save, x, nnzloc);

            // Return volume in cube format.
            if (mycoords[ROW] == 0 && mycoords[COL] != 0 ){
              srcoords[ROW] = 0;
              srcoords[COL] = 0;
              MPI_Cart_rank(comm_2d, srcoords, &srpid);

              MPI_Send(x, nnzloc, MPI_FLOAT, srpid, my2dpid, comm_2d);
            }

            if (mycoords[ROW] == 0 && mycoords[COL] == 0 ){  
              printf("Exit Code 2: Time = %11.3e\n", MPI_Wtime()-t0);
              xvol->set_size(nx, nx, nx);
              xvol->to_zero();
              float * voldata = xvol->get_data();
              float * xvol_sph = new float[nnz];

              for(int i=0; i< dims[COL]; i++){
                if (i==0){
                  for(int i=0; i<nnzloc; i++)
                    xvol_sph[i] = x[i];
                }else{
                  srcoords[ROW] = 0;
                  srcoords[COL] = i;
                  MPI_Cart_rank(comm_2d, srcoords, &srpid);

                  MPI_Recv(&xvol_sph[ptrs[ptrstart[i]]-1], nnzall[i], MPI_FLOAT, srpid, srpid, comm_2d, &mpistatus);
                }
              }
  // unpack the spherical volume back out into the original EMData object
              ierr = sph2cb(xvol_sph, volsize, nrays, radius, nnz, ptrs, cord, voldata);
              EMDeleteArray(xvol_sph);
            }
            return 0;
          } 
          else { 
            /* It was just a bump... keep going*/
            warning = 0;
            iterations_save = maxiter;
          }
        }
        /* ----- No warning: Check GCV function-----------------------*/
        else if ( ! warning ){
          if (GCV[i-2] < GCV[i-1]){
            warning = 1;
            matvec_mult(Vloc, f2,nnzloc, i+1, x_save);
            iterations_save = i;
          }
        }

        free_matrix(BB2);
        free_matrix(Ub2);
        free_matrix(Vb2);
        free_vector(Sb2);
        free_vector(vec2);
      }
    matvec_mult(Vloc, f2, nnzloc, i+1, x);
    }
    if (my2dpid ==0) printf("Iteration %d, alpha = %f, omega is %f\n", i+1, alpha, w);
  } //end HyBR main loop

  // Bring all parts of spherical volume together and turn it into cube format.
  if (mycoords[ROW] == 0 && mycoords[COL] != 0 ){
    srcoords[ROW] = 0;
    srcoords[COL] = 0;
    MPI_Cart_rank(comm_2d, srcoords, &srpid);

    MPI_Send(x, nnzloc, MPI_FLOAT, srpid, my2dpid, comm_2d);
  }

  if (mycoords[ROW] == 0 && mycoords[COL] == 0 ){
    printf("Exit Code 3: Time = %11.3e\n", MPI_Wtime()-t0);
  
    xvol->set_size(nx, nx, nx);
    xvol->to_zero();
    float * voldata = xvol->get_data();
    float * xvol_sph = new float[nnz];

    for(int i=0; i< dims[COL]; i++){
      if (i==0){
        for(int i=0; i<nnzloc; i++)
          xvol_sph[i] = x[i];
      }else{
        srcoords[ROW] = 0;
        srcoords[COL] = i;
        MPI_Cart_rank(comm_2d, srcoords, &srpid);

        MPI_Recv(&xvol_sph[ptrs[ptrstart[i]]-1], nnzall[i], MPI_FLOAT, srpid, srpid, comm_2d, &mpistatus);
      }
    }
    // unpack the spherical volume back out into the original EMData object
    ierr = sph2cb(xvol_sph, volsize, nrays, radius, nnz, ptrs, cord, voldata);
    EMDeleteArray(xvol_sph);
  }
 
  EMDeleteArray(ptrs);
  EMDeleteArray(cord);
  EMDeleteArray(psize);
  EMDeleteArray(nbase);
  EMDeleteArray(nnzpart);
  EMDeleteArray(nnzbase);
  EMDeleteArray(ptrstart);

  return 0; // recons3d_HyBR_mpi_Cart
}


/*------LBD_Cart------------------------------------------------------*/

void LBD_Cart(MPI_Comm comm_2d, MPI_Comm comm_row, MPI_Comm comm_col, float *angleshift, Vec3i volsize, int nraysloc, int nnzloc, Vec3i origin, int radius, int *ptrs, int myptrstart, int *cord, int nangloc, int nx, int ny, float *Ukloc, float **B, float **Vloc, int  m, int n, int k, int maxiter,std::string symmetry)
 /*
  *   void LBD_Cart(MPI_Comm comm_2d, MPI_Comm comm_row, MPI_Comm comm_col, 
  *     float *angleshift, Vec3i volsize, int nraysloc, int nnzloc,
  *     Vec3i origin, int radius, int *ptrs, int myptrstart, int *cord, 
  *     int nangloc, int nx, int ny, float *Ukloc, float **B, float **Vloc,
  *     int  m, int n, int k, int maxiter,std::string symmetry)
  *
  *  Perform one step of Lanczos bidiagonalization without
  *  reorthogonalization, no preconditioner here.
  */
{
  float * vc = myvector(nnzloc);
  float * vc_loc = myvector(nnzloc);
  float * uc = myvector(nangloc*nx*ny);
  float * uc_loc = myvector(nangloc*nx*ny);
  float alpha, alpha_loc, beta, beta_loc;
  float * Vk = myvector(nnzloc);
  int ierr;
  float phi, theta, psi;
  float dm[8];

  MPI_Status mpistatus;
  int ROW = 0, COL = 1;
  int  my2dpid;
  int mycoords[2], dims[2], periods[2];
  int srpid, srcoords[2]; // Send/receive info

  // Get dims and my coordinates
  MPI_Cart_get(comm_2d, 2, dims, periods, mycoords);
  MPI_Comm_rank(comm_2d, &my2dpid);

  int nsym = 0;
  Transform3D RA;
  Transform3D Tf;
  nsym = Tf.get_nsym(symmetry);
  Transform3D::EulerType EULER_SPIDER = Transform3D::SPIDER;
  Dict angdict;

  get_col(Vloc, nnzloc, k-1, Vk);
  if (k == 0){ // BACKPROJECTION CODE HERE!!!
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
  
        ierr = bckpj3_Cart(volsize, nraysloc, nnzloc, dm, origin, radius, ptrs, cord, myptrstart, &Ukloc[nx*nx*i], vc_loc);
      }
    }
    ierr = MPI_Allreduce(vc_loc, vc, nnzloc, MPI_FLOAT, MPI_SUM, comm_col);

  } 
  else {  //BACKPROJECTION CODE HERE!!!!
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

        ierr = bckpj3_Cart(volsize, nraysloc, nnzloc, dm, origin, radius, ptrs, cord, myptrstart, &Ukloc[nx*nx*i], vc_loc);
      }
    }
    ierr = MPI_Allreduce(vc_loc, vc, nnzloc, MPI_FLOAT, MPI_SUM, comm_col);

    for(int i=0; i<nnzloc; i++)
      vc[i] = vc[i] - B[k][k-1]*Vk[i];
  }

  //alpha = norm(vc, n);
  alpha_loc = 0.0;
  for(int i=0; i<nnzloc; i++)
    alpha_loc += vc[i]* vc[i];

  MPI_Allreduce (&alpha_loc, &alpha, 1, MPI_FLOAT, MPI_SUM, comm_row);
  alpha = sqrt(alpha);

  for (int i=0; i<nnzloc; i++)
    vc[i] = vc[i] / alpha;

  // FORWARDPROJECTION CODE HERE!!!!!
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

      ierr = fwdpj3_Cart(volsize, nraysloc, nnzloc, dm, origin, radius, ptrs, cord, myptrstart, vc, &uc_loc[nx*nx*i]);
    }
  }
  // Now an all reduce along the rows
  ierr = MPI_Allreduce(uc_loc, uc, nangloc*nx*nx, MPI_FLOAT, MPI_SUM, comm_row);

  for (int i=0; i<nangloc*nx*nx; i++)
    uc[i] = uc[i] - alpha*Ukloc[i];

  //beta = norm(uc, m);
  beta_loc = 0.0;
  for(int i=0; i<nangloc*nx*nx; i++)
    beta_loc += uc[i]* uc[i];

  MPI_Allreduce (&beta_loc, &beta, 1, MPI_FLOAT, MPI_SUM, comm_col);
  beta = sqrt(beta);

  for (int i=0; i<nangloc*nx*nx; i++)
    uc[i] = uc[i] / beta;

  copy_vector(uc, Ukloc,nangloc*nx*nx);
  put_col(Vloc, nnzloc, vc, k);

  B[k][k] = alpha;
  B[k+1][k] = beta;

  free_vector(vc);
  free_vector(vc_loc);
  free_vector(uc);
  free_vector(uc_loc);   
  free_vector(Vk);
} //end LBD_Cart

/*-----findomega------------------------------------------------*/

float findomega(float * bhat, float *s, int m, int n, int insolve)
 /*
  *   omega = findomega(bhat, s, m, n, insolve)
  *
  *  This function computes a value for the omega parameter.
  *
  *  The method: Assume the 'optimal' regularization parameter to be the
  *  smallest singular value.  Then we take the derivative of the GCV
  *  function with respect to alpha, evaluate it at alpha_opt, set the 
  *  derivative equal to zero and then solve for omega.
  *  
  *  Input:   bhat -  vector U'*b, where U = left singular vectors
  *              s -  vector containing the singular values
  *            m,n - dimensions of the (projected) problem
  *        insolve - inner solver
  *  Output:     omega - computed value for the omega parameter.
  */
{
  int i;
  float alpha, t0, t1, t3, t4, t5, v2, omega;
  float * s2 = myvector(n);
  float * tt = myvector(n);
  float * t2 = myvector(n);
  float * v1 = myvector(n);
  float * hold = myvector(m-n);
  float * hold2 = myvector(n);
  float * hold3 = myvector(n);
  float * hold4 = myvector(n);
  float * hold5 = myvector(n);
  float * hold6 = myvector(n);

  if(insolve ==0){
    int kopt = n;
    omega = (m*bhat[kopt-1]*bhat[kopt-1])/(kopt*bhat[kopt-1]*bhat[kopt-1] +2*bhat[kopt]*bhat[kopt]);
  } else {
    /*   First assume the 'optimal' regularization parameter to be the smallest
     *   singular value.
     */
    alpha = s[n-1];

    /*
     * Compute the needed elements for the function.
     */
    for (i = 0; i < m-n; i++)
      hold[i] = fabs(bhat[n+i])*fabs(bhat[n+i]);

    t0 = sum_vector(hold, m-n);

    for (i=0; i<n ; i++){
      s2[i] = fabs(s[i]) * fabs(s[i]);
      tt[i] = (float)1.0 / (s2[i] + alpha*alpha);
      hold2[i] = s2[i]* tt[i];
    }
    
    t1 = sum_vector(hold2,n);
    for (i=0; i<n; i++){
      t2[i] = fabs(bhat[i]*alpha*s[i]) * fabs(bhat[i]*alpha*s[i]);
      hold3[i] = t2[i] * fabs(tt[i]*tt[i]*tt[i]);
      hold4[i] = (s[i]*tt[i])*(s[i]*tt[i]);
      hold5[i] = fabs(alpha*alpha*bhat[i]*tt[i])*(fabs(alpha*alpha*bhat[i]*tt[i]));
      v1[i] = fabs(bhat[i]*s[i])*fabs(bhat[i]*s[i]);
      hold6[i]= v1[i] * fabs(tt[i]*tt[i]*tt[i]);
    }
    t3 = sum_vector(hold3,n);
    t4 = sum_vector(hold4,n);
    t5 = sum_vector(hold5,n);
    
    v2 = sum_vector(hold6,n);
    /*
     * Now compute the omega value.
     */
    omega = ((float)m*alpha*alpha*v2) / (t1*t3 + t4*(t5 + t0));
  }
  free_vector(s2);
  free_vector(tt);
  free_vector(t2);
  free_vector(v1);
  free_vector(hold);
  free_vector(hold2);
  free_vector(hold3);
  free_vector(hold4);
  free_vector(hold5);
  free_vector(hold6); 

  return omega;
}

/*------gcvstopfun------------------------------------------------------------*/

float gcvstopfun(float alpha, float * u, float *s, float beta, int k, int m, int n, int insolve)
 /*
  *  Evaluate the GCV function G(k, alpha) to determine a stopping index.
  *  
  *  Input:   alpha - regularization parameter at current iteration (k)
  *               u - vector Ub'*(vector), where Ub = left singular vectors of BB
  *               s - vector containing singular values of BB
  *            beta - norm of b(right hand side) 
  *               k - current iteration
  *               n - dimension of the original problem
  *         insolve - inner solver
  *  Output:      G - computed value of the GCV function
  */
{
  float beta2, alpha2, num, den, G;
  float *s2 = myvector(k);
  float *t1 = myvector(k);
  float *t2 = myvector(k+1);
  float *t3 = myvector(k);
  int i;

  beta2 = beta*beta;
  alpha2 = alpha*alpha;

  if(insolve == 0){
    for(i=(int) alpha; i<k+1; i++){
      t2[i] = fabs(u[i])*fabs(u[i]);
    }
    G = (n*beta2*sum_vector(t2,k+1))/((m-alpha)*(m-alpha));

  }else{
    for(i=0; i<k; i++)
      s2[i] = fabs(s[i]) * fabs(s[i]) ;

    for (i=0; i<k; i++){
      t1[i] = (float)1.0 / (s2[i] + alpha2);
      t2[i] = fabs(alpha2*u[i] * t1[i]) * fabs(alpha2*u[i] * t1[i]);
      t3[i] = s2[i] * t1[i];
    }

    num = n*beta2*(sum_vector(t2,k) + fabs(u[k])*fabs(u[k]));
    den = (((float)m - sum_vector(t3,k)))*(((float)m - sum_vector(t3, k)));
    
    G = num / den;
  }
  free_vector(s2);
  free_vector(t1);
  free_vector(t2);
  free_vector(t3);

  return G;
}

/*-------CGLS_mpi_Cart----------------------------------------------*/

void recons3d_CGLS_mpi_Cart(MPI_Comm comm_2d, MPI_Comm comm_row , MPI_Comm comm_col, 
                            EMData ** images, float * angleshift, EMData *& xvol   , 
                            int nangloc     , int radius        , int maxiter      ,
                            std::string symmetry)
 /*
  *   recons3d_CGLS_mpi_Cart(MPI_Comm comm_2d, MPI_Comm comm_row, MPI_Comm comm_col, EMData ** images, float * angleshift, EMData *& xvol, int nangloc, int radius, int maxiter, std::string symmetry)
  *
  * CGLS: Conjugate Gradient for Least Squares
  *  
  * Reference: A. Bjorck, "Numerical Methods for Least Squares Problems"       
  * SIAM, 1996, pg. 289.
  * A method for solving large-scale, ill-posed inverse problems of the form:
  *             b = A*x + noise
  *  
  * This code was written specifically for the cryo-EM project,to be used 
  *   with MPI parallel implementation on a Cartesian virtual topology 
  *  
  *  Input:  
  *   comm_2d, comm_row, comm_col : MPI communicators
  *                        images : 2D projection images (distributed)
  *                    angleshift : vector of angles (distributed)
  *                       nangloc : local number of 2D images
  *                         radius: used to determine spherical format
  *                       maxiter : maximum number of HyBR iterations
  *                      symmetry : symmetry parameter
  *
  *  Output:              xvol   : CGLS solution
  */
{
  MPI_Status mpistatus;
  int ROW = 0, COL = 1;
  int my2dpid, ierr;
  int mycoords[2], dims[2], periods[2];
  int srpid, srcoords[2]; // Send/receive info

  // Get dims and my coordinates
  MPI_Cart_get(comm_2d, 2, dims, periods, mycoords);
  MPI_Comm_rank(comm_2d, &my2dpid); //Get my pid in the new 2D topology

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

  // vector of symmetrized angles
  std::vector<float> symangles(3,0.0); 

  // Now distribute the volume (in spherical format) among columns 
  // of processors and use nnz to determine the splitting.
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

  int m = nx*nx*nangles;
  int n = nnz; 

  int k;
  float * trAb_loc = myvector(nnzloc);
  float * trAb = myvector(nnzloc);
  float * q_loc = myvector(nx*nx*nangloc); 
  float * q = myvector(nx*nx*nangloc); 
  float * s = myvector(nx*nx*nangloc);
  float * x = myvector(nnzloc); 
  float * r_loc = myvector(nnzloc);
  float * r = myvector(nnzloc);
  float * p = myvector(nnzloc);
  float nrm_trAb_loc, nrm_trAb, tol, gamma_loc, gamma, beta, oldgamma, nq_loc, nq, alpha;
  float phi, theta, psi;
  float dm[8];
  float * Rnrm = myvector(maxiter+1);

  for (int i=0; i< nnzloc; i++){
    trAb_loc [i] = 0.0;
    trAb [i] = 0.0;
    x[i] = 0.0;
    r_loc[i] = 0.0;
    r[i] = 0.0;
    p[i] = 0.0;
  }
  for (int i=0; i< nx*nx*nangloc; i++){
    q_loc[i] = 0.0;
    q[i] = 0.0;
    s[i] = 0.0;
  }

  EMData * current_image;
  Transform3D RA;
  Transform3D Tf;
  nsym = Tf.get_nsym(symmetry);
  Transform3D::EulerType EULER_SPIDER = Transform3D::SPIDER;
  Dict angdict;

  float * image_data;

  // trAb = A'*b;
  for (int i=0; i<nangloc; i++){
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
      ierr = bckpj3_Cart(volsize, nraysloc, nnzloc, dm, origin, radius, ptrs, cord, 
                         myptrstart, image_data, trAb_loc);
    }
  }
  // Now an all reduce along the columns
  ierr = MPI_Allreduce (trAb_loc, trAb, nnzloc, MPI_FLOAT, MPI_SUM, comm_col);

  nrm_trAb_loc = 0.0;
  for(int i=0; i<nnzloc; i++)
    nrm_trAb_loc += trAb[i]*trAb[i];

  ierr = MPI_Allreduce (&nrm_trAb_loc, &nrm_trAb, 1, MPI_FLOAT, MPI_SUM, comm_row);

  nrm_trAb = sqrt(nrm_trAb);
  tol = sqrt(FLT_EPSILON)*nrm_trAb;

  // s = A*x
  // The following code is needed if we have an initial guess for x.  
/*for (int i=0; i<nangs; i++){
    phi = angles[5*i+0] * PI/180.0;
    theta = angles[5*i+1] * PI/180.0;
    psi = angles[5*i+2] * PI/180.0;
    dm[6] = angles[5*i+3] * -1.0;
    dm[7] = angles[5*i+4] * -1.0;
    
    make_proj_mat(phi, theta, psi, dm);
		  
    ierr = fwdpj3(volsize, nrays, nnz, dm, origin, radius, ptrs, cord, x, &s[nx*ny*i]);
  }
*/
//s = b - s;
  for (int i=0; i<nangloc; i++){
    current_image = images[i];
    image_data = current_image->get_data(); 
    for (int j=0; j<nx*nx; j++)
      s[i*nx*nx+j] = image_data[j]- s[i*nx*nx+j];
  }

//r = A'*s;
  for (int i=0; i<nangloc; i++){
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
      ierr = bckpj3_Cart(volsize, nraysloc, nnzloc, dm, origin, radius, ptrs, cord, 
                         myptrstart, &s[nx*nx*i], r_loc);
    }
  }
  // Now an all reduce along the columns
  ierr = MPI_Allreduce (r_loc, r, nnzloc, MPI_FLOAT, MPI_SUM, comm_col);

  gamma_loc = 0.0;
  for(int i=0; i<nnzloc; i++)
    gamma_loc += r[i]*r[i];
  
  ierr = MPI_Allreduce (&gamma_loc, &gamma, 1, MPI_FLOAT, MPI_SUM, comm_row);
  Rnrm[0] = sqrt(gamma);
  if (my2dpid==0) printf("Iteration %d, Rnorm = %f\n", 1, Rnrm[0]/nrm_trAb);

  k = 0;
  /*
   * Main CGLS loop.
   */
  double t0;
  t0 = MPI_Wtime();
  while (Rnrm[k]>tol && k<= maxiter-1){
    k++;
    if(k==1){
      for(int i=0;i<nnzloc; i++)
        p[i] = r[i];
    }else{
      beta = gamma / oldgamma;
      for(int i=0;i<nnzloc; i++)
        p[i] = r[i]+beta*p[i];
    }
    // q = A*p;
    for (int i=0; i < nx*nx*nangloc; i++)
      q_loc[i] = 0.0;

    for (int i=0; i < nangloc; i++){
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
	
	ierr = fwdpj3_Cart(volsize, nraysloc, nnzloc, dm, origin, radius, ptrs, cord, 
                           myptrstart, p, &q_loc[nx*nx*i]);
      }
    }
    // Now an all reduce along the rows
    ierr = MPI_Allreduce(q_loc, q, nangloc*nx*nx, MPI_FLOAT, MPI_SUM, comm_row);

    nq_loc = 0.0;
    for(int i=0; i< nangloc*nx*nx; i++)
      nq_loc += q[i]*q[i];

    MPI_Allreduce (&nq_loc, &nq, 1, MPI_FLOAT, MPI_SUM, comm_col);

    alpha = gamma/nq;

    for(int i=0; i<nnzloc; i++){
      x[i] = x[i]+alpha*p[i];
    }

    for(int i=0; i<nangloc*nx*nx; i++)
      s[i] = s[i]- alpha*q[i];

    //r = A'*s;
    for(int i=0; i< nnzloc; i++)
      r_loc[i] = 0.0;

    for (int i=0; i<nangloc; i++){
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

        ierr = bckpj3_Cart(volsize, nraysloc, nnzloc, dm, origin, radius, ptrs, cord, 
                           myptrstart, &s[nx*nx*i], r_loc);
      }
    }
    ierr = MPI_Allreduce(r_loc, r, nnzloc, MPI_FLOAT, MPI_SUM, comm_col);

    oldgamma = gamma;
    gamma_loc = 0.0;
    for(int i=0; i< nnzloc; i++)
      gamma_loc += r[i]*r[i];

    MPI_Allreduce (&gamma_loc, &gamma, 1, MPI_FLOAT, MPI_SUM, comm_row);

    Rnrm[k] = sqrt(gamma);
    if(k!=1)
      if (my2dpid ==0) printf("Iteration %d, Rnorm = %f\n", k, Rnrm[k]/nrm_trAb);
  }
  if (my2dpid == 0)printf("Exit CGLS: time = %11.3e\n", MPI_Wtime()-t0);
  // Bring all parts of spherical volume together and turn it into cube format.
  if (mycoords[ROW] == 0 && mycoords[COL] != 0 ){
     srcoords[ROW] = 0;
     srcoords[COL] = 0;
     MPI_Cart_rank(comm_2d, srcoords, &srpid);
     MPI_Send(x, nnzloc, MPI_FLOAT, srpid, my2dpid, comm_2d);
  }

  if (mycoords[ROW] == 0 && mycoords[COL] == 0 ){
     xvol->set_size(nx, nx, nx);
     xvol->to_zero();
     float * voldata = xvol->get_data();
     float * xvol_sph = new float[nnz];

     for(int i=0; i< dims[COL]; i++){
       if (i==0){
          for(int i=0; i<nnzloc; i++) xvol_sph[i] = x[i];
       }
       else {
          srcoords[ROW] = 0;
          srcoords[COL] = i;
          MPI_Cart_rank(comm_2d, srcoords, &srpid);
          MPI_Recv(&xvol_sph[ptrs[ptrstart[i]]-1], nnzall[i], MPI_FLOAT, srpid, srpid, 
                   comm_2d, &mpistatus);
       }
     }
     // unpack the spherical volume back out into the original EMData object
     ierr = sph2cb(xvol_sph, volsize, nrays, radius, nnz, ptrs, cord, voldata);
     EMDeleteArray(xvol_sph);
  }

  EMDeleteArray(ptrs);
  EMDeleteArray(cord);
  EMDeleteArray(psize);
  EMDeleteArray(nbase);
  EMDeleteArray(nnzpart);
  EMDeleteArray(nnzbase);
  EMDeleteArray(ptrstart);

  free_vector(x);
  free_vector(trAb_loc);
  free_vector(trAb);
  free_vector(s);
  free_vector(q_loc);
  free_vector(q);
  free_vector(r_loc);
  free_vector(r);
  free_vector(p);
  free_vector(Rnrm);
}
