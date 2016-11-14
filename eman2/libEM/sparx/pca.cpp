/*
 * Author: Chao Yang
 * Copyright (c) 2000-2006
 *
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#include "emdata.h"
#include "util.h"
#include "emutil.h"

#include "pca.h"
#include "lapackblas.h"

using namespace EMAN;

// return all right singular vectors
int PCA::dopca(vector <EMData*> imgstack, EMData *mask)
{
   // performs PCA on a list of images (each under a mask)
   // returns a list of eigenimages

   int i, status = 0;

   vector<EMData*> img1dlst;
   vector<EMData*> eigvecs;

   int nimgs = imgstack.size();

   for (i=0; i<nimgs; i++) {
      img1dlst.push_back(Util::compress_image_mask(imgstack[i],mask));
   }

   // for right now, compute a full SVD
   eigvecs = Util::svdcmp(img1dlst, 0);

   img1dlst.clear();

   for (i=0; i<nimgs; i++) {
      eigenimages.push_back(Util::reconstitute_image_mask(eigvecs[i],mask));
   }

   eigvecs.clear();

   return status;
}

//------------------------------------------------------------------
// return a subset of right singular vectors
int PCA::dopca(vector <EMData*> imgstack, EMData *mask, int nvec)
{
   // performs PCA on a list of images (each under a mask)
   // returns a list of eigenimages

   vector<EMData*> img1dlst;
   vector<EMData*> eigvecs;
   int status = 0;

   int nimgs = imgstack.size();

   for (int i=0; i<nimgs; i++) {
      img1dlst.push_back(Util::compress_image_mask(imgstack[i],mask));
   }

   if ( nimgs < nvec || nvec == 0) nvec = nimgs;
   eigvecs = Util::svdcmp(img1dlst, nvec);
   img1dlst.clear();

   for (int i=0; i<nvec; i++) {
      eigenimages.push_back(Util::reconstitute_image_mask(eigvecs[i],mask));
   }

   eigvecs.clear();

   return status;
}

//------------------------------------------------------------------
// PCA by Lanczos
#define qmat(i,j) qmat[((j)-1)*kstep + (i) -1]
#define diag(i)   diag[(i)-1]

int PCA::dopca_lan(vector <EMData*> imgstack, EMData *mask, int nvec)
{
   // performs PCA on a list of images (each under a mask)
   // returns a list of eigenimages

   int ione = 1; 
   float one = 1.0, zero = 0.0;
   char trans;

   vector<EMData*> img1dlst;
   vector<EMData*> eigvecs;

   int status = 0;
   int nimgs = imgstack.size();
   for (int i=0; i<nimgs; i++) {
      img1dlst.push_back(Util::compress_image_mask(imgstack[i],mask));
   }

   float resnrm = 0.0;

   if ( nvec > nimgs || nvec ==0 ) nvec = nimgs;

   int nx = img1dlst[0]->get_xsize();
   // the definition of kstep is purely a heuristic for right now
   int kstep = nvec + 20;
   if (kstep > nimgs) kstep = nimgs;

   float *diag    = new float[kstep];
   float *subdiag = new float[kstep-1];
   float *vmat    = new float[nx*kstep];

   // run kstep-step Lanczos factorization
   status = Lanczos(img1dlst, &kstep, diag, subdiag, vmat, &resnrm);

   char jobz[2] = "V";
   float *qmat  = new float[kstep*kstep];
   // workspace size will be optimized later
   int   lwork  = 100 + 4*kstep + kstep*kstep;
   int   liwork = 3+5*kstep;

   float *work  = new float[lwork];
   int   *iwork = new int[liwork]; 
   int   info = 0;

   // call LAPACK tridiagonal eigensolver
   sstevd_(jobz, &kstep, diag, subdiag, qmat, &kstep, work, &lwork,
           iwork, &liwork, &info);

   // store singular values
   // eigenvalues of the cov matrix have been sorted in ascending order
   for (int j = kstep; j > kstep - nvec; j--) {
      singular_vals.push_back((float)sqrt(diag(j)));
   }

   img1dlst.clear();
   
   EMData *eigvec = new EMData();

   eigvec->set_size(nx, 1, 1);
   float *ritzvec = eigvec->get_data(); 

   // compute Ritz vectors (approximate eigenvectors) one at a time
   for (int j=0; j<nvec; j++) {
      trans = 'N';
      sgemv_(&trans, &nx, &kstep, &one, vmat, &nx, &qmat(1,kstep-j), &ione,
             &zero , ritzvec, &ione);  
      eigenimages.push_back(Util::reconstitute_image_mask(eigvec,mask));
   }

   eigvecs.clear();

   EMDeleteArray(diag);
   EMDeleteArray(subdiag);
   EMDeleteArray(vmat);
   EMDeleteArray(qmat);
   EMDeleteArray(work);
   EMDeleteArray(iwork);
   EMDeletePtr(eigvec);
   return status;
}

//------------------------------------------------------------------
// out of core version of PCA, not completed yet
int PCA::dopca_ooc(const string &filename_in, const string &filename_out, 
                   const string &lanscratch,  EMData *mask, int nvec)
{
   int status = 0, ione = 1;
   EMData *image_raw = new EMData();
   EMData *image_masked = NULL;

   int nimgs = EMUtil::get_image_count(filename_in);

   if (nimgs <= 0) {
      status = 2;
      fprintf(stderr,"dopca_ooc: no image in %s\n", filename_in.c_str());
   }
   for (int i=0; i<nimgs; i++) {
       image_raw->read_image(filename_in, i);      
       image_masked=Util::compress_image_mask(image_raw,mask);
       image_masked->write_image("temp_masked_images.img",i); 
   }

   int ndim = image_masked->get_ndim();
   if (ndim != 1) {
       fprintf(stderr,"dopca_ooc: masked image should be 1-D\n");
       status = 3; // images should all be 1-d
   }

   int nx = image_masked->get_xsize();
   
   float resnrm = 0.0;

   if ( nvec > nimgs || nvec ==0 ) nvec = nimgs;

   // the definition of kstep is purely a heuristic for right now
   int kstep = nvec + 20;
   if (kstep > nimgs) kstep = nimgs;

   float *diag    = new float[kstep];
   float *subdiag = new float[kstep-1];

   status = Lanczos_ooc("temp_masked_images.img", &kstep, diag, subdiag,
                        lanscratch , &resnrm);

   char jobz[2] = "V";
   float *qmat  = new float[kstep*kstep];
   // workspace size will be optimized later
   int   lwork  = 100 + 4*kstep + kstep*kstep;
   int   liwork = 3+5*kstep;

   float *work  = new float[lwork];
   int   *iwork = new int[liwork]; 
   int   info = 0;

   // call LAPACK tridiagonal eigensolver
   sstevd_(jobz, &kstep, diag, subdiag, qmat, &kstep, work, &lwork,
           iwork, &liwork, &info);

   // store singular values
   // eigenvalues of the cov matrix have been sorted in ascending order
   for (int j = kstep; j > kstep - nvec; j--) {
      singular_vals.push_back((float)sqrt(diag(j)));
   }

   EMData *eigvec = new EMData();
   eigvec->set_size(nx, 1, 1);
   float *ritzvec = eigvec->get_data(); 
   EMData *newimage = 0;

   // compute Ritz vectors (approximate eigenvectors) one at a time
   FILE *fp;
   float *vlan = new float[nx];
   fp = fopen(lanscratch.c_str(),"rb");
   for (int j=0; j<nvec; j++) {
      for (int i = 0; i<nx; i++) ritzvec[i]=0.0;
      for (int jlan = 1; jlan <= kstep; jlan++) {
           size_t nr = fread(vlan, sizeof(float), nx, fp); nr++;
           saxpy_(&nx, &qmat(jlan,kstep-j), vlan, &ione, ritzvec, &ione);
      }
      rewind(fp);
      newimage = Util::reconstitute_image_mask(eigvec,mask);
      newimage->write_image(filename_out,j);
   }
   fclose(fp);

   EMDeleteArray(diag);
   EMDeleteArray(subdiag);
   EMDeleteArray(qmat);
   EMDeleteArray(work);
   EMDeleteArray(iwork);
   EMDeleteArray(vlan);
   EMDeletePtr(eigvec);
   EMDeletePtr(newimage);

   return status;
}
#undef diag
#undef qmat

//------------------------------------------------------------------
#define TOL 1e-7
#define V(i,j)      V[((j)-1)*imgsize + (i) - 1]
#define v0(i)       v0[(i)-1]
#define Av(i)       Av[(i)-1]
#define subdiag(i)  subdiag[(i)-1]
#define diag(i)     diag[(i)-1]
#define hvec(i)     hvec[(i)-1]

int PCA::Lanczos(vector <EMData*> imgstack, int *kstep, 
                 float  *diag, float *subdiag, float *V, 
                 float  *beta)
{
    /*
        Purpose: Compute a kstep-step Lanczos factorization
                 on the covariant matrix X*trans(X), where 
                 X (imgstack) contains a set of images;

        Input: 
           imgstack (vector <EMData*>) a set of images on which PCA is 
                                       to be performed;
           
           kstep (int*) The maximum number of Lanczos iterations allowed.
                          If Lanczos terminates before kstep steps
                          is reached (an invariant subspace is found), 
                          kstep returns the number of steps taken;
      
        Output:
           diag (float *) The projection of the covariant matrix into a
                          Krylov subspace of dimension at most kstep.
                          The projection is a tridiagonal matrix. The
                          diagonal elements of this matrix is stored in 
                          the diag array.

           subdiag (float*) The subdiagonal elements of the projection
                            is stored here.

           V (float *)    an imgsize by kstep array that contains a 
                          set of orthonormal Lanczos basis vectors;

           beta (float *) the residual norm of the factorization;
    */
    int i, iter;
    
    float alpha;
    int   ione = 1;
    float zero = 0.0, one = 1.0, mone = -1.0;
    int   status = 0;
    
    char trans;
    int  nimgs = 0, imgsize = 0, ndim = 0; 
    float *v0, *Av, *hvec, *htmp, *imgdata;

    nimgs   = imgstack.size();
    if (nimgs <= 0) {
	status = 2; // no image in the stack
        goto EXIT; 
    }

    ndim = imgstack[0]->get_ndim();
    if (ndim != 1) {
        status = 3; // images should all be 1-d
        goto EXIT; 
    }

    imgsize = imgstack[0]->get_xsize();
     
    v0   = new float[imgsize];
    Av   = new float[imgsize];
    hvec = new float[*kstep];
    htmp = new float[*kstep];

    if (v0 == NULL || Av == NULL || hvec == NULL || htmp == NULL ) {
        fprintf(stderr, "Lanczos: failed to allocate v0,Av,hvec,htmp\n"); 
	status = -1;
        goto EXIT;
    }

    // may choose a random starting guess here     
    for ( i = 1; i <= imgsize; i++) v0(i) = 1.0;

    // normalize the starting vector
    *beta  = (float) snrm2_(&imgsize, v0, &ione);
    for (i = 1; i<=imgsize; i++)
	V(i,1) = v0(i) / (*beta);

    // do Av <-- A*v0, where A is a cov matrix
    for (i = 0; i < nimgs; i++) {
       imgdata = imgstack[i]->get_data();
       alpha = (float) sdot_(&imgsize, imgdata, &ione, V, &ione); 
       saxpy_(&imgsize, &alpha, imgdata, &ione, Av, &ione);
    }

    // Av <--- Av - V(:,1)*V(:,1)'*Av 
    diag(1) = (float) sdot_(&imgsize, V, &ione, Av, &ione); 
    alpha   = -diag(1);
    saxpy_(&imgsize, &alpha, V, &ione, Av, &ione);

    // main loop 
    for ( iter = 2 ; iter <= *kstep ; iter++ ) {
        *beta = (float) snrm2_(&imgsize, Av, &ione);

        if (*beta < TOL) {
	    // found an invariant subspace, exit
            *kstep = iter;
            break;
        }
 
        subdiag(iter-1) = *beta;
	for ( i = 1 ; i <= imgsize ; i++ ) {
	    V(i,iter) = Av(i) / (*beta);
	}	

        // do Av <-- A*V(:,iter), where A is a cov matrix
        for (i = 0; i < imgsize; i++) Av[i] = 0;
        for (i = 0; i < nimgs; i++) {
           imgdata = imgstack[i]->get_data();
           alpha = (float) sdot_(&imgsize, imgdata, &ione, &V(1,iter), &ione); 
           saxpy_(&imgsize, &alpha, imgdata, &ione, Av, &ione);
        }
	
        // f <--- Av - V(:,1:iter)*V(:,1:iter)'*Av
        trans = 'T';
        status = sgemv_(&trans, &imgsize, &iter, &one, V, &imgsize, Av, &ione,
                        &zero , hvec    , &ione); 
        trans = 'N';
        status = sgemv_(&trans, &imgsize, &iter, &mone, V, &imgsize, hvec, 
                        &ione , &one    , Av, &ione);

        // one step of reorthogonalization
        trans = 'T';
        status = sgemv_(&trans, &imgsize, &iter, &one, V, &imgsize, Av, &ione,
                        &zero , htmp    , &ione); 
        saxpy_(&iter, &one, htmp, &ione, hvec, &ione); 
        trans = 'N';
        status = sgemv_(&trans, &imgsize, &iter, &mone, V, &imgsize, htmp, 
                        &ione , &one    , Av, &ione);
        diag(iter) = hvec(iter);
    }

    EMDeleteArray(v0);
    EMDeleteArray(Av);
    EMDeleteArray(hvec);
    EMDeleteArray(htmp);

EXIT:
    return status;
}

#undef v0
#undef Av
#undef V
#undef hvec
#undef diag
#undef subdiag
#undef TOL

//------------------------------------------------------------------
#define TOL 1e-7
#define v(i)        v[(i) - 1]
#define resid(i)    resid[(i) - 1]
#define Av(i)       Av[(i)-1]
#define subdiag(i)  subdiag[(i)-1]
#define diag(i)     diag[(i)-1]

// currently reading one image at a time, can be improved by
// allocate a buffer to store a subset of images and Lanczos
// vectors

int PCA::Lanczos_ooc(string const& filename_in, int *kstep, 
    float  *diag, float *subdiag, string const& lanscratch,
    float  *beta)
{
    /*
        Purpose: Compute a kstep-step Lanczos factorization
                 on the covariant matrix X*trans(X), where 
                 X (imgstack) contains a set of images;

        Input: 
           filename_in  The name of the input file that contains
           (string)     masked images;
           
           kstep (int*) The maximum number of Lanczos iterations allowed.
                        If Lanczos terminates before kstep steps
                        is reached (an invariant subspace is found), 
                        kstep returns the number of steps taken;

           lanscratch   The name of the file that will be used to 
           (string)     store a set of orthonormal Lanczos basis vectors;

        Output:
           diag (float*)    The projection of the covariant matrix into a
                            Krylov subspace of dimension at most kstep.
                            The projection is a tridiagonal matrix. The
                            diagonal elements of this matrix is stored in 
                            the diag array;

           subdiag (float*) The subdiagonal elements of the projection
                            is stored here;

           beta (float *)   The residual norm of the factorization;
    */
    int i, iter;
    EMData *maskedimage;
    
    float alpha;
    int   ione = 1;
    int   status = 0;
    
    int   imgsize = 0, ndim = 0;
    float *v, *Av, *resid, *imgdata;
    float h = 0.0, htmp=0.0;

    int nimgs = EMUtil::get_image_count(filename_in);
    if (nimgs <= 0) {
	status = 2; // no image in the stack
        goto EXIT; 
    }

    maskedimage = new EMData();
    maskedimage->read_image(filename_in,0);
    // what happens when filename_in does not exist, or when read fails?

    ndim = maskedimage->get_ndim();
    if (ndim != 1) {
        status = 3; // images should all be 1-d
        goto EXIT; 
    }

    imgsize = maskedimage->get_xsize();
     
    v     = new float[imgsize];
    Av    = new float[imgsize];
    resid = new float[imgsize];

    if (v == NULL || Av == NULL ) {
        fprintf(stderr, "Lanczos: failed to allocate v,Av\n"); 
	status = -1;
        goto EXIT;
    }

    // may choose a random starting guess here     
    for ( i = 1; i <= imgsize; i++) v(i) = 1.0;

    // normalize the starting vector
    *beta  = (float) snrm2_(&imgsize, v, &ione);
    alpha = 1/(*beta);    
    sscal_(&imgsize, &alpha, v, &ione);
    // write v out to a scratch file
    FILE *fp;
    fp = fopen(lanscratch.c_str(), "wb");
    fwrite(v, sizeof(float), imgsize, fp);
    fclose(fp);

    // do Av <-- A*v, where A is a cov matrix
    for (i = 0; i < nimgs; i++) {
       maskedimage->read_image(filename_in,i);
       imgdata = maskedimage->get_data();
       alpha = (float) sdot_(&imgsize, imgdata, &ione, v, &ione); 
       saxpy_(&imgsize, &alpha, imgdata, &ione, Av, &ione);
    }

    // Av <--- Av - V(:,1)*V(:,1)'*Av 
    diag(1) = (float) sdot_(&imgsize, v, &ione, Av, &ione); 
    alpha   = -diag(1);
    scopy_(&imgsize, Av, &ione, resid, &ione);
    saxpy_(&imgsize, &alpha, v, &ione, resid, &ione);

    // main loop 
    for ( iter = 2 ; iter <= *kstep ; iter++ ) {
        *beta = (float) snrm2_(&imgsize, resid, &ione);

        if (*beta < TOL) {
	    // found an invariant subspace, exit
            *kstep = iter;
            break;
        }
 
        subdiag(iter-1) = *beta;
	for ( i = 1 ; i <= imgsize ; i++ ) {
	    v(i) = resid(i) / (*beta);
	}	

        // write v out to a scratch file at appropriate position
        fp = fopen(lanscratch.c_str(),"ab");
        fwrite(v, sizeof(float), imgsize, fp);
        fclose(fp);

        // do Av <-- A*V(:,iter), where A is a cov matrix
        for (i = 0; i < imgsize; i++) Av[i] = 0;
        for (i = 0; i < nimgs; i++) {
           maskedimage->read_image(filename_in,i);
           imgdata = maskedimage->get_data();
           alpha = (float) sdot_(&imgsize, imgdata, &ione, v, &ione); 
           saxpy_(&imgsize, &alpha, imgdata, &ione, Av, &ione);
        }

        // f <--- Av - V(:,1:iter)*V(:,1:iter)'*Av
        // the out-of-core version reads one Lanczos vector at a time
        scopy_(&imgsize, Av, &ione, resid, &ione);
        fp = fopen(lanscratch.c_str(),"rb");
        for (int jlan = 1; jlan <= iter; jlan++) {
           size_t nr = fread(v, sizeof(float), imgsize, fp); nr++;
           h     = (float) sdot_(&imgsize, v, &ione, Av, &ione);
           alpha = -h;
           saxpy_(&imgsize, &alpha, v, &ione, resid, &ione);
        }
        fclose(fp);

        // one step of reorthogonalization
        scopy_(&imgsize, resid, &ione, Av, &ione);
        fp = fopen(lanscratch.c_str(),"rb");
        for (int jlan = 1; jlan <= iter; jlan++) {
           size_t nr = fread(v, sizeof(float), imgsize, fp); nr++;
           htmp  = (float) sdot_(&imgsize, v, &ione, Av, &ione);
           alpha = -htmp;
           saxpy_(&imgsize, &alpha, v, &ione, resid, &ione);
           h += htmp;
        }
        fclose(fp);
        diag(iter) = h;
    }

    EMDeleteArray(v);
    EMDeleteArray(Av);
    EMDeleteArray(resid);
    EMDeletePtr(maskedimage);

EXIT:
    return status;
}

#undef Av
#undef v
#undef resid
#undef diag
#undef subdiag
#undef TOL

//------------------------------------------------------------------
vector<float> PCA::get_vals()
{
   // method for retrieving singular values
   return singular_vals;
}

vector<EMData*> PCA::get_vecs()
{
   // method for retrieving right singular vectors (eigenimages)
   return eigenimages;
}

void PCA::clear()
{
    // flush singular values and vectors
   singular_vals.clear();
   eigenimages.clear();
}
