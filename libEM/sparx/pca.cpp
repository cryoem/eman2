#include "emdata.h"
#include "util.h"
#include "emutil.h"

#include "pca.h"
#include "lapackblas.h"

using namespace EMAN;

// return all right singular vectors
vector <EMData*> PCA::dopca(vector <EMData*> imgstack, EMData *mask)
{
   // performs PCA on a list of images (each under a mask)
   // returns a list of eigenimages

   int i;

   vector<EMData*> img1dlst;
   vector<EMData*> eigvecs;
   vector<EMData*> eigenimages;

   int nimgs = imgstack.size();
   // printf("number of images in the stack = %d\n", nimgs);

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

   return eigenimages;
}

// return a subset of right singular vectors
vector <EMData*> PCA::dopca(vector <EMData*> imgstack, EMData *mask, int nvec)
{
   // performs PCA on a list of images (each under a mask)
   // returns a list of eigenimages

   int i;

   vector<EMData*> img1dlst;
   vector<EMData*> eigvecs;
   vector<EMData*> eigenimages;

   int nimg = imgstack.size();
   // printf("number of images in the stack = %d\n", nimg);

   for (i=0; i<nimg; i++) {
      img1dlst.push_back(Util::compress_image_mask(imgstack[i],mask));
   }

   // for right now, compute a full SVD
   eigvecs = Util::svdcmp(img1dlst, 0);

   img1dlst.clear();

   if (nimg < nvec) nvec = nimg;

   for (i=0; i<nvec; i++) {
      eigenimages.push_back(Util::reconstitute_image_mask(eigvecs[i],mask));
   }

   eigvecs.clear();

   return eigenimages;
}

// out of core version of PCA, not completed yet
char* PCA::dopca_ooc(const string &filename_in, EMData *mask, int nvec)
{
   char *filename_out = NULL;
   EMData *image_raw = new EMData();
   EMData *image_masked;

   int nimgs = EMUtil::get_image_count(filename_in);
   for (int i=0; i<nimgs; i++) {
       image_raw->read_image(filename_in, i);      
       image_masked=Util::compress_image_mask(image_raw,mask);
       image_masked->write_image("temp_masked_imaged.img",i); 
   }

   return filename_out;
}

#define TOL 1e-7
#define V(i,j) V[((j)-1)*imgsize + (i) - 1]
#define T(i,j) T[((j)-1)*(*maxiter) + (i) - 1]
#define v0(i)  v0[(i)-1]
#define Av(i)  Av[(i)-1]

int Lanczos(vector <EMData*> imgstack, int *maxiter, 
            float  *T, float *V, float *beta)
{
    /*
        Purpose: Compute a maxiter-step Lanczos factorization
                 on the covariant matrix X*trans(X), where 
                 X (imgstack) contains a set of images;

        Input: 
           imgstack (vector <EMData*>) a set of images on which PCA is 
                                       to be performed;
           
           maxiter (int*) The maximum number of Lanczos iterations allowed.
                          If Lanczos terminates before maxiter steps
                          is reached (an invariant subspace is found), 
                          maxiter returns the number of steps taken;
      
        Output:
           T (float *)    The projection of the covariant matrix into a
                          Krylov subspace of dimension at most maxiter.
                          The matrix T is declared as a maxiter*maxiter
                          1D float array. It should be a symmetric 
                          tridiagonal matrix in exact arithmetic;

           V (float *)    an imgsize by maxiter array that contains a 
                          set of orthonormal Lanczos basis vectors;

           beta (float *) the residual norm of the factorization;
    */
    int i, j;
    
    float alpha;
    int   ione = 1;
    float zero = 0.0, one = 1.0, mone = -1.0;
    int   status = 0;
    
    float *v0, *Av, *imgdata;
    char trans;
    int  nimgs = 0, imgsize = 0, ndim = 0; 

    nimgs   = imgstack.size();
    if (nimgs <= 0) {
	status = 2; // no image in the stack
        goto EXIT; 
    }

    ndim    = imgstack[0]->get_ndim();
    if (ndim != 1) {
        status = 3; // images should all be 1-d
        goto EXIT; 
    }

    imgsize = imgstack[0]->get_xsize();
     
    v0  = (float *) calloc(imgsize, sizeof(float));
    Av  = (float *) calloc(imgsize, sizeof(float));

    if (v0 == NULL || Av == NULL) {
	status = -1;
        goto EXIT;
    }

    // may choose a random starting guess here     
    for ( i = 1; i <= imgsize; i++) v0(i) = 1.0;

    // normalize the starting vector
    *beta  = snrm2_(&imgsize, v0, &ione);
    for (i = 1; i<=imgsize; i++)
	V(i,1) = v0(i) / (*beta);

    // do Av <-- A*v0, where A is a cov matrix
    for (i = 0; i < nimgs; i++) {
       imgdata = imgstack[i]->get_data();
       alpha = sdot_(&imgsize, imgdata, &ione, v0, &ione); 
       saxpy_(&imgsize, &alpha, imgdata, &ione, Av, &ione);
    }

    // Av <--- Av - V(:,1)*V(:,1)'*Av 
    T(1,1) = sdot_(&imgsize, V, &ione, Av, &ione); 
    alpha = -T(1,1);
    saxpy_(&imgsize, &alpha, V, &ione, Av, &ione);

    // main loop 
    for ( j = 2 ; j <= *maxiter ; j++ ) {
        *beta = snrm2_(&imgsize, Av, &ione);
        if (*beta < TOL) {
	    // found an invariant subspace, exit
            *maxiter = j;
            break;
        }
 
        T(j,j-1) = *beta;
	for ( i = 1 ; i <= imgsize ; i++ ) {
	    V(i,j) = Av(i) / T(j,j-1);
	}	

        // do Av <-- A*V(:,j), where A is a cov matrix
        for (i = 0; i < imgsize; i++) Av[i] = 0;
        for (i = 0; i < nimgs; i++) {
           imgdata = imgstack[i]->get_data();
           alpha = sdot_(&imgsize, imgdata, &ione, &V(1,j), &ione); 
           saxpy_(&imgsize, &alpha, imgdata, &ione, Av, &ione);
        }
	
        // f <--- Av - V(:,1:j)*V(:,1:j)'*Av
        trans = 'T';
        status = sgemv_(&trans, &imgsize, &j, &one, V, &imgsize, Av, &ione,
                        &zero , &T(1,j) , &ione); 
        trans = 'N';
        status = sgemv_(&trans, &imgsize, &j, &mone, V, &imgsize, &T(1,j), 
                        &ione , &one    , Av, &ione);
    }

    free(v0);
    free(Av);

EXIT:
    return status;
}

#undef v0
#undef Av
#undef V
#undef T
