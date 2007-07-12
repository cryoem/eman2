/**
 * $Id$
 * */


/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
 * Copyright (c) 2000-2006 Baylor College of Medicine
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
 * */

#include "emdata.h"
#include "analyzer.h"
#include "sparx/analyzer_sparx.h"
#include "util.h"
#include "sparx/lapackblas.h"
#include "sparx/varimax.h"

using namespace EMAN;

namespace EMAN {

	template <> Factory < Analyzer >::Factory()
	{
		force_add(&PCAsmall::NEW);
		force_add(&PCAlarge::NEW);
		force_add(&varimax::NEW);
	}

}

#define covmat(i,j) covmat[ ((j)-1)*nx + (i)-1 ]
#define imgdata(i)  imgdata[ (i)-1 ]
int PCAsmall::insert_image(EMData * image)
{
   EMData *maskedimage = Util::compress_image_mask(image,mask);

   int nx = maskedimage->get_xsize();
   float *imgdata = maskedimage->get_data();
   if (nx != ncov) {
      fprintf(stderr,"insert_image: something is wrong...\n");
      exit(1);
   }

   // there is a faster version of the following rank-1 update 
   nimages++;
   for (int j = 1; j <= nx; j++)
       for (int i = 1; i<=nx; i++) {
           covmat(i,j) += imgdata(i)*imgdata(j);
   }   

   EMDeletePtr(maskedimage);
   return 0;
}
#undef covmat

#define eigvec(i,j) eigvec[(j)*ncov + (i)]
vector<EMData*> PCAsmall::analyze()
{
        float *eigvec;
	int status = 0;
	printf("start analyzing..., ncov = %d\n", ncov);
        eigval = (float*)calloc(ncov,sizeof(float));
        eigvec = (float*)calloc(ncov*ncov,sizeof(float));
        status = Util::coveig(ncov, covmat, eigval, eigvec);

        for (int i=1; i<=nvec; i++) printf("eigval = %11.4e\n", 
            eigval[ncov-i]);

        // pack eigenvectors into the return imagelist
        EMData *eigenimage = new EMData();
        eigenimage->set_size(ncov,1,1);
        float *rdata = eigenimage->get_data();
        for (int j = 1; j<= nvec; j++) {
	    for (int i = 0; i < ncov; i++)
		rdata[i] = eigvec(i,ncov-j);
	    images.push_back(Util::reconstitute_image_mask(eigenimage,mask));
        }

        free(eigvec);
        EMDeletePtr(eigenimage); 

	return images;
}
#undef eigvec

void PCAsmall::set_params(const Dict & new_params)
{
	params = new_params;
	mask = params["mask"];
	nvec = params["nvec"];

        // count the number of pixels under the mask
        // (this is really ugly!!!)
        EMData *dummy = new EMData();

        int nx = mask->get_xsize();
        int ny = mask->get_ysize();
        int nz = mask->get_zsize();

        dummy->set_size(nx,ny,nz);

        EMData *dummy1d = Util::compress_image_mask(dummy,mask);
        ncov = dummy1d->get_xsize();
        EMDeletePtr(dummy);
        EMDeletePtr(dummy1d);

	// allocate and set up the covriance matrix
        nimages = 0;
	covmat = (float*)calloc(ncov*ncov,sizeof(float));
}

//------------------------------------------------------------------
// for large-scale PCA incore

int PCAlarge::insert_image(EMData * image)
{
   EMData *maskedimage = Util::compress_image_mask(image,mask);

   FILE *fp;
   string scratchfile = string("maskedimages.scratch");

   fp = fopen(scratchfile.c_str(),"a");

   int nx = maskedimage->get_xsize();
   float *imgdata = maskedimage->get_data();
   fwrite(imgdata, sizeof(float), nx, fp);
   nimages++;

   fclose(fp);

   EMDeletePtr(maskedimage);

   return 0;
}

void PCAlarge::set_params(const Dict & new_params)
{
	params = new_params;
	mask = params["mask"];
	nvec = params["nvec"];

        // count the number of pixels under the mask
        // (this is really ugly!!!)
        EMData *dummy = new EMData();

        int nx = mask->get_xsize();
        int ny = mask->get_ysize();
        int nz = mask->get_zsize();

        dummy->set_size(nx,ny,nz);

        EMData *dummy1d = Util::compress_image_mask(dummy,mask);

        ncov = dummy1d->get_xsize();

        EMDeletePtr(dummy);
        EMDeletePtr(dummy1d);
        // no need to allocate the covariance matrix
        nimages = 0;
}

#define qmat(i,j)   qmat[((j)-1)*kstep + (i) -1]
#define diag(i)     diag[(i)-1]
#define rdata(i)    rdata[(i)-1]
#define eigvec(i,j) eigvec[((j)-1)*ncov + (i)-1]
#define eigval(i)   eigval[(i)-1]

vector<EMData*> PCAlarge::analyze()
{
	int status = 0;
	int ione = 1; 
	float one = 1.0, zero = 0.0;
	char trans;
        float *eigvec;
        string scratchfile = string("maskedimages.scratch");
        char command[100];

	printf("start analyzing..., ncov = %d\n", ncov);

        float resnrm = 0.0;

        if ( nvec > nimages || nvec ==0 ) nvec = nimages;
        int nx = ncov;

        // the definition of kstep is purely a heuristic for right now
	int kstep = nvec + 20;
	if (kstep > nimages) kstep = nimages;

	float *diag    = new float[kstep];
	float *subdiag = new float[kstep-1];
	float *vmat    = new float[nx*kstep];

        // run kstep-step Lanczos factorization
	status = Lanczos(scratchfile, &kstep, diag, subdiag, 
                         vmat, &resnrm);

        // remove scratch file
        sprintf(command,"rm -f %s\n", scratchfile.c_str());
        status = system(command);
        if (status != 0) {
	    fprintf(stderr,"PCAlarge: cannot remove scratchfile\n");
        }

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

        // store eigenvalues
        eigval = (float*)calloc(ncov,sizeof(float));
        eigvec = (float*)calloc(ncov*nvec,sizeof(float));

	for (int j = 0; j < nvec; j++) {
            eigval[j] = diag(kstep-j);
        }

        for (int i=0; i<nvec; i++) printf("eigval = %11.4e\n", 
            eigval[i]);

        // compute eigenvectors
        for (int j=1; j<=nvec; j++) {
            trans = 'N';
            sgemv_(&trans, &nx,  &kstep, &one, vmat, &nx, &qmat(1,kstep-j+1), 
                   &ione, &zero, &eigvec(1,j), &ione);  
        }

        // pack eigenvectors into the return imagelist
        EMData *eigenimage = new EMData();
        eigenimage->set_size(ncov,1,1);
        float *rdata = eigenimage->get_data();
        for (int j = 1; j<= nvec; j++) {
	    for (int i = 1; i <= ncov; i++)
		rdata(i) = eigvec(i,j);
  
            EMData* recons_eigvec = Util::reconstitute_image_mask(eigenimage,mask);

            recons_eigvec->set_attr( "eigval", eigval[j-1] );

	    images.push_back( recons_eigvec );
        }

        free(eigvec);
        EMDeletePtr(eigenimage); 

	return images;
}
#undef qmat
#undef diag
#undef rdata
#undef eigvec
#undef eigval

#define TOL 1e-7
#define V(i,j)      V[((j)-1)*imgsize + (i) - 1]
#define v0(i)       v0[(i)-1]
#define Av(i)       Av[(i)-1]
#define subdiag(i)  subdiag[(i)-1]
#define diag(i)     diag[(i)-1]
#define hvec(i)     hvec[(i)-1]

int PCAlarge::Lanczos(const string &maskedimages, int *kstep, 
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
    int  imgsize = 0; 
    float *v0, *Av, *hvec, *htmp, *imgdata;
    FILE  *fp=NULL;

    if (nimages <= 0) {
	status = 2; // no image in the stack
        goto EXIT; 
    }

    imgsize = ncov;
    if (nimages <= 0) {
	status = 3; // no image in the stack
        goto EXIT; 
    }
     
    v0   = new float[imgsize];
    Av   = new float[imgsize];
    hvec = new float[*kstep];
    htmp = new float[*kstep];
    imgdata = new float[imgsize];

    if (v0 == NULL || Av == NULL || hvec == NULL || 
        htmp == NULL || imgdata == NULL) {
        fprintf(stderr, "Lanczos: failed to allocate v0,Av,hvec,htmp\n"); 
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
    fp = fopen(maskedimages.c_str(),"r");
    if (fp==NULL) {
	fprintf(stderr,"Lanczos: cannot open %s\n", maskedimages.c_str());
    }
    for (i = 0; i < nimages; i++) {
       fread(imgdata, sizeof(float), imgsize, fp);
       alpha = sdot_(&imgsize, imgdata, &ione, V, &ione); 
       saxpy_(&imgsize, &alpha, imgdata, &ione, Av, &ione);
    }
    fclose(fp);


    // Av <--- Av - V(:,1)*V(:,1)'*Av 
    diag(1) = sdot_(&imgsize, V, &ione, Av, &ione); 
    alpha   = -diag(1);
    saxpy_(&imgsize, &alpha, V, &ione, Av, &ione);

    // main loop 
    for ( iter = 2 ; iter <= *kstep ; iter++ ) {
        *beta = snrm2_(&imgsize, Av, &ione);

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
        fp = fopen(maskedimages.c_str(),"r");
        for (i = 0; i < nimages; i++) {
           fread(imgdata, sizeof(float), imgsize, fp);
           alpha = sdot_(&imgsize, imgdata, &ione, &V(1,iter), &ione); 
           saxpy_(&imgsize, &alpha, imgdata, &ione, Av, &ione);
        }
        fclose(fp); 
	
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
    EMDeleteArray(imgdata);

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

void varimax::set_params(const Dict & new_params)
{
	params = new_params;
	m_mask = params["mask"];

        // count the number of pixels under the mask
        // (this is really ugly!!!)
        EMData *dummy = new EMData();

        int nx = m_mask->get_xsize();
        int ny = m_mask->get_ysize();
        int nz = m_mask->get_zsize();

        dummy->set_size(nx,ny,nz);

        EMData *dummy1d = Util::compress_image_mask(dummy,m_mask);

        m_nlen = dummy1d->get_xsize();
        m_nfac = 0;

        EMDeletePtr(dummy);
        EMDeletePtr(dummy1d);
}

int varimax::insert_image(EMData* image)
{
    EMData* img1d = Util::compress_image_mask(image,m_mask);

    m_data.insert( m_data.end(), img1d->get_data(), img1d->get_data() + m_nlen );

    m_nfac++;

    Assert( (int)m_data.size() == m_nfac*m_nlen);
    
    return 0;
}

vector<EMData*> varimax::analyze()
{
    int itmax = 10000;
    float eps = 1e-4;
    int verbose = 1;
    float params[4];
    params[0] = 1.0;
    varmx( &m_data[0], m_nlen, m_nfac, IVARIMAX, params, NULL, itmax, eps, verbose);

    vector<EMData*> images;

    EMData* img1d = new EMData();
    img1d->set_size(m_nlen, 1, 1);
    for( int i=0; i < m_nfac; ++i )
    {
	float* imgdata = img1d->get_data();

	int offset = i * m_nlen;
	for( int i=0; i < m_nlen; ++i )
	{
	    imgdata[i] = m_data[offset+i];
        }

        EMData* img = Util::reconstitute_image_mask(img1d,m_mask);
	images.push_back(img);
    }

    EMDeletePtr(img1d);

    return images;
}

void EMAN::dump_analyzers()
{
	dump_factory < Analyzer > ();
}

map<string, vector<string> > EMAN::dump_analyzers_list()
{
	return dump_factory_list < Analyzer > ();
}







