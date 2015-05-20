/**
 * $Id$
 */


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

#include <ctime>
#include <memory>
#include "emdata.h"
#include "analyzer.h"
#include "sparx/analyzer_sparx.h"
#include "util.h"
#include "cmp.h"
#include "sparx/lapackblas.h"
#include "sparx/varimax.h"

using namespace EMAN;

namespace EMAN {

	const string PCAsmall::NAME = "pca";
	const string PCAlarge::NAME = "pca_large";
	const string varimax::NAME = "varimax";
	const string InertiaMatrixAnalyzer::NAME = "inertiamatrix";
	const string ShapeAnalyzer::NAME = "shape";
	const string KMeansAnalyzer::NAME = "kmeans";
	const string SVDAnalyzer::NAME = "svd_gsl";
	const string CircularAverageAnalyzer::NAME = "cir_avg";

	template <> Factory < Analyzer >::Factory()
	{
		force_add<PCAsmall>();
		force_add<PCAlarge>();
		force_add<varimax>();
 		force_add<InertiaMatrixAnalyzer>();
		force_add<ShapeAnalyzer>();
 		force_add<KMeansAnalyzer>();
 		force_add<SVDAnalyzer>();
 		force_add<CircularAverageAnalyzer>();
	}

}

int Analyzer::insert_images_list(vector<EMData *> image_list)
{
	vector<EMData *>::const_iterator iter;
		for(iter=image_list.begin(); iter!=image_list.end(); ++iter) {
			insert_image(*iter);
		}
	return 0;
}

vector<EMData *> InertiaMatrixAnalyzer::analyze() {
	int verbose = params.set_default("verbose",0);
	EMData *mx = new EMData(3,3);	// result is a 3x3 matrix
	mx->to_zero();
	ret.push_back(mx);

	if (images.size()!=1) throw ImageDimensionException("Inertia matrix computation accepts only a single volume as input");
	int nx=images[0]->get_xsize();
	int ny=images[0]->get_ysize();
	int nz=images[0]->get_zsize();
	if (nz==1 || ny==1 || nz==1) throw ImageDimensionException("Map must be 3-D");

	if (verbose>0) printf("Inertia volume size: %d %d %d\n",nx,ny,nz);

	for (int z=0; z<nz; z++) {
		for (int y=0; y<ny; y++) {
			for (int x=0; x<nx; x++) {
				int xx=x-nx/2;
				int yy=y-ny/2;
				int zz=z-nz/2;
				float v=images[0]->get_value_at(x,y,z);
				mx->set_value_at(0,0,mx->get_value_at(0,0)+v*(yy*yy+zz*zz));
				mx->set_value_at(0,1,mx->get_value_at(0,1)+v*(-xx*yy));
				mx->set_value_at(0,2,mx->get_value_at(0,2)+v*(-xx*zz));
				mx->set_value_at(1,0,mx->get_value_at(1,0)+v*(-xx*yy));
				mx->set_value_at(1,1,mx->get_value_at(1,1)+v*(zz*zz+xx*xx));
				mx->set_value_at(1,2,mx->get_value_at(1,2)+v*(-yy*zz));
				mx->set_value_at(2,0,mx->get_value_at(2,0)+v*(-xx*zz));
				mx->set_value_at(2,1,mx->get_value_at(2,1)+v*(-yy*zz));
				mx->set_value_at(2,2,mx->get_value_at(2,2)+v*(xx*xx+yy*yy));
			}
		}
	}
	mx->mult(1.0f/(nx*ny*nz));

	if (verbose>0) {
		printf("%1.3g\t%1.3g\t%1.3g\n",mx->get_value_at(0,0),mx->get_value_at(1,0),mx->get_value_at(2,0));
		printf("%1.3g\t%1.3g\t%1.3g\n",mx->get_value_at(0,1),mx->get_value_at(1,1),mx->get_value_at(2,1));
		printf("%1.3g\t%1.3g\t%1.3g\n",mx->get_value_at(0,2),mx->get_value_at(1,2),mx->get_value_at(2,2));
	}

	return ret;
}

vector<EMData *> ShapeAnalyzer::analyze() {
	int verbose = params.set_default("verbose",0);
	EMData *mx = new EMData(4,2,1);	// result is 4 values
	mx->to_zero();
	ret.push_back(mx);

	if (images.size()!=1) throw ImageDimensionException("Shape computation accepts only a single volume as input");
	int nx=images[0]->get_xsize();
	int ny=images[0]->get_ysize();
	int nz=images[0]->get_zsize();
	if (nz==1 || ny==1 || nz==1) throw ImageDimensionException("Map must be 3-D");

	if (verbose>0) printf("Shape size: %d %d %d\n",nx,ny,nz);

	for (int z=0; z<nz; z++) {
		for (int y=0; y<ny; y++) {
			for (int x=0; x<nx; x++) {
				int xx=x-nx/2;
				int yy=y-ny/2;
				int zz=z-nz/2;
				float v=images[0]->get_value_at(x,y,z);
				mx->set_value_at(0,0,mx->get_value_at(0,0)+v*(xx*xx));
				mx->set_value_at(1,0,mx->get_value_at(1,0)+v*(yy*yy));
				mx->set_value_at(2,0,mx->get_value_at(2,0)+v*(zz*zz));
				mx->set_value_at(3,0,mx->get_value_at(3,0)+v*(xx*xx+yy*yy+zz*zz)); 
				// sum(m*r^2), in which r is the distance to the center. Used for minicircle classification
				mx->set_value_at(0,1,mx->get_value_at(0,0)+v*abs(xx));
				mx->set_value_at(1,1,mx->get_value_at(1,0)+v*abs(yy));
				mx->set_value_at(2,1,mx->get_value_at(2,0)+v*abs(zz));
				mx->set_value_at(3,1,mx->get_value_at(3,0)+v*(float)sqrt((float)(xx*xx+yy*yy+zz*zz)));
			}
		}
	}
	mx->mult(1.0f/(nx*ny*nz));

	return ret;
}


void KMeansAnalyzer::set_params(const Dict & new_params)
{
	params = new_params;
	if (params.has_key("ncls")) ncls = params["ncls"];
	if (params.has_key("maxiter"))maxiter = params["maxiter"];
	if (params.has_key("minchange"))minchange = params["minchange"];
	if (params.has_key("mininclass"))mininclass = params["mininclass"];
	if (params.has_key("slowseed"))slowseed = params["slowseed"];
	if (params.has_key("verbose"))verbose = params["verbose"];
	if (params.has_key("calcsigmamean")) calcsigmamean=params["calcsigmamean"];

}

vector<EMData *> KMeansAnalyzer::analyze()
{
if (ncls<=1) return vector<EMData *>();
//srandom(time(0));

// These are the class centers, start each with a random image
int nptcl=images.size();
int nclstot=ncls;
if (calcsigmamean) centers.resize(nclstot*2);
else centers.resize(nclstot);
if (mininclass<1) mininclass=1;

for (int i=0; i<nptcl; i++) images[i]->set_attr("is_ok_center",(int)5);  // if an image becomes part of too small a set, it will (eventually) be marked as a bad center

if (slowseed) {
	if (ncls>25) slowseed=ncls/25+1;	// this becomes the number to seed in each step
//	if (maxiter<ncls*3+20) maxiter=ncls*3+20;	// We need to make sure we have enough iterations to seed all of the classes
//	ncls=2;
}

for (int i=0; i<ncls; i++) {
	// Fixed by d.woolford, Util.get_irand is inclusive (added a -1)
	centers[i]=images[Util::get_irand(0,nptcl-1)]->copy();

}

if (calcsigmamean) {
	for (int i=nclstot; i<nclstot*2; i++) centers[i]=new EMData(images[0]->get_xsize(),images[0]->get_ysize(),images[0]->get_zsize());
}


for (int i=0; i<maxiter; i++) {
	nchanged=0;
	reclassify();
	if (verbose) printf("iter %d>  %d (%d)\n",i,nchanged,ncls);
	if (nchanged<minchange && ncls==nclstot) break;
	update_centers();

	if (slowseed && i%3==2 && ncls<nclstot) {
		for (int j=0; j<slowseed && ncls<nclstot; j++) {
			centers[ncls]=0;
			ncls++;
		}
		reseed();
	}
}
update_centers(calcsigmamean);

return centers;
}

void KMeansAnalyzer::update_centers(int sigmas) {
int nptcl=images.size();
//int repr[ncls];
int * repr = new int[ncls];

for (int i=0; i<ncls; i++) {
	centers[i]->to_zero();
	if (sigmas) centers[i+ncls]->to_zero();
	repr[i]=0;
}

// compute new position for each center
for (int i=0; i<nptcl; i++) {
	int cid=images[i]->get_attr("class_id");
	if ((int)images[i]->get_attr("is_ok_center")>0) {
		centers[cid]->add(*images[i]);
		if (sigmas) centers[cid+ncls]->addsquare(*images[i]);
		repr[cid]++;
	}
}

for (int i=0; i<ncls; i++) {
	// If this class is too small
	if (repr[i]<mininclass) {
		// find all of the particles in the class, and decrement their "is_ok_center" counter.
		// when it reaches zero the particle will no longer participate in determining the location of a center
		for (int j=0; j<nptcl; j++) {
			if ((int)images[j]->get_attr("class_id")==i) images[i]->set_attr("is_ok_center",(int)images[i]->get_attr("is_ok_center")-1);
		}
		// Mark the center for reseeding
		delete centers[i];
		centers[i]=0;
		repr[i]=0;
	}
	// finishes off the statistics we started computing above
	else {
		centers[i]->mult((float)1.0/(float)(repr[i]));
		centers[i]->set_attr("ptcl_repr",repr[i]);
		if (sigmas) {
			centers[i+ncls]->mult((float)1.0/(float)(repr[i]));		// sum of squares over n
			centers[i+ncls]->subsquare(*centers[i]);					// subtract the mean value squared
			centers[i+ncls]->process("math.sqrt");					// square root
			centers[i+ncls]->mult((float)1.0/(float)sqrt((float)repr[i]));		// divide by sqrt(N) to get std. dev. of mean
		}

	}
	if (verbose>1) printf("%d(%d)\t",i,(int)repr[i]);
}

if (verbose>1) printf("\n");

reseed();

delete [] repr;
}

// This will look for any unassigned points and reseed each inside the class with the broadest distribution widely distributed
void KMeansAnalyzer::reseed() {
int nptcl=images.size();
int i,j;

// if no classes need reseeding just return
for (i=0; i<ncls; i++) {
	if (!centers[i]) break;
}
if (i==ncls) return;

// make a list of all particles which could be centers
vector<int> goodcen;
for (int i=0; i<nptcl; i++) if ((int)images[i]->get_attr("is_ok_center")>0) goodcen.push_back(i);

if (goodcen.size()==0) throw UnexpectedBehaviorException("Kmeans ran out of valid center particles with the provided parameters");

// pick a random particle for the new seed
for (i=0; i<ncls; i++) {
	if (centers[i]) continue;		// center doesn't need reseeding
	j=Util::get_irand(0,goodcen.size()-1);
	centers[i]=images[j]->copy();
	centers[i]->set_attr("ptcl_repr",1);
	printf("reseed %d -> %d\n",i,j);
}


}


void KMeansAnalyzer::reclassify() {
int nptcl=images.size();

Cmp *c = Factory < Cmp >::get("sqeuclidean");
for (int i=0; i<nptcl; i++) {
	float best=1.0e38f;
	int bestn=0;
	for (int j=0; j<ncls; j++) {
		float d=c->cmp(images[i],centers[j]);
//images[i]->cmp("sqeuclidean",centers[j]);
		if (d<best) { best=d; bestn=j; }
	}
	int oldn=images[i]->get_attr_default("class_id",0);
	if (oldn!=bestn) nchanged++;
	images[i]->set_attr("class_id",bestn);
}
delete c;
}

#define covmat(i,j) covmat[ ((j)-1)*nx + (i)-1 ]
#define imgdata(i)  imgdata[ (i)-1 ]
int PCAsmall::insert_image(EMData * image)
{
	if(mask==0)
		throw NullPointerException("Null mask image pointer, set_params() first");

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
//		printf("start analyzing..., ncov = %d\n", ncov);
        eigval = (float*)calloc(ncov,sizeof(float));
        eigvec = (float*)calloc(ncov*ncov,sizeof(float));
        status = Util::coveig(ncov, covmat, eigval, eigvec);
//       for (int i=1; i<=nvec; i++) printf("eigval = %11.4e\n",
//            eigval[ncov-i]);

        // pack eigenvectors into the return imagelist
        EMData *eigenimage = new EMData();
        eigenimage->set_size(ncov,1,1);
        float *rdata = eigenimage->get_data();
        for (int j = 1; j<= nvec; j++) {
	    for (int i = 0; i < ncov; i++) rdata[i] = eigvec(i,ncov-j);

		EMData* recons_eigvec = Util::reconstitute_image_mask(eigenimage,mask);
		recons_eigvec->set_attr( "eigval", eigval[j-1] );
	    images.push_back(recons_eigvec);
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
	if(mask==0)
		throw NullPointerException("Null mask image pointer, set_params() first");

   EMData *maskedimage = Util::compress_image_mask(image,mask);

   FILE *fp;
   string scratchfile = string("maskedimages.scratch");

   fp = fopen(scratchfile.c_str(),"ab");

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
#ifdef _WIN32
	if (_unlink(scratchfile.c_str()) == -1) {
		fprintf(stderr,"PCAlarge: cannot remove scratchfile\n");
	}
#else
	sprintf(command,"rm -f %s\n", scratchfile.c_str());
	status = system(command);
	if (status != 0) {
		fprintf(stderr,"PCAlarge: cannot remove scratchfile\n");
	}
#endif	//_WIN32

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

//         for (int i=0; i<nvec; i++) printf("eigval = %11.4e\n",
//             eigval[i]);

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
    for ( i = 1; i <= imgsize; i++)
    {
        v0(i) = 1.0;
        Av(i) = 0.0;
    }

    // normalize the starting vector
    *beta  = snrm2_(&imgsize, v0, &ione);
    for (i = 1; i<=imgsize; i++)
	V(i,1) = v0(i) / (*beta);

    // do Av <-- A*v0, where A is a cov matrix
    fp = fopen(maskedimages.c_str(),"rb");
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
        fp = fopen(maskedimages.c_str(),"rb");
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
	if(m_mask==0)
		throw NullPointerException("Null mask image pointer, set_params() first");

    EMData* img1d = Util::compress_image_mask(image,m_mask);

    m_data.insert( m_data.end(), img1d->get_data(), img1d->get_data() + m_nlen );

    m_nfac++;

    Assert( (int)m_data.size() == m_nfac*m_nlen);

    return 0;
}

vector<EMData*> varimax::analyze()
{
    int itmax = 10000;
    float eps = 1e-4f;
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

int SVDAnalyzer::insert_image(EMData * image)
{
	if (mask==0)
		throw NullPointerException("Null mask image pointer, set_params() first");

	// count pixels under mask
	size_t totpix=mask->get_xsize()*mask->get_ysize()*mask->get_zsize();
	float  *d=image->get_data();
	float *md=mask ->get_data();
	for (size_t i=0,j=0; i<totpix; ++i) {
		if (md[i]) {
			gsl_matrix_set(A,j,nsofar,d[i]);
			j++;
		}
	}
	nsofar++;

   return 0;
}
#undef covmat

#define eigvec(i,j) eigvec[(j)*ncov + (i)]
vector<EMData*> SVDAnalyzer::analyze()
{
// Allocate the working space
gsl_vector *work=gsl_vector_alloc(nimg);
gsl_vector *S=gsl_vector_alloc(nimg);
gsl_matrix *V=gsl_matrix_alloc(nimg,nimg);
gsl_matrix *X=gsl_matrix_alloc(nimg,nimg);


// Do the decomposition. All the real work is here
gsl_linalg_SV_decomp_mod (A,X, V, S, work);
//else gsl_linalg_SV_decomp_jacobi(A,V,S);

vector<EMData*> ret;
//unpack the results and write the output file
float *md=mask->get_data();
size_t totpix=mask->get_xsize()*mask->get_ysize()*mask->get_zsize();
for (int k=0; k<nvec; k++) {
	EMData *img = new EMData;
	img->set_size(mask->get_xsize(),mask->get_ysize(),mask->get_zsize());

	float  *d=img->get_data();
	for (size_t i=0,j=0; i<totpix; ++i) {
		if (md[i]) {
			d[i]=(float)gsl_matrix_get(A,j,k);
			j++;
		}
	}
	img->set_attr( "eigval", gsl_vector_get(S,k));
	ret.push_back(img);
}

gsl_vector_free(work);
gsl_vector_free(S);
gsl_matrix_free(V);
gsl_matrix_free(X);

gsl_matrix_free(A);
A=NULL;
mask=NULL;

return ret;
}

void SVDAnalyzer::set_params(const Dict & new_params)
{
	params = new_params;
	mask = params["mask"];
	nvec = params["nvec"];
	nimg = params["nimg"];

	// count pixels under mask
	pixels=0;
	size_t totpix=mask->get_xsize()*mask->get_ysize()*mask->get_zsize();
	float *d=mask->get_data();
	for (size_t i=0; i<totpix; ++i) if (d[i]) ++pixels;

	printf("%d,%d\n",pixels,nimg);
	A=gsl_matrix_alloc(pixels,nimg);
	nsofar=0;
}


void EMAN::dump_analyzers()
{
	dump_factory < Analyzer > ();
}

map<string, vector<string> > EMAN::dump_analyzers_list()
{
	return dump_factory_list < Analyzer > ();
}



vector<EMData *> CircularAverageAnalyzer::analyze() {
// 	for (int i=0; i<10; i++)
// 		avg->set_value_at(i,0,i);
	
	if (images.size()!=1) throw ImageDimensionException("Only takes a single image as input");
	int nx=images[0]->get_xsize();
	int ny=images[0]->get_ysize();
	int nz=images[0]->get_zsize();
	if (nz>1)	
		throw ImageDimensionException("Only takes 2D images.");
	int maxr=params.set_default("maxr",nx/2-1);
	int step=params.set_default("step",2);
	
	EMData *avg = new EMData(maxr/step+1,1);

	int ix,iy,it,count;
	for (it=0; it<maxr; it+=step){
		float mn=0;
		count=0;
		for (ix=-maxr-1; ix<=maxr+1; ix++){
			for (iy=-maxr-1; iy<=maxr+1; iy++){
				int d2=ix*ix+iy*iy;
				if (d2>=it*it && d2<(it+step)*(it+step)){
					count++;
					mn+=images[0]->sget_value_at(ix+nx/2,iy+ny/2);
					
				}
			}
		}
		
		mn/=count;
		if(verbose>0) printf("%d,%d,%f\n",it,count,mn);
		avg->set_value_at(it/step,0,mn);
	}
	
	
	ret.push_back(avg);
	return ret;
}



