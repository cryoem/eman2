#include "mpi.h"
#include "stdlib.h"

// include EMAN2
#include "emdata.h"
#include "emassert.h"
#include "projector.h"

#include "ali3d_unified_mpi.h"
#include "project3d.h"
#include "utilcomm.h"
#include "fgcalc.h"

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




