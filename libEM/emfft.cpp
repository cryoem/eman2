/**
 * $Id$
 */
#include "emfft.h"

using namespace EMAN;

static int get_rank(int ny, int nz)
{
	int rank = 3;
	if (ny == 1) {
		rank = 1;
	}
	else if (nz == 1) {
		rank = 2;
	}
	return rank;
}

#ifdef FFTW2

list < FftwPlan * >EMfft::fwplans = list < FftwPlan * >();
FftwPlan EMfft::planf_1d = FftwPlan();
FftwPlan EMfft::planr_1d = FftwPlan();

//rfftw_plan FftwPlan::plan_1d = 0;
//rfftwnd_plan FftwPlan::plan_nd = 0;

FftwPlan::FftwPlan()
:	rank(0), direction(UNKNOWN_DIRECTION)
{
	memset(dims, 0, NDIMS * sizeof(int));
}


FftwPlan::FftwPlan(int n, FftwDirection dir, int flags)
:rank(1), direction(dir)
{
	dims[NDIMS - 1] = n;
	for (int i = 0; i < NDIMS - 1; i++) {
		dims[i] = 1;
	}

	if (dir == REAL_TO_COMPLEX) {
		plan_1d = rfftw_create_plan(n, FFTW_REAL_TO_COMPLEX, flags);
	}
	else if (dir == COMPLEX_TO_REAL) {
		plan_1d = rfftw_create_plan(n, FFTW_COMPLEX_TO_REAL, flags);
	}
}


FftwPlan::FftwPlan(int ndim, int nx, int ny, int nz, FftwDirection dir, int flags)
:rank(ndim), direction(dir)
{
	dims[0] = nz;
	dims[1] = ny;
	dims[2] = nx;

	if (dir == REAL_TO_COMPLEX) {
		plan_nd = rfftwnd_create_plan(rank, dims + (NDIMS - rank), FFTW_REAL_TO_COMPLEX, flags);
	}
	else if (dir == COMPLEX_TO_REAL) {
		plan_nd = rfftwnd_create_plan(rank, dims + (NDIMS - rank), FFTW_COMPLEX_TO_REAL, flags);
	}
}



FftwPlan::~FftwPlan()
{
#if 0
	if (plan_nd) {
		rfftwnd_destroy_plan(plan_nd);
		plan_nd = 0;
	}

	if (plan_1d) {
		rfftw_destroy_plan(plan_1d);
		plan_1d = 0;
	}
#endif
}


bool EMAN::operator==(const FftwPlan & p1, const FftwPlan & p2)
{
	if (p1.rank == p2.rank && p1.direction == p2.direction) {
		for (int i = 0; i < FftwPlan::NDIMS; i++) {
			if (p1.dims[i] != p2.dims[i]) {
				return false;
			}
		}
		return true;
	}
	return false;
}


int EMfft::real_to_complex_1d(float *real_data, float *complex_data, int n)
{
	if (planf_1d.get_dim(0) != n) {
		planf_1d = FftwPlan(n, REAL_TO_COMPLEX, FFTW_ESTIMATE);
	}
	rfftw_one(planf_1d.get_plan_1d(), (fftw_real *) real_data, (fftw_real *) complex_data);
	return 0;
}

int EMfft::complex_to_real_1d(float *complex_data, float *real_data, int n)
{
	if (planr_1d.get_dim(0) != n) {
		planr_1d = FftwPlan(n, COMPLEX_TO_REAL, FFTW_ESTIMATE);
	}

	rfftw_one(planr_1d.get_plan_1d(), (fftw_real *) complex_data, (fftw_real *) real_data);
	return 0;
}



int EMfft::real_to_complex_nd(float *real_data, float *complex_data, int nx, int ny, int nz)
{
	// in-place or out-of-place?
	FFTPLACE fftplace = (real_data == complex_data) 
		? FFT_IN_PLACE
		: FFT_OUT_OF_PLACE;
	FftwPlan *fwplan = make_plan(nx, ny, nz, REAL_TO_COMPLEX, fftplace);
	rfftwnd_one_real_to_complex(fwplan->get_plan_nd(), (fftw_real *) real_data,
								(fftw_complex *) complex_data);
	return 0;
}

int EMfft::complex_to_real_nd(float *complex_data, float *real_data, int nx, int ny, int nz)
{
	FftwPlan *fwplan = make_plan(nx, ny, nz, COMPLEX_TO_REAL);
	rfftwnd_one_complex_to_real(fwplan->get_plan_nd(), (fftw_complex *) complex_data,
								(fftw_real *) real_data);
	return 0;
}

FftwPlan *EMfft::make_plan(int nx, int ny, int nz, 
		                   FftwDirection dir,
						   FFTPLACE fftplace)
{
	int rank = get_rank(ny, nz);
	int flags = FFTW_ESTIMATE | FFTW_USE_WISDOM | FFTW_THREADSAFE;
	if (dir == REAL_TO_COMPLEX) {
		if (FFT_OUT_OF_PLACE == fftplace) {
			flags |= FFTW_OUT_OF_PLACE;
		} else {
			flags |= FFTW_IN_PLACE;
		}
	}
	else if (dir == COMPLEX_TO_REAL) {
		// Always in place because the inverse transform
		// destroys the input array.
		flags |= FFTW_IN_PLACE;
	}


	// plan caching seems to fail if out-of-place transform
	// is followed by an in-place transform, since the in-place
	// transform then uses the out-of-place plan (which
	// has the wrong flags)
	list < FftwPlan * >::iterator qi = fwplans.begin();
	for (; qi != fwplans.end(); qi++) {
		if (**qi == FftwPlan(rank, nx, ny, nz, dir, 0)) {
			break;
		}
	}
	if (qi != fwplans.end()) {
		return *qi;
	}

	if (fwplans.size() == MAX_NUM_PLANS) {
		delete fwplans.front();
		fwplans.pop_front();
	}

	FftwPlan *new_plan = new FftwPlan(rank, nx, ny, nz, dir, flags);
	fwplans.push_back(new_plan);
	return new_plan;
}

#endif

#ifdef FFTW3

int EMfft::real_to_complex_1d(float *real_data, float *complex_data, int n)
{
	fftwf_plan plan = fftwf_plan_dft_r2c_1d(n, real_data, (fftwf_complex *) complex_data,
											FFTW_ESTIMATE);
	fftwf_execute(plan);
	fftwf_destroy_plan(plan);
	return 0;
}

int EMfft::complex_to_real_1d(float *complex_data, float *real_data, int n)
{
	fftwf_plan plan = fftwf_plan_dft_c2r_1d(n, (fftwf_complex *) complex_data, real_data,
											FFTW_ESTIMATE);
	fftwf_execute(plan);
	fftwf_destroy_plan(plan);
	return 0;
}

int EMfft::real_to_complex_nd(float *real_data, float *complex_data, int nx, int ny, int nz)
{
	int rank = get_rank(ny, nz);
	int dims[3];
	dims[0] = nz;
	dims[1] = ny;
	dims[2] = nx;

	fftwf_plan plan = fftwf_plan_dft_r2c(rank, dims + (3 - rank), real_data,
										 (fftwf_complex *) complex_data, FFTW_ESTIMATE);

	fftwf_execute(plan);
	fftwf_destroy_plan(plan);
	return 0;
}

int EMfft::complex_to_real_nd(float *complex_data, float *real_data, int nx, int ny, int nz)
{
	int rank = get_rank(ny, nz);
	int dims[3];
	dims[0] = nz;
	dims[1] = ny;
	dims[2] = nx;

	fftwf_plan plan = fftwf_plan_dft_c2r(rank, dims + (3 - rank), (fftwf_complex *) complex_data,
										 real_data, FFTW_ESTIMATE);

	fftwf_execute(plan);
	fftwf_destroy_plan(plan);
	return 0;
}

#endif
