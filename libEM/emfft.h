/**
 * $Id$
 */
#ifndef eman_emfft_h__
#define eman_emfft_h__


#ifdef FFTW2

#include <list>
#include <rfftw.h>

using std::list;

namespace EMAN
{
	enum FftwDirection
	{
		REAL_TO_COMPLEX = 0,
		COMPLEX_TO_REAL = 1,
		UNKNOWN_DIRECTION
	};

	/** FftwPlan wraps 1d/nd fftw plan operation.
     *
     *  Typical usage:
     *
     *  FftwPlan *plan_1d = new FftwPlan(12, REAL_TO_COMPLEX, FFTW_ESTIMATE);
     *  rfftw_plan plan1 = plan_1d->get_plan_1d();
     * 
     *  FftwPlan *plan_nd = new FftwPlan(rank, nx, ny, nz, REAL_TO_COMPLEX, FFTW_ESTIMATE);
     *  rfftwnd_plan *plan2 = plan_nd->get_plan_nd()
     */
	class FftwPlan
	{
	  public:
		FftwPlan();
		FftwPlan(int n, FftwDirection direction, int flags);
		FftwPlan(int ndim, int nx, int ny, int nz, FftwDirection direction, int flags);
		 ~FftwPlan();

		rfftwnd_plan get_plan_nd() const
		{
			return plan_nd;
		}

		rfftw_plan get_plan_1d() const
		{
			return plan_1d;
		}

		int get_dim(int i) const
		{
			return dims[NDIMS - i - 1];
		}

		friend bool operator==(const FftwPlan & plan1, const FftwPlan & plan2);
	  private:
		enum
		{ NDIMS = 3 };

		int rank;
		int dims[NDIMS];
		FftwDirection direction;
		static rfftwnd_plan plan_nd;
		static rfftw_plan plan_1d;
	};

	bool operator==(const FftwPlan & plan1, const FftwPlan & plan2);

	/** EMfft converts 1d/nd data from real to complex or from complex to real.
     */
	class EMfft
	{
	  public:
		static int real_to_complex_1d(float *real_data, float *complex_data, int n);
		static int complex_to_real_1d(float *complex_data, float *real_data, int n);

		static int real_to_complex_nd(float *real_data, float *complex_data, int nx, int ny,
									  int nz);
		static int complex_to_real_nd(float *complex_data, float *real_data, int nx, int ny,
									  int nz);

	  private:
		static FftwPlan *make_plan(int nx, int ny, int nz, FftwDirection dir);

		enum
		{ MAX_NUM_PLANS = 6 };
		static list < FftwPlan * >fwplans;
		static FftwPlan planf_1d;
		static FftwPlan planr_1d;
	};
}
#endif

#ifdef FFTW3
#include <fftw3.h>

namespace EMAN
{
	/** EMfft converts 1d/nd data from real to complex or from complex to real.
     */
	class EMfft
	{
	  public:
		static int real_to_complex_1d(float *real_data, float *complex_data, int n);
		static int complex_to_real_1d(float *complex_data, float *real_data, int n);

		static int real_to_complex_nd(float *real_data, float *complex_data, int nx, int ny,
									  int nz);
		static int complex_to_real_nd(float *complex_data, float *real_data, int nx, int ny,
									  int nz);
	};
}
#endif

#endif
