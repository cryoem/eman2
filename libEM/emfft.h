/**
 * $Id$
 */
#ifndef eman_emfft_h__
#define eman_emfft_h__


#ifdef FFTW2

#include <list>
#include <srfftw.h>

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
		rfftwnd_plan plan_nd;
		rfftw_plan plan_1d;
		int flags;
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
		enum FFTPLACE { FFT_OUT_OF_PLACE, FFT_IN_PLACE };
		static FftwPlan *make_plan(int nx, int ny, int nz, FftwDirection dir, FFTPLACE fftplace = FFT_IN_PLACE);

		enum
		{ MAX_NUM_PLANS = 6 };
		static list < FftwPlan * >fwplans;
		static FftwPlan planf_1d;
		static FftwPlan planr_1d;
	};
}
#endif	//FFTW2

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
#endif	//FFTW3

#ifdef NATIVE_FFT
namespace EMAN
{
	/** EMfft converts 1d/nd data from real to complex or from complex to real.
     */
	class EMfft
	{
	  public:
		static int real_to_complex_1d(float *real_data, float *complex_data, int n);
		static int complex_to_real_1d(float *complex_data, float *real_data, int n);

		static int real_to_complex_nd(float *real_data, float *complex_data, int nx, int ny, int nz);
		static int complex_to_real_nd(float *complex_data, float *real_data, int nx, int ny, int nz);
	};
}
#endif	//NATIVE_FFT

#ifdef	ACML
#include <acml.h>
#include <functional>

namespace EMAN
{
	/** EMfft converts 1d/nd data from real to complex or from complex to real.
     */
	class EMfft
	{
	  public:
		static int real_to_complex_1d(float *real_data, float *complex_data, int n);
		static int complex_to_real_1d(float *complex_data, float *real_data, int n);

		static int real_to_complex_nd(float *real_data, float *complex_data, int nx, int ny, int nz);
		static int complex_to_real_nd(float *complex_data, float *real_data, int nx, int ny, int nz);
	  
	  private:
		class time_sqrt_n : public std::unary_function<float, float> {
		  public:
			time_sqrt_n(int n) : n_(n), factor(sqrt(float(n_))) {}
			float operator()(float x) const {return x*factor;}	
		  private:
		    int n_;
		    float factor;
		};
		
		class divide_sqrt_n : public std::unary_function<float, float> {
		  public:
			divide_sqrt_n(int n) : n_(n), factor(sqrt(float(n_))) {}
			float operator()(float x) const {return x/factor;}
		  private:
			int n_;
			float factor;
		};
	};
}	
#endif	//ACML

#endif	//eman_emfft_h__
