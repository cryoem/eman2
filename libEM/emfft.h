#ifndef eman_emfft_h__
#define eman_emfft_h__


#ifdef FFTW2
    
#include <list>
#include <srfftw.h>
    
using std::list;

namespace EMAN
{
    enum FftwDirection {
	REAL_TO_COMPLEX = 0,
	COMPLEX_TO_REAL = 1,
	UNKNOWN_DIRECTION
    };

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
	enum { NDIMS = 3 };

	int rank;
	int dims[NDIMS];
	FftwDirection direction;
	rfftwnd_plan plan_nd;
	rfftw_plan plan_1d;
    };
    
    bool operator==(const FftwPlan & plan1, const FftwPlan & plan2);

    class EMfft
    {
    public:
	static int real_to_complex_1d(float *real_data, float *complex_data, int n);
	static int complex_to_real_1d(float *complex_data, float *real_data, int n);

	static int real_to_complex_nd(float *real_data, float *complex_data, int nx, int ny, int nz);
	static int complex_to_real_nd(float *complex_data, float *real_data, int nx, int ny, int nz);

    private:
	static FftwPlan *make_plan(int nx, int ny, int nz, FftwDirection dir);

	enum { MAX_NUM_PLANS = 6 };
	static list<FftwPlan *> fwplans;
	static FftwPlan planf_1d;
	static FftwPlan planr_1d;
    };
}
#endif

#ifdef FFTW3
#include <fftw3.h>

namespace EMAN
{
    class EMfft
    {
    public:
	static int real_to_complex_1d(float *real_data, float *complex_data, int n);
	static int complex_to_real_1d(float *complex_data, float *real_data, int n);

	static int real_to_complex_nd(float *real_data, float *complex_data, int nx, int ny, int nz);
	static int complex_to_real_nd(float *complex_data, float *real_data, int nx, int ny, int nz);
    };
}
#endif

#endif
