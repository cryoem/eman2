/**
 * $Id$
 */
#ifndef eman_emfft_h__
#define eman_emfft_h__


#ifdef FFTW2

#include <srfftw.h>
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
		static int real_to_complex_2d(float *real_data, float *complex_data, int nx, int ny);
		static int complex_to_real_2d(float *complex_data, float *real_data, int nx, int ny);
		static int real_to_complex_3d(float *real_data, float *complex_data, int nx, int ny, int nz);
		static int complex_to_real_3d(float *complex_data, float *real_data, int nx, int ny, int nz);
		
		class time_sqrt_n : public std::unary_function<float, float> {
		  public:
			time_sqrt_n(int n) : n_(n), factor(sqrt(float(n_))) {}
			float operator()(float x) const {return x*factor;}	
		  private:
		    const int n_;
		    const float factor;
		};
		
		class divide_sqrt_n : public std::unary_function<float, float> {
		  public:
			divide_sqrt_n(int n) : n_(n), factor(sqrt(float(n_))) {}
			float operator()(float x) const {return x/factor;}
		  private:
			const int n_;
			const float factor;
		};
	};
}	
#endif	//ACML

#endif	//eman_emfft_h__
