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
		private:
#ifdef  FFTW_PLAN_CACHING
#define EMFFTW2_ND_CACHE_SIZE 20
#define EMFFTW2_1D_CACHE_SIZE 20	
		static const int EMAN2_REAL_2_COMPLEX;
		static const int EMAN2_COMPLEX_2_REAL;
		static const int EMAN2_FFTW2_INPLACE;
		static const int EMAN2_FFTW2_OUT_OF_PLACE;
		class EMfftw2_cache_nd
		{
		public:
			EMfftw2_cache_nd();
			~EMfftw2_cache_nd();

			rfftwnd_plan get_plan(const int rank, const int x, const int y, const int z, const int r2c_flag, int ip_flag);
		private:
			// Prints useful debug information to standard out
			void debug_plans();
		
			// Store the number of plans currently contained in this EMfftw3_cache object
			int num_plans;

			// Store the rank of the plan
			int rank[EMFFTW2_ND_CACHE_SIZE];
			// Store the dimensions of the plan (always in 3D, if dimensions are "unused" they are taken to be 1)
			int plan_dims[EMFFTW2_ND_CACHE_SIZE][3];
			// store whether the plan was real to complex or vice versa
			int r2c[EMFFTW2_ND_CACHE_SIZE];
			// Store whether or not the plan was inplace
			int ip[EMFFTW2_ND_CACHE_SIZE];
			// Store the plans themselves
			rfftwnd_plan rfftwnd_plans[EMFFTW2_ND_CACHE_SIZE];
		};
		
		class EMfftw2_cache_1d
		{
		public:
			EMfftw2_cache_1d();
			~EMfftw2_cache_1d();

			rfftw_plan get_plan(const int x, const int r2c_flag, const int ip_flag);
		private:
			// Prints useful debug information to standard out
			void debug_plans();
	
			// Store the number of plans currently contained in this EMfftw3_cache object
			int num_plans;

			// Store the dimensions of the plan (always in 3D, if dimensions are "unused" they are taken to be 1)
			int plan_dims[EMFFTW2_1D_CACHE_SIZE];
			// store whether the plan was real to complex or vice versa
			int r2c[EMFFTW2_1D_CACHE_SIZE];
			// store whether the plan was real to complex or vice versa
			int ip[EMFFTW2_1D_CACHE_SIZE];
			// Store the plans themselves
			rfftw_plan rfftw1d_plans[EMFFTW2_1D_CACHE_SIZE];
		};
		
		static EMfftw2_cache_nd plan_nd_cache;
		static EMfftw2_cache_1d plan_1d_cache;
#endif 
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
	  private:
#ifdef FFTW_PLAN_CACHING
#define EMFFTW3_CACHE_SIZE 20
		static const int EMAN2_REAL_2_COMPLEX;
		static const int EMAN2_COMPLEX_2_REAL;
		static const int EMAN2_FFTW2_INPLACE;
		static const int EMAN2_FFTW2_OUT_OF_PLACE;
		class EMfftw3_cache
		{
		public:
			EMfftw3_cache();
			~EMfftw3_cache();

			fftwf_plan get_plan(const int rank, const int x, const int y, const int z, const int r2c_flag,const int ip_flag, fftwf_complex* complex_data, float* real_data);
		private:
			// Prints useful debug information to standard out
			void debug_plans();
			
			// Store the number of plans currently contained in this EMfftw3_cache object
			int num_plans;

			// Store the rank of the plan
			int rank[EMFFTW3_CACHE_SIZE];
			// Store the dimensions of the plan (always in 3D, if dimensions are "unused" they are taken to be 1)
			int plan_dims[EMFFTW3_CACHE_SIZE][3];
			// store whether the plan was real to complex or vice versa
			int r2c[EMFFTW3_CACHE_SIZE];
			// Store the plans themselves
			fftwf_plan fftwplans[EMFFTW3_CACHE_SIZE];
			// Store whether or not the plan was inplace
			int ip[EMFFTW3_CACHE_SIZE];
		
		};

		static EMfftw3_cache plan_cache;
#endif 
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
