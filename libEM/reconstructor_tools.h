/**
 * $Id$
 */

/*
 * Author: David Woolford, 07/25/2007 (woolford@bcm.edu)
 * Copyright (c) 2000-2007 Baylor College of Medicine
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
//#define RECONDEBUG		// This is used to add some data objects for debugging reconstructors, it should normally be commented out!

#ifndef eman_reconstructor_toosl_h__
#define eman_reconstructor_tools_h__ 1

#include "emdata.h"
#include "interp.h"

#include <string>
using std::string;

#include <cstdlib>

//debug
#include <iostream>
using std::cout;
using std::endl;


namespace EMAN
{

	/** FourierPixelInserter3D class defines a way a continuous pixel in 3D
	 * is inserted into the discrete 3D volume - there are various schemes for doing this
	 * including simply finding the nearest neighbor to more elaborate schemes that involve
	 * interpolation using the nearest 8 voxels and so on.
	 *
	 * FourierPixelInserter3D class is the base class for all pixel inserters.
	 * Each specific pixel inserter class has a unique ID name. This name
	 * is used to create a FourierPixelInserter3D instance or do a pixel insertion.
	 * Currently these IDs are (strings) 1,2,3,4,5,6,7 - where mode 2 is currently the
	 * most popularly used pixel insertion mode as it generally performs well, and was the
	 * default mode in EMAN1 - it interpolates to the nearest 8 voxels using a gaussian
	 * weighting based on Euclidian distance
	 *
	 * All FourierPixelInserter3D classes in EMAN2 are managed by a Factory
	 * pattern. So each FourierPixelInserter3D class must define:
	 *   - a unique name to idenfity itself in the factory.
	 *   - a static method to register itself in the factory.
	 *
     * Typical usages of FourierPixelInserter3D are as follows:
	 *
     *  - How to get all the FourierPixelInserter3D names:
	 *@code
	 *    vector<string> all_pixel_inserters = Factory<FourierPixeldtInserter3D>::get_list();
	@endcode
	 *
	 *  - How to use a FourierPixelInserter3D in your code
	 *@code
	 *  // First set up the params correctly this is essential - setting up the params requires
	 *  // these 5 parameters (only and always these 5). nx,ny, and nz are the dims of data and
	 *  // norm, data is the (Fourier) EMData object, and norm is the associated normalization matrix.
	 *	parms["data"] = image;
	 *	parms["norm"] = tmp_data->get_data();
	 *  // The two is the ID of the FourierInserter3DMode2 in this case.
	 *  FourierPixelInserter3D* r = Factory<FourierPixelInserter3D>::get("2", params);
	 *  // Then call init - this causes all of the internal data to be stored correctly
	 *  r->init()
	 *  // Then tell the inserter to insert the pixel with real and imaginary components dt[0] and dt[1] at floating point coords [xx,yy,zz]
	 *  // using the given weight - this is typically done for pixels in a set of many slices - see FourierReconstructor::insert_slice
	 *  r->insert_pixel(xx, yy, zz, dt, weight)
	@endcode
	 */
	class FourierPixelInserter3D : public FactoryBase
	{
		public:
		/** Construct a FourierPixelInserter3D
		 */
		FourierPixelInserter3D() : norm(0), data(0), nx(0), ny(0), nz(0), nxyz(0)
		{}

		/** Desctruct a FourierPixelInserter3D
		 */
		virtual ~FourierPixelInserter3D()
		{
			free_memory();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("data", EMObject::EMDATA);
			d.put("norm", EMObject::FLOAT_POINTER);
			return d;
		}

		/** Insert a complex pixel [dt[0]+dt[1]i] at (float) coordinate [xx,yy,zz] with weighting into a discrete 3D volume
		* @param xx the floating point x coordinate
		* @param yy the floating point y coordinate
		* @param zz the floating point z coordinate
		* @param dt the complex pixel value (dt[0] is real, dt[1] is imaginary)
		* @param weight the weight to given to this complex pixel
		* @return A boolean that indicates the pixel has been inserted (or not)
		 */
		virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const std::complex<float> dt, const float& weight=1.0) = 0;


		virtual void init();

#ifdef RECONDEBUG
		double *ddata;
		double *dnorm;
#endif

		protected:
			/// A pointer to the constructor argument normalize_values
			float * norm;
			/// A pointer to the constructor argument real_data
			EMData * data;

			/// Image volume data sizes a convenience variable used here and there
			int nx, ny, nz,nxyz;
			int nx2,ny2,nz2;
			int subx0,suby0,subz0,fullnx,fullny,fullnz;

		private:
		// Disallow copy and assignment by default
			FourierPixelInserter3D( const FourierPixelInserter3D& );
			FourierPixelInserter3D& operator=( const FourierPixelInserter3D& );

			void free_memory()
			{
			}
	};

	/** FourierPixelInserter3DMode1  - encapsulates "method 1" for inserting a 2D Fourier slice into a 3D volume
	 * See comments in FourierPixelInserter3D for explanations
	 */
	class FourierInserter3DMode1 : public FourierPixelInserter3D
	{
		public:
			FourierInserter3DMode1() {}
			virtual ~FourierInserter3DMode1() {}

			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const std::complex<float> dt, const float& weight=1.0);

			static FourierPixelInserter3D *NEW()
			{
				return new FourierInserter3DMode1();
			}

			virtual string get_name() const
			{
				return NAME;
			}

			virtual string get_desc() const
			{
				return "Fourier pixel insertion using nearest neighbor";
			}

			static const string NAME;

		private:
		// Disallow copy and assignment by default
			FourierInserter3DMode1( const FourierInserter3DMode1& );
			FourierInserter3DMode1& operator=( const FourierInserter3DMode1& );
	};

	/** FourierPixelInserter3DMode2  - encapsulates "method 2" for inserting a 2D Fourier slice into a 3D volume
	 * See comments in FourierPixelInserter3D for explanations
	 */
	class FourierInserter3DMode2 : public FourierPixelInserter3D
	{
		public:
			FourierInserter3DMode2() {}
			virtual ~FourierInserter3DMode2() {}

			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const std::complex<float> dt, const float& weight=1.0);

			static FourierPixelInserter3D *NEW()
			{
				return new FourierInserter3DMode2();
			}

			virtual string get_name() const
			{
				return NAME;
			}

			virtual string get_desc() const
			{
				return "Fourier pixel insertion using interpolation and the nearest 8 voxels";
			}

			static const string NAME;

		// Disallow copy and assignment by default
			FourierInserter3DMode2( const FourierInserter3DMode2& );
			FourierInserter3DMode2& operator=( const FourierInserter3DMode2& );
	};

	/** FourierPixelInserter3DMode3  - encapsulates "method 3" for inserting a 2D Fourier slice into a 3D volume
	 * See comments in FourierPixelInserter3D for explanations
	 */
	class FourierInserter3DMode3 : public FourierPixelInserter3D
	{
		public:
			FourierInserter3DMode3() {}
			virtual ~FourierInserter3DMode3() {}

			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const std::complex<float> dt, const float& weight=1.0);

			static FourierPixelInserter3D *NEW()
			{
				return new FourierInserter3DMode3();
			}

			virtual string get_name() const
			{
				return NAME;
			}

			virtual string get_desc() const
			{
				return "Fourier pixel insertion using a 3x3x3 Gaussian kernel";
			}

			static const string NAME;

		private:
		// Disallow copy and assignment by default
			FourierInserter3DMode3( const FourierInserter3DMode3& );
			FourierInserter3DMode3& operator=( const FourierInserter3DMode3& );
	};

	/** FourierPixelInserter3DMode5  - encapsulates "method 5" for inserting a 2D Fourier slice into a 3D volume
	 * See comments in FourierPixelInserter3D for explanations
	 */
	class FourierInserter3DMode5 : public FourierPixelInserter3D
	{
		public:
			FourierInserter3DMode5()
			{
//				gimx = EMAN::Interp::get_gimx();
			}

			virtual ~FourierInserter3DMode5()
			{
				// Don't delete gimx it causes a seg fault
// 				if ( gimx != 0 )
// 				{
// 					delete gimx;
// 					gimx = 0;
// 				}
			}

			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const std::complex<float> dt, const float& weight=1.0);

			static FourierPixelInserter3D *NEW()
			{
				return new FourierInserter3DMode5();
			}

			virtual string get_name() const
			{
				return NAME;
			}

			virtual string get_desc() const
			{
				return "Fourier pixel insertion mode 5";
			}

			static const string NAME;

		private:
		// Disallow copy and assignment by default
			FourierInserter3DMode5( const FourierInserter3DMode5& );
			FourierInserter3DMode5& operator=( const FourierInserter3DMode5& );

//			float * gimx;
	};

	/** FourierPixelInserter3DMode6  - encapsulates "method 6" for inserting a 2D Fourier slice into a 3D volume
	 * See comments in FourierPixelInserter3D for explanations
	 */
	class FourierInserter3DMode6 : public FourierPixelInserter3D
	{
		public:
			FourierInserter3DMode6() {}
			virtual ~FourierInserter3DMode6() {}

			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const std::complex<float> dt, const float& weight=1.0);

			static FourierPixelInserter3D *NEW()
			{
				return new FourierInserter3DMode6();
			}

			virtual string get_name() const
			{
				return NAME;
			}

			virtual string get_desc() const
			{
				return "More exact version of gauss_5";
			}

			static const string NAME;

		private:
		// Disallow copy and assignment by default
			FourierInserter3DMode6( const FourierInserter3DMode6& );
			FourierInserter3DMode6& operator=( const FourierInserter3DMode6& );
	};

	/** FourierPixelInserter3DMode7  - encapsulates "method 7" for inserting a 2D Fourier slice into a 3D volume
	 * See comments in FourierPixelInserter3D for explanations
	 */
	class FourierInserter3DMode7 : public FourierPixelInserter3D
	{
		public:
			FourierInserter3DMode7() {}
			virtual ~FourierInserter3DMode7() {}

			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const std::complex<float> dt, const float& weight=1.0);

			static FourierPixelInserter3D *NEW()
			{
				return new FourierInserter3DMode7();
			}

			virtual string get_name() const
			{
				return NAME;
			}

			virtual string get_desc() const
			{
				return "Hypergeometric kernel 5x5x5";
			}

			static const string NAME;

		private:
		// Disallow copy and assignment by default
			FourierInserter3DMode7( const FourierInserter3DMode7& );
			FourierInserter3DMode7& operator=( const FourierInserter3DMode7& );
	};

	/** FourierPixelInserter3DMode8  - encapsulates "method 8" for inserting a 2D Fourier slice into a 3D volume
	 * See comments in FourierPixelInserter3D for explanations
	 */
	class FourierInserter3DMode8 : public FourierPixelInserter3D
	{
		public:
			FourierInserter3DMode8() : W(0)
			{

			}
			virtual ~FourierInserter3DMode8()
			{
				if ( W != 0 )
					delete [] W;
			}

			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const std::complex<float> dt, const float& weight=1.0);

			static FourierPixelInserter3D *NEW()
			{
				return new FourierInserter3DMode8();
			}

			virtual string get_name() const
			{
				return NAME;
			}

			virtual string get_desc() const
			{
				return "Fourier pixel insertion mode 8";
			}

			virtual void init();

			static const string NAME;

		private:
			int mFreqCutoff;
			float mDFreq;
		// Disallow copy and assignment by default
			FourierInserter3DMode8( const FourierInserter3DMode8& );
			FourierInserter3DMode8& operator=( const FourierInserter3DMode8& );

			float* W;
	};

	/** FourierPixelInserter3DMode9  - encapsulates "method 9" for inserting a 2D Fourier slice into a 3D volume
	 * See comments in FourierPixelInserter3D for explanations
	 */
	class FourierInserter3DMode9 : public FourierPixelInserter3D
	{
		public:
			FourierInserter3DMode9() {}
			virtual ~FourierInserter3DMode9() {}

			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const std::complex<float> dt, const float& weight=1.0);

			static FourierPixelInserter3D *NEW()
			{
				return new FourierInserter3DMode9();
			}

			virtual string get_name() const
			{
				return NAME;
			}

			virtual string get_desc() const
			{
				return "Kaiser-bessel (KB) kernel 8x8x8";
			}

			static const string NAME;

		private:
		// Disallow copy and assignment by default
			FourierInserter3DMode9( const FourierInserter3DMode9& );
			FourierInserter3DMode9& operator=( const FourierInserter3DMode9& );
	};

	/** FourierPixelInserter3DMode10  - encapsulates "method 10" for inserting a 2D Fourier slice into a 3D volume
	 * See comments in FourierPixelInserter3D for explanations
	 */
	class FourierInserter3DMode10 : public FourierPixelInserter3D
	{
		public:
			FourierInserter3DMode10() {}
			virtual ~FourierInserter3DMode10() {}

			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const std::complex<float> dt, const float& weight=1.0);

			static FourierPixelInserter3D *NEW()
			{
				return new FourierInserter3DMode10();
			}

			virtual string get_name() const
			{
				return NAME;
			}

			virtual string get_desc() const
			{
				return "(imprecise) Kaiser-bessel derived (KBD) kernel 8x8x8";
			}

			static const string NAME;

		private:
		// Disallow copy and assignment by default
			FourierInserter3DMode10( const FourierInserter3DMode10& );
			FourierInserter3DMode10& operator=( const FourierInserter3DMode10& );
	};

	// Factory for FourierPixelInserter3D
	template <> Factory < FourierPixelInserter3D >::Factory();

} // namespace EMAN


#endif // eman_reconstructor_tools_h__
