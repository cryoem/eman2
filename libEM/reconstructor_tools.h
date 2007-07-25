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
	/** PixelOperation is an abstract class that PixelAddOperation and PixelMinusOperation derive from
	 * These classes are used as members FourierPixelInserter3D to enable the insertion operation
	 * to switch from subtraction to addition
	 */
	struct PixelOperation
	{
		virtual ~PixelOperation() {}

		/** Perform an operation on target memory using a given value
		 *@param target The memory location the data to operate on
		 *@param value The value to use as the basis of the operation
		 */
		virtual void operate( float* const target, const float& value ) = 0;

		/** Operator() is an alternative way to call operate
		 *@param target The memory location the data to operate on
		 *@param value The value to use as the basis of the operation
		 *@return The value stored in the target memory after the operation has been performed
		 */
		float& operator()(float* const target, const float& value) { operate(target,value); return *target; }
	};

	/** PixelAddOperation - encapsulates the the concept of +=
	 * see PixelOperation comments for details.
	 */
	struct PixelAddOperation : public PixelOperation
	{
		virtual ~PixelAddOperation() {};

		virtual void operate( float* const target, const float& value ) { *target += value; }
	};

	/** PixelAddOperation - encapsulates the the concept of -=
	 * see PixelOperation comments for details.
	 */
	struct PixelMinusOperation : public PixelOperation
	{
		virtual ~PixelMinusOperation() {};

		virtual void operate( float* const target, const float& value ) { *target -= value; }
	};

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
	 *    vector<string> all_pixel_inserters = Factory<FourierPixelInserter3D>::get_list();
	@endcode
	 *
	 *  - How to use a FourierPixelInserter3D in your code
	 *@code
	 *  // First set up the params correctly this is essential - setting up the params requires
	 *  // these 5 parameters (only and always these 5). nx,ny, and nz are the dims of rdata and
	 *  // norm, rdata is the real (Fourier) data, and norm is the associated normalization matrix.
	 *	parms["rdata"] = image->get_data();
	 *	parms["norm"] = tmp_data->get_data();
	 *	parms["nx"] = nx;
	 *	parms["ny"] = ny;
	 *	parms["nz"] = nz;
	 *  // The two is the ID of the FourierInserter3DMode2 in this case.
	 *  FourierPixelInserter3D* r = Factory<FourierPixelInserter3D>::get("2", params);
	 *  // Then call init - this causes all of the internal data to be stored correctly
	 *  r->init()
	 *  // Then tell the inserter to insert the pixel with real and imaginary components dt[0] and dt[1] at floating point coords [xx,yy,zz]
	 *  // using the given weight - this is typically done for pixels in a set of many slices - see FourierReconstructor::insert_slice
	 *  r->insert_pixel(xx, yy, zz, dt, weight)
	@endcode
	 */
	class FourierPixelInserter3D
	{
		public:
		/** Construct a FourierPixelInserter3D
		 */
			FourierPixelInserter3D() : norm(0), rdata(0), nx(0), ny(0), nz(0), nxy(0)
			{
				pixel_operation = new PixelAddOperation;
				other_pixel_operation = new PixelMinusOperation;
			}
		
		/** Desctruct a FourierPixelInserter3D
		 */
			virtual ~FourierPixelInserter3D()
			{
				free_memory();
			}
			
		/** Set the Reconstructor's parameters using a key/value dictionary.
			 * @param new_params A dictionary containing the new parameters.
		 */
			void set_params(const Dict & new_params)
			{
			// note but this is really inserting OR individually replacing...
			// the old data will be kept if it is not written over
			// This is different from the original design but has not yet been confirmed
			// as OK with Steve Ludtke
			// This shouldn't present any problems.
				TypeDict permissable_params = get_param_types();
				for ( Dict::const_iterator it = new_params.begin(); it != new_params.end(); ++it )
				{
				
					if ( !permissable_params.find_type(it->first) )
					{
						throw InvalidParameterException(it->first);
					}
					params[it->first] = it->second;
				}	
			}
		
			TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("nx", EMObject::INT);
				d.put("ny", EMObject::INT);
				d.put("nz", EMObject::INT);
				d.put("rdata", EMObject::FLOAT_POINTER);
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
			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight) = 0;
		
		/** set_pixel_minus_operation makes the default pixel insertion operation to mimic -=
			 * This function is useful when a pixel (more specifically all of the pixels in a slice) needs to be removed from a volume
			 * while maintaining the other contents of the volume. Originally added for testing purposes.
		 */	
			void set_pixel_minus_operation()
			{
				free_memory();
				pixel_operation = new PixelMinusOperation;
				other_pixel_operation = new PixelAddOperation;
			}

		/** set_pixel_add_operation makes the default pixel insertion operation to mimic +=
			 * The default behaviour of this class is to always add the incoming pixel to the volume,
			 * however sometimes it might be necessary to switch between adding and subtracting pixels
			 * from the volume, so if the client ever calls set_pixel_minus_operation they can switch
			 * back to additive behavior by calling this function
		 */	
			void set_pixel_add_operation()
			{
				free_memory();
				pixel_operation = new PixelAddOperation;
				other_pixel_operation = new PixelMinusOperation;
			}

		
		/** Get the FourierPixelInserter3D's name. Each FourierPixelInserter3D is
			 * identified by a unique name, them being (strings) 1,2,3,4,5,6 and 7
			 * @return The FourierPixelInserter3D's ID.
		 */
			virtual string get_name() const = 0;

		/** Get the FourierPixelInserter3D's desc.
			 * @return The FourierPixelInserter3D's description.
		 */
			virtual string get_desc() const = 0;

			virtual void init();
		
		/** A function added for testing purposes only
			 * Called with the given floating point coordinate [xx,yy,zz] checks to see if 
			 * the pixels that the FourierPixelInserter3D's insert_pixel function would
			 * have effected (or changed) are near enough to zero.
		 */
		
			virtual bool effected_pixels_are_zero(const float& xx, const float& yy, const float& zz) = 0;
		protected:
			mutable Dict params;
			/// A pointer to the constructor argument normalize_values
			float * norm;
			/// A pointer to the constructor argument real_data
			float * rdata;
		
			/// Image volume data sizes, and nxy a convenience variable used here and there
			int nx, ny, nz, nxy;

			/// A pixel operation, which is addition by default, but which can be changed to subtraction etc.
			PixelOperation* pixel_operation;

			/// Minus pixel operation is a hack needed to for inserter modes 5,6 and 7.
			PixelOperation* other_pixel_operation;
		
			/// A measure of tolerance used when testing to see if pixels are zero when calling effected_pixels_are_zero
			static float tolerance;
		private:
		// Disallow copy and assignment by default
			FourierPixelInserter3D( const FourierPixelInserter3D& );
			FourierPixelInserter3D& operator=( const FourierPixelInserter3D& );
		
			void free_memory()
			{
				if ( pixel_operation != 0 )
				{
					delete pixel_operation;
					pixel_operation = 0;
				}
				if ( other_pixel_operation != 0 )
				{
					delete other_pixel_operation;
					other_pixel_operation = 0;
				}
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
		
			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& = 1);

			static FourierPixelInserter3D *NEW()
			{
				return new FourierInserter3DMode1();
			}

			virtual string get_name() const
			{
				return "1";
			}
		
			virtual string get_desc() const
			{
				return "Fourier pixel insertion using nearest neighbor";
			}
		
			virtual bool effected_pixels_are_zero(const float& xx, const float& yy, const float& zz);

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
		
			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1);
		
			static FourierPixelInserter3D *NEW()
			{
				return new FourierInserter3DMode2();
			}

			virtual string get_name() const
			{
				return "2";
			}
		
			virtual string get_desc() const
			{
				return "Fourier pixel insertion using interpolation and the nearest 8 voxels";
			}
		
			virtual bool effected_pixels_are_zero(const float& xx, const float& yy, const float& zz);
		
			virtual void init();
		private:
			int off[8];
			float g[8];
		
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
		
			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1);
				
			static FourierPixelInserter3D *NEW()
			{
				return new FourierInserter3DMode3();
			}

			virtual string get_name() const
			{
				return "3";
			}
		
			virtual string get_desc() const
			{
				return "Fourier pixel insertion mode 3";
			}

			virtual bool effected_pixels_are_zero(const float& xx, const float& yy, const float& zz);
		
		private:
		// Disallow copy and assignment by default
			FourierInserter3DMode3( const FourierInserter3DMode3& );
			FourierInserter3DMode3& operator=( const FourierInserter3DMode3& );
	};
	
	/** FourierPixelInserter3DMode4  - encapsulates "method 4" for inserting a 2D Fourier slice into a 3D volume
	 * See comments in FourierPixelInserter3D for explanations
	 */
	class FourierInserter3DMode4 : public FourierPixelInserter3D
	{
		public:
			FourierInserter3DMode4() {}
			virtual ~FourierInserter3DMode4() {}
		
			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1);
		
			static FourierPixelInserter3D *NEW()
			{
				return new FourierInserter3DMode4();
			}

			virtual string get_name() const
			{
				return "4";
			}
		
			virtual string get_desc() const
			{
				return "Fourier pixel insertion mode 4";
			}

			virtual bool effected_pixels_are_zero(const float& xx, const float& yy, const float& zz);
		
		private:
		// Disallow copy and assignment by default
			FourierInserter3DMode4( const FourierInserter3DMode4& );
			FourierInserter3DMode4& operator=( const FourierInserter3DMode4& );
	};
	
	/** FourierPixelInserter3DMode5  - encapsulates "method 5" for inserting a 2D Fourier slice into a 3D volume
	 * See comments in FourierPixelInserter3D for explanations
	 */
	class FourierInserter3DMode5 : public FourierPixelInserter3D
	{
		public:
			FourierInserter3DMode5()
			{
				gimx = EMAN::Interp::get_gimx();
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
		
			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1);
		
			static FourierPixelInserter3D *NEW()
			{
				return new FourierInserter3DMode5();
			}

			virtual string get_name() const
			{
				return "5";
			}
		
			virtual string get_desc() const
			{
				return "Fourier pixel insertion mode 5";
			}
		
			virtual bool effected_pixels_are_zero(const float& xx, const float& yy, const float& zz);
		
		private:
		// Disallow copy and assignment by default
			FourierInserter3DMode5( const FourierInserter3DMode5& );
			FourierInserter3DMode5& operator=( const FourierInserter3DMode5& );
		
			float * gimx;
	};
	
	/** FourierPixelInserter3DMode6  - encapsulates "method 6" for inserting a 2D Fourier slice into a 3D volume
	 * See comments in FourierPixelInserter3D for explanations
	 */
	class FourierInserter3DMode6 : public FourierPixelInserter3D
	{
		public:
			FourierInserter3DMode6() {}
			virtual ~FourierInserter3DMode6() {}
		
			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1);
		
			static FourierPixelInserter3D *NEW()
			{
				return new FourierInserter3DMode6();
			}

			virtual string get_name() const
			{
				return "6";
			}
		
			virtual string get_desc() const
			{
				return "Fourier pixel insertion mode 6";
			}
		
			virtual bool effected_pixels_are_zero(const float& xx, const float& yy, const float& zz);
		
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
		
			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1);
		
			static FourierPixelInserter3D *NEW()
			{
				return new FourierInserter3DMode7();
			}

			virtual string get_name() const
			{
				return "7";
			}
		
			virtual string get_desc() const
			{
				return "Fourier pixel insertion mode 7";
			}
		
			virtual bool effected_pixels_are_zero(const float& xx, const float& yy, const float& zz);
		
		private:
		// Disallow copy and assignment by default
			FourierInserter3DMode7( const FourierInserter3DMode7& );
			FourierInserter3DMode7& operator=( const FourierInserter3DMode7& );
	};
	
	class QualityScores
	{
		public:
			QualityScores() : frc_integral(0), snr_normed_frc_intergral(0), normed_snr_integral(0), norm(0) {}
			QualityScores( const QualityScores& that ) : frc_integral(that.frc_integral), 
				snr_normed_frc_intergral(that.snr_normed_frc_intergral), normed_snr_integral(that.normed_snr_integral), norm(that.norm) {}
			QualityScores& operator=( const QualityScores& that ) 
			{
				frc_integral = that.frc_integral; 
				snr_normed_frc_intergral = that.snr_normed_frc_intergral;
				normed_snr_integral  = that.normed_snr_integral;
				norm = that.norm;
				return *this;
			}

			~QualityScores() {}

			float get_frc_integral() { return frc_integral; }
			float get_snr_normed_frc_integral() { return snr_normed_frc_intergral; }
			float get_normed_snr_integral() { return normed_snr_integral; }
			float get_norm() { return norm; }

			void set_frc_integral( const float& score ) { frc_integral = score; }
			void set_snr_normed_frc_integral(const float& score) { snr_normed_frc_intergral = score; }
			void set_normed_snr_integral(const float& score) { normed_snr_integral = score; }
			void set_norm( const float& score ) { norm = score; }

			void debug_print()
			{
				cout << "frc " << frc_integral << " nfrc " << snr_normed_frc_intergral << " nsnr " << normed_snr_integral << " norm " << norm << endl;
			}
		private:

			float frc_integral, snr_normed_frc_intergral, normed_snr_integral, norm;
		
	};

	class InterpolationFunctiod
	{
	public:
		InterpolationFunctiod() {}
		virtual ~InterpolationFunctiod() {}
		
		virtual float operate( const float radius ) const = 0;
	};
	
	class InterpolationFunctiodMode3 : public InterpolationFunctiod
	{
	public:
		InterpolationFunctiodMode3() {}
		virtual ~InterpolationFunctiodMode3() {}
		
		virtual float operate( const float radius ) const
		{
			return  exp(-radius / EMConsts::I3G);
		}
	};
	
	class InterpolationFunctiodMode4 : public InterpolationFunctiod
	{
	public:
		InterpolationFunctiodMode4() {}
		virtual ~InterpolationFunctiodMode4() {}
		
		virtual float operate( const float radius ) const
		{
			return  exp(-radius / EMConsts::I4G);
		}
	};
	
	class InterpolationFunctiodMode5
	{
		public:
			InterpolationFunctiodMode5() { gimx = EMAN::Interp::get_gimx(); }
			virtual ~InterpolationFunctiodMode5()
			{
				if ( gimx != 0 )
				{
// 					delete gimx;
// 					gimx = 0;
				}
			}
		
			virtual float operate( const int mmx, const int mmy, const int mmz ) const
			{
				return gimx[abs(mmx) + abs(mmy) * 100 + abs(mmz) * 10000];
			}
		private:
			float * gimx;
	};

	
	class InterpolationFunctiodMode6 : public InterpolationFunctiod
	{
	public:
		InterpolationFunctiodMode6() {}
		virtual ~InterpolationFunctiodMode6() {}
		
		virtual float operate( const float radius ) const
		{
			return  exp(-radius / EMConsts::I5G);
		}
	};
	
	class InterpolationFunctiodMode7 : public InterpolationFunctiod
	{
	public:
		InterpolationFunctiodMode7() {}
		virtual ~InterpolationFunctiodMode7() {}
		
		virtual float operate( const float radius ) const
		{
			return  EMAN::Interp::hyperg(radius);
		}
	};
	
	
	
	class InterpolatedFRC
	{
		public:
			InterpolatedFRC() : threed_rdata(0), frc(0), frc_norm_rdata(0), frc_norm_dt(0), size(0), pixel_radius_max(0), r(0), rn(0) {}

			InterpolatedFRC(float* const rdata, const int xsize, const int ysize, const int zsize, const float& sampling=1.0 );
			~InterpolatedFRC()
			{
				free_memory();
			}

		// Copy and assignment
			InterpolatedFRC( const InterpolatedFRC& that );
			InterpolatedFRC& operator=( const InterpolatedFRC& that);
	
			bool continue_frc_calc1(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1.0);
			bool continue_frc_calc2(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1.0);
			bool continue_frc_calc3(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1.0);
			bool continue_frc_calc4(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1.0);
			bool continue_frc_calc5(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1.0);
			bool continue_frc_calc6(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1.0);
			bool continue_frc_calc7(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1.0);
	
			bool continue_frc_calc_functoid(const float& xx, const float& yy, const float& zz, const float dt[], const InterpolationFunctiod& functoid, const float& weight = 1.0 );
			
			unsigned int get_size() { return size; }
	
			float operator[](const unsigned int& idx) { return frc[idx]; }

			QualityScores finish(const unsigned int  num_particles);
	
			void reset();
		private:
			void free_memory()
			{
				if ( frc != 0 )
				{
					delete [] frc;
					frc = 0;
				}
				if ( frc_norm_rdata != 0 )
				{
					delete [] frc_norm_rdata;
					frc_norm_rdata = 0;
				}
				if ( frc_norm_dt != 0 )
				{
					delete [] frc_norm_dt;
					frc_norm_dt = 0;
				}
			}
			// Pointers to the 3D (complex) data 
			float* threed_rdata;

			// I wish I could make these unsigned but everything else is ints, so these are too.
			int nx, ny, nz, nxy;
	
			float bin;

			float* frc;
			float* frc_norm_rdata;
			float* frc_norm_dt;
	
			// The maximum dimension (radius) of a considered pixel
			int size;
			int pixel_radius_max;
			int pixel_radius_max_square;
			
			
			// values used for calculating the normalization value
			float r, rn;
			
			int off[8];
			float g[8];
	};
	
	template <> Factory < FourierPixelInserter3D >::Factory();
	
} // namespace EMAN


#endif // eman_reconstructor_tools_h__
