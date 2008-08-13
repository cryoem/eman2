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
	class FourierPixelInserter3D : public FactoryBase
	{
		public:
		/** Construct a FourierPixelInserter3D
		 */
		FourierPixelInserter3D() : norm(0), rdata(0), nx(0), ny(0), nz(0), nxy(0)
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
		virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight=1.0) = 0;
	

		virtual void init();
		
		/** A function added for testing purposes only
		* Called with the given floating point coordinate [xx,yy,zz] checks to see if 
		* the pixels that the FourierPixelInserter3D's insert_pixel function would
		* have effected (or changed) are near enough to zero.
		 */
		virtual bool effected_pixels_are_zero(const float& xx, const float& yy, const float& zz) = 0;
		
		protected:
			/// A pointer to the constructor argument normalize_values
			float * norm;
			/// A pointer to the constructor argument real_data
			float * rdata;
		
			/// Image volume data sizes, and nxy a convenience variable used here and there
			int nx, ny, nz, nxy;
		
			/// A measure of tolerance used when testing to see if pixels are zero when calling effected_pixels_are_zero
			static float tolerance;
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
		
			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight=1.0);

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
		
			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight=1.0);
		
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
		
			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight=1.0);
				
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
		
			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight=1.0);
		
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
		
			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight=1.0);
		
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
		
			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight=1.0);
		
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
		
			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight=1.0);
		
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
		
			virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight=1.0);
		
			static FourierPixelInserter3D *NEW()
			{
				return new FourierInserter3DMode8();
			}

			virtual string get_name() const
			{
				return "8";
			}
		
			virtual string get_desc() const
			{
				return "Fourier pixel insertion mode 8";
			}
		
			virtual bool effected_pixels_are_zero(const float&, const float&, const float&) { throw; }
			
			virtual void init();
		private:
			int mFreqCutoff;
			float mDFreq;
		// Disallow copy and assignment by default
			FourierInserter3DMode8( const FourierInserter3DMode8& );
			FourierInserter3DMode8& operator=( const FourierInserter3DMode8& );
			
			float* W;
	};
	
	/** QualityScores class is used by the FourierReconstructor and InterpolatedFRC for storing useful quality information.
	 * It's basically a data storage object that has a whole heap of getter and setter methods, and nothing more.
     */
	class QualityScores
	{
		public:
			/** Default constructor
			*/
			QualityScores() : frc_integral(0), snr_normed_frc_intergral(0), normed_snr_integral(0), norm(0) {}
			
			/** Copy constructor
			* @param that the quality scores to be constructed from
			*/
			QualityScores( const QualityScores& that ) : frc_integral(that.frc_integral), 
				snr_normed_frc_intergral(that.snr_normed_frc_intergral), normed_snr_integral(that.normed_snr_integral), norm(that.norm) {}
			
			/** Assignment operator 
			* @param that the quality scores to be constructed from
			*/
			QualityScores& operator=( const QualityScores& that ) 
			{
				if ( &that != this )
				{
					frc_integral = that.frc_integral; 
// 					snr_normed_frc_intergral = that.snr_normed_frc_intergral;
					normed_snr_integral  = that.normed_snr_integral;
					norm = that.norm;
				}
				return *this;
			}

			/** Deconstructor
			* no memory is alloced by this object, but 4 private variables fall out of scope
			*/
			~QualityScores() {}

			/// Various setter and getter methods are below
			
			float get_frc_integral() const { return frc_integral; }
			float get_snr_normed_frc_integral() const { return snr_normed_frc_intergral; }
			float get_normed_snr_integral() const { return normed_snr_integral; }
			float get_norm() const { return norm; }

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

	/** InterpolationFunctoid is an abstract base class, having basically one function which is "operate(float radius)"
	* It simplifies the implementation of InterpolatedFRC::continue_frc_calc? (where ? = 3,4,6 or 7)
	* The other cases (1,2 and 5) must be handled case by case. Note InterpolationFunctoidMode5 is declared 
	* below but it does not derive from InterpolationFunctoid, i.e. it is a special case.
	* Each of these modes corresponds tightly with the FourierInserter3DMode? method of the same number and this
	* is intentional
	*/
	class InterpolationFunctoid
	{
	public:
		InterpolationFunctoid() {}
		virtual ~InterpolationFunctoid() {}
		
		virtual float operate( const float radius ) const = 0;
	};
	
	/** InterpolationFunctoidMode3
	 * see comments for abstract base class InterpolationFunctoid
	*/
	class InterpolationFunctoidMode3 : public InterpolationFunctoid
	{
	public:
		InterpolationFunctoidMode3() {}
		virtual ~InterpolationFunctoidMode3() {}
		
		virtual float operate( const float radius ) const
		{
			return  exp(-radius / EMConsts::I3G);
		}
	};
	
	/** InterpolationFunctoidMode4
	* see comments for abstract base class InterpolationFunctoid
	*/
	class InterpolationFunctoidMode4 : public InterpolationFunctoid
	{
	public:
		InterpolationFunctoidMode4() {}
		virtual ~InterpolationFunctoidMode4() {}
		
		virtual float operate( const float radius ) const
		{
			return  exp(-radius / EMConsts::I4G);
		}
	};
	
	/** InterpolationFunctoidMode5
	* Handles the special case of mode5 interpolation - see FourierInserter3DMode5
	*/
	class InterpolationFunctoidMode5
	{
		public:
			InterpolationFunctoidMode5() { gimx = EMAN::Interp::get_gimx(); }
			virtual ~InterpolationFunctoidMode5()
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

	/** InterpolationFunctoidMode6
	* see comments for abstract base class InterpolationFunctoid
	*/
	class InterpolationFunctoidMode6 : public InterpolationFunctoid
	{
	public:
		InterpolationFunctoidMode6() {}
		virtual ~InterpolationFunctoidMode6() {}
		
		virtual float operate( const float radius ) const
		{
			return  exp(-radius / EMConsts::I5G);
		}
	};
	
	/** InterpolationFunctoidMode7
	* see comments for abstract base class InterpolationFunctoid
	*/
	class InterpolationFunctoidMode7 : public InterpolationFunctoid
	{
	public:
		InterpolationFunctoidMode7() {}
		virtual ~InterpolationFunctoidMode7() {}
		
		virtual float operate( const float radius ) const
		{
			return  EMAN::Interp::hyperg(radius);
		}
	};
	
	
	/** Interpolated FRC - oversees calculation of the FRC and normalization values in Fourier Reconstruction (compares a slice (in some orientation) to a volume)
	* This class works in a similar fashion to FourierPixelInserter3D objects in that the class is first initialized,
	* all of the pixels are "inserted" iteratively, and finally a "finish" is called which calculates the quality scores
	* and returns them in an approprate data object (QualityScores).
	* Typical programmatic usage of this class is as follows
	*
	* InterpatedFRC ifrc(rdata,norm, xsize,ysize,zsize)
	* for pixels in slice
	*		// Figure out x,y and z of each pixel using its Euler angle (InterpatedFRC doesn't do this for you)
	*		// Store the complex pixel value in dt[2]
	*		ifrc.continue_frc_calc2(x,y,z,dt)
	* end
	* QualityScores qualityScores = ifrc.finish(slice->get_attr["ptcl_repr"])
	*
	* Note that the qualityScores objects contains the FRC integral, the SNR weighted FRC integral (which is the new EMAN2 similarity metric), the SNR integral,
	* and the slice normalization score. SNR weighted FRC integral is elaborated upon in the FourierReconstructor comments (reconstructor.h) and will be described and tested
	* online in the EMAN2 wiki (http://blake.bcm.edu/emanwiki/e2make3d). Further, the slice normalization score is the (inverse of the) constant to be applied
	* to the image slice (comprised of all the pixels that were part of the calculation) to achieve normalization with respect to the pixels in the volume it intersects.
	*
	* As you can see from the above example, in order for this approach to work you must initialize the InterpatedFRC
	* object with pointers to 3D volumes containing the true pixel data (rdata) and the associated normalization volume (norm).
	* Note that the norm volume should be half the size of the rdata volume - this is because InterpolatedFRC assumes the real data
	* is complex, and that for each complex pixel there is only one normalization value.
	*
	* Note also that I used continue_frc_calc2 in the above example but could have substituted the "2" with any number
	* from 1 through 7, and in doing so a different technique for approximating the FRC (and normalization) would be used. Each of these 7
	* techniques is meant to be used in conjunction with the associated FourierInserter3DMode? that was used to insert
	* the pixels into the 3D volume in the first place. Because mode 2 is quite often the best method for inserting pixels,
	* it is the mode I used in the example above. Testing showed that, for modes 1 and 2, the FRC scores were in the high 0.8's 
	* to low 0.9s when inserting a noisey image into a 3D volume and then testing it's own FRC. For the other 
	* methods the results weren't as high, varying from 0.7 to 0.5, where mode 7 seemed to give the next best results. Also note that when a 
	* a slice is inserted into a 3D volume and, following this, if the normalization constant is calculated it is usually greater than 1, and this is
	* a consequence of voxels being contributed to by more than one pixel (from the slice) in the 3D volume.
	*/
	class InterpolatedFRC
	{
		public:
			/** Normal constructor
			* @param new_params 
			*/
			InterpolatedFRC(const Dict & new_params );
			
			/** Destructor
			*/
			~InterpolatedFRC()
			{
				free_memory();
			}
	
			void set_params(const Dict & new_params)
			{
		// note but this is really inserting OR individually replacing...
		// the old data will be kept if it is not written over
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
				d.put("z_scale", EMObject::FLOAT_POINTER);
				d.put("y_scale", EMObject::FLOAT_POINTER);
				d.put("x_scale", EMObject::FLOAT_POINTER);
				d.put("sampling", EMObject::FLOAT_POINTER);
				return d;
			}
			
			/** continue_frc_calc1 - function for including an additional pixel in the calculation of the FRC
			* FRC calculated using nearest neighbor. Meant for use in conjunction with FourierInserter3DMode1
			* @param xx the floating point x location of the incoming pixel 
			* @param yy the floating point y location of the incoming pixel
			* @param zz the floating point z location of the incoming pixel
			* @param dt the complex pixel value stored in a float array of length 2
			* @param weight the weight that this pixel has - corresponds to the weight with which the pixel was original inserted into the 3D volume
			*/
			bool continue_frc_calc1(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1.0);
			
			/** continue_frc_calc2 - function for including an additional pixel in the calculation of the FRC
			 * FRC calculated using weighted average of 8 nearest neighbors. Meant for use in conjunction with FourierInserter3DMode2.
			 * See comments in continue_frc_calc1 for parameter details.
			 */
			bool continue_frc_calc2(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1.0);
			
			/** continue_frc_calc3 - function for including an additional pixel in the calculation of the FRC
			 * Meant for use in conjunction with FourierInserter3DMode3.
			 * See comments in continue_frc_calc1 for parameter details.
			 */
			bool continue_frc_calc3(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1.0);
			
			/** continue_frc_calc4 - function for including an additional pixel in the calculation of the FRC
			 * Meant for use in conjunction with FourierInserter3DMode4.
			 * See comments in continue_frc_calc1 for parameter details.
			 */
			bool continue_frc_calc4(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1.0);
			
			/** continue_frc_calc5 - function for including an additional pixel in the calculation of the FRC
			 * Meant for use in conjunction with FourierInserter3DMode5.
			 * See comments in continue_frc_calc1 for parameter details.
			 */
			bool continue_frc_calc5(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1.0);
			
			/** continue_frc_calc6 - function for including an additional pixel in the calculation of the FRC
			 * Meant for use in conjunction with FourierInserter3DMode6.
			 * See comments in continue_frc_calc1 for parameter details.
			 */
			bool continue_frc_calc6(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1.0);
			
			/** continue_frc_calc7 - function for including an additional pixel in the calculation of the FRC
			 * Meant for use in conjunction with FourierInserter3DMode7.
			 * See comments in continue_frc_calc1 for parameter details.
			 */
			bool continue_frc_calc7(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1.0);
			
			unsigned int get_size() { return size; }
	
			float operator[](const unsigned int& idx) { return frc[idx]; }

			QualityScores finish(const unsigned int num_particles);
	
			void reset();
		
		protected:
			mutable Dict params;
			
		private:
			// Disallow copy construction
			InterpolatedFRC( const InterpolatedFRC& that );
			// Disallow assigment
			InterpolatedFRC& operator=( const InterpolatedFRC& that);
			
			/** continue_frc_calc_functoid
		 	* Meant for convenience in continue_frc_calc3, continue_frc_calc4, continue_frc_calc6, and continue_frc_calc7
		 	* See comments in continue_frc_calc1 for parameter details.
			*/
			bool continue_frc_calc_functoid(const float& xx, const float& yy, const float& zz, const float dt[], const InterpolationFunctoid& functoid, const float& weight = 1.0 );
			
			/** free_memory - frees all associated memory
			*/
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
			
			/** Loads all information from parms into private variables
			* Called in the constructor only.
			*/
			void init();
			
			// Pointers to the 3D (complex) data 
			float* threed_rdata, *norm_data;

			// I wish I could make these unsigned but everything else is ints, so these are too.
			int nx, ny, nz, nxy;
			
			// scale factors are stored when dimensions are being stretched or shrunken in the Fourier reconstructor
			float x_scale, y_scale, z_scale;
	
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
	
	// Factory for FourierPixelInserter3D
	template <> Factory < FourierPixelInserter3D >::Factory();
	
} // namespace EMAN


#endif // eman_reconstructor_tools_h__
