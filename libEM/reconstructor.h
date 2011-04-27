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

#ifndef eman_reconstructor_h__
#define eman_reconstructor_h__ 1
#include <fstream>
#include <boost/shared_ptr.hpp>
#include "emdata.h"
#include "exception.h"
#include "emobject.h"
#include "interp.h"

using std::vector;
using std::map;
using std::string;
using boost::shared_ptr;

using std::cout;
using std::cerr;
using std::endl;

#include <utility>
using std::pair;

#include "reconstructor_tools.h"

namespace EMAN
{

	class Transform;
	class EMData;

	/** Reconstructor class defines a way to do 3D recontruction.
	 * A reconstruction is done by 3 steps:
	 *   - set up.
	 *   - insert all image slices.
	 *   - finish up. The last step will return the result.
	 *
	 * Reconstructor class is the base class for all reconstructors.
     * Each specific Reconstructor class has a unique ID name. This name
     * is used to create a Reconstructor instance or do a reconstruction.
     *
	 * All Reconstructor classes in EMAN are managed by a Factory
	 * pattern. So each Reconstructor class must define:
	 *   - a unique name to idenfity itself in the factory.
	 *   - a static method to register itself in the factory.
	 *
     * Typical usages of Reconstructors are as follows:
     *
     *  - How to get all the Reconstructor names:
     *@code
     *    vector<string> all_reconstructors = Factory<Reconstructor>::get_list();
     @endcode
	 *
     *  - How to use a Reconstructor
     *@code
     *    Reconstructor* r = Factory<Reconstructor>::get("fourier");
     *    r->setup();
     *    r->determine_slice_agreement(slice,euler,weight,true);
     *    ...
     *    r->insert_slice(slice, euler, weight);
     *    ...
     *    EMData* result = r->finish();
     @endcode
	 *
     *  - How to define a new Reconstructor type \n
     *    A new XYZReconstructor class must implement the following functions:
     *    (Please replace 'XYZ' with your own class name).
	 @code
     *        void setup();
     *        int insert_slice(const EMData* const slice, const Transform & t);
     *        EMData * finish();
     *        string get_name() const { return "xyz"; }
     *        static Reconstructor *NEW() { return new XYZReconstructor(); }
     *        TypeDict get_param_types() const;
	 @endcode
	*/
	class Reconstructor : public FactoryBase
	{
	  public:
		Reconstructor() {}
		virtual ~Reconstructor() {}
		/** Initialize the reconstructor.
		 */
		virtual void setup() = 0;
		
		/** Initialize the reconstructor with a seed volume. This can be used to provide some 'default' value
		when there is missing data in Fourier space. The passed 'seed' must be of the appropriate padded size, must be
		in Fourier space, and the same EMData* object will be returned by finish(), meaning the Reconstructor is 
		implicitly taking ownership of the object. However, in Python, this means the seed may be passed in without
		copying, as the same EMData will be coming back out at the end. The seed_weight determines how 'strong' the
		seed volume should be as compared to other inserted slices in Fourier space. Raises an exception if not
		supported by the Reconstructor, or if there is an error with the size. **/
		virtual void setup_seed(EMData* seed,float seed_weight) {throw;}

	  	/** While you can just insert unprocessed slices, if you call preprocess_slice yourself, and insert the returned
		 * slice instead, repeatedly, it can save a fair bit of computation. The default operation just returns a copy
		 of the image, as the preprocessing is reconstructor-specific.
	  	 * @return the processed slice
	  	 * @param slice the slice to be prepocessed
	  	 * @param t transform
	  	 * @exception InvalidValueException when the specified padding value is less than the size of the images
		 */
		virtual EMData* preprocess_slice( const EMData* const slice, const Transform& t = Transform() ) { EMData *ret=slice->copy(); ret->set_attr("reconstruct_preproc",(int)1); return ret; }

		/** Insert an image slice to the reconstructor. To insert multiple
		 * image slices, call this function multiple times.
		 *
		 * @param slice Image slice.
		 * @param euler Euler angle of this image slice.
	  	 * @param weight A weighting factor for this slice, generally the number of particles in a class-average. May be ignored by some reconstructors
		 * @return 0 if OK. 1 if error.
		 */
		virtual int insert_slice(const EMData* const slice, const Transform & euler,const float weight=1.0) {throw;}

		/** Compares a slice to the current reconstruction volume and computes a normalization factor and
		 * quality. Normalization and quality are returned via attributes set in the passed slice. You may freely mix calls
		 * to determine_slice_agreement with calls to insert_slice, but note that determine_slice_agreement can only use information
		 * from slices that have already been inserted.
		 * reconstruct_norm contains the relative normalization factor which should be applied before inserting the slice
	  	 * reconstruct_qual contains a quality factor (larger better) for this slice as compared to the existing reconstruction
	  	 * @param input_slice The EMData slice to be compared
	  	 * @param euler The orientation of the slice as a Transform object
	  	 * @param weight A weighting factor for this slice, generally the number of particles in a class-average. May be ignored by some reconstructors
		 * @param sub Flag indicating whether to subtract the slice from the volume before comparing. May be ignored by some reconstructors
		 * @return 0 if OK. 1 if error.
	  	 * @exception
		 */
		virtual int determine_slice_agreement(EMData* slice, const Transform &euler, const float weight=1.0, bool sub=true ) { throw; }

		/** Finish reconstruction and return the complete model.
		 * @param doift A flag indicating whether the returned object should be guaranteed to be in real-space (true) or should be left in whatever space the reconstructor generated
		 * @return The result 3D model.
		 */
		virtual EMData *finish(bool doift=true) { throw; }
		
		/** set the volume and tmp_volume data to zero, for use in Monte Carlo reconstructors
		*/
		virtual void clear() {throw; }

		/** Print the current parameters to std::out
		 */
		void print_params() const
		{
			std::cout << "Printing reconstructor params" << std::endl;
			for ( Dict::const_iterator it = params.begin(); it != params.end(); ++it )
			{
				std::cout << (it->first) << " " << (it->second).to_str() << std::endl;
			}
			std::cout << "Done printing reconstructor params" << std::endl;
		}


		EMObject& operator[]( const string& key ) { return params[key]; }

	  private:
		// Disallow copy construction
		Reconstructor(const Reconstructor& that);
		Reconstructor& operator=(const Reconstructor& );

	};

	/** This is a Mixin class 
	 *  A class object encapsulating the volume data required by Reconstructors
	 *  It basically stores two (pointers) to EMData objectsd stores the dimensions of the image volume.
	 *  One EMData object basically stores the real pixel data, the other is used for storing normalization values.
	 *  This class was originally added simply to encapsulate the
	 *  the things common to FourierReconstructor, WienerFourierReconstructor
	 *  and BackProjectionReconstructor. It was never expected to instantiated on its own,
	 *  and is intended to be a parent of the Reconstructor class.
	 *  d.woolford May 2007
	 */

	class ReconstructorVolumeData
	{
		public:
			/** Only constructor
			 * All member variables are zeroed
			 */
			inline ReconstructorVolumeData() : image(0), tmp_data(0), nx(0), ny(0), nz(0), subnx(0), subny(0), subnz(0), subx0(0), suby0(0), subz0(0) {}
			
			/** Destructor safely frees memory
			 */
			virtual ~ReconstructorVolumeData() { free_memory(); }

			/** Get the main image pointer, probably redundant (not used)
			 */
			const EMData* get_emdata() { return image; }
		protected:
			//These EMData pointers will most probably be allocated in setup() and released in finish()
			/// Inheriting class allocates this, probably in setup().
			EMData* image;
			/// Inheriting class may allocate this, probably in setup()
			EMData* tmp_data;

			// nx,ny,nz generally will store the dimensions of image
			int nx,nx2;
			int ny,ny2;
			int nz,nz2;
			
			int subnx;
			int subny;
			int subnz;
			
			int subx0;
			int suby0;
			int subz0;

		protected:
			/** Free allocated memorys 
			 * The inherited class may have allocated image of tmp_data
			 * In either case you can safely call this function to delete
			 * either of those pointers, even if they bdb:refine_03#threed_00are NULL
			 */
			void free_memory()
			{
				if (image != 0)  {delete image; image = 0;}
				if ( tmp_data != 0 ) { delete tmp_data; tmp_data = 0; }
			}

			/** Normalize on the assumption that image is a Fourier volume
			 * and that tmp_data is a volume of weights corresponding in size
			 * to this Fourier volume. This means tmp_data is assumed to 
			 * have have as many x pixels as image.
			 */
			virtual void normalize_threed(const bool sqrt_damp=false,const bool wiener=false);

			/** Sends the pixels in tmp_data and image to zero
			 * Convenience only
			 */
			virtual void zero_memory()
			{
				if (tmp_data != 0 ) tmp_data->to_zero();
				if (image != 0 ) image->to_zero();
			}

		private:
		/** Disallow copy construction */
		ReconstructorVolumeData(const ReconstructorVolumeData& that);
		/** Disallow  assignment */
		ReconstructorVolumeData& operator=(const ReconstructorVolumeData& );

	};

	/** This class originally added for 2D experimentation and prototying.
	 * It is basically a replica of the FourierReconstructor, but works in 2D
	 * @author David Woolford and Phil Baldwin
	 * @date early 2008
	 */
	class FourierReconstructorSimple2D : public Reconstructor, public ReconstructorVolumeData
	{
		public:
			FourierReconstructorSimple2D() {}

			virtual ~FourierReconstructorSimple2D() { }

			virtual void setup();

			virtual int insert_slice(const EMData* const slice, const Transform & euler,const float weight=1.0);

			virtual EMData *finish(bool doift=true);

			virtual string get_name() const { return NAME; }

			virtual string get_desc() const { return "performs 2D reconstruction"; }

			static Reconstructor *NEW()
			{
				return new FourierReconstructorSimple2D();
			}


			virtual TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("nx", EMObject::INT, "Necessary. The x dimension of the input images.");
// 				d.put("sym", EMObject::STRING, "Symmetry - assumed to be C1 if not specified");
				return d;
			}
			
			static const string NAME;
	};



	/** Fourier space 3D reconstruction
	 * The Fourier reconstructor is designed to work in an iterative fashion, where similarity ("quality") metrics
	 * are used to determine if a slice should be inserted into the 3D in each subsequent iteration.
	 * The client creates a Fourier reconstructor to insert real images into a 3D volume. The return image is a real space image
	 *
	 *
	 * This reconstructor is based on EMAN1's Fourier reconstructor with a handful of modifications including
	 * 1. - Fourier ring correlation (FRC) as opposed to the mean phase residual is used to estimate slice quality.
	 *		The FRC of the slice in the 3D volume is determined - but the slice is removed from the 3D volume before
	 * 		doing this so the score reflects the extent to which the slice agrees with the contribution of the other
	 *		slices in the 3D volume. The FRC is converted to SNR using the relationship described by Penczek
	 *		( Three-dimensional spectral signal to noise ratio for a class of reconstruction algorithms, JSB 2002 138 (24-46) )
	 *		FRC = S/(sqrt(S+N1)sqrt(S+N2))
	 *		Where N1 is the noise in the slice of the 3D volume and N2 is the noise in the image slice being inserted.
	 *		We make the assumption that the noise in the 3D volume is 0 (N1=0) to get
	 *		FRC^2 = SNR/(1+SNR)
	 *		which gives a spectral SNR plot - we then divide each SNR value by the number of particles in the class average
	 *		(seeing as SNR should scale linearly with the number of particles) to get the estimated SNR per contributing particle
	 * 		in this class average. If the particles that have been averaged are not homogenous this score should be low etc.
	 *		The scaled SNR curve is then converted back to a FRC curve and integrated. This integral is the similarity metric,
	 *		and depends on how far information extends to in Fourier space - typical values range from 0.05 to 0.2, but can vary
	 *		substantially depending on the data.
	 *
	 * 2 - 	Uses half of the memory used by EMAN1's equivalent reconstruction algorithm
	 *
	 *
	 * - Fourier reconstructor usage
	 * @ingroup CUDA_ENABLED
	 *@code
	 *	Reconstructor* r = Factory<Reconstructor>::get("fourier", params);
	 *	r->setup();
	 *	for k in 0:num_iterations-1
	 *		// First do a round of slice quality (metric) determination - only possible if a 3D volume has
	 *		// already been generated (k>0)
	 *		if ( k  > 0 )
	 *			// Determine the agreement of the slices with the previous reconstructed volume (in memory)
	 *			for i in 0:num_slices-1
	 *				r->determine_slice_agreement(image[i], image[i].euler_orientation);
	 *
	 *		// Insert the slices into the 3D volume
	 *		// Will decide not to insert the slice if the its "quality" is not good enough
	 *		for i in 0:num_slices-1
	 *			int failure = r->insert_slice(image[i], image[i].euler_orientation);
	 *			if ( failure ) cout << "Slice was not inserted due to poor quality" << endl;
	 *
	 *	// Get the resulting volume
	 *	EMData* result = r->finish();
	 *	result->write_image("threed.mrc");
	@endcode
	 */
	class FourierReconstructor : public Reconstructor, public ReconstructorVolumeData
	{
	  public:
		/** Default constructor
		* calls load_default_settings()
		*/
		FourierReconstructor() { load_default_settings(); }

		/** Deconstructor
		* calls free_memory()
		*/
		virtual ~FourierReconstructor() { free_memory(); }

		/** Setup the Fourier reconstructor
		* @exception InvalidValueException When one of the input parameters is invalid
		*/
		virtual void setup();

		/** Initialize the reconstructor with a seed volume. This can be used to provide some 'default' value
		when there is missing data in Fourier space. The passed 'seed' must be of the appropriate padded size, must be
		in Fourier space, and the same EMData* object will be returned by finish(), meaning the Reconstructor is 
		implicitly taking ownership of the object. However, in Python, this means the seed may be passed in without
		copying, as the same EMData will be coming back out at the end. The seed_weight determines how 'strong' the
		seed volume should be as compared to other inserted slices in Fourier space.
		* @exception InvalidValueException When one of the input parameters is invalid
		*/
		virtual void setup_seed(EMData* seed,float seed_weight);

	  	/** Preprocess the slice prior to insertion into the 3D volume
		 * this Fourier tranforms the slice and make sure all the pixels are in the right position
		 * it always returns a copy of the provided slice, so it should be deleted by someone eventually.
	  	 * @return the processed slice
	  	 * @param slice the slice to be prepocessed
	  	 * @param t transform to be used for insertion
	  	 * @exception InvalidValueException when the specified padding value is less than the size of the images
		 */
		virtual EMData* preprocess_slice( const EMData* const slice, const Transform& t = Transform() );

		/** Insert a slice into a 3D volume, in a given orientation
		* @return 0 if successful, 1 otherwise
		* @param slice the image slice to be inserted into the 3D volume
		* @param euler Euler angle of this image slice.
		* @param weight A weighting factor for this slice, generally the number of particles in a class-average. May be ignored by some reconstructors
		* @return 0 if OK. 1 if error.
		* @exception NullPointerException if the input EMData pointer is null
		* @exception ImageFormatException if the image is complex as opposed to real
		*/
		virtual int insert_slice(const EMData* const slice, const Transform & euler,const float weight=1.0);


		/** Compares a slice to the current reconstruction volume and computes a normalization factor and
		 * quality. Normalization and quality are returned via attributes set in the passed slice. You may freely mix calls
		 * to determine_slice_agreement with calls to insert_slice, but note that determine_slice_agreement can only use information
		 * from slices that have already been inserted. Attributes set in the slice are:
		 * reconstruct_norm    the relative normalization factor which should be applied before inserting the slice
	  	 * reconstruct_qual    a scaled quality factor (larger better) for this slice as compared to the existing reconstruction
		 * reconstruct_absqual the absolute (not scaled based on weight) quality factor comparing this slice to the existing reconstruction
		 * reconstruct_weight  the summed weights from all voxels intersecting with the inserted slice, larger -> more overlap with other slices
	  	 * @param input_slice The EMData slice to be compared
	  	 * @param euler The orientation of the slice as a Transform object
	  	 * @param weight A weighting factor for this slice, generally the number of particles in a class-average. May be ignored by some reconstructors
		 * @param sub Flag indicating whether to subtract the slice from the volume before comparing. May be ignored by some reconstructors
		 * @return 0 if OK. 1 if error.
		* @exception NullPointerException if the input EMData pointer is null
		* @exception ImageFormatException if the image is complex as opposed to real
		*/
		virtual int determine_slice_agreement(EMData* slice, const Transform &euler, const float weight=1.0, bool sub=true );

		/** Get the reconstructed volume
		* Normally will return the volume in real-space with the requested size. The calling application is responsible for 
		* removing any padding.
		* @param doift A flag indicating whether the returned object should be guaranteed to be in real-space (true) or should be left in whatever space the reconstructor generated
		* @return The real space reconstructed volume
		*/
		virtual EMData *finish(bool doift=true);
		
		/** clear the volume and tmp_data for use in Monte Carlo reconstructions
		*/
		virtual void clear();

		/** Get the unique name of the reconstructor
		*/
		virtual string get_name() const
		{
			return NAME;
		}

		/** Get the one line description of the reconstructor
		*/
		virtual string get_desc() const
		{
			return "Reconstruction via direct Fourier methods using one of a variety of different kernels, most of which are Gaussian based";
		}

		/** Factory incorporation uses the pointer of this function
		* @return a Reconstructor pointer to a newly allocated FourierReconstructor
		*/
		static Reconstructor *NEW()
		{
			return new FourierReconstructor();
		}

		/** Get the parameter types of this object
		* @return a TypeDict detailing all of the acceptable (and necessary) parameters
		*/
		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("size", EMObject::INTARRAY, "Required. The dimensions of the real-space output volume, including any padding (must be handled by the calling application). Assumed that apix x/y/z identical.");
			d.put("sym", EMObject::STRING, "Optional. The symmetry of the reconstructed volume, c?, d?, oct, tet, icos, h?. Default is c1, ie - an asymmetric object");
			d.put("mode", EMObject::STRING, "Optional. Fourier pixel insertion mode name (nearest_neighbor, gauss_2, gauss_3, gauss_5, gauss_5_slow, gypergeom_5, experimental) gauss_2 is the default.");
			d.put("sqrtnorm", EMObject::BOOL, "Optional. When normalizing, additionally divides by the sqrt of the normalization factor to damp exaggerated features. Is this justifyable ? No idea (yet). Default is false.");
			d.put("verbose", EMObject::BOOL, "Optional. Toggles writing useful information to standard out. Default is false.");
			d.put("quiet", EMObject::BOOL, "Optional. If false, print verbose information.");
			d.put("subvolume",EMObject::INTARRAY, "Optional. (xorigin,yorigin,zorigin,xsize,ysize,zsize) all in Fourier pixels. Useful for parallelism.");
			d.put("savenorm",EMObject::STRING, "Debug. Will cause the normalization volume to be written directly to the specified file when finish() is called.");
			return d;
		}
		
		static const string NAME;

	  protected:
		/** Load default settings
		*/
		virtual void load_default_settings();

		/** Frees the memory owned by this object (but not parent objects)
		 * Deletes the FourierPixelInserter3D pointer
		 */
		virtual void free_memory();

		/** Load the pixel inserter based on the information in params
		 */
		virtual void load_inserter();

		/** A function to perform the nuts and bolts of inserting an image slice
		 * @param input_slice the slice to insert into the 3D volume
		 * @param euler a transform storing the slice euler angle
		 * @param weight weighting factor for this slice (usually number of particles in a class-average)
		 */
		virtual void do_insert_slice_work(const EMData* const input_slice, const Transform & euler,const float weight);

		/** A function to perform the nuts and bolts of comparing an image slice
		 * @param input_slice the slice to insert into the 3D volume
		 * @param euler a transform storing the slice euler angle
		 */
		virtual void do_compare_slice_work(EMData* input_slice, const Transform & euler,float weight);

		/** This is a mode-2 pixel extractor
		 * @param xx,yy,zz voxel coordinates (need not be integers)
		 * @param dt float pointer with 3 floats allocated for returned complex value and weight sum
		 */
		virtual bool pixel_at(const float& xx, const float& yy, const float& zz, float *dt);

		/// A pixel inserter pointer which inserts pixels into the 3D volume using one of a variety of insertion methods
		FourierPixelInserter3D* inserter;

	  private:
		 /** Disallow copy construction
  		 */
  		FourierReconstructor( const FourierReconstructor& that );
  		/**Disallow assignment
  		 */
  		FourierReconstructor& operator=( const FourierReconstructor& );

	};



	/** Fourier space 3D reconstruction
	 * This is a modified version of the normal FourierReconstructor which is aware of the SSNR information stored
	 * in individual class-average headers as "ctf_snr_total" and "ctf_wiener_filtered". It will perform a reconstruction
	 * with a nonisotropic Wiener filter applied to the final reconstruction, and will 'undo' the Wiener filter on the
	 * individual class-averages if ctf_wiener_filtered is set. This represents something which was not possible to accomplish
	 * in EMAN1, and should produce superior results, with proper anisotropic filtering. Still, the filtration makes the assumption
	 * that the original SNR estimates were accurate, and that the data will average completely coherently, which is not
	 * truly the case. This may produce models which are somewhat underfiltered in the Wiener sense, but since B-factor
	 * corrections are not applied in the ctf.auto averager, this effect is likely already more than compensated for.
	*/
	class WienerFourierReconstructor : public FourierReconstructor
	{
	  public:
		/** Default constructor
		* calls load_default_settings()
		*/
		WienerFourierReconstructor() {};

		/** Deconstructor
		* calls free_memory()
		*/
		virtual ~WienerFourierReconstructor() { }


		/** Insert a slice into a 3D volume, in a given orientation
		* @return 0 if successful, 1 otherwise
		* @param slice the image slice to be inserted into the 3D volume
		* @param euler Euler angle of this image slice.
		* @param weight This is ignored in this reconstructor, since the SSNR from the particle header is used instead
		* @return 0 if OK. 1 if error.
		* @exception NullPointerException if the input EMData pointer is null
		* @exception ImageFormatException if the image is complex as opposed to real
		*/
		virtual int insert_slice(const EMData* const slice, const Transform & euler,const float weight=1.0);


		/** Compares a slice to the current reconstruction volume and computes a normalization factor and
		 * quality. Normalization and quality are returned via attributes set in the passed slice. You may freely mix calls
		 * to determine_slice_agreement with calls to insert_slice, but note that determine_slice_agreement can only use information
		 * from slices that have already been inserted. Attributes set in the slice are:
		 * reconstruct_norm    the relative normalization factor which should be applied before inserting the slice
	  	 * reconstruct_qual    a scaled quality factor (larger better) for this slice as compared to the existing reconstruction
		 * reconstruct_absqual the absolute (not scaled based on weight) quality factor comparing this slice to the existing reconstruction
		 * reconstruct_weight  the summed weights from all voxels intersecting with the inserted slice, larger -> more overlap with other slices
	  	 * @param input_slice The EMData slice to be compared
	  	 * @param euler The orientation of the slice as a Transform object
	  	 * @param weight This is ignored except for it's sign, since the SSNR from the particle header is used instead
		 * @param sub Flag indicating whether to subtract the slice from the volume before comparing. May be ignored by some reconstructors
		 * @return 0 if OK. 1 if error.
		* @exception NullPointerException if the input EMData pointer is null
		* @exception ImageFormatException if the image is complex as opposed to real
		*/
		virtual int determine_slice_agreement(EMData* slice, const Transform &euler, const float weight=1.0, bool sub=true );

		/** Get the reconstructed volume
		* Normally will return the volume in real-space with the requested size. The calling application is responsible for 
		* removing any padding.
		* @param doift A flag indicating whether the returned object should be guaranteed to be in real-space (true) or should be left in whatever space the reconstructor generated
		* @return The real space reconstructed volume
		*/
		virtual EMData *finish(bool doift=true);

		/** Get the unique name of the reconstructor
		*/
		virtual string get_name() const
		{
			return NAME;
		}

		/** Get the one line description of the reconstructor
		*/
		virtual string get_desc() const
		{
			return "Reconstruction via direct Fourier methods using one of a variety of different kernels, most of which are Gaussian based. This version also incorporates a nonisotropic Wiener filter based on SNR estimates stored in the class-average headers by the ctf.auto averager.";
		}

		/** Factory incorporation uses the pointer of this function
		* @return a Reconstructor pointer to a newly allocated WienerFourierReconstructor
		*/
		static Reconstructor *NEW()
		{
			return new WienerFourierReconstructor();
		}
		
		static const string NAME;

	  protected:
		
		virtual void do_insert_slice_work(const EMData* const input_slice, const Transform & euler,const float weight);

		/** A function to perform the nuts and bolts of comparing an image slice
		 * @param input_slice the slice to insert into the 3D volume
		 * @param euler a transform storing the slice euler angle
		 */
		virtual void do_compare_slice_work(EMData* input_slice, const Transform & euler,float weight);

		/** This is a mode-2 pixel extractor
		 * @param xx,yy,zz voxel coordinates (need not be integers)
		 * @param dt float pointer with 3 floats allocated for returned complex value and weight sum
		 */
		virtual bool pixel_at(const float& xx, const float& yy, const float& zz, float *dt);

		/// A pixel inserter pointer which inserts pixels into the 3D volume using one of a variety of insertion methods
//		FourierPixelInserter3D* inserter;

	  private:
		 /** Disallow copy construction
  		 */
  		WienerFourierReconstructor( const WienerFourierReconstructor& that );
  		/**Disallow assignment
  		 */
  		WienerFourierReconstructor& operator=( const WienerFourierReconstructor& );

	};

	/** Real space 3D reconstruction using back projection.
     *
     * Back-projection is a method of 3D reconstruction from 2D
     * projections. It is based on superposing 3D functions
     * ("back-projection bodies") obtained by translating the
     * 2D projections along the directions of projection.
     */
	class BackProjectionReconstructor:public Reconstructor, public ReconstructorVolumeData
	{
	  public:
		BackProjectionReconstructor() { load_default_settings();  }

		virtual ~BackProjectionReconstructor() {}

		virtual void setup();

		/** Insert an image slice to the reconstructor. To insert multiple
		 * image slices, call this function multiple times.
		 *
		 * @param slice Image slice.
		 * @param euler Euler angle of this image slice.
	  	 * @param weight A weighting factor for this slice, generally the number of particles in a class-average. May be ignored by some reconstructors
		 * @return 0 if OK. 1 if error.
		 */
		virtual int insert_slice(const EMData* const slice, const Transform & euler,const float weight=1.0);

		virtual EMData *finish(bool doift=true);

		virtual string get_name() const
		{
			return NAME;
		}

		virtual string get_desc() const
		{
			return "Simple (unfiltered) back-projection reconstruction. Weighting by contributing particles in the class average is optional and default behaviour";
		}

		static Reconstructor *NEW()
		{
			return new BackProjectionReconstructor();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("size", EMObject::INT, "Necessary. The x and y dimensions of the input images.");
			d.put("weight", EMObject::FLOAT, "Optional. A temporary value set prior to slice insertion, indicative of the inserted slice's weight. Default sis 1.");
			d.put("sym", EMObject::STRING, "Optional. The symmetry to impose on the final reconstruction. Default is c1");
			d.put("zsample", EMObject::INT, "Optional. The z dimensions of the reconstructed volume.");
			return d;
		}
		
		static const string NAME;
		
	  private:
		// Disallow copy construction
		BackProjectionReconstructor( const BackProjectionReconstructor& that);
		// Disallow assignment
		BackProjectionReconstructor& operator=( const BackProjectionReconstructor& );

		void load_default_settings()
		{
			params["weight"] = 1.0;
			params["use_weights"] = true;
			params["size"] = 0;
			params["sym"] = "c1";
			params["zsample"] = 0;
		}

		EMData* preprocess_slice(const EMData* const slice, const Transform& t);
	};


	/** Direct Fourier inversion Reconstructor
     *
     */
	EMData* padfft_slice( const EMData* const slice, const Transform& t, int npad );

	class nn4Reconstructor:public Reconstructor
	{
	  public:
		nn4Reconstructor();

		nn4Reconstructor( const string& symmetry, int size, int npad );

		virtual ~nn4Reconstructor();

		virtual void setup();

		/** Insert an image slice to the reconstructor. To insert multiple
		 * image slices, call this function multiple times.
		 *
		 * @param slice Image slice.
		 * @param euler Euler angle of this image slice.
	  	 * @param weight A weighting factor for this slice, generally the number of particles in a class-average. May be ignored by some reconstructors
		 * @return 0 if OK. 1 if error.
		 */
		virtual int insert_slice(const EMData* const slice, const Transform & euler,const float weight=1.0);

		virtual EMData *finish(bool doift=true);

		virtual string get_name() const
		{
			return NAME;
		}

		virtual string get_desc() const
		{
			return "Direct Fourier inversion routine";
		}

		static Reconstructor *NEW()
		{
			return new nn4Reconstructor();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("size",		EMObject::INT);
			d.put("npad",		EMObject::INT);
			d.put("sign",		EMObject::INT);
			d.put("ndim",		EMObject::INT);
			d.put("snr",		EMObject::FLOAT);
			d.put("symmetry",	EMObject::STRING);
			d.put("snr",		EMObject::FLOAT);
			d.put("fftvol",		EMObject::EMDATA);
			d.put("weight",		EMObject::EMDATA);
			d.put("weighting",      EMObject::INT);
			return d;
		}

		void setup( const string& symmetry, int size, int npad );

		int insert_padfft_slice( EMData* padded, const Transform& trans, int mult=1 );

		static const string NAME;

	  private:
		EMData* m_volume;
		EMData* m_wptr;
		string  m_symmetry;
		int m_weighting;
		int m_vnx, m_vny, m_vnz;
		int m_npad;
		int m_nsym;
		int m_ndim;
		int m_vnzp, m_vnyp, m_vnxp;
		int m_vnzc, m_vnyc, m_vnxc;
		void buildFFTVolume();
		void buildNormVolume();
		float m_wghta;
		float m_wghtb;
		float m_osnr;
		void load_default_settings()
		{
			//params["use_weights"] = false;
		}
	};


	/** Direct Fourier inversion Reconstructor for extremly rectangular object
     *
     */
	

	class nn4_rectReconstructor:public Reconstructor
	{
	  public:
		nn4_rectReconstructor();

		nn4_rectReconstructor( const string& symmetry, int size, int npad );

		virtual ~nn4_rectReconstructor();

		virtual void setup();

		/** Insert an image slice to the reconstructor. To insert multiple
		 * image slices, call this function multiple times.
		 *
		 * @param slice Image slice.
		 * @param euler Euler angle of this image slice.
	  	 * @param weight A weighting factor for this slice, generally the number of particles in a class-average. May be ignored by some reconstructors
		 * @return 0 if OK. 1 if error.
		 */
		virtual int insert_slice(const EMData* const slice, const Transform & euler,const float weight=1.0);

		virtual EMData *finish(bool doift=true);

		virtual string get_name() const
		{
			return NAME;
		}

		virtual string get_desc() const
		{
			return "Direct Fourier inversion routine";
		}

		static Reconstructor *NEW()
		{
			return new nn4_rectReconstructor();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("sizeprojection", EMObject::INT);
			d.put("sizex",		EMObject::INT);
			d.put("sizey",		EMObject::INT);
			d.put("sizez",		EMObject::INT);
			d.put("xratio",		EMObject::FLOAT);
			d.put("yratio", 	EMObject::FLOAT);
			d.put("npad",		EMObject::INT);
			d.put("sign",		EMObject::INT);
			d.put("ndim",		EMObject::INT);
			d.put("snr",		EMObject::FLOAT);
			d.put("symmetry",	EMObject::STRING);
			d.put("snr",		EMObject::FLOAT);
			d.put("fftvol",		EMObject::EMDATA);
			d.put("weight",		EMObject::EMDATA);
			d.put("weighting",      EMObject::INT);
			return d;
		}

		void setup( const string& symmetry, int size, int npad );

		int insert_padfft_slice( EMData* padded, const Transform& trans, int mult=1 );

		static const string NAME;

	  private:
		EMData* m_volume;
		EMData* m_wptr;
		string  m_symmetry;
		int m_weighting;
		int m_vnx, m_vny, m_vnz;
		int m_npad;
		int m_nsym;
		int m_ndim;
		int m_vnzp, m_vnyp, m_vnxp;
		int m_vnzc, m_vnyc, m_vnxc;
		int m_count;
		float m_xratio,m_yratio,m_zratio;//ratio of x,y,z direction in the 3d volume comparing to the cubic case
		float m_xscale,m_yscale;//ratior of x,y direction of 2D FFT after scaling and roatating operations
		int m_sizeofprojection;
		void buildFFTVolume();
		void buildNormVolume();
		float m_wghta;
		float m_wghtb;
		float m_osnr;
		void load_default_settings()
		{
			//params["use_weights"] = false;
		}
	};


     /* Fourier Reconstruction by nearest neighbor with 3D SSNR
        Added by Zhengfan Yang on 03/16/07
     */

	class nnSSNR_Reconstructor:public Reconstructor
	{

	  public:
		nnSSNR_Reconstructor();

		nnSSNR_Reconstructor( const string& symmetry, int size, int npad);

		~nnSSNR_Reconstructor();

		virtual void setup();

		/** Insert an image slice to the reconstructor. To insert multiple
		 * image slices, call this function multiple times.
		 *
		 * @param slice Image slice.
		 * @param euler Euler angle of this image slice.
	  	 * @param weight A weighting factor for this slice, generally the number of particles in a class-average. May be ignored by some reconstructors
		 * @return 0 if OK. 1 if error.
		 */
		virtual int insert_slice(const EMData* const slice, const Transform & euler,const float weight=1.0);

		virtual EMData *finish(bool doift=true);

		virtual string get_name() const
		{
			return NAME;
		}

		virtual string get_desc() const
		{
			return "Reconstruction by nearest neighbor with 3D SSNR";
		}

		static Reconstructor *NEW()
		{
			return new nnSSNR_Reconstructor();
		}

		virtual	TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("size", EMObject::INT);
			d.put("npad", EMObject::INT);
			d.put("symmetry", EMObject::STRING);
			d.put("fftvol", EMObject::EMDATA);
			d.put("weight", EMObject::EMDATA);
			d.put("weight2", EMObject::EMDATA);
			d.put("SSNR", EMObject::EMDATA);
			d.put("vol_ssnr", EMObject::EMDATA);
			d.put("w", EMObject::FLOAT);
			return d;
		}

		void setup( const string& symmetry, int size, int npad);

		int insert_padfft_slice( EMData* padded, const Transform& trans, int mult=1 );

		static const string NAME;
		
	  private:
		EMData* m_volume;
		EMData* m_wptr;
		EMData* m_wptr2;
		string  m_symmetry;
		int m_weighting;
		int m_vnx, m_vny, m_vnz;
		int m_npad;
		int m_nsym;
		int m_vnzp, m_vnyp, m_vnxp;
		int m_vnzc, m_vnyc, m_vnxc;
		void buildFFTVolume();
		void buildNormVolume();
		void buildNorm2Volume();
		float m_wghta;
		float m_wghtb;
	};


	/** nn4_ctf Direct Fourier Inversion Reconstructor
     *
     */
	class nn4_ctfReconstructor:public Reconstructor
	{
	  public:
		nn4_ctfReconstructor();

		nn4_ctfReconstructor( const string& symmetry, int size, int npad, float snr, int sign );

		virtual ~nn4_ctfReconstructor();

		virtual void setup();

		/** Insert a slice into a 3D volume, in a given orientation
		* @return 0 if successful, 1 otherwise
		* @param slice the image slice to be inserted into the 3D volume
		* @param euler Euler angle of this image slice.
		* @param weight A weighting factor for this slice, generally the number of particles in a class-average. May be ignored by some reconstructors
		* @return 0 if OK. 1 if error.
		* @exception NullPointerException if the input EMData pointer is null
		* @exception ImageFormatException if the image is complex as opposed to real
		*/
		virtual int insert_slice(const EMData* const slice, const Transform & euler,const float weight=1.0);

		virtual EMData *finish(bool doift=true);

		virtual string get_name() const
		{
			return NAME;
		}

		virtual string get_desc() const
		{
			return "Direct Fourier inversion reconstruction routine";
		}

		static Reconstructor *NEW()
		{
			return new nn4_ctfReconstructor();
		}


		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("size",		EMObject::INT);
			d.put("npad",		EMObject::INT);
			d.put("sign",		EMObject::INT);
			d.put("symmetry",	EMObject::STRING);
			d.put("snr",		EMObject::FLOAT);
			d.put("fftvol",		EMObject::EMDATA);
			d.put("weight",		EMObject::EMDATA);
			d.put("weighting",      EMObject::INT);
			d.put("varsnr",         EMObject::INT);
			return d;
		}

		void setup( const string& symmetry, int size, int npad, float snr, int sign );

		int insert_padfft_slice( EMData* padfft, const Transform& trans, int mult=1);

		int insert_buffed_slice( const EMData* buffer, int mult );
		
		static const string NAME;
		
	  private:
		EMData* m_volume;
		EMData* m_wptr;
		int m_vnx, m_vny, m_vnz;
		int m_vnzp, m_vnyp, m_vnxp;
		int m_vnxc, m_vnyc, m_vnzc;
		int m_npad;
		int m_sign;
		int m_varsnr;
		int m_weighting;
		float m_wghta, m_wghtb;
		float m_snr;
		string m_symmetry;
		int m_nsym;

		void buildFFTVolume();
		void buildNormVolume();

	};


	/** nn4_ctf_rectDirect Fourier Inversion Reconstructor
     *
     */
	class nn4_ctf_rectReconstructor:public Reconstructor
	{
	  public:
		nn4_ctf_rectReconstructor();

		nn4_ctf_rectReconstructor( const string& symmetry, int size, int npad, float snr, int sign );

		virtual ~nn4_ctf_rectReconstructor();

		virtual void setup();

		/** Insert a slice into a 3D volume, in a given orientation
		* @return 0 if successful, 1 otherwise
		* @param slice the image slice to be inserted into the 3D volume
		* @param euler Euler angle of this image slice.
		* @param weight A weighting factor for this slice, generally the number of particles in a class-average. May be ignored by some reconstructors
		* @return 0 if OK. 1 if error.
		* @exception NullPointerException if the input EMData pointer is null
		* @exception ImageFormatException if the image is complex as opposed to real
		*/
		virtual int insert_slice(const EMData* const slice, const Transform & euler, const float weight=1.0);

		virtual EMData *finish(bool doift=true);

		virtual string get_name() const
		{
			return NAME;
		}

		virtual string get_desc() const
		{
			return "Direct Fourier inversion reconstruction routine";
		}

		static Reconstructor *NEW()
		{
			return new nn4_ctf_rectReconstructor();
		}


		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("sizeprojection", EMObject::INT);
			d.put("sizex",		EMObject::INT);
			d.put("sizey",		EMObject::INT);
			d.put("sizez",		EMObject::INT);
			d.put("xratio",		EMObject::FLOAT);
			d.put("yratio", 	EMObject::FLOAT);
			d.put("size",		EMObject::INT);
			d.put("npad",		EMObject::INT);
			d.put("sign",		EMObject::INT);
			d.put("symmetry",	EMObject::STRING);
			d.put("snr",		EMObject::FLOAT);
			d.put("fftvol",		EMObject::EMDATA);
			d.put("weight",		EMObject::EMDATA);
			d.put("weighting",  EMObject::INT);
			d.put("varsnr",     EMObject::INT);
			return d;
		}

		void setup( const string& symmetry, int size, int npad, float snr, int sign );

		int insert_padfft_slice( EMData* padfft, const Transform& trans, int mult=1);

		int insert_buffed_slice( const EMData* buffer, int mult );
		
		static const string NAME;
		
	  private:
		EMData* m_volume;
		EMData* m_wptr;
		int m_vnx, m_vny, m_vnz;
		int m_vnzp, m_vnyp, m_vnxp;
		int m_vnxc, m_vnyc, m_vnzc;
		int m_count;
		float m_xratio, m_yratio, m_zratio;//ratio of x,y,z direction in the 3d volume comparing to the cubic case
		float m_xscale, m_yscale;//ratior of x,y direction of 2D FFT after scaling and roatating operations
		int m_sizeofprojection;
		int m_npad;
		int m_sign;
	        int m_varsnr;
		int m_weighting;
		float m_wghta, m_wghtb;
		float m_snr;
		string m_symmetry;
		int m_nsym;

		void buildFFTVolume();
		void buildNormVolume();
	};



     /* Fourier Reconstruction by nearest neighbor with 3D SSNR and CTF
        Added by Zhengfan Yang on 04/11/07
     */

	class nnSSNR_ctfReconstructor:public Reconstructor
	{

	  public:
		nnSSNR_ctfReconstructor();

		nnSSNR_ctfReconstructor( const string& symmetry, int size, int npad, float snr, int sign);

		~nnSSNR_ctfReconstructor();

		virtual void setup();

		/** Insert a slice into a 3D volume, in a given orientation
		* @return 0 if successful, 1 otherwise
		* @param slice the image slice to be inserted into the 3D volume
		* @param euler Euler angle of this image slice.
		* @param weight A weighting factor for this slice, generally the number of particles in a class-average. May be ignored by some reconstructors
		* @return 0 if OK. 1 if error.
		* @exception NullPointerException if the input EMData pointer is null
		* @exception ImageFormatException if the image is complex as opposed to real
		*/
		virtual int insert_slice(const EMData* const slice, const Transform & euler,const float weight=1.0);


		virtual EMData *finish(bool doift=true);

		virtual string get_name() const
		{
			return NAME;
		}

		virtual string get_desc() const
		{
			return "Reconstruction by nearest neighbor with 3D SSNR with CTF";
		}

		static Reconstructor *NEW()
		{
			return new nnSSNR_ctfReconstructor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("size",     EMObject::INT);
			d.put("npad",     EMObject::INT);
			d.put("symmetry", EMObject::STRING);
			d.put("fftvol",   EMObject::EMDATA);
			d.put("fftwvol",  EMObject::EMDATA);
			d.put("weight",   EMObject::EMDATA);
			d.put("weight2",  EMObject::EMDATA);
			d.put("weight3",  EMObject::EMDATA);
			d.put("SSNR",     EMObject::EMDATA);
			d.put("vol_ssnr", EMObject::EMDATA);
			d.put("w",        EMObject::FLOAT);
			d.put("sign",     EMObject::INT);
			d.put("snr",      EMObject::FLOAT);
			return d;
		}
		void setup( const string& symmetry, int size, int npad, float snr, int sign);

		int insert_padfft_slice( EMData* padded, const Transform& trans, int mult=1 );

		static const string NAME;     
                
	  private:
		EMData* m_volume;
		EMData* m_wptr;
		EMData* m_wptr2;
		EMData* m_wptr3;
	        string  m_symmetry;
		int m_weighting;
		int m_vnx, m_vny, m_vnz;
		int m_npad;
		int m_nsym;
		int m_vnzp, m_vnyp, m_vnxp;
		int m_vnzc, m_vnyc, m_vnxc;
		void buildFFTVolume();
		void buildNormVolume();
		void buildNorm2Volume();
		void buildNorm3Volume();
		float m_wghta;
		float m_wghtb;
		int   m_sign;
		float m_snr;
		int wiener;
	};

	template <> Factory < Reconstructor >::Factory();

	void dump_reconstructors();
	map<string, vector<string> > dump_reconstructors_list();


    struct point_t
    {
        int pos2;
        float real;
        float imag;
        float ctf2;
    };


        class newfile_store
        {
        public:
		newfile_store( const string& prefix, int npad, bool ctf );

		virtual ~newfile_store();

		void add_image( EMData* data, const Transform& tf );

		void add_tovol( EMData* fftvol, EMData* wgtvol, const vector<int>& mults, int pbegin, int pend );

		void get_image( int id, EMData* buf );

        void read( int nprj );

		void restart( );

	private:
		int m_npad;

		bool m_ctf;

		string m_bin_file;
		string m_txt_file;

		shared_ptr<std::ofstream> m_bin_of;
		shared_ptr<std::ofstream> m_txt_of;
		shared_ptr<std::ifstream> m_bin_if;
		vector< std::istream::off_type > m_offsets;

        vector< point_t > m_points;
	};

	class file_store
	{
	  public:
		file_store(const string& filename, int npad, int write, bool CTF);

		virtual ~file_store();

		void add_image(EMData* data, const Transform& tf);

		void get_image(int id, EMData* padfft);

		void restart();
	  private:
		shared_ptr<std::ifstream> m_ihandle;
		shared_ptr<std::ofstream> m_bin_ohandle;
		shared_ptr<std::ofstream> m_txt_ohandle;
		string m_bin_file;
		string m_txt_file;
		int m_ctf;
		int m_npad;
		int m_prev;
		int m_x_out;
		int m_y_out;
		int m_z_out;
		int m_write;
		std::istream::off_type m_totsize;
		float m_Cs;
		float m_pixel;
		float m_voltage;
		float m_ctf_applied;
		float m_amp_contrast;
		vector< float > m_defocuses;
		vector< float > m_phis;
		vector< float > m_thetas;
		vector< float > m_psis;
	};

}

#endif

/* vim: set ts=4 noet: */
