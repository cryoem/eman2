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

	class Transform3D;
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
     *    r->insert_slice(slice1, euler1);
     *    insert more
     *    EMData* result = r->finish();
     @endcode
	 *
     *  - How to define a new Reconstructor type \n
     *    A new XYZReconstructor class must implement the following functions:
     *    (Please replace 'XYZ' with your own class name).
	 @code
     *        void setup();
     *        int insert_slice(const EMData* const slice, const Transform3D & t);
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

		/** Insert an image slice to the reconstructor. To insert multiple
		 * image slices, call this function multiple times.
		 *
		 * @param slice Image slice.
		 * @param euler Euler angle of this image slice.
		 * @return 0 if OK. 1 if error.
		 */
		virtual int insert_slice(const EMData* const slice, const Transform3D & euler) {throw;};

		// This is pure virtual yet, but it should be... changed by Wei
		virtual int insert_slice(const EMData* const slice, const Transform & euler) {throw;};

		/**
	  	 * @return
	  	 * @param input_slice
	  	 * @param arg
	  	 * @param num_particles_in_slice
	  	 * @exception
		 */
		virtual int determine_slice_agreement(const EMData* const, const Transform3D &, const unsigned int)
		{
			cout << "You called determine slice agreement but nothing happened - there is no functionality for determing slice agreement using this " << get_name() << " reconstructor" << endl;
			return 0;
		}

		virtual int determine_slice_agreement(const EMData* const, const Transform &, const unsigned int)
		{
			throw;
		}

		/** Finish reconstruction and return the complete model.
		 * @return The result 3D model.
		 */
		virtual EMData *finish() = 0;

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

		// A function used only by the Fourier reconstructor for testing and writing to std out purposes in e2make3d.py
		virtual float get_score(const unsigned int) { return 0; }
		// A function used only by the Fourier reconstructor for testing and writing to std out purposes in e2make3d.py
		virtual float get_norm(const unsigned int) { return 0; }


		virtual int insert_slice_weights(const Transform&) {return 0;}

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
			inline ReconstructorVolumeData() : image(0), tmp_data(0), nx(0), ny(0), nz(0) {}
			
			/** Destructor safely frees memory
			 */
			virtual ~ReconstructorVolumeData() { free_memory(); }

			/** Get the main image pointer, probably redundant (not used)
			 */
			const EMData* const get_emdata() { return image; }
		protected:
			//These EMData pointers will most probably be allocated in setup() and released in finish()
			/// Inheriting class allocates this, probably in setup().
			EMData* image;
			/// Inheriting class may allocate this, probably in setup()
			EMData* tmp_data;

			// nx,ny,nz generally will store the dimensions of image
			int nx;
			int ny;
			int nz;

		protected:
			/** Free allocated memorys 
			 * The inherited class may have allocated image of tmp_data
			 * In either case you can safely call this function to delete
			 * either of those pointers, even if they are NULL
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
			virtual void normalize_threed();

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

			virtual int insert_slice(const EMData* const slice, const Transform3D & euler) {
				throw UnexpectedBehaviorException("Transform3D interface is redundant");
			}

			virtual int insert_slice(const EMData* const slice, const Transform & euler);

			virtual EMData *finish();

			virtual string get_name() const { return "fouriersimple2D"; }

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
		FourierReconstructor() : image_idx(0), inserter(0), interp_FRC_calculator(0), slice_insertion_flag(true),
		slice_agreement_flag(false), x_scale_factor(0.0), y_scale_factor(0.0), z_scale_factor(0.0),
		max_input_dim(0), output_x(0), output_y(0), output_z(0) { load_default_settings(); }

		/** Deconstructor
		* calls free_memory()
		*/
		virtual ~FourierReconstructor() { free_memory(); }

		/** Setup the Fourier reconstructor
		* relies on the correct parameters
		* throws a variety of exceptions if the parameters are unworkable
		* @exception InvalidValueException when the x_in parameter is odd or less than zero
		* @exception InvalidValueException when the y_in parameter is odd or less than zero
		* @exception InvalidValueException when the x_out parameter is odd or less than zero
		* @exception InvalidValueException when the y_out parameter is odd or less than zero
		* @exception InvalidValueException when the z_out parameter is odd or less than zero
		* @exception InvalidValueException when the pad parameter is less than the greatest dimension of the input images (the greater of x_in and y_in)
		* Note the restraint that the Fourier volumes should be even will be lifted once debugging of the the xform.fourierorigin processor is performed,
		* but however, when this happens a rigorous testing should be performed because no testing has been done in this regard.
		*/
		virtual void setup();

		/** Insert a slice into a 3D volume, in a given orientation
		* @return 0 if successful, 1 otherwise
		* @param slice the image slice to be inserted into the 3D volume
		* @param t3d the Transform that stores the image Euler angles
		* @exception NullPointerException if the input EMData pointer is null
		* @exception ImageFormatException if the image is complex as opposed to real
		*/
		virtual int insert_slice(const EMData* const slice, const Transform & t3d = Transform());


		/** Determine slice agreement with the current reconstruction
		* @return 0 if successful, 1 otherwise
		* @param input_slice the image slice to compared against the 3D volume
		* @param t3d the Transform that the image Euler angles
		* @param num_particles_in_slice the number of particles in the slice - used to determine the SNR normalized FSC, defaults to 1.
		* @exception NullPointerException if the input EMData pointer is null
		* @exception ImageFormatException if the image is complex as opposed to real
		*/
		virtual int determine_slice_agreement(const EMData* const input_slice, const Transform & t3d, const unsigned int  num_particles_in_slice = 1);


		virtual int insert_slice(const EMData* const slice, const Transform3D & t3d) {
			throw UnexpectedBehaviorException("No support for Transform3D in insert_slice, use Transform instead");
		};

		virtual int determine_slice_agreement(const EMData* const input_slice, const Transform3D & t3d, const unsigned int  num_particles_in_slice = 1) {
			throw UnexpectedBehaviorException("No support for Transform3D in determine_slice_agreement, use Transform instead");
		};

		/** Get the reconstructed volume
		* Peforms Fourier inversion on a potentially large volume and may take
		* several minutes.
		* @return The real space reconstructed volume
		*/
		virtual	EMData *finish();

		/** Get the unique name of the reconstructor
		*/
		virtual string get_name() const
		{
			return "fourier";
		}

		/** Get the one line description of the reconstructor
		*/
		virtual string get_desc() const
		{
			return "Reconstruction via direct Fourier methods using a combination of nearest neighbour and Gaussian kernels";
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
			d.put("x_in", EMObject::INT, "Necessary. The x dimension of the input images.");
			d.put("y_in", EMObject::INT, "Necessary. The y dimension of the input images.");
			d.put("zsample", EMObject::INT, "Optional. The z dimension (Fourier sampling) of the reconstructed volume, very useful for tomographic reconstruction. Works for general volumes.");
			d.put("ysample", EMObject::INT, "Optional. The y dimension (Fourier sampling) of the reconstructed volume, works for general volumes. Not commonly specified.");
			d.put("xsample", EMObject::INT, "Optional. The x dimension (Fourier sampling) of the reconstructed volume, works for general volumes. Not commonly specified.");
			d.put("mode", EMObject::INT, "Optional. Fourier pixel insertion mode [1-8] - mode 2 is default.");
			d.put("hard", EMObject::FLOAT, "Optional. The quality metric threshold. Default is 0 (off).");
			d.put("sym", EMObject::STRING, "Optional. The symmetry of the reconstructed volume, c?, d?, oct, tet, icos, h?. Default is c1");
			d.put("quiet", EMObject::BOOL, "Optional. Toggles writing useful information to standard out. Default is false.");
			d.put("3damp", EMObject::BOOL, "Optional. Toggles writing the 3D FFT amplitude image. Default is false.");
			d.put("weight", EMObject::FLOAT, "Weight of the slice that is being inserted. Default is 1.0.");
			d.put("start_model", EMObject::EMDATA, "Start model");
			d.put("start_model_weight", EMObject::FLOAT, "start_model_weight");
			return d;
		}

		/** Get the quality score that has been determined for this slice
		 * this is used in e2make3d.py to print information to std out
		 * @return the snr normalized Fourier ring correlation score
		 * @param idx the index of the slice in the quality_scores vector
		 * @exception GenericException(just throw) when the idx is beyond the range of the quality_scores vector
		 */
		virtual float get_score(const unsigned int idx) { if ( quality_scores.size() > idx ) return quality_scores[idx].get_snr_normed_frc_integral(); else throw UnexpectedBehaviorException("The requested index was beyond the length of the quality scores vector."); }

		/** Get the normalization value that has been determined for this slice
		 * this is used in e2make3d.py to print information to std out
		 * @return the normalization value
		 * @param idx the index of the slice in the quality_scores vector
		 * @exception GenericException(throw) when the idx is beyond the range of the quality_scores vector
		 */
		virtual float get_norm(const unsigned int idx) { if ( quality_scores.size() > idx ) return quality_scores[idx].get_norm();  else throw UnexpectedBehaviorException("The requested index was beyond the length of the quality scores vector."); }

		virtual void zero_memory();

	  protected:
	  	/** Preprocess the slice prior to insertion into the 3D volume
		 * this Fourier tranforms the slice and make sure all the pixels are in the right position
	  	 * @return the processed slice
	  	 * @param slice the slice to be prepocessed
	  	 * @param t transform
	  	 * @exception InvalidValueException when the specified padding value is less than the size of the images
		 */
		EMData* preprocess_slice( const EMData* const slice, const Transform& t = Transform() );
		
		/** Load default settings
		*/
		void load_default_settings()
		{
// 			params.put("x_in", = 0;
// 			params.put("y_in", = 0;
// 			params["zsample"] = 0;
// 			params["ysample"] = 0;
// 			params["xsample"] = 0;
// 			params["mode"] = 2;
// 			params["hard"] = 0.05;
// 			params["sym"] = "c1";
// 			params["quiet"] = true;
// 			params["3damp"] = false;
		}

		/** Frees the memory owned by this object (but not parent objects)
		 * Deletes the FourierPixelInserter3D pointer
		 */
		void free_memory();

		/** Load the pixel inserter based on the information in params
		 */
		void load_inserter();

		/** Load the pixel inserter based on the information in params
		 */
		void load_interp_FRC_calculator();

		/** A function to perform the nuts and bolts of inserting an image slice
		 * @param input_slice the slice to insert into the 3D volume
		 * @param euler a transform3D storing the slice euler angle
		 */
		void do_insert_slice_work(const EMData* const input_slice, const Transform & euler);

		/** print stats is called internally at various points in the reconstruction routine and is for the benefit of people using e2make3d.py
		 */
		void print_stats(const vector<InterpolatedFRC::QualityScores>& scores);

		/// Quality scores vectors for storing normalization constants, and SNR normalized Fourier ring correlation scores
		vector<InterpolatedFRC::QualityScores> quality_scores, prev_quality_scores;

		/// Index used internally to index into the quality scores vectors at different points ( in insert_slice and determine_slice_agreement)
		unsigned int image_idx;

		/// A pixel inserter pointer which inserts pixels into the 3D volume using one of a variety of insertion methods
		FourierPixelInserter3D* inserter;

		/// A pixel inserter pointer which inserts pixels into the 3D volume using one of a variety of insertion methods
		InterpolatedFRC* interp_FRC_calculator;

		/// Internal flags used to perform memory zeroing and normalization transparently
		bool slice_insertion_flag;
		bool slice_agreement_flag;

		/// Used for scaling frequency axes when any of the x_out, y_out or z_out parameters are specified
		float x_scale_factor, y_scale_factor, z_scale_factor;

		/// Keeps track of the maximum dimension size
		int max_input_dim;

		/// Keeps track of the eventual size of the output real space volume
		int output_x, output_y, output_z;
	  private:
		 /** Disallow copy construction
  		 */
  		FourierReconstructor( const FourierReconstructor& that );
  		/**Disallow assignment
  		 */
  		FourierReconstructor& operator=( const FourierReconstructor& );

	};

	/** A more mathematically rigorou Fourier reconstructor
	 * That did not seem to outperform the FourierReconstructor
	 * @author David Woolford (programming) and Phil Baldwin (mathematics)
	 * @date early 2008  
	 */
	class BaldwinWoolfordReconstructor : public FourierReconstructor
	{
		public:
		BaldwinWoolfordReconstructor() : W(0) {}

		/** Deconstructor
		 */
		virtual ~BaldwinWoolfordReconstructor() {
			if ( W != 0 )
				delete [] W;
		}

		/** Get the unique name of the reconstructor
		 */
		virtual string get_name() const
		{
			return "baldwinwoolford";
		}

		/** Get the one line description of the reconstructor
		 */
		virtual string get_desc() const
		{
			return "Reconstruction via direct Fourier inversion using gridding and delta function based weights";
		}

		/** Factory themed method allocating a new FourierReconstructor
		 * @return a Reconstructor pointer to a newly allocated FourierReconstructor
		 */
		static Reconstructor *NEW()
		{
			return new BaldwinWoolfordReconstructor();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("mode", EMObject::INT, "Optional. Fourier pixel insertion mode [1-8] - mode 2 is default.");
			d.put("x_in", EMObject::INT, "Necessary. The x dimension of the input images.");
			d.put("y_in", EMObject::INT, "Necessary. The y dimension of the input images.");
			d.put("zsample", EMObject::INT, "Optional. The z dimension (Fourier sampling) of the reconstructed volume, very useful for tomographic reconstruction. Works for general volumes.");
			d.put("ysample", EMObject::INT, "Optional. The y dimension (Fourier sampling) of the reconstructed volume, works for general volumes. Not commonly specified.");
			d.put("xsample", EMObject::INT, "Optional. The x dimension (Fourier sampling) of the reconstructed volume, works for general volumes. Not commonly specified.");
			d.put("sym", EMObject::STRING, "Optional. The symmetry of the reconstructed volume, c?, d?, oct, tet, icos, h?. Default is c1");
			d.put("maskwidth", EMObject::INT, "The width of the Fourier space kernel used to interpolate data to grid points" );
			d.put("postmultiply", EMObject::BOOL, "A flag that controls whether or not the reconstructed volume is post multiplied in real space by the IFT of the weighting function. Default is on");
			// Currently redundant
			d.put("3damp", EMObject::BOOL, "this doesn't work, fixme dsaw");
			d.put("hard", EMObject::FLOAT, "Optional. The quality metric threshold. Default is 0 (off).");
			d.put("quiet", EMObject::BOOL, "Optional. Toggles writing useful information to standard out. Default is false.");
			return d;
		}
		/** Finish reconstruction and return the complete model.
		 * @return The result 3D model.
		 */
		virtual EMData *finish();

		virtual int insert_slice_weights(const Transform& t3d);

		virtual int insert_slice(const EMData* const image, const Transform3D& t3d) {
			throw UnexpectedBehaviorException("Transform3D interface is redundant");
		}
		virtual int insert_slice(const EMData* const image, const Transform& t3d);
		void insert_pixel(const float& x, const float& y, const float& z, const float dt[2]);

		void insert_density_at(const float& x, const float& y, const float& z);

		virtual void setup();

		protected:
		/** Load default settings
		*/
		void load_default_settings()
		{
			params["mode"] = 1;
			params["x_in"] = 0;
			params["y_in"] = 0;
			params["zsample"] = 0;
			params["ysample"] = 0;
			params["xsample"] = 0;
			params["sym"] = "c1";
			params["maskwidth"] = 3;

			// Currently redundant
			params["3damp"] = false;
			params["hard"] = 0.05;
			params["quiet"] = false;
		}

		private:
		/** Disallow copy construction
  		 */
		BaldwinWoolfordReconstructor( const BaldwinWoolfordReconstructor& that );
  		/**Disallow assignment
  		 */
		BaldwinWoolfordReconstructor& operator=( const BaldwinWoolfordReconstructor& );

		float* W;
		float dfreq;

	};

	/** Fourier space 3D reconstruction with slices already Wiener filter processed.
     */
	class WienerFourierReconstructor:public Reconstructor, public ReconstructorVolumeData
	{
	  public:
		WienerFourierReconstructor() { load_default_settings(); };
		virtual ~WienerFourierReconstructor() {};

		virtual void setup();

		virtual int insert_slice(const EMData* const slice, const Transform3D & euler);

		virtual EMData *finish();

		virtual string get_name() const
		{
			return "wiener_fourier";
		}

		virtual string get_desc() const
		{
			return "Experimental - Direct Fourier reconstruction taking into account the Wiener filtration of the individual images.";
		}

		static Reconstructor *NEW()
		{
			return new WienerFourierReconstructor();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			// FIXME:: double check all of thes are need, expecially dlog and weight
			d.put("size", EMObject::INT);
			d.put("mode", EMObject::INT);
			d.put("weight", EMObject::FLOAT);
			d.put("use_weights", EMObject::BOOL);
			d.put("dlog", EMObject::BOOL);
			d.put("padratio", EMObject::FLOAT);
			d.put("snr", EMObject::FLOATARRAY);
			return d;
		}
	  private:
		// Disallow copy construction
		WienerFourierReconstructor( const WienerFourierReconstructor& that);
		// Disallow assignment
		WienerFourierReconstructor& operator=( const WienerFourierReconstructor& );

		void load_default_settings()
		{
			params["size"] = 0;
			params["mode"] = 2;
			params["weight"] = 1.0;
			params["use_weights"] = true;
			params["dlog"] = false;
			params["padratio"] = 1.0;
		}
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

		virtual int insert_slice(const EMData* const slice, const Transform3D & euler) {
			throw UnexpectedBehaviorException("Use of Tranform3D with BackProjection is deprecated"); }

		virtual int insert_slice(const EMData* const slice, const Transform & euler);

		virtual EMData *finish();

		virtual string get_name() const
		{
			return "back_projection";
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

		virtual int insert_slice(const EMData* const slice, const Transform & euler);

 	       virtual EMData *finish();

		virtual string get_name() const
		{
			return "nn4";
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


	  private:
		EMData* m_volume;
		EMData* m_wptr;
		EMData* m_result;
		bool m_delete_volume;
		bool m_delete_weight;
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

		virtual int insert_slice(const EMData* const slice, const Transform & euler);

		virtual EMData* finish();

		virtual string get_name() const
		{
			return "nnSSNR";
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
			d.put("w", EMObject::FLOAT);
			return d;
		}

		void setup( const string& symmetry, int size, int npad);

		int insert_padfft_slice( EMData* padded, const Transform& trans, int mult=1 );

	  private:
		EMData* m_volume;
		EMData* m_wptr;
		EMData* m_wptr2;
		EMData* m_result;
		bool m_delete_volume;
		bool m_delete_weight;
		bool m_delete_weight2;
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

		virtual int insert_slice(const EMData* const slice, const Transform& euler);

		virtual EMData *finish();

		virtual string get_name() const
		{
			return "nn4_ctf";
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
            d.put("weighting",  EMObject::INT);
            d.put("varsnr",     EMObject::INT);
			return d;
		}

		void setup( const string& symmetry, int size, int npad, float snr, int sign );

		int insert_padfft_slice( EMData* padfft, const Transform& trans, int mult=1);

		int insert_buffed_slice( const EMData* buffer, int mult );
	  private:
		EMData* m_volume;
		EMData* m_result;
		EMData* m_wptr;
		bool m_delete_volume;
		bool m_delete_weight;
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

		virtual int insert_slice(const EMData *const slice, const Transform& euler);


		virtual EMData* finish();

		virtual string get_name() const
		{
			return "nnSSNR_ctf";
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
			d.put("w",        EMObject::FLOAT);
			d.put("sign",     EMObject::INT);
			d.put("snr",      EMObject::FLOAT);
			return d;
		}
		void setup( const string& symmetry, int size, int npad, float snr, int sign);

                int insert_padfft_slice( EMData* padded, const Transform& trans, int mult=1 );

	  private:
		EMData* m_volume;
		EMData* m_wptr;
		EMData* m_wptr2;
		EMData* m_wptr3;
		EMData* m_result;
		bool m_delete_volume;
		bool m_delete_weight;
		bool m_delete_weight2;
		bool m_delete_weight3;
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
