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

	/** A class object encapsulating the volume data required by Reconstructors
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
		inline ReconstructorVolumeData() : image(0), tmp_data(0), nx(0), ny(0), nz(0) {}
		virtual ~ReconstructorVolumeData() { free_memory(); }
		
		const EMData* const get_emdata() { return image; }
	  protected: 
		EMData* image;
		//tmp_data is the substitute of misused parent in reconstruction
		//the memory will be allocated in setup() and released in finish()
		EMData* tmp_data;
		
		int nx;
		int ny;
		int nz;
		
	  protected:
		void free_memory()
		{
			if (image != 0)  {delete image; image = 0;} 
			if ( tmp_data != 0 ) { delete tmp_data; tmp_data = 0; }
		}

		virtual void normalize_threed();
		
		void zero_memory() 
		{
			if (tmp_data != 0 ) tmp_data->to_zero();
			if (image != 0 ) image->to_zero();
		}
		
		private:
		// Disallow copy construction
		ReconstructorVolumeData(const ReconstructorVolumeData& that);
		// Disallow  assignment
		ReconstructorVolumeData& operator=(const ReconstructorVolumeData& );

	};

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
	class Reconstructor : public ReconstructorVolumeData
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
		virtual int insert_slice(const EMData* const slice, const Transform3D & euler) = 0;
		
		/** 
	  	 * @return 
	  	 * @param input_slice
	  	 * @param arg
	  	 * @param num_particles_in_slice
	  	 * @exception 
		 */
		virtual int determine_slice_agreement(const EMData* const input_slice, const Transform3D & arg, const unsigned int  num_particles_in_slice = 1)
		{  
			cout << "You called determine slice agreement but nothing happened - there is no functionality for determing slice agreement using this " << get_name() << " reconstructor" << endl;
			return 0;
		}


		/** Finish reconstruction and return the complete model.
		 * @return The result 3D model.
		 */
		virtual EMData *finish() = 0;

		/** Get the Reconstructor's name. Each Reconstructor is
		 * identified by a unique name.
		 * @return The Reconstructor's name.
		 */
		virtual string get_name() const = 0;

		
		virtual string get_desc() const = 0;

		/** Get the Reconstructor's parameters in a key/value dictionary.
		 * @return A key/value pair dictionary containing the parameters.
		 */
		Dict get_params() const
		{
			return params;
		}

		/** Set the Reconstructor's parameters using a key/value dictionary.
		 * @param new_params A dictionary containing the new parameters.
		 */
		void set_params(const Dict & new_params)
		{
			// This commented out by d.woolford. 
			//params.clear();
			
			insert_params(new_params);
		}
		
		/** Insert the Reconstructor's parameters using a key/value dictionary.
		 * @param new_params A dictionary containing the new parameters.
		 * @exception InvalidParameterException thrown when the parameter being inserted is not in the list of permissable params
		 */
		void insert_params(const Dict & new_params)
		{
			// this is really inserting OR individually replacing...
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
		/** Get reconstructor parameter information in a dictionary. 
		 * Each parameter has one record in the dictionary. Each 
		 * record contains its name, data-type, and description.
		 *
		 * @return A dictionary containing the parameter info.
		 */	 
		virtual TypeDict get_param_types() const = 0;

		EMObject& operator[]( const string& key ) { return params[key]; }

		// A function used only by the Fourier reconstructor for testing and writing to std out purposes in e2make3d.py
		virtual float get_score(const unsigned int idx) { return 0; }
		// A function used only by the Fourier reconstructor for testing and writing to std out purposes in e2make3d.py
		virtual float get_norm(const unsigned int idx) { return 0; }

	  protected:
		mutable Dict params;
		
	  private:
		// Disallow copy construction
		Reconstructor(const Reconstructor& that);
		Reconstructor& operator=(const Reconstructor& );
	
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
	class FourierReconstructor : public Reconstructor
	{
	  public:
		/** Default constructor
		* calls load_default_settings()
		*/
		FourierReconstructor() : image_idx(0), inserter(0), interpFRC_calculator(0), slice_insertion_flag(true),
		slice_agreement_flag(false), x_scale_factor(0.0), y_scale_factor(0.0), z_scale_factor(0.0), 
		max_padded_dim(0), output_x(0), output_y(0), output_z(0) { load_default_settings(); }
		
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
		* @param t3d the Transform3D that stores the image Euler angles
		* @exception NullPointerException if the input EMData pointer is null
		* @exception ImageFormatException if the image is complex as opposed to real
		*/
		virtual int insert_slice(const EMData* const slice, const Transform3D & t3d);
		
		/** Determine slice agreement with the current reconstruction
		* @return 0 if successful, 1 otherwise
		* @param input_slice the image slice to compared against the 3D volume
		* @param t3d the Transform3D that the image Euler angles
		* @param num_particles_in_slice the number of particles in the slice - used to determine the SNR normalized FSC, defaults to 1.
		* @exception NullPointerException if the input EMData pointer is null
		* @exception ImageFormatException if the image is complex as opposed to real
		*/
		virtual int determine_slice_agreement(const EMData* const input_slice, const Transform3D & t3d, const unsigned int  num_particles_in_slice = 1);

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

		/** Factory themed method allocating a new FourierReconstructor
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
			d.put("pad", EMObject::INT, "Optional. The amount to pad the input images to - should be greater than the image size.");
			d.put("x_pad", EMObject::INT, "Optional. The amount to pad the input images to in the x direction - should be greater than the image x size.");
			d.put("y_pad", EMObject::INT, "Optional. The amount to pad the input images to in the y direction - should be greater than the image y size.");
			d.put("mode", EMObject::INT, "Optional. Fourier pixel insertion mode [1-7] - mode 2 is default.");
			d.put("weight", EMObject::FLOAT, "Optional. A temporary weight variable, used to weight slices as they are inserted. Default is 1");
			d.put("hard", EMObject::FLOAT, "Optional. The quality metric threshold. Default is 0 (off).");
			d.put("edgenorm", EMObject::BOOL, "Optional. Whether or not to perform edge normalization on the inserted slices before Fourier transforming them. Default is true." );
			d.put("sym", EMObject::STRING, "Optional. The symmetry of the reconstructed volume, c?, d?, oct, tet, icos, h?. Default is c1");
			
			d.put("apix", EMObject::FLOAT, "Optional. Angstrom per pixel of the input images, default is 1.0.");
			d.put("tomo_weight", EMObject::BOOL, "Optional. A tomographic reconstruction flag that causes inserted slices to be weighted by 1/cos(alt) - alt is the tilt angle. Default is false.");
			d.put("tomo_mask", EMObject::BOOL, "Optional. A tomographic reconstruction flag that causes inserted slices to have their edge pixels masked according to tilt angle, ensuring that each projection image depicts approximately the same volume, default is false." );
			d.put("t_emm", EMObject::BOOL, "Optional. Read as tomo edge mean mask - experimental, default false");
			d.put("t_emm_gauss", EMObject::INT, "Optional. An optional gaussain fall off for tomo_mask" );
			d.put("dlog", EMObject::BOOL, "Optional. This is a residual parameter from EMAN1 that has not yet been addressed in the EMAN2 implementation.");
			d.put("quiet", EMObject::BOOL, "Optional. Toggles writing useful information to standard out. Default is false.");
			d.put("3damp", EMObject::BOOL, "Optional. Toggles writing the 3D FFT amplitude image. Default is false.");
			return d;
		}

		/** Get the quality score that has been determined for this slice
		 * this is used in e2make3d.py to print information to std out
		 * @return the snr normalized Fourier ring correlation score
		 * @param idx the index of the slice in the quality_scores vector
		 * @exception GenericException(just throw) when the idx is beyond the range of the quality_scores vector
		 */
		virtual float get_score(const unsigned int idx) { if ( quality_scores.size() > idx ) return quality_scores[idx].get_snr_normed_frc_integral(); else throw; }
		
		/** Get the normalization value that has been determined for this slice
		 * this is used in e2make3d.py to print information to std out
		 * @return the normalization value
		 * @param idx the index of the slice in the quality_scores vector
		 * @exception GenericException(throw) when the idx is beyond the range of the quality_scores vector
		 */
		virtual float get_norm(const unsigned int idx) { if ( quality_scores.size() > idx ) return quality_scores[idx].get_norm();  else throw; }
		
	  protected:
	  	/** Preprocess the slice prior to insertion into the 3D volume
		 * this Fourier tranforms the slice and make sure all the pixels are in the right position
	  	 * @return the processed slice
	  	 * @param slice the slice to be prepocessed
	  	 * @param transform
	  	 * @exception InvalidValueException when the specified padding value is less than the size of the images
		 */
		EMData* preprocess_slice( const EMData* const slice, const Transform3D transform = Transform3D() );
	  private:
		// Disallow copy construction
		FourierReconstructor( const FourierReconstructor& that );
		// Disallow assignment
		FourierReconstructor& operator=( const FourierReconstructor& );
		
		/** Load default settings
		*/
		void load_default_settings()
		{
			params["x_in"] = 0;
			params["y_in"] = 0;
			params["zsample"] = 0;
			params["ysample"] = 0;
			params["xsample"] = 0;
			params["x_pad"] = 0;
			params["y_pad"] = 0;
			params["mode"] = 2;
			params["weight"] = 1.0;
			params["dlog"] = false;
			params["hard"] = 0.05;
			params["sym"] = "c1";
			params["apix"] = 1.0;
			params["tomo_weight"] = false;
			params["tomo_mask"] = false;
			params["t_emm"] = false;
			params["edgenorm"] = true;
			params["t_emm_gauss"] = 10;
			params["quiet"] = false;
			params["3damp"] = false;
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
		void load_interpFRC_calculator();
		
		/** A function to perform the nuts and bolts of inserting an image slice
		 * @param input_slice the slice to insert into the 3D volume
		 * @param euler a transform3D storing the slice euler angle
		 */
		void do_insert_slice_work(const EMData* const input_slice, const Transform3D & euler);
		
		/** print stats is called internally at various points in the reconstruction routine and is for the benefit of people using e2make3d.py
		 */
		void print_stats(const vector<QualityScores>& scores);
		
		/// Quality scores vectors for storing normalization constants, and SNR normalized Fourier ring correlation scores
		vector<QualityScores> quality_scores, prev_quality_scores;
		
		/// Index used internally to index into the quality scores vectors at different points ( in insert_slice and determine_slice_agreement)
		unsigned int image_idx;
		
		/// A pixel inserter pointer which inserts pixels into the 3D volume using one of a variety of insertion methods
		FourierPixelInserter3D* inserter;
		
		/// A pixel inserter pointer which inserts pixels into the 3D volume using one of a variety of insertion methods
		InterpolatedFRC* interpFRC_calculator;

		/// Internal flags used to perform memory zeroing and normalization transparently
		bool slice_insertion_flag;
		bool slice_agreement_flag;
		
		/// Used for scaling frequency axes when any of the x_out, y_out or z_out parameters are specified
		float x_scale_factor, y_scale_factor, z_scale_factor;
		
		/// Keeps track of the maximum dimension size
		int max_padded_dim;
		
		/// Keeps track of the eventual size of the output real space volume
		int output_x, output_y, output_z;
	};

	/** Fourier space 3D reconstruction with slices already Wiener filter processed.
     */
	class WienerFourierReconstructor:public Reconstructor
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
	class BackProjectionReconstructor:public Reconstructor
	{
	  public:
		BackProjectionReconstructor() { load_default_settings();  }
	
		virtual ~BackProjectionReconstructor() {}

		virtual void setup();
		
		virtual int insert_slice(const EMData* const slice, const Transform3D & euler);
		
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
		
		EMData* preprocess_slice(const EMData* const slice);
	};


	/** Direct Fourier inversion Reconstructor
     * 
     */
	EMData* padfft_slice( const EMData* const slice, int npad );

	class nn4Reconstructor:public Reconstructor
	{
	  public:
		nn4Reconstructor();

		nn4Reconstructor( const string& symmetry, int size, int npad );

		virtual ~nn4Reconstructor();

		virtual void setup();

		virtual int insert_slice(const EMData* const slice, const Transform3D & euler);

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
			d.put("symmetry",	EMObject::STRING);
			d.put("snr",		EMObject::FLOAT);
			d.put("fftvol",		EMObject::EMDATA);
			d.put("weight",		EMObject::EMDATA);
            d.put("weighting",      EMObject::INT);
			return d;
		}

		void setup( const string& symmetry, int size, int npad );

		int insert_padfft_slice( EMData* padded, const Transform3D& trans, int mult=1 );


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
		int m_vnzp, m_vnyp, m_vnxp;
		int m_vnzc, m_vnyc, m_vnxc;
		void buildFFTVolume();
		void buildNormVolume();
		float m_wghta;
		float m_wghtb;

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

		virtual int insert_slice(const EMData* const slice, const Transform3D & euler);
		
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

		int insert_padfft_slice( EMData* padded, const Transform3D& trans, int mult=1 );


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


	class bootstrap_nnReconstructor:public Reconstructor
	{
	  public:
		bootstrap_nnReconstructor();

		~bootstrap_nnReconstructor();

		virtual void setup();

		virtual int insert_slice(const EMData* const slice, const Transform3D & euler);
		
		virtual EMData *finish();

		virtual string get_name() const
		{
			return "bootstrap_nn";
		}
		
		string get_desc() const
		{
			return "Bootstrap nearest neighbour constructor";
		}

		static Reconstructor *NEW()
		{
			return new bootstrap_nnReconstructor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("size", EMObject::INT);
			d.put("npad", EMObject::INT);
			d.put("mult", EMObject::INTARRAY);
			d.put("media", EMObject::STRING);
			d.put("symmetry", EMObject::STRING);
			return d;
		}

	  private:
		string m_ctf;
		string m_media;
		string m_symmetry;
		float m_snr;
		int m_size;
		int m_npad;
		int m_nsym;
		vector< EMData* > m_padffts;
		vector< Transform3D* > m_transes;
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

		virtual int insert_slice(const EMData* const slice, const Transform3D & euler);

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
                        d.put("weighting",      EMObject::INT);
			return d;
		}

		void setup( const string& symmetry, int size, int npad, float snr, int sign );

		int insert_padfft_slice( EMData* padfft, const Transform3D& trans, int mult=1);

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

		virtual int insert_slice(const EMData *const slice, const Transform3D & euler);
		

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

                int insert_padfft_slice( EMData* padded, const Transform3D& trans, int mult=1 );

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

	class bootstrap_nnctfReconstructor:public Reconstructor
	{
	  public:
		bootstrap_nnctfReconstructor();

		~bootstrap_nnctfReconstructor();

		virtual void setup();

		virtual int insert_slice(const EMData* const slice, const Transform3D & euler);
		
		virtual EMData *finish();

		virtual string get_name() const
		{
			return "bootstrap_nnctf";
		}
		
		string get_desc() const
		{
			return "Bootstrap nearest neighbour CTF constructor";
		}

		static Reconstructor *NEW()
		{
			return new bootstrap_nnctfReconstructor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("snr", EMObject::FLOAT);
			d.put("size", EMObject::INT);
			d.put("npad", EMObject::INT);
			d.put("sign", EMObject::INT);
			d.put("mult", EMObject::INTARRAY);
			d.put("media", EMObject::STRING);
			d.put("symmetry", EMObject::STRING);
			return d;
		}

	  private:
		string m_ctf;
		string m_media;
		string m_symmetry;
		float m_snr;
		int m_size;
		int m_npad;
		int m_nsym;
		int m_sign;
		vector< EMData* > m_padffts;
		vector< Transform3D* > m_transes;
	};

	template <> Factory < Reconstructor >::Factory();

	void dump_reconstructors();
	map<string, vector<string> > dump_reconstructors_list();

	class file_store
	{
	  public: 
		file_store(const string& filename, int npad, int write);

		virtual ~file_store();

		void add_image(EMData* data);

		void get_image(int id, EMData* padfft);

		void restart();
	  private:
		shared_ptr<std::ifstream> m_ihandle;
		shared_ptr<std::ofstream> m_bin_ohandle;
		shared_ptr<std::ofstream> m_txt_ohandle;
		string m_bin_file;
		string m_txt_file;
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
