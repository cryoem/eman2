/**
 * $Id$
 */

/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
 * Contributing Author: Dave Woolford, 06/01/2007 (woolford@bcm.edu)
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

namespace EMAN
{

	class Transform3D;
	class EMData;
	
	class QualityScores
	{
	public:
		QualityScores() : frc_integral(0), snr_normed_frc_intergral(0), normed_snr_integral(0) {}
		QualityScores( const QualityScores& that ) : frc_integral(that.frc_integral), 
		snr_normed_frc_intergral(that.snr_normed_frc_intergral), normed_snr_integral(that.normed_snr_integral) {}
		QualityScores& operator=( const QualityScores& that ) 
		{
			frc_integral = that.frc_integral; 
			snr_normed_frc_intergral = that.snr_normed_frc_intergral;
			normed_snr_integral  = that.normed_snr_integral;
			return *this;
		}
		
		~QualityScores() {}

		float get_frc_integral() { return frc_integral; }
		float get_snr_normed_frc_integral() { return snr_normed_frc_intergral; }
		float get_normed_snr_integral() { return normed_snr_integral; }

		void set_frc_integral( const float& score ) { frc_integral = score; }
		void set_snr_normed_frc_integral(const float& score) { snr_normed_frc_intergral = score; }
		void set_normed_snr_integral(const float& score) { normed_snr_integral = score; }

		void debug_print()
		{
			cout << "frc " << frc_integral << " nfrc " << snr_normed_frc_intergral << " nsnr " << normed_snr_integral << endl;
		}
	private:

		float frc_integral, snr_normed_frc_intergral, normed_snr_integral;
		
	};

	class InterpolatedFRC
	{
	public:
		InterpolatedFRC() : threed_rdata(0), frc(0), frc_norm_rdata(0), frc_norm_dt(0), size(0), pixel_radius_max(0) {}

		InterpolatedFRC(float* const rdata, const int xsize, const int ysize, const int zsize, const float& sampling=1.0 );
		~InterpolatedFRC()
		{
			free_memory();
		}

		// Copy and assignment
		InterpolatedFRC( const InterpolatedFRC& that );
		InterpolatedFRC& operator=( const InterpolatedFRC& that);
	
		bool continue_frc_calc(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1.0);
	
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

		int off[8];
		float g[8];
	};

	/** A Reconstructor containing an EMData pointer and x,y,z dimensions
	 *  This is a pseudo abstract base class 
     *  because it derives from Reconstructor but does not implement 
	 *  any of its pure virtual functions. It basically stores 
	 *  a pointer to an image object and stores the dimensions of the image volume.
	 *  This class was originally added simply to encapsulate the 
	 *  the things common to FourierReconstructor, WienerFourierReconstructor
	 *  and BackProjectionReconstructor
	 *  d.woolford May 2007
     */

	class ReconstructorWithVolumeData
	{
	  public:
		inline ReconstructorWithVolumeData() : image(0), tmp_data(0), nx(0), ny(0), nz(0) {}
		inline virtual ~ReconstructorWithVolumeData() { free_memory(); }

		/** Copy constructor
		 */
		ReconstructorWithVolumeData(const ReconstructorWithVolumeData& that) { copyData(that); }

		/** Assignment operator
		 */
		ReconstructorWithVolumeData& operator=(const ReconstructorWithVolumeData& );

		void zero_memory() 
		{
			if (tmp_data != 0 ) tmp_data->to_zero();
			if (image != 0 ) image->to_zero();
		}
	  protected:
		EMData * image;
		//tmp_data is the substitute of misused parent in reconstruction
		//the memory will be allocated in setup() and released in finish()
		EMData* tmp_data;
		
		int nx;
		int ny;
		int nz;

	  private:
		/** a convenience function for copying data, called by the assigment operator and the copy constructor,
		 *  deep and complete complete.
		 */
		void copyData( const ReconstructorWithVolumeData& that );

		void free_memory()
		{
			if (image)  {delete image; image = 0;} 
			if ( tmp_data != 0 ) { delete tmp_data; tmp_data = 0; }
		}
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
	class Reconstructor : public ReconstructorWithVolumeData
	{
	  public:
		Reconstructor() {}
		inline virtual ~Reconstructor()
		{
		}

		/** Copy constructor
		 */
		Reconstructor(const Reconstructor&);

		/** Assignment operator
		 */
		Reconstructor& operator=(const Reconstructor& );


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
		virtual int determine_slice_agreement(const EMData* const input_slice, const Transform3D & arg, const unsigned int  num_particles_in_slice = 1) { return 0;}


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

		// A function used only in fourier reconstructor atm for reseting information in between iterations
		virtual void iteration_reset() {}
		// A function used only by the Fourier reconstructor for testing purposes 
		virtual float get_score(const unsigned int idx) { return 0; }

	  protected:
		mutable Dict params;

	  private:
		/** a convenience function for copying data, called by the assigment operator and the copy constructor
		 */
		void copyData( const Reconstructor& );
	
	};


	/** Fourier space 3D reconstruction
     */
	class FourierReconstructor:public Reconstructor
	{
	  public:
		FourierReconstructor() { load_default_settings(); }
		virtual ~FourierReconstructor() {}
	
		/** Copy constructor
		 */
		FourierReconstructor( const FourierReconstructor& that ) : Reconstructor(that)
		 { load_default_settings(); }

		/** Assignment operator
		 */
		FourierReconstructor& operator=( const FourierReconstructor& );

		/** Assignment operator
		 * @exception InvalidValueException
		 */
		virtual void setup();

		virtual int insert_slice(const EMData* const slice, const Transform3D & euler);

		/** 
	  	 * @return 
	  	 * @param 
	  	 * @exception 
		 */
		virtual int determine_slice_agreement(const EMData* const input_slice, const Transform3D & arg, const unsigned int  num_particles_in_slice);

		virtual	EMData *finish();

		virtual string get_name() const
		{
			return "fourier";
		}
		
		virtual string get_desc() const
		{
			return "Reconstruction via direct Fourier methods using a Gaussian kernel";
		}

		static Reconstructor *NEW()
		{
			return new FourierReconstructor();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("size", EMObject::INT);
			d.put("mode", EMObject::INT);
			d.put("weight", EMObject::FLOAT);
			d.put("hard", EMObject::FLOAT);
			d.put("use_weights", EMObject::BOOL);
			d.put("dlog", EMObject::BOOL);
			d.put("sym", EMObject::STRING);
			d.put("pad", EMObject::INT);
			d.put("apix", EMObject::FLOAT);
			return d;
		}
		
		virtual void iteration_reset() { quality_scores.clear(); idx = 0;}

		virtual float get_score(const unsigned int idx) { if ( quality_scores.size() > idx ) return quality_scores[idx].get_frc_integral(); else {cout << "foo" << endl; return 2;}}
	  protected:
	  	/** Preprocess the slice prior to insertion into the 3D volume
		 * this Fourier tranforms the slice and make sure all the pixels are in the right positions
	  	 * @return the processed slice
	  	 * @param slice the slice to be prepocessed
	  	 * @exception InvalidValueException
		 */
	  	EMData* preprocess_slice( const EMData* const slice );


	  private:
		void load_default_settings()
		{
			params["size"] = 0;
			params["mode"] = 2;
			params["weight"] = 1.0;
			params["use_weights"] = true;
			params["dlog"] = false;
			params["hard"] = 0.07;
			params["sym"] = "unknown";
		}

		vector<QualityScores> quality_scores;
		unsigned int idx;
	};

	/** Fourier space 3D reconstruction with slices already Wiener filter processed.
     */
	class WienerFourierReconstructor:public Reconstructor
	{
	  public:
		WienerFourierReconstructor() { load_default_settings(); };
		virtual ~WienerFourierReconstructor() {};
	
		/** Copy constructor
		 */
		WienerFourierReconstructor( const WienerFourierReconstructor& that) : Reconstructor(that) {}

		/** Assignment operator
		 */
		WienerFourierReconstructor& operator=( const WienerFourierReconstructor& );


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

		/** Copy constructor
		 */
		BackProjectionReconstructor( const BackProjectionReconstructor& that) : Reconstructor(that) {}

		/** Assignment operator
		 */
		BackProjectionReconstructor& operator=( const BackProjectionReconstructor& );

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
			d.put("size", EMObject::INT);
			d.put("weight", EMObject::FLOAT);
			d.put("use_weights", EMObject::BOOL);
			return d;
		}
	  private:
		void load_default_settings()
		{
			params["weight"] = 1.0; 
			params["use_weights"] = true;
			params["size"] = 0;
		}
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
			d.put("size", EMObject::INT);
			d.put("npad", EMObject::INT);
			d.put("symmetry", EMObject::STRING);
			d.put("fftvol", EMObject::EMDATA);
			d.put("weight", EMObject::EMDATA);
			d.put("use_weights", EMObject::BOOL);
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

		~nn4_ctfReconstructor();

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
			d.put("weight4",  EMObject::EMDATA);
			d.put("weight5",  EMObject::EMDATA);
			d.put("SSNR",     EMObject::EMDATA);
			d.put("w",        EMObject::FLOAT);
			d.put("sign",     EMObject::INT);
			d.put("snr",      EMObject::FLOAT);
			d.put("wiener",   EMObject::INT);
			return d;
		}
		void setup( const string& symmetry, int size, int npad, float snr, int sign);

                int insert_padfft_slice( EMData* padded, const Transform3D& trans, int mult=1 );

	  private:
		EMData* m_volume;
		EMData* m_wvolume;
		EMData* m_wptr;
		EMData* m_wptr2;
		EMData* m_wptr3;
		EMData* m_wptr4;
		EMData* m_wptr5;
		EMData* m_result;
		bool m_delete_volume;
		bool m_delete_wvolume;
		bool m_delete_weight;
		bool m_delete_weight2;
		bool m_delete_weight3;
		bool m_delete_weight4;
		bool m_delete_weight5;
	        string  m_symmetry;
		int m_weighting;
		int m_vnx, m_vny, m_vnz;		
		int m_npad;
		int m_nsym;
		int m_vnzp, m_vnyp, m_vnxp;
		int m_vnzc, m_vnyc, m_vnxc;
		void buildFFTVolume();
		void buildWFFTVolume();
		void buildNormVolume();
		void buildNorm2Volume();
		void buildNorm3Volume();
		void buildNorm4Volume();
		void buildNorm5Volume();
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

	/** FourierPixelInserter3D - encapsulates a method for inserting a (continuous) 2D Fourier pixel into a discrete 3D volume
	 * This is an abstract base class
	 * Deriving classes must over ride the pure virtual insert_pixel function.
	 * 
	 * At the moment this is only useful in conjunction with FourierReconstructor
	 * 
	 * David Woolford, June 2007.
	*/
	class FourierPixelInserter3D
	{
	public:
		/** Constructor a FourierSliceInserter
		 * @param normalize_values a block of memory equal in size to memory associated with real_data
		 * @param real_data a pointer to the memory containing the discrete 3D volume pixel data
		 * @param xsize the xsize of the discrete 3D volume
		 * @param ysize the ysize of the discrete 3D volume
		 * @param zsize the zsize of the discrete 3D volume
		 */
		FourierPixelInserter3D( float * const normalize_values, float * const real_data, const unsigned int xsize, const unsigned int ysize, const unsigned int zsize ) :
			norm(normalize_values), rdata(real_data), nx(xsize), ny(ysize), nz(zsize) {nxy = nx * ny;}
		virtual ~FourierPixelInserter3D() {}
		
		/** Insert a complex pixel [dt[0]+dt[1]i] at (float) coordinate [xx,yy,zz] with weighting into a discrete 3D volume
		 * @param xx the floating point x coordinate
		 * @param yy the floating point y coordinate
		 * @param zz the floating point z coordinate
		 * @param dt the complex pixel value (dt[0] is real, dt[1] is imaginary)
		 * @param weight the weight to given to this complex pixel
		 * @return A boolean that indicates the pixel has been inserted (or not)
		 */
		virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight) = 0;
		
	protected:
		// A pointer to the constructor argument normalize_values
		float * const norm;
		// A pointer to the constructor argument real_data
		float * const rdata;
		
		// Image volume data sizes, and nxy a convenience variable used here and there
		int nx, ny, nz, nxy;
		
	private:
		// Disallow copy and assignment by default
		FourierPixelInserter3D( const FourierPixelInserter3D& );
		FourierPixelInserter3D& operator=( const FourierPixelInserter3D& );
	};
	
	/** FourierPixelInserter3DMode1  - encapsulates "method 1" for inserting a 2D Fourier slice into a 3D volume
	 * See comments in FourierPixelInserter3D for explanations
	*/
	class FourierInserter3DMode1 : public FourierPixelInserter3D
	{
	public:
		FourierInserter3DMode1(float * const normalize_values, float * const real_data, const unsigned int xsize, const unsigned int ysize, const unsigned int zsize) :
			FourierPixelInserter3D(normalize_values, real_data, xsize, ysize, zsize) {}
		virtual ~FourierInserter3DMode1() {}
		
		virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& = 1);
		
		/** get_mode_number get the unique mode number
		 * Is static because it makes the implementation of the associated FourierPixelInserterMaker constructor
		 * independent of the this inserters mode number. Should be virtual, but static virtual functions
		 * are not allowable.
		 * @return a unique mode number
		 */
		static /* virtual */ unsigned int get_mode_number() { return 1; }
		
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
		FourierInserter3DMode2(float * const normalize_values, float * const real_data, const unsigned int xsize, const unsigned int ysize, const unsigned int zsize);
		virtual ~FourierInserter3DMode2() {}
		
		virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1);
		
		/** get_mode_number get the unique mode number
		 * Is static because it makes the implementation of the associated FourierPixelInserterMaker constructor
		 * independent of the this inserters mode number. Should be virtual, but static virtual functions
		 * are not allowable.
		 * @return a unique mode number
		 */
		static /* virtual */ unsigned int  get_mode_number() { return 2; }
		
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
		FourierInserter3DMode3(float * const normalize_values, float * const real_data, const unsigned int xsize, const unsigned int ysize, const unsigned int zsize) :
			FourierPixelInserter3D(normalize_values, real_data, xsize, ysize, zsize) {}
		virtual ~FourierInserter3DMode3() {}
		
		virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1);
		
		/** get_mode_number get the unique mode number
		 * Is static because it makes the implementation of the associated FourierPixelInserterMaker constructor
		 * independent of the this inserters mode number. Should be virtual, but static virtual functions
		 * are not allowable.
		 * @return a unique mode number
		 */
		static /* virtual */ unsigned int get_mode_number() { return 3; }
		
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
		FourierInserter3DMode4(float * const normalize_values, float * const real_data, const unsigned int xsize, const unsigned int ysize, const unsigned int zsize) :
			FourierPixelInserter3D(normalize_values, real_data, xsize, ysize, zsize) {}
		virtual ~FourierInserter3DMode4() {}
		
		virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1);
		
		/** get_mode_number get the unique mode number
		 * Is static because it makes the implementation of the associated FourierPixelInserterMaker constructor
		 * independent of the this inserters mode number. Should be virtual, but static virtual functions
		 * are not allowable.
		 * @return a unique mode number
		 */
		static /* virtual */ unsigned int get_mode_number() { return 4; }
		
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
		FourierInserter3DMode5(float * const normalize_values, float * const real_data, const unsigned int xsize, const unsigned int ysize, const unsigned int zsize) :
			FourierPixelInserter3D(normalize_values, real_data, xsize, ysize, zsize)
		{
			gimx = Interp::get_gimx();
		}
		virtual ~FourierInserter3DMode5() {}
		
		virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1);
		
		/** get_mode_number get the unique mode number
		 * Is static because it makes the implementation of the associated FourierPixelInserterMaker constructor
		 * independent of the this inserters mode number. Should be virtual, but static virtual functions
		 * are not allowable.
		 * @return a unique mode number
		 */
		static /* virtual */ unsigned int get_mode_number() { return 5; }
		
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
		FourierInserter3DMode6(float * const normalize_values, float * const real_data, const unsigned int xsize, const unsigned int ysize, const unsigned int zsize) :
			FourierPixelInserter3D(normalize_values, real_data, xsize, ysize, zsize) {}
		virtual ~FourierInserter3DMode6() {}
		
		virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1);
		
		/** get_mode_number get the unique mode number
		 * Is static because it makes the implementation of the associated FourierPixelInserterMaker constructor
		 * independent of the this inserters mode number. Should be virtual, but static virtual functions
		 * are not allowable.
		 * @return a unique mode number
		 */
		static /* virtual */ unsigned int get_mode_number() { return 6; }
		
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
		FourierInserter3DMode7(float * const normalize_values, float * const real_data, const unsigned int xsize, const unsigned int ysize, const unsigned int zsize) :
			FourierPixelInserter3D(normalize_values, real_data, xsize, ysize, zsize) {}
		virtual ~FourierInserter3DMode7() {}
		
		/** get_mode_number get the unique mode number
		 * Is static because it makes the implementation of the associated FourierPixelInserterMaker constructor
		 * independent of the this inserters mode number. Should be virtual, but static virtual functions
		 * are not allowable.
		 * @return a unique mode number
		 */
		virtual bool insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight = 1);
		
		static /* virtual */ unsigned int get_mode_number() { return 7; }
		
	private:
		// Disallow copy and assignment by default
		FourierInserter3DMode7( const FourierInserter3DMode7& );
		FourierInserter3DMode7& operator=( const FourierInserter3DMode7& );
	};
	
	/** FourierPixelInserterMaker - a factory that returns an instance of a FourierPixelInserter3D
	 * This implementation is an example of an abstract pluggable factory.
	 * It provides a static interface for performing the factory functionality (get_inserter).
	 * It also acts as an abstract base class for all FourierPixelInserterMakers inserted into maker_registry.
	 * FourierPixelInserterMakers are inserted into the maker_registry when they are constructed (see the
	 * constructor in this (the abstract base) class. Hence, any concrete instance of a FourierPixelInserterMaker
	 * should only be instantiated once, and this is achieved using static initialisation. Each concrete
	 * FourierPixelInserterMaker that derives from this class must contain a static member variable of its own
	 * type. When these static members are initialised they automatically insert themselves into maker_registry
	 * in this, the base class.
	 * 
	 * David Woolford, June 2007.
	 */
	class FourierPixelInserterMaker
	{
	public:
		/** FourierPixelInserterMaker constructor
		 * @param key the key of the inserted (child) FourierPixelInserterMaker
		 */
		FourierPixelInserterMaker( const unsigned int key )
		{
			if (maker_registry.find(key) == maker_registry.end() )
			{
				maker_registry[key] = this;
			}
			else
			{
				LOGERR("Attempted to insert a key/value pair in the FourierPixelInserterMaker registry that already existed");
				throw InvalidCallException("Attempted to insert a key/value pair in the FourierPixelInserterMaker registry that already existed");
			}
		}
		virtual ~FourierPixelInserterMaker() {}
		
		/** make_inserter get an instance of a FourierPixelInserter3D*
		 * @param normalize_values a block of memory equal in size to memory associated with real_data
		 * @param real_data a pointer to the memory containing the discrete 3D volume pixel data
		 * @param xsize the xsize of the discrete 3D volume
		 * @param ysize the ysize of the discrete 3D volume
		 * @param zsize the zsize of the discrete 3D volume
		 */
		virtual FourierPixelInserter3D* make_inserter(float * const normalize_values, float * const real_data, const unsigned int xsize, const unsigned int ysize, const unsigned int zsize) = 0;
		
		/** make_inserter get an instance of a FourierPixelInserter3D*
		 * @param key the key to get the associated maker from the maker registry
		 * @param normalize_values a block of memory equal in size to memory associated with real_data to construct the maker with 
		 * @param real_data a pointer to the memory containing the discrete 3D volume pixel data to construct the maker with
		 * @param xsize the xsize of the discrete 3D volume to construct the maker with
		 * @param ysize the ysize of the discrete 3D volume to construct the maker with
		 * @param zsize the zsize of the discrete 3D volume to construct the maker with
		 */
		static FourierPixelInserter3D* get_inserter( const unsigned int key, float * const normalize_values, float * const real_data, const unsigned int xsize,
			 const unsigned int ysize, const unsigned int zsize )
		{
			if ( maker_registry.find(key) != maker_registry.end() )
				return maker_registry[key]->make_inserter(normalize_values, real_data, xsize, ysize, zsize);
			else
				return 0;	
		}
		
	private:
		// The maker registry to underly the Factory implementation
		static map<unsigned int, FourierPixelInserterMaker* > maker_registry;
		
		// Disallow copy and assignment by default
		FourierPixelInserterMaker( const FourierPixelInserterMaker& );
		FourierPixelInserterMaker& operator=( const FourierPixelInserterMaker& );
	};
	
	/** FourierInserter3DMode1Maker  - returns a newly made (pointer to) a FourierInserter3DMode1
	 * See comments in FourierPixelInserterMaker for explanations of functionality
	 */
	class FourierInserter3DMode1Maker : public FourierPixelInserterMaker
	{
	public:
		FourierInserter3DMode1Maker(const unsigned int key = FourierInserter3DMode1::get_mode_number() ) : FourierPixelInserterMaker(key) {}
		virtual ~FourierInserter3DMode1Maker() {}
		
		virtual FourierPixelInserter3D* make_inserter(float * const normalize_values, float * const real_data, const unsigned int xsize, const unsigned int ysize, const unsigned int zsize)
		{
			return new FourierInserter3DMode1(normalize_values, real_data, xsize, ysize, zsize);
		}
	private:
		static const FourierInserter3DMode1Maker registerMode1Maker;
		
		// Disallow copy and assignment by default
		FourierInserter3DMode1Maker( const FourierInserter3DMode1Maker& );
		FourierInserter3DMode1Maker& operator=( const FourierInserter3DMode1Maker& );
	};
	
	/** FourierInserter3DMode2Maker  - returns a newly made (pointer to) a FourierInserter3DMode2
	 * See comments in FourierPixelInserterMaker for explanations of functionality
	 */
	class FourierInserter3DMode2Maker : public FourierPixelInserterMaker
	{
	public:
		FourierInserter3DMode2Maker(const unsigned int key = FourierInserter3DMode2::get_mode_number() ): FourierPixelInserterMaker(key) {}
		virtual ~FourierInserter3DMode2Maker() {}
		
		virtual FourierPixelInserter3D* make_inserter(float * const normalize_values, float * const real_data, const unsigned int xsize, const unsigned int ysize, const unsigned int zsize)
		{
			return new FourierInserter3DMode2(normalize_values, real_data, xsize, ysize, zsize);
		}
	private:
		static const FourierInserter3DMode2Maker registerMode2Maker;
		
		// Disallow copy and assignment by default
		FourierInserter3DMode2Maker( const FourierInserter3DMode2Maker& );
		FourierInserter3DMode2Maker& operator=( const FourierInserter3DMode2Maker& );
	};
	
	/** FourierInserter3DMode3Maker  - returns a newly made (pointer to) a FourierInserter3DMode3
	 * See comments in FourierPixelInserterMaker for explanations of functionality
	 */
	class FourierInserter3DMode3Maker : public FourierPixelInserterMaker
	{
	public:
		FourierInserter3DMode3Maker(const unsigned int key = FourierInserter3DMode3::get_mode_number() ): FourierPixelInserterMaker(key) {}
		virtual ~FourierInserter3DMode3Maker() {}
		
		virtual FourierPixelInserter3D* make_inserter(float * const normalize_values, float * const real_data, const unsigned int xsize, const unsigned int ysize, const unsigned int zsize)
		{
			return new FourierInserter3DMode3(normalize_values, real_data, xsize, ysize, zsize);
		}
	private:
		static const FourierInserter3DMode3Maker registerMode3Maker;
		
		// Disallow copy and assignment by default
		FourierInserter3DMode3Maker( const FourierInserter3DMode3Maker& );
		FourierInserter3DMode3Maker& operator=( const FourierInserter3DMode3Maker& );
	};
	
	/** FourierInserter3DMode4Maker  - returns a newly made (pointer to) a FourierInserter3DMode4
	 * See comments in FourierPixelInserterMaker for explanations of functionality
	 */
	class FourierInserter3DMode4Maker : public FourierPixelInserterMaker
	{
	public:
		FourierInserter3DMode4Maker(const unsigned int key = FourierInserter3DMode4::get_mode_number() ): FourierPixelInserterMaker(key) {}
		virtual ~FourierInserter3DMode4Maker() {}
		
		virtual FourierPixelInserter3D* make_inserter(float * const normalize_values, float * const real_data, const unsigned int xsize, const unsigned int ysize, const unsigned int zsize)
		{
			return new FourierInserter3DMode4(normalize_values, real_data, xsize, ysize, zsize);
		}
	private:
		static const FourierInserter3DMode4Maker registerMode4Maker;
		
		// Disallow copy and assignment by default
		FourierInserter3DMode4Maker( const FourierInserter3DMode4Maker& );
		FourierInserter3DMode4Maker& operator=( const FourierInserter3DMode4Maker& );
	};
	
	/** FourierInserter3DMode5Maker  - returns a newly made (pointer to) a FourierInserter3DMode5
	 * See comments in FourierPixelInserterMaker for explanations of functionality
	 */
	class FourierInserter3DMode5Maker : public FourierPixelInserterMaker
	{
	public:
		FourierInserter3DMode5Maker(const unsigned int key = FourierInserter3DMode5::get_mode_number() ): FourierPixelInserterMaker(key) {}
		virtual ~FourierInserter3DMode5Maker() {}
		
		virtual FourierPixelInserter3D* make_inserter(float * const normalize_values, float * const real_data, const unsigned int xsize, const unsigned int ysize, const unsigned int zsize)
		{
			return new FourierInserter3DMode5(normalize_values, real_data, xsize, ysize, zsize);
		}
	private:
		static const FourierInserter3DMode5Maker registerMode5Maker;
		
		// Disallow copy and assignment by default
		FourierInserter3DMode5Maker( const FourierInserter3DMode5Maker& );
		FourierInserter3DMode5Maker& operator=( const FourierInserter3DMode5Maker& );
	};
	
	/** FourierInserter3DMode6Maker  - returns a newly made (pointer to) a FourierInserter3DMode6
	 * See comments in FourierPixelInserterMaker for explanations of functionality
	 */
	class FourierInserter3DMode6Maker : public FourierPixelInserterMaker
	{
	public:
		FourierInserter3DMode6Maker(const unsigned int key = FourierInserter3DMode6::get_mode_number() ): FourierPixelInserterMaker(key) {}
		virtual ~FourierInserter3DMode6Maker() {}
		
		virtual FourierPixelInserter3D* make_inserter(float * const normalize_values, float * const real_data, const unsigned int xsize, const unsigned int ysize, const unsigned int zsize)
		{
			return new FourierInserter3DMode6(normalize_values, real_data, xsize, ysize, zsize);
		}
	private:
		static const FourierInserter3DMode6Maker registerMode6Maker;
		
		// Disallow copy and assignment by default
		FourierInserter3DMode6Maker( const FourierInserter3DMode6Maker& );
		FourierInserter3DMode6Maker& operator=( const FourierInserter3DMode6Maker& );
	};
	
	/** FourierInserter3DMode7Maker  - returns a newly made (pointer to) a FourierInserter3DMode7
	 * See comments in FourierPixelInserterMaker for explanations of functionality
	 */
	class FourierInserter3DMode7Maker : public FourierPixelInserterMaker
	{
	public:
		FourierInserter3DMode7Maker(const unsigned int key = FourierInserter3DMode7::get_mode_number() ): FourierPixelInserterMaker(key) {}
		virtual ~FourierInserter3DMode7Maker() {}
		
		virtual FourierPixelInserter3D* make_inserter(float * const normalize_values, float * const real_data, const unsigned int xsize, const unsigned int ysize, const unsigned int zsize)
		{
			return new FourierInserter3DMode7(normalize_values, real_data, xsize, ysize, zsize);
		}
	private:
		static const FourierInserter3DMode7Maker registerMode7Maker;
		
		// Disallow copy and assignment by default
		FourierInserter3DMode7Maker( const FourierInserter3DMode7Maker& );
		FourierInserter3DMode7Maker& operator=( const FourierInserter3DMode7Maker& );
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
		int m_xsize;
		int m_ysize;
		int m_zsize;
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
