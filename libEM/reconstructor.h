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

	/** A Reconstructor containing an EMData pointer and x,y,z dimensions
	 *  It basically stores a pointer to an image object and stores the dimensions of the image volume.
	 *  This class was originally added simply to encapsulate the 
	 *  the things common to FourierReconstructor, WienerFourierReconstructor
	 *  and BackProjectionReconstructor. It was never expected to instantiated on its own,
	 *  and is intended to be a parent of the Reconstructor class.
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
		FourierReconstructor() : inserter(0) { load_default_settings(); }
		virtual ~FourierReconstructor() { free_memory(); }
	
		/** Copy constructor
		 */
		FourierReconstructor( const FourierReconstructor& that ) : Reconstructor(that)
		{
			load_inserter();
		}

		/** Assignment operator
		 */
		FourierReconstructor& operator=( const FourierReconstructor& );

		/** Assignment operator
		 * @exception InvalidValueException
		 */
		virtual void setup();

		virtual int insert_slice(const EMData* const slice, const Transform3D & euler);
		// these functions are for testing purposes
		bool do_insert_remove_test(const EMData* const input_slice, const Transform3D & arg);
		int remove_slice(const EMData* const input_slice, const Transform3D & arg);
		bool test_pixel_wise_zero(const EMData* const input_slice, const Transform3D & arg);
		
		/** 
	  	 * @return 
	  	 * @param input_slice
	  	 * @param arg
	  	 * @param num_particles_in_slice
	  	 * @exception 
		 */
		virtual int determine_slice_agreement(const EMData* const input_slice, const Transform3D & arg, const unsigned int  num_particles_in_slice = 1);

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
		
		virtual void iteration_reset()
		{
			static bool first_time = true;
			if ( !first_time )
			{
				//cout << "Normalizing" << endl;
// 				float* norm = tmp_data->get_data();
// 				float* rdata = image->get_data();
// 				for (int i = 0; i < nx * ny * nz; i += 2) {
// 					float d = norm[i];
// 					if (d == 0) {
// 						rdata[i] = 0;
// 						rdata[i + 1] = 0;
// 					}
// 					else {
// 					rdata[i] /= d;
// 					rdata[i + 1] /= d;
// 					}
// 				}	
			}
			//else cout << "Not normalizing" << endl;
			
			first_time = false;
			quality_scores.clear(); idx = 0;
		}

		virtual float get_score(const unsigned int idx) { if ( quality_scores.size() > idx ) return quality_scores[idx].get_frc_integral(); else {cout << "foo" << endl; return 2;}}
		virtual float get_norm(const unsigned int idx) { if ( quality_scores.size() > idx ) return quality_scores[idx].get_norm(); else {cout << "foo" << endl; return 2;}}
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

		FourierPixelInserter3D* inserter;

		void free_memory()
		{
			if ( inserter != 0 )
			{
				delete inserter;
				inserter = 0;
			}
		}
	
		void load_inserter()
		{
			if ( inserter != 0 )
			{
				delete inserter;
				inserter = 0;
			}
			
			string mode = (string)params["mode"];
			Dict parms;
			parms["rdata"] = image->get_data();
			parms["norm"] = tmp_data->get_data();
			parms["nx"] = nx;
			parms["ny"] = ny;
			parms["nz"] = nz;
			
			inserter = Factory<FourierPixelInserter3D>::get(mode, parms);
			inserter->init();
		}

		void do_insert_slice_work(const EMData* const input_slice, const Transform3D & arg);
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
