/**
 * $Id$
 */

/*
 * Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
 * Copyright (c) 2000-2006 The University of Texas - Houston Medical School
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
 */

#ifndef eman_processor_sparx_h__
#define eman_processor_sparx_h__ 1

#include "emdata.h"

namespace EMAN
{


	/** mirror an image
	 * @param axis  ''x', 'y', or 'z' axis, means mirror by changing the sign of the respective axis;
	 */
	class MirrorProcessor:public Processor
	{

	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new MirrorProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("axis", EMObject::STRING, "'x', 'y', or 'z' axis");
			return d;
		}

		string get_desc() const
		{
			return "Mirrors an image along the specified axis. This will shift the image center for even box sizes. Use the 'xform.flip' processor to preserve center.";
		}
		
		static const string NAME;
	};


	/** Base class for Fourier processors
	 *@param sigma Gaussian sigma (0-.5)
	 *@param cutoff_abs Processor radius in terms of Nyquist (0-.5)
	 *@param cutoff_pixels Width in Fourier pixels (0 - size()/2
	 *@param cutoff_freq Resolution in 1/A (0 - 1 / size*apix)
	 *@param apix Override A/pix in the image header (changes x,y and z)
	 */
	class NewFourierProcessor:public Processor
	{
	  public:
		//virtual void process_inplace(EMData * image);

		static string get_group_desc()
		{
			return "Fourier Filter Processors are frequency domain processors. The input image can be either real or Fourier, and the output processed image format corresponds to that of the input file. FourierFilter class is the base class of fourier space processors. The processors can be either low-pass, high-pass, band-pass, or homomorphic. The processor parameters are in absolute frequency units, valid range is ]0,0.5], where 0.5 is Nyquist freqeuncy. ";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
//			d.put("cutoff_frequency", EMObject::FLOAT, "Absolute [0,0.5] cut-off frequency.");	//use cutoff_abs
			d.put("sigma", EMObject::FLOAT, "Gaussian sigma (0-.5)");							//use cutoff_abs
			d.put("cutoff_abs", EMObject::FLOAT, "Processor radius in terms of Nyquist (0-.5)");
			d.put("cutoff_pixels", EMObject::FLOAT, "Width in Fourier pixels (0 - size()/2");
			d.put("cutoff_freq", EMObject::FLOAT, "Resolution in 1/A (0 - 1 / size*apix)");
			d.put("apix", EMObject::FLOAT, " Override A/pix in the image header (changes x,y and z)");
			return d;
		}

	  protected:
		virtual void preprocess(EMData * image) {
			if(params.has_key("apix")) {
				image->set_attr("apix_x", (float)params["apix"]);
				image->set_attr("apix_y", (float)params["apix"]);
				image->set_attr("apix_z", (float)params["apix"]);
			}

			const Dict dict = image->get_attr_dict();
			if( params.has_key("sigma")) {
				params["cutoff_abs"] = (float)params["sigma"];
			}
			else if( params.has_key("cutoff_abs") ) {
				params["sigma"] = (float)params["cutoff_abs"];
			}
			else if( params.has_key("cutoff_freq") ) {
				float val =  (float)params["cutoff_freq"] * (float)dict["apix_x"]; 
				params["cutoff_abs"] = val;
				params["sigma"] = val;
			}
			else if( params.has_key("cutoff_pixels") ) {
				float val = (0.5f*(float)params["cutoff_pixels"] / (float)dict["nx"]);
				params["cutoff_abs"] = val;
				params["sigma"] = val;
			}
		}
		virtual void preprocessandconvertpars(EMData * image)
		{
			if(params.has_key("apix")) {
				image->set_attr("apix_x", (float)params["apix"]);
				image->set_attr("apix_y", (float)params["apix"]);
				image->set_attr("apix_z", (float)params["apix"]);
			}

			const Dict dict = image->get_attr_dict();
			if( params.has_key("sigma")) {
				params["cutoff_abs"] = (float)params["sigma"];
			}
			else if( params.has_key("cutoff_abs") ) {
			        // Here I have added a patch 1/sqrt(2) to compensate for the different Gaussian used in EMAN1 vs EMAN2 (John Flanagan)
				float val = (float)params["cutoff_abs"] / sqrt(2.0f);
				params["cutoff_abs"] = val;
				params["sigma"] = val;
				
			}
			else if( params.has_key("cutoff_freq") ) {
			        // Here I have added a patch 1/sqrt(2) to compensate for the different Gaussian used in EMAN1 vs EMAN2 (John Flanagan)
				float val =  (float)params["cutoff_freq"] * (float)dict["apix_x"] / sqrt(2.0f); 
				params["cutoff_abs"] = val;
				params["sigma"] = val;
			}
			else if( params.has_key("cutoff_pixels") ) {
			        // Here I have added a patch 1/sqrt(2) to compensate for the different Gaussian used in EMAN1 vs EMAN2 (John Flanagan)
				float val = (0.5f*(float)params["cutoff_pixels"] / (float)dict["nx"]) / sqrt(2.0f);
				params["cutoff_abs"] = val;
				params["sigma"] = val;
			}
		}
		virtual void setbutterworthdefaults(EMData * image)
		{
		        float highlowratio = 0.15f;
			const Dict dict = image->get_attr_dict();
			
			if(params.has_key("cutoff_abs") && !params.has_key("low_cutoff_frequency"))
			{
			        params["low_cutoff_frequency"] = (float)params["cutoff_abs"];
				
			        float val = (float)params["low_cutoff_frequency"];
			        params["high_cutoff_frequency"] = highlowratio*log10(1.0f/val) + val;
			}
			
                        else if(params.has_key("cutoff_freq") && !params.has_key("low_cutoff_frequency"))
			{
			        params["low_cutoff_frequency"] = (float)params["cutoff_freq"] * (float)dict["apix_x"];
				
			        float val = (float)params["low_cutoff_frequency"];
			        params["high_cutoff_frequency"] = highlowratio*log10(1.0f/val) + val;  
			}
			
			else if(params.has_key("cutoff_pixels") && !params.has_key("low_cutoff_frequency"))
			{
			        params["low_cutoff_frequency"] = (0.5f*(float)params["cutoff_pixels"] / (float)dict["nx"]);
				
			        float val = (float)params["low_cutoff_frequency"];
			        params["high_cutoff_frequency"] = highlowratio*log10(1.0f/val) + val;  
			}
			
		}
	};

	/**Lowpass top-hat filter processor applied in Fourier space.
	 * @param sigma Absolute [0,0.5] cut-off frequency.
	 */
	class NewLowpassTopHatProcessor:public NewFourierProcessor
	{
	  public:
		string get_name() const
		{ return NAME; }
		static Processor *NEW()
		{ return new NewLowpassTopHatProcessor(); }
		string get_desc() const
		{
			return "Lowpass top-hat filter processor applied in Fourier space.";
		}
		void process_inplace(EMData* image) {
			params["filter_type"] = TOP_HAT_LOW_PASS;
			preprocess(image);
			EMFourierFilterInPlace(image, params);
		}
		
		static const string NAME;
	};

	/** Highpass top-hat filter applied in Fourier NewLowpassGaussProcessorspace.
	 * @param sigma Absolute [0,0.5] cut-off frequency.
	 */
	class NewHighpassTopHatProcessor:public NewFourierProcessor
	{
	  public:
		string get_name() const
		{ return NAME; }
		static Processor *NEW()
		{ return new NewHighpassTopHatProcessor(); }
		string get_desc() const
		{
			return "Highpass top-hat filter applied in Fourier space.";
		}
		void process_inplace(EMData* image) {
			params["filter_type"] = TOP_HAT_HIGH_PASS;
			preprocess(image);
			EMFourierFilterInPlace(image, params);
		}
		
		static const string NAME;
	};

	/**Bandpass top-hat filter processor applied in Fourier space.
	 *@param low_cutoff_frequency Absolute [0,0.5] low cut-off frequency.
	 *@param high_cutoff_frequency Absolute [0,0.5] high cut-off frequency.
	 */
	class NewBandpassTopHatProcessor:public NewFourierProcessor
	{
	  public:
		string get_name() const
		{ return NAME; }
		static Processor *NEW()
		{ return new NewBandpassTopHatProcessor(); }
		string get_desc() const
		{
			return "Bandpass top-hat filter processor applied in Fourier space.";
		}
		void process_inplace(EMData* image) {
			params["filter_type"] = TOP_HAT_BAND_PASS;
			EMFourierFilterInPlace(image, params);
		}
		TypeDict get_param_types() const
		{
			TypeDict d = NewFourierProcessor::get_param_types();
			d.put("low_cutoff_frequency", EMObject::FLOAT, "Absolute [0,0.5] low cut-off frequency.");
			d.put("high_cutoff_frequency", EMObject::FLOAT, "Absolute [0,0.5] high cut-off frequency.");
			return d;
		}
		
		static const string NAME;
	};

	/**Homomorphic top-hat filter processor applied in Fourier space.
	 *@param low_cutoff_frequency Absolute [0,0.5] low cut-off frequency.
	 *@param high_cutoff_frequency Absolute [0,0.5] high cut-off frequency.
	 *@param value_at_zero_frequency Value at zero frequency.
	 */
	class NewHomomorphicTopHatProcessor:public NewFourierProcessor
	{
	  public:
		string get_name() const
		{ return NAME; }
		static Processor *NEW()
		{ return new NewHomomorphicTopHatProcessor(); }
		string get_desc() const
		{
			return "Homomorphic top-hat filter processor applied in Fourier space.";
		}
		void process_inplace(EMData* image) {
			params["filter_type"] = TOP_HOMOMORPHIC;
			EMFourierFilterInPlace(image, params);
		}
		TypeDict get_param_types() const
		{
			TypeDict d = NewFourierProcessor::get_param_types();
			d.put("low_cutoff_frequency", EMObject::FLOAT, "Absolute [0,0.5] low cut-off frequency.");
			d.put("high_cutoff_frequency", EMObject::FLOAT, "Absolute [0,0.5] high cut-off frequency.");
			d.put("value_at_zero_frequency", EMObject::FLOAT, "Value at zero frequency.");
			return d;
		}
		
		static const string NAME;
	};

	/**Lowpass Gauss filter processor applied in Fourier space.
	 * @param cutoff_abs Gaussian sigma.
	 */
	class NewLowpassGaussProcessor:public NewFourierProcessor
	{
	  public:
		string get_name() const
		{ return NAME; }
		static Processor *NEW()
		{ return new NewLowpassGaussProcessor(); }
		string get_desc() const
		{
			return "Lowpass Gauss filter processor applied in Fourier space.";
		}
		void process_inplace(EMData* image) {
			params["filter_type"] = GAUSS_LOW_PASS;
			preprocessandconvertpars(image);
			
			if(params.has_key("cutoff_resolv")){
			  
			        const Dict dict = image->get_attr_dict();
			
			        // here I have added a little function to filter to a resolvability (John Flanagan 20/09/2010)
				float R = 1/((float)params["cutoff_resolv"]*(float)dict["apix_x"]);    // convert to pixels
				float rsigma = sqrt((-4*log(0.36f))/(pow(M_PI,2)*pow(R,2)));        // find optimal sigma
				params["cutoff_abs"] = rsigma / sqrt(2.0f);                               // patch to fix the 2s^2 problem
				params["sigma"] = rsigma / sqrt(2.0f);                                    // patch to fix the 2s^2 problem
			}
			
			EMFourierFilterInPlace(image, params);
		}
		TypeDict get_param_types() const
		{
			TypeDict d = NewFourierProcessor::get_param_types();
			d.put("cutoff_resolv", EMObject::FLOAT, "Resolvibility in 1/A, applied using filter.lowpass.gauss where cutoff_freq = sqrt(-4ln(0.36)/pi^2R^2) & R = 1/cf*apix");
			return d;
		}
		
		static const string NAME;
	};

	/**Highpass Gauss filter processor applied in Fourier space.
	 * @param cutoff_abs Gaussian sigma.
	 */
	class NewHighpassGaussProcessor:public NewFourierProcessor
	{
	  public:
		string get_name() const
		{ return NAME; }
		static Processor *NEW()
		{ return new NewHighpassGaussProcessor(); }
		string get_desc() const
		{
			return "Highpass Gauss filter processor applied in Fourier space.";
		}
		void process_inplace(EMData* image) {
			params["filter_type"] = GAUSS_HIGH_PASS;
			preprocessandconvertpars(image);
			EMFourierFilterInPlace(image, params);
		}
		
		static const string NAME;
	};

	/**Bandpass Gauss filter processor applied in Fourier space.
	 * @param cutoff_abs Gaussian sigma.
	 * @param center Gaussian center.
	 */
	class NewBandpassGaussProcessor:public NewFourierProcessor
	{
	  public:
		string get_name() const
		{ return NAME; }
		static Processor *NEW()
		{ return new NewBandpassGaussProcessor(); }
		string get_desc() const
		{
			return "Bandpass Gauss filter processor applied in Fourier space.";
		}
		void process_inplace(EMData* image) {
			params["filter_type"] = GAUSS_BAND_PASS;
			preprocess(image);
			EMFourierFilterInPlace(image, params);
		}
		TypeDict get_param_types() const
		{
			TypeDict d = NewFourierProcessor::get_param_types();
			d.put("center", EMObject::FLOAT, "Gaussian center.");
			return d;
		}
		
		static const string NAME;
	};

	/**Homomorphic Gauss filter processor applied in Fourier space.
	 * @param cutoff_abs Gaussian sigma.
	 * @param value_at_zero_frequency Value at zero frequency.
	 */
	class NewHomomorphicGaussProcessor:public NewFourierProcessor
	{
	  public:
		string get_name() const
		{ return NAME; }
		static Processor *NEW()
		{ return new NewHomomorphicGaussProcessor(); }
		string get_desc() const
		{
			return "Homomorphic Gauss filter processor applied in Fourier space.";
		}
		void process_inplace(EMData* image) {
			params["filter_type"] = GAUSS_HOMOMORPHIC;
			preprocess(image);
			EMFourierFilterInPlace(image, params);
		}
		TypeDict get_param_types() const
		{
			TypeDict d = NewFourierProcessor::get_param_types();
			d.put("value_at_zero_frequency", EMObject::FLOAT, "Value at zero frequency.");
			return d;
		}
		
		static const string NAME;
	};

	/**Divide by a Gaussian in Fourier space.
	 * @param cutoff_abs Gaussian sigma.
	 */
	class NewInverseGaussProcessor:public NewFourierProcessor
	{
	  public:
		string get_name() const
		{ return NAME; }
		static Processor *NEW()
		{ return new NewInverseGaussProcessor(); }
		string get_desc() const
		{
			return "Divide by a Gaussian in Fourier space.";
		}
		void process_inplace(EMData* image) {
			params["filter_type"] = GAUSS_INVERSE;
			preprocess(image);
			EMFourierFilterInPlace(image, params);
		}
		
		static const string NAME;
	};

	/**Shift by phase multiplication in Fourier space.
	 */
	class SHIFTProcessor:public NewFourierProcessor
	{
	  public:
		string get_name() const
		{ return NAME; }
		static Processor *NEW()
		{ return new SHIFTProcessor(); }
		string get_desc() const
		{
			return "Shift by phase multiplication in Fourier space.";
		}
		void process_inplace(EMData* image) {
			params["filter_type"] = SHIFT;
			EMFourierFilterInPlace(image, params);
		}
		TypeDict get_param_types() const
		{
			TypeDict d = NewFourierProcessor::get_param_types();
			d.put("x_shift", EMObject::FLOAT, "Shift x");
			d.put("y_shift", EMObject::FLOAT, "Shift y");
			d.put("z_shift", EMObject::FLOAT, "Shift z");
			return d;
		}
		
		static const string NAME;
	};

	/**Divide by a Kaiser-Bessel I0 func in Fourier space.
	 */
	class InverseKaiserI0Processor:public NewFourierProcessor
	{
	  public:
		string get_name() const
		{ return NAME; }
		static Processor *NEW()
		{ return new InverseKaiserI0Processor(); }
		string get_desc() const
		{
			return "Divide by a Kaiser-Bessel I0 func in Fourier space.";
		}
		void process_inplace(EMData* image) {
			params["filter_type"] = KAISER_I0_INVERSE;
			EMFourierFilterInPlace(image, params);
		}
		TypeDict get_param_types() const
		{
			TypeDict d = NewFourierProcessor::get_param_types();
			return d;
		}
		
		static const string NAME;
	};

	/**Divide by a Kaiser-Bessel Sinh func in Fourier space.
	 */
	class InverseKaiserSinhProcessor:public NewFourierProcessor
	{
	  public:
		string get_name() const
		{ return NAME; }
		static Processor *NEW()
		{ return new InverseKaiserSinhProcessor(); }
		string get_desc() const
		{
			return "Divide by a Kaiser-Bessel Sinh func in Fourier space.";
		}
		void process_inplace(EMData* image) {
			params["filter_type"] = KAISER_SINH_INVERSE;
			EMFourierFilterInPlace(image, params);
		}
		TypeDict get_param_types() const
		{
			TypeDict d = NewFourierProcessor::get_param_types();
			return d;
		}
		
		static const string NAME;
	};

	/**Filter with tabulated data in Fourier space.
	 *@param table Tabulated filter data.
	 */
	class NewRadialTableProcessor:public NewFourierProcessor
	{
	  public:
		string get_name() const
		{ return NAME; }
		
		static Processor *NEW()
		{ return new NewRadialTableProcessor(); }
		
		string get_desc() const
		{
			return "Filter with tabulated data in EMFourierFilterFuncFourier space. 1 value per Fourier pixel, extending to corners. Missing value assumed to be 0.";
		}
		void process_inplace(EMData* image) {
			params["filter_type"] = RADIAL_TABLE;
			EMFourierFilterInPlace(image, params);
		}
		TypeDict get_param_types() const
		{
//			TypeDict d = NewFourierProcessor::get_param_types();
			TypeDict d;
			d.put("table", EMObject::FLOATARRAY, "Radial data array. 1 value per Fourier image pixel.");
			return d;
		}
		
		static const string NAME;
	};

	/**Lowpass Butterworth filter processor applied in Fourier space.
	 *@param low_cutoff_frequency Absolute [0,0.5] low cut-off frequency.
	 *@param high_cutoff_frequency Absolute [0,0.5] high cut-off frequency.
	 */
	class NewLowpassButterworthProcessor:public NewFourierProcessor
	{
	  public:
		string get_name() const
		{ return NAME; }
		static Processor *NEW()
		{ return new NewLowpassButterworthProcessor(); }
		string get_desc() const
		{
			return "Lowpass Butterworth filter processor applied in Fourier space.";
		}
		void process_inplace(EMData* image) {
		        setbutterworthdefaults(image);
			params["filter_type"] = BUTTERWORTH_LOW_PASS;
			EMFourierFilterInPlace(image, params);
		}
		TypeDict get_param_types() const
		{
			TypeDict d = NewFourierProcessor::get_param_types();
			d.put("low_cutoff_frequency", EMObject::FLOAT, "Absolute [0,0.5] low cut-off frequency.");
			d.put("high_cutoff_frequency", EMObject::FLOAT, "Absolute [0,0.5] high cut-off frequency.");
			return d;
		}
		
		static const string NAME;
	};

	/**Highpass Butterworth filter processor applied in Fourier space.
	 *@param low_cutoff_frequency Absolute [0,0.5] low cut-off frequency.
	 *@param high_cutoff_frequency Absolute [0,0.5] high cut-off frequency.
	 */
	class NewHighpassButterworthProcessor:public NewFourierProcessor
	{
	  public:
		string get_name() const
		{ return NAME; }
		static Processor *NEW()
		{ return new NewHighpassButterworthProcessor(); }
		string get_desc() const
		{
			return "Highpass Butterworth filter processor applied in Fourier space.";
		}
		void process_inplace(EMData* image) {
			params["filter_type"] = BUTTERWORTH_HIGH_PASS;
			setbutterworthdefaults(image);
			EMFourierFilterInPlace(image, params);
		}
		TypeDict get_param_types() const
		{
			TypeDict d = NewFourierProcessor::get_param_types();
			d.put("low_cutoff_frequency", EMObject::FLOAT, "Absolute [0,0.5] low cut-off frequency.");
			d.put("high_cutoff_frequency", EMObject::FLOAT, "Absolute [0,0.5] high cut-off frequency.");
			return d;
		}
		
		static const string NAME;
	};

	/**Homomorphic Butterworth filter processor applied in Fourier space.
	 *@param low_cutoff_frequency Absolute [0,0.5] low cut-off frequency.
	 *@param high_cutoff_frequency Absolute [0,0.5] high cut-off frequency.
	 *@param value_at_zero_frequency Value at zero frequency.
	 */
	class NewHomomorphicButterworthProcessor:public NewFourierProcessor
	{
	  public:
		string get_name() const
		{ return NAME; }
		static Processor *NEW()
		{ return new NewHomomorphicButterworthProcessor(); }
		string get_desc() const
		{
			return "Homomorphic Butterworth filter processor applied in Fourier space.";
		}
		void process_inplace(EMData* image) {
			params["filter_type"] = BUTTERWORTH_HOMOMORPHIC;
			EMFourierFilterInPlace(image, params);
		}
		TypeDict get_param_types() const
		{
			TypeDict d = NewFourierProcessor::get_param_types();
			d.put("low_cutoff_frequency", EMObject::FLOAT, "Absolute [0,0.5] low cut-off frequency.");
			d.put("high_cutoff_frequency", EMObject::FLOAT, "Absolute [0,0.5] high cut-off frequency.");
			d.put("value_at_zero_frequency", EMObject::FLOAT, "Value at zero frequency.");
			return d;
		}
		
		static const string NAME;
	};

	/**Lowpass tanh filter processor applied in Fourier space.
	 *@param cutoff_abs Absolute [0,0.5] cut-off frequency.
	 *@param fall_off Tanh decay rate.
	 */
	class NewLowpassTanhProcessor:public NewFourierProcessor
	{
	  public:
		string get_name() const
		{ return NAME; }
		static Processor *NEW()
		{ return new NewLowpassTanhProcessor(); }
		string get_desc() const
		{
			return "Lowpass tanh filter processor applied in Fourier space.";
		}
		void process_inplace(EMData* image) {
			params["filter_type"] = TANH_LOW_PASS;
			preprocess(image);
			params.set_default("fall_off",.5f); // Only sets it if is not already set
			EMFourierFilterInPlace(image, params);
		}
		TypeDict get_param_types() const
		{
			TypeDict d = NewFourierProcessor::get_param_types();
			d.put("fall_off", EMObject::FLOAT, "Tanh decay rate.");
			return d;
		}
		
		static const string NAME;
	};

	/**Highpass tanh filter processor applied in Fourier space.
	 *@param cutoff_abs Absolute [0,0.5] cut-off frequency.
	 *@param fall_off Tanh decay rate.
	 */
	class NewHighpassTanhProcessor:public NewFourierProcessor
	{
	  public:
		string get_name() const
		{ return NAME; }
		static Processor *NEW()
		{ return new NewHighpassTanhProcessor(); }
		string get_desc() const
		{
			return "Highpass tanh filter processor applied in Fourier space.";
		}
		void process_inplace(EMData* image) {
			params["filter_type"] = TANH_HIGH_PASS;
			params.set_default("fall_off",.5f); // Only sets it if is not already set
			preprocess(image);
			EMFourierFilterInPlace(image, params);
		}
		TypeDict get_param_types() const
		{
			TypeDict d = NewFourierProcessor::get_param_types();
			d.put("fall_off", EMObject::FLOAT, "Tanh decay rate.");
			return d;
		}
		
		static const string NAME;
	};

	/**Homomorphic Tanh processor applied in Fourier space
	 *@param cutoff_abs Absolute [0,0.5] cut-off frequency.
	 *@param fall_off Tanh decay rate.
	 *@param value_at_zero_frequency Value at zero frequency.
	 */
	class NewHomomorphicTanhProcessor:public NewFourierProcessor
	{
	  public:
		string get_name() const
		{ return NAME; }
		static Processor *NEW()
		{ return new NewHomomorphicTanhProcessor(); }
		string get_desc() const
		{
			return "Homomorphic Tanh processor applied in Fourier space";
		}
		void process_inplace(EMData* image) {
			params["filter_type"] = TANH_HOMOMORPHIC;
			params.set_default("fall_off",.5f); // Only sets it if is not already set
			preprocess(image);
			EMFourierFilterInPlace(image, params);
		}
		TypeDict get_param_types() const
		{
			TypeDict d = NewFourierProcessor::get_param_types();
			d.put("fall_off", EMObject::FLOAT, "Tanh decay rate.");
			d.put("value_at_zero_frequency", EMObject::FLOAT, "Value at zero frequency.");
			return d;
		}
		
		static const string NAME;
	};

	/**Bandpass tanh processor applied in Fourier space.
	 *@param low_cutoff_frequency Absolute [0,0.5] low cut-off frequency.
	 *@param Low_fall_off Tanh low decay rate.
	 *@param high_cutoff_frequency Absolute [0,0.5] high cut-off frequency.
	 *@param high_fall_off Tanh high decay rate.
	 *@param fall_off Tanh decay rate.
	 */
	class NewBandpassTanhProcessor:public NewFourierProcessor
	{
	  public:
		string get_name() const
		{ return NAME; }
		static Processor *NEW()
		{ return new NewBandpassTanhProcessor(); }
		string get_desc() const
		{
			return "Bandpass tanh processor applied in Fourier space.";
		}
		void process_inplace(EMData* image) {
			params["filter_type"] = TANH_BAND_PASS;
			EMFourierFilterInPlace(image, params);
		}
		TypeDict get_param_types() const
		{
			TypeDict d = NewFourierProcessor::get_param_types();
			d.put("low_cutoff_frequency", EMObject::FLOAT, "Absolute [0,0.5] low cut-off frequency.");
			d.put("Low_fall_off", EMObject::FLOAT, "Tanh low decay rate.");
			d.put("high_cutoff_frequency", EMObject::FLOAT, "Absolute [0,0.5] high cut-off frequency.");
			d.put("high_fall_off", EMObject::FLOAT, "Tanh high decay rate.");
			d.put("fall_off", EMObject::FLOAT, "Tanh decay rate.");
			return d;
		}
		
		static const string NAME;
	};

	class CTF_Processor:public NewFourierProcessor
	{
	  public:
		string get_name() const
		{ return NAME; }
		static Processor *NEW()
		{ return new CTF_Processor(); }
		string get_desc() const
		{
			return "CTF_ is applied in Fourier image.";
		}
		void process_inplace(EMData* image) {
			params["filter_type"] = CTF_;
			EMFourierFilterInPlace(image, params);
		}
		TypeDict get_param_types() const
		{
			TypeDict d = NewFourierProcessor::get_param_types();
			d.put("defocus", EMObject::FLOAT, "defocus value in Angstrom.");
			d.put("cs", EMObject::FLOAT, "cs in CM.");
			d.put("voltage", EMObject::FLOAT, "voltage in Kv.");
			d.put("ps", EMObject::FLOAT, "pixel size.");
			d.put("b_factor", EMObject::FLOAT, "Gaussian like evelope function (b_factor).");
			d.put("wgh", EMObject::FLOAT, "Amplitude contrast ratio.");
			d.put("sign", EMObject::FLOAT, "Sign of Contrast transfer function,and use -1 to compensate.");
			d.put("npow", EMObject::FLOAT, "power of transfer function.");
			return d;
		}
		
		static const string NAME;
	};
}

#endif	//eman_processor_sparx_h__
