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

#ifndef eman_processor_h__
#define eman_processor_h__ 1

#include "emobject.h"
#include "util.h"
#include "geometry.h"

#include <cfloat>
#include <climits>
#include <cstring>

using std::vector;
using std::map;
using std::string;

namespace EMAN
{
	class EMData;
	class Transform3D;

	/** Typical usage of Processors are as follows:
     *
     *   - How to get all the processor names
     *@code
     *      vector<string> all_processors = Factory<Processor>::get_list();
     @endcode
     *   - How to use a processor
     *@code
     *      EMData *img = ...;
     *      img->process_inplace("PROCESSORNAME", Dict("sigma", 12));
     @endcode
     *   - How to define a new XYZProcessor \n
     *      XYZProcessor should either extend the base class 'Processor' or a
     *      subclass of 'Processor'. At a minimum, it should define:
	 *      (Please replace 'XYZ' with your own class name).
	 *@code
     *          string get_name() const { return "processorname"; }
     *          static Processor *NEW() { return XYZProcessor(); }
	 @endcode
	 *      If XYZProcessor is a parent class, it should define:
	 *@code
	 *          static string get_group_desc();
	 @endcode
	 *      Otherwise, it should define:
	 *@code
	 *          string get_desc() const;
	 @endcode
     *      If XYZProcessor need parameters not defined by its parent
     *      class, it should define:
	 *@code
	 *          Dict get_params() const;
	 *          void set_params(const Dict & new_params);
     *          TypeDict get_param_types() const;
     @endcode
     */
	class Processor
	{
	  public:
		virtual ~ Processor()
		{
		}

		/** To process an image in-place. 
		 * For those processors which can only be processed out-of-place, override this function 
		 * to just print out some error message to remind user call the out-of-place version.
		 * @param image The image to be processed.
		 */
		virtual void process_inplace(EMData *image) = 0;
		
		/** To proccess an image out-of-place. 
		 * For those processors which can only be processed out-of-place, override this function 
		 * to give the right behavior.
		 * @param image The image will be copied, actual process happen on copy of image.
		 * @return the image processing result, may or may not be the same size of the input image  
		 * */
		virtual EMData* process(const EMData * const image);

		/** To process multiple images using the same algorithm.
		 * @param images Multiple images to be processed.
		 */
		virtual void process_list_inplace(vector < EMData * > & images)
		{
			for (size_t i = 0; i < images.size(); i++) {
				process_inplace(images[i]);
			}
		}

		/** Get the processor's name. Each processor is identified by a unique name.
		 * @return The processor's name.
		 */
		virtual string get_name() const = 0;
		
		/** Get the processor parameters in a key/value dictionary.
		 * @return A key/value pair dictionary containing the parameters.
		 */
		virtual Dict get_params() const
		{
			return params;
		}
		
		/** Set the processor parameters using a key/value dictionary.
		 * @param new_params A dictionary containing the new parameters.
		 */
		virtual void set_params(const Dict & new_params)
		{
			params = new_params;
		}
		
		/** Get processor parameter information in a dictionary. Each
		 * parameter has one record in the dictionary. Each record
		 * contains its name, data-type, and description.
		 *
		 * @return A dictionary containing the parameter info.
		 */	 
		virtual TypeDict get_param_types() const
		{
			return TypeDict();
		}

		/** Get the description of this group of processors. This
		 * function is defined in a parent class. It gives a
		 * introduction to a group of processors.
		 *
		 * @return The description of this group of processors.
		 */
		static string get_group_desc()
		{
			return "EMAN processors are in-place image processors. You may apply a processor to process a single image or process multiple images. Processor class is the base class for all processor. <br> \
The basic design of EMAN Processors: <br>\
    1) Each Processor class defines an image-processinging algorithm. <br>\
    2) All the Processor classes in EMAN are managed by a Factory pattern. So each Processor class must define: <br> a) a unique name to idenfity itself in the factory. <br>b) a static method to register itself in the factory.<br>\
    3) Each Processor class defines its own parameter set.<br>\
    4) Each Processor class defines functions to return its documentation including parameter information, and processor description. These functions enable EMAN to generate processor manuals dynamically.";
		}

		/** Get the descrition of this specific processor. This function
		 * must be overwritten by a subclass.
		 *
		 * @return The description of this processor.
		 */
		virtual string get_desc() const = 0;

		/** Fourier filter Processor type enum.
		 *  New Fourier filter processors are computed in a single function,
		 *  EMFourierFilterFunc, that uses a large switch statement to 
		 *  apply the correct filter processor.  This enum specifies the
		 *  filter processor to be applied.
		 */
		enum fourier_filter_types {
			TOP_HAT_LOW_PASS,
			TOP_HAT_HIGH_PASS,
			TOP_HAT_BAND_PASS,
			TOP_HOMOMORPHIC,
			GAUSS_LOW_PASS,
			GAUSS_HIGH_PASS,
			GAUSS_BAND_PASS,
			GAUSS_INVERSE,
			GAUSS_HOMOMORPHIC,
			BUTTERWORTH_LOW_PASS,
			BUTTERWORTH_HIGH_PASS,
			BUTTERWORTH_HOMOMORPHIC,
			KAISER_I0,
			KAISER_SINH,
			KAISER_I0_INVERSE,
			KAISER_SINH_INVERSE,
			SHIFT,
			TANH_LOW_PASS,
			TANH_HIGH_PASS,
			TANH_HOMOMORPHIC,
			TANH_BAND_PASS,	
			RADIAL_TABLE,
		        CTF_,
		};

		/** Compute a Fourier-filter processed image in place.
		 *
		 *  @par Purpose: Apply selected Fourier space processor to 1-,2-, or 3-D image.
		 *  @par Method: 
		 *  
		 *  @param     fimage  Input image object to be processed, either
		 *                     a real-space image or a Fourier-space image.
		 *                     Image may be 1-, 2-, or 3-dimensional.  The
		 *                     original input image is not touched by 
		 *                     this routine.
		 *
		 *  @param[in] params  Processor parameters.  Different processors require
		 *                     different parameters, so we this routine accepts
		 *                     a dictionary of parameters and looks up the
		 *                     appropriate params for the chosen processor at
		 *                     run time.  All processors use the "dopad" 
		 *                     parameter to determine whether the 
		 *                     Fourier workspace array should be zero-
		 *                     padded to twice the original length
		 *                     (dopad == 1) or not zero-padded at all
		 *                     (dopad == 0).
		 *  @return No explicit return.  The image fimage is modified
		 *  in place.
		 */
		static void
		EMFourierFilterInPlace(EMData* fimage, Dict params) {
			bool doInPlace = true;
			EMFourierFilterFunc(fimage, params, doInPlace);
		}
		
		/** Compute a Fourier-processor processed image without altering the original image.
		 *
		 *  @par Purpose: Apply selected Fourier space processor to 1-,2-, or 3-D image.
		 *  @par Method: 
		 *  
		 *  @param     fimage  Input image object to be processeded, either
		 *                     a real-space image or a Fourier-space image.
		 *                     Image may be 1-, 2-, or 3-dimensional.  
		 *
		 *  @param[in] params  Processor parameters.  Different processors require
		 *                     different parameters, so we this routine accepts
		 *                     a dictionary of parameters and looks up the
		 *                     appropriate params for the chosen processor processor at
		 *                     run time.  All processors use the "dopad" 
		 *                     parameter to determine whether the 
		 *                     Fourier workspace array should be zero-
		 *                     padded to twice the original length
		 *                     (dopad == 1) or not zero-padded at all
		 *                     (dopad == 0).
		 *  @return 1-, 2-, or 3-dimensional filter processed image.  If the
		 *          input image is a real-space image, then the returned
		 *          output image will also be a real-space image.  
		 *          Similarly, if the input image is already a Fourier image,
		 *          then the output image will be a Fourier image.
		 */
		static EMData*
		EMFourierFilter(EMData* fimage, Dict params) {
			bool doInPlace = false;
			return EMFourierFilterFunc(fimage, params, doInPlace);
		}
		
	  private:
		/** Compute a Fourier-filter processed image.
		 *  This function is called by either of the convience functions
		 *  EMFourierFilter or EMFourierFilterInPlace.
		 *
		 *  @par Purpose: Apply selected Fourier space processor to 1-,2-, or 3-D image.
		 *  @par Method: 
		 *  
		 *  @param     fimage  Input image object to be processed, either
		 *                     a real-space image or a Fourier-space image.
		 *                     Image may be 1-, 2-, or 3-dimensional.  Image
		 *                     fimage will not be changed unless 
		 *                     inplace == true.
		 *  @param[in] params  Processor parameters.  Different processor processors require
		 *                     different parameters, so we this routine accepts
		 *                     a dictionary of parameters and looks up the
		 *                     appropriate params for the chosen processor processor at
		 *                     run time.  All processors use the "dopad" 
		 *                     parameter to determine whether the 
		 *                     Fourier workspace array should be zero-
		 *                     padded to twice the original length
		 *                     (dopad == 1) or not zero-padded at all
		 *                     (dopad == 0).
		 *  @param[in] doInPlace Inplace flag.  If this flag is true then 
		 *                     fimage will contain the processeded image 
		 *                     when this function returns.
		 *
		 *  @return 1-, 2-, or 3-dimensional filter processed image.  If the
		 *          input image is a real-space image, then the returned
		 *          output image will also be a real-space image.  
		 *          Similarly, if the input image is already a Fourier image,
		 *          then the output image will be a Fourier image.
		 *          In either case, if inplace == true then the output
		 *          image (pointer) will be the same as the input image
		 *          (pointer).
		 */
		static EMData* 
		EMFourierFilterFunc(EMData* fimage, Dict params, bool doInPlace=true);
	
	  protected:
		mutable Dict params;
	};

	class ImageProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);
		
		static string get_group_desc()
		{
			return "An Image Processor defines a way to create a processor image. The processor image is used to multiply the input-image in the fourier space. ImageFilter class is the base class. Each specific ImageFilter class must define function create_processor_image(). ";
		}
		
	  protected:
		virtual EMData * create_processor_image() const = 0;
	};

	/**
	 */
	class FourierProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);
		
		static string get_group_desc()
		{
			return "Fourier Filter processors are a group of processor in the frequency domain. Before using such processors on an image, the image must be transformed from real space to the fourier space. FourierProcessor class is the base class of fourier space processors. Each specific processor is either a lowpass filter processor, or a highpass filter processor, or neighter. The unit of lowpass and highpass parameters are in terms of Nyquist, valid range is [0,0.5]. ";
		}
		
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("cutoff_abs", EMObject::FLOAT, "Processor radius in terms of Nyquist (0-.5)");
			d.put("cutoff_pixels", EMObject::FLOAT, "0 - get_xsize()/2");
			d.put("cutoff_freq", EMObject::FLOAT, "0 - 1.0/(get_xsize()*apix_x) same for y and z");
			return d;
		}
		
	  protected:
		virtual void create_radial_func(vector < float >&radial_mask) const = 0;
	};

	class AmpweightFourierProcessor:public Processor
	{
	  public:
		string get_name() const
		{
			return "filter.ampweight";
		}

		void process_inplace(EMData * image);
		
		void set_params(const Dict & new_params)
		{
			params = new_params;
			sum = params["sum"];
			dosqrt = params["sqrt"];
//			printf("%s %f\n",params.keys()[0].c_str(),lowpass);
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("sum", EMObject::EMDATA, "Adds the weights to sum for normalization");
			d.put("sqrt", EMObject::INT, "Weights using sqrt of the amplitude if set");
			return d;
		}

		static Processor *NEW()
		{
			return new AmpweightFourierProcessor();
		}

		string get_desc() const
		{
			return "Multiplies each Fourier pixel by its amplitude";
		}

		protected:
		EMData *sum;
		int dosqrt;
	};
	/**
	 * This processor performs fast convolution in Fourier space
	 *@author David Woolford
	 *@date 2007/12/04
	 */
	class ConvolutionProcessor : public Processor 
	{
		public:
			ConvolutionProcessor() {}
		  
			string get_name() const
			{
				return "convolution";
			}
		
			void process_inplace(EMData *image);

			static Processor *NEW()
			{
				return new ConvolutionProcessor();
			}

			string get_desc() const
			{
				return "Performs Fourier space convolution. Maintains the space that the image is in - i.e. if image is real, the result is real and vice versa.";
			}

			TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("with", EMObject::EMDATA, "The image that will convolute the other image");
				return d;
			}
	};
	
	/** Determines the partial derivatives in the x direction
	 * Does this by constructing edge kernels in real space but convoluting in Fourier space
	 *
	 *@author David Woolford
	 *@date 2007/12/04
	 */
	class XGradientProcessor : public Processor 
	{
	 public:
		XGradientProcessor() {}
		  
		string get_name() const
		{
			return "edge.xgradient";
		}
		
		void process_inplace(EMData *image);

		static Processor *NEW()
		{
			return new XGradientProcessor();
		}

		string get_desc() const
		{
			return "Determines the image gradient in the x direction";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			return d;
		}
	};
	
	class YGradientProcessor : public Processor 
	{
		public:
			YGradientProcessor() {}
		  
			string get_name() const
			{
				return "edge.ygradient";
			}
		
			void process_inplace(EMData *image);

			static Processor *NEW()
			{
				return new YGradientProcessor();
			}

			string get_desc() const
			{
				return "Determines the image gradient in the y direction";
			}
	

			TypeDict get_param_types() const
			{
				TypeDict d;
				return d;
			}
	};
	
	class ZGradientProcessor : public Processor 
	{
		public:
			ZGradientProcessor() {}
		  
			string get_name() const
			{
				return "edge.zgradient";
			}
		
			void process_inplace(EMData *image);

			static Processor *NEW()
			{
				return new ZGradientProcessor();
			}

			string get_desc() const
			{
				return "Determines the image gradient in the z direction";
			}
		
			TypeDict get_param_types() const
			{
				TypeDict d;
				return d;
			}
	};
	
	/** Automatically determines the background for the image then uses this to perform
	 * Wiener filters on overlapping subregions of the image, which are then 
	 * combined using linear interpolation 
	 *@param size[in] size in pixels of the boxes to chop the image into during processing
	 *
	 *@author Steve Ludtke
	 *@date 2007/11/06
	 */
	class Wiener2DAutoAreaProcessor:public Processor
	{
	  public:
		string get_name() const
		{
			return "filter.wiener2dauto";
		}

		EMData* process(const EMData * const image);
		
		void process_inplace(EMData *image);

		void set_params(const Dict & new_params)
		{
			params = new_params;
			bgsize = params["size"];
//			printf("%s %f\n",params.keys()[0].c_str(),lowpass);
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("size", EMObject::INT, "Size in pixels of the boxes to chop the image into");
			return d;
		}

		static Processor *NEW()
		{
			return new Wiener2DAutoAreaProcessor();
		}

		string get_desc() const
		{
			return "Automatically detrmines the background for the image then uses this to perform Wiener filters on overlapping subregions of the image, which are then combined using linear interpolation";
		}

		protected:
		int bgsize;
	};


	/**Wiener filter based on externally provided background image
	 *@param background[in] Fourier intensity image containing the background to use in SNR calculation
	 *
	 *@author Steve Ludtke
	 *@date 2007/11/06
	 */
	class Wiener2DFourierProcessor:public Processor
	{
	  public:
		string get_name() const
		{
			return "filter.wiener2d";
		}

		EMData* process(const EMData * const image);
		
		void process_inplace(EMData *image);

		void set_params(const Dict & new_params)
		{
			params = new_params;
			background = params["background"];
//			printf("%s %f\n",params.keys()[0].c_str(),lowpass);
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("background", EMObject::EMDATA, "Fourier image containing averaged background from the frame");
			return d;
		}

		static Processor *NEW()
		{
			return new Wiener2DFourierProcessor();
		}

		string get_desc() const
		{
			return "Applies a 2-D Wiener filter to an image based on the provided background image, which must be the same size as the image being filtered.";
		}

		protected:
		EMData *background;
	};

	/**Low-pass processor attenuates amplitudes at high spatial frequencies. It has the result of blurring the image, and of eliminating sharp edges and noise. The base class for all low pass fourier processors.
	 *@param lowpass Processor radius in terms of Nyquist (0-.5)
	 */
	class LowpassFourierProcessor:public FourierProcessor
	{
	  public:
		LowpassFourierProcessor():lowpass(0)
		{
		}

		void set_params(const Dict & new_params)
		{
			params = new_params;
			if( params.has_key("cutoff_abs") ) {
				lowpass = params["cutoff_abs"];
			}
			else if( params.has_key("lowpass") ) {
				lowpass = params["lowpass"];
			}
//			printf("%s %f\n",params.keys()[0].c_str(),lowpass);
		}

		TypeDict get_param_types() const
		{
			TypeDict d = FourierProcessor::get_param_types();
			d.put("lowpass", EMObject::FLOAT, "Processor radius in terms of Nyquist (0-.5)");
			return d;
		}
		
		static string get_group_desc() 
		{
			return "Low-pass processor attenuates amplitudes at high spatial frequencies. It has the result of blurring the image, and of eliminating sharp edges and noise. The base class for all low pass fourier processors.";
		}
		
	  protected:
		float lowpass;
	};

	/**High-pass processor is rotationally symmetric 2D function. It attenuates amplitudes at low spatial frequencies, and increases amplitudes for high spatial frequencies. It has the result of enhancing the edges in the image while suppressing all slow-moving variations.	<br> HighpassFourierProcessor class is the base class for all high pass fourier processors.
	 *@param highpass Processor radius in terms of Nyquist (0-.5)
	 */
	class HighpassFourierProcessor:public FourierProcessor
	{
	  public:
		HighpassFourierProcessor():highpass(0)
		{
		}

		void set_params(const Dict & new_params)
		{
			params = new_params;
			if( params.has_key("cutoff_abs") ) {
				highpass = params["cutoff_abs"];
			}
			else if( params.has_key("highpass") ) {
				highpass = params["highpass"];
			}
		}

		TypeDict get_param_types() const
		{
			TypeDict d = FourierProcessor::get_param_types();
			d.put("highpass", EMObject::FLOAT, "Processor radius in terms of Nyquist (0-.5)");
			return d;
		}
		
		static string get_group_desc()
		{
			return "High-pass processor is rotationally symmetric 2D function. It attenuates amplitudes at low spatial frequencies, and increases amplitudes for high spatial frequencies. It has the result of enhancing the edges in the image while suppressing all slow-moving variations.	<br> HighpassFourierProcessor class is the base class for all high pass fourier processors.";
		}
		
	  protected:
		float highpass;
	};

	/**processor radial function: if x <= lowpass, f(x) = 1; else f(x) = 0;
	 */
	class LowpassSharpCutoffProcessor:public LowpassFourierProcessor
	{
	  public:
		string get_name() const
		{
			return "eman1.filter.lowpass.sharp";
		}

		static Processor *NEW()
		{
			return new LowpassSharpCutoffProcessor();
		}

		string get_desc() const
		{
			return "processor radial function: if x <= lowpass, f(x) = 1; else f(x) = 0;";
		}
		
	  protected:
		void create_radial_func(vector < float >&radial_mask) const;
	};

	/**processor radial function: if x >= highpass, f(x) = 1; else f(x) = 0;
	 */
	class HighpassSharpCutoffProcessor:public HighpassFourierProcessor
	{
	  public:
		string get_name() const
		{
			return "eman1.filter.highpass.sharp";
		}

		static Processor *NEW()
		{
			return new HighpassSharpCutoffProcessor();
		}

		string get_desc() const
		{
			return "processor radial function: if x >= highpass, f(x) = 1; else f(x) = 0;";
		}

		
	  protected:
		void create_radial_func(vector < float >&radial_mask) const;
	};

	/**processor radial function: if lowpass > 0, f(x) = exp(-x*x/(lowpass*lowpass)); else f(x) = exp(x*x/(lowpass*lowpass))
	 */
	class LowpassGaussProcessor:public LowpassFourierProcessor
	{
	  public:
		string get_name() const
		{
			return "eman1.filter.lowpass.gaussian";
		}

		static Processor *NEW()
		{
			return new LowpassGaussProcessor();
		}

		string get_desc() const
		{
			return "processor radial function: if lowpass > 0, f(x) = exp(-x*x/(lowpass*lowpass)); else f(x) = exp(x*x/(lowpass*lowpass));";
		}

	  protected:
		void create_radial_func(vector < float >&radial_mask) const;
	};

	/**processor radial function: f(x) = 1.0-exp(-x*x/(highpass*highpass))
	 */
	class HighpassGaussProcessor:public HighpassFourierProcessor
	{
	  public:
		string get_name() const
		{
			return "eman1.filter.highpass.gaussian";
		}
		static Processor *NEW()
		{
			return new HighpassGaussProcessor();
		}

		string get_desc() const
		{
			return "processor radial function: f(x) = 1.0-exp(-x*x/(highpass*highpass);";
		}

	  protected:
		void create_radial_func(vector < float >&radial_mask) const;
	};

	/**processor radial function: f(x)=tanh(lowpass-x)/2.0 + 0.5;
	 */
	class LowpassTanhProcessor:public LowpassFourierProcessor
	{
	  public:
		string get_name() const
		{
			return "eman1.filter.lowpass.tanh";
		}
		static Processor *NEW()
		{
			return new LowpassTanhProcessor();
		}

		string get_desc() const
		{
			return "processor radial function: f(x)=tanh(lowpass-x)/2.0 + 0.5;";
		}

	  protected:
		void create_radial_func(vector < float >&radial_mask) const;
	};


	/**processor radial function: f(x)=tanh(x-highpass)/2.0+0.5;
	 */
	class HighpassTanhProcessor:public HighpassFourierProcessor
	{
	  public:
		string get_name() const
		{
			return "eman1.filter.highpass.tanh";
		}
		static Processor *NEW()
		{
			return new HighpassTanhProcessor();
		}

		string get_desc() const
		{
			return "processor radial function: f(x)=tanh(x-highpass)/2.0+0.5;";
		}

	  protected:
		void create_radial_func(vector < float >&radial_mask) const;
	};

	/**processor radial function: f(x) = 1/(1+t*t)
	 */
	class HighpassButterworthProcessor:public HighpassFourierProcessor
	{
	  public:
		string get_name() const
		{
			return "eman1.filter.highpass.butterworth";
		}
		static Processor *NEW()
		{
			return new HighpassButterworthProcessor();
		}

		string get_desc() const
		{
			return "processor radial function: f(x) = 1/(1+t*t);";
		}

	  protected:
		void create_radial_func(vector < float >&radial_mask) const;
	};

	/**processor radial function: f(x) = slope * x + intercept
	 */
	class LinearRampProcessor:public FourierProcessor
	{
	  public:
		LinearRampProcessor():intercept(0), slope(0)
		{
		}

		string get_name() const
		{
			return "eman1.filter.ramp";
		}
		static Processor *NEW()
		{
			return new LinearRampProcessor();
		}

		string get_desc() const
		{
			return "processor radial function: f(x) = slope * x + intercept;";
		}

		void set_params(const Dict & new_params)
		{
			params = new_params;
			intercept = params["intercept"];
			slope = params["slope"];
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("intercept", EMObject::FLOAT);
			d.put("slope", EMObject::FLOAT);
			return d;
		}

	  protected:
		void create_radial_func(vector < float >&radial_mask) const;

	  private:
		float intercept;
		float slope;
	};

	/**The base class for real space processor working on individual pixels. The processor won't consider the pixel's coordinates and neighbors.
	 */
	class RealPixelProcessor:public Processor
	{
	  public:
		RealPixelProcessor():value(0), maxval(1), mean(0), sigma(0)
		{
		}
		void process_inplace(EMData * image);

		void set_params(const Dict & new_params)
		{
			params = new_params;
			if (params.size() == 1) {
				vector < EMObject > dict_values = params.values();
				value = dict_values[0];
			}
		}

		static string get_group_desc()
		{
			return "The base class for real space processor working on individual pixels. The processor won't consider the pixel's coordinates and neighbors.";
		}
		
	  protected:
		virtual void process_pixel(float *x) const = 0;
		virtual void calc_locals(EMData *)
		{
		}
		virtual void normalize(EMData *) const
		{
		}

		float value;
		float maxval;
		float mean;
		float sigma;
	};

	/**f(x) = |x|
	 */
	class AbsoluateValueProcessor:public RealPixelProcessor
	{
	  public:
		string get_name() const
		{
			return "math.absvalue";
		}
		static Processor *NEW()
		{
			return new AbsoluateValueProcessor();
		}
	  protected:
		void process_pixel(float *x) const
		{
			*x = fabs(*x);
		}

		string get_desc() const
		{
			return "f(x) = |x|";
		}
		
	};
	

	/**f(x) = 0 if x = 0; f(x) = 1 if x != 0
	 */
	class BooleanProcessor:public RealPixelProcessor
	{
	  public:
		string get_name() const
		{
			return "threshold.notzero";
		}
		static Processor *NEW()
		{
			return new BooleanProcessor();
		}

		string get_desc() const
		{
			return "f(x) = 0 if x = 0; f(x) = 1 if x != 0;";
		}
		
	  protected:
		void process_pixel(float *x) const
		{
			if (*x != 0)
			{
				*x = 1.0;
			}
		}
	};

	class ValuePowProcessor:public RealPixelProcessor
	{
	  public:
		string get_name() const
		{
			return "math.pow";
		}
		static Processor *NEW()
		{
			return new ValuePowProcessor();
		}

		void set_params(const Dict & new_params)
		{
			params = new_params;
			pwr = params["pow"];
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("pow", EMObject::FLOAT, "Each pixel is raised to this power");
			return d;
		}

		string get_desc() const
		{
			return "f(x) = x ^ pow;";
		}
		
	  protected:
		void process_pixel(float *x) const
		{
			if (*x<0 && pwr!=(int)pwr) *x=0;
			else (*x) = pow(*x,pwr);
		}
	  private:
		float pwr;
	};


	class ValueSquaredProcessor:public RealPixelProcessor
	{
	  public:
		string get_name() const
		{
			return "math.squared";
		}
		static Processor *NEW()
		{
			return new ValueSquaredProcessor();
		}


		string get_desc() const
		{
			return "f(x) = x * x;";
		}
		
	  protected:
		void process_pixel(float *x) const
		{
			(*x) *= (*x);
		}
	};

	/**f(x) = sqrt(x)
	 */
	class ValueSqrtProcessor:public RealPixelProcessor
	{
	  public:
		string get_name() const
		{
			return "math.sqrt";
		}
		static Processor *NEW()
		{
			return new ValueSqrtProcessor();
		}

		string get_desc() const
		{
			return "f(x) = sqrt(x)";
		}
		
	  protected:
		void process_pixel(float *x) const
		{
			*x = sqrt(*x);
		}
	};
	
	/**f(x) = x if x >= minval; f(x) = 0 if x < minval
	*@param minval
	 */
	class ToZeroProcessor:public RealPixelProcessor
	{
		public:
			string get_name() const
			{
				return "threshold.belowtozero";
			}
			static Processor *NEW()
			{
				return new ToZeroProcessor();
			}
			TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("minval", EMObject::FLOAT, "Everything below this value is set to zero");
				return d;
			}

			string get_desc() const
			{
				return "f(x) = x if x >= minval; f(x) = 0 if x < minval.";
			}
		
		protected:
			inline void process_pixel(float *x) const
			{
				if (*x < value) {
					*x = 0;
				}
			}
	};

	/**f(x) = maxval if f(x) > maxval;
	  * f(x) = minval if f(x) < minval
	  * @param minval the minimum value to clamp to
	  * @param maxval the maximum value to clamp to
	 */
	class ClampingProcessor :public Processor
	{
	  public:
		ClampingProcessor() : default_max(1.0), default_min(0.0) {}
		
		string get_name() const
		{
			return "threshold.clampminmax";
		}
		static Processor *NEW()
		{
			return new ClampingProcessor();
		}
		
		void process_inplace(EMData *image);
		
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("minval", EMObject::FLOAT, "The pixel values that bounds the smallest pixel value in the output image" );
			d.put("maxval", EMObject::FLOAT, "The pixel values that bounds the largest pixel value in the output image" );
			return d;
		}

		string get_desc() const
		{
			return "This function clamps the min and max vals in the image at minval and maxval, respectively. In a sense this a bi-truncation of the data.";
		}
		
	  protected:
		float default_max, default_min;
	};

	class NSigmaClampingProcessor : public ClampingProcessor
	{
		public:
			NSigmaClampingProcessor() : default_sigma(2.0) {}
			
			string get_name() const
			{
				return "threshold.clampminmax.nsigma";
			}
			
			static Processor *NEW()
			{
				return new NSigmaClampingProcessor();
			}
			
			TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("nsigma", EMObject::FLOAT, "The number (n) of sigmas to clamp min and max vals at, so that the clamped boundaries are mean-n*sigma and mean+n*sigma" );
				return d;
			}
		
			void process_inplace(EMData *image);
			
			string get_desc() const
			{
				return "This function clamps the min and max vals in the image at minval and maxval at mean-n*sigma and mean+n*sigma, respectively. The parameter specified by the user is n, the default value of n is 2.";
			}
			
		protected:
			float default_sigma;
	};

	/**f(x) = x if x >= minval; f(x) = minval if x < minval
	 *@param minval
	 */
	class ToMinvalProcessor:public RealPixelProcessor
	{
	  public:
		string get_name() const
		{
			return "threshold.belowtominval";
		}
		static Processor *NEW()
		{
			return new ToMinvalProcessor();
		}
		
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("minval", EMObject::FLOAT, "Everything below this value is set to this value");
			return d;
		}

		string get_desc() const
		{
			return "f(x) = x if x >= minval; f(x) = minval if x < minval.";
		}
		
	protected:
		inline void process_pixel(float *x) const
		{
			if (*x < value) {
				*x = value;
			}
		}
	};


	/**f(x) = x-minval if x >= minval; f(x) = 0 if x < minval
	 *@param minval
	 */
	class CutToZeroProcessor:public RealPixelProcessor
	{
	  public:
		string get_name() const
		{
			return "threshold.belowtozero_cut";
		}
		static Processor *NEW()
		{
			return new CutToZeroProcessor();
		}
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("minval", EMObject::FLOAT, "the value that will be set to zero - all values below will also be set to zero. Values above get minval subtracted from them" );
			return d;
		}

		string get_desc() const
		{
			return "f(x) = x-minval if x >= minval; f(x) = 0 if x < minval.";
		}
		
	  protected:
		void process_pixel(float *x) const
		{
		        *x = *x - value;
			if (*x < 0) {
				*x = 0;
			}
		}
	};
	
	/**f(x) = 0 if x < value; f(x) = 1 if x >= value.
	 *@param value
	 */
	class BinarizeProcessor:public RealPixelProcessor
	{
	  public:
		string get_name() const
		{
			return "threshold.binary";
		}
		static Processor *NEW()
		{
			return new BinarizeProcessor();
		}
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("value", EMObject::FLOAT, "The thresholding value. If a pixel value is equal to or above the threshold it is set to 1. If it is below it is set to 0" );
			return d;
		}

		string get_desc() const
		{
			return "f(x) = 0 if x < value; f(x) = 1 if x >= value.";
		}
		
	  protected:
		void process_pixel(float *x) const
		{
			if (*x < value)
			{
				*x = 0;
			}
			else
			{
				*x = 1;
			}
		}
	};

	/**f(x): if v-r<x<v+r -> v; if x>v+r -> x-r; if x<v-r -> x+r
	 *@param range
	 *@param value
	 */
	class CollapseProcessor:public RealPixelProcessor
	{
	  public:
		string get_name() const
		{
			return "threshold.compress";
		}
		static Processor *NEW()
		{
			return new CollapseProcessor();
		}

		void set_params(const Dict & new_params)
		{
			params = new_params;
			range = params["range"];
			value = params["value"];
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("range", EMObject::FLOAT, "The range about 'value' which will be collapsed to 'value'");
			d.put("value", EMObject::FLOAT, "The pixel value where the focus of the collapse operation is");
			return d;
		}

		string get_desc() const
		{
			return "f(x): if v-r<x<v+r -> v; if x>v+r -> x-r; if x<v-r -> x+r";
		}
		
	  protected:
		void process_pixel(float *x) const
		{
			if (*x>value+range) *x-=range;
			else if (*x<value-range) *x+=range;
			else *x=value;
		}
		float range;
	};

	/**linear transform processor: f(x) = x * scale + shift
	 *@param shift
	 *@param scale
	 */
	class LinearXformProcessor:public RealPixelProcessor
	{
	  public:
		LinearXformProcessor():shift(0), scale(0)
		{
		}

		string get_name() const
		{
			return "math.linear";
		}
		static Processor *NEW()
		{
			return new LinearXformProcessor();
		}

		void set_params(const Dict & new_params)
		{
			params = new_params;
			shift = params.get("shift");
			scale = params.get("scale");
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("shift", EMObject::FLOAT, "The amount to shift pixel values by before scaling");
			d.put("scale", EMObject::FLOAT, "The scaling factor to be applied to pixel values");
			return d;
		}

		string get_desc() const
		{
			return "linear transform processor: f(x) = x * scale + shift. This is equivalent to a regular contrast stretching operation";
		}
		
	  protected:
		void process_pixel(float *x) const
		{
			*x = (*x) * scale + shift;
		}

	  private:
		float shift;
		float scale;
	};

	/**f(x) = exp( x / low - high)
	 *@param low
	 *@param high
	 */
	class ExpProcessor:public RealPixelProcessor
	{
	  public:
		ExpProcessor():low(0), high(0)
		{
		}

		string get_name() const
		{
			return "math.exp";
		}

		static Processor *NEW()
		{
			return new ExpProcessor();
		}

		void set_params(const Dict & new_params)
		{
			params = new_params;
			low = params.get("low");
			high = params.get("high");
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("low", EMObject::FLOAT, "Pixels are divided by (low - high) prior to the exponential operation");
			d.put("high", EMObject::FLOAT, "Pixels are divided by (low - high) prior to the exponential operation");
			return d;
		}

		string get_desc() const
		{
			return "f(x) = exp( x / low - high)";
		}

	  protected:
	/**
	 * '40' is used to avoid floating number overflow.
	 */
		void process_pixel(float *x) const
		{
			float v = *x / low - high;
			if (v > 40) {
				v = 40;
			}
			*x = exp(v);
		}

	  private:
		float low;
		float high;
	};

	/**f(x) = 1 if (low <= x <= high); else f(x) = 0
	 *@param low
	 *@param high
	 */
	class RangeThresholdProcessor:public RealPixelProcessor
	{
	  public:
		RangeThresholdProcessor():low(0), high(0)
		{
		}

		string get_name() const
		{
			return "threshold.binaryrange";
		}
		static Processor *NEW()
		{
			return new RangeThresholdProcessor();
		}

		void set_params(const Dict & new_params)
		{
			params = new_params;
			low = params.get("low");
			high = params.get("high");
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("low", EMObject::FLOAT, "The lower limit of the range that will be set to 1");
			d.put("high", EMObject::FLOAT, "The upper limit of the range that will be set to 1");
			return d;
		}

		string get_desc() const
		{
			return "Range thresholding. A range of values is set to 1, all else is set to 0. f(x) = 1 if (low <= x <= high); else f(x) = 0";
		}
		
	  protected:
		void process_pixel(float *x) const
		{
			if (*x >= low && *x <= high) {
				*x = 1;
			}
			else {
				*x = 0;
			}
		}
	  private:
		float low;
		float high;

	};

	/**f(x) = mean if x<(mean-v2*sigma) or x>(mean+v1*sigma); else f(x) = x;
	 *@param value1
	 *@param value2
	 */
	class SigmaProcessor:public RealPixelProcessor
	{
	  public:
		string get_name() const
		{
			return "math.sigma";
		}
		static Processor *NEW()
		{
			return new SigmaProcessor();
		}

		void set_params(const Dict & new_params)
		{
			params = new_params;
			value1 = params.get("value1");
			value2 = params.get("value2");
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("value1", EMObject::FLOAT, "A number reflecting total standard deviations in the right direction");
			d.put("value2", EMObject::FLOAT, "A number reflecting total standard deviations in the left direction");
			return d;
		}

		string get_desc() const
		{
			return "f(x) = mean if x<(mean-v2*sigma) or x>(mean+v1*sigma); else f(x) = x;";
		}

	  protected:
		void process_pixel(float *x) const
		{
			if (*x < (mean - value2 * sigma) || *x > (mean + value1 * sigma))
			{
				*x = mean;
			}
		}

	  private:
		float value1;
		float value2;
	};

	/**f(x) = log10(x/max_pixel) if x > 0; else f(x) = 0
	 */
	class LogProcessor:public RealPixelProcessor
	{
	  public:
		string get_name() const
		{
			return "math.log";
		}
		static Processor *NEW()
		{
			return new LogProcessor();
		}

		string get_desc() const
		{
			return "f(x) = log10(x/max_pixel) if x > 0; else f(x) = 0;";
		}

	  protected:
		void process_pixel(float *x) const
		{
			if (*x > 0)
			{
				*x = log10(*x / maxval);
			}
			else
			{
				*x = 0;
			}
		}
	};

	/**CoordinateProcessor applies processing based on a pixel's value and it coordinates. This is the base class. Specific coordinate processor should implement process_pixel().
	 */
	class CoordinateProcessor:public Processor
	{
	  public:
		CoordinateProcessor():nx(0), ny(0), nz(0), mean(0), sigma(0), maxval(0), is_complex(false)
		{
		}
		void process_inplace(EMData * image);

		static string get_group_desc()
		{
			return "CoordinateProcessor applies processing based on a pixel's value and it coordinates. This is the base class. Specific coordinate processor should implement process_pixel().";
		}

	  protected:
		virtual void process_pixel(float *pixel, int xi, int yi, int zi) const = 0;
		virtual void calc_locals(EMData *)
		{
		}
		virtual bool is_valid() const
		{
			return true;
		}

		int nx;
		int ny;
		int nz;
		float mean;
		float sigma;
		float maxval;

		bool is_complex;
	};

	/**CircularMaskProcessor applies a circular mask to the data.This is the base class for specific circular mask processors.Its subclass must implement process_dist_pixel().
	 *@param inner_radius inner mask radius. optional, default=-1
	 *@param outer_radius outer mask radius
	 *@param dx Modify mask center by dx relative to the default center nx/2
	 *@param dy Modify mask center by dy relative to the default center ny/2
	 *@param dz Modify mask center by dz relative to the default center nz/2
	 */
	class CircularMaskProcessor:public CoordinateProcessor
	{
	  public:
		CircularMaskProcessor():inner_radius(0), outer_radius(0), inner_radius_square(0),
			outer_radius_square(0), dx(0), dy(0), dz(0), xc(0), yc(0), zc(0)
		{
		}

		void set_params(const Dict & new_params)
		{
			params = new_params;

			if (params.has_key("inner_radius")) {
				inner_radius = params["inner_radius"];
				inner_radius_square = inner_radius * inner_radius;
			}
			else {
				inner_radius = -1;
				inner_radius_square = -1;
			}

			if (params.has_key("outer_radius")) {
				outer_radius = params["outer_radius"];
				outer_radius_square = outer_radius * outer_radius;
			}
			else {
				outer_radius = INT_MAX;
				outer_radius_square = INT_MAX;
			}

			if (params.has_key("xc")) xc = params["xc"];
			if (params.has_key("yc")) yc = params["yc"];
			if (params.has_key("zc")) zc = params["zc"];
			if (params.has_key("dx")) dx = params["dx"];
			if (params.has_key("dy")) dy = params["dy"];
			if (params.has_key("dz")) dz = params["dz"];
		}
		
		string get_desc() const
		{
			return "CircularMaskProcessor applies a circular mask to the data.This is the base class for specific circular mask processors.Its subclass must implement process_dist_pixel().";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;

			d.put("inner_radius", EMObject::INT, "inner mask radius. optional, default=-1");
			d.put("outer_radius", EMObject::INT, "outer mask radius");

			d.put("dx", EMObject::FLOAT,
				  "Modify mask center by dx relative to the default center nx/2");
			d.put("dy", EMObject::FLOAT,
				  "Modify mask center by dy relative to the default center ny/2");
			d.put("dz", EMObject::FLOAT,
				  "Modify mask center by dz relative to the default center nz/2");

			return d;
		}
	  protected:
		void calc_locals(EMData * image);

		bool is_valid() const
		{
			return (!is_complex);
		}

		void process_pixel(float *pixel, int xi, int yi, int zi) const
		{
			float dist = (xi - xc) * (xi - xc) + (yi - yc) * (yi - yc) + (zi - zc) * (zi - zc);
			process_dist_pixel(pixel, dist);
		}

		virtual void process_dist_pixel(float *pixel, float dist) const = 0;

		int inner_radius;
		int outer_radius;
		int inner_radius_square;
		int outer_radius_square;
		float dx, dy, dz;
		float xc, yc, zc;
	};

	/**step cutoff to a user-given value in both inner and outer circles.
	 *@param value step cutoff to this value
	 */
	class MaskSharpProcessor:public CircularMaskProcessor
	{
	  public:
		MaskSharpProcessor():value(0)
		{
		}

		string get_name() const
		{
			return "mask.sharp";
		}
		static Processor *NEW()
		{
			return new MaskSharpProcessor();
		}

		void set_params(const Dict & new_params)
		{
			CircularMaskProcessor::set_params(new_params);
			value = params["value"];
		}

		TypeDict get_param_types() const
		{
			TypeDict d = CircularMaskProcessor::get_param_types();
			d.put("value", EMObject::FLOAT, "step cutoff to this value.");
			return d;
		}
		
		string get_desc() const
		{
			return "step cutoff to a user-given value in both inner and outer circles.";
		}

	  protected:
		void process_dist_pixel(float *pixel, float dist) const
		{
			if (dist >= outer_radius_square || dist < inner_radius_square)
			{
				*pixel = value;
			}
		}

		float value;
	};

	
	/**A step cutoff to the the mean value in a ring centered on the outer radius
	 *@param ring_width The width of the mask ring.
	 */
	class MaskEdgeMeanProcessor:public CircularMaskProcessor
	{							// 6
	  public:
		string get_name() const
		{
			return "mask.ringmean";
		}
		static Processor *NEW()
		{
			return new MaskEdgeMeanProcessor();
		}

		void set_params(const Dict & new_params)
		{
			CircularMaskProcessor::set_params(new_params);
			ring_width = params["ring_width"];
			if (ring_width == 0) {
				ring_width = 1;
			}
		}

		TypeDict get_param_types() const
		{
			TypeDict d = CircularMaskProcessor::get_param_types();
			d.put("ring_width", EMObject::INT, "The width of the mask ring.");
			return d;
		}

		string get_desc() const
		{
			return "A step cutoff to the the mean value in a ring centered on the outer radius";
		}

	  protected:
		void calc_locals(EMData * image);


		void process_dist_pixel(float *pixel, float dist) const
		{
			if (dist >= outer_radius_square){
				*pixel = ring_avg;
			}
		}

	  private:
		int ring_width;
		float ring_avg;
	};

	/**fills masked region
	 */
	class MaskNoiseProcessor:public CircularMaskProcessor
	{
	  public:
		string get_name() const
		{
			return "mask.noise";
		}
		static Processor *NEW()
		{
			return new MaskNoiseProcessor();
		}

		string get_desc() const
		{
			return "fills masked region";
		}

	  protected:
		void process_dist_pixel(float *pixel, float dist) const
		{
			if (dist >= outer_radius_square || dist < inner_radius_square)
			{
				*pixel = Util::get_gauss_rand(mean, sigma);
			}
		}
	};

	/**a gaussian falloff to zero, radius is the 1/e of the width.
	 */
	class MaskGaussProcessor:public CircularMaskProcessor
	{
	  public:
		string get_name() const
		{
			return "mask.gaussian";
		}
		static Processor *NEW()
		{
			return new MaskGaussProcessor();
		}

		string get_desc() const
		{
			return "a gaussian falloff to zero, radius is the 1/e of the width.";
		}

	  protected:
		void process_dist_pixel(float *pixel, float dist) const
		{
			(*pixel) *= exp(-dist / outer_radius_square);
		}
	};

	/**f(x) = f(x) / exp(-radius*radius * gauss_width / (ny*ny))
	 *@param gauss_width
	 */
	class MaskGaussInvProcessor:public CircularMaskProcessor
	{
	  public:
		TypeDict get_param_types() const
		{
			TypeDict d = CircularMaskProcessor::get_param_types();
			d.put("gauss_width", EMObject::FLOAT, "Used to calculate the constant factor - gauss_width / (ny*ny)" );
			d.put("ring_width", EMObject::INT, "The width of the mask ring.");
			return d;
		}

		string get_name() const
		{
			return "math.gausskernelfix";
		}

		static Processor *NEW()
		{
			return new MaskGaussInvProcessor();
		}

		string get_desc() const
		{
			return "f(x) = f(x) / exp(-radius*radius * gauss_width / (ny*ny))";
		}
		
	  protected:
		void calc_locals(EMData *)
		{
			float gauss_width = params["gauss_width"];
			slice_value = gauss_width / (ny * ny);
		}

		void process_dist_pixel(float *pixel, float dist) const
		{
			(*pixel) /= exp(-dist * slice_value);
		}
	  private:
		float slice_value;
	};


	/**Multiplies image by a 'linear pyramid'
       1-(|x-xsize/2|*|y-ysize/2|*4/(xsize*ysize))
       This is useful in averaging together boxed out regions with 50% overlap.
	 */
	class LinearPyramidProcessor:public Processor
	{
	  public:
		string get_name() const
		{
			return "math.linearpyramid";
		}

		void process_inplace(EMData *image);

		static Processor *NEW()
		{
			return new LinearPyramidProcessor();
		}

		string get_desc() const
		{
			return "Multiplies image by a 'linear pyramid', 1-(|x-xsize/2|*|y-ysize/2|*4/(xsize*ysize))";
		}

	};


	/**overwrites input, f(x) = radius * radius
	 */
	class MakeRadiusSquaredProcessor:public CircularMaskProcessor
	{
	  public:
		string get_name() const
		{
			return "math.toradiussqr";
		}
		static Processor *NEW()
		{
			return new MakeRadiusSquaredProcessor();
		}

		string get_desc() const
		{
			return "overwrites input, f(x) = radius * radius";
		}

	  protected:
		void process_dist_pixel(float *pixel, float dist) const
		{
			*pixel = dist;
		}
	};

	/**overwrites input, f(x) = radius
	 */
	class MakeRadiusProcessor:public CircularMaskProcessor
	{
	  public:
		string get_name() const
		{
			return "math.toradius";
		}
		static Processor *NEW()
		{
			return new MakeRadiusProcessor();
		}

		string get_desc() const
		{
			return "overwrites input, f(x) = radius;";
		}

	  protected:
		void process_dist_pixel(float *pixel, float dist) const
		{
			*pixel = sqrt(dist);
		}
	};

	/**The base class for fourier space processor working on individual pixels. ri2ap() is called before processing, so individual pixels will be A/P rather than R/I. The processor won't consider the pixel's coordinates and neighbors.
	 */
	class ComplexPixelProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);
		
		static string get_group_desc()
		{
			return "The base class for fourier space processor working on individual pixels. ri2ap() is called before processing, so individual pixels will be A/P rather than R/I. The processor won't consider the pixel's coordinates and neighbors.";
		}
		
	  protected:
		virtual void process_pixel(float *x) const = 0;
	};

	/**Each Fourier pixel will be normalized. ie - amp=1, phase=unmodified. Useful for performing phase-residual-like computations with dot products.
	 */
	class ComplexNormPixel:public ComplexPixelProcessor
	{
	  public:
		string get_name() const
		{
			return "complex.normpixels";
		}
		static Processor *NEW()
		{
			return new ComplexNormPixel();
		}

		string get_desc() const
		{
			return "Each Fourier pixel will be normalized. ie - amp=1, phase=unmodified. Useful for performing phase-residual-like computations with dot products.";
		}

	  protected:
		void process_pixel(float *x) const
		{
			*x=1.0;
		}
	};

	/**AreaProcessor use pixel values and coordinates of a real-space square area. This is the base class. Specific AreaProcessor needs to implement function create_kernel().
	 *@param areasize
	 */
	class AreaProcessor:public Processor
	{
	  public:
		AreaProcessor():areasize(0), kernel(0), nx(0), ny(0), nz(0)
		{
		}

		void process_inplace(EMData * image);

		void set_params(const Dict & new_params)
		{
			params = new_params;
			areasize = params["areasize"];
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("areasize", EMObject::INT, "The width of the area to process (not radius)");
			return d;
		}

		string get_desc() const
		{
			return "AreaProcessor use pixel values and coordinates of a real-space square area. This is the base class. Specific AreaProcessor needs to implement function create_kernel().";
		}
		
	  protected:
		virtual void process_pixel(float *pixel, float, float, float, float *area_matrix) const
		{
			for (int i = 0; i < matrix_size; i++)
			{
				*pixel += area_matrix[i] * kernel[i];
			}
		}

		virtual void create_kernel() const = 0;

		int areasize;
		int matrix_size;
		float *kernel;
		int nx;
		int ny;
		int nz;
	};

	/**Discrete approximation to Laplacian. Edge enchancement, but works poorly in the presence of noise. Laplacian processor (x -> d^2/dx^2 + d^2/dy^2 + d^2/dz^2).
	 */
	class LaplacianProcessor:public AreaProcessor
	{
	  public:
		string get_name() const
		{
			return "math.laplacian";
		}
		static Processor *NEW()
		{
			return new LaplacianProcessor();
		}

		string get_desc() const
		{
			return "Discrete approximation to Laplacian. Edge enchancement, but works poorly in the presence of noise. Laplacian processor (x -> d^2/dx^2 + d^2/dy^2 + d^2/dz^2).";
		}

	  protected:
		void create_kernel() const;

	};

	/**Contraction of data, if any nearest neighbor is 0, value -> 0, generally used iteratively
	 */
	class ZeroConstantProcessor:public AreaProcessor
	{
	  public:
		string get_name() const
		{
			return "mask.contract";
		}
		static Processor *NEW()
		{
			return new ZeroConstantProcessor();
		}

		string get_desc() const
		{
			return "Contraction of data, if any nearest neighbor is 0, value -> 0, generally used iteratively";
		}
		
	  protected:
		void process_pixel(float *pixel, float, float, float, float *matrix) const
		{
			if (*pixel != 0)
			{
				if (*pixel == matrix[1] || *pixel == matrix[3] || *pixel == matrix[5] ||
					*pixel == matrix[7] || matrix[1] == 0 || matrix[3] == 0 ||
					matrix[5] == 0 || matrix[7] == 0) {
					*pixel = 0;
				}
			}
		}

		void create_kernel() const
		{
		}
	};

	/**BoxStatProcessor files are a kind of neighborhood processors. These processors compute every output pixel using information from a reduced region on the neighborhood of the input pixel. The classical form are the 3x3 processors. BoxStatProcessors could perform diverse tasks ranging from noise reduction, to differential , to mathematical morphology. BoxStatProcessor class is the base class. Specific BoxStatProcessor needs to define process_pixel(float *pixel, const float *array, int n).
	 */
	class BoxStatProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		static string get_group_desc()
		{
			return "BoxStatProcessor files are a kind of neighborhood processors. These processors compute every output pixel using information from a reduced region on the neighborhood of the input pixel. The classical form are the 3x3 processors. BoxStatProcessors could perform diverse tasks ranging from noise reduction, to differential , to mathematical morphology. BoxStatProcessor class is the base class. Specific BoxStatProcessor needs to define process_pixel(float *pixel, const float *array, int n).";
		}
		
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("radius", EMObject::INT, "The radius of the search box, default is 1 which results in a 3x3 box (3 = 2xradius + 1)");
			return d;
		}
		
	  protected:
		virtual void process_pixel(float *pixel, const float *array, int n) const = 0;
	};


	/**A processor for noise reduction. pixel = median of values surrounding pixel.
	 */
	class BoxMedianProcessor:public BoxStatProcessor
	{
	  public:
		string get_name() const
		{
			return "eman1.filter.median";
		}
		static Processor *NEW()
		{
			return new BoxMedianProcessor();
		}

		string get_desc() const
		{
			return "A processor for noise reduction. pixel = median of values surrounding pixel.";
		}


	  protected:
		void process_pixel(float *pixel, const float *array, int n) const
		{
			float *data = new float[n];
			memcpy(data, array, sizeof(float) * n);

			for (int i = 0; i <= n / 2; i++)
			{
				for (int j = i + 1; j < n; j++)
				{
					if (data[i] < data[j]) {
						float t = data[i];
						data[i] = data[j];
						data[j] = t;
					}
				}
			}

			if (n % 2 != 0)
			{
				*pixel = data[n / 2];
			}
			else {
				*pixel = (data[n / 2] + data[n / 2 - 1]) / 2;
			}
			if( data )
			{
				delete[]data;
				data = 0;
			}
		}
	};

	/**pixel = standard deviation of values surrounding pixel.
	 */
	class BoxSigmaProcessor:public BoxStatProcessor
	{
	  public:
		string get_name() const
		{
			return "math.localsigma";
		}
		static Processor *NEW()
		{
			return new BoxSigmaProcessor();
		}

		string get_desc() const
		{
			return "pixel = standard deviation of values surrounding pixel.";
		}

	  protected:
		void process_pixel(float *pixel, const float *data, int n) const
		{
			float sum = 0;
			float square_sum = 0;
			for (int i = 0; i < n; i++)
			{
				sum += data[i];
				square_sum += data[i] * data[i];
			}

			float mean = sum / n;
			*pixel = sqrt(square_sum / n - mean * mean);
		}
	};

	/**peak processor: pixel = max of values surrounding pixel.
	 */
	class BoxMaxProcessor:public BoxStatProcessor
	{
	  public:
		string get_name() const
		{
			return "math.localmax";
		}
		static Processor *NEW()
		{
			return new BoxMaxProcessor();
		}

		string get_desc() const
		{
			return "peak processor: pixel = max of values surrounding pixel.";
		}

	  protected:
		void process_pixel(float *pixel, const float *data, int n) const
		{
			float maxval = -FLT_MAX;
			for (int i = 0; i < n; i++)
			{
				if (data[i] > maxval) {
					maxval = data[i];
				}
			}
			 *pixel = maxval;
		}
	};

	/**peak processor: pixel = pixel - max of values surrounding pixel. This is a sort of positive peak-finding algorithm.
	 */
	class MinusPeakProcessor:public BoxStatProcessor
	{
	  public:
		string get_name() const
		{
			return "math.submax";
		}
		static Processor *NEW()
		{
			return new MinusPeakProcessor();
		}

		string get_desc() const
		{
			return "peak processor: pixel = pixel - max of values surrounding pixel. This is a sort of positive peak-finding algorithm.";
		}

	  protected:
		void process_pixel(float *pixel, const float *data, int n) const
		{
			float maxval = -FLT_MAX;
			for (int i = 0; i < n; i++)
			{
				if (data[i] > maxval) {
					maxval = data[i];
				}
			}
			 *pixel -= maxval;
		}
	};

	/**peak processor -> if npeaks or more surrounding values >= value, value->0
	 *@param npeaks
	 */
	class PeakOnlyProcessor:public BoxStatProcessor
	{
	  public:
		string get_name() const
		{
			return "mask.onlypeaks";
		}
		static Processor *NEW()
		{
			return new PeakOnlyProcessor();
		}
		void set_params(const Dict & new_params)
		{
			params = new_params;
			npeaks = params["npeaks"];
			if (npeaks == 0) {
				npeaks = 1;
			}
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("npeaks", EMObject::INT);
			return d;
		}

		string get_desc() const
		{
			return "peak processor -> if npeaks or more surrounding values >= value, value->0";
		}

	  protected:
		void process_pixel(float *pixel, const float *data, int n) const
		{
			int r = 0;

			for (int i = 0; i < n; i++)
			{
				if (data[i] >= *pixel) {
					r++;
				}
			}

			if (r > npeaks)
			{
				*pixel = 0;
			}
		}
	  private:
		int npeaks;
	};

	/**peak processor -> if npeaks or more surrounding values >= value, value->0
	 *@param npeaks
	 */
	class XrayPixelProcessor:public BoxStatProcessor
	{
	  public:
		string get_name() const
		{
			return "mask.xraypixel";
		}
		static Processor *NEW()
		{
			return new XrayPixelProcessor();
		}
		void set_params(const Dict & new_params)
		{
			params = new_params;
			max_peak_ratio = params["max_peak_ratio"];
			if (max_peak_ratio == 0) {
				max_peak_ratio=0.1f;
			}
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("max_peak_ratio", EMObject::FLOAT,"Max ratio of neighbor amp/peak amp, default 0.1");
			return d;
		}

		string get_desc() const
		{
			return "Removes isolated pixels where all neighbors are less than a limiting fraction of the pixel amplitude. Will not work after standard normalization.";
		}

	  protected:
		void process_pixel(float *pixel, const float *data, int n) const
		{
			
			for (int i = 0; i < n; i++)
			{
				if (data[i] >= *pixel*max_peak_ratio) return;
			}

		}
	  private:
		float max_peak_ratio;
	};

	/**
	 */
	class DiffBlockProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);
		
		string get_name() const
		{
			return "eman1.filter.blockrange";
		}
		static Processor *NEW()
		{
			return new DiffBlockProcessor();
		}

		string get_desc() const
		{
			return "averages over cal_half_width, then sets the value in a local block";
		}
		
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("cal_half_width", EMObject::FLOAT);
			d.put("fill_half_width", EMObject::FLOAT);
			return d;
		}
	};

	/**Block processor, val1 is dx/dy, val2 is lp freq cutoff in pixels. Mystery processor.
	 *@param value1
	 *@param value2
	 */
	class CutoffBlockProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);
		
		string get_name() const
		{
			return "eman1.filter.blockcutoff";
		}
		static Processor *NEW()
		{
			return new CutoffBlockProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("value1", EMObject::FLOAT);
			d.put("value2", EMObject::FLOAT);
			return d;
		}

		string get_desc() const
		{
			return "Block processor, val1 is dx/dy, val2 is lp freq cutoff in pixels. Mystery processor.";
		}
		
	};

	/** BooleanShrinkProcessor encapsulates code common to
	* MaxShrinkProcessor and MinShrinkProcessor - the processors use more or less
	* identical code, the main difference being the logical operator.
	* Both of these instances are written at compile time using templates.
	*/
	class BooleanShrinkProcessor: public Processor
	{
		protected:
			/** Boolean shrink an image, returning the processed image
			* @param image the image to operate on
			* @param LogicOp a logical operator class
			* @exception ImageFormatException if the image is complex
			* @exception NullPointerException if the image pointer is null
			* @return the image that results from the operation
			 */
			template<class LogicOp>
			EMData* process(const EMData *const image);
			
			/** Boolean shrink an image inplace
			 * @param image the image to operate on
			 * @param LogicOp a logical operator class
			 * @exception ImageFormatException if the image is complex
			 * @exception NullPointerException if the image pointer is null
			 */
			template<class LogicOp>
			void process_inplace(EMData * image);
			
	};
	
	/** MaxShrinkProcessors shrinks an image by in an integer amount, 
	 * keeping the maximum pixel value - useful when constructing binary search
	 * trees in the marching cubes algorithm
	 * @author David Woolford
	 * @date September 2007
	 */
	class MaxShrinkProcessor:public BooleanShrinkProcessor
	{
		public:
			/** The max shrink processor has its own process function
			* to minise memory usage - if this function was not over written
			* the base Processor class would create copy of the input image
			* and hand it to the process_inplace function. This latter approach
			* mallocs more memory than necessary
			*/
			virtual EMData* process(const EMData *const image)
			{
				return BooleanShrinkProcessor::process<GreaterThan>(image);
			}
			
			// resizes the image
			virtual void process_inplace(EMData * image)
			{
				BooleanShrinkProcessor::process_inplace<GreaterThan>(image);
			}
			
			string get_desc() const
			{
				return "Shrink an image by a given amount (default 2), using the maximum value found in the pixel neighborhood.";
			}
			
			string get_name() const
			{
				return "math.maxshrink";
			}
			static Processor *NEW()
			{
				return new MaxShrinkProcessor();
			}
			
			TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("n", EMObject::INT, "The shrink factor");
				d.put("search", EMObject::INT, "The search area (cubic volume width, usually the same as shrink)");
				return d;
			}
		
		private:
			struct GreaterThan
			{
				inline bool operator()(float left,float right) const { return left > right; }
				inline float get_start_val() { return -10000000; }
			};
	};
	
	/** MinShrinkProcessor shrinks an image by in an integer amount, 
	 * keeping the minimum pixel value - useful when constructing binary search
	 * trees in the marching cubes algorithm
	 * @author David Woolford
	 * @date September 2007
	 */
	class MinShrinkProcessor:public BooleanShrinkProcessor
	{
		public:
			/** The min shrink processor has its own process function
			* to minise memory usage - if this function was not over written
			* the base Processor class would create copy of the input image
			* and hand it to the process_inplace function. This latter approach
			* mallocs and copies more memory than necessary
			*/
			virtual EMData* process(const EMData *const image)
			{
				return BooleanShrinkProcessor::process<LessThan>(image);
			}
			
			// resizes the image
			virtual void process_inplace(EMData * image)
			{
				BooleanShrinkProcessor::process_inplace<LessThan>(image);
			}
			string get_desc() const
			{
				return "Shrink an image by a given amount (default 2), using the minimum value found in the pixel neighborhood.";
			}
			
			string get_name() const
			{
				return "math.minshrink";
			}
			static Processor *NEW()
			{
				return new MinShrinkProcessor();
			}
			
			TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("n", EMObject::INT, "The shrink factor");
				d.put("search", EMObject::INT, "The search area (cubic volume width, usually the same as shrink)");
				return d;
			}
		
		private:
		struct LessThan
		{
			inline bool operator()(float left,float right) const { return left < right; }
			inline float get_start_val() { return 9999999999.0f; }
		};
	};
	/** MeanShrinkProcessor shrinks an image by in an integer amount (and optionally by 1.5) 
	 * taking the mean of the pixel neighbourhood
	 * @author David Woolford (But is basically a copy of the old EMData::mean_shrink, probably written by Steven Ludtke )
	 * @date May 2008
	 */
	class MeanShrinkProcessor : public Processor
	{
		public:
			/** The meanshrink processor has its own process function
			* to minise memory usage - if this function was not over written
			* the base Processor class would create copy of the input image
			* and hand it to the process_inplace function. This latter approach
	 		* mallocs and copies more memory than necessary
			* @param image the image that will be used to generate a 'mean shrunken' image
			* @exception ImageFormatException if the image is complex
			* @exception ImageDimensionException if the image is 1D
			* @exception InvalidValueException if the shrink amount is a nonzero integer, unless it is 1.5, which is an exceptional circumstance
			 */
			virtual EMData* process(const EMData *const image);
			
			/** Mean shrink inplace
			* @param image the image that will be 'mean shrunken' inplace
			* @exception ImageFormatException if the image is complex
			* @exception ImageDimensionException if the image is 1D
			* @exception InvalidValueException if the shrink amount is a nonzero integer, unless it is 1.5, which is an exceptional circumstance
			*/
			virtual void process_inplace(EMData * image);
			
			string get_desc() const
			{
				return "Shrink an image by a given amount (default 2), using the mean value found in the pixel neighborhood.";
			}
			
			string get_name() const
			{
				return "math.meanshrink";
			}
			static Processor *NEW()
			{
				return new MeanShrinkProcessor();
			}
			
			TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("n", EMObject::FLOAT, "The shrink factor");
				return d;
			}
		private:
			/** Accrue the local mean in the image 'from' to the image 'to' using the given shrinkfactor
			 * An internal function that encapsulates a routine common to both process and 
			 * process inplace
			 * @param to the smaller image that will store the mean values
			 * @param from the larger image that will be used to calculate the mean values
			 * @param shrinkfactor the shrink amount
			*/
			void accrue_mean(EMData* to, const EMData *const from, const int shrinkfactor);
			
			/** Accrue the local mean in the image 'from' to the image 'to' using the the special case shrink factor of 1.5
			 * This is an internal function that encapsulates a routine common to both process and 
			 * process inplace
			 * @param to the smaller image that will store the mean values
			 * @param from the larger image that will be used to calculate the mean values
			 */
			void accrue_mean_one_p_five(EMData* to, const EMData * const from);
	};
	

		/** MeanShrinkProcessor shrinks an image by in an integer amount
		* taking the median of the pixel neighbourhood
		* @author David Woolford (But is basically a copy of the old EMData::median_shrink, probably written by Steven Ludtke )
		* @date May 2008
		*/
		class MedianShrinkProcessor : public Processor
		{
		public:
			/** The medianshrink processor has its own process function
			* to minise memory usage - if this function was not over written
			* the base Processor class would create copy of the input image
			* and hand it to the process_inplace function. This latter approach
			* mallocs and copies more memory than necessary
			* @param image the image that will be used to generate a 'median shrunken' image
		 	* @exception ImageFormatException  if the image is complex
			* @exception InvalidValueException if the shrink amount is not a non zero, positive integer
			* @exception InvalidValueException if any of the image dimensions are not divisible by the the shrink amount 
			*/
			virtual EMData* process(const EMData *const image);
				
			/** Median shrink the image
			* @param image the image the image that will be 'median shrunken' inplace
			* @exception ImageFormatException  if the image is complex
			* @exception InvalidValueException if the shrink amount is not a non zero, positive integer
			* @exception InvalidValueException if any of the image dimensions are not divisible by the the shrink amount 
			*/
			virtual void process_inplace(EMData * image);
				
			string get_desc() const
			{
				return "Shrink an image by a given amount (default 2), using the median value found in the pixel neighborhood.";
			}
				
			string get_name() const
			{
				return "math.medianshrink";
			}
			static Processor *NEW()
			{
				return new MedianShrinkProcessor();
			}
				
			TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("n", EMObject::INT, "The shrink factor");
				return d;
			}
		private:
			/** Accrue the local median in the image 'from' to the image 'to' using the given shrinkfactor
			* An internal function that encapsulates a routine common to both process and 
			* process inplace
			* @param to the smaller image that will store the calculated median values
			* @param from the larger image that will be used to calculate the median values
			* @param shrinkfactor the shrink amount
			 */
			void accrue_median(EMData* to, const EMData* const from,const int shrink_factor);
		};
	
	
	/**Gradient remover, does a rough plane fit to find linear gradients.
	 */
	class GradientRemoverProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);
		
		string get_name() const
		{
			return "math.lineargradientfix";
		}
		static Processor *NEW()
		{
			return new GradientRemoverProcessor();
		}

		string get_desc() const
		{
			return "Gradient remover, does a rough plane fit to find linear gradients.";
		}
		
	};

	/**Gradient removed by least square plane fit
	 *@param mask[in] optional EMData object to mask the pixels used to fit the plane
	 *@param changeZero[in] optional bool to specify if the zero value pixels are modified
	 *@param planeParam[out] optional return parameters [nx, ny, nz, cx, cy, cz] for the fitted plane defined as (x-cx)*nx+(y-cy)*ny+(z-cz)*nz=0
	 *
	 *@author Wen Jiang
	 *@date 2006/7/18
	 */
	class GradientPlaneRemoverProcessor:public Processor
    {
	  public:
		void process_inplace(EMData * image);
		
		string get_name() const
		{
			return "filter.gradientPlaneRemover";
		}
		static Processor *NEW()
		{
			return new GradientPlaneRemoverProcessor();
		}

		string get_desc() const
		{
			return "Remove gradient by least square plane fit";
		}
		
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("mask", EMObject::EMDATA, "mask object: nonzero pixel positions will be used to fit plane. default = 0");
			d.put("changeZero", EMObject::INT, "if zero pixels are modified when removing gradient. default = 0");
			d.put("planeParam", EMObject::FLOATARRAY, "fitted plane parameters output");
			return d;
		}	
	};
	
	
	/** Flattens the background by subtracting the local mean
	 *@param map an EMData object that defines the extent of the local neighbourhood - will be used for convolution
	 *@param radius exclusive of the mask parameter, this defines the radius of a circle/sphere that will be used for local mean subtraction
	 *@author David Woolford 
	 *@date April 2008
	 */
	class FlattenBackgroundProcessor:public Processor
	{
		public:
			void process_inplace(EMData * image);
		
			string get_name() const
			{
				return "filter.flattenbackground";
			}
			
			static Processor *NEW()
			{
				return new FlattenBackgroundProcessor();
			}

			string get_desc() const
			{
				return "Flattens the background by subtracting the local mean";
			}
		
			TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("mask", EMObject::EMDATA, "A mask the defines the local neighborhood that will be used to find the local mean. Exclusive of the radius argument");
				d.put("radius", EMObject::INT, "The radius of circle/sphere that defines the local neighborhood. Exclusive of the mask argument");
				return d;
			}	
	};
	
	
	/**Ramp processor -- Fits a least-squares plane to the picture, and subtracts the plane from the picture.  A wedge-shaped overall density profile can thus be removed from the picture.
	 */
	class RampProcessor:public Processor
    {
	  public:
		void process_inplace(EMData * image);
		
		string get_name() const
		{
			return "filter.ramp";
		}
		static Processor *NEW()
		{
			return new RampProcessor();
		}

		string get_desc() const
		{
			return "Ramp processor -- Fits a least-squares plane "
				   "to the picture, and subtracts the plane from "
				   "the picture.  A wedge-shaped overall density "
				   "profile can thus be removed from the picture.";
		}
		
	};

	/**Tries to fix images scanned on the zeiss for poor ccd normalization.
	 */
	class VerticalStripeProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "math.verticalstripefix";
		}

		static Processor *NEW()
		{
			return new VerticalStripeProcessor();
		}

		string get_desc() const
		{
			return "Tries to fix images scanned on the zeiss for poor ccd normalization.";
		}
		
	};
	
	/**This will replace the image with a full-circle 2D fft amplitude rendering.
	 */
	class RealToFFTProcessor:public Processor
	{
		public:
		void process_inplace(EMData *image);
	
		string get_name() const
		{
			return "math.realtofft";
		}

		static Processor *NEW()
		{
			return new RealToFFTProcessor();
		}

		string get_desc() const
		{
			return "This will replace the image with a full-circle 2D fft amplitude rendering. Note that this renders amplitude, when intensity is more common.";
		}

	};
			

	/**Fill zeroes at edges with nearest horizontal/vertical value.
	 */
	class SigmaZeroEdgeProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "mask.zeroedgefill";
		}
		static Processor *NEW()
		{
			return new SigmaZeroEdgeProcessor();
		}

		string get_desc() const
		{
			return "Fill zeroes at edges with nearest horizontal/vertical value.";
		}
		
	};

	/**Try to eliminate beamstop in electron diffraction patterns. If value1<0 also does radial subtract.
	 *@param value1 sig multiplier
	 *@param value2 x of center
	 *@param value3 y of center
	 */
	class BeamstopProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "mask.beamstop";
		}

		static Processor *NEW()
		{
			return new BeamstopProcessor();
		}

		string get_desc() const
		{
			return "Try to eliminate beamstop in electron diffraction patterns. value1=sig multiplier; value2,value3 are x,y of center, if value1<0 also does radial subtract.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("value1", EMObject::FLOAT);
			d.put("value2", EMObject::FLOAT);
			d.put("value3", EMObject::FLOAT);
			return d;
		}
	};

	/**Fill zeroes at edges with nearest horizontal/vertical value damped towards Mean2.
	 */
	class MeanZeroEdgeProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "mask.dampedzeroedgefill";
		}

		static Processor *NEW()
		{
			return new MeanZeroEdgeProcessor();
		}

		string get_desc() const
		{
			return "Fill zeroes at edges with nearest horizontal/vertical value damped towards Mean2.";
		}
		
	};


	/**Average along Y and replace with average
	 */
	class AverageXProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "math.averageovery";
		}

		static Processor *NEW()
		{
			return new AverageXProcessor();
		}

		string get_desc() const
		{
			return "Average along Y and replace with average";
		}
		
	};

	/**zero edges of image on top and bottom, and on left and right.
	 *@param x0
	 *@param x1
	 *@param y0
	 *@param y1
	 */
	class ZeroEdgeRowProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);
		string get_name() const
		{
			return "mask.zeroedge2d";
		}

		static Processor *NEW()
		{
			return new ZeroEdgeRowProcessor();
		}

		string get_desc() const
		{
			return "zero edges of image on top and bottom, and on left and right.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("x0", EMObject::INT, "The number of columns to zero from left");
			d.put("x1", EMObject::INT, "The number of columns to zero from right");
			d.put("y0", EMObject::INT, "The number of rows to zero from the bottom");
			d.put("y1", EMObject::INT, "The number of rows to zero from the top");
			return d;
		}
	};

	/**zero edges of volume on all sides
	 *@param x0
	 *@param x1
	 *@param y0
	 *@param y1
	 *@param z0
	 *@param z1
	 */
	class ZeroEdgePlaneProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);
		string get_name() const
		{
			return "mask.zeroedge3d";
		}

		static Processor *NEW()
		{
			return new ZeroEdgePlaneProcessor();
		}

		string get_desc() const
		{
			return "zero edges of volume on all sides";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("x0", EMObject::INT, "The number of columns to zero from left");
			d.put("x1", EMObject::INT, "The number of columns to zero from right");
			d.put("y0", EMObject::INT, "The number of rows to zero from the bottom");
			d.put("y1", EMObject::INT, "The number of rows to zero from the top");
			d.put("z0", EMObject::INT, "The number of slices to zero from the bottom");
			d.put("z1", EMObject::INT, "The number of slices to zero from the top");
			return d;
		}
	};


	/**Bilateral processing on 2D or 3D volume data. Bilateral processing does non-linear weighted averaging processing within a certain window.
	 *@param distance_sigma means how large the voxel has impact on its neighbors in spatial domain. The larger it is, the more blurry the resulting image.
	 *@param value_sigma eans how large the voxel has impact on its in  range domain. The larger it is, the more blurry the resulting image.
	 *@param niter how many times to apply this processing on your data.
	 *@param half_width processing window size = (2 * half_widthh + 1) ^ 3.
	 */
	class BilateralProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);
		string get_name() const
		{
			return "bilateral";
		}

		string get_desc() const
		{
			return "Bilateral processing on 2D or 3D volume data. Bilateral processing does non-linear weighted averaging processing within a certain window. ";
		}

		static Processor *NEW()
		{
			return new BilateralProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("distance_sigma", EMObject::FLOAT, "means how large the voxel has impact on its neighbors in spatial domain. The larger it is, the more blurry the resulting image.");
			d.put("value_sigma", EMObject::FLOAT, "means how large the voxel has impact on its in  range domain. The larger it is, the more blurry the resulting image.");
			d.put("niter", EMObject::INT, "how many times to apply this processing on your data.");
			d.put("half_width", EMObject::INT, "processing window size = (2 * half_widthh + 1) ^ 3.");
			return d;
		}
	};

	/**Base class for normalization processors. Each specific normalization processor needs to define how to calculate mean and how to calculate sigma.
	 */
	class NormalizeProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		static string get_group_desc()
		{
			return "Base class for normalization processors. Each specific normalization processor needs to define how to calculate mean and how to calculate sigma.";
		}
		
	  protected:
		virtual float calc_sigma(EMData * image) const;
		virtual float calc_mean(EMData * image) const = 0;
	};

	/**Normalize an image so its vector length is 1.0.
	 */
	class NormalizeUnitProcessor:public NormalizeProcessor
	{
	  public:
		string get_name() const
		{
			return "normalize.unitlen";
		}

		static Processor *NEW()
		{
			return new NormalizeUnitProcessor();
		}

		string get_desc() const
		{
			return "Normalize an image so its vector length is 1.0.";
		}

	  protected:
		float calc_sigma(EMData * image) const;
		float calc_mean(EMData * image) const;
	};
	
	inline float NormalizeUnitProcessor::calc_mean(EMData *) const { return 0; }

	/**Normalize an image so its elements sum to 1.0 (fails if mean=0)
	 */
	class NormalizeUnitSumProcessor:public NormalizeProcessor
	{
	  public:
		string get_name() const
		{
			return "normalize.unitsum";
		}

		static Processor *NEW()
		{
			return new NormalizeUnitSumProcessor();
		}

		string get_desc() const
		{
			return "Normalize an image so its elements sum to 1.0 (fails if mean=0)";
		}

	  protected:
		float calc_sigma(EMData * image) const;
		float calc_mean(EMData * image) const;
	};
	
	inline float NormalizeUnitSumProcessor::calc_mean(EMData *) const { return 0; }

		
	/**do a standard normalization on an image.
	 */
	class NormalizeStdProcessor:public NormalizeProcessor
	{
	  public:
		string get_name() const
		{
			return "normalize";
		}

		static Processor *NEW()
		{
			return new NormalizeStdProcessor();
		}

		string get_desc() const
		{
			return "do a standard normalization on an image.";
		}

	  protected:
		float calc_mean(EMData * image) const;
	};

	/**Uses a 1/0 mask defining a region to use for the zero-normalization.if no_sigma is 1, standard deviation not modified.
	 *@param mask
	 *@param no_sigma
	 */
	class NormalizeMaskProcessor:public NormalizeProcessor
	{
	  public:
		string get_name() const
		{
			return "normalize.mask";
		}

		string get_desc() const
		{
			return "Uses a 1/0 mask defining a region to use for the zero-normalization.if no_sigma is 1, standard deviation not modified.";
		}

		static Processor *NEW()
		{
			return new NormalizeMaskProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("mask", EMObject::EMDATA, "the 1/0 mask defining a region to use for the zero-normalization");
			d.put("no_sigma", EMObject::INT, "if this flag is zero, only average under the mask will be substracted. set this flag to 1, standard deviation not modified");
			return d;
		}

	  protected:
		float calc_sigma(EMData * image) const;
		float calc_mean(EMData * image) const;
	};

	/**normalizes an image, mean value equals to edge mean.
	 */
	class NormalizeEdgeMeanProcessor:public NormalizeProcessor
	{
	  public:
		string get_name() const
		{
			return "normalize.edgemean";
		}

		static Processor *NEW()
		{
			return new NormalizeEdgeMeanProcessor();
		}

		string get_desc() const
		{
			return "normalizes an image, mean value equals to edge mean.";
		}

	  protected:
		float calc_mean(EMData * image) const;
	};

	/**normalizes an image, mean value equals to mean of 2 pixel circular border.
	 */
	class NormalizeCircleMeanProcessor:public NormalizeProcessor
	{
	  public:
		string get_name() const
		{
			return "normalize.circlemean";
		}

		static Processor *NEW()
		{
			return new NormalizeCircleMeanProcessor();
		}

		string get_desc() const
		{
			return "normalizes an image, mean value equals to mean of 2 pixel circular border.";
		}

	  protected:
		float calc_mean(EMData * image) const;
	};

	/**normalizes an image, uses 2 pixels on left and right edge
	 */
	class NormalizeLREdgeMeanProcessor:public NormalizeProcessor
	{
	  public:
		string get_name() const
		{
			return "normalize.lredge";
		}

		static Processor *NEW()
		{
			return new NormalizeLREdgeMeanProcessor();
		}

		string get_desc() const
		{
			return "normalizes an image, uses 2 pixels on left and right edge";
		}

	  protected:
		float calc_mean(EMData * image) const;
	};

	/**normalizes an image. mean -> (maxval-minval)/2; std dev = (maxval+minval)/2;
	 */
	class NormalizeMaxMinProcessor:public NormalizeProcessor
	{
	  public:
		string get_name() const
		{
			return "normalize.maxmin";
		}

		static Processor *NEW()
		{
			return new NormalizeMaxMinProcessor();
		}

		string get_desc() const
		{
			return "normalizes an image. mean -> (maxval-minval)/2; std dev = (maxval+minval)/2;";
		}

	  protected:
		float calc_sigma(EMData * image) const;
		float calc_mean(EMData * image) const;
	};

	/**normalizes each row in the image individually
	 */
	class NormalizeRowProcessor:public Processor
	{
	  public:
		string get_name() const
		{
			return "normalize.rows";
		}

		static Processor *NEW()
		{
			return new NormalizeRowProcessor();
		}

		string get_desc() const
		{
			return "normalizes each row in the image individually";
		}

		void process_inplace(EMData * image);
	};

	/**multiply 'this' by a constant so it is scaled to the signal in 'to'.keepzero will exclude zero values, and keep them at zero in the result.
	 *@param noisy
	 *@param keepzero
	 *@param invert
	 *@param mult
	 *@param add
	 */
	class NormalizeToStdProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "normalize.toimage";
		}

		static Processor *NEW()
		{
			return new NormalizeToStdProcessor();
		}

		string get_desc() const
		{
			return "multiply 'this' by a constant so it is scaled to the signal in 'to'.keepzero will exclude zero values, and keep them at zero in the result.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("noisy", EMObject::EMDATA,"Image to normalize to");
			d.put("keepzero", EMObject::INT,"Ignore zeroes");
			d.put("invert", EMObject::INT,"Noisy is the less noisy image");
			d.put("sigmax",EMObject::FLOAT,"Values greater or less than sigma from zero will be excluded from the normalization");
			d.put("mult", EMObject::FLOAT);
			d.put("add", EMObject::FLOAT);
			return d;
		}
	};

	/**Multiply this image by a constant so it is scaled to the signal in 'noisyfile'
	 * @param noisyfile 
	 * @param keepzero  exclude zero values
	 *@param invert
	 *@param mult
	 *@param add dfgd
	 */
	class NormalizeToFileProcessor:public NormalizeToStdProcessor
	{
	  public:
		string get_name() const
		{
			return "normalize.tofile";
		}

		static Processor *NEW()
		{
			return new NormalizeToFileProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("noisyfile", EMObject::STRING);
			d.put("keepzero", EMObject::INT, "exclude zero values");
			d.put("invert", EMObject::INT);
			d.put("mult", EMObject::FLOAT);
			d.put("add", EMObject::FLOAT);
			return d;
		}

		string get_desc() const
		{
			return "Multiply this image by a constant so it is scaled to the signal in 'noisyfile'";
		}
		
	};

	/**use least square method to normalize
	 * @param to 
	 * @param low_threshold 
	 * @param high_threshold 
	 */
	class NormalizeToLeastSquareProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "normalize.toimage.lsq";
		}

		static Processor *NEW()
		{
			return new NormalizeToLeastSquareProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("to", EMObject::EMDATA);
			d.put("low_threshold", EMObject::FLOAT);
			d.put("high_threshold", EMObject::FLOAT);
			return d;
		}

		string get_desc() const
		{
			return "use least square method to normalize";
		}
		
	};

	/**makes image circularly symmetric.
	 */
	class RadialAverageProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "math.radialaverage";
		}

		static Processor *NEW()
		{
			return new RadialAverageProcessor();
		}

		string get_desc() const
		{
			return "makes image circularly symmetric.";
		}
		
	};

	/**subtracts circularly symmetric part of an image.
	 */
	class RadialSubstractProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "math.radialsubtract";
		}

		static Processor *NEW()
		{
			return new RadialSubstractProcessor();
		}

		string get_desc() const
		{
			return "subtracts circularly symmetric part of an image.";
		}
		
	};


	/** flip an image around an axis
	 * @param axis  'x', 'y', or 'z' axis. 'x' means horizonal flip; 'y' means vertical flip;
	 */
	class FlipProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "xform.flip";
		}

		static Processor *NEW()
		{
			return new FlipProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("axis", EMObject::STRING, "'x', 'y', or 'z' axis. 'x' means horizonal flip; 'y' means vertical flip;");
			return d;
		}

		string get_desc() const
		{
			return "flip an image around an axis.";
		}
		
	};
	
	/*class FlipProcessor1:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "xform.flip1";
		}

		static Processor *NEW()
		{
			return new FlipProcessor1();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("axis", EMObject::STRING, "'x', 'y', or 'z' axis. 'x' means horizonal flip; 'y' means vertical flip;");
			return d;
		}

		string get_desc() const
		{
			return "flip an image around an axis.";
		}
		
	};*/

	/** add noise to an image
	 * @param noise 
	 */
	class AddNoiseProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "math.addnoise";
		}

		static Processor *NEW()
		{
			return new AddNoiseProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("noise", EMObject::FLOAT);
			d.put("seed", EMObject::INT, "seed for random number generator");
			return d;
		}

		string get_desc() const
		{
			return "add noise to an image";
		}

	  protected:
		virtual float get_sigma(EMData *)
		{
			return 1.0;
		}
	};

	/** add sigma noise
	 */
	class AddSigmaNoiseProcessor:public AddNoiseProcessor
	{
	  public:
		string get_name() const
		{
			return "math.addsignoise";
		}

		static Processor *NEW()
		{
			return new AddSigmaNoiseProcessor();
		}

		string get_desc() const
		{
			return "add sigma noise.";
		}
		
	  protected:
		float get_sigma(EMData * image);
	};

	/**add random noise
	 *@param n 
	 *@param x0 
	 *@param dx 
	 *@param y 
	 *@param interpolation 
	 */
	class AddRandomNoiseProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "addspectralnoise";
		}

		static Processor *NEW()
		{
			return new AddRandomNoiseProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("n", EMObject::INT);
			d.put("x0", EMObject::FLOAT);
			d.put("dx", EMObject::FLOAT);
			d.put("y", EMObject::FLOATARRAY);
			d.put("interpolation", EMObject::INT);
			d.put("seed", EMObject::INT, "seed for random number generator");
			return d;
		}

		string get_desc() const
		{
			return "add random noise.";
		}

	};
	
	/** Undo the effects of the FourierToCenterProcessor
	 * @author David Woolford <woolford@bcm.edu> 
	 * @date October 2007
	 * @ingroup tested3c
	 */
	class FourierToCornerProcessor:public Processor
	{
		public:
			/** Fourier origin shift the image in the backwards direction
			* Should only be called after the application of FourierToCenterProcessor
			* @param image the image to operate on
			* @exception ImageFormatException if the image is not complex
			 */
			void process_inplace(EMData * image);

			string get_name() const
			{
				return "xform.fourierorigin.tocorner";
			}

			static Processor *NEW()
			{
				return new FourierToCornerProcessor();
			}

			string get_desc() const
			{
				return "Undoes the xform.fourierorigin.tocenter processor";
			}
	};
			

	/** Translates the origin in Fourier space from the corner to the center in y and z
	 * Handles 2D and 3D, and handles all combinations of even and oddness
	 * Typically you call this function after Fourier transforming a real space image.
	 * After this you operate on the Fourier image in convenient format, then
	 * you call FourierToCornerProcessor (above) and then inverse FT to get to the
	 * original image
	 *
	 * @author David Woolford <woolford@bcm.edu> 
	 * @date October 2007
	 * @ingroup tested3c
	 */
	class FourierToCenterProcessor:public Processor
	{
		public:
			/** Fourier origin shift the image in the forward direction
			*
			* @param image the image to operate on
			* @exception ImageFormatException if the image is not complex
			*/
			void process_inplace(EMData * image);

			string get_name() const
			{
				return "xform.fourierorigin.tocenter";
			}

			static Processor *NEW()
			{
				return new FourierToCenterProcessor();
			}

			string get_desc() const
			{
				return "Translates the origin in Fourier space from the corner to the center in y and z - works in 2D and 3D";
			}

	};
	
	/**
	 * This class is abstract. It contains functionality common to the PhaseToCenterProcessor and PhaseToCornerProcessor
	 * processors
	 *
	 * @author David Woolford <woolford@bcm.edu> 
	 * @date October 2007
	 * @ingroup tested3c
	 * though the testing of this processor is really implicit
	 */
	class Phase180Processor:public Processor
	{
		protected:
			/** swap_corners_180 - works on 2D and 3D images
			* 
			* Implements the conventional 180 degree phase shift required to put the center of the image
			* at the bottom left of the image - is used in conjunction with swap_central_slices_180
			* if any of the image dimensions are odd, but by itself will perform the entire operation on even 
			* images. This functions is never called by anyone except for the PhaseToCenterProcessor and
			* PhaseToCornerProcessor classes.
			* Highly specialised function to handle all cases of even and oddness
			* @param image the image to be operated upon
			* @exception ImageDimensionException if the image is 1D
			* @exception NullPointerException if the image is null
			*/
			void swap_corners_180(EMData * image);
			
			/** swap_central_slices_180 - works on 2D and 3D images
			 * 
			 * swaps pixels values in central slices, only required if the image has one or more
			 * odd dimensions. Should be used striclty in conjunction with swap_central_slices_180 function
			 * and never called by anyone except for PhaseToCenterProcessor and PhaseToCornerProcessor
			 * classes.
			 * Highly specialised function to handle all cases of even and oddness
			 * @param image the image to be operated upon
			 * @exception ImageDimensionException if the image is 1D
			 * @exception NullPointerException if the image is null
			 */
			void swap_central_slices_180(EMData * image);
			
			/** fourier_phaseshift180 - fourier phase shift by 180
			 * this function is called internally if the argument to the process_inplace function
			 * is complex.
			 * @param image the image to be operated upon
			 * @exception ImageFormatException if the image is not in complex format
			 */
			void fourier_phaseshift180(EMData * image);

	};
	
	/** Translates a cornered image to the center
	 * Undoes the PhaseToCornerProcessor
	 *
	 * works for 1D, 2D and 3D images, for all combinations of even and oddness
	 *
	 * @author David Woolford <woolford@bcm.edu> 
	 * @date October 2007
	 * @ingroup tested3c
	 */
	class PhaseToCenterProcessor:public Phase180Processor
	{
		public:
			void process_inplace(EMData * image);

			string get_name() const
			{
				return "xform.phaseorigin.tocenter";
			}

			static Processor *NEW()
			{
				return new PhaseToCenterProcessor();
			}
		
			string get_desc() const
			{
				return "Undoes the effect of the xform.phaseorigin.tocorner processor";
			}

	};
	
	/** Translates a centered image to the corner
	 * works for 1D, 2D and 3D images, for all combinations of even and oddness
	 *
	 * @author David Woolford <woolford@bcm.edu> 
	 * @date October 2007
	 * @ingroup tested3c
	 */
	class PhaseToCornerProcessor:public Phase180Processor
	{
		public:
			void process_inplace(EMData * image);

			string get_name() const
			{
				return "xform.phaseorigin.tocorner";
			}

			static Processor *NEW()
			{
				return new PhaseToCornerProcessor();
			}
		
			string get_desc() const
			{
				return "Translates a centered image to the corner in a forward fashion";
			}

	};

	/**Attempts to automatically mask out the particle, excluding other particles in the box, etc.
	 * @param threshold  runs from ~ -2 to 2, negative numbers for dark protein and positive numbers for light protein (stain).
	 * @param filter  is expressed as a fraction of the fourier radius.
	 */
	class AutoMask2DProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "mask.auto2d";
		}

		static Processor *NEW()
		{
			return new AutoMask2DProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("threshold", EMObject::FLOAT, "runs from ~ -2 to 2, negative numbers for dark protein and positive numbers for light protein (stain).");
			d.put("filter", EMObject::FLOAT, "is expressed as a fraction of the fourier radius.");
			return d;
		}
		
		string get_desc() const
		{
			return "Attempts to automatically mask out the particle, excluding other particles in the box, etc.";
		}

	};

	/**Tries to mask out only interesting density
	 * @param threshold1 
	 * @param threshold2 
	 */
	class AutoMask3DProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "mask.auto3d.thresh";
		}

		static Processor *NEW()
		{
			return new AutoMask3DProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("threshold1", EMObject::FLOAT);
			d.put("threshold2", EMObject::FLOAT);
			return d;
		}

		string get_desc() const
		{
			return "Tries to mask out only interesting density";
		}

		static void search_nearby(float *dat, float *dat2, int nx, int ny, int nz, float thr);
		static void fill_nearby(float *dat2, int nx, int ny, int nz);
	};

	/** Tries to mask out only interesting density
	 * @param radius 
	 * @param threshold 
	 * @param nshells 
	 */
	class AutoMask3D2Processor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "mask.auto3d";
		}

		static Processor *NEW()
		{
			return new AutoMask3D2Processor();
		}

		string get_desc() const
		{
			return "Tries to mask out only interesting density";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("radius", EMObject::INT);
			d.put("threshold", EMObject::FLOAT);
			d.put("nshells", EMObject::INT);
			return d;
		}
	};

	/**Add additional shells/rings to an existing 1/0 mask image
	 * @param nshells   number of shells to add
	*/
	class AddMaskShellProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "mask.addshells";
		}

		string get_desc() const
		{
			return "Add additional shells/rings to an existing 1/0 mask image";
		}

		static Processor *NEW()
		{
			return new AddMaskShellProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("nshells", EMObject::INT, "number of shells to add");
			return d;
		}
	};

	/**ToMassCenterProcessor centers image at center of mass, ignores old dx, dy.
	 * @param int_shift_only 
	 */
	class ToMassCenterProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "xform.centerofmass";
		}

		static Processor *NEW()
		{
			return new ToMassCenterProcessor();
		}

		string get_desc() const
		{
			return "ToMassCenterProcessor centers image at center of mass, ignores old dx, dy.";
		}
		
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("int_shift_only", EMObject::INT);
			return d;
		}
	};

	/**Center image using CCF with 180 degree rotation.
	 * @param is3d 
	 */
	class ACFCenterProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "xform.centeracf";
		}

		static Processor *NEW()
		{
			return new ACFCenterProcessor();
		}

		string get_desc() const
		{
			return "Center image using CCF with 180 degree rotation.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("is3d", EMObject::INT);
			return d;
		}
	};

	/** Processor the images by the estimated SNR in each image.if parameter 'wiener' is 1, then wiener processor the images using the estimated SNR with CTF amplitude correction.
	 * @param wiener 
	 * @param snrfile 
	 */
	class SNRProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "eman1.filter.snr";
		}

		static Processor *NEW()
		{
			return new SNRProcessor();
		}

		string get_desc() const
		{
			return "Processor the images by the estimated SNR in each image.if parameter 'wiener' is 1, then wiener processor the images using the estimated SNR with CTF amplitude correction.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("wiener", EMObject::INT);
			d.put("snrfile", EMObject::STRING);
			return d;
		}
	};

	/** A fourier processor specified in a 2 column text file.
	 * @param filename 
	*/
	class FileFourierProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "eman1.filter.byfile";
		}

		string get_desc() const
		{
			return "A fourier processor specified in a 2 column text file.";
		}

		static Processor *NEW()
		{
			return new FileFourierProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("filename", EMObject::STRING);
			return d;
		}
	};

	/** Identifiy the best symmetry in the given symmetry list for each pixel and then apply the best symmetry to each pixel
	 *
	 * @author Wen Jiang <wjiang@bcm.tmc.edu> 
	 * @date 2005-1-9
	 * @param sym[in] the list of symmetries to search
	 * @param thresh[in] the minimal level of symmetry to be accepted (0-1)
	 * @param output_symlabel[in] if output the symmetry label map in which the pixel value is the index of symmetry in the symmetry list
	 * @param symlabel_map[out] the optional return map when output_symlabel=1
	 */

	class SymSearchProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "misc.symsearch";
		}

		string get_desc() const
		{
			return "Identifiy the best symmetry in the given symmetry list for each pixel and then apply the best symmetry to each pixel.";
		}

		static Processor *NEW()
		{
			return new SymSearchProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("sym", EMObject::STRINGARRAY, "the list of symmetries to search");
			d.put("thresh", EMObject::FLOAT, "the minimal level of symmetry to be accepted (0-1)");
			d.put("output_symlabel", EMObject::INT, "if output the symmetry label map in which the pixel value is the index of symmetry in the symmetry list");
			d.put("symlabel_map", EMObject::EMDATA, "the optional return map when output_symlabel=1");
			return d;
		}
	};

	/**This processor attempts to perform a 'local normalization' so low density and high density features will be on a more even playing field in an isosurface display. threshold is an isosurface threshold at which all desired features are visible, radius is a normalization size similar to an lp= value.
	 *@param threshold
	 *@param radius
	 *@param apix
	 */
	class LocalNormProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "misc.localnorm";
		}

		static Processor *NEW()
		{
			return new LocalNormProcessor();
		}

		string get_desc() const
		{
			return "This processor attempts to perform a 'local normalization' so low density and high density features will be on a more even playing field in an isosurface display. threshold is an isosurface threshold at which all desired features are visible, radius is a normalization size similar to an lp= value.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("threshold", EMObject::FLOAT);
			d.put("radius", EMObject::FLOAT);
			d.put("apix", EMObject::FLOAT);
			return d;
		}
	};

	/**Multiplies the image by the specified file using pixel indices. The images must be same size. If 'ismaskset=' is 1, it will take a file containing a set of masks and apply the first mask to the image.
	 *@param filename
	 *@param ismaskset
	 */
	class IndexMaskFileProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "mask.fromfile";
		}

		static Processor *NEW()
		{
			return new IndexMaskFileProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("filename", EMObject::STRING);
			d.put("ismaskset", EMObject::INT);
			return d;
		}

		string get_desc() const
		{
			return "Multiplies the image by the specified file using pixel indices. The images must be same size. If 'ismaskset=' is 1, it will take a file containing a set of masks and apply the first mask to the image.";
		}

	};

	/**Multiplies the image by the specified file using pixel coordinates instead of pixel indices. The images can be different size.
	 *@param filename
	 */
	class CoordinateMaskFileProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "mask.fromfile.sizediff";
		}

		static Processor *NEW()
		{
			return new CoordinateMaskFileProcessor();
		}

		string get_desc() const
		{
			return "Multiplies the image by the specified file using pixel coordinates instead of pixel indices. The images can be different size.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("filename", EMObject::STRING);
			return d;
		}
	};

	/**'paints' a circle into the image at x,y,z with values inside r1 set to v1, values between r1 and r2 will be set to a
		value between v1 and v2, and values outside r2 will be unchanged
	 */
	class PaintProcessor:public Processor
	{
	  public:
		PaintProcessor():x(0), y(0), z(0),r1(0), v1(0.0), r2(0), v2(0.0)
		{
		}

		string get_name() const
		{
			return "mask.paint";
		}

		static Processor *NEW()
		{
			return new PaintProcessor();
		}

		string get_desc() const
		{
			return "Paints a circle with a decaying edge into the image. r<r1 -> v1, r1<r<r2 -> (v1,v2), r>r2 unchanged";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("x", EMObject::INT,"Center of circle");
			d.put("y", EMObject::INT);
			d.put("z", EMObject::INT);
			d.put("r1", EMObject::INT,"First radius");
			d.put("v1", EMObject::FLOAT,"Inner value");
			d.put("r2", EMObject::INT,"Second radius");
			d.put("v2", EMObject::FLOAT,"Outer Value");
			return d;
		}

		void set_params(const Dict & new_params)
		{
			params = new_params;

			if (params.has_key("x")) x = params["x"];
			if (params.has_key("y")) y = params["y"];
			if (params.has_key("z")) z = params["z"];
			if (params.has_key("r1")) r1 = params["r1"];
			if (params.has_key("r2")) r2 = params["r2"];
			if (params.has_key("v1")) v1 = params["v1"];
			if (params.has_key("v2")) v2 = params["v2"];
		}

		protected:
		void process_inplace(EMData *image);

		int x,y,z,r1;
		float v1;
		int r2;
		float v2;

	};


	/**Sets the structure factor based on a 1D x/y text file.
	 *@param filename
	 */
	class SetSFProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "misc.setpowspec";
		}

		string get_desc() const
		{
			return "Sets the structure factor based on a 1D x/y text file.";
		}

		static Processor *NEW()
		{
			return new SetSFProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("filename", EMObject::STRING);
			return d;
		}
	};

	/**Smart mask processor
	 *@param mask
	 */
	class SmartMaskProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "mask.smart";
		}

		static Processor *NEW()
		{
			return new SmartMaskProcessor();
		}

		string get_desc() const
		{
			return "Smart mask processor.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("mask", EMObject::FLOAT);
			return d;
		}
	};

	/**Iterative expansion of a binary mask, val1 is number of pixels to expand, if val2!=0 will make a soft Gaussian edge starting after val2 pixels.
	 */
	class IterBinMaskProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "mask.addshells.gauss";
		}

		string get_desc() const
		{
			return "Iterative expansion of a binary mask, val1 is number of pixels to expand, if val2!=0 will make a soft Gaussian edge starting after val2 pixels.";
		}

		static Processor *NEW()
		{
			return new IterBinMaskProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("val1", EMObject::FLOAT, "number of pixels to expand");
			d.put("val2", EMObject::FLOAT, "if val2!=0 will make a soft Gaussian edge starting after val2 pixels");
			return d;
		}
	};

	/**This is a group of 'processor' used to create test image.
	 */
	class TestImageProcessor : public Processor
	{
	public:
		static string get_group_desc()
		{
			return "This is a group of 'processors' used to create test image.";
		}
	
	protected:
		void preprocess(const EMData * const image);
		int nx, ny, nz; //this is the size of the source image
	};
	
	/**Replace a source image as a strict Gaussian 
	 *@param sigma sigma value for this Gaussian blob
	 */
	class TestImagePureGaussian : public TestImageProcessor
	{
	public:
		void process_inplace(EMData * image);
		
		string get_name() const
		{
			return "testimage.puregaussian";
		}
		
		string get_desc() const
		{
			return "Replace a source image as a strict Gaussian ";
		}
		
		static Processor * NEW()
		{
			return new TestImagePureGaussian();
		}
		
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("x_sigma", EMObject::FLOAT, "sigma value for this Gaussian blob on x direction");
			d.put("y_sigma", EMObject::FLOAT, "sigma value for this Gaussian blob on y direction");
			d.put("z_sigma", EMObject::FLOAT, "sigma value for this Gaussian blob on z direction");
			d.put("x_center", EMObject::FLOAT, "center for this Gaussian blob on x direction" );
			d.put("y_center", EMObject::FLOAT, "center for this Gaussian blob on y direction" );
			d.put("z_center", EMObject::FLOAT, "center for this Gaussian blob on z direction" );
			return d;
		}
	};
	
	/**Make an image useful for tomographic reconstruction testing
	 * this is a 3D phantom image based on the 2D phantom described in
	 * Delaney and Bresler, "Globally convergent edge-preserving regularized reconstruction: An application to limited-angle tomography". IEEE
	 * Transactions on Image Processing, 7(2), Feb 1998, 204-221.
	 * @author David Woolford
	 * @date November 2007
	 */
	class TestTomoImage : public TestImageProcessor
	{
		public:
			/** Make a useful tomographic phantom image 
			 * @param image the image to operate upon
			 */
			void process_inplace(EMData * image);
		
			string get_name() const
			{
				return "testimage.tomo.objects";
			}
		
			string get_desc() const
			{
				return "Make an image consisting various objects, useful for tomographic testing";
			}
		
			static Processor * NEW()
			{
				return new TestTomoImage();
			}
		private:
			void insert_solid_ellipse( EMData* image, const Region& region, const float& value, const Transform3D* const t3d = NULL );
			void insert_hollow_ellipse( EMData* image, const Region& region, const float& value, const int& radius, const Transform3D* const t3d = NULL );
			void insert_rectangle( EMData* image, const Region& region, const float& value, const Transform3D* const t3d = NULL );
	};
	
	/** Put a gradient in the image of the form y = mx+b : "x" is a string indicating any of the image axes, i.e., x,y or z.
	 * @author David Woolford <woolford@bcm.edu>
	 * @date 01/10/2008
	 * @param axis The axis the will be used to determine pixel values. Must be x,y or z. Default is x
	 * @param m m in the equation m*axis+b. Default is 1.0
	 * @param b b in the equation m*axis+b. Default is 0.0
	*/
	class TestImageGradient : public TestImageProcessor
	{
		public:
			void process_inplace(EMData * image);
		
			string get_name() const
			{
				return "testimage.gradient";
			}
		
			string get_desc() const
			{
				return "Make a gradient image of the form y=mx+b, where x is any of the image axes.";
			}
		
			static Processor * NEW()
			{
				return new TestImageGradient();
			}
		
			TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("axis", EMObject::STRING, "The axis the will be used to determine pixel values. Must be x,y or z");
				d.put("m", EMObject::FLOAT, "m in the equation m*axis+b. Default is 1.0");
				d.put("b", EMObject::FLOAT, "b in the equation m*axis+b. Default is 0.0");
				return d;
			}
	};
	
	/**Make an image consisting of a single cross, with lines
	 * going in the axial directions, intersecting at the origin.
	 * @author David Woolford <woolford@bcm.edu>
	 * @date October 2007
	 *@param radius the radial length of the lines from the origin
	 *@param fill the value to assign to pixels made non zero
	 */
	class TestImageAxes : public TestImageProcessor
	{
		public:
			/** Make an image where the axes (where x,y and z=0) are some
			 * nono zero value
			* @param image the image to operate upon
			 */
			void process_inplace(EMData * image);
		
			string get_name() const
			{
				return "testimage.axes";
			}
		
			string get_desc() const
			{
				return "Make an image consisting of a single cross";
			}
		
			static Processor * NEW()
			{
				return new TestImageAxes();
			}
		
			TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("int", EMObject::FLOAT, "radius of the lines emanating from the origin");
				d.put("fill", EMObject::FLOAT, "value to make non-zero pixels");
				return d;
			}
	};
	
	/** Replace a source image as a Gaussian Blob
	 *@param sigma sigma value for this Gaussian blob
	 *@param axis specify a major axis for asymmetric features
	 *@param c distance between focus and the center of an ellipse
	 */
	class TestImageGaussian : public TestImageProcessor
	{
	public:
		void process_inplace(EMData * image);
		
		string get_name() const
		{
			return "testimage.gaussian";
		}
		
		string get_desc() const
		{
			return "Replace a source image as a Gaussian Blob";
		}
		
		static Processor * NEW()
		{
			return new TestImageGaussian();
		}
		
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("sigma", EMObject::FLOAT, "sigma value for this Gaussian blob");
			d.put("axis", EMObject::STRING, "specify a major axis for asymmetric features");
			d.put("c", EMObject::FLOAT, "distance between focus and the center of an ellipse");
			return d;
		}
	};
	
	/**Replace a source image with a lumpy S-curve used for alignment testing
	 */
	class TestImageScurve : public TestImageProcessor
	{
	public:
		void process_inplace(EMData * image);
		
		string get_name() const
		{
			return "testimage.scurve";
		}
		
		string get_desc() const
		{
			return "Replace a source image with a lumpy S-curve used for alignment testing";
		}
		
		static Processor * NEW()
		{
			return new TestImageScurve();
		}
		
		TypeDict get_param_types() const
		{
			TypeDict d;
			return d;
		}
	};
	

	/**Replace a source image as a sine wave in specified wave length
	 *@param wavelength this value is the d in function |sin(x/d)|
	 *@param axis specify a major axis for asymmetric features
	 *@param c distance between focus and the center of an ellipse
	 *@param phase (optional)phase for sine wave, default is 0
	 */
	class TestImageSinewave : public TestImageProcessor
	{
	public:
		void process_inplace(EMData * image);
		
		string get_name() const
		{
			return "testimage.sinewave";
		}
		
		string get_desc() const
		{
			return "Replace a source image as a sine wave in specified wave length";
		}
		
		static Processor * NEW()
		{
			return new TestImageSinewave();
		}
		
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("wave_length", EMObject::FLOAT, "wave_length in equation sin(x*2*PI/wave_length - phase*180/PI)");
			d.put("axis", EMObject::STRING, "(optional) specify a major axis for asymmetric features, default x axis");
			d.put("phase", EMObject::FLOAT, "(optional) the phase in equation sin(x*2*PI/wave_length - phase*180/PI)");
			d.put("az", EMObject::FLOAT, "(optional) angle in degree. for 2D image, this is the rotated angle of the image, \
												in 3D image, it's az for euler angle. default is zero");
			d.put("alt", EMObject::FLOAT, "(optional) angle in degree. only in 3D case, alt for euler angle, default is zero"); 
			d.put("phi", EMObject::FLOAT, "(optional) angle in degree. only in 3D case, phi for euler angle, default is zero"); 
			return d;
		}
	};
	
	/**Replace a source image as a circular sine wave in specified wave length
	 *@param wave_length float this value is the d in function |sin(x/d)|
	 *@param axis char (optional) specify a major axis for asymmetric features
	 *@param c distance (optional) float between focus and the center of an ellipse
	 *@param phase degree (optional) phase for sine wave, default is 0
	 */
	class TestImageSinewaveCircular : public TestImageProcessor
	{
	public:
		void process_inplace(EMData * image);
		
		string get_name() const
		{
			return "testimage.sinewave.circular";
		}
		
		string get_desc() const
		{
			return "Replace a source image as a circular sine wave in specified wave length";
		}
		
		static Processor * NEW()
		{
			return new TestImageSinewaveCircular();
		}
		
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("wave_length", EMObject::FLOAT, "(required)this value is the d in function |sin(x/d)|, unit: pixel");
			d.put("axis", EMObject::STRING, "specify a major axis for asymmetric features");
			d.put("c", EMObject::FLOAT, "distance between focus and the center of an ellipse");
			d.put("phase", EMObject::FLOAT, "(optional)phase for sine wave, default is 0");
			return d;
		}
	};
	
	/**Replace a source image as a square or cube depends on 2D or 3D of the source image
	 *@param edge_length edge length of the square or cube
	 *@param axis specify a major axis for asymmetric features
	 *@param odd_edge edge length for the asymmetric axis
	 *@param fill answer 'yes' or 'no' to specify if it's filled or hollow, default filled
	 */
	class TestImageSquarecube : public TestImageProcessor
	{
	public:
		void process_inplace(EMData * image);
		
		string get_name() const
		{
			return "testimage.squarecube";
		}
		
		string get_desc() const
		{
			return "Replace a source image as a square or cube depends on 2D or 3D of the source image";
		}
		
		static Processor * NEW()
		{
			return new TestImageSquarecube();
		}
		
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("edge_length", EMObject::FLOAT, "edge length of the square or cube, unit: pixel");
			d.put("axis", EMObject::STRING, "specify a major axis for asymmetric features");
			d.put("odd_edge", EMObject::FLOAT, "edge length for the asymmetric axis");
			d.put("fill", EMObject::INT, "Flag indicating if image is filled, default filled, 1 for filled, 0 for blank");
			return d;
		}
	};
	
	/**Replace a source image as a circle or sphere depends on 2D or 3D of the source image
	 *@param radius radius of circle or sphere
	 *@param axis specify a major axis for asymmetric features
	 *@param c distance between focus and the center of an ellipse
	 *@param fill answer 'yes' or 'no' to specify if it's filled or hollow, default filled
	 */
	class TestImageCirclesphere : public TestImageProcessor
	{
	public:
		void process_inplace(EMData * image);
		
		string get_name() const
		{
			return "testimage.circlesphere";
		}
		
		string get_desc() const
		{
			return "Replace a source image as a circle or sphere depends on 2D or 3D of the source image";
		}
		
		static Processor * NEW()
		{
			return new TestImageCirclesphere();
		}
		
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("radius", EMObject::FLOAT, "radius of circle or sphere, unit: pixel");
			d.put("axis", EMObject::STRING, "specify a major axis for asymmetric features");
			d.put("c", EMObject::FLOAT, "distance between focus and the center of an ellipse");
			d.put("fill", EMObject::INT, "Flag indicating if image is filled, default filled, 1 for filled, 0 for blank.");
			return d;
		}
	};
	
	/**Replace a source image as a uniform random noise, random number generated from gsl_rng_mt19937, 
	 * the pixel value is from 0 to 1, [0, 1)
	 */
	class TestImageNoiseUniformRand : public TestImageProcessor
	{
	public:
		void process_inplace(EMData * image);
		
		string get_name() const
		{
			return "testimage.noise.uniform.rand";
		}
		
		string get_desc() const
		{
			return "Replace a source image as a uniform random noise, random number generated from gsl_rng_mt19937, the pixel value is [0, 1)";
		}
		
		static Processor * NEW()
		{
			return new TestImageNoiseUniformRand();
		}
		
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("seed", EMObject::INT, "seed for random number generator");
			return d;
		}
	};
	
	/**Replace a source image with gaussian distributed random noise
	 * If you don't provide a seed at all, it should be seeded using the best available
	 * source of randomness( time(0) in this implementation).
	 * The testimage classes using random numbers should take an int 'seed'
	 * parameter. If this parameter is provided, it will be cast into an unsigned int. 
	 * This will permit initialization to a known state if desired.
	 */
	class TestImageNoiseGauss : public TestImageProcessor
	{
	public:
		void process_inplace(EMData * image);
		
		string get_name() const
		{
			return "testimage.noise.gauss";
		}
		
		string get_desc() const
		{
			return "Replace a source image as a random noise, the random value is gaussian distributed";
		}
		
		static Processor * NEW()
		{
			return new TestImageNoiseGauss();
		}
		
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("sigma", EMObject::FLOAT, "sigma value of gausian distributed noise, default is 0.5");
			d.put("mean", EMObject::FLOAT, "mean value of gausian distributed noise, default is zero.");
			d.put("seed", EMObject::INT, "the seed for random number generator, default is not to reseed. Uses random().");

			return d;
		}
	};
	
	/** Replace a source image with a cylinder */
	class TestImageCylinder : public TestImageProcessor
	{
	public:
		void process_inplace(EMData * image);
		
		string get_name() const
		{
			return "testimage.cylinder";
		}	
		
		string get_desc() const
		{
			return "Replace a source image as a cylinder";
		}
		
		static Processor * NEW()
		{
			return new TestImageCylinder();
		}
		
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("radius", EMObject::FLOAT, "radius for the cylinder");
			d.put("height", EMObject::FLOAT, "height for the cylinder, by default it's the nz");
			
			return d;
		}
	};

	/** Try to normalize the 4 quadrants of a CCD image
	 * @author Deepy Mann <dsmann@bcm.tmc.edu> 
	 * @date 9-2005
	 * @param width number of pixels on either side of the seam to sample
	 */
	class CCDNormProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "filter.ccdnorm";
		}

		static Processor *NEW()
		{
			return new CCDNormProcessor();
		}

		string get_desc() const
		{
			return "normalize the 4 quadrants of a CCD image";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("width", EMObject::INT, "number of pixels on either side of the seam to sample");
			return d;
		}
	};

 	/** Perform a Wavelet transform using GSL
	 * @author Steve Ludtke <sludtke@bcm.edu>
	 * @date 10/15/2006
	 * @param type  "daub", "harr", or "bspl"
	 * @param dir  1 for forward transform, -1 for inverse transform
	 * @param ord  for Daubechies (4,6,8,...,20), for Harr (2), for B-Splines (103, 105, 202, 204, 206, 208, 301, 303, 305 307, 309)
	 */
	class WaveletProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "basis.wavelet";
		}

		static Processor *NEW()
		{
			return new WaveletProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("type", EMObject::STRING, "'daub', 'harr' or 'bspl'");
			d.put("dir", EMObject::INT, "1 for forward transform, -1 for inverse transform");
			d.put("ord", EMObject::INT, "Daubechies (4,6,8,...,20), for Harr (2), for B-Splines (103, 105, 202, 204, 206, 208, 301, 303, 305 307, 309)");
			return d;
		}

		string get_desc() const
		{
			return "Computes the DWT (discrete wavelet transform) of an image in one of 3 possible bases";
		}
		
	};
	
	/** A processor designed specifically for tomographic tilt series data.
	 * This processors masks out 'mass' in tilted images that is not present in the zero-tilt (0 degrees) image.
	 * It does this based on the tilt angle. The tilt angle can be extracted from the image metadata (stored as the euler_alt attribute),
	 * or it may be specified explicitly (specifying the angle is the default behavior). The masked out regions at both sides of the image are set to 0 by default, 
	 * but can  also be set to the mean of the nearest non-masked data edge (in the y direction), or similarly the mean of both non-masked data
	 * edges on either side of the image. A gaussian fall-off is optional (but off by default).
	 * @author David Woolford <woolford@bcm.edu>
	 * @date 01/10/2008
	 * @param biedgemean Mutually  exclusive of edgemean. Experimental. Causes the pixels in the masked out areas to take the average value of both the left and right edge pixel strips
	 * @param edgemean Mutually  exclusive of biedgemean. Masked pixels values assume the mean edge pixel value, independently, for both sides of the image
	 * @param angle The angle that the image is, with respect to the zero tilt image
	 * @param angle_fim Read fim as 'from image metadata' - this causes the altitude angle stored in by the image object (i.e. as extracted from the header, as currently stored in memory) to be used as the angle. This overrides the angle argument
	 * @param gauss_falloff Causes the edge masking to have a smooth Gaussian fall-off - this parameter specifies how many pixels the fall-off will proceed over. Default is 0
	 * @param gauss_sigma The sigma of the Gaussian function used to smooth the edge fall-off (functional form is exp(-(pixel distance)^2/sigma^2)
	 */
	class TomoTiltEdgeMaskProcessor : public Processor
	{
	public:
		void process_inplace(EMData* image);
			
		string get_name() const
		{
			return "tomo.tiltedgemask";
		}
		
		static Processor *NEW()
		{
			return new TomoTiltEdgeMaskProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("biedgemean", EMObject::BOOL, "Mutually  exclusive of edgemean. Experimental. Causes the pixels in the masked out areas to take the average value of both the left and right edge pixel strips");
			d.put("edgemean", EMObject::BOOL, "Mutually  exclusive of biedgemean. Masked pixels values assume the mean edge pixel value, independently, for both sides of the image.");
			d.put("angle", EMObject::INT, "The angle that the image is, with respect to the zero tilt image");
			d.put("gauss_falloff",EMObject::INT, "Causes the edge masking to have a smooth Gaussian fall-off - this parameter specifies how many pixels the fall-off will proceed over. Default is 0.");
			d.put("gauss_sigma",EMObject::FLOAT,"The sigma of the Gaussian function used to smooth the edge fall-off (functional form is exp(-(pixel distance)^2/sigma^2)");
			d.put("angle_fim",EMObject::BOOL,"Read fim as 'from image metadata' - this causes the altitude angle stored in by the image object (i.e. as extracted from the header, as currently stored in memory) to be used as the angle. This overrides the angle argument");
			return d;
		}

		string get_desc() const
		{
			return "Masks the part of the image which is not present in the 0-tilt image. Masked areas can be 0 or set to the edgemean (of the nearest or both edges). Masked areas can also have a Gaussian fall-off to make the appearance smooth.";
		}

	};
	
	/** A processor that can be used to weight an image by 1/cos(angle)
	 * This processor evolved originally as an experimental tool for weighting tomographic data
	 * by the width of its cross section relative to the electron beam. The relative width
	 * can be derived using elementary trigonometry to be 1/cos(tiltangle). This processor 
	 * should hence probably be called OneOverCosineWeightingProcessor. You can specify the
	 * angle explicitly (which is the default behavior), or you can force the angle
	 * to be the altitude angle as derived from the EMData metadata. The processor could obviously
	 * be made more robust if the angle derived from the EMData header could be specified...
	 *
	 * @author David Woolford <woolford@bcm.edu>
	 * @date 02/11/2008
	 * @param angle The angle that the image is, with respect to the zero tilt image
	 * @param angle_fim Read fim as 'from image metadata' - this causes the altitude angle stored in by the image object (i.e. as extracted from the header, as currently stored in memory) to be used as the angle. This overrides the angle argument
	 */
	class TomoTiltAngleWeightProcessor : public Processor
	{
		public:
			void process_inplace(EMData* image);
			
			string get_name() const
			{
				return "tomo.tiltangleweight";
			}
		
			static Processor *NEW()
			{
				return new TomoTiltAngleWeightProcessor();
			}

			TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("angle", EMObject::INT, "The angle that the image is, with respect to the zero tilt image");
				d.put("angle_fim",EMObject::BOOL,"Read fim as 'from image metadata' - this causes the altitude angle stored in by the image object (i.e. as extracted from the header, as currently stored in memory) to be used as the angle. This overrides the angle argument");
				return d;
			}

			string get_desc() const
			{
				return "Weights the image by 1/cos(angle)";
			}

	};

	/** Perform a FFT transform by calling EMData::do_fft() and EMData::do_ift()
	 * @param dir 1 for forward transform, -1 for inverse transform, forward transform by default
	 */
	class FFTProcessor : public Processor
	{
	  public:
		void process_inplace(EMData * image);
		
		string get_name() const
		{
			return "basis.fft";
		}	
		
		static Processor *NEW()
		{
			return new FFTProcessor();
		}
		
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("dir", EMObject::INT, "1 for forward transform, -1 for inverse transform");
			return d;
		}
		
		string get_desc() const
		{
			return "Computes the DFFT (Discrete Fast Fourier Transform) of an image";
		}
	};
	
	/** Perform a multiplication of real image with a radial table
	 * @param table a radial table for multiplication
	 * @exception ImageFormatException this filter only apply to real image
	 * */
	class RadialProcessor : public Processor
	{
	public:
		void process_inplace(EMData * image);
		
		string get_name() const
		{
			return "filter.radialtable";
		}
		
		static Processor *NEW()
		{
			return new RadialProcessor();
		}
		
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("table", EMObject::FLOATARRAY, "radial float array for multiplication");
			return d;
		}
		
		string get_desc() const
		{
			return "multiply an image in real-space by a radial function";
		}
	};
	
	class HistogramBin : public Processor
	{
		public:
			HistogramBin() : default_bins(1024) {}
			
			void process_inplace(EMData * image);
		
			string get_name() const
			{
				return "histogram.bin";
			}
		
			static Processor *NEW()
			{
				return new HistogramBin();
			}
		
			TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("nbins", EMObject::INT, "The number of bins the pixel values will be compressed into");
				d.put("debug", EMObject::BOOL, "Outputs debugging information (number of pixels per bin)");
				return d;
			}
		
			string get_desc() const
			{
				return "Bins pixel values, similar to calculating a histogram. The histogram is comprised of 'nbins' bins, and the value assigned to each pixel in the bin is the midpoint of the bin's upper and lower limits. Defaults to 256 bins";
			}
		protected:
			int default_bins;
	};

#if 0

	class XYZProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return "XYZ";
		}

		static Processor *NEW()
		{
			return new XYZProcessor();
		}

		string get_desc() const
		{
			return "N/A";
		}
		
		TypeDict get_param_types() const
		{
			TypeDict d;
			return d;
		}
	};


#endif


	int multi_processors(EMData * image, vector < string > processornames);
	void dump_processors();
	map<string, vector<string> > dump_processors_list();
	map<string, vector<string> > group_processors();

	template <> Factory < Processor >::Factory();
}

#endif	//eman_filter_h__

/* vim: set ts=4 noet: */
