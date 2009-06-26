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
#include "transform.h"

#include <cfloat>
#include <climits>
#include <cstring>

using std::vector;
using std::map;
using std::string;

namespace EMAN
{
	class EMData;

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

	/** base class for Fourier filters
	 * @param cutoff_abs Processor radius in terms of Nyquist (0-.5).
	 * @param cutoff_pixels Width in Fourier pixels (0 - size()/2).
	 * @param cutoff_freq Resolution in 1/A (0 - 1 / size*apix).
	 * @param apix Override A/pix in the image header (changes x,y and z).
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
			d.put("cutoff_pixels", EMObject::FLOAT, " Width in Fourier pixels (0 - size()/2)");
			d.put("cutoff_freq", EMObject::FLOAT, "Resolution in 1/A (0 - 1 / size*apix)");
			d.put("apix", EMObject::FLOAT, " Override A/pix in the image header (changes x,y and z)");
			return d;
		}

	  protected:
		  virtual void preprocess(EMData * image) {}
		  virtual void create_radial_func(vector < float >&radial_mask) const = 0;
	};

	/**Multiplies each Fourier pixel by its amplitude
	 *@param sum Adds the weights to sum for normalization
	 *@param sqrt Weights using sqrt of the amplitude if set
	 */
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
	 *@param with The image that will convolute the other image
	 */
	class ConvolutionProcessor : public Processor
	{
		public:
			ConvolutionProcessor() {}

			string get_name() const
			{
				return "math.convolution";
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
			return "math.edge.xgradient";
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
				return "math.edge.ygradient";
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
				return "math.edge.zgradient";
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

		virtual EMData* process(const EMData * const image);

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


	/**Wiener filter based on a Ctf object either in the image header
	 *@param ctf[in] A Ctf object to use
	 *
	 *@author Steve Ludtke
	 *@date 2008/11/03
	 */
	class Wiener2DFourierProcessor:public Processor
	{
	  public:
		string get_name() const
		{
			return "filter.wiener2d";
		}

		virtual EMData* process(const EMData * const image);

		void process_inplace(EMData *image);

		void set_params(const Dict & new_params)
		{
			params = new_params;
			ctf = params["ctf"];
//			printf("%s %f\n",params.keys()[0].c_str(),lowpass);
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("ctf", EMObject::EMDATA, "Ctf object to use for Wiener filter parameters");
			return d;
		}

		static Processor *NEW()
		{
			return new Wiener2DFourierProcessor();
		}

		string get_desc() const
		{
			return "Applies a 2-D Wiener filter to an image based on its Ctf parameters";
		}

		protected:
		Ctf *ctf;
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
			if( params.has_key("lowpass") ) {
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
		  virtual void preprocess(EMData * image);
		  float lowpass;
	};

	class LinearRampFourierProcessor:public FourierProcessor
	{
		public:
			virtual string get_name() const { return "filter.linearfourier"; }

			virtual string get_desc() const
			{
				return "";
			}

			static Processor *NEW()
			{
				return new LinearRampFourierProcessor();
			}

		protected:
			virtual void create_radial_func(vector < float >&radial_mask) const ;
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
			if( params.has_key("highpass") ) {
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
		  virtual void preprocess(EMData * image);
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
	 *@param intercept intercept in 'f(x) = slope * x + intercept
	 *@param slope slope in 'f(x) = slope * x + intercept
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
			d.put("intercept", EMObject::FLOAT, "intercept in 'f(x) = slope * x + intercept'");
			d.put("slope", EMObject::FLOAT, "slope in 'f(x) = slope * x + intercept'");
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

	/**Invert image as if f(x) != 0: f(x) = 1/f(x) else: f(x) = zero_to
	 *@param zero_to  Inverted zero values are set to this value, default is 0
	 */
	class InvertCarefullyProcessor:public RealPixelProcessor
	{
		public:
			string get_name() const
			{
				return "math.invert.carefully";
			}
			static Processor *NEW()
			{
				return new InvertCarefullyProcessor();
			}

			void set_params(const Dict & new_params)
			{
				params = new_params;
				zero_to = params.set_default("zero_to",0.0f);
			}

			TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("zero_to", EMObject::FLOAT, "Inverted zero values are set to this value, default is 0.");
				return d;
			}

			string get_desc() const
			{
				return "if f(x) != 0: f(x) = 1/f(x) else: f(x) = zero_to";
			}

		protected:
			void process_pixel(float *x) const
			{
				if (*x == 0) *x = zero_to;
				else *x = 1/(*x);
			}
		private:
			float zero_to;
	};

	/**Do a math power operation on image, f(x) = x ^ pow;
	 *@param pow Each pixel is raised to this power
	 */
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

	/**Do a square operation on image, f(x) = x * x;
	 */
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

	/**Rotate by 180 using pixel swapping, works for 2D only
	 * @author David Woolford
	 * @date July 29th 2008
	 * @ingroup CUDA_ENABLED
	 */
	class Rotate180Processor:public Processor
	{
		public:
			string get_name() const
			{
				return "math.rotate.180";
			}
			static Processor *NEW()
			{
				return new Rotate180Processor();
			}

			/**
			 * @exception ImageDimensionException if the image dimensions are not 2D
			 */
			void process_inplace(EMData* image);

			string get_desc() const
			{
				return "The 2D image is rotated by 180 degree by carefully swapping image pixel values. No explicit matrix multiplication is performed. If image dimensions are even will change pixels along x=0 and y=0. Works for all combinations of even and oddness.";
			}
	};

	/** Transform the image using a Transform object
	 *@author David Woolford
	 *@date September 2008
	 *@param transform The Transform object that will be applied to the image
	 *@ingroup CUDA_ENABLED
	 */
	class TransformProcessor:public Processor
	{
		public:
			virtual string get_name() const
			{
				return "math.transform";
			}
			static Processor *NEW()
			{
				return new TransformProcessor();
			}

			/**
			 * @exception ImageDimensionException if the image is not 2D or 3D
			 * @exception InvalidParameterException if the Transform parameter is not specified
			 */
			virtual void process_inplace(EMData* image);

			/**
			 * @exception ImageDimensionException if the image is not 2D or 3D
			 * @exception InvalidParameterException if the Transform parameter is not specified
			 */
			virtual EMData* process(const EMData* const image);

			virtual TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("transform", EMObject::TRANSFORM, "The Transform object that will be applied to the image" );
				return d;
			}

			virtual string get_desc() const
			{
				return "The image is transformed using Transform parameter.";
			}
		private:
			float* transform(const EMData* const image, const Transform& t) const;

			void update_emdata_attributes(EMData* const image, const Dict& attr_dict, const float& scale) const;


			void assert_valid_aspect(const EMData* const image) const;
	};

	/** Translate the image an integer amount
	 * Uses EMData::clip_inplace (inplace) and EMData::get_clip (out of place) to do the translation
	 *@ingroup CUDA_ENABLED
	 *@author David Woolford
	 *@date March 2009
	 *@param trans The displacement array, can be length 1-3
	 */
	class IntTranslateProcessor:public Processor
	{
		public:
			virtual string get_name() const
			{
				return "math.translate.int";
			}

			static Processor *NEW()
			{
				return new IntTranslateProcessor();
			}

			/**
			 * @exception ImageDimensionException if the image is not 1,2 or 3D
			 * @exception InvalidParameterException if the Transform parameter is not specified
			 */
			virtual void process_inplace(EMData* image);

			/**
			 * @exception ImageDimensionException if the image is not 1,2 or 3D
			 * @exception InvalidParameterException if the Transform parameter is not specified
			 */
			virtual EMData* process(const EMData* const image);

			virtual TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("trans", EMObject::INTARRAY, "The displacement array, can be length 1-3" );
				return d;
			}

			virtual string get_desc() const
			{
				return "The image is translated an integer amount";
			}
		private:
			/** Check that the particular aspect is valid
		 	* @exception ImageDimensionException if the image is not 1,2 or 3D
			*/
			void assert_valid_aspect(const vector<int>& translation, const EMData* const image) const;

			/** Get the clip region that will achieve the desired translation
			 * @exception ImageDimensionException if the image is not 1,2 or 3D
			 * @param translation the amount by which to translate
			 * @param image the image that will eventually used for the translation operation
			 */
			Region get_clip_region(vector<int>& translation, const EMData* const image) const;
	};

		/** Scale the image with control over the output dimensions
		 *@author David Woolford
		 *@date June 2009
		 *@param scale The amount by which to scale
		 *@param cip The length of each output dimension. Non sophisticated, output dimensions can't be different
		 */
		class ScaleTransformProcessor:public Processor
		{
			public:
				virtual string get_name() const
				{
					return "math.transform.scale";
				}
				static Processor *NEW()
				{
					return new ScaleTransformProcessor();
				}

				/**
				 * @exception ImageDimensionException if the image is not 2D or 3D
				 * @exception InvalidParameterException if the Transform parameter is not specified
				 */
				virtual void process_inplace(EMData* image);

				/**
				 * @exception ImageDimensionException if the image is not 2D or 3D
				 * @exception InvalidParameterException if the Transform parameter is not specified
				 */
				virtual EMData* process(const EMData* const image);

				virtual TypeDict get_param_types() const
				{
					TypeDict d;
					d.put("scale", EMObject::FLOAT, "The amount by which to scale" );
					d.put("clip", EMObject::INT, "The length of each output dimension. Non sophisticated, output dimensions can't be different" );
					return d;
				}

				virtual string get_desc() const
				{
					return "The image is scaled with the clip variable in mind, being sure to preserve as much pixel information as possible.";
				}
		};

	/**f(x) = maxval if f(x) > maxval;
	  * f(x) = minval if f(x) < minval
	  *@param minval the minimum value to clamp to
	  *@param maxval the maximum value to clamp to
	  *@tomean Replace outlying pixels values with the mean pixel value instead
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
			d.put("tomean", EMObject::BOOL, "Replace outlying pixels values with the mean pixel value instead" );
			return d;
		}

		string get_desc() const
		{
			return "This function clamps the min and max vals in the image at minval and maxval, respectively. In a sense this a bi-truncation of the data.";
		}

	  protected:
		float default_max, default_min;
	};

	/**This function clamps the min and max vals in the image at minval and maxval at mean-n*sigma
	 * and mean+n*sigma, respectively. The parameter specified by the user is n, the default value of n is 2.
	 *@param nsigma The number (n) of sigmas to clamp min and max vals at, so that the clamped boundaries are mean-n*sigma and mean+n*sigma
	 *@param tomean Replace outlying pixels values with the mean pixel value instead
	 */
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
				d.put("tomean", EMObject::BOOL, "Replace outlying pixels values with the mean pixel value instead" );
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
	 *@param minval Everything below this value is set to this value
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
	 *@param minval the value that will be set to zero - all values below will also be set to zero. Values above get minval subtracted from them
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
	 *@param value The thresholding value. If a pixel value is equal to or above the threshold it is set to 1. If it is below it is set to 0
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

	/** A thresholding processor for Fourier images based on the amplitude component.
	 * Useful in tomography when you want to count large complex pixels. Not the fastest approach,
	 * if you were going for efficiency it would probably be better just to iterate through the pixels
	 * and count. But if you do it this way you can just get the mean of the resulting image (and multiplying by 2). So it's
	 * basically easier, but lazy.
	 * Originally added for use by e2tomohunter.py
	 *@author David Woolford
	 *@date April 29th 2009
	 *@param value The Fourier amplitude threshold cutoff
	 */
	class BinarizeFourierProcessor:public Processor
		{
		  public:
			virtual string get_name() const
			{
				return "threshold.binary.fourier";
			}
			static Processor *NEW()
			{
				return new BinarizeFourierProcessor();
			}

			/**
			 * @exception ImageFormatException if the input image is not complex
			 * @exception InvalidParameterException if the threshold is less than 0
			 * Note result is always in real-imaginary format
			 * Note input can be real-imaginary or amplitude-phase
			 */
			virtual void process_inplace(EMData* image);

			virtual TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("value", EMObject::FLOAT, "The Fourier amplitude threshold cutoff" );
				return d;
			}

			virtual string get_desc() const
			{
				return "f(k) = 0 + 0i if ||f(k)|| < value; f(k) = 1 + 0i if ||f(k)|| >= value.";
			}
		};

	/**f(x): if v-r<x<v+r -> v; if x>v+r -> x-r; if x<v-r -> x+r
	 *@param range The range about 'value' which will be collapsed to 'value'
	 *@param value The pixel value where the focus of the collapse operation is
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
	 *@param shift The amount to shift pixel values by before scaling
	 *@param scale The scaling factor to be applied to pixel values
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
	 *@param low Pixels are divided by (low - high) prior to the exponential operation
	 *@param high Pixels are divided by (low - high) prior to the exponential operation
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
	
	/**f(x) = f(x) if f(x) is finite | to if f(x) is not finite
	 *@param to Pixels which are not finite will be set to this value
	 */
	class FiniteProcessor:public RealPixelProcessor
	{
		public:
			FiniteProcessor():to(0)
			{
			}

			string get_name() const
			{
				return "math.finite";
			}

			static Processor *NEW()
			{
				return new FiniteProcessor();
			}

			void set_params(const Dict & new_params)
			{
				if (new_params.has_key("to") )
					to = params["to"];
			}

			TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("to", EMObject::FLOAT, "Pixels which are not finite will be set to this value");
				return d;
			}

			string get_desc() const
			{
				return "f(x) = f(x) if f(x) is finite | to if f(x) is not finite";
			}

		protected:
			/**
			* 
			*/
			void process_pixel(float *x) const; 
		private:
			float to;
	};

	/**f(x) = 1 if (low <= x <= high); else f(x) = 0
	 *@param low The lower limit of the range that will be set to 1
	 *@param high The upper limit of the range that will be set to 1
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
	 *@param value1 A number reflecting total standard deviations in the right direction
	 *@param value2 A number reflecting total standard deviations in the left direction
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

	/**f(x) = log10(x) if x > 0; else f(x) = 0
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
			return "f(x) = log10(x) if x > 0; else f(x) = 0;";
		}

	  protected:
		void process_pixel(float *x) const
		{
			if (*x > 0)
			{
				*x = log10(*x);
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
			value = params.set_default("value",0.0f);
		}

		TypeDict get_param_types() const
		{
			TypeDict d = CircularMaskProcessor::get_param_types();
			d.put("value", EMObject::FLOAT, "step cutoff to this value. Default is 0.");
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
			if (dist >= outer_radius_square || dist < inner_radius_square){
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
			return "a gaussian falloff to zero, radius is the 1/e of the width. If inner_radius>0, then \
outer radius specifies width of Gaussian starting at inner_radius rather than total radius.";
		}

	  protected:
		void process_dist_pixel(float *pixel, float dist) const
		{
			if (inner_radius_square>0) {
				if (dist>inner_radius_square)
					(*pixel) *= exp(-pow(sqrt(dist)-inner_radius,2.0f) / outer_radius_square);
			}
			else (*pixel) *= exp(-dist / outer_radius_square);
		}
	};

	/**a gaussian falloff to zero, with nonisotropic widths along x,y,z
	 *@param radius_x x-axis radius
	 *@param radius_y y-axis radius
	 *@param radius_z z-axis radius
	 *@param gauss_width Gaussian falloff width, relative to each radius, default 0.05
	 */
	class MaskGaussNonuniformProcessor:public CoordinateProcessor
	{
	  public:
		MaskGaussNonuniformProcessor():radius_x(0), radius_y(0), radius_z(0), gauss_width(0)
		{
		}

		void set_params(const Dict & new_params)
		{
			params = new_params;

			if (params.has_key("radius_x")) radius_x=params["radius_x"];
			else radius_x=5.0;

			if (params.has_key("radius_y")) radius_y=params["radius_y"];
			else radius_y=5.0;

			if (params.has_key("radius_z")) radius_z=params["radius_z"];
			else radius_z=5.0;

			if (params.has_key("gauss_width")) gauss_width=params["gauss_width"];
			else gauss_width=0.05f;
		}

		TypeDict get_param_types() const
		{
			TypeDict d;

			d.put("radius_x", EMObject::INT, "x-axis radius");
			d.put("radius_y", EMObject::INT, "y-axis radius");
			d.put("radius_z", EMObject::INT, "z-axis radius");
			d.put("gauss_width", EMObject::FLOAT, "Gaussian falloff width, relative to each radius, default 0.05");

			return d;
		}

		string get_name() const
		{
			return "mask.gaussian.nonuniform";
		}
		static Processor *NEW()
		{
			return new MaskGaussNonuniformProcessor();
		}

		string get_desc() const
		{
			return "A Gaussian falloff to zero. Nonisotropic, specify inner radius for x,y,z and Gaussian falloff width. Falloff \
width is also nonisotropic and relative to the radii, with 1 being equal to the radius on that axis.";
		}

	  protected:
		void process_pixel(float *pixel, int xi, int yi, int zi) const
		{
			float dist = pow((xi - nx/2)/radius_x,2.0f) + pow((yi - ny/2)/radius_y,2.0f) + pow((zi - nz/2)/radius_z,2.0f);
			if (dist>1.0) (*pixel)*=exp(-pow((sqrt(dist)-1.0f)/gauss_width,2.0f));
		}

		float radius_x,radius_y,radius_z,gauss_width;
	};

	/**f(x) = f(x) / exp(-radius*radius * gauss_width / (ny*ny))
	 *@param gauss_width Used to calculate the constant factor - gauss_width / (ny*ny)
	 *@param ring_width The width of the mask ring
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
	 *@param areasize The width of the area to process (not radius)
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

	/**BoxStatProcessor files are a kind of neighborhood processors. These processors
	 * compute every output pixel using information from a reduced region on the neighborhood
	 * of the input pixel. The classical form are the 3x3 processors. BoxStatProcessors could
	 * perform diverse tasks ranging from noise reduction, to differential , to mathematical
	 * morphology. BoxStatProcessor class is the base class. Specific BoxStatProcessor needs
	 * to define process_pixel(float *pixel, const float *array, int n).
	 *@param radius The radius of the search box, default is 1 which results in a 3x3 box (3 = 2xradius + 1)
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
	 *@param npeaks the number of surrounding peaks allow to >= pixel values
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
			d.put("npeaks", EMObject::INT, "the number of surrounding peaks allow to >= pixel values");
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

	/**averages over cal_half_width, then sets the value in a local block
	 *@param cal_half_width cal_half_width is dx/dy for calculating an average
	 *@param fill_half_width fill_half_width is dx/dy for fill/step
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
			d.put("cal_half_width", EMObject::FLOAT, "cal_half_width is dx/dy for calculating an average");
			d.put("fill_half_width", EMObject::FLOAT, "fill_half_width is dx/dy for fill/step");
			return d;
		}
	};

	/**Block processor, val1 is dx/dy, val2 is lp freq cutoff in pixels. Mystery processor.
	 *@param value1 val1 is dx/dy
	 *@param value2 val2 is lowpass freq cutoff in pixels
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
			d.put("value1", EMObject::FLOAT, "val1 is dx/dy");
			d.put("value2", EMObject::FLOAT, "val2 is lowpass freq cutoff in pixels");
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
	class BooleanShrinkProcessor
	{
		protected:
			/** Boolean shrink an image, returning the processed image
			* @param image the image to operate on
			* @param params parameter dictionary
			* @exception ImageFormatException if the image is complex
			* @exception NullPointerException if the image pointer is null
			* @return the image that results from the operation
			 */
			template<class LogicOp>
			EMData* process(const EMData *const image, Dict& params);

			/** Boolean shrink an image inplace
			 * @param image the image to operate on
			 * @param params parameter dictionary
			 * @exception ImageFormatException if the image is complex
			 * @exception NullPointerException if the image pointer is null
			 */
			template<class LogicOp>
			void process_inplace(EMData * image, Dict& params);

	};

	/** MaxShrinkProcessors shrinks an image by in an integer amount,
	 * keeping the maximum pixel value - useful when constructing binary search
	 * trees in the marching cubes algorithm
	 *@author David Woolford
	 *@date September 2007
	 *@param n The shrink factor
	 *@param search The search area (cubic volume width, usually the same as shrink)
	 */
	class MaxShrinkProcessor:public BooleanShrinkProcessor, public Processor
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
				return BooleanShrinkProcessor::process<GreaterThan>(image, params);
			}

			// resizes the image
			virtual void process_inplace(EMData * image)
			{
				BooleanShrinkProcessor::process_inplace<GreaterThan>(image, params);
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
	 *@author David Woolford
	 *@date September 2007
	 *@param n The shrink factor
	 *@param search The search area (cubic volume width, usually the same as shrink)
	 */
	class MinShrinkProcessor:public BooleanShrinkProcessor, public Processor
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
				return BooleanShrinkProcessor::process<LessThan>(image, params);
			}

			// resizes the image
			virtual void process_inplace(EMData * image)
			{
				BooleanShrinkProcessor::process_inplace<LessThan>(image, params);
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
	 *@author David Woolford (But is basically a copy of the old EMData::mean_shrink, probably written by Steven Ludtke )
	 *@date May 2008
	 *@param n The shrink factor
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
				return "Shrink an image by a given amount , using the mean value found in the pixel neighborhood.";
			}

			virtual string get_name() const
			{
				return "math.meanshrink";
			}
			static Processor *NEW()
			{
				return new MeanShrinkProcessor();
			}

			virtual TypeDict get_param_types() const
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
	*@author David Woolford (But is basically a copy of the old EMData::median_shrink, probably written by Steven Ludtke )
	*@date May 2008
	*@param n The shrink factor
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
			return "Shrink an image by a given amount , using the median value found in the pixel neighborhood.";
		}

		virtual string get_name() const
		{
			return "math.medianshrink";
		}
		static Processor *NEW()
		{
			return new MedianShrinkProcessor();
		}

		virtual TypeDict get_param_types() const
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
		* @param shrink_factor the shrink amount
			*/
		void accrue_median(EMData* to, const EMData* const from,const int shrink_factor);
	};


	/** FFTResampleProcessor resamples an image by clipping the Fourier Transform
	 * This is the same as multipyling the FT by a box window, in real space this is a convolution that will induce rippling.
	 * Hence it should probably be combined with a damping function near the edge
	 * Works for even/odd, 1D, 2D and 3D images (completely robust, tested)
	 *@author David Woolford
	 *@date June 2009
	 *@param n sampe_rate
	 */
	class FFTResampleProcessor : public Processor
	{
		public:
			virtual EMData* process(const EMData *const image);

			virtual void process_inplace(EMData * image);

			string get_desc() const
			{
				return "Robust resampling of an image by clipping its Fourier transform.";
			}

			string get_name() const
			{
				return "math.fft.resample";
			}
			static Processor *NEW()
			{
				return new FFTResampleProcessor();
			}

			TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("n", EMObject::FLOAT, "The sample rate. Less than one enlarges the image, greater than one shrinks it.");
				return d;
			}
		private:
			/** An internal function that encapsulates a routine common to both process and
			* process inplace
			* @param to the smaller image that will store the shrunken values
			* @param from the larger image that will be used to calculate the shrunken values
			* @param shrinkfactor the shrink amount
			 */
			void fft_resample(EMData* to, const EMData *const from, const float& sample_rate);

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
			d.put("value1", EMObject::FLOAT, "sig multiplier");
			d.put("value2", EMObject::FLOAT, "x of center");
			d.put("value3", EMObject::FLOAT, "y of center");
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

	/**Decay edges of image to zero
	 *@param width Width of the decay region around the edge of the image in pixels
	 */
	class DecayEdgeProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);
		string get_name() const
		{
			return "mask.decayedge2d";
		}

		static Processor *NEW()
		{
			return new DecayEdgeProcessor();
		}

		string get_desc() const
		{
			return "Decay edges of image to zero";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("width", EMObject::INT, "Width of the decay region around the edge of the image in pixels");
			return d;
		}
	};

	/**zero edges of image on top and bottom, and on left and right.
	 *@param x0 The number of columns to zero from left
	 *@param x1 The number of columns to zero from right
	 *@param y0 The number of rows to zero from the bottom
	 *@param y1 The number of rows to zero from the top
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
	 *@param x0 The number of columns to zero from left
	 *@param x1 The number of columns to zero from right
	 *@param y0 The number of rows to zero from the bottom
	 *@param y1 The number of rows to zero from the top
	 *@param z0 The number of slices to zero from the bottom
	 *@param z1 The number of slices to zero from the top
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
	 *@param mask the 1/0 mask defining a region to use for the zero-normalization
	 *@param no_sigma if this flag is zero, only average under the mask will be substracted. set this flag to 1, standard deviation not modified
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

	/**Normalize the image whilst also removing any ramps. Ramps are removed first, then mean and sigma becomes 0 and 1 respectively
	* This is essential Pawel Penczek's preferred method of particle normalization
	* @author David Woolford
	* @date Mid 2008
	*/
	class NormalizeRampNormVar: public Processor
	{
		public:
			string get_name() const
			{
				return "normalize.ramp.normvar";
			}

			static Processor *NEW()
			{
				return new NormalizeRampNormVar();
			}

			string get_desc() const
			{
				return "First call filter.ramp on the image, then make the mean 0 and norm 1";
			}

			void process_inplace(EMData * image);
	};

	/** Normalize the mass of the image assuming a density of 1.35 g/ml (0.81 Da/A^3).
	 * Only works for 3D images. Essentially a replica of Volume.C in EMAN1.
	 *@author David Woolford (a direct port of Steve Ludtke's code)
	 *@date 01/17/09
	 *@param apix Angstrom per pixel of the image. If not set will use the apix_x attribute of the image
	 *@param mass The approximate mass of protein/structure in kilodaltons
	 *@param thr The isosurface threshold which encapsulates the structure
	 */
	class NormalizeByMassProcessor: public Processor
	{
		public:
			string get_name() const
			{
				return "normalize.bymass";
			}

			static Processor *NEW()
			{
				return new NormalizeByMassProcessor();
			}

			string get_desc() const
			{
				return "Normalize the mass of the image assuming a density of 1.35 g/ml (0.81 Da/A^3) (3D only)";
			}

			TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("apix", EMObject::FLOAT,"Angstrom per pixel of the image. If not set will use the apix_x attribute of the image");
				d.put("mass", EMObject::FLOAT,"The approximate mass of protein/structure in kilodaltons");
				d.put("thr", EMObject::FLOAT,"The isosurface threshold which encapsulates the structure");
				return d;
			}

			void process_inplace(EMData * image);
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
	 *@param noisy Image to normalize to
	 *@param keepzero set to 1 to ignore zeroes
	 *@param invert set to one to invert image
	 *@param sigmax Values greater or less than sigma from zero will be excluded from the normalization
	 *@param mult multiply by this factor
	 *@param add add by this factor
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
			d.put("keepzero", EMObject::INT,"set to 1 to ignore zeroes");
			d.put("invert", EMObject::INT,"set to one to invert image");
			d.put("sigmax",EMObject::FLOAT,"Values greater or less than sigma from zero will be excluded from the normalization");
			d.put("mult", EMObject::FLOAT, "multiply by this factor");
			d.put("add", EMObject::FLOAT, "add by this factor");
			return d;
		}
	};

	/**Multiply this image by a constant so it is scaled to the signal in 'noisyfile'
	 *@param noisyfile Image to normalize to
	 *@param keepzero  set to 1 to exclude zero values
	 *@param invert set to 1 to invert image
	 *@param mult multiply by this factor
	 *@param add add by this factor
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
			d.put("noisyfile", EMObject::STRING, "Image file name to normalize to");
			d.put("keepzero", EMObject::INT, "set to 1 to ignore zeroes");
			d.put("invert", EMObject::INT, "set to 1 to invert image");
			d.put("mult", EMObject::FLOAT, "multiply by this factor");
			d.put("add", EMObject::FLOAT, "add by this factor");
			return d;
		}

		string get_desc() const
		{
			return "Multiply this image by a constant so it is scaled to the signal in 'noisyfile'";
		}

	};

	/**use least square method to normalize
	 * @param to reference image normalize to
	 * @param low_threshold only take into account the reference image's pixel value between high and low threshold (zero is ignored)
	 * @param high_threshold only take into account the reference image's pixel value between high and low threshold (zero is ignored)
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
			d.put("to", EMObject::EMDATA, "reference image normalize to");
			d.put("low_threshold", EMObject::FLOAT, "only take into account the reference image's pixel value between high and low threshold (zero is ignored)");
			d.put("high_threshold", EMObject::FLOAT, "only take into account the reference image's pixel value between high and low threshold (zero is ignored)");
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
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "math.radialsubtract";
		}

		static Processor *NEW()
		{
			return new RadialSubstractProcessor();
		}

		virtual string get_desc() const
		{
			return "subtracts circularly symmetric part of an image.";
		}

	};

	/** Transpose a 2D image
	 * @author David Woolford
	 * @date April 27th 2009
	 * @ingroup tested3c
	 */
	class TransposeProcessor:public Processor
	{
	  public:

		/** See Processor comments for more details
		 * @exception UnexpectedBehaviorException if the image is not 2D
		 * @exception UnexpectedBehaviorException if the image is complex
		 */
		virtual void process_inplace(EMData * image);

		/** See Processor comments for more details
		 * @exception UnexpectedBehaviorException if the image is not 2D
		 * @exception UnexpectedBehaviorException if the image is complex
		 */
		virtual EMData* process(const EMData * const image);

		virtual string get_name() const
		{
			return "xform.transpose";
		}

		static Processor *NEW()
		{
			return new TransposeProcessor();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			return d;
		}

		virtual string get_desc() const
		{
			return "Get the transpose of an image. Works for 2D only";
		}

	};


	/** flip an image around an axis
	 * @param axis  'x', 'y', or 'z' axis. 'x' means horizonal flip; 'y' means vertical flip;
	 */
	class FlipProcessor:public Processor
	{
	  public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "xform.flip";
		}

		static Processor *NEW()
		{
			return new FlipProcessor();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("axis", EMObject::STRING, "'x', 'y', or 'z' axis. 'x' means horizonal flip; 'y' means vertical flip;");
			return d;
		}

		virtual string get_desc() const
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
	 * @param noise noise factor used to generate Gaussian distribution random noise
	 * @param seed seed for random number generator
	 */
	class AddNoiseProcessor:public Processor
	{
	  public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "math.addnoise";
		}

		static Processor *NEW()
		{
			return new AddNoiseProcessor();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("noise", EMObject::FLOAT, "noise factor used to generate Gaussian distribution random noise");
			d.put("seed", EMObject::INT, "seed for random number generator");
			return d;
		}

		virtual string get_desc() const
		{
			return "add noise to an image, image multiply by noise then add a random value";
		}

	  protected:
		virtual float get_sigma(EMData *)
		{
			return 1.0;
		}
	};

	/** add sigma noise, multiply image's sigma value to noise
	 */
	class AddSigmaNoiseProcessor:public AddNoiseProcessor
	{
	  public:
		virtual string get_name() const
		{
			return "math.addsignoise";
		}

		static Processor *NEW()
		{
			return new AddSigmaNoiseProcessor();
		}

		virtual string get_desc() const
		{
			return "add sigma noise.";
		}

	  protected:
		float get_sigma(EMData * image);
	};

	/**add spectral noise to a complex image
	 *@param n
	 *@param x0
	 *@param dx
	 *@param y
	 *@param interpolation
	 *@param seed seed for random number generator
	 */
	class AddRandomNoiseProcessor:public Processor
	{
	  public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "addspectralnoise";
		}

		static Processor *NEW()
		{
			return new AddRandomNoiseProcessor();
		}

		virtual TypeDict get_param_types() const
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

		virtual string get_desc() const
		{
			return "add spectral noise to a complex image.";
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
			virtual void process_inplace(EMData * image);

			virtual string get_name() const
			{
				return "xform.fourierorigin.tocorner";
			}

			static Processor *NEW()
			{
				return new FourierToCornerProcessor();
			}

			virtual string get_desc() const
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
			virtual void process_inplace(EMData * image);

			virtual string get_name() const
			{
				return "xform.fourierorigin.tocenter";
			}

			static Processor *NEW()
			{
				return new FourierToCenterProcessor();
			}

			virtual string get_desc() const
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
			virtual void process_inplace(EMData * image);

			virtual string get_name() const
			{
				return "xform.phaseorigin.tocenter";
			}

			static Processor *NEW()
			{
				return new PhaseToCenterProcessor();
			}

			virtual string get_desc() const
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
			virtual void process_inplace(EMData * image);

			virtual string get_name() const
			{
				return "xform.phaseorigin.tocorner";
			}

			static Processor *NEW()
			{
				return new PhaseToCornerProcessor();
			}

			virtual string get_desc() const
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
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "mask.auto2d";
		}

		static Processor *NEW()
		{
			return new AutoMask2DProcessor();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("threshold", EMObject::FLOAT, "runs from ~ -2 to 2, negative numbers for dark protein and positive numbers for light protein (stain).");
			d.put("filter", EMObject::FLOAT, "is expressed as a fraction of the fourier radius.");
			return d;
		}

		virtual string get_desc() const
		{
			return "Attempts to automatically mask out the particle, excluding other particles in the box, etc.";
		}

	};


	/**Tries to mask out only interesting density
	 * @param au  The asymmetric unit to mask, if -1 masks all asymmetric units but assigns each a unique value
	 * @param sym The symmetry type, for example, "tet"
	 * @author David Woolford
	 * @date February 2009
	*/
	class AutoMaskAsymUnit:public Processor
	{
		public:
			virtual void process_inplace(EMData * image);

			virtual string get_name() const
			{
				return "mask.asymunit";
			}

			static Processor *NEW()
			{
				return new AutoMaskAsymUnit();
			}

			virtual TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("au", EMObject::INT, "The asymmetric unit to mask out. If this is -1 will mask all asymmetric units, giving each a unique number.");
				d.put("sym", EMObject::STRING, "The symmetry, for example, d7");
				return d;
			}

			virtual string get_desc() const
			{
				return "Masks out a specific asymmetric unit of the given symmetry. If the au parameter is -1 will mask all asymmetric units, assigning the asymetric unit number to the masked area.";
			}

	};

	/**Tries to mask out only interesting density
	 * @param threshold1
	 * @param threshold2
	 */
	class AutoMask3DProcessor:public Processor
	{
	  public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "mask.auto3d.thresh";
		}

		static Processor *NEW()
		{
			return new AutoMask3DProcessor();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("threshold1", EMObject::FLOAT);
			d.put("threshold2", EMObject::FLOAT);
			return d;
		}

		virtual string get_desc() const
		{
			return "Tries to mask out only interesting density";
		}

		static void search_nearby(float *dat, float *dat2, int nx, int ny, int nz, float thr);
		static void fill_nearby(float *dat2, int nx, int ny, int nz);
	};

	/** Tries to mask out only interesting density
	 * @param radius Pixel radius of a ball which is used to seed the flood filling operation
	 * @param threshold An isosurface threshold that suitably encases the mass
	 * @param nshells The number of dilation operations
	 * @param nshellsgauss number of Gaussian pixels to expand, following the dilation operations
	 * @return_mask If true the result of the operation will produce the mask, not the masked volume
	 */
	class AutoMask3D2Processor:public Processor
	{
	  public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "mask.auto3d";
		}

		static Processor *NEW()
		{
			return new AutoMask3D2Processor();
		}

		virtual string get_desc() const
		{
			return "Tries to mask out only interesting density using something akin to a flood file approach.";
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("radius", EMObject::INT,"Pixel radius of a ball which is used to seed the flood filling operation.");
			d.put("threshold", EMObject::FLOAT, "An isosurface threshold that suitably encases the mass.");
			d.put("nshells", EMObject::INT, "The number of dilation operations");
			d.put("nshellsgauss", EMObject::INT, "number of Gaussian pixels to expand, following the dilation operations");
			d.put("return_mask", EMObject::BOOL, "If true the result of the operation will produce the mask, not the masked volume.");
			return d;

		}
	};

	/**Add additional shells/rings to an existing 1/0 mask image
	 * @param nshells   number of shells to add
	*/
	class AddMaskShellProcessor:public Processor
	{
	  public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "mask.addshells";
		}

		virtual string get_desc() const
		{
			return "Add additional shells/rings to an existing 1/0 mask image";
		}

		static Processor *NEW()
		{
			return new AddMaskShellProcessor();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("nshells", EMObject::INT, "number of shells to add");
			return d;
		}
	};

	/**ToMassCenterProcessor centers image at center of mass, ignores old dx, dy.
	 * @param int_shift_only set to 1 only shift by integer, no interpolation
	 * @ingroup tested3c
	 */
	class PhaseToMassCenterProcessor:public Processor
	{
		public:
			virtual void process_inplace(EMData * image);

			virtual string get_name() const
			{
				return "xform.phasecenterofmass";
			}

			static Processor *NEW()
			{
				return new PhaseToMassCenterProcessor();
			}

			virtual string get_desc() const
			{
				return "centers the image the center of mass, which is calculated using Fourier phases, ignores old dx, dy.";
			}

			virtual TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("int_shift_only", EMObject::INT, "set to 1 only shift by integer, no interpolation");
				return d;
			}
	};

	/**ToMassCenterProcessor centers image at center of mass, ignores old dx, dy.
	 * @param int_shift_only set to 1 only shift by integer, no interpolation
	 * @ingroup tested3c
	 */
	class ToMassCenterProcessor:public Processor
	{
	  public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "xform.centerofmass";
		}

		static Processor *NEW()
		{
			return new ToMassCenterProcessor();
		}

		virtual string get_desc() const
		{
			return "ToMassCenterProcessor centers image at center of mass, ignores old dx, dy.";
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("int_shift_only", EMObject::INT, "set to 1 only shift by integer, no interpolation");
			return d;
		}
	};

	/**Center image using auto convolution with 180 degree rotation.
	 * @ingroup tested3c
	 */
	class ACFCenterProcessor:public Processor
	{
	  public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "xform.centeracf";
		}

		static Processor *NEW()
		{
			return new ACFCenterProcessor();
		}

		virtual string get_desc() const
		{
			return "Center image using self-convolution.";
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			return d;
		}
	};

	/** Processor the images by the estimated SNR in each image.if parameter 'wiener' is 1, then wiener processor the images using the estimated SNR with CTF amplitude correction.
	 * @param wiener if set to 1,  then use wiener processor to process the images using the estimated SNR with CTF amplitude correction
	 * @param snrfile structure factor file name
	 */
	class SNRProcessor:public Processor
	{
	  public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "eman1.filter.snr";
		}

		static Processor *NEW()
		{
			return new SNRProcessor();
		}

		virtual string get_desc() const
		{
			return "Processor the images by the estimated SNR in each image.if parameter 'wiener' is 1, then wiener processor the images using the estimated SNR with CTF amplitude correction.";
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("wiener", EMObject::INT, "if set to 1,  then use wiener processor to process the images using the estimated SNR with CTF amplitude correction");
			d.put("snrfile", EMObject::STRING, "structure factor file name");
			return d;
		}
	};

	/** A fourier processor specified in a 2 column text file.
	 * @param filename file name for a 2 column text file which specified a radial function data array
	*/
	class FileFourierProcessor:public Processor
	{
	  public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "eman1.filter.byfile";
		}

		virtual string get_desc() const
		{
			return "A fourier processor specified in a 2 column text file.";
		}

		static Processor *NEW()
		{
			return new FileFourierProcessor();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("filename", EMObject::STRING, "file name for a 2 column text file which specified a radial function data array.");
			return d;
		}
	};

	/** Identifiy the best symmetry in the given symmetry list for each pixel and then apply the best symmetry to each pixel
	 *
	 *@author Wen Jiang <wjiang@bcm.tmc.edu>
	 *@date 2005-1-9
	 *@param sym[in] the list of symmetries to search
	 *@param thresh[in] the minimal level of symmetry to be accepted (0-1)
	 *@param output_symlabel[in] if output the symmetry label map in which the pixel value is the index of symmetry in the symmetry list
	 *@param symlabel_map[out] the optional return map when output_symlabel=1
	 */

	class SymSearchProcessor:public Processor
	{
	  public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "misc.symsearch";
		}

		virtual string get_desc() const
		{
			return "Identifiy the best symmetry in the given symmetry list for each pixel and then apply the best symmetry to each pixel.";
		}

		static Processor *NEW()
		{
			return new SymSearchProcessor();
		}

		virtual TypeDict get_param_types() const
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
	 *@param threshold an isosurface threshold at which all desired features are visible
	 *@param radius a normalization size similar to an lp= value
	 *@param apix Angstrom per pixel ratio
	 */
	class LocalNormProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "misc.localnorm";
		}

		static Processor *NEW()
		{
			return new LocalNormProcessor();
		}

		virtual string get_desc() const
		{
			return "This processor attempts to perform a 'local normalization' so low density and high density features will be on a more even playing field in an isosurface display. threshold is an isosurface threshold at which all desired features are visible, radius is a normalization size similar to an lp= value.";
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("threshold", EMObject::FLOAT, "an isosurface threshold at which all desired features are visible");
			d.put("radius", EMObject::FLOAT, "a normalization size similar to an lp= value");
			d.put("apix", EMObject::FLOAT, "Angstrom per pixel ratio");
			return d;
		}
	};

	/**Multiplies the image by the specified file using pixel indices. The images must be same size. If 'ismaskset=' is 1, it will take a file containing a set of masks and apply the first mask to the image.
	 *@param filename mask image file name
	 *@param ismaskset If set to 1, it will take a file containing a set of masks and apply the first mask to the image
	 */
	class IndexMaskFileProcessor:public Processor
	{
	  public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "mask.fromfile";
		}

		static Processor *NEW()
		{
			return new IndexMaskFileProcessor();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("filename", EMObject::STRING, "mask image file name");
			d.put("ismaskset", EMObject::INT, "If set to 1, it will take a file containing a set of masks and apply the first mask to the image");
			return d;
		}

		virtual string get_desc() const
		{
			return "Multiplies the image by the specified file using pixel indices. The images must be same size. If 'ismaskset=' is 1, it will take a file containing a set of masks and apply the first mask to the image.";
		}

	};

	/**Multiplies the image by the specified file using pixel coordinates instead of pixel indices. The images can be different size.
	 *@param filename mask image file name
	 */
	class CoordinateMaskFileProcessor:public Processor
	{
	  public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "mask.fromfile.sizediff";
		}

		static Processor *NEW()
		{
			return new CoordinateMaskFileProcessor();
		}

		virtual string get_desc() const
		{
			return "Multiplies the image by the specified file using pixel coordinates instead of pixel indices. The images can be different size.";
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("filename", EMObject::STRING, "mask image file name");
			return d;
		}
	};

	/**'paints' a circle into the image at x,y,z with values inside r1 set to v1, values between r1 and r2 will be set to a
	 * value between v1 and v2, and values outside r2 will be unchanged
	 *@param x x coordinate for Center of circle
	 *@param y y coordinate for Center of circle
	 *@param z z coordinate for Center of circle
	 *@param r1 Inner radius
	 *@param v1 Inner value
	 *@param r2 Outter radius
	 *@param v2 Outer Value
	 */
	class PaintProcessor:public Processor
	{
	  public:
		PaintProcessor():x(0), y(0), z(0),r1(0), v1(0.0), r2(0), v2(0.0)
		{
		}

		virtual string get_name() const
		{
			return "mask.paint";
		}

		static Processor *NEW()
		{
			return new PaintProcessor();
		}

		virtual string get_desc() const
		{
			return "Paints a circle with a decaying edge into the image. r<r1 -> v1, r1<r<r2 -> (v1,v2), r>r2 unchanged";
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("x", EMObject::INT, "x coordinate for Center of circle");
			d.put("y", EMObject::INT, "y coordinate for Center of circle");
			d.put("z", EMObject::INT, "z coordinate for Center of circle");
			d.put("r1", EMObject::INT, "Inner radius");
			d.put("v1", EMObject::FLOAT, "Inner value");
			d.put("r2", EMObject::INT, "Outter radius");
			d.put("v2", EMObject::FLOAT, "Outer Value");
			return d;
		}

		virtual void set_params(const Dict & new_params)
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
		virtual void process_inplace(EMData *image);

		int x,y,z,r1;
		float v1;
		int r2;
		float v2;

	};


	/**Does a projection in one the axial directions
	 * Doesn't support process_inplace (because the output has potentially been compressed in one dimension)
	 * @param direction The direction of the sum, either x,y or z
	 */
	class DirectionalSumProcessor : public Processor
	{
	  public:
		virtual string get_name() const
		{
			return "misc.directional_sum";
		}

		static Processor *NEW()
		{
			return new DirectionalSumProcessor();
		}

		/**
		 * @exception InvalidParameterException raised if the direction parameter is not "x", "y" or "z"
		 */
		virtual EMData* process(const EMData* const image);

		/**
		 * @exception InvalidCallException raised if this function is called
		 */
		virtual void process_inplace(EMData* image ) {
			throw InvalidCallException("The directional sum processor does not work inplace");
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("direction", EMObject::STRING,"The direction of the sum, either x,y or z");
			return d;
		}

		string get_desc() const
		{
			return "Calculates the projection of the image along one of the axial directions, either x, y or z";
		}

	};

	/**'paints' a circle into the image at x,y,z with values inside r1 set to v1, values between r1 and r2 will be set to a
	 *	value between v1 and v2, and values outside r2 will be unchanged
	 *	@param xpoints x coordinates
	 *	@param ypoints y coordinates
	 *	@param zpoints z coordinates
	 *	@param minval min value
	 */
	class WatershedProcessor:public Processor
	{
	  public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "watershed";
		}

		static Processor *NEW()
		{
			return new WatershedProcessor();
		}

		virtual string get_desc() const
		{
			return "Does a watershed";
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("xpoints", EMObject::FLOATARRAY,"x coordinates");
			d.put("ypoints", EMObject::FLOATARRAY,"y coordinates");
			d.put("zpoints", EMObject::FLOATARRAY,"z coordinates");
			d.put("minval", EMObject::FLOAT,"min value");
			return d;
		}
	  private:
		  vector<Vec3i > watershed(EMData* mask, EMData* image, const float& threshold, const Vec3i& cordinate, const int mask_value);
		  vector<Vec3i > find_region(EMData* mask,const vector<Vec3i >& coords, const int mask_value, vector<Vec3i >& region);

	};

	/**Sets the structure factor based on a 1D x/y text file.
	 *@param filename file name for structure factor
	 */
	class SetSFProcessor:public Processor
	{
	  public:
		  virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "misc.setpowspec";
		}

		virtual string get_desc() const
		{
			return "Sets the structure factor based on a 1D x/y text file.";
		}

		static Processor *NEW()
		{
			return new SetSFProcessor();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("filename", EMObject::STRING, "file name for structure factor");
			return d;
		}
	};

	/**Smart mask processor
	 *@param mask mask value
	 */
	class SmartMaskProcessor:public Processor
	{
	  public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "mask.smart";
		}

		static Processor *NEW()
		{
			return new SmartMaskProcessor();
		}

		virtual string get_desc() const
		{
			return "Smart mask processor.";
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("mask", EMObject::FLOAT, "mask value");
			return d;
		}
	};

	/**Iterative expansion of a binary mask, val1 is number of pixels to expand, if val2!=0 will make a soft Gaussian edge starting after val2 pixels.
	 * @param val1 number of pixels to expand
	 * @param val2 number of Gaussian pixels to expand, following the first expansion
	 */
	class IterBinMaskProcessor:public Processor
	{
	  public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "mask.addshells.gauss";
		}

		virtual string get_desc() const
		{
			return "Iterative expansion of a binary mask, val1 is number of pixels to expand, if val2!=0 will make a soft Gaussian edge starting after val2 pixels.";
		}

		static Processor *NEW()
		{
			return new IterBinMaskProcessor();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("val1", EMObject::FLOAT, "number of pixels to expand");
			d.put("val2", EMObject::FLOAT, "number of Gaussian pixels to expand, following the first expansion");
			return d;
		}
	};

	/**Base class for a group of 'processor' used to create test image.
	 */
	class TestImageProcessor : public Processor
	{
	public:
		static string get_group_desc()
		{
			return "Base class for a group of 'processors' used to create test image.";
		}

	protected:
		void preprocess(EMData * image);
		int nx, ny, nz; //this is the size of the source image
	};

	/**Replace a source image as a strict Gaussian
	 *@param x_sigma sigma value for this Gaussian blob on x direction
	 *@param y_sigma sigma value for this Gaussian blob on y direction
	 *@param z_sigma sigma value for this Gaussian blob on z direction
	 *@param x_center center for this Gaussian blob on x direction
	 *@param y_center center for this Gaussian blob on y direction
	 *@param z_center center for this Gaussian blob on z direction
	 */
	class TestImagePureGaussian : public TestImageProcessor
	{
	public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "testimage.puregaussian";
		}

		virtual string get_desc() const
		{
			return "Replace a source image as a strict Gaussian ";
		}

		static Processor * NEW()
		{
			return new TestImagePureGaussian();
		}

		virtual TypeDict get_param_types() const
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

	/**Replace a source image as a strict Gaussian
	 *@param sigma sigma value for this Gaussian blob
	 */
	class TestImageFourierNoiseGaussian : public TestImageProcessor
	{
	public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "testimage.noise.fourier.gaussian";
		}

		virtual string get_desc() const
		{
			return "Replace a source image with pink Fourier noise, based on a Gaussian. Random phase.";
		}

		static Processor * NEW()
		{
			return new TestImageFourierNoiseGaussian();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("sigma", EMObject::FLOAT, "sigma value");
			return d;
		}
	};
	
	/**
	 * @author David Woolford
	 * @date June 15th 2009
	 */
	class TestImageFourierNoiseProfile : public TestImageProcessor
	{
	public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "testimage.noise.fourier.profile";
		}

		virtual string get_desc() const
		{
			return "Replace a source image with Fourier noise using amplitude information that is stored in a profile.";
		}

		static Processor * NEW()
		{
			return new TestImageFourierNoiseProfile();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("profile", EMObject::FLOATARRAY, "The noise profile, squared amplitude. As in, what is the EMAN2CTF.background attribute");
			return d;
		}
	};
	
	
	/**
	* @author David Woolford
	* @date June 15th 2009
	*/
	class CTFSNRWeightProcessor : public TestImageProcessor
	{
		public:
			virtual void process_inplace(EMData * image);

			virtual string get_name() const
			{
				return "ctf.snr.weight";
			}

			virtual string get_desc() const
			{
				return "Weight the amplitudes of an image based on radial noise and snr curves ";
			}

			static Processor * NEW()
			{
				return new CTFSNRWeightProcessor();
			}

			virtual TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("noise", EMObject::FLOATARRAY, "The noise profile, squared amplitude. As in, what is the EMAN2CTF.background attribute");
				d.put("snr", EMObject::FLOATARRAY, "Squared amplitude divided by squared noise amplitude. As in, what is the EMAN2CTF.snr attribute");
				d.put("boost", EMObject::FLOAT, "Multiplicative signal boost");
				return d;
			}
	};



	/** Treats the pixels as though they are 1D (even if the image is 2D or 3D),
	 * inserting a sine wave of pixel period extracted from the parameters (default is 10)
	 *@param period the period of the sine wave
	 */
	class TestImageLineWave : public TestImageProcessor
	{
		public:
			virtual void process_inplace(EMData * image);

			virtual string get_name() const
			{
				return "testimage.linewave";
			}

			virtual string get_desc() const
			{
				return "Insert an oscillating sine wave into the pixel data";
			}

			static Processor * NEW()
			{
				return new TestImageLineWave();
			}

			virtual TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("period", EMObject::FLOAT, "The period of the oscillating sine wave. Default 10.");
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
			virtual void process_inplace(EMData * image);

			virtual string get_name() const
			{
				return "testimage.tomo.objects";
			}

			virtual string get_desc() const
			{
				return "Make an image consisting various objects, useful for tomographic testing";
			}

			static Processor * NEW()
			{
				return new TestTomoImage();
			}
		private:
			void insert_solid_ellipse( EMData* image, const Region& region, const float& value, const Transform& t3d  = Transform() );
			void insert_hollow_ellipse( EMData* image, const Region& region, const float& value, const int& radius, const Transform& t3d = Transform() );
			void insert_rectangle( EMData* image, const Region& region, const float& value, const Transform& t3d = Transform() );
	};

	/** Put a gradient in the image of the form y = mx+b : "x" is a string indicating any of the image axes, i.e., x,y or z.
	 *@author David Woolford <woolford@bcm.edu>
	 *@date 01/10/2008
	 *@param axis The axis the will be used to determine pixel values. Must be x,y or z. Default is x
	 *@param m m in the equation m*axis+b. Default is 1.0
	 *@param b b in the equation m*axis+b. Default is 0.0
	*/
	class TestImageGradient : public TestImageProcessor
	{
		public:
			virtual void process_inplace(EMData * image);

			virtual string get_name() const
			{
				return "testimage.gradient";
			}

			virtual string get_desc() const
			{
				return "Make a gradient image of the form y=mx+b, where x is any of the image axes.";
			}

			static Processor * NEW()
			{
				return new TestImageGradient();
			}

			virtual TypeDict get_param_types() const
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
	 *@author David Woolford <woolford@bcm.edu>
	 *@date October 2007
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
			virtual void process_inplace(EMData * image);

			virtual string get_name() const
			{
				return "testimage.axes";
			}

			virtual string get_desc() const
			{
				return "Make an image consisting of a single cross";
			}

			static Processor * NEW()
			{
				return new TestImageAxes();
			}

			virtual TypeDict get_param_types() const
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
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "testimage.gaussian";
		}

		virtual string get_desc() const
		{
			return "Replace a source image as a Gaussian Blob";
		}

		static Processor * NEW()
		{
			return new TestImageGaussian();
		}

		virtual TypeDict get_param_types() const
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
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "testimage.scurve";
		}

		virtual string get_desc() const
		{
			return "Replace a source image with a lumpy S-curve used for alignment testing";
		}

		static Processor * NEW()
		{
			return new TestImageScurve();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			return d;
		}
	};

	/**Replace a source image as a sine wave in specified wave length
	 *@param wavelength cos(2*pi*r/wavelength+phase)
	 *@param phase in radians
	 *@param x center of the spherical wave
	 *@param y center of the spherical wave
	 *@param z center of the spherical wave
	 */
	class TestImageSphericalWave : public TestImageProcessor
	{
	public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "testimage.sphericalwave";
		}

		virtual string get_desc() const
		{
			return "Replace a source image in 2d or 3d with a spherical wave cos(2*pi*r/wavelength+phase) also 1/r (2d) or 1/r^2 (3d)";
		}

		static Processor * NEW()
		{
			return new TestImageSphericalWave();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("wavelength", EMObject::FLOAT, "cos(2*pi*r/wavelength+phase)");
			d.put("phase", EMObject::FLOAT, "in radians");
			d.put("x", EMObject::FLOAT, "center of the spherical wave");
			d.put("y", EMObject::FLOAT, "center of the spherical wave");
			d.put("z", EMObject::FLOAT, "center of the spherical wave");
			return d;
		}
	};


	/**Replace a source image as a sine wave in specified wave length
	 *@param wavelength wavelength in equation sin(x*2*PI/wavelength - phase*180/PI)
	 *@param axis (optional) specify a major axis for asymmetric features, default x axis
	 *@param phase (optional) the phase in equation sin(x*2*PI/wavelength - phase*180/PI)
	 *@param az (optional) angle in degree. for 2D image, this is the rotated angle of the image, in 3D image, it's az for euler angle. default is zero
	 *@param alt (optional) angle in degree. only in 3D case, alt for euler angle, default is zero
	 *@param phi (optional) angle in degree. only in 3D case, phi for euler angle, default is zero
	 */
	class TestImageSinewave : public TestImageProcessor
	{
	public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "testimage.sinewave";
		}

		virtual string get_desc() const
		{
			return "Replace a source image as a sine wave in specified wave length";
		}

		static Processor * NEW()
		{
			return new TestImageSinewave();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("wavelength", EMObject::FLOAT, "wavelength in equation sin(x*2*PI/wavelength - phase*180/PI)");
			d.put("axis", EMObject::STRING, "(optional) specify a major axis for asymmetric features, default x axis");
			d.put("phase", EMObject::FLOAT, "(optional) the phase in equation sin(x*2*PI/wavelength - phase*180/PI)");
			d.put("az", EMObject::FLOAT, "(optional) angle in degree. for 2D image, this is the rotated angle of the image, \
												in 3D image, it's az for euler angle. default is zero");
			d.put("alt", EMObject::FLOAT, "(optional) angle in degree. only in 3D case, alt for euler angle, default is zero");
			d.put("phi", EMObject::FLOAT, "(optional) angle in degree. only in 3D case, phi for euler angle, default is zero");
			return d;
		}
	};

	/**Replace a source image as a circular sine wave in specified wave length
	 *@param wavelength float this value is the d in function |sin(x/d)|
	 *@param axis char (optional) specify a major axis for asymmetric features
	 *@param c distance (optional) float between focus and the center of an ellipse
	 *@param phase degree (optional) phase for sine wave, default is 0
	 */
	class TestImageSinewaveCircular : public TestImageProcessor
	{
	public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "testimage.sinewave.circular";
		}

		virtual string get_desc() const
		{
			return "Replace a source image as a circular sine wave in specified wave length";
		}

		static Processor * NEW()
		{
			return new TestImageSinewaveCircular();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("wavelength", EMObject::FLOAT, "(required)this value is the d in function |sin(x/d)|, unit: pixel");
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
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "testimage.squarecube";
		}

		virtual string get_desc() const
		{
			return "Replace a source image as a square or cube depends on 2D or 3D of the source image";
		}

		static Processor * NEW()
		{
			return new TestImageSquarecube();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("edge_length", EMObject::FLOAT, "edge length of the square or cube, unit: pixel");
			d.put("axis", EMObject::STRING, "specify a major axis for asymmetric features");
			d.put("odd_edge", EMObject::FLOAT, "edge length for the asymmetric axis");
			d.put("fill", EMObject::INT, "Flag indicating if image is filled, default filled, 1 for filled, 0 for blank");
			return d;
		}
	};

	/**Generate an ellipse or ellipsoid image
	 *@param a equatorial radii along x axes
	 *@param b equatorial radii along y axes
	 *@param c polar radius
	 *@param transform Optionally transform the ellipse
	 *@param fill value you want to fill in ellipse, default to 1.0
	 */
	class TestImageEllipse : public TestImageProcessor
	{
	public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "testimage.ellipsoid";
		}

		virtual string get_desc() const
		{
			return "Insert an ellipse into the image.";
		}

		static Processor * NEW()
		{
			return new TestImageEllipse();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("a", EMObject::FLOAT, "equatorial radii along x axes");
			d.put("b", EMObject::FLOAT, "equatorial radii along y axes");
			d.put("c", EMObject::FLOAT, "polar radius");
			d.put("transform", EMObject::TRANSFORM, "Optionally transform the ellipse");
			d.put("fill", EMObject::FLOAT, "value you want to fill in ellipse, default to 1.0");
			return d;
		}
	};

	/**Generate an ellipse/ellipsoid image with an inner hollow ellipse/ellipsoid
	 *@param xwidth inner equatorial radii along x axes
	 *@param ywidth inner equatorial radii along y axes
	 *@param zwidth inner polar radius
	 *@param a outter equatorial radii along x axes
	 *@param b outter equatorial radii along y axes
	 *@param c outter polar radius
	 *@param width specify the width or specify each width explicitly - xwidth, ywidth, zwidth
	 *@param transform Optionally transform the ellipse
	 *@param fill value you want to fill in hollow ellipse, default to 1.0
	 */
	class TestImageHollowEllipse : public TestImageProcessor
		{
		public:
			virtual void process_inplace(EMData * image);

			virtual string get_name() const
			{
				return "testimage.ellipsoid.hollow";
			}

			virtual string get_desc() const
			{
				return "Insert a hollow ellipse into the image.";
			}

			static Processor * NEW()
			{
				return new TestImageHollowEllipse();
			}

			virtual TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("xwidth", EMObject::FLOAT, "inner equatorial radii along x axes");
				d.put("ywidth", EMObject::FLOAT, "inner equatorial radii along y axes");
				d.put("zwidth", EMObject::FLOAT, "inner polar radius");
				d.put("a", EMObject::FLOAT, "outter equatorial radii along x axes");
				d.put("b", EMObject::FLOAT, "outter equatorial radii along y axes");
				d.put("c", EMObject::FLOAT, "outter polar radius");
				d.put("width",EMObject::FLOAT, "width - specify the width or specify each width explicitly - xwidth, ywidth, zwidth");
				d.put("transform", EMObject::TRANSFORM, "Optionally transform the ellipse");
				d.put("fill", EMObject::FLOAT, "value you want to fill in hollow ellipse, default to 1.0");
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
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "testimage.circlesphere";
		}

		virtual string get_desc() const
		{
			return "Replace a source image as a circle or sphere depends on 2D or 3D of the source image";
		}

		static Processor * NEW()
		{
			return new TestImageCirclesphere();
		}

		virtual TypeDict get_param_types() const
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
	 *@param seed seed for random number generator
	 */
	class TestImageNoiseUniformRand : public TestImageProcessor
	{
	public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "testimage.noise.uniform.rand";
		}

		virtual string get_desc() const
		{
			return "Replace a source image as a uniform random noise, random number generated from gsl_rng_mt19937, the pixel value is [0, 1)";
		}

		static Processor * NEW()
		{
			return new TestImageNoiseUniformRand();
		}

		virtual TypeDict get_param_types() const
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
	 *@param sigma sigma value of gausian distributed noise, default is 0.5
	 *@param mean mean value of gausian distributed noise, default is zero
	 *@param seed mean value of gausian distributed noise, default is zero
	 */
	class TestImageNoiseGauss : public TestImageProcessor
	{
	public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "testimage.noise.gauss";
		}

		virtual string get_desc() const
		{
			return "Replace a source image as a random noise, the random value is gaussian distributed";
		}

		static Processor * NEW()
		{
			return new TestImageNoiseGauss();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("sigma", EMObject::FLOAT, "sigma value of gausian distributed noise, default is 0.5");
			d.put("mean", EMObject::FLOAT, "mean value of gausian distributed noise, default is zero.");
			d.put("seed", EMObject::INT, "the seed for random number generator, default is not to reseed.");

			return d;
		}
	};

	/** Replace a source image with a cylinder
	 *@param radius radius for the cylinder
	 *@param height height for the cylinder, by default it's the nz
	 * */
	class TestImageCylinder : public TestImageProcessor
	{
	public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "testimage.cylinder";
		}

		virtual string get_desc() const
		{
			return "Replace a source image as a cylinder";
		}

		static Processor * NEW()
		{
			return new TestImageCylinder();
		}

		virtual TypeDict get_param_types() const
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
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "filter.ccdnorm";
		}

		static Processor *NEW()
		{
			return new CCDNormProcessor();
		}

		virtual 	string get_desc() const
		{
			return "normalize the 4 quadrants of a CCD image";
		}

		virtual TypeDict get_param_types() const
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
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return "basis.wavelet";
		}

		static Processor *NEW()
		{
			return new WaveletProcessor();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("type", EMObject::STRING, "'daub', 'harr' or 'bspl'");
			d.put("dir", EMObject::INT, "1 for forward transform, -1 for inverse transform");
			d.put("ord", EMObject::INT, "Daubechies (4,6,8,...,20), for Harr (2), for B-Splines (103, 105, 202, 204, 206, 208, 301, 303, 305 307, 309)");
			return d;
		}

		virtual string get_desc() const
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
	 *@author David Woolford <woolford@bcm.edu>
	 *@date 01/10/2008
	 *@param biedgemean Mutually  exclusive of edgemean. Experimental. Causes the pixels in the masked out areas to take the average value of both the left and right edge pixel strips
	 *@param edgemean Mutually  exclusive of biedgemean. Masked pixels values assume the mean edge pixel value, independently, for both sides of the image
	 *@param angle The angle that the image is, with respect to the zero tilt image
	 *@param angle_fim Read fim as 'from image metadata' - this causes the altitude angle stored in by the image object (i.e. as extracted from the header, as currently stored in memory) to be used as the angle. This overrides the angle argument
	 *@param gauss_falloff Causes the edge masking to have a smooth Gaussian fall-off - this parameter specifies how many pixels the fall-off will proceed over. Default is 0
	 *@param gauss_sigma The sigma of the Gaussian function used to smooth the edge fall-off (functional form is exp(-(pixel distance)^2/sigma^2)
	 */
	class TomoTiltEdgeMaskProcessor : public Processor
	{
	public:
		virtual void process_inplace(EMData* image);

		virtual string get_name() const
		{
			return "tomo.tiltedgemask";
		}

		static Processor *NEW()
		{
			return new TomoTiltEdgeMaskProcessor();
		}

		virtual TypeDict get_param_types() const
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

		virtual string get_desc() const
		{
			return "Masks the part of the image which is not present in the 0-tilt image. Masked areas can be 0 or set to the edgemean (of the nearest or both edges). Masked areas can also have a Gaussian fall-off to make the appearance smooth.";
		}

	private:
		class GaussianFunctoid
		{
			public:
				GaussianFunctoid(const float sigma, const float mean = 0.0) : m_mean(mean), m_sigma_squared(sigma*sigma) {}
				~GaussianFunctoid() {}

				float operator()(const float distance )
				{
					return exp( -(distance-m_mean)*(distance-m_mean)/ (m_sigma_squared ));
				}
			private:
				float m_mean, m_sigma_squared;
		};

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
	 *@author David Woolford <woolford@bcm.edu>
	 *@date 02/11/2008
	 *@param angle The angle that the image is, with respect to the zero tilt image
	 *@param angle_fim Read fim as 'from image metadata' - this causes the altitude angle stored in by the image object (i.e. as extracted from the header, as currently stored in memory) to be used as the angle. This overrides the angle argument
	 */
	class TomoTiltAngleWeightProcessor : public Processor
	{
		public:
			virtual void process_inplace(EMData* image);

			virtual string get_name() const
			{
				return "tomo.tiltangleweight";
			}

			static Processor *NEW()
			{
				return new TomoTiltAngleWeightProcessor();
			}

			virtual TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("angle", EMObject::INT, "The angle that the image is, with respect to the zero tilt image");
				d.put("angle_fim",EMObject::BOOL,"Read fim as 'from image metadata' - this causes the altitude angle stored in by the image object (i.e. as extracted from the header, as currently stored in memory) to be used as the angle. This overrides the angle argument");
				return d;
			}

			virtual string get_desc() const
			{
				return "Weights the image by 1/cos(angle)";
			}

	};

	/** Perform a FFT transform by calling EMData::do_fft() and EMData::do_ift()
	 *@param dir 1 for forward transform, -1 for inverse transform, forward transform by default
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
	 *@param table a radial table for multiplication
	 *@exception ImageFormatException this filter only apply to real image
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

	/**Bins pixel values, similar to calculating a histogram. The histogram is comprised
	 * of 'nbins' bins, and the value assigned to each pixel in the bin is the midpoint
	 * of the bin's upper and lower limits. Defaults to 256 bins
	 *@param nbins The number of bins the pixel values will be compressed into
	 *@param debug Outputs debugging information (number of pixels per bin)
	 */
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


#ifdef EMAN2_USING_CUDA
	/** Cuda based constant multiplication processor
	 * @author David Woolford
	 * @date Feb 24 2009
	 * @param scale The amount to multiply each pixel by
	*/
	class CudaMultProcessor: public Processor
	{
		public:

			virtual void process_inplace(EMData * image);

			virtual string get_name() const
			{
				return "cuda.math.mult";
			}

			static Processor *NEW()
			{
				return new CudaMultProcessor();
			}

			virtual TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("scale", EMObject::FLOAT, "The amount to multiply each pixel by");
				return d;
			}

			virtual string get_desc() const
			{
				return "Multiplies each pixel by a constant value";
			}
		protected:
	};

	/** Cuda based auto correlation processor
	 * @author David Woolford
	 * @date Feb 24 2009
	 * @param with That which to perform the cross correlation with.
	 */
	class CudaCorrelationProcessor: public Processor
	{
		public:

			virtual void process_inplace(EMData * image);

			string get_name() const
			{
				return "cuda.correlate";
			}

			static Processor *NEW()
			{
				return new CudaCorrelationProcessor();
			}

			virtual TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("with", EMObject::EMDATA, "That which to perform the cross correlation with.");
				return d;
			}

			virtual string get_desc() const
			{
				return "Performs Fourier based cross correlation on the GPU";
			}
		protected:
	};


	/* class CUDA kmeans processor
	* 02/13/2009 JB
	* python wrap for cuda_kmeans.cu
	*/
	class CUDA_kmeans {
	public:
		CUDA_kmeans();
		~CUDA_kmeans();
		int setup(int  extm, int extN, int extK, float extF, float extT0, int extmaxite, int extrnd);
		void append_flat_image(EMData* im, int pos);
		int kmeans();
		vector<EMData*> get_averages();
		vector<int> get_partition();
		void set_K(int valK);
		void set_rnd(int valrnd);
		Dict get_info();
	private:
		// params
		int m;
		int N;
		int K;
		int maxite;
		int nb_part;
		float F;
		float T0;
		int rnd;
		// host memory
		float* h_IM;
		float* h_INFO;
		float* h_AVE;
		unsigned short int* h_ASG;
	};
#endif //EMAN2_USING_CUDA

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
