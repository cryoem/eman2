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
 * it under the terms of the GNU General Public License as published by.edge
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
#include "emdata.h"
#include "gorgon/skeletonizer.h"

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
			d.put("cutoff_freq", EMObject::FLOAT, "1/Resolution in 1/A (0 - 1 / 2*apix). eg - a 20 A filter is cutoff_freq=0.05");
			d.put("apix", EMObject::FLOAT, " Override A/pix in the image header (changes x,y and z)");
			d.put("return_radial", EMObject::BOOL, "Return the radial filter function as an attribute (filter_curve)");
			return d;
		}

	  protected:
		  virtual void preprocess(EMData *image) {
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
		  virtual void create_radial_func(vector < float >&radial_mask) const = 0;
	};

	/** Same as FourierProcessor, except first computes the current image radial power spectrum and passes it to the processor
	 * (no radial oversampling, number of elements = radius)
	 * @param cutoff_abs Processor radius in terms of Nyquist (0-.5).
	 * @param cutoff_pixels Width in Fourier pixels (0 - size()/2).
	 * @param cutoff_freq Resolution in 1/A (0 - 1 / size*apix).
	 * @param apix Override A/pix in the image header (changes x,y and z).
	 */
	class FourierAnlProcessor:public Processor
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
			d.put("cutoff_freq", EMObject::FLOAT, "1/Resolution in 1/A (0 - 1 / 2*apix). eg - a 20 A filter is cutoff_freq=0.05");
			d.put("apix", EMObject::FLOAT, " Override A/pix in the image header (changes x,y and z)");
			d.put("return_radial", EMObject::BOOL, "Return the radial filter function as an attribute (filter_curve)");
			d.put("interpolate", EMObject::BOOL, "Whether or not to interpolate the radial scaling function. Default=false. Prb should be true.");
			return d;
		}

	  protected:
		  virtual void preprocess(EMData *) {}
		  virtual void create_radial_func(vector < float >&radial_mask,EMData *image) const = 0;
	};

	/** Similar to FourierProcessor, but enhances or compresses azimuthal contrast rather than the
	 * typical radial linear filter
	 * @param az_scale Scale factor, >1 enhances contrast <1 decreases
	 */
	class AzSharpProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);
		
		string get_name() const
		{
			return NAME;
		}
		
		static Processor *NEW()
		{
			return new AzSharpProcessor();
		}

		string get_desc() const
		{
			return "Typical linear filters are radial, but certain operations can lead to inhomogeneities so the balance between radial and azimuthal power is off. This processor permits enhancing or suppressing azimuthal contrast in Fourier space.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("az_scale", EMObject::FLOAT, "Scale factor, >1 enhances contrast <1 decreases");
			return d;
		}
		
		static const string NAME;

	};

	
	/** Evaluate individual particle images using a tenchique similar to that used for CTF evaluation
	 */
	class SNREvalProcessor:public Processor
	{
		public:
		string get_name() const
		{
			return NAME;
		}

		void process_inplace(EMData * image);

		void set_params(const Dict & new_params)
		{
			params = new_params;
//			sum = params["sum"];
//			dosqrt = params["sqrt"];
//			printf("%s %f\n",params.keys()[0].c_str(),lowpass);
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
//			d.put("sum", EMObject::EMDATA, "Adds the weights to sum for normalization");
//			d.put("sqrt", EMObject::INT, "Weights using sqrt of the amplitude if set");
			return d;
		}

		static Processor *NEW()
		{
			return new SNREvalProcessor();
		}

		string get_desc() const
		{
			return "Evaluates the SNR of the particle using a masking method similar to that used in the CTF analysis process. The image is not changed. The resulting value is placed in the particle dictionary as eval_maskedsnr";
		}

		static const string NAME;

		protected:
		EMData *sum;
		int dosqrt;
	};

	/**Zeroes the values on the X=0 and y=0 Fourier axes (except x=y=0)
	 */
	class Axis0FourierProcessor:public Processor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}

		void process_inplace(EMData * image);

		void set_params(const Dict & new_params)
		{
			params = new_params;
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("x", EMObject::INT, "If set, zeroes along X axis. Default True.");
			d.put("y", EMObject::INT, "If set, zeroes along Y axis. Default True.");
			d.put("neighbor", EMObject::INT, "If set, interpolates neighbor values rather than zeroing");
			d.put("neighbornorm", EMObject::INT, "In neighbor mode it sums the neighboring 2 (complex) pixels, then divides by this factor. default = sqrt(2), which is good for very noisy images");
			return d;
		}

		static Processor *NEW()
		{
			return new Axis0FourierProcessor();
		}

		string get_desc() const
		{
			return "Sets values along X/Y Fourier axes to 0, except origin";
		}

		static const string NAME;

		protected:
	};

	/**Zeroes the values on the X=0 and y=0 Fourier axes (except x=y=0)
	 */
	class GaussZFourierProcessor:public Processor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}

		void process_inplace(EMData * image);

		void set_params(const Dict & new_params)
		{
			params = new_params;
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("cutoff_abs", EMObject::FLOAT, "Processor radius in terms of Nyquist (0-.5)");
			d.put("cutoff_pixels", EMObject::FLOAT, " Width in Fourier pixels (0 - size()/2)");
			d.put("cutoff_freq", EMObject::FLOAT, "1/Resolution in 1/A (0 - 1 / 2*apix). eg - a 20 A filter is cutoff_freq=0.05");
			d.put("apix", EMObject::FLOAT, " Override A/pix in the image header (changes x,y and z)");
			return d;
		}

		static Processor *NEW()
		{
			return new GaussZFourierProcessor();
		}

		string get_desc() const
		{
			return "Applies a Gaussian lowpass filter (or its inverse), but only along the Z axis. May be useful in anisotropic filtering of tomograms.";
		}

		static const string NAME;

		protected:
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
			return NAME;
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

		static const string NAME;

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
				return NAME;
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

			static const string NAME;
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
			return NAME;
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

		static const string NAME;
	};

	class YGradientProcessor : public Processor
	{
		public:
			YGradientProcessor() {}

			string get_name() const
			{
				return NAME;
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

			static const string NAME;
	};

	class ZGradientProcessor : public Processor
	{
		public:
			ZGradientProcessor() {}

			string get_name() const
			{
				return NAME;
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

			static const string NAME;
	};

	/** Computes the image divergence using David's partial derivative processors
	 *
	 *@author James Michael Bell
	 *@date 06/26/2015
	 */
	class ImageDivergenceProcessor : public Processor
	{
	 public:
		ImageDivergenceProcessor() {}

		string get_name() const
		{
			return NAME;
		}

		void process_inplace(EMData *image);

		static Processor *NEW()
		{
			return new ImageDivergenceProcessor();
		}

		string get_desc() const
		{
			return "Determines the divergence of a 2D image.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			return d;
		}

		static const string NAME;
	};

	/** Determines the magnitude of an approximate image gradient using David's image gradient processors
	 *
	 *@author James Michael Bell
	 *@date 06/26/2015
	 */
	class GradientMagnitudeProcessor : public Processor
	{
	 public:
		GradientMagnitudeProcessor() {}

		string get_name() const
		{
			return NAME;
		}

		void process_inplace(EMData *image);

		static Processor *NEW()
		{
			return new GradientMagnitudeProcessor();
		}

		string get_desc() const
		{
			return "Determines the magnitude of the gradient of a 2D image.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			return d;
		}

		static const string NAME;
	};

	/** Determines the direction of an approximate image gradient using David's image gradient processors
	 *
	 *@author James Michael Bell
	 *@date 06/26/2015
	 */
	class GradientDirectionProcessor : public Processor
	{
	 public:
		GradientDirectionProcessor() {}

		string get_name() const
		{
			return NAME;
		}

		void process_inplace(EMData *image);

		static Processor *NEW()
		{
			return new GradientDirectionProcessor();
		}

		string get_desc() const
		{
			return "Determines the direction of the gradient of a 2D image.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			return d;
		}

		static const string NAME;
	};

	/** Determines the direction of an approximate image laplacian using David's image gradient processors
	 *
	 *@author James Michael Bell
	 *@date 06/26/2015
	 */
	class LaplacianMagnitudeProcessor : public Processor
	{
	 public:
		LaplacianMagnitudeProcessor() {}

		string get_name() const
		{
			return NAME;
		}

		void process_inplace(EMData *image);

		static Processor *NEW()
		{
			return new LaplacianMagnitudeProcessor();
		}

		string get_desc() const
		{
			return "Determines the magnitude of the laplacian of a 2D image.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			return d;
		}

		static const string NAME;
	};

	/** Determines the direction of an approximate image laplacian using David's image gradient processors
	 *
	 *@author James Michael Bell
	 *@date 06/26/2015
	 */
	class LaplacianDirectionProcessor : public Processor
	{
	 public:
		LaplacianDirectionProcessor() {}

		string get_name() const
		{
			return NAME;
		}

		void process_inplace(EMData *image);

		static Processor *NEW()
		{
			return new LaplacianDirectionProcessor();
		}

		string get_desc() const
		{
			return "Determines the direction of the laplacian of a 2D image.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			return d;
		}

		static const string NAME;
	};

	/** Determines the second derivative in the gradient direction using David's image gradient processors
	 *
	 *@author James Michael Bell
	 *@date 06/26/2015
	 */
	class SDGDProcessor: public Processor
	{
	 public:
		SDGDProcessor() {}

		string get_name() const
		{
			return NAME;
		}

		void process_inplace(EMData *image);

		static Processor *NEW()
		{
			return new SDGDProcessor();
		}

		string get_desc() const
		{
			return "Determines the second derivative of a 2D image in the gradient direction.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			return d;
		}

		static const string NAME;
	};

	/** Sets pixel values in a binary image equal to their element wise manhattan distance.
	 * Credit for this code goes to Stephen Ostermiller (http://blog.ostermiller.org/dilate-and-erode).
	 *
	 *@author James Michael Bell
	 *@date 06/27/2015
	 */
	class ManhattanDistanceProcessor: public Processor
	{
	 public:
		ManhattanDistanceProcessor() {}

		string get_name() const
		{
			return NAME;
		}

		virtual EMData* process(const EMData * const image);
		
		void process_inplace(EMData * image);
		
		static Processor *NEW()
		{
			return new ManhattanDistanceProcessor();
		}

		string get_desc() const
		{
			return "Sets pixel values in a binary image equal to their element wise manhattan distance.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			return d;
		}

		static const string NAME;
	};


	/** Performs a morphological dilation operation on an image.
	 *
	 *@author James Michael Bell
	 *@date 06/27/2015
	 */
	class BinaryDilationProcessor: public Processor
	{
	 public:
		BinaryDilationProcessor() {}

		string get_name() const
		{
			return NAME;
		}

		virtual EMData* process(const EMData * const image);
		
		void process_inplace(EMData *image);

		static Processor *NEW()
		{
			return new BinaryDilationProcessor();
		}

		string get_desc() const
		{
			return "Performs a morphological dilation of a (binary) 2D image.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("k", EMObject::INT, "The number of pixels to dilate the input image.");
			d.put("selem",EMObject::EMDATA, "The structuring element with which you want to dilate.");
			return d;
		}

		static const string NAME;
	};

	/** Performs a morphological erosion operation on an image.
	 *
	 *@author James Michael Bell
	 *@date 06/27/2015
	 */
	class BinaryErosionProcessor: public Processor
	{
	 public:
		BinaryErosionProcessor() {}

		string get_name() const
		{
			return NAME;
		}

		virtual EMData* process(const EMData * const image);
		
		void process_inplace(EMData *image);

		static Processor *NEW()
		{
			return new BinaryErosionProcessor();
		}

		string get_desc() const
		{
			return "Performs a morphological k-pixel erosion of a (binary) 2D image.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("k", EMObject::INT, "The number of pixels to close the input image.");
			return d;
		}

		static const string NAME;
	};


	/** Performs a morphological closing operation on an image.
	 *
	 *@author James Michael Bell
	 *@date 06/27/2015
	 */
	class BinaryClosingProcessor: public Processor
	{
	 public:
		BinaryClosingProcessor() {}

		string get_name() const
		{
			return NAME;
		}

		virtual EMData* process(const EMData * const image);
		
		void process_inplace(EMData *image);

		static Processor *NEW()
		{
			return new BinaryClosingProcessor();
		}

		string get_desc() const
		{
			return "Performs a morphological k-pixel opening of a (binary) 2D image.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("k", EMObject::INT, "The number of pixels to close the input image.");
			return d;
		}

		static const string NAME;
	};


	/** Performs a morphological opening operation on an image.
	 *
	 *@author James Michael Bell
	 *@date 06/27/2015
	 */
	class BinaryOpeningProcessor: public Processor
	{
	 public:
		BinaryOpeningProcessor() {}

		string get_name() const
		{
			return NAME;
		}

		virtual EMData* process(const EMData * const image);
		
		void process_inplace(EMData *image);

		static Processor *NEW()
		{
			return new BinaryOpeningProcessor();
		}

		string get_desc() const
		{
			return "Performs a morphological k-pixel closing of a (binary) 2D image.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("k", EMObject::INT, "The number of pixels to close the input image.");
			return d;
		}

		static const string NAME;
	};


	/** Computes an internal morphological gradient operation on an image.
	 *
	 *@author James Michael Bell
	 *@date 06/27/2015
	 */
	class BinaryInternalGradientProcessor: public Processor
	{
	 public:
		BinaryInternalGradientProcessor() {}

		string get_name() const
		{
			return NAME;
		}

		virtual EMData* process(const EMData * const image);
		
		void process_inplace(EMData *image);

		static Processor *NEW()
		{
			return new BinaryInternalGradientProcessor();
		}

		string get_desc() const
		{
			return "Computes an internal morphological graduent using k-pixel-width operations on a (binary) 2D image.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("k", EMObject::INT, "The number of pixels to close the input image.");
			return d;
		}

		static const string NAME;
	};
	
	
	/** Computes an external morphological gradient operation on an image.
	 *
	 *@author James Michael Bell
	 *@date 06/27/2015
	 */
	class BinaryExternalGradientProcessor: public Processor
	{
	 public:
		BinaryExternalGradientProcessor() {}

		string get_name() const
		{
			return NAME;
		}

		virtual EMData* process(const EMData * const image);
		
		void process_inplace(EMData *image);

		static Processor *NEW()
		{
			return new BinaryExternalGradientProcessor();
		}

		string get_desc() const
		{
			return "Computes an external morphological graduent using k-pixel-width operations on a (binary) 2D image.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("k", EMObject::INT, "The number of pixels to close the input image.");
			return d;
		}

		static const string NAME;
	};
	
	
	/** Computes the morphological gradient operation on an image.
	 *
	 *@author James Michael Bell
	 *@date 06/27/2015
	 */
	class BinaryMorphGradientProcessor: public Processor
	{
	 public:
		BinaryMorphGradientProcessor() {}

		string get_name() const
		{
			return NAME;
		}

		virtual EMData* process(const EMData * const image);
		
		void process_inplace(EMData *image);

		static Processor *NEW()
		{
			return new BinaryMorphGradientProcessor();
		}

		string get_desc() const
		{
			return "Computes the morphological graduent using k-pixel-width operations on a (binary) 2D image.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("k", EMObject::INT, "The number of pixels to close the input image.");
			return d;
		}

		static const string NAME;
	};
	
	
	/** Performs a morphological top hat operation on an image.
	 *
	 *@author James Michael Bell
	 *@date 06/27/2015
	 */
	class BinaryTopHatProcessor: public Processor
	{
	 public:
		BinaryTopHatProcessor() {}

		string get_name() const
		{
			return NAME;
		}

		virtual EMData* process(const EMData * const image);
		
		void process_inplace(EMData *image);

		static Processor *NEW()
		{
			return new BinaryTopHatProcessor();
		}

		string get_desc() const
		{
			return "Performs a morphological top hat operation on a (binary) 2D image.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("k", EMObject::INT, "The number of pixels to close the input image.");
			return d;
		}

		static const string NAME;
	};


	/** Performs a morphological black hat operation on an image.
	 *
	 *@author James Michael Bell
	 *@date 06/27/2015
	 */
	class BinaryBlackHatProcessor: public Processor
	{
	 public:
		BinaryBlackHatProcessor() {}

		string get_name() const
		{
			return NAME;
		}

		virtual EMData* process(const EMData * const image);
		
		void process_inplace(EMData *image);

		static Processor *NEW()
		{
			return new BinaryBlackHatProcessor();
		}

		string get_desc() const
		{
			return "Performs a morphological black hat operation on a (binary) 2D image.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("k", EMObject::INT, "The number of pixels to close the input image.");
			return d;
		}

		static const string NAME;
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
			return NAME;
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

		static const string NAME;

		protected:
		int bgsize;
	};

	/** Segment a volume about:homeinto subvolumes based on a center separation value. For linear densities
	 * such as skeletons this should fill linear regions with uniformly separated points
	 *
	 *@author Steve Ludtke
	 *@date 2010/07/14
	 */
	class DistanceSegmentProcessor:public Processor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}

		virtual EMData* process(const EMData * const image);
		void process_inplace( EMData * image);

		TypeDict get_param_types() const
		{
			TypeDict d ;
			d.put("thr",EMObject::FLOAT,"Optional : Isosurface threshold value. Pixels below this will not be segment centers (default = 0.9)");
			d.put("minsegsep",EMObject::FLOAT,"Required: Minimum segment separation in pixels. Segments too close will trigger a reseed");
			d.put("maxsegsep",EMObject::FLOAT,"Required: Maximum segment separation in pixels. Segments too close will trigger a reseed");
			d.put("verbose",EMObject::INT,"Be verbose while running");
			return d;
		}

		static Processor *NEW()
		{
			return new DistanceSegmentProcessor();
		}

		string get_desc() const
		{
			return "Segments a volume into pieces separated by distances in the specified range.";
		}

		static const string NAME;

	};

	/** Segment a volume into ~n subvolumes using K-means classification
	 *
	 *@author Steve Ludtke
	 *@date 2008/11/03
	 *@param ctf[in] A Ctf object to use
	 */
	class KmeansSegmentProcessor:public Processor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}

		virtual EMData* process(const EMData * const image);
		void process_inplace( EMData * image);

		TypeDict get_param_types() const
		{
			TypeDict d ;
			d.put("nseg", EMObject::INT, "Number of segments to divide the image into. default=12" );
			d.put("thr",EMObject::FLOAT,"Isosurface threshold value. Pixels below this will not be segmented");
			d.put("ampweight",EMObject::INT,"If set, will weight centers by voxel amplitude. default = 1");
			d.put("maxsegsize",EMObject::FLOAT,"Maximum radial distance from segment center to member voxel. Default=10000");
			d.put("minsegsep",EMObject::FLOAT,"Minimum segment separation. Segments too close will trigger a reseed");
			d.put("maxiter",EMObject::FLOAT,"Maximum number of iterations to run before stopping. Default=100");
			d.put("maxvoxmove",EMObject::FLOAT,"Maximum number of voxels that can move before quitting. Default=25");
			d.put("verbose",EMObject::INT,"Be verbose while running");
			/**
			 *An option for pseudoatom generation in pathwalker. Instead of random seeding, seed on the gird initially.
			 *@author Muyuan Chen
			 *@date 2014/06/05
	         */
			d.put("pseudoatom",EMObject::BOOL,"Doing pseudoatom generation");
			d.put("sep",EMObject::FLOAT,"Separation distance, used only in pseudoatom generation. Default=3.78");
			return d;
		}

		static Processor *NEW()
		{
			return new KmeansSegmentProcessor();
		}

		string get_desc() const
		{
			return "Performs K-means segmentation on a volume. Note that this method uses random seeds, and thus will return different results each time it is run. Returned map contains number of segment for each voxel (or 0 for unsegmented voxels). Segmentation centers are stored in 'segmentcenters' attribute, consisting of a list of 3n floats in x,y,z triples.";
		}

		static const string NAME;

	};


	/** CTF simulation processor. Takes individual CTF parameters, suitable for use with programs like
	 * e2filtertool.py. Can use an internal noise profile or an external profile from a text file.
	 *@param defocus[in]	Defocus in microns (underfocus positive)
	 *@param ampcont[in]	% amplitude contrast (0-100)
	 *@param bfactor[in]	B-factor in A^2, uses MRC convention rather than EMAN1 convention
	 *@param noiseamp[in]	Amplitude of the added empirical pink noise
	 *@param noiseampwhite[in]	Amplitude of added white noise
	 *@param voltage[in]	Microscope voltage in KV
	 *@param cs[in]	Cs of microscope in mm
	 *@param apix[in]	A/pix of data
	 *@author Steve Ludtke
	 *@date 2011/11/03
	 */
	class CtfSimProcessor:public Processor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}

		virtual EMData* process(const EMData * const image);

		void process_inplace(EMData *image);

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("defocus", EMObject::FLOAT, "Defocus in microns (underfocus positive)");
			d.put("ampcont", EMObject::FLOAT, "% amplitude contrast (0-100)");
			d.put("bfactor", EMObject::FLOAT, "B-factor in A^2, uses MRC convention rather than EMAN1 convention");
			d.put("noiseamp", EMObject::FLOAT, "Amplitude of the added empirical pink noise");
			d.put("noiseampwhite", EMObject::FLOAT, "Amplitude of added white noise");
			d.put("voltage", EMObject::FLOAT, "Microscope voltage in KV");
			d.put("cs", EMObject::FLOAT, "Cs of microscope in mm");
			d.put("apix", EMObject::FLOAT, "A/pix of data");
			return d;
		}

		static Processor *NEW()
		{
			return new CtfSimProcessor();
		}

		string get_desc() const
		{
			return "Applies a simulated CTF with noise to an image. The added noise is either white or based on an empirical curve generated from cryoEM data. ";
		}

		static const string NAME;

//		protected:
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
			return NAME;
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

		static const string NAME;

		protected:
		Ctf *ctf;
	};

	class LinearRampFourierProcessor:public FourierProcessor
	{
		public:
			virtual string get_name() const
			{
				return NAME;
			}

			virtual string get_desc() const
			{
				return "";
			}

			static Processor *NEW()
			{
				return new LinearRampFourierProcessor();
			}

			static const string NAME;

		protected:
			virtual void create_radial_func(vector < float >&radial_mask) const ;
	};


	   /**Lowpass Phase Randomization processor applied in Fourier space.
         */
        class LowpassRandomPhaseProcessor:public FourierProcessor
        {
          public:
                string get_name() const
                { return NAME; }
                static Processor *NEW() { return new LowpassRandomPhaseProcessor(); }
                string get_desc() const
                {
                        return "Above the cutoff frequency, phases will be completely randomized, but amplitudes will be unchanged. Used for testing for noise bias. If you can reconstruct information that isn't there, then you have noise bias problems.";
                }
				void process_inplace(EMData * image);
		  		void create_radial_func(vector < float >&radial_mask) const;

                static const string NAME;
        };


	/**processor radial function: if lowpass > 0, f(x) = exp(-x*x/(lowpass*lowpass)); else f(x) = exp(x*x/(lowpass*lowpass))
	 */
	class LowpassAutoBProcessor:public FourierAnlProcessor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new LowpassAutoBProcessor();
		}

		string get_desc() const
		{
			return "This processor sharpens a map based on the concept that the power spectrum should be roughly flat over the ~15 A-Nyquist resolution range, then combines this inverse B-factor with the specified low-pass Gaussian filter parameters to produce a single aggregate Gaussian filter.";
		}

		static const string NAME;

		virtual TypeDict get_param_types() const
		{
			TypeDict d = FourierAnlProcessor::get_param_types();
			d.put("bfactor", EMObject::FLOAT, "B-factor in terms of e^-(B s^2/4)");
			d.put("noisecutoff", EMObject::FLOAT, "Spatial frequency past which inverse-B will not be applied, default=1/6A");
//			d.put("adaptnoise", EMObject::INT, "Dual linear fit separating lower resolution signal from higher resolution noise. Noise region not upweighted.");
			d.put("return_radial", EMObject::BOOL, "Return the radial filter function as an attribute (filter_curve)");
			d.put("verbose", EMObject::INT, "Print information about the determined B-factor");
			return d;
		}

	  protected:
		virtual void preprocess(EMData * image) {
			if(params.has_key("apix")) {
				image->set_attr("apix_x", (float)params["apix"]);
				image->set_attr("apix_y", (float)params["apix"]);
				image->set_attr("apix_z", (float)params["apix"]);
			}
			float apix=(float)image->get_attr("apix_x");

			const Dict dict = image->get_attr_dict();
			if (params.has_key("cutoff_abs")) {
				params["bfactor"] = pow(apix/(float)params["cutoff_abs"],2.0f);
			}
			else if( params.has_key("cutoff_freq") ) {
				float val =  (float)params["cutoff_freq"] * apix;
				params["cutoff_abs"] = val;
				params["bfactor"] = pow(apix/(float)params["cutoff_abs"],2.0f);
			}
			else if( params.has_key("cutoff_pixels") ) {
				float val = 0.5f*(float)params["cutoff_pixels"] / (float)dict["nx"];
				params["cutoff_abs"] = val;
				params["bfactor"] = pow(apix/(float)params["cutoff_abs"],2.0f);
			}
		}

		void create_radial_func(vector < float >&radial_mask,EMData *image) const;
	};

	/** This processor attempts to remove the low resolution peak present in all cryoEM data
	 */
	class HighpassAutoPeakProcessor:public FourierAnlProcessor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}
		static Processor *NEW()
		{
			return new HighpassAutoPeakProcessor();
		}

		string get_desc() const
		{
			return "Attempts to automatically remove the low resolution peak present in virtually all cryoEM data.";
		}

		static const string NAME;

	  protected:
		void create_radial_func(vector < float >&radial_mask, EMData *image) const;
		  virtual void preprocess(EMData * image);
		  float highpass;
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
			return NAME;
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

		static const string NAME;

	  protected:
		void create_radial_func(vector < float >&radial_mask) const;

	  private:
		float intercept;
		float slope;
	};

	/**processor radial function: f(x) = ((x^2 - s^2)/s^4)e^-(x^2/2s^2)
	 *@param sigma LoG sigma
	 */
	class LoGFourierProcessor:public FourierProcessor
	{
	  public:
		LoGFourierProcessor():sigma(0)
		{
		}

		string get_name() const
		{
			return NAME;
		}
		static Processor *NEW()
		{
			return new LoGFourierProcessor();
		}

		string get_desc() const
		{
			return "processor radial function: f(x) = ((x^2 - s^2)/s^4)e^-(x^2/2s^2)";
		}

		void set_params(const Dict & new_params)
		{
			params = new_params;
			sigma = params["sigma"];
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("sigma", EMObject::FLOAT, "LoG sigma");
			return d;
		}

		static const string NAME;

	  protected:
		void create_radial_func(vector < float >&radial_mask) const;

	  private:
		float sigma;
	};

	/**processor radial function: f(x) = 1/sqrt(2*pi)*[1/sigma1*exp-(x^2/2*sigma1^2) - 1/sigma2*exp-(x^2/2*sigma2^2)]
	 *@param sigma1 DoG sigma1
	 *@param sigma1 DoG sigma2
	 */
	class DoGFourierProcessor:public FourierProcessor
	{
	  public:
		DoGFourierProcessor():sigma1(0), sigma2(0)
		{
		}

		string get_name() const
		{
			return NAME;
		}
		static Processor *NEW()
		{
			return new DoGFourierProcessor();
		}

		string get_desc() const
		{
			return "processor radial function: f(x) = 1/sqrt(2*pi)*[1/sigma1*exp-(x^2/2*sigma1^2) - 1/sigma2*exp-(x^2/2*sigma2^2)]";
		}

		void set_params(const Dict & new_params)
		{
			params = new_params;
			sigma1 = params["sigma1"];
			sigma2 = params["sigma2"];
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("sigma1", EMObject::FLOAT, "DoG sigma1");
			d.put("sigma2", EMObject::FLOAT, "DoG sigma2");
			return d;
		}

		static const string NAME;

	  protected:
		void create_radial_func(vector < float >&radial_mask) const;

	  private:
		float sigma1;
		float sigma2;
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

		virtual void set_params(const Dict & new_params)
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
			return NAME;
		}
		static Processor *NEW()
		{
			return new AbsoluateValueProcessor();
		}

		static const string NAME;

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

	/**f(x) = floor(x)
	 */
	class FloorValueProcessor:public RealPixelProcessor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}
		static Processor *NEW()
		{
			return new FloorValueProcessor();
		}

		static const string NAME;

	  protected:
		void process_pixel(float *x) const
		{
			*x = floor(*x);
		}

		string get_desc() const
		{
			return "f(x) = floor(x)";
		}
	};

	/** This processor can be used to correct errors when reading signed data as unsigned and vice-versa
	 */
	class FixSignProcessor:public RealPixelProcessor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}
		static Processor *NEW()
		{
			return new FixSignProcessor();
		}

		static const string NAME;

		void set_params(const Dict & new_params)
			{
				params = new_params;
				if (params.set_default("byte_stou",0)) mode=1;
				else if (params.set_default("byte_utos",1)) mode=2;
				else mode=0;
			}
		
		TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("byte_stou", EMObject::BOOL, "8 bit data signed read -> unsigned");
				d.put("byte_utos", EMObject::BOOL, "8 bit data unsigned read -> signed");
				return d;
			}

	  protected:
		void process_pixel(float *x) const
		{
			switch (mode) {
			case 1:
				if (*x<0) *x+=256;
				break;
			case 2:
				if (*x>127) *x-=256;
				break;
			}
				
		}

		string get_desc() const
		{
			return "Fixes errors with reading signed/unsigned data. Need to specify the correct mode.";
		}
		
		int mode;
	};


	/**f(x) = 0 if x = 0; f(x) = 1 if x != 0
	 */
	class BooleanProcessor:public RealPixelProcessor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}
		static Processor *NEW()
		{
			return new BooleanProcessor();
		}

		string get_desc() const
		{
			return "f(x) = 0 if x = 0; f(x) = 1 if x != 0;";
		}

		static const string NAME;

	  protected:
		void process_pixel(float *x) const
		{
			if (*x != 0)
			{
				*x = 1.0;
			}
		}
	};

	/**Reciprocal image as if f(x) != 0: f(x) = 1/f(x) else: f(x) = zero_to
	 *@param zero_to  Inverted zero values are set to this value, default is 0
	 */
	class RecipCarefullyProcessor:public RealPixelProcessor
	{
		public:
			string get_name() const
			{
				return NAME;
			}
			static Processor *NEW()
			{
				return new RecipCarefullyProcessor();
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

			static const string NAME;

		protected:
			void process_pixel(float *x) const
			{
				if (*x == 0.0) *x = zero_to;
				else *x = 1.0f/(*x);
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
			return NAME;
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

		static const string NAME;

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
			return NAME;
		}
		static Processor *NEW()
		{
			return new ValueSquaredProcessor();
		}


		string get_desc() const
		{
			return "f(x) = x * x;";
		}

		static const string NAME;

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
			return NAME;
		}
		static Processor *NEW()
		{
			return new ValueSqrtProcessor();
		}

		string get_desc() const
		{
			return "f(x) = sqrt(x)";
		}

		static const string NAME;

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
				return NAME;
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

			static const string NAME;

		protected:
			inline void process_pixel(float *x) const
			{
				if (*x < value) {
					*x = 0;
				}
			}
	};

	/**f(x) = x if x <= maxval; f(x) = 0 if x > maxval
	 * @param maxval
	 */
	class AboveToZeroProcessor:public RealPixelProcessor
	{
	public:
		string get_name() const
		{
			return NAME;
		}
		static Processor *NEW()
		{
			return new AboveToZeroProcessor();
		}
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("maxval", EMObject::FLOAT, "Everything above this value is set to zero");
			return d;
		}

		string get_desc() const
		{
			return "f(x) = x if x <= maxval; f(x) = 0 if x > maxval.";
		}

		static const string NAME;

	protected:
		inline void process_pixel(float *x) const
		{
			if (*x > value) {
				*x = 0;
			}
		}
	};

	/**Rotate by 180 using pixel swapping, works for 2D only
	 * @author David Woolford
	 * @date March 21, 2014
	 */
	class AddShapeProcessor:public Processor
	{
		public:
			string get_name() const
			{
				return NAME;
			}

			static Processor *NEW()
			{
				return new AddShapeProcessor();
			}

			void process_inplace(EMData* image);

			TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("shape", EMObject::STRING, "Name of the shape to add (");
				d.put("x", EMObject::FLOAT, "X coordinate of object center");
				d.put("y", EMObject::FLOAT, "Y coordinate of object center");
				d.put("z", EMObject::FLOAT, "Z coordinate of object center");
				d.put("size1", EMObject::FLOAT, "Size of the object. Meaning varies by shape.");
				d.put("size2", EMObject::FLOAT, "2nd axis size of the object. Meaning varies by shape.");
				d.put("size3", EMObject::FLOAT, "3rd axis size of the object. Meaning varies by shape.");
				d.put("val1", EMObject::FLOAT, "First pixel value. Meaning varies by shape.");
				d.put("val2", EMObject::FLOAT, "2nd pixel value. Meaning varies with shape");

				return d;
			}

			string get_desc() const
			{
				return "Adds a specified shape to a volume.";
			}

			static const string NAME;
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
				return NAME;
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

			static const string NAME;
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
				return NAME;
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
				d.put("zerocorners",EMObject::INT,"If set, corners (anything beyond radius/2-1) may be zeroed out in real or Fourier space. This will produce a considerable speedup in Fourier rotations. ");
				return d;
			}

			virtual string get_desc() const
			{
				return "The image is transformed using Transform parameter.";
			}

			static const string NAME;

			float* transform(const EMData* const image, const Transform& t) const;
		private:
			// This function became redundant if favor of EMData::scale_pixel
			//void update_emdata_attributes(EMData* const image, const Dict& attr_dict, const float& scale) const;


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
				return NAME;
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

			static const string NAME;

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

	 /** Applies a symmetry to a 3D model. The model must be on aligned to its symmetry axis(via align3d or other mechanism)
	 *@author Steve Ludtke and John Flanagan
	 *@date June 2011
	 *@param sym A string specifying the symmetry under which to do the alignment
	 */
	class ApplySymProcessor:public Processor
	{
		public:
			virtual string get_name() const
			{
				return NAME;
			}

			static Processor *NEW()
			{
				return new ApplySymProcessor();
			}

			virtual void process_inplace(EMData* image);

			virtual EMData* process(const EMData* const image);

			virtual TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("sym", EMObject::STRING, "The symmetry under which to do the alignment, Default=c1" );
				d.put("averager", EMObject::STRING, "Name of an Averager to use. default=mean" );
				return d;
			}

			virtual string get_desc() const
			{
				return "Symmetry is imposed on a 2-D image (Cn only) or 3-D volume";
			}

			static const string NAME;

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
				return NAME;
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
				/**d.put("clipx", EMObject::INT, "The length of the output x dimension. Exclusive of the clip.");
				*/
				return d;
			}

			virtual string get_desc() const
			{
				return "The image is scaled with the clip variable in mind, being sure to preserve as much pixel information as possible.";
			}

			static const string NAME;
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
			return NAME;
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
			d.put("tozero", EMObject::BOOL, "Replace outlying pixels values with zero" );
			return d;
		}

		string get_desc() const
		{
			return "This function clamps the min and max vals in the image at minval and maxval, respectively. In a sense this a bi-truncation of the data.";
		}

		static const string NAME;

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
				return NAME;
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
				d.put("tozero", EMObject::BOOL, "Replace outlying pixels values with zero" );
				return d;
			}

			void process_inplace(EMData *image);

			string get_desc() const
			{
				return "This function clamps the min and max vals in the image at minval and maxval at mean-n*sigma and mean+n*sigma, respectively. The parameter specified by the user is n, the default value of n is 2.";
			}

			static const string NAME;

		protected:
			float default_sigma;
	};

	/**f(x) = x if x >= minval; f(x) = minval if x < minval
	 *@param minval Everything below this value is set to this value
	 */
	class ToMinvalProcessor:public Processor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}
		static Processor *NEW()
		{
			return new ToMinvalProcessor();
		}

		void process_inplace(EMData *image);

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("minval", EMObject::FLOAT, "Everything below this value is set to this value");
			d.put("newval", EMObject::FLOAT, "If set, values below minval will be set to newval instead of minval ");
			return d;
		}

		string get_desc() const
		{
			return "f(x) = x if x >= minval; f(x) = minval|newval if x < minval.";
		}

		static const string NAME;

	protected:

	};



	/**f(x) = x-minval if x >= minval; f(x) = 0 if x < minval
	 *@param minval the value that will be set to zero - all values below will also be set to zero. Values above get minval subtracted from them
	 */
	class CutToZeroProcessor:public RealPixelProcessor
	{
	  public:
		string get_name() const
		{
			return NAME;
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

		static const string NAME;

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
			return NAME;
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

		static const string NAME;

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
	 * All maps set below value are set to zero
	 * Useful in tomography when you want to toss complex components with low amplitides
	 *@author John Flanagan
	 *@date Oct 25th 2010
	 *@param value The Fourier amplitude threshold cutoff
	 */
	class BinarizeFourierProcessor:public Processor
		{
		  public:
			virtual string get_name() const
			{
				return NAME;
			}
			static Processor *NEW()
			{
				return new BinarizeFourierProcessor();
			}

			/**
			 * @exception ImageFormatException if the input image is not complex
			 * Note result is always in real-imaginary format
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
				return "f(k) = 0 + 0i if ||f(k)|| < value; f(k) = a + bi if ||f(k)|| >= value.";
			}

			static const string NAME;
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
			return NAME;
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

		static const string NAME;

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
			return NAME;
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
			d.put("shift", EMObject::FLOAT, "The amount to shift pixel values by after scaling");
			d.put("scale", EMObject::FLOAT, "The scaling factor to be applied to pixel values");
			return d;
		}

		string get_desc() const
		{
			return "linear transform processor: f(x) = x * scale + shift. This is equivalent to a regular contrast stretching operation";
		}

		static const string NAME;

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
			return NAME;
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

		static const string NAME;

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
				return NAME;
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

			static const string NAME;

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
			return NAME;
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

		static const string NAME;

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
			return NAME;
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

		static const string NAME;

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
			return NAME;
		}
		static Processor *NEW()
		{
			return new LogProcessor();
		}

		string get_desc() const
		{
			return "f(x) = log10(x) if x > 0; else f(x) = 0;";
		}

		static const string NAME;

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

	/**MaskAzProcessor masks out pixels within a specified cylindrical (or circular) arc.
	 *@param phi1 angle in degrees ccw from the x-axis. Starting angle to be set to 0.
	 *@param phi2 Ending angle to be set to 0
	 *@param dx Modify mask center by dx relative to the default center nx/2
	 *@param dy Modify mask center by dy relative to the default center ny/2
	 */
	class MaskAzProcessor:public Processor
	{
	  public:

		void process_inplace(EMData * image);
		
		static Processor *NEW()
		{
			return new MaskAzProcessor();
		}
		
		string get_name() const
		{
			return NAME;
		}
		
		static const string NAME;

		string get_desc() const
		{
			return "Masks out an angular arc in circular/cylindrical coordinates with a sharp edge.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;

			d.put("phicen", EMObject::FLOAT,"Angle in degrees ccw from the x-axis. Center of the region to NOT set to zero.");
			d.put("phirange", EMObject::FLOAT,"Angle in degrees. Region phicen+-phirange will not be zeroed");
			d.put("phitriangle", EMObject::BOOL, "If set mask will fall from 1 at phicen to 0 at phicen+-phirange");
			d.put("cx", EMObject::FLOAT,"Mask X center. Default nx/2");
			d.put("cy", EMObject::FLOAT,"Mask Y center. Default ny/2");
			d.put("inner_radius", EMObject::INT, "inner mask radius. optional. Default 0");
			d.put("outer_radius", EMObject::INT, "outer mask radius. optional. Default nx+ny. Negative value -> box radius + outer_radius +1");

			return d;
		}
	  protected:

		
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

			d.put("inner_radius", EMObject::INT, "inner mask radius. optional");
			d.put("outer_radius", EMObject::INT, "outer mask radius. Negative value -> box radius + outer_radius +1");

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

		virtual void process_dist_pixel(float *pixel, float dist) const = 0;		// note that this function gets the distance SQUARED !

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
			return NAME;
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

		static const string NAME;

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

	/**step cutoff to a user-given value in both inner and outer circles.
	 *@param value step cutoff to this value
	 */
	class MaskSoftProcessor:public CircularMaskProcessor
	{
	  public:
		MaskSoftProcessor():value(0)
		{
		}

		string get_name() const
		{
			return NAME;
		}
		static Processor *NEW()
		{
			return new MaskSoftProcessor();
		}

		void set_params(const Dict & new_params)
		{
			CircularMaskProcessor::set_params(new_params);
			value = params.set_default("value",0.0f);
			width = params.set_default("width",4.0f);
		}

		TypeDict get_param_types() const
		{
			TypeDict d = CircularMaskProcessor::get_param_types();
			d.put("value", EMObject::FLOAT, "cutoff to this value. default=0");
			d.put("width", EMObject::FLOAT, "1/e width of Gaussian falloff (both inner and outer) in pixels. default=4");
			return d;
		}

		string get_desc() const
		{
			return "Outer (optionally also inner) mask to a user-provided value with a soft Gaussian edge.";
		}

		static const string NAME;

	  protected:
		void process_dist_pixel(float *pixel, float dist) const
		{
			if (dist>=inner_radius_square && dist<=outer_radius_square) return;

			if (dist<inner_radius_square) *pixel=value+(*pixel-value)*exp(-pow((inner_radius-sqrt(dist))/width,2.0f));
			else *pixel=value+(*pixel-value)*exp(-pow((sqrt(dist)-outer_radius)/width,2.0f));
		}

		float value,width;
	};

	/**A step cutoff to the the mean value in a ring centered on the outer radius
	 *@param ring_width The width of the mask ring.
	 */
	class MaskEdgeMeanProcessor:public CircularMaskProcessor
	{							// 6
	  public:
		string get_name() const
		{
			return NAME;
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

		static const string NAME;

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
			return NAME;
		}
		static Processor *NEW()
		{
			return new MaskNoiseProcessor();
		}

		string get_desc() const
		{
			return "fills masked region";
		}

		static const string NAME;

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
			return NAME;
		}
		static Processor *NEW()
		{
			return new MaskGaussProcessor();
		}

		void set_params(const Dict & new_params)
		{
			CircularMaskProcessor::set_params(new_params);
			exponent = params["exponent"];
			if (exponent <= 0.0) {
				exponent = 2.0;
			}
		}

		TypeDict get_param_types() const
		{
			TypeDict d = CircularMaskProcessor::get_param_types();
			d.put("exponent", EMObject::FLOAT, "The exponent, f in e^-Bs^f. default 2.0, producing a Gaussian");
			return d;
		}

		string get_desc() const
		{
			return "a gaussian falloff to zero, radius is the 1/e of the width. If inner_radius>0, then \
outer radius specifies width of Gaussian starting at inner_radius rather than total radius.";
		}

		static const string NAME;

	  protected:
		float exponent;
		void process_dist_pixel(float *pixel, float dist) const
		{
			if (inner_radius_square>0) {
				if (dist>inner_radius_square) {
					if (exponent==2.0f) (*pixel) *= exp(-pow(sqrt(dist)-inner_radius,2.0f) / outer_radius_square);
					else (*pixel) *= exp(-pow(sqrt(dist)-inner_radius,exponent) / pow((float)outer_radius_square,exponent/2.0f));
				}
			}
			else {
				if (exponent==2.0f) (*pixel) *= exp(-dist / outer_radius_square);
				else (*pixel) *= exp(-pow(dist,exponent/2.0f) / pow((float)outer_radius_square,exponent/2.0f));
			}
		}
	};

	/**a gaussian falloff to zero, with anisotropic widths along x,y,z
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
			return NAME;
		}
		static Processor *NEW()
		{
			return new MaskGaussNonuniformProcessor();
		}

		string get_desc() const
		{
			return "A Gaussian falloff to zero. Anisotropic, specify inner radius for x,y,z and Gaussian falloff width. Falloff \
width is also anisotropic and relative to the radii, with 1 being equal to the radius on that axis.";
		}

		static const string NAME;

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
//			TypeDict d = CircularMaskProcessor::get_param_types();
			TypeDict d;
			d.put("gauss_width", EMObject::FLOAT, "Used to calculate the constant factor - gauss_width / (ny*ny)" );
//			d.put("ring_width", EMObject::INT, "The width of the mask ring.");
			return d;
		}

		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new MaskGaussInvProcessor();
		}

		string get_desc() const
		{
			return "f(x) = f(x) / exp(-radius*radius * gauss_width / (ny*ny))";
		}

		static const string NAME;

	  protected:
		void calc_locals(EMData *)
		{
			xc = Util::fast_floor(nx/2.0f) + dx;
			yc = Util::fast_floor(ny/2.0f) + dy;
			zc = Util::fast_floor(nz/2.0f) + dz;

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
			return NAME;
		}

		void process_inplace(EMData *image);

		static Processor *NEW()
		{
			return new LinearPyramidProcessor();
		}
		
		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("x0", EMObject::FLOAT);
			d.put("y0", EMObject::FLOAT);
			d.put("z0", EMObject::FLOAT);
			d.put("xwidth", EMObject::FLOAT);
			d.put("ywidth", EMObject::FLOAT);
			d.put("zwidth", EMObject::FLOAT);
			return d;
		}

		string get_desc() const
		{
			return "Multiplies image by a 'linear pyramid' in 1-3 dimensions. The origin and total base width of the pyramid can be specified. Default is centered with the total image size.";
		}

		static const string NAME;
	};


	/**overwrites input, f(x) = radius * radius
	 */
	class MakeRadiusSquaredProcessor:public CircularMaskProcessor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}
		static Processor *NEW()
		{
			return new MakeRadiusSquaredProcessor();
		}

		string get_desc() const
		{
			return "overwrites input, f(x) = radius * radius";
		}

		static const string NAME;

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
			return NAME;
		}
		static Processor *NEW()
		{
			return new MakeRadiusProcessor();
		}

		string get_desc() const
		{
			return "overwrites input, f(x) = radius;";
		}

		static const string NAME;

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
			return NAME;
		}
		static Processor *NEW()
		{
			return new ComplexNormPixel();
		}

		string get_desc() const
		{
			return "Each Fourier pixel will be normalized. ie - amp=1, phase=unmodified. Useful for performing phase-residual-like computations with dot products.";
		}

		static const string NAME;

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
			return NAME;
		}
		
		void process_inplace(EMData *image);
		
		static Processor *NEW()
		{
			return new LaplacianProcessor();
		}

		string get_desc() const
		{
			return "Discrete approximation to Laplacian. Edge enchancement, but works poorly in the presence of noise. Laplacian processor (x -> d^2/dx^2 + d^2/dy^2 + d^2/dz^2).";
		}

		static const string NAME;
		

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
			return NAME;
		}
		static Processor *NEW()
		{
			return new ZeroConstantProcessor();
		}

		string get_desc() const
		{
			return "Contraction of data, if any nearest neighbor is 0, value -> 0, generally used iteratively";
		}

		static const string NAME;

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
			return NAME;
		}
		static Processor *NEW()
		{
			return new BoxMedianProcessor();
		}

		string get_desc() const
		{
			return "A processor for noise reduction. pixel = median of values surrounding pixel.";
		}

		static const string NAME;

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
			return NAME;
		}
		static Processor *NEW()
		{
			return new BoxSigmaProcessor();
		}

		string get_desc() const
		{
			return "pixel = standard deviation of values surrounding pixel.";
		}

		static const string NAME;

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
			return NAME;
		}
		static Processor *NEW()
		{
			return new BoxMaxProcessor();
		}

		string get_desc() const
		{
			return "peak processor: pixel = max of values surrounding pixel.";
		}

		static const string NAME;

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
			return NAME;
		}
		static Processor *NEW()
		{
			return new MinusPeakProcessor();
		}

		string get_desc() const
		{
			return "peak processor: pixel = pixel - max of values surrounding pixel. This is a sort of positive peak-finding algorithm.";
		}

		static const string NAME;

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
			return NAME;
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
			d.put("npeaks", EMObject::INT, "The number of pixels adjacent to the pixel under consideration which may be higher and still be a valid peak. If 0, finds pure peaks");
			d.put("usemean", EMObject::BOOL, "Count all pixels with value higher than the mean of adjacent pixels as peaks. Overwrite npeaks.");
			return d;
		}

		string get_desc() const
		{
			return "Zeros all pixels with adjacent pixels >= the value being considered. That is, it leaves behind only local maxima.";
		}

		static const string NAME;

	  protected:
		void process_pixel(float *pixel, const float *data, int n) const
		{
			if (params["usemean"]){
				float mean=0;
				for (int i = 0; i < n; i++)
				{
					mean+=data[i];
				}

				if (*pixel < mean/float(n))
				{
					*pixel = 0;
				}
			}
			else{
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
			return NAME;
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

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
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
				return NAME;
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

			static const string NAME;

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
				return NAME;
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

			static const string NAME;

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
				return NAME;
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

			static const string NAME;

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
			return NAME;
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

		static const string NAME;

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
				return NAME;
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

			static const string NAME;

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
			return NAME;
		}
		static Processor *NEW()
		{
			return new GradientRemoverProcessor();
		}

		string get_desc() const
		{
			return "Gradient remover, does a rough plane fit to find linear gradients.";
		}

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
	};

	/** Make a curve or surface non-convex (planar or concave), iteratively.
	 *
	 *@author Steve Ludtke
	 *@date 2011/08/11
	 */
	class NonConvexProcessor:public Processor
    {
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return NAME;
		}
		static Processor *NEW()
		{
			return new NonConvexProcessor();
		}

		string get_desc() const
		{
			return "Makes a curve or plane monotonically decreasing and non-convex. Useful in generating background curves from power spectra. Anchored at edges and (in 2d) at the center. If local value > mean(surrounding values) => mean(surrounding values).";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
/*			d.put("mask", EMObject::EMDATA, "mask object: nonzero pixel positions will be used to fit plane. default = 0");
			d.put("changeZero", EMObject::INT, "if zero pixels are modified when removing gradient. default = 0");
			d.put("planeParam", EMObject::FLOATARRAY, "fitted plane parameters output");*/
			return d;
		}

		static const string NAME;
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
				return NAME;
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

			static const string NAME;
	};


	/**Ramp processor -- Fits a least-squares plane to the picture, and subtracts the plane from the picture.  A wedge-shaped overall density profile can thus be removed from the picture.
	 */
	class RampProcessor:public Processor
    {
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return NAME;
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

		static const string NAME;
	};

	/**Tries to fix images scanned on the zeiss for poor ccd normalization.
	 */
	class VerticalStripeProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new VerticalStripeProcessor();
		}

		string get_desc() const
		{
			return "Tries to fix images scanned on the zeiss for poor ccd normalization.";
		}

		static const string NAME;
	};

	/**This will replace the image with a full-circle 2D fft amplitude rendering.
	 */
	class RealToFFTProcessor:public Processor
	{
		public:
		void process_inplace(EMData *image);

		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new RealToFFTProcessor();
		}

		string get_desc() const
		{
			return "This will replace the image with a full-circle 2D fft amplitude rendering. Note that this renders amplitude, when intensity is more common.";
		}

		static const string NAME;
	};

	/** Fill missing wedge with information from another image
	 */
	class WedgeFillProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return NAME;
		}
		static Processor *NEW()
		{
			return new WedgeFillProcessor();
		}

		string get_desc() const
		{
			return "Identifies missing wedge voxels and fills them with data extracted from another image";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("fillsource", EMObject::EMDATA, "The image from which to draw the missing values");
			return d;
		}

		static const string NAME;
	};

	/**Fill zeroes at edges with nearest horizontal/vertical value.
	 */
	class SigmaZeroEdgeProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return NAME;
		}
		static Processor *NEW()
		{
			return new SigmaZeroEdgeProcessor();
		}

		string get_desc() const
		{
			return "Fill zeroes at edges with nearest horizontal/vertical value.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("nonzero", EMObject::BOOL, "If set, will look for constant non-zero values to fill");
			return d;
		}

		static const string NAME;
	};

	/** This processor will try and remove outliers (and optionally exactly zero values), replacing any identified values with the local mean value
	 * @param sigma outliers are defined as mean+-x*sigma where x is the specified value
	 * @param fix_zero if set, any values that are exactly zero will be treated as outliers)
	 */
	class OutlierProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new OutlierProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("fix_zero", EMObject::BOOL, "If set, any pixels that are exactly zero are considered to be outliers, default=false");
			d.put("sigma", EMObject::FLOAT, "outliers are defined as mean+-x*sigma where x is the specified value, default=3.0");
			return d;
		}

		string get_desc() const
		{
			return "Identifies any pixels that are outliers and adjusts them to be the average of any nearby non-outlier pixels. Operates iteratively when required, so large 'outlier' areas can be corrected.";
		}

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
	};

	/**Fill zeroes at edges with nearest horizontal/vertical value damped towards Mean2.
	 */
	class MeanZeroEdgeProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new MeanZeroEdgeProcessor();
		}

		string get_desc() const
		{
			return "Fill zeroes at edges with nearest horizontal/vertical value damped towards Mean2.";
		}

		static const string NAME;
	};


	/**Average along Y and replace with average
	 */
	class AverageXProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new AverageXProcessor();
		}

		string get_desc() const
		{
			return "Average along Y and replace with average";
		}

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
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
			return NAME;
		}

		static Processor *NEW()
		{
			return new NormalizeUnitProcessor();
		}

		string get_desc() const
		{
			return "Normalize an image so its vector length is 1.0.";
		}

		static const string NAME;

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
			return NAME;
		}

		static Processor *NEW()
		{
			return new NormalizeUnitSumProcessor();
		}

		string get_desc() const
		{
			return "Normalize an image so its elements sum to 1.0 (fails if mean=0)";
		}

		static const string NAME;

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
			return NAME;
		}

		static Processor *NEW()
		{
			return new NormalizeStdProcessor();
		}

		string get_desc() const
		{
			return "do a standard normalization on an image.";
		}

		static const string NAME;

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
			return NAME;
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

		static const string NAME;

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
				return NAME;
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

			static const string NAME;
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
				return NAME;
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
				d.put("verbose", EMObject::INT,"If set will give details about the normalization");
				return d;
			}

			void process_inplace(EMData * image);

			static const string NAME;
	};


	/**normalizes an image, mean value equals to edge mean.
	 */
	class NormalizeEdgeMeanProcessor:public NormalizeProcessor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new NormalizeEdgeMeanProcessor();
		}

		string get_desc() const
		{
			return "normalizes an image, mean value equals to edge mean.";
		}

		static const string NAME;

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
			return NAME;
		}

		static Processor *NEW()
		{
			return new NormalizeCircleMeanProcessor();
		}

		string get_desc() const
		{
			return "normalizes an image, mean value equals to mean of 2 pixel circular radius or of the circular border if no radius is set.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("radius", EMObject::FLOAT,"Radius of 2 pixel circular border");
			return d;
		}

		static const string NAME;

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
			return NAME;
		}

		static Processor *NEW()
		{
			return new NormalizeLREdgeMeanProcessor();
		}

		string get_desc() const
		{
			return "normalizes an image, uses 2 pixels on left and right edge";
		}

		static const string NAME;

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
			return NAME;
		}

		static Processor *NEW()
		{
			return new NormalizeMaxMinProcessor();
		}

		string get_desc() const
		{
			return "normalizes an image. mean -> (maxval-minval)/2; std dev = (maxval+minval)/2;";
		}

		static const string NAME;

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
			return NAME;
		}

		static Processor *NEW()
		{
			return new NormalizeRowProcessor();
		}

		string get_desc() const
		{
			return "normalizes each row in the image individually";
		}

		static const string NAME;

		void process_inplace(EMData * image);
	};

	/**Sorry for the pun. This processor will take a second image and try to filter/scale it to optimally subtract it
	from the original image. The idea here is that if you have an image with noise plus a linear-filter modified projection,
	that a good measure of the similarity of the image to the projection would be to try and remove the projection from
	the image as optimally as possible, then compute the standard deviation of what's left.

	Now you might say that if the total energy in the noisy image is normalized then this should be equivalent to just
	integrating the FSC, which is what we use to do the optimal subtraction in the first place. This would be true, but
	this "optimal subtraction" has other purposes as well, such as the e2extractsubparticles program.
	 * @param ref Reference image to subtract
	 * @param return_radial Will return the radial filter function applied to ref as filter_curve
	 */
	class SubtractOptProcessor:public Processor
	{
	  public:
		virtual void process_inplace(EMData *image);
		virtual EMData* process(const EMData * const image);

		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new SubtractOptProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("ref", EMObject::EMDATA, "Reference image to subtract");
			d.put("actual", EMObject::EMDATA, "If specified, ref is used for normalization, but actual is subtracted.");
			d.put("low_cutoff_frequency", EMObject::FLOAT, "Absolute [0,0.5] low cut-off frequency.");
			d.put("high_cutoff_frequency", EMObject::FLOAT, "Absolute [0,0.5] high cut-off frequency.");
			d.put("ctfweight",EMObject::BOOL, "Filter the image by CTF before subtraction");
			d.put("return_fft",EMObject::BOOL, "Skips the final IFT, and returns the FFT of the subtracted image");
			d.put("return_subim", EMObject::BOOL, "Instead of returning the image after subtraction, returns the filtered image which would have been subtracted from the image.");
			d.put("return_radial", EMObject::BOOL, "Return the radial filter function as an attribute (filter_curve)");
			d.put("return_presigma", EMObject::BOOL, "Return the sigma of the pre-subtracted image in real-space with the specified filter applied as sigma_presub. This is an expensive option.");
			return d;
		}

		string get_desc() const
		{
			return "This will filter/scale 'ref' optimally and subtract it from image using ring dot products in Fourier space for normalization. Cutoff frequencies apply a bandpass tophat filter to the output.";
		}

		static const string NAME;
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
			return NAME;
		}

		static Processor *NEW()
		{
			return new NormalizeToLeastSquareProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("to", EMObject::EMDATA, "reference image normalize to");
			d.put("ignore_zero", EMObject::BOOL, "If set, ignores any pixels which are exactly zero in either image. Defaut = True.");
			d.put("ignore_lowsig", EMObject::FLOAT, "If >0, then any pixels closer to the mean than val*sigma in either image excluded");
			d.put("low_threshold", EMObject::FLOAT, "only take into account the reference image's pixel value between high and low threshold (zero is always ignored)");
			d.put("high_threshold", EMObject::FLOAT, "only take into account the reference image's pixel value between high and low threshold (zero is always ignored)");
			return d;
		}

		string get_desc() const
		{
			return "use least square method to normalize";
		}

		static const string NAME;
	};

	/**makes image circularly symmetric.
	 */
	class RotationalAverageProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new RotationalAverageProcessor();
		}

		string get_desc() const
		{
			return "Makes image circularly/spherically symmetric.";
		}

		static const string NAME;
	};

	/**subtracts circularly symmetric part of an image.
	 */
	class RotationalSubstractProcessor:public Processor
	{
	  public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new RotationalSubstractProcessor();
		}

		virtual string get_desc() const
		{
			return "subtracts circularly/spherically symmetric part of an image.";
		}

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
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
			return NAME;
		}

		static Processor *NEW()
		{
			return new FlipProcessor();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("axis", EMObject::STRING, "'x', 'y', or 'z' axis.");
			return d;
		}

		virtual string get_desc() const
		{
			return "Mirrors an image along the specified axis, preserving the center. This will introduce a plane of 0's for even box sizes. Use 'xform.mirror' processor to avoid the zero plane, but not preserve the center.";
		}

		static const string NAME;
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
			return NAME;
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

		static const string NAME;

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
			return NAME;
		}

		static Processor *NEW()
		{
			return new AddSigmaNoiseProcessor();
		}

		virtual string get_desc() const
		{
			return "add sigma noise.";
		}

		static const string NAME;

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
			return NAME;
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

		static const string NAME;
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
				return NAME;
			}

			static Processor *NEW()
			{
				return new FourierToCornerProcessor();
			}

			virtual string get_desc() const
			{
				return "Undoes the xform.fourierorigin.tocenter processor";
			}

			static const string NAME;
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
				return NAME;
			}

			static Processor *NEW()
			{
				return new FourierToCenterProcessor();
			}

			virtual string get_desc() const
			{
				return "Translates the origin in Fourier space from the corner to the center in y and z - works in 2D and 3D";
			}

			static const string NAME;
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
				return NAME;
			}

			static Processor *NEW()
			{
				return new PhaseToCenterProcessor();
			}

			virtual string get_desc() const
			{
				return "Undoes the effect of the xform.phaseorigin.tocorner processor";
			}

			static const string NAME;
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
				return NAME;
			}

			static Processor *NEW()
			{
				return new PhaseToCornerProcessor();
			}

			virtual string get_desc() const
			{
				return "Translates a centered image to the corner in a forward fashion";
			}

			static const string NAME;
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
			return NAME;
		}

		static Processor *NEW()
		{
			return new AutoMask2DProcessor();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("radius", EMObject::INT,"Pixel radius of a ball which is used to seed the flood filling operation. ");
			d.put("nmaxseed",EMObject::INT,"Use the n highest valued pixels in the map as a seed. Alternative to radius. Useful for viruses.");
			d.put("threshold", EMObject::FLOAT, "An isosurface threshold that suitably encases the mass.");
			d.put("sigma", EMObject::FLOAT, "Alternative to threshold based on mean + x*sigma");
			d.put("nshells", EMObject::INT, "The number of dilation operations");
			d.put("nshellsgauss", EMObject::INT, "number of Gaussian pixels to expand, following the dilation operations");
			d.put("return_mask", EMObject::BOOL, "If true the result of the operation will produce the mask, not the masked volume.");
			d.put("verbose", EMObject::INT, "How verbose to be (stdout)");
			return d;
		}

		virtual string get_desc() const
		{
			return "2D version of mask.auto3d";
		}

		static const string NAME;
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
				return NAME;
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

			static const string NAME;
	};

	/** A "dust removal" filter which will remove above threshold densities smaller than a given size
	 * @param voxels
	 * @param threshold
	 */
	class AutoMaskDustProcessor:public Processor
	{
	  public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new AutoMaskDustProcessor();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("threshold", EMObject::FLOAT,"Only considers densities above the threshold");
			d.put("voxels", EMObject::INT,"If a connected mass is smaller than this many voxels it is removed");
			d.put("verbose", EMObject::INT, "Level of verbosity, 0 default. 1 will print each non-excluded zone");
			return d;
		}

		virtual string get_desc() const
		{
			return "A dust removal filter which will remove above threshold densities smaller than a given size";
		}

		static const string NAME;

		protected:
		EMData *mask;
		EMData *image;
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
			return NAME;
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

		static const string NAME;
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
			return NAME;
		}

		static Processor *NEW()
		{
			return new AutoMask3D2Processor();
		}

		virtual string get_desc() const
		{
			return "This will mask a 3-D volume using a 'flood filling' approach. It begins with a seed generated either as a sphere with \
specified 'radius' or with the 'nmaxseed' highest values. It then includes any mass connected to the seed with value higher than 'threshold'.\
Next, the mask is expanded by 'nshells'+'nshellsgauss'/2 voxels. Finally a gaussian low-pass filter is applied with a width of 'nshellsgauss'.";
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("radius", EMObject::INT,"Pixel radius of a ball which is used to seed the flood filling operation. ");
			d.put("nmaxseed",EMObject::INT,"Use the n highest valued pixels in the map as a seed. Alternative to radius. Useful for viruses.");
			d.put("threshold", EMObject::FLOAT, "An isosurface threshold that suitably encases the mass.");
			d.put("sigma", EMObject::FLOAT, "Alternative to threshold based on mean + x*sigma");
			d.put("nshells", EMObject::INT, "Number of 1-voxel shells to expand the mask by.");
			d.put("nshellsgauss", EMObject::INT, "Width in voxels of a Gaussian decay at the edge of the mask.");
			d.put("return_mask", EMObject::BOOL, "If true the result of the operation will produce the mask, not the masked volume.");
			d.put("verbose", EMObject::INT, "How verbose to be (stdout)");
			return d;
		}

		static const string NAME;
	};

	/** This expands a multilevel mask volume so inter-mask boundaries are preserved
	 * @param nshells   number of shells to add
	*/
	class IterMultiMaskProcessor:public Processor
	{
	  public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return NAME;
		}

		virtual string get_desc() const
		{
			return "A multilevel mask has an integer value at each pixel location. -1 indicates unmasked regions. 0-n-1 are individual masks. Expands the masked regions into unmasked areas by nshells.";
		}

		static Processor *NEW()
		{
			return new IterMultiMaskProcessor();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("nshells", EMObject::INT, "number of shells to add");
			return d;
		}

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
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
				return NAME;
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

			static const string NAME;
	};
	
	/**ToCenterProcessor centers image, ignores old dx, dy.
	 * @ingroup tested3c
	 */
	class ToCenterProcessor:public Processor
	{
	  public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new ToCenterProcessor();
		}

		virtual string get_desc() const
		{
			return "Centers similar to the way a human would, by identifying the shape of the object and centering its sillouette. May be inaccurate if sillouette cannot be clearly identified. ";
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
//			d.put("int_shift_only", EMObject::INT, "set to 1 only shift by integer, no interpolation");
//			d.put("threshold", EMObject::FLOAT, "Only values larger than the threshold are included in the center of mass computation. Default is 0.");
//			d.put("positive", EMObject::INT, "uses only densities >0 for the calculatton");
			return d;
		}

		static const string NAME;
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
			return NAME;
		}

		static Processor *NEW()
		{
			return new ToMassCenterProcessor();
		}

		virtual string get_desc() const
		{
			return "ToMassCenterProcessor centers image at center of mass, with a threshold. Only values higher than the threshold are considered.";
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("int_shift_only", EMObject::INT, "set to 1 only shift by integer, no interpolation");
			d.put("threshold", EMObject::FLOAT, "Only values larger than the threshold are included in the center of mass computation. Default is 0.");
//			d.put("positive", EMObject::INT, "uses only densities >0 for the calculatton");
			return d;
		}

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
	};


	/** This processor will apply a Wiener filter to a volume based on a provided FSC curve. The assumption is that the FSC curve
	 * represents a "gold standard" FSC between two 1/2 sets, and that the filter is being applied to the combined average.
	 * @param snrmult This multiplier is applied to the computed SNR before Wiener filtration. This permits the filter to be applied to 1/2 images, etc. Default=1.0
	 * @param fscfile File containing the FSC curve to use
	 * @param sscale This rescales the S axis to produce empirical under/overfiltration. sscale=1.1 for example will extend the resolution (underfilter) by 10%. Default=1.0
	 */
	class FSCFourierProcessor:public Processor
	{
	  public:
		virtual EMData* process(EMData const *image);
		virtual void process_inplace(EMData *image);

		virtual string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new FSCFourierProcessor();
		}

		virtual string get_desc() const
		{
			return "This processor will apply a Wiener filter to a volume based on a provided FSC curve. The assumption is that the FSC curve represents \
a gold standard FSC between two 1/2 sets, and that the filter is being applied to the combined average. Hence the default fscmult of 2, \
since the SSNR is being computed as FSC/(1-FSC). Ie - the SSNR of the combined halves is twice as high.";
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("snrmult", EMObject::FLOAT, "This multiplier is applied to the computed SNR before Wiener filtration. This permits the filter to be applied to 1/2 images, etc. Default=2.0");
			d.put("sscale", EMObject::FLOAT, "This rescales the S axis to produce empirical under/overfiltration. sscale=1.1 for example will extend the resolution (underfilter) by 10%. Default=1.0");
			d.put("maxfreq", EMObject::FLOAT, "This acts as a high resolution limit to prevent FSC artifacts from iteratively reinforcing themselves. Above this spatial frequency, the FSC is forced to decrease monotonically. Default=1.0");
			d.put("fscfile", EMObject::STRING, "filename of a file containing the FSC curve to use for the SNR computation");
			return d;
		}

		static const string NAME;
	};

	/** Processor the images by the estimated SNR in each image.if parameter 'wiener' is 1, then wiener processor the images using the estimated SNR with CTF amplitude correction.
	 * @param defocus mean defocus in microns
	 * @param voltage microscope voltage in Kv
	 * @param ac amplitude contrast in %
	 */
	class CTFCorrProcessor:public Processor
	{
	  public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new CTFCorrProcessor();
		}

		virtual string get_desc() const
		{
			return "One of the strongest visual impacts of CTF on a structure is the low resolution high-pass filter effect caused by \
phase contrast. This Processor performs a simple linear filter to roughly correct for this. This is not a substitution or replacement \
for the full CTF correction routine available for single particle work in EMAN, but if you are in a situation where accurate CTF \
correction is not possible, this will allow you to approximate the correction to relieve some of the visual artifacts.";
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("defocus", EMObject::FLOAT, "Mean defocus to correct for in microns");
			d.put("ac", EMObject::FLOAT, "Amplitude contrast in % (default 10%)");
			d.put("voltage", EMObject::FLOAT, "Microscope Voltage in Kv (default 300)");
			d.put("apix", EMObject::FLOAT, "A/pix (default value from image header)");
			return d;
		}

		static const string NAME;
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
			return NAME;
		}

		static Processor *NEW()
		{
			return new SNRProcessor();
		}

		virtual string get_desc() const
		{
			return "Process the images by the estimated SNR in each image.if parameter 'wiener' is 1, then wiener processor the images using the estimated SNR with CTF amplitude correction.";
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("wiener", EMObject::INT, "if set to 1,  then use wiener processor to process the images using the estimated SNR with CTF amplitude correction");
			d.put("snrfile", EMObject::STRING, "structure factor file name");
			return d;
		}

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
	};

	/** This processor will remove specific bad lines from CCD images, generally due to faulty lines/rows in the detector.
	 * Specify only one of xloc or yloc
	 *@param xloc x location of a bad vertical line
	 *@param yloc y location of a bad horizontal line
	 */
	class BadLineXYProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new BadLineXYProcessor();
		}

		virtual string get_desc() const
		{
			return "This processor will remove localized 'striping' along the x/y axes, caused by issues with CCD/CMOS readout. In theory this should be done by dark/gain correction, but in many cases, there are residual effects that this will help eliminate. This can produce high-pass filter-like effects, so generally large length values are suggested. Integration covers +-xlen/ylen. Y and X axes are corrected sequentially, not simultaneously, Y first";
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("xloc", EMObject::INT, "X coordinate of a bad vertical line");
			d.put("yloc", EMObject::INT, "Y coordinate of a bad vertical line. Specify only one of xloc or yloc");
			return d;
		}

		static const string NAME;
	};

	
	/** This processor will remove localized 'striping' along the x/y axes, caused by issues with CCD/CMOS readout. In theory this should be done by dark/gain correction, but in many cases, there are residual effects that this will help eliminate.
	 *@param threshold an isosurface threshold at which all desired features are visible
	 *@param radius a normalization size similar to an lp= value
	 *@param apix Angstrom per pixel ratio
	 */
	class StripeXYProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new StripeXYProcessor();
		}

		virtual string get_desc() const
		{
			return "This processor will remove localized 'striping' along the x/y axes, caused by issues with CCD/CMOS readout. In theory this should be done by dark/gain correction, but in many cases, there are residual effects that this will help eliminate. This can produce high-pass filter-like effects, so generally large length values are suggested. Integration covers +-xlen/ylen. Y and X axes are corrected sequentially, not simultaneously, Y first";
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("xlen", EMObject::INT, "Integration 1/2 length on x axis in pixels. Default=10");
			d.put("ylen", EMObject::INT, "Integration 1/2 length on y axis in pixels. Default=10");
			return d;
		}

		static const string NAME;
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
			return NAME;
		}

		static Processor *NEW()
		{
			return new LocalNormProcessor();
		}

		virtual string get_desc() const
		{
			return "This processor attempts to perform a 'local normalization' so low density and high density features will be on a more even playing field in an isosurface display. threshold is an isosurface threshold at which all desired features are visible, radius is a feature size over which to equalize.";
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("threshold", EMObject::FLOAT, "Only values above the threshold will be used to compute the normalization. Generally a good isosurface value.");
			d.put("radius", EMObject::FLOAT, "Fourier filter radius expressed in pixels in Fourier space. cutoff_pixels in filter.lowpass.gauss");
			d.put("apix", EMObject::FLOAT, "Angstroms per pixel");
			return d;
		}

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
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
			return NAME;
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

		static const string NAME;

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
			return NAME;
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
		virtual void process_inplace(EMData*) {
			throw InvalidCallException("The directional sum processor does not work inplace");
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("axis", EMObject::STRING,"The direction of the sum, either x,y or z. Returned axes are xy, xz or zy.");
			d.put("first", EMObject::INT,"The first position along the speficied axis to use in the sum. Neg val -> nx/y/z+first (default=0)");
			d.put("last", EMObject::INT,"The last position along the speficied axis to use in the sum. Neg val -> nx/y/z+last (default=-1)");
			return d;
		}

		string get_desc() const
		{
			return "Calculates the projection of the image along one of the axial directions, either x, y or z";
		}

		static const string NAME;
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
		virtual EMData* process(const EMData* const image);
		virtual void process_inplace(EMData*);

		virtual string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new WatershedProcessor();
		}

		virtual string get_desc() const
		{
			return "Watershed segmentation. Warning: uses up to 2.5x the map size in RAM. This will segment all voxels above threshold except for a 1-voxel wide border on all edges.";
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("nseg", EMObject::INT, "Number of segments to (attempt) to produce. The actual number may be fewer. (default=12)" );
			d.put("thr",EMObject::FLOAT,"Isosurface threshold value. Pixels below this value will not be segmented. All voxels above this value will be segmented. (default=0.5)");
			d.put("segbymerge", EMObject::INT, "If set, will achieve the specified number of segments by progressively merging the most connected segments. Can produce very different results." );
			d.put("verbose", EMObject::INT, "If set, will print console output while running" );
			return d;
		}

		static const string NAME;

//	  private:

	};

	/** Operates on two images, returning an image containing the maximum/minimum/multiplied pixel (etc, you choose) at each location
	 * The actual operation depends on what template argument you use. Currently the MaxPixelOperator and MinPixelOperator
	 * are the only ones utilized - see processor.cpp where the Processor Factory constructor is defined to get an idea of how to add another one
	 * Initially added at the request of Dr Matthew Baker
	 * NOTE: binary is meant in the sense of the standard algorithms, NOT in the sense of a black and white image
	 * @param with the other image that will be used generate the image with the maximum (or minimum, etc) pixel values
	 * @author David Woolford
	 * @date June 2009
	 */
	template<class Type>
	class BinaryOperateProcessor : public Processor{
		public:
			/**
			* @exception InvalidParameterException if with is not specified
			* @exception ImageDimensionException if image dimensions do not match
			*/
			virtual void process_inplace(EMData * image) {
				if ( ! params.has_key("with") ) throw InvalidParameterException("You must supply the \"with\" parameter");
				EMData* with = params["with"];

				if ( with->get_xsize() != image->get_xsize() || with->get_ysize() != image->get_ysize() || with->get_zsize() != image->get_zsize() )
					throw ImageDimensionException("The images you are operating on do not have the same dimensions");

				float* image_data = image->get_data();
				float* with_data = with->get_data();

				std::transform(image_data,image_data+image->get_size(),with_data,image_data,Type::binary_operate);
				image->update();
			}

			virtual string get_name() const
			{
				return op.get_name();
			}

			virtual string get_desc() const
			{
				return op.get_desc();
			}

			static Processor *NEW()
			{
				return new BinaryOperateProcessor<Type>();
			}

			virtual TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("with", EMObject::EMDATA,"The second image");
				return d;
			}

			static const string NAME;
		private:
			Type op;
	};

	class MaxPixelOperator {
		public:
		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Compares pixels in two images, returning an image with the maximum pixel value in each pixel location";
		}

		static float binary_operate(const float& left, const float& right) {
			if (left > right) return left;
			return right;
		}

		static const string NAME;
	};

	class MinPixelOperator {
		public:
			string get_name() const
			{
				return NAME;
			}

			string get_desc() const
			{
				return "Compares pixels in two images, returning an image with the minimum pixel value in each pixel location";
			}

			static float binary_operate(const float& left, const float& right) {
				if (left < right) return left;
				return right;
			}

			static const string NAME;
	};

	/**Sets the structure factor To match a second provided image/volume
	 *@param to EMData object to match to. Make sure apix values are set properly
	 */
	class MatchSFProcessor:public FourierAnlProcessor
	{
	  public:

		virtual string get_name() const
		{
			return NAME;
		}

		virtual string get_desc() const
		{
			return "Filters the image so its 1-D power spectrum matches a second image. Optionally can incorporate a dot product to better match noise.";
		}

		static Processor *NEW()
		{
			return new MatchSFProcessor();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("to", EMObject::EMDATA, "The image to match with. Make sure apix values are correct.");
			d.put("bydot", EMObject::BOOL, "Rather than matching the intensity profile, uses the complex dot product as a function of resolution to match only the portion that agrees.");
			d.put("keephires", EMObject::BOOL, "If the reference being matched is heavily filtered, total information loss may occur at some resolutions. This insures that some information is kept at all resolutions.");
			d.put("return_radial", EMObject::BOOL, "Return the radial filter function as an attribute (filter_curve)");
			d.put("interpolate", EMObject::BOOL, "Whether or not to interpolate the radial scaling function. Default=true");
			return d;
		}

		static const string NAME;

	  protected:
		void create_radial_func(vector < float >&radial_mask, EMData *image) const;
	};


	/**Sets the structure factor based on a 1D s/intensity curve as an XYData object.
	 *@param strucfac XYData object with the curve to be imposed as intensity as a function of s
	 *@param apix A/pix value. Overrides and replaces apix_x/y/z in image
	 */
	class SetSFProcessor:public FourierAnlProcessor
	{
	  public:

		virtual string get_name() const
		{
			return NAME;
		}

		virtual string get_desc() const
		{
			return "Filters the image so its 1-D power spectrum matches a supplied S,Y curve. If the S axis does not extend to Nyquist, only a uniform scaling will be applied beyond the end of the supplied curve. ";
		}

		static Processor *NEW()
		{
			return new SetSFProcessor();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("strucfac", EMObject::XYDATA, "An XYData object contaning the curve to be imposed as a function of S");
			d.put("apix", EMObject::FLOAT, " Override A/pix in the image header (changes x,y and z)");
			return d;
		}

		static const string NAME;

	  protected:
		void create_radial_func(vector < float >&radial_mask, EMData *image) const;
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
			return NAME;
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

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
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
				return NAME;
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

			static const string NAME;
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
				return NAME;
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

			static const string NAME;
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
				return NAME;
			}

			virtual string get_desc() const
			{
				return "Make an image consisting various objects, useful for tomographic testing";
			}

			static Processor * NEW()
			{
				return new TestTomoImage();
			}

			static const string NAME;

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
				return NAME;
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

			static const string NAME;
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
				return NAME;
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

			static const string NAME;
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
			return NAME;
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

		static const string NAME;
	};

	/**Replace a source image with a lumpy S-curve used for alignment testing
	 */
	class TestImageScurve : public TestImageProcessor
	{
	public:
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return NAME;
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

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
	};


	/**Replace a source image as a sine wave in specified wave length
	 *@param wavelength wavelength in equation sin(x*2*PI/wavelength - phase*180/PI)
	 *@param axis (optional) specify a major axis for asymmetric features, default x axis
	 *@param phase (optional) the phase in radians
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
			return NAME;
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
			d.put("phase", EMObject::FLOAT, "(optional) the phase in radians");
			d.put("az", EMObject::FLOAT, "(optional) angle in degree. for 2D image, this is the rotated angle of the image, \
												in 3D image, it's az for euler angle. default is zero");
			d.put("alt", EMObject::FLOAT, "(optional) angle in degree. only in 3D case, alt for euler angle, default is zero");
			d.put("phi", EMObject::FLOAT, "(optional) angle in degree. only in 3D case, phi for euler angle, default is zero");
			return d;
		}

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
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
			return NAME;
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
			d.put("a", EMObject::FLOAT, "equatorial radius along x axes (major semiaxes)");
			d.put("b", EMObject::FLOAT, "equatorial radius along y axes (minor semiaxes)");
			d.put("c", EMObject::FLOAT, "polar radius for ellipsoid (x^2/a^2+y^2/b^2+z^2/c^2=1)");
			d.put("transform", EMObject::TRANSFORM, "Optionally transform the ellipse");
			d.put("fill", EMObject::FLOAT, "value you want to fill in ellipse, default to 1.0");
			return d;
		}

		static const string NAME;
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
				return NAME;
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

			static const string NAME;
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
			return NAME;
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

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
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
			return NAME;
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

		static const string NAME;
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
			return NAME;
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

		static const string NAME;

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
				return NAME;
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

			static const string NAME;

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
			return NAME;
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

		static const string NAME;
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
			return NAME;
		}

		static Processor *NEW()
		{
			return new RadialProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("table", EMObject::FLOATARRAY, "Radial array of floats, 1 float/pixel");
			return d;
		}

		string get_desc() const
		{
			return "Multiply a real-space image by a radial function. 1 value / pixel, extending to corner. Missing values -> 0.";
		}

		static const string NAME;
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
				return NAME;
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

			static const string NAME;

		protected:
			int default_bins;
	};

	class ModelHelixProcessor : public Processor
	{
	  protected:
		float radprofile(float r, int type);
	};

	class ModelEMCylinderProcessor : public ModelHelixProcessor
	{
	  public:
		void process_inplace(EMData * in);

		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new ModelEMCylinderProcessor();
		}

		string get_desc() const
		{
			return "Adds a cylinder with a radial density profile similar to that of an alpha helix.";
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("type", EMObject::INT, "Radial profile of density method, defaults to 2: 0 = pure Gaussian falloff; 1 = Gaussian falloff + dip, so mean is zero; 2 = polynomial fitting of real helix density");
			d.put("length", EMObject::FLOAT, "cylinder length in angstroms, defaults to 3 turns (16.2 Angstroms)");
			d.put("x0", EMObject::INT, "x coordinate in pixels for the midpoint of the cylinder's axis, defaults to center of map");
			d.put("y0", EMObject::INT, "y coordinate in pixels for the midpoint of the cylinder's axis, defaults to center of map");
			d.put("z0", EMObject::INT, "z coordinate in pixels for the midpoint of the cylinder's axis, defaults to center of map");
			//TODO: Check with Matt Baker about description strings
			return d;
		}

		static const string NAME;
	};

	class ApplyPolynomialProfileToHelix : public ModelHelixProcessor
	{
	public:
		void process_inplace(EMData * in);

		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new ApplyPolynomialProfileToHelix();
		}

		string get_desc() const
		{
			return "Finds the CM of each z-axis slice and applies a polynomial radial profile about it.";
		}
		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("length", EMObject::FLOAT, "Helix length in angstroms.");
			d.put("z0", EMObject::INT, "z coordinate in pixels for the midpoint of the cylinder's axis, defaults to center of map");
			return d;
		}

		static const string NAME;
	};

	class BinarySkeletonizerProcessor : public Processor
	{
	public:
		virtual EMData* process(EMData * image);
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return NAME;
//			return "gorgon.binary_skel";
		}
		static Processor *NEW()
		{
			return new BinarySkeletonizerProcessor();
		}
		string get_desc() const
		{
			return "Creates a skeleton of the 3D image by considering whether density is above or below a threshold value.";
		}
		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("threshold", EMObject::FLOAT, "Threshold value.");
			d.put("min_curve_width", EMObject::INT, "Minimum curve width.");
			d.put("min_surface_width", EMObject::INT, "Minimum surface width.");
			d.put("mark_surfaces", EMObject::BOOL, "Mark surfaces with a value of 2.0f, whereas curves are 1.0f.");
			return d;
		}
		static const string NAME;
	};

	class ConvolutionKernelProcessor : public Processor
	{
	public:
		virtual EMData* process(const EMData* const image);
		virtual void process_inplace(EMData * image);

		virtual string get_name() const
		{
			return NAME;
		}
		static Processor *NEW()
		{
			return new ConvolutionKernelProcessor();
		}
		string get_desc() const
		{
			return "Filters an image with a convolution kernel in real space.";
		}
		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("kernel", EMObject::FLOATARRAY, "the convolution kernel");
			d.put("selem", EMObject::EMDATA, "the structuring element");
			return d;
		}
		static const string NAME;
	};

	class RotateInFSProcessor : public Processor
		{
		public:
			//virtual EMData* process(const EMData* const image);
			virtual void process_inplace(EMData * image);
			virtual EMData* process(const EMData* const image);

			virtual string get_name() const
			{
				return NAME;
			}
			static Processor *NEW()
			{
				return new RotateInFSProcessor( );
			}
			string get_desc() const
			{
				return "Rotates a Fourier object using a kernel.";
			}
			virtual TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("transform", EMObject::TRANSFORM, "transform");
				d.put("interpCutoff", EMObject::FLOAT, "cutoff for interpolation");
//				d.put("offset", EMObject::FLOAT, "offset for FT centering");
//				d.put("angle", EMObject::FLOAT, "angle");
				return d;
			}
			static const string NAME;
		};
	
	/*
	 * The base class for morphological processors. 
	 */
	class MorphologicalProcessor:public Processor
	{
	  public:
		MorphologicalProcessor(): value(0), maxval(1), mean(0), sigma(0)
		{
		}
		void process_inplace(EMData * image);
		
		virtual void set_params(const Dict & new_params)
		{
			params = new_params;
			if (params.size() == 1) {
				vector < EMObject > dict_values = params.values();
				value = dict_values[0];
			}
		}
		
		static string get_group_desc()
		{
			return "The base class for morphological image processors.";
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
		
	/**  Binarize an image based on the circular average around each pixel in real space.
	 *   The pixel is set to 1 when the ring average around the pixel keeps decreasing for n 
	 *   pixels as the radius of the ring increases. Here n is the threshold. This essentially 
	 *   picks out the local maximum pixels with some noise tolerance.
	 *   @author: Muyuan Chen
	 *   @date: 03/2015
	 */
	class CircularAverageBinarizeProcessor:public Processor
	{
	public:
		virtual void process_inplace(EMData * image);
		virtual EMData* process(const EMData* const image);

		virtual string get_name() const
		{
			return NAME;
		}
		static Processor *NEW()
		{
			return new CircularAverageBinarizeProcessor( );
		}
		string get_desc() const
		{
			return "Binarize an image based on the circular average around each pixel in real space.";
		}
		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("thresh", EMObject::INT, "The radius threshold that the circular average of density keep dropping.");
			return d;
		}
		static const string NAME;
	};
	
	/**  Replace the value of each pixel with the sum of density of the object it belongs to. 
	 *   Objects are defined by continius density above the given threshold.
	 *   @author: Muyuan Chen
	 *   @date: 03/2015
	 */	
	class ObjDensityProcessor:public Processor
	{
	public:
		virtual void process_inplace(EMData * image);
		virtual EMData* process(const EMData* const image);

		virtual string get_name() const
		{
			return NAME;
		}
		static Processor *NEW()
		{
			return new ObjDensityProcessor();
		}
		string get_desc() const
		{
			return "Sum of density of each object above threshold. Treats a 3D volume as 2D slices.";
		}
		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("thresh", EMObject::FLOAT, "The threshold to seperate objects.");
			return d;
		}
		static const string NAME;
	};
	
	/**  Label each object in a black-white image. Also return the center of each object.
	 *   @author: Muyuan Chen
	 *   @date: 03/2015
	 */	
	class ObjLabelProcessor:public Processor
	{
	public:
		virtual void process_inplace(EMData * image);
		virtual EMData* process(const EMData* const image);

		virtual string get_name() const
		{
			return NAME;
		}
		static Processor *NEW()
		{
			return new ObjLabelProcessor();
		}
		string get_desc() const
		{
			return "Label each object above threshold. Also return the center of each object. Treats a 3D volume as 2D slices.";
		}
		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("write_centers", EMObject::BOOL, "Write the center of each object in the attribute obj_centers.");
			return d;
		}
		static const string NAME;
	};
	
	/**  Thinning a binary map to skelton using the Zhang-Suen thinning algorithm. (1984, ACM)
	 *   @author: Muyuan Chen
	 *   @date: 03/2015
	 */	
	class BwThinningProcessor:public Processor
	{
	public:
		virtual void process_inplace(EMData * image);
		virtual EMData* process(const EMData* const image);

		virtual string get_name() const
		{
			return NAME;
		}
		static Processor *NEW()
		{
			return new BwThinningProcessor();
		}
		string get_desc() const
		{
			return "Thinning a binary map to skelton using the Zhang-Suen thinning algorithm.";
		}
		int process_pixel(float* data, float* array, int step);
		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("thresh", EMObject::FLOAT, "The threshold to binarize the map.");
			d.put("verbose", EMObject::INT, "Verbose");
			d.put("preserve_value", EMObject::BOOL, "The value on the skeleton is the same as the original map.");
			d.put("ntimes", EMObject::INT, "Number of iterations in the thinning process. Default: -1 (perform thinning until the image is skeltonized");
			return d;
		}
		static const string NAME;
	};

	/**  Set a pixel to white when >= N neighbors are white.
	 *   @author: Muyuan Chen
	 *   @date: 04/2015
	 */	
	class BwMajorityProcessor:public BoxStatProcessor
	{
	public:
		virtual string get_name() const
		{
			return NAME;
		}
		static Processor *NEW()
		{
			return new BwMajorityProcessor();
		}
		string get_desc() const
		{
			return "Set a pixel to white when >= N neighbors are white.";
		}
		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("thresh", EMObject::FLOAT, "The threshold to binarize the map.");
			d.put("nmaj", EMObject::INT, "Number of neighbors needed to set to white.");
			return d;
		}
		static const string NAME;
		
	protected:
		void process_pixel(float *pixel, const float *array, int n) const
		{
			float thresh=params.set_default("thresh",0);
			int nmaj=params.set_default("nmaj",n/2+1);
			for (int i=0; i<n; i++){
				if (array[i]>thresh)
					nmaj--;
			}
			*pixel=nmaj<=0?1:0;			
		}
	};
	
	
	/**  Prune branches from the skeleton. Remove a piece when the minimum distance through density 
	 *   from an endpoint to the nearest branch point is shorter than a given value.
	 *   @author: Muyuan Chen
	 *   @date: 04/2015
	 */	
	class PruneSkeletonProcessor:public Processor
	{
	public:
		virtual void process_inplace(EMData * image);
		virtual EMData* process(const EMData* const image);

		virtual string get_name() const
		{
			return NAME;
		}
		static Processor *NEW()
		{
			return new PruneSkeletonProcessor();
		}
		string get_desc() const
		{
			return "Prune branches from the skeleton. Remove a piece when the minimum distance through density from an endpoint to the nearest branch point is shorter than a given value.";
		}
		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("thresh", EMObject::FLOAT, "The threshold to binarize the map.");
			d.put("verbose", EMObject::INT, "Verbose");
			d.put("maxdist", EMObject::INT, "Maximum distance from the endpoint to branchpoint");
			return d;
		}
		static const string NAME;
	};
	
	/**  Grow a skeleton map toward a local direction. Image should be binarized.
	 *   @author: Muyuan Chen
	 *   @date: 06/2015
	 */	
	class GrowSkeletonProcessor:public Processor
	{
	public:
		virtual void process_inplace(EMData * image);
		virtual EMData* process(const EMData* const image);

		virtual string get_name() const
		{
			return NAME;
		}
		static Processor *NEW()
		{
			return new GrowSkeletonProcessor();
		}
		string get_desc() const
		{
			return "Grow a skeleton map toward a local direction.";
		}
		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("verbose", EMObject::INT, "Verbose");
			d.put("radius", EMObject::INT, "Half of the box size to determine the local direction.");
			return d;
		}
		static const string NAME;
	};
	
	/**  Calculate the z thickness of each pixel in a binarized 3d image
	 *   @author: Muyuan Chen
	 *   @date: 02/2016
	 */	
	class ZThicknessProcessor:public Processor
	{
	public:
		virtual void process_inplace(EMData * image);
		virtual EMData* process(const EMData* const image);

		virtual string get_name() const
		{
			return NAME;
		}
		static Processor *NEW()
		{
			return new ZThicknessProcessor();
		}
		string get_desc() const
		{
			return "Calculate the z thickness of each pixel in a binarized 3d image.";
		}
		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("thresh", EMObject::FLOAT, "Threshold for binarization");
			return d;
		}
		static const string NAME;
	};
	
	/** Replace the value of each pixel with a value in a given array.
	 * i.e. given an array of [3,7,9], pixels with value of 0 will become 3, 1 becomes 7, 2 becomes 9.
	 *   @author: Muyuan Chen
	 *   @date: 02/2016
	 */
	class ReplaceValuefromListProcessor:public RealPixelProcessor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}
		static Processor *NEW()
		{
			return new ReplaceValuefromListProcessor();
		}
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("num", EMObject::INT, "Length of the array.");
			d.put("colorlst", EMObject::FLOATARRAY, "Array of values to replace.");
			return d;
		}
		static const string NAME;

	  protected:
		void process_pixel(float *x) const
		{
			int num=params["num"];
			vector < float >lst =params["colorlst"];
			if (*x<num){
				*x=lst[int(*x)];
			}
			else{
				*x=0;
			}
		}

		string get_desc() const
		{
			return "Replace the value of each pixel with a value in a given array, i.e. given an array of [3,7,9], pixels with value of 0 will become 3, 1 becomes 7, 2 becomes 9. The input image has to be int, or it will be round down. Values exceed the length of array are set to zero. Designed for labeled image coloring.";
		}
	};
	
#ifdef SPARX_USING_CUDA
	/* class MPI CUDA kmeans processor
	 * 2009-02-13 17:34:45 JB first version
	 * 2009-09-02 11:19:10 JB for MPI version
	 * python wrap for GPU cluster
	 */
	class MPICUDA_kmeans {
	public:
		MPICUDA_kmeans();
		~MPICUDA_kmeans();
		int setup(int extm, int extN, int extn, int extK, int extn_start);
		void append_flat_image(EMData* im, int pos);
		int init_mem(int numdev);
		void compute_im2();
		int random_ASG(long int rnd);
		vector<int> get_ASG();
		vector<int> get_asg();
		void compute_NC();
		vector<int> get_NC();
		void set_ASG(const vector <int>& ASG);
		void set_NC(const vector <int>& NC);
		int get_ct_im_mv();
		void set_T(float extT);
		float get_T();
		void compute_AVE();
		void set_AVE(EMData* im, int pos);
		vector<EMData*> get_AVE();
		int one_iter();
		//int one_iter_SSE();
		//int AVE_to_host();
		int one_iter_SA();
		vector<float> compute_ji();
		vector<float> compute_criterion(const vector <float>& Ji);
		int shutdown();
	private:
		// params
		int m;
		int N;
		int n;
		int K;
		int nb_part;
		int n_start;
		int size_im;
		int size_IM;
		int size_AVE;
		int size_dist;
		int BLOCK_SIZE;
		int NB;
		int ins_BLOCK;
		int ite;
		float T;
		// debug
		int ct_im_mv;
		// host memory
		float* h_IM;
		float* h_im;
		float* h_AVE;
		float* h_dist;
		float* h_AVE2;
		float* h_im2;
		unsigned short int* h_ASG;
		unsigned short int* h_asg;
		unsigned int* h_NC;
		int* params;
		float ttt;
		// device memory
		float* d_im;
		float* d_AVE;
		float* d_dist;
		//int init_dist(); // intial h_dist and d_dist for SSE
                float compute_tt();
	};

#endif //EMAN2_USING_CUDA

#if 0

	class XYZProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return NAME;
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

		static const string NAME;
	};


#endif


#if 0

	class XYZProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return NAME;
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

		static const string NAME;
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

