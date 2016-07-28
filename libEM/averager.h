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

#ifndef eman_averager_h__
#define eman_averager_h__ 1

#include "emobject.h"
#include "emdata.h"

#include <vector>
using std::vector;

namespace EMAN
{
	class EMData;
	class XYData;

	/** Averager class defines a way to do averaging on a set
     * of images. A user may add one or more images to the Averager at
     * one time. The images are averaged at the time of adding to the
     * Averager. After all images are added, Average will return the
     * averaged result.
	 *
     * Averager class is the base class for all averager classes. Each
     * specific averager has a unique ID name. This name is used to
     * call a averager.
	 *
	 * All Averager classes in EMAN are managed by a Factory
	 * pattern. So each Averager class must define:
	 *   - a unique name to idenfity itself in the factory.
	 *   - a static method to register itself in the factory.
	 *
     * Typical usages of Averager:
     *
     *  - How to get all Averager types
     @code
     *    vector<string> all_averagers = Factory<Averager>::get_list();
     @endcode
	 *
     *  - How to use an Averager
     @code
     *    Averager *imgavg = Factory<Averager>::get("image");
     *    vector<EMData *> images(2);
     *    EMData *image1 = ...;
     *    EMData *image2 = ...;
     *    images[0] = image1;
     *    images[1] = image2;
	 *    imgavg->add_image(image1);
	 *    imgavg->add_image(image2);
     *    EMData *result = imgavg->finish();
	 @endcode
     *
     *  - How to define a new XYZAverager \n
     *    XYZAverager should extend Averager and implement the
     *    following functions:
     *    (Please replace 'XYZ' with your own class name).
	 @code
     *        void add_image(EMData * image);
	 *        EMData * finish();
     *        string get_name() const { return "XYZ"; }
     *        static Averager *NEW() { return new XYZAverager(); }
	 @endcode
     */
	class Averager
	{
	  public:
		Averager() : result(0) {}

		virtual ~ Averager()
		{
		}

		/** To add an image to the Averager. This image will be
		 * averaged in this function.
		 * @param image The image to be averaged.
		 */
		virtual void add_image(EMData * image) ;

		/** To add multiple images to the Averager. All the
		 * newly-added images are averaged in this function.
		 * @param images The images to be averaged.
		 */
		virtual void add_image_list(const vector<EMData*> & images);

		/** Finish up the averaging and return the result.
		 *
		 * @return The averaged image.
		 */
		virtual EMData * finish() = 0;

		/** Get the Averager's name. Each Averager is identified by a unique name.
		 * @return The Averager's name.
		 */
		virtual string get_name() const = 0;

		virtual string get_desc() const = 0;

		/** Set the Averager parameters using a key/value dictionary.
		 * @param new_params A dictionary containing the new parameters.
		 */
		virtual void set_params(const Dict & new_params)
		{
			params = new_params;
		}

		/** Multiply the result image by some floating point constant
		 * This is useful when weighting the input images prior to
		 * calling add_image - a situation where it is likely you
		 * want to divide by the sum of the weights. Hence call mult
		 * after all of the weighted images have been added.
		 * @param s the scaling factor.
		 * @exception NullPointerException if the EMData pointer (result) is NULL
		 */
		virtual void mult(const float& s);

		/** Get Averager  parameter information in a dictionary. Each
		 * parameter has one record in the dictionary. Each record
		 * contains its name, data-type, and description.
		 *
		 * @return A dictionary containing the parameter info.
		 */
		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			return d;
		}
		
	  protected:
		mutable Dict params;
		EMData *result;
	};

	/** ImageAverager averages a list of images. It optionally makes
     * a sigma image.
     *@param sigma sigma value
     *@param ignore0 if set, ignore zero value pixels
     */
	class ImageAverager:public Averager
	{
	  public:
		ImageAverager();

		void add_image( EMData * image);
		EMData * finish();

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Simple mean average of images";
		}

		static Averager *NEW()
		{
			return new ImageAverager();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("sigma", EMObject::EMDATA, "sigma value");
			d.put("normimage", EMObject::EMDATA, "In conjunction with ignore0, the number of non zero values for each pixel will be stored in this image.");
			d.put("ignore0", EMObject::INT, "if set, ignore zero value pixels");
			return d;
		}

		virtual void mult(const float&) { }

		static const string NAME;
		
	private:
		EMData *sigma_image,*normimage;
		int ignore0;
		int nimg;
		int freenorm;
	};

	/** FourierWeightAverager makes an average of a set of images in Fourier space using a per-image radial weight. The provided XYData object for each inserted
	 * image should range from x=0 - 0.5*sqrt(2), and contains the radial weights from 0 - Nyquist at the point. If the x range is insufficient, values will be
	 * clamped at the ends of the available x-range. 2-D Images only, but will work with rectangular images.
     *@param normimage	After finish() will contain the sum of the weights in each Fourier location. Size must be ((nx+1)/2,y)
     */
	class FourierWeightAverager:public Averager
	{
	  public:
		FourierWeightAverager();

		void add_image( EMData * image);
		EMData * finish();

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Weighted mean of images in Fourier space. Each image must have weighting curve in its header, an XYData object called 'avg_weight'.";
		}

		static Averager *NEW()
		{
			return new FourierWeightAverager();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
//			d.put("weight", EMObject::XYDATA, "Radial weight. X: 0 - 0.5*sqrt(2). Y contains weights.");
			d.put("normimage", EMObject::EMDATA, "After finish() will contain the sum of the weights in each Fourier location. Size must be ((nx+1)/2,y)");
			return d;
		}

		static const string NAME;
		
	private:
		EMData *normimage;
		int freenorm;
		int nimg;
	};

	
	/** TomoAverager averages a list of volumes in Fourier space. It excludes values near zero
	 * from the average, assuming they are part of the missing-cone/wedge. A threshold
     *@param thresh_sigma a number, multiplied by the standard deviation of the image, below-which values are considered zero
     */
	class TomoAverager:public Averager
	{
	  public:
		TomoAverager();

		void add_image( EMData * image);
		EMData * finish();

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Average of volumes in Fourier space, excluding any pixels with near 0 intensity.";
		}

		static Averager *NEW()
		{
			return new TomoAverager();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("thresh_sigma", EMObject::FLOAT, "multiplied by the standard deviation of the image, below-which values are considered zero. Default = .01");
			d.put("save_norm", EMObject::INT, "If set, will save the normalization volume as norm.hdf. Mainly for debugging purposes.");
			return d;
		}

		virtual void mult(const float&) { }

		static const string NAME;
		
	private:
		EMData *norm_image;
		float thresh_sigma;
	};

	
	/** ImageAverager averages a list of images. It optionally makes
     * a sigma image.
     *@param max If set, will find the max value, otherwise finds min
     */
	class MinMaxAverager:public Averager
	{
	  public:
		MinMaxAverager();

		void add_image( EMData * image);
		EMData * finish();

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Finds the minimum or maximum value in each pixel";
		}

		static Averager *NEW()
		{
			return new MinMaxAverager();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("max", EMObject::INT, "If set, will find the max value, otherwise finds min");
			d.put("owner", EMObject::EMDATA, "Contains the number of the input image which 'owns' the max/min value. Value will be insertion sequence number unless 'ortid' is set in each image being averaged.");
			return d;
		}

		virtual void mult(const float&) { }
		
		static const string NAME;
		
	private:
		int max;
		int nimg;
	};

	/** AbsMaxMinAverager averages a list of images to the maximum(or minimum of the absolute pixel value)
	 *  It optionally makes a sigma image.
     *@param min If set, will find the min value, otherwise finds max
     */
	class AbsMaxMinAverager:public Averager
	{
	public:
		AbsMaxMinAverager();

		void add_image( EMData * image);
		EMData * finish();

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Average to maximum(or minimum if set parameter 'min' to non-zero) absolute value in each pixel";
		}

		static Averager *NEW()
		{
			return new AbsMaxMinAverager();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("min", EMObject::INT, "If set, will average to minimum absolute value, by default average to max");
			return d;
		}

		static const string NAME;

	private:
		int min;
		int nimg;
	};

	/** IterationAverager averages images by doing the smoothing iteration.
     */
	class IterationAverager:public Averager
	{
	  public:
		IterationAverager();
		void add_image( EMData * image);
		EMData * finish();

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Unknown";
		}

		static Averager *NEW()
		{
			return new IterationAverager();
		}
		
		static const string NAME;
		
	private:
		EMData * sigma_image;
		int nimg;
	};

	/** CtfWtAverager
     */
	class CtfWtAverager:public Averager
	{
	  public:
	    CtfWtAverager();

		void add_image( EMData * image);
		EMData * finish();

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Average without CTF correction but with CTF weighting. Smoothed SNR can still have large uncertainty, so weighting by envelope-free CTF may provide more uniform results.";
		}

		static Averager *NEW()
		{
			return new CtfWtAverager();
		}

		void set_params(const Dict & new_params)
		{
			params = new_params;
//			outfile = params["outfile"];
		}
		
		static const string NAME;
		
	  protected:
		EMData *ctfsum;   // contains the summed SNR for the average
		int nimg;
	};

		/** CtfWtAverager
     */
	class CtfWtFiltAverager:public Averager
	{
	  public:
	    CtfWtFiltAverager();

		void add_image( EMData * image);
		EMData * finish();

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Average without CTF correction but with CTF weighting and automatic filter estimated from the data. Smoothed SNR can still have large uncertainty, so weighting by envelope-free CTF may provide more uniform results.";
		}

		static Averager *NEW()
		{
			return new CtfWtFiltAverager();
		}

		void set_params(const Dict & new_params)
		{
			params = new_params;
//			outfile = params["outfile"];
		}
		
		static const string NAME;
		
	  protected:
		EMData *results[2];		// even/odd split for filter estimate
		EMData *ctfsum[2];		// contains the summed SNR for the average
		int nimg[2],eo;
	};


	/** CtfCWautoAverager averages the images with CTF correction with a Wiener filter.
     *  The Weiner filter is estimated directly from the data.
     */
	class CtfCAutoAverager:public Averager
	{
	  public:
	    CtfCAutoAverager();

		void add_image( EMData * image);
		EMData * finish();

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Averaging with automatic CTF correction and SNR weight. No B-factor correction (as this is best done in 3-D). Bases estimated SSNR on CTF parameters, so requires EMAN2 CTF parameters.";
		}

		static Averager *NEW()
		{
			return new CtfCAutoAverager();
		}

		void set_params(const Dict & new_params)
		{
			params = new_params;
//			outfile = params["outfile"];
		}
		
		static const string NAME;
		
	  protected:
		EMData *snrsum;   // contains the summed SNR for the average
		int nimg;
	};

	/** CtfCWautoAverager averages the images with CTF correction with a Wiener filter.
     *  The Weiner filter is estimated directly from the data.
     */
	class CtfCWautoAverager:public Averager
	{
	  public:
	    CtfCWautoAverager();

		void add_image( EMData * image);
		EMData * finish();

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Averaging with autmatic CTF correction. Does not require a structure factor, but only works with EMAN2's CTF model";
		}

		static Averager *NEW()
		{
			return new CtfCWautoAverager();
		}

		void set_params(const Dict & new_params)
		{
			params = new_params;
//			outfile = params["outfile"];
		}
		
		static const string NAME;
		
	  protected:
		EMData *snrsum;   // contains the summed SNR for the average
		int nimg;
	};

	/** 
	* VarianceAverager computes the pixel-wise variance of a list of images.
        */
        class VarianceAverager:public Averager
	{
	  public:
		VarianceAverager();

		void add_image( EMData * image);
		EMData * finish();

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Pixel-wise variance of images";
		}

		static Averager *NEW()
		{
			return new VarianceAverager();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			return d;
		}

		virtual void mult(const float&) { }

		static const string NAME;
		
	private:
		EMData *mean;
		int nimg;
	};

	template <> Factory < Averager >::Factory();

	void dump_averagers();
	map<string, vector<string> > dump_averagers_list();
}

#endif
