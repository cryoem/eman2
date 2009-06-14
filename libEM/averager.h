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
		virtual void add_image(EMData * image) = 0;

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
			return "mean";
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
			d.put("ignore0", EMObject::INT, "if set, ignore zero value pixels");
			return d;
		}

		virtual void mult(const float&) { }

	private:
		EMData *sigma_image;
		int *nimg_n0;
		int ignore0;
		int nimg;
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
			return "minmax";
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
			return d;
		}

		virtual void mult(const float&) { }

	private:
		int max;
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
			return "iteration";
		}

		string get_desc() const
		{
			return "Unknown";
		}

		static Averager *NEW()
		{
			return new IterationAverager();
		}
	private:
		EMData * sigma_image;
		int nimg;
	};

	/** CtfAverager is the base Averager class for CTF correction or SNR weighting.
    */
	class CtfAverager:public Averager
	{
	  public:
		CtfAverager();
		void add_image( EMData * image);
		EMData * finish();

		vector < float >get_snr() const
		{
			return snr;
		}

	  protected:
		XYData *sf;
		EMData *curves;
		bool need_snr;
		const char *outfile;
	  private:
		mutable vector < float >snr;
		EMData * image0_fft;
		EMData * image0_copy;

		vector<vector<float> > ctf;
		vector<vector<float> > ctfn;

		float *snri;
		float *snrn;
		float *tdr;
		float *tdi;
		float *tn;

		float filter;
		int nimg;
		int nx;
		int ny;
		int nz;
	};

	/** WeightingAverager averages the images with SNR weighting, but no CTF correction.
	 *@param curves
	 *@param sf
     */
	class WeightingAverager:public CtfAverager
	{
	  public:
		string get_name() const
		{
			return "snr_weight";
		}

		string get_desc() const
		{
			return "SNR Weighted average without CTF amplitude correction";
		}

		static Averager *NEW()
		{
			return new WeightingAverager();
		}

		void set_params(const Dict & new_params)
		{
			params = new_params;
			curves = params["curves"];
			sf = params["sf"];
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("curves", EMObject::EMDATA);
			d.put("sf", EMObject::XYDATA);
			return d;
		}
	};

	/** CtfCAverager averages the images with CTF correction.
     */
	class CtfCAverager:public CtfAverager
	{
	  public:
		string get_name() const
		{
			return "ctfc";
		}

		string get_desc() const
		{
			return "CTF amplitude corrected average, including SNR weight, but result is unprocessed.";
		}

		static Averager *NEW()
		{
			return new CtfCAverager();
		}
	};

	/** CtfCWAverager averages the images with CTF correction.
     */
	class CtfCWAverager:public CtfAverager
	{
	  public:
		string get_name() const
		{
			return "ctfcw";
		}

		string get_desc() const
		{
			return "CTF amplitude corrected average, including SNR weight and Wiener filtration";
		}

		static Averager *NEW()
		{
			return new CtfCWAverager();
		}

		void set_params(const Dict & new_params)
		{
			params = new_params;
			if ((int) params["need_snr"]) {
				need_snr = true;
			}
			else {
				need_snr = false;
			}
		}
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
			return "ctf.auto";
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
	  protected:
		EMData *snrsum;   // contains the summed SNR for the average
		int nimg;
	};

	template <> Factory < Averager >::Factory();

	void dump_averagers();
	map<string, vector<string> > dump_averagers_list();
}


#endif
