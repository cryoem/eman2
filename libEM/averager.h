/**
 * $Id$
 */
#ifndef eman_averager_h__
#define eman_averager_h__ 1

#include <vector>
using std::vector;

#include "emobject.h"

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
     *    Averager *imgavg = Factory<Averager>::get("Image");
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

		/** Set the Averager parameters using a key/value dictionary.
		 * @param new_params A dictionary containing the new parameters.
		 */
		virtual void set_params(const Dict & new_params)
		{
			params = new_params;
		}

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
	};

	/** ImageAverager averages a list of images. It optionally makes
     * a sigma image.
     */
	class ImageAverager:public Averager
	{
	  public:
		ImageAverager();
		
		void add_image( EMData * image);
		EMData * finish();
		
		string get_name() const
		{
			return "Image";
		}

		static Averager *NEW()
		{
			return new ImageAverager();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("sigma", EMObject::EMDATA);
			d.put("ignore0", EMObject::INT);
			return d;
		}

	private:
		EMData *result;
		EMData *sigma_image;
		int *nimg_n0;
		int ignore0;
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
			return "Iteration";
		}

		static Averager *NEW()
		{
			return new IterationAverager();
		}
	private:
		EMData * result;
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
		EMData * result;
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
     */
	class WeightingAverager:public CtfAverager
	{
	  public:
		string get_name() const
		{
			return "Weighting";
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
			return "CtfC";
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
			return "CtfCW";
		}

		static Averager *NEW()
		{
			return new CtfCWAverager();
		}

		void set_params(const Dict & new_params)
		{
			params = new_params;
			need_snr = (bool) ((int) params["need_snr"]);
		}
	};

	/** CtfCWautoAverager averages the images with CTF correction with a Wiener filter.
     *  The Weiner filter is estimated directly from the data.
     */
	class CtfCWautoAverager:public CtfAverager
	{
	  public:
		string get_name() const
		{
			return "CtfCWauto";
		}

		static Averager *NEW()
		{
			return new CtfCWautoAverager();
		}

		void set_params(const Dict & new_params)
		{
			params = new_params;
			outfile = params["outfile"];
		}
	};

	template <> Factory < Averager >::Factory();

	void dump_averagers();

}


#endif
