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

    /** Averager is the base class for all averager classes.
     * Each subclass averager defines a way to do averaging on a set
     * of images. Each specific averager has a unique ID name. This name
     * is used to call a averager.
     *
     * Typical usages of Averager:
     *
     * 1. How to get all Averager types
     *
     *    vector<string> all_averagers = Factory<Averager>.get_list();
     *
     * 2. How to use an Averager
     *
     *    Averager *imgavg = Factory<Averager>.get("Image");
     *    vector<EMData *> images(2);
     *    EMData *image1 = ...;
     *    EMData *image2 = ...;
     *    images[0] = image1;
     *    images[1] = image2;
     *    EMData *result = imgavg->average(images);
     *
     * 3. How to define a new XYZAverager
     *
     *    XYZAverager should extend Averager and implement the
     *    following functions:
     *
     *        EMData *average(const vector<EMData *> & image_list) const;
     *        string get_name() const { return "XYZ"; }
     *        static Averager *NEW() { return new XYZAverager(); }
     */
    class Averager
    {
    public:
	virtual ~Averager() { }
	virtual EMData *average(const vector<EMData *> & image_list) const = 0;

	virtual string get_name() const = 0;
	
	virtual void set_params(const Dict & new_params)
	{
	    params = new_params;
	}
	
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
    class ImageAverager : public Averager
    {
    public:
	EMData * average(const vector<EMData *> & image_list) const;

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
    };

    /** IterationAverager averages images by doing the smoothing iteration.
     */
    class IterationAverager : public Averager
    {
    public:
	EMData * average(const vector<EMData *> & image_list) const;

	string get_name() const
	{
	    return "Iteration";
	}

	static Averager *NEW()
	{
	    return new IterationAverager();
	}
    };

    /** CtfAverager is the base Averager class for CTF correction or SNR weighting.
    */
    class CtfAverager : public Averager
    {
    public:
	EMData * average(const vector<EMData *> & image_list) const;

	vector<float> get_snr() const
	{
	    return snr;
	}
	
    protected:
	CtfAverager() : sf(0), curves(0), need_snr(false) 
	{
	}
	
	XYData *sf;
	EMData *curves;
	bool need_snr;
	const char *outfile;
    private:
	mutable vector<float> snr;
    };

    /** WeightingAverager averages the images with SNR weighting, but no CTF correction.
     */
    class WeightingAverager : public CtfAverager
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
    class CtfCAverager : public CtfAverager
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
    class CtfCWAverager : public CtfAverager
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
	    need_snr = (bool) (int)params["need_snr"];
	}
    };

    /** CtfCWautoAverager averages the images with CTF correction with a Wiener filter.
     *  The Weiner filter is estimated directly from the data.
     */
    class CtfCWautoAverager : public CtfAverager
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

    template<> Factory<Averager>::Factory();

    void dump_averagers();
    
}


#endif
