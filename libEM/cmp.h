/**
 * $Id$
 */
#ifndef eman_cmp__h__
#define eman_cmp__h__ 1


#include "emobject.h"

namespace EMAN
{
    
    class EMData;
    class Transform;
    
    /** Cmp is the base class for image Comparison methods.
     * The bigger the comparison result, the better. Each specific
     * coomparison class has a unique name. This name is used to
     * create a new Cmp instance or call a Cmp.
     *
     * Typical usage of Cmp:
     *
     * 1. How to get all Cmp types
     *
     *      vector<string> all_cmps = Factory<Cmp>.instance()->get_list();
     *
     * 2. How to use a Cmp
     *
     *      EMData *img = ...;
     *      EMData *with = ...;
     *      float result = img->cmp("CMP_NAME", Dict("with", with));
     *
     * 3. How to define a new Cmp class
     *
     *    A new XYZCmp class should implement the following functions:
     *        float cmp(EMData * image, Transform * transform = 0) const;
     *        TypeDict get_param_types() const;
     *        string get_name() const { return "XYZ"; }
     *        static Cmp *NEW() { return XYZCmp(); }
     */
    
    class Cmp
    {
    public:
	virtual ~Cmp() { }
	virtual float cmp(EMData * image, Transform * transform = 0) const = 0;
	virtual TypeDict get_param_types() const = 0;
	virtual string get_name() const = 0;
	
	virtual Dict get_params() const
	{
	    return params;
	}

	virtual void set_params(const Dict & new_params)
	{
	    params = new_params;
	}
	
    protected:
	mutable Dict params;
    };

    /** Use dot product of 2 same-size images to do the comparison.
     * complex does not check r/i vs a/p.
    */
    class DotCmp : public Cmp
    {
    public:
	float cmp(EMData * image, Transform * transform = 0) const;

	string get_name() const
	{
	    return "Dot";
	}

	static Cmp *NEW() 
	{
	    return new DotCmp();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);
	    d.put("evenonly", EMObject::INT);
	    return d;
	}

    };

    /** Linear comparison of 2 data sets. 'image' should be noisy and
     * 'with' should be less noisy. Scaling of 'this' is determined to
     * make the density histogram of the difference between the data
     * sets as symmetric as possible scale will optionally return
     * the scale factor which would be multiplied by 'this' to achieve
     * this normalization shift will return the corresponding shift.
     * If modifying 'this', scale should be applied first, then b
     * should be added
     */
    class VarianceCmp : public Cmp
    {
    public:
	float cmp(EMData * image, Transform * transform = 0) const;

	string get_name() const
	{
	    return "Variance";
	}

	static Cmp *NEW()
	{
	    return new VarianceCmp();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);
	    d.put("keepzero", EMObject::INT);
	    return d;
	}
    };
    
    /** Amplitude weighted mean phase difference (radians) with optional
     * SNR weight. SNR should be an array as returned by ctfcurve()
     * 'data' should be the less noisy image, since it's amplitudes 
     * will be used to weight the phase residual. 2D only.

     * Use Phase Residual as a measure of similarity
     * Differential phase residual (DPR) is a measure of statistical
     * dependency between two averages, computed over rings in Fourier
     * space as a function of ring radius (= spatial frequency, or resolution) 
     */
    class PhaseCmp : public Cmp
    {
    public:
	float cmp(EMData * image, Transform * transform = 0) const;

	string get_name() const
	{
	    return "Phase";
	}

	static Cmp *NEW()
	{
	    return new PhaseCmp();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);

	    return d;
	}
    };

    /** returns a quality factor based on FRC between images.
     *  Fourier ring correlation (FRC) is a measure of statistical
     * dependency between two averages, computed by comparison of
     * rings in Fourier space. 1 means prefect agreement. 0 means no
     * correlation.    
     */
    class FRCCmp : public Cmp
    {
    public:
	float cmp(EMData * image, Transform * transform = 0) const;

	string get_name() const
	{
	    return "FRC";
	}

	static Cmp *NEW()
	{
	    return new FRCCmp();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);
	    d.put("snr", EMObject::FLOATARRAY);
	    return d;
	}
    };

    template<> Factory<Cmp>::Factory();
}


#endif
