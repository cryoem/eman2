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

	/** Cmp class defines image comparison method. Before doing the
	 * comparison, an optional transformation may be used to
	 * transform the 2 images. The bigger the comparison result, the
	 * more similar of the 2 images.
	 *
	 * Cmp class is the base comparison class. Each specific
     * comparison class is a subclass of Cmp, and must have a unique
     * name. The name is used to  create a new Cmp instance or call a Cmp.
	 *
	 * All Cmp classes in EMAN are managed by a Factory pattern. So 
	 * each Cmp class must define:
	 *   a) a unique name to idenfity itself in the factory.
	 *   b) a static method to register itself in the factory.
	 *
	 *
     * Typical usage of Cmp:
     *
     * 1. How to get all Cmp names
     *
     *      vector<string> all_cmps = Factory<Cmp>::get_list();
     *
     * 2. How to use a Cmp
     *
     *      EMData *image1 = ...;
     *      EMData *image2 = ...;
     *      float result = image1->cmp("CMP_NAME", Dict("with", image2));
     *
     * 3. How to define a new Cmp class
     *
     *    A new XYZCmp class should implement the following functions:
	 *    (Please replace 'XYZ' with your own class name).
	 *
     *        float cmp(EMData * image, Transform * transform = 0) const;
     *        TypeDict get_param_types() const;
     *        string get_name() const { return "XYZ"; }
     *        static Cmp *NEW() { return XYZCmp(); }
     */

	class Cmp
	{
	  public:
		virtual ~ Cmp()
		{
		}
		
		/** To compare 'image' with another image passed in through
		 * its parameters. An optional transformation may be used
		 * to transform the 2 images.
		 *
		 * @param image The image to be compared.
		 * @param transform Defines the transformation.
		 * @return The comparison result. The bigger, the better.
		 */			
		virtual float cmp(EMData * image, Transform * transform = 0) const = 0;
		
		/** Get the Cmp's name. Each Cmp is identified by a unique name.
		 * @return The Cmp's name.
		 */
		virtual string get_name() const = 0;

		/** Get the Cmp parameters in a key/value dictionary.
		 * @return A key/value pair dictionary containing the parameters.
		 */
		virtual Dict get_params() const
		{
			return params;
		}

		/** Set the Cmp parameters using a key/value dictionary.
		 * @param new_params A dictionary containing the new parameters.
		 */
		virtual void set_params(const Dict & new_params)
		{
			params = new_params;
		}

		/** Get Cmp parameter information in a dictionary. Each
		 * parameter has one record in the dictionary. Each record
		 * contains its name, data-type, and description.
		 *
		 * @return A dictionary containing the parameter info.
		 */	 
		virtual TypeDict get_param_types() const = 0;
		
	  protected:
		  mutable Dict params;
	};

	/** Use dot product of 2 same-size images to do the comparison.
     * For complex images, it does not check r/i vs a/p.
    */
	class DotCmp:public Cmp
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
	class VarianceCmp:public Cmp
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
	 *
     * Use Phase Residual as a measure of similarity
     * Differential phase residual (DPR) is a measure of statistical
     * dependency between two averages, computed over rings in Fourier
     * space as a function of ring radius (= spatial frequency, or resolution) 
     */
	class PhaseCmp:public Cmp
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

	/** FRCCmp returns a quality factor based on FRC between images.
     *  Fourier ring correlation (FRC) is a measure of statistical
     * dependency between two averages, computed by comparison of
     * rings in Fourier space. 1 means prefect agreement. 0 means no
     * correlation.    
     */
	class FRCCmp:public Cmp
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

	template <> Factory < Cmp >::Factory();

	void dump_cmps();
}


#endif
