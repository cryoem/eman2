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

#ifndef eman_cmp__h__
#define eman_cmp__h__ 1


#include "emobject.h"

namespace EMAN
{

	class EMData;
	/** Cmp class defines image comparison method. Using default
	 * arguments, smaller values indicate more similar images.
	 *
	 * Cmp class is the base comparison class. Each specific
     * comparison class is a subclass of Cmp, and must have a unique
     * name. The name is used to  create a new Cmp instance or call a Cmp.
	 *
	 * All Cmp classes in EMAN are managed by a Factory pattern. So
	 * each Cmp class must define:
	 *   - a unique name to idenfity itself in the factory.
	 *   - a static method to register itself in the factory.
	 *
	 *
     * Typical usage of Cmp:
     *
     *  - How to get all Cmp names
     @code
     *      vector<string> all_cmps = Factory<Cmp>::get_list();
     @endcode
	 *
     *  - How to use a Cmp
     @code
     *      EMData *image1 = ...;
     *      EMData *image2 = ...;
	 *      Dict params = ...;
     *      float result = image1->cmp("CMP_NAME", image2, params);
     @endcode
	 *
     *  - How to define a new Cmp class \n
     *    A new XYZCmp class should implement the following functions:
	 *    (Please replace 'XYZ' with your own class name).
	 @code
     *        float cmp(EMData * image, EMData * with) const;
     *        TypeDict get_param_types() const;
     *        string get_name() const { return "XYZ"; }
     *        static Cmp *NEW() { return XYZCmp(); }
	 @endcode
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
		 * @param image The first image to be compared.
		 * @param with The second image to be comppared.
		 * @return The comparison result. Smaller better by default
		 */
		virtual float cmp(EMData * image, EMData * with) const = 0;

		/** Get the Cmp's name. Each Cmp is identified by a unique name.
		 * @return The Cmp's name.
		 */
		virtual string get_name() const = 0;

		virtual string get_desc() const = 0;

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
		void validate_input_args(const EMData * image, const EMData *with) const;

		mutable Dict params;
	};

	/** Compute the cross-correlation coefficient between two images.
	 *
	 * The cross-correlation coefficient is defined as:
	 *       <AB> - <A><B>
	 * CCC = -------------
	 *       sig(A)sig(B)
	 *
	 * where the angle brackets denote averages and "sig" is the
	 * standard deviation.  In the case of a mask, only pixels under
	 * the mask are included in the calculation of averages.
	 *
	 * For complex images, this routine currently bails.
	 * @ingroup CUDA_ENABLED
	 * @author Grant Goodyear (grant.goodyear@uth.tmc.edu)
	 * @date 2005-10-03
	 * @param negative Returns -1 * ccc, default true
	 */
	class CccCmp:public Cmp
	{
	  public:
		float cmp(EMData * image, EMData * with) const;

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Cross-correlation coefficient (default -1 * ccc)";
		}

		static Cmp *NEW()
		{
			return new CccCmp();
		}

		//param mask Image mask
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("negative", EMObject::INT, "If set, returns -1 * ccc product. Set by default so smaller is better");
			d.put("mask", EMObject::EMDATA, "image mask");
			return d;
		}

		static const string NAME;
	};


	/* Compues the L^1 difference, after normalization.*/
	class LodCmp:public Cmp
	{
	  public:
		float cmp(EMData * image, EMData * with) const;

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "L^1 normalized difference (positive by default)";
		}

		static Cmp *NEW()
		{
			return new LodCmp();
		}

		//param mask Image mask
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("negative", EMObject::INT, "If set (which is the default), returns Lod. (The smaller the better)");
			d.put("normalize", EMObject::INT, "If set, normalizes maps prior to computing the difference. Default=0 (no normalization)");
			d.put("mask", EMObject::EMDATA, "image mask");
			return d;
		}

		static const string NAME;
	};


	/** Squared Euclidean distance normalized by n between 'this' and 'with'*/
	//  I corrected this as there is no such thing as "variance between two images"
	//  I corrected naive coding to avoid square
	//  Also, the equation in return statement was incorrect, grrrr!!!
	//  Finally, I added mask option  PAP 04/23/06
	class SqEuclideanCmp:public Cmp
	{
	  public:
		SqEuclideanCmp() {}

		float cmp(EMData * image, EMData * with) const;

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Squared Euclidean distance (sum(a - b)^2)/n.";
		}

		static Cmp *NEW()
		{
			return new SqEuclideanCmp();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("mask", EMObject::EMDATA, "image mask");
			d.put("zeromask", EMObject::INT, "If set, zero pixels in either image will be excluded from the statistics");
			d.put("normto",EMObject::INT,"If set, 'with' is normalized to 'this' before computing the distance");
			return d;
		}

		static const string NAME;
	};


	/** Use dot product of 2 same-size images to do the comparison.  // Added mask option PAP 04/23/06
	* For complex images, it does not check r/i vs a/p.
	* @ingroup CUDA_ENABLED
	* @author Steve Ludtke (sludtke@bcm.tmc.edu)
	* @date 2005-07-13
	* @param negative Returns -1 * dot product, default true
	* @param normalize Returns normalized dot product -1.0 - 1.0
    */
	class DotCmp:public Cmp
	{
	  public:
		float cmp(EMData * image, EMData * with) const;

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Dot product (default -1 * dot product)";
		}

		static Cmp *NEW()
		{
			return new DotCmp();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("negative", EMObject::INT, "If set, returns -1 * dot product. Set by default so smaller is better");
			d.put("normalize", EMObject::INT, "If set, returns normalized dot product (cosine of the angle) -1.0 - 1.0.");
			d.put("mask", EMObject::EMDATA, "image mask");
			return d;
		}
		
		static const string NAME;
	};

	/** This implements the technique of Mike Schmid where by the cross correlation is normalized
        * in an effort to remove the effects of the missing wedge. Somewhat of a heuristic solution, but it seems
        * to work. Basiaclly it relies on the observation that 'good' matchs will conentrate the correlation
        * signal in the peak, wheras 'bad' correlations will distribute the signal.
	* @ingroup CUDA_ENABLED
	* @author John Flanagan (a port of Mike Schmid's code - Mike Schmid is the intellectual author)
	* @date 2010-10-18
	* @param norm Normalize the ccf (you need to do this if a missing wedge is present)
	* @param ccf Supply your own ccf function (can be used for speedups)
	* @param searchx The maximum acceptable tx from the origin
	* @param searchy The maximum acceptable ty from the origin
	* @param searchz The maximum acceptable tz from the origin
	*/
	class TomoCccCmp:public Cmp
	{
	  public:
		virtual float cmp(EMData * image, EMData * with) const;

		virtual string get_name() const 
		{
			return NAME;
		}

		virtual string get_desc() const
		{
			return "Ccc with consideration given for the missing wedge";
		}

		static Cmp *NEW()
		{
			return new TomoCccCmp();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("norm", EMObject::BOOL,"Whether the cross correlation image should be normalized (should be for normalized images). Default is true.");
			d.put("ccf", EMObject::EMDATA,"The ccf image, can be provided if it already exists to avoid recalculating it");
			d.put("normalize", EMObject::EMDATA,"Return the negative value (which is EMAN2 convention), Defalut is true(1)");
			d.put("searchx", EMObject::INT, "The maximum range of the peak location in the x direction. Default is sizex/4");
			d.put("searchy", EMObject::INT, "The maximum range of the peak location in the y direction. Default is sizey/4");
			d.put("searchz", EMObject::INT, "The maximum range of the peak location in the z direction. Default is sizez/4");
			return d;
		}
		
		static const string NAME;
	};

	/** Both images should be FFTs. Computes a CCC excluding pixels with very small value (assumes they are missing wedge)
    * Returns overlap amount as "fft_overlap" in "image"
	* @date 2010-10-18
	*/
	class TomoWedgeCccCmp:public Cmp
	{
	  public:
		virtual float cmp(EMData * image, EMData * with) const;

		virtual string get_name() const 
		{
			return NAME;
		}

		virtual string get_desc() const
		{
			return "Ccc of two FFTs with missing wedge taken into account explicitly. Also stores overlap as 'fft_overlap' in \
first input image. Individual voxels are include only if they exceed the threshold value in both images. fft_overlap expresses \
the fractional number of voxels included in the comparison. Actual returned correlations already account for this effect, but \
the overlap can be used to estimate the relative uncertainty in the resulting correlation coefficient.";
		}

		static Cmp *NEW()
		{
			return new TomoWedgeCccCmp();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
// 			d.put("norm", EMObject::BOOL,"Whether the cross correlation image should be normalized (should be for normalized images). Default is true.");
// 			d.put("ccf", EMObject::EMDATA,"The ccf image, can be provided if it already exists to avoid recalculating it");
 			d.put("sigmaimg", EMObject::FLOAT, "Sigma coefficient for thresholding values included in the dot product. default = 0.5");
 			d.put("sigmawith", EMObject::FLOAT, "Sigma coefficient for thresholding values included in the dot product in the 'with' image. Default = 0.5");
			d.put("negative", EMObject::INT, "If set, returns -1 * ccc product. Set by default so smaller is better");
			return d;
		}
		
		static const string NAME;
	};

	/** Both images should be FFTs. Computes a CCC excluding pixels with very small value (assumes they are missing wedge)
    * Returns overlap amount as "fft_overlap" in "image"
	* @date 2010-10-18
	*/
	class TomoWedgeFscCmp:public Cmp
	{
	  public:
		virtual float cmp(EMData * image, EMData * with) const;

		virtual string get_name() const 
		{
			return NAME;
		}

		virtual string get_desc() const
		{
			return "FSC of two FFTs with missing wedge taken into account via a dynamic threshold. Also stores overlap as 'fft_overlap' in \
first input image on return. This comparator can take the missing wedge into account automatically without any precomputation of wedge \
orientation. If sigmaimg and/or sigmawith are provided, they must contain floating point arrays with at least as many elements as the \
Fourier radius in pixels, and these curves represent the thresholds for each image for including the values in the FSC, or considering them \
to be missing values. If not provided, these values will be computed automatically. These options are primarily provided for speed.";
		}

		static Cmp *NEW()
		{
			return new TomoWedgeFscCmp();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
// 			d.put("norm", EMObject::BOOL,"Whether the cross correlation image should be normalized (should be for normalized images). Default is true.");
// 			d.put("ccf", EMObject::EMDATA,"The ccf image, can be provided if it already exists to avoid recalculating it");
 			d.put("sigmaimg", EMObject::FLOATARRAY, "Res dependent amplitude threshold. If provided this must contain a floating point array with as many elements as the Fourier radius of the image. ");
 			d.put("sigmawith", EMObject::FLOATARRAY, "Resolution dependent coefficient for thresholding values included in the dot product in the 'with' image. Default = 0.1 and is normally fine");
			d.put("negative", EMObject::INT, "If set, returns -1 * ccc product. Set by default so smaller is better");
			return d;
		}
		
		static const string NAME;
	};

	
	/** This is a FSC comparitor for tomography. I only counts voxels above a threshold, which is obtained from subtomogram metadata
	* @ingroup CUDA_ENABLED
	* @author John Flanagan
	* @date 2012-2-02
	* @param sigmas The number of times the standard deviation of Fourier amplitudes to accept
	* @param minres The minimum resolution to accept
	* @param maxes The maximum resolution to accept
	* @param apix The angstroms per pixel to use
	*/
	class TomoFscCmp:public Cmp
	{
	  public:
		virtual float cmp(EMData * image, EMData * with) const;

		virtual string get_name() const 
		{
			return NAME;
		}

		virtual string get_desc() const
		{
			return "EXPERIMENTAL - Fsc with consideration given for the missing wedge. Not ready for routine use.";
		}

		static Cmp *NEW()
		{
			return new TomoFscCmp();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("normalize", EMObject::EMDATA,"Return the negative value (which is EMAN2 convention), Defalut is true(1)");
			d.put("sigmas", EMObject::FLOAT, "The number of times the standard deviation of Fourier amplitudes to accept");
			d.put("minres", EMObject::FLOAT, "The minimum resolution to accept (1/A) Default is inf");
			d.put("maxres", EMObject::FLOAT, "The maximum resolution to accept (1/A) Default=0.0");
			d.put("apix", EMObject::FLOAT, "The angstroms per pixel to use. Default = apix_x(1.0 if not present)");
			return d;
		}
		
		static const string NAME;
	};
	
	/** This will calculate the dot product for each quadrant of the image and
	* return the worst value
	* @author Steve Ludtke (sludtke@bcm.tmc.edu)
	* @date 2005-07-13
	* @param negative Returns -1 * dot product, default true
	* @param normalize Returns normalized dot product -1.0 - 1.0
	*/
	class QuadMinDotCmp:public Cmp
	{
	  public:
		float cmp(EMData * image, EMData * with) const;

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Calculates dot product for each quadrant and returns worst value (default -1 * dot product)";
		}

		static Cmp *NEW()
		{
			return new QuadMinDotCmp();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("negative", EMObject::INT, "If set, returns -1 * dot product. Default = true (smaller is better)");
			d.put("normalize", EMObject::INT, "If set, returns normalized dot product -1.0 - 1.0.");
			return d;
		}
		
		static const string NAME;
	};

	/** Uses math.sub.optimal to remove the density of the reference from the image as much as possible.
	 * What is left should be pure noise if there is a good match. If the match is not as good there should
	 * be additional residual power. This comparator returns the average residual power after subtraction
	 * weighted in various ways
	*/
	class OptSubCmp:public Cmp
	{
	  public:
		OptSubCmp() {}

		float cmp(EMData * image, EMData * with) const;

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Residual power left in the image after optimally subtracting the reference with math.sub.optimal. Smaller values indicate a closer match.";
		}

		static Cmp *NEW()
		{
			return new OptSubCmp();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("minres", EMObject::FLOAT, "Lowest resolution to use in comparison (soft cutoff). Requires accurate A/pix in image. <0 disables. Default=200");
			d.put("maxres", EMObject::FLOAT, "Highest resolution to use in comparison (soft cutoff). Requires accurate A/pix in image. <0 disables.  Default=10");
			d.put("zeromask", EMObject::BOOL, "Treat zero regions in 'with' as a mask");
			d.put("ctfweight", EMObject::BOOL, "Weight frequencies using CTF");
//			d.put("swap", EMObject::BOOL, "Swaps the identity of this and with, since A->B comparisons are not equivalent to B->A");
			d.put("mask", EMObject::EMDATA, "Real space mask. Only computes the residual power under the mask. Significant speed penalty if specified. Default=None.");
			return d;
		}
		
		static const string NAME;

	};

	
	/** Variance between two data sets after various modifications.
	* Generally, 'this' should be noisy and 'with' should be less noisy.
	* linear scaling (mx + b) of the densities in 'this' is performed
	* to produce the smallest possible variance between images.
	*
	* If keepzero is set, then zero pixels are left at zero (b is not added).
	* if matchfilt is set, then 'with' is filtered so its radial power spectrum matches 'this'
	* If invert is set, then (y-b)/m is applied to the second image rather than mx+b to the first.
	*
	* To modify 'this' to match the operation performed here, scale
	* should be applied first, then b should be added
	*/
	class OptVarianceCmp:public Cmp
	{
	  public:
		OptVarianceCmp() : scale(0), shift(0) {}

		float cmp(EMData * image, EMData * with) const;

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Real-space variance after density optimization, self should be noisy and target less noisy. Linear transform applied to density to minimize variance.";
		}

		static Cmp *NEW()
		{
			return new OptVarianceCmp();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("invert", EMObject::INT, "If set, 'with' is rescaled rather than 'this'. 'this' should still be the noisier image. (default=0)");
			d.put("keepzero", EMObject::INT, "If set, zero pixels will not be adjusted in the linear density optimization. (default=1)");
			d.put("matchfilt", EMObject::INT, "If set, with will be filtered so its radial power spectrum matches 'this' before density optimization of this. (default=1)");
			d.put("matchamp", EMObject::INT, "Takes per-pixel Fourier amplitudes from self and imposes them on the target, but leaves the phases alone. (default=0)");
			d.put("radweight", EMObject::INT, "Upweight variances closer to the edge of the image. (default=0)");
			d.put("debug", EMObject::INT, "Performs various debugging actions if set.");
			return d;
		}

		float get_scale() const
		{
			return scale;
		}

		float get_shift() const
		{
			return shift;
		}
		
		static const string NAME;

	private:
		mutable float scale;
		mutable float shift;
	};
	
	/** Amplitude weighted mean phase difference (radians) with optional
     * SNR weight. SNR should be an array as returned by ctfcurve()
     * 'with' should be the less noisy image, since it's amplitudes
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
		float cmp(EMData * image, EMData * with) const;

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Mean phase difference";
		}

		static Cmp *NEW()
		{
			return new PhaseCmp();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("snrweight", EMObject::INT, "If set, the SNR of 'this' will be used to weight the result. If 'this' lacks CTF info, it will check 'with'. (default=0)");
			d.put("snrfn", EMObject::INT, "If nonzero, an empirical function will be used as a radial weight rather than the true SNR. (1 - exp decay)'. (default=0)");
			d.put("ampweight", EMObject::INT, "If set, the amplitude of 'with' will be used as a weight in the averaging'. (default=0)");
			d.put("zeromask", EMObject::INT, "Treat regions in either image that are zero as a mask");
			d.put("minres", EMObject::FLOAT, "Lowest resolution to use in comparison (soft cutoff). Requires accurate A/pix in image. <0 disables. Default=500");
			d.put("maxres", EMObject::FLOAT, "Highest resolution to use in comparison (soft cutoff). Requires accurate A/pix in image. <0 disables.  Default=10");
			return d;
		}
		
		static const string NAME;

//#ifdef EMAN2_USING_CUDA
//		 float cuda_cmp(EMData * image, EMData *with) const;
//#endif //EMAN2_USING_CUDA
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
		float cmp(EMData * image, EMData * with) const;

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Computes the mean Fourier Ring Correlation between the image and reference (with optional weighting factors).";
		}

		static Cmp *NEW()
		{
			return new FRCCmp();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("snrweight", EMObject::INT, "If set, the SNR of 'this' will be used to weight the result. If 'this' lacks CTF info, it will check 'with'. (default=0)");
			d.put("ampweight", EMObject::INT, "If set, the amplitude of 'this' will be used to weight the result (default=0)");
			d.put("sweight", EMObject::INT, "If set, weight the (1-D) average by the number of pixels in each ring (default=1)");
			d.put("nweight", EMObject::INT, "Downweight similarity based on number of particles in reference (default=0)");
			d.put("zeromask", EMObject::INT, "Treat regions in either image that are zero as a mask");
			d.put("minres", EMObject::FLOAT, "Lowest resolution to use in comparison (soft cutoff). Requires accurate A/pix in image. <0 disables. Default=500");
			d.put("maxres", EMObject::FLOAT, "Highest resolution to use in comparison (soft cutoff). Requires accurate A/pix in image. <0 disables.  Default=10");
			return d;
		}
		
		static const string NAME;
	};
	
	
	
	/** 
    */
	class VerticalCmp:public Cmp
	{
	  public:
		float cmp(EMData * image, EMData * with) const;

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Align the 2d image vertically. Reference image is ignored";
		}

		static Cmp *NEW()
		{
			return new VerticalCmp();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			return d;
		}
		
		static const string NAME;
	};
	
		
	

	template <> Factory < Cmp >::Factory();

	void dump_cmps();
	map<string, vector<string> > dump_cmps_list();
}


#endif

/* vim: set ts=4 noet: */
