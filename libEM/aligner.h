/**
 * $Id$
 */
#ifndef eman__aligner_h__
#define eman__aligner_h__ 1


#include "emobject.h"


namespace EMAN
{
	class EMData;
	class Cmp;

	/** Aligner class defines image alignment method. It aligns 2
	 * images based on a user-given comparison method. 
	 * Aligner class is the base class for all aligners. Each
     * specific Aligner class has a unique name. This name is used to
     * create a new Aligner instance or call an Aligner.
     *
	 * All Aligner classes in EMAN are managed by a Factory
	 * pattern. So each Aligner class must define:
	 *   - a unique name to idenfity itself in the factory.
	 *   - a static method to register itself in the factory.
	 *
	 *
     * Typical usage of Aligners:
     *
     *  - How to get all the Aligner types
     @code
     *    vector<string> all_aligners = Factory<Aligner>::get_list();
     @endcode
	 *
     *  - How to use an Aligner
     @code
     *    EMData *image1 = ...;
     *    EMData *image2 = ...;
     *    image1->align("ALIGNER_NAME", image2);
     @endcode
	 
     *  - How to define a new Aligner class \n
     *    A new XYZAligner class should implement the following functions:
     *    (Please replace 'XYZ' with your own class name).
	 @code
     *        EMData *align(EMData * this_img, EMData * to_img, string cmp_name = "") const;
     *        TypeDict get_param_types() const;
     *        string get_name() const { return "XYZ"; }
     *        static Aligner* NEW() { return new XYZAligner(); }
	 @endcode
     */
	class Aligner
	{
	  public:
		virtual ~ Aligner()
		{
		}
		
		virtual EMData *align(EMData * this_img, EMData * to_img) const = 0;
		
		/** To align 'this_img' with another image passed in through
		 * its parameters. The alignment uses a user-given comparison
		 * method to compare the two images. If none is given, a default
		 * one is used.
		 *
		 * @param this_img The image to be compared.
		 * @param to_img 'this_img" is aligned with 'to_img'.
		 * @param cmp_name The comparison method to compare the two images.
		 * @param cmp_params The parameter dictionary for comparison method.
		 * @return The aligned image.
		 */
		virtual EMData *align(EMData * this_img, EMData * to_img, 
							  const string & cmp_name, const Dict& cmp_params) const = 0;
		
		/** Get the Aligner's name. Each Aligner is identified by a unique name.
		 * @return The Aligner's name.
		 */
		virtual string get_name() const = 0;

		virtual string get_desc() const = 0;
		/** Get the Aligner parameters in a key/value dictionary.
		 * @return A key/value pair dictionary containing the parameters.
		 */
		virtual Dict get_params() const
		{
			return params;
		}

		/** Set the Aligner parameters using a key/value dictionary.
		 * @param new_params A dictionary containing the new parameters.
		 */
		virtual void set_params(const Dict & new_params)
		{
			params = new_params;
		}

		virtual TypeDict get_param_types() const = 0;

	  protected:
		mutable Dict params;
	};

	/** Translational 2D Alignment using cross correlation.
     * It calculates the shift for a translational alignment, then
     * do the translation.
     */
	class TranslationalAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img, 
						const string & cmp_name, const Dict& cmp_params) const;

		EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "", Dict());
		}
		
		string get_name() const
		{
			return "translational";
		}
		
		string get_desc() const
		{
			return "Translational 2D alignment by cross-correlation";
		}

		static Aligner *NEW()
		{
			return new TranslationalAligner();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("intonly", EMObject::INT);
			d.put("maxshift", EMObject::INT);
			return d;
		}
	};

	/** Translational 3D Alignment using  cross correlation
     * It calculates the shift for a translational alignment, then
     * do the translation.
     */
	class Translational3DAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img, 
						const string & cmp_name, const Dict& cmp_params) const;

		EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "", Dict());
		}
		
		string get_name() const
		{
			return "translational3d";
		}
		
		string get_desc() const
		{
			return "Translational 3D alignment by cross-correlation";
		}

		static Aligner *NEW()
		{
			return new Translational3DAligner();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("intonly", EMObject::INT);
			return d;
		}

	};

	/** rotational alignment using angular correlation
     */
	class RotationalAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img, 
						const string & cmp_name = "", const Dict& cmp_params = Dict()) const;

		EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "", Dict());
		}
		
		string get_name() const
		{
			return "rotational";
		}
		
		string get_desc() const
		{
			return "Rotational 2D alignment by angular correlation using make_rotational_footprint(), ~translationally independent";
		}

		static Aligner *NEW()
		{
			return new RotationalAligner();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			return d;
		}
	};

	/** rotational alignment assuming centers are correct
     */
	class RotatePrecenterAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img, 
						const string & cmp_name = "", const Dict& cmp_params = Dict()) const;

		EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "", Dict());
		}
		
		string get_name() const
		{
			return "rotate_precenter";
		}
		
		string get_desc() const
		{
			return "Rotational 2D alignment by angular correlation using unwrap(4,ny*7/16,...), precentered images only";
		}

		static Aligner *NEW()
		{
			return new RotatePrecenterAligner();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			return d;
		}
	};

	/** rotational alignment via circular harmonic
     */
	class RotateCHAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img, 
						const string & cmp_name, const Dict& cmp_params = Dict()) const;

		EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "", Dict());
		}
		
		string get_name() const
		{
			return "rotate_ch";
		}
		
		string get_desc() const
		{
			return "Rotational 2D alignment using circular harmonics";
		}

		static Aligner *NEW()
		{
			return new RotateCHAligner();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			
			d.put("irad", EMObject::INT);
			d.put("orad", EMObject::INT);
			return d;
		}
	};

	/** rotational, translational alignment
     */
	class RotateTranslateAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img, 
						const string & cmp_name, const Dict& cmp_params) const;

		EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "SqEuclidean", Dict());
		}
		
		string get_name() const
		{
			return "rotate_translate";
		}
		
		string get_desc() const
		{
			return "Full 2D alignment using 'Rotational' and 'Translational'.";
		}

		static Aligner *NEW()
		{
			return new RotateTranslateAligner();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			//d.put("usedot", EMObject::INT);
			d.put("maxshift", EMObject::INT);
			return d;
		}
	};

	/** rotational, translational alignment
     */
	class RotateTranslateBestAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img, 
						const string & cmp_name, const Dict& cmp_params) const;

		EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "frc", Dict());
		}

		string get_name() const
		{
			return "rotate_translate_best";
		}
		
		string get_desc() const
		{
			return "Full 2D alignment using 'Rotational' and 'Translational', also incorporates 2D 'Refine' alignments.";
		}

		static Aligner *NEW()
		{
			return new RotateTranslateBestAligner();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;			
			d.put("maxshift", EMObject::INT);
			d.put("snr", EMObject::FLOATARRAY);
			return d;
		}
	};

	/** rotational, translational alignment with Radon transforms
     */
	class RotateTranslateRadonAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img, 
						const string & cmp_name, const Dict& cmp_params) const;
		EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "", Dict());
		}
		string get_name() const
		{
			return "rotate_translate_radon";
		}

		string get_desc() const
		{
			return "Experimental 2D alignment using Radon transforms.";
		}

		static Aligner *NEW()
		{
			return new RotateTranslateRadonAligner();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;

			d.put("maxshift", EMObject::INT);
			d.put("radonwith", EMObject::EMDATA);
			d.put("radonthis", EMObject::EMDATA);
			return d;
		}
	};

	/** rotational and flip alignment
     */
	class RotateFlipAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img, 
						const string & cmp_name, const Dict& cmp_params) const;
		EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "", Dict());
		}
		string get_name() const
		{
			return "rotate_flip";
		}

		string get_desc() const
		{
			return "'Rotational' alignment with possible handedness change";
		}

		static Aligner *NEW()
		{
			return new RotateFlipAligner();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			
			d.put("flip", EMObject::EMDATA);
			d.put("imask", EMObject::INT);
			return d;
		}
	};

	/** rotational, translational and flip alignment
     */
	class RotateTranslateFlipAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img, 
						const string & cmp_name, const Dict& cmp_params) const;
		EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "SqEuclidean", Dict());
		}
		
		string get_name() const
		{
			return "rotate_translate_flip";
		}

		string get_desc() const
		{
			return "'RotatateTranslate' alignment with possible handedness change";
		}

		static Aligner *NEW()
		{
			return new RotateTranslateFlipAligner();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			
			d.put("flip", EMObject::EMDATA);
			d.put("usedot", EMObject::INT);
			d.put("maxshift", EMObject::INT);
			return d;
		}
	};

	/** rotational, translational and flip alignment using real-space methods. slow
    */
	class RTFSlowAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img, 
						const string & cmp_name, const Dict& cmp_params) const;
		EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "SqEuclidean", Dict());
		}
		
		string get_name() const
		{
			return "rtf_slow";
		}

		string get_desc() const
		{
			return "Experimental full 2D alignment with handedness check using semi-exhaustive search (not necessarily better than RTFBest)";
		}
		
		static Aligner *NEW()
		{
			return new RTFSlowAligner();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;

			d.put("flip", EMObject::EMDATA);
			d.put("maxshift", EMObject::INT);
			return d;
		}
	};
	/** rotational, translational and flip alignment using exhaustive search. VERY SLOW
     */
	class RTFSlowestAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img, 
						const string & cmp_name, const Dict& cmp_params) const;
		EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "SqEuclidean", Dict());
		}
		string get_name() const
		{
			return "rtf_slowest";
		}

		string get_desc() const
		{
			return "Experimental full 2D alignment with handedness check using more exhaustive search (not necessarily better than RTFBest)";
		}
		
		static Aligner *NEW()
		{
			return new RTFSlowestAligner();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			
			d.put("flip", EMObject::EMDATA);
			d.put("maxshift", EMObject::INT);
			return d;
		}
	};

	/** rotational, translational and flip alignment using fscmp at multiple locations, slow
     * but this routine probably produces the best results
     */
	class RTFBestAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img, 
						const string & cmp_name, const Dict& cmp_params) const;

		EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "SqEuclidean", Dict());
		}

		string get_name() const
		{
			return "rtf_best";
		}

		string get_desc() const
		{
			return "Best full 2D alignment with handedness check currently in EMAN.";
		}
		
		static Aligner *NEW()
		{
			return new RTFBestAligner();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			
			d.put("flip", EMObject::EMDATA);
			d.put("maxshift", EMObject::INT);
			d.put("snr", EMObject::FLOATARRAY);
			return d;
		}
	};

	/** rotational, translational and flip alignment with Radon transforms
     */
	class RTFRadonAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img, 
						const string & cmp_name, const Dict& cmp_params) const;
		EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "SqEuclidean", Dict());
		}
		string get_name() const
		{
			return "rtf_radon";
		}

		string get_desc() const
		{
			return "Experimental 2D alignment using Radon transforms with handedness check.";
		}

		static Aligner *NEW()
		{
			return new RTFRadonAligner();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			
			d.put("maxshift", EMObject::INT);
			d.put("thisf", EMObject::EMDATA);
			d.put("radonwith", EMObject::EMDATA);
			d.put("radonthis", EMObject::EMDATA);
			d.put("radonthisf", EMObject::EMDATA);
			return d;
		}
	};

	/** refine alignment
     */
	class RefineAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img, 
						const string & cmp_name, const Dict& cmp_params) const;

		EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "SqEuclidean", Dict());
		}
		
		string get_name() const
		{
			return "refine";
		}

		string get_desc() const
		{
			return "Refines a preliminary 2D alignment using a simplex algorithm. Subpixel accuracy.";
		}

		static Aligner *NEW()
		{
			return new RefineAligner();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;

			d.put("mode", EMObject::INT);
			d.put("snr", EMObject::FLOATARRAY);
			d.put("alot", EMObject::FLOAT);
			d.put("az", EMObject::FLOAT);
			d.put("phi", EMObject::FLOAT);
			d.put("dx", EMObject::FLOAT);
			d.put("dy", EMObject::FLOAT);
			d.put("dz", EMObject::FLOAT);
			d.put("alt", EMObject::FLOAT);

			return d;
		}
	};

	template <> Factory < Aligner >::Factory();

	void dump_aligners();
	map<string, vector<string> > dump_aligners_list();
}

#endif
