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

//		/** Get a Transform pointer that is currently stored in the image header corresponding to the given key.
//		 * If non existant then initialize a new Transform pointer, set it as the image attribute, and return it.
//		 * The point being that if you have to access the xform.align2d or the xform.align3d
//		 * attributes from multiple aligners then you always want to get the same one.
//		 * @param key the alignment key such as "xform.align2d" or "xform.align3d"
//		 * @param image the image from which the Transform pointer will be extracted (and potentially assigned to if non existant)
//		 * @return the Transform pointer that is currently stored in the image header corresponding to the given key.
//		 */
//		static Transform* get_align_attr(const string& key, EMData* const image );
//
//		static Transform* get_set_align_attr(const string& key, EMData* const to_image, const EMData* const from_image  );
	};

	/** Translational 2D Alignment using cross correlation.
     * It calculates the shift for a translational alignment, then
     * do the translation.
	 * @ingroup CUDA_ENABLED
	 * @param intonly Integer pixel translations only
	 * @param maxshift Maximum translation in pixels
	 * @param nozero Zero translation not permitted (useful for CCD images)
     */
	class TranslationalAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img,
						const string & cmp_name="dot", const Dict& cmp_params = Dict()) const;

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
			return "Translational 2D and 3D alignment by cross-correlation";
		}

		static Aligner *NEW()
		{
			return new TranslationalAligner();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("intonly", EMObject::INT,"Integer pixel translations only");
			d.put("maxshift", EMObject::INT,"Maximum translation in pixels");
			d.put("nozero", EMObject::INT,"Zero translation not permitted (useful for CCD images)");
			return d;
		}
	};

	/** rotational alignment using angular correlation
	 * @param rfp_mode Either 0,1 or 2. A temporary flag for testing the rotational foot print. O is the original eman1 way. 1 is just using calc_ccf without padding. 2 is using calc_mutual_correlation without padding
	 * @ingroup CUDA_ENABLED
     */
	class RotationalAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img,
						const string & cmp_name = "dot", const Dict& cmp_params = Dict()) const;

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
			return "Performs rotational alignment,works accurately if the image is precentered, normally called internally in combination with translational and flip alignment";
		}

		static Aligner *NEW()
		{
			return new RotationalAligner();
		}

		static EMData * align_180_ambiguous(EMData * this_img, EMData * to_img, int rfp_mode = 0);

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("rfp_mode", EMObject::INT,"Either 0,1 or 2. A temporary flag for testing the rotational foot print. O is the original eman1 way. 1 is just using calc_ccf without padding. 2 is using calc_mutual_correlation without padding.");
			return d;
		}
	};

	/** rotational alignment assuming centers are correct
     */
	class RotatePrecenterAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img,
						const string & cmp_name = "dot", const Dict& cmp_params = Dict()) const;

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
			return "Performs rotational alignment and works accurately if the image is precentered";
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

	/** rotational, translational alignment
	 * @param maxshift Maximum translation in pixels
	 * @param nozero Zero translation not permitted (useful for CCD images)
	 * @param rfp_mode Either 0,1 or 2. A temporary flag for testing the rotational foot print
	 * @ingroup CUDA_ENABLED
     */
	class RotateTranslateAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img,
					   const string & cmp_name="dot", const Dict& cmp_params = Dict()) const;

		EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "sqeuclidean", Dict());
		}

		string get_name() const
		{
			return "rotate_translate";
		}

		string get_desc() const
		{
			return "Performs rotational alignment and follows this with translational alignment.";
		}

		static Aligner *NEW()
		{
			return new RotateTranslateAligner();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			//d.put("usedot", EMObject::INT);
			d.put("maxshift", EMObject::INT, "Maximum translation in pixels");
			d.put("nozero", EMObject::INT,"Zero translation not permitted (useful for CCD images)");
			d.put("rfp_mode", EMObject::INT,"Either 0,1 or 2. A temporary flag for testing the rotational foot print");
			return d;
		}
	};

	/** rotational, translational alignment
	 * @param maxshift Maximum translation in pixels
	 * @param snr signal to noise ratio array
     */
	class RotateTranslateBestAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img,
					   const string & cmp_name="dot", const Dict& cmp_params = Dict()) const;

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
			d.put("maxshift", EMObject::INT, "Maximum translation in pixels");
			d.put("snr", EMObject::FLOATARRAY, "signal to noise ratio array");
			return d;
		}


	};

	/** rotational and flip alignment
	 * @param imask
	 * @param rfp_mode Either 0,1 or 2. A temporary flag for testing the rotational foot print
	 * @ingroup CUDA_ENABLED
     */
	class RotateFlipAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img,
					   const string & cmp_name="dot", const Dict& cmp_params = Dict()) const;
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
			return "Performs two rotational alignments, one using the original image and one using the hand-flipped image. Decides which alignment is better using a comparitor and returns it";
		}

		static Aligner *NEW()
		{
			return new RotateFlipAligner();
		}

		TypeDict get_param_types() const
		{
			return static_get_param_types();
		}

		static TypeDict static_get_param_types() {
			TypeDict d;

			d.put("imask", EMObject::INT);
			d.put("rfp_mode", EMObject::INT,"Either 0,1 or 2. A temporary flag for testing the rotational foot print");
			return d;
		}

	};

	/** rotational, translational and flip alignment
	 * @param flip
	 * @param usedot
	 * @param maxshift Maximum translation in pixels
	 * @param rfp_mode Either 0,1 or 2. A temporary flag for testing the rotational foot print
	 * @ingroup CUDA_ENABLED
     */
	class RotateTranslateFlipAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img,
					   const string & cmp_name="dot", const Dict& cmp_params = Dict()) const;
		EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "sqeuclidean", Dict());
		}

		string get_name() const
		{
			return "rotate_translate_flip";
		}

		string get_desc() const
		{
			return " Does two 'rotate_translate' alignments, one to accommodate for possible handedness change. Decided which alignment is better using a comparitor and returns the aligned image as the solution";
		}

		static Aligner *NEW()
		{
			return new RotateTranslateFlipAligner();
		}

		TypeDict get_param_types() const
		{
			return static_get_param_types();
		}

		static TypeDict static_get_param_types() {
			TypeDict d;

			d.put("flip", EMObject::EMDATA);
			d.put("usedot", EMObject::INT);
			d.put("maxshift", EMObject::INT, "Maximum translation in pixels");
			d.put("rfp_mode", EMObject::INT,"Either 0,1 or 2. A temporary flag for testing the rotational foot print");
			return d;
		}
	};

	/** rotational, translational and flip alignment using real-space methods. slow
	 * @param flip
	 * @param maxshift Maximum translation in pixels
	 * */
	class RTFExhaustiveAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img,
					   const string & cmp_name="dot", const Dict& cmp_params = Dict()) const;
		EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "sqeuclidean", Dict());
		}

		string get_name() const
		{
			return "rtf_exhaustive";
		}

		string get_desc() const
		{
			return "Experimental full 2D alignment with handedness check using semi-exhaustive search (not necessarily better than RTFBest)";
		}

		static Aligner *NEW()
		{
			return new RTFExhaustiveAligner();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;

			d.put("flip", EMObject::EMDATA);
			d.put("maxshift", EMObject::INT, "Maximum translation in pixels");
			return d;
		}
	};

	/** rotational, translational and flip alignment using exhaustive search. This is very slow
	 * but can ensure localization of the global maximum
	 * @param flip Optional. This is the flipped version of the images that is being aligned. If specified it will be used for the handedness check, if not a flipped copy of the image will be made
	 * @param maxshift The maximum length of the detectable translational shift
	 * @param transtep The translation step to take when honing the alignment, which occurs after coarse alignment
	 * @param angstep The angular step (in degrees) to take in the exhaustive search for the solution angle. Typically very small i.e. 3 or smaller
     */
	class RTFSlowExhaustiveAligner:public Aligner
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
			return "rtf_slow_exhaustive";
		}

		string get_desc() const
		{
			return "Experimental full 2D alignment with handedness check using more exhaustive search (not necessarily better than RTFBest)";
		}

		static Aligner *NEW()
		{
			return new RTFSlowExhaustiveAligner();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;

			d.put("flip", EMObject::EMDATA,"Optional. This is the flipped version of the images that is being aligned. If specified it will be used for the handedness check, if not a flipped copy of the image will be made");
			d.put("maxshift", EMObject::INT,"The maximum length of the detectable translational shift");
			d.put("transtep", EMObject::FLOAT,"The translation step to take when honing the alignment, which occurs after coarse alignment");
			d.put("angstep", EMObject::FLOAT,"The angular step (in degrees) to take in the exhaustive search for the solution angle. Typically very small i.e. 3 or smaller.");
			return d;
		}
	};

	/** refine alignment. Refines a preliminary 2D alignment using a simplex algorithm. Subpixel precision.
     */
	class RefineAligner:public Aligner
	{
	  public:
		EMData * align(EMData * this_img, EMData * to_img,
					   const string & cmp_name="dot", const Dict& cmp_params = Dict()) const;

		EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "sqeuclidean", Dict());
		}

		string get_name() const
		{
			return "refine";
		}

		string get_desc() const
		{
			return "Refines a preliminary 2D alignment using a simplex algorithm. Subpixel precision.";
		}

		static Aligner *NEW()
		{
			return new RefineAligner();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;

			d.put("mode", EMObject::INT);
			d.put("xform.align2d", EMObject::TRANSFORM);
			d.put("stepx", EMObject::FLOAT);
			d.put("stepy", EMObject::FLOAT);
			d.put("stepaz", EMObject::FLOAT);
			d.put("precision", EMObject::FLOAT);
			d.put("maxiter", EMObject::INT,"The maximum number of iterations that can be performed by the Simplex minimizer");
			d.put("maxshift", EMObject::INT,"Maximum translation in pixels in any direction");
			return d;
		}
	};


	/** refine alignment. Refines a preliminary 3D alignment using a simplex algorithm. Subpixel precision.
	 */
	class Refine3DAligner:public Aligner
	{
		public:
			EMData * align(EMData * this_img, EMData * to_img,
						   const string & cmp_name="dot", const Dict& cmp_params = Dict()) const;

			EMData * align(EMData * this_img, EMData * to_img) const
			{
				return align(this_img, to_img, "sqeuclidean", Dict());
			}

			string get_name() const
			{
				return "refine3d";
			}

			string get_desc() const
			{
				return "Refines a preliminary 3D alignment using a simplex algorithm. Subpixel precision.";
			}

			static Aligner *NEW()
			{
				return new Refine3DAligner();
			}

			TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("xform.align3d", EMObject::TRANSFORM);
				d.put("stepx", EMObject::FLOAT);
				d.put("stepy", EMObject::FLOAT);
				d.put("stepz", EMObject::FLOAT);
				d.put("stepphi", EMObject::FLOAT);
				d.put("stepalt", EMObject::FLOAT);
				d.put("stepaz", EMObject::FLOAT);
				d.put("precision", EMObject::FLOAT);
				d.put("maxiter", EMObject::INT, "The maximum number of iterations that can be performed by the Simplex minimizer");
				d.put("maxshift", EMObject::INT,"Maximum translation in pixels in any direction");
				return d;
			}
	};

	
	class CUDA_Aligner
	{
	  public:
	  	CUDA_Aligner();

		~CUDA_Aligner();

		void setup(int nima, int nx, int ny, int ring_length, int nring, float step, int kx, int ky);

		void insert_image(EMData *image, int num);

		vector<float> alignment_2d(EMData *ref_image);

	  private:
	        float *image_stack;
		float *ccf;
		int NIMA, NX, NY, RING_LENGTH, NRING, KX, KY;
		float STEP;
	};


	template <> Factory < Aligner >::Factory();

	void dump_aligners();
	map<string, vector<string> > dump_aligners_list();
}

#endif
