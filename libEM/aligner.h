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

		/** This function first added in the context of the 3D aligners used by e2tomohunter:
		 * which wants the n best solutions, as opposed to just the best. Return value is an
		 * ordered vector of Dicts of length nsoln. The data with idx 0 has the best solution in
		 * it. The Dicts in the vector have two keys, "score" (which is a floating point score,
		 * probably correlation score), and "xform.align3d", which is a Transform containing
		 * the alignment parameters
		 * @param this_img the image that will be  aligned (transformed) and compared to to_img as part of the process of finding the best alignment
		 * @param to_img the image that will stay still as this_img is transformed and compared to it
		 * @param nsoln the number of solutions you want to receive in the return vector.
		 * @param cmp_name the name of a comparator - may be unused
		 * @param cmp_params the params of the comparator - may be unused
		 * @return an ordered vector of Dicts of length nsoln. The Dicts in the vector have keys "score" (i.e. correlation score) and "xform.align3d" (Transform containing the alignment)
		 */
		virtual vector<Dict> xform_align_nbest(EMData * this_img, EMData * to_img, const unsigned int nsoln, const string & cmp_name, const Dict& cmp_params) const {
			vector<Dict> solns;
			return solns;
		}

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
		virtual EMData * align(EMData * this_img, EMData * to_img,
						const string & cmp_name="dot", const Dict& cmp_params = Dict()) const;

		virtual EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "", Dict());
		}

		virtual string get_name() const
		{
			return NAME;
		}

		virtual string get_desc() const
		{
			return "Translational 2D and 3D alignment by cross-correlation";
		}

		static Aligner *NEW()
		{
			return new TranslationalAligner();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("intonly", EMObject::INT,"Integer pixel translations only");
			d.put("useflcf", EMObject::INT,"Use Fast Local Correlation Function rather than CCF");
			d.put("maxshift", EMObject::INT,"Maximum translation in pixels");
			d.put("masked", EMObject::INT,"Treat zero pixels in 'this' as a mask for normalization (default false)");
			d.put("nozero", EMObject::INT,"Zero translation not permitted (useful for CCD images)");
			return d;
		}
		
		static const string NAME;
	};

	/** rotational alignment using angular correlation
	 * @param rfp_mode Either 0,1 or 2. A temporary flag for testing the rotational foot print. O is the original eman1 way. 1 is just using calc_ccf without padding. 2 is using calc_mutual_correlation without padding
	 * @ingroup CUDA_ENABLED
     */
	class RotationalAligner:public Aligner
	{
	  public:
		virtual EMData * align(EMData * this_img, EMData * to_img,
						const string & cmp_name = "dot", const Dict& cmp_params = Dict()) const;

		virtual EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "", Dict());
		}

		virtual string get_name() const
		{
			return NAME;
		}

		virtual string get_desc() const
		{
			return "Performs rotational alignment,works accurately if the image is precentered, normally called internally in combination with translational and flip alignment";
		}

		static Aligner *NEW()
		{
			return new RotationalAligner();
		}

		static EMData * align_180_ambiguous(EMData * this_img, EMData * to_img, int rfp_mode = 0);

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("rfp_mode", EMObject::INT,"Either 0,1 or 2. A temporary flag for testing the rotational foot print. O is the original eman1 way. 1 is just using calc_ccf without padding. 2 is using calc_mutual_correlation without padding.");
			return d;
		}
		
		static const string NAME;
	};

	/** rotational alignment assuming centers are correct
     */
	class RotatePrecenterAligner:public Aligner
	{
	  public:
		virtual EMData * align(EMData * this_img, EMData * to_img,
						const string & cmp_name = "dot", const Dict& cmp_params = Dict()) const;

		virtual EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "", Dict());
		}

		virtual string get_name() const
		{
			return NAME;
		}

		virtual string get_desc() const
		{
			return "Performs rotational alignment and works accurately if the image is precentered";
		}

		static Aligner *NEW()
		{
			return new RotatePrecenterAligner();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			return d;
		}
		
		static const string NAME;
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
		  virtual EMData * align(EMData * this_img, EMData * to_img,
					   const string & cmp_name="dot", const Dict& cmp_params = Dict()) const;

		virtual EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "sqeuclidean", Dict());
		}

		virtual string get_name() const
		{
			return NAME;
		}

		virtual string get_desc() const
		{
			return "Performs rotational alignment and follows this with translational alignment.";
		}

		static Aligner *NEW()
		{
			return new RotateTranslateAligner();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			//d.put("usedot", EMObject::INT);
			d.put("maxshift", EMObject::INT, "Maximum translation in pixels");
			d.put("nozero", EMObject::INT,"Zero translation not permitted (useful for CCD images)");
			d.put("rfp_mode", EMObject::INT,"Either 0,1 or 2. A temporary flag for testing the rotational foot print");
			d.put("useflcf", EMObject::INT,"Use Fast Local Correlation Function rather than CCF for translational alignment");
			return d;
		}
		
		static const string NAME;
	};

	/** rotational, translational alignment
	 * @param maxshift Maximum translation in pixels
	 * @param snr signal to noise ratio array
     */
	class RotateTranslateBestAligner:public Aligner
	{
	  public:
		virtual EMData * align(EMData * this_img, EMData * to_img,
					   const string & cmp_name="dot", const Dict& cmp_params = Dict()) const;

		virtual EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "frc", Dict());
		}

		virtual string get_name() const
		{
			return NAME;
		}

		virtual string get_desc() const
		{
			return "Full 2D alignment using 'Rotational' and 'Translational', also incorporates 2D 'Refine' alignments.";
		}

		static Aligner *NEW()
		{
			return new RotateTranslateBestAligner();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("maxshift", EMObject::INT, "Maximum translation in pixels");
			d.put("snr", EMObject::FLOATARRAY, "signal to noise ratio array");
			return d;
		}

		static const string NAME;
	};

	/** rotational and flip alignment
	 * @param imask
	 * @param rfp_mode Either 0,1 or 2. A temporary flag for testing the rotational foot print
	 * @ingroup CUDA_ENABLED
     */
	class RotateFlipAligner:public Aligner
	{
	  public:
		virtual EMData * align(EMData * this_img, EMData * to_img,
					   const string & cmp_name="dot", const Dict& cmp_params = Dict()) const;
		virtual EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "", Dict());
		}
		virtual string get_name() const
		{
			return NAME;
		}

		virtual string get_desc() const
		{
			return "Performs two rotational alignments, one using the original image and one using the hand-flipped image. Decides which alignment is better using a comparitor and returns it";
		}

		static Aligner *NEW()
		{
			return new RotateFlipAligner();
		}

		virtual TypeDict get_param_types() const
		{
			return static_get_param_types();
		}

		static TypeDict static_get_param_types() {
			TypeDict d;

			d.put("imask", EMObject::INT);
			d.put("rfp_mode", EMObject::INT,"Either 0,1 or 2. A temporary flag for testing the rotational foot print");
			return d;
		}

		static const string NAME;
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
		virtual EMData * align(EMData * this_img, EMData * to_img,
					   const string & cmp_name="dot", const Dict& cmp_params = Dict()) const;
		virtual EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "sqeuclidean", Dict());
		}

		virtual string get_name() const
		{
			return NAME;
		}

		virtual string get_desc() const
		{
			return " Does two 'rotate_translate' alignments, one to accommodate for possible handedness change. Decided which alignment is better using a comparitor and returns the aligned image as the solution";
		}

		static Aligner *NEW()
		{
			return new RotateTranslateFlipAligner();
		}

		virtual TypeDict get_param_types() const
		{
			return static_get_param_types();
		}

		static TypeDict static_get_param_types() {
			TypeDict d;

			d.put("flip", EMObject::EMDATA);
			d.put("usedot", EMObject::INT);
			d.put("maxshift", EMObject::INT, "Maximum translation in pixels");
			d.put("rfp_mode", EMObject::INT,"Either 0,1 or 2. A temporary flag for testing the rotational foot print");
			d.put("useflcf", EMObject::INT,"Use Fast Local Correlation Function rather than CCF for translational alignment");
			return d;
		}
		
		static const string NAME;
	};

	/** rotational, translational and flip alignment using real-space methods. slow
	 * @param flip
	 * @param maxshift Maximum translation in pixels
	 * */
	class RTFExhaustiveAligner:public Aligner
	{
	  public:
		virtual EMData * align(EMData * this_img, EMData * to_img,
					   const string & cmp_name="dot", const Dict& cmp_params = Dict()) const;
		virtual EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "sqeuclidean", Dict());
		}

		virtual string get_name() const
		{
			return NAME;
		}

		virtual string get_desc() const
		{
			return "Experimental full 2D alignment with handedness check using semi-exhaustive search (not necessarily better than RTFBest)";
		}

		static Aligner *NEW()
		{
			return new RTFExhaustiveAligner();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;

			d.put("flip", EMObject::EMDATA);
			d.put("maxshift", EMObject::INT, "Maximum translation in pixels");
			return d;
		}
		
		static const string NAME;
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
		virtual EMData * align(EMData * this_img, EMData * to_img,
						const string & cmp_name, const Dict& cmp_params) const;
		virtual EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "sqeuclidean", Dict());
		}
		virtual string get_name() const
		{
			return NAME;
		}

		virtual string get_desc() const
		{
			return "Experimental full 2D alignment with handedness check using more exhaustive search (not necessarily better than RTFBest)";
		}

		static Aligner *NEW()
		{
			return new RTFSlowExhaustiveAligner();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("flip", EMObject::EMDATA,"Optional. This is the flipped version of the images that is being aligned. If specified it will be used for the handedness check, if not a flipped copy of the image will be made");
			d.put("maxshift", EMObject::INT,"The maximum length of the detectable translational shift");
			d.put("transtep", EMObject::FLOAT,"The translation step to take when honing the alignment, which occurs after coarse alignment");
			d.put("angstep", EMObject::FLOAT,"The angular step (in degrees) to take in the exhaustive search for the solution angle. Typically very small i.e. 3 or smaller.");
			return d;
		}
		
		static const string NAME;
	};

	/** refine alignment. Refines a preliminary 2D alignment using a simplex algorithm. Subpixel precision.
     */
	class RefineAligner:public Aligner
	{
	  public:
		virtual EMData * align(EMData * this_img, EMData * to_img,
					   const string & cmp_name="dot", const Dict& cmp_params = Dict()) const;

		virtual EMData * align(EMData * this_img, EMData * to_img) const
		{
			return align(this_img, to_img, "sqeuclidean", Dict());
		}

		virtual string get_name() const
		{
			return NAME;
		}

		virtual string get_desc() const
		{
			return "Refines a preliminary 2D alignment using a simplex algorithm. Subpixel precision.";
		}

		static Aligner *NEW()
		{
			return new RefineAligner();
		}

		virtual TypeDict get_param_types() const
		{
			TypeDict d;

			d.put("mode", EMObject::INT, "Currently unused");
			d.put("xform.align2d", EMObject::TRANSFORM, "The Transform storing the starting guess. If unspecified the identity matrix is used");
			d.put("stepx", EMObject::FLOAT, "The x increment used to create the starting simplex. Default is 1");
			d.put("stepy", EMObject::FLOAT, "The y increment used to create the starting simplex. Default is 1");
			d.put("stepaz", EMObject::FLOAT, "The rotational increment used to create the starting simplex. Default is 5");
			d.put("precision", EMObject::FLOAT, "The precision which, if achieved, can stop the iterative refinement before reaching the maximum iterations. Default is 0.04.");
			d.put("maxiter", EMObject::INT,"The maximum number of iterations that can be performed by the Simplex minimizer");
			d.put("maxshift", EMObject::INT,"Maximum translation in pixels in any direction. If the solution yields a shift beyond this value in any direction, then the refinement is judged a failure and the original alignment is used as the solution.");
			return d;
		}
		
		static const string NAME;
	};

	//No longer does this (doesn't work so well):
	/*
	   Uses quaternions extensively, separates the in plane (phi) rotation and the 2 rotations that define
	 * a point on the sphere (POTS), manipulating them independently. The POTS is"jiggled" in a local circular
	 * sub manifold of the sphere (of radius stepdelta). Phi is allowed to vary as a normal variable. The result
	 * is combined into a single transform and then a similarity is computed, etc.
	
	*/
	/** Refine alignment. Refines a preliminary 3D alignment using a simplex algorithm. Subpixel precision.
	 * Target function for the simplex  algorithm is a transformation of the 3D model by az, alt, phi, tx, ty, tz
	 * The simplex algorithm downs the function downhill in a ameboa like fasion, hence it may get stuck in a local 
	 * minima if the two 3D models are already roughly aligned.
	 * @author David Woolford and John Flanagan
	 * @date June 23 2009 and Oct 8th 2010
	 */
	class Refine3DAligner:public Aligner
	{
		public:
			virtual EMData * align(EMData * this_img, EMData * to_img,
						   const string & cmp_name="sqeuclidean", const Dict& cmp_params = Dict()) const;

			virtual EMData * align(EMData * this_img, EMData * to_img) const
			{
				return align(this_img, to_img, "sqeuclidean", Dict());
			}

			virtual string get_name() const
			{
				return NAME;
			}

			virtual string get_desc() const
			{
				return "Refines a preliminary 3D alignment using a simplex algorithm. Subpixel precision.";
			}

			static Aligner *NEW()
			{
				return new Refine3DAligner();
			}

			virtual TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("xform.align3d", EMObject::TRANSFORM,"The Transform storing the starting guess. If unspecified the identity matrix is used");
				d.put("stepx", EMObject::FLOAT, "The x increment used to create the starting simplex. Default is 1");
				d.put("stepy", EMObject::FLOAT, "The y increment used to create the starting simplex. Default is 1");
				d.put("stepz", EMObject::FLOAT, "The z increment used to create the starting simplex. Default is 1." );
				d.put("stepaz", EMObject::FLOAT, "The az increment used to create the starting simplex. Default is 5." );
				d.put("stepalt", EMObject::FLOAT, "The alt increment used to create the starting simplex. Default is 5." );
				d.put("stepphi", EMObject::FLOAT, "The phi incremenent used to create the starting simplex. Default is 5." );
				//d.put("stepdelta", EMObject::FLOAT,"The angular increment which represents a good initial step along the sphere, thinking in terms of quaternions. Default is 5.");
				d.put("precision", EMObject::FLOAT, "The precision which, if achieved, can stop the iterative refinement before reaching the maximum iterations. Default is 0.04." );
				d.put("maxiter", EMObject::INT, "The maximum number of iterations that can be performed by the Simplex minimizer. Default is 60.");
				d.put("maxshift", EMObject::INT,"Maximum translation in pixels in any direction. If the solution yields a shift beyond this value in any direction, then the refinement is judged a failure and the original alignment is used as the solution.");
				return d;
			}
			
			static const string NAME;
	};

	/** rotational and translational alignment using a square qrid of Altitude and Azimuth values (the phi range is specifiable)
	 * This aligner is ported from the original tomohunter.py - it is less efficient than searching on the sphere (RT3DSphereAligner),
	 * but very useful  if you want to search in a specific, small, local area.
	 * @author David Woolford (ported from Mike Schmid's e2tomohuntThis is the increment applied to the inplane rotationer code - Mike Schmid is the intellectual author)
	 * @date June 23 2009
	 */
	class RT3DGridAligner:public Aligner
	{
		public:
			/** See Aligner comments for more details
			 */
			virtual EMData * align(EMData * this_img, EMData * to_img,
						   const string & cmp_name, const Dict& cmp_params) const;
			/** See Aligner comments for more details
			 */
			virtual EMData * align(EMData * this_img, EMData * to_img) const
			{
				return align(this_img, to_img, "ccc.tomo", Dict());
			}


			/** See Aligner comments for more details
			 */
			virtual vector<Dict> xform_align_nbest(EMData * this_img, EMData * to_img, const unsigned int nsoln, const string & cmp_name, const Dict& cmp_params) const;

			virtual string get_name() const
			{
				return NAME;
			}

			virtual string get_desc() const
			{
				return "3D rotational and translational alignment using specified ranges and maximum shifts";
			}

			static Aligner *NEW()
			{
				return new RT3DGridAligner();
			}

			virtual TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("daz", EMObject::FLOAT,"The angle increment in the azimuth direction. Default is 10");
				d.put("raz", EMObject::FLOAT,"The range of angles to sample in the azimuth direction. Default is 360.");
				d.put("dphi", EMObject::FLOAT,"The angle increment in the phi direction. Default is 10.");
				d.put("rphi", EMObject::FLOAT,"The range of angles to sample in the phi direction. Default is 180.");
				d.put("dalt", EMObject::FLOAT,"The angle increment in the altitude direction. Default is 10.");
				d.put("ralt", EMObject::FLOAT,"The range of angles to sample in the altitude direction. Default is 180.");
				d.put("dotrans", EMObject::BOOL,"Do a translational search. Default is True(1)");
				d.put("search", EMObject::INT,"The maximum length of the detectable translational shift - if you supply this parameter you can not supply the maxshiftx, maxshifty or maxshiftz parameters. Each approach is mutually exclusive.");
				d.put("searchx", EMObject::INT,"The maximum length of the detectable translational shift in the x direction- if you supply this parameter you can not supply the maxshift parameters. Default is 3.");
				d.put("searchy", EMObject::INT,"The maximum length of the detectable translational shift in the y direction- if you supply this parameter you can not supply the maxshift parameters. Default is 3.");
				d.put("searchz", EMObject::INT,"The maximum length of the detectable translational shift in the z direction- if you supply this parameter you can not supply the maxshift parameters. Default is 3");
				d.put("verbose", EMObject::BOOL,"Turn this on to have useful information printed to standard out.");
				return d;
			}
			
			static const string NAME;
	};

	/** 3D rotational and translational alignment using spherical sampling, can reduce the search space based on symmetry.
	 * can also make use of different OrientationGenerators (random, for example)
	 * 160-180% more efficient than the RT3DGridAligner
	 * @author David Woolford
	 * @date June 23 2009
	 */
	class RT3DSphereAligner:public Aligner
	{
		public:
			/** See Aligner comments for more details
			 */
			virtual EMData * align(EMData * this_img, EMData * to_img,
								   const string & cmp_name, const Dict& cmp_params) const;
			/** See Aligner comments for more details
			 */
			virtual EMData * align(EMData * this_img, EMData * to_img) const
			{
				return align(this_img, to_img, "sqeuclidean", Dict());
			}


			/** See Aligner comments for more details
			 */
			virtual vector<Dict> xform_align_nbest(EMData * this_img, EMData * to_img, const unsigned int nsoln, const string & cmp_name, const Dict& cmp_params) const;

			virtual string get_name() const
			{
				return NAME;
			}

			virtual string get_desc() const
			{
				return "3D rotational and translational alignment using spherical sampling. Can reduce the search space if symmetry is supplied";
			}

			static Aligner *NEW()
			{
				return new RT3DSphereAligner();
			}

			virtual TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("sym", EMObject::STRING,"The symmtery to use as the basis of the spherical sampling. Default is c1 (asymmetry).");
				d.put("orientgen", EMObject::STRING,"Advanced. The orientation generation strategy. Default is eman");
				d.put("delta", EMObject::FLOAT,"Angle the separates points on the sphere. This is exclusive of the \'n\' paramater. Default is 10");
				d.put("n", EMObject::INT,"An alternative to the delta argument, this is the number of points you want generated on the sphere. Default is OFF");
				d.put("dphi", EMObject::FLOAT,"The angle increment in the phi direction. Default is 10.");
				d.put("rphi", EMObject::FLOAT,"The range of angles to sample in the phi direction. Default is 180.");
				d.put("dotrans", EMObject::BOOL,"Do a translational search. Default is True(1)");
				d.put("search", EMObject::INT,"The maximum length of the detectable translational shift - if you supply this parameter you can not supply the maxshiftx, maxshifty or maxshiftz parameters. Each approach is mutually exclusive.");
				d.put("searchx", EMObject::INT,"The maximum length of the detectable translational shift in the x direction- if you supply this parameter you can not supply the maxshift parameters. Default is 3.");
				d.put("searchy", EMObject::INT,"The maximum length of the detectable translational shift in the y direction- if you supply this parameter you can not supply the maxshift parameters. Default is 3.");
				d.put("searchz", EMObject::INT,"The maximum length of the detectable translational shift in the z direction- if you supply this parameter you can not supply the maxshift parameters. Default is 3");
				d.put("verbose", EMObject::BOOL,"Turn this on to have useful information printed to standard out.");
				return d;
			}
			
			static const string NAME;
	};
	

	class FRM2DAligner:public Aligner
			{
				public:
					virtual EMData * align(EMData * this_img, EMData * to_img,
							const string& cmp_name, const Dict& cmp_params=Dict()) const; //ming add ="frc"

					virtual EMData * align(EMData * this_img, EMData * to_img) const
					{
						return align(this_img, to_img, "frc", Dict());
					}

					string get_name() const
					{
						return NAME;
					}

					string get_desc() const
					{
						return "FRM2D uses two rotational parameters and one translational parameter";
					}

					static Aligner *NEW()
					{
						return new FRM2DAligner();
					}
					virtual	TypeDict get_param_types() const
					{
							TypeDict d;
							d.put("maxshift", EMObject::INT,"Maximum translation in pixels in any direction. If the solution yields a shift beyond this value in any direction, then the refinement is judged a failure and the original alignment is used as the solution.");

							//d.put("p_max", EMObject::FLOAT,"p_max is");
							return d;
					}

					static const string NAME;
		};

	
	class CUDA_Aligner
	{
	  public:
	  	CUDA_Aligner();

		void finish();

		void setup(int nima, int nx, int ny, int ring_length, int nring, int ou, float step, int kx, int ky, bool ctf);

		void insert_image(EMData *image, int num);

		void filter_stack(vector<float> ctf_params, int id);
		
		void sum_oe(vector<float> ctf_params, vector<float> ali_params, EMData* ave1, EMData *ave2, int id);

		vector<float> alignment_2d(EMData *ref_image, vector<float> sx, vector<float> sy, int id, int silent);

		vector<float> ali2d_single_iter(EMData *ref_image, vector<float> ali_params, float csx, float csy, int id, int silent, float delta);

	  private:
	        float *image_stack, *image_stack_filtered;
		float *ccf;
		int NIMA, NX, NY, RING_LENGTH, NRING, OU, KX, KY;
		bool CTF;
		float STEP;
	};


	template <> Factory < Aligner >::Factory();

	void dump_aligners();
	map<string, vector<string> > dump_aligners_list();
}

#endif
