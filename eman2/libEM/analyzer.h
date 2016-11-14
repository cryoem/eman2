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

#ifndef eman_analyzer_h__
#define eman_analyzer_h__

#include "emobject.h"
#include <gsl/gsl_linalg.h>
using std::vector;

namespace EMAN
{
	class EMData;

	/** Analyzer class defines a way to take a List of images as input,
	 * and returns a new List of images.
	 *
	 * Analyzer class is the base class for all analyzer classes. Each
     * specific analyzer has a unique ID name. This name is used to
     * call a analyzer.
     *
     * All Analyzer classes in EMAN are managed by a Factory
	 * pattern. So each Analyzer class must define:
	 *   - a unique name to idenfity itself in the factory.
	 *   - a static method to register itself in the factory.
	 *
	 */
	class Analyzer
	{
	  public:
	  	Analyzer() {}

		virtual ~Analyzer()
		{}

		/** insert a image to the list of input images
		 * @param image
		 * @return int 0 for success, <0 for fail
		 * */
		virtual int insert_image(EMData * image) = 0;

		/** insert a list of images to the list of input images
		 * @param image_list
		 * @return int 0 for success, <0 for fail
		 * */
		virtual int insert_images_list(vector<EMData *> image_list);

		/** main function for Analyzer, analyze input images and create output images
		 * @return vector<EMData *> result os images analysis
		 * */
		virtual vector<EMData*> analyze() = 0;

		/** Get the Analyzer's name. Each Analyzer is identified by a unique name.
		 * @return The Analyzer's name.
		 */
		virtual string get_name() const = 0;

		/** Get the Analyzer's description.
		 * @return The Analyzer's description.
		 */
		virtual string get_desc() const = 0;

		/** Set the Analyzer parameters using a key/value dictionary.
		 * @param new_params A dictionary containing the new parameters.
		 */
		virtual void set_params(const Dict & new_params)
		{
			params = new_params;
		}

		/** Get the Reconstructor's parameters in a key/value dictionary.
		 * @return A key/value pair dictionary containing the parameters.
		 */
		virtual Dict get_params() const
		{
			return params;
		}

		/** Get Analyzer parameter information in a dictionary. Each
		 * parameter has one record in the dictionary. Each record
		 * contains its name, data-type, and description.
		 *
		 * @return A dictionary containing the parameter info.
		 */
		virtual TypeDict get_param_types() const = 0;

	  protected:
		mutable Dict params;
		vector<EMData *> images;
	};

	/** Inertia Matrix computer
	 * Computes the Inertia Matrix for a 3-D volume
	 * @author Steve Ludtke
	 * @date 12/29/2013
	 * @param verbose Display progress if set, more detail with larger numbers
	 *
	 */
	class InertiaMatrixAnalyzer:public Analyzer
	{
	  public:
		InertiaMatrixAnalyzer() : verbose(0) {}

		virtual int insert_image(EMData *image) {
			images.push_back(image);
			if (images.size()>1) { printf("InertiaMatrixAnalyzer takes only a single image\n"); return 1; }
			return 0;
		}

		virtual vector<EMData*> analyze();

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Compute Inertia matrix for a volume";
		}

		static Analyzer * NEW()
		{
			return new InertiaMatrixAnalyzer();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("verbose", EMObject::INT, "Display progress if set, more detail with larger numbers");
			return d;
		}

		static const string NAME;

	  protected:
		int verbose;
		vector<EMData *> ret;		// This will contain only a single image
	};

	/** Shape characterization
	 * Computes a set of values characteristic of the shape of a volume.
	 * The first 3 values are distributions along X, Y and Z axes respectively, and are not orientation independent.
	 * @author Steve Ludtke
	 * @date 12/29/2013
	 * @param verbose Display progress if set, more detail with larger numbers
	 *
	 */
	class ShapeAnalyzer:public Analyzer
	{
	  public:
		ShapeAnalyzer() : verbose(0) {}

		virtual int insert_image(EMData *image) {
			images.push_back(image);
			if (images.size()>1) { printf("ShapeAnalyzer takes only a single image\n"); return 1; }
			return 0;
		}

		virtual vector<EMData*> analyze();

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Experimental. Computes a set of values characterizing a 3-D volume. Returns a 3x2x1 image containing X, Y and Z axial distributions using axis squared and axis linear weighting.";
		}

		static Analyzer * NEW()
		{
			return new ShapeAnalyzer();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("verbose", EMObject::INT, "Display progress if set, more detail with larger numbers");
			return d;
		}

		static const string NAME;

	  protected:
		int verbose;
		vector<EMData *> ret;		// This will contain only a single image
	};


	
	/** KMeansAnalyzer
	 * Performs k-means classification on a set of input images (shape/size arbitrary)
	 * returned result is a set of classification vectors
	 * @author Steve Ludtke
	 * @date 03/02/2008
	 * @param verbose Display progress if set, more detail with larger numbers (9 max)
	 * @param ncls number of desired classes
	 * @param maxiter maximum number of iterations
	 * @param minchange Terminate if fewer than minchange members move in an iteration
	 * @param mininclass Minumum number of particles to keep a class as good (not enforced at termination
	 * @param slowseed Instead of seeding all classes at once, it will gradually increase the number of classes by adding new seeds in groups with large standard deviations
	 * @param calcsigmamean Computes standard deviation of the mean image for each class-average (center), and returns them at the end of the list of centers
	 *
	 */
	class KMeansAnalyzer:public Analyzer
	{
	  public:
		KMeansAnalyzer() : ncls(0),verbose(0),minchange(0),maxiter(100),mininclass(2),slowseed(0) {}

		virtual int insert_image(EMData *image) {
			images.push_back(image);
			return 0;
		}

		virtual vector<EMData*> analyze();

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "k-means classification";
		}

		static Analyzer * NEW()
		{
			return new KMeansAnalyzer();
		}

		void set_params(const Dict & new_params);

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("verbose", EMObject::INT, "Display progress if set, more detail with larger numbers (9 max)");
			d.put("ncls", EMObject::INT, "number of desired classes");
			d.put("maxiter", EMObject::INT, "maximum number of iterations");
			d.put("minchange", EMObject::INT, "Terminate if fewer than minchange members move in an iteration");
			d.put("mininclass", EMObject::INT, "Minumum number of particles to keep a class as good (not enforced at termination");
			d.put("slowseed",EMObject::INT, "Instead of seeding all classes at once, it will gradually increase the number of classes by adding new seeds in groups with large standard deviations");
			d.put("calcsigmamean",EMObject::INT, "Computes standard deviation of the mean image for each class-average (center), and returns them at the end of the list of centers");
			return d;
		}

		static const string NAME;

	  protected:
		void update_centers(int sigmas=0);
		void reclassify();
		void reseed();

		vector<EMData *> centers;
		int ncls;	//number of desired classes
		int verbose;
		int minchange;
		int maxiter;
		int mininclass;
		int nchanged;
		int slowseed;
		int calcsigmamean;

	};

	/**Singular Value Decomposition from GSL. Comparable to pca
	 *@param mask mask image
	 *@param nvec number of desired basis vectors
	 *@param nimg total number of input images, required even with insert_image()
	 */
	class SVDAnalyzer : public Analyzer
	{
	  public:
		SVDAnalyzer() : mask(0), nvec(0), nimg(0), A(NULL) {}

		virtual int insert_image(EMData * image);

		virtual int insert_images_list(vector<EMData *> image_list) {
			vector<EMData*>::const_iterator iter;
			for(iter=image_list.begin(); iter!=image_list.end(); ++iter) {
				images.push_back(*iter);
			}
			return 0;
		}

		virtual vector<EMData*> analyze();

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Singular Value Decomposition from GSL. Comparable to pca";
		}

		static Analyzer * NEW()
		{
			return new SVDAnalyzer();
		}

		void set_params(const Dict & new_params);

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("mask", EMObject::EMDATA, "mask image");
			d.put("nvec", EMObject::INT, "number of desired basis vectors");
			d.put("nimg", EMObject::INT, "total number of input images, required even with insert_image()");
			return d;
		}

		static const string NAME;

	  protected:
		EMData * mask;
		int nvec;	//number of desired principal components
		int pixels;	// pixels under the mask
		int nimg; // number of input images

		private:
		int nsofar;
		gsl_matrix *A;
	};

		
	/**  Calculate the circular average around the center in real space.
	 *   @author: Muyuan Chen
	 *   @date: 04/2015
	 */
	class CircularAverageAnalyzer:public Analyzer
	{
	  public:
		CircularAverageAnalyzer() : verbose(0) {}

		virtual int insert_image(EMData *image) {
			images.push_back(image);
			return 0;
		}

		virtual vector<EMData*> analyze();

		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Calculate the circular average around the center in real space";
		}

		static Analyzer * NEW()
		{
			return new CircularAverageAnalyzer();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("verbose", EMObject::INT, "Display progress if set, more detail with larger numbers");
			d.put("maxr", EMObject::INT, "Maximum radius.");
			d.put("step", EMObject::INT, "Thickness of the ring.");
			return d;
		}

		static const string NAME;

	  protected:
		int verbose;
		vector<EMData *> ret;		// This will contain only a single image
	};
	
	template <> Factory < Analyzer >::Factory();

	void dump_analyzers();
	map<string, vector<string> > dump_analyzers_list();
}

#endif	//eman_analyzer_h__
