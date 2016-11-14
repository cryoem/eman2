/**
 * $Id$
 */

/*
 * Author: Chao Yang
 * Copyright (c) 2000-2006
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
 */

#ifndef eman_analyzer_sparx_h__
#define eman_analyzer_sparx_h__

#include "analyzer.h"

namespace EMAN
{
	/** Principal component analysis
	 */
	class PCAsmall : public Analyzer
	{
	  public:
		PCAsmall() : mask(0), nvec(0) {}
		
//		virtual int insert_image(EMData * image) {
//		images.push_back(image);
//		return 0;
//
//}
		virtual int insert_image(EMData * image);
		
		virtual vector<EMData*> analyze();
		
		string get_name() const
		{
			return NAME;
		}	  	
		
		string get_desc() const
		{
			return "Principal component analysis";
		}
		
		static Analyzer * NEW()
		{
			return new PCAsmall();
		}
		
		void set_params(const Dict & new_params);

//		void set_params(const Dict & new_params)
//		{
//			params = new_params;
//			mask = params["mask"];
//			nvec = params["nvec"];
//		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("mask", EMObject::EMDATA, "mask image");
			d.put("nvec", EMObject::INT, "number of desired principal components");
			return d;
		}
		
		static const string NAME;
		
	  protected:
		EMData * mask;
		int nvec;	//number of desired principal components

          private:
                float *covmat; // covariance matrix 
                int   ncov;    // dimension of the covariance matrix
                int   nimages; // number of images
                float *eigval; // array for storing computed eigvalues
	}; 

        //-------------------------------------------------------------

	class PCAlarge : public Analyzer
	{
	  public:
		PCAlarge() : mask(0), nvec(0) {}
		
		virtual int insert_image(EMData * image);
		
		virtual vector<EMData*> analyze();
		
		string get_name() const
		{
			return NAME;
		}	  	
		
		string get_desc() const
		{
			return "Principal component analysis";
		}
		
		static Analyzer * NEW()
		{
			return new PCAlarge();
		}
		
		void set_params(const Dict & new_params);

                int Lanczos(const string &maskedimages, int *kstep, 
                            float  *diag, float *subdiag, float *V, 
                            float  *beta);


		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("mask", EMObject::EMDATA, "mask image");
			d.put("nvec", EMObject::INT, "number of desired principal components");
			d.put("tmpfile", EMObject::STRING, "Name of temporary file during processing");
			return d;
		}
		
		static const string NAME;
		
	  protected:
		EMData * mask;
		int nvec;	//number of desired principal components

          private:
                int   ncov;    // dimension of the covariance matrix
                int   nimages; // number of images used in the analysis
                float *eigval; // array for storing computed eigvalues
	}; 

	class varimax : public Analyzer
	{
	  public:
		varimax() : m_mask(NULL) {}
		
		virtual int insert_image(EMData * image);
		
		virtual vector<EMData*> analyze();
		
		string get_name() const
		{
			return NAME;
		}	  	
		
		string get_desc() const
		{
			return "varimax rotation of PCA results";
		}
		
		static Analyzer * NEW()
		{
			return new varimax();
		}
		
		virtual void set_params(const Dict & new_params);

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("mask", EMObject::EMDATA, "mask image");
			return d;
		}
		
		static const string NAME;
		
	  private:
		int m_nlen;
		int m_nfac;
		EMData *m_mask;
		vector<float> m_data;
	}; 
}

#endif	//eman_analyzer_sparx_h__ 
