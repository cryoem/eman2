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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#ifndef eman_analyzer_h__
#define eman_analyzer_h__ 

#include "emobject.h"

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
		virtual ~Analyzer()
		{
		}
		
		virtual vector<EMData*> analyze(const vector<EMData*> & images) = 0;
		
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
		
		/** Get Analyzer parameter information in a dictionary. Each
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
	
	
	
}

#endif	//eman_analyzer_h__
