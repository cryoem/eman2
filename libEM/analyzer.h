/**
 * $Id$
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
