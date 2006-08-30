/**
 * $Id$
 */

#ifndef eman_analyzer_sparx_h__
#define eman_analyzer_sparx_h__

#include "analyzer.h"

namespace EMAN
{
	/** 
	 */
	class PcaAnalyzer : public Analyzer
	{
	  public:
		PcaAnalyzer() : mask(0), nvec(0) {}
/*		virtual ~PcaAnalyzer() {
			if(mask) {
				delete mask;
				mask = 0;
			}
		}
*/		
		virtual vector<EMData*> analyze(const vector<EMData*> & images);
		
		string get_name() const
		{
			return "pca";
		}	  	
		
		string get_desc() const
		{
			return "";
		}
		
		static Analyzer * NEW()
		{
			return new PcaAnalyzer();
		}
		
		void set_params(const Dict & new_params)
		{
			params = new_params;
			mask = params["mask"];
			nvec = params["nvec"];
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("mask", EMObject::EMDATA);
			d.put("nvec", EMObject::INT);
			return d;
		}
		
	  protected:
		EMData * mask;
		int nvec;
	}; 
}

#endif	//eman_analyzer_sparx_h__ 
