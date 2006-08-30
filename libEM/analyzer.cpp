/**
 * $Id$
 */

#include "analyzer.h"
#include "sparx/analyzer_sparx.h"

using namespace EMAN;

namespace EMAN {

	template <> Factory < Analyzer >::Factory()
	{
		force_add(&PcaAnalyzer::NEW);
	}

}

vector<EMData*> PcaAnalyzer::analyze(const vector<EMData*> & images)
{
	
}

