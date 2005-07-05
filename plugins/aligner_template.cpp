#include "aligner_template.h"

using namespace EMAN;

EMData *XYZAligner::align(EMData * this_img, EMData * to_img, 
			const string & cmp_name, const Dict& cmp_params) const
{
	if (!this_img) {
		return 0;
	}
#if 0

	int param1 = params["param1"];
	float param2 = params["param2"];
#endif

	return 0;
}
