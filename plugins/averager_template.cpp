#include "averager_template.h"

using namespace EMAN;

EMData *XYZAverager::average(const vector<EMData*>& image_list) const
{
    if (image_list.size() == 0) {
	return 0;
    }
    
    return 0;
}
