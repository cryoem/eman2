#include "averager_template.h"
#include "emdata.h"

using namespace EMAN;

XYZAverager::XYZAverager() : result(0)
{
}

void XYZAverager::add_image(EMData * image)
{
	if (!image) {
		return;
	}
	result = new EMData();
}


EMData *XYZAverager::finish()
{
	return result;
}
