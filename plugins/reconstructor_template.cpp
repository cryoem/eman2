#include "reconstructor_template.h"

using namespace EMAN;

XYZReconstructor::XYZReconstructor()
{
}

XYZReconstructor::~XYZReconstructor()
{
}

int XYZReconstructor::setup()
{
	return 0;
}

int XYZReconstructor::insert_slice(EMData * slice, const Rotation & euler)
{
	return 0;
}

EMData *XYZReconstructor::finish()
{
	return image;
}
