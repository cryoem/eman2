#include "sparx.h"

using namespace EMAN;

scitbx::af::shared<float> em2flex(const EMData& image);
{
	size_t size = image.get_xsize() * image.get_ysize() * image.get_zsize();
	scitbx::af::shared<float> result;
	result.reserve(size);
	float *data = image.get_data();
	
	for(size_t i = 0; i < size; i++) {
		result.push_back(data[i]);
	}
	
	return result;
}

