#include "sparx.h"

using namespace EMAN;

void em2flex(EMData& image, Dict & header_dict, scitbx::af::shared<float> & flex_array)
{
	header_dict = image.get_attr_dict();
	size_t size = image.get_xsize() * image.get_ysize() * image.get_zsize();
	
	flex_array.reserve(size);
	float *data = image.get_data();
	
	for(size_t i = 0; i < size; i++) {
		flex_array.push_back(data[i]);
	}
}

void flex2em(const scitbx::af::shared<float> & flex_array, Dict & header_dict, EMData & image)
{
	int nx = header_dict["nx"];
	int ny = header_dict["ny"];
	int nz = header_dict["nz"];
	image.set_size(nx, ny, nz);

	size_t size = nx * ny * nz;
	float * data = image.get_data();
	
	for(size_t i = 0; i < size; i++) {
		data[i] = flex_array[i];
	}

	image.set_attr_dict(header_dict);
}
