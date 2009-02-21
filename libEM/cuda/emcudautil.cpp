

#include "../emdata.h"

using namespace EMAN;

void set_emdata_cuda_array_handle(const int handle, void* emdata_pointer)
{
	EMData* e = static_cast<EMData*>(emdata_pointer);
	e->set_cuda_array_handle(handle);
}

