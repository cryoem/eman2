#include "emdata.h"
#include <cstdio>

using namespace EMAN;

int main()
{
    int err = 0;
    EMData e1;

    char filename[1024];
    sprintf(filename, "%s/images/search.dm3", getenv("HOME"));

	try {
		e1.read_image(filename); 
 
		EMData e2;
		e2.set_size(e1.get_xsize(), e1.get_ysize(), e1.get_zsize());
		float* data = e2.get_data();
		int size = e1.get_xsize() * e1.get_ysize() * e1.get_zsize();
    
		for (int i = 0; i < size; i++) {
			data[i] = 0;
		}

		e2.done_data();
		e2.write_image("search.mrc", 0, EMUtil::IMAGE_MRC);
		Dict d = e2.get_attr_dict();
		EMUtil::dump_dict(d);
    }
	catch(...) {
		err = 1;
		throw;
	}
    return err;
}
