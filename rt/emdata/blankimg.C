#include "emdata.h"
#include <stdio.h>

using namespace EMAN;

int main()
{
    int err = 0;
    EMData e1;
    err = e1.read_image("/home/lpeng/images/search.dm3");

    if (!err) {    
	EMData e2;
	e2.set_size(e1.get_x(), e1.get_y(), e1.get_z());
	float* data = e2.get_data();
	int size = e1.get_x() * e1.get_y() * e1.get_z();
    
	for (int i = 0; i < size; i++) {
	    data[i] = 0;
	}

	e2.done_data();
	err = e2.write_image("search.mrc", 0, EMUtil::IMAGE_MRC);
	Dict d = e2.get_attr_dict();
	EMUtil::dump_dict(d);
    }
    
    return err;
}
