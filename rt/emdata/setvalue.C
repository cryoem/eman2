#include "emdata.h"
#include <stdio.h>
#include <stdlib.h>

using namespace EMAN;

int main(int argc, char* argv[])
{
    if (argc != 2) {
	printf("usage: %s n\n", argv[0]);
	printf("   n: 0 or 1\n");
	return 1;
    }

    int method = atoi(argv[1]);

    int err = 0;
    EMData e1;
    err = e1.read_image("/home/lpeng/images/jj0880f.mrc");

    
    int nx = e1.get_x();
    int ny = e1.get_y();
    int nz = e1.get_z();
    
    if (!err) {
	if (method == 0) {
	    float* data = e1.get_data();
	    int size = nx*ny*nz;
	    for (int i = 0; i < size; i++) {
		data[i] = 1;
	    }
	    
	    e1.done_data();
	}
	else {
	    for (int i = 0; i < nz; i++) {
		for (int j = 0; j < ny; j++) {
		    for (int k = 0; k < nx; k++) {
			e1.set_value_at(k, j, i, 1);
		    }
		}
	    }
	}
    }
    
    return 0;
}
	
