#include "emdata.h"
#include "util.h"
#include <stdio.h>

using namespace EMAN;

void create_lattice_image(const char* filename)
{
    int size = 256;
    int nslices = 4;
    int width = 5;
    
    EMData e;
    e.set_size(size, size, 1);

    float * data = e.get_data();
	
    for (int i = 2; i < nslices; i++) {
	int start = size/nslices*i - width;
	int end = size/nslices*i + width;
	
	for (int j = start; j < end; j++) {
	    for (int k = 0; k < size; k++) {
		data[j*size+k]  = 1;
	    }
	}
    }
    
    for (int l = 0; l < size; l++) {
	for (int i = 1; i < nslices; i++) {
	    int start = size/nslices*i - width;
	    int end = size/nslices*i + width;
	    
	    for (int j = start; j < end; j++) {
		data[l*size+j] = 1;
	    }
	}
    }
    
    e.done_data();
    e.write_image(filename, 0, EMUtil::IMAGE_MRC);
    
    return;
}
    


int main(int argc, char* argv[])
{
    int err = 0;
    Util::set_log_level(argc, argv);

    const char* filename = "lattice.mrc";
    create_lattice_image(filename);
    
    EMData e1;
    err = e1.read_image(filename);

    if (!err) {
	e1.do_fft();
	err = e1.write_image("lattice_fft.mrc", 0, EMUtil::IMAGE_MRC);
    }
        
    return err;
}
