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
			for (int k = 0; k < size/2; k++) {
				data[j*size+k]  = 1;
			}
		}
    }
    
    for (int l = 0; l < size; l++) {
		for (int i = 1; i < nslices-1; i++) {
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

	try {
		EMData e1;
		e1.read_image(filename);

		EMData* fft1 = e1.do_fft();
		fft1->write_image("lattice_fft.mrc", 0, EMUtil::IMAGE_MRC);

		EMData* fft1_copy = fft1->copy();
		
		EMData* ift1 = fft1->do_ift();	    
		ift1->write_image("lattice_fft_ift.mrc", 0, EMUtil::IMAGE_MRC);

		EMData* ift2 = fft1_copy->do_ift();	    
		ift2->write_image("lattice_fft_ift2.mrc", 0, EMUtil::IMAGE_MRC);

	}
	catch(...) {
		err = 1;
	}
    
    return err;
}
