#include "emdata.h"

using namespace EMAN;


int main()
{
	EMData * image = new EMData();
	const char* test_imagename = Util::get_debug_image("groel2d.mrc");

	int err = 0;
	try {
		image->read_image(test_imagename);
		image->set_size(2*image->get_xsize(), image->get_ysize(), image->get_zsize());
		image->write_image(test_imagename, 0, EMUtil::IMAGE_MRC);
	}
	catch(E2Exception & e) {
		err = 1;
		printf("%s\n", e.what());
	} 
		 
	return err;
}
	
