#include "emdata.h"

using namespace EMAN;

const char* get_test_image()
{
	static char filename[256];
	static bool done = false;
	if (!done) {
		sprintf(filename, "%s/images/groel2d.mrc", getenv("HOME"));
		done = true;
	}
	return filename;
}

int main()
{
	EMData * image = new EMData();
	const char* test_imagename = get_test_image();

	int err = 0;
	try {
		image->read_image(test_imagename);
		image->set_size(2*image->get_xsize(), image->get_ysize(), image->get_zsize());
		image->write_image(test_imagename, 0, EMUtil::IMAGE_MRC);
	}
	catch(Exception *e) {
		err = 1;
		e->dump();
	} 
		 
	return err;
}
	
