#include "emdata.h"
#include <math.h>

using namespace EMAN;

void rotate(EMData * image, float alt, float phi, const char * imagename)
{
	char outfile[128];
	
	float f = (float)M_PI / 180;
	const char* imagefile = Util::get_debug_image(imagename);
	image->read_image(imagefile, 0, false, 0, true);
	image->rotate(alt*f, 0, phi*f);
	
	sprintf(outfile, "%s_%d_%d.mrc", imagename, (int)alt, (int)phi);
	image->write_image(outfile, 0, EMUtil::IMAGE_MRC);
}

void test_rotate(EMData * image, const char * imagename)
{
	rotate(image, 0, 45, imagename);
#if 1
	rotate(image, 0, 90, imagename);
	rotate(image, 0, 180, imagename);
	
	rotate(image, 45, 0, imagename);
	rotate(image, 90, 0, imagename);
	rotate(image, 180, 0, imagename);
#endif
}


int main()
{
	EMData *image = new EMData();

	//test_rotate(image, "lattice.mrc");
	//test_rotate(image, "3d.mrc");
	test_rotate(image, "3d99.hed");
	delete image;
	image = 0;

	return 0;
}
