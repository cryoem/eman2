#include "EMData.h"
#include "Euler.h"

const char* get_test_image()
{
	static char filename[256];
	static bool done = false;
	if (!done) {
		sprintf(filename, "%s/images/3d99.hed", getenv("HOME"));
		done = true;
	}
	return filename;
}

void rotate(EMData * image, float phi)
{
	char outfile[128];
	sprintf(outfile, "test_%d.hed", (int)phi);
	float f = (float)M_PI / 180;
	const char* imagefile = get_test_image();
	image->readImage(imagefile, -1);
	image->setRAlign(0,0,phi*f);
	image->rotateAndTranslate();
	image->writeIMAGIC3D(outfile);
}

int main()
{
	int err = 0;
	
	EMData *image = new EMData();
#if 0
	for (int i = 0; i <= 180; i += 45) {
		for (int j = 0; j < 360; j += 60) {
			for (int z = 0; z < 360; z += 45) {
				rotate(image, i, j, z);
			}
		}
	}
#endif

	rotate(image, 45);

	if( image )
	{
		delete image;
		image = 0;
	}
	
	return  err;
}
