#include "EMData.h"

const char* get_test_image(const char* img)
{
	static char filename[256];
	static bool done = false;
	if (!done) {
		sprintf(filename, "%s/images/%s", getenv("HOME"), img);
		done = true;
	}
	return filename;
}

int main()
{
	EMData *a = new EMData();
	const char * img1 = get_test_image("3d86_1.mrc");
	const char * img2 = get_test_image("3d86_2.mrc");
	
	a->readImage(img1, -1);

	EMData *b = new EMData();
	b->readImage(img2, -1);

	EMData *c = a->calcCCF(b, 0);
	
	printf("nx = %d, ny = %d, nz = %d\n", a->xSize(), a->ySize(), a->zSize());

	printf("nx = %d, ny = %d, nz = %d\n", c->xSize(), c->ySize(), c->zSize());
	
	c->writeImage("aaatest1.mrc", -1);

	if( a )
	{
		delete a;
		a = 0;
	}

	if( b )
	{
		delete b;
		b = 0;
	}
	
	return 0;
}
