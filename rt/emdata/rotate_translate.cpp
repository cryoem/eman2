#include "emdata.h"

using namespace EMAN;

const char* get_test_image()
{
	static char filename[256];
	static bool done = false;
	if (!done) {
		sprintf(filename, "%s/images/monomer.mrc", getenv("HOME"));
		done = true;
	}
	return filename;
}

void r1()
{
	EMData *image = new EMData();

	image->read_image(get_test_image());
	image->rotate(1.0329837512591338,3.7260642381912579,5.7671541529246966);
	image->write_image("r1.mrc");
	delete image;
	image = 0;
}

void t1()
{
	EMData *image = new EMData();

	image->read_image(get_test_image());
	image->rotate_translate(0,0,0,16,16,16);
	image->write_image("t1.mrc");
	delete image;
	image = 0;
}

void rt1()
{
	EMData *image = new EMData();

	image->read_image(get_test_image());
	image->rotate_translate(1.0329837512591338,3.7260642381912579,5.7671541529246966,16,16,16);
	image->write_image("rt1.mrc");
	delete image;
	image = 0;
}

void rt2()
{
	EMData *image = new EMData();

	image->read_image(get_test_image());
	Rotation r = Rotation(1.0329837512591338,3.7260642381912579,5.7671541529246966, Rotation::EMAN);
	image->rotate_translate(r,Vec3<float>(16,16,16));
	image->write_image("rt2.mrc");
	delete image;
	image = 0;
}

int main()
{
	try {
		r1();
		t1();
		rt1();
		rt2();
	}
	catch (Exception & e) {
		e.what();
	}
	
	
	return 0;
}
