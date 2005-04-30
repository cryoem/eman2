#include "emdata.h"
#include "testutil.h"

using namespace EMAN;

string get_test_image()
{
	return TestUtil::get_debug_image("monomer.mrc");
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
	Transform3D r = Transform3D(Vec3f(16,16,16), 1.0329837512591338,3.7260642381912579,
							5.7671541529246966);
	image->rotate_translate(r);
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
	catch (E2Exception & e) {
		e.what();
	}
	
	
	return 0;
}
