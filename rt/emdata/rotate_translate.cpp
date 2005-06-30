#include "emdata.h"
#include "testutil.h"
#include <iostream>

using namespace EMAN;
using namespace std;

string get_test_image()
{
	return TestUtil::get_debug_image("tg.mrc");
}

void r1()
{
	EMData *image = new EMData();

	image->read_image(get_test_image());
	image->rotate(1.0329837512591338,3.7260642381912579,5.7671541529246966);
	image->write_image("r1.mrc");
	if( image )
	{
		delete image;
		image = 0;
	}
}

void t1()
{
	EMData *image = new EMData();

	image->read_image(get_test_image());
	image->rotate_translate(0,0,0,16,16,16);
	image->write_image("t1.mrc");
	if( image )
	{
		delete image;
		image = 0;
	}
}

void t2()
{
	EMData *image = new EMData();

	image->read_image(get_test_image());
	image->translate(-0.5f,-0.5f, 0.0f);
	image->write_image("t2.mrc");
	if( image )
	{
		delete image;
		image = 0;
	}
}

void rt1()
{
	EMData *image = new EMData();

	image->read_image(get_test_image());
	image->rotate_translate(1.0329837512591338,3.7260642381912579,5.7671541529246966,16,16,16);
	image->write_image("rt1.mrc");
	if( image ) 
	{
		delete image;
		image = 0;
	}
}

void rt2()
{
	EMData *image = new EMData();

	image->read_image(get_test_image());
	Transform3D r = Transform3D(Vec3f(16,16,16), 1.0329837512591338,3.7260642381912579,
							5.7671541529246966);
	image->rotate_translate(r);
	image->write_image("rt2.mrc");
	if( image )
	{
		delete image;
		image = 0;
	}
}

void compare_image()
{
	EMData *image1 = new EMData();
	EMData *image2 = new EMData();
	
	image1->read_image(get_test_image());
	image2->read_image(get_test_image());
	
	image1->translate(Vec3f(100.5, 100.5, 0));
	image1->write_image("tran1.mrc");
	
	Transform3D r = Transform3D(Vec3f(100.5, 100.5, 0), 0,0,0);
	image2->rotate_translate(r);
	image2->write_image("tran2.mrc");
	
	float *data1 = image1->get_data();
	float *data2 = image2->get_data();
	int size = sizeof(float) * image1->get_xsize() * image1->get_ysize() * image1->get_zsize();
	int cmp = memcmp( (void*)data1, (void*)data2, size );
	cout << "Compare reault is: " << cmp << endl;
	
	if( image1 ) {
		delete image1;
		image1 = 0;
	}
	if( image2 ) {
		delete image2;
		image2 = 0;
	}
}

int main()
{
	cout << "Starting to test rotate_translate..." << endl;
	
	try {
		//r1();
		//t1();
		//rt1();
		//rt2();
		//t2();
		compare_image();
	}
	catch (E2Exception & e) {
		e.what();
	}
	
	
	return 0;
}
