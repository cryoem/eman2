#include "emdata.h"

using namespace EMAN;

void test_2d()
{
	EMData * image = new EMData();
	image->set_size(64,64,1);
	image->to_one();
	
	image->write_image("test2d1.img", 0, EMUtil::IMAGE_IMAGIC);
	EMData * image2 = image->get_clip(Region(-32,-32,128,128));
	image2->write_image("test2d2.img", 0, EMUtil::IMAGE_IMAGIC);

	delete image;
	image = 0;
	delete image2;
	image2 = 0;

}

void test_3d()
{
	EMData * image = new EMData();
	image->set_size(64,64,64);
	image->to_one();
	
	image->write_image("test3d1.img", 0, EMUtil::IMAGE_IMAGIC);
	EMData * image2 = image->get_clip(Region(-32,-32,-32,128,128,128));
	image2->write_image("test3d2.img", 0, EMUtil::IMAGE_IMAGIC);

	delete image;
	image = 0;
	delete image2;
	image2 = 0;
}



int main()
{
	test_2d();
	test_3d();
	return 0;
}
