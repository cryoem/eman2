#include "emdata.h"

using namespace EMAN;


void test_clip(const IntSize & size0, const Region & region, const string & id)
{
	string old_file = "old" + id + ".img";
	string new_file = "new" + id + ".img";
	
	EMData * image = new EMData();
	image->set_size(size0[0], size0[1], size0[2]);
	image->to_one();
	image->write_image(old_file, 0, EMUtil::IMAGE_IMAGIC);
	
	EMData * image2 = image->get_clip(region);
	image2->write_image(new_file, 0, EMUtil::IMAGE_IMAGIC);

	delete image;
	image = 0;
	delete image2;
	image2 = 0;
}



int main()
{
	test_clip(IntSize(64,64,1), Region(-32,-32,128,128), "2do");
	test_clip(IntSize(64,64,64), Region(-32,-32,-32,128,128,128), "3do");
	test_clip(IntSize(164,164,1), Region(32,32,68,68), "2di");
	test_clip(IntSize(164,164,164), Region(32,32,32,68,68,68), "3di");

	return 0;
}






