#include "emdata.h"
#include "testutil.h"

using namespace EMAN;

int main()
{
	EMData e;
	vector<EMData*> all_imgs = e.read_images(TestUtil::get_debug_image("start.hed"));

	for (size_t i = 0; i < all_imgs.size(); i++) {
		all_imgs[i]->append_image("start_append.hed");
		all_imgs[i]->append_image("start_append.spi", EMUtil::IMAGE_SPIDER);
	}

	return 0;
}
