#include "emdata.h"
#include "testutil.h"

using namespace EMAN;

int main()
{
	EMData *a = new EMData();
	a->read_image(TestUtil::get_debug_image("3d86_1.mrc"));

	EMData *b = new EMData();
	b->read_image(TestUtil::get_debug_image("3d86_2.mrc"));
#if 1
	EMData *c = a->calc_ccf(b, 0);
	delete b;
	b = 0;
#endif
	
#if 0
	EMData *fft = a->do_fft();
	EMData *c = fft->copy(false);
#endif
	
	c->write_image("aatest2.mrc");

	delete a;
	a = 0;

	
	return 0;
}
