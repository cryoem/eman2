#include "emdata.h"
#include "testutil.h"

using namespace EMAN;

int main()
{
	EMData *a = new EMData();
	a->read_image(TestUtil::get_debug_image("square.mrc"));

	EMData *fft = a->do_fft();
	fft->write_image("aatest2.mrc");

	delete fft;
	fft = 0;
	
	delete a;
	a = 0;

	
	return 0;
}
