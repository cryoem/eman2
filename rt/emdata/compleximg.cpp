#include "emdata.h"

using namespace EMAN;

int main()
{
	const char* img1 = Util::get_debug_image("square.mrc");
	
	EMData *a = new EMData();
	a->read_image(img1);

	EMData *fft = a->do_fft();
	fft->write_image("aatest2.mrc");

	delete fft;
	fft = 0;
	
	delete a;
	a = 0;

	
	return 0;
}
