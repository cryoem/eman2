#include "emdata.h"
#include "log.h"
#include "testutil.h"
using namespace EMAN;

void test1()
{
	try {
		EMData * e1 = new EMData();
		e1->read_image(TestUtil::get_debug_image("square288.mrc"));
		e1->calc_fourier_shell_correlation(0);
	}
	catch (_NullPointerException& e) {
		LOGERR("%s", e.what());
	}
}

void test2()
{
	try {
		EMData e1;
		e1.set_size(10,10,10);
		e1.dot(0, false);
	}
	catch (E2Exception & e) {
		LOGERR("%s", e.what());
	}
}


int main()
{
	test1();
	test2();
	return 0;
}

