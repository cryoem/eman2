#include "emdata.h"

using namespace EMAN;

int main(int argc, char *argv[])
{
	try {
		EMData e;
		e.set_size(100,100,1);
		e.write_image("test.mrc");
		e.write_image("test.mrc");
		e.write_image("test.mrc", 0, EMUtil::IMAGE_MRC, true);
	}
	catch(Exception &e) {
		e.what();
	}
	return 0;
}
