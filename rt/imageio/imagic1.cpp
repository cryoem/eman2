#include "emdata.h"

using namespace EMAN;

void append_to_newfile()
{
	EMData * e = new EMData();
	const char * img1 = Util::get_debug_image("search.dm3");
	e->read_image(img1);
	e->write_image("hello.img", 4, EMUtil::IMAGE_IMAGIC);
	delete e;
	e = 0;
}

void append_to_existing_file()
{
	EMData * e = new EMData();
	const char * img1 = Util::get_debug_image("x.hed");
	e->read_image(img1, 0, false, 0, true);
	e->write_image("x2.hed", 0, EMUtil::IMAGE_IMAGIC);
	delete e;
	e = 0;
#if 1
	EMData * e2 = new EMData();
	e2->read_image(img1, 0);
	e2->write_image("x2.img", 9, EMUtil::IMAGE_IMAGIC);
	delete e2;
	e2 = 0;
#endif
}


int main(int argc, char *argv[])
{
	try {
		//append_to_newfile();
		append_to_existing_file();
	}
	catch(Exception &e) {
		printf("%s", e.what());
	}
	return 0;
}
