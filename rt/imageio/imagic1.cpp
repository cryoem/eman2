#include "emdata.h"

using namespace EMAN;

void append_to_newfile()
{
	EMData * e = new EMData();
	const char * img1 = Util::get_debug_image("search.dm3");
	e->read_image(img1);
	e->append_image("appendnew.img", EMUtil::IMAGE_IMAGIC);
	delete e;
	e = 0;
}

void append_to_existing_file()
{
	EMData * e = new EMData();
	const char * img1 = Util::get_debug_image("x.hed");
	e->read_image(img1, 0, false, 0, true);
	e->write_image("appendexisting.hed", 0, EMUtil::IMAGE_IMAGIC);
	delete e;
	e = 0;

	EMData * e2 = new EMData();
	const char * img2 = Util::get_debug_image("square288.mrc");
	e2->read_image(img2, 0);
	e2->append_image("appendexisting.img", EMUtil::IMAGE_IMAGIC);
	delete e2;
	e2 = 0;

}

void insert_to_newfile()
{
	EMData * e = new EMData();
	const char * img1 = Util::get_debug_image("square288.mrc");
	e->read_image(img1);
	e->write_image("insertnew.img", 4, EMUtil::IMAGE_IMAGIC);
	delete e;
	e = 0;
}

void insert_beyond_existing_file()
{
	const char * outfile = "insert_beyond_existing.hed";
	EMData * e = new EMData();
	const char * img1 = Util::get_debug_image("x.hed");
	e->read_image(img1, 0, false, 0, true);
	e->write_image(outfile, 0, EMUtil::IMAGE_IMAGIC);
	delete e;
	e = 0;

	EMData * e2 = new EMData();
	const char * img2 = Util::get_debug_image("square288.mrc");
	e2->read_image(img2, 0);
	e2->write_image(outfile, 9, EMUtil::IMAGE_IMAGIC);

	e2->write_image(outfile, 14, EMUtil::IMAGE_IMAGIC);
	
	delete e2;
	e2 = 0;
}

void insert_inside_existing_file()
{
	const char * outfile = "insert_in_existing.img";
	EMData * e = new EMData();
	const char * img1 = Util::get_debug_image("x.hed");
	e->read_image(img1, 0, false, 0, true);
	e->write_image(outfile, 0, EMUtil::IMAGE_IMAGIC);
	delete e;
	e = 0;

	EMData * e2 = new EMData();
	const char * img2 = Util::get_debug_image("square288.mrc");
	e2->read_image(img2, 0);
	e2->write_image(outfile, 2, EMUtil::IMAGE_IMAGIC);
	delete e2;
	e2 = 0;
}




int main(int argc, char *argv[])
{
	try {
		append_to_newfile();
		append_to_existing_file();
		insert_to_newfile();
		insert_beyond_existing_file();
		insert_inside_existing_file();
	}
	catch(Exception &e) {
		printf("%s", e.what());
	}
	return 0;
}
