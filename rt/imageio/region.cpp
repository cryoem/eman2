#include "emdata.h"

using namespace EMAN;

const char* get_test_image()
{
        static char filename[256];
        static bool done = false;
        if (!done) {
                sprintf(filename, "%s/images/tablet.mrc", getenv("HOME"));
                done = true;
        }
        return filename;
}

void test_em()
{
	EMData e;
	e.set_size(100,100,1);
	e.to_one();
	e.write_image("test.em", 0, EMUtil::IMAGE_EM);
	Region area(0, 0, 20, 40);
	e.to_zero();
	e.write_image("test.em", 0, EMUtil::IMAGE_EM, false, &area);
}


void test_mrc()
{
	EMData e;
	Region region(20, 50, 200, 400);
	const char * imgfile = get_test_image();
	e.read_image(imgfile, 0, false, &region);
	e.write_image("tablet_region.mrc");
}

int main(int argc, char *argv[])
{
	try {

		test_mrc();
		
					  
	}
	catch(Exception &e) {
		e.what();
	}
	return 0;
}
