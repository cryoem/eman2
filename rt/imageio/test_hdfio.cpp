#include <stdio.h>
#include <stdlib.h>
#include "emdata.h"
#include "util.h"
#include "testutil.h"

using namespace EMAN;

void test_read()
{
	EMData * e = new EMData();
	char filename[64];
	for (int i = 0; i < 1000; i++) {
		sprintf(filename, "test_hdfio_%d.h5", i);
		e->read_image(TestUtil::get_debug_image("samesize1.mrc"));
		e->write_image(filename, 0, EMUtil::IMAGE_HDF);
	}
	delete e;
	e = 0;
}

void test_write()
{
	EMData * e = new EMData();
	e->read_image(TestUtil::get_debug_image("3d.mrc"));
	e->write_image("test_hdfio_write_2.h5", 0, EMUtil::IMAGE_HDF);
	delete e;
	e = 0;
}


int main()
{
	test_read();
	test_write();
	return 0;
}

		
	
