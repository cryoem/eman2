#include "emdata.h"
#include "testutil.h"
using namespace EMAN;

int main()
{
	EMData * image = new EMData();
	EMData * image2 = new EMData();
	
	string test_imagename = TestUtil::get_debug_image("search.dm3");
	int err = 0;
	
	try {
		image->read_image(test_imagename);
		float * data = image->get_data();
		int nx = image->get_xsize();
		int ny = image->get_ysize();
		int nz = image->get_zsize();

		image2->set_shared_data(nx,ny,nz, data);
		image2->write_image("search2.mrc");

	}
	catch(E2Exception & e) {
		err = 1;
		printf("%s\n", e.what());
	}

	if( image )
	{
		delete image;
		image = 0;
	}
	if( image2 )
	{
		delete image2;
		image2 = 0;
	}
	
	return err;
}
	

	
		
		
