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
		image2->read_image(test_imagename);

		image->align("rotate_translate_flip", image2, Dict());

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
	

	
		
		
