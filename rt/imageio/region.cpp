#include "emdata.h"
#include "testutil.h"
#include <assert.h>

using namespace EMAN;


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


void test_region(EMUtil::ImageType imgtype, int ndims)
{
	const char * testfile = "";
	bool is_3d = false;
	
	if (imgtype == EMUtil::IMAGE_MRC) {
		testfile = "groel3d.mrc";
	}
	else if (imgtype == EMUtil::IMAGE_IMAGIC) {
		testfile = "start.hed";
		is_3d = true;
	}
		
		
	const char * imgfile = TestUtil::get_debug_image(testfile);
	string imgbase = Util::remove_filename_ext(testfile);
	string ext = Util::get_filename_ext(testfile);
	string imgfile2 = imgbase + "_write_region." + ext;
	string imgfile3 = imgbase + "_read_region." + ext;
	
	EMData e;
	e.read_image(imgfile, 0, false, 0, is_3d);
	e.write_image(imgfile2, 0, imgtype);
	
	int nx = e.get_xsize();
	int ny = e.get_ysize();
	int nz = e.get_zsize();
	
	int x0 = nx/4;
	int y0 = ny/4;
	int z0 = nz/4;
	int xsize = nx/2;
	int ysize = ny/2;
	int zsize = nz/2;
	if (zsize == 0) {
		zsize = 1;
	}

	int image_index = 0;
	Region region;
	if (ndims == 2) {
		region = Region(x0, y0, xsize, ysize);
		image_index = nz / 2;
	}
	else if (ndims == 3) {
		region = Region(x0, y0, z0, xsize, ysize, zsize);
	}
	
	EMData e2;
	e2.read_image(imgfile, image_index, false, &region);
	e2.write_image(imgfile3, 0, imgtype);

	EMData e3;
	e3.set_size(xsize, ysize, zsize);
	e3.to_zero();
	
	e3.write_image(imgfile2, image_index, imgtype, false, &region);
}

	

int main(int argc, char *argv[])
{
	try {
		test_region(EMUtil::IMAGE_MRC, 2);
		test_region(EMUtil::IMAGE_IMAGIC, 2);
					  
	}
	catch(E2Exception &e) {
		e.what();
	}
	return 0;
}
