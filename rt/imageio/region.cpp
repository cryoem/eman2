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


void test_region(EMUtil::ImageType imgtype, const char * testfile)
{	
	bool is_3d = false;
	if (imgtype == EMUtil::IMAGE_IMAGIC) {
		is_3d = true;
	}
		
	const char * imgfile = TestUtil::get_debug_image(testfile);
	string imgbase = Util::remove_filename_ext(testfile);
	string ext = Util::get_filename_ext(testfile);
	string writefile_2d = imgbase + "_write_region_2d." + ext;
	string writefile_3d = imgbase + "_write_region_3d." + ext;
	string readfile_2d = imgbase + "_read_region_2d." + ext;
	string readfile_3d = imgbase + "_read_region_3d." + ext;
	
	EMData e;
	e.read_image(imgfile, 0, false, 0, is_3d);
	e.write_image(writefile_2d, 0, imgtype);
	e.write_image(writefile_3d, 0, imgtype);
	
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
	int ndims = e.get_ndim();
	
	Region region_2d;
	Region region_3d;
	
	if (ndims == 2) {
		region_2d = Region(x0, y0, xsize, ysize);
	}
	else if (ndims == 3) {
		region_2d = Region(x0, y0, z0, xsize, ysize, 1);
		region_3d = Region(x0, y0, z0, xsize, ysize, zsize);
	}
	
	
	EMData e2;
	e2.read_image(imgfile, 0, false, &region_2d, is_3d);
	e2.write_image(readfile_2d, 0, imgtype);

	if (ndims == 3) {
		EMData e4;
		e4.read_image(imgfile, 0, false, &region_3d, is_3d);
		e4.write_image(readfile_3d, 0, imgtype);
	}

	EMData e3;
	e3.set_size(xsize, ysize, zsize);
	e3.to_zero();
	
	e3.write_image(writefile_2d, 0, imgtype, false, &region_2d);
	if (ndims == 3) {
		e3.write_image(writefile_3d, 0, imgtype, false, &region_3d);
	}
}

	

int main(int argc, char *argv[])
{
	try {
#if 0
		test_region(EMUtil::IMAGE_MRC, "groel3d.mrc");
		test_region(EMUtil::IMAGE_MRC, "samesize1.mrc");
		test_region(EMUtil::IMAGE_MRC, "tablet.mrc");
#endif
#if 1
		test_region(EMUtil::IMAGE_IMAGIC, "start.hed");
#endif
		
	}
	catch(E2Exception &e) {
		printf("%s\n", e.what());
	}
	return 0;
}
