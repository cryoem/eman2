#include "emdata.h"
#include "testutil.h"

#include <assert.h>

using namespace EMAN;

void test_get_view_raw(EMData * image)
{
	const int ndims = 3;
	MArray3D marray = image->get_3dview();

	float * data = image->get_data();

	assert((int)marray.num_dimensions() == ndims);
	
	const size_t* shape = marray.shape();

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	assert((int)shape[0] == nx);
	assert((int)shape[1] == ny);
	assert((int)shape[2] == nz);
	
	int l = 0;
	for (int i = 0; i < nz; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nx; k++) {
				assert(marray[i][j][k] == data[l]);
				l++;
			}
		}
	}
}

void test_get_view_trans_2d(EMData * image)
{
	const int ndims = 2;
	const int x0 = 50;
	const int y0 = 150;
	
	MArray2D marray = image->get_2dview(x0,y0);

	float * data = image->get_data();

	assert((int)marray.num_dimensions() == ndims);
	
	const size_t* shape = marray.shape();

	int nx = image->get_xsize();
	int ny = image->get_ysize();

	assert((int)shape[0] == nx);
	assert((int)shape[1] == ny);
	
	int l = 0;
	for (int j = 0; j < ny; j++) {
		for (int k = 0; k < nx; k++) {
			assert(marray[x0+j][y0+k] == data[l]);
			l++;
		}
	}
}




int main()
{
	EMData *image3d = new EMData();
	image3d->read_image(TestUtil::get_debug_image("monomer.mrc"));
	test_get_view_raw(image3d);

	
	EMData *image2d = new EMData();
	image2d->read_image(TestUtil::get_debug_image("samesize1.mrc"));
	test_get_view_trans_2d(image2d);
	
	return 0;
}
	
