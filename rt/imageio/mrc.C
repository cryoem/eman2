#include "mrcio.h"
#include <stdlib.h>
#include "transform.h"
#include <string>
#include "util.h"

using namespace EMAN;
using std::string;


void test_read_mrc(const char* filename, int image_index, Region* area, bool is_3d, const char* outfile)
{
    MrcIO mrc(filename, ImageIO::READ_ONLY);
    map<string, EMObject> dict;
    
    int err = mrc.read_header(dict, image_index, area, is_3d);

    if (err) {
	fprintf(stderr, "mrc read_header() failed.\n");
	return;
    }

    Util::dump_dict(dict);
    
    int size = dict["nx"].get_int() * dict["ny"].get_int() * dict["nz"].get_int();

    float* data = new float[size];
    err = mrc.read_data(data, image_index, area, is_3d);

    if (err) {
	fprintf(stderr, "mrc read_data() failed.\n");
	return;
    }

    
    FILE* out = fopen(outfile, "wb");

    int nsteps =  dict["ny"].get_int() * dict["nz"].get_int();
    int nx = dict["nx"].get_int();
    int stepsize = nx*sizeof(float);
    
    for (int i = 0; i < nsteps; i++) {
	fwrite(&data[i+nx], stepsize, 1, out);
    }
    
    
    delete [] data;
    fclose(out);
}


int main()
{
    char file1d[256];
    char file3d[256];
#if 1
    sprintf(file1d, "/home/lpeng/images/tablet.mrc");
    Region good_2d1(10, 11.4, 30, 40);
    Region bad_2d1(-2, -4, 10, 20);
    Region bad_2d2(1, 2, 3000, 400);
    
    test_read_mrc(file1d, 0, 0, false, "tablet.out");
    test_read_mrc(file1d, 0, &good_2d1, false, "tablet_good1.out");
    
    test_read_mrc(file1d, -3, 0, false, "tablet.out");
    test_read_mrc(file1d, 0, &bad_2d1, false, "tablet_bad1.out");
    test_read_mrc(file1d, 0, &bad_2d2, false, "tablet_bad2.out");
#endif
#if 1
    sprintf(file3d, "/home/lpeng/images/3d.mrc");
    
    Region good_3d1(0, 0, 0, 20, 20, 20);
    Region bad_3d1(0, 1, 2, 400, 230, 5);
    Region bad_3d2(0, -3, -5, 3, 5, 9);

    // positive tests
    test_read_mrc(file3d, 40, 0, false, "3d-1.out");
    test_read_mrc(file3d, 0, 0, true, "3d-all-1.out");
    
    test_read_mrc(file3d, 10, 0, true, "3d-all-2.out");
    test_read_mrc(file3d, 0, &good_3d1, true, "3d_good1.out");

    // negative tests
    test_read_mrc(file3d, 120, 0, false, "3d.out");
    test_read_mrc(file3d, 0, &bad_3d1, true, "3d_bad1.out");
    test_read_mrc(file3d, 0, &bad_3d2, true, "3d_bad2.out");
#endif
    
    return 0;
}
