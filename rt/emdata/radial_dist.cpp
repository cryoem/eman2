#include "emdata.h"

using namespace EMAN;

int main(int argc, char *argv[])
{
    if (argc != 2) {
	printf("usage: %s imagefile\n", argv[0]);
	exit(1);
    }
    
    EMData e;
    e.read_image(argv[1]);
    int ny = e.get_ysize();
    e.calc_radial_dist(ny, 0, 0.5);

    return 0;
}
