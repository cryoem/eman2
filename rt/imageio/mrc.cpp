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
    e.set_attr_dict("MRC.label1", "Liwei Peng");
    int nx = e.get_xsize();
    int ny = e.get_ysize();
    EMData *c = e.get_clip(Region(0, 0, nx/2, ny/2));
    
    c->write_image("a.mrc", 0, EMUtil::IMAGE_MRC);

    delete c;
    c = 0;

    EMData *c2 = e.get_clip(Region(0, 0, nx*2, ny*2));
    c2->write_image("a2.mrc", 0, EMUtil::IMAGE_MRC);

    delete c2;
    c2 = 0;
    
    return 0;
}


