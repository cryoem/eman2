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
    e.write_image("a.mrc", 0, EMUtil::IMAGE_MRC);
    return 0;
}
