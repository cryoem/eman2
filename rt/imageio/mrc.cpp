#include "emdata.h"

using namespace EMAN;

int main(int argc, char *argv[])
{
    EMData e;
    e.read_image(argv[1]);
    e.set_attr_dict("MRC.label1", "Liwei Peng");
    e.write_image("a.mrc", 0, EMUtil::IMAGE_MRC);
    return 0;
}
