#include "emdata.h"
#include "ctf.h"

using namespace EMAN;

int main(int argc, char *argv[])
{
    if (argc != 2) {
	printf("usage: %s imagefile\n", argv[0]);
	exit(1);
    }
    
    EMData e;
    e.read_image(argv[1]);
    Ctf *ctf = e.get_ctf();
    if (ctf) {
	printf("%s CTF: %s\n", argv[1], ctf->to_string().c_str());
    }
    e.write_image("a.h5", 0, EMUtil::IMAGE_HDF);
#if 0
    SimpleCtf *ctf2 = (SimpleCtf*)ctf;
    ctf2->cs = 100;
    ctf2->apix = 10;
    e.set_ctf(ctf2);

    e.write_image("b.h5", 0, EMUtil::IMAGE_HDF);
#endif
    
    return 0;
}
