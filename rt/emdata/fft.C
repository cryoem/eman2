#include "emdata.h"
#include "util.h"
#include <stdio.h>

using namespace EMAN;

int main(int argc, char* argv[])
{
    int err = 0;
    Util::set_log_level(argc, argv);
    
    char filename[1024];
    sprintf(filename, "%s/images/square.mrc", getenv("HOME"));    
    EMData e1;
    
    err = e1.read_image(filename);

    if (!err) {
	e1.do_fft();
	err = e1.write_image("square_fft.mrc", 0, EMUtil::IMAGE_MRC);
    }

    return err;
}
