#include "emdata.h"
#include <stdio.h>
#include "log.h"
#include "ctf.h"

using namespace EMAN;

void usage(const char *progname)
{
    printf("\n%s [-H] [-vN] [-stat] <image file>\n", progname);
    printf("    -H: show all header information.\n");
    printf("   -vN: set verbosity level. N=0,1,2,3. large N means more verbose.\n");
    printf(" -stat: show statistical information about the image(s).\n");
    printf("\n");
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
	usage(argv[0]);
	exit(1);
    }

    bool show_all_header = false;
    bool stat = false;

    Util::set_log_level(argc, argv);
    for (int i = 1; i < argc - 1; i++) {
	if (strcmp(argv[i], "-H") == 0) {
	    show_all_header = true;
	}
	else if (strcmp(argv[i], "-stat") == 0) {
	    stat = true;
	}
    }
    
    const char* imagefile = argv[argc-1];
    
    int nimg = EMUtil::get_image_count(imagefile);
    const char* imgtype = EMUtil::get_imagetype_name(EMUtil::get_image_type(imagefile));
    
    printf("\n%20s: %d\n", "Number of Images", nimg);
    printf("%20s: %s\n", "Image Format", imgtype);

    EMData* d = new EMData();
    if (!stat) {
	d->read_image(imagefile, 0, EMData::HEADER_ONLY);
    }
    else {
	d->read_image(imagefile, 0, EMData::HEADER_AND_DATA);
    }
    
    printf("%20s: %d x %d x %d\n", "Image Dimensions",
	   d->get_xsize(), d->get_ysize(), d->get_zsize());

    if (stat) {
	printf("mean=%1.3g sigma=%1.3g skewness=%1.3g kurtosis=%1.3g\n",
	       d->get_mean(), d->get_sigma(), d->get_skewness(), d->get_kurtosis());
    }
    
    Ctf* ctf = d->get_ctf();
    if (ctf) {
	printf("CTF: %s\n", ctf->to_string().c_str());
    }
    
    if (show_all_header) {
	Dict dict = d->get_attr_dict();
	printf("\nDetailed Header Information:\n");
	EMUtil::dump_dict(dict);
    }
    printf("\n");
    
    delete d;
    d = 0;
}


    
