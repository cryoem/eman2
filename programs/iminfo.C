#include "emdata.h"
#include <stdio.h>
#include "log.h"
#include "ctf.h"

using namespace EMAN;

void usage(const char* progname)
{
    printf("\n%s [-H] [-vN] <image file>\n", progname);
    printf("    -H: to display all header information\n");
    printf("   -vN: set verbosity level. N=0,1,2,3. large N means more verbose.\n");
    printf("\n");
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
	usage(argv[0]);
	exit(1);
    }

    Log::LogLevel log_level = Log::ERROR_LOG;
    bool show_all_header = false;
    
    for (int i = 1; i < argc - 1; i++) {
	if (strncmp(argv[i], "-v", 2) == 0) {
	    char str1[32];
	    strcpy(str1, argv[1]+2);
	    log_level = (Log::LogLevel) atoi(str1);
	}
	else if (strcmp(argv[i], "-H") == 0) {
	    show_all_header = true;
	}
    }
    
    Log::logger()->set_level(log_level);
    
    
    const char* imagefile = argv[argc-1];
    
    int nimg = EMUtil::get_image_count(imagefile);
    const char* imgtype = EMUtil::get_imagetype_name(EMUtil::get_image_type(imagefile));
    
    printf("\n%20s: %d\n", "Number of Images", nimg);
    printf("%20s: %s\n", "Image Format", imgtype);

    EMData* d = new EMData();
    d->read_image(imagefile, 0, EMData::HEADER_ONLY);
    Dict dict = d->get_attr_dict();
    
    int nx = dict.get("nx").get_int();
    int ny = dict.get("ny").get_int();
    int nz = dict.get("nz").get_int();
    
    printf("%20s: %dx%dx%d\n", "Image Dimensions", nx, ny, nz);

    Ctf* ctf = d->get_ctf();
    if (ctf) {
	printf("CTF: %s\n", ctf->to_string().c_str());
    }
    
    if (show_all_header) {
	printf("\nDetailed Header Information:\n");
	EMUtil::dump_dict(dict);
    }
    printf("\n");
    
    delete d;
    d = 0;
}


    
