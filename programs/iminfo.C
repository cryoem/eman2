#include "emdata.h"
#include <stdio.h>
#include "log.h"

using namespace EMAN;

void usage(const char* progname)
{
    printf("\n%s [-H] [-vN] <image file>\n", progname);
    printf("    -H: to display all header information\n");
    printf("   -vN: set log level. N=0,1,2,3, from most detailed to least detailed\n");
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
    
    printf("\nNumber of Images:  %d\n", nimg);
    printf("Image Format: %s\n", imgtype);

    EMData* d = new EMData();
    d->read_image(imagefile, 0, true);
    map<string, EMObject> dict = d->get_attr_dict();
    
    int nx = dict["nx"].get_int();
    int ny = dict["ny"].get_int();
    int nz = dict["nz"].get_int();
    
    printf("Image Dimensions: %dx%dx%d\n\n", nx, ny, nz);

    if (show_all_header) {
	printf("\nDetailed Header Information:\n");
	EMUtil::dump_dict(dict);
    }
    
    delete d;
    d = 0;
}


    
