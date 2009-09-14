/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
 * Copyright (c) 2000-2006 Baylor College of Medicine
 *
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * */

#include <cstring>
#include "emdata.h"
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
    EMData* d = 0;
	try {
		int nimg = EMUtil::get_image_count(imagefile);
        EMUtil::ImageType imgtype = EMUtil::get_image_type(imagefile);
		const char* imgtypename = EMUtil::get_imagetype_name(imgtype);
		int image_index = 0;
        if (imgtype == EMUtil::IMAGE_SPIDER && !stat) {
            image_index = -1;
        }
		printf("\n%20s: %d\n", "Number of Images", nimg);
		printf("%20s: %s\n", "Image Format", imgtypename);

		d = new EMData();

		if (!stat) {
			d->read_image(imagefile, image_index, true);
		}
		else {
			d->read_image(imagefile, image_index, false);
		}

		printf("%20s: %d x %d x %d\n", "Image Dimensions",
			   d->get_xsize(), d->get_ysize(), d->get_zsize());

		printf("%20s: %s\n", "Image Data Type",
				EMUtil::get_datatype_string((EMUtil::EMDataType)((int)(d->get_attr("datatype")))));

		if (stat) {
			printf("mean=%1.3g sigma=%1.3g skewness=%1.3g kurtosis=%1.3g\n",
				   (float) d->get_attr("mean"),
				   (float) d->get_attr("sigma"),
				   (float) d->get_attr("skewness"),
				   (float) d->get_attr("kurtosis"));
		}

		Ctf* ctf = d->get_ctf();
		if (ctf) {
			printf("CTF: %s\n", ctf->to_string().c_str());
			delete ctf;
			ctf = 0;
		}

		if (show_all_header) {
			Dict dict = d->get_attr_dict();
			printf("\nDetailed Header Information:\n");
			EMUtil::dump_dict(dict);
		}
		printf("\n");



    	if( d )
    	{
			delete d;
			d = 0;
    	}
	}
	catch(E2Exception  & e) {
		if (d) {
			delete d;
			d = 0;
		}
		printf("%s\n", e.what());
	}


	return 0;
}



