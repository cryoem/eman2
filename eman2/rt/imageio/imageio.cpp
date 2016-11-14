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

#include <cstdio>
#include <cstdlib>
#include "emdata.h"
#include "util.h"
#include "transform.h"
#include "emutil.h"
#include "log.h"

using namespace EMAN;

static int total_ntests = 0;
static int fail_ntests = 0;
static int err_code = 0;

int write_image(EMData* em, const char* infile, char* outfile,
				int r_image_index, EMUtil::ImageType image_type,
				int w_image_index = 0)
{
    const char* imgext = EMUtil::get_imagetype_name(image_type);
    bool is_new_file = false;
    int err = 0;
    
    if (outfile == 0) {
		is_new_file = true;
		outfile = new char[256];
		strcpy(outfile, infile);
		const char* ext = strrchr(infile, '.');
		outfile[strlen(infile) - strlen(ext)] = '\0';
		sprintf(outfile, "%s_%d.%s", outfile, r_image_index, imgext);
    }

    //em->dump_data(outfile);
	try {
		em->write_image(outfile, w_image_index, image_type);
	}
	catch(...) {
		err = 1;
	}
	
    if (is_new_file) {
    	if( outfile )
    	{
			delete [] outfile;
			outfile = 0;
    	}
    }
    
    return err;
}


int test_image(const char* base_filename, int r_image_index = 0,
			   Region* area=0, bool is_3d = false,
			   int expected_err = 0, char* outfile = 0, 
			   EMUtil::ImageType image_type = EMUtil::IMAGE_UNKNOWN,
			   int w_image_index = 0)
{
    char home[128];
    sprintf(home, getenv("HOME"));
    char filename[256];
    sprintf(filename, "%s/images/%s", home, base_filename);
    int err = 0;

	try {
		EMData* em = new EMData();
		em->read_image(filename, r_image_index, false, area, is_3d);

		if (image_type ==  EMUtil::IMAGE_UNKNOWN) {
			image_type = EMUtil::IMAGE_MRC;
		}
		
		write_image(em, base_filename, outfile, r_image_index, image_type, w_image_index);    
		
		if( em )
		{
			delete em;
			em = 0;
		}
	}
	catch(...) {
		err = 1;
		LOGERR("image read/write error: %s\n", base_filename);
	}
	
    total_ntests++;
    if (err != expected_err) {
		fail_ntests++;
		err_code = 1;
    }
    
    return err;
}

int fail_test(const char* base_filename, int r_image_index = 0,
			  Region* area=0, bool is_3d = false, char* outfile = 0, 
			  EMUtil::ImageType image_type = EMUtil::IMAGE_UNKNOWN,
			  int w_image_index = 0)
{
    return test_image(base_filename, r_image_index, area, is_3d, 1,
					  outfile, image_type, w_image_index);
}

int pass_test(const char* base_filename, int r_image_index = 0,
			  Region* area=0, bool is_3d = false, char* outfile = 0, 
			  EMUtil::ImageType image_type = EMUtil::IMAGE_UNKNOWN,
			  int w_image_index = 0)
{
    return test_image(base_filename, r_image_index, area, is_3d, 0,
					  outfile, image_type, w_image_index);
}


int test_mrc()
{
    const char* file3d = "3d.mrc";
    Region good_3d1(0, 0, 0, 100, 100, 10);
    pass_test(file3d, 0, &good_3d1, false, "3d-all-1.mrc");

    Region good_3d2(0, 0, 24, 100, 100, 10);
    pass_test(file3d, 0, &good_3d2, false, "3d-all-2.mrc");

    char filename1[32];
    
    const int ids[] = {0, 20, 99};
    int n_ids = sizeof(ids) / sizeof(int);

    for (int i = 0; i < n_ids; i++) {
		Region d3_d2(0, 0, ids[i], 100, 100, 1);
		sprintf(filename1, "3d_each_%d.mrc", ids[i]);
		pass_test(file3d, 0, &d3_d2, false, filename1);
    }

    const char* file1d = "tablet.mrc";
    Region good_2d1(0, 1500, 300, 400);
    pass_test(file1d, 0, 0, false, "tablet_all1.mrc");
    pass_test(file1d, 0, &good_2d1, false, "tablet_good1.mrc");

    pass_test(file1d, 0, 0, false, "tablet_all1.mrc");
    pass_test(file3d, 0, 0, false, "3d-all-1.mrc");
    pass_test(file1d, 0, &good_2d1, false, "tablet_good1.mrc");

    Region good_3d3(0, 0, 0, 20, 20, 20);
    
    pass_test(file3d, 0, &good_3d3, false, "3d_good1.mrc");
    pass_test(file3d, 0, 0, false, "3d-all-1.mrc");
    
    Region bad_2d1(-2, -4, 10, 20);
    Region bad_2d2(1, 2, 3000, 400);
    Region bad_3d1(0, 1, 2, 400, 230, 5);
    Region bad_3d2(0, -3, -5, 3, 5, 9);
    
    fail_test(file1d, -3, 0, false, "tablet_no.mrc");
    fail_test(file1d, 0, &bad_2d1, false, "tablet_bad1.mrc");
    fail_test(file1d, 0, &bad_2d2, false, "tablet_bad2.mrc");
    
    fail_test(file3d, 40, 0, false, "3d-1.mrc");
    fail_test(file3d, 10, 0, false, "3d-all-2.mrc");

    fail_test(file3d, 120, 0, false, "3d_bad0.mrc");
    fail_test(file3d, 0, &bad_3d1, false, "3d_bad1.mrc");
    fail_test(file3d, 0, &bad_3d2, false, "3d_bad2.mrc");

    return 0;
}

int test_spider()
{
    pass_test("spider-single.spi");
    pass_test("spider-stack.spi", 0, 0, true);
    pass_test("spider-stack.spi", 27, 0, false, "spider-3d_27.mrc");

    Region r1(10, 10, 400, 800);
    pass_test("spider-single.spi", 0, &r1, false, "spider-single_r1.mrc");

    Region r2(0, 0, 40, 100);
    pass_test("spider-stack.spi", 27, &r2, false, "spider-3d_r2.mrc");
    
    pass_test("search.dm3", 0, 0, false, 0, EMUtil::IMAGE_SPIDER);
    pass_test("search.dm3", 0, 0, false, 0, EMUtil::IMAGE_SINGLE_SPIDER);

    pass_test("tablet.mrc", 0, 0, false, 0, EMUtil::IMAGE_SPIDER);
    pass_test("tablet.mrc", 0, 0, false, 0, EMUtil::IMAGE_SINGLE_SPIDER);

    for (int i = 0; i < 20; i++) {
		Region d32(0, 0, i, 100, 100, 1);
		pass_test("3d.mrc", 0, &d32, false, "3d_all.spi", EMUtil::IMAGE_SPIDER, i);
    }

    return err_code;
}


int test_dm3()
{
    Region good1(0, 0, 400, 400);
    Region good2(120, 230, 500, 600);
    Region bad1(0, 0, 1400, 400);
    pass_test("ccd.dm3");
    pass_test("search.dm3");
    
    pass_test("search.dm3", 0, &good1, false, "search_good1.mrc");
    pass_test("search.dm3", 0, &good2, false, "search_good2.mrc");
    fail_test("search.dm3", 0, &bad1, false, "search_bad1.mrc");
    
    pass_test("ccd.dm3", 0, &good1, false, "ccd_good1.mrc");
    pass_test("ccd.dm3", 0, &good2, false, "ccd_good2.mrc");

    return err_code;
}

int test_icos()
{
    pass_test("icos3f.map", 12);
    pass_test("icos2f.map", 12);
    pass_test("icos2f.map", 220);    
    pass_test("icos2f.map", 2, 0, true);
    fail_test("icos3f.map", 111);
    pass_test("icos3f.map", 0, 0, true);    
    pass_test("tablet.mrc", 0, 0, false, 0, EMUtil::IMAGE_ICOS);
    
    return err_code;
}

int test_tiff()
{
    pass_test("ducky-16bits.tif");
    pass_test("ducky-8bits.tif");
    pass_test("test2.tif");

    Region good1(1, 2, 12, 23);
    Region good2(1, 2, 220, 212);
    Region good3(1001, 1, 412, 8974);

    pass_test("ducky-16bits.tif", 0, &good1, false, "ducky16-good1.mrc");
    pass_test("ducky-8bits.tif", 0, &good2, false, "ducky8-good1.mrc");
    pass_test("test2.tif", 0, &good3, false, "test2_good3.mrc");
    
    Region bad1(1234, 0, 200, 300);
    Region bad2(123, 234, 1220, 212);

    fail_test("test2.tif", 1);
    fail_test("ducky-16bits.tif", 0, &bad1, false, "ducky16-bad1.mrc");
    fail_test("ducky-8bits.tif", 0, &bad2, false, "ducky8-bad2.mrc");
    
    return err_code;
}

int test_hdf()
{
    Region r1(0, 0, 800, 800);
    pass_test("search.h5", 0, &r1, false, "search_r1.mrc");

    Region r2(0, 0, 20, 100, 100, 1);
    pass_test("3d.h5", 0, &r2, false, "3d_r1.mrc");

    pass_test("t1.h5");
    pass_test("m1.h5");

    Region good1(0, 0, 0, 255, 255, 1);
    Region good2(0, 0, 0, 255, 255, 10);
    Region good3(0, 0, 40, 255, 255, 20);
    Region good4(10, 10, 0, 111, 121, 50);

    pass_test("m1.h5", 0, &good1, false, "m1_good1.mrc");
    pass_test("m1.h5", 0, &good2, false, "m1_good2.mrc");
    pass_test("m1.h5", 0, &good3, false, "m1_good3.mrc");
    pass_test("m1.h5", 0, &good4, false, "m1_good4.mrc");

    pass_test("bob_ctf.mrc", 0, 0, false, 0, EMUtil::IMAGE_HDF);
    pass_test("tablet.mrc", 0, 0, false, 0, EMUtil::IMAGE_HDF);
    pass_test("search.dm3", 0, 0, false, 0, EMUtil::IMAGE_HDF);
    pass_test("3d.mrc", 0, 0, false, 0, EMUtil::IMAGE_HDF);

    
    return err_code;
}

int test_pgm()
{
    pass_test("clinton.pgm", 0, 0, false, 0, EMUtil::IMAGE_IMAGIC);
    pass_test("bob.pgm");

    Region r1(20, 20, 200, 200);
    
    pass_test("clinton.pgm", 0, &r1, false, "clinton_1.mrc");
    pass_test("bob.pgm", 0, &r1, false, "bob_1.mrc");
    
    return err_code;
}

int test_lst()
{
    for (int k = 0; k < 5; k++) {
		pass_test("lst1.lst", k);
    }

    pass_test("cls0000.lst", 0);
    pass_test("cls0000.lst", 1);
    pass_test("cls0000.lst", 19);
    pass_test("cls0000.lst", 48);
    fail_test("cls0000.lst", -1);
    fail_test("cls0000.lst", 49);
    fail_test("cls0000.lst", 100);

    Region r1(20, 15, 55, 56);
    pass_test("cls0000.lst", 4);

    pass_test("cls0000.lst", 0, &r1, false, "cls_r1.mrc");

    return err_code;
}

int test_png()
{
    pass_test("window.png");
    pass_test("poker.png");
    pass_test("tablet.mrc",  0, 0, false, 0, EMUtil::IMAGE_PNG);
    pass_test("search.dm3",  0, 0, false, 0, EMUtil::IMAGE_PNG);
    
    return err_code;
}

int test_sal()
{
    pass_test("pds_tablet.hdr");

    Region r1(1, 1, 11850, 1200);
    Region r2(100, 200, 1500, 800);
    Region r3(100, 200, 11000, 200);
    
    pass_test("pds_tablet.hdr", 0, &r1, false, "pds_tablet_r1.mrc");
    pass_test("pds_tablet.hdr", 0, &r2, false, "pds_tablet_r2.mrc");
    pass_test("pds_tablet.hdr", 0, &r3, false, "pds_tablet_r3.mrc");
    
    return err_code;
}


int test_imagic()
{
#if 0
    pass_test("micam.hed", 12, 0, true);
    pass_test("ali.hed", 111);
    pass_test("classes.to.hed", 271);

    Region r1(3, 5, 50, 75);
    Region r2(12, 13, 31, 100, 100, 100);

    pass_test("classes.to.hed", 271, &r1, false, "classes_r1.mrc");
    pass_test("micam.hed", 0, &r2, true, "micam_r2.mrc");
    
    pass_test("start.hed", 100);
    pass_test("start.hed", 0, 0, true);

    pass_test("start.hed", 100);
    fail_test("start.hed", 500);
    fail_test("start.hed", 600);
    pass_test("tablet.mrc", 0, 0, false, 0, EMUtil::IMAGE_IMAGIC);
    pass_test("start.hed", 0, 0, true);

    pass_test("start.hed", 0, 0, true, 0, EMUtil::IMAGE_IMAGIC);
    pass_test("3d.mrc", 0, 0, false, 0, EMUtil::IMAGE_IMAGIC);
#endif
    
    EMData* e = new EMData();
    e->read_image("/home/lpeng/images/tablet.mrc");

    for (int i = 0; i < 2; i++) {
		e->write_image("tablet1.img", -1, EMUtil::IMAGE_IMAGIC);
    }
    
    if( e )
    {
    	delete  e;
    	e = 0;
    }
    
    return err_code;
}

int test_pif()
{    
    pass_test("sv-2d.pif");
    pass_test("sv-2d.pif", 19);
    fail_test("sv-2d.pif", -1);
    fail_test("sv-2d.pif", 158);
    pass_test("sv-3d.pif", 0);
    fail_test("sv-3d.pif", 1);
    fail_test("sv-3d.pif", -1);

    return err_code;
}

int test_performance()
{
    const char* imagefile = "/home/lpeng/raw_images/stress/start.lst";
    
    int nimg = EMUtil::get_image_count(imagefile);
    EMData* d = new EMData();
    double x_sum = 0;
    
    for (int i = 0; i < nimg; i++) {
		d->read_image(imagefile, i, true);
		Dict dict = d->get_attr_dict();
		int nx = dict.get("nx");
		x_sum += nx;
    }
    printf("nimg = %d, nx sum = %f\n", nimg, x_sum);
    
    if( d )
    {
    	delete d;
    	d = 0;
    }

    return 0;
}

    
    
void usage()
{
    printf("usage: imageio -vN dm3|tiff|hdf|pif|mrc|spider|pgm|lst|icos|png|sal|amira|gatan2|imagic\n");
}

int main(int argc, char* argv[])
{
    if (argc == 1) {
		usage();
		exit(1);
    }

    Util::set_log_level(argc, argv);
			     
    const char* imageformat = argv[argc-1];

    printf("Testing '%s' imageio\n\n", imageformat);
    
    if (strcmp(imageformat, "dm3") == 0) {
		test_dm3();
    }
    else if (strcmp(imageformat, "tiff") == 0) {
		test_tiff();
    }
    else if (strcmp(imageformat, "hdf") == 0) {
		test_hdf();
    }
    else if (strcmp(imageformat, "pif") == 0) {
		test_pif();
    }
    else if (strcmp(imageformat, "mrc") == 0) {
		test_mrc();
    }
    else if (strcmp(imageformat, "spider") == 0) {
		test_spider();
    }
    else if (strcmp(imageformat, "pgm") == 0) {
		test_pgm();
    }
    else if (strcmp(imageformat, "lst") == 0) {
		test_lst();
    }
    else if (strcmp(imageformat, "icos") == 0) {
		test_icos();
    }
    else if (strcmp(imageformat, "png") == 0) {
		test_png();
    }
    else if (strcmp(imageformat, "sal") == 0) {
		test_sal();
    }
    else if (strcmp(imageformat, "amira") == 0) {
		pass_test("tablet.mrc", 0, 0, false, 0, EMUtil::IMAGE_AMIRA);
    }
    else if (strcmp(imageformat, "gatan2") == 0) {
		pass_test("gatan2.dm2");
    }
    else if (strcmp(imageformat, "imagic") == 0) {
		test_imagic();
    }
    else if (strcmp(imageformat, "perf") == 0) {
		test_performance();
    }
    else {
		usage();
		exit(1);
    }
    
    printf("Total # Tests: %d. Failed # Tests: %d\n", total_ntests, fail_ntests);
    
    return err_code;
}
