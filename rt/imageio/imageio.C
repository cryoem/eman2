#include <stdio.h>
#include <stdlib.h>
#include "emdata.h"
#include "util.h"
#include "transform.h"
#include "emutil.h"
#include "log.h"


using namespace EMAN;

int write_image(EMData* em, const char* infile, char* outfile,
		int r_image_index, EMUtil::ImageType image_type,
		int w_image_index = 0)
{
    const char* imgext = EMUtil::get_imagetype_name(image_type);
    bool is_new_file = false;
    
    if (outfile == 0) {
	is_new_file = true;
	outfile = new char[256];
	strcpy(outfile, infile);
	char* ext = strrchr(infile, '.');
	outfile[strlen(infile) - strlen(ext)] = '\0';
	sprintf(outfile, "%s_%d.%s", outfile, r_image_index, imgext);
    }

    em->dump_data(outfile);
    em->write_image(outfile, w_image_index, image_type);

    if (is_new_file) {
	delete [] outfile;
	outfile = 0;
    }
    
    return 0;
}


int test_image(const char* base_filename, int r_image_index = 0,
	       Region* area=0, bool is_3d = false, char* outfile = 0,
	       EMUtil::ImageType image_type = EMUtil::IMAGE_UNKNOWN,
	       int w_image_index = 0)
{
    char home[128];
    sprintf(home, getenv("HOME"));
    char filename[256];
    sprintf(filename, "%s/images/%s", home, base_filename);

    EMData* em = new EMData();
    if (em->read_image(filename, r_image_index, false, area, is_3d) != 0) {
	fprintf(stderr, "EMData read_image() failed\n");
	delete em;
	return 1;
    }

    
    map<string, EMObject> dict = em->get_attr_dict();
    
    if (image_type ==  EMUtil::IMAGE_UNKNOWN) {
	image_type = EMUtil::IMAGE_MRC;
    }
    
    write_image(em, base_filename, outfile, r_image_index, image_type, w_image_index);    
    
    delete em;
    em = 0;

    return 0;
}

int test_mrc()
{

    test_image("bob_0.MRC");
    
#if 0
    const char* file3d = "3d.mrc";
    Region good_3d1(0, 0, 0, 100, 100, 10);
    test_image(file3d, 0, &good_3d1, false, "3d-all-1.mrc");

    Region good_3d2(0, 0, 24, 100, 100, 10);
    test_image(file3d, 0, &good_3d2, false, "3d-all-2.mrc");

    char filename1[32];
    
    for (int i = 0; i < 100; i++) {
	Region d3_d2(0, 0, i, 100, 100, 1);
	sprintf(filename1, "3d_each_%d.mrc", i);
	test_image(file3d, 0, &d3_d2, false, filename1);
    }
    

    const char* file1d = "tablet.mrc";
    Region good_2d1(0, 1500, 300, 400);
    test_image(file1d, 0, 0, false, "tablet_all1.mrc");
    test_image(file1d, 0, &good_2d1, false, "tablet_good1.mrc");


    Region bad_2d1(-2, -4, 10, 20);
    Region bad_2d2(1, 2, 3000, 400);
    
    test_image(file1d, 0, 0, false, "tablet_all1.mrc");
    test_image(file3d, 0, 0, true, "3d-all-1.mrc");
    
    test_image(file1d, 0, &good_2d1, false, "tablet_good1.mrc");
    
    test_image(file1d, -3, 0, false, "tablet_no.mrc");
    test_image(file1d, 0, &bad_2d1, false, "tablet_bad1.mrc");
    test_image(file1d, 0, &bad_2d2, false, "tablet_bad2.mrc");
    
    Region good_3d1(0, 0, 0, 20, 20, 20);
    Region bad_3d1(0, 1, 2, 400, 230, 5);
    Region bad_3d2(0, -3, -5, 3, 5, 9);

    // positive tests
    test_image(file3d, 40, 0, false, "3d-1.mrc");
    test_image(file3d, 0, 0, true, "3d-all-1.mrc");
    
    test_image(file3d, 10, 0, true, "3d-all-2.mrc");
    test_image(file3d, 0, &good_3d1, true, "3d_good1.mrc");

    // negative tests
    test_image(file3d, 120, 0, false, "3d_bad0.mrc");
    test_image(file3d, 0, &bad_3d1, true, "3d_bad1.mrc");
    test_image(file3d, 0, &bad_3d2, true, "3d_bad2.mrc");
#endif
    return 0;
}

int test_spider()
{
    int err = 0;

    err = test_image("tablet_0.SPIDER");
    
#if 0
    err = test_image("spider-single.spi");
    err = test_image("spider-stack.spi", 0, 0, true);
    err = test_image("spider-stack.spi", 27, 0, false, "spider-3d_27.mrc");

    Region r1(10, 10, 400, 800);
    err = test_image("spider-single.spi", 0, &r1, false, "spider-single_r1.mrc");

    Region r2(0, 0, 40, 100);
    err = test_image("spider-stack.spi", 27, &r2, false, "spider-3d_r2.mrc");
    
    err = test_image("search.dm3", 0, 0, false, 0, EMUtil::IMAGE_SPIDER);
    err = test_image("search.dm3", 0, 0, false, 0, EMUtil::IMAGE_SINGLE_SPIDER);
#endif
#if 0
    err = test_image("tablet.mrc", 0, 0, false, 0, EMUtil::IMAGE_SPIDER);
    err = test_image("tablet.mrc", 0, 0, false, 0, EMUtil::IMAGE_SINGLE_SPIDER);

    for (int i = 0; i < 100; i++) {
	Region d32(0, 0, i, 100, 100, 1);
	err = test_image("3d.mrc", 0, &d32, false, "3d_all.spi", EMUtil::IMAGE_SPIDER, i);
    }
#endif
    return err;
}


int test_dm3()
{
    int err = 0;

    Region good1(0, 0, 400, 400);
    Region good2(120, 230, 500, 600);
    
    Region bad1(0, 0, 1400, 400);
    
    err = test_image("search.dm3");
    err = test_image("search.dm3", 0, &good1, false, "search_good1.mrc");
    err = test_image("search.dm3", 0, &good2, false, "search_good2.mrc");
    err = test_image("search.dm3", 0, &bad1, false, "search_bad1.mrc");
    
    err = test_image("ccd.dm3");
    err = test_image("ccd.dm3", 0, &good1, false, "ccd_good1.mrc");
    err = test_image("ccd.dm3", 0, &good2, false, "ccd_good2.mrc");
    
    return err;
}

int test_icos()
{
    int err = 0;
#if 0
    err = test_image("icos2f.map", 220);
    err = test_image("icos2f.map", 215);    
    err = test_image("icos2f.map", 2, 0, true);
#endif
    
    err = test_image("tablet.mrc", 0, 0, false, 0, EMUtil::IMAGE_ICOS);

    //    err = test_image("3d.mrc", 0, 0, false, 0, EMUtil::IMAGE_ICOS);

    
    return err;
}

int test_tiff()
{
    int err = 0;

    err = test_image("ducky-16bits.tif");
    err = test_image("ducky-8bits.tif");
    err = test_image("test2.tif");

    Region good1(1, 2, 12, 23);
    Region good2(1, 2, 220, 212);
    Region good3(1001, 1, 412, 8974);

    err = test_image("ducky-16bits.tif", 0, &good1, false, "ducky16-good1.mrc");
    err = test_image("ducky-8bits.tif", 0, &good2, false, "ducky8-good1.mrc");
    err = test_image("test2.tif", 0, &good3, false, "test2_good3.mrc");
    
    Region bad1(1234, 0, 200, 300);
    Region bad2(123, 234, 1220, 212);

    err = test_image("test2.tif", 1);
    err = test_image("ducky-16bits.tif", 0, &bad1, false, "ducky16-bad1.mrc");
    err = test_image("ducky-8bits.tif", 0, &bad2, false, "ducky8-bad2.mrc");
    
    return err;
}

int test_hdf()
{
    int err = 0;

    Region r1(0, 0, 800, 800);
    err = test_image("search.h5", 0, &r1, false, "search_r1.mrc");

    Region r2(0, 0, 20, 100, 100, 1);
    err = test_image("3d.h5", 0, &r2, false, "3d_r1.mrc");
#if 1
    err = test_image("t1.h5");
    err = test_image("m1.h5");

    Region good1(0, 0, 0, 255, 255, 1);
    Region good2(0, 0, 0, 255, 255, 10);
    Region good3(0, 0, 40, 255, 255, 20);
    Region good4(10, 10, 0, 111, 121, 50);

    err = test_image("m1.h5", 0, &good1, false, "m1_good1.mrc");
    err = test_image("m1.h5", 0, &good2, false, "m1_good2.mrc");
    err = test_image("m1.h5", 0, &good3, false, "m1_good3.mrc");
    err = test_image("m1.h5", 0, &good4, false, "m1_good4.mrc");

    err = test_image("tablet.mrc", 0, 0, false, 0, EMUtil::IMAGE_HDF);
    err = test_image("search.dm3", 0, 0, false, 0, EMUtil::IMAGE_HDF);
    err = test_image("3d.mrc", 0, 0, false, 0, EMUtil::IMAGE_HDF);
#endif
    
    return err;
}

int test_pgm()
{
    int err = 0;
#if 1
    err = test_image("clinton.pgm", 0, 0, false, 0, EMUtil::IMAGE_IMAGIC);
    err = test_image("bob.pgm");
#endif
    Region r1(20, 20, 200, 200);
    
    err = test_image("clinton.pgm", 0, &r1, false, "clinton_1.mrc");
    err = test_image("bob.pgm", 0, &r1, false, "bob_1.mrc");

    
    return err;
}

int test_lst()
{
    int err = 0;
#if 1
    const char* imagefile = "/home/lpeng/raw_images/stress/lst2.lst";
    
    int nimg = EMUtil::get_image_count(imagefile);
    EMData* d = new EMData();
    double x_sum = 0;
    
    for (int i = 0; i < nimg; i++) {
	d->read_image(imagefile, i, true);
	map<string, EMObject> dict = d->get_attr_dict();
	int nx = dict["nx"].get_int();
	x_sum += nx;
	printf("%i ",i);
    }
    printf("nimg = %d, nx sum = %f\n", nimg, x_sum);
    
    delete d;
    d = 0;
#endif
#if 0
    for (int i = 0; i < 5; i++) {
	err = test_image("lst1.lst", i);
    }

    err = test_image("cls0000.lst", 0);

    err = test_image("cls0000.lst", -1);
    err = test_image("cls0000.lst", 1);
    err = test_image("cls0000.lst", 19);
    err = test_image("cls0000.lst", 48);
    err = test_image("cls0000.lst", 49);
    err = test_image("cls0000.lst", 100);

    Region r1(20, 15, 55, 56);
    
    err = test_image("cls0000.lst", 0, &r1, false, "cls_r1.mrc");
#endif
    return err;
}

int test_png()
{
    int err = 0;
    
    err = test_image("window.png");
    err = test_image("poker.png");
    err = test_image("tablet.mrc",  0, 0, false, 0, EMUtil::IMAGE_PNG);
    err = test_image("search.dm3",  0, 0, false, 0, EMUtil::IMAGE_PNG);

    return err;
}

int test_sal()
{
    int err = 0;
    err = test_image("pds_tablet.hdr");

    Region r1(1, 1, 11850, 1200);
    Region r2(100, 200, 1500, 800);
    Region r3(100, 200, 11000, 200);
    
    err = test_image("pds_tablet.hdr", 0, &r1, false, "pds_tablet_r1.mrc");
    err = test_image("pds_tablet.hdr", 0, &r2, false, "pds_tablet_r2.mrc");
    err = test_image("pds_tablet.hdr", 0, &r3, false, "pds_tablet_r3.mrc");
    
    return err;
}


int test_imagic()
{
    int err = 0;

    err = test_image("micam.hed", 12, 0, true);
    err = test_image("ali.hed", 111);
    err = test_image("classes.to.hed", 271);

    Region r1(3, 5, 50, 75);
    Region r2(12, 13, 31, 100, 100, 100);

    err = test_image("classes.to.hed", 271, &r1, false, "classes_r1.mrc");
    err = test_image("micam.hed", 0, &r2, true, "micam_r2.mrc");
    
    err = test_image("start.hed", 100);
    err = test_image("start.hed", 0, 0, true);

    err = test_image("start.hed", 100);
    err = test_image("start.hed", 500);
    err = test_image("start.hed", 600);
    err = test_image("tablet.mrc", 0, 0, false, 0, EMUtil::IMAGE_IMAGIC);
    err = test_image("start.hed", 0, 0, true);

    err = test_image("start.hed", 0, 0, true, 0, EMUtil::IMAGE_IMAGIC);
    err = test_image("3d.mrc", 0, 0, false, 0, EMUtil::IMAGE_IMAGIC);

    return err;
}

int test_performance()
{
    const char* imagefile = "/home/lpeng/raw_images/stress/start.lst";
    
    int nimg = EMUtil::get_image_count(imagefile);
    EMData* d = new EMData();
    double x_sum = 0;
    
    for (int i = 0; i < nimg; i++) {
	d->read_image(imagefile, i, true);
	map<string, EMObject> dict = d->get_attr_dict();
	int nx = dict["nx"].get_int();
	x_sum += nx;
	//printf("%i ",i);
    }
    printf("nimg = %d, nx sum = %f\n", nimg, x_sum);
    
    delete d;
    d = 0;

    return 0;
}

    
    
void usage()
{
    printf("usage: imageio -vN dm3|tiff|hdf|pif|mrc|spi|pgm|lst|icos|png|sal|amira|gatan2|imagic\n");
}

int main(int argc, char* argv[])
{
    if (argc == 1) {
	usage();
	exit(1);
    }

    Log::LogLevel log_level = Log::WARNING_LOG;
    
    if (strncmp(argv[1], "-v", 2) == 0) {
	char str1[32];
	strcpy(str1, argv[1]+2);
	log_level = (Log::LogLevel) atoi(str1);
    }
    
    Log::logger()->set_level(log_level);
			     
    const char* imageformat = argv[argc-1];
    int err = 0;
    
    if (strcmp(imageformat, "dm3") == 0) {
	err = test_dm3();
    }
    else if (strcmp(imageformat, "tiff") == 0) {
	err = test_tiff();
    }
    else if (strcmp(imageformat, "hdf") == 0) {
	err = test_hdf();
    }
    else if (strcmp(imageformat, "pif") == 0) {
	err = test_image("sv-2d.pif");
	err = test_image("sv-2d.pif", 19);

	err = test_image("sv-2d.pif", -1);
	
	err = test_image("sv-2d.pif", 158);
	
	err = test_image("sv-3d.pif", 0);
	err = test_image("sv-3d.pif", 1);
	err = test_image("sv-3d.pif", -1);
    }
    else if (strcmp(imageformat, "mrc") == 0) {
	err = test_mrc();
    }
    else if (strcmp(imageformat, "spi") == 0) {
	err = test_spider();
    }
    else if (strcmp(imageformat, "pgm") == 0) {
	err = test_pgm();
    }
    else if (strcmp(imageformat, "lst") == 0) {
	err = test_lst();
    }
    else if (strcmp(imageformat, "icos") == 0) {
	err = test_icos();
	
    }
    else if (strcmp(imageformat, "png") == 0) {
	err = test_png();
    }
    else if (strcmp(imageformat, "sal") == 0) {
	err = test_sal();
    }
    else if (strcmp(imageformat, "amira") == 0) {
	err = test_image("tablet.mrc", 0, 0, false, 0, EMUtil::IMAGE_AMIRA);
    }
    else if (strcmp(imageformat, "gatan2") == 0) {
	err = test_image("gatan2.dm2");
    }
    else if (strcmp(imageformat, "imagic") == 0) {
	err = test_imagic();
    }
    else if (strcmp(imageformat, "perf") == 0) {
	err = test_performance();
    }
    else {
	usage();
	exit(1);
    }
    
    return err;
}



