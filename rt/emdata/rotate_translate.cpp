#include "emdata.h"
#include "testutil.h"
//#include "exception.h"
#include <iostream>

using namespace EMAN;
using namespace std;

string get_test_image()
{
	return TestUtil::get_debug_image("tg.mrc");
}

void r1()
{
	EMData *image = new EMData();

	image->read_image(get_test_image());
	image->rotate(1.0329837512591338,3.7260642381912579,5.7671541529246966);
	image->write_image("r1.mrc");
	if( image )
	{
		delete image;
		image = 0;
	}
}

void t1()
{
	EMData *image = new EMData();

	image->read_image(get_test_image());
	image->rotate_translate(0,0,0,16,16,16);
	image->write_image("t1.mrc");
	if( image )
	{
		delete image;
		image = 0;
	}
}

void t2()
{
	EMData *image = new EMData();

	image->read_image(get_test_image());
	image->translate(-0.5f,-0.5f, 0.0f);
	image->write_image("t2.mrc");
	if( image )
	{
		delete image;
		image = 0;
	}
}

void rt1()
{
	EMData *image = new EMData();

	image->read_image(get_test_image());
	image->rotate_translate(1.0329837512591338,3.7260642381912579,5.7671541529246966,16,16,16);
	image->write_image("rt1.mrc");
	if( image ) 
	{
		delete image;
		image = 0;
	}
}

void rt2()
{
	EMData *image = new EMData();

	image->read_image(get_test_image());
	Transform3D r = Transform3D(Vec3f(16,16,16), 1.0329837512591338,3.7260642381912579,
							5.7671541529246966);
	image->rotate_translate(r);
	image->write_image("rt2.mrc");
	if( image )
	{
		delete image;
		image = 0;
	}
}

void compare_image()
{
	EMData *image1 = new EMData();
	EMData *image2 = new EMData();
	
	image1->read_image(get_test_image());
	image2->read_image(get_test_image());
	
	image1->translate(Vec3f(100.5, 100.5, 0));
	image1->write_image("tran1.mrc");
	
	Transform3D r = Transform3D(Vec3f(100.5, 100.5, 0), 0,0,0);
	image2->rotate_translate(r);
	image2->write_image("tran2.mrc");
	
	float *data1 = image1->get_data();
	float *data2 = image2->get_data();
	int size = sizeof(float) * image1->get_xsize() * image1->get_ysize() * image1->get_zsize();
	int cmp = memcmp( (void*)data1, (void*)data2, size );
	cout << "Compare reault is: " << cmp << endl;
	
	if( image1 ) {
		delete image1;
		image1 = 0;
	}
	if( image2 ) {
		delete image2;
		image2 = 0;
	}
}

void debug_align()
{
	cout << "Enter the function debug_align()... " << endl;
	
	EMData *a = new EMData();
	
	cout << "can we go there?" << endl;
	a->set_size(128, 128, 1);
	a->process("testimage.noise.uniform.rand");
	//a->write_image("scurve.mrc");
/*	
	try {
		printf("Dic of a is: \n");
		Dict dic3 = a->get_attr_dict();
		EMUtil::dump_dict(dic3);
	}
	catch(_NotExistingObjectException& e) {
		cout << "catch an _NotExistingObjectException in dump_dict()..." << endl;
		cout << e.what();
	}
	catch(E2Exception& e) {
		cout << "catch an E2Exception in dump_dict()..." << endl;
		cout << e.what();
	}
	catch(...) {
		cout << "catch unknown exception in dump_dict()..." << endl;
	}
*/	
	EMData* b=a->rot_scale_trans2D(0.0f,1.0f,1.0f,0.0f);
/*	
	try {
		printf("Dic of b is: \n");
		Dict dic2 = b->get_attr_dict();
		EMUtil::dump_dict(dic2);
		
		printf("Dic of a is: \n");
		Dict dic3 = a->get_attr_dict();
		EMUtil::dump_dict(dic3);
	}
	catch(_NotExistingObjectException& e) {
		cout << "catch an _NotExistingObjectException in dump_dict()..." << endl;
		cout << e.what();
	}
	catch(E2Exception& e) {
		cout << "catch an E2Exception in dump_dict()..." << endl;
		cout << e.what();
	}
	catch(...) {
		cout << "catch unknown exception in dump_dict()..." << endl;
	}
*/	
	
	EMData* d=b->align("rotate_translate",a,Dict(),"variance",Dict());
	
	printf("Dic of d is: \n");
	Dict dic1 = d->get_attr_dict();
	EMUtil::dump_dict(dic1);
	
	try {
		printf("Dic of b is: \n");
		Dict dic2 = b->get_attr_dict();
		EMUtil::dump_dict(dic2);
		
		printf("Dic of a is: \n");
		Dict dic3 = a->get_attr_dict();
		EMUtil::dump_dict(dic3);
	}
	catch(_NotExistingObjectException& e) {
		cout << "catch an _NotExistingObjectException in dump_dict()..." << endl;
		cout << e.what();
	}
	catch(E2Exception& e) {
		cout << "catch an E2Exception in dump_dict()..." << endl;
		cout << e.what();
	}
	catch(...) {
		cout << "catch unknown exception in dump_dict()..." << endl;
	}
	
	try {
		delete d;
		delete b;
		delete a;
	}
	catch( E2Exception & e ) {
		cout << "catch an E2Exception." << endl;
		e.what();
	}
	catch(...) {
		cout << "catch unknown exception in delete..." << endl;
	}
	
	cout << "Finish debug_align()... " << endl;
}

void debug_log()
{
	cout << "try to debug memory issue in log." << endl;
	
	LOGERR("Simple log test.");
}

//no memory leak
void debug_set_size()
{
	cout << "strating degugging set_size()." << endl;
	EMData *a = new EMData();
	a->set_size(1024,1024,1);
	delete a;
	
	LOGERR("test finished.");
}

void debug_footprint()
{
	cout << "begin debug_footprint()" << endl;
	
	EMData *a = new EMData();
	a->set_size(128, 128, 1);
	a->process("testimage.noise.uniform.rand");
	
	a->make_rotational_footprint(true);
	
	printf("Dic of a is: \n");
	Dict dic3 = a->get_attr_dict();
	EMUtil::dump_dict(dic3);
	
	try {
		delete a;
	}
	catch( E2Exception & e ) {
		cout << e.what();
	}
	catch(...) {
		cout << "catch unknown exception..." << endl;
	}
	
	cout << "end debug_footprint()" << endl;
}

void debug_get_clip()
{
	cout << "begin debug_get_clip()" << endl;
	
	EMData *a = new EMData();
	a->set_size(128, 128, 1);
	a->process("testimage.noise.uniform.rand");
//	a->write_image("rand3.mrc");
/*	
	try {
		cout << "At beginning" << endl;
		printf("Dic of a is: \n");
		Dict dic2 = a->get_attr_dict();
		EMUtil::dump_dict(dic2);
	}
	catch(_NotExistingObjectException& e) {
		std::cout << "catch an _NotExistingObjectException in dump_dict()..." << std::endl;
		std::cout << e.what();
	}
	catch(E2Exception& e) {
		std::cout << "catch an E2Exception in dump_dict()..." << std::endl;
		std::cout << e.what();
	}
	catch(...) {
		std::cout << "catch unknown exception in dump_dict()..." << std::endl;
	}
*/	
	int nx = a->get_xsize();
	int ny = a->get_ysize();
	int cs = (((nx * 7 / 4) & 0xfffff8) - nx) / 2;
	Region r1;
	r1 = Region(-cs, -cs, nx + 2 * cs, ny + 2 * cs);
	std::cout << "The region r1 is: " << r1.get_string() << std::endl;
	
	try {
		cout << "Before get_clip()" << endl;
		printf("Dic of a is: \n");
		Dict dic2 = a->get_attr_dict();
		EMUtil::dump_dict(dic2);
	}
	catch(_NotExistingObjectException& e) {
		std::cout << "catch an _NotExistingObjectException in dump_dict()..." << std::endl;
		std::cout << e.what();
	}
	catch(E2Exception& e) {
		std::cout << "catch an E2Exception in dump_dict()..." << std::endl;
		std::cout << e.what();
	}
	catch(...) {
		std::cout << "catch unknown exception in dump_dict()..." << std::endl;
	}
	
	
	EMData *tmp2 = 0;
	tmp2 = a->get_clip(r1);
	
	try {
		cout << "After get_clip()" << endl;
		printf("Dic of a is: \n");
		Dict dic2 = a->get_attr_dict();
		EMUtil::dump_dict(dic2);
	}
	catch(_NotExistingObjectException& e) {
		std::cout << "catch an _NotExistingObjectException in dump_dict()..." << std::endl;
		std::cout << e.what();
	}
	catch(E2Exception& e) {
		std::cout << "catch an E2Exception in dump_dict()..." << std::endl;
		std::cout << e.what();
	}
	catch(...) {
		std::cout << "catch unknown exception in dump_dict()..." << std::endl;
	}
	
//	tmp2->write_image("rand_enlarge.mrc");
	
	try {
		delete a;
		a = 0;
		delete tmp2;
		tmp2 = 0;
	}
	catch( E2Exception & e ) {
		cout << e.what();
	}
	catch(...) {
		cout << "catch unknown exception..." << endl;
	}

	cout << "end debug_get_clip()" << endl;
}

int main()
{
	cout << "Starting to test rotate_translate..." << endl;
	
	try {
		//r1();
		//t1();
		//rt1();
		//rt2();
		//t2();
		//compare_image();
		debug_align();
		//debug_log();
		//debug_set_size();
		//debug_footprint();
		//debug_get_clip();
	}
	catch (E2Exception & e) {
		cout << e.what();
	}
	
	
	return 0;
}
