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

//void rt2()
//{
//	EMData *image = new EMData();

//	image->read_image(get_test_image());
//	Transform3D r = Transform3D(Vec3f(16,16,16), 1.0329837512591338,3.7260642381912579,
//							5.7671541529246966);
//	image->rotate_translate(r);
//	image->write_image("rt2.mrc");
//	if( image )
//	{
//		delete image;
//		image = 0;
//	}
//}

//void compare_image()
//{
//	EMData *image1 = new EMData();
//	EMData *image2 = new EMData();
	
//	image1->read_image(get_test_image());
//	image2->read_image(get_test_image());
	
//	image1->translate(Vec3f(100.5, 100.5, 0));
//	image1->write_image("tran1.mrc");
	
//	Transform3D r = Transform3D(Vec3f(100.5, 100.5, 0), 0,0,0);
//	image2->rotate_translate(r);
//	image2->write_image("tran2.mrc");
	
//	float *data1 = image1->get_data();
//	float *data2 = image2->get_data();
//	int size = sizeof(float) * image1->get_xsize() * image1->get_ysize() * image1->get_zsize();
//	int cmp = memcmp( (void*)data1, (void*)data2, size );
//	cout << "Compare reault is: " << cmp << endl;
	
//	if( image1 ) {
//		delete image1;
//		image1 = 0;
//	}
//	if( image2 ) {
//		delete image2;
//		image2 = 0;
//	}
//}

void debug_align()
{
	cout << "Enter the function debug_align()... " << endl;
	
	EMData *a = new EMData();
	
	cout << "can we go there?" << endl;
	a->set_size(128, 128, 1);
	a->process_inplace("testimage.noise.uniform.rand");
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
	
	EMData* d=b->align("rotate_translate",a,Dict(),"sqeuclidean",Dict());
	
	printf("Dic of d is: \n");
	Dict dic1 = d->get_attr_dict();
	EMUtil::dump_dict(dic1);
	
	try {
		printf("Dic of b is: \n");
		Dict dic2 = b->get_attr_dict();
		cout << "prepare to dump" << endl;
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
	a->process_inplace("testimage.noise.uniform.rand");
	
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
	a->process_inplace("testimage.noise.uniform.rand");
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

void debug_write_image()
{
	EMData *a = new EMData();
	a->set_size(128, 128, 1);
	a->process_inplace("testimage.noise.uniform.rand");
	a->write_image("rand3.mrc");
	Dict dic = a->get_attr_dict();
	
	try {
		printf("Dic of a is: \n");
		EMUtil::dump_dict(dic);
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
}

void debeg_dict()
{
	cout << "Enter debug_dict()..." << endl;
	
	Dict dict;
	cout << "try LHS []..." << endl;
	dict["first"] = "one";
	
	cout << "try RHS []..." << endl;
	float tt = (float)(EMObject)dict["second"];
	
	cout << "Undifined dict['second'] = " << tt << endl; 
	
	const Dict dict2("first", "one", "second", 1.1f);
	float ttt = dict2["second"];
	dict2["forth"] = 4;
	cout << "const dict2['second'] = " << ttt << endl;
	
	float tttt = dict2["third"];
	cout << "undifined const dict2['third'] = " << tttt << endl;
	
	cout << "Leave debug_dict()..." << endl;
}

void debug_complex_image_arithmetic()
{
	cout << "start to test debug_complex_image_arithmetic()..." << endl;
	
	EMData * e = new EMData();
	e->set_size(75,75,1);
	e->process_inplace("testimage.circlesphere", Dict("radius",20,"fill","yes"));
	e->write_image("origin.mrc");
	EMData * eft = e->do_fft();
	
	float * data1 = eft->get_data();
	cout << "------------------------------------------------------------------" << endl;
	cout << "eft's data is: " << endl;
	for(int i=0; i<10; i++)
	{
 		cout << data1[i] << "\t";
	}
	cout << endl << "------------------------------------------------------------------" << endl;
	
	EMData * delta = new EMData();
	delta->set_size(75,75,1);
	delta->to_zero();
	delta->set_value_at(0,0,1.0);
	EMData * deltaft = delta->do_fft();
	
	float * data = deltaft->get_data();
	cout << "------------------------------------------------------------------" << endl;
	cout << "deltaft's data is: " << endl;
	for(int i=0; i<10; i++)
	{
 		cout << data[i] << "\t";
	}
	cout << endl << "------------------------------------------------------------------" << endl;
	
	EMData * divft = (*eft) / (*deltaft);
	
	float * data2 = divft->get_data();
	cout << "------------------------------------------------------------------" << endl;
	cout << "div's data is: " << endl;
	for(int i=0; i<10; i++)
	{
 		cout << data2[i] << "\t";
	}
	cout << endl << "------------------------------------------------------------------" << endl;
	
	EMData * div = divft->do_ift();
	div->write_image("div.mrc");
	
	if( e != 0 )
	{
		cout << "Delete e ..." << endl;
		delete e;
		e = 0;
	}
	if( eft != 0 )
	{
		cout << "Delete eft ..." << endl;
		delete eft;
		eft = 0;
	}
	if( delta != 0 )
	{
		cout << "Delete delta ..." << endl;
		delete delta;
		delta = 0;
	}
	if( deltaft != 0 )
	{
		cout << "Delete deltaft ..." << endl;
		delete deltaft;
		deltaft = 0;
	}
	
	cout << "Leave test debug_complex_image_arithmetic()..." << endl;
}

void debug_big_file()
{
	cout << "Start debug_big_file()..." << endl << endl;
	
	const string filename = "big_file.mrc";
	const Region* reg = new Region(0,0,0,10,10,10);
	
	EMData * e = new EMData();
	e->read_image(filename, 0, true, reg);
	
	cout << "End debug_big_file()..." << endl << endl;
}

void debug_insert_clip()
{
	cout << "Starting debug_insert_clip() ... " << endl;
	
	EMData * e = new EMData();
	e->set_size(4,4,4);
	e->to_zero();
	
	EMData * e2 = new EMData();
	e2->set_size(2,2,2);
	e2->to_one();
	
	e->insert_clip(e2, IntPoint(1,1,1));
	MArray3D data3 = e->get_3dview();
	for(int k=0; k<4; ++k)
	{
		for(int j=0; j<4; ++j)
		{
			for(int i=0; i<4; ++i)
			{
				cout << data3[i][j][k] << " ";
			}
			cout << endl;
		}
		cout << endl << endl;
	}
}

void debug_do_fft()
{
	cout << "Starting debug_do_fft() ... " << endl;
	
	EMData * e = new EMData();
	e->set_size(32, 32, 32);
	e->process_inplace("testimage.noise.uniform.rand");
	
	EMData * e2 = e->do_fft();
	EMData * e3 = e2->do_fft();
	EMData * e4 = e3->do_fft();
	cout << "e2 = " << e2 << endl
		 << "e3 = " << e3 << endl
		 << "e4 = " << e4 << endl;
}

void debug_do_fft_inplace()
{
	cout << "Starting debug_do_fft_inplace() ... " << endl;

	EMData * e = new EMData();
	e->set_size(32, 32, 32);
	e->process_inplace("testimage.noise.uniform.rand");

	e->do_fft_inplace();
	
	if(e) {
		delete e;
		e = 0;
	}
	
}
/*
void debug_calc_radial_dist()
{
	cout << "Starting debug_calc_radial_dist() ... " << endl;
	
	EMData * e = new EMData();
	e->set_size(32, 32, 1);
	e->process_inplace("testimage.noise.uniform.rand");
	
	int ny = e->get_ysize();
	vector<float> v;
	v = e->calc_radial_dist(ny, 0, 0.5);
	
	vector<float>::iterator iter;
	for(iter=v.begin(); iter!=v.end(); iter++)
	{
		cout << *iter << endl;
	}
	
	cout << "Finish debug_calc_radial_dist() ... " << endl;
}
*/
void debug_calc_hist()
{
	cout << "Starting debug_calc_hist() ... " << endl;
	
	EMData * e = new EMData();
	e->set_size(32, 32, 1);
	e->process_inplace("testimage.noise.uniform.rand");
	
	vector<float> v;
	v = e->calc_hist();	
}

void debug_set_value_at()
{
	cout << "Starting debug_set_value_at() ... " << endl;
	
	EMData * e = new EMData();
	e->set_size(32, 32, 32);
	e->process_inplace("testimage.noise.uniform.rand");
	
	e->set_value_at(1, 1, 1, 1.0);
	float f = e->get_value_at(1, 1, 1);
	cout << "f = " << f << endl;
	
	e->set_value_at(1000000,2,2,122.0);
	cout << e->get_value_at(1000000, 2, 2) << endl;
}

void debug_common_lines()
{
	cout << "Starting debug_common_lines() ..." << endl;
	
	EMData * e = new EMData();
	e->set_size(32,32,1);
	e->process_inplace("testimage.noise.uniform.rand");
	
	EMData * e2 = new EMData();
	e2->set_size(32,32,1);
	e2->process_inplace("testimage.noise.uniform.rand");
	
	EMData * e3 = new EMData();
	e3->set_size(32,32,1);
	e3->process_inplace("testimage.noise.uniform.rand");
	
//	EMData * e4 = e2->do_fft();
//	EMData * e5 = e3->do_fft();
	
//	cout << "before common_lines() call ..." << endl;
//	e->common_lines(e4, e5);
//	cout << "after common_lines() call ..." << endl;
	
	cout << "before common_lines_real() call ..." << endl;
	e->common_lines_real(e2, e3);
	cout << "after common_lines_real() call ..." << endl;
}

void debug_sigma_processor()
{
	cout << "Start debug_sigma_processor function ..." << endl;
	
	EMData * e = new EMData();
	e->set_size(2,2,1);
	e->set_value_at(0, 0, 1.0);
	e->set_value_at(0, 1, 2.0);
	e->set_value_at(1, 0, 3.0);
	e->set_value_at(1, 1, 4.0);
	
	float mean = e->get_attr("mean");
	float sigma = e->get_attr("sigma");
	
	cout << "mean = " << mean << endl;
	cout << "sigma = " << sigma << endl;
	
	MArray2D data0 = e->get_2dview();
	for(int i=0; i<2; i++) {
		for(int j=0; j<2; j++) {
			cout << data0[i][j] << " ";
		}
		cout << endl;
	}	
	
	e->process_inplace("math.sigma", Dict("value1", 1.0, "value2", 1.0));
	MArray2D data = e->get_2dview();
	for(int i=0; i<2; i++) {
		for(int j=0; j<2; j++) {
			cout << data[i][j] << " ";
		}
		cout << endl;
	}	
}
/*
void debug_peak_search()
{
	cout << "Start debug_peak_search function ..." << endl;
	
	EMData * e = new EMData();
	e->set_size(1024, 1024, 1);
	e->process_inplace("testimage.noise.uniform.rand");
	
	SparxUtil * su = new SparxUtil();
	vector<Peak> vp = su->peak_search(e, 1, 0.5);
	
	cout << "End debug_peak_search function ..." << endl;
}
*/

void emobject_cast_test(EMObject& obj)
{
	cout << "Enter emobject_cast_test() ..." << endl;
	
	//int dims = ((vector <float>)obj).size();
	//float dim = (float)obj;
	//double dim = (double)obj;
	//const char * dim = (const char *)obj;
	//float dim = (float)obj;
	//EMData * em = (EMData *)obj;
	//XYData * dim = (XYData *)obj;
	//vector<string> dim = obj;
	vector<float> dim = obj.operator vector<float>();
	
	cout << "Leave emobject_cast_test() ..." << endl;
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
		//debug_write_image();
		//debeg_dict();
		//debug_complex_image_arithmetic();
		//debug_big_file();
		//debug_insert_clip();
		//debug_do_fft();
		//debug_do_fft_inplace();
		//debug_calc_radial_dist();
		//debug_calc_hist();
		//debug_set_value_at();
		//debug_common_lines();
		//debug_sigma_processor();
		//debug_peak_search();
		//EMObject obj;
		//emobject_cast_test(obj);
	}
	catch (E2Exception & e) {
		cout << e.what();
	}
	catch (...) {
		cout << "Unknown exception ???" << endl;
	}
	
	
	return 0;
}
