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

#include <iostream>
#include "emdata.h"
#include <vector>
#include "randnum.h"

#include <boost/shared_array.hpp>
using boost::shared_array;
#include <boost/shared_ptr.hpp>
using boost::shared_ptr;

using namespace std;
using namespace EMAN;

void testfunc1()
{
	char t;

	const int SIZE = 20000;
	Randnum *r = Randnum::Instance();

	EMData * pEM[SIZE];
	for(int j=0; j<SIZE; ++j) {
		EMData *volume = new EMData(64,64); // initial volume
		volume->process_inplace("testimage.noise.uniform.rand");
		pEM[j] = volume;
	}

	cout << "Before set float array attribute...type a character to continue:" << endl;
	cin>>t;

	for(int i=0; i<SIZE; ++i) {
		vector<float> vf;
		for (int j=0; j<SIZE; ++j) {
			vf.push_back(r->get_frand());
		}

//		cout << "Before actually set float array attribute...xsize = " << (*itvEM)->get_xsize() << endl;
		pEM[i]->set_attr("farr", vf);
	}

	cout << "Before release memory...type a character to continue:" << endl;
	cin>>t;

	vector<float> vff =  pEM[1]->get_attr("farr");
	for(int i=0; i<100; ++i) {
		cout << vff[i] << endl;
	}

	for(int i=0; i<SIZE; ++i) {
		if (pEM[i]) delete pEM[i];
	}

	cout << "After release EMData memory...type a character to continue:" << endl;
	cin>>t;

//	vector<float> vff2 =  pEM[1]->get_attr("farr");
//	for(int i=0; i<100; ++i) {
//		cout << vff2[i] << endl;
//	}
}

void testfunc2()
{
	char t;

	const int SIZE = 20000;
	Randnum *r = Randnum::Instance();

	EMData * pEM[SIZE];
	for(int j=0; j<SIZE; ++j) {
		EMData *volume = new EMData(64,64); // initial volume
		volume->process_inplace("testimage.noise.uniform.rand");
		pEM[j] = volume;
	}

	cout << "Before set float array attribute...type a character to continue:" << endl;
	cin>>t;

	for(int i=0; i<SIZE; ++i) {
		vector<float> vf;
		for (int j=0; j<SIZE; ++j) {
			vf.push_back(r->get_frand());
		}

//		cout << "Before actually set float array attribute...xsize = " << (*itvEM)->get_xsize() << endl;
		pEM[i]->set_attr("farr", vf);
	}

	cout << "Before release memory...type a character to continue:" << endl;
	cin>>t;

	for(int i=0; i<SIZE; ++i) {
		if (pEM[i]) delete pEM[i];
	}

	cout << "After release EMData memory...type a character to continue:" << endl;
	cin>>t;
}

void testfunc3()
{
	char t;
	const int SIZE = 100000000;
	Randnum *r = Randnum::Instance();

	EMData *volume = new EMData(64,64); // initial volume
	volume->process_inplace("testimage.noise.uniform.rand");

	cout << "Before set float array attribute...type a character to continue:" << endl;
	cin>>t;

	vector<float> vf;
	for (int j=0; j<SIZE; ++j) {
		vf.push_back(r->get_frand());
	}

	volume->set_attr("farr", vf);


	cout << "Before release memory...type a character to continue:" << endl;
	cin>>t;

	delete volume;

	cout << "After release EMData memory...type a character to continue:" << endl;
	cin>>t;
}

void testfunc4()
{
	char t;
	const int SIZE = 100000000;
	Randnum *r = Randnum::Instance();

	EMData *volume = new EMData(64,64); // initial volume
	volume->process_inplace("testimage.noise.uniform.rand");

	cout << "Before set float array attribute...type a character to continue:" << endl;
	cin>>t;

	vector<float> vf;
	for (int j=0; j<SIZE; ++j) {
		vf.push_back(r->get_frand());
	}

	cout << "After create float vector...type a character to continue:" << endl;
	cin>>t;

	volume->set_attr("farr", vf);

	cout << "After set float vector as attribute...type a character to continue:" << endl;
	cin>>t;

	vector<float> vf2 = volume->get_attr("farr");
	for(int i=SIZE-1; i>SIZE-10; --i) {
		cout << vf2[i] << endl;
	}

	cout << "After retrieve float vector attribute, Before release memory...type a character to continue:" << endl;
	cin>>t;

	delete volume;

	cout << "After release EMData memory...type a character to continue:" << endl;
	cin>>t;
}

//test shared_array<float>
void testfunc5()
{
	char t;
	const int SIZE = 100000000;
	Randnum *r = Randnum::Instance();

	EMData *volume = new EMData(64,64); // initial volume
	volume->process_inplace("testimage.noise.uniform.rand");

	cout << "Before set float array attribute...type a character to continue:" << endl;
	cin>>t;

	shared_array<float> saf(new float[SIZE]);
	for (int j=0; j<SIZE; ++j) {
		saf[j] = r->get_frand();
	}

	for(int i=SIZE-1; i>SIZE-10; --i) {
		cout << saf[i] << '\t';
	}
	cout << endl;

//	shared_array<float> ssff = saf;

	cout << "After create float array...type a character to continue: use_count = " << saf.use_count() << ", raw pointer = " << (void*)saf.get() << endl;
	cin>>t;

	volume->set_attr("farr", saf);

	cout << "After set float vector as attribute...type a character to continue: use_count = " << saf.use_count() << ", raw pointer = " << (void*)saf.get() << endl;
	cin>>t;

	shared_array<float> saf2 = volume->get_attr("farr");

	cout << "after get attr, use_count = " << saf2.use_count() << endl;

	for(int i=SIZE-1; i>SIZE-10; --i) {
		cout << saf2[i] << '\t';
	}
	cout << endl;

	cout << "After retrieve float vector attribute, Before release memory...type a character to continue: use_count = " << saf2.use_count() << ", raw pointer = " << (void*)saf.get() << endl;
	cin>>t;

	delete volume;

	cout << "After release EMData memory...type a character to continue: use_count = " << saf.use_count() << ", raw pointer = " << (void*)saf.get() << endl;
	cin>>t;
}

//test shared_ptr< vector<float> >
void testfunc6()
{
	char t;
	const int SIZE = 100000000;
	Randnum *r = Randnum::Instance();

	EMData *volume = new EMData(64,64); // initial volume
	volume->process_inplace("testimage.noise.uniform.rand");

	cout << "Before set float array attribute...type a character to continue:" << endl;
	cin>>t;

	shared_ptr< vector<float> > spvf(new vector<float>());
//	shared_array<float> saf(new float[SIZE]);
	for (int j=0; j<SIZE; ++j) {
		spvf->push_back(r->get_frand());
	}

	for(int i=SIZE-1; i>SIZE-10; --i) {
		cout << spvf->operator[](i) << '\t';
	}
	cout << endl;

//	shared_array<float> ssff = saf;

	cout << "After create float array...type a character to continue: use_count = " << spvf.use_count() << ", raw pointer = " << (void*)spvf.get() << endl;
	cin>>t;

	volume->set_attr("farr", spvf);

	cout << "After set float vector as attribute...type a character to continue: use_count = " << spvf.use_count() << ", raw pointer = " << (void*)spvf.get() << endl;
	cin>>t;

	shared_ptr< vector<float> > spvf2 = volume->get_attr("farr");

	cout << "after get attr, use_count = " << spvf2.use_count() << ", raw pointer = " << (void*)spvf2.get() << endl;

	for(int i=SIZE-1; i>SIZE-10; --i) {
		cout << spvf2->operator[](i) << '\t';
	}
	cout << endl;

	cout << "After retrieve float vector attribute, Before release memory...type a character to continue: use_count = " << spvf.use_count() << ", raw pointer = " << (void*)spvf.get() << endl;
	cin>>t;

	delete volume;

	cout << "After release EMData memory...type a character to continue: use_count = " << spvf.use_count() << ", raw pointer = " << (void*)spvf.get() << endl;
	cin>>t;
}

int main()
{
	cout << "Hello, memory leak valgrind test!" << endl;
	char t;

//	cout << "testfunc1: access vector attribute before delete image..." << endl;
//	testfunc1();

//	cout << "testfunc2: NO access to vector attribute before delete image..." << endl;
//	testfunc2();

//	cout << "testfunc3: Single image test..." << endl;
//	testfunc3();
//	cout << "After testfunc3...type a character to continue:" << endl;
//	cin>>t;

//	cout << "testfunc4: Single image test, access attribute after set..." << endl;
//	testfunc4();
//	cout << "After testfunc4...type a character to continue:" << endl;
//	cin>>t;

//	cout << "testfunc5: Single image test, access attribute after set..." << endl;
//	testfunc5();
//	cout << "After testfunc5...type a character to continue:" << endl;
//	cin>>t;

	cout << "testfunc6: Single image test, access attribute after set..." << endl;
	testfunc6();
	cout << "After testfunc6...type a character to continue:" << endl;
	cin>>t;

	return 0;
}
