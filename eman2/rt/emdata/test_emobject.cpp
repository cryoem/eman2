/*
 * Author: David Woolford, 06/12/2007 (woolford@bcm.edu)
 * Copyright (c) 2000-2007 Baylor College of Medicine
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

 /** Code to test the functioning of EMObject and Dict
  *
  * For EMObject the things tested are
  * 	construction, copy construction and assigment behave as expected for each of the supported types
  * 	operator== behaves as expected for all types
  * 	all of the conversion operators are accurate and true
  * 	get_object_type_name, get_type and get_type_string are all internally consistent, for all types
  * 	the bevavior of is_null is as expected
  *
  *
  * For Dict things tested are
  * 	const and non const square brackets operators behave as expected
  * 	operator== behaves as expected (the correct behaviour of operator!= is therefore implicit)
  * 	the behavior of find is as expected
  * 	copy construction and assignment behave as expected
  * 	construction from map works as expected
  * 	the testing of the Dict iterator is implicit but not explicit
  * Things not tested are
  * 	keys, values, has_key, get, put, erase, set_default
  * 	Not that get and put work implicitly if the operator[] functions work.
  *
  * For EMObject there is also some curiosity tests to see what unsupported types are cast to, and some
  * commented out code that shows some un-compilable types
  */

#include "emobject.h"
#include "emdata.h"
#include "transform.h"
#include "xydata.h"
#include "ctf.h"

#include <iostream>
using std::cout;
using std::endl;

#include <vector>
using std::vector;

#include <string>
using std::string;

#include <map>
using std::map;

using namespace EMAN;

template<typename Type>
void test_emobject_specific_conversion( const Type& type )
{
	if ( static_cast<Type>(EMObject(type)) != type )
	{
		cout << "FAILED" << endl;
	}
	else
		cout << "passed" << endl;
}

template<typename Type>
void test_emobject_specific_conversion_vector( const vector<Type>& type )
{
	vector<Type> converted = EMObject(type);

	bool is_the_same = true;

	if ( converted.size() != type.size() )
	{
		is_the_same = false;
	}
	else
	{
		typename vector<Type>::const_iterator it1 = type.begin(), it2 = converted.begin();
		for( ;it1 != type.end(); ++it1 )
		{
			if ( *it1 != *it2 )
			{
				is_the_same = false;
				break;
			}
		}
	}

	if ( !is_the_same )
		cout << "FAILED" << endl;
	else
		cout << "passed" << endl;
}

// argh this is a horrible mess I wish all code would use strings exclusively, but oh well...
void test_emobject_specific_conversion( const char* type )
{
	if ( static_cast<string>(EMObject(type)) != string(type) )
		cout << "FAILED" << endl;
	else
		cout << "passed" << endl;
}


void test_emobject_assignment_and_equality( const EMObject& object1 )
{
	EMObject object2(object1);

	cout << "EMObject testing assignment and equality using " << object1.get_type_string() << " ..... ";
	// Test both the inequality and equality operators in one go
	// Though in reality it shouldn't matter, they are both tightly coupled.
	if ( object2 != object1 || !(object2 == object1) )
	{
		cout << "FAILED";
		if ( object2 != object1 )
			cout << "The inequality operator returned true" << endl;
		if ( !(object1 == object2) )
			cout << "The equality operator returned false" << endl;
	}
	else
	{
		cout << "passed";
	}
	cout << endl;
}

void test_emobject_isnull( const vector<EMObject>& objects )
{
	// I don't know why I put this test in it is hard to conceive a situation where it fails...
	// i.e. am I testing against impossible failure and being stupid?
	cout << "Testing the behavior of EMObject.is_null() for all types ..... ";

	bool this_works = true;

	for( vector<EMObject>::const_iterator it = objects.begin(); it != objects.end(); ++it )
	{
		if ( it->get_type() == EMObject::UNKNOWN && !it->is_null() )
		{
			this_works = false;
			break;
		}
		else
		{
			if ( it->get_type() != EMObject::UNKNOWN && it->is_null() )
			{
			 	this_works = false;
				break;
			}
		}
	}

	if ( !this_works )
		cout << "FAILED" << endl;
	else
		cout << "passed" << endl;
}

void test_emobject_object_types_bevaviour( const vector<EMObject>& objects )
{
	// Again, I don't know why I put this test in it is hard to conceive of a situation where it fails...
	// i.e. am I testing against impossible failure and being stupid?

	cout << "Testing the behavior of EMObject.get_type() and EMObject.get_type_string(): " << endl;

	for( vector<EMObject>::const_iterator it = objects.begin(); it != objects.end(); ++it )
	{
		cout << "Testing internal consisteny for type " << it->get_type_string() << " ..... ";

		if  ( EMObject::get_object_type_name(it->get_type()) == it->get_type_string() )
			cout << "passed";
		else
			cout << "FAILED";

		cout << endl;
	}

}

void test_emobject_unsupported_types()
{
	// test cases where the type is not supported by the EMObject
	char c1 = 'a';
	cout << "Testing for character behaviour even though it is not supported, this was turned into type ..... " << EMObject(c1).get_type_string() << endl;

	unsigned char c2 = 125;
	cout << "Testing for unsigned character behaviour even though it is not supported, this was turned into type ..... ";
	cout << EMObject(c2).get_type_string() << endl;

	// TYPES such as these do not compile
//	unsigned int ui = 0;
//	cout << "Testing for unsigned int behaviour even though it is not supported, this was turned into type ..... ";
//	cout << EMObject(ui).get_type_string() << endl;

//	long l1 = 0;
//	cout << "Testing for long (int) behaviour even though it is not supported, this was turned into type ..... ";
//	cout << EMObject(l1).get_type_string() << endl;
//
//	long double l2 = 0;
//	cout << "Testing for long double behaviour even though it is not supported, this was turned into type ..... ";
//	cout << EMObject(l2).get_type_string() << endl;

	// TYPES such as these (pointers) are recorded as bools in EMObjects on my OS (Fedora Core 6)
	int i1 = 0;
	cout << "Testing for int * behaviour even though it is not supported, this was turned into type ..... ";
	cout << EMObject(&i1).get_type_string() << endl;

	float f1 = 0;
	cout << "Testing for float * behaviour even though it is not supported, this was turned into type ..... ";
	cout << EMObject(&f1).get_type_string() << endl;

	double d1 = 0;
	cout << "Testing for double * behaviour even though it is not supported, this was turned into type ..... ";
	cout << EMObject(&d1).get_type_string() << endl;

	cout << "Currently unsigned ints and longs are not supported... a more rigorous list of unsupported types is probably needed " << endl;
}

void test_emobject_conversion()
{
	cout << "Testing boolean conversion operator ..... ";
	test_emobject_specific_conversion( (bool) false );

	cout << "Testing integer conversion operator ..... ";
	test_emobject_specific_conversion( (int) 1235412 );

	cout << "Testing float conversion operator ..... ";
	test_emobject_specific_conversion( (float) 12.312 );

	cout << "Testing double conversion operator ..... ";
	test_emobject_specific_conversion( (double) 0.00201 );

	cout << "Testing const char* conversion operator..... ";
	test_emobject_specific_conversion( "hello world" );

	cout << "Testing float pointer conversion operator..... ";
	float* pFloat = new float;
	test_emobject_specific_conversion( pFloat );
	delete pFloat;

	cout << "Testing xydata pointer conversion operator..... ";
	XYData* xydata = new XYData;
	test_emobject_specific_conversion( xydata );
	delete xydata;

	cout << "Testing Transform pointer conversion operator..... ";
	Transform* pTransform = new Transform;
	test_emobject_specific_conversion( pTransform );
	delete pTransform;

	cout << "Testing Ctf pointer conversion operator..... ";
	Ctf* pCtf = new EMAN2Ctf();
	test_emobject_specific_conversion( pCtf );
	delete pCtf;

	cout << "Testing EMData pointer conversion operator..... ";
	EMData* pEMData = new EMData;
	test_emobject_specific_conversion( pEMData );
	delete pEMData;

	cout << "Testing vector<int> conversion operator..... ";
	test_emobject_specific_conversion_vector( vector<int>(1000,123) );

	cout << "Testing vector<float> conversion operator..... ";
	test_emobject_specific_conversion_vector( vector<float>(20,23.23421f) );

	cout << "Testing vector<string> conversion operator..... ";
	test_emobject_specific_conversion_vector( vector<string>(100,"empty") );
}

// Returns a vector containing an instance of each type of EMObject - is not dynamic, so if changes are made to
// EMObject this function would need to be updated
// no memory management is happening here, there are memory leaks in this program, but oh well it's only for testing
// purposes. This is not a trivial problem to fix, because in the EMAN2 body of code much copying of dictionaries occurs,
// and it is often imperative that the pointers in an EMObject are not destroyed - because ownership is never assumed.
// A solution might be the use of policies, or smart pointers etc...
vector<EMObject> get_test_emobjects()
{
	vector<EMObject> objects;

	objects.push_back(EMObject());

	bool b = false;
	objects.push_back(EMObject(b));

	int i = 0;
	objects.push_back(EMObject(i));

	float f = 12345.2124f;
	objects.push_back(EMObject(f));

	double d = .00002;
	objects.push_back(EMObject(d));

	const char* c = "hello from test code";
	objects.push_back(EMObject(c));

	string s("a string in test code");
	objects.push_back(EMObject(s));

	XYData* xydata = new XYData;
	objects.push_back(EMObject(xydata));

	EMData *a = new EMData();
	objects.push_back(EMObject(a));

	float *fp = new float;
	objects.push_back(EMObject(fp));

	Transform* pTransform = new Transform;
	objects.push_back(EMObject(pTransform));
	delete pTransform;

	Ctf* ctf = new EMAN1Ctf();
	objects.push_back(EMObject(ctf));
	delete ctf;

	vector<int> iv(100,2);
	objects.push_back(EMObject(iv));

	vector<float> fv(1234,2.00);
	objects.push_back(EMObject(fv));

	vector<string> sv(100000,"empty");
	objects.push_back(EMObject(sv));

	return objects;
}

// Used primarily for testing DICT
map<string, EMObject> get_test_dict_keys()
{
	vector<EMObject> objects = get_test_emobjects();
	map<string, EMObject> return_objects;

	for( vector<EMObject>::const_iterator it = objects.begin(); it != objects.end(); ++it )
	{
		return_objects[it->get_type_string()] = *it;
	}

	return return_objects;
}

void test_emobject()
{

	vector<EMObject> objects = get_test_emobjects();

	cout << "-----------------------------------------------------------------------------" << endl;

	for( vector<EMObject>::const_iterator it = objects.begin(); it != objects.end(); ++it )
	{
		test_emobject_assignment_and_equality(*it);
	}

	cout << "-----------------------------------------------------------------------------" << endl;
	// Now check to see if copy works irrespective of the current tyoe
	bool can_copy_irrespective_of_type = true;
	for( vector<EMObject>::const_iterator it = objects.begin(); it != objects.end(); ++it )
	{
		for( vector<EMObject>::const_iterator it2 = objects.begin(); it2 != objects.end(); ++it2 )
		{
			EMObject a = *it, b=*it2;
			b = a;
			if ( b != a)
			{
				cout <<	"Copying type " << it2->get_type_string() << " to an object of type " << it->get_type_string();
				cout << " FAILED " << endl;
				can_copy_irrespective_of_type = false;
			}
		}
	}
	cout << "Testing to make sure objects that already have types are accuracately copied ..... ";
	if ( can_copy_irrespective_of_type )
		cout << "passed" << endl;
	else
		cout << "FAILED" << endl;

	cout << "-----------------------------------------------------------------------------" << endl;
	test_emobject_unsupported_types();
	cout << "-----------------------------------------------------------------------------" << endl;
	test_emobject_conversion();
	cout << "-----------------------------------------------------------------------------" << endl;
	test_emobject_isnull(objects);
	cout << "-----------------------------------------------------------------------------" << endl;
	test_emobject_object_types_bevaviour(objects);
}

void print_success_message( const bool success, const string& message = "Overal the test ..... " )
{
	cout << endl << message;
	if ( success )
		cout << "passed";
	else
		cout << "FAILED";
	cout << endl;
}

bool test_dict_map_constructor()
{
	cout << "Testing construction from map" << endl << endl;;
	map<string, EMObject> objects = get_test_dict_keys();
	Dict dict(objects);

	bool success = true;

	map<string, EMObject>::const_iterator mapIt = objects.begin();
	for( ; mapIt != objects.end(); ++mapIt )
	{
		bool found = false;
		cout << "Trying to find map key " << mapIt->first << " in the Dict constructed from it ..... ";
		// Could use Dict::find but that is a different function to test
		for( Dict::const_iterator it = dict.begin(); it != dict.end(); ++it )
		{
			if ( it->first == mapIt-> first && it->second == mapIt->second )
			{
				found = true;
				break;
			}
		}

		if ( found )
			cout << "passed";
		else
		{
			cout << "FAILED";
			success = false;
		}

		cout << endl;
	}

	print_success_message(success);

	return success;
}

bool test_copy_construction_and_assigment()
{
	map<string, EMObject> objects = get_test_dict_keys();
	Dict dict(objects);
	Dict dict2 = dict;

	bool success = true;

	cout << "Ierating through two maps, one containing an instance of every EMObject, the other copy constructed from the first," << endl;
	cout << "to test pair wise equality:" << endl << endl;

	Dict::const_iterator it = dict.begin(), it2 = dict2.begin();
	for(; it != dict.end(); ++it, ++it2 )
	{
		cout << "Iter1 type is " << it->first << " iter2 type is " << it2->first << ", testing equality ...... ";
		if ( it->first != it2->first || it->second != it2->second )
		{
			success = false;
			cout << "FAILED";
			if ( it->first != it2->first )
				cout << " ..... pair.first members were unequal (strings)";
			if ( it->second != it2->second )
				cout << " ..... pair.second members were unequal (EMObjects)";
		}
		else
			cout << "passed";

		cout << endl;
	}

	print_success_message( success );

	cout << "-----------------------------------------------------------------------------" << endl;
	cout << "Testing assigment and equality operator==..... ";
	success = true;
	if ( dict != dict2 )
	{
		success = false;
		cout << "FAILED";
	}
	else
		cout << "passed";
	cout << endl;

	print_success_message( success );
	return success;
}

bool test_dict_non_const_square_brackets_operator()
{
	cout << "Testing non const operator[] to ensure it is inserting values as expected" << endl << endl;
	map<string, EMObject> objects = get_test_dict_keys();

	Dict dict;

	// Here is where the insertion occurs
	map<string, EMObject>::const_iterator mapIt = objects.begin();
	for( ; mapIt != objects.end(); ++mapIt )
	{
		dict[mapIt->first] = mapIt->second;
	}

	bool success = true;


	for( mapIt = objects.begin(); mapIt != objects.end(); ++mapIt )
	{
		bool found = false;
		cout << "Testing the insertion of map key/value pair with key " << mapIt->first << " (inserted into Dict with operator[]) ..... ";
		for( Dict::const_iterator it = dict.begin(); it != dict.end(); ++it )
		{
			if ( it->first == mapIt-> first && it->second == mapIt->second )
			{
				found = true;
				break;
			}
		}

		if ( found )
			cout << "passed";
		else
		{
			cout << "FAILED";
			success = false;
		}

		cout << endl;
	}

	print_success_message(success);
	return success;
}

bool test_dict_find()
{
	cout << "Testing find functionality by iterating through a dict containing all of the EMObjects and ensuring they can be found within itself" << endl << endl;
	Dict dict(get_test_dict_keys());

	bool success = true;
	for(Dict::const_iterator it = dict.begin();it != dict.end();++it )
	{
		cout << "Searching self for object with key (" << it->first << ") ..... ";
		if ( dict.find(it->first) == dict.end() )
		{
			success = false;
			cout << "FAILED";
		}
		else
			if ( dict.find(it->first)->second != it->second )
			{
				success = false;
				cout << " was found but the pair.second was not the same ..... FAILED";
			}
			else
				cout << "passed";

		cout << endl;
	}

	print_success_message(success);

	return success;
}

bool test_square_brackets_operator()
{
	cout << "Testing operator[] functionality by iterating through the Dict and asking if it->second != dict[it->first] " << endl << endl;

	Dict dict(get_test_dict_keys());

	bool success = true;

	for(Dict::const_iterator it = dict.begin();it != dict.end();++it )
	{
		cout << "For type " << it->first << " operator[] ..... ";
		if ( it->second != dict[it->first] )
		{
			success = false;
			cout << "FAILED" << endl;
		}
		else
			cout << "passed" << endl;
	}

	print_success_message(success);

	return success;
}

void test_dict()
{
	cout << "-----------------------------------------------------------------------------" << endl;
	if ( test_dict_map_constructor() )
	{
		cout << "-----------------------------------------------------------------------------" << endl;
		test_copy_construction_and_assigment();
		cout << "-----------------------------------------------------------------------------" << endl;
		test_dict_find();
		cout << "-----------------------------------------------------------------------------" << endl;
		test_square_brackets_operator();
	}
	else
	{
		cout << "FOOBAR - WARNING: because construction of a Dict from a map failed, many other tests will fail, aborted" << endl;
		return;
	}

	cout << "-----------------------------------------------------------------------------" << endl;
	test_dict_non_const_square_brackets_operator();

}

void test_attr()
{
	cout << "-----------------------------------------------------------------------------" << endl;
	bool success = true;
	EMData * img = new EMData(32,32);
	float arr[11] = {0.1f, 1.1f, 2.2f, 3.3f, 4.4f, 5.5f, 6.6f, 7.7f, 8.8f, 9.9f, 10.1f};
	vector<float> v1(arr, arr+11);

	Ctf* ctf1 = new EMAN1Ctf();
	ctf1->from_vector(v1);
	img->set_attr("ctf1", ctf1);
	delete ctf1;
	img->write_image("test_image1.hdf");

	EMData* img2 = new EMData("test_image1.hdf");
	Ctf* ctf11 = img2->get_attr("ctf1");
	string s1 = ctf11->to_string();
	cout << "ctf1 string: " << s1 << endl;
	if(success && s1 != "O0.1 1.1 2.2 3.3 4.4 5.5 6.6 7.7 8.8 9.9 10.1") {
		success = false;
	}
	delete ctf11;
	delete img2;

//	Ctf* ctf2  = new EMAN2Ctf();
//	ctf2->from_vector(v1);
//	img->set_attr("ctf2", ctf2);
//	delete ctf2;
//	img->write_image("test_image2.hdf");
//
//	EMData* img3 = new EMData("test_image2.hdf");
//	Ctf* ctf22 = img3->get_attr("ctf2");
//	string s2 = ctf22->to_string();
//	cout << "ctf2 string: " << s2 << endl;
//	if(success && s2 != "E0.1 1.1 2.2 3.3 4.4 5.5") {
//		success = false;
//	}
//	delete ctf22;

	delete img;
	if(success) {
		cout << "Testing set Ctf object as image attribute ............................. passed" << endl;
	}
	else {
		cout << "Testing set Ctf object as image attribute ............................. failed" << endl;
	}


}

int main(int, char**)
{
	cout << "-----------------------------------------------------------------------------" << endl;
	cout << "TESTING EMOBJECT" << endl << endl;
	test_emobject();
	cout << endl << endl;
	cout << "-----------------------------------------------------------------------------" << endl;
	cout << "TESTING DICT" << endl << endl;
	test_dict();
	cout << endl;
	cout << "-----------------------------------------------------------------------------" << endl;
	cout << "TESTING Attribute" << endl << endl;
	test_attr();
	cout << endl;

	return 0;
}

