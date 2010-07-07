/**
 * $Id$
 */

/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
 * Probable contributor: Liwei Peng (what dates?)
 * Contributing author: David Woolford 06/11/2007
 *
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

#include "emobject.h"
#include <cmath>
#ifdef WIN32
#define M_PI 3.14159265358979323846f
#endif

#include <algorithm>
// using copy

#include "util.h"

using namespace EMAN;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

const float EMConsts::I2G = (float) (4.0 / (M_PI*M_PI));
const float EMConsts::I3G = (float) (6.4 / (M_PI*M_PI));
const float EMConsts::I4G = (float) (8.8 / (M_PI*M_PI));
const float EMConsts::I5G = (float) (10.4 / (M_PI*M_PI));
// This stolen from wikipedia.org
const double EMConsts::pi = 3.141592653589793238462643383279502884197169399;
const double EMConsts::deg2rad = pi/180.0;
const double EMConsts::rad2deg = 180.0/pi;


#include <sstream>
using std::stringstream;

#include "transform.h"
#include "ctf.h"

// Static init
map< EMObject::ObjectType, string>  EMObject::type_registry = init();;

//-------------------------------EMObjectTypes-----------------------------------------
map< EMObject::ObjectType, string> EMObject::init()
{
	map< EMObject::ObjectType, string> mymap;
	static bool first_construction = true;
	if ( first_construction )
	{
		// Initialize the the type registry once and for all
		mymap[BOOL] = "BOOL";
		mymap[INT] = "INT";
		mymap[UNSIGNEDINT] = "UNSIGNEDINT";
		mymap[FLOAT] = "FLOAT";
		mymap[DOUBLE] = "DOUBLE";
		mymap[STRING] = "STRING";
		mymap[EMDATA] = "EMDATA";
		mymap[XYDATA] = "XYDATA";
		mymap[INTARRAY] = "INTARRAY";
		mymap[FLOATARRAY] = "FLOATARRAY";
		mymap[STRINGARRAY] = "STRINGARRAY";
		mymap[TRANSFORM] = "TRANFORM";
		mymap[CTF] = "CTF";
		mymap[FLOAT_POINTER] = "FLOAT_POINTER";
		mymap[INT_POINTER] = "INT_POINTER";
		mymap[UNKNOWN] = "UNKNOWN";
		mymap[VOID_POINTER] = "VOID_POINTER";
		first_construction = false;
	}
	
	return mymap;
}

//-------------------------------EMObject--------------------------------------------

void EMObject::printInfo() const
{
	cout << "The address of my type is " << &type << endl;
	cout << " Now printing the enumerated values in type_registry " << endl;
	for( map< ObjectType, string>::const_iterator it = type_registry.begin(); it != type_registry.end(); ++it )
	{
		cout << it->first << " " << it->second << endl;
	}
	cout << "My type is " << to_str(type) << " and its enumerated value is " << type << endl;
	cout << "The address of the static type registry is " << &type_registry <<", it should be same for all EMObjects" << endl;
}

EMObject::~EMObject() {
#ifdef MEMDEBUG
	allobjlist.erase(this);
	printf("  -(%6d) %p\n",allobjlist.count(),this);
#endif
}

EMObject::EMObject() :
	d(0), type(UNKNOWN)
{
#ifdef MEMDEBUG
	allobjlist.insert(this);
	printf("  +(%6d) %p\n",allobjlist.count(),this);
#endif
}

EMObject::EMObject(bool boolean) :
	b(boolean), type(BOOL)
{
#ifdef MEMDEBUG
	allobjlist.insert(this);
	printf("  +(%6d) %p\n",allobjlist.count(),this);
#endif
}

EMObject::EMObject(int num) :
	n(num), type(INT)
{
#ifdef MEMDEBUG
	allobjlist.insert(this);
	printf("  +(%6d) %p\n",allobjlist.count(),this);
#endif
}

EMObject::EMObject(unsigned int num) :
	ui(num), type(UNSIGNEDINT)
{
#ifdef MEMDEBUG
	allobjlist.insert(this);
	printf("  +(%6d) %p\n",allobjlist.count(),this);
#endif
}

EMObject::EMObject(float ff) :
	type(FLOAT)
{

#ifdef MEMDEBUG
	allobjlist.insert(this);
	printf("  +(%6d) %p\n",allobjlist.count(),this);
#endif
	if(Util::goodf(&ff)) {
		f = ff;
	}
	else{
		f = 0.0f;
	}
}

EMObject::EMObject(double dd) :
	type(DOUBLE)
{

#ifdef MEMDEBUG
	allobjlist.insert(this);
	printf("  +(%6d) %p\n",allobjlist.count(),this);
#endif
	if(Util::goodf(&dd)) {
		d = dd;
	}
	else{
		d = 0.0;
	}
}

EMObject::EMObject(const char *s) :
	str(s), type(STRING)
{
#ifdef MEMDEBUG
	allobjlist.insert(this);
	printf("  +(%6d) %p\n",allobjlist.count(),this);
#endif
}

EMObject::EMObject(const string & s) :
	str(s), type(STRING)
{
#ifdef MEMDEBUG
	allobjlist.insert(this);
	printf("  +(%6d) %p\n",allobjlist.count(),this);
#endif
}

EMObject::EMObject(float *f) :
	fp(f), type(FLOAT_POINTER)
{
#ifdef MEMDEBUG
	allobjlist.insert(this);
	printf("  +(%6d) %p\n",allobjlist.count(),this);
#endif
}

EMObject::EMObject(int *i) :
	ip(i), type(INT_POINTER)
{
#ifdef MEMDEBUG
	allobjlist.insert(this);
	printf("  +(%6d) %p\n",allobjlist.count(),this);
#endif
}

EMObject::EMObject(void *v) :
	vp(v), type(VOID_POINTER)
{
#ifdef MEMDEBUG
	allobjlist.insert(this);
	printf("  +(%6d) %p\n",allobjlist.count(),this);
#endif
}

EMObject::EMObject(EMData * em)	:
	emdata(em), type(EMDATA)
{
#ifdef MEMDEBUG
	allobjlist.insert(this);
	printf("  +(%6d) %p\n",allobjlist.count(),this);
#endif
}

EMObject::EMObject(XYData * xy) :
	xydata(xy), type(XYDATA)
{
#ifdef MEMDEBUG
	allobjlist.insert(this);
	printf("  +(%6d) %p\n",allobjlist.count(),this);
#endif
}

EMObject::EMObject(Transform* t) :
	farray(t->get_matrix()), type(TRANSFORM)
{
#ifdef MEMDEBUG
	allobjlist.insert(this);
	printf("  +(%6d) %p\n",allobjlist.count(),this);
#endif
}

EMObject::EMObject(Ctf * ctf) :
	str(ctf->to_string()), type(CTF)
{
#ifdef MEMDEBUG
	allobjlist.insert(this);
	printf("  +(%6d) %p\n",allobjlist.count(),this);
#endif
}

EMObject::EMObject(const vector< int >& v ) :
	iarray(v), type(INTARRAY)
{
#ifdef MEMDEBUG
	allobjlist.insert(this);
	printf("  +(%6d) %p\n",allobjlist.count(),this);
#endif
}

EMObject::EMObject(const vector < float >&v) :
	farray(v), type(FLOATARRAY)
{
#ifdef MEMDEBUG
	allobjlist.insert(this);
	printf("  +(%6d) %p\n",allobjlist.count(),this);
#endif
}

EMObject:: EMObject(const vector <string>& sarray) :
	strarray(sarray), type(STRINGARRAY)
{
#ifdef MEMDEBUG
	allobjlist.insert(this);
	printf("  +(%6d) %p\n",allobjlist.count(),this);
#endif
}

EMObject::operator bool () const
{
	if (type == BOOL) {
		return b;
	}
	else if (type == INT) {
		return n != 0;
	}
	else if (type == UNSIGNEDINT) {
		return ui != 0;
	}
	else if (type == FLOAT) {
		return f != 0;
	}
	else if (type == DOUBLE) {
		return d != 0;
	}
	else if (type == EMDATA) {
		return emdata != 0;
	}
	else if (type == XYDATA) {
		return xydata != 0;
	}
	else if (type == FLOAT_POINTER) {
		return fp != 0;
	}
	else if (type == INT_POINTER) {
		return ip != 0;
	}
	else if (type == VOID_POINTER) {
		return vp != 0;
	}
//	else if (type == TRANSFORM) {
//		return transform != 0;
//	}
	// It seemed unconventional to return a boolean for the stl objects
	else {
		if (type != UNKNOWN) {
			throw TypeException("Cannot convert to bool this data type ",
								get_object_type_name(type));
		}
	}
	return 0;
}

EMObject::operator int () const
{
	if (type == INT) {
		return n;
	}
	else if (type == UNSIGNEDINT) {
		return (int) ui;
	}
	else if (type == FLOAT) {
		return (int) f;
	}
	else if (type == DOUBLE) {
		return (int) d;
	}
	else if (type == BOOL) {
		return b?1:0;
	}
	else {
		if (type != UNKNOWN) {
			throw TypeException("Cannot convert to int this data type ",
								get_object_type_name(type));
		}
	}
	return 0;
}

EMObject::operator unsigned int () const
{
	if (type == UNSIGNEDINT) {
		return (unsigned int) ui;
	}
	else {
		if (type != UNKNOWN) {
			throw TypeException("Cannot convert to int this data type ",
								get_object_type_name(type));
		}
	}
	return 0;
}

EMObject::operator float () const
{
	if (type == BOOL) {
		return b?1.0f:0.0f;
	}
	else if (type == FLOAT) {
		return f;
	}
	else if (type == INT) {
		return (float) n;
	}
	else if (type == UNSIGNEDINT) {
		return (float) ui;
	}
	else if (type == DOUBLE) {
		return (float) d;
	}
	else {
		if (type != UNKNOWN) {
			throw TypeException("Cannot convert to float from this data type",
								get_object_type_name(type));
		}
	}

	return 0;
}

EMObject::operator double () const
{
	if (type == BOOL) {
	return b?1.0:0.0;
	}
	else if (type == DOUBLE) {
		return d;
	}
	else if (type == INT) {
		return (double) n;
	}
	else if (type == UNSIGNEDINT) {
		return (double) ui;
	}
	else if (type == FLOAT) {
		return (double) f;
	}
	else {
		if (type != UNKNOWN) {
			throw TypeException("Cannot convert to double from this data type",
								get_object_type_name(type));
		}
	}
	return 0;
}

EMObject::operator int * () const
{
	if (type != INT_POINTER)
	{
		if (type != UNKNOWN)
			throw TypeException("Cannot convert to float pointer from this data type",
								get_object_type_name(type));

		return 0;
	}

	return ip;
}


EMObject::operator float * () const
{
	if (type != FLOAT_POINTER)
	{
		if (type != UNKNOWN)
			throw TypeException("Cannot convert to float pointer from this data type",
								get_object_type_name(type));

		return 0;
	}

	return fp;
}

// MAYBE REMOVE? FIXME
EMObject::operator void * () const
{
	if (type == VOID_POINTER) return vp;
	else if (type == FLOAT_POINTER) return (void *)fp;
	else if (type == INT_POINTER) return (void *)ip;
	else if (type == EMDATA) return (void *) emdata;
	else if (type == XYDATA) return (void *) xydata;
// 	else if (type == TRANSFORM) return (void *) transform;
	else throw TypeException("Cannot convert to void pointer from this data type", get_object_type_name(type));
}

EMObject::operator const char * () const
{
	if (type != STRING && type != CTF) {
		stringstream ss;
		string return_string;
		if ( type == INT )
		{
			ss << n;
			ss >> return_string;
			return return_string.c_str();
		}
		if ( type == UNSIGNEDINT )
		{
			ss << ui;
			ss >> return_string;
			return return_string.c_str();
		}
		else
		if ( type == FLOAT )
		{
			ss << f;
			ss >> return_string;
			return return_string.c_str();
		}
		else
		if ( type == DOUBLE )
		{
			ss << d;
			ss >> return_string;
			return return_string.c_str();
		}
		else if (type != UNKNOWN) {
			throw TypeException("Cannot convert to string from this data type",
								get_object_type_name(type));
		}

		return "";
	}
	return str.c_str();
}

EMObject::operator EMData * () const
{
	if (type != EMDATA) {
		if (type != UNKNOWN) {
			throw TypeException("Cannot convert to EMData* from this data type",
				   get_object_type_name(type));
		}
		return 0;
	}
	return emdata;
}

EMObject::operator XYData * () const
{
	if (type != XYDATA) {
		if (type != UNKNOWN) {
			throw TypeException("Cannot convert to XYData* from this data type",
				   get_object_type_name(type));
		}
		return 0;
	}
	return xydata;
}

EMObject::operator Transform* () const
{
	if(type != TRANSFORM) {
		if(type != UNKNOWN) {
			throw TypeException("Cannot convert to TRANSFORM* from this data type",
								get_object_type_name(type));
		}
	}
	Transform * transform = new Transform();
	transform->set_matrix(farray);
	return transform;
}

EMObject::operator Ctf* () const
{
/*	if(type != CTF) {
		if(type != CTF) {
			throw TypeException("Cannot convert to TRANSFORM* from this data type",
								get_object_type_name(type));
		}
	}*/
	Ctf * ctf = 0;
	if(str[0] == 'O') {
		ctf = new EMAN1Ctf();
		ctf->from_string(str);
	}
	else if(str[0] == 'E') {
		ctf = new EMAN2Ctf();
		ctf->from_string(str);
	}
	return ctf;
}

EMObject::operator vector<int>() const
{
    if( type != INTARRAY )
    {
        if( type != UNKNOWN ) {
	    	throw TypeException("Cannot convert to vector<int> from this data type", get_object_type_name(type) );
		}
		return vector<int>();
    }
    return iarray;
}

EMObject::operator vector < float > () const
{
	if (type != FLOATARRAY) {
		if (type != UNKNOWN) {
			throw TypeException("Cannot convert to vector<float> from this data type",
								get_object_type_name(type));
		}
		return vector < float >();
	}
	return farray;
}

EMObject::operator vector<string> () const
{
	if (type != STRINGARRAY) {
		if (type != UNKNOWN) {
			throw TypeException("Cannot convert to vector<string> from this data type",
								get_object_type_name(type));
		}
		return vector<string>();
	}
	return strarray;
}

bool EMObject::is_null() const
{
	return (type == UNKNOWN);
}

string EMObject::to_str() const
{
	return to_str(type);
}

string EMObject::to_str(ObjectType argtype) const
{
	if (argtype == STRING) {
		return str;
	}
	else {
		char tmp_str[32];
		if (argtype == BOOL) {
			if (b)
				sprintf(tmp_str, "true");
			else
				sprintf(tmp_str, "false");
		}
		else if (argtype == INT) {
			sprintf(tmp_str, "%d", n);
		}
		else if (argtype == UNSIGNEDINT) {
			sprintf(tmp_str, "%d", ui);
		}
		else if (argtype == FLOAT) {
			sprintf(tmp_str, "%f", f);
		}
		else if (argtype == DOUBLE) {
			sprintf(tmp_str, "%f", d);
		}
		else if (argtype == EMDATA) {
			sprintf(tmp_str, "EMDATA");
		}
		else if (argtype == FLOAT_POINTER) {
			sprintf(tmp_str, "FLOAT_POINTER");
		}
		else if (argtype == INT) {
			sprintf(tmp_str, "INT_POINTER");
		}
		else if (argtype == VOID_POINTER) {
			sprintf(tmp_str, "VOID_POINTER");
		}
		else if (argtype == XYDATA) {
			sprintf(tmp_str, "XYDATA");
		}
		else if (argtype == INTARRAY) {
			sprintf(tmp_str, "INTARRAY");
		}
		else if (argtype == FLOATARRAY) {
			sprintf(tmp_str, "FLOATARRAY");
		}
		else if (argtype == STRINGARRAY) {
			sprintf(tmp_str, "STRINGARRAY");
		}
		else if (argtype == TRANSFORM) {
			sprintf(tmp_str, "TRANSFORMD");
		}
		else if (argtype == CTF) {
			sprintf(tmp_str, "CTF");
		}
		else if (argtype == UNKNOWN) {
			sprintf(tmp_str, "UNKNOWN");
		}
		else {
			LOGERR("No such EMObject defined");
			throw NotExistingObjectException("EMObject", "unknown type");
		}
		return string(tmp_str);
	}
}

string EMObject::get_object_type_name(ObjectType t)
{
#ifdef _WIN32
	if (t == BOOL) {
		return "BOOL";
	}else
	if ( t == INT){
		return "INT";
	}else
	if ( t == UNSIGNEDINT){
		return "UNSIGNEDINT";
	} else
	if ( t == FLOAT){
		return "FLOAT";
	} else
	if ( t == DOUBLE){
		return "DOUBLE";
	}else
	if ( t == STRING){
		return "STRING";
	}else
	if ( t == EMDATA){
		return "EMDATA";
	}
	else
	if ( t == XYDATA){
		return "XYDATA";
	}else
	if ( t == INTARRAY){
		return "INTARRAY";
	}else
	if ( t == FLOATARRAY){
		return "FLOATARRAY";
	} else
	if ( t == STRINGARRAY){
		return "STRINGARRAY";
	}else
	if ( t == TRANSFORM){
		return "TRANSFORM";
	}else
	if ( t == CTF){
		return "CTF";
	}else
	if ( t == FLOAT_POINTER){
		return "FLOAT_POINTER";
	}else
	if ( t == INT_POINTER){
		return "INT_POINTER";
	}else
	if ( t == UNKNOWN){
		return "UNKNOWN";
	} else
	if ( t == VOID_POINTER){
		return "VOID_POINTER";
	}
	else {
		LOGERR("No such EMObject defined");
		throw NotExistingObjectException("EMObject", "unknown type");
	}

#else

	if  ( type_registry.find(t) != type_registry.end() )
		return type_registry[t];
	else
		LOGERR("No such EMObject defined");
		throw NotExistingObjectException("EMObject", "unknown type");
#endif	//_WIN32
}

bool EMAN::operator==(const EMObject &e1, const EMObject & e2)
{

	if (e1.type != e2.type) {
		return false;
	}

	switch (e1.type) {
	case  EMObject::BOOL:
		return (e1.b == e2.b);
	break;
	case  EMObject::INT:
		return (e1.n == e2.n);
	break;
	case  EMObject::UNSIGNEDINT:
		return (e1.ui == e2.ui);
	break;
	case  EMObject::FLOAT:
		return (e1.f == e2.f);
	break;
	case  EMObject::DOUBLE:
		return (e1.d == e2.d);
	break;
	case EMObject::CTF:
	case  EMObject::STRING:
		return (e1.str == e2.str);
	break;
	case  EMObject::FLOAT_POINTER:
		return (e1.fp == e2.fp);
	break;
	case  EMObject::INT_POINTER:
		return (e1.ip == e2.ip);
	break;
	case  EMObject::VOID_POINTER:
		return (e1.vp == e2.vp);
	break;
	case  EMObject::EMDATA:
		return (e1.emdata == e2.emdata);
	break;
	case  EMObject::XYDATA:
		return (e1.xydata == e2.xydata);
	break;
	case  EMObject::TRANSFORM:
	case  EMObject::FLOATARRAY:
		if (e1.farray.size() == e2.farray.size()) {
			for (size_t i = 0; i < e1.farray.size(); i++) {
				if (e1.farray[i] != e2.farray[i]) {
					return false;
				}
			}
			return true;
		}
		else {
			return false;
		}
	break;
	case  EMObject::INTARRAY:
		if (e1.iarray.size() == e2.iarray.size()) {
			for (size_t i = 0; i < e1.iarray.size(); i++) {
				if (e1.iarray[i] != e2.iarray[i]) {
					return false;
				}
			}
			return true;
		}
	break;
	case  EMObject::STRINGARRAY:
		if (e1.strarray.size() == e2.strarray.size()) {
			for (size_t i = 0; i < e1.strarray.size(); i++) {
				if (e1.strarray[i] != e2.strarray[i]) {
					return false;
				}
			}
			return true;
		}
		else {
			return false;
		}
	break;
	case  EMObject::UNKNOWN:
		// UNKNOWN really means "no type" and if two objects both have
		// type UNKNOWN they really are the same
		return (e1.type == e2.type);
	break;
	default:
		return false;
	break;
	}
	return false;
}

bool EMAN::operator!=(const EMObject &e1, const EMObject & e2)
{
	return !(e1 == e2);
}

// Copy constructor
EMObject::EMObject(const EMObject& that)
{
	// init isn't necessary because that must have already called it!
	*this = that;
}


// Assignment operator -  - copies only the variable associated with the type of the argument.
// It would be possible just to do a dumb copy of everything, but that seems opposed to
// the concept of an EMObject, which is always of a single type.
EMObject& EMObject::operator=( const EMObject& that )
{

// This test breaks assignment when either the current or assigned values are (float)nan
// it's also somewhat inherently stupid, since testing for equivalence may cost as much
// as assignment.
//	if ( *this != that )
	{
		// First store the type of the input, At first I forgot to do this and it was a very
		// difficult bug to track down
		type = that.type;

// 		if ( type != this_type ) throw


		switch (type)
		{
		case BOOL:
			b = that.b;
		break;
		case INT:
			n = that.n;
		break;
		case UNSIGNEDINT:
			ui = that.ui;
		break;
		case FLOAT:
			f = that.f;
		break;
		case DOUBLE:
			d = that.d;
		break;
		case CTF:
		case STRING:
			str = that.str;
		break;
		case FLOAT_POINTER:
			// Warning - Pointer address copy.
			fp = that.fp;
		break;
		case INT_POINTER:
		// Warning - Pointer address copy.
			ip = that.ip;
		break;
		case VOID_POINTER:
			// Warning - Pointer address copy.
			vp = that.vp;
		break;
		case EMDATA:
			// Warning - Pointer address copy.
			emdata = that.emdata;
		break;
		case XYDATA:
			// Warning - Pointer address copy.
			xydata = that.xydata;
		break;
		case TRANSFORM:
		case FLOATARRAY:
			farray = that.farray;
		break;
		case INTARRAY:
			iarray = that.iarray;
		break;
		case STRINGARRAY:
			strarray = that.strarray;
		break;
		case UNKNOWN:
			// This is possible, nothing should happen
			// The EMObject's default constructor has been called and
			// as yet has no type - doing nothing is exactly as the
			// the assignment operator should work.
		break;
		default:
			LOGERR("No such EMObject defined");
			throw NotExistingObjectException("EMObject", "unknown type");
		break;
		}
	}
//	else
//	{
//		cerr << "Warning - attempt to assign EMObject onto itself. No action taken" << endl;
//		cerr << "My type is " << get_object_type_name(type) << endl;
//	}

	return *this;
}

//-------------------------------TypeDict--------------------------------------------

void TypeDict::dump()
{
	map < string, string >::iterator p;
	for (p = type_dict.begin(); p != type_dict.end(); p++) {
		printf("\t%s    %s  %s\n",
			   p->first.c_str(), p->second.c_str(), desc_dict[p->first].c_str());
	}
}

//-------------------------------Dict--------------------------------------------

Dict::Dict(const Dict& that)
{
	*this = that;
}

Dict& Dict::operator=(const Dict& that)
{
	if ( this != &that )
	{
		dict.clear();
		copy(that.begin(), that.end(), inserter(dict, dict.begin()));
		// or use this
		// dict.insert( that.begin(), that.end());
	}
	else
	{
		cerr << "Warning - attempted to assign a Dict object to itself. No action taken" << endl;
	}

	return *this;
}

bool EMAN::operator==(const Dict& d1, const Dict& d2)
{
	// Just make use of map's version of operator==
	return (d1.dict == d2.dict);
}

bool EMAN::operator!=(const Dict& d1, const Dict& d2)
{
	return !(d1 == d2);
}


// Iterator support
// This is just a wrapper, everything is inherited from the map<string,EMObject>::iterator
// so the interface is the same as you would expect
// iterator support added by d.woolford May 2007

Dict::iterator Dict::begin( void )
{
	return iterator( dict.begin() );
}

Dict::const_iterator Dict::begin( void ) const
{
	return const_iterator( (map < string, EMObject >::const_iterator) dict.begin() );
}

// Wraps map.find(const string& key)
Dict::iterator Dict::find( const string& key )
{
	return iterator( dict.find(key) );
}

Dict::iterator Dict::end( void )
{
	return iterator( dict.end() );
}

Dict::const_iterator Dict::end( void ) const
{
	return const_iterator( (map < string, EMObject >::const_iterator)dict.end() );
}

Dict::const_iterator Dict::find( const string& key ) const
{
	return const_iterator( (map < string, EMObject >::const_iterator)dict.find(key) );
}

//
// iterator
//
Dict::iterator::iterator( map< string, EMObject >::iterator parent_it  ) :
	map< string, EMObject >::iterator( parent_it )
{
}


Dict::iterator::iterator( const iterator& that ) :
	map < string, EMObject >::iterator( that )
{
}


Dict::iterator& Dict::iterator::operator=( const iterator& that )
{
	if( this != &that )
	{
		map < string, EMObject >::iterator::operator=( that );
	}
	return *this;
}

//
// const_iterator
//

Dict::const_iterator::const_iterator( const map < string, EMObject >::const_iterator parent_it  ) :
	map< string, EMObject >::const_iterator( parent_it )
{
}

Dict::const_iterator::const_iterator( const Dict::iterator& it ) :
	map< string, EMObject >::const_iterator(it)
{
}

Dict::const_iterator::const_iterator( const const_iterator& it ) :
	map< string, EMObject >::const_iterator(it)
{
}

Dict::const_iterator& Dict::const_iterator::operator=( const const_iterator& that )
{
	if( this != &that )
	{
		map < string, EMObject >::const_iterator::operator=( that );
	}
	return *this;
}

EMObject Dict::get_ci(const string & key) const
{
	string lower_key = Util::str_to_lower(key);

	for (map < string, EMObject >::const_iterator it = dict.begin(); it != dict.end(); ++it ) {
		string lower = Util::str_to_lower(it->first);
		if (lower == lower_key) return it->second;
	}

	throw NotExistingObjectException("EMObject", "Nonexisting key (" + key + ") in Dict");
}

bool Dict::has_key_ci(const string & key) const
{
	string lower_key = Util::str_to_lower(key);

	for (map < string, EMObject >::const_iterator it = dict.begin(); it != dict.end(); ++it ) {
		string lower = Util::str_to_lower(it->first);
		if (lower == lower_key) return true;
	}
	return false;
}
