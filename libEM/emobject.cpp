/**
 * $Id$
 */

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

#include "emobject.h"
#include <cmath>
#ifdef WIN32
#define M_PI 3.14159265358979323846f
#endif

using namespace EMAN;


const float EMConsts::I2G = (float) (4.0 / (M_PI*M_PI));  
const float EMConsts::I3G = (float) (6.4 / (M_PI*M_PI));  
const float EMConsts::I4G = (float) (8.8 / (M_PI*M_PI));  
const float EMConsts::I5G = (float) (10.4 / (M_PI*M_PI)); 

EMObject::operator int () const
{
	if (type == INT) {
		return n;
	}
	else if (type == FLOAT) {
		return (int) f;
	}
	else if (type == DOUBLE) {
		return (int) d;
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
	if (type == FLOAT) {
		return f;
	}
	else if (type == INT) {
		return (float) n;
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
	if (type == DOUBLE) {
		return d;
	}
	else if (type == INT) {
		return (double) n;
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

EMObject::operator  const char *() const
{
	if (type != STRING) {
		if (type != UNKNOWN) {
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

EMObject::operator  XYData * () const
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

EMObject::operator  Transform3D *() const
{
	if(type != TRANSFORM3D) {
		if(type != UNKNOWN) {
			throw TypeException("Cannot convert to TRANSFORM3D* from this data type",
				   get_object_type_name(type));
		}
	}
	return transform3d;
}

EMObject::operator shared_ptr< vector<int> >() const
{
    if( type != INTARRAY )
    {
        if( type != UNKNOWN )
	{
	    throw TypeException("Cannot convert to int array from ", get_object_type_name(type) );
	}

	return shared_ptr< vector<int> >();
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
	if (type == STRING) {
		return str;
	}
	else {
		char tmp_str[32];

		if (type == INT) {
			sprintf(tmp_str, "%d", n);
		}
		else if (type == FLOAT) {
			sprintf(tmp_str, "%f", f);
		}
		else if (type == DOUBLE) {
			sprintf(tmp_str, "%f", d);
		}
		else if (type == EMDATA) {
			sprintf(tmp_str, "EMDATA");
		}
		else if (type == XYDATA) {
			sprintf(tmp_str, "XYDATA");
		}
		else {
			LOGERR("No such EMObject defined");
			throw NotExistingObjectException("EMObject", "unknown type");
		}
		return string(tmp_str);
	}
}

EMObject::ObjectType EMObject::get_type() const
{
	return type;
}


const char *EMObject::get_object_type_name(ObjectType t)
{
	switch (t) {
	case INT:
		return "INT";
	case FLOAT:
		return "FLOAT";
	case DOUBLE:
		return "DOUBLE";
	case STRING:
		return "STRING";
	case EMDATA:
		return "EMDATA";
	case XYDATA:
		return "XYDATA";
	case TRANSFORM3D:
		return "TRANSFORM3D";
	case INTARRAY:
		return "INTARRAY";
	case FLOATARRAY:
		return "FLOATARRAY";
	case STRINGARRAY:
		return "STRINGARRAY";
	case UNKNOWN:
		LOGERR("No such EMObject defined");
		throw NotExistingObjectException("EMObject", "unknown type");
	}

	return "UNKNOWN";
}


bool EMAN::operator==(const EMObject &e1, const EMObject & e2)
{
#if 0
	if (e1.type != e2.type) {
		return false;
	}
#endif
	switch (e1.type) {
	case EMObject::INT:
		return (e1.n == e2.n);
	case EMObject::FLOAT:
		return (e1.f == e2.f);
	case EMObject::DOUBLE:
		return (e1.d == e2.d);
	case EMObject::STRING:
		return (e1.str == e2.str);
	case EMObject::EMDATA:
		return (e1.emdata == e2.emdata);
	case EMObject::XYDATA:
		return (e1.xydata == e2.xydata);
	case EMObject::FLOATARRAY:
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
	case EMObject::STRINGARRAY:
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
	default:
		return false;
	}
	return false;
}

bool EMAN::operator!=(const EMObject &e1, const EMObject & e2)
{
	return !(e1 == e2);
}


void TypeDict::dump() 
{
	map < string, string >::iterator p;
	for (p = type_dict.begin(); p != type_dict.end(); p++) {
		printf("\t%s    %s  %s\n",
			   p->first.c_str(), p->second.c_str(), desc_dict[p->first].c_str());
	}
}
