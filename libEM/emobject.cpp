#include "emobject.h"
#include "exception.h"
#include <math.h>
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
			LOGERR("type error. Cannot convert to int from data type '%s'",
				   get_object_type_name(type));
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
			LOGERR("type error. Cannot convert to float from data type '%s'",
					  get_object_type_name(type));
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
			LOGERR("type error. Cannot convert to double from data type '%s'",
				   get_object_type_name(type));
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
			LOGERR("type error. Cannot convert to string from data type '%s'",
				   get_object_type_name(type));
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
			LOGERR("type error. Cannot convert to EMData* from data type '%s'",
				   get_object_type_name(type));
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
			LOGERR("type error. Cannot convert to XYData* data type '%s'",
				   get_object_type_name(type));
			throw TypeException("Cannot convert to XYData* from this data type",
				   get_object_type_name(type));
		}
		return 0;
	}
	return xydata;
}

vector < float >EMObject::get_farray() const
{
	if (type != FLOATARRAY) {
		if (type != UNKNOWN) {
			LOGERR("type error. Cannot call get_farray for data type '%s'",
				   get_object_type_name(type));
			throw TypeException("Cannot call get_farray for this data type",
								get_object_type_name(type));
		}
		return vector < float >();
	}
	return farray;
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
	case FLOATARRAY:
		return "FLOATARRAY";
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
		printf("%20s    %s  %s\n",
			   p->first.c_str(), p->second.c_str(), desc_dict[p->first].c_str());
	}
}
