#include "emobject.h"

using namespace EMAN;

EMObject::operator  int ()
	 const
	 {
		 if (type == INT)
{
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
				 Log::logger()->error("type error. Cannot convert to int from data type '%s'",
									  get_object_type_name(type));
			 }
		 }
		 return 0;
	 }

EMObject::operator  float ()
	 const
	 {
		 if (type == FLOAT)
{
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
				 Log::logger()->
					 error("type error. Cannot convert to float from data with type '%s'",
						   get_object_type_name(type));
			 }
		 }

		 return 0;
	 }

EMObject::operator  double ()
	 const
	 {
		 if (type == DOUBLE)
{
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
				 Log::logger()->error("type error. Cannot convert to double from data type '%s'",
									  get_object_type_name(type));
			 }
		 }
		 return 0;
	 }

EMObject::operator  const char *()
	 const
	 {
		 if (type != STRING)
		 {
			 if (type != UNKNOWN) {
				 Log::logger()->error("type error. Cannot convert to string from data type '%s'",
									  get_object_type_name(type));
			 }
			 return "";
		 }
		 return str.c_str();
	 }

EMObject::operator  EMData * ()
	 const
	 {
		 if (type != EMDATA)
		 {
			 if (type != UNKNOWN) {
				 Log::logger()->error("type error. Cannot convert to EMData* from data type '%s'",
									  get_object_type_name(type));
			 }
			 return 0;
		 }
		 return emdata;
	 }

EMObject::operator  XYData * ()
	 const
	 {
		 if (type != XYDATA)
		 {
			 if (type != UNKNOWN) {
				 Log::logger()->error("type error. Cannot convert to XYData* data type '%s'",
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
			Log::logger()->error("type error. Cannot call get_farray for data type '%s'",
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
			sprintf(tmp_str, "Unknown");
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
		return "UNKNOWN";
	}

	return "UNKNOWN";
}


void TypeDict::dump() const
{
	map < string, string >::const_iterator p;
	for (p = dict.begin(); p != dict.end(); p++) {
		printf("%20s    %s\n", p->first.c_str(), p->second.c_str());
	}
}
