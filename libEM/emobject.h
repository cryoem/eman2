#ifndef eman__obejct__em__
#define eman__obejct__em__ 1

#include <string>
#include <map>
#include <vector>
#include "log.h"

using std::vector;
using std::string;
using std::map;

namespace EMAN {

    enum MapInfoType {
	NORMAL,
	ICOS2F_FIRST_OCTANT,
	ICOS2F_FULL,
	ICOS2F_HALF,
	ICOS3F_HALF,
	ICOS3F_FULL,
	ICOS5F_HALF,
	ICOS5F_FULL,
	ICOS_UNKNOWN
    };

    class EMData;
    
    class EMObject {
    public:
	EMObject() : type(UNKNOWN_OBJECT) {}
	EMObject(int num) : n(num), type(INT_OBJECT) {}
	EMObject(float ff) : f(ff), type(FLOAT_OBJECT) {}
	EMObject(double dd) : d(dd), type(DOUBLE_OBJECT) {}
	EMObject(string s) : str(s), type(STRING_OBJECT) {}
	EMObject(EMData* em) : emdata(em), type(EMDATA_OBJECT) {}
	
	~EMObject() { }

	int get_int() const
	{
	    if (type == INT_OBJECT) {
		return n;
	    }
	    else if (type == FLOAT_OBJECT) {
		return (int) f;
	    }
	    else if (type == DOUBLE_OBJECT) {
		return (int) d;
	    }
	    else {
		if (type != UNKNOWN_OBJECT) {
		    Log::logger()->error("type error. Cannot call get_int() for data type '%s'", get_object_type_name(type));
		}
	    }
	    return 0;
	}
	
	float get_float() const
	{
	    if (type == FLOAT_OBJECT) {
		return f;
	    }
	    else if (type == INT_OBJECT) {
		return (float) n;
	    }
	    else if (type == DOUBLE_OBJECT) {
		return (float) d;
	    }
	    else {
		if (type != UNKNOWN_OBJECT) {
		    Log::logger()->error("type error. Cannot call get_float() for data type '%s'", get_object_type_name(type));
		}
	    }
    
	    return 0;
	}
	
	double get_double() const
	{
	    if (type == DOUBLE_OBJECT) {
		return d;
	    }
	    else if (type == INT_OBJECT) {
		return (double)n;
	    }
	    else if (type == FLOAT_OBJECT) {
		return (double) f;
	    }
	    else {
		if (type != UNKNOWN_OBJECT) {
		    Log::logger()->error("type error. Cannot call get_double() for data type '%s'", get_object_type_name(type));
		}
	    }    
	    return 0;
	}
	
	string get_string() const
	{
	    if (type != STRING_OBJECT) {
		if (type != UNKNOWN_OBJECT) {
		    Log::logger()->error("type error. Cannot call get_string() for data type '%s'", get_object_type_name(type));
		}
		return "";
	    }    
	    return str;
	}

	    
	EMData* get_EMData() const
	{
	    if (type != EMDATA_OBJECT) {
		if (type != UNKNOWN_OBJECT) {
		    Log::logger()->error("type error. Cannot call get_EMData() for data type '%s'", get_object_type_name(type));
		}
		return 0;
	    }    
	    return emdata;
	}

	bool is_null() const { return (type == UNKNOWN_OBJECT); }
	
	string to_str() const;
   
    private:
	enum ObjectType {
	    INT_OBJECT,
	    FLOAT_OBJECT,
	    DOUBLE_OBJECT,
	    STRING_OBJECT,
	    EMDATA_OBJECT,
	    UNKNOWN_OBJECT
	};

	const char* get_object_type_name(ObjectType t) const
	{
	    switch(type) {
	    case INT_OBJECT:
		return "INT";
	    case FLOAT_OBJECT:
		return "FLOAT";
	    case DOUBLE_OBJECT:
		return "DOUBLE";
	    case STRING_OBJECT:
		return "STRING";
	    case EMDATA_OBJECT:
		return "EMDATA";
	    case UNKNOWN_OBJECT:
		return "UNKNOWN";
	    }
    
	    return "UNKNOWN";
	}
	
    private:
	union {
	    int n;
	    float f;
	    double d;
	};
	
	EMData* emdata;
	string str;
	ObjectType type;
    };


    class Dict {
    public:
	Dict() {}
	
	Dict(const map<string, EMObject>& d)
	{
	    map<string, EMObject>::const_iterator p;
	    for (p = d.begin(); p != d.end(); p++) {
		dict[p->first] = p->second;
	    }
	}
	
	~Dict() {}

	vector<string> keys() const
	{
	    vector<string> result;
    
	    map<string, EMObject>::const_iterator p = 0;
	    for (p = dict.begin(); p != dict.end(); p++) {
		result.push_back(p->first);
	    }
    
	    return result;
	}
	vector<EMObject> values() const
	{
	    vector<EMObject> result;
    
	    map<string, EMObject>::const_iterator p = 0;
	    for (p = dict.begin(); p != dict.end(); p++) {
		result.push_back(p->second);
	    }
    
	    return result;
	}

	bool has_key(const string& key) const
	{
	    map<string, EMObject>::const_iterator p = dict.find(key);
	    if (p != dict.end()) {
		return true;
	    }
	    return false;
	}

	int size() const { return dict.size(); }
	
	EMObject get(const string& key) { return dict[key]; }
	
	void put(string key, EMObject val) { dict[key] = val; }
	
	map<string, EMObject>& get_dict() { return dict; }
	
	map<string, EMObject> get_dict() const { return dict; }

	EMObject& operator[](const string& key) { return dict[key]; }
	
    private:
	map<string, EMObject> dict;
    };
    
}

#endif
