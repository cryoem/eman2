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
	enum ObjectType {
	    INT,
	    FLOAT,
	    DOUBLE,
	    STRING,
	    EMDATA,
	    FLOATARRAY,
	    UNKNOWN
	};

    public:
	EMObject() : type(UNKNOWN) {}
	EMObject(int num) : n(num), type(INT) {}
	EMObject(float ff) : f(ff), type(FLOAT) {}
	EMObject(double dd) : d(dd), type(DOUBLE) {}
	EMObject(string s) : str(s), type(STRING) {}
	EMObject(EMData* em) : emdata(em), type(EMDATA) {}
	EMObject(const vector<float>& v) : farray(v), type(FLOATARRAY) {}

	~EMObject() { }

	int get_int() const
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
		    Log::logger()->error("type error. Cannot call get_int() for data type '%s'", get_object_type_name(type));
		}
	    }
	    return 0;
	}
	
	float get_float() const
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
		    Log::logger()->error("type error. Cannot call get_float() for data type '%s'", get_object_type_name(type));
		}
	    }
    
	    return 0;
	}
	
	double get_double() const
	{
	    if (type == DOUBLE) {
		return d;
	    }
	    else if (type == INT) {
		return (double)n;
	    }
	    else if (type == FLOAT) {
		return (double) f;
	    }
	    else {
		if (type != UNKNOWN) {
		    Log::logger()->error("type error. Cannot call get_double() for data type '%s'", get_object_type_name(type));
		}
	    }    
	    return 0;
	}
	
	string get_string() const
	{
	    if (type != STRING) {
		if (type != UNKNOWN) {
		    Log::logger()->error("type error. Cannot call get_string() for data type '%s'", get_object_type_name(type));
		}
		return "";
	    }    
	    return str;
	}

	    
	EMData* get_EMData() const
	{
	    if (type != EMDATA) {
		if (type != UNKNOWN) {
		    Log::logger()->error("type error. Cannot call get_EMData() for data type '%s'", get_object_type_name(type));
		}
		return 0;
	    }    
	    return emdata;
	}

	vector<float> get_farray() const
	{
	    if (type != FLOATARRAY) {
		if (type != UNKNOWN) {
		    Log::logger()->error("type error. Cannot call get_farray for data type '%s'", get_object_type_name(type));
		}
		return vector<float>();
	    }
	    return farray;
	}
	
	bool is_null() const { return (type == UNKNOWN); }
	
	string to_str() const;
 
	static const char* get_object_type_name(ObjectType t) 
	{
	    switch(t) {
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
	    case FLOATARRAY:
		return "FLOATARRAY";
	    case UNKNOWN:
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
	vector<float> farray;
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

    class TypeDict {
    public:
	TypeDict() {}
	~TypeDict() {}

	vector<string> keys() const
	{
	    vector<string> result;
    
	    map<string, string>::const_iterator p = 0;
	    for (p = dict.begin(); p != dict.end(); p++) {
		result.push_back(p->first);
	    }
    
	    return result;
	}
	
	int size() const { return dict.size(); }

	void put(string key, EMObject::ObjectType o)
	{
	    dict[key] = EMObject::get_object_type_name(o);
	}
	
	string get(string key) { return dict[key]; }
	
	string operator[](const string& key) { return dict[key]; }

    private:
	map<string, string> dict;
    };
    
}

#endif
