/**
 * $Id$
 */
#ifndef eman__obejct__em__
#define eman__obejct__em__ 1

#include <string>
#include <map>
#include <vector>
#include "log.h"

using std::vector;
using std::string;
using std::map;

namespace EMAN
{    
    class EMData;
    class XYData;
    class Aligner;
    class Averager;
    class Cmp;
    class Filter;
    class Projector;
    class Reconstructor;
    
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

    class EMObject
    {
    public:
	enum ObjectType {
	    INT,
	    FLOAT,
	    DOUBLE,
	    STRING,
	    EMDATA,
	    XYDATA,
	    FLOATARRAY,
	    UNKNOWN
	};

    public:
	EMObject() : type(UNKNOWN)
	{
	}
	EMObject(int num) : n(num), type(INT)
	{
	}
	EMObject(float ff) : f(ff), type(FLOAT)
	{
	}
	EMObject(double dd) : d(dd), type(DOUBLE)
	{
	}
	EMObject(string s) : str(s), type(STRING)
	{
	}
	EMObject(EMData * em) : emdata(em), type(EMDATA)
	{
	}
	EMObject(XYData * xy) : xydata(xy), type(XYDATA)
	{
	}
	EMObject(const vector<float> & v) : farray(v), type(FLOATARRAY)
	{
	}

	~EMObject() {
	}

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
		    Log::logger()->error("type error. Cannot call get_int() for data type '%s'",
					 get_object_type_name(type));
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
		    Log::logger()->error("type error. Cannot call get_float() for data type '%s'",
					 get_object_type_name(type));
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
		return (double) n;
	    }
	    else if (type == FLOAT) {
		return (double) f;
	    }
	    else {
		if (type != UNKNOWN) {
		    Log::logger()->error("type error. Cannot call get_double() for data type '%s'",
					 get_object_type_name(type));
		}
	    }
	    return 0;
	}

	string get_string() const
	{
	    if (type != STRING) {
		if (type != UNKNOWN) {
		    Log::logger()->error("type error. Cannot call get_string() for data type '%s'",
					 get_object_type_name(type));
		}
		return "";
	    }
	    return str;
	}


	EMData *get_emdata() const
	{
	    if (type != EMDATA) {
		if (type != UNKNOWN) {
		    Log::logger()->error("type error. Cannot call get_emdata() for data type '%s'",
					 get_object_type_name(type));
		}
		return 0;
	    }
	    return emdata;
	}

	XYData *get_xydata() const
	{
	    if (type != XYDATA) {
		if (type != UNKNOWN) {
		    Log::logger()->error("type error. Cannot call get_xydata() for data type '%s'",
					 get_object_type_name(type));
		}
		return 0;
	    }
	    return xydata;
	}

	vector<float> get_farray() const
	{
	    if (type != FLOATARRAY) {
		if (type != UNKNOWN) {
		    Log::logger()->error("type error. Cannot call get_farray for data type '%s'",
					 get_object_type_name(type));
		}
		return vector<float> ();
	    }
	    return farray;
	}

	bool is_null() const
	{
	    return (type == UNKNOWN);
	}

	string to_str() const
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
	
	static const char *get_object_type_name(ObjectType t)
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

    private:
	union
	{
	    int n;
	    float f;
	    double d;
	};

	EMData *emdata;
	XYData *xydata;
	string str;
	vector<float> farray;
	ObjectType type;
    };


    class Dict
    {
    public:
	Dict()
	{
	}

	Dict(const map<string, EMObject> & d)
	{
	    map<string, EMObject>::const_iterator p;
	    for (p = d.begin(); p != d.end(); p++) {
		dict[p->first] = p->second;
	    }
	}

	~Dict() {
	}

	vector<string> keys() const
	{
	    vector<string> result;

	    map<string, EMObject>::const_iterator p;
	    for (p = dict.begin(); p != dict.end(); p++) {
		result.push_back(p->first);
	    }

	    return result;
	}

	vector<EMObject> values() const
	{
	    vector<EMObject> result;

	    map<string, EMObject>::const_iterator p;
	    for (p = dict.begin(); p != dict.end(); p++) {
		result.push_back(p->second);
	    }

	    return result;
	}

	bool has_key(const string & key) const
	{
	    map<string, EMObject>::const_iterator p = dict.find(key);
	    if (p != dict.end()) {
		return true;
	    }
	    return false;
	}

	size_t size() const
	{
	    return dict.size();
	}

	EMObject get(const string & key)
	{
	    return dict[key];
	}

	void put(string key, EMObject val)
	{
	    dict[key] = val;
	}

	map<string, EMObject> & get_dict() {
	    return dict;
	}

	map<string, EMObject> get_dict()const
	{
	    return dict;
	}

	EMObject & operator[] (const string & key)
	{
	    return dict[key];
	}

	EMObject operator[] (const string & key) const
	{
	    return dict[key];
	}

    private:
	mutable map<string, EMObject> dict;
    };

    class TypeDict
    {
    public:
	TypeDict()
	{
	}
	
	~TypeDict()
	{
	}

	vector<string> keys()const
	{
	    vector<string> result;
	    map<string, string>::const_iterator p;
	    
	    for (p = dict.begin(); p != dict.end(); p++) {
		result.push_back(p->first);
	    }

	    return result;
	}

	size_t size() const
	{
	    return dict.size();
	}

	void put(string key, EMObject::ObjectType o)
	{
	    dict[key] = EMObject::get_object_type_name(o);
	}

	string get(string key)
	{
	    return dict[key];
	}

	string operator[] (const string & key)
	{
	    return dict[key];
	}

    private:
	map<string, string> dict;
    };

    template <class T> class Factory {
    public:
	typedef T *(*InstanceType)();
	
	static Factory<T> *instance();

	void add(InstanceType i);
	T *get(string instance_name);
	T *get(string instance_name, const Dict & params);
	
	vector<string> get_list();

    private:
	Factory();
	Factory(const Factory<T> &);
	~Factory();
	
	void force_add(InstanceType i);

	static Factory<T> *my_instance;
	map<string, InstanceType> my_dict;
    };
    
    template <class T> Factory<T> *Factory<T>::my_instance = 0;

    template <class T> Factory<T> *Factory<T>::instance()
    {
	if (!my_instance) {
	    my_instance = new Factory<T>();
	}
	return my_instance;
    }

    template <class T> void Factory<T>::force_add(InstanceType new_instance)
    {
	T *i = new_instance();
	string name = i->get_name();
	my_dict[name] = new_instance;
	delete i;
	i = 0;
    }


    template <class T> void Factory<T>::add(InstanceType new_instance)
    {
	T *i = new_instance();
	string name = i->get_name();
    
	typename map<string, InstanceType>::iterator fi = my_dict.find(name);

	if (fi != my_dict.end()) {
	    my_dict[name] = new_instance;
	}
	delete i;
	i = 0;
    }

    template <class T> T *Factory<T>::get(string instancename)
    {
	typename map<string, InstanceType>::iterator fi = my_dict.find(instancename);
	if (fi != my_dict.end()) {
	    return my_dict[instancename]();
	}
	Log::logger()->error("No such an instance existing: %s", instancename.c_str());
	return 0;
    }

    template <class T> T *Factory<T>::get(string instancename, const Dict & params)
    {
	typename map<string, InstanceType>::iterator fi = my_dict.find(instancename);
    
	if (fi != my_dict.end()) {
	    T *i = my_dict[instancename]();
	    i->set_params(params);
	    return i;
	}
	Log::logger()->error("No such an instance existing: %s", instancename.c_str());
	return 0;
    }

    template <class T> vector<string> Factory<T>::get_list()
    {
	vector<string> result;
	typename map<string, InstanceType>::const_iterator p;
	for (p = my_dict.begin(); p != my_dict.end(); p++) {
	    result.push_back(p->first);
	}

	return result;
    }


}

#endif
