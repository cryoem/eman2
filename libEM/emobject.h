/**
 * $Id$
 */
#ifndef eman__object__h__
#define eman__object__h__ 1

#include <string>
#include <map>
#include <vector>
#include <stdio.h>
#include "log.h"
#include "exception.h"

using std::vector;
using std::string;
using std::map;

namespace EMAN
{
	class EMConsts {
	public:
		static const float I2G; // 2 interpolation
		static const float I3G; // used for 3 and 5x5x5 interpolation
		static const float I4G; // used for 4 interpolation
		static const float I5G; // used for 5x5x5 interpolation
	};
	
	class EMData;
	class XYData;
	class Aligner;
	class Averager;
	class Cmp;
	class Processor;
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

	/** EMObject is a wrapper class for types including int, float,
     * double, etc as defined in ObjectType. Each type is typically used 
     * as follows ('int' is the example):
     *
     *    int a = 12;
     *    EMObject o(a);
     *    EMObject o2 = a; // implicit converter from int to EMObject. 
     *    int a1 = o;      // implicit converter from EMObject to int.
     */
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
			STRINGARRAY,
			UNKNOWN
		};

	public:
		EMObject():type(UNKNOWN)
		{
			n = 0;
			f = 0;
			d = 0;
			emdata = 0;
			xydata = 0;
		}

		EMObject(int num):n(num), emdata(0), xydata(0), type(INT)
		{
		}
		EMObject(float ff):f(ff), emdata(0), xydata(0), type(FLOAT)
		{
		}
		EMObject(double dd):d(dd), emdata(0), xydata(0), type(DOUBLE)
		{
		}
		EMObject(const char *s): n(0), emdata(0), xydata(0), str(string(s)), type(STRING)
		{
		}
		EMObject(const string & s):n(0), emdata(0), xydata(0), str(s), type(STRING)
		{
		}
		EMObject(EMData * em):n(0), emdata(em), xydata(0), type(EMDATA)
		{
		}
		EMObject(XYData * xy):n(0),  emdata(0), xydata(xy),type(XYDATA)
		{
		}
		EMObject(const vector < float >&v)
			:n(0), emdata(0), xydata(0), farray(v),type(FLOATARRAY)
		{
		}

		EMObject(const vector <string>& sarray)
			:n(0),emdata(0),xydata(0),strarray(sarray),type(STRINGARRAY)
		{
		}
		
		~EMObject() {
		}

		operator  int () const;
		operator  float () const;
		operator  double () const;
		operator  const char *() const;
		operator  EMData *() const;
		operator  XYData *() const;

		operator vector < float > () const;
		operator vector<string> () const;
		
		bool is_null() const;
		string to_str() const;
		ObjectType get_type() const;
		static const char *get_object_type_name(ObjectType t);

		friend bool operator==(const EMObject &e1, const EMObject & e2);
		friend bool operator!=(const EMObject &e1, const EMObject & e2);
		
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
		vector < float >farray;
		vector < string> strarray;
		ObjectType type;
	};

	bool operator==(const EMObject &e1, const EMObject & e2);
	bool operator!=(const EMObject &e1, const EMObject & e2);

	
	/** Dict is a dictionary to store <string, EMObject> pair.
     * Typical ways to construct a Dict:
     *
     *      Dict d;
     *      d["lowpass"] = 12.23;
     *      float lowpass1 = d["lowpass"];
     *
     *      Dict d2("lowpass", 12.23);
     */
	class Dict
	{
	public:
		Dict()
		{
		}

		Dict(const string & key1, EMObject val1)
		{
			dict[key1] = val1;
		}

		Dict(const string & key1, EMObject val1,
			 const string & key2, EMObject val2)
		{
			dict[key1] = val1;
			dict[key2] = val2;
		}

		Dict(const string & key1, EMObject val1,
			 const string & key2, EMObject val2,
			 const string & key3, EMObject val3)
		{
			dict[key1] = val1;
			dict[key2] = val2;
			dict[key3] = val3;
		}

		Dict(const string & key1, EMObject val1,
			 const string & key2, EMObject val2,
			 const string & key3, EMObject val3,
			 const string & key4, EMObject val4)
		{			
			dict[key1] = val1;
			dict[key2] = val2;
			dict[key3] = val3;
			dict[key4] = val4;
		}


		Dict(const map < string, EMObject > &d)
		{
			map < string, EMObject >::const_iterator p;
			for (p = d.begin(); p != d.end(); p++) {
				dict[p->first] = p->second;
			}
		}

		~Dict() {
		}

		vector < string > keys()const
		{
			vector < string > result;

			map < string, EMObject >::const_iterator p;
			for (p = dict.begin(); p != dict.end(); p++) {
				result.push_back(p->first);
			}

			return result;
		}

		vector < EMObject > values()const
		{
			vector < EMObject > result;

			map < string, EMObject >::const_iterator p;
			for (p = dict.begin(); p != dict.end(); p++) {
				result.push_back(p->second);
			}

			return result;
		}

		bool has_key(const string & key) const
		{
			map < string, EMObject >::const_iterator p = dict.find(key);
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

		void put(const string & key, EMObject val)
		{
			dict[key] = val;
		}

		EMData *set_default(const string & key, EMData * val)
		{
			if (!has_key(key)) {
				dict[key] = EMObject(val);
			}
			return dict[key];
		}

		int set_default(const string & key, int val)
		{
			if (!has_key(key)) {
				dict[key] = val;
			}
			return dict[key];
		}

		float set_default(const string & key, float val)
		{
			if (!has_key(key)) {
				dict[key] = val;
			}
			return dict[key];
		}

		double set_default(const string & key, double val)
		{
			if (!has_key(key)) {
				dict[key] = val;
			}
			return dict[key];
		}

		map < string, EMObject > &get_dict() {
			return dict;
		}

		map < string, EMObject > get_dict()const
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
		mutable map < string, EMObject > dict;
	};

	/** TypeDict is a dictionary to store <string, EMObject::ObjectType> pair.
     * It is mainly used to store processor-like class's parameter
     * information: <parameter-name, parameter-type>.
     * Typical usage of this class:
     *
     *    TypeDict d;
     *    d.put("with", EMObject::EMDATA);
     *    d.put("lowpass", EMObject::FLOAT);
     *
     *    string lowpass_type = d["lowpass"];
     */
	class TypeDict
	{
	public:
		TypeDict()
		{
		}

		~TypeDict()
		{
		}

		vector < string > keys() const
		{
			vector < string > result;
			map < string, string >::const_iterator p;

			for (p = type_dict.begin(); p != type_dict.end(); p++) {
				result.push_back(p->first);
			}

			return result;
		}

		size_t size() const
		{
			return type_dict.size();
		}

		void put(const string& key, EMObject::ObjectType o, const string& desc = "")
		{
			type_dict[key] = EMObject::get_object_type_name(o);
			desc_dict[key] = desc;
		}

		string get_type(const string& key)
		{
			return type_dict[key];
		}

		string get_desc(const string& key)
		{
			return desc_dict[key];
		}
		
		string operator[] (const string & key)
		{
			return type_dict[key];
		}

		void dump();

	private:
		map < string, string > type_dict;
		map < string, string > desc_dict;
	};

	/** Factory is used to store objects to create new instances.
     * It is a singleton template. Typical usages are as follows:
     *
     *   1. How to define a new factory (e.g. Processor Factory):
     *   
     *      template<> Factory<Processor>::Factory()
     *      {
     *         force_add(&AbsoluateValueProcessor::NEW);
     *         force_add(&BooleanProcessor::NEW);
     *      }
     *
     *
     *   2. How to use a Factory (e.g. Processor Factory):
     *
     *	    Filter *f1 = Factory<Processor>::get("AbsoluateValue");
     *      Filter *f2 = Factory<Processor>::get("LowpassGauss", Dict("lowpass", EMObject(12));
     */
	template < class T > class Factory
	{
	public:
		typedef T *(*InstanceType) ();

		static void add(InstanceType i);
		static T *get(const string & instance_name);
		static T *get(const string & instance_name, const Dict & params);
		static vector < string > get_list();

	private:
		Factory();
		Factory(const Factory < T > &);
		~Factory();
		static void init();
		void force_add(InstanceType i);

		static Factory < T > *my_instance;
		map < string, InstanceType > my_dict;
	};

	template < class T > Factory < T > *Factory < T >::my_instance = 0;

	template < class T > void Factory < T >::init()
	{
		if (!my_instance) {
			my_instance = new Factory < T > ();
		}
	}

	template < class T > void Factory < T >::force_add(InstanceType new_instance)
	{
		T *i = new_instance();
		string name = i->get_name();
		my_dict[name] = new_instance;
		if( i )
		{
			delete i;
			i = 0;
		}
	}


	template < class T > void Factory < T >::add(InstanceType new_instance)
	{
		init();

		T *i = new_instance();
		string name = i->get_name();
		typename map < string, InstanceType >::iterator fi =
			my_instance->my_dict.find(name);

		if (fi == my_instance->my_dict.end()) {
			my_instance->my_dict[name] = new_instance;
		}
		if( i )
		{
			delete i;
			i = 0;
		}
	}

	template < class T > T * Factory < T >::get(const string & instancename)
	{
		init();
		typename map < string, InstanceType >::iterator fi =
			my_instance->my_dict.find(instancename);
		if (fi != my_instance->my_dict.end()) {
			return my_instance->my_dict[instancename] ();
		}
		
		throw NotExistingObjectException(instancename, "No such an instance existing");
	}

	template < class T > T * Factory < T >::get(const string & instancename,
												const Dict & params)
	{
		init();

		typename map < string, InstanceType >::iterator fi =
			my_instance->my_dict.find(instancename);

		if (fi != my_instance->my_dict.end()) {
			T *i = my_instance->my_dict[instancename] ();
			i->set_params(params);
			return i;
		}		

		throw NotExistingObjectException(instancename, "No such an instance existing");
	}

	template < class T > vector < string > Factory < T >::get_list() {
		init();
		vector < string > result;
		typename map < string, InstanceType >::const_iterator p;
		for (p = my_instance->my_dict.begin(); p != my_instance->my_dict.end(); p++) {
			result.push_back(p->first);
		}

		return result;
	}

	template < class T > void dump_factory()
	{
		vector < string > item_names = Factory < T >::get_list();

		for (size_t i = 0; i < item_names.size(); i++) {
			T *item = Factory < T >::get(item_names[i]);
			printf("%s :  %s\n", item->get_name().c_str(),item->get_desc().c_str());
			TypeDict td = item->get_param_types();
			td.dump();
		}
	}
	
}

#endif
