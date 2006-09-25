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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * 
 * */
 
#ifndef eman__object__h__
#define eman__object__h__ 1

#include <map>
#include <vector>
#include <iostream>
#include <boost/shared_ptr.hpp>

#include "log.h"
#include "exception.h"

using std::vector;
using std::string;
using std::map;

using boost::shared_ptr;

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
	class Analyzer;

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
			INTARRAY,
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

		EMObject(const vector< int >& v )
		    : n(0), emdata(0), xydata(0), iarray( new vector<int>(v) ), type(INTARRAY)
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

                operator shared_ptr< vector<int> >() const;
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
		shared_ptr< vector<int> > iarray;
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
			if( has_key(key) ) {
				return dict[key];
			}
			else {
				LOGERR("No such key exist in this Dict");
				throw NotExistingObjectException("EMObject", "Nonexisting key (" + key + ") in Dict");
			}
		}

		void put(const string & key, EMObject val)
		{
			dict[key] = val;
		}

		void erase(const string & key)
		{
			dict.erase(key);
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
//			static EMObject nullreturn;
//			if( has_key(key) )  return dict[key];
//			else return nullreturn;

//			if( has_key(key) ) {
				return dict[key];
//			}
//			else {
//				LOGERR("No such key exist in this Dict");
//				throw NotExistingObjectException("EMObject", "Nonexisting key (" + key + ") in Dict");
//			}
		}

		EMObject operator[] (const string & key) const
		{
//			if( has_key(key) )  return dict[key];
//			else return EMObject();
			return dict[key];

//			else {
//				LOGERR("No such key exist in this Dict");
//				throw NotExistingObjectException("EMObject", "Nonexisting key (" + key + ") in Dict");
//			}
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
     *	    Processor *f1 = Factory<Processor>::get("eman1.math.absvalue");
     *      Processor *f2 = Factory<Processor>::get("eman1.filter.lowpass.gaussian", Dict("lowpass", EMObject(12));
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
			
			const vector<string> para_keys = params.keys();
			const vector<string> valid_keys = i->get_param_types().keys();
			typename vector<string>::const_iterator it;
			for(it=para_keys.begin(); it!=para_keys.end(); ++it) {
				if( find(valid_keys.begin(), valid_keys.end(), *it) == valid_keys.end() ) {
					throw InvalidParameterException(*it);
				}
			}
			
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
	
	template < class T > map<string, vector<string> > dump_factory_list()
	{
		vector < string > item_names = Factory < T >::get_list();
		map<string, vector<string> >	factory_list;
		
		typename vector<string>::const_iterator p;
		for(p = item_names.begin(); p !=item_names.end(); ++p) {
			T *item = Factory<T>::get(*p);
			
			string name = item->get_name();
			
			vector<string> content;
			content.push_back(item->get_desc());
			TypeDict td = item->get_param_types();
			vector<string> keys = td.keys();
			for(unsigned int i=0; i<td.size(); ++i) {
				content.push_back(keys[i]);
				content.push_back( td.get_type(keys[i]) );
				content.push_back( td.get_desc(keys[i]) );
			}
			factory_list[name] = content;
		}
		
		return factory_list;
	}
}

#endif
