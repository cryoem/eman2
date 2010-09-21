/**
 * $Id$
 */

/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
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

#ifndef eman__object__h__
#define eman__object__h__ 1

#include <map>
using std::map;

#include <set>
using std::set;

#include <vector>
using std::vector;

#include <string>
using std::string;

#include <utility>
using std::pair;

#include <algorithm>
// using copy

#include <iterator>

#include "log.h"
#include "exception.h"

// #include "transform.h" // Trnasform3D::EulerType

// debug
#include <iostream>
using std::cout;
using std::endl;

#include <cctype> // tolower
#include <algorithm> //tolower
namespace EMAN
{
	class EMConsts {
	public:
		static const float I2G; // 2 interpolation
		static const float I3G; // used for 3 and 5x5x5 interpolation
		static const float I4G; // used for 4 interpolation
		static const float I5G; // used for 5x5x5 interpolation

		static const double rad2deg; // radians to degree constant factor
		static const double deg2rad; // degrees to radians constant factor
		static const double pi; // degrees to radians constant factor
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
	class Transform;
	class Ctf;

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
     *
     *  EMObjects may store pointers but they currently do not assume ownership - that
     *  is, the memory associated with a pointer is never freed by an EMObject.
     *
     * This type of class design is sometimes referred to as the Variant pattern.
     *
     * See the testing code in rt/emdata/test_emobject.cpp for prewritten testing code
     */
	class EMObject
	{
	public:
		enum ObjectType {
			UNKNOWN,
			BOOL,
			UNSIGNEDINT,
			INT,
			FLOAT,
			DOUBLE,
			STRING,
			EMDATA,
			XYDATA,
			INTARRAY,
			FLOATARRAY,
			STRINGARRAY,
			TRANSFORM,
			CTF,
			FLOAT_POINTER,
			INT_POINTER,
			VOID_POINTER
		};


		/** Constructors for each type
		 * More types could be added, but all of the accompanying functions would have
		 * to be altered to ensure correct functioning
		 */
		EMObject();
		EMObject(bool boolean);
		EMObject(int num);
		EMObject(unsigned int num);
		EMObject(float ff);
		EMObject(double dd);
		EMObject(const char *s);
		EMObject(const string & s);
		EMObject(float * fp);
		EMObject(int * ip);
		EMObject(void * vp);
		EMObject(EMData * em);
		EMObject(XYData * xy);
		EMObject(Transform * t);
		EMObject(Ctf * ctf);
		EMObject(const vector< int >& v );
		EMObject(const vector < float >&v);
		EMObject(const vector <string>& sarray);

		/** Copy constructor.
		 * copies pointer locations - does not take ownership
		 * deep copies all non pointer objects
		 */
		EMObject(const EMObject& that);

		/** Assigment operator
		 * copies pointer locations (emdata, xydata, transform3d) - does not take ownership
		 * deep copies all non pointer objects
		 */
		EMObject& operator=(const EMObject& that);

		/** Desctructor
		 * Does not free pointers.
		 */
		~EMObject();
		/** Conversion operators
		 */
		operator bool () const;
		operator int () const;
		operator unsigned int () const;
		operator float () const;
		operator double () const;
		operator const char *() const;
		operator float * () const;
		operator int * () const;
		operator void * () const;
		operator EMData *() const;
		operator XYData *() const;
		operator Transform *() const;
		operator Ctf *() const;
		operator vector < int > () const;
		operator vector < float > () const;
		operator vector<string> () const;

		/** Checks to see if the EMObject is interpretable
		 * This basically equates to checking to see if the
		 * type is UNKNOWN
		 */
		bool is_null() const;

		/** Calls to_str( this->type)
		 */
		string to_str() const;

		/** Get the ObjectType
		 * This is an enumerated type first declared in the class
		 * EMObjectTypes
		 */
		ObjectType get_type() const { return type; }

		/** Get the ObjectType as a string
		 * This is an enumerated type first declared in the class
		 * EMObjectTypes
		 */
		string get_type_string() const { return get_object_type_name(type); }


		/** Write the EMObject's value to a string
		 * Literally copies into a string, except for the case
		 * where the type is boolen where it writes true of false,
		 * as opposed to 1 or 0.
		 */
		string to_str(ObjectType type) const;

		/** Get an ObjectType as a string statically
		 * Can be accessed without the instantiation of a class object
		 */
		static string get_object_type_name(ObjectType t);

		/** Friend declaration operator==
		 * namespace EMAN2 operator== accesses private variables
		 */
		friend bool operator==(const EMObject &e1, const EMObject & e2);

		/** Friend declaration operator!=
		 * namespace EMAN2 operator!= accesses private variables
		 */
		friend bool operator!=(const EMObject &e1, const EMObject & e2);

	private:
		union
		{
			bool b;
			int n;
			unsigned int ui;
			float f;
			double d;
			float * fp;
			int * ip;
			void * vp;
			EMData *emdata;
			XYData *xydata;
		};

		string str;
		vector < int > iarray;
		vector < float >farray;
		vector < string> strarray;
		ObjectType type;

		/** A debug function that prints as much information as possibe to cout
		 */
		void printInfo() const;

//		void init();
		
		static map< ObjectType, string> init();
		static map< ObjectType, string> type_registry;
		
	};

	bool operator==(const EMObject &e1, const EMObject & e2);
	bool operator!=(const EMObject &e1, const EMObject & e2);

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

			inline bool find_type( const string& type ) {  if ( type_dict.find(type) != type_dict.end() ) return true; return false; }

		private:
			map < string, string > type_dict;
			map < string, string > desc_dict;
	};


	/** Dict is a dictionary to store <string, EMObject> pair.
     * Typical ways to construct a Dict:
     *
     *      Dict d;
     *      d["lowpass"] = 12.23;
     *      float lowpass1 = d["lowpass"];
     *
     *      Dict d2("lowpass", 12.23);
     *
     * You can iterate through a dict:
     *	for ( Dict::const_iterator it = params.begin(); it != params.end(); ++it ) { //do things to it }
     * And similary use the Dict iterator as arguments to the generic algorithms that are feasible, such as copy.
     *
     * You can find things in the iterator style:
     * 	if(	d.find("lowpass") != d.end() ) cout << "D has a lowpass key" << endl;\
     * Or like this
     * 	if( d.has_key("lowpass") ) ...
     *
     * A Dict has copy and assignment operators.
     *
     *
     * See the testing code in rt/emdata/test_emobject.cpp for prewritten testing code
     */
	class Dict
	{
	public:
		Dict()
		{
		}

		/** Construct a Dict object from 1 key/value pair
		 * It's probably more conventional to intialize key/value pairs
		 * using operator[], but either approach is fine.
		 */
		Dict(const string & key1, EMObject val1)
		{
			dict[key1] = val1;
		}

		/** Construct a Dict object from 2 key/value pairs
		 */
		Dict(const string & key1, EMObject val1,
			 const string & key2, EMObject val2)
		{
			dict[key1] = val1;
			dict[key2] = val2;
		}

		/** Construct a Dict object from 3 key/value pairs
		 */
		Dict(const string & key1, EMObject val1,
			 const string & key2, EMObject val2,
			 const string & key3, EMObject val3)
		{
			dict[key1] = val1;
			dict[key2] = val2;
			dict[key3] = val3;
		}

		/** Construct a Dict object from 4 key/value pairs
		 */
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

		/** Construct a Dict object from a map object
		 * Calls the generic algorithm "copy".
		 */
		Dict(const map < string, EMObject > &d)
		{
			copy(d.begin(), d.end(), inserter(dict, dict.begin()));
			// Or use
			// dict.insert(d.begin(), d.end());
		}

		/** Destructor
		 * Performs no explicit action besides what the compiler automatically does.
		 */
		~Dict() {}

		/** Copy constructor
		 * Copies all elements in dict
		 */
		Dict( const Dict& that);

		/** Assignment operator
		 * Copies all elements in dict
		 */
		Dict& operator=(const Dict& that);

		/**	Get a vector containing all of the (string) keys in this dictionary.
		 */
		vector < string > keys()const
		{
			vector < string > result;

			map < string, EMObject >::const_iterator p;
			for (p = dict.begin(); p != dict.end(); p++) {
				result.push_back(p->first);
			}

			return result;
		}

		/** Get a vector containing copies of each of the EMObjects in this dictionary.
		 */
		vector < EMObject > values()const
		{
			vector < EMObject > result;

			map < string, EMObject >::const_iterator p;
			for (p = dict.begin(); p != dict.end(); p++) {
				result.push_back(p->second);
			}

			return result;
		}

		/** Ask the Dictionary if it as a particular key in a case insensitive way
		 * @param key the (string) key to find
		 */
		bool has_key_ci(const string & key) const;

		/** Ask the Dictionary if it as a particular key
		 * @param key the (string) key to find
		 */
		bool has_key(const string & key) const
		{
			map < string, EMObject >::const_iterator p = dict.find(key);
			if (p != dict.end()) {
				return true;
			}
			return false;
		}

		/** Ask the Dictionary for its size
		 */
		size_t size() const
		{
			return dict.size();
		}

		/** Get the EMObject corresponding to the particular key
		 * Probably better to just use operator[]
		 */
		EMObject get(const string & key) const
		{
			if( has_key(key) ) {
				return dict[key];
			}
			else {
				LOGERR("No such key exist in this Dict");
				throw NotExistingObjectException("EMObject", "Nonexisting key (" + key + ") in Dict");
			}
		}

		/** Get the EMObject corresponding to the particular key using case insensitivity
		 * @param key the key you want to check for in a case insensitive way
		 */
		EMObject get_ci(const string & key) const;
		/** Put the value/key pair into the dictionary
		 * probably better to just use operator[]
		 */
		void put(const string & key, EMObject val)
		{
			dict[key] = val;
		}

		/** Remove a particular key
		 */
		void erase(const string & key)
		{
			dict.erase(key);
		}

		/** Clear all keys
		 * wraps map.clear()
		 */
		void clear()
		{
			dict.clear();
		}

		/** Default setting behavior
		 * This can be achieved using a template - d.woolford Jan 2008 (before there was a function being written for every type)
		 */
		template<typename type>
		type set_default(const string & key, type val)
		{
			if (!has_key(key)) {
				dict[key] = val;
			}
			return dict[key];
		}

		Dict copy_exclude_keys(const vector<string>& excluded_keys) const
		{
			Dict ret(*this);

			for ( vector<string>::const_iterator it = excluded_keys.begin(); it != excluded_keys.end(); ++it ) {
				if (ret.has_key(*it)) ret.erase(*it);
			}

			return ret;
		}

		Dict copy_exclusive_keys(const vector<string>& exclusive_keys) const
		{
			Dict ret;
			for ( vector<string>::const_iterator it = exclusive_keys.begin(); it != exclusive_keys.end(); ++it ) {
				if (has_key(*it)) ret[*it] = (*this)[*it];
			}

			return ret;
		}

		Dict copy_keys_in( const TypeDict& tdict ) const {
			vector<string> keys = tdict.keys();
			return copy_exclusive_keys(keys);
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

		/** Friend declaration operator==
		 * namespace EMAN2 operator== accesses private variables
		 */
		friend bool operator==(const Dict& d1, const Dict& d2);

		/** Friend declaration operator!=
		 * namespace EMAN2 operator!= accesses private variables
		 */
		friend bool operator!=(const Dict& d1, const Dict& d2);

	private:
		mutable map < string, EMObject > dict;

	public:
		/** Non const iterator support for the Dict object
		* This is just a wrapper, everything is inherited from the map<string,EMObject>::iterator
		* so the interface is the same as you would expect
		* i.e for ( Dict::iterator it = params.begin(); it != params.end(); ++it )
		* @author David Woolford
		* @date Mid 2007
		*/
		class iterator : public map < string, EMObject >::iterator
		{
		public:
			typedef std::bidirectional_iterator_tag iterator_category;
 			typedef pair<string, EMObject> value_type;

		public:
			iterator( map < string, EMObject >::iterator parent_it );
			virtual ~iterator(){}

			iterator( const iterator& that );
			iterator& operator=( const iterator& that );
		};

		/** Const iterator support for the Dict object
		 * This is just a wrapper, everything is inherited from the map<string,EMObject>::cons_iterator
		 * so the interface is the same as you would expect
		 * i.e for ( Dict::const_iterator it = params.begin(); it != params.end(); ++it )
		 * @author David Woolford
		 * @date Mid 2007
		 */
		class const_iterator :  public map < string, EMObject >::const_iterator
		{
		public:
			typedef std::bidirectional_iterator_tag iterator_category;
			typedef pair<string, EMObject> value_type; // Note that value_type should NOT be const even though the container elements are const
		public:
			const_iterator( const map < string, EMObject >::const_iterator parent_it);
			virtual ~const_iterator(){}
			const_iterator( const Dict::iterator& it );

			const_iterator( const const_iterator& that );
			const_iterator& operator=( const const_iterator& that );
		};

		// Iterator support
		iterator begin( void );
		const_iterator begin( void ) const;

		iterator end( void );
		const_iterator end( void ) const;

		// Wraps map.find(const string& key)
		iterator find( const string& key );
		const_iterator find( const string& key ) const;
	};

	// These operators were originally added for the purposes of making testing code but might come in handy for other things
	// operator== simply wraps map<string, EMObject>::operator==
	bool operator==(const Dict &d1, const Dict& d2);
	bool operator!=(const Dict &d1, const Dict& d2);


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
     *	    Processor *f1 = Factory<Processor>::get("math.absvalue");
     *      Processor *f2 = Factory<Processor>::get("filter.lowpass.gauss", Dict("cufoff_freq", EMObject(12));
	 * @author Liwei Peng
     */
	template < class T > class Factory
	{
	public:
		typedef T *(*InstanceType) ();

		template <class ClassType> static void add();
		static T *get(const string & instance_name);
		static T *get(const string & instance_name, const Dict & params);
		static vector < string > get_list();

	private:
		Factory();
		Factory(const Factory < T > &);
		~Factory();
		static void init();
		template <class ClassType> void force_add();

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

	template < class T > 
	template < class ClassType > 
	void Factory < T >::force_add()
	{
		string name = ClassType::NAME;
		my_dict[name] = &ClassType::NEW;
	}


	template < class T > 
	template < class ClassType >
	void Factory < T >::add()
	{
		init();

		string name = ClassType::NAME;
		typename map < string, InstanceType >::iterator fi =
			my_instance->my_dict.find(name);

		if (fi == my_instance->my_dict.end()) {
			my_instance->my_dict[name] = &ClassType::NEW;
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

		string lower = instancename;
		for (unsigned int i=0; i<lower.length(); i++) lower[i]=tolower(lower[i]);

		fi = my_instance->my_dict.find(lower);
		if (fi != my_instance->my_dict.end()) {
			return my_instance->my_dict[lower] ();
		}

		throw NotExistingObjectException(instancename, "The named object doesn't exist");
	}

	template < class T > T * Factory < T >::get(const string & instancename,
												const Dict & params)
	{
		init();

		typename map < string, InstanceType >::iterator fi =
			my_instance->my_dict.find(instancename);

		string lower = instancename;
		if (fi == my_instance->my_dict.end()) {
			for (unsigned int i=0; i<lower.length(); i++) lower[i]=tolower(lower[i]);
			fi = my_instance->my_dict.find(lower);
		}

		if (fi != my_instance->my_dict.end()) {
			T *i = my_instance->my_dict[lower] ();

			const vector<string> para_keys = params.keys();
//			std::cout << "the number of keys is " << para_keys.size() << std::endl; // PRB May 19th
			const vector<string> valid_keys = i->get_param_types().keys();
			typename vector<string>::const_iterator it;
			for(it=para_keys.begin(); it!=para_keys.end(); ++it) {
// 				std::cout << "the iterator  is " << *it << std::endl; // PRB May 19th
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

	/** A class one may inherit from to ensure that the responsibilities of
	* being incorporated into an EMAN2::Factory are met.
	* This class is abstract.
	* @author David Woolford
	* @date Feb 2008
	*/
	class FactoryBase
	{
	public:
		FactoryBase() {}
		virtual ~FactoryBase() {};

		/** Get the unique name of this class (especially for factory based instantiation access)
		 * @return the unique name of this class
		*/
		virtual string get_name() const = 0;

		/** Get a clear, concise description of this class
		 * @return a clear, concise description of this class
		 */
		virtual string get_desc() const = 0;

		/** get a copy of the parameters of this class
		 * @return a copy of the parameters of this class
		 */
		Dict get_params() const	{ return params; }

		/** Set new parameters. Old parameters are cleared
		 * @param new_params the new parameters
		 */
		void set_params(const Dict & new_params)
		{
			params.clear();
			insert_params(new_params);
		}
		
		inline void set_param(const string key,const EMObject val) { params[key]=val; }
		
		/** @return a TypeDict defining and describing the feasible parameters of this class 
		 */
		virtual TypeDict get_param_types() const = 0;

		/** Insert parameters. Previously present parameters are replaced, new ones are inserted.
		 * @param new_params the parameters to insert
		 */
		void insert_params(const Dict & new_params)
		{
		// this is really inserting OR individually replacing...
		// the old data will be kept if it is not written over
			TypeDict permissable_params = get_param_types();
			for ( Dict::const_iterator it = new_params.begin(); it != new_params.end(); ++it )
			{

				if ( !permissable_params.find_type(it->first) )
				{
					throw InvalidParameterException(it->first);
				}
				params[it->first] = it->second;
			}
		}

		Dict copy_relevant_params(const FactoryBase* const that) const
		{
			return params.copy_keys_in(that->get_param_types());

		}

		protected:
		/// This is the dictionary the stores the parameters of the object
		mutable Dict params;
	};
}

#endif
