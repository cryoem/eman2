#ifndef eman__pylist_h__
#define eman__pylist_h__

#include <boost/python.hpp>
#include <boost/python/to_python_converter.hpp>
#include <boost/python/detail/api_placeholder.hpp>

#include <vector>
#include <map>
#include <string>

namespace python = boost::python;
using std::vector;
using std::map;

namespace EMAN {
    class PyList {
    public:

		template <class T>
		static vector<T> list2vector(const python::list& l)
		{
			vector<T> v;
			for (int i = 0; i < python::len(l); i++) {
				v.push_back(python::extract<T>(l[i]));
		
			}
			return v;
		}

		template <class T>
		static python::list vector2list(const vector<T>& v)
		{
			python::list l;
			for (unsigned int i = 0; i < v.size(); i++) {
				l.append(v[i]);
			}
			return l;
		}
    };
    
    template <class T>
    struct vector_to_python : python::to_python_converter<vector<T>, vector_to_python<T> >
    {
		static PyObject* convert(vector<T> const& v)
		{
			python::list result;
	    
			for (size_t i = 0; i < v.size(); i++) {
				result.append(v[i]);
			}
	    
			return python::incref(python::list(result).ptr());
		}
    };

    template <class T>
    struct map_to_python : python::to_python_converter<map<std::string, T>, map_to_python<T> >
    {
		static PyObject* convert(map<std::string, T> const& d)
		{
			python::dict result;

			typedef typename map<std::string, T>::const_iterator MI;
			for (MI p = d.begin(); p != d.end(); p++) {
				result[p->first] = p->second;
			}
	
			return python::incref(python::dict(result).ptr());
		}
    };

    struct Dict_to_python : python::to_python_converter<Dict, Dict_to_python>
    {
		static PyObject* convert(Dict const& dd)
		{
			python::dict result;
			vector<std::string> keys = dd.keys();
			vector<EMObject> values = dd.values();
			for (unsigned int i = 0; i < keys.size(); i++) {
				result[keys[i]] = values[i];
			}
	
			return python::incref(python::dict(result).ptr());
		}
    };
    
    template <class T>
    struct vector_from_python
    {
		vector_from_python()
		{
			python::converter::registry::push_back(&convertible, &construct,
												   python::type_id<vector<T> >());
		}
    
		static void* convertible(PyObject* obj_ptr)
		{
			if (!(PyList_Check(obj_ptr) || PyTuple_Check(obj_ptr)
				  || PyIter_Check(obj_ptr)  || PyRange_Check(obj_ptr))) {
				return 0;
			}
	
			return obj_ptr;
		}

    
		static void construct(PyObject* obj_ptr,
							  python::converter::rvalue_from_python_stage1_data* data)
		{
			void* storage = ((python::converter::rvalue_from_python_storage<vector<T> >*) data)->storage.bytes;
			new (storage) vector<T>();

			data->convertible = storage;

			vector<T>& result = *((vector<T>*) storage);
	
			python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
	
			while(1) {
				python::handle<> py_elem_hdl(python::allow_null(PyIter_Next(obj_iter.get())));
				if (PyErr_Occurred()) {
					python::throw_error_already_set();
				}
	    
				if (!py_elem_hdl.get()) {
					break;
				}
	    
				python::object py_elem_obj(py_elem_hdl);
				python::extract<T> elem_proxy(py_elem_obj);
				result.push_back(elem_proxy());
			}
		}
    };
    
    template <class T>
    struct map_from_python
    {
		map_from_python()
		{
			python::converter::registry::push_back(&convertible, &construct,
												   python::type_id<map<std::string, T> >());
		}
    
		static void* convertible(PyObject* obj_ptr)
		{
			if (!(PyDict_Check(obj_ptr))) {
				return 0;
			}
	
			return obj_ptr;
		}

    
		static void construct(PyObject* obj_ptr,
							  python::converter::rvalue_from_python_stage1_data* data)
		{
			void* storage = ((python::converter::rvalue_from_python_storage<map<std::string, T> >*) data)->storage.bytes;
			new (storage) map<std::string, T>();
			data->convertible = storage;
			map<std::string, T>& result = *((map<std::string, T>*) storage);

			python::handle<> obj_handle(obj_ptr);
			python::object dict_obj(obj_handle);
	    
			python::dict d(dict_obj);
		       
			python::list k = d.keys();
			python::list v = d.values();
			long l = python::len(k);
	
			for(long i = 0; i < l; i++) {
				std::string key = python::extract<std::string>(k[i]);
				T val = python::extract<T>(v[i]);
				result[key] = val;
			}

		}
    };

    struct Dict_from_python
    {
		Dict_from_python()
		{
			python::converter::registry::push_back(&convertible, &construct,
												   python::type_id<Dict>());
		}
    
		static void* convertible(PyObject* obj_ptr)
		{
			if (!(PyDict_Check(obj_ptr))) {
				return 0;
			}
	
			return obj_ptr;
		}

    
		static void construct(PyObject* obj_ptr,
							  python::converter::rvalue_from_python_stage1_data* data)
		{
			void* storage = ((python::converter::rvalue_from_python_storage<Dict>*) data)->storage.bytes;
			new (storage) Dict();
			data->convertible = storage;
			Dict& result = *((Dict*) storage);

			python::handle<> obj_handle(obj_ptr);
			python::object dict_obj(obj_handle);
	    
			python::dict d(dict_obj);
	    
			python::list k = d.keys();
			python::list v = d.values();
			long l = python::len(k);
	
			for(long i = 0; i < l; i++) {
				std::string key = python::extract<std::string>(k[i]);
				EMObject val = python::extract<EMObject>(v[i]);
				result.put(key, val);
			}

		}
    };
    
}


#endif
