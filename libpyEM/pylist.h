#ifndef __pylist_h__
#define __pylist_h__

#include <boost/python.hpp>
#include <vector>

namespace python = boost::python;

namespace EMAN {
    class PyList {
    public:
	static void list2array(const python::list& l, int* array);
	static void list2array(const python::list& l, float* array);
	static void list2array(const python::list& l, const char** array);
	static void array2list(const int* array, python::list& l, int nitems = 0);
	static void array2list(const float* array, python::list& l, int nitems = 0);
	
	template <class T>
	static python::list vector2list(const std::vector<T> & v)
	{
	    python::list res;
	    for (unsigned int i = 0; i < v.size(); i++) {
		res.append(v[i]);
	    }
	    
	    Py_XINCREF(res.ptr());
	    return res;
	}
    };
}


#endif
