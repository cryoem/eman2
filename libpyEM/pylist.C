#include "pylist.h"
#include <boost/python/detail/api_placeholder.hpp>

using namespace EMAN;

void PyList::list2array(const python::list& l, float* array)
{
    if (!array) {
	return;
    }
    for (int i = 0; i < python::len(l); i++) {
	array[i] = python::extract<float>(l[i]);
	
    }
}

void PyList::list2array(const python::list& l, int* array)
{
    if (!array) {
	return;
    }
    for (int i = 0; i < python::len(l); i++) {
	array[i] = python::extract<int>(l[i]);
	
    }
}


void PyList::list2array(const python::list& l, const char** array)
{
    if (!array) {
	return;
    }
    for (int i = 0; i < python::len(l); i++) {
	array[i] = python::extract<const char*>(l[i]);
	
    }
}


void PyList::array2list(const float* array, python::list& l, int nitems)
{
    if (!array) {
	return;
    }

    if (nitems != 0) {
	for (int i = 0; i < nitems; i++) {
	    l.append(array[i]);
	}
    }
    else {    
	for (int i = 0; i < python::len(l); i++) {
	    l[i] = array[i];
	}
    }
}

void PyList::array2list(const int* array, python::list& l, int nitems)
{
    if (!array) {
	return;
    }

    if (nitems != 0) {
	for (int i = 0; i < nitems; i++) {
	    l.append(array[i]);
	}
    }
    else {    
	for (int i = 0; i < python::len(l); i++) {
	    l[i] = array[i];
	}
    }
}
