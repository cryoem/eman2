#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/exception_translator.hpp>

#include "pyexception.h"
#include "exception.h"

//using namespace EMAN;

void translate1(const EMAN::_ImageFormatException& e)
{
	PyErr_SetString(PyExc_RuntimeError, "liwei peng");
}

