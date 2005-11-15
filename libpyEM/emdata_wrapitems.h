#include "emdata.h"
#include <boost/python.hpp>

#ifndef EMDATA_WRAPPER_H_
#define EMDATA_WRAPPER_H_

using boost::python::object;

object emdata_getitem(object self, object key);
void emdata_setitem(object self, object key, object val);

#endif // EMDATA_WRAPPER_H_
