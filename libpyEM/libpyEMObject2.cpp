
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <emdata.h>
#include <emobject.h>
#include <xydata.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_TypeDict_put_overloads_2_3, put, 2, 3)


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyEMObject2)
{
    class_< EMAN::TypeDict >("TypeDict", init<  >())
        .def(init< const EMAN::TypeDict& >())
        .def("keys", &EMAN::TypeDict::keys)
        .def("size", &EMAN::TypeDict::size)
        .def("put", &EMAN::TypeDict::put, EMAN_TypeDict_put_overloads_2_3())
        .def("get_type", &EMAN::TypeDict::get_type)
        .def("get_desc", &EMAN::TypeDict::get_desc)
        .def("dump", &EMAN::TypeDict::dump)
    ;

}

