
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <emdata.h>
#include <emobject.h>
#include <reconstructor.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
BOOST_PYTHON_MODULE(libpyReconstructor2)
{
    def("dump_reconstructors", &EMAN::dump_reconstructors);
    class_< EMAN::Factory<EMAN::Reconstructor>, boost::noncopyable >("Reconstructors", no_init)
        .def("get", (EMAN::Reconstructor* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >))&EMAN::Factory<EMAN::Reconstructor>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Reconstructor* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&))&EMAN::Factory<EMAN::Reconstructor>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Reconstructor>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

}

