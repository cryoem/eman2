
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <emdata.h>
#include <emobject.h>
#include <projector.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
BOOST_PYTHON_MODULE(libpyProjector2)
{
    def("dump_projectors", &EMAN::dump_projectors);
    class_< EMAN::Factory<EMAN::Projector>, boost::noncopyable >("Projectors", no_init)
        .def("get", (EMAN::Projector* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >))&EMAN::Factory<EMAN::Projector>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Projector* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&))&EMAN::Factory<EMAN::Projector>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Projector>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

}

