
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <averager.h>
#include <emdata.h>
#include <emobject.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
BOOST_PYTHON_MODULE(libpyAverager2)
{
    class_< EMAN::Factory<EMAN::Averager>, boost::noncopyable >("Averagers", no_init)
        .def("get", (EMAN::Averager* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >))&EMAN::Factory<EMAN::Averager>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Averager* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&))&EMAN::Factory<EMAN::Averager>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Averager>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

    def("dump_averagers", &EMAN::dump_averagers);
}

