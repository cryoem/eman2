
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <cmp.h>
#include <emdata.h>
#include <emobject.h>
#include <log.h>
#include <transform.h>
#include <xydata.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_Log_end_overloads_1_3, end, 1, 3)


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyCmp2)
{
    def("dump_cmps", &EMAN::dump_cmps);
    class_< EMAN::Factory<EMAN::Cmp>, boost::noncopyable >("Cmps", no_init)
        .def("get", (EMAN::Cmp* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >))&EMAN::Factory<EMAN::Cmp>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Cmp* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&))&EMAN::Factory<EMAN::Cmp>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Cmp>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

    scope* EMAN_Log_scope = new scope(
    class_< EMAN::Log, boost::noncopyable >("Log", no_init)
        .def("logger", &EMAN::Log::logger, return_value_policy< reference_existing_object >())
        .def("begin", &EMAN::Log::begin)
        .def("end", &EMAN::Log::end, EMAN_Log_end_overloads_1_3())
        .def("set_level", &EMAN::Log::set_level)
        .def("set_logfile", &EMAN::Log::set_logfile)
        .def("loc", &EMAN::Log::loc)
        .staticmethod("logger")
    );

    enum_< EMAN::Log::LogLevel >("LogLevel")
        .value("ERROR_LOG", EMAN::Log::ERROR_LOG)
        .value("VARIABLE_LOG", EMAN::Log::VARIABLE_LOG)
        .value("WARNING_LOG", EMAN::Log::WARNING_LOG)
        .value("DEBUG_LOG", EMAN::Log::DEBUG_LOG)
    ;

    delete EMAN_Log_scope;

    scope* EMAN_XYData_scope = new scope(
    class_< EMAN::XYData >("XYData", init<  >())
        .def(init< const EMAN::XYData& >())
        .def("read_file", &EMAN::XYData::read_file)
        .def("write_file", &EMAN::XYData::write_file)
        .def("calc_correlation", &EMAN::XYData::calc_correlation)
        .def("update", &EMAN::XYData::update)
        .def("get_yatx", &EMAN::XYData::get_yatx)
        .def("get_x", &EMAN::XYData::get_x)
        .def("set_x", &EMAN::XYData::set_x)
        .def("get_y", &EMAN::XYData::get_y)
        .def("set_y", &EMAN::XYData::set_y)
        .def("get_size", &EMAN::XYData::get_size)
        .def("get_miny", &EMAN::XYData::get_miny)
        .def("get_maxy", &EMAN::XYData::get_maxy)
        .def("is_validx", &EMAN::XYData::is_validx)
    );

    class_< EMAN::XYData::Pair >("Pair", init< const EMAN::XYData::Pair& >())
        .def(init< float, float >())
        .def_readwrite("x", &EMAN::XYData::Pair::x)
        .def_readwrite("y", &EMAN::XYData::Pair::y)
        .def( self < self )
    ;

    delete EMAN_XYData_scope;

}

