
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

struct EMAN_Cmp_Wrapper: EMAN::Cmp
{
    EMAN_Cmp_Wrapper(PyObject* self_, const EMAN::Cmp& p0):
        EMAN::Cmp(p0), self(self_) {}

    EMAN_Cmp_Wrapper(PyObject* self_):
        EMAN::Cmp(), self(self_) {}

    float cmp(EMAN::EMData* p0, EMAN::EMData* p1) const {
        return call_method< float >(self, "cmp", p0, p1);
    }

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
    }

    std::string get_desc() const {
        return call_method< std::string >(self, "get_desc");
    }

    EMAN::Dict get_params() const {
        return call_method< EMAN::Dict >(self, "get_params");
    }

    EMAN::Dict default_get_params() const {
        return EMAN::Cmp::get_params();
    }

    void set_params(const EMAN::Dict& p0) {
        call_method< void >(self, "set_params", p0);
    }

    void default_set_params(const EMAN::Dict& p0) {
        EMAN::Cmp::set_params(p0);
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(self, "get_param_types");
    }

    PyObject* self;
};

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_Log_end_overloads_1_3, end, 1, 3)


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyCmp2)
{
    def("dump_cmps", &EMAN::dump_cmps);
    class_< EMAN::Factory<EMAN::Cmp>, boost::noncopyable >("Cmps", no_init)
        .def("get", (EMAN::Cmp* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&))&EMAN::Factory<EMAN::Cmp>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Cmp* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&, const EMAN::Dict&))&EMAN::Factory<EMAN::Cmp>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Cmp>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

    class_< EMAN::Cmp, boost::noncopyable, EMAN_Cmp_Wrapper >("__Cmp", init<  >())
        .def("cmp", pure_virtual(&EMAN::Cmp::cmp))
        .def("get_name", pure_virtual(&EMAN::Cmp::get_name))
        .def("get_desc", pure_virtual(&EMAN::Cmp::get_desc))
        .def("get_params", &EMAN::Cmp::get_params, &EMAN_Cmp_Wrapper::default_get_params)
        .def("set_params", &EMAN::Cmp::set_params, &EMAN_Cmp_Wrapper::default_set_params)
        .def("get_param_types", pure_virtual(&EMAN::Cmp::get_param_types))
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

