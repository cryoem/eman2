
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <emdata.h>
#include <emobject.h>
#include <filter.h>
#include <floatstat.h>
#include <log.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct EMAN_Filter_Wrapper: EMAN::Filter
{
    EMAN_Filter_Wrapper(PyObject* self_, const EMAN::Filter& p0):
        EMAN::Filter(p0), self(self_) {}

    EMAN_Filter_Wrapper(PyObject* self_):
        EMAN::Filter(), self(self_) {}

    void process(EMAN::EMData* p0) {
        call_method< void >(self, "process", p0);
    }

    void default_process(EMAN::EMData* p0) {
        EMAN::Filter::process(p0);
    }

    void process_list(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >& p0) {
        call_method< void >(self, "process_list", p0);
    }

    void default_process_list(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >& p0) {
        EMAN::Filter::process_list(p0);
    }

    EMAN::Dict get_params() const {
        return call_method< EMAN::Dict >(self, "get_params");
    }

    EMAN::Dict default_get_params() const {
        return EMAN::Filter::get_params();
    }

    void set_params(const EMAN::Dict& p0) {
        call_method< void >(self, "set_params", p0);
    }

    void default_set_params(const EMAN::Dict& p0) {
        EMAN::Filter::set_params(p0);
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(self, "get_param_types");
    }

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
    }

    PyObject* self;
};


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyFilter)
{
    class_< EMAN::Filter, boost::noncopyable, EMAN_Filter_Wrapper >("Filter", init<  >())
        .def("process", &EMAN::Filter::process, &EMAN_Filter_Wrapper::default_process)
        .def("process_list", &EMAN::Filter::process_list, &EMAN_Filter_Wrapper::default_process_list)
        .def("get_params", &EMAN::Filter::get_params, &EMAN_Filter_Wrapper::default_get_params)
        .def("set_params", &EMAN::Filter::set_params, &EMAN_Filter_Wrapper::default_set_params)
        .def("get_param_types", pure_virtual(&EMAN::Filter::get_param_types))
        .def("get_name", pure_virtual(&EMAN::Filter::get_name))
    ;

    class_< EMAN::Factory<EMAN::Filter>, boost::noncopyable >("FilterFactory", no_init)
        .def("instance", &EMAN::Factory<EMAN::Filter>::instance, return_value_policy< reference_existing_object >())
        .def("get", (EMAN::Filter* (EMAN::Factory<EMAN::Filter>::*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >) )&EMAN::Factory<EMAN::Filter>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Filter* (EMAN::Factory<EMAN::Filter>::*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&) )&EMAN::Factory<EMAN::Filter>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Filter>::get_list)
        .staticmethod("instance")
    ;

}

