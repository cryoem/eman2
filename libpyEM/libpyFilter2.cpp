
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <emdata.h>
#include <emobject.h>
#include <filter.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct EMAN_Processor_Wrapper: EMAN::Processor
{
    void process(EMAN::EMData* p0) {
        call_method< void >(py_self, "process", p0);
    }

    void default_process(EMAN::EMData* p0) {
        EMAN::Processor::process(p0);
    }

    void process_list(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >& p0) {
        call_method< void >(py_self, "process_list", p0);
    }

    void default_process_list(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >& p0) {
        EMAN::Processor::process_list(p0);
    }

    std::string get_name() const {
        return call_method< std::string >(py_self, "get_name");
    }

    EMAN::Dict get_params() const {
        return call_method< EMAN::Dict >(py_self, "get_params");
    }

    EMAN::Dict default_get_params() const {
        return EMAN::Processor::get_params();
    }

    void set_params(const EMAN::Dict& p0) {
        call_method< void >(py_self, "set_params", p0);
    }

    void default_set_params(const EMAN::Dict& p0) {
        EMAN::Processor::set_params(p0);
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(py_self, "get_param_types");
    }

    EMAN::TypeDict default_get_param_types() const {
        return EMAN::Processor::get_param_types();
    }

    std::string get_desc() const {
        return call_method< std::string >(py_self, "get_desc");
    }

    PyObject* py_self;
};


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyFilter2)
{
    scope* EMAN_Processor_scope = new scope(
    class_< EMAN::Processor, boost::noncopyable, EMAN_Processor_Wrapper >("__Processor", no_init)
        .def("process", &EMAN::Processor::process, &EMAN_Processor_Wrapper::default_process)
        .def("process_list", &EMAN::Processor::process_list, &EMAN_Processor_Wrapper::default_process_list)
        .def("get_name", pure_virtual(&EMAN::Processor::get_name))
        .def("get_params", &EMAN::Processor::get_params, &EMAN_Processor_Wrapper::default_get_params)
        .def("set_params", &EMAN::Processor::set_params, &EMAN_Processor_Wrapper::default_set_params)
        .def("get_param_types", &EMAN::Processor::get_param_types, &EMAN_Processor_Wrapper::default_get_param_types)
        .def("get_desc", pure_virtual(&EMAN::Processor::get_desc))
        .def("get_group_desc", &EMAN::Processor::get_group_desc)
        .def("EMFourierFilterInPlace", &EMAN::Processor::EMFourierFilterInPlace)
        .def("EMFourierFilter", &EMAN::Processor::EMFourierFilter, return_value_policy< manage_new_object >())
        .staticmethod("EMFourierFilterInPlace")
        .staticmethod("get_group_desc")
        .staticmethod("EMFourierFilter")
    );

    enum_< EMAN::Processor::fourier_filter_types >("fourier_filter_types")
        .value("GAUSS_HIGH_PASS", EMAN::Processor::GAUSS_HIGH_PASS)
        .value("TANH_LOW_PASS", EMAN::Processor::TANH_LOW_PASS)
        .value("GAUSS_BAND_PASS", EMAN::Processor::GAUSS_BAND_PASS)
        .value("BUTTERWORTH_LOW_PASS", EMAN::Processor::BUTTERWORTH_LOW_PASS)
        .value("TOP_HAT_LOW_PASS", EMAN::Processor::TOP_HAT_LOW_PASS)
        .value("GAUSS_HOMOMORPHIC", EMAN::Processor::GAUSS_HOMOMORPHIC)
        .value("GAUSS_LOW_PASS", EMAN::Processor::GAUSS_LOW_PASS)
        .value("TOP_HAT_BAND_PASS", EMAN::Processor::TOP_HAT_BAND_PASS)
        .value("BUTTERWORTH_HOMOMORPHIC", EMAN::Processor::BUTTERWORTH_HOMOMORPHIC)
        .value("TOP_HOMOMORPHIC", EMAN::Processor::TOP_HOMOMORPHIC)
        .value("TANH_HOMOMORPHIC", EMAN::Processor::TANH_HOMOMORPHIC)
        .value("TOP_HAT_HIGH_PASS", EMAN::Processor::TOP_HAT_HIGH_PASS)
        .value("TANH_HIGH_PASS", EMAN::Processor::TANH_HIGH_PASS)
        .value("BUTTERWORTH_HIGH_PASS", EMAN::Processor::BUTTERWORTH_HIGH_PASS)
        .value("TANH_BAND_PASS", EMAN::Processor::TANH_BAND_PASS)
    ;

    delete EMAN_Processor_scope;

    def("dump_filters", &EMAN::dump_filters);
    def("multi_filters", &EMAN::multi_filters);
    def("group_filters", &EMAN::group_filters);
    class_< EMAN::Factory<EMAN::Processor>, boost::noncopyable >("Processors", no_init)
        .def("get", (EMAN::Processor* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&))&EMAN::Factory<EMAN::Processor>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Processor* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&, const EMAN::Dict&))&EMAN::Factory<EMAN::Processor>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Processor>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

}

