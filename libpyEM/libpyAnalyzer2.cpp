
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <analyzer.h>
#include <emdata.h>
#include <emobject.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct EMAN_Analyzer_Wrapper: EMAN::Analyzer
{
    EMAN_Analyzer_Wrapper(PyObject* py_self_, const EMAN::Analyzer& p0):
        EMAN::Analyzer(p0), py_self(py_self_) {}

    EMAN_Analyzer_Wrapper(PyObject* py_self_):
        EMAN::Analyzer(), py_self(py_self_) {}

    std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> > analyze(const std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >& p0) {
        return call_method< std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> > >(py_self, "analyze", p0);
    }

    std::string get_name() const {
        return call_method< std::string >(py_self, "get_name");
    }

    std::string get_desc() const {
        return call_method< std::string >(py_self, "get_desc");
    }

    void set_params(const EMAN::Dict& p0) {
        call_method< void >(py_self, "set_params", p0);
    }

    void default_set_params(const EMAN::Dict& p0) {
        EMAN::Analyzer::set_params(p0);
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(py_self, "get_param_types");
    }

    EMAN::TypeDict default_get_param_types() const {
        return EMAN::Analyzer::get_param_types();
    }

    PyObject* py_self;
};


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyAnalyzer2)
{
    class_< EMAN::Analyzer, boost::noncopyable, EMAN_Analyzer_Wrapper >("__Analyzer", init<  >())
        .def("analyze", pure_virtual(&EMAN::Analyzer::analyze))
        .def("get_name", pure_virtual(&EMAN::Analyzer::get_name))
        .def("get_desc", pure_virtual(&EMAN::Analyzer::get_desc))
        .def("set_params", &EMAN::Analyzer::set_params, &EMAN_Analyzer_Wrapper::default_set_params)
        .def("get_param_types", &EMAN::Analyzer::get_param_types, &EMAN_Analyzer_Wrapper::default_get_param_types)
    ;

    class_< EMAN::Factory<EMAN::Analyzer>, boost::noncopyable >("Analyzers", no_init)
        .def("get", (EMAN::Analyzer* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&))&EMAN::Factory<EMAN::Analyzer>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Analyzer* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&, const EMAN::Dict&))&EMAN::Factory<EMAN::Analyzer>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Analyzer>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

}

