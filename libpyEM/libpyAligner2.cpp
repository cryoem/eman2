
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <aligner.h>
#include <emdata.h>
#include <emobject.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct EMAN_Aligner_Wrapper: EMAN::Aligner
{
    EMAN_Aligner_Wrapper(PyObject* self_, const EMAN::Aligner& p0):
        EMAN::Aligner(p0), self(self_) {}

    EMAN_Aligner_Wrapper(PyObject* self_):
        EMAN::Aligner(), self(self_) {}

    EMAN::EMData* align(EMAN::EMData* p0, std::string p1) const {
        return call_method< EMAN::EMData* >(self, "align", p0, p1);
    }

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
    }

    EMAN::Dict get_params() const {
        return call_method< EMAN::Dict >(self, "get_params");
    }

    EMAN::Dict default_get_params() const {
        return EMAN::Aligner::get_params();
    }

    void set_params(const EMAN::Dict& p0) {
        call_method< void >(self, "set_params", p0);
    }

    void default_set_params(const EMAN::Dict& p0) {
        EMAN::Aligner::set_params(p0);
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(self, "get_param_types");
    }

    PyObject* self;
};


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyAligner2)
{
    class_< EMAN::Aligner, boost::noncopyable, EMAN_Aligner_Wrapper >("Aligner", init<  >())
        .def("align", pure_virtual(&EMAN::Aligner::align), return_value_policy< manage_new_object >())
        .def("get_name", pure_virtual(&EMAN::Aligner::get_name))
        .def("get_params", &EMAN::Aligner::get_params, &EMAN_Aligner_Wrapper::default_get_params)
        .def("set_params", &EMAN::Aligner::set_params, &EMAN_Aligner_Wrapper::default_set_params)
        .def("get_param_types", pure_virtual(&EMAN::Aligner::get_param_types))
    ;

    def("dump_aligners", &EMAN::dump_aligners);
    class_< EMAN::Factory<EMAN::Aligner>, boost::noncopyable >("Aligners", no_init)
        .def("get", (EMAN::Aligner* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >))&EMAN::Factory<EMAN::Aligner>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Aligner* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&))&EMAN::Factory<EMAN::Aligner>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Aligner>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

}

