
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <emdata.h>
#include <emobject.h>
#include <reconstructor.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct EMAN_Reconstructor_Wrapper: EMAN::Reconstructor
{
    EMAN_Reconstructor_Wrapper(PyObject* self_, const EMAN::Reconstructor& p0):
        EMAN::Reconstructor(p0), self(self_) {}

    EMAN_Reconstructor_Wrapper(PyObject* self_):
        EMAN::Reconstructor(), self(self_) {}

    void setup() {
        call_method< void >(self, "setup");
    }

    int insert_slice(EMAN::EMData* p0, const EMAN::Transform& p1) {
        return call_method< int >(self, "insert_slice", p0, p1);
    }

    EMAN::EMData* finish() {
        return call_method< EMAN::EMData* >(self, "finish");
    }

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
    }

    EMAN::Dict get_params() const {
        return call_method< EMAN::Dict >(self, "get_params");
    }

    EMAN::Dict default_get_params() const {
        return EMAN::Reconstructor::get_params();
    }

    void set_params(const EMAN::Dict& p0) {
        call_method< void >(self, "set_params", p0);
    }

    void default_set_params(const EMAN::Dict& p0) {
        EMAN::Reconstructor::set_params(p0);
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(self, "get_param_types");
    }

    PyObject* self;
};


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyReconstructor2)
{
    def("dump_reconstructors", &EMAN::dump_reconstructors);
    class_< EMAN::Reconstructor, boost::noncopyable, EMAN_Reconstructor_Wrapper >("__Reconstructor", init<  >())
        .def("setup", pure_virtual(&EMAN::Reconstructor::setup))
        .def("insert_slice", pure_virtual(&EMAN::Reconstructor::insert_slice))
        .def("finish", pure_virtual(&EMAN::Reconstructor::finish), return_value_policy< manage_new_object >())
        .def("get_name", pure_virtual(&EMAN::Reconstructor::get_name))
        .def("get_params", &EMAN::Reconstructor::get_params, &EMAN_Reconstructor_Wrapper::default_get_params)
        .def("set_params", &EMAN::Reconstructor::set_params, &EMAN_Reconstructor_Wrapper::default_set_params)
        .def("get_param_types", pure_virtual(&EMAN::Reconstructor::get_param_types))
    ;

    class_< EMAN::Factory<EMAN::Reconstructor>, boost::noncopyable >("Reconstructors", no_init)
        .def("get", (EMAN::Reconstructor* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >))&EMAN::Factory<EMAN::Reconstructor>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Reconstructor* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&))&EMAN::Factory<EMAN::Reconstructor>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Reconstructor>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

}

