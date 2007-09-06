
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <emdata.h>
#include <emobject.h>
#include <reconstructor.h>

// Using =======================================================================
using namespace boost::python;

// debug
#include <iostream>
using std::cout;
using std::endl;

// Declarations ================================================================
namespace  {

struct EMAN_Reconstructor_Wrapper: EMAN::Reconstructor
{
//     EMAN_Reconstructor_Wrapper(PyObject* py_self_, const EMAN::Reconstructor& p0):
//         EMAN::Reconstructor(p0), py_self(py_self_) {}

    EMAN_Reconstructor_Wrapper(PyObject* py_self_):
        EMAN::Reconstructor(), py_self(py_self_) {}

    void setup() {
        call_method< void >(py_self, "setup");
    }

    int insert_slice(const EMAN::EMData* const p0, const EMAN::Transform3D& p1) {
        return call_method< int >(py_self, "insert_slice", p0, p1);
    }

    EMAN::EMData* finish() {
        return call_method< EMAN::EMData* >(py_self, "finish");
    }

    std::string get_name() const {
        return call_method< std::string >(py_self, "get_name");
    }

    std::string get_desc() const {
        return call_method< std::string >(py_self, "get_desc");
    }

    EMAN::Dict get_params() const {
		printf("call goes here!\n");
        return call_method< EMAN::Dict >(py_self, "get_params");
    }
/*
	void print_params() const {
        call_method< void >(py_self, "print_params");
	}*/

    EMAN::Dict default_get_params() const {
        return EMAN::Reconstructor::get_params();
    }

    void set_params(const EMAN::Dict& p0) {
        call_method< void >(py_self, "set_params", p0);
    }

    void default_set_params(const EMAN::Dict& p0) {
        EMAN::Reconstructor::set_params(p0);
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(py_self, "get_param_types");
    }

    PyObject* py_self;
};


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyReconstructor2)
{
    def("dump_reconstructors", &EMAN::dump_reconstructors);
    def("dump_reconstructors_list", &EMAN::dump_reconstructors_list);
    class_< EMAN::Reconstructor, boost::noncopyable, EMAN_Reconstructor_Wrapper >("__Reconstructor", init<  >())
        .def("setup", pure_virtual(&EMAN::Reconstructor::setup))
        .def("insert_slice", pure_virtual(&EMAN::Reconstructor::insert_slice))
		.def("determine_slice_agreement", &EMAN::Reconstructor::determine_slice_agreement)
        .def("finish", pure_virtual(&EMAN::Reconstructor::finish), return_internal_reference< 1 >())
        .def("get_name", pure_virtual(&EMAN::Reconstructor::get_name))
		.def("get_score", pure_virtual(&EMAN::Reconstructor::get_score))
		.def("get_norm", pure_virtual(&EMAN::Reconstructor::get_norm))
        .def("get_desc", pure_virtual(&EMAN::Reconstructor::get_desc))
		.def("get_emdata", (&EMAN::Reconstructor::get_emdata),  return_internal_reference< 1 >())
        //.def("get_params", &EMAN::Reconstructor::get_params, &EMAN_Reconstructor_Wrapper::default_get_params) 
		.def("get_params", &EMAN::Reconstructor::get_params)
		.def("insert_params", &EMAN::Reconstructor::insert_params)
        .def("set_params", &EMAN::Reconstructor::set_params, &EMAN_Reconstructor_Wrapper::default_set_params)
		.def("print_params",  &EMAN::Reconstructor::print_params) // Why is this different to set_params and get_params? Why is the wrapper needed? d.woolford May 2007
        .def("get_param_types", pure_virtual(&EMAN::Reconstructor::get_param_types))
    ;

    class_< EMAN::Factory<EMAN::Reconstructor>, boost::noncopyable >("Reconstructors", no_init)
        .def("get", (EMAN::Reconstructor* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&))&EMAN::Factory<EMAN::Reconstructor>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Reconstructor* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&, const EMAN::Dict&))&EMAN::Factory<EMAN::Reconstructor>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Reconstructor>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

}

