
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <emdata.h>
#include <emobject.h>
#include <projector.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct EMAN_Projector_Wrapper: EMAN::Projector
{
    EMAN_Projector_Wrapper(PyObject* py_self_, const EMAN::Projector& p0):
        EMAN::Projector(p0), py_self(py_self_) {}

    EMAN_Projector_Wrapper(PyObject* py_self_):
        EMAN::Projector(), py_self(py_self_) {}

    EMAN::EMData* project3d(EMAN::EMData* p0) const {
        return call_method< EMAN::EMData* >(py_self, "project3d", p0);
    }

    std::string get_name() const {
        return call_method< std::string >(py_self, "get_name");
    }

    std::string get_desc() const {
        return call_method< std::string >(py_self, "get_desc");
    }

    EMAN::Dict get_params() const {
        return call_method< EMAN::Dict >(py_self, "get_params");
    }

    EMAN::Dict default_get_params() const {
        return EMAN::Projector::get_params();
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(py_self, "get_param_types");
    }

    EMAN::TypeDict default_get_param_types() const {
        return EMAN::Projector::get_param_types();
    }

    PyObject* py_self;
};


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyProjector2)
{
    def("dump_projectors", &EMAN::dump_projectors);
    class_< EMAN::Projector, boost::noncopyable, EMAN_Projector_Wrapper >("__Projector", init<  >())
        .def("project3d", pure_virtual(&EMAN::Projector::project3d), return_value_policy< manage_new_object >())
        .def("get_name", pure_virtual(&EMAN::Projector::get_name))
        .def("get_desc", pure_virtual(&EMAN::Projector::get_desc))
        .def("get_params", &EMAN::Projector::get_params, &EMAN_Projector_Wrapper::default_get_params)
        .def("get_param_types", &EMAN::Projector::get_param_types, &EMAN_Projector_Wrapper::default_get_param_types)
        .def("set_params", &EMAN::Projector::set_params)
    ;

    class_< EMAN::Factory<EMAN::Projector>, boost::noncopyable >("Projectors", no_init)
        .def("get", (EMAN::Projector* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&))&EMAN::Factory<EMAN::Projector>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Projector* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&, const EMAN::Dict&))&EMAN::Factory<EMAN::Projector>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Projector>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

}

