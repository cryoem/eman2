
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <cmp.h>
#include <emdata.h>
#include <emobject.h>

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

    float cmp(EMAN::EMData* p0, EMAN::Transform* p1) const {
        return call_method< float >(self, "cmp", p0, p1);
    }

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
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


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyCmp2)
{
    class_< EMAN::Cmp, boost::noncopyable, EMAN_Cmp_Wrapper >("Cmp", init<  >())
        .def("cmp", pure_virtual(&EMAN::Cmp::cmp))
        .def("get_name", pure_virtual(&EMAN::Cmp::get_name))
        .def("get_params", &EMAN::Cmp::get_params, &EMAN_Cmp_Wrapper::default_get_params)
        .def("set_params", &EMAN::Cmp::set_params, &EMAN_Cmp_Wrapper::default_set_params)
        .def("get_param_types", pure_virtual(&EMAN::Cmp::get_param_types))
    ;

    def("dump_cmps", &EMAN::dump_cmps);
    class_< EMAN::Factory<EMAN::Cmp>, boost::noncopyable >("Cmps", no_init)
        .def("get", (EMAN::Cmp* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >))&EMAN::Factory<EMAN::Cmp>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Cmp* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&))&EMAN::Factory<EMAN::Cmp>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Cmp>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

}

