
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <averager.h>
#include <emdata.h>
#include <emobject.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct EMAN_Averager_Wrapper: EMAN::Averager
{
    EMAN_Averager_Wrapper(PyObject* py_self_, const EMAN::Averager& p0):
        EMAN::Averager(p0), py_self(py_self_) {}

    EMAN_Averager_Wrapper(PyObject* py_self_):
        EMAN::Averager(), py_self(py_self_) {}

    void add_image(EMAN::EMData* p0) {
        call_method< void >(py_self, "add_image", p0);
    }

    void add_image_list(const std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >& p0) {
        call_method< void >(py_self, "add_image_list", p0);
    }

    void default_add_image_list(const std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >& p0) {
        EMAN::Averager::add_image_list(p0);
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

    void set_params(const EMAN::Dict& p0) {
        call_method< void >(py_self, "set_params", p0);
    }

    void default_set_params(const EMAN::Dict& p0) {
        EMAN::Averager::set_params(p0);
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(py_self, "get_param_types");
    }

    EMAN::TypeDict default_get_param_types() const {
        return EMAN::Averager::get_param_types();
    }

    PyObject* py_self;
};


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyAverager2)
{
    class_< EMAN::Averager, boost::noncopyable, EMAN_Averager_Wrapper >("__Averager", init<  >())
        .def("add_image", pure_virtual(&EMAN::Averager::add_image))
        .def("add_image_list", &EMAN::Averager::add_image_list, &EMAN_Averager_Wrapper::default_add_image_list)
        .def("finish", pure_virtual(&EMAN::Averager::finish), return_value_policy< manage_new_object >())
        .def("get_name", pure_virtual(&EMAN::Averager::get_name))
        .def("get_desc", pure_virtual(&EMAN::Averager::get_desc))
        .def("set_params", &EMAN::Averager::set_params, &EMAN_Averager_Wrapper::default_set_params)
        .def("get_param_types", &EMAN::Averager::get_param_types, &EMAN_Averager_Wrapper::default_get_param_types)
    ;

    class_< EMAN::Factory<EMAN::Averager>, boost::noncopyable >("Averagers", no_init)
        .def("get", (EMAN::Averager* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&))&EMAN::Factory<EMAN::Averager>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Averager* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&, const EMAN::Dict&))&EMAN::Factory<EMAN::Averager>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Averager>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

    def("dump_averagers", &EMAN::dump_averagers);
    def("dump_averagers_list", &EMAN::dump_averagers_list);
}

