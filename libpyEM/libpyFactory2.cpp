
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <aligner.h>
#include <averager.h>
#include <cmp.h>
#include <ctf.h>
#include <emdata.h>
#include <emobject.h>
#include <filter.h>
#include <interp.h>
#include <log.h>
#include <projector.h>
#include <reconstructor.h>
#include <util.h>
#include <xydata.h>

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

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
    }

    PyObject* self;
};

struct EMAN_Cmp_Wrapper: EMAN::Cmp
{
    EMAN_Cmp_Wrapper(PyObject* self_, const EMAN::Cmp& p0):
        EMAN::Cmp(p0), self(self_) {}

    EMAN_Cmp_Wrapper(PyObject* self_):
        EMAN::Cmp(), self(self_) {}

    float cmp(EMAN::EMData* p0, EMAN::Transform* p1) const {
        return call_method< float >(self, "cmp", p0, p1);
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(self, "get_param_types");
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

    PyObject* self;
};

struct EMAN_Averager_Wrapper: EMAN::Averager
{
    EMAN_Averager_Wrapper(PyObject* self_, const EMAN::Averager& p0):
        EMAN::Averager(p0), self(self_) {}

    EMAN_Averager_Wrapper(PyObject* self_):
        EMAN::Averager(), self(self_) {}

    EMAN::EMData* average(const std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >& p0) const {
        return call_method< EMAN::EMData* >(self, "average", p0);
    }

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
    }

    void set_params(const EMAN::Dict& p0) {
        call_method< void >(self, "set_params", p0);
    }

    void default_set_params(const EMAN::Dict& p0) {
        EMAN::Averager::set_params(p0);
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(self, "get_param_types");
    }

    EMAN::TypeDict default_get_param_types() const {
        return EMAN::Averager::get_param_types();
    }

    PyObject* self;
};

struct EMAN_Projector_Wrapper: EMAN::Projector
{
    EMAN_Projector_Wrapper(PyObject* self_, const EMAN::Projector& p0):
        EMAN::Projector(p0), self(self_) {}

    EMAN_Projector_Wrapper(PyObject* self_):
        EMAN::Projector(), self(self_) {}

    EMAN::EMData* project3d(EMAN::EMData* p0) const {
        return call_method< EMAN::EMData* >(self, "project3d", p0);
    }

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
    }

    EMAN::Dict get_params() const {
        return call_method< EMAN::Dict >(self, "get_params");
    }

    EMAN::Dict default_get_params() const {
        return EMAN::Projector::get_params();
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(self, "get_param_types");
    }

    EMAN::TypeDict default_get_param_types() const {
        return EMAN::Projector::get_param_types();
    }

    PyObject* self;
};

struct EMAN_Reconstructor_Wrapper: EMAN::Reconstructor
{
    EMAN_Reconstructor_Wrapper(PyObject* self_, const EMAN::Reconstructor& p0):
        EMAN::Reconstructor(p0), self(self_) {}

    EMAN_Reconstructor_Wrapper(PyObject* self_):
        EMAN::Reconstructor(), self(self_) {}

    int setup() {
        return call_method< int >(self, "setup");
    }

    int insert_slice(EMAN::EMData* p0, const EMAN::Rotation& p1) {
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
BOOST_PYTHON_MODULE(libpyFactory2)
{
    scope* EMAN_EMObject_scope = new scope(
    class_< EMAN::EMObject >("EMObject", init<  >())
        .def(init< const EMAN::EMObject& >())
        .def(init< int >())
        .def(init< float >())
        .def(init< double >())
        .def(init< std::string >())
        .def(init< EMAN::EMData* >())
        .def(init< EMAN::XYData* >())
        .def(init< const std::vector<float,std::allocator<float> >& >())
        .def("get_int", &EMAN::EMObject::get_int)
        .def("get_float", &EMAN::EMObject::get_float)
        .def("get_double", &EMAN::EMObject::get_double)
        .def("get_string", &EMAN::EMObject::get_string)
        .def("get_emdata", &EMAN::EMObject::get_emdata, return_internal_reference< 1 >())
        .def("get_xydata", &EMAN::EMObject::get_xydata, return_internal_reference< 1 >())
        .def("get_farray", &EMAN::EMObject::get_farray)
        .def("is_null", &EMAN::EMObject::is_null)
        .def("to_str", &EMAN::EMObject::to_str)
        .def("get_object_type_name", &EMAN::EMObject::get_object_type_name)
        .staticmethod("get_object_type_name")
    );

    enum_< EMAN::EMObject::ObjectType >("ObjectType")
        .value("EMDATA", EMAN::EMObject::EMDATA)
        .value("STRING", EMAN::EMObject::STRING)
        .value("INT", EMAN::EMObject::INT)
        .value("DOUBLE", EMAN::EMObject::DOUBLE)
        .value("FLOAT", EMAN::EMObject::FLOAT)
        .value("XYDATA", EMAN::EMObject::XYDATA)
        .value("FLOATARRAY", EMAN::EMObject::FLOATARRAY)
        .value("UNKNOWN", EMAN::EMObject::UNKNOWN)
    ;

    delete EMAN_EMObject_scope;

    class_< EMAN::Aligner, boost::noncopyable, EMAN_Aligner_Wrapper >("Aligner", init<  >())
        .def("align", pure_virtual(&EMAN::Aligner::align), return_value_policy< manage_new_object >())
        .def("get_params", &EMAN::Aligner::get_params, &EMAN_Aligner_Wrapper::default_get_params)
        .def("set_params", &EMAN::Aligner::set_params, &EMAN_Aligner_Wrapper::default_set_params)
        .def("get_param_types", pure_virtual(&EMAN::Aligner::get_param_types))
        .def("get_name", pure_virtual(&EMAN::Aligner::get_name))
    ;

    class_< EMAN::Factory<EMAN::Aligner>, boost::noncopyable >("AlignerFactory", no_init)
        .def("instance", &EMAN::Factory<EMAN::Aligner>::instance, return_value_policy< reference_existing_object >())
        .def("get", (EMAN::Aligner* (EMAN::Factory<EMAN::Aligner>::*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >) )&EMAN::Factory<EMAN::Aligner>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Aligner* (EMAN::Factory<EMAN::Aligner>::*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&) )&EMAN::Factory<EMAN::Aligner>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Aligner>::get_list)
        .staticmethod("instance")
    ;

    class_< EMAN::Cmp, boost::noncopyable, EMAN_Cmp_Wrapper >("Cmp", init<  >())
        .def("cmp", pure_virtual(&EMAN::Cmp::cmp))
        .def("get_param_types", pure_virtual(&EMAN::Cmp::get_param_types))
        .def("get_name", pure_virtual(&EMAN::Cmp::get_name))
        .def("get_params", &EMAN::Cmp::get_params, &EMAN_Cmp_Wrapper::default_get_params)
        .def("set_params", &EMAN::Cmp::set_params, &EMAN_Cmp_Wrapper::default_set_params)
    ;

    class_< EMAN::Factory<EMAN::Cmp>, boost::noncopyable >("CmpFactory", no_init)
        .def("instance", &EMAN::Factory<EMAN::Cmp>::instance, return_value_policy< reference_existing_object >())
        .def("get", (EMAN::Cmp* (EMAN::Factory<EMAN::Cmp>::*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >) )&EMAN::Factory<EMAN::Cmp>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Cmp* (EMAN::Factory<EMAN::Cmp>::*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&) )&EMAN::Factory<EMAN::Cmp>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Cmp>::get_list)
        .staticmethod("instance")
    ;

    class_< EMAN::Averager, boost::noncopyable, EMAN_Averager_Wrapper >("Averager", init<  >())
        .def("average", pure_virtual(&EMAN::Averager::average), return_value_policy< manage_new_object >())
        .def("get_name", pure_virtual(&EMAN::Averager::get_name))
        .def("set_params", &EMAN::Averager::set_params, &EMAN_Averager_Wrapper::default_set_params)
        .def("get_param_types", &EMAN::Averager::get_param_types, &EMAN_Averager_Wrapper::default_get_param_types)
    ;

    class_< EMAN::Factory<EMAN::Averager>, boost::noncopyable >("AveragerFactory", no_init)
        .def("instance", &EMAN::Factory<EMAN::Averager>::instance, return_value_policy< reference_existing_object >())
        .def("get", (EMAN::Averager* (EMAN::Factory<EMAN::Averager>::*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >) )&EMAN::Factory<EMAN::Averager>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Averager* (EMAN::Factory<EMAN::Averager>::*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&) )&EMAN::Factory<EMAN::Averager>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Averager>::get_list)
        .staticmethod("instance")
    ;

    class_< EMAN::Projector, boost::noncopyable, EMAN_Projector_Wrapper >("Projector", init<  >())
        .def("project3d", pure_virtual(&EMAN::Projector::project3d), return_value_policy< manage_new_object >())
        .def("get_name", pure_virtual(&EMAN::Projector::get_name))
        .def("get_params", &EMAN::Projector::get_params, &EMAN_Projector_Wrapper::default_get_params)
        .def("get_param_types", &EMAN::Projector::get_param_types, &EMAN_Projector_Wrapper::default_get_param_types)
        .def("set_params", &EMAN::Projector::set_params)
    ;

    class_< EMAN::Factory<EMAN::Projector>, boost::noncopyable >("ProjectorFactory", no_init)
        .def("instance", &EMAN::Factory<EMAN::Projector>::instance, return_value_policy< reference_existing_object >())
        .def("get", (EMAN::Projector* (EMAN::Factory<EMAN::Projector>::*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >) )&EMAN::Factory<EMAN::Projector>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Projector* (EMAN::Factory<EMAN::Projector>::*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&) )&EMAN::Factory<EMAN::Projector>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Projector>::get_list)
        .staticmethod("instance")
    ;

    class_< EMAN::Reconstructor, boost::noncopyable, EMAN_Reconstructor_Wrapper >("Reconstructor", init<  >())
        .def("setup", pure_virtual(&EMAN::Reconstructor::setup))
        .def("insert_slice", pure_virtual(&EMAN::Reconstructor::insert_slice))
        .def("finish", pure_virtual(&EMAN::Reconstructor::finish), return_value_policy< manage_new_object >())
        .def("get_name", pure_virtual(&EMAN::Reconstructor::get_name))
        .def("get_params", &EMAN::Reconstructor::get_params, &EMAN_Reconstructor_Wrapper::default_get_params)
        .def("set_params", &EMAN::Reconstructor::set_params, &EMAN_Reconstructor_Wrapper::default_set_params)
        .def("get_param_types", pure_virtual(&EMAN::Reconstructor::get_param_types))
    ;

    class_< EMAN::Factory<EMAN::Reconstructor>, boost::noncopyable >("ReconstructorFactory", no_init)
        .def("instance", &EMAN::Factory<EMAN::Reconstructor>::instance, return_value_policy< reference_existing_object >())
        .def("get", (EMAN::Reconstructor* (EMAN::Factory<EMAN::Reconstructor>::*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >) )&EMAN::Factory<EMAN::Reconstructor>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Reconstructor* (EMAN::Factory<EMAN::Reconstructor>::*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&) )&EMAN::Factory<EMAN::Reconstructor>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Reconstructor>::get_list)
        .staticmethod("instance")
    ;

}

