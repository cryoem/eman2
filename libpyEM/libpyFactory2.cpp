
// Boost Includes ==============================================================
#include <Python.h>
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <aligner.h>
#include <averager.h>
#include <cmp.h>
#include <ctf.h>
#include <emdata.h>
#include <emobject.h>
#include <exception.h>
#include <filter.h>
#include <interp.h>
#include <projector.h>
#include <reconstructor.h>
#include <util.h>
#include <xydata.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct EMAN_Exception_Wrapper: EMAN::Exception
{
    EMAN_Exception_Wrapper(PyObject* self_, const EMAN::Exception& p0):
        EMAN::Exception(p0), self(self_) {}

    EMAN_Exception_Wrapper(PyObject* self_):
        EMAN::Exception(), self(self_) {}

    EMAN_Exception_Wrapper(PyObject* self_, const std::string& p0):
        EMAN::Exception(p0), self(self_) {}

    EMAN_Exception_Wrapper(PyObject* self_, const std::string& p0, int p1):
        EMAN::Exception(p0, p1), self(self_) {}

    EMAN_Exception_Wrapper(PyObject* self_, const std::string& p0, int p1, const std::string& p2):
        EMAN::Exception(p0, p1, p2), self(self_) {}

    EMAN_Exception_Wrapper(PyObject* self_, const std::string& p0, int p1, const std::string& p2, const std::string& p3):
        EMAN::Exception(p0, p1, p2, p3), self(self_) {}

    void set_desc(const std::string& p0) {
        call_method< void >(self, "set_desc", p0);
    }

    void default_set_desc(const std::string& p0) {
        EMAN::Exception::set_desc(p0);
    }

    const char* get_file() const {
        return call_method< const char* >(self, "get_file");
    }

    const char* default_get_file() const {
        return EMAN::Exception::get_file();
    }

    const char* get_desc() const {
        return call_method< const char* >(self, "get_desc");
    }

    const char* default_get_desc() const {
        return EMAN::Exception::get_desc();
    }

    int get_line_num() const {
        return call_method< int >(self, "get_line_num");
    }

    int default_get_line_num() const {
        return EMAN::Exception::get_line_num();
    }

    const char* get_objname() const {
        return call_method< const char* >(self, "get_objname");
    }

    const char* default_get_objname() const {
        return EMAN::Exception::get_objname();
    }

    const char* what() const throw() {
        return call_method< const char* >(self, "what");
    }

    const char* default_what() const {
        return std::exception::what();
    }

    PyObject* self;
};

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

    void add_image(EMAN::EMData* p0) {
        call_method< void >(self, "add_image", p0);
    }

    EMAN::EMData* finish() {
        return call_method< EMAN::EMData* >(self, "finish");
    }

    void add_image_list(const std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >& p0) {
        call_method< void >(self, "add_image_list", p0);
    }

    void default_add_image_list(const std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >& p0) {
        EMAN::Averager::add_image_list(p0);
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

    EMAN::TypeDict default_get_param_types() const {
        return EMAN::Filter::get_param_types();
    }

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
    }

    PyObject* self;
};


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyFactory2)
{
    class_< EMAN::Exception, EMAN_Exception_Wrapper >("Exception", init< const EMAN::Exception& >())
        .def(init< optional< const std::string&, int, const std::string&, const std::string& > >())
        .def("set_desc", &EMAN::Exception::set_desc, &EMAN_Exception_Wrapper::default_set_desc)
        .def("get_file", &EMAN::Exception::get_file, &EMAN_Exception_Wrapper::default_get_file)
        .def("get_desc", &EMAN::Exception::get_desc, &EMAN_Exception_Wrapper::default_get_desc)
        .def("get_line_num", &EMAN::Exception::get_line_num, &EMAN_Exception_Wrapper::default_get_line_num)
        .def("get_objname", &EMAN::Exception::get_objname, &EMAN_Exception_Wrapper::default_get_objname)
        .def("what", (const char* (std::exception::*)() const throw())&std::exception::what, (const char* (EMAN_Exception_Wrapper::*)() const)&EMAN_Exception_Wrapper::default_what)
        .def("dump", &EMAN::Exception::dump)
    ;

    scope* EMAN_EMObject_scope = new scope(
    class_< EMAN::EMObject >("EMObject", init<  >())
        .def(init< const EMAN::EMObject& >())
        .def(init< int >())
        .def(init< float >())
        .def(init< double >())
        .def(init< const char* >())
        .def(init< std::string >())
        .def(init< EMAN::EMData* >())
        .def(init< EMAN::XYData* >())
        .def(init< const std::vector<float,std::allocator<float> >& >())
        .def("get_farray", &EMAN::EMObject::get_farray)
        .def("is_null", &EMAN::EMObject::is_null)
        .def("to_str", &EMAN::EMObject::to_str)
        .def("get_object_type_name", &EMAN::EMObject::get_object_type_name)
        .staticmethod("get_object_type_name")
        .def("__int__", &EMAN::EMObject::operator int)
        .def("__float__", &EMAN::EMObject::operator float)
        .def("__float__", &EMAN::EMObject::operator double)
        .def("to_EMAN_EMData", &EMAN::EMObject::operator EMAN::EMData*, return_internal_reference< 1 >())
        .def("to_EMAN_XYData", &EMAN::EMObject::operator EMAN::XYData*, return_internal_reference< 1 >())
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

    class_< EMAN::NotExistingObjectError, bases< EMAN::Exception >  >("NotExistingObjectError", init< const EMAN::NotExistingObjectError& >())
        .def(init< const std::string&, optional< const std::string&, int, const std::string& > >())
    ;

    class_< EMAN::Aligner, boost::noncopyable, EMAN_Aligner_Wrapper >("Aligner", init<  >())
        .def("align", pure_virtual(&EMAN::Aligner::align), return_value_policy< manage_new_object >())
        .def("get_params", &EMAN::Aligner::get_params, &EMAN_Aligner_Wrapper::default_get_params)
        .def("set_params", &EMAN::Aligner::set_params, &EMAN_Aligner_Wrapper::default_set_params)
        .def("get_param_types", pure_virtual(&EMAN::Aligner::get_param_types))
        .def("get_name", pure_virtual(&EMAN::Aligner::get_name))
    ;

    class_< EMAN::Factory<EMAN::Aligner>, boost::noncopyable >("Aligners", no_init)
        .def("get", (EMAN::Aligner* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >))&EMAN::Factory<EMAN::Aligner>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Aligner* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&))&EMAN::Factory<EMAN::Aligner>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Aligner>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

    class_< EMAN::Cmp, boost::noncopyable, EMAN_Cmp_Wrapper >("Cmp", init<  >())
        .def("cmp", pure_virtual(&EMAN::Cmp::cmp))
        .def("get_param_types", pure_virtual(&EMAN::Cmp::get_param_types))
        .def("get_name", pure_virtual(&EMAN::Cmp::get_name))
        .def("get_params", &EMAN::Cmp::get_params, &EMAN_Cmp_Wrapper::default_get_params)
        .def("set_params", &EMAN::Cmp::set_params, &EMAN_Cmp_Wrapper::default_set_params)
    ;

    class_< EMAN::Factory<EMAN::Cmp>, boost::noncopyable >("Cmps", no_init)
        .def("get", (EMAN::Cmp* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >))&EMAN::Factory<EMAN::Cmp>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Cmp* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&))&EMAN::Factory<EMAN::Cmp>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Cmp>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

    class_< EMAN::Averager, boost::noncopyable, EMAN_Averager_Wrapper >("Averager", init<  >())
        .def("add_image", pure_virtual(&EMAN::Averager::add_image))
        .def("finish", pure_virtual(&EMAN::Averager::finish), return_value_policy< manage_new_object >())
        .def("add_image_list", &EMAN::Averager::add_image_list, &EMAN_Averager_Wrapper::default_add_image_list)
        .def("get_name", pure_virtual(&EMAN::Averager::get_name))
        .def("set_params", &EMAN::Averager::set_params, &EMAN_Averager_Wrapper::default_set_params)
        .def("get_param_types", &EMAN::Averager::get_param_types, &EMAN_Averager_Wrapper::default_get_param_types)
    ;

    class_< EMAN::Factory<EMAN::Averager>, boost::noncopyable >("Averagers", no_init)
        .def("get", (EMAN::Averager* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >))&EMAN::Factory<EMAN::Averager>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Averager* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&))&EMAN::Factory<EMAN::Averager>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Averager>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

    class_< EMAN::Projector, boost::noncopyable, EMAN_Projector_Wrapper >("Projector", init<  >())
        .def("project3d", pure_virtual(&EMAN::Projector::project3d), return_value_policy< manage_new_object >())
        .def("get_name", pure_virtual(&EMAN::Projector::get_name))
        .def("get_params", &EMAN::Projector::get_params, &EMAN_Projector_Wrapper::default_get_params)
        .def("get_param_types", &EMAN::Projector::get_param_types, &EMAN_Projector_Wrapper::default_get_param_types)
        .def("set_params", &EMAN::Projector::set_params)
    ;

    class_< EMAN::Factory<EMAN::Projector>, boost::noncopyable >("Projectors", no_init)
        .def("get", (EMAN::Projector* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >))&EMAN::Factory<EMAN::Projector>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Projector* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&))&EMAN::Factory<EMAN::Projector>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Projector>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
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

    class_< EMAN::Factory<EMAN::Reconstructor>, boost::noncopyable >("Reconstructors", no_init)
        .def("get", (EMAN::Reconstructor* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >))&EMAN::Factory<EMAN::Reconstructor>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Reconstructor* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&))&EMAN::Factory<EMAN::Reconstructor>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Reconstructor>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

    class_< EMAN::Filter, boost::noncopyable, EMAN_Filter_Wrapper >("Filter", init<  >())
        .def("process", &EMAN::Filter::process, &EMAN_Filter_Wrapper::default_process)
        .def("process_list", &EMAN::Filter::process_list, &EMAN_Filter_Wrapper::default_process_list)
        .def("get_params", &EMAN::Filter::get_params, &EMAN_Filter_Wrapper::default_get_params)
        .def("set_params", &EMAN::Filter::set_params, &EMAN_Filter_Wrapper::default_set_params)
        .def("get_param_types", &EMAN::Filter::get_param_types, &EMAN_Filter_Wrapper::default_get_param_types)
        .def("get_name", pure_virtual(&EMAN::Filter::get_name))
    ;

    class_< EMAN::Factory<EMAN::Filter>, boost::noncopyable >("Filters", no_init)
        .def("get", (EMAN::Filter* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >))&EMAN::Factory<EMAN::Filter>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Filter* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&))&EMAN::Factory<EMAN::Filter>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Filter>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

}

