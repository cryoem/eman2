
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <emdata.h>
#include <emobject.h>
#include <exception.h>
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

    const char* what() const throw() {
        return call_method< const char* >(self, "what");
    }

    const char* default_what() const {
        return EMAN::Exception::what();
    }

    PyObject* self;
};

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_TypeDict_put_overloads_2_3, put, 2, 3)


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyEMObject2)
{
    class_< EMAN::Exception, EMAN_Exception_Wrapper >("Exception", init< const EMAN::Exception& >())
        .def(init< optional< const std::string&, int, const std::string&, const std::string& > >())
        .def_readwrite("filename", &EMAN::Exception::filename)
        .def_readwrite("line", &EMAN::Exception::line)
        .def_readwrite("desc", &EMAN::Exception::desc)
        .def_readwrite("objname", &EMAN::Exception::objname)
        .def("what", (const char* (EMAN::Exception::*)() const throw())&EMAN::Exception::what, (const char* (EMAN_Exception_Wrapper::*)() const)&EMAN_Exception_Wrapper::default_what)
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

    class_< EMAN::TypeDict >("TypeDict", init<  >())
        .def(init< const EMAN::TypeDict& >())
        .def("keys", &EMAN::TypeDict::keys)
        .def("size", &EMAN::TypeDict::size)
        .def("put", &EMAN::TypeDict::put, EMAN_TypeDict_put_overloads_2_3())
        .def("get_type", &EMAN::TypeDict::get_type)
        .def("get_desc", &EMAN::TypeDict::get_desc)
        .def("dump", &EMAN::TypeDict::dump)
    ;

}

