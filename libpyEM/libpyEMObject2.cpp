
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <emdata.h>
#include <emobject.h>
#include <xydata.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_TypeDict_put_overloads_2_3, put, 2, 3)


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyEMObject2)
{
    scope* EMAN_EMObject_scope = new scope(
    class_< EMAN::EMObject >("EMObject", init<  >())
        .def(init< const EMAN::EMObject& >())
        .def(init< int >())
        .def(init< float >())
        .def(init< double >())
        .def(init< const char* >())
        .def(init< const std::string& >())
        .def(init< EMAN::EMData* >())
        .def(init< EMAN::XYData* >())
        .def(init< const std::vector<float,std::allocator<float> >& >())
        .def(init< const std::vector<std::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::allocator<std::basic_string<char, std::char_traits<char>, std::allocator<char> > > >& >())
        .def("is_null", &EMAN::EMObject::is_null)
        .def("to_str", &EMAN::EMObject::to_str)
        .def("get_object_type_name", &EMAN::EMObject::get_object_type_name)
        .staticmethod("get_object_type_name")
        .def( self == self )
        .def( self != self )
        .def("__int__", &EMAN::EMObject::operator int)
        .def("__float__", &EMAN::EMObject::operator float)
        .def("__float__", &EMAN::EMObject::operator double)
        .def("to_EMAN_EMData", &EMAN::EMObject::operator EMAN::EMData*, return_internal_reference< 1 >())
        .def("to_EMAN_XYData", &EMAN::EMObject::operator EMAN::XYData*, return_internal_reference< 1 >())
        .def("to_std_vector_float_std_allocator_float", &EMAN::EMObject::operator std::vector<float,std::allocator<float> >)
        .def("__str__", &EMAN::EMObject::operator std::vector<std::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::allocator<std::basic_string<char, std::char_traits<char>, std::allocator<char> > > >)
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
        .value("STRINGARRAY", EMAN::EMObject::STRINGARRAY)
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

