
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <quaternion.h>
#include <transform.h>
#include <vec3.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_Transform_set_rotation_overloads_4_5, set_rotation, 4, 5)


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyTransform2)
{
    class_< EMAN::Vec3f >("Vec3f", init<  >())
        .def(init< const EMAN::Vec3f& >())
        .def(init< float, float, optional< float > >())
        .def(init< const std::vector<float,std::allocator<float> >& >())
        .def(init< const EMAN::Vec3i& >())
        .def("normalize", &EMAN::Vec3f::normalize)
        .def("length", &EMAN::Vec3f::length)
        .def("dot", &EMAN::Vec3f::dot)
        .def("cross", &EMAN::Vec3f::cross)
        .def("as_list", &EMAN::Vec3f::as_list)
        .def("set_value", (void (EMAN::Vec3f::*)(const std::vector<float,std::allocator<float> >&) )&EMAN::Vec3f::set_value)
        .def("set_value", (void (EMAN::Vec3f::*)(float, float, float) )&EMAN::Vec3f::set_value)
        .def("at", &EMAN::Vec3f::at)
        .def( other< float >() - self )
        .def( self - other< float >() )
        .def( self * other< float >() )
        .def( self / other< float >() )
        .def( self == self )
        .def( self != self )
        .def( other< float >() + self )
        .def( self + other< float >() )
        .def( other< float >() * self )
        .def( other< EMAN::Vec3i >() * self )
        .def( self * other< EMAN::Vec3i >() )
        .def( -self )
        .def( self - other< EMAN::Vec3i >() )
        .def( other< EMAN::Vec3i >() - self )
        .def( self + other< EMAN::Vec3i >() )
        .def( other< EMAN::Vec3i >() + self )
        .def( self * self )
        .def( self - self )
        .def( self + self )
        .def( self += self )
        .def( self += other< EMAN::Vec3i >() )
        .def( self += other< float >() )
        .def( self -= self )
        .def( self -= other< EMAN::Vec3i >() )
        .def( self -= other< float >() )
        .def( self *= other< float >() )
        .def( self /= other< float >() )
    ;

    class_< EMAN::Vec3i >("Vec3i", init<  >())
        .def(init< const EMAN::Vec3i& >())
        .def(init< int, int, optional< int > >())
        .def(init< const std::vector<int,std::allocator<int> >& >())
        .def("normalize", &EMAN::Vec3i::normalize)
        .def("length", &EMAN::Vec3i::length)
        .def("dot", &EMAN::Vec3i::dot)
        .def("cross", &EMAN::Vec3i::cross)
        .def("as_list", &EMAN::Vec3i::as_list)
        .def("set_value", (void (EMAN::Vec3i::*)(const std::vector<int,std::allocator<int> >&) )&EMAN::Vec3i::set_value)
        .def("set_value", (void (EMAN::Vec3i::*)(int, int, int) )&EMAN::Vec3i::set_value)
        .def( self + self )
        .def( other< int >() + self )
        .def( self + other< int >() )
        .def( other< int >() * self )
        .def( self * other< EMAN::Vec3f >() )
        .def( other< EMAN::Vec3f >() * self )
        .def( other< EMAN::Vec3f >() - self )
        .def( self - other< EMAN::Vec3f >() )
        .def( other< EMAN::Vec3f >() + self )
        .def( self + other< EMAN::Vec3f >() )
        .def( self * other< int >() )
        .def( other< int >() - self )
        .def( self * self )
        .def( self != self )
        .def( self / other< int >() )
        .def( self == self )
        .def( self - self )
        .def( self - other< int >() )
        .def( self += self )
        .def( self += other< int >() )
        .def( self -= self )
        .def( self -= other< int >() )
        .def( self *= other< int >() )
        .def( self /= other< int >() )
    ;

    class_< EMAN::Quaternion >("Quaternion", init<  >())
        .def(init< const EMAN::Quaternion& >())
        .def(init< float, float, float, float >())
        .def(init< float, const EMAN::Vec3f& >())
        .def(init< const EMAN::Vec3f&, float >())
        .def(init< const std::vector<float,std::allocator<float> >& >())
        .def(init< const EMAN::Matrix3f& >())
        .def("norm", &EMAN::Quaternion::norm)
        .def("conj", &EMAN::Quaternion::conj)
        .def("abs", &EMAN::Quaternion::abs)
        .def("normalize", &EMAN::Quaternion::normalize)
        .def("inverse", &EMAN::Quaternion::inverse, return_internal_reference< 1 >())
        .def("create_inverse", &EMAN::Quaternion::create_inverse)
        .def("rotate", &EMAN::Quaternion::rotate)
        .def("to_angle", &EMAN::Quaternion::to_angle)
        .def("to_axis", &EMAN::Quaternion::to_axis)
        .def("to_matrix3", &EMAN::Quaternion::to_matrix3)
        .def("real", &EMAN::Quaternion::real)
        .def("unreal", &EMAN::Quaternion::unreal)
        .def("as_list", &EMAN::Quaternion::as_list)
        .def("interpolate", &EMAN::Quaternion::interpolate)
        .staticmethod("interpolate")
        .def( self * other< float >() )
        .def( other< float >() * self )
        .def( self + self )
        .def( self - self )
        .def( self * self )
        .def( self / self )
        .def( self == self )
        .def( self != self )
        .def( self += self )
        .def( self -= self )
        .def( self *= self )
        .def( self *= other< float >() )
        .def( self /= self )
        .def( self /= other< float >() )
    ;

    scope* EMAN_Transform_scope = new scope(
    class_< EMAN::Transform >("Transform", init<  >())
        .def(init< const EMAN::Transform& >())
        .def(init< EMAN::Transform::EulerType, float, float, float, optional< float > >())
        .def(init< const EMAN::Vec3f&, EMAN::Transform::EulerType, float, float, float, optional< float > >())
        .def(init< const EMAN::Vec3f&, const EMAN::Vec3f&, EMAN::Transform::EulerType, float, float, float, optional< float > >())
        .def(init< EMAN::Transform::EulerType, std::map<std::basic_string<char, std::char_traits<char>, std::allocator<char> >,float,std::less<std::basic_string<char, std::char_traits<char>, std::allocator<char> > >,std::allocator<std::pair<const std::basic_string<char, std::char_traits<char>, std::allocator<char> >, float> > >& >())
        .def(init< const EMAN::Vec3f&, EMAN::Transform::EulerType, std::map<std::basic_string<char, std::char_traits<char>, std::allocator<char> >,float,std::less<std::basic_string<char, std::char_traits<char>, std::allocator<char> > >,std::allocator<std::pair<const std::basic_string<char, std::char_traits<char>, std::allocator<char> >, float> > >& >())
        .def(init< const EMAN::Vec3f&, const EMAN::Vec3f&, EMAN::Transform::EulerType, std::map<std::basic_string<char, std::char_traits<char>, std::allocator<char> >,float,std::less<std::basic_string<char, std::char_traits<char>, std::allocator<char> > >,std::allocator<std::pair<const std::basic_string<char, std::char_traits<char>, std::allocator<char> >, float> > >& >())
        .def_readonly("ERR_LIMIT", &EMAN::Transform::ERR_LIMIT)
        .def("to_identity", &EMAN::Transform::to_identity)
        .def("is_identity", &EMAN::Transform::is_identity)
        .def("orthogonalize", &EMAN::Transform::orthogonalize)
        .def("inverse", &EMAN::Transform::inverse)
        .def("set_pretrans", &EMAN::Transform::set_pretrans)
        .def("set_posttrans", &EMAN::Transform::set_posttrans)
        .def("set_center", &EMAN::Transform::set_center)
        .def("set_rotation", (void (EMAN::Transform::*)(EMAN::Transform::EulerType, float, float, float, float) )&EMAN::Transform::set_rotation, EMAN_Transform_set_rotation_overloads_4_5())
        .def("set_rotation", (void (EMAN::Transform::*)(EMAN::Transform::EulerType, std::map<std::basic_string<char, std::char_traits<char>, std::allocator<char> >,float,std::less<std::basic_string<char, std::char_traits<char>, std::allocator<char> > >,std::allocator<std::pair<const std::basic_string<char, std::char_traits<char>, std::allocator<char> >, float> > >&) )&EMAN::Transform::set_rotation)
        .def("set_scale", &EMAN::Transform::set_scale)
        .def("get_pretrans", &EMAN::Transform::get_pretrans)
        .def("get_posttrans", &EMAN::Transform::get_posttrans)
        .def("get_center", &EMAN::Transform::get_center)
        .def("get_rotation", &EMAN::Transform::get_rotation)
        .def("get_matrix3_col", &EMAN::Transform::get_matrix3_col)
        .def("get_matrix3_row", &EMAN::Transform::get_matrix3_row)
        .def("get_scale", &EMAN::Transform::get_scale)
        .def("orthogonality", &EMAN::Transform::orthogonality)
        .def("get_nsym", &EMAN::Transform::get_nsym)
        .def("get_sym", &EMAN::Transform::get_sym)
        .def("at", &EMAN::Transform::at)
        .staticmethod("get_nsym")
        .def( self * self )
        .def( other< EMAN::Vec3f >() * self )
        .def( self *= self )
    );

    enum_< EMAN::Transform::EulerType >("EulerType")
        .value("UNKNOWN", EMAN::Transform::UNKNOWN)
        .value("IMAGIC", EMAN::Transform::IMAGIC)
        .value("SPIDER", EMAN::Transform::SPIDER)
        .value("QUATERNION", EMAN::Transform::QUATERNION)
        .value("SGIROT", EMAN::Transform::SGIROT)
        .value("MRC", EMAN::Transform::MRC)
        .value("SPIN", EMAN::Transform::SPIN)
        .value("EMAN", EMAN::Transform::EMAN)
    ;

    delete EMAN_Transform_scope;

}

