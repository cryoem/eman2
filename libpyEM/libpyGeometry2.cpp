
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <geometry.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
BOOST_PYTHON_MODULE(libpyGeometry2)
{
    class_< EMAN::IntSize >("IntSize", init<  >())
        .def(init< const EMAN::IntSize& >())
        .def(init< int, int >())
        .def(init< int, int, int >())
        .def_readwrite("x", &EMAN::IntSize::x)
        .def_readwrite("y", &EMAN::IntSize::y)
        .def_readwrite("z", &EMAN::IntSize::z)
        .def("get_ndim", &EMAN::IntSize::get_ndim)
    ;

    class_< EMAN::FloatSize >("FloatSize", init<  >())
        .def(init< const EMAN::FloatSize& >())
        .def(init< float, float >())
        .def(init< float, float, float >())
        .def_readwrite("x", &EMAN::FloatSize::x)
        .def_readwrite("y", &EMAN::FloatSize::y)
        .def_readwrite("z", &EMAN::FloatSize::z)
        .def("get_ndim", &EMAN::FloatSize::get_ndim)
    ;

    class_< EMAN::IntPoint >("IntPoint", init<  >())
        .def(init< const EMAN::IntPoint& >())
        .def(init< int, int >())
        .def(init< int, int, int >())
        .def_readwrite("x", &EMAN::IntPoint::x)
        .def_readwrite("y", &EMAN::IntPoint::y)
        .def_readwrite("z", &EMAN::IntPoint::z)
        .def("get_ndim", &EMAN::IntPoint::get_ndim)
    ;

    class_< EMAN::FloatPoint >("FloatPoint", init<  >())
        .def(init< const EMAN::FloatPoint& >())
        .def(init< float, float >())
        .def(init< float, float, float >())
        .def_readwrite("x", &EMAN::FloatPoint::x)
        .def_readwrite("y", &EMAN::FloatPoint::y)
        .def_readwrite("z", &EMAN::FloatPoint::z)
        .def("get_ndim", &EMAN::FloatPoint::get_ndim)
    ;

    class_< EMAN::Region >("Region", init<  >())
        .def(init< const EMAN::Region& >())
        .def(init< float, float, float, float >())
        .def(init< float, float, float, float, float, float >())
        .def(init< const EMAN::FloatPoint&, const EMAN::FloatSize& >())
        .def_readwrite("origin", &EMAN::Region::origin)
        .def_readwrite("size", &EMAN::Region::size)
        .def("inside_region", (bool (EMAN::Region::*)() const)&EMAN::Region::inside_region)
        .def("inside_region", (bool (EMAN::Region::*)(const EMAN::FloatPoint&) const)&EMAN::Region::inside_region)
        .def("inside_region", (bool (EMAN::Region::*)(float, float) const)&EMAN::Region::inside_region)
        .def("inside_region", (bool (EMAN::Region::*)(float, float, float) const)&EMAN::Region::inside_region)
        .def("get_ndim", &EMAN::Region::get_ndim)
        .def("get_string", &EMAN::Region::get_string)
    ;

    class_< EMAN::Pixel >("Pixel", init< const EMAN::Pixel& >())
        .def(init< int, int, int, float >())
        .def_readwrite("x", &EMAN::Pixel::x)
        .def_readwrite("y", &EMAN::Pixel::y)
        .def_readwrite("z", &EMAN::Pixel::z)
        .def_readwrite("value", &EMAN::Pixel::value)
        .def("get_point", &EMAN::Pixel::get_point)
        .def("get_value", &EMAN::Pixel::get_value)
        .def( self < self )
    ;

}

