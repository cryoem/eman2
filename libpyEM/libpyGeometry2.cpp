
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
    class_< EMAN::Region >("Region", init<  >())
        .def(init< const EMAN::Region& >())
        .def(init< int, int, int, int >())
        .def(init< int, int, int, int, int, int >())
        .def(init< float, float, float, float >())
        .def(init< float, float, float, float, float, float >())
        .def(init< double, double, double, double >())
        .def(init< double, double, double, double, double, double >())
        .def(init< const EMAN::FloatPoint&, const EMAN::FloatSize& >())
        .def_readwrite("origin", &EMAN::Region::origin)
        .def_readwrite("size", &EMAN::Region::size)
        .def("inside_region", (bool (EMAN::Region::*)() const)&EMAN::Region::inside_region)
        .def("inside_region", (bool (EMAN::Region::*)(const EMAN::FloatPoint&) const)&EMAN::Region::inside_region)
        .def("inside_region", (bool (EMAN::Region::*)(float, float) const)&EMAN::Region::inside_region)
        .def("inside_region", (bool (EMAN::Region::*)(float, float, float) const)&EMAN::Region::inside_region)
        .def("is_region_in_box", &EMAN::Region::is_region_in_box)
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
        .def( self != self )
        .def( self == self )
        .def( self < self )
    ;

}

