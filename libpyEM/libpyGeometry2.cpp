
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
		.def(init< int, int >())
		.def(init< int, int, int, int >())
		.def(init< int, int, int, int, int, int >())
		.def(init< float, float >())
		.def(init< float, float, float, float >())
		.def(init< float, float, float, float, float, float >())
		.def(init< double, double >())
		.def(init< double, double, double, double >())
		.def(init< double, double, double, double, double, double >())
		.def(init< const EMAN::FloatPoint&, const EMAN::FloatSize& >())
		.def(init< const EMAN::Region& >())
		.def_readwrite("origin", &EMAN::Region::origin)
		.def_readwrite("size", &EMAN::Region::size)
		.def("inside_region", (bool (EMAN::Region::*)() const)&EMAN::Region::inside_region)
		.def("inside_region", (bool (EMAN::Region::*)(const EMAN::FloatPoint&) const)&EMAN::Region::inside_region)
		.def("inside_region", (bool (EMAN::Region::*)(float) const)&EMAN::Region::inside_region)
		.def("inside_region", (bool (EMAN::Region::*)(float, float) const)&EMAN::Region::inside_region)
		.def("inside_region", (bool (EMAN::Region::*)(float, float, float) const)&EMAN::Region::inside_region)
		.def("is_region_in_box", &EMAN::Region::is_region_in_box)
		.def("get_ndim", &EMAN::Region::get_ndim)
		.def("get_string", &EMAN::Region::get_string)
		.def("get_width", &EMAN::Region::get_width)
		.def("get_height", &EMAN::Region::get_height)
		.def("get_depth", &EMAN::Region::get_depth)
		.def("set_width", &EMAN::Region::set_width)
		.def("set_height", &EMAN::Region::set_height)
		.def("set_depth", &EMAN::Region::set_depth)
		.def("x_origin", &EMAN::Region::x_origin)
		.def("y_origin", &EMAN::Region::y_origin)
		.def("z_origin", &EMAN::Region::z_origin)
		.def("get_origin", &EMAN::Region::get_origin)
		.def("set_origin", &EMAN::Region::set_origin)
		.def("get_size", &EMAN::Region::get_size)
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
		.def( self != self )
		.def( self == self )
	;

}

