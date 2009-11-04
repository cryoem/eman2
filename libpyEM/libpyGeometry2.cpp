/*
 * Copyright (c) 2000-2006 Baylor College of Medicine
 *
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * */

#ifdef _WIN32
	#pragma warning(disable:4819)
#endif	//_WIN32

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
	class_< EMAN::Region >("Region",
			"Region defines a 2D or 3D rectangular region specified by its\n"
			"origin coordinates and all edges' sizes. The coordinates and"
			"edge sizes can be integer or floating numbers.",
			init<  >())
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

	class_< EMAN::Pixel >("Pixel",
			"Pixel describes a 3D pixel's coordinates and its intensity value.",
			init< const EMAN::Pixel& >())
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

