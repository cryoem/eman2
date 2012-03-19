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
#include <boost/shared_ptr.hpp>

// Includes ====================================================================
#include <typeconverter.h>
#include <xydata.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
BOOST_PYTHON_MODULE(libpyTypeConverter2)
{
    class_< EMAN::EMNumPy >("EMNumPy", init<  >())
        .def(init< const EMAN::EMNumPy& >())
        .def("em2numpy", &EMAN::EMNumPy::em2numpy)
        .def("numpy2em", &EMAN::EMNumPy::numpy2em, return_value_policy< manage_new_object >())
        .staticmethod("em2numpy")
        .staticmethod("numpy2em")
    ;


	import_array();
	python::numeric::array::set_module_and_type("numpy", "ndarray");


	EMAN::vector_to_python<int>();
	EMAN::vector_to_python<long>();
	EMAN::vector_to_python<float>();
	EMAN::vector_to_python<double>();
//	EMAN::vector_to_python<EMAN::Transform3D>();
	EMAN::vector_to_python<EMAN::Transform>();
	EMAN::vector_to_python<EMAN::Ctf*>();
	EMAN::vector_to_python<EMAN::EMAN1Ctf>();
	EMAN::vector_to_python<EMAN::EMAN2Ctf>();
	EMAN::vector_to_python<std::string>();
	EMAN::vector_to_python<EMAN::EMData*>();
	EMAN::vector_to_python< boost::shared_ptr<EMAN::EMData> >();
	EMAN::vector_to_python<EMAN::Pixel>();
	EMAN::vector_to_python<EMAN::EMObject>();
	EMAN::vector_to_python<EMAN::Vec3f>();
	EMAN::vector_to_python<EMAN::Vec3i>();
	EMAN::vector_to_python<EMAN::IntPoint>();
	EMAN::vector_to_python< std::vector<EMAN::Vec3f> >();
	EMAN::vector_to_python<EMAN::Dict>();
	EMAN::vector_from_python<int>();
	EMAN::vector_from_python<long>();
	EMAN::vector_from_python<float>();
//	EMAN::vector_from_python<EMAN::Transform3D>();
	EMAN::vector_from_python<EMAN::Ctf*>();
	EMAN::vector_from_python<EMAN::Transform>();
	//EMAN::vector_from_python<EMAN::Transform*>();
	EMAN::vector_from_python<EMAN::EMAN1Ctf>();
	EMAN::vector_from_python<EMAN::EMAN2Ctf>();
	EMAN::vector_from_python<std::string>();
	EMAN::vector_from_python<EMAN::EMData*>();
	EMAN::vector_from_python<EMAN::Pixel>();
	EMAN::vector_from_python<EMAN::EMObject>();
	EMAN::vector_from_python<EMAN::Vec3f>();
	EMAN::vector_from_python<std::vector<float> >();
	EMAN::map_to_python_2<unsigned int, unsigned int>();
	EMAN::map_to_python<int>();
	EMAN::map_to_python<long>();
	EMAN::map_to_python<float>();
	EMAN::map_to_python<double>();
	EMAN::map_to_python<std::string>();
	EMAN::map_to_python<vector<string> >();

	EMAN::map_from_python<int>();
	EMAN::map_from_python<long>();
	EMAN::map_from_python<float>();
	EMAN::map_from_python<std::string>();
	EMAN::map_from_python<vector<string> >();

    EMAN::EMObject_to_python();
	EMAN::Dict_to_python();
	EMAN::Dict_from_python();

	EMAN::tuple3_to_python<EMAN::IntPoint>();
	EMAN::tuple3_to_python<EMAN::FloatPoint>();

	EMAN::tuple3_to_python<EMAN::IntSize>();
	EMAN::tuple3_to_python<EMAN::FloatSize>();

	EMAN::tuple3_from_python<EMAN::IntPoint, int>();
	EMAN::tuple3_from_python<EMAN::FloatPoint, float>();

	EMAN::tuple3_from_python<EMAN::IntSize, int>();
	EMAN::tuple3_from_python<EMAN::FloatSize, float>();

	EMAN::tuple3_from_python<EMAN::Vec3i, int>();
	EMAN::tuple3_from_python<EMAN::Vec3f, float>();

// 	EMAN::tuple2_from_python<EMAN::Vec2i, int>();
// 	EMAN::tuple2_from_python<EMAN::Vec2f, float>();

	EMAN::emobject_array_from_python();
	EMAN::emobject_emdata_from_python();
	EMAN::emobject_string_from_python();
	EMAN::emobject_xydata_from_python();
// 	EMAN::emobject_transform3d_from_python();
	EMAN::emobject_transform_from_python();
	EMAN::emobject_ctf_from_python();
	EMAN::emobject_eman1ctf_from_python();
	EMAN::emobject_eman2ctf_from_python();
	EMAN::emobject_null_from_python();

	register_ptr_to_python< boost::shared_ptr<EMAN::EMData> >();

	implicitly_convertible<int, EMAN::EMObject>();
	//implicitly_convertible<float, EMAN::EMObject>();
	implicitly_convertible<double, EMAN::EMObject>();
	implicitly_convertible<const char*, EMAN::EMObject>();
	implicitly_convertible<EMAN::EMData*, EMAN::EMObject>();
	implicitly_convertible<EMAN::XYData*, EMAN::EMObject>();
// 	implicitly_convertible<EMAN::Transform3D*, EMAN::EMObject>();
	implicitly_convertible<EMAN::Transform*, EMAN::EMObject>();
	implicitly_convertible<EMAN::Ctf*, EMAN::EMObject>();
	implicitly_convertible<EMAN::EMAN1Ctf*, EMAN::Ctf*>();
	implicitly_convertible<EMAN::EMAN2Ctf*, EMAN::Ctf*>();

	EMAN::MArrayND_to_python<2>();
	EMAN::MArrayND_to_python<3>();
	EMAN::MCArrayND_to_python<2>();
	EMAN::MCArrayND_to_python<3>();

}

