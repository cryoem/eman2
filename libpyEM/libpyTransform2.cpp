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
#include <transform.h>
#include <symmetry.h>
#include <emdata.h>
#include <emdata_pickle.h>
#include <quaternion.h>
#include <vec3.h>

#include <string>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct EMAN_Symmetry3D_Wrapper : public EMAN::Symmetry3D
{
    EMAN_Symmetry3D_Wrapper(PyObject* py_self_):
			EMAN::Symmetry3D(), py_self(py_self_) {}

	int get_max_csym() const {
		return call_method< int >(py_self, "get_max_csym");
	}

	int get_nsym() const {
		return call_method< int >(py_self, "get_nsym");
	}

	EMAN::Transform get_sym(const int n) const {
		return call_method< EMAN::Transform >(py_self, "get_sym", n);
	}

	EMAN::Transform get_sym_proj(const string s) const {
		return call_method< EMAN::Transform >(py_self, "get_sym_proj", s);
	}

	EMAN::Dict get_delimiters(const bool b) const {
		return call_method< EMAN::Dict >(py_self, "get_delimiters",b);
	}

	std::string get_name() const {
		return call_method< std::string >(py_self, "get_name");
	}
	std::string get_desc() const {
		return call_method< std::string >(py_self, "get_desc");
	}

	EMAN::TypeDict get_param_types() const {
		return call_method< EMAN::TypeDict >(py_self, "get_param_types");
	}

	std::vector<EMAN::Vec3f > get_asym_unit_points(bool b) const {
		return call_method< std::vector<EMAN::Vec3f > >(py_self, "get_asym_unit_points", b);
	}

	std::vector<std::vector<EMAN::Vec3f > > get_asym_unit_triangles(bool b) const {
		return call_method< std::vector<std::vector<EMAN::Vec3f > > >(py_self, "get_asym_unit_triangles", b);
	}

	void insert_params(const EMAN::Dict& d) {
		return call_method< void >(py_self, "insert_params", d);
	}

	bool is_in_asym_unit(const float& altitude, const float& azimuth, const bool inc_mirror) const {
		return call_method< bool >(py_self, "is_in_asymm_init",altitude,azimuth,inc_mirror);
	}

	PyObject* py_self;
};

struct EMAN_OrientationGenerator_Wrapper : public EMAN::OrientationGenerator
{
    EMAN_OrientationGenerator_Wrapper(PyObject* py_self_):
			EMAN::OrientationGenerator(), py_self(py_self_) {}

	std::vector<EMAN::Transform> gen_orientations(const EMAN::Symmetry3D* const sym) const {
		return call_method< std::vector<EMAN::Transform> >(py_self, "gen_orientations", sym);
	}

	std::string get_name() const {
		return call_method< std::string >(py_self, "get_name");
	}
	std::string get_desc() const {
		return call_method< std::string >(py_self, "get_desc");
	}

	EMAN::TypeDict get_param_types() const {
		return call_method< EMAN::TypeDict >(py_self, "get_param_types");
	}

	int get_orientations_tally(const EMAN::Symmetry3D* const sym, const float& delta ) const {
		return call_method< int >(py_self, "get_orientations_tally",sym,delta);
	}

	PyObject* py_self;
};

//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_Transform3D_get_rotation_overloads_0_1, get_rotation, 0, 1)
//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_Transform3D_get_pretrans_overloads_0_1, get_pretrans, 0, 1)
//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_Transform3D_get_posttrans_overloads_0_1, get_posttrans, 0, 1)
//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_Transform3D_set_pretrans_overloads_2_3, set_pretrans, 2, 3)
//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_Transform3D_set_posttrans_overloads_2_3, set_posttrans, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_Transform_settrans_overloads_2_3, set_trans, 2, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_Transform_translate_overloads_2_3, translate,2, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_Transform_translate_newbasis_overloads_2_3, translate_newBasis,3, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_Transform_get_rotation_overloads_0_1, get_rotation, 0, 1)

}// namespace


// Module ======================================================================

BOOST_PYTHON_MODULE(libpyTransform2)
{
	def("dump_orientgens", &EMAN::dump_orientgens);
	def("dump_orientgens_list", &EMAN::dump_orientgens_list);

	def("dump_symmetries", &EMAN::dump_symmetries);
	def("dump_symmetries_list", &EMAN::dump_symmetries_list);

	class_< EMAN::Symmetry3D, boost::noncopyable, EMAN_Symmetry3D_Wrapper >("Symmetry3D",
			"A base class for 3D Symmetry objects.\n"
			"Objects of this type must provide delimiters for the asymmetric unit (get_delimiters), and\n"
			"must also provide all of the rotational symmetric operations (get_sym(const int n)). They must also\n"
			"provide the total number of unique symmetric operations with get_nsym (except in helical symmetry).\n"
			"get_delimiter returns a dictionary with 'alt_max' and 'az_max' keys, which correspond to the\n"
			"encompassing azimuth and altitude angles of the asymmetric unit. These can be interpreted in a\n"
			"relatively straight forward fashion when dealing with C and D symmetries to demarcate the asymmetric\n"
			"unit, however when dealing with Platonic symmetries the asymmetric unit is not so trivial.\n"
			"see http://blake.bcm.edu/emanwiki/EMAN2/Symmetry for figures and description of what we're doing\n"
			"here, for all the symmetries, and look in the comments of the PlatonicSym classes themselves.\n"
			"It inherits from a factory base, making it amenable to incorporation in EMAN2 style factories\n",
			init<  >())
		.def("get_delimiters", pure_virtual(&EMAN::Symmetry3D::get_delimiters), args("inc_mirror"), "Every Symmetry3D object must return a dictionary containing the delimiters\nthat define its asymmetric unit (this is not strictly true in the case of the PlatonicSym class)\n \ninc_mirror - whether or not the mirror part of the asymmetric unit should be included in the consideration (default = False)\n \nreturn a dictionary containing atleast 'alt_max' and 'az_max'\n")
		.def("get_sym", pure_virtual(&EMAN::Symmetry3D::get_sym), args("n"), "Every Symmetry3D object must provide access to the full set of\nits symmetry operators via this function\n \nn - the symmetry operator number\n \nreturn a Transform object describing the symmetry operation\n")
		.def("is_in_asym_unit", pure_virtual(&EMAN::Symmetry3D::is_in_asym_unit), args("altitude", "azimuth", "inc_mirror"), "A function to be used when generating orientations over portion of the unit sphere\ndefined by parameters returned by get_delimiters. In platonic symmetry altitude and azimuth\nalone are not enough to correctly demarcate the asymmetric unit. See the get_delimiters comments.\n \naltitude - the EMAN style altitude of the 3D orientation in degrees\nazimuth - the EMAN style azimuth of the 3D orientation in degrees\ninc_mirror - whether or not to include orientations if they are in the mirror portion of the asymmetric unit\n \nreturn true or false, depending on whether or not the orientation is within the asymmetric unit\n")
		.def("get_asym_unit_triangles", pure_virtual(&EMAN::Symmetry3D::get_asym_unit_triangles), args("get_asym_unit_triangles"), "Get triangles that precisely occlude the projection area of the default asymmetric unit. This will be used\nfor collision detection in Symmetry3D::reduce\n \ninc_mirror - whether to include the mirror portion of the asymmetric unit\n")
		.def("get_nsym",pure_virtual(&EMAN::Symmetry3D::get_nsym), "Gets the total number of unique roational symmetry operations associated with this symmetry\nFor C symmetry, this is simply nsym\n \nreturn the degree of of cyclic symmetry (nsym)\n")
		.def("get_name",pure_virtual(&EMAN::Symmetry3D::get_name), "Return DSym::NAME\n \nreturn the unique name of this class\n")
		.def("is_platonic_sym",pure_virtual(&EMAN::Symmetry3D::is_platonic_sym), "Determines whether or not this Symmetry3D is the platonic type - returns true\n \nreturn true - indicating that this is a platonic symmetry object\n")
		.def("get_az_alignment_offset", pure_virtual(&EMAN::Symmetry3D::get_az_alignment_offset), "Get the azimuth alignment offset required to ensure that orientations align correctly\nwith symmetric axes of the tetrahedron. This offset is directly related to the way\nthe symmetric operations are generated by get_sym. All orientations generated as a\nresult of using the delimiters supplied by this class should by offset by this azimuth\nto ensure proper alignment with tetrahedral objects in EMAN2\n")
		.def("is_h_sym", pure_virtual(&EMAN::Symmetry3D::is_h_sym), "A function that is used to determine if this is a Helical symmetry object\nThis function is only virtually overidden by the HSym symmetry, which returns true, not false\n \nreturn false - indicating that this is not a helical symmetry object\n")
		.def("get_params", &EMAN::Symmetry3D::get_params, "get a copy of the parameters of this class\n \nreturn a copy of the parameters of this class")
		.def("gen_orientations", &EMAN::Symmetry3D::gen_orientations, args("generatorname", "parms"), "Ask the Symmetry3D object to generate a set of orientations in its asymmetric unit\nusing an OrientationGenerator constructed from the given parameters (using a Factory).\nThis is reminiscent of the strategy design pattern\n \ngeneratorname - the string name of the OrientationGenerator, as accessed for the OrientationGenerator factory\nparms - the parameters handed to OrientationGenerator::set_params after initial construction\n \nreturn a set of orientations in the unit sphere\n")
		.def("get_asym_unit_points", pure_virtual(&EMAN::Symmetry3D::get_asym_unit_points), args("inc_mirror"), "to demarcate the asymmetric unit. The last should may be connected to the first.\n \ninc_mirror - whether or not to include the mirror portion of the asymmetric unit (default=False)\n \nreturn a cyclic set of points which can be connected using great arcs on the unit sphere\n")
		.def("insert_params", pure_virtual(&EMAN::Symmetry3D::insert_params), args("new_params"), "Insert parameters. Previously present parameters are replaced, new ones are inserted.\n \nnew_params - the parameters to insert\n")
		.def("reduce", &EMAN::Symmetry3D::reduce, args("t3d", "n"), "A function that will reduce an orientation, as characterized by Euler anges, into a specific asymmetric unit.\nDefault behavior is to map the given orientation into the default asymmetric unit of the symmetry (n=0). This\nis a concrete implementation that works for all symmetries, relying on a concrete instance\nof the get_asym_unit_triangles function\n \nt3d - a Transform characterizing an orientation\nn - the number of the asymmetric unit you wish to map the given orientation into. There is a strong relationship between n and to Symmetry3D::get_sym (default = 0)\n \nreturn the orientation the specified asymmetric unit (by default this is the default asymmetric unit of the symmetry)")
		.def("in_which_asym_unit", &EMAN::Symmetry3D::in_which_asym_unit, args("t3d"), "A function that will determine in which asymmetric unit a given orientation resides\nThe asymmetric unit 'number' will depend entirely on the order in which different symmetry\noperations are return by the Symmetry3D::get_sym function\n \nt3d - a Transform characterizing an orientation\n \nreturn the asymmetric unit number the the orientation is in\n")
		.def("point_in_which_asym_unit", &EMAN::Symmetry3D::point_in_which_asym_unit, args("v"), "A function that will determine in which asymmetric unit a given vector resides\nThe asymmetric unit 'number' will depend entirely on the order in which different\nsymmetry operations are return by the Symmetry3D::get_sym function\n \nv a Vec3f characterizing a point\n \nreturn the asymmetric unit number the the orientation is in\n")
		.def("get_touching_au_transforms",&EMAN::Symmetry3D::get_touching_au_transforms, args("inc_mirror"), "Gets a vector of Transform objects that define the set of asymmetric units that touch the default\nasymmetric unit. The 'default asymmetric unit' is defined by the results of Symmetry3d::get_asym_unit_points\nand is sensitive to whether or not you want to include the mirror part of the asymmetric unit.\nThis function is useful when used in conjunction with Symmetry3D::reduce, and particularly when finding\nthe angular deviation of particles through different stages of iterative Single Particle Reconstruction\nThis function could be expanded to work for an asymmetric unit number supplied by the user.\n \ninc_mirror - whether or not to include the mirror portion of the asymmetric unit\n \nreturn a vector of Transform objects that map the default asymmetric unit to the neighboring asymmetric unit\n")
		.def("get_syms", &EMAN::Symmetry3D::get_syms, "")
		.def("get_symmetries", &EMAN::Symmetry3D::get_symmetries, "")
		.staticmethod("get_symmetries")
		;


	class_< EMAN::Factory<EMAN::Symmetry3D>, boost::noncopyable >("Symmetries", no_init)
		.def("get", (EMAN::Symmetry3D* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&))&EMAN::Factory<EMAN::Symmetry3D>::get, return_value_policy< manage_new_object >())
		.def("get", (EMAN::Symmetry3D* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&, const EMAN::Dict&))&EMAN::Factory<EMAN::Symmetry3D>::get, return_value_policy< manage_new_object >())
		.def("get_list", &EMAN::Factory<EMAN::Symmetry3D>::get_list)
		.staticmethod("get_list")
		.staticmethod("get")
		;

	class_< EMAN::OrientationGenerator, boost::noncopyable, EMAN_OrientationGenerator_Wrapper >("__OrientationGenerator", init<  >())
			.def("gen_orientations", pure_virtual(&EMAN::OrientationGenerator::gen_orientations))
			.def("get_orientations_tally", pure_virtual(&EMAN::OrientationGenerator::get_orientations_tally))
		;

	class_< EMAN::Factory<EMAN::OrientationGenerator>, boost::noncopyable >("OrientGens", no_init)
		.def("get", (EMAN::OrientationGenerator* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&))&EMAN::Factory<EMAN::OrientationGenerator>::get, return_value_policy< manage_new_object >())
		.def("get", (EMAN::OrientationGenerator* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&, const EMAN::Dict&))&EMAN::Factory<EMAN::OrientationGenerator>::get, return_value_policy< manage_new_object >())
		.def("get_list", &EMAN::Factory<EMAN::OrientationGenerator>::get_list)
		.staticmethod("get_list")
		.staticmethod("get")
		;
	typedef void (EMAN::Vec4f::*vec4f_set_value_at_float)(int, const float&);
	typedef void (EMAN::Vec4f::*vec4f_set_value_float)(const float&, const float&, const float&, const float&);
	typedef float (EMAN::Vec4f::*vec4f_at_float)(int) const;	
	class_< EMAN::Vec4f >("Vec4f", "typedef Vec4<float> Vec4f;", init<  >())
		.def(init< float, float, float, float >())
		.def(init< const std::vector<float,std::allocator<float> >& >())
		.def(init< const EMAN::Vec4i& >())
		.def(init< const EMAN::Vec4f& >())
		.def("normalize", &EMAN::Vec4f::normalize)
		.def("length", &EMAN::Vec4f::length)
		.def("set_value", (void (EMAN::Vec4f::*)(const std::vector<float,std::allocator<float> >&) )&EMAN::Vec4f::set_value)
		.def("set_value_at", vec4f_set_value_at_float(&EMAN::Vec4f::set_value_at))
		.def("set_value", (void (EMAN::Vec4f::*)(float, float, float, float) )vec4f_set_value_float(&EMAN::Vec4f::set_value))
	.def("at", vec4f_at_float(&EMAN::Vec4f::at))
	.def("__getitem__",vec4f_at_float(&EMAN::Vec4f::at))
		.def("__setitem__",vec4f_set_value_at_float(&EMAN::Vec4f::set_value_at))
		.def("__len__",&EMAN::Vec4f::number_of_element)
		.def("__iter__", range(&EMAN::Vec4f::begin, &EMAN::Vec4f::end))
    ;
	typedef float (EMAN::Vec3f::*dot_float)(const EMAN::Vec3f&) const;
	typedef EMAN::Vec3f (EMAN::Vec3f::*cross_float)(const EMAN::Vec3f&) const;
	typedef void (EMAN::Vec3f::*set_value_at_float)(int, const float&);
	typedef void (EMAN::Vec3f::*set_value_float)(const float&, const float&,const float&);
	typedef float (EMAN::Vec3f::*at_float)(int) const;
	
    class_< EMAN::Vec3f >("Vec3f", "typedef Vec3<float> Vec3f;", init<  >())
    	.def_pickle(Vec3f_pickle_suite())
        .def(init< float, float, float >())
        .def(init< const std::vector<float,std::allocator<float> >& >())
        .def(init< const EMAN::Vec3i& >())
        .def(init< const EMAN::Vec3f& >())
        .def("normalize", &EMAN::Vec3f::normalize)
        .def("length", &EMAN::Vec3f::length)
        .def("dot", dot_float(&EMAN::Vec3f::dot))
		.def("cross", cross_float(&EMAN::Vec3f::cross))
        .def("as_list", &EMAN::Vec3f::as_list)
        .def("set_value", (void (EMAN::Vec3f::*)(const std::vector<float,std::allocator<float> >&) )&EMAN::Vec3f::set_value)
		.def("set_value_at", set_value_at_float(&EMAN::Vec3f::set_value_at))
		.def("set_value", (void (EMAN::Vec3f::*)(float, float, float) )set_value_float(&EMAN::Vec3f::set_value))
        .def("at", at_float(&EMAN::Vec3f::at))
        .def("__getitem__",at_float(&EMAN::Vec3f::at))
		.def("__setitem__",set_value_at_float(&EMAN::Vec3f::set_value_at))
		.def("__len__",&EMAN::Vec3f::number_of_element)
		.def("__iter__", range(&EMAN::Vec3f::begin, &EMAN::Vec3f::end))
        .def( self + other< EMAN::Vec3i >() )
        .def( other< EMAN::Vec3i >() - self )
        .def( self - other< EMAN::Vec3i >() )
        .def( -self )
        .def( self * other< EMAN::Vec3i >() )
        .def( other< EMAN::Vec3i >() * self )
//         .def( other< float >() + self )
        .def( self + other< float >() )
//         .def( other< float >() - self )
        .def( self - other< float >() )
        .def( self * other< float >() )
        .def( other< float >() * self )
        .def( self != self )
        .def( self == self )
        .def( self / other< float >() )
        .def( self - self )
        .def( self + self )
        .def( other< EMAN::Vec3i >() + self )
        .def( self * self )
        .def( self += self )
        .def( self += other< EMAN::Vec3i >() )
        .def( self += other< float >() )
        .def( self -= self )
        .def( self -= other< EMAN::Vec3i >() )
        .def( self -= other< float >() )
        .def( self *= other< float >() )
        .def( self /= other< float >() )
    ;

	typedef int (EMAN::Vec3i::*dot_int)(const EMAN::Vec3i&) const;
	typedef EMAN::Vec3i (EMAN::Vec3i::*cross_int)(const EMAN::Vec3i&) const;
	typedef void (EMAN::Vec3i::*set_value_at_int)(int, const int&);
	typedef void (EMAN::Vec3i::*set_value_int)(const int&, const int&,const int&);
	typedef int (EMAN::Vec3i::*at_int)(int) const;
    class_< EMAN::Vec3i >("Vec3i", "typedef Vec3<int> Vec3i;", init<  >())
        .def(init< int, int,  int  >())
        .def(init< const std::vector<int,std::allocator<int> >& >())
        .def(init< const EMAN::Vec3i& >())
        .def("normalize", &EMAN::Vec3i::normalize)
        .def("length", &EMAN::Vec3i::length)
        .def("dot", dot_int(&EMAN::Vec3i::dot))
        .def("cross", cross_int(&EMAN::Vec3i::cross))
        .def("as_list", &EMAN::Vec3i::as_list)
        .def("set_value", (void (EMAN::Vec3i::*)(const std::vector<int,std::allocator<int> >&) )&EMAN::Vec3i::set_value)
		.def("set_value_at",set_value_at_int( &EMAN::Vec3i::set_value_at))
        .def("set_value", (void (EMAN::Vec3i::*)(int, int, int) ) set_value_int(&EMAN::Vec3i::set_value))
		.def("at", at_int(&EMAN::Vec3i::at))
		.def("__getitem__",at_int(&EMAN::Vec3i::at))
		.def("__setitem__",set_value_at_int(&EMAN::Vec3i::set_value_at))
		.def("__len__",&EMAN::Vec3i::number_of_element)
		.def("__iter__", range(&EMAN::Vec3i::begin, &EMAN::Vec3i::end))
        .def( other< EMAN::Vec3f >() + self )
        .def( self - other< EMAN::Vec3f >() )
        .def( other< EMAN::Vec3f >() - self )
        .def( other< EMAN::Vec3f >() * self )
        .def( self * other< EMAN::Vec3f >() )
        .def( self + other< int >() )
//         .def( other< int >() + self )
        .def( self + self )
//         .def( other< int >() - self )
        .def( self * self )
        .def( self - self )
        .def( self - other< int >() )
        .def( self * other< int >() )
        .def( other< int >() * self )
        .def( self == self )
        .def( self / other< int >() )
        .def( self != self )
        .def( self + other< EMAN::Vec3f >() )
        .def( self += self )
        .def( self += other< int >() )
        .def( self -= self )
        .def( self -= other< int >() )
        .def( self *= other< int >() )
        .def( self /= other< int >() )
    ;

	typedef float (EMAN::Vec2f::*vec2f_dot_float)(const EMAN::Vec2f&) const;
	typedef void (EMAN::Vec2f::*vec2f_set_value_at_float)(int, const float&);
	typedef void (EMAN::Vec2f::*vec2f_set_value_float)(const float&, const float&);
	typedef float (EMAN::Vec2f::*vec2f_at_float)(int) const;
	class_< EMAN::Vec2f >("Vec2f", "typedef Vec2<float> Vec2f;", init<  >())
		.def(init< float, float >())
		.def(init< const std::vector<float,std::allocator<float> >& >())
		.def(init< const EMAN::Vec2i& >())
		.def(init< const EMAN::Vec2f& >())
		.def("normalize", &EMAN::Vec2f::normalize)
		.def("length", &EMAN::Vec2f::length)
		.def("dot", vec2f_dot_float(&EMAN::Vec2f::dot))
		.def("as_list", &EMAN::Vec2f::as_list)
		.def("set_value", (void (EMAN::Vec2f::*)(const std::vector<float,std::allocator<float> >&) )&EMAN::Vec2f::set_value)
		.def("set_value_at", vec2f_set_value_at_float(&EMAN::Vec2f::set_value_at))
		.def("set_value", (void (EMAN::Vec2f::*)(float, float) )vec2f_set_value_float(&EMAN::Vec2f::set_value))
		.def("at", vec2f_at_float(&EMAN::Vec2f::at))
		.def("__getitem__",vec2f_at_float(&EMAN::Vec2f::at))
		.def("__setitem__",vec2f_set_value_at_float(&EMAN::Vec2f::set_value_at))
		.def("__len__",&EMAN::Vec2f::number_of_element)
		.def("__iter__", range(&EMAN::Vec2f::begin, &EMAN::Vec2f::end))
		.def( self + other< EMAN::Vec2i >() )
		.def( other< EMAN::Vec2i >() - self )
		.def( self - other< EMAN::Vec2i >() )
		.def( -self )
		.def( self * other< EMAN::Vec2i >() )
		.def( other< EMAN::Vec2i >() * self )
		.def( self + other< float >() )
		.def( self - other< float >() )
		.def( self * other< float >() )
		.def( other< float >() * self )
		.def( self != self )
		.def( self == self )
		.def( self / other< float >() )
		.def( self - self )
		.def( self + self )
		.def( other< EMAN::Vec2i >() + self )
		.def( self * self )
		.def( self += self )
		.def( self += other< EMAN::Vec2i >() )
		.def( self += other< float >() )
		.def( self -= self )
		.def( self -= other< EMAN::Vec2i >() )
		.def( self -= other< float >() )
		.def( self *= other< float >() )
		.def( self /= other< float >() )
	;

    class_< EMAN::Quaternion >("Quaternion",
    		"Quaternion is used in Rotation and Transformation to replace Euler angles.\n\n"
    		"Quaternions extend the concept of rotation in three dimensions to\n"
    		"rotation in four dimensions. This avoids the problem of \"gimbal-lock\"\n"
    		"and allows for the implementation of smooth and continuous rotation.\n\n"
    		"Euler angles have the disadvantage of being\n"
    		"susceptible to \"Gimbal lock\" where attempts to rotate an\n"
    		"object fail due to the order in which the rotations are performed.\n\n"
    		"Quaternions are a solution to this problem. Instead of rotating an\n"
    		"object through a series of successive rotations, a quaternion allows\n"
    		"the programmer to rotate an object through a single arbitary rotation axis.\n\n"
    		"Because the rotation axis is specifed as a unit direction vector,\n"
    		"it may be calculated through vector mathematics or from spherical\n"
    		"coordinates ie (longitude/latitude).\n\n"
    		"Quaternions offer another advantage in that they be interpolated.\n"
    		"This allows for smooth and predictable rotation effects.",
    		init<  >())
        .def(init< const EMAN::Quaternion& >())
        .def(init< float, float, float, float >())
        .def(init< float, const EMAN::Vec3f& >())
        .def(init< const EMAN::Vec3f&, float >())
        .def(init< const std::vector<float,std::allocator<float> >& >())
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
        .def( self + self )
        .def( self - self )
        .def( self * self )
        .def( other< float >() * self )
        .def( self * other< float >() )
        .def( self != self )
        .def( self == self )
        .def( self / self )
        .def( self += self )
        .def( self -= self )
        .def( self *= self )
        .def( self *= other< float >() )
        .def( self /= self )
        .def( self /= other< float >() )
    ;

//	scope* EMAN_Transform3D_scope = new scope(
//	class_< EMAN::Transform3D >("Transform3D",
//			"These are  a collection of transformation tools: rotation, translation,\n"
//			"and construction of symmetric objects\n"
//			"Transform defines a transformation, which can be rotation,\n"
//			"translation, scale, and their combinations.\n\n"
//			"Internally a transformation is stored in a 4x4 matrix.\n"
//			"        a b c d\n"
//			"        e f g h           R        v\n"
//			" M=     j k m n    =      vpre     1    , where R is 3by3, v is 3by1\n"
//			"        p q r 1\n"
//			"The standard computer graphics convention is identical to ours after setting vpre to\n"
//			"zero and can be found in many references including Frackowiak et al; Human Brain Function\n\n"
//			"The left-top 3x3 submatrix\n\n"
//			"       a b c\n"
//			"  R =  e f g\n"
//			"       j k m\n\n"
//			"provides rotation, scaling and skewing (not yet implimented).\n\n"
//			"The cumulative translation is stored in (d, h, n).\n"
//			"We put the post-translation into (p, q, r), since it is convenient\n"
//			"to carry along at times. When matrices are multiplied or printed, these\n"
//			"are hidden to the user. They can only be found by applying the post_translation\n"
//			"method, and these elements are non-zero. Otherwise the post-translation method returns\n"
//			"the cumulative translationmlb\n\n"
//			"If rotations need to be found\n"
//			"around alternate origins, then brief calculations need to be performed\n"
//			"Pre and Post Translations should be kept as separate vectors\n\n"
//			"a matrix  R is called orthogonal if\n"
//			"          R * transpose(R) = 1.\n"
//			"All Real Orthogonal Matrices have eigenvalues with unit modulus and determinant"
//			"therefore equal to  \"\\pm 1\"",
//			init<  >())
//		.def(init< const EMAN::Transform3D& >())
//		.def(init< const float&, const float&, const float& >())
//// 		.def(init< const EMAN::Dict&, const string&, optional<const EMAN::Transform3D::EulerType> >())
//		.def(init< const float&, const float&, const float&, const EMAN::Vec3f& >())
//		.def(init< EMAN::Transform3D::EulerType, const float&, const float&, const float& >())
//		.def(init< EMAN::Transform3D::EulerType, const float&, const float&, const float&, const float&>())
//		.def(init< EMAN::Transform3D::EulerType, const EMAN::Dict& >())
//		.def(init< const EMAN::Vec3f&, const float&, const float&, const float&, const EMAN::Vec3f& >())
//		.def(init< const float&, const float&, const float&, const float&, const float&, const float&, const float&, const float&, const float& >())
//		.add_static_property("ERR_LIMIT", make_getter(EMAN::Transform3D::ERR_LIMIT))
//		.def("set_posttrans", (void (EMAN::Transform3D::*)(const float&, const float&, const float&) )&EMAN::Transform3D::set_posttrans,EMAN_Transform3D_set_posttrans_overloads_2_3())
//		.def("set_posttrans", (void (EMAN::Transform3D::*)(const EMAN::Vec3f&) )&EMAN::Transform3D::set_posttrans)
//		.def("set_posttrans", (void (EMAN::Transform3D::*)(const EMAN::Vec2f&) )&EMAN::Transform3D::set_posttrans)
//		.def("apply_scale", &EMAN::Transform3D::apply_scale)
//		.def("set_scale", &EMAN::Transform3D::set_scale)
//		.def("orthogonalize", &EMAN::Transform3D::orthogonalize)
//		.def("transpose", &EMAN::Transform3D::transpose)
//		.def("set_rotation", (void (EMAN::Transform3D::*)(const float&, const float&, const float&) )&EMAN::Transform3D::set_rotation)
//		.def("set_rotation", (void (EMAN::Transform3D::*)(EMAN::Transform3D::EulerType, const float&, const float&, const float&) )&EMAN::Transform3D::set_rotation)
//		.def("set_rotation", (void (EMAN::Transform3D::*)(EMAN::Transform3D::EulerType, const float&, const float&, const float&, const float&) )&EMAN::Transform3D::set_rotation)
//		.def("set_rotation", (void (EMAN::Transform3D::*)(const float&, const float&, const float&, const float&, const float&, const float&, const float&, const float&, const float&) )&EMAN::Transform3D::set_rotation)
//		.def("set_rotation", (void (EMAN::Transform3D::*)(EMAN::Transform3D::EulerType, const EMAN::Dict&) )&EMAN::Transform3D::set_rotation)
//		.def("set_rotation", (void (EMAN::Transform3D::*)(const EMAN::Vec3f&, const EMAN::Vec3f&, const EMAN::Vec3f&, const EMAN::Vec3f&) )&EMAN::Transform3D::set_rotation)
//		.def("get_post_x_mirror",&EMAN::Transform3D::get_post_x_mirror)
//		.def("set_post_x_mirror",&EMAN::Transform3D::set_post_x_mirror)
//		.def("get_mag", &EMAN::Transform3D::get_mag)
//		.def("get_finger", &EMAN::Transform3D::get_finger)
//		.def("get_pretrans", &EMAN::Transform3D::get_pretrans, EMAN_Transform3D_get_pretrans_overloads_0_1())
//		.def("get_posttrans", &EMAN::Transform3D::get_posttrans, EMAN_Transform3D_get_posttrans_overloads_0_1())
//		.def("get_total_posttrans", &EMAN::Transform3D::get_total_posttrans)
//		.def("get_total_pretrans", &EMAN::Transform3D::get_total_pretrans)
//		.def("get_center", &EMAN::Transform3D::get_center)
//		.def("get_matrix3_col", &EMAN::Transform3D::get_matrix3_col)
//		.def("get_matrix3_row", &EMAN::Transform3D::get_matrix3_row)
//		.def("transform", &EMAN::Transform3D::transform)
//		.def("rotate", &EMAN::Transform3D::rotate)
//		.def("inverse", &EMAN::Transform3D::inverse)
//		.def("inverseUsingAngs", &EMAN::Transform3D::inverseUsingAngs)
//		.def("get_rotation", &EMAN::Transform3D::get_rotation, EMAN_Transform3D_get_rotation_overloads_0_1())
//		.def("printme", &EMAN::Transform3D::printme)
//		.def("at", &EMAN::Transform3D::at)
//		.def("get_nsym", &EMAN::Transform3D::get_nsym)
//		.def("get_sym", &EMAN::Transform3D::get_sym)
//		.def("set_center", &EMAN::Transform3D::set_center)
//		.def("set_pretrans", (void (EMAN::Transform3D::*)(const float&,const float&,const float&))&EMAN::Transform3D::set_pretrans, EMAN_Transform3D_set_pretrans_overloads_2_3())
//		.def("set_pretrans", (void (EMAN::Transform3D::*)(const EMAN::Vec3f&))&EMAN::Transform3D::set_pretrans)
//		.def("set_pretrans", (void (EMAN::Transform3D::*)(const EMAN::Vec2f&))&EMAN::Transform3D::set_pretrans)
//		.def("get_scale", &EMAN::Transform3D::get_scale)
//		.def("to_identity", &EMAN::Transform3D::to_identity)
//		.def("is_identity", &EMAN::Transform3D::is_identity)
//		.def("angles2tfvec", &EMAN::Transform3D::angles2tfvec)
//		.staticmethod("angles2tfvec")
//		.staticmethod("get_nsym")
//		.def( other< EMAN::Vec3f >() * self )
//		.def( self * self )
//		.def( self * other< EMAN::Vec3f >() )
//		.def( self * other< EMAN::Vec2f >() )
//// 		.def( self + self )
//// 		.def( self - self )
//    );
//
//    enum_< EMAN::Transform3D::EulerType >("EulerType")
//        .value("MATRIX", EMAN::Transform3D::MATRIX)
//        .value("UNKNOWN", EMAN::Transform3D::UNKNOWN)
//        .value("XYZ", EMAN::Transform3D::XYZ)
//        .value("IMAGIC", EMAN::Transform3D::IMAGIC)
//        .value("SPIDER", EMAN::Transform3D::SPIDER)
//        .value("QUATERNION", EMAN::Transform3D::QUATERNION)
//        .value("SGIROT", EMAN::Transform3D::SGIROT)
//        .value("MRC", EMAN::Transform3D::MRC)
//        .value("SPIN", EMAN::Transform3D::SPIN)
//        .value("EMAN", EMAN::Transform3D::EMAN)
//    ;
//
//	delete EMAN_Transform3D_scope;

	class_< EMAN::Transform >("Transform",
			"A Transform object is a somewhat specialized object designed specifically for EMAN2/Sparx storage of\n"
			"alignment parameters and euler orientations.\n"
			"It's designed to store four transformations in a specific order, namely\n"
			"Transform = MTSR\n"
			"Where M is a mirroring operation (about the x-axis) or the identity\n"
			"T is a Translation matrix\n"
			"S is a uniform scaling matrix\n"
			"R is a rotation matrix\n"
			"This means you can call set_scale, set_trans, set_rotation in any order but still have the operations arranged\n"
			"internally in the order of MTSR. This is somewhat restrictive, for example in the context of how OpenGL handles\n"
			"transformations, but in practice is nicely suited to the situations that arise in EMAN2 - namely, alignment and\n"
			"projection orientation characterization.\n\n"
			"Note that you can fool the Transform object into storing any matrix by using the constructors that take array arguments.\n"
			"This can useful, for example, for shearing your image.\n\n"
			"See http://blake.bcm.tmc.edu/emanwiki/Eman2TransformInPython for using it from Python and detailed discussion\n"
			"See test_transform.py for examples of the way it is unit tested\n"
			"See http://blake.bcm.tmc.edu/emanwiki/EMAN2/Tutorials/RotateTranslate for examples showing how to transform EMDatas with it.\n",
			init<  >())
//	class_< EMAN::Transform, std::auto_ptr<EMAN::Transform>  >("Transform", init<  >())
		.def(init< const EMAN::Transform& >())
		.def(init<const EMAN::Dict& >())
		.def(init<const std::vector<float>& >())
		.def_pickle(Transform_pickle_suite())
		.add_static_property("ERR_LIMIT", make_getter(EMAN::Transform::ERR_LIMIT))
		.def("set_trans", (void (EMAN::Transform::*)(const float&, const float&, const float&) )&EMAN::Transform::set_trans, EMAN_Transform_settrans_overloads_2_3(args("x", "y", "z"), "Set the post translation component\n \nx - the x translation\ny - the y translation\nz - the z translation (default = 0)\n"))
		.def("set_trans", (void (EMAN::Transform::*)(const EMAN::Vec3f&) )&EMAN::Transform::set_trans, args("v"), "Set the post translation component using a Vec3f\n \nv - the 3D translation vector\n")
		.def("set_trans", (void (EMAN::Transform::*)(const EMAN::Vec2f&) )&EMAN::Transform::set_trans, args("v"), "v - the 2D translation vector\n")
		.def("translate", (void (EMAN::Transform::*)(const float&, const float&, const float&) )&EMAN::Transform::translate, EMAN_Transform_translate_overloads_2_3(args("x", "y", "z"), "Increment the post translation component\n \nx - the x translation\ny - the y translation\nz - the z translation (default = 0)\n"))
		.def("translate", (void (EMAN::Transform::*)(const EMAN::Vec3f&) )&EMAN::Transform::translate, args("v"), "Increment the post translation component using a Vec3f\n \nv - the 3D translation vector\n")
		.def("translate", (void (EMAN::Transform::*)(const EMAN::Vec2f&) )&EMAN::Transform::translate, args("v"), "v - the 2D translation vector\n")
		.def("translate_newbasis", (void (EMAN::Transform::*)(const EMAN::Transform&, const float&, const float&, const float&) )&EMAN::Transform::translate_newBasis, EMAN_Transform_translate_newbasis_overloads_2_3(args("x", "y", "z"), "Increment the post translation component\n \nx - the x translation\ny - the y translation\nz - the z translation (default = 0)\n"))
		.def("translate_newbasis", (void (EMAN::Transform::*)(const EMAN::Transform&, const EMAN::Vec3f&) )&EMAN::Transform::translate_newBasis, args("v"), "Increment the post translation component using a Vec3f\n \nv - the 3D translation vector\n")
		.def("set_scale", &EMAN::Transform::set_scale, args("scale"), "Set the scale\n \nscale - the amount to scale by\n")
		.def("set_rotation", (void (EMAN::Transform::*)(const EMAN::Dict&) )&EMAN::Transform::set_rotation, args("rotation"), "Set a rotation using a specific Euler type and the dictionary interface\nWorks for all Euler types\n \nrotation - a dictionary containing all key-entry pair required of the associated Euler type\n")
		.def("set_rotation", (void (EMAN::Transform::*)(const EMAN::Vec3f&) )&EMAN::Transform::set_rotation, args("v"), "Determine the rotation that would transform a vector pointing in the Z direction\nso that it points in the direction of the argument vector\nAutomatically normalizes the vector\n \nv - the direction you want to solve for\n")
		.def("rotate_origin", (void (EMAN::Transform::*)(const EMAN::Transform&) )&EMAN::Transform::rotate_origin, args("rotation"), "Increment the rotation using the rotation bit of a tranfromation matrix.")
		.def("rotate_origin_newbasis", (void (EMAN::Transform::*)(const EMAN::Transform&, const float&, const float&, const float&, const float&) )&EMAN::Transform::rotate_origin_newBasis, args("rotation, omega, n1, n2, n3"), "Increment the rotation using the rotation bit of a tranfromation matrix.")
		.def("rotate", (void (EMAN::Transform::*)(const EMAN::Transform&) )&EMAN::Transform::rotate, args("rotation"), "Increment the rotation using a tranfromation matrix.")
		.def("get_hflip_transform", (EMAN::Transform (EMAN::Transform::*)() const )&EMAN::Transform::get_hflip_transform, "How do I get the transform that will yield the horizontally flipped projection?\n \nreturn the transform that will yield the horizontally flipped projection\n")
		.def("get_vflip_transform", (EMAN::Transform (EMAN::Transform::*)() const )&EMAN::Transform::get_vflip_transform, "How do I get the transform that will yield the vertically flipped projection?\n \nreturn the transform that will yield the vertically flipped projection\n")
		.def("get_mirror",&EMAN::Transform::get_mirror, "Query whether x_mirroring is occuring\n \nreturn whether x_mirroring is occuring\n")
		.def("set_mirror",&EMAN::Transform::set_mirror, args("x_mirror"), "Set whether or not x_mirroring is occuring\n \nx_mirror - whether x_mirroring should be applied\n")
		.def("set_pre_trans", (void (EMAN::Transform::*)(const EMAN::Vec2f&) ) &EMAN::Transform::set_pre_trans, args("v"), "Set the translational component of the matrix as though it was MSRT_ not MTSR, where\nT_ is the pre translation. Internally the correct form of MTSR is computed.\n \nv - the vector Vec2f that is the pre trans\n")
		.def("set_pre_trans", (void (EMAN::Transform::*)(const EMAN::Vec3f&) ) &EMAN::Transform::set_pre_trans, args("v"), "Set the translational component of the matrix as though it was MSRT_ not MTSR, where\nT_ is the pre translation. Internally the correct form of MTSR is computed.\n \nv - the vector Vec3f that is the pre trans\n")
		.def("get_pre_trans", &EMAN::Transform::get_pre_trans, "Get the translation vector as though this object was MSRT_ not MTSR, where T_ is what you want\nNote M means post x mirror, T means translation, S means scale, and R means rotaiton\n \nreturn the pre translation vector\n")
		.def("get_pre_trans_2d", &EMAN::Transform::get_pre_trans_2d, "2D version of getting the translation vector as though this object was MSRT_ not MTSR, where T_ is what you want\nNote M means post x mirror, T means translation, S means scale, and R means rotation\n \nreturn the pre translation vector\n")
		.def("get_pre_trans_2D", &EMAN::Transform::get_pre_trans_2d, "2D version of getting the translation vector as though this object was MSRT_ not MTSR, where T_ is what you want\nNote M means post x mirror, T means translation, S means scale, and R means rotation\n \nreturn the pre translation vector\n")
		.def("get_trans", &EMAN::Transform::get_trans, "Get the post trans as a vec3f\n \nreturn the  translation vector\n")
		.def("get_trans_2d", &EMAN::Transform::get_trans_2d, "Get the degenerant 2D post trans as a vec2f\n \nreturn the 2D translation vector\n")
		.def("get_trans_2D", &EMAN::Transform::get_trans_2d, "Get the degenerant 2D post trans as a vec2f\n \nreturn the 2D translation vector\n")
		.def("get_rotation", &EMAN::Transform::get_rotation, EMAN_Transform_get_rotation_overloads_0_1(args("euler_type"), "Get a rotation in any Euler format\n \neuler_type - the requested Euler type(default='eman')\n \nreturn a dictionary containing the key-entry pairs describing the rotations in terms of the requested Euler type"))
		.def("negate", &EMAN::Transform::negate, "Negates the Transform - a useful way, for example, for getting an orientation on the opposite side\nof the sphere\n \nreturn a transform that has been negated\n")
		.def("set_params", &EMAN::Transform::set_params, args("d"), "Set the parameters of the entire transform.\nkeys acted upon are 'type' - if this exists then the correct euler angles need to be included -\nalso 'tx','ty','tz', 'scale', and 'mirror'\n \nd - the dictionary containing the parameters\n")
		.def("get_params", &EMAN::Transform::get_params, args("euler_type"), "Get the parameters of the entire transform, using a specific euler convention\n \neuler_type - the euler type of the retrieved rotation\n \nreturn a dictionary containing the parameters\n")
		.def("transform", (EMAN::Vec2f (EMAN::Transform::*)(const EMAN::Vec2f&) const )&EMAN::Transform::transform, args("v"), "Transform a 2D vector using the internal transformation matrix\n \nv - a two dimensional vector to be transformed\n \nreturn the transformed vector")
		.def("transform", (EMAN::Vec3f (EMAN::Transform::*)(const EMAN::Vec3f&) const )&EMAN::Transform::transform, args("v"), "Transform a 3D vector using the internal transformation matrix\n \nv a three dimensional vector to be transformed\n \nreturn the transformed vector\n")
		.def("transform", (EMAN::Vec2f (EMAN::Transform::*)(const float&, const float&) const )&EMAN::Transform::transform, args("x", "y"), "Transform 2D coordinates using the internal transformation matrix\n \nx - the x coordinate of the transformed point\n \ny - the y coordinate of the transformed point\n \nreturn the transformed vector\n")
		.def("transform", (EMAN::Vec3f (EMAN::Transform::*)(const float&, const float&, const float&) const )&EMAN::Transform::transform, args("x", "y", "z"), "Transform 3D coordinates using the internal transformation matrix\n \nx - the x coordinate of the transformed point\ny - the y coordinate of the transformed point\nz - the z coordinate of the transformed point\n \nreturn the transformed vector\n")
		.def("inverse", &EMAN::Transform::inverse, "Get the inverse of this transformation matrix\n \nreturn the inverse of this transformation matrix\n")
		.def("invert", &EMAN::Transform::invert, "Get the inverse of this transformation matrix\n \nreturn the inverse of this transformation matrix\n")
		.def("transpose", &EMAN::Transform::transpose, "Get the transpose of this transformation matrix\n \nreturn the transpose of this transformation matrix\n")
		.def("orthogonalize", &EMAN::Transform::orthogonalize, "Reorthogonalize the rotation part of the matrix in place.\nDoes this by performing the SVD decomposition of the rotation matrix R\nsuch that R = USV^T - since the eigenvalues of a rotation matrix are all 1 we\nenforce that S should be the identity and produce a corrected matrix R' = UV^T\n")
		.def("printme", &EMAN::Transform::printme, "Print the contents of the internal matrix verbatim to standard out\n")
		.def("at", &EMAN::Transform::at, "Get the value stored in the internal transformation matrix at at coordinate (r,c)\n")
		.def("get_nsym", &EMAN::Transform::get_nsym, args("sym"), "get the number of symmetries associated with the given symmetry name\n")
		.def("get_sym", &EMAN::Transform::get_sym, args("sym", "n"), "Apply the symmetry deduced from the function arguments to this Transform and\nreturn the result\n")
		.def("get_sym_proj", &EMAN::Transform::get_sym_proj, args("sym", "s"), "Who knows  Apply the symmetry deduced from the function arguments to this Transform and\nreturn the result\n")
		.def("get_scale", &EMAN::Transform::get_scale, "Get the scale that was applied\n \nreturn the scale factor\n")
		.def("scale", &EMAN::Transform::scale, args("scale"), "Increment the scale\n \nscale - the amount to scale by\n")
		.def("to_identity", &EMAN::Transform::to_identity, "Force the internal matrix to become the identity\n")
		.def("is_identity", &EMAN::Transform::is_identity, "Returns whethers or this matrix is the identity\n")
		.def("is_rot_identity", &EMAN::Transform::is_rot_identity, "Returns whethers or this matrix rotation is the identity\n")
		.def("get_determinant", &EMAN::Transform::get_determinant, "Get the determinant of the matrix\n \nreturn the determinant\n")
		.def("get_params_inverse",&EMAN::Transform::get_params_inverse, args("euler_type"), "Get the parameters of the inverse of the transform as though it were in RSMT order not MTSR\n \neuler_type - the euler type of the retrieved rotation\n \nreturn a dictionary containing the parameters\n")
		.def("set_params_inverse",&EMAN::Transform::set_params_inverse, args("d"), "Set the parameters of the entire transform as though they there in the inverse format\nin other words, calling set_params_inverse(get_params_inverse()) should essentially leave\nthe object unchanged.\n \nd - the dictionary containing the inverse parameters\n")
		.def("get_matrix",&EMAN::Transform::get_matrix, "Get the transformation matrix using a vector.\n \nreturn a vector - 3 rows of 4 - that stores the values of the transformation matrix\n")
		.def("get_matrix_4x4",&EMAN::Transform::get_matrix_4x4, "Get the 4x4 transformation matrix using a vector.\n \nreturn a vector - 4 rows of 4 - that stores the values of the transformation matrix\n")
		.def("set_matrix",&EMAN::Transform::set_matrix, "Set the transformation matrix using a vector. Must be of length 12.\n \nv - the transformation matrix stored as a vector - 3 rows of 4.\n")
		.def( self * self )
		.def( other< EMAN::Vec3f >() * self )
		.def( self * other< EMAN::Vec3f >() )
		.def( self * other< EMAN::Vec2f >() )
		.staticmethod("get_nsym")
		.def("icos_5_to_2", &EMAN::Transform::icos_5_to_2)
		.staticmethod("icos_5_to_2")
		.def("tet_3_to_2", &EMAN::Transform::tet_3_to_2)
		.staticmethod("tet_3_to_2")
		.def("__eq__", (bool (EMAN::Transform::*)(const EMAN::Transform&) const)&EMAN::Transform::operator==)
		.def("__ne__", (bool (EMAN::Transform::*)(const EMAN::Transform&) const)&EMAN::Transform::operator!=)
	;

}


