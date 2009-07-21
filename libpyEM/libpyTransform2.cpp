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

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_Transform3D_get_rotation_overloads_0_1, get_rotation, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_Transform3D_get_pretrans_overloads_0_1, get_pretrans, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_Transform3D_get_posttrans_overloads_0_1, get_posttrans, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_Transform3D_set_pretrans_overloads_2_3, set_pretrans, 2, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_Transform3D_set_posttrans_overloads_2_3, set_posttrans, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_Transform_settrans_overloads_2_3, set_trans, 2, 3)

}// namespace


// Module ======================================================================

BOOST_PYTHON_MODULE(libpyTransform2)
{
	def("dump_orientgens", &EMAN::dump_orientgens);
	def("dump_orientgens_list", &EMAN::dump_orientgens_list);

	def("dump_symmetries", &EMAN::dump_symmetries);
	def("dump_symmetries_list", &EMAN::dump_symmetries_list);

	class_< EMAN::Symmetry3D, boost::noncopyable, EMAN_Symmetry3D_Wrapper >("Symmetry3D", init<  >())
		.def("get_delimiters", pure_virtual(&EMAN::Symmetry3D::get_delimiters))
		.def("get_sym", pure_virtual(&EMAN::Symmetry3D::get_sym))
		.def("is_in_asym_unit", pure_virtual(&EMAN::Symmetry3D::is_in_asym_unit))
		.def("get_asym_unit_triangles", pure_virtual(&EMAN::Symmetry3D::get_asym_unit_triangles))
		.def("get_nsym",pure_virtual(&EMAN::Symmetry3D::get_nsym))
		.def("get_name",pure_virtual(&EMAN::Symmetry3D::get_name))
		.def("is_platonic_sym",pure_virtual(&EMAN::Symmetry3D::is_platonic_sym))
		.def("get_az_alignment_offset", pure_virtual(&EMAN::Symmetry3D::get_az_alignment_offset))
		.def("is_h_sym", pure_virtual(&EMAN::Symmetry3D::is_h_sym))
		.def("get_params", &EMAN::Symmetry3D::get_params)
		.def("gen_orientations", &EMAN::Symmetry3D::gen_orientations)
		.def("get_asym_unit_points", pure_virtual(&EMAN::Symmetry3D::get_asym_unit_points))
		.def("insert_params", pure_virtual(&EMAN::Symmetry3D::insert_params))
		.def("get_params", &EMAN::Symmetry3D::get_params)
		.def("reduce", &EMAN::Symmetry3D::reduce)
		.def("in_which_asym_unit", &EMAN::Symmetry3D::in_which_asym_unit)
		.def("point_in_which_asym_unit", &EMAN::Symmetry3D::point_in_which_asym_unit)
		.def("get_touching_au_transforms",&EMAN::Symmetry3D::get_touching_au_transforms)
		.def("get_syms", &EMAN::Symmetry3D::get_syms)
		.def("get_symmetries", &EMAN::Symmetry3D::get_symmetries)
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
	typedef float (EMAN::Vec3f::*dot_float)(const EMAN::Vec3f&) const;
	typedef EMAN::Vec3f (EMAN::Vec3f::*cross_float)(const EMAN::Vec3f&) const;
	typedef void (EMAN::Vec3f::*set_value_at_float)(int, const float&);
	typedef void (EMAN::Vec3f::*set_value_float)(const float&, const float&,const float&);
	typedef float (EMAN::Vec3f::*at_float)(int) const;
    class_< EMAN::Vec3f >("Vec3f", init<  >())
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
    class_< EMAN::Vec3i >("Vec3i", init<  >())
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
	class_< EMAN::Vec2f >("Vec2f", init<  >())
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

    class_< EMAN::Quaternion >("Quaternion", init<  >())
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

	scope* EMAN_Transform3D_scope = new scope(
	class_< EMAN::Transform3D >("Transform3D", init<  >())
		.def(init< const EMAN::Transform3D& >())
		.def(init< const float&, const float&, const float& >())
// 		.def(init< const EMAN::Dict&, const string&, optional<const EMAN::Transform3D::EulerType> >())
		.def(init< const float&, const float&, const float&, const EMAN::Vec3f& >())
		.def(init< EMAN::Transform3D::EulerType, const float&, const float&, const float& >())
		.def(init< EMAN::Transform3D::EulerType, const float&, const float&, const float&, const float&>())
		.def(init< EMAN::Transform3D::EulerType, const EMAN::Dict& >())
		.def(init< const EMAN::Vec3f&, const float&, const float&, const float&, const EMAN::Vec3f& >())
		.def(init< const float&, const float&, const float&, const float&, const float&, const float&, const float&, const float&, const float& >())
		.def_readonly("ERR_LIMIT", &EMAN::Transform3D::ERR_LIMIT)
		.def("set_posttrans", (void (EMAN::Transform3D::*)(const float&, const float&, const float&) )&EMAN::Transform3D::set_posttrans,EMAN_Transform3D_set_posttrans_overloads_2_3())
		.def("set_posttrans", (void (EMAN::Transform3D::*)(const EMAN::Vec3f&) )&EMAN::Transform3D::set_posttrans)
		.def("set_posttrans", (void (EMAN::Transform3D::*)(const EMAN::Vec2f&) )&EMAN::Transform3D::set_posttrans)
		.def("apply_scale", &EMAN::Transform3D::apply_scale)
		.def("set_scale", &EMAN::Transform3D::set_scale)
		.def("orthogonalize", &EMAN::Transform3D::orthogonalize)
		.def("transpose", &EMAN::Transform3D::transpose)
		.def("set_rotation", (void (EMAN::Transform3D::*)(const float&, const float&, const float&) )&EMAN::Transform3D::set_rotation)
		.def("set_rotation", (void (EMAN::Transform3D::*)(EMAN::Transform3D::EulerType, const float&, const float&, const float&) )&EMAN::Transform3D::set_rotation)
		.def("set_rotation", (void (EMAN::Transform3D::*)(EMAN::Transform3D::EulerType, const float&, const float&, const float&, const float&) )&EMAN::Transform3D::set_rotation)
		.def("set_rotation", (void (EMAN::Transform3D::*)(const float&, const float&, const float&, const float&, const float&, const float&, const float&, const float&, const float&) )&EMAN::Transform3D::set_rotation)
		.def("set_rotation", (void (EMAN::Transform3D::*)(EMAN::Transform3D::EulerType, const EMAN::Dict&) )&EMAN::Transform3D::set_rotation)
		.def("set_rotation", (void (EMAN::Transform3D::*)(const EMAN::Vec3f&, const EMAN::Vec3f&, const EMAN::Vec3f&, const EMAN::Vec3f&) )&EMAN::Transform3D::set_rotation)
		.def("get_post_x_mirror",&EMAN::Transform3D::get_post_x_mirror)
		.def("set_post_x_mirror",&EMAN::Transform3D::set_post_x_mirror)
		.def("get_mag", &EMAN::Transform3D::get_mag)
		.def("get_finger", &EMAN::Transform3D::get_finger)
		.def("get_pretrans", &EMAN::Transform3D::get_pretrans, EMAN_Transform3D_get_pretrans_overloads_0_1())
		.def("get_posttrans", &EMAN::Transform3D::get_posttrans, EMAN_Transform3D_get_posttrans_overloads_0_1())
		.def("get_total_posttrans", &EMAN::Transform3D::get_total_posttrans)
		.def("get_total_pretrans", &EMAN::Transform3D::get_total_pretrans)
		.def("get_center", &EMAN::Transform3D::get_center)
		.def("get_matrix3_col", &EMAN::Transform3D::get_matrix3_col)
		.def("get_matrix3_row", &EMAN::Transform3D::get_matrix3_row)
		.def("transform", &EMAN::Transform3D::transform)
		.def("rotate", &EMAN::Transform3D::rotate)
		.def("inverse", &EMAN::Transform3D::inverse)
		.def("inverseUsingAngs", &EMAN::Transform3D::inverseUsingAngs)
		.def("get_rotation", &EMAN::Transform3D::get_rotation, EMAN_Transform3D_get_rotation_overloads_0_1())
		.def("printme", &EMAN::Transform3D::printme)
		.def("at", &EMAN::Transform3D::at)
		.def("get_nsym", &EMAN::Transform3D::get_nsym)
		.def("get_sym", &EMAN::Transform3D::get_sym)
		.def("set_center", &EMAN::Transform3D::set_center)
		.def("set_pretrans", (void (EMAN::Transform3D::*)(const float&,const float&,const float&))&EMAN::Transform3D::set_pretrans, EMAN_Transform3D_set_pretrans_overloads_2_3())
		.def("set_pretrans", (void (EMAN::Transform3D::*)(const EMAN::Vec3f&))&EMAN::Transform3D::set_pretrans)
		.def("set_pretrans", (void (EMAN::Transform3D::*)(const EMAN::Vec2f&))&EMAN::Transform3D::set_pretrans)
		.def("get_scale", &EMAN::Transform3D::get_scale)
		.def("to_identity", &EMAN::Transform3D::to_identity)
		.def("is_identity", &EMAN::Transform3D::is_identity)
		.def("angles2tfvec", &EMAN::Transform3D::angles2tfvec)
		.staticmethod("angles2tfvec")
		.staticmethod("get_nsym")
		.def( other< EMAN::Vec3f >() * self )
		.def( self * self )
		.def( self * other< EMAN::Vec3f >() )
		.def( self * other< EMAN::Vec2f >() )
// 		.def( self + self )
// 		.def( self - self )
    );

    enum_< EMAN::Transform3D::EulerType >("EulerType")
        .value("MATRIX", EMAN::Transform3D::MATRIX)
        .value("UNKNOWN", EMAN::Transform3D::UNKNOWN)
        .value("XYZ", EMAN::Transform3D::XYZ)
        .value("IMAGIC", EMAN::Transform3D::IMAGIC)
        .value("SPIDER", EMAN::Transform3D::SPIDER)
        .value("QUATERNION", EMAN::Transform3D::QUATERNION)
        .value("SGIROT", EMAN::Transform3D::SGIROT)
        .value("MRC", EMAN::Transform3D::MRC)
        .value("SPIN", EMAN::Transform3D::SPIN)
        .value("EMAN", EMAN::Transform3D::EMAN)
    ;

	delete EMAN_Transform3D_scope;

	class_< EMAN::Transform >("Transform", init<  >())
//	class_< EMAN::Transform, std::auto_ptr<EMAN::Transform>  >("Transform", init<  >())
		.def(init< const EMAN::Transform& >())
		.def(init<const EMAN::Dict& >())
		.def(init<const std::vector<float>& >())
		.def_pickle(Transform_pickle_suite())
		.def_readonly("ERR_LIMIT", &EMAN::Transform::ERR_LIMIT)
		.def("set_trans", (void (EMAN::Transform::*)(const float&, const float&, const float&) )&EMAN::Transform::set_trans, EMAN_Transform_settrans_overloads_2_3())
		.def("set_trans", (void (EMAN::Transform::*)(const EMAN::Vec3f&) )&EMAN::Transform::set_trans)
		.def("set_trans", (void (EMAN::Transform::*)(const EMAN::Vec2f&) )&EMAN::Transform::set_trans)
		.def("set_scale", &EMAN::Transform::set_scale)
		.def("set_rotation", (void (EMAN::Transform::*)(const EMAN::Dict&) )&EMAN::Transform::set_rotation)
		.def("set_rotation", (void (EMAN::Transform::*)(const EMAN::Vec3f&) )&EMAN::Transform::set_rotation)
		.def("get_hflip_transform", (EMAN::Transform (EMAN::Transform::*)() const )&EMAN::Transform::get_hflip_transform)
		.def("get_vflip_transform", (EMAN::Transform (EMAN::Transform::*)() const )&EMAN::Transform::get_vflip_transform)
		.def("get_mirror",&EMAN::Transform::get_mirror)
		.def("set_mirror",&EMAN::Transform::set_mirror)
		.def("set_pre_trans", (void (EMAN::Transform::*)(const EMAN::Vec2f&) ) &EMAN::Transform::set_pre_trans)
		.def("set_pre_trans", (void (EMAN::Transform::*)(const EMAN::Vec3f&) ) &EMAN::Transform::set_pre_trans)
		.def("get_pre_trans", &EMAN::Transform::get_pre_trans)
		.def("get_pre_trans_2d", &EMAN::Transform::get_pre_trans_2d)
		.def("get_pre_trans_2D", &EMAN::Transform::get_pre_trans_2d)
		.def("get_trans", &EMAN::Transform::get_trans)
		.def("get_trans_2d", &EMAN::Transform::get_trans_2d)
		.def("get_trans_2D", &EMAN::Transform::get_trans_2d)
		.def("get_rotation", &EMAN::Transform::get_rotation)
		.def("spherical_opposite", &EMAN::Transform::spherical_opposite)
		.def("set_params", &EMAN::Transform::set_params)
		.def("get_params", &EMAN::Transform::get_params)
		.def("transform", (EMAN::Vec2f (EMAN::Transform::*)(const EMAN::Vec2f&) const )&EMAN::Transform::transform)
		.def("transform", (EMAN::Vec3f (EMAN::Transform::*)(const EMAN::Vec3f&) const )&EMAN::Transform::transform)
		.def("transform", (EMAN::Vec2f (EMAN::Transform::*)(const float&, const float&) const )&EMAN::Transform::transform)
		.def("transform", (EMAN::Vec3f (EMAN::Transform::*)(const float&, const float&, const float&) const )&EMAN::Transform::transform)
		.def("inverse", &EMAN::Transform::inverse)
		.def("invert", &EMAN::Transform::invert)
	        .def("transpose", &EMAN::Transform::transpose)
		.def("orthogonalize", &EMAN::Transform::orthogonalize)
		.def("printme", &EMAN::Transform::printme)
		.def("at", &EMAN::Transform::at)
		.def("get_nsym", &EMAN::Transform::get_nsym)
		.def("get_sym", &EMAN::Transform::get_sym)
		.def("get_scale", &EMAN::Transform::get_scale)
		.def("to_identity", &EMAN::Transform::to_identity)
		.def("is_identity", &EMAN::Transform::is_identity)
		.def("get_determinant", &EMAN::Transform::get_determinant)
		.def("get_params_inverse",&EMAN::Transform::get_params_inverse)
		.def("set_params_inverse",&EMAN::Transform::set_params_inverse)
		.def("get_matrix",&EMAN::Transform::get_matrix)
		.def("set_matrix",&EMAN::Transform::set_matrix)
		.def( self * self )
		.def( other< EMAN::Vec3f >() * self )
		.def( self * other< EMAN::Vec3f >() )
		.def( self * other< EMAN::Vec2f >() )
		.staticmethod("get_nsym")
		.def("icos_5_to_2", &EMAN::Transform::icos_5_to_2)
		.staticmethod("icos_5_to_2")
		.def("tet_3_to_2", &EMAN::Transform::tet_3_to_2)
		.staticmethod("tet_3_to_2")
	;

}


