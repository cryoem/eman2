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
#include <pointarray.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_PointArray_distmx_overloads_0_1, distmx, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_PointArray_mask_overloads_1_2, mask, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_PointArray_set_from_overloads_1_3, set_from, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_PointArray_set_from_overloads_2_4, set_from, 2, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_PointArray_set_from_density_map_overloads_4_5, set_from_density_map, 4, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_PointArray_sort_by_axis_overloads_0_1, sort_by_axis, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_PointArray_projection_by_nfft_overloads_2_3, projection_by_nfft, 2, 3)


}// namespace


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyPointArray2)
{
    scope* EMAN_PointArray_scope = new scope(
    class_< EMAN::PointArray >("PointArray",
    		"PointArray defines a double array of points with values in a 3D space.",
    		init<  >())
        .def(init< const EMAN::PointArray& >())
        .def(init< unsigned int >())
        .def("zero", &EMAN::PointArray::zero)
        .def("copy", &EMAN::PointArray::copy, return_value_policy< manage_new_object >())
        .def("get_number_points", &EMAN::PointArray::get_number_points)
        .def("__len__", &EMAN::PointArray::get_number_points)
        .def("set_number_points", &EMAN::PointArray::set_number_points)
         .def("get_points", &EMAN::PointArray::get_points)
        .def("read_from_pdb", &EMAN::PointArray::read_from_pdb)
        .def("save_to_pdb", &EMAN::PointArray::save_to_pdb)
        .def("get_center", &EMAN::PointArray::get_center)
        .def("center_to_zero", &EMAN::PointArray::center_to_zero)
        .def("get_bounding_box", &EMAN::PointArray::get_bounding_box)
        .def("get_vector_at", &EMAN::PointArray::get_vector_at)
        .def("__getitem__", &EMAN::PointArray::get_vector_at)
        .def("get_value_at", &EMAN::PointArray::get_value_at)
        .def("set_vector_at", (void (EMAN::PointArray::*)(int, EMAN::Vec3f, double) )&EMAN::PointArray::set_vector_at)
        .def("set_vector_at", (void (EMAN::PointArray::*)(int, std::vector<double,std::allocator<double> >) )&EMAN::PointArray::set_vector_at)
        .def("set_points_array", &EMAN::PointArray::set_points_array)
        .def("distmx", &EMAN::PointArray::distmx, return_value_policy< manage_new_object >(), EMAN_PointArray_distmx_overloads_0_1())
        .def("match_points", &EMAN::PointArray::match_points)
        .def("transform", &EMAN::PointArray::transform)
		.def("right_transform", &EMAN::PointArray::right_transform)
        .def("align_2d", &EMAN::PointArray::align_2d,return_value_policy< manage_new_object >())
        .def("align_trans_2d", &EMAN::PointArray::align_trans_2d)
        .def("mask", &EMAN::PointArray::mask, EMAN_PointArray_mask_overloads_1_2())
        .def("mask_asymmetric_unit", &EMAN::PointArray::mask_asymmetric_unit)
        .def("set_from", (void (EMAN::PointArray::*)(EMAN::PointArray*, const std::string&, EMAN::Transform*) )&EMAN::PointArray::set_from, EMAN_PointArray_set_from_overloads_1_3())
        .def("set_from", (void (EMAN::PointArray::*)(double*, int, const std::string&, EMAN::Transform*) )&EMAN::PointArray::set_from, EMAN_PointArray_set_from_overloads_2_4())
        .def("set_from", (void (EMAN::PointArray::*)(vector<float>) )&EMAN::PointArray::set_from)
        .def("set_from_density_map", &EMAN::PointArray::set_from_density_map, EMAN_PointArray_set_from_density_map_overloads_4_5())
        .def("sort_by_axis", &EMAN::PointArray::sort_by_axis, EMAN_PointArray_sort_by_axis_overloads_0_1())
        .def("pdb2mrc_by_nfft", &EMAN::PointArray::pdb2mrc_by_nfft, return_value_policy< manage_new_object >())
        .def("pdb2mrc_by_summation", &EMAN::PointArray::pdb2mrc_by_summation, return_value_policy< manage_new_object >())
        .def("projection_by_nfft", &EMAN::PointArray::projection_by_nfft, EMAN_PointArray_projection_by_nfft_overloads_2_3()[ return_value_policy< manage_new_object >() ])
        .def("projection_by_summation", &EMAN::PointArray::projection_by_summation, return_value_policy< manage_new_object >())
        .def("replace_by_summation", &EMAN::PointArray::replace_by_summation)
        .def("opt_from_proj", &EMAN::PointArray::opt_from_proj)
        .def("sim_set_pot_parms", &EMAN::PointArray::sim_set_pot_parms)
        .def("sim_minstep", &EMAN::PointArray::sim_minstep)
        .def("sim_minstep_seq", &EMAN::PointArray::sim_minstep_seq)
        .def("sim_rescale", &EMAN::PointArray::sim_rescale)
        .def("sim_potential", &EMAN::PointArray::sim_potential)
        .def("sim_printstat", &EMAN::PointArray::sim_printstat)
        .def("sim_add_point_double", &EMAN::PointArray::sim_add_point_double)
        .def("sim_add_point_one", &EMAN::PointArray::sim_add_point_one)
        .def("calc_total_length", &EMAN::PointArray::calc_total_length)
        .def("fit_helix", &EMAN::PointArray::fit_helix)
        .def("reverse_chain", &EMAN::PointArray::reverse_chain)
        .def("save_pdb_with_helix", &EMAN::PointArray::save_pdb_with_helix)
        .def("remove_helix_from_map", &EMAN::PointArray::remove_helix_from_map)
        .def("merge_to", &EMAN::PointArray::merge_to)
        .def("delete_point", &EMAN::PointArray::delete_point)
        .def("read_ca_from_pdb", &EMAN::PointArray::read_ca_from_pdb)
        .def("calc_transform", &EMAN::PointArray::calc_transform)
		
    );

    enum_< EMAN::PointArray::Density2PointsArrayAlgorithm >("Density2PointsArrayAlgorithm")
        .value("PEAKS_DIV", EMAN::PointArray::PEAKS_DIV)
        .value("KMEANS", EMAN::PointArray::KMEANS)
        .value("PEAKS_SUB", EMAN::PointArray::PEAKS_SUB)
    ;

    delete EMAN_PointArray_scope;

}

