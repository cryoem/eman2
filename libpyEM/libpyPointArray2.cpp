
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <pointarray.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

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
    class_< EMAN::PointArray >("PointArray", init<  >())
        .def(init< const EMAN::PointArray& >())
        .def(init< unsigned int >())
        .def("zero", &EMAN::PointArray::zero)
        .def("copy", &EMAN::PointArray::copy, return_value_policy< manage_new_object >())
        .def("get_number_points", &EMAN::PointArray::get_number_points)
        .def("set_number_points", &EMAN::PointArray::set_number_points)
        .def("read_from_pdb", &EMAN::PointArray::read_from_pdb)
        .def("save_to_pdb", &EMAN::PointArray::save_to_pdb)
        .def("get_center", &EMAN::PointArray::get_center)
        .def("center_to_zero", &EMAN::PointArray::center_to_zero)
        .def("get_bounding_box", &EMAN::PointArray::get_bounding_box)
        .def("get_vector_at", &EMAN::PointArray::get_vector_at)
        .def("get_value_at", &EMAN::PointArray::get_value_at)
        .def("set_vector_at", (void (EMAN::PointArray::*)(int, EMAN::Vec3f, double) )&EMAN::PointArray::set_vector_at)
        .def("set_vector_at", (void (EMAN::PointArray::*)(int, std::vector<double,std::allocator<double> >) )&EMAN::PointArray::set_vector_at)
        .def("set_points_array", &EMAN::PointArray::set_points_array)
        .def("mask", &EMAN::PointArray::mask, EMAN_PointArray_mask_overloads_1_2())
        .def("mask_asymmetric_unit", &EMAN::PointArray::mask_asymmetric_unit)
        .def("set_from", (void (EMAN::PointArray::*)(EMAN::PointArray*, std::string, EMAN::Transform*) )&EMAN::PointArray::set_from, EMAN_PointArray_set_from_overloads_1_3())
        .def("set_from", (void (EMAN::PointArray::*)(double*, unsigned int, std::string, EMAN::Transform*) )&EMAN::PointArray::set_from, EMAN_PointArray_set_from_overloads_2_4())
        .def("set_from_density_map", &EMAN::PointArray::set_from_density_map, EMAN_PointArray_set_from_density_map_overloads_4_5())
        .def("sort_by_axis", &EMAN::PointArray::sort_by_axis, EMAN_PointArray_sort_by_axis_overloads_0_1())
        .def("pdb2mrc_by_nfft", &EMAN::PointArray::pdb2mrc_by_nfft, return_value_policy< manage_new_object >())
        .def("pdb2mrc_by_summation", &EMAN::PointArray::pdb2mrc_by_summation, return_value_policy< manage_new_object >())
        .def("projection_by_nfft", &EMAN::PointArray::projection_by_nfft, return_value_policy< manage_new_object >(), EMAN_PointArray_projection_by_nfft_overloads_2_3())
        .def("projection_by_summation", &EMAN::PointArray::projection_by_summation, return_value_policy< manage_new_object >())
    );

    enum_< EMAN::PointArray::Density2PointsArrayAlgorithm >("Density2PointsArrayAlgorithm")
        .value("PEAKS_DIV", EMAN::PointArray::PEAKS_DIV)
        .value("KMEANS", EMAN::PointArray::KMEANS)
        .value("PEAKS_SUB", EMAN::PointArray::PEAKS_SUB)
    ;

    delete EMAN_PointArray_scope;

}

