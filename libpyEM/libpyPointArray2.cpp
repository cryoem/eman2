
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <pointarray.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_PointArray_set_from_overloads_2_3, set_from, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_PointArray_set_from_overloads_3_4, set_from, 3, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_PointArray_sort_by_axis_overloads_0_1, sort_by_axis, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_PointArray_projection_by_nfft_overloads_2_3, projection_by_nfft, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_PointArray_refine_overloads_1_4, refine, 1, 4)


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyPointArray2)
{
    scope* EMAN_PointArray_scope = new scope(
    class_< EMAN::PointArray >("PointArray", init<  >())
        .def(init< const EMAN::PointArray& >())
        .def(init< unsigned int >())
        .def("get_number_points", &EMAN::PointArray::get_number_points)
        .def("set_number_points", &EMAN::PointArray::set_number_points)
        .def("read_from_pdb", &EMAN::PointArray::read_from_pdb)
        .def("get_center", &EMAN::PointArray::get_center)
        .def("center_to_zero", &EMAN::PointArray::center_to_zero)
        .def("get_bounding_box", &EMAN::PointArray::get_bounding_box)
        .def("set_points_array", &EMAN::PointArray::set_points_array)
        .def("set_from", (void (EMAN::PointArray::*)(EMAN::PointArray*, EMAN::Transform*, std::string) )&EMAN::PointArray::set_from, EMAN_PointArray_set_from_overloads_2_3())
        .def("set_from", (void (EMAN::PointArray::*)(double*, unsigned int, EMAN::Transform*, std::string) )&EMAN::PointArray::set_from, EMAN_PointArray_set_from_overloads_3_4())
        .def("sort_by_axis", &EMAN::PointArray::sort_by_axis, EMAN_PointArray_sort_by_axis_overloads_0_1())
        .def("pdb2mrc_by_nfft", &EMAN::PointArray::pdb2mrc_by_nfft, return_value_policy< manage_new_object >())
        .def("pdb2mrc_by_summation", &EMAN::PointArray::pdb2mrc_by_summation, return_value_policy< manage_new_object >())
        .def("projection_by_nfft", &EMAN::PointArray::projection_by_nfft, return_value_policy< manage_new_object >(), EMAN_PointArray_projection_by_nfft_overloads_2_3())
        .def("projection_by_summation", &EMAN::PointArray::projection_by_summation, return_value_policy< manage_new_object >())
        .def("refine", &EMAN::PointArray::refine, EMAN_PointArray_refine_overloads_1_4())
    );

    enum_< EMAN::PointArray::Optimizer >("Optimizer")
        .value("ConjugateGradientFletcherReeves", EMAN::PointArray::ConjugateGradientFletcherReeves)
        .value("SimplexNelderMead", EMAN::PointArray::SimplexNelderMead)
        .value("ConjugateGradientBFGS", EMAN::PointArray::ConjugateGradientBFGS)
        .value("ConjugateGradientPolakRibiere", EMAN::PointArray::ConjugateGradientPolakRibiere)
    ;


    enum_< EMAN::PointArray::OptimizedParameters >("OptimizedParameters")
        .value("BeamTilt", EMAN::PointArray::BeamTilt)
        .value("Map", EMAN::PointArray::Map)
        .value("Scale", EMAN::PointArray::Scale)
        .value("Orientation", EMAN::PointArray::Orientation)
        .value("Defocus", EMAN::PointArray::Defocus)
        .value("Drift", EMAN::PointArray::Drift)
        .value("BFactor", EMAN::PointArray::BFactor)
        .value("Distortion", EMAN::PointArray::Distortion)
        .value("Astigmatism", EMAN::PointArray::Astigmatism)
        .value("DepthOfView", EMAN::PointArray::DepthOfView)
        .value("Center", EMAN::PointArray::Center)
    ;

    delete EMAN_PointArray_scope;

}

