
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <aligner.h>
#include <cmp.h>
#include <ctf.h>
#include <emdata.h>
#include <emfft.h>
#include <filter.h>
#include <io.h>
#include <pylist.h>
#include <transform.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_read_image_overloads_1_5, read_image, 1, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_write_image_overloads_1_4, write_image, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_append_image_overloads_1_3, append_image, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_filter_overloads_1_2, filter, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_align_overloads_2_3, align, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_copy_overloads_0_2, copy, 0, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_to_mass_center_overloads_0_1, to_mass_center, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_fast_translate_overloads_0_1, fast_translate, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_rotate_translate_overloads_0_5, rotate_translate, 0, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_fast_rotate_translate_overloads_0_1, fast_rotate_translate, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_little_big_dot_overloads_1_2, little_big_dot, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_ccf_overloads_1_3, calc_ccf, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_make_rotational_footprint_overloads_0_2, make_rotational_footprint, 0, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_ccfx_overloads_1_4, calc_ccfx, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_mutual_correlation_overloads_1_3, calc_mutual_correlation, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_unwrap_overloads_0_6, unwrap, 0, 6)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_apply_radial_func_overloads_3_4, apply_radial_func, 3, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_add_random_noise_overloads_4_5, add_random_noise, 4, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_hist_overloads_1_4, calc_hist, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_dist_overloads_1_2, calc_dist, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_flcf_overloads_1_3, calc_flcf, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_dot_overloads_1_2, dot, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_common_lines_overloads_2_5, common_lines, 2, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_common_lines_real_overloads_2_4, common_lines_real, 2, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_cut_slice_overloads_2_6, cut_slice, 2, 6)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_uncut_slice_overloads_2_5, uncut_slice, 2, 5)


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyEMData2)
{
    class_< EMAN::EMData >("EMData", init<  >())
        .def(init< const EMAN::EMData& >())
        .def_readwrite("HEADER_ONLY", &EMAN::EMData::HEADER_ONLY)
        .def_readwrite("HEADER_AND_DATA", &EMAN::EMData::HEADER_AND_DATA)
        .def_readwrite("IS_3D", &EMAN::EMData::IS_3D)
        .def_readwrite("NOT_3D", &EMAN::EMData::NOT_3D)
        .def_readwrite("DATA_READ_ONLY", &EMAN::EMData::DATA_READ_ONLY)
        .def_readwrite("DATA_READ_WRITE", &EMAN::EMData::DATA_READ_WRITE)
        .def("read_image", &EMAN::EMData::read_image, EMAN_EMData_read_image_overloads_1_5())
        .def("write_image", &EMAN::EMData::write_image, EMAN_EMData_write_image_overloads_1_4())
        .def("append_image", &EMAN::EMData::append_image, EMAN_EMData_append_image_overloads_1_3())
        .def("filter", &EMAN::EMData::filter, EMAN_EMData_filter_overloads_1_2())
        .def("cmp", &EMAN::EMData::cmp)
        .def("align", &EMAN::EMData::align, return_value_policy< manage_new_object >(), EMAN_EMData_align_overloads_2_3())
        .def("project", &EMAN::EMData::project, return_value_policy< manage_new_object >())
        .def("copy", &EMAN::EMData::copy, return_value_policy< manage_new_object >(), EMAN_EMData_copy_overloads_0_2())
        .def("copy_head", &EMAN::EMData::copy_head, return_value_policy< manage_new_object >())
        .def("get_clip", &EMAN::EMData::get_clip, return_value_policy< manage_new_object >())
        .def("insert_clip", &EMAN::EMData::insert_clip)
        .def("do_fft", &EMAN::EMData::do_fft, return_value_policy< manage_new_object >())
        .def("do_ift", &EMAN::EMData::do_ift, return_value_policy< manage_new_object >())
        .def("gimme_fft", &EMAN::EMData::gimme_fft)
        .def("normalize_slice", &EMAN::EMData::normalize_slice)
        .def("ri2ap", &EMAN::EMData::ri2ap)
        .def("ap2ri", &EMAN::EMData::ap2ri)
        .def("to_mass_center", &EMAN::EMData::to_mass_center, EMAN_EMData_to_mass_center_overloads_0_1())
        .def("rotate_x", &EMAN::EMData::rotate_x)
        .def("rotate_180", &EMAN::EMData::rotate_180)
        .def("fast_translate", &EMAN::EMData::fast_translate, EMAN_EMData_fast_translate_overloads_0_1())
        .def("rotate_translate", &EMAN::EMData::rotate_translate, EMAN_EMData_rotate_translate_overloads_0_5())
        .def("fast_rotate_translate", &EMAN::EMData::fast_rotate_translate, EMAN_EMData_fast_rotate_translate_overloads_0_1())
        .def("dot_rotate_translate", &EMAN::EMData::dot_rotate_translate)
        .def("little_big_dot", &EMAN::EMData::little_big_dot, return_value_policy< manage_new_object >(), EMAN_EMData_little_big_dot_overloads_1_2())
        .def("do_radon", &EMAN::EMData::do_radon, return_value_policy< manage_new_object >())
        .def("calc_ccf", &EMAN::EMData::calc_ccf, return_value_policy< manage_new_object >(), EMAN_EMData_calc_ccf_overloads_1_3())
        .def("make_rotational_footprint", &EMAN::EMData::make_rotational_footprint, return_value_policy< manage_new_object >(), EMAN_EMData_make_rotational_footprint_overloads_0_2())
        .def("calc_ccfx", &EMAN::EMData::calc_ccfx, return_value_policy< manage_new_object >(), EMAN_EMData_calc_ccfx_overloads_1_4())
        .def("calc_mutual_correlation", &EMAN::EMData::calc_mutual_correlation, return_value_policy< manage_new_object >(), EMAN_EMData_calc_mutual_correlation_overloads_1_3())
        .def("unwrap", &EMAN::EMData::unwrap, return_value_policy< manage_new_object >(), EMAN_EMData_unwrap_overloads_0_6())
        .def("mean_shrink", &EMAN::EMData::mean_shrink)
        .def("median_shrink", &EMAN::EMData::median_shrink)
        .def("apply_radial_func", &EMAN::EMData::apply_radial_func, EMAN_EMData_apply_radial_func_overloads_3_4())
        .def("calc_radial_dist", (std::vector<float,std::allocator<float> > (EMAN::EMData::*)(int, float, float) )&EMAN::EMData::calc_radial_dist)
        .def("calc_radial_dist", (std::vector<float,std::allocator<float> > (EMAN::EMData::*)(int, float, float, float, float) )&EMAN::EMData::calc_radial_dist)
        .def("add", (int (EMAN::EMData::*)(float) )&EMAN::EMData::add)
        .def("add", (int (EMAN::EMData::*)(const EMAN::EMData&) )&EMAN::EMData::add)
        .def("sub", (int (EMAN::EMData::*)(float) )&EMAN::EMData::sub)
        .def("sub", (int (EMAN::EMData::*)(const EMAN::EMData&) )&EMAN::EMData::sub)
        .def("mult", (int (EMAN::EMData::*)(float) )&EMAN::EMData::mult)
        .def("mult", (int (EMAN::EMData::*)(const EMAN::EMData&) )&EMAN::EMData::mult)
        .def("div", (int (EMAN::EMData::*)(float) )&EMAN::EMData::div)
        .def("div", (int (EMAN::EMData::*)(const EMAN::EMData&) )&EMAN::EMData::div)
        .def("done_data", &EMAN::EMData::done_data)
        .def("update", &EMAN::EMData::update)
        .def("to_zero", &EMAN::EMData::to_zero)
        .def("to_one", &EMAN::EMData::to_one)
        .def("dump_data", &EMAN::EMData::dump_data)
        .def("add_incoherent", &EMAN::EMData::add_incoherent)
        .def("add_mask_shell", &EMAN::EMData::add_mask_shell)
        .def("add_random_noise", &EMAN::EMData::add_random_noise, EMAN_EMData_add_random_noise_overloads_4_5())
        .def("calc_fourier_shell_correlation", &EMAN::EMData::calc_fourier_shell_correlation)
        .def("calc_hist", &EMAN::EMData::calc_hist, EMAN_EMData_calc_hist_overloads_1_4())
        .def("calc_az_dist", &EMAN::EMData::calc_az_dist)
        .def("calc_dist", &EMAN::EMData::calc_dist, EMAN_EMData_calc_dist_overloads_1_2())
        .def("calc_flcf", &EMAN::EMData::calc_flcf, return_value_policy< manage_new_object >(), EMAN_EMData_calc_flcf_overloads_1_3())
        .def("convolute", &EMAN::EMData::convolute, return_value_policy< manage_new_object >())
        .def("has_ctff", &EMAN::EMData::has_ctff)
        .def("dot", &EMAN::EMData::dot, EMAN_EMData_dot_overloads_1_2())
        .def("common_lines", &EMAN::EMData::common_lines, EMAN_EMData_common_lines_overloads_2_5())
        .def("common_lines_real", &EMAN::EMData::common_lines_real, EMAN_EMData_common_lines_real_overloads_2_4())
        .def("cut_slice", &EMAN::EMData::cut_slice, EMAN_EMData_cut_slice_overloads_2_6())
        .def("uncut_slice", &EMAN::EMData::uncut_slice, EMAN_EMData_uncut_slice_overloads_2_5())
        .def("get_edge_mean", &EMAN::EMData::get_edge_mean)
        .def("get_circle_mean", &EMAN::EMData::get_circle_mean)
        .def("setup_insert_slice", &EMAN::EMData::setup_insert_slice)
        .def("get_ctf", &EMAN::EMData::get_ctf, return_internal_reference< 1 >())
        .def("set_ctf", &EMAN::EMData::set_ctf)
        .def("get_translation", &EMAN::EMData::get_translation)
        .def("set_translation", &EMAN::EMData::set_translation)
        .def("get_rotation", &EMAN::EMData::get_rotation)
        .def("get_trans_align", &EMAN::EMData::get_trans_align)
        .def("set_size", &EMAN::EMData::set_size)
        .def("set_path", &EMAN::EMData::set_path)
        .def("set_pathnum", &EMAN::EMData::set_pathnum)
        .def("get_row", &EMAN::EMData::get_row, return_value_policy< manage_new_object >())
        .def("set_row", &EMAN::EMData::set_row)
        .def("get_col", &EMAN::EMData::get_col, return_value_policy< manage_new_object >())
        .def("set_col", &EMAN::EMData::set_col)
        .def("set_talign_params", (void (EMAN::EMData::*)(float, float) )&EMAN::EMData::set_talign_params)
        .def("set_talign_params", (void (EMAN::EMData::*)(float, float, float) )&EMAN::EMData::set_talign_params)
        .def("set_ralign_params", (void (EMAN::EMData::*)(float, float, float) )&EMAN::EMData::set_ralign_params)
        .def("set_ralign_params", (void (EMAN::EMData::*)(const EMAN::Rotation&) )&EMAN::EMData::set_ralign_params)
        .def("get_align_score", &EMAN::EMData::get_align_score)
        .def("set_align_score", &EMAN::EMData::set_align_score)
        .def("get_density_center", &EMAN::EMData::get_density_center)
        .def("get_sigma_diff", &EMAN::EMData::get_sigma_diff)
        .def("get_attr_dict", &EMAN::EMData::get_attr_dict)
        .def("set_attr_dict", &EMAN::EMData::set_attr_dict)
        .def("get_max", &EMAN::EMData::get_max)
        .def("get_min", &EMAN::EMData::get_min)
        .def("get_mean", &EMAN::EMData::get_mean)
        .def("get_sigma", &EMAN::EMData::get_sigma)
        .def("get_skewness", &EMAN::EMData::get_skewness)
        .def("get_kurtosis", &EMAN::EMData::get_kurtosis)
        .def("get_min_location", &EMAN::EMData::get_min_location)
        .def("get_max_location", &EMAN::EMData::get_max_location)
        .def("get_min_index", &EMAN::EMData::get_min_index)
        .def("get_max_index", &EMAN::EMData::get_max_index)
        .def("get_xsize", &EMAN::EMData::get_xsize)
        .def("get_ysize", &EMAN::EMData::get_ysize)
        .def("get_zsize", &EMAN::EMData::get_zsize)
        .def("get_ndim", &EMAN::EMData::get_ndim)
        .def("get_parent", &EMAN::EMData::get_parent, return_value_policy< manage_new_object >())
        .def("set_parent", &EMAN::EMData::set_parent)
        .def("get_name", &EMAN::EMData::get_name)
        .def("set_name", &EMAN::EMData::set_name)
        .def("get_value_at", (float (EMAN::EMData::*)(int, int, int) const)&EMAN::EMData::get_value_at)
        .def("get_value_at", (float (EMAN::EMData::*)(int, int) const)&EMAN::EMData::get_value_at)
        .def("sget_value_at", (float (EMAN::EMData::*)(int, int, int) const)&EMAN::EMData::sget_value_at)
        .def("sget_value_at", (float (EMAN::EMData::*)(int, int) const)&EMAN::EMData::sget_value_at)
        .def("get_value_at_interp", (float (EMAN::EMData::*)(float, float) const)&EMAN::EMData::get_value_at_interp)
        .def("get_value_at_interp", (float (EMAN::EMData::*)(float, float, float) const)&EMAN::EMData::get_value_at_interp)
        .def("set_value_at", (void (EMAN::EMData::*)(int, int, int, float) )&EMAN::EMData::set_value_at)
        .def("set_value_at", (void (EMAN::EMData::*)(int, int, float) )&EMAN::EMData::set_value_at)
        .def("is_complex", &EMAN::EMData::is_complex)
        .def("set_complex", &EMAN::EMData::set_complex)
        .def("is_complex_x", &EMAN::EMData::is_complex_x)
        .def("set_complex_x", &EMAN::EMData::set_complex_x)
        .def("is_flipped", &EMAN::EMData::is_flipped)
        .def("set_flipped", &EMAN::EMData::set_flipped)
        .def("is_ri", &EMAN::EMData::is_ri)
        .def("set_ri", &EMAN::EMData::set_ri)
        .def( self / other< float >() )
        .def( self * other< float >() )
        .def( self - other< float >() )
        .def( self + other< float >() )
        .def( other< float >() / self )
        .def( other< float >() + self )
        .def( other< float >() - self )
        .def( self + self )
        .def( self - self )
        .def( other< float >() * self )
        .def( self * self )
        .def( self / self )
        .def( self += other< float >() )
        .def( self -= other< float >() )
        .def( self *= other< float >() )
        .def( self /= other< float >() )
        .def( self += self )
        .def( self -= self )
        .def( self *= self )
        .def( self /= self )
    ;


    EMAN::vector_to_python<int>();
    EMAN::vector_to_python<float>();
    EMAN::vector_to_python<std::string>();
    EMAN::vector_to_python<EMAN::EMData*>();

    EMAN::vector_from_python<int>();
    EMAN::vector_from_python<float>();
    EMAN::vector_from_python<std::string>();
    EMAN::vector_from_python<EMAN::EMData*>();

    EMAN::map_to_python<EMAN::EMObject>();
    EMAN::map_from_python<EMAN::EMObject>();

    EMAN::Dict_to_python();
    EMAN::Dict_from_python();
}

