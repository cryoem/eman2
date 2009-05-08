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
#include <aligner.h>
#include <cmp.h>
#include <ctf.h>
#include <emdata_pickle.h>
#include <emdata_wrapitems.h>
#include <emfft.h>
#include <processor.h>
#include <transform.h>
#include <xydata.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_read_image_overloads_1_5, read_image, 1, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_write_image_overloads_1_7, write_image, 1, 7)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_append_image_overloads_1_3, append_image, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_write_lst_overloads_1_4, write_lst, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_print_image_overloads_0_2, print_image, 0, 2)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_EMData_read_images_overloads_1_3, EMAN::EMData::read_images, 1, 3)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_EMData_read_images_ext_overloads_3_5, EMAN::EMData::read_images_ext, 3, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_set_size_overloads_1_3, set_size, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_set_complex_size_overloads_1_3, set_complex_size, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_get_attr_default_overloads_1_2, get_attr_default, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_clip_inplace_overloads_1_2,clip_inplace, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_process_inplace_overloads_1_2, process_inplace, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_process_overloads_1_2, process, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_cmp_overloads_2_3, cmp, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_align_overloads_2_5, align, 2, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_project_overloads_1_2, project, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_backproject_overloads_1_2, backproject, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_insert_scaled_sum_overloads_2_4, insert_scaled_sum, 2, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_add_overloads_1_2, add, 1, 2)

#ifndef	_WIN32
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_set_array_offsets_overloads_0_3, set_array_offsets, 0, 3)
#endif

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_real2complex_overloads_0_1, real2complex, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_FH2F_overloads_2_3, FH2F, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_FH2Real_overloads_2_3, FH2Real, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_fourier_shell_correlation_overloads_1_2, calc_fourier_shell_correlation, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_nn_overloads_3_4, nn, 3, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_nn_SSNR_overloads_4_5, nn_SSNR, 4, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_nn_SSNR_ctf_overloads_5_6, nn_SSNR_ctf, 5, 6)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_rot_scale_trans2D_overloads_1_4, rot_scale_trans2D, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_rot_scale_conv_overloads_4_5, rot_scale_conv, 4, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_rot_scale_conv_new_overloads_4_5, rot_scale_conv_new, 4, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_downsample_overloads_1_2, downsample, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_getconvpt2d_kbi0_overloads_3_4, getconvpt2d_kbi0, 3, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_pad_fft_overloads_0_1, pad_fft, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_FourInterpol_overloads_1_4, FourInterpol, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_FourTruncate_overloads_1_4, FourTruncate, 1, 4)

//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_FourInterpol_i_overloads_1_4, FourInterpol_i, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_Four_ds_overloads_1_4, Four_ds, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_Four_shuf_ds_cen_us_overloads_1_4, Four_shuf_ds_cen_us, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_filter_by_image_overloads_1_2, filter_by_image, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_replace_amplitudes_overloads_1_2, replace_amplitudes, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_get_rotated_clip_overloads_2_3, get_rotated_clip, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_little_big_dot_overloads_1_2, little_big_dot, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_ccf_overloads_1_3, calc_ccf, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_ccfx_overloads_1_4, calc_ccfx, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_make_rotational_footprint_overloads_0_1, make_rotational_footprint, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_make_rotational_footprint_e1_overloads_0_1, make_rotational_footprint_e1, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_make_rotational_footprint_cmc_overloads_0_1, make_rotational_footprint_cmc, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_mutual_correlation_overloads_1_3, calc_mutual_correlation, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_unwrap_overloads_0_7, unwrap, 0, 7)

#ifdef EMAN2_USING_CUDA
//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_unwrap_cuda_overloads_0_6, unwrap_cuda, 0, 6)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_get_clip_cuda_overloads_1_2, get_clip_cuda, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_make_rotational_footprint_cuda_overloads_0_1, make_rotational_footprint_cuda, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_do_ift_cuda_overloads_0_1, do_ift_cuda, 0, 1)
#endif //EMAN2_USING_CUDA
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_apply_radial_func_overloads_3_4, apply_radial_func, 3, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_hist_overloads_0_3, calc_hist, 0, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_dist_overloads_1_2, calc_dist, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_common_lines_overloads_2_5, common_lines, 2, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_common_lines_real_overloads_2_4, common_lines_real, 2, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_cut_slice_overloads_2_3, cut_slice, 2, 3)

// BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_uncut_slice_overloads_1_2, uncut_slice, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_get_clip_overloads_1_2, get_clip, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_mult_overloads_1_2, mult, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_norm_pad_overloads_2_3, EMAN::EMData::norm_pad, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_read_data_2_3, read_data, 2, 6)

}// namespace


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyEMData2)
{
    scope* EMAN_EMData_scope = new scope(
    class_< EMAN::EMData >("EMData", init<  >())
//    class_< EMAN::EMData, std::auto_ptr<EMAN::EMData> >("EMData", init<  >())
	.def_pickle(EMData_pickle_suite())
	.def(init< const EMAN::EMData& >())
	.def(init< const std::string&, optional< int > >())
	.def(init< int, int, optional< int, bool > >())
	.def_readwrite("totalalloc", &EMAN::EMData::totalalloc)
	.def("read_image", &EMAN::EMData::read_image, EMAN_EMData_read_image_overloads_1_5())
	.def("write_image", &EMAN::EMData::write_image, EMAN_EMData_write_image_overloads_1_7())
	.def("append_image", &EMAN::EMData::append_image, EMAN_EMData_append_image_overloads_1_3())
	.def("write_lst", &EMAN::EMData::write_lst, EMAN_EMData_write_lst_overloads_1_4())
	.def("print_image", &EMAN::EMData::print_image, EMAN_EMData_print_image_overloads_0_2())
	.def("read_images", &EMAN::EMData::read_images, EMAN_EMData_read_images_overloads_1_3())
	.def("read_images_ext", &EMAN::EMData::read_images_ext, EMAN_EMData_read_images_ext_overloads_3_5())
	.def("get_fft_amplitude", &EMAN::EMData::get_fft_amplitude, return_value_policy< manage_new_object >())
	.def("get_fft_amplitude2D", &EMAN::EMData::get_fft_amplitude2D, return_value_policy< manage_new_object >())
	.def("get_fft_phase", &EMAN::EMData::get_fft_phase, return_value_policy< manage_new_object >())
	.def("update", &EMAN::EMData::update)
	.def("has_ctff", &EMAN::EMData::has_ctff)
	.def("calc_center_density", &EMAN::EMData::calc_center_density)
	.def("calc_sigma_diff", &EMAN::EMData::calc_sigma_diff)
	.def("calc_min_location", &EMAN::EMData::calc_min_location)
	.def("calc_max_location", &EMAN::EMData::calc_max_location)
	.def("calc_max_location_wrap", &EMAN::EMData::calc_max_location_wrap)
	.def("calc_center_of_mass", &EMAN::EMData::calc_center_of_mass)
	.def("calc_min_index", &EMAN::EMData::calc_min_index)
	.def("calc_max_index", &EMAN::EMData::calc_max_index)
	.def("calc_highest_locations", &EMAN::EMData::calc_highest_locations)
	.def("get_edge_mean", &EMAN::EMData::get_edge_mean)
	.def("get_circle_mean", &EMAN::EMData::get_circle_mean)
	.def("get_ctf", &EMAN::EMData::get_ctf, return_value_policy< manage_new_object >())
	.def("set_ctf", &EMAN::EMData::set_ctf)
	.def("get_translation", &EMAN::EMData::get_translation)
	.def("set_translation", (void (EMAN::EMData::*)(const EMAN::Vec3f&) )&EMAN::EMData::set_translation)
	.def("set_translation", (void (EMAN::EMData::*)(float, float, float) )&EMAN::EMData::set_translation)
	.def("get_transform", &EMAN::EMData::get_transform)
	.def("set_rotation", (void (EMAN::EMData::*)(float, float, float) )&EMAN::EMData::set_rotation)
	.def("set_rotation", (void (EMAN::EMData::*)(const EMAN::Transform3D&) )&EMAN::EMData::set_rotation)
	.def("set_size", &EMAN::EMData::set_size, EMAN_EMData_set_size_overloads_1_3())
	.def("set_complex_size", &EMAN::EMData::set_complex_size, EMAN_EMData_set_complex_size_overloads_1_3())
	.def("set_path", &EMAN::EMData::set_path)
	.def("set_pathnum", &EMAN::EMData::set_pathnum)
#ifdef EMAN2_USING_OPENGL
	.def("gen_glu_mipmaps", (unsigned int (EMAN::EMData::*)() const)&EMAN::EMData::gen_glu_mipmaps)
	.def("gen_gl_texture", (unsigned int (EMAN::EMData::*)() const)&EMAN::EMData::gen_gl_texture)
// 	.def("get_data_void_pointer", (EMAN::Dict (EMAN::EMData::*)() const)&EMAN::EMData::get_data_void_pointer) removed by d.woolford - we never got the idea off the ground
	.def("render_amp8_gl_texture", &EMAN::EMData::render_amp8_gl_texture)
#endif //EMAN2_USING_OPENGL
	.def("get_2dview", (EMAN::MArray2D (EMAN::EMData::*)() const)&EMAN::EMData::get_2dview)
	.def("get_3dview", (EMAN::MArray3D (EMAN::EMData::*)() const)&EMAN::EMData::get_3dview)
	.def("get_2dcview", (EMAN::MCArray2D (EMAN::EMData::*)() const)&EMAN::EMData::get_2dcview)
	.def("get_3dcview", (EMAN::MCArray3D (EMAN::EMData::*)() const)&EMAN::EMData::get_3dcview)
	.def("get_3dcviewptr", &EMAN::EMData::get_3dcviewptr, return_value_policy< reference_existing_object >())
	.def("get_2dview", (EMAN::MArray2D (EMAN::EMData::*)(int, int) const)&EMAN::EMData::get_2dview)
	.def("get_3dview", (EMAN::MArray3D (EMAN::EMData::*)(int, int, int) const)&EMAN::EMData::get_3dview)
	.def("get_2dcview", (EMAN::MCArray2D (EMAN::EMData::*)(int, int) const)&EMAN::EMData::get_2dcview)
	.def("get_3dcview", (EMAN::MCArray3D (EMAN::EMData::*)(int, int, int) const)&EMAN::EMData::get_3dcview)
	.def("get_attr", &EMAN::EMData::get_attr)
	.def("get_attr_default", &EMAN::EMData::get_attr_default, EMAN_EMData_get_attr_default_overloads_1_2())
	.def("set_attr", &EMAN::EMData::set_attr_python)
//	.def("set_attr", &EMAN::EMData::set_attr)
	.def("get_attr_dict", &EMAN::EMData::get_attr_dict)
	.def("set_attr_dict", &EMAN::EMData::set_attr_dict)
	.def("del_attr", &EMAN::EMData::del_attr)
	.def("has_attr", &EMAN::EMData::has_attr)
	.def("del_attr_dict", &EMAN::EMData::del_attr_dict)
	.def("debug_print_parms", &EMAN::EMData::debug_print_parms)
	.def("get_xsize", &EMAN::EMData::get_xsize)
	.def("get_ysize", &EMAN::EMData::get_ysize)
	.def("get_zsize", &EMAN::EMData::get_zsize)
	.def("get_size", &EMAN::EMData::get_size)
	.def("get_data_as_vector", &EMAN::EMData::get_data_as_vector)
	.def("get_ndim", &EMAN::EMData::get_ndim)
	.def("is_shuffled", &EMAN::EMData::is_shuffled)
	.def("is_FH", &EMAN::EMData::is_FH)
	.def("is_complex", &EMAN::EMData::is_complex)
	.def("is_real", &EMAN::EMData::is_real)
	.def("set_shuffled", &EMAN::EMData::set_shuffled)
	.def("set_FH", &EMAN::EMData::set_FH)
	.def("set_complex", &EMAN::EMData::set_complex)
	.def("is_complex_x", &EMAN::EMData::is_complex_x)
	.def("set_complex_x", &EMAN::EMData::set_complex_x)
	.def("is_flipped", &EMAN::EMData::is_flipped)
	.def("set_flipped", &EMAN::EMData::set_flipped)
	.def("is_ri", &EMAN::EMData::is_ri)
	.def("set_ri", &EMAN::EMData::set_ri)
	.def("is_fftpadded", &EMAN::EMData::is_fftpadded)
	.def("set_fftpad", &EMAN::EMData::set_fftpad)
	.def("is_fftodd", &EMAN::EMData::is_fftodd)
	.def("set_fftodd", &EMAN::EMData::set_fftodd)
	.def("set_nxc", &EMAN::EMData::set_nxc)
	.def("get_flags", &EMAN::EMData::get_flags)
	.def("set_flags", &EMAN::EMData::set_flags)
	.def("get_changecount", &EMAN::EMData::get_changecount)
	.def("set_changecount", &EMAN::EMData::set_changecount)
	.def("get_xoff", &EMAN::EMData::get_xoff)
	.def("get_yoff", &EMAN::EMData::get_yoff)
	.def("get_zoff", &EMAN::EMData::get_zoff)
	.def("set_xyzoff", &EMAN::EMData::set_xyzoff)
	.def("get_path", &EMAN::EMData::get_path)
	.def("get_pathnum", &EMAN::EMData::get_pathnum)
	.def("get_data_pickle", &EMAN::EMData::get_data_pickle)
	.def("set_data_pickle", &EMAN::EMData::set_data_pickle)
	.def("get_supp_pickle", &EMAN::EMData::get_supp_pickle)
	.def("set_supp_pickle", &EMAN::EMData::set_supp_pickle)
	.def("write_data",&EMAN::EMData::write_data)
	.def("read_data",&EMAN::EMData::read_data,EMAN_EMData_read_data_2_3())
	.def("process_inplace", (void (EMAN::EMData::*)(const std::string&, const EMAN::Dict&) )&EMAN::EMData::process_inplace, EMAN_EMData_process_inplace_overloads_1_2())
	.def("process_inplace", (void (EMAN::EMData::*)(EMAN::Processor*) )&EMAN::EMData::process_inplace)
	.def("process", (EMAN::EMData* (EMAN::EMData::*)(const std::string&, const EMAN::Dict&) const )&EMAN::EMData::process, EMAN_EMData_process_overloads_1_2()[ return_value_policy< manage_new_object >() ])
	.def("process", (EMAN::EMData* (EMAN::EMData::*)(EMAN::Processor*) const )&EMAN::EMData::process, return_value_policy< manage_new_object >())
	.def("cmp", &EMAN::EMData::cmp, EMAN_EMData_cmp_overloads_2_3())
	.def("align", &EMAN::EMData::align, EMAN_EMData_align_overloads_2_5()[ return_value_policy< manage_new_object >() ])
	.def("project", (EMAN::EMData* (EMAN::EMData::*)(const std::string&, const EMAN::Dict&) )&EMAN::EMData::project, EMAN_EMData_project_overloads_1_2()[ return_value_policy< manage_new_object >() ])
	.def("project", (EMAN::EMData* (EMAN::EMData::*)(const std::string&, const EMAN::Transform&) )&EMAN::EMData::project, return_value_policy< manage_new_object >() )
	.def("backproject", &EMAN::EMData::backproject, EMAN_EMData_backproject_overloads_1_2()[ return_value_policy< manage_new_object >() ])
	.def("do_fft", &EMAN::EMData::do_fft, return_value_policy< manage_new_object >())
	.def("do_fft_inplace", &EMAN::EMData::do_fft_inplace, return_value_policy< reference_existing_object >())
	//.def("do_fft_inplace", &EMAN::EMData::do_fft_inplace)
	.def("do_ift", &EMAN::EMData::do_ift, return_value_policy< manage_new_object >())
	.def("do_ift_inplace", &EMAN::EMData::do_ift_inplace, return_value_policy< reference_existing_object >())
	.def("bispecRotTransInvN", &EMAN::EMData::bispecRotTransInvN, return_value_policy< reference_existing_object >())
	.def("bispecRotTransInvDirect", &EMAN::EMData::bispecRotTransInvDirect, return_value_policy< reference_existing_object >())
#ifdef EMAN2_USING_CUDA
	.def("do_fft_cuda", &EMAN::EMData::do_fft_cuda, return_value_policy< manage_new_object >())
	.def("do_ift_cuda", &EMAN::EMData::do_ift_cuda, EMAN_EMData_do_ift_cuda_overloads_0_1()[ return_value_policy< manage_new_object >() ])
	.def("calc_ccf_cuda", &EMAN::EMData::calc_ccf_cuda, return_value_policy< manage_new_object >())
	.def("cut_slice_cuda", &EMAN::EMData::cut_slice_cuda, return_value_policy< manage_new_object >())
	.def("set_gpu_rw_current", &EMAN::EMData::set_gpu_rw_current)
	.def("column_sum_cuda",&EMAN::EMData::column_sum_cuda,return_value_policy< manage_new_object >() )
	.def("make_rotational_footprint_cuda", &EMAN::EMData::make_rotational_footprint_cuda, EMAN_EMData_make_rotational_footprint_cuda_overloads_0_1()[ return_value_policy< manage_new_object >() ])
	// These ones are currently meant mainly for testing purposes
	.def("_copy_gpu_rw_to_cpu", &EMAN::EMData::copy_gpu_rw_to_cpu)
	.def("_copy_cpu_to_gpu_rw", &EMAN::EMData::copy_cpu_to_gpu_rw)
	.def("_copy_cpu_to_gpu_ro", &EMAN::EMData::copy_cpu_to_gpu_ro)
	.def("_copy_gpu_rw_to_gpu_ro", &EMAN::EMData::copy_gpu_rw_to_gpu_ro)
	.def("_copy_gpu_ro_to_gpu_rw", &EMAN::EMData::copy_gpu_ro_to_gpu_rw)
	.def("_copy_gpu_ro_to_cpu", &EMAN::EMData::copy_gpu_ro_to_cpu)
	.def("print_this", &EMAN::EMData::print_this)
	.def("cuda_lock", &EMAN::EMData::cuda_lock)
	.def("cuda_unlock", &EMAN::EMData::cuda_unlock)
#endif // EMAN2_USING_CUDA
	.def("render_amp8", &EMAN::EMData::render_amp8)
	.def("render_ap24", &EMAN::EMData::render_ap24)
	.def("ri2ap", &EMAN::EMData::ri2ap)
	.def("ap2ri", &EMAN::EMData::ap2ri)
	.def("ri2inten", &EMAN::EMData::ri2inten)
	.def("insert_clip", &EMAN::EMData::insert_clip)
	.def("insert_scaled_sum", &EMAN::EMData::insert_scaled_sum, EMAN_EMData_insert_scaled_sum_overloads_2_4())
	.def("copy", &EMAN::EMData::copy, return_value_policy< manage_new_object >())
	.def("copy_head", &EMAN::EMData::copy_head, return_value_policy< manage_new_object >())
	.def("add", (void (EMAN::EMData::*)(float, int) )&EMAN::EMData::add, EMAN_EMData_add_overloads_1_2())
	.def("add", (void (EMAN::EMData::*)(const EMAN::EMData&) )&EMAN::EMData::add)
	.def("addsquare", (void (EMAN::EMData::*)(const EMAN::EMData&) )&EMAN::EMData::addsquare)
	.def("sub", (void (EMAN::EMData::*)(float) )&EMAN::EMData::sub)
	.def("sub", (void (EMAN::EMData::*)(const EMAN::EMData&) )&EMAN::EMData::sub)
	.def("subsquare", (void (EMAN::EMData::*)(const EMAN::EMData&) )&EMAN::EMData::subsquare)
	.def("mult", (void (EMAN::EMData::*)(const EMAN::EMData&, bool) )&EMAN::EMData::mult, EMAN_EMData_mult_overloads_1_2())
	.def("mult", (void (EMAN::EMData::*)(int) )&EMAN::EMData::mult)
	.def("mult", (void (EMAN::EMData::*)(float) )&EMAN::EMData::mult)
	.def("div", (void (EMAN::EMData::*)(float) )&EMAN::EMData::div)
	.def("div", (void (EMAN::EMData::*)(const EMAN::EMData&) )&EMAN::EMData::div)
	.def("to_zero", &EMAN::EMData::to_zero)
	.def("to_one", &EMAN::EMData::to_one)
	.def("dot", &EMAN::EMData::dot)
	.def("get_row", &EMAN::EMData::get_row, return_value_policy< manage_new_object >())
	.def("set_row", &EMAN::EMData::set_row)
	.def("get_col", &EMAN::EMData::get_col, return_value_policy< manage_new_object >())
	.def("set_col", &EMAN::EMData::set_col)
	.def("get_value_at", (float (EMAN::EMData::*)(int, int, int) const)&EMAN::EMData::get_value_at)
	.def("get_value_at", (float (EMAN::EMData::*)(int, int) const)&EMAN::EMData::get_value_at)
	.def("get_value_at", (float (EMAN::EMData::*)(size_t) const)&EMAN::EMData::get_value_at)
	.def("sget_value_at", (float (EMAN::EMData::*)(int, int, int) const)&EMAN::EMData::sget_value_at)
	.def("sget_value_at", (float (EMAN::EMData::*)(int, int) const)&EMAN::EMData::sget_value_at)
	.def("sget_value_at", (float (EMAN::EMData::*)(size_t) const)&EMAN::EMData::sget_value_at)
	.def("sget_value_at_interp", (float (EMAN::EMData::*)(float, float) const)&EMAN::EMData::sget_value_at_interp)
	.def("sget_value_at_interp", (float (EMAN::EMData::*)(float, float, float) const)&EMAN::EMData::sget_value_at_interp)
	.def("set_value_at", (void (EMAN::EMData::*)(int, int, int, float) )&EMAN::EMData::set_value_at)
	.def("set_value_at_fast", (void (EMAN::EMData::*)(int, int, int, float) )&EMAN::EMData::set_value_at_fast)
	.def("set_value_at", (void (EMAN::EMData::*)(int, int, float) )&EMAN::EMData::set_value_at)
	.def("set_value_at_fast", (void (EMAN::EMData::*)(int, int, float) )&EMAN::EMData::set_value_at_fast)
	.def("set_value_at", (void (EMAN::EMData::*)(int, float) )&EMAN::EMData::set_value_at)
#ifndef	_WIN32
	.def("set_array_offsets", (void (EMAN::EMData::*)(const int, const int, const int) )&EMAN::EMData::set_array_offsets, EMAN_EMData_set_array_offsets_overloads_0_3())
	.def("set_array_offsets", (void (EMAN::EMData::*)(std::vector<int,std::allocator<int> >) )&EMAN::EMData::set_array_offsets)
#endif
	.def("get_array_offsets", &EMAN::EMData::get_array_offsets)
	.def("power", &EMAN::EMData::power, return_value_policy< manage_new_object >())
	.def("sqrt", &EMAN::EMData::sqrt, return_value_policy< manage_new_object >())
	.def("log", &EMAN::EMData::log, return_value_policy< manage_new_object >())
	.def("log10", &EMAN::EMData::log10, return_value_policy< manage_new_object >())
	.def("real", &EMAN::EMData::real, return_value_policy< manage_new_object >())
	.def("imag", &EMAN::EMData::imag, return_value_policy< manage_new_object >())
	.def("absi", &EMAN::EMData::absi, return_value_policy< manage_new_object >())
	.def("amplitude", &EMAN::EMData::amplitude, return_value_policy< manage_new_object >())
	.def("phase", &EMAN::EMData::phase, return_value_policy< manage_new_object >())
	.def("real2complex", &EMAN::EMData::real2complex, EMAN_EMData_real2complex_overloads_0_1()[ return_value_policy< manage_new_object >() ])
	.def("real2FH", &EMAN::EMData::real2FH, return_value_policy< manage_new_object >())
	.def("FH2F", &EMAN::EMData::FH2F, EMAN_EMData_FH2F_overloads_2_3()[ return_value_policy< manage_new_object >() ])
	.def("FH2Real", &EMAN::EMData::FH2Real, EMAN_EMData_FH2Real_overloads_2_3()[ return_value_policy< manage_new_object >() ])
	.def("rotavg", &EMAN::EMData::rotavg, return_value_policy< manage_new_object >())
	.def("rotavg_i", &EMAN::EMData::rotavg_i, return_value_policy< manage_new_object >())
	.def("mult_radial", &EMAN::EMData::mult_radial, return_value_policy< manage_new_object >())
	.def("cog", &EMAN::EMData::cog)
	.def("calc_fourier_shell_correlation", &EMAN::EMData::calc_fourier_shell_correlation, EMAN_EMData_calc_fourier_shell_correlation_overloads_1_2())
	.def("average_circ_sub", &EMAN::EMData::average_circ_sub, return_value_policy< manage_new_object >())
	.def("onelinenn", &EMAN::EMData::onelinenn)
	.def("onelinenn_mult", &EMAN::EMData::onelinenn_mult)
	.def("nn", &EMAN::EMData::nn, EMAN_EMData_nn_overloads_3_4())
	.def("nn_SSNR", &EMAN::EMData::nn_SSNR, EMAN_EMData_nn_SSNR_overloads_4_5())
	.def("nn_SSNR_ctf", &EMAN::EMData::nn_SSNR_ctf, EMAN_EMData_nn_SSNR_ctf_overloads_5_6())
	.def("symplane0", &EMAN::EMData::symplane0)
	.def("symplane1", &EMAN::EMData::symplane1)
	.def("symplane2", &EMAN::EMData::symplane2)
	.def("onelinenn_ctf", &EMAN::EMData::onelinenn_ctf)
	.def("nn_ctf", &EMAN::EMData::nn_ctf)
	.def("onelinenn_ctf_applied", &EMAN::EMData::onelinenn_ctf_applied)
	.def("nn_ctf_applied", &EMAN::EMData::nn_ctf_applied)
	.def("symplane0_ctf", &EMAN::EMData::symplane0_ctf)
	.def("symvol", &EMAN::EMData::symvol, return_value_policy< manage_new_object >())
	.def("rot_scale_trans2D", &EMAN::EMData::rot_scale_trans2D, EMAN_EMData_rot_scale_trans2D_overloads_1_4()[ return_value_policy< manage_new_object >() ])
	.def("rot_scale_trans", &EMAN::EMData::rot_scale_trans, return_value_policy< manage_new_object >())
	.def("cm_euc", &EMAN::EMData::cm_euc)
	.def("rot_scale_conv", &EMAN::EMData::rot_scale_conv, EMAN_EMData_rot_scale_conv_overloads_4_5()[ return_value_policy< manage_new_object >() ])
	.def("rot_scale_conv7", &EMAN::EMData::rot_scale_conv7, return_value_policy< manage_new_object >())
	.def("rot_scale_conv_new", &EMAN::EMData::rot_scale_conv_new, EMAN_EMData_rot_scale_conv_new_overloads_4_5()[ return_value_policy< manage_new_object >() ])
	.def("downsample", &EMAN::EMData::downsample, EMAN_EMData_downsample_overloads_1_2()[ return_value_policy< manage_new_object >() ])
	.def("get_pixel_conv", &EMAN::EMData::get_pixel_conv)
	.def("get_pixel_conv7", &EMAN::EMData::get_pixel_conv7)
	.def("getconvpt2d_kbi0", &EMAN::EMData::getconvpt2d_kbi0, EMAN_EMData_getconvpt2d_kbi0_overloads_3_4())
	.def("fft_shuffle", &EMAN::EMData::fft_shuffle)
	.def("extractpoint", &EMAN::EMData::extractpoint)
	.def("extract_plane", &EMAN::EMData::extract_plane, return_value_policy< manage_new_object >())
	.def("fouriergridrot2d", &EMAN::EMData::fouriergridrot2d, return_value_policy< manage_new_object >())
	.def("fouriergridrot_shift2d", &EMAN::EMData::fouriergridrot_shift2d, return_value_policy< manage_new_object >())
	.def("delete_disconnected_regions", &EMAN::EMData::delete_disconnected_regions, return_value_policy< manage_new_object >())
	.def("helicise", &EMAN::EMData::helicise, return_value_policy< manage_new_object >())
	.def("divkbsinh", &EMAN::EMData::divkbsinh)
	.def("peak_search", &EMAN::EMData::peak_search)
	.def("phase_cog", &EMAN::EMData::phase_cog)
	.def("find_3d_threshold", &EMAN::EMData::find_3d_threshold)
	.def("peak_ccf", &EMAN::EMData::peak_ccf)
	.def("debug_print_params", &EMAN::EMData::debug_print_parms)
	.def("get_pow", &EMAN::EMData::get_pow, return_value_policy< manage_new_object >())
	.def("conjg", &EMAN::EMData::conjg, return_value_policy< manage_new_object >())
	.def("extractline", &EMAN::EMData::extractline, return_value_policy< manage_new_object >())
	.def("center_origin", &EMAN::EMData::center_origin)
	.def("center_origin_yz", &EMAN::EMData::center_origin_yz)
	.def("center_origin_fft", &EMAN::EMData::center_origin_fft)
	.def("depad", &EMAN::EMData::depad)
	.def("depad_corner", &EMAN::EMData::depad_corner)
	.def("FourInterpol", &EMAN::EMData::FourInterpol, EMAN_EMData_FourInterpol_overloads_1_4()[ return_value_policy< manage_new_object >() ])
	.def("FourTruncate", &EMAN::EMData::FourTruncate, EMAN_EMData_FourTruncate_overloads_1_4()[ return_value_policy< manage_new_object >() ])
//       .def("FourInterpol_i", &EMAN::EMData::FourInterpol_i, EMAN_EMData_FourInterpol_i_overloads_1_4()[ return_value_policy< manage_new_object >() ])
	.def("Four_ds", &EMAN::EMData::Four_ds, EMAN_EMData_Four_ds_overloads_1_4()[ return_value_policy< manage_new_object >() ])
	.def("Four_shuf_ds_cen_us", &EMAN::EMData::Four_shuf_ds_cen_us, EMAN_EMData_Four_shuf_ds_cen_us_overloads_1_4()[ return_value_policy< manage_new_object >() ])
	.def("filter_by_image", &EMAN::EMData::filter_by_image, EMAN_EMData_filter_by_image_overloads_1_2()[ return_value_policy< manage_new_object >() ])
	.def("replace_amplitudes", &EMAN::EMData::replace_amplitudes, EMAN_EMData_replace_amplitudes_overloads_1_2()[ return_value_policy< manage_new_object >() ])
	.def("norm_pad", &EMAN::EMData::norm_pad, EMAN_EMData_norm_pad_overloads_2_3() [ return_value_policy< manage_new_object >()])
	.def("get_clip", &EMAN::EMData::get_clip, EMAN_EMData_get_clip_overloads_1_2() [ return_value_policy< manage_new_object >()])
	.def("clip_inplace", &EMAN::EMData::clip_inplace, EMAN_EMData_clip_inplace_overloads_1_2()[return_value_policy< reference_existing_object >()])
	.def("get_top_half", &EMAN::EMData::get_top_half, return_value_policy< manage_new_object >())
	.def("get_rotated_clip", &EMAN::EMData::get_rotated_clip, EMAN_EMData_get_rotated_clip_overloads_2_3()[ return_value_policy< manage_new_object >() ])
	.def("window_center", &EMAN::EMData::window_center, return_value_policy< manage_new_object >())
	.def("scale", &EMAN::EMData::scale)
	.def("zero_corner_circulant", &EMAN::EMData::zero_corner_circulant)
	.def("translate", (void (EMAN::EMData::*)(float, float, float) )&EMAN::EMData::translate)
	.def("translate", (void (EMAN::EMData::*)(const EMAN::Vec3f&) )&EMAN::EMData::translate)
	.def("translate", (void (EMAN::EMData::*)(int, int, int) )&EMAN::EMData::translate)
	.def("translate", (void (EMAN::EMData::*)(const EMAN::Vec3i&) )&EMAN::EMData::translate)
	.def("rotate", (void (EMAN::EMData::*)(const EMAN::Transform3D&) )&EMAN::EMData::rotate)
	.def("rotate", (void (EMAN::EMData::*)(float, float, float) )&EMAN::EMData::rotate)
	.def("rotate_translate", (void (EMAN::EMData::*)(const EMAN::Transform3D&) )&EMAN::EMData::rotate_translate)
	.def("transform", &EMAN::EMData::transform)
	.def("rotate_translate", (void (EMAN::EMData::*)(const EMAN::Transform&) )&EMAN::EMData::rotate_translate)
	.def("rotate_translate", (void (EMAN::EMData::*)(float, float, float, float, float, float) )&EMAN::EMData::rotate_translate)
	.def("rotate_translate", (void (EMAN::EMData::*)(float, float, float, float, float, float, float, float, float) )&EMAN::EMData::rotate_translate)
	.def("rotate_x", &EMAN::EMData::rotate_x)
	.def("rotate_180", &EMAN::EMData::rotate_180)
	.def("dot_rotate_translate", &EMAN::EMData::dot_rotate_translate)
	.def("little_big_dot", &EMAN::EMData::little_big_dot, EMAN_EMData_little_big_dot_overloads_1_2()[ return_value_policy< manage_new_object >() ])
	.def("do_radon", &EMAN::EMData::do_radon, return_value_policy< manage_new_object >())
	.def("calc_ccf", &EMAN::EMData::calc_ccf, EMAN_EMData_calc_ccf_overloads_1_3()[ return_value_policy< manage_new_object >() ])
	.def("calc_ccfx", &EMAN::EMData::calc_ccfx, EMAN_EMData_calc_ccfx_overloads_1_4()[ return_value_policy< manage_new_object >() ])
	.def("calc_fast_sigma_image",&EMAN::EMData::calc_fast_sigma_image, return_value_policy< manage_new_object >())
	.def("make_rotational_footprint", &EMAN::EMData::make_rotational_footprint, EMAN_EMData_make_rotational_footprint_overloads_0_1()[ return_value_policy< manage_new_object >() ])
	.def("make_rotational_footprint_e1", &EMAN::EMData::make_rotational_footprint_e1, EMAN_EMData_make_rotational_footprint_e1_overloads_0_1()[ return_value_policy< manage_new_object >() ])
	.def("make_rotational_footprint_cmc", &EMAN::EMData::make_rotational_footprint_cmc, EMAN_EMData_make_rotational_footprint_cmc_overloads_0_1()[ return_value_policy< manage_new_object >() ])
	.def("make_footprint", &EMAN::EMData::make_footprint, return_value_policy< manage_new_object >())
	.def("calc_mutual_correlation", &EMAN::EMData::calc_mutual_correlation, EMAN_EMData_calc_mutual_correlation_overloads_1_3()[ return_value_policy< manage_new_object >() ])
	.def("unwrap", &EMAN::EMData::unwrap, EMAN_EMData_unwrap_overloads_0_7()[ return_value_policy< manage_new_object >() ])
	.def("apply_radial_func", &EMAN::EMData::apply_radial_func, EMAN_EMData_apply_radial_func_overloads_3_4())
	.def("calc_radial_dist", (std::vector<float,std::allocator<float> > (EMAN::EMData::*)(int, float, float, bool) )&EMAN::EMData::calc_radial_dist)
	.def("calc_radial_dist", (std::vector<float,std::allocator<float> > (EMAN::EMData::*)(int, float, float, int, bool) )&EMAN::EMData::calc_radial_dist)
	.def("cconj", &EMAN::EMData::cconj)
	.def("add_incoherent", &EMAN::EMData::add_incoherent)
	.def("calc_hist", &EMAN::EMData::calc_hist, EMAN_EMData_calc_hist_overloads_0_3())
	.def("calc_az_dist", &EMAN::EMData::calc_az_dist)
	.def("calc_dist", &EMAN::EMData::calc_dist, EMAN_EMData_calc_dist_overloads_1_2())
	.def("calc_flcf", &EMAN::EMData::calc_flcf, return_value_policy< manage_new_object >())
	.def("convolute", &EMAN::EMData::convolute, return_value_policy< manage_new_object >())
	.def("common_lines", &EMAN::EMData::common_lines, EMAN_EMData_common_lines_overloads_2_5())
	.def("common_lines_real", &EMAN::EMData::common_lines_real, EMAN_EMData_common_lines_real_overloads_2_4())
	.def("cut_slice", &EMAN::EMData::cut_slice, EMAN_EMData_cut_slice_overloads_2_3())
	.def("mask_contig_region",&EMAN::EMData::mask_contig_region)
	.def("uncut_slice", &EMAN::EMData::uncut_slice)
	.def("set_xyz_origin", &EMAN::EMData::set_xyz_origin)
	.def("__getitem__", &emdata_getitem)
	.def("__setitem__", &emdata_setitem)
	.staticmethod("read_images_ext")
	.staticmethod("read_images")
	.def("__add__", (EMAN::EMData* (*)(const EMAN::EMData&, const EMAN::EMData&) )&EMAN::operator+, return_value_policy< manage_new_object >() )
	.def("__sub__", (EMAN::EMData* (*)(const EMAN::EMData&, const EMAN::EMData&) )&EMAN::operator-, return_value_policy< manage_new_object >() )
	.def("__mul__", (EMAN::EMData* (*)(const EMAN::EMData&, const EMAN::EMData&) )&EMAN::operator*, return_value_policy< manage_new_object >() )
	.def("__div__", (EMAN::EMData* (*)(const EMAN::EMData&, const EMAN::EMData&) )&EMAN::operator/, return_value_policy< manage_new_object >() )
	.def("__add__", (EMAN::EMData* (*)(const EMAN::EMData&, float) )&EMAN::operator+, return_value_policy< manage_new_object >() )
	.def("__sub__", (EMAN::EMData* (*)(const EMAN::EMData&, float) )&EMAN::operator-, return_value_policy< manage_new_object >() )
	.def("__mul__", (EMAN::EMData* (*)(const EMAN::EMData&, float) )&EMAN::operator*, return_value_policy< manage_new_object >() )
	.def("__div__", (EMAN::EMData* (*)(const EMAN::EMData&, float) )&EMAN::operator/, return_value_policy< manage_new_object >() )
        .def("__radd__", (EMAN::EMData* (*)(const EMAN::EMData&, float) )&EMAN::operator+, return_value_policy< manage_new_object >() )
	.def("__rsub__", (EMAN::EMData* (*)(const EMAN::EMData&, float) )&EMAN::rsub, return_value_policy< manage_new_object >() )
        .def("__rmul__", (EMAN::EMData* (*)(const EMAN::EMData&, float) )&EMAN::operator*, return_value_policy< manage_new_object >() )
        .def("__rdiv__", (EMAN::EMData* (*)(const EMAN::EMData&, float) )&EMAN::rdiv, return_value_policy< manage_new_object >() )
        .def( self += other< float >() )
        .def( self -= other< float >() )
        .def( self *= other< float >() )
        .def( self /= other< float >() )
        .def( self += self )
        .def( self -= self )
		.def( self == self )
        .def( self *= self )
        .def( self /= self )
        .def("__call__", (float& (EMAN::EMData::*)(const int, const int, const int) const )&EMAN::EMData::operator (), return_value_policy< copy_non_const_reference >())
#ifndef	_WIN32
        .def("__call__", (float& (EMAN::EMData::*)(const int, const int) const )&EMAN::EMData::operator (), return_value_policy< copy_non_const_reference >())
        .def("__call__", (float& (EMAN::EMData::*)(const int) const )&EMAN::EMData::operator (), return_value_policy< copy_non_const_reference >())
#endif	//_WIN32
	);

	enum_< EMAN::EMData::FFTPLACE >("FFTPLACE")
	    .value("FFT_IN_PLACE", EMAN::EMData::FFT_IN_PLACE)
	    .value("FFT_OUT_OF_PLACE", EMAN::EMData::FFT_OUT_OF_PLACE)
	;


	enum_< EMAN::EMData::WINDOWPLACE >("WINDOWPLACE")
	    .value("WINDOW_OUT_OF_PLACE", EMAN::EMData::WINDOW_OUT_OF_PLACE)
	    .value("WINDOW_IN_PLACE", EMAN::EMData::WINDOW_IN_PLACE)
	;

	delete EMAN_EMData_scope;

}

