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
#include <emdata.h>
#include <emdata_pickle.h>
#include <emdata_wrapitems.h>
#include <emfft.h>
#include <processor.h>
#include <transform.h>
#include <xydata.h>/** return the FFT amplitude which is greater than thres %
 *
 * @exception ImageFormatException If the image is not a complex image.
 * @return The FFT amplitude which is greater than thres %.
 */
float get_amplitude_thres(float thres);

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_read_image_overloads_1_5, read_image, 1, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_read_binedimage_overloads_1_5, read_binedimage, 1, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_write_image_overloads_1_7, write_image, 1, 7)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_append_image_overloads_1_3, append_image, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_write_lst_overloads_1_4, write_lst, 1, 4)

//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_print_image_overloads_0_2, EMAN::EMData::print_image, 0, 2)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_EMData_read_images_overloads_1_3, EMAN::EMData::read_images, 1, 3)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_EMData_read_images_ext_overloads_3_5, EMAN::EMData::read_images_ext, 3, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_set_size_overloads_1_4, EMAN::EMData::set_size, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_set_complex_size_overloads_1_3, EMAN::EMData::set_complex_size, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_get_attr_default_overloads_1_2, EMAN::EMData::get_attr_default, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_clip_inplace_overloads_1_2, EMAN::EMData::clip_inplace, 1, 2)

//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_process_inplace_overloads_1_2, EMAN::EMData::process_inplace, 1, 2)

//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_process_overloads_1_2, EMAN::EMData::process, 1, 2)

//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_cmp_overloads_2_3, EMAN::EMData::cmp, 2, 3)

//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_align_overloads_2_5, EMAN::EMData::align, 2, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_xform_align_nbest_overloads_2_6, EMAN::EMData::xform_align_nbest, 2, 6)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_project_overloads_1_2, EMAN::EMData::project, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_backproject_overloads_1_2, EMAN::EMData::backproject, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_insert_scaled_sum_overloads_2_4, EMAN::EMData::insert_scaled_sum, 2, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_add_overloads_1_2, EMAN::EMData::add, 1, 2)

#ifndef	_WIN32
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_set_array_offsets_overloads_0_3, EMAN::EMData::set_array_offsets, 0, 3)
#endif

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_real2complex_overloads_0_1, EMAN::EMData::real2complex, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_FH2F_overloads_2_3, EMAN::EMData::FH2F, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_FH2Real_overloads_2_3, EMAN::EMData::FH2Real, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_fourier_shell_correlation_overloads_1_2, EMAN::EMData::calc_fourier_shell_correlation, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_nn_overloads_3_4, EMAN::EMData::nn, 3, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_nn_SSNR_overloads_4_5, EMAN::EMData::nn_SSNR, 4, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_nn_SSNR_ctf_overloads_5_6, EMAN::EMData::nn_SSNR_ctf, 5, 6)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_rot_scale_trans2D_overloads_1_4, EMAN::EMData::rot_scale_trans2D, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_rot_scale_conv_overloads_4_5, EMAN::EMData::rot_scale_conv, 4, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_rot_scale_conv_new_overloads_4_5, EMAN::EMData::rot_scale_conv_new, 4, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_rot_scale_conv_new_3D_overloads_7_9, EMAN::EMData::rot_scale_conv_new_3D, 7, 9)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_downsample_overloads_1_2, EMAN::EMData::downsample, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_getconvpt2d_kbi0_overloads_3_4, EMAN::EMData::getconvpt2d_kbi0, 3, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_FourInterpol_overloads_1_4, EMAN::EMData::FourInterpol, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_FourTruncate_overloads_1_4, EMAN::EMData::FourTruncate, 1, 4)

//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_FourInterpol_i_overloads_1_4, FourInterpol_i, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_Four_ds_overloads_1_4, EMAN::EMData::Four_ds, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_Four_shuf_ds_cen_us_overloads_1_4, EMAN::EMData::Four_shuf_ds_cen_us, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_filter_by_image_overloads_1_2, EMAN::EMData::filter_by_image, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_replace_amplitudes_overloads_1_2, EMAN::EMData::replace_amplitudes, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_get_rotated_clip_overloads_2_3, EMAN::EMData::get_rotated_clip, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_little_big_dot_overloads_1_2, EMAN::EMData::little_big_dot, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_ccf_overloads_0_3, EMAN::EMData::calc_ccf, 0, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_ccfx_overloads_1_4, EMAN::EMData::calc_ccfx, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_make_rotational_footprint_overloads_0_1, EMAN::EMData::make_rotational_footprint, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_make_rotational_footprint_e1_overloads_0_1, EMAN::EMData::make_rotational_footprint_e1, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_make_rotational_footprint_cmc_overloads_0_1, EMAN::EMData::make_rotational_footprint_cmc, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_mutual_correlation_overloads_1_3, EMAN::EMData::calc_mutual_correlation, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_unwrap_overloads_0_7, EMAN::EMData::unwrap, 0, 7)


//#ifdef EMAN2_USING_CUDA
//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_unwrap_cuda_overloads_0_6, unwrap_cuda, 0, 6)
//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_get_clip_cuda_overloads_1_2, get_clip_cuda, 1, 2)
//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_make_rotational_footprint_cuda_overloads_0_1, EMAN::EMData::make_rotational_footprint_cuda, 0, 1)
//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_do_ift_cuda_overloads_0_1, do_ift_cuda, 0, 1)
//#endif //EMAN2_USING_CUDA
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_apply_radial_func_overloads_3_4, EMAN::EMData::apply_radial_func, 3, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_hist_overloads_0_5, EMAN::EMData::calc_hist, 0, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_dist_overloads_1_2, EMAN::EMData::calc_dist, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_common_lines_overloads_2_5, EMAN::EMData::common_lines, 2, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_common_lines_real_overloads_2_4, EMAN::EMData::common_lines_real, 2, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_cut_slice_overloads_2_3, EMAN::EMData::cut_slice, 2, 3)

// BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_uncut_slice_overloads_1_2, uncut_slice, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_get_clip_overloads_1_2, EMAN::EMData::get_clip, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_mult_overloads_1_2, EMAN::EMData::mult, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_norm_pad_overloads_2_3, EMAN::EMData::norm_pad, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_read_data_overloads_2_6, EMAN::EMData::read_data, 2, 6)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_write_data_overloads_2_6, EMAN::EMData::write_data, 2, 6)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_rot_scale_trans2D_background_overloads_1_4, EMAN::EMData::rot_scale_trans2D_background, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_rot_scale_conv_new_background_overloads_4_5, EMAN::EMData::rot_scale_conv_new_background, 4, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_rot_scale_conv_new_background_twice_overloads_4_5, EMAN::EMData::rot_scale_conv_new_background_twice, 4, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_rot_scale_conv_new_background_3D_overloads_7_9, EMAN::EMData::rot_scale_conv_new_background_3D, 7, 9)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_delete_disconnected_regions_overloads_0_3, EMAN::EMData::delete_disconnected_regions, 0, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_helicise_overloads_3_6, EMAN::EMData::helicise, 3, 6)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_helicise_grid_overloads_4_7, EMAN::EMData::helicise_grid, 4, 7)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_zero_corner_circulant_overloads_0_1, EMAN::EMData::zero_corner_circulant, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_dot_rotate_translate_overloads_4_5, EMAN::EMData::dot_rotate_translate, 4, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_make_footprint_overloads_0_1, EMAN::EMData::make_footprint, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_compute_missingwedge_overloads_1_3, EMAN::EMData::compute_missingwedge, 1, 3)

}// namespace

using namespace EMAN;

//// These give us threadsafety. Couldn't find a more elegant way to do it with overloading :^/
EMData *EMData_align_wrapper2(EMData &ths, const string & aligner_name, EMData * to_img) {
	EMData *ret;
	PyThreadState *_save = PyEval_SaveThread();

	try {
		ret = ths.align(aligner_name,to_img);
	}
	catch (std::exception &e) {
		PyEval_RestoreThread(_save);
		cerr << e.what() << endl;
		throw e;
	}
	PyEval_RestoreThread(_save);
	return ret;
}
EMData *EMData_align_wrapper3(EMData &ths, const string & aligner_name, EMData * to_img,const Dict & params) {
	EMData *ret;
	PyThreadState *_save = PyEval_SaveThread();

	try {
		ret=ths.align(aligner_name,to_img,params);
	}
	catch (std::exception &e) {
		PyEval_RestoreThread(_save);
		cerr << e.what() << endl;
		throw e;
	}
	PyEval_RestoreThread(_save);
	return ret;
}
EMData *EMData_align_wrapper4(EMData &ths, const string & aligner_name, EMData * to_img,const Dict & params, const string & cmp_name) {
	EMData *ret;
	PyThreadState *_save = PyEval_SaveThread();

	try {
		ret=ths.align(aligner_name,to_img,params,cmp_name);
	}
	catch (std::exception &e) {
		PyEval_RestoreThread(_save);
		cerr << e.what() << endl;
		throw e;
	}
	PyEval_RestoreThread(_save);
	return ret;
}
EMData *EMData_align_wrapper5(EMData &ths, const string & aligner_name, EMData * to_img,const Dict & params, const string & cmp_name, const Dict& cmp_params) {
	EMData *ret;
	PyThreadState *_save = PyEval_SaveThread();

	try {
		ret=ths.align(aligner_name,to_img,params,cmp_name,cmp_params);
	}
	catch (std::exception &e) {
		PyEval_RestoreThread(_save);
		cerr << e.what() << endl;
		throw e;
	}
	PyEval_RestoreThread(_save);
	return ret;
}

void EMData_process_inplace_wrapper1(EMData &ths,const string & processorname) {
	PyThreadState *_save = PyEval_SaveThread();

	try {
		ths.process_inplace(processorname);
	}
	catch (std::exception &e) {
		PyEval_RestoreThread(_save);
		cerr << e.what() << endl;
		throw e;
	}
	PyEval_RestoreThread(_save);
}

void EMData_process_inplace_wrapper2(EMData &ths,const string & processorname, const Dict & params) {
	PyThreadState *_save = PyEval_SaveThread();

	try {
		ths.process_inplace(processorname,params);
	}
	catch (std::exception &e) {
		PyEval_RestoreThread(_save);
		cerr << e.what() << endl;
		throw e;
	}
	PyEval_RestoreThread(_save);
}

EMData *EMData_process_wrapper1(EMData &ths,const string & processorname) {
	EMData *ret;
	PyThreadState *_save = PyEval_SaveThread();

	try {
		ret=ths.process(processorname);
	}
	catch (std::exception &e) {
		PyEval_RestoreThread(_save);
		cerr << e.what() << endl;
		throw e;
	}
	PyEval_RestoreThread(_save);
	return ret;
}
EMData *EMData_process_wrapper2(EMData &ths,const string & processorname, const Dict & params) {
	EMData *ret;
	PyThreadState *_save = PyEval_SaveThread();

	try {
		ret=ths.process(processorname,params);
	}
	catch (std::exception &e) {
		PyEval_RestoreThread(_save);
		cerr << e.what() << endl;
		throw e;
	}
	PyEval_RestoreThread(_save);
	return ret;
}

float EMData_cmp_wrapper2(EMData &ths,const string & cmpname, EMData * with) {
	float ret;
	PyThreadState *_save = PyEval_SaveThread();

	try {
		ret=ths.cmp(cmpname,with);
	}
	catch (std::exception &e) {
		PyEval_RestoreThread(_save);
		cerr << e.what() << endl;
		throw e;
	}
	PyEval_RestoreThread(_save);
	return ret;
}
float EMData_cmp_wrapper3(EMData &ths,const string & cmpname, EMData * with, const Dict & params) {
	float ret;
	PyThreadState *_save = PyEval_SaveThread();

	try {
		ret=ths.cmp(cmpname,with,params);
	}
	catch (std::exception &e) {
		PyEval_RestoreThread(_save);
		cerr << e.what() << endl;
		throw e;
	}
	PyEval_RestoreThread(_save);
	return ret;
}


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyEMData2)
{
    scope* EMAN_EMData_scope = new scope(
    class_< EMAN::EMData >("EMData",
    		"EMData stores an image's data and defines core image processing routines.\n"
    		"The image is 1D, 2D or 3D, in real space or fourier space (complex image).\n"
    		"Data are ordered with x increasing fastest, then y, then z.",
    		init<  >())
//    class_< EMAN::EMData, std::auto_ptr<EMAN::EMData> >("EMData", init<  >())
	.def_pickle(EMData_pickle_suite())
	.def(init< const EMAN::EMData& >(args("that"), "Construct from an EMData (copy constructor).\nPerforms a deep copy.\n \nthat - the EMData to copy"))
	.def(init< const std::string&, optional< int > >(args("filename", "image_index"), "Construct from an image file.\n \nfilename - the image file name\nimage_index the image index for stack image file(default = 0)"))
	.def(init< int, int, optional< int, bool > >(args("nx", "ny", "nz", "is_real"), "makes an image of the specified size, either real or complex.\nFor complex image, the user would specify the real-space dimensions.\n \nnx - size for x dimension\nny - size for y dimension\nnz size for z dimension(default=1)\nis_real - boolean to specify real(true) or complex(false) image(default=True)"))
	.add_static_property("totalalloc", make_getter(EMAN::EMData::totalalloc), make_setter(EMAN::EMData::totalalloc))
	.def("read_image", &EMAN::EMData::read_image, EMAN_EMData_read_image_overloads_1_5(args("filename", "img_index", "header_only", "region", "is_3d"), "read an image file and stores its information to this EMData object.\n\nIf a region is given, then only read a\nregion of the image file. The region will be this\nEMData object. The given region must be inside the given\nimage file. Otherwise, an error will be created.\n\nfilename The image file name.\nimg_index The nth image you want to read.\nheader_only To read only the header or both header and data.\nregion To read only a region of the image.\nis_3d  Whether to treat the image as a single 3D or a set of 2Ds. This is a hint for certain image formats which has no difference between 3D image and set of 2Ds.\nexception ImageFormatException\nexception ImageReadException"))
	.def("read_binedimage", &EMAN::EMData::read_binedimage, EMAN_EMData_read_binedimage_overloads_1_5(args("filename", "img_index", "binfactor", "fast", "is_3d"), "read an image file and stores its information to this EMData object.\nfilename The image file name.\nimg_index The nth image you want to read.\nbinfactor The amount by which to bin by. Must be an integer\nfast bin very binfactor xy slice otherwise meanshrink z slice\nis_3d  Whether to treat the image as a single 3D or a set of 2Ds. This is a hint for certain image formats which has no difference between 3D image and set of 2Ds.\nexception ImageFormatException\nexception ImageReadException"))
	.def("write_image", &EMAN::EMData::write_image, EMAN_EMData_write_image_overloads_1_7(args("filename", "img_index", "imgtype", "header_only", "region", "filestoragetype", "use_host_endian"), "write the header and data out to an image.\n\nIf the img_index = -1, append the image to the given image file.\n\nIf the given image file already exists, this image\nformat only stores 1 image, and no region is given, then\ntruncate the image file  to  zero length before writing\ndata out. For header writing only, no truncation happens.\n\nIf a region is given, then write a region only.\n\nfilename - The image file name.\nimg_index - The nth image to write as.\nimgtype - Write to the given image format type. if not specified, use the 'filename' extension to decide.\nheader_only - To write only the header or both header and data.\nregion - Define the region to write to.\nfilestoragetype - The image data type used in the output file.\nuse_host_endian - To write in the host computer byte order.\n\nexception - ImageFormatException\nexception ImageWriteException"))
	.def("append_image", &EMAN::EMData::append_image, EMAN_EMData_append_image_overloads_1_3(args("filename", "imgtype", "header_only"), "append to an image file; If the file doesn't exist, create one.\nfilename - The image file name.\nimgtype - Write to the given image format type. if not specified, use the 'filename' extension to decide.\nheader_only - To write only the header or both header and data."))
	.def("write_lst", &EMAN::EMData::write_lst, EMAN_EMData_write_lst_overloads_1_4(args("filename", "reffile", "refn", "comment"), "Append data to a LST image file.\nfilename - The LST image file name.\nreffile - Reference file name.\nrefn The reference file number.\ncomment - The comment to the added reference file."))
//	.def("print_image", &EMAN::EMData::print_image, EMAN_EMData_print_image_overloads_0_2(args("filename", "output_stream"), "Print the image data to a file stream (standard out by default).\nfilename - image file to be printed.\noutput_stream - Output stream; cout by default."))
	.def("read_images", &EMAN::EMData::read_images, EMAN_EMData_read_images_overloads_1_3(args("filename", "img_indices", "header_only"),"Read a set of images from file specified by 'filename'.\nWhich images are read is set by 'img_indices'.\nfilename The image file name.\nimg_indices Which images are read. If it is empty, all images are read. If it is not empty, only those in this array are read.\nheader_only If true, only read image header. If false, read both data and header.\nreturn The set of images read from filename."))
	.def("read_images_ext", &EMAN::EMData::read_images_ext, EMAN_EMData_read_images_ext_overloads_3_5(args("filename", "img_index_start", "img_index_end", "header_only", "ext"), "Read a set of images from file specified by 'filename'. If\nthe given 'ext' is not empty, replace 'filename's extension it.\nImages with index from img_index_start to img_index_end are read.\n \nfilename - The image file name.\nimg_index_start Starting image index.\nimg_index_end - Ending image index.\nheader_only - If true, only read image header. If false, read both data and header.\next - The new image filename extension.\n \nreturn The set of images read from filename."))
	.def("get_fft_amplitude", &EMAN::EMData::get_fft_amplitude, return_value_policy< manage_new_object >(), "return the amplitudes of the FFT including the left half\n \nreturn The current FFT image's amplitude image.\nexception - ImageFormatException If the image is not a complex image.")
	.def("get_fft_amplitude2D", &EMAN::EMData::get_fft_amplitude2D, return_value_policy< manage_new_object >(), "return the amplitudes of the 2D FFT including the left half, PRB\n \nreturn The current FFT image's amplitude image.\nexception - ImageFormatException If the image is not a complex image.")
	.def("get_fft_phase", &EMAN::EMData::get_fft_phase, return_value_policy< manage_new_object >(), "return the phases of the FFT including the left half\n \nreturn The current FFT image's phase image.\nexception - ImageFormatException If the image is not a complex image.")
	.def("update", &EMAN::EMData::update, "Mark EMData as changed, statistics, etc will be updated at need.")
	.def("has_ctff", &EMAN::EMData::has_ctff, "check whether the image physical file has the CTF info or not.\n \nreturn True if it has the CTF information. Otherwise, false.")
	.def("calc_center_density", &EMAN::EMData::calc_center_density, "Calculates the density value at the peak of the\nimage histogram, sort of like the mode of the density.\n \nreturn The density value at the peak of the image histogram.")
	.def("calc_sigma_diff", &EMAN::EMData::calc_sigma_diff, "Calculates sigma above and below the mean and returns the\ndifference between them.\n \nreturn The difference between sigma above and below the mean.")
	.def("calc_min_location", &EMAN::EMData::calc_min_location, "Calculates the coordinates of the minimum-value pixel.\n \nreturn The coordinates of the minimum-value pixel.")
	.def("calc_max_location", &EMAN::EMData::calc_max_location, "Calculates the coordinates of the maximum-value pixel.\n \nreturn The coordinates of the maximum-value pixel.")
	.def("calc_max_location_wrap", &EMAN::EMData::calc_max_location_wrap, "Calculates the wrapped coordinates of the maximum value.\nThis function is useful in the context of Fourier correlation.\nyou can call this function to find the correct translational shift when using calc_ccf etc.\n \nreturn the wrapped coordinates of the maximum")
	.def("calc_center_of_mass", &EMAN::EMData::calc_center_of_mass, "Calculate the center of mass with a threshold (eg - 0 means use only positive values)")
	.def("calc_min_index", &EMAN::EMData::calc_min_index, "Calculates the index of minimum-value pixel when assuming\nall pixels are in a 1D array.\n \nreturn Index of the minimum-value pixel.")
	.def("calc_max_index", &EMAN::EMData::calc_max_index, "Calculates the index of maximum-value pixel when assuming\nall pixels are in a 1D array.\n \nreturn Index of the maximum-value pixel.")
	.def("calc_highest_locations", &EMAN::EMData::calc_highest_locations, "Calculate and return a sorted list of pixels whose values\nare above a specified threshold. The pixels are sorted\nfrom high to low.\n \nthreshold - The specified pixel value. Returned pixels should have higher values than it.\n \nreturn A sorted list of pixels with their values, and locations. Their values are higher than threshold.")
	.def("calc_n_highest_locations", &EMAN::EMData::calc_n_highest_locations, "Calculate and return a sorted list of N highest pixels. The pixels are sorted\nfrom high to low.\n \nn - number of pixels to return\n \nreturn A sorted list of N pixels with their values, and locations.")
	.def("find_pixels_with_value", &EMAN::EMData::find_pixels_with_value, "Finds pixels with exactly the specified value, and returns a list of them.")
	.def("get_edge_mean", &EMAN::EMData::get_edge_mean, "Calculates the mean pixel values around the (1 pixel) edge of the image.\n \nreturn The mean pixel values around the (1 pixel) edge.")
	.def("to_value", &EMAN::EMData::to_value, "sets all the pixels to a given value")
	.def("get_circle_mean", &EMAN::EMData::get_circle_mean, "Calculates the circular edge mean by applying a circular\nmask on 'this' image.\n \nreturn The circular edge mean.")
	.def("get_ctf", &EMAN::EMData::get_ctf, return_value_policy< manage_new_object >(), "Get ctf parameter of this image.\n \nreturn The ctf parameter.")
	.def("set_ctf", &EMAN::EMData::set_ctf, args("ctf"), "Set the CTF parameter of this image.\n \nctf - The CTF parameter object.")
	.def("get_translation", &EMAN::EMData::get_translation, "Get 'this' image's translation vector from the original location.\n \nreturn 'this' image's translation vector from the original location.")
	.def("set_translation", (void (EMAN::EMData::*)(const EMAN::Vec3f&) )&EMAN::EMData::set_translation, args("t"), "Set 'this' images' translation vector from the original location.\n \nt - The new translation vector.")
	.def("set_translation", (void (EMAN::EMData::*)(float, float, float) )&EMAN::EMData::set_translation, args("dx", "dy", "dz"), "Set 'this' images' translation vector from the original location.\n \ndx - The translation distance in x direction.\ndy - The translation distance in y direction.\ndz - The translation distance in z direction.")
	.def("get_transform", &EMAN::EMData::get_transform, "Get the 3D orientation of 'this' image.\n \nreturn The 3D orientation of 'this' image.")
	.def("set_rotation", (void (EMAN::EMData::*)(float, float, float) )&EMAN::EMData::set_rotation, args("az", "alt", "phi"), "Define the 3D orientation of this particle, also\nused to indicate relative rotations for reconstructions.\n \naz - 'az' Euler angle in EMAN convention.\nalt - 'alt' Euler angle in EMAN convention.\nphi - 'phi' Euler angle in EMAN convention.")
	.def("set_rotation", (void (EMAN::EMData::*)(const EMAN::Transform&) )&EMAN::EMData::set_rotation, args("t3d"), "Define the 3D orientation of this particle Orientation\ninformation is extracted from a Transform object and\nstored internally in EMAN (az,alt,phi) format.\n \nt3d - a Transform object containing the particle orientation.")
	.def("set_size", &EMAN::EMData::set_size, EMAN_EMData_set_size_overloads_1_4(args("nx", "ny", "nz","noalloc"), "Resize this EMData's main board memory pointer.\n \nnx - x size of this image.\nny - y size of this image. default to 1.\nnz - z size of this image. default to 1.\nnoalloc - if true, memory will not be allocated \n\nBadAllocException if memory allocation returns a null pointer."))
	.def("set_complex_size", &EMAN::EMData::set_complex_size, EMAN_EMData_set_complex_size_overloads_1_3(args("nx", "ny", "nz"), "Resize 'this' complex image.\n \nnx - x size of this image.\nny - y size of this image. deault to 1.\nnz - z size of this image. default to 1."))
//	.def("set_path", &EMAN::EMData::set_path, args("new_path"), "Set the path.\n \nnew_path - The new path.")
//	.def("set_pathnum", &EMAN::EMData::set_pathnum, args("n"), "Set the number of paths.\n \nn - The number of paths.")
	.def("get_2dview", (EMAN::MArray2D (EMAN::EMData::*)() const)&EMAN::EMData::get_2dview, "Get image raw pixel data in a 2D multi-array format.\nThe array shares the memory space with the image data.\nNotice: the subscription order is d[y][x] in Python, it's d[x][y] in C++\nIt should be used on 2D image only.\n \nreturn 2D multi-array format of the raw data.")
	.def("get_3dview", (EMAN::MArray3D (EMAN::EMData::*)() const)&EMAN::EMData::get_3dview, "Get image raw pixel data in a 3D multi-array format.\nThe array shares the memory space with the image data.\nNotice: the subscription order is d[z][y][x] in Python, it's d[x][y][z] in C++\nIt should be used on 3D image only.\n \nreturn 3D multi-array format of the raw data.")
	.def("get_2dcview", (EMAN::MCArray2D (EMAN::EMData::*)() const)&EMAN::EMData::get_2dcview, "Get complex image raw pixel data in a 2D multi-array format.\nThe array shares the memory space with the image data.\nIt should be used on 2D image only.\n \nreturn 2D multi-array format of the raw data.")
	.def("get_3dcview", (EMAN::MCArray3D (EMAN::EMData::*)() const)&EMAN::EMData::get_3dcview, "Get complex image raw pixel data in a 3D multi-array format.\nThe array shares the memory space with the image data.\nIt should be used on 3D image only.\n \nreturn 3D multi-array format of the raw data.")
	.def("get_3dcviewptr", &EMAN::EMData::get_3dcviewptr, return_value_policy< reference_existing_object >(), "Get pointer to a complex image raw pixel data in a 3D multi-array format.\nThe array shares the memory space with the image data.\nIt should be used on 3D image only.\n \nreturn Pointer to a 3D multi-array format of the raw data.")
	.def("get_2dview", (EMAN::MArray2D (EMAN::EMData::*)(int, int) const)&EMAN::EMData::get_2dview, args("x0, y0"), "Get image raw pixel data in a 2D multi-array format.\nThe data coordinates is translated by (x0,y0) such that\narray[y0][x0] points to the pixel at the origin location.\nthe data coordiates translated by (x0,y0). The\narray shares the memory space with the image data.\nIt should be used on 2D image only.\n \nx0 - X-axis translation amount.\ny0 - Y-axis translation amount.\n \nreturn 2D multi-array format of the raw data.")
	.def("get_3dview", (EMAN::MArray3D (EMAN::EMData::*)(int, int, int) const)&EMAN::EMData::get_3dview, args("x0, y0", "z0"), "Get image raw pixel data in a 3D multi-array format. The\ndata coordinates is translated by (x0,y0,z0) such that\narray[z0][y0][x0] points to the pixel at the origin location.\nthe data coordiates translated by (x0,y0,z0). The\narray shares the memory space with the image data.\nIt should be used on 3D image only.\n \nx0 - X-axis translation amount.\ny0 - Y-axis translation amount.\nz0 - Z-axis translation amount.\n \nreturn 3D multi-array format of the raw data.")
	.def("get_2dcview", (EMAN::MCArray2D (EMAN::EMData::*)(int, int) const)&EMAN::EMData::get_2dcview, args("x0, y0"), "Get complex image raw pixel data in a 2D multi-array format. The\ndata coordinates is translated by (x0,y0) such that\narray[y0][x0] points to the pixel at the origin location.\nthe data coordiates translated by (x0,y0). The\narray shares the memory space with the image data.\nIt should be used on 2D image only.\n \nx0 - X-axis translation amount.\ny0 - Y-axis translation amount.\n \nreturn 2D multi-array format of the raw data.")
	.def("get_3dcview", (EMAN::MCArray3D (EMAN::EMData::*)(int, int, int) const)&EMAN::EMData::get_3dcview, args("x0, y0", "z0"), "Get complex image raw pixel data in a 3D multi-array format. The\ndata coordinates is translated by (x0,y0,z0) such that\narray[z0][y0][x0] points to the pixel at the origin location.\nthe data coordiates translated by (x0,y0,z0). The\narray shares the memory space with the image data.\nIt should be used on 3D image only.\n \nx0 - X-axis translation amount.\ny0 - Y-axis translation amount.\nz0 - Z-axis translation amount.\n \nreturn 3D multi-array format of the raw data.")
	.def("get_attr", &EMAN::EMData::get_attr, args("attr_name"), "The generic way to get any image header information\ngiven a header attribute name. If the attribute does not exist,\nit will raise an exception.\n \nattr_name - The header attribute name.\n \nreturn The attribute value.\nexception - NotExistingObjectException when attribute not exist")
	.def("get_attr_default", &EMAN::EMData::get_attr_default, EMAN_EMData_get_attr_default_overloads_1_2(args("attr_name", "em_obj"), "The generic way to get any image header information\ngiven a header attribute name. If the attribute does not exist,\nit will return a default EMObject() object, which will be converted\nto None in Python. Or return any object user submit.\n \nattr_name - The header attribute name.\nem_obj - the default attribute to return when this attr_name not exist in attr_dict. default to None."))
	.def("set_attr", &EMAN::EMData::set_attr, args("key", "val"), "Set a header attribute's value from Python.\n \nkey - The header attribute name.\nval - The attribute value.")
	.def("get_attr_dict", &EMAN::EMData::get_attr_dict, "Get the image attribute dictionary containing all the\nimage attribute names and attribute values.\n \nreturn The image attribute dictionary containing all attribute names and values.")
	.def("set_attr_dict", &EMAN::EMData::set_attr_dict, args("new_dict"), "Merge the new values with the existing dictionary.\n \nnew_dict - The new attribute dictionary.")
	.def("del_attr", &EMAN::EMData::del_attr, args("attr_name"), "Delete the attribute from dictionary.\n \nattr_name - the attribute name to be removed.")
	.def("has_attr", &EMAN::EMData::has_attr, args("key"), "Ask if the header has a particular attribute.\n \nkey - the header attribute name.\n \nreturn whether or not the header has the name as a key/value entry")
	.def("del_attr_dict", &EMAN::EMData::del_attr_dict, args("del_keys"), "Delete the attributes from the dictionary.\n \ndel_keys - the attrutes' names to be removed.")
#ifdef DEBUG
	.def("debug_print_parms", &EMAN::EMData::debug_print_parms)
#endif	//DEBUG
	.def("get_xsize", &EMAN::EMData::get_xsize, "Get the image x-dimensional size.\n \nreturn Image x-dimensional size.")
	.def("get_ysize", &EMAN::EMData::get_ysize, "Get the image y-dimensional size.\n \nreturn Image y-dimensional size.")
	.def("get_zsize", &EMAN::EMData::get_zsize, "Get the image z-dimensional size.\n \nreturn Image z-dimensional size.")
	.def("get_size", &EMAN::EMData::get_size, "Get the number of allocated floats in the image (nx*ny*nz)\n \nreturn nx*ny*nz")
// 	.def("get_data_char", &EMAN::EMData::get_data_char, return_value_policy<reference_existing_object>(),"Get the binary pixel data as a raw uchar pointer")
	.def("get_data_as_vector", &EMAN::EMData::get_data_as_vector, "Get the pixel data as a vector\n \nreturn a vector containing the pixel data.")
	.def("get_data_string",&EMAN::EMData::get_data_pickle,"Returns a string representation of the floating point data in the image")
	.def("set_data_string",&EMAN::EMData::set_data_pickle, args("data_string"), "Sets the floating point data array from a string of binary data. Must be exactly the correct length.")
	.def("get_ndim", &EMAN::EMData::get_ndim, "Get image dimension.\n \nreturn image dimension.")
	.def("is_shuffled", &EMAN::EMData::is_shuffled, "Has this image been shuffled?\n \nreturn Whether this image has been shuffled to put origin in the center.")
	.def("is_FH", &EMAN::EMData::is_FH, "Is this a FH image?\n \nreturn Whether this is a FH image or not.")
	.def("is_complex", &EMAN::EMData::is_complex, "Is this a complex image?\n \nreturn Whether this is a complex image or not.")
	.def("is_real", &EMAN::EMData::is_real, "Is this a real image?\n \nreturn Whether this is image is real (not complex) or not.")
	.def("set_shuffled", &EMAN::EMData::set_shuffled, args("is_shuffled"), "Mark this image as a shuffled image.\n \nis_shuffled - If true, a shuffled image. If false, not a shuffled image.")
	.def("set_FH", &EMAN::EMData::set_FH, args("is_FH"), "Mark this complex image as a FH image.\n \nis_FH - If true, a FH image. If false, not a FH image.")
	.def("set_complex", &EMAN::EMData::set_complex, args("is_complex"), "Mark this image as a complex image.\n \nis_complex - If true, a complex image. If false, a real image.")
	.def("is_complex_x", &EMAN::EMData::is_complex_x, "Is this image a 1D FFT image in X direction?\n \nreturn Whether this image is a 1D FFT image in X direction.")
	.def("set_complex_x", &EMAN::EMData::set_complex_x, args("is_complex_x"), "Marks this image a 1D FFT image in X direction.\n \nis_complex_x - If true, a 1D FFT image in X direction; if false, not such an image.")
	.def("is_flipped", &EMAN::EMData::is_flipped, "Is this image flipped?\n \nreturn Whether this image is flipped or not.")
	.def("set_flipped", &EMAN::EMData::set_flipped, args("is_flipped"), "Mark this image as flipped.\n \nis_flipped - If true, mark this image as flipped;\nIf false, mark this image as not flipped.")
	.def("is_ri", &EMAN::EMData::is_ri, "Is this image a real/imaginary format complex image?\n \nreturn Whether this image is real/imaginary format complex image.")
	.def("set_ri", &EMAN::EMData::set_ri, args("is_ri"), "Mark this image as a real/imaginary format complex image.\n \nis_ri - If true, mark as real/imaginary format; If false, mark as amp/phase format.")
	.def("is_fftpadded", &EMAN::EMData::is_fftpadded, "Is this image already extended along x for ffts?\n \nreturn Whether this image is extended along x for ffts.")
	.def("set_fftpad", &EMAN::EMData::set_fftpad, args("is_fftpadded"), "Mark this image as already extended along x for ffts.\n \nis_fftpadded - If true, mark as padded along x; If\nfalse, mark as not padded along x.")
	.def("is_fftodd", &EMAN::EMData::is_fftodd, "Does this image correspond to a (real-space) odd nx?\n \nreturn Whether this image has a (real-space) odd nx.")
	.def("set_fftodd", &EMAN::EMData::set_fftodd, args("is_fftodd"), "Mark this image as having (real-space) odd nx.\n \nis_fftodd If true, mark as nx odd; If false, mark as nx not odd.")
	.def("set_nxc", &EMAN::EMData::set_nxc, args("nxc"), "Set the number of complex elements along x.\n \nnxc - is the number of complex elements along x.")
//	.def("get_path", &EMAN::EMData::get_path)
//	.def("get_pathnum", &EMAN::EMData::get_pathnum)
	.def("write_data",&EMAN::EMData::write_data,EMAN_EMData_write_data_overloads_2_6(args("fsp", "loc", "area", "file_nx", "file_ny", "file_nz"), "Dump the image pixel data in native byte order to a disk file.\n \nfsp - The filename to read the image data from\nloc - Location to seek to in the file before writing (size_t)\narea - The image region you want to read, default 0 means read the whole image(default=Null)\nfile_nx - Image x size.(default=0)\nfile_ny - Image y size.(default=0)\nfile_nz Image z size.(default=0)"))
	.def("read_data",&EMAN::EMData::read_data,EMAN_EMData_read_data_overloads_2_6(args("fsp", "loc", "area", "file_nx", "file_ny", "file_nz"), "Read the image pixel data in native byte order from a disk file.\nThe image should already have the correct dimensions.\n \nfsp - The filename to read the image data from\nloc - Location to seek to in the file before writing (size_t)\narea - The image region you want to read, default 0 means read the whole image(default=Null)\nfile_nx - Image x size.(default=0)\nfile_ny - Image y size.(default=0)\nfile_nz Image z size.(default=0)"))
	.def("process_inplace", &EMData_process_inplace_wrapper1,args("processorname"),return_value_policy< manage_new_object >(), "Apply a processor with its parameters on this image.\n \nprocessorname - Processor Name.\nparams - Processor parameters in a keyed dictionary. default to None.\n \nNotExistingObjectError If the processor doesn't exist.")
	.def("process_inplace", &EMData_process_inplace_wrapper2,args("processorname", "params"),return_value_policy< manage_new_object >(), "Apply a processor with its parameters on this image.\n \nprocessorname - Processor Name.\nparams - Processor parameters in a keyed dictionary. default to None.\n \nNotExistingObjectError If the processor doesn't exist.")
	.def("process", &EMData_process_wrapper1,args("processorname"),return_value_policy< manage_new_object >(), "Apply a processor with its parameters on a copy of this image, return result\nas a a new image. The returned image may or may not be the same size as this image.\n \nprocessorname - Processor Name.\nparams - Processor parameters in a keyed dictionary.\n \nreturn the processed result, a new image\n \nexception - NotExistingObjectError If the processor doesn't exist.")
	.def("process", &EMData_process_wrapper2,args("processorname", "params"),return_value_policy< manage_new_object >(), "Apply a processor with its parameters on a copy of this image, return result\nas a a new image. The returned image may or may not be the same size as this image.\n \nprocessorname - Processor Name.\nparams - Processor parameters in a keyed dictionary.\n \nreturn the processed result, a new image\n \nexception - NotExistingObjectError If the processor doesn't exist.")
	.def("process", (EMAN::EMData* (EMAN::EMData::*)(EMAN::Processor*) const )&EMAN::EMData::process, args("p"), "Call the process with an instance od Processor, usually this instance can\nbe get by (in Python) Processors.get('name', {'k':v, 'k':v})\n \np - the processor object", return_value_policy< manage_new_object >())
	.def("process_inplace", (void (EMAN::EMData::*)(EMAN::Processor*) )&EMAN::EMData::process_inplace, args("p"), "Call the process_inplace with an instance od Processor, usually this instancecan\nbe get by (in Python) Processors.get('name', {'k':v, 'k':v}).\n \np - the processor object")
	.def("cmp", &EMData_cmp_wrapper2, args("cmpname", "with"), "Compare this image with another image.\n \ncmpname - Comparison algorithm name.\nwith - The image you want to compare to.\nparams - Comparison parameters in a keyed dictionary, default to Null.\n \nreturn comparison score. The bigger, the better.\nexception - NotExistingObjectError If the comparison algorithm doesn't exist.")
	.def("cmp", &EMData_cmp_wrapper3, args("cmpname", "with", "params"), "Compare this image with another image.\n \ncmpname - Comparison algorithm name.\nwith - The image you want to compare to.\nparams - Comparison parameters in a keyed dictionary, default to Null.\n \nreturn comparison score. The bigger, the better.\nexception - NotExistingObjectError If the comparison algorithm doesn't exist.")
	.def("xform_align_nbest", &EMAN::EMData::xform_align_nbest, EMAN_EMData_xform_align_nbest_overloads_2_6(args("aligner_name", "to_img", "params", "nsoln", "cmp_name", "cmp_params"), "Align this image with another image, return the parameters of the \"n best\" solutions.\nThis function first added in the context of the 3D aligners used by e2tomohunter:\nwhich wants the n best solutions, as opposed to just the best. Return value is an\nordered vector of Dicts of length nsoln. The data with idx 0 has the best solution in it.\n \naligner_name - Alignment algorithm name.\nto_img - The image 'this' image aligns to.\nparams - Alignment algorithm parameters in a keyed dictionary, default to Null.\nnsoln - the number of solutions you want to receive in the return vector, default to 1.\ncmp_name - Comparison algorithm used in alignment, default to 'dot'.\ncmp_params - Parameter dictionary for comparison algorithm, default to NUll.\n \nreturn an ordered vector of Dicts of length nsoln. The Dicts in the vector have keys \"score\" (i.e. correlation score) and \"xform.align3d\" (Transform containing the alignment)\nexception - NotExistingObjectError If the alignment algorithm doesn't exist."))
	.def("align", &EMData_align_wrapper2,args("aligner_name", "to_img"),return_value_policy< manage_new_object >(), "Align this image with another image and return the result image.\n \naligner_name - Alignment algorithm name.\nto_img - The image 'this' image aligns to.\nparams - Alignment algorithm parameters in a keyed dictionary, default to Null.\ncmp_name - Comparison algorithm used in alignment, default to 'dot'.\ncmp_params - Parameter dictionary for comparison algorithm, default to Null.\n \nreturn The result image.\nexception - NotExistingObjectError If the alignment algorithm doesn't exist.")
	.def("align", &EMData_align_wrapper3,args("aligner_name", "to_img", "params"),return_value_policy< manage_new_object >(), "Align this image with another image and return the result image.\n \naligner_name - Alignment algorithm name.\nto_img - The image 'this' image aligns to.\nparams - Alignment algorithm parameters in a keyed dictionary, default to Null.\ncmp_name - Comparison algorithm used in alignment, default to 'dot'.\ncmp_params - Parameter dictionary for comparison algorithm, default to Null.\n \nreturn The result image.\nexception - NotExistingObjectError If the alignment algorithm doesn't exist.")
	.def("align", &EMData_align_wrapper4,args("aligner_name", "to_img", "params", "cmp_name"),return_value_policy< manage_new_object >(), "Align this image with another image and return the result image.\n \naligner_name - Alignment algorithm name.\nto_img - The image 'this' image aligns to.\nparams - Alignment algorithm parameters in a keyed dictionary, default to Null.\ncmp_name - Comparison algorithm used in alignment, default to 'dot'.\ncmp_params - Parameter dictionary for comparison algorithm, default to Null.\n \nreturn The result image.\nexception - NotExistingObjectError If the alignment algorithm doesn't exist.")
	.def("align", &EMData_align_wrapper5,args("aligner_name", "to_img", "params", "cmp_name", "cmp_params"),return_value_policy< manage_new_object >(), "Align this image with another image and return the result image.\n \naligner_name - Alignment algorithm name.\nto_img - The image 'this' image aligns to.\nparams - Alignment algorithm parameters in a keyed dictionary, default to Null.\ncmp_name - Comparison algorithm used in alignment, default to 'dot'.\ncmp_params - Parameter dictionary for comparison algorithm, default to Null.\n \nreturn The result image.\nexception - NotExistingObjectError If the alignment algorithm doesn't exist.")
	.def("project", (EMAN::EMData* (EMAN::EMData::*)(const std::string&, const EMAN::Dict&) )&EMAN::EMData::project, EMAN_EMData_project_overloads_1_2(args("projector_name", "params"), "Calculate the projection of this image and return the result.\n \nprojector_name - Projection algorithm name.\nparams - Projection Algorithm parameters, default to Null.\n \nreturn The result image.\nexception - NotExistingObjectError If the projection algorithm doesn't exist.")[ return_value_policy< manage_new_object >() ])
	.def("project", (EMAN::EMData* (EMAN::EMData::*)(const std::string&, const EMAN::Transform&) )&EMAN::EMData::project, args("projector_name", "t3d"), "Calculate the projection of this image and return the result.\n \nprojector_name - Projection algorithm name.\nt3d - Transform object used to do projection.\n \nreturn The result image.\nexception - NotExistingObjectError If the projection algorithm doesn't exist.", return_value_policy< manage_new_object >() )
	.def("backproject", &EMAN::EMData::backproject, EMAN_EMData_backproject_overloads_1_2(args("peojector_name", "params"), "Calculate the backprojection of this image (stack) and return the result.\n \nprojector_name - Projection algorithm name. Only \"pawel\" and \"chao\" have been implemented now.\nparams - Projection Algorithm parameters, default to Null.\n \nreturn The result image.\nexception - NotExistingObjectError If the projection algorithm doesn't exist.")[ return_value_policy< manage_new_object >() ])
	.def("do_fft", &EMAN::EMData::do_fft, return_value_policy< manage_new_object >(), "return the fast fourier transform (FFT) image of the current\nimage. the current image is not changed. The result is in\nreal/imaginary format.\n \nreturn The FFT of the current image in real/imaginary format.")
	.def("do_fft_inplace", &EMAN::EMData::do_fft_inplace, return_value_policy< reference_existing_object >(), "Do FFT inplace. And return the FFT image.\n \nreturn The FFT of the current image in real/imaginary format.")
	.def("do_ift", &EMAN::EMData::do_ift, return_value_policy< manage_new_object >(), "return the inverse fourier transform (IFT) image of the current\nimage. the current image may be changed if it is in amplitude/phase\nformat as opposed to real/imaginary format - if this change is\nperformed it is not undone.\n \nreturn The current image's inverse fourier transform image.\nexception - ImageFormatException If the image is not a complex image.")
	.def("do_ift_inplace", &EMAN::EMData::do_ift_inplace, return_value_policy< reference_existing_object >(), "Do IFT inplace. And return the IFT image.\n \nreturn The IFT image.")
	.def("bispecRotTransInvN", &EMAN::EMData::bispecRotTransInvN, return_value_policy< reference_existing_object >(), args("N", "NK"), "This computes the rotational and translational bispectral\ninvariants of an image. The invariants are labelled by the Fourier\nHarmonic label given by N.\nNK is the number of Fourier components one wishes to use in calculating this bispectrum.\nthe output is a single 2D image whose x,y labels are lengths, corresponding to the two lengths of sides of a triangle.")
	.def("bispecRotTransInvDirect", &EMAN::EMData::bispecRotTransInvDirect, return_value_policy< reference_existing_object >(), args("type"), "This computes the rotational and translational bispectral\ninvariants of an image.\nthe output is a single 3d Volume whose x,y labels are lengths,\ncorresponding to the two lengths of sides of a triangle.\nthe z label is for the angle.")
#ifdef EMAN2_USING_CUDA
	.def("do_fft_cuda", &EMAN::EMData::do_fft_cuda, return_value_policy< manage_new_object >(), "return the fast fourier transform (FFT) image of the current\nimage inplace. The result is in\nreal/imaginary format and exists only on the GPU.\n \nreturn The FFT of the current image in real/imaginary format, existing on the GPU.")
	.def("do_fft_inplace_cuda", &EMAN::EMData::do_fft_inplace_cuda, return_value_policy< manage_new_object >(), "return the fast fourier transform (FFT) image of the current\nimage. the current image is not changed. The result is in\nreal/imaginary format and exists only on the GPU.\n \nreturn The FFT of the current image in real/imaginary format, existing on the GPU.")
	.def("do_ift_cuda", &EMAN::EMData::do_ift_cuda, return_value_policy< manage_new_object >(), "return the inverse fourier transform (IFT) image of the current\nimage. The result exists only on the GPU.\n \npreserve_input - whether or not this EMData object should be preserved. If this is unecessary than we can avoid a copy and run faster(default=True).\n \nreturn The FFT of the current image in real/imaginary format, existing on the GPU.")
	.def("cuda_cleanup", &EMAN::EMData::cuda_cleanup, "clean up CUDA")
	.def("cuda_initialize", &EMAN::EMData::cuda_initialize, "initialize cuda")
	.staticmethod("cuda_cleanup")
	.staticmethod("cuda_initialize")
	.def("copy_to_cudaro", &EMAN::EMData::copy_to_cudaro, "Copy RO data to the CUDA device\n Can be useful to force CUDA operations")
	.def("copy_to_cuda", &EMAN::EMData::copy_to_cuda, "Copy RW data to the CUDA device\n Can be useful to force CUDA operations")
	.def("copy_rw_to_ro", &EMAN::EMData::copy_rw_to_ro, "Copy RW data to the CUDA device RO \n Can be useful to force CUDA operations")
	.def("switchoncuda", &EMAN::EMData::switchoncuda, "Turn on CUDA")
	.def("switchoffcuda", &EMAN::EMData::switchoffcuda, "Turn off CUDA")
	.def("getcudalock", &EMAN::EMData::getcudalock, "Return where the /tmp/cuda locks are stored")
	.staticmethod("switchoncuda")
	.staticmethod("switchoffcuda")
	.staticmethod("getcudalock")
#endif // EMAN2_USING_CUDA
	.def("render_ap24", &EMAN::EMData::render_ap24, args("x", "y", "xsize", "ysize", "bpl", "scale", "min_gray", "max_gray", "min_render", "max_render", "gamma", "flags"), "Render the image into an 8-bit image. 2D images only.\nflags provide a way to do unusual things with this function, such\nas calculating a histogram of the rendered area.\n \nx - origin of the area to render\ny - origin of the area to render\nxsize - size of the area to render in output pixels\nysize - size of the area to render in output pixels\nbpl - bytes per line, if asrgb remember *3\nscale - scale factor for rendering\nmin_gray - minimum gray value to render (0-255)\nmax_gray - maximum gray value to render (0-255)\nmin_render - float image density corresponding to min_gray\nmax_render - float image density corresponding to max_gray\ngamma - \nflags	- 1-duplicate each output pixel 3x for RGB rendering,2-add a 256 int greyscale histogram to the end of the image array,4-invert y axis,8-render 32 bit 0xffRRGGBB\n \nexception - ImageDimensionException If the image is not 2D.")
	.def("ri2ap", &EMAN::EMData::ri2ap, "convert the complex image from real/imaginary to amplitude/phase")
	.def("ap2ri", &EMAN::EMData::ap2ri, "convert the complex image from amplitude/phase to real/imaginary")
	.def("ri2inten", &EMAN::EMData::ri2inten, "convert the complex image from real/imaginary to Intensity/0.\nThis conversion cannot be reversed, and the image remains marked as R/I")
	.def("insert_clip", &EMAN::EMData::insert_clip, args("block", "orogin"), "Insert a clip into this image.\nVery robust clip insertion code works in all way you might think possible.\n \nblock - An image block.\norigin - The origin location to insert the clip.")
	.def("insert_scaled_sum", &EMAN::EMData::insert_scaled_sum, EMAN_EMData_insert_scaled_sum_overloads_2_4(args("block", "center", "scale", "mult_factor"), "Add a scaled image into another image at a specified location.\nThis is used, for example, to accumulate gaussians in\nprograms like pdb2mrc.py. The center of 'block' will be positioned at\n'center' with scale factor 'scale'. Densities will be interpolated in\n'block' and multiplied by 'mult'.\n \nblock - The image to inserted.\ncenter - The center of the inserted block in 'this'.\nscale - Scale factor, default to 1.0.\nmult_factor - Number used to multiply the block's densities, default to 1.0.\n \nexception - ImageDimensionException If 'this' image is not 2D/3D."))
	.def("copy", &EMAN::EMData::copy, return_value_policy< manage_new_object >(), "Make a copy of this image including both data and header.\n \nreturn A copy of this image including both data and header.")
	.def("copy_head", &EMAN::EMData::copy_head, return_value_policy< manage_new_object >(), "Make an image with a copy of the current image's header.\n \nreturn An image with a copy of the current image's header.")
	.def("add", (void (EMAN::EMData::*)(float, int) )&EMAN::EMData::add, EMAN_EMData_add_overloads_1_2(args("f", "keepzero"), "add a number to each pixel value of the image. Image may be real or complex.\n \nf - The number added to 'this' image.\nkeepzero - If set will not modify pixels that are exactly zero, default to 0."))
	.def("add", (void (EMAN::EMData::*)(const EMAN::EMData&) )&EMAN::EMData::add, args("image"), "add a same-size image to this image pixel by pixel.\n \nimage - The image added to 'this' image.\n \nexception - ImageFormatException If the 2 images are not same size.")
	.def("addsquare", (void (EMAN::EMData::*)(const EMAN::EMData&) )&EMAN::EMData::addsquare, args("image"), "add the squared value of each pixel from a same-size image to this image.\n \nimage - The image whose square is added to 'this' image.\n \nexception ImageFormatException If the 2 images are not same size.")
	.def("sub", (void (EMAN::EMData::*)(float) )&EMAN::EMData::sub, args("f"), "subtract a float number to each pixel value of the image.\n \nf - The float number subtracted from 'this' image.")
	.def("sub", (void (EMAN::EMData::*)(const EMAN::EMData&) )&EMAN::EMData::sub, args("image"), "subtract a same-size image from this image pixel by pixel.\n \nimage - The image subtracted  from 'this' image.\n \nexception - ImageFormatException If the 2 images are not same size.")
	.def("subsquare", (void (EMAN::EMData::*)(const EMAN::EMData&) )&EMAN::EMData::subsquare, args("image"), "subtract the squared value of each pixel from a same-size image to this image.\n \nimage - The image whose square is subtracted from 'this' image.\n \nexception - ImageFormatException If the 2 images are not same size.")
	.def("mult", (void (EMAN::EMData::*)(const EMAN::EMData&, bool) )&EMAN::EMData::mult, EMAN_EMData_mult_overloads_1_2(args("image", "prevent_complex_multiplication"), "multiply each pixel of this image with each pixel of some other same-size image.\n \nimage - The image multiplied to 'this' image.\nprevent_complex_multiplication - If the image is complex, this flag will override complex multiplication and just multiply each pixel by the other.(default=False)\n \nexception - ImageFormatException If the 2 images are not same size."))
	.def("mult", (void (EMAN::EMData::*)(int) )&EMAN::EMData::mult, args("n"), "multiply an integer number to each pixel value of the image.\n \nn - The integer multiplied to 'this' image.")
	.def("mult", (void (EMAN::EMData::*)(float) )&EMAN::EMData::mult, args("f"), "multiply a float number to each pixel value of the image.\n \nf - The float multiplied to 'this' image.")
	.def("div", (void (EMAN::EMData::*)(float) )&EMAN::EMData::div, args("f"), "make each pixel value divided by a float number.\n \nf - The float number 'this' image divided by.")
	.def("div", (void (EMAN::EMData::*)(const EMAN::EMData&) )&EMAN::EMData::div, args("image"), "make each pixel value divided by pixel value of another same size image.\n \nimage - The image 'this' image divided by.\n \nexception - ImageFormatException If the 2 images are not same size.")
	.def("to_zero", &EMAN::EMData::to_zero, "Set all the pixel value = 0.")
	.def("to_one", &EMAN::EMData::to_one, "Set all the pixel value = 1.")
	.def("dot", &EMAN::EMData::dot, args("with"), "Dot product 2 images. The 2 images must be of same size.\nShortcut for cmp('dot').\n \nwith - The image to do dot product with.\n \nreturn The dot product result.\nexception - NullPointerException if with is a NULL image.")
	.def("get_row", &EMAN::EMData::get_row, return_value_policy< manage_new_object >(), args("row_index"), "Get one row of a 1D/2D image.\n \nrow_index - Index of the row.\n \nreturn A 1D image with the row data.\nexception - ImageDimensionException If this image is 3D.")
	.def("set_row", &EMAN::EMData::set_row, args("data", "row_index"), "Set one row of a 1D/2D image.\n \ndata - The row image data.\nrow_index - Index of the row.\n \nexception - ImageDimensionException If this image is 3D.")
	.def("get_col", &EMAN::EMData::get_col, return_value_policy< manage_new_object >(), args("col_index"), "Get one column of a 2D images.\n \ncol_index - Index of the column.\n \nreturn A 1D image with the column data.\nexception - ImageDimensionException If this image is not 2D.")
	.def("set_col", &EMAN::EMData::set_col, args("data", "col_index"), "Set one column of a 2D image.\n \ndata - The column image data.\ncol_index - Index of the column.\n \nexception - ImageDimensionException If this image is not 2D.")
	.def("get_value_at", (float (EMAN::EMData::*)(int, int, int) const)&EMAN::EMData::get_value_at, args("x", "y", "z"), "Get the pixel density value at coordinates (x,y,z).\nThe validity of x, y, and z is not checked.\n \nx - The x cooridinate.\ny - The y cooridinate.\nz - The z cooridinate.\n \nreturn The pixel density value at coordinates (x,y,z).")
	.def("get_value_at", (float (EMAN::EMData::*)(int, int) const)&EMAN::EMData::get_value_at, args("x", "y"), "Get the pixel density value at coordinates (x,y). 2D only.\nThe validity of x, y is not checked.\n \nx - The x cooridinate.\ny The y cooridinate.\n\nreturn The pixel density value at coordinates (x,y).")
	.def("get_value_at", (float (EMAN::EMData::*)(size_t) const)&EMAN::EMData::get_value_at, args("i"), "Get the pixel density value given an index 'i' assuming\nthe pixles are stored in a 1D array. The validity of i\nis not checked.\n \ni - 1D data array index.\n \nreturn The pixel density value")
	.def("get_value_at_wrap", (float (EMAN::EMData::*)(int, int, int) const)&EMAN::EMData::get_value_at_wrap, args("x", "y", "z"), "Get the pixel density value at coordinates (x,y,z).\nShould only be called on 3D images - no errors are thrown\nWraps pixel values if they are negative - i.e. is circulant\nFor example, if x = -1, then the pixel at nx-1  is returned\n \nx - The x cooridinate.\ny - The y cooridinate.\nz - The z cooridinate.\n \nreturn The pixel density value at circulant coordinates (x,y,z).")
	.def("get_value_at_wrap", (float (EMAN::EMData::*)(int, int) const)&EMAN::EMData::get_value_at_wrap, args("x", "y"), "Get the pixel density value at coordinates (x,y).\nShould only be called on 2D images - no errors are thrown\nWraps pixel values if they are negative - i.e. is circulant\nFor example, if x = -1, then the pixel at nx-1  is returned\n \nx - The x cooridinate.\ny - The y cooridinate.\n \nreturn The pixel density value at circulant coordinates (x,y).")
	.def("get_value_at_wrap", (float (EMAN::EMData::*)(int) const)&EMAN::EMData::get_value_at_wrap, args("x"), "Get the pixel density value at coordinates (x).\nShould only be called on 1D images - no errors are thrown\nWraps pixel values if they are negative - i.e. is circulant\nFor example, if x = -1, then the pixel at nx-1  is returned.\n \nx - The x cooridinate.\n \nreturn The pixel density value at circulant coordinates (x).")
	.def("get_complex_at", (std::complex<float> (EMAN::EMData::*)(const int &, const int &) const)&EMAN::EMData::get_complex_at)
	.def("get_complex_at", (std::complex<float> (EMAN::EMData::*)(const int &, const int &, const int &) const)&EMAN::EMData::get_complex_at)
	.def("get_complex_index", (size_t (EMAN::EMData::*)(const int &, const int &, const int &) const)&EMAN::EMData::get_complex_index)
	.def("set_complex_at", (void (EMAN::EMData::*)(const int &, const int &, const std::complex<float> &) )&EMAN::EMData::set_complex_at)
	.def("set_complex_at", (void (EMAN::EMData::*)(const int &, const int &, const int &, const std::complex<float> &) )&EMAN::EMData::set_complex_at)
	.def("add_complex_at", (size_t (EMAN::EMData::*)(const int &, const int &, const int &,const std::complex<float> &) )&EMAN::EMData::add_complex_at)
	.def("sget_value_at", (float (EMAN::EMData::*)(int, int, int) const)&EMAN::EMData::sget_value_at, args("x", "y", "z"), "A safer, slower way to get the pixel density value at\ncoordinates (x,y,z). The validity of x, y, and z is checked.\nIf the coordinates are out of range, return 0.\n \nx - The x cooridinate.\ny The y cooridinate.\nz The z cooridinate.\n \nreturn The pixel density value at coordinates (x,y,z).")
	.def("sget_value_at", (float (EMAN::EMData::*)(int, int) const)&EMAN::EMData::sget_value_at, args("x", "y"), "A safer, slower way to get the pixel density value at\ncoordinates (x,y). 2D only. The validity of x, y is checked.\nIf the coordinates are out of range, return 0.\n \nx - The x cooridinate.\ny - The y cooridinate.\n \nreturn The pixel density value at coordinates (x,y).")
	.def("sget_value_at", (float (EMAN::EMData::*)(size_t) const)&EMAN::EMData::sget_value_at, args("i"), "A safer, slower way to get the pixel density value\ngiven an index 'i' assuming\nthe pixles are stored in a 1D array. The validity of i\nis checked. If i is out of range, return 0.\n \ni - 1D data array index.\n \nreturn The pixel density value")
	.def("sget_value_at_interp", (float (EMAN::EMData::*)(float, float) const)&EMAN::EMData::sget_value_at_interp,args("x", "y"), "Get pixel density value at interpolation of (x,y).\nThe validity of x, y is checked.2D image only.\n \nx - The x cooridinate.\ny - The y cooridinate.\n \nreturn The pixel density value at coordinates (x,y).")
	.def("sget_value_at_interp", (float (EMAN::EMData::*)(float, float, float) const)&EMAN::EMData::sget_value_at_interp, args("x", "y", "z"), "Get the pixel density value at interpolation of (x,y,z).\nThe validity of x, y, and z is checked.\n \nx - The x cooridinate.\ny - The y cooridinate.\nz - The z cooridinate.\n \nreturn The pixel density value at coordinates (x,y,z).")
	.def("set_value_at", (void (EMAN::EMData::*)(int, int, int, float) )&EMAN::EMData::set_value_at, args("x", "y", "z", "v"), "Set the pixel density value at coordinates (x,y,z).\nThis implementation does bounds checking.\n \nx - The x cooridinate.\ny - The y cooridinate.\nz - The z cooridinate.\nv - The pixel density value at coordinates (x,y,z).\n \nexception - OutofRangeException wehn index out of image data's range.")
	.def("set_value_at_fast", (void (EMAN::EMData::*)(int, int, int, float) )&EMAN::EMData::set_value_at_fast, args("x", "y", "z", "v"), "Set the pixel density value at coordinates (x,y,z).\nThe validity of x, y, and z is not checked.\nThis implementation has no bounds checking.\n \nx - The x cooridinate.\ny - The y cooridinate.\nz - The z cooridinate.\nv - The pixel density value at coordinates (x,y,z).")
	.def("set_value_at", (void (EMAN::EMData::*)(int, int, float) )&EMAN::EMData::set_value_at, args("x", "y", "v"), "Set the pixel density value at coordinates (x,y).\n2D image only.\n \nx - The x cooridinate.\ny - The y cooridinate.\nv - The pixel density value at coordinates (x,y).\n \nexception - OutofRangeException wehn index out of image data's range.")
	.def("set_value_at_fast", (void (EMAN::EMData::*)(int, int, float) )&EMAN::EMData::set_value_at_fast, args("x", "y", "v"), "Set the pixel density value at coordinates (x,y).\n2D image only. The validity of x, y, is not checked.\n \nx - The x cooridinate.\ny - The y cooridinate.\nv - The pixel density value at coordinates (x,y).")
	.def("set_value_at", (void (EMAN::EMData::*)(int, float) )&EMAN::EMData::set_value_at, args("x", "v"), "Set the pixel density value at coordinate (x).\n1D image only.\n \nx - The x cooridinate.\nv The pixel density value at coordinate (x).\n \nexception - OutofRangeException wehn index out of image data's range.")
//#ifndef	_WIN32
//	.def("set_array_offsets", (void (EMAN::EMData::*)(const int, const int, const int) )&EMAN::EMData::set_array_offsets, EMAN_EMData_set_array_offsets_overloads_0_3(args("xoff", "yoff", "zoff"), "Set the array offset by three integers."))
//	.def("set_array_offsets", (void (EMAN::EMData::*)(std::vector<int,std::allocator<int> >) )&EMAN::EMData::set_array_offsets, args("offsets"), "Set the array offset by a vector of three integers.")
//#endif
	//.def("get_array_offsets", &EMAN::EMData::get_array_offsets)
	.def("power", &EMAN::EMData::power, return_value_policy< manage_new_object >(), args("n"), "return a image to the power of n.\n \nn	- the power of this image\n \nreturn a image which is the nth power of this image\nexception - InvalidValueException n must be >= 0")
	.def("sqrt", &EMAN::EMData::sqrt, return_value_policy< manage_new_object >(), "return square root of current image\n \nreturn a image which is the square root of this image\nexception - ImageFormatException real image only")
	.def("log", &EMAN::EMData::log, return_value_policy< manage_new_object >(), "return natural logarithm image for a image\n \nreturn a image which is the natural logarithm of this image\nexception - InvalidValueException pixel value must be >= 0\n")
	.def("log10", &EMAN::EMData::log10, return_value_policy< manage_new_object >(), "return base 10 logarithm image for a image\n \nreturn a image which is the base 10 logarithm of this image\nexception - InvalidValueException pixel value must be >= 0\nexception - ImageFormatException real image only")
	.def("real", &EMAN::EMData::real, return_value_policy< manage_new_object >(), "return real part of a complex image as a real image format,\nif this image is a real image, return a copy of this image.\n \nreturn a real image which is the real part of this image.")
	.def("imag", &EMAN::EMData::imag, return_value_policy< manage_new_object >(), "return imaginary part of a complex image as a real image format.\n \nreturn a real image which is the imaginary part of this image.\nexception - InvalidCallException if this image is a real image\nexception - InvalidCallException if this image is a complex image in amplitude/phase format")
	.def("absi", &EMAN::EMData::absi, return_value_policy< manage_new_object >(), "For a real image, it returns a same size image with abs() of each pixel.\nFor a complex image, it returns a image in size (nx/2,ny,nz),\nthe pixel value output[i]=sqrt(input[i]*input[i]+input[i+1]*input[i+1])\n \nexception InvalidCallException this function call require a complex image in real/imaginary format.")
	.def("amplitude", &EMAN::EMData::amplitude, return_value_policy< manage_new_object >(), "return amplitude part of a complex image as a real image format\n \nreturn EMData * a real image which is the amplitude part of this image\nexception - InvalidCallException if this image is a real image or is in real/imaginary format")
	.def("phase", &EMAN::EMData::phase, return_value_policy< manage_new_object >(), "return phase part of a complex image as a real image format\n \nreturn EMData * a real image which is the phase part of this image\nexception - InvalidCallException if this image is a real image or is in real/imaginary format")
	.def("real2complex", &EMAN::EMData::real2complex, EMAN_EMData_real2complex_overloads_0_1(args("img"), "create a complex image from a real image, this complex image is in real/imaginary format\n \nimg - give an artificial imaginary part, default to 0.0.\nreturn a complex image which is generated from a real image\nexception  - InvalidCallException this function can not be called by complex image")[ return_value_policy< manage_new_object >() ])
	.def("real2FH", &EMAN::EMData::real2FH, return_value_policy< manage_new_object >(), args("OverSamplekB"), "returns the fourier harmonic transform (FH) image of the current\nimage (in real space). The current image is not changed. The result is in\nreal/imaginary format. The FH switch is set on.\n \nOverSamplekB - is a parameter controlling the fineness of the Fourier sampling\n \nreturn the Fourier Harmonic image\n")
	.def("FH2F", &EMAN::EMData::FH2F, EMAN_EMData_FH2F_overloads_2_3(args("Size", "OverSamplekB", "IntensityFlag"), "returns the fourier version of the image from the FH version.\nThe current image is not changed. The result is in real/imaginary format.\nThe FH switch is set off.\n \nSize - is the size of the image to be returned\nOverSamplekB - is a parameter controlling the fineness of the Fourier sampling\nIntensityFlag - =0 is the usual; =1 means that the input was an intensity, default to 0.\n \nreturn the shuffled version of the FFT")[ return_value_policy< manage_new_object >() ])
	.def("FH2Real", &EMAN::EMData::FH2Real, EMAN_EMData_FH2Real_overloads_2_3(args("Size", "OverSamplekB", "IntensityFlag"), "returns the real version of the image from the FH version.\nThe current image is not changed.\nThe result is in real format.\n \nSize - is the size of the image to be returned\nOverSamplekB - is a parameter controlling the fineness of the Fourier sampling\nIntensityFlag - =0 is the usual; =1 means that the input was an intensity, default to 0.\n \nreturn the real version of the data")[ return_value_policy< manage_new_object >() ])
	.def("rotavg", &EMAN::EMData::rotavg, return_value_policy< manage_new_object >(), "Create a (1-D) rotationally averaged image.\n \nreturn 1-D rotationally-averaged image\n \nexception - ImageDimensionException If 'this' image is 1D.")
	.def("rotavg_i", &EMAN::EMData::rotavg_i, return_value_policy< manage_new_object >(), "Create a 2-D or 3-D rotationally averaged image.\n \nreturn 2-D or 3-D rotationally-averaged image\nexception - ImageDimensionException If 'this' image is 1D.")
	.def("mult_radial", &EMAN::EMData::mult_radial, return_value_policy< manage_new_object >(), args("radial"), "Multiply radially a 2-D or 3-D image by a 1-D image.\n \nradial - the 1-D image multiply to\n \nreturn 2-D or 3-D radially multiplied image\nexception - ImageDimensionException If 'this' image is 1D.")
	.def("cog", &EMAN::EMData::cog, "Calculates the Center of Gravity\nand the Radius of Gyration of the image.\n \nreturns the mass and the radius as vectors.")
	.def("calc_fourier_shell_correlation", &EMAN::EMData::calc_fourier_shell_correlation, EMAN_EMData_calc_fourier_shell_correlation_overloads_1_2(args("with", "width"), "Calculate CCF in Fourier space as a function of spatial frequency\nbetween a pair of 2-3D images (corners not included).\nThe input image 'with' must have the same size to 'this' image.\nInput images can be either real or Fourier in arbitrary combination.\n \nwith -The image used to caculate the fourier shell.\nwidth - Ring/shell width in Fourier space, default to 1.0.\n \nexception - ImageFormatException If the 2 images are not same size.\nexception - NullPointerException if the input image is null\nexception - Cannot calculate FSC for 1D images\nreturn  Vector of 3*k FSC results (frequencies, FSC values, error)\nk - length of FSC curve, depends on dimensions of the image and ring width\n1 column - normalized frequency [0,0.5]\n2 column - FSC,\n3 column - error of the FSC = 1/sqrt(n), where n is the number of Fourier coefficients within given shell."))
	.def("scale_factors", &EMAN::EMData::scale_factors, "Calculate scale_factors in Fourier space as a function of spatial frequency\nbetween a pair of 2-3D images (corners not included).\nThe input image 'with' must have the same size to 'this' image.\nInput images can be either real or Fourier in arbitrary combination.\n \nwith -The image used to caculate the fourier shell.\beg - beginning Ring/shell in Fourier space.\end - ending Ring/shell in Fourier space.\n \nexception - ImageFormatException If the 2 images are not same size.\nexception - NullPointerException if the input image is null\nexception - Cannot calculate FSC for 1D images\nreturn  Vector of 3*k FSC results (frequencies, FSC values, error)\nk - length of FSC curve, depends on dimensions of the image and ring width\n1 column - normalized frequency [0,0.5]\n2 column - FSC,\n3 column - error of the FSC = 1/sqrt(n), where n is the number of Fourier coefficients within given shell.")
	.def("average_circ_sub", &EMAN::EMData::average_circ_sub, return_value_policy< manage_new_object >(), "Subtract average outside of a circle\n \nreturn image with sbtracted average outside of a circle.")
//	.def("onelinenn", &EMAN::EMData::onelinenn)
//	.def("onelinenn_mult", &EMAN::EMData::onelinenn_mult)
	.def("nn", &EMAN::EMData::nn, EMAN_EMData_nn_overloads_3_4(args("wptr", "myfft", "tf", "mult"), "Nearest Neighbor interpolation.\nModifies the current object.\n \nwptr - Normalization data.\nmyfft - FFT data.\ntf - Transform reference\nmult - default to 1."))
	.def("nn_SSNR", &EMAN::EMData::nn_SSNR, EMAN_EMData_nn_SSNR_overloads_4_5(args("wptr", "wptr2", "myfft", "tf", "mult"), "Nearest Neighbor interpolation, meanwhile return necessary data such as\nKn, sum_k(F_k^n) ans sum_k(|F_k^n|^2)\nModifies the current object.\n \nwptr - Normalization data.\nwptr2 - \nmyfft - FFT data.\ntf - Transform reference\nmult - default to 1."))
	.def("nn_SSNR_ctf", &EMAN::EMData::nn_SSNR_ctf, EMAN_EMData_nn_SSNR_ctf_overloads_5_6(args("wptr", "wptr2", "wptr3", "myfft", "tf", "mult"), "Nearest Neighbor interpolation, meanwhile return necessary data such as\nKn, sum_k(F_k^n) ans sum_k(|F_k^n|^2)\nModifies the current object.\n \nwptr - Normalization data.\nwptr2 - \nwptr3 - \nmyfft - FFT data.\ntf - Transform reference\nmult - default to 1."))
	.def("symplane0", &EMAN::EMData::symplane0, args("norm"), "Calculate Wiener summation from the inserted 2D slice\nput the summation into 3D grids using nearest neighbour approximation\na. Map the 2D coordinates of the interted slice into 3D grid using 3D transformation\nb. calculate 2D CTF_K^2  and CTF_K*F_K, and put them on the voxel of 3D volume\nc. count the number of images entering each boxel wptr3")
	.def("symplane1", &EMAN::EMData::symplane1, args("norm", "norm2"), "Symmetrize plane 0\nModifies the current object.\n \nnorm - Normalization data.\nnorm2 -")
	.def("symplane2", &EMAN::EMData::symplane2, args("norm", "norm2", "norm3"), "Symmetrize plane 0\nModifies the current object.\n \nnorm - Normalization data.\nnorm2 - \nnorm3 -")
//	.def("onelinenn_ctf", &EMAN::EMData::onelinenn_ctf)
	.def("nn_ctf", &EMAN::EMData::nn_ctf, args("w", "myfft", "tf", "mult"), "Nearest Neighbor interpolation.\nModifies the current object.\n \nw - Normalization data.\nmyfft - FFT data.\ntf - Transform reference\nmult - ")
//	.def("onelinenn_ctf_applied", &EMAN::EMData::onelinenn_ctf_applied)
//	.def("nn_ctf_applied", &EMAN::EMData::nn_ctf_applied)
//	.def("symplane0_ctf", &EMAN::EMData::symplane0_ctf)
	.def("symvol", &EMAN::EMData::symvol, return_value_policy< manage_new_object >(), args("symmetry"), "Symmetrize volume in real space.\n \nsymmetry - Point group of the target volume.\n \nreturn New symmetrized volume object.")
	.def("rot_scale_trans2D", &EMAN::EMData::rot_scale_trans2D, EMAN_EMData_rot_scale_trans2D_overloads_1_4(args("ang", "delx", "dely", "scale"), "Rotate-Shift-Scale-Circulantly image.\nIf the image is a volume, then all slices are rotated/translated/scaled.\n \nang -Rotation angle in degrees.\ndelx - Amount to shift rotation origin along x\ndely - Amount to shift rotation origin along y\nscale - Scaling factor (default=1.0)\n \nreturn New rotated/shifted/scaled image\nexception ImageDimensionException can not rotate 1 D image\nexception ImageDimensionException can not rotate 3 D image")[ return_value_policy< manage_new_object >() ])
	.def("rot_scale_trans2D_background", &EMAN::EMData::rot_scale_trans2D_background, EMAN_EMData_rot_scale_trans2D_background_overloads_1_4(args("ang", "delx", "dely", "scale"), "Rotate-Shift-Scale image\nIn contrast to rot_scale_trans2D, wrap aroud is not done circulantly so as to\nprevent artifacts from occurring.\nIf the image is a volume, then all slices are\nrotated/translated/scaled.\n \nang - Rotation angle in degrees.\ndelx - Amount to shift rotation origin along x(default=0.0)\ndely - Amount to shift rotation origin along y(default=0.0)\nscale - Scaling factor (default=1.0)\n \nreturn New rotated/shifted/scaled image\nexception - ImageDimensionException can not rotate 1 D image\nexception - ImageDimensionException can not rotate 3 D image")[ return_value_policy< manage_new_object >() ])
	.def("rot_scale_trans", &EMAN::EMData::rot_scale_trans, return_value_policy< manage_new_object >(), args("RA"), "Rotate-Shift-Scale-Circulantly image\nIf the image is a volume, then all slices are\nrotated/translated/scaled.\n \nRA - Transform object\n \nreturn New rotated/shifted/scaled image\nexception - ImageDimensionException can not rotate 1 D image")
	.def("rot_scale_trans_background", &EMAN::EMData::rot_scale_trans_background, return_value_policy< manage_new_object >(), args("RA"), "Rotate-Shift-Scale image\nIn contrast to rot_scale_trans, wrap around is not done circulantly\nso as to prevent artifacts occurring during rotation.\nIf the image is a volume, then all slices are\nrotated/translated/scaled.\n \nRA - Transform object\n \nreturn New rotated/shifted/scaled image\nexception - ImageDimensionException can not rotate 1 D image")
	.def("cm_euc", &EMAN::EMData::cm_euc, args("sinoj", "n1", "n2"), "euclidean distance between two line\n \nsinoj - \nn1 - \nn2 - ")
	.def("rot_scale_conv", &EMAN::EMData::rot_scale_conv, EMAN_EMData_rot_scale_conv_overloads_4_5(args("ang", "delx", "dely", "kb", "scale"), "Rotate-Shift-Scale-Circulantly image using convolution\nIf the image is a volume, then all slices are rotated/translated/scaled.\n \nang - Rotation angle in degrees.\ndelx - Amount to shift rotation origin along x\ndely - Amount to shift rotation origin along y\nkb - convolution kernel\nscale - Scaling factor (default=1.0)\n \nreturn New rotated/shifted/scaled image\nexception - ImageDimensionException can not rotate 1 D image\nexception - ImageDimensionException can not rotate 3 D image")[ return_value_policy< manage_new_object >() ])
	.def("rot_scale_conv7", &EMAN::EMData::rot_scale_conv7, return_value_policy< manage_new_object >(), args("ang", "delx", "dely", "kb", "scale_input"), " ")
	.def("rot_scale_conv_new", &EMAN::EMData::rot_scale_conv_new, EMAN_EMData_rot_scale_conv_new_overloads_4_5(args("ang", "delx", "dely", "kb", "scale"), " ")[ return_value_policy< manage_new_object >() ])
	.def("rot_scale_conv_new_3D", &EMAN::EMData::rot_scale_conv_new_3D,	EMAN_EMData_rot_scale_conv_new_3D_overloads_7_9(args("phi", "theta", "psi", "delx", "dely", "delz", "kb", "scale", "wrap"), " ")[ return_value_policy< manage_new_object >() ])
	.def("rot_scale_conv_new_background", &EMAN::EMData::rot_scale_conv_new_background, EMAN_EMData_rot_scale_conv_new_background_overloads_4_5(args("ang", "delx", "dely", "kb", "scale"), "")[return_value_policy< manage_new_object >()])
	.def("rot_scale_conv_new_background_twice", &EMAN::EMData::rot_scale_conv_new_background_twice, EMAN_EMData_rot_scale_conv_new_background_twice_overloads_4_5(args("ang", "delx", "dely", "kb", "scale"), "")[return_value_policy< manage_new_object >()])
	.def("rot_scale_conv_new_background_3D", &EMAN::EMData::rot_scale_conv_new_background_3D, EMAN_EMData_rot_scale_conv_new_background_3D_overloads_7_9(args("phi", "theta", "psi", "delx", "dely", "delz", "kb", "scale", "wrap"), "")[return_value_policy< manage_new_object >()])
	.def("downsample", &EMAN::EMData::downsample, EMAN_EMData_downsample_overloads_1_2(args("kb", "scale"), " ")[ return_value_policy< manage_new_object >() ])
	.def("get_pixel_conv", &EMAN::EMData::get_pixel_conv, args("delx", "dely", "delz", "kb"), "Get pixel value image using convolution\nIf the image is a volume, then all slices are\nrotated/translated/scaled.\n \ndelx - Amount to shift rotation origin along x\ndely - Amount to shift rotation origin along y\ndelz - Amount to shift rotation origin along z\nkb - convolution kernel\n \nreturn New rotated/shifted/scaled image\nexception - ImageDimensionException can not rotate 1 D image")
	.def("get_pixel_conv7", &EMAN::EMData::get_pixel_conv7, args("delx", "dely", "delz", "kb"), " ")
	.def("getconvpt2d_kbi0", &EMAN::EMData::getconvpt2d_kbi0, EMAN_EMData_getconvpt2d_kbi0_overloads_3_4(args("x", "y", "win", "size"), "Value of 2-D analytic masking (or 2-D convolution) at off-grid point.\nThe only requirement for the window function object is that\nit overload operator()(const float) and return a float.\n \nx - x-value of the desired (potentially off-grid) point\ny - y-value of the desired (potentially off-grid) point\nwin - Window (mask/kernel) function object.\nsize - Size of real-space kernel/mask.\n \nreturn Value of masked/convolved image at (x,y)"))
	.def("fft_shuffle", &EMAN::EMData::fft_shuffle, "fft_shuffle -- Shuffle a Fourier image to put the origin at (0,ny/2)\nOur usual FFT convention puts the origin at (0,0), but then\ngrid points corresponding to iy > ny/2 correspond to\n(unnormalized) frequencies iy-ny.  This routine rearranges\nthe columns of the Fourier image so that iy varies from\n-ny/2 to ny/2 (or ny/2 - 1 for ny even).  This method acts\nas a toggle, so to unshuffle a Fourier image just call\nthis method a second time.")
	.def("extractpoint", &EMAN::EMData::extractpoint, args("xin", "yin", "kb"), "extractpoint -- Gridding convolution\nNote: Expected to be used in combination with fouriergridrot2d.\nsee P.A. Penczek, R. Renka, and H. Schomberg, J. Opt. Soc. Am. A _21_, (2004)\n \nxin - x-position\nyin - y-position\nkb - Kaiser-Bessel window\n \nreturn Complex gridding result")
	.def("extract_plane", &EMAN::EMData::extract_plane, return_value_policy< manage_new_object >(), args("tf", "kb"), "extractplane -- Gridding convolution in 3D along a plane\nNote: Expected to be used in combination with fourier gridding projections.\nsee P.A. Penczek, R. Renka, and H. Schomberg, J. Opt. Soc. Am. A _21_, 499-509 (2004)\n \ntf - transform matrix defining the intended plane.\nkb - Kaiser-Bessel window\n \nreturn Complex gridding plane")
	.def("extract_plane_rect", &EMAN::EMData::extract_plane_rect, return_value_policy< manage_new_object >(), args("tf", "kbx","kby","kbz"), "extractplane square fft plane from 3d rectangualr fft volume -- Gridding convolution in 3D along a plane\nNote: Expected to be used in combination with fourier gridding projections.\nsee P.A. Penczek, R. Renka, and H. Schomberg, J. Opt. Soc. Am. A _21_, 499-509 (2004)\n \ntf - transform matrix defining the intended plane.\nkb - Kaiser-Bessel window\n \nreturn Complex gridding plane")
	.def("extract_plane_rect_fast", &EMAN::EMData::extract_plane_rect_fast, return_value_policy< manage_new_object >(), args("tf", "kbx","kby","kbz"), "extractplane rectangular fft plane from 3d rectangualr fft volume and pad to square after ifft -- Gridding convolution in 3D along a plane\nNote: Expected to be used in combination with fourier gridding projections.\nsee P.A. Penczek, R. Renka, and H. Schomberg, J. Opt. Soc. Am. A _21_, 499-509 (2004)\n \ntf - transform matrix defining the intended plane.\nkb - Kaiser-Bessel window\n \nreturn Complex gridding plane")
	.def("fouriergridrot2d", &EMAN::EMData::fouriergridrot2d, return_value_policy< manage_new_object >(), args("ang", "scale", "kb"), " ")
	.def("fouriergridrot_shift2d", &EMAN::EMData::fouriergridrot_shift2d, return_value_policy< manage_new_object >(), args("ang", "scale", "kb"), " ")
	.def("delete_disconnected_regions", &EMAN::EMData::delete_disconnected_regions, EMAN_EMData_delete_disconnected_regions_overloads_0_3(args("ix","iy", "iz"), "Delete disconnected regions in a binary image\nWorks only for a volume.\n \nix - x coordinate (with respect to the center) from which the search of the compact region begins(default=0).\niy - y coordinate (with respect to the center) from which the search of the compact region begins(default=0).\niz - z coordinate (with respect to the center) from which the search of the compact region begins(default=0).\n \nreturn New binary image")[return_value_policy< manage_new_object >()])
	.def("helicise", &EMAN::EMData::helicise, EMAN_EMData_helicise_overloads_3_6(args("pixel_size", "dp", "dphi", "section_use", "radius"), "Apply helical symmetry\nWorks only for a volume.\n \npixel_size - pixel size in Angstroms.\ndp - repeat in z direction in Angstroms.\ndphi - angular repeat in degrees.\nsection_use - how much of z section to use for symmetrization (between zero and one)(default=1.0).\nradius - radius of the structure (default nx/2-1)(default=-1.0).\ninrad - minimum radius of the structure (default 0.0)(default=-1.0).\n \nreturn New image")[return_value_policy< manage_new_object >()])
	.def("helicise_grid", &EMAN::EMData::helicise_grid, EMAN_EMData_helicise_grid_overloads_4_7(args("pixel_size", "dp", "dphi", "kb", "section_use", "rmax", "rmin"), "Apply helical symmetry\nWorks only for a volume.\n \npixel_size - pixel size in Angstroms.\ndp - repeat in z direction in Angstroms.\ndphi - angular repeat in degrees.\nsection_use - how much of z section to use for symmetrization (between zero and one)(default=1.0).\nradius - radius of the structure (default nx/2-1)(default=-1.0).\ninrad - minimum radius of the structure (default 0.0)(default=-1.0).\n \nreturn New image")[return_value_policy< manage_new_object >()])
	.def("divkbsinh", &EMAN::EMData::divkbsinh, args("kb"), "Divide image by a Kaiser-Bessel sinh window.\nNote: Ideally this method really should be a 'processor'\ninstead, but at the moment a KaiserBessel object\ncannot be passed as part of a Dict, making the usual\nEMData::project() interface rather awkward here.\n \nkb - Kaiser-Bessel window object")
	.def("divkbsinh_rect", &EMAN::EMData::divkbsinh_rect, args("kbx","kby","kbz"), "Divide image by a Kaiser-Bessel sinh in the rect case window.\nNote: Ideally this method really should be a 'processor'\ninstead, but at the moment a KaiserBessel object\ncannot be passed as part of a Dict, making the usual\nEMData::project() interface rather awkward here.\n \nkb - Kaiser-Bessel window object")
	.def("peak_search", &EMAN::EMData::peak_search, args("ml", "invert"), "Search specified number peaks in 1D, 2D, or 3D real images.\nand output the peaks in descendent order:\nThe numbers coming out are: image dimension, then\n1D: pixel value, x coord, relative peak value, x coord( NX/2 center),...\n2D: pixel value, x coord, y coord, realative peak value, x coord(NX/2 center) y coord(NY/2 center)...\n3D  pixel value, x coord, y coord, z coord, realative peak value, x coord(NX/2 center) y coord(NY/2 center) z coord(NZ/2 center)...\nThe function is supposed to return 0 dimension and first pixel value (0,0,0) when the image is constant.\n \nml - \ninvert - ")
	.def("phase_cog", &EMAN::EMData::phase_cog, "Calculate the Phase approximation to center of gravity\nThis operations works for 1-2-3-d images.\n \nreturns both the center of gravity and the phase approximated center of gravity values.")
	.def("find_3d_threshold", &EMAN::EMData::find_3d_threshold, args("mass", "pixel_size"), " ")
	.def("peak_ccf", &EMAN::EMData::peak_ccf,args("hf_p"), "Peak (with a radius of hf_p) search for particle picking\n \nhf_p - ")
#ifdef	DEBUG
	.def("debug_print_params", &EMAN::EMData::debug_print_parms, "Printing EMData params for debugging purpose.")
#endif	//DEBUG
	.def("get_pow", &EMAN::EMData::get_pow, return_value_policy< manage_new_object >(), args("n_pow"), "pixel power operation function\n \nn_pow - ")
	.def("conjg", &EMAN::EMData::conjg, return_value_policy< manage_new_object >(), "pixel conjugate operation function")
	.def("extractline", &EMAN::EMData::extractline, return_value_policy< manage_new_object >(), args("kb", "nuxnew", "nuynew"), " ")
	.def("center_origin", &EMAN::EMData::center_origin)
	.def("center_origin_yz", &EMAN::EMData::center_origin_yz)
	.def("center_origin_fft", &EMAN::EMData::center_origin_fft, "Multiply a Fourier image by (-1)**(ix+iy+iz) to center it.")
	.def("depad", &EMAN::EMData::depad, "De-pad, and and remove Fourier extension convenience function.\nPurpose: De-pad, and and remove Fourier extension from a real image.\nMethod: Remove padding and extension along x for fft, and return the new  image.\n \nreturn depadded input image.")
	.def("depad_corner", &EMAN::EMData::depad_corner, "De-pad, and and remove Fourier extension convenience function.\nPurpose: De-pad, and and remove Fourier extension from a real image.\nMethod: Remove padding and extension along x for fft, and return the new  image.\n \nreturn depadded input image.")
	.def("FourInterpol", &EMAN::EMData::FourInterpol, EMAN_EMData_FourInterpol_overloads_1_4(args("nxni", "nyni", "nzni", "RetReal"), " ")[ return_value_policy< manage_new_object >() ])
	.def("FourTruncate", &EMAN::EMData::FourTruncate, EMAN_EMData_FourTruncate_overloads_1_4(args("nxni", "nyni", "nzni", "RetReal"), "Truncate Fourier transform of an image, it will reduce its size.  (It is a form of decimation).\n \nnxni - new x size (has to be larger/equal than the original x size)\nnyni - new y size (has to be larger/equal than the original y size)(default=0)\nnzni new z size (has to be larger/equal than the original z size)(default=0)\nRetReal - (default=True)\n \nreturn New truncated up image.")[ return_value_policy< manage_new_object >() ])
//       .def("FourInterpol_i", &EMAN::EMData::FourInterpol_i, EMAN_EMData_FourInterpol_i_overloads_1_4()[ return_value_policy< manage_new_object >() ])
	.def("Four_ds", &EMAN::EMData::Four_ds, EMAN_EMData_Four_ds_overloads_1_4(args("nxni", "nyni", "nzni", "RetReal"), "nxni - new x size (has to be larger/equal than the original x size)\nnyni - new y size (has to be larger/equal than the original y size)(default=0)\nnzni - new z size (has to be larger/equal than the original z size)(default=0)\nRetReal - (default=True)")[ return_value_policy< manage_new_object >() ])
	.def("Four_shuf_ds_cen_us", &EMAN::EMData::Four_shuf_ds_cen_us, EMAN_EMData_Four_shuf_ds_cen_us_overloads_1_4(args("nxni", "nyni", "nzni", "RetReal"), "nxni - new x size (has to be larger/equal than the original x size)\nnyni - new y size (has to be larger/equal than the original y size)(default=0)\nnzni new z size (has to be larger/equal than the original z size)(default=0)\nRetReal - (default=True)")[ return_value_policy< manage_new_object >() ])
	.def("filter_by_image", &EMAN::EMData::filter_by_image, EMAN_EMData_filter_by_image_overloads_1_2(args("image", "RetReal"), " ")[ return_value_policy< manage_new_object >() ])
	.def("replace_amplitudes", &EMAN::EMData::replace_amplitudes, EMAN_EMData_replace_amplitudes_overloads_1_2(args("image", "RetReal"), " ")[ return_value_policy< manage_new_object >() ])
	.def("norm_pad", &EMAN::EMData::norm_pad, EMAN_EMData_norm_pad_overloads_2_3(args("do_norm", "npad", "valtype"), "Normalize, pad, and Fourier extend convenience function.\nPurpose: Create a new [normalized] [zero-padded] Fourier image.\nMethod: Normalize (if requested), pad with zeros (if\nrequested), extend along x for fft, and return the new  image.\n \ndo_norm - If true then perform normalization.\nnpad - Amount of zero-padding to use (defaults to 2 if do_pad is true).(default=1)\nvaltype - (default=0)\n \nreturn [normalized,] [zero-padded,] [ft-extended] input image.") [ return_value_policy< manage_new_object >()])
	.def("get_clip", &EMAN::EMData::get_clip, EMAN_EMData_get_clip_overloads_1_2(args("area", "fill"), "Get an inclusive clip. Pads to fill if larger than this image.\n \narea - The clip area, can be 2D/3D.\nfill - the value to assign new pixels outside the area of the original image.(default=0)\n \nreturn The clip image.\nexception ImageDimensionException if any of the dimensions of the argument region are negative") [ return_value_policy< manage_new_object >()])
	.def("clip_inplace", &EMAN::EMData::clip_inplace, EMAN_EMData_clip_inplace_overloads_1_2(args("region", "fill_value"), "Clip the image inplace - clipping region must be smaller than the current region.\ninternally memory is reallocated.\n \nregion - The clip area, can be 2D/3D.\nfill_value - the value fill that new region. default to 0.")[return_value_policy< reference_existing_object >()])
	.def("get_top_half", &EMAN::EMData::get_top_half, return_value_policy< manage_new_object >(), "Get the top half of this 3D image.\n \nreturn The top half of this image.\nexception - ImageDimensionException If this image is not 3D.")
	.def("get_rotated_clip", &EMAN::EMData::get_rotated_clip, EMAN_EMData_get_rotated_clip_overloads_2_3(args("xform", "size", "scale"), "This will extract an arbitrarily oriented and sized region from the image.\n \nxform - The transformation of the region.\nsize - Size of the clip.\nscale - Scaling put on the returned image(default=1.0).\n \nreturn The clip image.")[ return_value_policy< manage_new_object >() ])
	.def("window_center", &EMAN::EMData::window_center, return_value_policy< manage_new_object >(), args("l"), "Window the center of an image.\nOften an image is padded with zeros for fourier interpolation.  In\nthat case the desired lxlxl volume (or lxl area) lies in the center\nof a larger volume (or area).  This routine creates a new object\nthat contains only the desired window.  (This routine is a thin\nwrapper around get_clip.)\n \nl - Length of the window.\n \nreturn An image object that has been windowed.")
	.def("scale", &EMAN::EMData::scale, args("scale_factor"), "scale the image by a factor.\n \nscale_factor - scale factor.")
//	.def("zero_corner_circulant", &EMAN::EMData::zero_corner_circulant, EMAN_EMData_zero_corner_circulant_overloads_0_1(args("radius"), "Zero the pixels in the bottom left corner of the image\nIf radius is greater than 1, than circulant zeroing occurs\nassuming that the center of operation starts in the bottom left\ncorner and proceed outwards to the NE and backwards in a circulant\nfashion towards the SW.\nIntended to zero the area corresponding to the middle of the image,\nas generated by calc_ccf\n \nradius - the radius to the zeroing operation(default=0)\n \nexception - ImageDimensionException if nx > 1 and nx < 2*radius + 1\nexception - ImageDimensionException if ny > 1 and ny < 2*radius + 1\nexception - ImageDimensionException if nz > 1 and nz < 2*radius + 1"))
	.def("translate", (void (EMAN::EMData::*)(float, float, float) )&EMAN::EMData::translate, args("dx", "dy", "dz"), "Translate this image.\n \ndx - Translation distance in x direction.\ndy - Translation distance in y direction.\ndz - Translation distance in z direction.")
	.def("translate", (void (EMAN::EMData::*)(const EMAN::Vec3f&) )&EMAN::EMData::translate, args("translation"), "Translate this image.\n \ntranslation - The translation distance vector.")
	.def("translate", (void (EMAN::EMData::*)(int, int, int) )&EMAN::EMData::translate, args("dx", "dy", "dz"), "Translate this image. integer only translation\ncould be done faster, without interpolation.\n \ndx - Translation distance in x direction.\ndy - Translation distance in y direction.\ndz - Translation distance in z direction.")
	.def("translate", (void (EMAN::EMData::*)(const EMAN::Vec3i&) )&EMAN::EMData::translate, args("translation"), "Translate this image. integer only translation\ncould be done faster, without interpolation.\n \ntranslation - The translation distance vector.")
	.def("max_3D_pixel_error", &EMAN::EMData::max_3D_pixel_error, args("t1", "t2", "r"), "")
//	.def("rotate", (void (EMAN::EMData::*)(const EMAN::Transform3D&) )&EMAN::EMData::rotate, args("t"), "Rotate this image.\nDEPRECATED USE EMData::Transform\n \nt - Transformation rotation.")
	.def("rotate", (void (EMAN::EMData::*)(float, float, float) )&EMAN::EMData::rotate, args("az", "alt", "phi"), "Rotate this image.\nDEPRECATED USE EMData::Transform\n \naz - Rotation euler angle az  in EMAN convention.\nalt - Rotation euler angle alt in EMAN convention.\nphi - Rotation euler angle phi in EMAN convention.")
//	.def("rotate_translate", (void (EMAN::EMData::*)(const EMAN::Transform3D&) )&EMAN::EMData::rotate_translate, args("t"), "Rotate then translate the image.\nDEPRECATED USE EMData::Transform\n \nt - The rotation and translation transformation to be done.")
	.def("transform", &EMAN::EMData::transform, args("t"), "Transform the image\n \nt - the transform object that describes the transformation to be applied to the image.")
	.def("rotate_translate", (void (EMAN::EMData::*)(const EMAN::Transform&) )&EMAN::EMData::rotate_translate, args("t"), "Apply a transformation to the image.\nDEPRECATED USE EMData::Transform\n \nt - transform object that describes the transformation to be applied to the image.")
	.def("rotate_translate", (void (EMAN::EMData::*)(float, float, float, float, float, float) )&EMAN::EMData::rotate_translate, args("az", "alt", "phi", "dx", "dy", "dz"), "Rotate then translate the image.\nDEPRECATED USE EMData::Transform\n \naz - Rotation euler angle az  in EMAN convention.\nalt - Rotation euler angle alt in EMAN convention.\nphi - Rotation euler angle phi in EMAN convention.\ndx - Translation distance in x direction.\ndy - Translation distance in y direction.\ndz - Translation distance in z direction.")
	.def("rotate_translate", (void (EMAN::EMData::*)(float, float, float, float, float, float, float, float, float) )&EMAN::EMData::rotate_translate, args("az", "alt", "phi", "dx", "dy", "dz", "pdx", "pdy", "pdz"), "Rotate then translate the image.\nDEPRECATED USE EMData::Transform\n \naz - Rotation euler angle az  in EMAN convention.\nalt - Rotation euler angle alt in EMAN convention.\nphi - Rotation euler angle phi in EMAN convention.\ndx - Translation distance in x direction.\ndy - Translation distance in y direction.\ndz - Translation distance in z direction.\npdx - Pretranslation distance in x direction.\npdy - Pretranslation distance in y direction.\npdz - Pretranslation distance in z direction.")
	.def("rotate_x", &EMAN::EMData::rotate_x, args("dx"), "This performs a translation of each line along x with wraparound.\nThis is equivalent to a rotation when performed on 'unwrapped' maps.\n \ndx - Translation distance align x direction.\n \nexception - ImageDimensionException If the image is 3D.")
	.def("rotate_180", &EMAN::EMData::rotate_180, "Fast rotation by 180 degrees. Square 2D image only.\n \nexception - ImageFormatException If the image is not square.\nexception - ImageDimensionException If the image is not 2D.")
	.def("dot_rotate_translate", &EMAN::EMData::dot_rotate_translate, EMAN_EMData_dot_rotate_translate_overloads_4_5(args("with", "dx", "dy", "da", "mirror"), "dot product of 2 images. Then 'this' image is rotated/translated.\nIt is much faster than Rotate/Translate then dot product.\n2D images only.\n \nwith - The image used to do the dot product.\ndx - Translation distance in x direction.\ndy - Translation distance in y direction.\nda - Rotation euler angle in degrees\nmirror - (default=False)\n \nexception - ImageFormatException If the 2 images are not the same size.\nexception - ImageDimensionException If the image is 3D."))
	.def("little_big_dot", &EMAN::EMData::little_big_dot, EMAN_EMData_little_big_dot_overloads_1_2(args("little_img", "do_sigma"), "This does a normalized dot product of a little image with a big image\nusing real-space methods. The result is the same size as 'this',\nbut a border 1/2 the size of 'little_img' will be zero.\nThis routine is only efficient when 'little_img' is fairly small.\n \nlittle_img - A small image.\ndo_sigma - Calculate sigma or not(default=False).\n \nreturn normalized dot product image.\nexception - ImageDimensionException If the image is not 1D/2D.")[ return_value_policy< manage_new_object >() ])
	.def("do_radon", &EMAN::EMData::do_radon, return_value_policy< manage_new_object >(), "Radon Transform: an algorithm that transforms an original\nimage into a series of equiangular projections. When\napplied to a 2D object, the output of the Radon transform is a\nseries of 1D lines.\nDo radon transformation on this image. This image must be 2D square.\n \nreturn Radon transform image in square.\nexception - ImageFormatException If the image is not square.\nexception - ImageDimensionException If the image is not 2D.")
	.def("calc_ccf", &EMAN::EMData::calc_ccf, EMAN_EMData_calc_ccf_overloads_0_3(args("with", "fpflag", "center"), "Calculate Cross-Correlation Function (CCF).\nCalculate the correlation of two 1-, 2-, or 3-dimensional images.\nNote: this method internally just calls the correlation function from fundamentals.h.\n \nwith - The image used to calculate the CCF. If 'with' is NULL, the autocorrelation function is computed instead.\nfpflag - Specify how periodicity (or normalization) should be handled. See fundamentals.h for specific flags.(default = 'CIRCULANT').\ncenter - whether or not to center the image (bring bottom left corner to center)(default=False)\n \nreturn Real-space image.\nexception - ImageDimensionException if nx > 1 and nx < 2*radius + 1")[ return_value_policy< manage_new_object >() ])
	.def("calc_ccfx", &EMAN::EMData::calc_ccfx, EMAN_EMData_calc_ccfx_overloads_1_4(args("with", "y0", "y1", "nosum"), "Calculate Cross-Correlation Function (CCF) in the x-direction and adds them up,\nresult in 1D.\nWARNING: this routine will modify the 'this' and 'with' to contain\n1D fft's without setting some flags. This is an optimization\nfor rotational alignment.\nsee calc_ccf()\n \nwith - The image used to calculate CCF.\ny0 - Starting position in x-direction(default=0).\ny1 - Ending position in x-direction. '-1' means the end of the row.(default=-1)\nnosum - If true, returns an image y1-y0+1 pixels high.(default=False)\n \nreturn The result image containing the CCF.\nexception - NullPointerException If input image 'with' is NULL.\nexception - ImageFormatException If 'with' and 'this' are not same size.\nexception - ImageDimensionException If 'this' image is 3D.")[ return_value_policy< manage_new_object >() ])
	.def("calc_fast_sigma_image",&EMAN::EMData::calc_fast_sigma_image, return_value_policy< manage_new_object >(), args("mask"), "Calculates the local standard deviation (sigma) image using the given\nmask image. The mask image is typically much smaller than this image,\nand consists of ones, or is a small circle consisting of ones. The extent\nof the non zero neighborhood explicitly defines the range over which\nthe local standard deviation is determined.\nFourier convolution is used to do the math, ala Roseman (2003, Ultramicroscopy)\nHowever, Roseman was just working on methods Van Heel had presented earlier.\nThe normalize flag causes the mask image to be processed so that it has a unit sum.\nWorks in 1,2 and 3D\n \nmask - the image that will be used to define the neighborhood for determine the local standard deviation\n \nreturn the sigma image, the phase origin is at the corner (not the center)\nexception - ImageDimensionException if the dimensions of with do not match those of this\nexception - ImageDimensionException if any of the dimensions sizes of with exceed of this image's.")
	.def("make_rotational_footprint", &EMAN::EMData::make_rotational_footprint, EMAN_EMData_make_rotational_footprint_overloads_0_1(args("unwrap"), "Makes a 'rotational footprint', which is an 'unwound'\nautocorrelation function. generally the image should be\nedge-normalized and masked before using this.\n \nunwrap - RFP undergoes polar->cartesian x-form,(default=True)\n \nreturn The rotaional footprint image.\nexception - ImageFormatException If image size is not even.")[ return_value_policy< manage_new_object >() ])
	.def("make_rotational_footprint_e1", &EMAN::EMData::make_rotational_footprint_e1, EMAN_EMData_make_rotational_footprint_e1_overloads_0_1(args("unwrap"), "unwrap - RFP undergoes polar->cartesian x-form,(default=True)")[ return_value_policy< manage_new_object >() ])
	.def("make_rotational_footprint_cmc", &EMAN::EMData::make_rotational_footprint_cmc, EMAN_EMData_make_rotational_footprint_cmc_overloads_0_1(args("unwrap"), "unwrap - RFP undergoes polar->cartesian x-form,(default=True)")[ return_value_policy< manage_new_object >() ])
	.def("make_footprint", &EMAN::EMData::make_footprint, EMAN_EMData_make_footprint_overloads_0_1(args("type"), "Makes a 'footprint' for the current image. This is image containing\na rotational & translational invariant of the parent image. The size of the\nresulting image depends on the selected type.\ntype 0- The original, default footprint derived from the rotational footprint\ntypes 1-6 - bispectrum-based\ntypes 1,3,5 - returns Fouier-like images\ntypes 2,4,6 - returns real-space-like images\ntype 1,2 - simple r1,r2, 2-D footprints\ntype 3,4 - r1,r2,anle 3D footprints\ntype 5,6 - same as 1,2 but with the cube root of the final products used\n \ntype - Select one of several possible algorithms for producing the invariants\n \nreturn The footprint image.\nexception - ImageFormatException If image size is not even.")[return_value_policy< manage_new_object >()])
	.def("calc_mutual_correlation", &EMAN::EMData::calc_mutual_correlation, EMAN_EMData_calc_mutual_correlation_overloads_1_3(args("with", "tocorner", "filter"), "Calculates mutual correlation function (MCF) between 2 images.\nIf 'with' is NULL, this does mirror ACF.\n \nwith - The image used to calculate MCF.\ntocorner - Set whether to translate the result image to the corner.(default=False)\nfilter - The filter image used in calculating MCF.(default=Null)\n \nreturn Mutual correlation function image.\nexception - ImageFormatException If 'with' is not NULL and it doesn't have the same size to 'this' image.\nexception NullPointerException If FFT returns NULL image.")[ return_value_policy< manage_new_object >() ])
	.def("unwrap", &EMAN::EMData::unwrap, EMAN_EMData_unwrap_overloads_0_7(args("r1", "r2", "xs", "dx", "dy", "do360", "weight_radial"), "Maps to polar coordinates from Cartesian coordinates. Optionaly radially weighted.\nWhen used with RFP, this provides 1 pixel accuracy at 75% radius.\n2D only.\n \nr1 - (default=-1)\nr2 - (default=-1)\nxs - (deffault=-1)\ndx - (default=0)\ndy - (default=0)\ndo360 - (default=False)\nweight_redial - (default=True)\n \nreturn The image in Cartesian coordinates.\nxception - ImageDimensionException If 'this' image is not 2D.\nexception - UnexpectedBehaviorException if the dimension of this image and the function arguments are incompatibale - i.e. the return image is less than 0 in some dimension.")[ return_value_policy< manage_new_object >() ])
	.def("apply_radial_func", &EMAN::EMData::apply_radial_func, EMAN_EMData_apply_radial_func_overloads_3_4(args("x0", "dx", "array", "interp"), "multiplies by a radial function in fourier space.\n \nx0 - starting point x coordinate.\ndx - step of x.\narray - radial function data array.\ninterp Do the interpolation or not.(default=True)"))
	.def("calc_radial_dist", (std::vector<float,std::allocator<float> > (EMAN::EMData::*)(int, float, float, int) )&EMAN::EMData::calc_radial_dist, args("n", "x0", "dx", "inten"), "calculates radial distribution. works for real and imaginary images.\ninten=0->mean amp, 1->mean inten (amp^2), 2->min, 3->max, 4->sigma. Note that the complex\norigin is at (0,0), with periodic boundaries. Note that the inten option is NOT\nequivalent to returning amplitude and squaring the result.\n \nn - number of points.\nx0 - starting point x coordinate.\ndx - step of x.\ninten returns intensity (amp^2) rather than amplitude if set\n \nreturn The radial distribution in an array.")
	.def("calc_radial_dist", (std::vector<float,std::allocator<float> > (EMAN::EMData::*)(int, float, float, int, float, bool) )&EMAN::EMData::calc_radial_dist, args("n", "x0", "dx", "nwedge", "offset", "inten"), "calculates radial distribution subdivided by angle. works for real and imaginary images.\n2-D only. The first returns a single vector of n*nwedge points, with radius varying first.\nThat is, the first n points represent the radial profile in the first wedge.\n \nn - number of points.\nx0 - starting x coordinate.\ndx - step of x.\nnwedge - int number of wedges to divide the circle into\noffset - angular offset in radians for start of first bin\ninten - returns intensity (amp^2) rather than amplitude if set\n \nreturn nwedge radial distributions packed into a single vector<float>\nexception - ImageDimensionException If 'this' image is not 2D.")
	.def("cconj", &EMAN::EMData::cconj, "Replace the image its complex conjugate.\n \nexception - ImageFormatException Image must be complex (and RI)")
	.def("add_incoherent", &EMAN::EMData::add_incoherent, args("obj"), "Adds 'obj' to 'this' incoherently. 'obj' and 'this' should\nbe same size. Both images should be complex images.\n \nobj - The image added to 'this' image.\n \nexception ImageFormatException If the 2 images are not same size; or if the 2 images are not complex images.")
	.def("calc_hist", &EMAN::EMData::calc_hist, EMAN_EMData_calc_hist_overloads_0_5(args("hist_size", "hist_min", "hist_max", "brt", "cont"), "Calculates the histogram of 'this' image. The result is\nstored in float array 'hist'. If hist_min = hist_max, use\nimage data min as hist_min; use image data max as hist_max.\n \nhist_size - Histogram array's size.(default=128)\nhist_min - Minimum histogram value.(default=0)\nhist_max - Maximum histogram value.(default=0)\nbrt - (default=0.0f)\ncont - (default=1.0f)\n \nreturn histogram array of this image."))
	.def("calc_az_dist", &EMAN::EMData::calc_az_dist, args("n", "a0", "da", "rmin", "rmax"), "Caculates the azimuthal distributions.\nworks for real or complex images, 2D only.\n \nn - Number of elements.\na0 - Starting angle.\nda - Angle step.\nrmin - Minimum radius.\nrmax  Maximum radius.\n \nreturn Float array to store the data.\nexception - ImageDimensionException If image is 3D.")
	.def("calc_dist", &EMAN::EMData::calc_dist, EMAN_EMData_calc_dist_overloads_1_2(args("second_img", "y_index"), "Calculates the distance between 2 vectors. 'this' image is\n1D, which contains a vector; 'second_img' may be nD. One of\nits row is used as the second vector. 'second_img' and\n'this' must have the same x size.\n \nsecond_img The image used to caculate the distance.\ny_index Specifies which row in 'second_img' is used to do the caculation.(default=0)\n \nreturn The distance between 2 vectors.\nexception ImageDimensionException If 'this' image is not 1D.\nexception ImageFormatException If the 2 images don't have same xsize."))
	.def("calc_flcf", &EMAN::EMData::calc_flcf, return_value_policy< manage_new_object >(), args("with"), "Calculates the cross correlation with local normalization\nbetween 2 images. This is a faster version of local correlation\nthat make use of Fourier convolution and correlation.\nWith is the template - the thing that you wish to find in the this image.\nIt should not be unecessarily padded. The function uses the size of with\nto determine the extent of the local neighborhood used in the local\nnormalization (for technical details, see calc_fast_sigma_image).\nNote that this function circularly masks the template at its radius\nso the calling function need not do this beforehand.\nWorks in 1,2 and 3D.\n \nwith - The image used to calculate cross correlation (the template)\n \nreturn the local normalized cross correlation image - the phase origin is at the corner of the image")
	.def("convolute", &EMAN::EMData::convolute, return_value_policy< manage_new_object >(), args("with"), "Convolutes 2 data sets. The 2 images must be of the same size.\n \nwith - One data set. 'this' image is the other data set.\n \nreturn The result image.\nexception - NullPointerException If FFT resturns NULL image.")
	.def("common_lines", &EMAN::EMData::common_lines, EMAN_EMData_common_lines_overloads_2_5(args("image1", "image2", "mode", "steps", "horizontal"), "Finds common lines between 2 complex images.\nThis function does not assume any symmetry, just blindly\ncompute the \"sinogram\" and the user has to take care how\nto interpret the returned \"sinogram\". it only considers\ninplane rotation and assumes prefect centering and identical scale.\n \nimage1 - The first complex image.\nimage2 - The second complex image.\nmode Either 0 or 1 or 2. mode 0 is a summed dot-product, larger value means better match; mode 1 is weighted phase residual, lower value means better match.(default=0)\nsteps - 1/2 of the resolution of the map.(default=180)\nhorizontal - In horizontal way or not.(default=False)\n \nexception - NullPointerException If 'image1' or 'image2' is NULL.\nexception - OutofRangeException If 'mode' is invalid.\nexception - ImageFormatException If 'image1' 'image2' are not same size."))
	.def("common_lines_real", &EMAN::EMData::common_lines_real, EMAN_EMData_common_lines_real_overloads_2_4(args("image1", "image2", "steps", "horizontal"), "Finds common lines between 2 real images.\n \nimage1 - The first image.\nimage2 - The second image.\nsteps - 1/2 of the resolution of the map.(default=180)\nhorizontal In horizontal way or not.(default=False)\n \nexception - NullPointerException If 'image1' or 'image2' is NULL.\nexception - ImageFormatException If 'image1' 'image2' are not same size."))
	.def("cut_slice", &EMAN::EMData::cut_slice, EMAN_EMData_cut_slice_overloads_2_3(args("map", "tr", "interpolate"), "cut a 2D slice out of a real 3D map. Put slice into 'this' image.\n \nmap - The real 3D map.\ntr - orientation of the slice as encapsulated in a Transform object.\ninterpolate - Do interpolation or not.(default=True)\n \nexception - NullPointerException If map is NULL.\nexception - ImageDimensionException If this image is not 2D or 3D.\nexception - ImageFormatException If this image is complex\nexception - ImageFormatException If map is complex"))
	.def("extract_box", &EMAN::EMData::extract_box, return_value_policy< manage_new_object >(), args("transform", "region"), "Extract a box from a subtomogram in a arbitary coord system")
	.def("get_amplitude_thres", &EMAN::EMData::get_amplitude_thres, args("threshold"), "Compute the Fourier Amp where thres % of the amps are lower than this functions return value")
	.def("mask_contig_region",&EMAN::EMData::mask_contig_region, args("val", "seed"), " ")
	.def("uncut_slice", &EMAN::EMData::uncut_slice, args("map", "tr"), "Opposite of the cut_slice(). It will take a slice and insert\nthe data into a real 3D map. It does not interpolate, it uses\nthe nearest neighbor.\n \nmap - The real 3D map.\ntr - Orientation of the slice.\n \nexception - NullPointerException If map is NULL.\nexception - ImageDimensionException If this image is not 2D.\nexception - ImageDimensionException If map image is not 3D.\nexception - ImageFormatException If this image is complex.\nexception - ImageFormatException If map is complex")
	.def("set_xyz_origin", &EMAN::EMData::set_xyz_origin, args("origin_x", "origin_y", "origin_z"), "Set the x,y,z origin of the image\n \norigin_x - the x origin\norigin_y - the y origin\norigin_z - the z origin")
	.def("equal", &EMAN::EMData::equal, args("that"), "compare the equality of two EMData object based on their pixel values")
	.def("getwedge", &EMAN::EMData::compute_missingwedge, EMAN_EMData_compute_missingwedge_overloads_1_3(args("wedgeangle", "start", "stop"), "get the missingt wedge mean ansd std")[ return_value_policy< manage_new_object >() ])
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
	.def("__call__", (float& (EMAN::EMData::*)(const size_t) const )&EMAN::EMData::operator (), return_value_policy< copy_non_const_reference >())
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

