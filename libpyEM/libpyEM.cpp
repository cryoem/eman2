
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <aligner.h>
#include <ctf.h>
#include <emdata.h>
#include <emobject.h>
#include <emutil.h>
#include <filter.h>
#include <imageio.h>
#include <log.h>
#include <projector.h>
#include <pyem.h>
#include <pylist.h>
#include <transform.h>
#include <xydata.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_EMUtil_get_imageio_overloads_2_3, EMAN::EMUtil::get_imageio, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_read_image_overloads_1_5, read_image, 1, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_write_image_overloads_1_4, write_image, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_append_image_overloads_1_3, append_image, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_filter_overloads_1_2, filter, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_align_overloads_2_3, align, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_copy_overloads_0_2, copy, 0, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_mask_normalize_overloads_1_2, mask_normalize, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_edge_normalize_overloads_0_1, edge_normalize, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_normalize_to_overloads_1_5, normalize_to, 1, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_least_square_normalize_to_overloads_1_3, least_square_normalize_to, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_hist_overloads_1_4, calc_hist, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_little_big_dot_overloads_1_2, little_big_dot, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_fast_translate_overloads_0_1, fast_translate, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_rotate_translate_overloads_0_5, rotate_translate, 0, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_fast_rotate_translate_overloads_0_1, fast_rotate_translate, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_ccf_overloads_1_3, calc_ccf, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_make_rotational_footprint_overloads_0_2, make_rotational_footprint, 0, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_ccfx_overloads_1_4, calc_ccfx, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_mutual_correlation_overloads_1_3, calc_mutual_correlation, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_unwrap_overloads_0_6, unwrap, 0, 6)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_add_random_noise_overloads_4_5, add_random_noise, 4, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_auto_mask_overloads_1_2, auto_mask, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_dist_overloads_1_2, calc_dist, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_flcf_overloads_1_3, calc_flcf, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_dot_overloads_1_2, dot, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_common_lines_overloads_2_5, common_lines, 2, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_common_lines_real_overloads_2_4, common_lines_real, 2, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_cut_slice_overloads_2_6, cut_slice, 2, 6)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_uncut_slice_overloads_2_5, uncut_slice, 2, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_to_mass_center_overloads_0_1, to_mass_center, 0, 1)

struct EMAN_Filter_Wrapper: EMAN::Filter
{
    EMAN_Filter_Wrapper(PyObject* self_, const EMAN::Filter& p0):
        EMAN::Filter(p0), self(self_) {}

    EMAN_Filter_Wrapper(PyObject* self_):
        EMAN::Filter(), self(self_) {}

    void process(EMAN::EMData* p0) {
        call_method< void >(self, "process", p0);
    }

    void default_process(EMAN::EMData* p0) {
        EMAN::Filter::process(p0);
    }

    void process_list(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >& p0) {
        call_method< void >(self, "process_list", p0);
    }

    void default_process_list(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >& p0) {
        EMAN::Filter::process_list(p0);
    }

    EMAN::Dict get_params() const {
        return call_method< EMAN::Dict >(self, "get_params");
    }

    EMAN::Dict default_get_params() const {
        return EMAN::Filter::get_params();
    }

    void set_params(const EMAN::Dict& p0) {
        call_method< void >(self, "set_params", p0);
    }

    void default_set_params(const EMAN::Dict& p0) {
        EMAN::Filter::set_params(p0);
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(self, "get_param_types");
    }

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
    }

    PyObject* self;
};

struct EMAN_Aligner_Wrapper: EMAN::Aligner
{
    EMAN_Aligner_Wrapper(PyObject* self_, const EMAN::Aligner& p0):
        EMAN::Aligner(p0), self(self_) {}

    EMAN_Aligner_Wrapper(PyObject* self_):
        EMAN::Aligner(), self(self_) {}

    EMAN::EMData* align(EMAN::EMData* p0, std::string p1) const {
        return call_method< EMAN::EMData* >(self, "align", p0, p1);
    }

    EMAN::Dict get_params() const {
        return call_method< EMAN::Dict >(self, "get_params");
    }

    EMAN::Dict default_get_params() const {
        return EMAN::Aligner::get_params();
    }

    void set_params(const EMAN::Dict& p0) {
        call_method< void >(self, "set_params", p0);
    }

    void default_set_params(const EMAN::Dict& p0) {
        EMAN::Aligner::set_params(p0);
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(self, "get_param_types");
    }

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
    }

    PyObject* self;
};

struct EMAN_Projector_Wrapper: EMAN::Projector
{
    EMAN_Projector_Wrapper(PyObject* self_, const EMAN::Projector& p0):
        EMAN::Projector(p0), self(self_) {}

    EMAN_Projector_Wrapper(PyObject* self_):
        EMAN::Projector(), self(self_) {}

    EMAN::Dict get_params() const {
        return call_method< EMAN::Dict >(self, "get_params");
    }

    EMAN::Dict default_get_params() const {
        return EMAN::Projector::get_params();
    }

    EMAN::EMData* project3d(EMAN::EMData* p0) const {
        return call_method< EMAN::EMData* >(self, "project3d", p0);
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(self, "get_param_types");
    }

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
    }

    PyObject* self;
};


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyEM)
{
    EMAN::vector_to_python<EMAN::EMData*>();
    EMAN::vector_from_python<int>();
    EMAN::vector_from_python<float>();
    EMAN::vector_to_python<std::string>();
    
    EMAN::map_to_python<EMAN::EMObject>();
    EMAN::map_from_python<EMAN::EMObject>();

    EMAN::Dict_to_python();
    EMAN::Dict_from_python();


    scope* EMAN_Log_scope = new scope(
    class_< EMAN::Log, boost::noncopyable >("Log", no_init)
        .def("logger", &EMAN::Log::logger, return_value_policy< reference_existing_object >())
        .def("set_level", &EMAN::Log::set_level)
        .def("set_logfile", &EMAN::Log::set_logfile)
        .staticmethod("logger")
    );

    enum_< EMAN::Log::LogLevel >("LogLevel")
        .value("ERROR_LOG", EMAN::Log::ERROR_LOG)
        .value("NORMAL_LOG", EMAN::Log::NORMAL_LOG)
        .value("WARNING_LOG", EMAN::Log::WARNING_LOG)
        .value("VARIABLE_LOG", EMAN::Log::VARIABLE_LOG)
    ;

    delete EMAN_Log_scope;

    scope* EMAN_EMObject_scope = new scope(
    class_< EMAN::EMObject >("EMObject", init<  >())
        .def(init< const EMAN::EMObject& >())
        .def(init< int >())
        .def(init< float >())
        .def(init< double >())
        .def(init< std::string >())
        .def(init< EMAN::EMData* >())
        .def(init< EMAN::XYData* >())
        .def(init< const std::vector<float,std::allocator<float> >& >())
        .def("get_int", &EMAN::EMObject::get_int)
        .def("get_float", &EMAN::EMObject::get_float)
        .def("get_double", &EMAN::EMObject::get_double)
        .def("get_string", &EMAN::EMObject::get_string)
        .def("get_emdata", &EMAN::EMObject::get_emdata, return_internal_reference< 1 >())
        .def("get_xydata", &EMAN::EMObject::get_xydata, return_internal_reference< 1 >())
        .def("get_farray", &EMAN::EMObject::get_farray)
        .def("is_null", &EMAN::EMObject::is_null)
        .def("to_str", &EMAN::EMObject::to_str)
        .def("get_object_type_name", &EMAN::EMObject::get_object_type_name)
        .staticmethod("get_object_type_name")
    );

    enum_< EMAN::EMObject::ObjectType >("ObjectType")
        .value("EMDATA", EMAN::EMObject::EMDATA)
        .value("STRING", EMAN::EMObject::STRING)
        .value("INT", EMAN::EMObject::INT)
        .value("DOUBLE", EMAN::EMObject::DOUBLE)
        .value("FLOAT", EMAN::EMObject::FLOAT)
        .value("XYDATA", EMAN::EMObject::XYDATA)
        .value("FLOATARRAY", EMAN::EMObject::FLOATARRAY)
        .value("UNKNOWN", EMAN::EMObject::UNKNOWN)
    ;

    delete EMAN_EMObject_scope;

    scope* EMAN_EMUtil_scope = new scope(
    class_< EMAN::EMUtil >("EMUtil", init<  >())
        .def(init< const EMAN::EMUtil& >())
        .def("vertical_acf", &EMAN::EMUtil::vertical_acf, return_value_policy< manage_new_object >())
        .def("make_image_median", &EMAN::EMUtil::make_image_median, return_value_policy< manage_new_object >())
        .def("get_image_type", &EMAN::EMUtil::get_image_type)
        .def("get_image_count", &EMAN::EMUtil::get_image_count)
        .def("get_imageio", &EMAN::EMUtil::get_imageio, return_internal_reference< 1 >(), EMAN_EMUtil_get_imageio_overloads_2_3())
        .def("get_imagetype_name", &EMAN::EMUtil::get_imagetype_name)
        .def("get_datatype_string", &EMAN::EMUtil::get_datatype_string)
        .def("dump_dict", &EMAN::EMUtil::dump_dict)
        .def("is_same_size", &EMAN::EMUtil::is_same_size)
        .staticmethod("vertical_acf")
        .staticmethod("get_datatype_string")
        .staticmethod("dump_dict")
        .staticmethod("get_imageio")
        .staticmethod("get_image_count")
        .staticmethod("get_imagetype_name")
        .staticmethod("get_image_type")
        .staticmethod("is_same_size")
        .staticmethod("make_image_median")
    );

    enum_< EMAN::EMUtil::EMDataType >("EMDataType")
        .value("EM_SHORT_COMPLEX", EMAN::EMUtil::EM_SHORT_COMPLEX)
        .value("EM_SHORT", EMAN::EMUtil::EM_SHORT)
        .value("EM_UCHAR", EMAN::EMUtil::EM_UCHAR)
        .value("EM_FLOAT_COMPLEX", EMAN::EMUtil::EM_FLOAT_COMPLEX)
        .value("EM_CHAR", EMAN::EMUtil::EM_CHAR)
        .value("EM_INT", EMAN::EMUtil::EM_INT)
        .value("EM_USHORT", EMAN::EMUtil::EM_USHORT)
        .value("EM_USHORT_COMPLEX", EMAN::EMUtil::EM_USHORT_COMPLEX)
        .value("EM_UNKNOWN", EMAN::EMUtil::EM_UNKNOWN)
        .value("EM_UINT", EMAN::EMUtil::EM_UINT)
        .value("EM_DOUBLE", EMAN::EMUtil::EM_DOUBLE)
        .value("EM_FLOAT", EMAN::EMUtil::EM_FLOAT)
    ;


    enum_< EMAN::EMUtil::ImageType >("ImageType")
        .value("IMAGE_XPLOR", EMAN::EMUtil::IMAGE_XPLOR)
        .value("IMAGE_MRC", EMAN::EMUtil::IMAGE_MRC)
        .value("IMAGE_GATAN2", EMAN::EMUtil::IMAGE_GATAN2)
        .value("IMAGE_ICOS", EMAN::EMUtil::IMAGE_ICOS)
        .value("IMAGE_AMIRA", EMAN::EMUtil::IMAGE_AMIRA)
        .value("IMAGE_LST", EMAN::EMUtil::IMAGE_LST)
        .value("IMAGE_DM3", EMAN::EMUtil::IMAGE_DM3)
        .value("IMAGE_TIFF", EMAN::EMUtil::IMAGE_TIFF)
        .value("IMAGE_SAL", EMAN::EMUtil::IMAGE_SAL)
        .value("IMAGE_IMAGIC", EMAN::EMUtil::IMAGE_IMAGIC)
        .value("IMAGE_VTK", EMAN::EMUtil::IMAGE_VTK)
        .value("IMAGE_HDF", EMAN::EMUtil::IMAGE_HDF)
        .value("IMAGE_SINGLE_SPIDER", EMAN::EMUtil::IMAGE_SINGLE_SPIDER)
        .value("IMAGE_EMIM", EMAN::EMUtil::IMAGE_EMIM)
        .value("IMAGE_SPIDER", EMAN::EMUtil::IMAGE_SPIDER)
        .value("IMAGE_PNG", EMAN::EMUtil::IMAGE_PNG)
        .value("IMAGE_PGM", EMAN::EMUtil::IMAGE_PGM)
        .value("IMAGE_EM", EMAN::EMUtil::IMAGE_EM)
        .value("IMAGE_PIF", EMAN::EMUtil::IMAGE_PIF)
        .value("IMAGE_UNKNOWN", EMAN::EMUtil::IMAGE_UNKNOWN)
    ;

    delete EMAN_EMUtil_scope;

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
        .def("copy", &EMAN::EMData::copy, return_value_policy< manage_new_object >(), EMAN_EMData_copy_overloads_0_2())
        .def("copy_head", &EMAN::EMData::copy_head, return_value_policy< manage_new_object >())
        .def("get_clip", &EMAN::EMData::get_clip, return_value_policy< manage_new_object >())
        .def("insert_clip", &EMAN::EMData::insert_clip)
        .def("do_fft", &EMAN::EMData::do_fft, return_value_policy< manage_new_object >())
        .def("do_ift", &EMAN::EMData::do_ift, return_value_policy< manage_new_object >())
        .def("gimme_fft", &EMAN::EMData::gimme_fft)
        .def("normalize", &EMAN::EMData::normalize)
        .def("mask_normalize", &EMAN::EMData::mask_normalize, EMAN_EMData_mask_normalize_overloads_1_2())
        .def("edge_normalize", &EMAN::EMData::edge_normalize, EMAN_EMData_edge_normalize_overloads_0_1())
        .def("row_normalize", &EMAN::EMData::row_normalize)
        .def("normalize_max", &EMAN::EMData::normalize_max)
        .def("normalize_slice", &EMAN::EMData::normalize_slice)
        .def("normalize_to", &EMAN::EMData::normalize_to, EMAN_EMData_normalize_to_overloads_1_5())
        .def("least_square_normalize_to", &EMAN::EMData::least_square_normalize_to, EMAN_EMData_least_square_normalize_to_overloads_1_3())
        .def("calc_hist", &EMAN::EMData::calc_hist, EMAN_EMData_calc_hist_overloads_1_4())
        .def("little_big_dot", &EMAN::EMData::little_big_dot, return_value_policy< manage_new_object >(), EMAN_EMData_little_big_dot_overloads_1_2())
        .def("calc_az_dist", &EMAN::EMData::calc_az_dist)
        .def("ri2ap", &EMAN::EMData::ri2ap)
        .def("ap2ri", &EMAN::EMData::ap2ri)
        .def("to_corner", &EMAN::EMData::to_corner)
        .def("rotate_x", &EMAN::EMData::rotate_x)
        .def("rotate_180", &EMAN::EMData::rotate_180)
        .def("fast_translate", &EMAN::EMData::fast_translate, EMAN_EMData_fast_translate_overloads_0_1())
        .def("rotate_translate", &EMAN::EMData::rotate_translate, EMAN_EMData_rotate_translate_overloads_0_5())
        .def("fast_rotate_translate", &EMAN::EMData::fast_rotate_translate, EMAN_EMData_fast_rotate_translate_overloads_0_1())
        .def("dot_rotate_translate", &EMAN::EMData::dot_rotate_translate)
        .def("do_radon", &EMAN::EMData::do_radon, return_value_policy< manage_new_object >())
        .def("vertical_flip", &EMAN::EMData::vertical_flip)
        .def("horizontal_flip", &EMAN::EMData::horizontal_flip)
        .def("calc_ccf", &EMAN::EMData::calc_ccf, return_value_policy< manage_new_object >(), EMAN_EMData_calc_ccf_overloads_1_3())
        .def("make_rotational_footprint", &EMAN::EMData::make_rotational_footprint, return_value_policy< manage_new_object >(), EMAN_EMData_make_rotational_footprint_overloads_0_2())
        .def("calc_ccfx", &EMAN::EMData::calc_ccfx, return_value_policy< manage_new_object >(), EMAN_EMData_calc_ccfx_overloads_1_4())
        .def("calc_mutual_correlation", &EMAN::EMData::calc_mutual_correlation, return_value_policy< manage_new_object >(), EMAN_EMData_calc_mutual_correlation_overloads_1_3())
        .def("unwrap", &EMAN::EMData::unwrap, return_value_policy< manage_new_object >(), EMAN_EMData_unwrap_overloads_0_6())
        .def("calc_fourier_shell_correlation", &EMAN::EMData::calc_fourier_shell_correlation)
        .def("mean_shrink", &EMAN::EMData::mean_shrink)
        .def("median_shrink", &EMAN::EMData::median_shrink)
        .def("apply_radial_func", &EMAN::EMData::apply_radial_func)
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
        .def("auto_mask", &EMAN::EMData::auto_mask, EMAN_EMData_auto_mask_overloads_1_2())
        .def("calc_dist", &EMAN::EMData::calc_dist, EMAN_EMData_calc_dist_overloads_1_2())
        .def("calc_flcf", &EMAN::EMData::calc_flcf, return_value_policy< manage_new_object >(), EMAN_EMData_calc_flcf_overloads_1_3())
        .def("calc_radial_dist", (void (EMAN::EMData::*)(int, float, float, float*) )&EMAN::EMData::calc_radial_dist)
        .def("calc_radial_dist", (void (EMAN::EMData::*)(int, float, float, float*, float, float) )&EMAN::EMData::calc_radial_dist)
        .def("convolute", &EMAN::EMData::convolute, return_value_policy< manage_new_object >())
        .def("has_ctff", &EMAN::EMData::has_ctff)
        .def("dot", &EMAN::EMData::dot, EMAN_EMData_dot_overloads_1_2())
        .def("common_lines", &EMAN::EMData::common_lines, EMAN_EMData_common_lines_overloads_2_5())
        .def("common_lines_real", &EMAN::EMData::common_lines_real, EMAN_EMData_common_lines_real_overloads_2_4())
        .def("cut_slice", &EMAN::EMData::cut_slice, EMAN_EMData_cut_slice_overloads_2_6())
        .def("uncut_slice", &EMAN::EMData::uncut_slice, EMAN_EMData_uncut_slice_overloads_2_5())
        .def("ift_slice", &EMAN::EMData::ift_slice, return_value_policy< manage_new_object >())
        .def("get_edge_mean", &EMAN::EMData::get_edge_mean)
        .def("get_circle_mean", &EMAN::EMData::get_circle_mean)
        .def("radial_average", &EMAN::EMData::radial_average)
        .def("radial_subtract", &EMAN::EMData::radial_subtract)
        .def("sub_noise", &EMAN::EMData::sub_noise)
        .def("setup_insert_slice", &EMAN::EMData::setup_insert_slice)
        .def("to_mass_center", &EMAN::EMData::to_mass_center, EMAN_EMData_to_mass_center_overloads_0_1())
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
        .def("get_mean", &EMAN::EMData::get_mean)
        .def("get_std", &EMAN::EMData::get_std)
        .def("get_min_location", &EMAN::EMData::get_min_location)
        .def("get_max_location", &EMAN::EMData::get_max_location)
        .def("get_min_index", &EMAN::EMData::get_min_index)
        .def("get_max_index", &EMAN::EMData::get_max_index)
        .def("get_x", &EMAN::EMData::get_x)
        .def("get_y", &EMAN::EMData::get_y)
        .def("get_z", &EMAN::EMData::get_z)
        .def("get_dim", &EMAN::EMData::get_dim)
        .def("get_parent", &EMAN::EMData::get_parent, return_value_policy< manage_new_object >())
        .def("set_parent", &EMAN::EMData::set_parent)
        .def("get_name", &EMAN::EMData::get_name)
        .def("set_name", &EMAN::EMData::set_name)
        .def("set_pixel_size", &EMAN::EMData::set_pixel_size)
        .def("get_pixel_size", &EMAN::EMData::get_pixel_size)
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
        .def( self - other< float >() )
        .def( self + other< float >() )
        .def( self - self )
        .def( self + self )
        .def( other< float >() / self )
        .def( other< float >() * self )
        .def( other< float >() - self )
        .def( other< float >() + self )
        .def( self / other< float >() )
        .def( self * other< float >() )
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

    class_< EMAN::FilterFactory, boost::noncopyable >("FilterFactory", no_init)
        .def("instance", &EMAN::FilterFactory::instance, return_value_policy< reference_existing_object >())
        .def("get", (EMAN::Filter* (EMAN::FilterFactory::*)(std::string) )&EMAN::FilterFactory::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Filter* (EMAN::FilterFactory::*)(std::string, const EMAN::Dict&) )&EMAN::FilterFactory::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::FilterFactory::get_list)
        .staticmethod("instance")
    ;

    class_< EMAN::Filter, boost::noncopyable, EMAN_Filter_Wrapper >("Filter", init<  >())
        .def("process", &EMAN::Filter::process, &EMAN_Filter_Wrapper::default_process)
        .def("process_list", &EMAN::Filter::process_list, &EMAN_Filter_Wrapper::default_process_list)
        .def("get_params", &EMAN::Filter::get_params, &EMAN_Filter_Wrapper::default_get_params)
        .def("set_params", &EMAN::Filter::set_params, &EMAN_Filter_Wrapper::default_set_params)
        .def("get_param_types", pure_virtual(&EMAN::Filter::get_param_types))
        .def("get_name", pure_virtual(&EMAN::Filter::get_name))
    ;

    class_< EMAN::AlignerFactory, boost::noncopyable >("AlignerFactory", no_init)
        .def("instance", &EMAN::AlignerFactory::instance, return_value_policy< reference_existing_object >())
        .def("get", (EMAN::Aligner* (EMAN::AlignerFactory::*)(std::string) )&EMAN::AlignerFactory::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Aligner* (EMAN::AlignerFactory::*)(std::string, const EMAN::Dict&) )&EMAN::AlignerFactory::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::AlignerFactory::get_list)
        .staticmethod("instance")
    ;

    class_< EMAN::Aligner, boost::noncopyable, EMAN_Aligner_Wrapper >("Aligner", init<  >())
        .def("align", pure_virtual(&EMAN::Aligner::align), return_value_policy< manage_new_object >())
        .def("get_params", &EMAN::Aligner::get_params, &EMAN_Aligner_Wrapper::default_get_params)
        .def("set_params", &EMAN::Aligner::set_params, &EMAN_Aligner_Wrapper::default_set_params)
        .def("get_param_types", pure_virtual(&EMAN::Aligner::get_param_types))
        .def("get_name", pure_virtual(&EMAN::Aligner::get_name))
    ;

    class_< EMAN::ProjectorFactory, boost::noncopyable >("ProjectorFactory", no_init)
        .def("instance", &EMAN::ProjectorFactory::instance, return_value_policy< reference_existing_object >())
        .def("get", (EMAN::Projector* (EMAN::ProjectorFactory::*)(std::string) )&EMAN::ProjectorFactory::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Projector* (EMAN::ProjectorFactory::*)(std::string, const EMAN::Dict&) )&EMAN::ProjectorFactory::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::ProjectorFactory::get_list)
        .staticmethod("instance")
    ;

    class_< EMAN::Projector, boost::noncopyable, EMAN_Projector_Wrapper >("Projector", init<  >())
        .def("get_params", &EMAN::Projector::get_params, &EMAN_Projector_Wrapper::default_get_params)
        .def("project3d", pure_virtual(&EMAN::Projector::project3d), return_value_policy< manage_new_object >())
        .def("get_param_types", pure_virtual(&EMAN::Projector::get_param_types))
        .def("get_name", pure_virtual(&EMAN::Projector::get_name))
        .def("set_params", &EMAN::Projector::set_params)
    ;

}

