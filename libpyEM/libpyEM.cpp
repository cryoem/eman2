
// Includes ====================================================================
#include <boost/python.hpp>
#include <filter.h>
#include <emobject.h>
#include <log.h>
#include <imageio.h>
#include <emdata.h>
#include <emutil.h>
#include <aligner.h>
#include <ctf.h>
#include <pyem.h>
#include <pylist.h>
#include <transform.h>
#include <projector.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {


struct EMAN_Aligner_Wrapper: EMAN::Aligner
{
    EMAN_Aligner_Wrapper(PyObject* self_, const EMAN::Aligner & p0):
        EMAN::Aligner(p0), self(self_) {}

    EMAN_Aligner_Wrapper(PyObject* self_):
        EMAN::Aligner(), self(self_) {}

    EMAN::EMData * align(EMAN::EMData * p0, std::basic_string<char,std::char_traits<char>,std::allocator<char> > p1) const {
        return call_method< EMAN::EMData * >(self, "align", p0, p1);
    }

    EMAN::Dict get_params() const {
        return call_method< EMAN::Dict >(self, "get_params");
    }

    EMAN::Dict default_get_params() const {
        return EMAN::Aligner::get_params();
    }

    void set_params(const EMAN::Dict & p0) {
        call_method< void >(self, "set_params", p0);
    }

    void default_set_params(const EMAN::Dict & p0) {
        EMAN::Aligner::set_params(p0);
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(self, "get_param_types");
    }

    std::basic_string<char,std::char_traits<char>,std::allocator<char> > get_name() const {
        return call_method< std::basic_string<char,std::char_traits<char>,std::allocator<char> > >(self, "get_name");
    }

    PyObject* self;
};

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMUtil_get_imageio_overloads_2_3, get_imageio, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_copy_overloads_0_2, copy, 0, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_read_image_overloads_1_5, read_image, 1, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_write_image_overloads_1_4, write_image, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_align_overloads_2_3, align, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_setup4slice_overloads_0_1, setup4slice, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_rotate_translate_overloads_0_5, rotate_translate, 0, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_fast_rotate_translate_overloads_0_1, fast_rotate_translate, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_fast_translate_overloads_0_1, fast_translate, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_dot_overloads_1_2, dot, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_ccf_overloads_1_3, calc_ccf, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_make_rotational_footprint_overloads_0_2, make_rotational_footprint, 0, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_ccfx_overloads_1_4, calc_ccfx, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_unwrap_overloads_0_6, unwrap, 0, 6)

struct EMAN_Projector_Wrapper: EMAN::Projector
{
    EMAN_Projector_Wrapper(PyObject* self_, const EMAN::Projector & p0):
        EMAN::Projector(p0), self(self_) {}

    EMAN_Projector_Wrapper(PyObject* self_):
        EMAN::Projector(), self(self_) {}

    EMAN::Dict get_params() const {
        return call_method< EMAN::Dict >(self, "get_params");
    }

    EMAN::Dict default_get_params() const {
        return EMAN::Projector::get_params();
    }

    EMAN::EMData * project3d(EMAN::EMData * p0) const {
        return call_method< EMAN::EMData * >(self, "project3d", p0);
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(self, "get_param_types");
    }

    std::basic_string<char,std::char_traits<char>,std::allocator<char> > get_name() const {
        return call_method< std::basic_string<char,std::char_traits<char>,std::allocator<char> > >(self, "get_name");
    }

    PyObject* self;
};

struct EMAN_Filter_Wrapper: EMAN::Filter
{
    EMAN_Filter_Wrapper(PyObject* self_, const EMAN::Filter & p0):
        EMAN::Filter(p0), self(self_) {}

    EMAN_Filter_Wrapper(PyObject* self_):
        EMAN::Filter(), self(self_) {}

    void process(EMAN::EMData * p0) {
        call_method< void >(self, "process", p0);
    }

    void default_process(EMAN::EMData * p0) {
        EMAN::Filter::process(p0);
    }

    void process_list(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> > & p0) {
        call_method< void >(self, "process_list", p0);
    }

    void default_process_list(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> > & p0) {
        EMAN::Filter::process_list(p0);
    }

    EMAN::Dict get_params() const {
        return call_method< EMAN::Dict >(self, "get_params");
    }

    EMAN::Dict default_get_params() const {
        return EMAN::Filter::get_params();
    }

    void set_params(const EMAN::Dict & p0) {
        call_method< void >(self, "set_params", p0);
    }

    void default_set_params(const EMAN::Dict & p0) {
        EMAN::Filter::set_params(p0);
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(self, "get_param_types");
    }

    std::basic_string<char,std::char_traits<char>,std::allocator<char> > get_name() const {
        return call_method< std::basic_string<char,std::char_traits<char>,std::allocator<char> > >(self, "get_name");
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



    
    class_< EMAN::Aligner, boost::noncopyable, EMAN_Aligner_Wrapper >("Aligner", init<  >())
        .def("get_params", &EMAN::Aligner::get_params, &EMAN_Aligner_Wrapper::default_get_params)
        .def("set_params", &EMAN::Aligner::set_params, &EMAN_Aligner_Wrapper::default_set_params)
    ;

    scope* EMAN_EMObject_scope = new scope(
    class_< EMAN::EMObject >("EMObject", init<  >())
        .def(init< const EMAN::EMObject & >())
        .def(init< int >())
        .def(init< float >())
        .def(init< double >())
        .def(init< std::basic_string<char,std::char_traits<char>,std::allocator<char> > >())
        .def(init< EMAN::EMData * >())
        .def(init< const std::vector<float,std::allocator<float> > & >())
        .def("get_int", &EMAN::EMObject::get_int)
        .def("get_float", &EMAN::EMObject::get_float)
        .def("get_double", &EMAN::EMObject::get_double)
        .def("get_string", &EMAN::EMObject::get_string)
        .def("get_EMData", &EMAN::EMObject::get_EMData, return_internal_reference< 1 >())
        .def("get_farray", &EMAN::EMObject::get_farray)
        .def("is_null", &EMAN::EMObject::is_null)
        .def("to_str", &EMAN::EMObject::to_str)
        .def("get_object_type_name", &EMAN::EMObject::get_object_type_name)
        .staticmethod("get_object_type_name")
    );

    enum_< EMAN::EMObject::ObjectType >("ObjectType")
        .value("FLOATARRAY", EMAN::EMObject::FLOATARRAY)
        .value("EMDATA", EMAN::EMObject::EMDATA)
        .value("STRING", EMAN::EMObject::STRING)
        .value("UNKNOWN", EMAN::EMObject::UNKNOWN)
        .value("INT", EMAN::EMObject::INT)
        .value("DOUBLE", EMAN::EMObject::DOUBLE)
        .value("FLOAT", EMAN::EMObject::FLOAT)
    ;

    delete EMAN_EMObject_scope;

    scope* EMAN_EMUtil_scope = new scope(
    class_< EMAN::EMUtil >("EMUtil", init<  >())
        .def(init< const EMAN::EMUtil & >())
        .def("get_image_type", &EMAN::EMUtil::get_image_type)
        .staticmethod("get_image_type")
        .def("get_image_count", &EMAN::EMUtil::get_image_count)
        .staticmethod("get_image_count")
        .def("get_imageio", &EMAN::EMUtil::get_imageio, return_internal_reference< 1 >(), EMAN_EMUtil_get_imageio_overloads_2_3())
        .staticmethod("get_imageio")
        .def("get_imagetype_name", &EMAN::EMUtil::get_imagetype_name)
        .staticmethod("get_imagetype_name")
        .def("get_datatype_string", &EMAN::EMUtil::get_datatype_string)
        .staticmethod("get_datatype_string")
        .def("dump_dict", &EMAN::EMUtil::dump_dict)
        .staticmethod("dump_dict")
        .def("is_same_size", &EMAN::EMUtil::is_same_size)
        .staticmethod("is_same_size")
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
        .value("IMAGE_PIF", EMAN::EMUtil::IMAGE_PIF)
        .value("IMAGE_AMIRA", EMAN::EMUtil::IMAGE_AMIRA)
        .value("IMAGE_LST", EMAN::EMUtil::IMAGE_LST)
        .value("IMAGE_DM3", EMAN::EMUtil::IMAGE_DM3)
        .value("IMAGE_TIFF", EMAN::EMUtil::IMAGE_TIFF)
        .value("IMAGE_SAL", EMAN::EMUtil::IMAGE_SAL)
        .value("IMAGE_IMAGIC", EMAN::EMUtil::IMAGE_IMAGIC)
        .value("IMAGE_HDF", EMAN::EMUtil::IMAGE_HDF)
        .value("IMAGE_SPIDER", EMAN::EMUtil::IMAGE_SPIDER)
        .value("IMAGE_EMIM", EMAN::EMUtil::IMAGE_EMIM)
        .value("IMAGE_SINGLE_SPIDER", EMAN::EMUtil::IMAGE_SINGLE_SPIDER)
        .value("IMAGE_PNG", EMAN::EMUtil::IMAGE_PNG)
        .value("IMAGE_PGM", EMAN::EMUtil::IMAGE_PGM)
        .value("IMAGE_EM", EMAN::EMUtil::IMAGE_EM)
        .value("IMAGE_VTK", EMAN::EMUtil::IMAGE_VTK)
        .value("IMAGE_UNKNOWN", EMAN::EMUtil::IMAGE_UNKNOWN)
    ;

    delete EMAN_EMUtil_scope;

    scope* EMAN_EMData_scope = new scope(
    class_< EMAN::EMData >("EMData", init<  >())
        .def(init< const EMAN::EMData & >())
        .def("copy", &EMAN::EMData::copy, return_value_policy< manage_new_object >(), EMAN_EMData_copy_overloads_0_2())
        .def("get_clip", &EMAN::EMData::get_clip, return_value_policy< manage_new_object >())
        .def("insert_clip", &EMAN::EMData::insert_clip)
        .def("read_image", &EMAN::EMData::read_image, EMAN_EMData_read_image_overloads_1_5())
        .def("write_image", &EMAN::EMData::write_image, EMAN_EMData_write_image_overloads_1_4())
        .def("filter", &EMAN::EMData::filter)
        .def("cmp", &EMAN::EMData::cmp)
        .def("align", &EMAN::EMData::align, return_value_policy< manage_new_object >(), EMAN_EMData_align_overloads_2_3())
        .def("normalize", &EMAN::EMData::normalize)
        .def("is_complex", &EMAN::EMData::is_complex)
        .def("set_complex", &EMAN::EMData::set_complex)
        .def("set_ri", &EMAN::EMData::set_ri)
        .def("ri2ap", &EMAN::EMData::ri2ap)
        .def("ap2ri", &EMAN::EMData::ap2ri)
        .def("get_parent", &EMAN::EMData::get_parent, return_value_policy< manage_new_object >())
        .def("set_parent", &EMAN::EMData::set_parent)
        .def("setup4slice", &EMAN::EMData::setup4slice, EMAN_EMData_setup4slice_overloads_0_1())
        .def("to_corner", &EMAN::EMData::to_corner)
        .def("do_fft", &EMAN::EMData::do_fft, return_value_policy< manage_new_object >())
        .def("do_ift", &EMAN::EMData::do_ift, return_value_policy< manage_new_object >())
        .def("gimme_fft", &EMAN::EMData::gimme_fft)
        .def("rotate_x", &EMAN::EMData::rotate_x)
        .def("rotate_translate", &EMAN::EMData::rotate_translate, EMAN_EMData_rotate_translate_overloads_0_5())
        .def("fast_rotate_translate", &EMAN::EMData::fast_rotate_translate, EMAN_EMData_fast_rotate_translate_overloads_0_1())
        .def("dot_rotate_translate", &EMAN::EMData::dot_rotate_translate)
        .def("fast_translate", &EMAN::EMData::fast_translate, EMAN_EMData_fast_translate_overloads_0_1())
        .def("rotate_180", &EMAN::EMData::rotate_180)
        .def("do_radon", &EMAN::EMData::do_radon, return_value_policy< manage_new_object >())
        .def("vertical_acf", &EMAN::EMData::vertical_acf, return_value_policy< manage_new_object >())
        .def("vertical_flip", &EMAN::EMData::vertical_flip)
        .def("horizontal_flip", &EMAN::EMData::horizontal_flip)
        .def("set_flipped", &EMAN::EMData::set_flipped)
        .def("dot", &EMAN::EMData::dot, EMAN_EMData_dot_overloads_1_2())
        .def("calc_ccf", &EMAN::EMData::calc_ccf, return_value_policy< manage_new_object >(), EMAN_EMData_calc_ccf_overloads_1_3())
        .def("make_rotational_footprint", &EMAN::EMData::make_rotational_footprint, return_value_policy< manage_new_object >(), EMAN_EMData_make_rotational_footprint_overloads_0_2())
        .def("calc_ccfx", &EMAN::EMData::calc_ccfx, return_value_policy< manage_new_object >(), EMAN_EMData_calc_ccfx_overloads_1_4())
        .def("unwrap", &EMAN::EMData::unwrap, return_value_policy< manage_new_object >(), EMAN_EMData_unwrap_overloads_0_6())
        .def("calc_fourier_shell_correlation", &EMAN::EMData::calc_fourier_shell_correlation)
        .def("mean_shrink", &EMAN::EMData::mean_shrink)
        .def("median_shrink", &EMAN::EMData::median_shrink)
        .def("apply_radial_func", &EMAN::EMData::apply_radial_func)
        .def("set_talign_params", (void (EMAN::EMData::*)(float, float) )&EMAN::EMData::set_talign_params)
        .def("set_talign_params", (void (EMAN::EMData::*)(float, float, float) )&EMAN::EMData::set_talign_params)
        .def("set_ralign_params", (void (EMAN::EMData::*)(float, float, float) )&EMAN::EMData::set_ralign_params)
        .def("set_ralign_params", (void (EMAN::EMData::*)(const EMAN::Rotation &) )&EMAN::EMData::set_ralign_params)
        .def("set_align_score", &EMAN::EMData::set_align_score)
        .def("get_align_score", &EMAN::EMData::get_align_score)
        .def("add", (int (EMAN::EMData::*)(float) )&EMAN::EMData::add)
        .def("add", (int (EMAN::EMData::*)(const EMAN::EMData &) )&EMAN::EMData::add)
        .def("sub", &EMAN::EMData::sub)
        .def("mult", (int (EMAN::EMData::*)(float) )&EMAN::EMData::mult)
        .def("mult", (int (EMAN::EMData::*)(const EMAN::EMData &) )&EMAN::EMData::mult)
        .def("div", (int (EMAN::EMData::*)(float) )&EMAN::EMData::div)
        .def("div", (int (EMAN::EMData::*)(const EMAN::EMData &) )&EMAN::EMData::div)
        .def("done_data", &EMAN::EMData::done_data)
        .def("update", &EMAN::EMData::update)
        .def("get_ctf", &EMAN::EMData::get_ctf, return_internal_reference< 1 >())
        .def("set_ctf", &EMAN::EMData::set_ctf)
        .def("set_size", &EMAN::EMData::set_size)
        .def("get_min_location", &EMAN::EMData::get_min_location)
        .def("get_max_location", &EMAN::EMData::get_max_location)
        .def("get_max_index", &EMAN::EMData::get_max_index)
        .def("get_attr_dict", &EMAN::EMData::get_attr_dict)
        .def("get_value_at", (float (EMAN::EMData::*)(int, int, int) const)&EMAN::EMData::get_value_at)
        .def("get_value_at", (float (EMAN::EMData::*)(int, int) const)&EMAN::EMData::get_value_at)
        .def("sget_value_at", (float (EMAN::EMData::*)(int, int, int) const)&EMAN::EMData::sget_value_at)
        .def("sget_value_at", (float (EMAN::EMData::*)(int, int) const)&EMAN::EMData::sget_value_at)
        .def("get_value_at_interp", (float (EMAN::EMData::*)(float, float) const)&EMAN::EMData::get_value_at_interp)
        .def("get_value_at_interp", (float (EMAN::EMData::*)(float, float, float) const)&EMAN::EMData::get_value_at_interp)
        .def("set_value_at", (void (EMAN::EMData::*)(int, int, int, float) )&EMAN::EMData::set_value_at)
        .def("set_value_at", (void (EMAN::EMData::*)(int, int, float) )&EMAN::EMData::set_value_at)
        .def("get_x", &EMAN::EMData::get_x)
        .def("get_y", &EMAN::EMData::get_y)
        .def("get_z", &EMAN::EMData::get_z)
        .def("get_translation", &EMAN::EMData::get_translation)
        .def("set_translation", &EMAN::EMData::set_translation)
        .def("get_rotation", &EMAN::EMData::get_rotation)
        .def("get_trans_align", &EMAN::EMData::get_trans_align)
        .def("dump_data", &EMAN::EMData::dump_data)
        .def( self / self )
        .def( self * self )
        .def( self - self )
        .def( self + self )
        .def( other< float >() / self )
        .def( other< float >() * self )
        .def( other< float >() - self )
        .def( self - other< float >() )
        .def( self + other< float >() )
        .def( other< float >() + self )
        .def( self / other< float >() )
        .def( self * other< float >() )
        .def( self += other< float >() )
        .def( self /= other< float >() )
        .def( self += self )
        .def( self /= self )
    );
    EMAN_EMData_scope->attr("HEADER_ONLY") = EMAN::EMData::HEADER_ONLY;
    EMAN_EMData_scope->attr("HEADER_AND_DATA") = EMAN::EMData::HEADER_AND_DATA;
    EMAN_EMData_scope->attr("IS_3D") = EMAN::EMData::IS_3D;
    EMAN_EMData_scope->attr("NOT_3D") = EMAN::EMData::NOT_3D;
    EMAN_EMData_scope->attr("DATA_READ_ONLY") = EMAN::EMData::DATA_READ_ONLY;
    EMAN_EMData_scope->attr("DATA_READ_WRITE") = EMAN::EMData::DATA_READ_WRITE;
    delete EMAN_EMData_scope;

    class_< EMAN::Projector, boost::noncopyable, EMAN_Projector_Wrapper >("Projector", init<  >())
        .def("set_params", &EMAN::Projector::set_params)
        .def("get_params", &EMAN::Projector::get_params, &EMAN_Projector_Wrapper::default_get_params)
    ;

    class_< EMAN::AlignerFactory, boost::noncopyable >("AlignerFactory", no_init)
        .def("instance", &EMAN::AlignerFactory::instance, return_value_policy< reference_existing_object >())
        .staticmethod("instance")
        .def("get", (EMAN::Aligner * (EMAN::AlignerFactory::*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >) )&EMAN::AlignerFactory::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Aligner * (EMAN::AlignerFactory::*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict &) )&EMAN::AlignerFactory::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::AlignerFactory::get_list)
    ;

    class_< EMAN::Filter, boost::noncopyable, EMAN_Filter_Wrapper >("Filter", init<  >())
        .def("process", &EMAN::Filter::process, &EMAN_Filter_Wrapper::default_process)
        .def("process_list", &EMAN::Filter::process_list, &EMAN_Filter_Wrapper::default_process_list)
        .def("get_params", &EMAN::Filter::get_params, &EMAN_Filter_Wrapper::default_get_params)
        .def("set_params", &EMAN::Filter::set_params, &EMAN_Filter_Wrapper::default_set_params)
    ;

    class_< EMAN::FilterFactory, boost::noncopyable >("FilterFactory", no_init)
        .def("instance", &EMAN::FilterFactory::instance, return_value_policy< reference_existing_object >())
        .staticmethod("instance")
        .def("get", (EMAN::Filter * (EMAN::FilterFactory::*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >) )&EMAN::FilterFactory::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Filter * (EMAN::FilterFactory::*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict &) )&EMAN::FilterFactory::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::FilterFactory::get_list)
    ;

    class_< EMAN::ProjectorFactory, boost::noncopyable >("ProjectorFactory", no_init)
        .def("instance", &EMAN::ProjectorFactory::instance, return_value_policy< reference_existing_object >())
        .staticmethod("instance")
        .def("get", (EMAN::Projector * (EMAN::ProjectorFactory::*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >) )&EMAN::ProjectorFactory::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Projector * (EMAN::ProjectorFactory::*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict &) )&EMAN::ProjectorFactory::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::ProjectorFactory::get_list)
    ;

    scope* EMAN_Log_scope = new scope(
    class_< EMAN::Log, boost::noncopyable >("Log", no_init)
        .def("logger", &EMAN::Log::logger, return_value_policy< reference_existing_object >())
        .staticmethod("logger")
        .def("set_level", &EMAN::Log::set_level)
        .def("set_logfile", &EMAN::Log::set_logfile)
    );

    enum_< EMAN::Log::LogLevel >("LogLevel")
        .value("ERROR_LOG", EMAN::Log::ERROR_LOG)
        .value("NORMAL_LOG", EMAN::Log::NORMAL_LOG)
        .value("WARNING_LOG", EMAN::Log::WARNING_LOG)
        .value("VARIABLE_LOG", EMAN::Log::VARIABLE_LOG)
    ;

    delete EMAN_Log_scope;

}
