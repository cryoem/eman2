
// Includes ====================================================================
#include <boost/python.hpp>
#include <emobject.h>
#include <ctf.h>
#include <pyem.h>
#include <filter.h>
#include <pylist.h>
#include <log.h>
#include <imageio.h>
#include <emdata.h>
#include <emutil.h>
#include <transform.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {


BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMUtil_get_imageio_overloads_2_3, get_imageio, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_copy_overloads_0_2, copy, 0, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_read_image_overloads_1_5, read_image, 1, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_write_image_overloads_1_4, write_image, 1, 4)

struct EMAN_Filter_Wrapper: EMAN::Filter
{
    EMAN_Filter_Wrapper(PyObject* self_, const EMAN::Filter & p0):
        EMAN::Filter(p0), self(self_) {}

    EMAN_Filter_Wrapper(PyObject* self_):
        EMAN::Filter(), self(self_) {}

    EMAN_Filter_Wrapper(PyObject* self_, const EMAN::Dict & p0):
        EMAN::Filter(p0), self(self_) {}

    void process(EMAN::EMData * p0) const {
        call_method< void >(self, "process", p0);
    }

    void default_process(EMAN::EMData * p0) const {
        EMAN::Filter::process(p0);
    }

    void process_list(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> > & p0) const {
        call_method< void >(self, "process_list", p0);
    }

    void default_process_list(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> > & p0) const {
        EMAN::Filter::process_list(p0);
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
    EMAN::vector_to_python<std::string>();
 
    EMAN::map_to_python<EMAN::EMObject>();
    EMAN::map_from_python<EMAN::EMObject>();

    EMAN::Dict_to_python();
    EMAN::Dict_from_python();

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

    class_< EMAN::EMObject >("EMObject", init<  >())
        .def(init< const EMAN::EMObject & >())
        .def(init< int >())
        .def(init< float >())
        .def(init< double >())
        .def(init< std::basic_string<char,std::char_traits<char>,std::allocator<char> > >())
        .def(init< EMAN::EMData * >())
        .def("get_int", &EMAN::EMObject::get_int)
        .def("get_float", &EMAN::EMObject::get_float)
        .def("get_double", &EMAN::EMObject::get_double)
        .def("get_string", &EMAN::EMObject::get_string)
        .def("get_EMData", &EMAN::EMObject::get_EMData, return_internal_reference< 1 >())
        .def("is_null", &EMAN::EMObject::is_null)
        .def("to_str", &EMAN::EMObject::to_str)
    ;

    scope* EMAN_EMData_scope = new scope(
    class_< EMAN::EMData >("EMData", init<  >())
        .def(init< const EMAN::EMData & >())
        .def("copy", &EMAN::EMData::copy, return_value_policy< manage_new_object >(), EMAN_EMData_copy_overloads_0_2())
        .def("get_clip", &EMAN::EMData::get_clip, return_value_policy< manage_new_object >())
        .def("insert_clip", &EMAN::EMData::insert_clip)
        .def("read_image", &EMAN::EMData::read_image, EMAN_EMData_read_image_overloads_1_5())
        .def("write_image", &EMAN::EMData::write_image, EMAN_EMData_write_image_overloads_1_4())
        .def("normalize", &EMAN::EMData::normalize)
        .def("is_complex", &EMAN::EMData::is_complex)
        .def("ri2ap", &EMAN::EMData::ri2ap)
        .def("ap2ri", &EMAN::EMData::ap2ri)
        .def("do_fft", &EMAN::EMData::do_fft, return_value_policy< manage_new_object >())
        .def("add", (int (EMAN::EMData::*)(float) )&EMAN::EMData::add)
        .def("add", (int (EMAN::EMData::*)(const EMAN::EMData &) )&EMAN::EMData::add)
        .def("sub", &EMAN::EMData::sub)
        .def("mult", (int (EMAN::EMData::*)(float) )&EMAN::EMData::mult)
        .def("mult", (int (EMAN::EMData::*)(const EMAN::EMData &) )&EMAN::EMData::mult)
        .def("div", (int (EMAN::EMData::*)(float) )&EMAN::EMData::div)
        .def("div", (int (EMAN::EMData::*)(const EMAN::EMData &) )&EMAN::EMData::div)
        .def("done_data", &EMAN::EMData::done_data)
        .def("get_ctf", &EMAN::EMData::get_ctf, return_internal_reference< 1 >())
        .def("set_ctf", &EMAN::EMData::set_ctf)
        .def("set_size", &EMAN::EMData::set_size)
        .def("get_attr_dict", &EMAN::EMData::get_attr_dict)
        .def("get_value_at", (float (EMAN::EMData::*)(int, int, int) const)&EMAN::EMData::get_value_at)
        .def("get_value_at", (float (EMAN::EMData::*)(int, int) const)&EMAN::EMData::get_value_at)
        .def("sget_value_at", (float (EMAN::EMData::*)(int, int, int) const)&EMAN::EMData::sget_value_at)
        .def("sget_value_at", (float (EMAN::EMData::*)(int, int) const)&EMAN::EMData::sget_value_at)
        .def("get_value_at_interp", &EMAN::EMData::get_value_at_interp)
        .def("set_value_at", (void (EMAN::EMData::*)(int, int, int, float) )&EMAN::EMData::set_value_at)
        .def("set_value_at", (void (EMAN::EMData::*)(int, int, float) )&EMAN::EMData::set_value_at)
        .def("get_x", &EMAN::EMData::get_x)
        .def("get_y", &EMAN::EMData::get_y)
        .def("get_z", &EMAN::EMData::get_z)
        .def("dump_data", &EMAN::EMData::dump_data)
        .def( self * self )
        .def( self / self )
        .def( other< float >() * self )
        .def( other< float >() - self )
        .def( self + self )
        .def( other< float >() / self )
        .def( self - other< float >() )
        .def( self * other< float >() )
        .def( self + other< float >() )
        .def( self - self )
        .def( other< float >() + self )
        .def( self / other< float >() )
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

    class_< EMAN::FilterFactory, boost::noncopyable >("FilterFactory", no_init)
        .def("instance", &EMAN::FilterFactory::instance, return_value_policy< reference_existing_object >())
        .staticmethod("instance")
        .def("get", (EMAN::Filter * (EMAN::FilterFactory::*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >) )&EMAN::FilterFactory::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Filter * (EMAN::FilterFactory::*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict &) )&EMAN::FilterFactory::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::FilterFactory::get_list)
    ;

    class_< EMAN::Filter, boost::noncopyable, EMAN_Filter_Wrapper >("Filter", init<  >())
        .def(init< const EMAN::Dict & >())
        .def("get_params", &EMAN::Filter::get_params)
        .def("set_params", &EMAN::Filter::set_params)
        .def("process", &EMAN::Filter::process, &EMAN_Filter_Wrapper::default_process)
        .def("process_list", &EMAN::Filter::process_list, &EMAN_Filter_Wrapper::default_process_list)
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
