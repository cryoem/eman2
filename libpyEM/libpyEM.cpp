
// Includes ====================================================================
#include <boost/python.hpp>
#include <ctf.h>
#include <emobject.h>
#include <pylist.h>
#include <pyem.h>
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

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_read_image_overloads_1_5, read_image, 1, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_write_image_overloads_1_4, write_image, 1, 4)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_EMData_read_images_by_index_overloads_2_3, py_read_images_by_index, 2, 3)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_EMData_read_images_by_ext_overloads_3_5, py_read_images_by_ext, 3, 5)



}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyEM)
{
    EMAN::vector_to_python<EMAN::EMData*>();
    EMAN::vector_from_python<int>();
    
    EMAN::map_to_python<EMAN::EMObject>();
    EMAN::map_from_python<EMAN::EMObject>();
    
    EMAN::Dict_to_python();
    EMAN::Dict_from_python();
    
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

    scope* EMAN_EMData_scope = new scope(
    class_< EMAN::EMData >("EMData", init<  >())
        .def(init< const EMAN::EMData & >())
        .def("read_image", &EMAN::EMData::read_image, EMAN_EMData_read_image_overloads_1_5())
        .def("write_image", &EMAN::EMData::write_image, EMAN_EMData_write_image_overloads_1_4())
        .def("get_attr_dict", &EMAN::EMData::get_attr_dict)
        .def("get_ctf", &EMAN::EMData::get_ctf, return_internal_reference< 1 >())
        .def("set_ctf", &EMAN::EMData::set_ctf)
        .def("get_value_at", (float (EMAN::EMData::*)(int, int, int) )&EMAN::EMData::get_value_at)
        .def("get_value_at", (float (EMAN::EMData::*)(int, int) )&EMAN::EMData::get_value_at)
        .def("sget_value_at", (float (EMAN::EMData::*)(int, int, int) )&EMAN::EMData::sget_value_at)
        .def("sget_value_at", (float (EMAN::EMData::*)(int, int) )&EMAN::EMData::sget_value_at)
        .def("get_value_at_interp", &EMAN::EMData::get_value_at_interp)
        .def("set_value_at", (void (EMAN::EMData::*)(int, int, int, float) )&EMAN::EMData::set_value_at)
        .def("set_value_at", (void (EMAN::EMData::*)(int, int, float) )&EMAN::EMData::set_value_at)
        .def("dump_data", &EMAN::EMData::dump_data)
        .def("read_images_by_index", &py_read_images_by_index, EMAN_EMData_read_images_by_index_overloads_2_3())
        .staticmethod("read_images_by_index")
        .def("read_images_by_ext", &py_read_images_by_ext, EMAN_EMData_read_images_by_ext_overloads_3_5())
        .staticmethod("read_images_by_ext")
    );
    EMAN_EMData_scope->attr("HEADER_ONLY") = EMAN::EMData::HEADER_ONLY;
    EMAN_EMData_scope->attr("HEADER_AND_DATA") = EMAN::EMData::HEADER_AND_DATA;
    EMAN_EMData_scope->attr("IS_3D") = EMAN::EMData::IS_3D;
    EMAN_EMData_scope->attr("NOT_3D") = EMAN::EMData::NOT_3D;
    delete EMAN_EMData_scope;

    scope* EMAN_Log_scope = new scope(
    class_< EMAN::Log, boost::noncopyable >("Log", no_init)
        .def("logger", &EMAN::Log::logger, return_internal_reference< 1 >())
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
