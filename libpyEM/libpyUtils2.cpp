
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <ctf.h>
#include <emdata.h>
#include <emutil.h>
#include <imageio.h>
#include <log.h>
#include <pylist.h>
#include <util.h>
#include <xydata.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_Log_end_overloads_1_3, end, 1, 3)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_EMUtil_get_imageio_overloads_2_3, EMAN::EMUtil::get_imageio, 2, 3)

struct EMAN_Ctf_Wrapper: EMAN::Ctf
{
    EMAN_Ctf_Wrapper(PyObject* self_, const EMAN::Ctf& p0):
        EMAN::Ctf(p0), self(self_) {}

    EMAN_Ctf_Wrapper(PyObject* self_):
        EMAN::Ctf(), self(self_) {}

    int from_string(std::string p0) {
        return call_method< int >(self, "from_string", p0);
    }

    std::string to_string() const {
        return call_method< std::string >(self, "to_string");
    }

    void from_dict(const EMAN::Dict& p0) {
        call_method< void >(self, "from_dict", p0);
    }

    EMAN::Dict to_dict() const {
        return call_method< EMAN::Dict >(self, "to_dict");
    }

    std::vector<float,std::allocator<float> > compute_1d(EMAN::EMData* p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        return call_method< std::vector<float,std::allocator<float> > >(self, "compute_1d", p0, p1, p2);
    }

    void compute_2d_real(EMAN::EMData* p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        call_method< void >(self, "compute_2d_real", p0, p1, p2);
    }

    void compute_2d_complex(EMAN::EMData* p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        call_method< void >(self, "compute_2d_complex", p0, p1, p2);
    }

    void copy_from(const EMAN::Ctf* p0) {
        call_method< void >(self, "copy_from", p0);
    }

    bool equal(const EMAN::Ctf* p0) const {
        return call_method< bool >(self, "equal", p0);
    }

    float get_defocus() const {
        return call_method< float >(self, "get_defocus");
    }

    float get_bfactor() const {
        return call_method< float >(self, "get_bfactor");
    }

    PyObject* self;
};
// Unique type for unnamed enums
#ifndef PYSTE_UNIQUE_INT_DEFINED
#define PYSTE_UNIQUE_INT_DEFINED
template<int num>
struct UniqueInt {
   int v;
   enum { value=num };
   UniqueInt(int v_):
       v(v_)
   {}
   operator int() const
   { return v; }
};
#endif // PYSTE_UNIQUE_INT_DEFINED 

struct EMAN_SimpleCtf_Wrapper: EMAN::SimpleCtf
{
    EMAN_SimpleCtf_Wrapper(PyObject* self_, const EMAN::SimpleCtf& p0):
        EMAN::SimpleCtf(p0), self(self_) {}

    EMAN_SimpleCtf_Wrapper(PyObject* self_):
        EMAN::SimpleCtf(), self(self_) {}

    std::vector<float,std::allocator<float> > compute_1d(EMAN::EMData* p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        return call_method< std::vector<float,std::allocator<float> > >(self, "compute_1d", p0, p1, p2);
    }

    std::vector<float,std::allocator<float> > default_compute_1d_2(EMAN::EMData* p0, EMAN::Ctf::CtfType p1) {
        return EMAN::SimpleCtf::compute_1d(p0, p1);
    }

    std::vector<float,std::allocator<float> > default_compute_1d_3(EMAN::EMData* p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        return EMAN::SimpleCtf::compute_1d(p0, p1, p2);
    }

    void compute_2d_real(EMAN::EMData* p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        call_method< void >(self, "compute_2d_real", p0, p1, p2);
    }

    void default_compute_2d_real_2(EMAN::EMData* p0, EMAN::Ctf::CtfType p1) {
        EMAN::SimpleCtf::compute_2d_real(p0, p1);
    }

    void default_compute_2d_real_3(EMAN::EMData* p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        EMAN::SimpleCtf::compute_2d_real(p0, p1, p2);
    }

    void compute_2d_complex(EMAN::EMData* p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        call_method< void >(self, "compute_2d_complex", p0, p1, p2);
    }

    void default_compute_2d_complex_2(EMAN::EMData* p0, EMAN::Ctf::CtfType p1) {
        EMAN::SimpleCtf::compute_2d_complex(p0, p1);
    }

    void default_compute_2d_complex_3(EMAN::EMData* p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        EMAN::SimpleCtf::compute_2d_complex(p0, p1, p2);
    }

    int from_string(std::string p0) {
        return call_method< int >(self, "from_string", p0);
    }

    int default_from_string(std::string p0) {
        return EMAN::SimpleCtf::from_string(p0);
    }

    std::string to_string() const {
        return call_method< std::string >(self, "to_string");
    }

    std::string default_to_string() const {
        return EMAN::SimpleCtf::to_string();
    }

    void from_dict(const EMAN::Dict& p0) {
        call_method< void >(self, "from_dict", p0);
    }

    void default_from_dict(const EMAN::Dict& p0) {
        EMAN::SimpleCtf::from_dict(p0);
    }

    EMAN::Dict to_dict() const {
        return call_method< EMAN::Dict >(self, "to_dict");
    }

    EMAN::Dict default_to_dict() const {
        return EMAN::SimpleCtf::to_dict();
    }

    void copy_from(const EMAN::Ctf* p0) {
        call_method< void >(self, "copy_from", p0);
    }

    void default_copy_from(const EMAN::Ctf* p0) {
        EMAN::SimpleCtf::copy_from(p0);
    }

    bool equal(const EMAN::Ctf* p0) const {
        return call_method< bool >(self, "equal", p0);
    }

    bool default_equal(const EMAN::Ctf* p0) const {
        return EMAN::SimpleCtf::equal(p0);
    }

    float get_defocus() const {
        return call_method< float >(self, "get_defocus");
    }

    float default_get_defocus() const {
        return EMAN::SimpleCtf::get_defocus();
    }

    float get_bfactor() const {
        return call_method< float >(self, "get_bfactor");
    }

    float default_get_bfactor() const {
        return EMAN::SimpleCtf::get_bfactor();
    }

    PyObject* self;
};

struct EMAN_ImageIO_Wrapper: EMAN::ImageIO
{
    EMAN_ImageIO_Wrapper(PyObject* self_, const EMAN::ImageIO& p0):
        EMAN::ImageIO(p0), self(self_) {}

    EMAN_ImageIO_Wrapper(PyObject* self_):
        EMAN::ImageIO(), self(self_) {}

    int read_header(EMAN::Dict& p0, int p1, const EMAN::Region* p2, bool p3) {
        return call_method< int >(self, "read_header", p0, p1, p2, p3);
    }

    int write_header(const EMAN::Dict& p0, int p1) {
        return call_method< int >(self, "write_header", p0, p1);
    }

    int read_data(float* p0, int p1, const EMAN::Region* p2, bool p3) {
        return call_method< int >(self, "read_data", p0, p1, p2, p3);
    }

    int write_data(float* p0, int p1) {
        return call_method< int >(self, "write_data", p0, p1);
    }

    int read_ctf(EMAN::Ctf& p0, int p1) {
        return call_method< int >(self, "read_ctf", p0, p1);
    }

    int default_read_ctf_1(EMAN::Ctf& p0) {
        return EMAN::ImageIO::read_ctf(p0);
    }

    int default_read_ctf_2(EMAN::Ctf& p0, int p1) {
        return EMAN::ImageIO::read_ctf(p0, p1);
    }

    int write_ctf(const EMAN::Ctf& p0, int p1) {
        return call_method< int >(self, "write_ctf", p0, p1);
    }

    int default_write_ctf_1(const EMAN::Ctf& p0) {
        return EMAN::ImageIO::write_ctf(p0);
    }

    int default_write_ctf_2(const EMAN::Ctf& p0, int p1) {
        return EMAN::ImageIO::write_ctf(p0, p1);
    }

    int get_nimg() {
        return call_method< int >(self, "get_nimg");
    }

    bool is_complex_mode() {
        return call_method< bool >(self, "is_complex_mode");
    }

    bool is_image_big_endian() {
        return call_method< bool >(self, "is_image_big_endian");
    }

    int init() {
        return call_method< int >(self, "init");
    }

    PyObject* self;
};

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_find_max_overloads_3_4, EMAN::Util::find_max, 3, 4)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_find_min_and_max_overloads_4_6, EMAN::Util::find_min_and_max, 4, 6)


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyUtils2)
{
    scope* EMAN_Log_scope = new scope(
    class_< EMAN::Log, boost::noncopyable >("Log", no_init)
        .def("logger", &EMAN::Log::logger, return_value_policy< reference_existing_object >())
        .def("begin", &EMAN::Log::begin)
        .def("end", &EMAN::Log::end, EMAN_Log_end_overloads_1_3())
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

    scope* EMAN_EMUtil_scope = new scope(
    class_< EMAN::EMUtil >("EMUtil", init<  >())
        .def(init< const EMAN::EMUtil& >())
        .def("vertical_acf", &EMAN::EMUtil::vertical_acf, return_value_policy< manage_new_object >())
        .def("make_image_median", &EMAN::EMUtil::make_image_median, return_value_policy< manage_new_object >())
        .def("get_image_ext_type", &EMAN::EMUtil::get_image_ext_type)
        .def("get_image_type", &EMAN::EMUtil::get_image_type)
        .def("get_image_count", &EMAN::EMUtil::get_image_count)
        .def("get_imageio", &EMAN::EMUtil::get_imageio, return_internal_reference< 1 >(), EMAN_EMUtil_get_imageio_overloads_2_3())
        .def("get_imagetype_name", &EMAN::EMUtil::get_imagetype_name)
        .def("get_datatype_string", &EMAN::EMUtil::get_datatype_string)
        .def("dump_dict", &EMAN::EMUtil::dump_dict)
        .def("is_same_size", &EMAN::EMUtil::is_same_size)
        .def("is_same_ctf", &EMAN::EMUtil::is_same_ctf)
        .def("dump_filters", &EMAN::EMUtil::dump_filters)
        .def("dump_aligners", &EMAN::EMUtil::dump_aligners)
        .def("dump_cmps", &EMAN::EMUtil::dump_cmps)
        .def("dump_averagers", &EMAN::EMUtil::dump_averagers)
        .def("dump_projectors", &EMAN::EMUtil::dump_projectors)
        .def("dump_reconstructors", &EMAN::EMUtil::dump_reconstructors)
        .staticmethod("vertical_acf")
        .staticmethod("get_datatype_string")
        .staticmethod("dump_dict")
        .staticmethod("dump_cmps")
        .staticmethod("get_imageio")
        .staticmethod("get_image_count")
        .staticmethod("get_imagetype_name")
        .staticmethod("get_image_type")
        .staticmethod("is_same_size")
        .staticmethod("dump_reconstructors")
        .staticmethod("make_image_median")
        .staticmethod("dump_averagers")
        .staticmethod("is_same_ctf")
        .staticmethod("dump_aligners")
        .staticmethod("dump_projectors")
        .staticmethod("dump_filters")
        .staticmethod("get_image_ext_type")
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

    scope* EMAN_Ctf_scope = new scope(
    class_< EMAN::Ctf, boost::noncopyable, EMAN_Ctf_Wrapper >("Ctf", init<  >())
        .def("from_string", pure_virtual(&EMAN::Ctf::from_string))
        .def("to_string", pure_virtual(&EMAN::Ctf::to_string))
        .def("from_dict", pure_virtual(&EMAN::Ctf::from_dict))
        .def("to_dict", pure_virtual(&EMAN::Ctf::to_dict))
        .def("compute_1d", pure_virtual(&EMAN::Ctf::compute_1d))
        .def("compute_2d_real", pure_virtual(&EMAN::Ctf::compute_2d_real))
        .def("compute_2d_complex", pure_virtual(&EMAN::Ctf::compute_2d_complex))
        .def("copy_from", pure_virtual(&EMAN::Ctf::copy_from))
        .def("equal", pure_virtual(&EMAN::Ctf::equal))
        .def("get_defocus", pure_virtual(&EMAN::Ctf::get_defocus))
        .def("get_bfactor", pure_virtual(&EMAN::Ctf::get_bfactor))
    );

    enum_< EMAN::Ctf::CtfType >("CtfType")
        .value("CTF_NOISE_S", EMAN::Ctf::CTF_NOISE_S)
        .value("CTF_WIENER_CTF_CORRECTION2", EMAN::Ctf::CTF_WIENER_CTF_CORRECTION2)
        .value("CTF_WIENER_CTF_CORRECTION1", EMAN::Ctf::CTF_WIENER_CTF_CORRECTION1)
        .value("CTF_ABS_AMP_S", EMAN::Ctf::CTF_ABS_AMP_S)
        .value("CTF_SNR_SIGN", EMAN::Ctf::CTF_SNR_SIGN)
        .value("CTF_RELATIVE_SNR", EMAN::Ctf::CTF_RELATIVE_SNR)
        .value("CTF_AMP_S", EMAN::Ctf::CTF_AMP_S)
        .value("CTF_WIENER_CTF_CORRECTION", EMAN::Ctf::CTF_WIENER_CTF_CORRECTION)
        .value("CTF_SIGN", EMAN::Ctf::CTF_SIGN)
        .value("CTF_SNR_WIENER", EMAN::Ctf::CTF_SNR_WIENER)
        .value("CTF_AMP", EMAN::Ctf::CTF_AMP)
        .value("CTF_SNR", EMAN::Ctf::CTF_SNR)
        .value("CTF_CTF_NO_BDECAY", EMAN::Ctf::CTF_CTF_NO_BDECAY)
        .value("CTF_ABS_SNR", EMAN::Ctf::CTF_ABS_SNR)
        .value("CTF_WIENER_FILTER", EMAN::Ctf::CTF_WIENER_FILTER)
        .value("CTF_AMP_NO_BDECAY", EMAN::Ctf::CTF_AMP_NO_BDECAY)
        .value("CTF_NOISE", EMAN::Ctf::CTF_NOISE)
        .value("CTF_BFACTOR", EMAN::Ctf::CTF_BFACTOR)
        .value("CTF_BACKGROUND", EMAN::Ctf::CTF_BACKGROUND)
        .value("CTF_TOTAL_CURVE", EMAN::Ctf::CTF_TOTAL_CURVE)
    ;


    enum_< UniqueInt<0> >("unnamed")
        .value("CTFOS", EMAN::Ctf::CTFOS)
        .export_values()
    ;

    delete EMAN_Ctf_scope;

    class_< EMAN::SimpleCtf, bases< EMAN::Ctf > , EMAN_SimpleCtf_Wrapper >("SimpleCtf", init<  >())
        .def(init< const EMAN::SimpleCtf& >())
        .def_readwrite("defocus", &EMAN::SimpleCtf::defocus)
        .def_readwrite("bfactor", &EMAN::SimpleCtf::bfactor)
        .def_readwrite("amplitude", &EMAN::SimpleCtf::amplitude)
        .def_readwrite("ampcont", &EMAN::SimpleCtf::ampcont)
        .def_readwrite("noise1", &EMAN::SimpleCtf::noise1)
        .def_readwrite("noise2", &EMAN::SimpleCtf::noise2)
        .def_readwrite("noise3", &EMAN::SimpleCtf::noise3)
        .def_readwrite("noise4", &EMAN::SimpleCtf::noise4)
        .def_readwrite("voltage", &EMAN::SimpleCtf::voltage)
        .def_readwrite("cs", &EMAN::SimpleCtf::cs)
        .def_readwrite("apix", &EMAN::SimpleCtf::apix)
        .def("compute_1d", (std::vector<float,std::allocator<float> > (EMAN::SimpleCtf::*)(EMAN::EMData*, EMAN::Ctf::CtfType, EMAN::XYData*) )&EMAN::SimpleCtf::compute_1d, (std::vector<float,std::allocator<float> > (EMAN_SimpleCtf_Wrapper::*)(EMAN::EMData*, EMAN::Ctf::CtfType, EMAN::XYData*))&EMAN_SimpleCtf_Wrapper::default_compute_1d_3)
        .def("compute_1d", (std::vector<float,std::allocator<float> > (EMAN_SimpleCtf_Wrapper::*)(EMAN::EMData*, EMAN::Ctf::CtfType))&EMAN_SimpleCtf_Wrapper::default_compute_1d_2)
        .def("compute_2d_real", (void (EMAN::SimpleCtf::*)(EMAN::EMData*, EMAN::Ctf::CtfType, EMAN::XYData*) )&EMAN::SimpleCtf::compute_2d_real, (void (EMAN_SimpleCtf_Wrapper::*)(EMAN::EMData*, EMAN::Ctf::CtfType, EMAN::XYData*))&EMAN_SimpleCtf_Wrapper::default_compute_2d_real_3)
        .def("compute_2d_real", (void (EMAN_SimpleCtf_Wrapper::*)(EMAN::EMData*, EMAN::Ctf::CtfType))&EMAN_SimpleCtf_Wrapper::default_compute_2d_real_2)
        .def("compute_2d_complex", (void (EMAN::SimpleCtf::*)(EMAN::EMData*, EMAN::Ctf::CtfType, EMAN::XYData*) )&EMAN::SimpleCtf::compute_2d_complex, (void (EMAN_SimpleCtf_Wrapper::*)(EMAN::EMData*, EMAN::Ctf::CtfType, EMAN::XYData*))&EMAN_SimpleCtf_Wrapper::default_compute_2d_complex_3)
        .def("compute_2d_complex", (void (EMAN_SimpleCtf_Wrapper::*)(EMAN::EMData*, EMAN::Ctf::CtfType))&EMAN_SimpleCtf_Wrapper::default_compute_2d_complex_2)
        .def("from_string", (int (EMAN::SimpleCtf::*)(std::string) )&EMAN::SimpleCtf::from_string, (int (EMAN_SimpleCtf_Wrapper::*)(std::string))&EMAN_SimpleCtf_Wrapper::default_from_string)
        .def("to_string", (std::string (EMAN::SimpleCtf::*)() const)&EMAN::SimpleCtf::to_string, (std::string (EMAN_SimpleCtf_Wrapper::*)() const)&EMAN_SimpleCtf_Wrapper::default_to_string)
        .def("from_dict", (void (EMAN::SimpleCtf::*)(const EMAN::Dict&) )&EMAN::SimpleCtf::from_dict, (void (EMAN_SimpleCtf_Wrapper::*)(const EMAN::Dict&))&EMAN_SimpleCtf_Wrapper::default_from_dict)
        .def("to_dict", (EMAN::Dict (EMAN::SimpleCtf::*)() const)&EMAN::SimpleCtf::to_dict, (EMAN::Dict (EMAN_SimpleCtf_Wrapper::*)() const)&EMAN_SimpleCtf_Wrapper::default_to_dict)
        .def("copy_from", (void (EMAN::SimpleCtf::*)(const EMAN::Ctf*) )&EMAN::SimpleCtf::copy_from, (void (EMAN_SimpleCtf_Wrapper::*)(const EMAN::Ctf*))&EMAN_SimpleCtf_Wrapper::default_copy_from)
        .def("equal", (bool (EMAN::SimpleCtf::*)(const EMAN::Ctf*) const)&EMAN::SimpleCtf::equal, (bool (EMAN_SimpleCtf_Wrapper::*)(const EMAN::Ctf*) const)&EMAN_SimpleCtf_Wrapper::default_equal)
        .def("get_defocus", (float (EMAN::SimpleCtf::*)() const)&EMAN::SimpleCtf::get_defocus, (float (EMAN_SimpleCtf_Wrapper::*)() const)&EMAN_SimpleCtf_Wrapper::default_get_defocus)
        .def("get_bfactor", (float (EMAN::SimpleCtf::*)() const)&EMAN::SimpleCtf::get_bfactor, (float (EMAN_SimpleCtf_Wrapper::*)() const)&EMAN_SimpleCtf_Wrapper::default_get_bfactor)
    ;

    scope* EMAN_ImageIO_scope = new scope(
    class_< EMAN::ImageIO, boost::noncopyable, EMAN_ImageIO_Wrapper >("ImageIO", init<  >())
        .def("read_header", pure_virtual(&EMAN::ImageIO::read_header))
        .def("write_header", pure_virtual(&EMAN::ImageIO::write_header))
        .def("read_data", pure_virtual(&EMAN::ImageIO::read_data))
        .def("write_data", pure_virtual(&EMAN::ImageIO::write_data))
        .def("read_ctf", &EMAN::ImageIO::read_ctf, &EMAN_ImageIO_Wrapper::default_read_ctf_2)
        .def("read_ctf", &EMAN_ImageIO_Wrapper::default_read_ctf_1)
        .def("write_ctf", &EMAN::ImageIO::write_ctf, &EMAN_ImageIO_Wrapper::default_write_ctf_2)
        .def("write_ctf", &EMAN_ImageIO_Wrapper::default_write_ctf_1)
        .def("get_nimg", pure_virtual(&EMAN::ImageIO::get_nimg))
        .def("is_complex_mode", pure_virtual(&EMAN::ImageIO::is_complex_mode))
        .def("is_image_big_endian", pure_virtual(&EMAN::ImageIO::is_image_big_endian))
    );

    enum_< EMAN::ImageIO::IOMode >("IOMode")
        .value("READ_ONLY", EMAN::ImageIO::READ_ONLY)
        .value("READ_WRITE", EMAN::ImageIO::READ_WRITE)
    ;

    delete EMAN_ImageIO_scope;

    scope* EMAN_XYData_scope = new scope(
    class_< EMAN::XYData >("XYData", init<  >())
        .def(init< const EMAN::XYData& >())
        .def("read_file", &EMAN::XYData::read_file)
        .def("write_file", &EMAN::XYData::write_file)
        .def("calc_correlation", &EMAN::XYData::calc_correlation)
        .def("update", &EMAN::XYData::update)
        .def("get_yatx", &EMAN::XYData::get_yatx)
        .def("get_x", &EMAN::XYData::get_x)
        .def("set_x", &EMAN::XYData::set_x)
        .def("get_y", &EMAN::XYData::get_y)
        .def("set_y", &EMAN::XYData::set_y)
        .def("get_size", &EMAN::XYData::get_size)
        .def("get_miny", &EMAN::XYData::get_miny)
        .def("get_maxy", &EMAN::XYData::get_maxy)
        .def("is_validx", &EMAN::XYData::is_validx)
    );

    class_< EMAN::XYData::Pair >("Pair", init< const EMAN::XYData::Pair& >())
        .def(init< float, float >())
        .def_readwrite("x", &EMAN::XYData::Pair::x)
        .def_readwrite("y", &EMAN::XYData::Pair::y)
        .def( self < self )
    ;

    delete EMAN_XYData_scope;

    class_< EMAN::Util >("Util", init<  >())
        .def(init< const EMAN::Util& >())
        .def("ap2ri", &EMAN::Util::ap2ri)
        .def("flip_complex_phase", &EMAN::Util::flip_complex_phase)
        .def("file_lock_wait", &EMAN::Util::file_lock_wait)
        .def("file_unlock", &EMAN::Util::file_unlock)
        .def("generate_machine_stamp", &EMAN::Util::generate_machine_stamp)
        .def("check_file_by_magic", &EMAN::Util::check_file_by_magic)
        .def("flip_image", &EMAN::Util::flip_image)
        .def("sstrncmp", &EMAN::Util::sstrncmp)
        .def("get_str_float", (bool (*)(const char*, const char*, float*))&EMAN::Util::get_str_float)
        .def("get_str_float", (bool (*)(const char*, const char*, float*, float*))&EMAN::Util::get_str_float)
        .def("get_str_float", (bool (*)(const char*, const char*, int*, float*, float*))&EMAN::Util::get_str_float)
        .def("get_str_int", (bool (*)(const char*, const char*, int*))&EMAN::Util::get_str_int)
        .def("get_str_int", (bool (*)(const char*, const char*, int*, int*))&EMAN::Util::get_str_int)
        .def("get_str_int", (bool (*)(const char*, const char*, int*, int*, int*))&EMAN::Util::get_str_int)
        .def("get_filename_by_ext", &EMAN::Util::get_filename_by_ext)
        .def("calc_least_square_fit", &EMAN::Util::calc_least_square_fit)
        .def("save_data", (void (*)(const std::vector<float,std::allocator<float> >&, const std::vector<float,std::allocator<float> >&, std::string))&EMAN::Util::save_data)
        .def("save_data", (void (*)(float, float, const std::vector<float,std::allocator<float> >&, std::string))&EMAN::Util::save_data)
        .def("save_data", (void (*)(float, float, float*, size_t, std::string))&EMAN::Util::save_data)
        .def("get_frand", &EMAN::Util::get_frand)
        .def("get_gauss_rand", &EMAN::Util::get_gauss_rand)
        .def("round", &EMAN::Util::round)
        .def("bilinear_interpolate", &EMAN::Util::bilinear_interpolate)
        .def("trilinear_interpolate", &EMAN::Util::trilinear_interpolate)
        .def("find_max", &EMAN::Util::find_max, EMAN_Util_find_max_overloads_3_4())
        .def("find_min_and_max", &EMAN::Util::find_min_and_max, EMAN_Util_find_min_and_max_overloads_4_6())
        .def("calc_best_fft_size", &EMAN::Util::calc_best_fft_size)
        .def("square", (int (*)(int))&EMAN::Util::square)
        .def("square", (float (*)(float))&EMAN::Util::square)
        .def("square", (double (*)(double))&EMAN::Util::square)
        .def("square_sum", &EMAN::Util::square_sum)
        .def("hypot3", &EMAN::Util::hypot3)
        .def("fast_floor", &EMAN::Util::fast_floor)
        .def("agauss", &EMAN::Util::agauss)
        .def("min", (float (*)(float, float))&EMAN::Util::min)
        .def("min", (float (*)(float, float, float))&EMAN::Util::min)
        .def("min", (float (*)(float, float, float, float))&EMAN::Util::min)
        .def("max", (float (*)(float, float))&EMAN::Util::max)
        .def("max", (float (*)(float, float, float))&EMAN::Util::max)
        .def("max", (float (*)(float, float, float, float))&EMAN::Util::max)
        .def("angle_sub_2pi", &EMAN::Util::angle_sub_2pi)
        .def("angle_sub_pi", &EMAN::Util::angle_sub_pi)
        .def("goodf", &EMAN::Util::goodf)
        .def("get_time_label", &EMAN::Util::get_time_label)
        .def("set_log_level", &EMAN::Util::set_log_level)
        .def("eman_copysign", &EMAN::Util::eman_copysign)
        .def("eman_erfc", &EMAN::Util::eman_erfc)
        .staticmethod("flip_complex_phase")
        .staticmethod("square")
        .staticmethod("get_str_int")
        .staticmethod("calc_best_fft_size")
        .staticmethod("get_frand")
        .staticmethod("angle_sub_pi")
        .staticmethod("get_time_label")
        .staticmethod("trilinear_interpolate")
        .staticmethod("file_lock_wait")
        .staticmethod("min")
        .staticmethod("save_data")
        .staticmethod("set_log_level")
        .staticmethod("get_str_float")
        .staticmethod("angle_sub_2pi")
        .staticmethod("get_gauss_rand")
        .staticmethod("find_max")
        .staticmethod("max")
        .staticmethod("file_unlock")
        .staticmethod("fast_floor")
        .staticmethod("ap2ri")
        .staticmethod("sstrncmp")
        .staticmethod("check_file_by_magic")
        .staticmethod("generate_machine_stamp")
        .staticmethod("agauss")
        .staticmethod("eman_erfc")
        .staticmethod("bilinear_interpolate")
        .staticmethod("calc_least_square_fit")
        .staticmethod("find_min_and_max")
        .staticmethod("get_filename_by_ext")
        .staticmethod("square_sum")
        .staticmethod("goodf")
        .staticmethod("hypot3")
        .staticmethod("flip_image")
        .staticmethod("eman_copysign")
        .staticmethod("round")
    ;

}

