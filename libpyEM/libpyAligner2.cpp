
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <aligner.h>
#include <ctf.h>
#include <emdata.h>
#include <emobject.h>
#include <imageio.h>
#include <xydata.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

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

    std::vector<float,std::allocator<float> > compute_1d(int p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
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

    std::vector<float,std::allocator<float> > compute_1d(int p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        return call_method< std::vector<float,std::allocator<float> > >(self, "compute_1d", p0, p1, p2);
    }

    std::vector<float,std::allocator<float> > default_compute_1d_2(int p0, EMAN::Ctf::CtfType p1) {
        return EMAN::SimpleCtf::compute_1d(p0, p1);
    }

    std::vector<float,std::allocator<float> > default_compute_1d_3(int p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
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

    int write_header(const EMAN::Dict& p0, int p1, const EMAN::Region* p2, bool p3) {
        return call_method< int >(self, "write_header", p0, p1, p2, p3);
    }

    int read_data(float* p0, int p1, const EMAN::Region* p2, bool p3) {
        return call_method< int >(self, "read_data", p0, p1, p2, p3);
    }

    int write_data(float* p0, int p1, const EMAN::Region* p2, bool p3) {
        return call_method< int >(self, "write_data", p0, p1, p2, p3);
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

    void flush() {
        call_method< void >(self, "flush");
    }

    int get_nimg() {
        return call_method< int >(self, "get_nimg");
    }

    int default_get_nimg() {
        return EMAN::ImageIO::get_nimg();
    }

    bool is_complex_mode() {
        return call_method< bool >(self, "is_complex_mode");
    }

    bool is_image_big_endian() {
        return call_method< bool >(self, "is_image_big_endian");
    }

    bool is_single_image_format() const {
        return call_method< bool >(self, "is_single_image_format");
    }

    bool default_is_single_image_format() const {
        return EMAN::ImageIO::is_single_image_format();
    }

    void init() {
        call_method< void >(self, "init");
    }

    PyObject* self;
};


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyAligner2)
{
    def("dump_aligners", &EMAN::dump_aligners);
    class_< EMAN::Factory<EMAN::Aligner>, boost::noncopyable >("Aligners", no_init)
        .def("get", (EMAN::Aligner* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >))&EMAN::Factory<EMAN::Aligner>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Aligner* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&))&EMAN::Factory<EMAN::Aligner>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Aligner>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

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
        .def("compute_1d", (std::vector<float,std::allocator<float> > (EMAN::SimpleCtf::*)(int, EMAN::Ctf::CtfType, EMAN::XYData*) )&EMAN::SimpleCtf::compute_1d, (std::vector<float,std::allocator<float> > (EMAN_SimpleCtf_Wrapper::*)(int, EMAN::Ctf::CtfType, EMAN::XYData*))&EMAN_SimpleCtf_Wrapper::default_compute_1d_3)
        .def("compute_1d", (std::vector<float,std::allocator<float> > (EMAN_SimpleCtf_Wrapper::*)(int, EMAN::Ctf::CtfType))&EMAN_SimpleCtf_Wrapper::default_compute_1d_2)
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
        .def("flush", pure_virtual(&EMAN::ImageIO::flush))
        .def("get_nimg", &EMAN::ImageIO::get_nimg, &EMAN_ImageIO_Wrapper::default_get_nimg)
        .def("is_complex_mode", pure_virtual(&EMAN::ImageIO::is_complex_mode))
        .def("is_image_big_endian", pure_virtual(&EMAN::ImageIO::is_image_big_endian))
        .def("is_single_image_format", &EMAN::ImageIO::is_single_image_format, &EMAN_ImageIO_Wrapper::default_is_single_image_format)
    );

    enum_< EMAN::ImageIO::IOMode >("IOMode")
        .value("READ_ONLY", EMAN::ImageIO::READ_ONLY)
        .value("READ_WRITE", EMAN::ImageIO::READ_WRITE)
        .value("WRITE_ONLY", EMAN::ImageIO::WRITE_ONLY)
    ;

    delete EMAN_ImageIO_scope;

}

