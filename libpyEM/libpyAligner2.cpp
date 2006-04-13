
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <aligner.h>
#include <ctf.h>
#include <emdata.h>
#include <emobject.h>
#include <xydata.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct EMAN_Aligner_Wrapper: EMAN::Aligner
{
    EMAN::EMData* align(EMAN::EMData* p0, EMAN::EMData* p1) const {
        return call_method< EMAN::EMData* >(py_self, "align", p0, p1);
    }

    EMAN::EMData* align(EMAN::EMData* p0, EMAN::EMData* p1, const std::string& p2, const EMAN::Dict& p3) const {
        return call_method< EMAN::EMData* >(py_self, "align", p0, p1, p2, p3);
    }

    std::string get_name() const {
        return call_method< std::string >(py_self, "get_name");
    }

    std::string get_desc() const {
        return call_method< std::string >(py_self, "get_desc");
    }

    EMAN::Dict get_params() const {
        return call_method< EMAN::Dict >(py_self, "get_params");
    }

    EMAN::Dict default_get_params() const {
        return EMAN::Aligner::get_params();
    }

    void set_params(const EMAN::Dict& p0) {
        call_method< void >(py_self, "set_params", p0);
    }

    void default_set_params(const EMAN::Dict& p0) {
        EMAN::Aligner::set_params(p0);
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(py_self, "get_param_types");
    }

    PyObject* py_self;
};

struct EMAN_Ctf_Wrapper: EMAN::Ctf
{
    int from_string(const std::string& p0) {
        return call_method< int >(py_self, "from_string", p0);
    }

    std::string to_string() const {
        return call_method< std::string >(py_self, "to_string");
    }

    void from_dict(const EMAN::Dict& p0) {
        call_method< void >(py_self, "from_dict", p0);
    }

    EMAN::Dict to_dict() const {
        return call_method< EMAN::Dict >(py_self, "to_dict");
    }

    std::vector<float,std::allocator<float> > compute_1d(int p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        return call_method< std::vector<float,std::allocator<float> > >(py_self, "compute_1d", p0, p1, p2);
    }

    void compute_2d_real(EMAN::EMData* p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        call_method< void >(py_self, "compute_2d_real", p0, p1, p2);
    }

    void compute_2d_complex(EMAN::EMData* p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        call_method< void >(py_self, "compute_2d_complex", p0, p1, p2);
    }

    void copy_from(const EMAN::Ctf* p0) {
        call_method< void >(py_self, "copy_from", p0);
    }

    bool equal(const EMAN::Ctf* p0) const {
        return call_method< bool >(py_self, "equal", p0);
    }

    float get_defocus() const {
        return call_method< float >(py_self, "get_defocus");
    }

    float get_bfactor() const {
        return call_method< float >(py_self, "get_bfactor");
    }

    PyObject* py_self;
};

struct EMAN_SimpleCtf_Wrapper: EMAN::SimpleCtf
{
    EMAN_SimpleCtf_Wrapper(PyObject* py_self_):
        EMAN::SimpleCtf(), py_self(py_self_) {}

    std::vector<float,std::allocator<float> > compute_1d(int p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        return call_method< std::vector<float,std::allocator<float> > >(py_self, "compute_1d", p0, p1, p2);
    }

    std::vector<float,std::allocator<float> > default_compute_1d_2(int p0, EMAN::Ctf::CtfType p1) {
        return EMAN::SimpleCtf::compute_1d(p0, p1);
    }

    std::vector<float,std::allocator<float> > default_compute_1d_3(int p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        return EMAN::SimpleCtf::compute_1d(p0, p1, p2);
    }

    void compute_2d_real(EMAN::EMData* p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        call_method< void >(py_self, "compute_2d_real", p0, p1, p2);
    }

    void default_compute_2d_real_2(EMAN::EMData* p0, EMAN::Ctf::CtfType p1) {
        EMAN::SimpleCtf::compute_2d_real(p0, p1);
    }

    void default_compute_2d_real_3(EMAN::EMData* p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        EMAN::SimpleCtf::compute_2d_real(p0, p1, p2);
    }

    void compute_2d_complex(EMAN::EMData* p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        call_method< void >(py_self, "compute_2d_complex", p0, p1, p2);
    }

    void default_compute_2d_complex_2(EMAN::EMData* p0, EMAN::Ctf::CtfType p1) {
        EMAN::SimpleCtf::compute_2d_complex(p0, p1);
    }

    void default_compute_2d_complex_3(EMAN::EMData* p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        EMAN::SimpleCtf::compute_2d_complex(p0, p1, p2);
    }

    int from_string(const std::string& p0) {
        return call_method< int >(py_self, "from_string", p0);
    }

    int default_from_string(const std::string& p0) {
        return EMAN::SimpleCtf::from_string(p0);
    }

    std::string to_string() const {
        return call_method< std::string >(py_self, "to_string");
    }

    std::string default_to_string() const {
        return EMAN::SimpleCtf::to_string();
    }

    void from_dict(const EMAN::Dict& p0) {
        call_method< void >(py_self, "from_dict", p0);
    }

    void default_from_dict(const EMAN::Dict& p0) {
        EMAN::SimpleCtf::from_dict(p0);
    }

    EMAN::Dict to_dict() const {
        return call_method< EMAN::Dict >(py_self, "to_dict");
    }

    EMAN::Dict default_to_dict() const {
        return EMAN::SimpleCtf::to_dict();
    }

    void copy_from(const EMAN::Ctf* p0) {
        call_method< void >(py_self, "copy_from", p0);
    }

    void default_copy_from(const EMAN::Ctf* p0) {
        EMAN::SimpleCtf::copy_from(p0);
    }

    bool equal(const EMAN::Ctf* p0) const {
        return call_method< bool >(py_self, "equal", p0);
    }

    bool default_equal(const EMAN::Ctf* p0) const {
        return EMAN::SimpleCtf::equal(p0);
    }

    float get_defocus() const {
        return call_method< float >(py_self, "get_defocus");
    }

    float default_get_defocus() const {
        return EMAN::SimpleCtf::get_defocus();
    }

    float get_bfactor() const {
        return call_method< float >(py_self, "get_bfactor");
    }

    float default_get_bfactor() const {
        return EMAN::SimpleCtf::get_bfactor();
    }

    PyObject* py_self;
};


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyAligner2)
{
    def("dump_aligners", &EMAN::dump_aligners);
    def("dump_aligners_list", &EMAN::dump_aligners_list);
    class_< EMAN::Aligner, boost::noncopyable, EMAN_Aligner_Wrapper >("__Aligner", no_init)
        .def("align", pure_virtual((EMAN::EMData* (EMAN::Aligner::*)(EMAN::EMData*, EMAN::EMData*) const)&EMAN::Aligner::align), return_value_policy< manage_new_object >())
        .def("align", pure_virtual((EMAN::EMData* (EMAN::Aligner::*)(EMAN::EMData*, EMAN::EMData*, const std::string&, const EMAN::Dict&) const)&EMAN::Aligner::align), return_value_policy< manage_new_object >())
        .def("get_name", pure_virtual(&EMAN::Aligner::get_name))
        .def("get_desc", pure_virtual(&EMAN::Aligner::get_desc))
        .def("get_params", &EMAN::Aligner::get_params, &EMAN_Aligner_Wrapper::default_get_params)
        .def("set_params", &EMAN::Aligner::set_params, &EMAN_Aligner_Wrapper::default_set_params)
        .def("get_param_types", pure_virtual(&EMAN::Aligner::get_param_types))
    ;

    class_< EMAN::Factory<EMAN::Aligner>, boost::noncopyable >("Aligners", no_init)
        .def("get", (EMAN::Aligner* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&))&EMAN::Factory<EMAN::Aligner>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Aligner* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&, const EMAN::Dict&))&EMAN::Factory<EMAN::Aligner>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Aligner>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

    scope* EMAN_Ctf_scope = new scope(
    class_< EMAN::Ctf, boost::noncopyable, EMAN_Ctf_Wrapper >("Ctf", no_init)
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


    scope().attr("CTFOS") = (int)EMAN::Ctf::CTFOS;

    delete EMAN_Ctf_scope;

    class_< EMAN::SimpleCtf, bases< EMAN::Ctf > , boost::noncopyable, EMAN_SimpleCtf_Wrapper >("SimpleCtf", init<  >())
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
        .def("from_string", (int (EMAN::SimpleCtf::*)(const std::string&) )&EMAN::SimpleCtf::from_string, (int (EMAN_SimpleCtf_Wrapper::*)(const std::string&))&EMAN_SimpleCtf_Wrapper::default_from_string)
        .def("to_string", (std::string (EMAN::SimpleCtf::*)() const)&EMAN::SimpleCtf::to_string, (std::string (EMAN_SimpleCtf_Wrapper::*)() const)&EMAN_SimpleCtf_Wrapper::default_to_string)
        .def("from_dict", (void (EMAN::SimpleCtf::*)(const EMAN::Dict&) )&EMAN::SimpleCtf::from_dict, (void (EMAN_SimpleCtf_Wrapper::*)(const EMAN::Dict&))&EMAN_SimpleCtf_Wrapper::default_from_dict)
        .def("to_dict", (EMAN::Dict (EMAN::SimpleCtf::*)() const)&EMAN::SimpleCtf::to_dict, (EMAN::Dict (EMAN_SimpleCtf_Wrapper::*)() const)&EMAN_SimpleCtf_Wrapper::default_to_dict)
        .def("copy_from", (void (EMAN::SimpleCtf::*)(const EMAN::Ctf*) )&EMAN::SimpleCtf::copy_from, (void (EMAN_SimpleCtf_Wrapper::*)(const EMAN::Ctf*))&EMAN_SimpleCtf_Wrapper::default_copy_from)
        .def("equal", (bool (EMAN::SimpleCtf::*)(const EMAN::Ctf*) const)&EMAN::SimpleCtf::equal, (bool (EMAN_SimpleCtf_Wrapper::*)(const EMAN::Ctf*) const)&EMAN_SimpleCtf_Wrapper::default_equal)
        .def("get_defocus", (float (EMAN::SimpleCtf::*)() const)&EMAN::SimpleCtf::get_defocus, (float (EMAN_SimpleCtf_Wrapper::*)() const)&EMAN_SimpleCtf_Wrapper::default_get_defocus)
        .def("get_bfactor", (float (EMAN::SimpleCtf::*)() const)&EMAN::SimpleCtf::get_bfactor, (float (EMAN_SimpleCtf_Wrapper::*)() const)&EMAN_SimpleCtf_Wrapper::default_get_bfactor)
    ;

}

