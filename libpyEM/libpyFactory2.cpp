
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <aligner.h>
#include <averager.h>
#include <cmp.h>
#include <ctf.h>
#include <emdata.h>
#include <emobject.h>
#include <exception.h>
#include <filter.h>
#include <interp.h>
#include <projector.h>
#include <reconstructor.h>
#include <util.h>
#include <xydata.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct EMAN_Exception_Wrapper: EMAN::Exception
{
    EMAN_Exception_Wrapper(PyObject* self_, const EMAN::Exception& p0):
        EMAN::Exception(p0), self(self_) {}

    EMAN_Exception_Wrapper(PyObject* self_):
        EMAN::Exception(), self(self_) {}

    EMAN_Exception_Wrapper(PyObject* self_, const std::string& p0):
        EMAN::Exception(p0), self(self_) {}

    EMAN_Exception_Wrapper(PyObject* self_, const std::string& p0, int p1):
        EMAN::Exception(p0, p1), self(self_) {}

    EMAN_Exception_Wrapper(PyObject* self_, const std::string& p0, int p1, const std::string& p2):
        EMAN::Exception(p0, p1, p2), self(self_) {}

    EMAN_Exception_Wrapper(PyObject* self_, const std::string& p0, int p1, const std::string& p2, const std::string& p3):
        EMAN::Exception(p0, p1, p2, p3), self(self_) {}

    void set_desc(const std::string& p0) {
        call_method< void >(self, "set_desc", p0);
    }

    void default_set_desc(const std::string& p0) {
        EMAN::Exception::set_desc(p0);
    }

    const char* get_file() const {
        return call_method< const char* >(self, "get_file");
    }

    const char* default_get_file() const {
        return EMAN::Exception::get_file();
    }

    const char* get_desc() const {
        return call_method< const char* >(self, "get_desc");
    }

    const char* default_get_desc() const {
        return EMAN::Exception::get_desc();
    }

    int get_line_num() const {
        return call_method< int >(self, "get_line_num");
    }

    int default_get_line_num() const {
        return EMAN::Exception::get_line_num();
    }

    void set_objname(const std::string& p0) {
        call_method< void >(self, "set_objname", p0);
    }

    void default_set_objname(const std::string& p0) {
        EMAN::Exception::set_objname(p0);
    }

    const char* get_objname() const {
        return call_method< const char* >(self, "get_objname");
    }

    const char* default_get_objname() const {
        return EMAN::Exception::get_objname();
    }

    const char* what() const throw() {
        return call_method< const char* >(self, "what");
    }

    const char* default_what() const {
        return std::exception::what();
    }

    PyObject* self;
};

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_TypeDict_put_overloads_2_3, put, 2, 3)

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

struct EMAN_Cmp_Wrapper: EMAN::Cmp
{
    EMAN_Cmp_Wrapper(PyObject* self_, const EMAN::Cmp& p0):
        EMAN::Cmp(p0), self(self_) {}

    EMAN_Cmp_Wrapper(PyObject* self_):
        EMAN::Cmp(), self(self_) {}

    float cmp(EMAN::EMData* p0, EMAN::Transform* p1) const {
        return call_method< float >(self, "cmp", p0, p1);
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(self, "get_param_types");
    }

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
    }

    EMAN::Dict get_params() const {
        return call_method< EMAN::Dict >(self, "get_params");
    }

    EMAN::Dict default_get_params() const {
        return EMAN::Cmp::get_params();
    }

    void set_params(const EMAN::Dict& p0) {
        call_method< void >(self, "set_params", p0);
    }

    void default_set_params(const EMAN::Dict& p0) {
        EMAN::Cmp::set_params(p0);
    }

    PyObject* self;
};

struct EMAN_Averager_Wrapper: EMAN::Averager
{
    EMAN_Averager_Wrapper(PyObject* self_, const EMAN::Averager& p0):
        EMAN::Averager(p0), self(self_) {}

    EMAN_Averager_Wrapper(PyObject* self_):
        EMAN::Averager(), self(self_) {}

    void add_image(EMAN::EMData* p0) {
        call_method< void >(self, "add_image", p0);
    }

    EMAN::EMData* finish() {
        return call_method< EMAN::EMData* >(self, "finish");
    }

    void add_image_list(const std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >& p0) {
        call_method< void >(self, "add_image_list", p0);
    }

    void default_add_image_list(const std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >& p0) {
        EMAN::Averager::add_image_list(p0);
    }

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
    }

    void set_params(const EMAN::Dict& p0) {
        call_method< void >(self, "set_params", p0);
    }

    void default_set_params(const EMAN::Dict& p0) {
        EMAN::Averager::set_params(p0);
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(self, "get_param_types");
    }

    EMAN::TypeDict default_get_param_types() const {
        return EMAN::Averager::get_param_types();
    }

    PyObject* self;
};

struct EMAN_Projector_Wrapper: EMAN::Projector
{
    EMAN_Projector_Wrapper(PyObject* self_, const EMAN::Projector& p0):
        EMAN::Projector(p0), self(self_) {}

    EMAN_Projector_Wrapper(PyObject* self_):
        EMAN::Projector(), self(self_) {}

    EMAN::EMData* project3d(EMAN::EMData* p0) const {
        return call_method< EMAN::EMData* >(self, "project3d", p0);
    }

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
    }

    EMAN::Dict get_params() const {
        return call_method< EMAN::Dict >(self, "get_params");
    }

    EMAN::Dict default_get_params() const {
        return EMAN::Projector::get_params();
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(self, "get_param_types");
    }

    EMAN::TypeDict default_get_param_types() const {
        return EMAN::Projector::get_param_types();
    }

    PyObject* self;
};

struct EMAN_Reconstructor_Wrapper: EMAN::Reconstructor
{
    EMAN_Reconstructor_Wrapper(PyObject* self_, const EMAN::Reconstructor& p0):
        EMAN::Reconstructor(p0), self(self_) {}

    EMAN_Reconstructor_Wrapper(PyObject* self_):
        EMAN::Reconstructor(), self(self_) {}

    int setup() {
        return call_method< int >(self, "setup");
    }

    int insert_slice(EMAN::EMData* p0, const EMAN::Rotation& p1) {
        return call_method< int >(self, "insert_slice", p0, p1);
    }

    EMAN::EMData* finish() {
        return call_method< EMAN::EMData* >(self, "finish");
    }

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
    }

    EMAN::Dict get_params() const {
        return call_method< EMAN::Dict >(self, "get_params");
    }

    EMAN::Dict default_get_params() const {
        return EMAN::Reconstructor::get_params();
    }

    void set_params(const EMAN::Dict& p0) {
        call_method< void >(self, "set_params", p0);
    }

    void default_set_params(const EMAN::Dict& p0) {
        EMAN::Reconstructor::set_params(p0);
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(self, "get_param_types");
    }

    PyObject* self;
};

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

    EMAN::TypeDict default_get_param_types() const {
        return EMAN::Filter::get_param_types();
    }

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
    }

    PyObject* self;
};

struct EMAN_RealPixelFilter_Wrapper: EMAN::RealPixelFilter
{
    EMAN_RealPixelFilter_Wrapper(PyObject* self_, const EMAN::RealPixelFilter& p0):
        EMAN::RealPixelFilter(p0), self(self_) {}

    EMAN_RealPixelFilter_Wrapper(PyObject* self_):
        EMAN::RealPixelFilter(), self(self_) {}

    void process(EMAN::EMData* p0) {
        call_method< void >(self, "process", p0);
    }

    void default_process(EMAN::EMData* p0) {
        EMAN::RealPixelFilter::process(p0);
    }

    void set_params(const EMAN::Dict& p0) {
        call_method< void >(self, "set_params", p0);
    }

    void default_set_params(const EMAN::Dict& p0) {
        EMAN::RealPixelFilter::set_params(p0);
    }

    void process_pixel(float* p0) const {
        call_method< void >(self, "process_pixel", p0);
    }

    void calc_locals(EMAN::EMData* p0) {
        call_method< void >(self, "calc_locals", p0);
    }

    void default_calc_locals(EMAN::EMData* p0) {
        EMAN::RealPixelFilter::calc_locals(p0);
    }

    void normalize(EMAN::EMData* p0) const {
        call_method< void >(self, "normalize", p0);
    }

    void default_normalize(EMAN::EMData* p0) const {
        EMAN::RealPixelFilter::normalize(p0);
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

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(self, "get_param_types");
    }

    EMAN::TypeDict default_get_param_types() const {
        return EMAN::Filter::get_param_types();
    }

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
    }

    PyObject* self;
};

struct EMAN_BoxStatFilter_Wrapper: EMAN::BoxStatFilter
{
    EMAN_BoxStatFilter_Wrapper(PyObject* self_, const EMAN::BoxStatFilter& p0):
        EMAN::BoxStatFilter(p0), self(self_) {}

    EMAN_BoxStatFilter_Wrapper(PyObject* self_):
        EMAN::BoxStatFilter(), self(self_) {}

    void process(EMAN::EMData* p0) {
        call_method< void >(self, "process", p0);
    }

    void default_process(EMAN::EMData* p0) {
        EMAN::BoxStatFilter::process(p0);
    }

    void process_pixel(float* p0, const float* p1, int p2) const {
        call_method< void >(self, "process_pixel", p0, p1, p2);
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

    EMAN::TypeDict default_get_param_types() const {
        return EMAN::Filter::get_param_types();
    }

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
    }

    PyObject* self;
};

struct EMAN_ComplexPixelFilter_Wrapper: EMAN::ComplexPixelFilter
{
    EMAN_ComplexPixelFilter_Wrapper(PyObject* self_, const EMAN::ComplexPixelFilter& p0):
        EMAN::ComplexPixelFilter(p0), self(self_) {}

    EMAN_ComplexPixelFilter_Wrapper(PyObject* self_):
        EMAN::ComplexPixelFilter(), self(self_) {}

    void process(EMAN::EMData* p0) {
        call_method< void >(self, "process", p0);
    }

    void default_process(EMAN::EMData* p0) {
        EMAN::ComplexPixelFilter::process(p0);
    }

    void process_pixel(float* p0) const {
        call_method< void >(self, "process_pixel", p0);
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

    EMAN::TypeDict default_get_param_types() const {
        return EMAN::Filter::get_param_types();
    }

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
    }

    PyObject* self;
};

struct EMAN_CoordinateFilter_Wrapper: EMAN::CoordinateFilter
{
    EMAN_CoordinateFilter_Wrapper(PyObject* self_, const EMAN::CoordinateFilter& p0):
        EMAN::CoordinateFilter(p0), self(self_) {}

    EMAN_CoordinateFilter_Wrapper(PyObject* self_):
        EMAN::CoordinateFilter(), self(self_) {}

    void process(EMAN::EMData* p0) {
        call_method< void >(self, "process", p0);
    }

    void default_process(EMAN::EMData* p0) {
        EMAN::CoordinateFilter::process(p0);
    }

    void process_pixel(float* p0, int p1, int p2, int p3) const {
        call_method< void >(self, "process_pixel", p0, p1, p2, p3);
    }

    void calc_locals(EMAN::EMData* p0) {
        call_method< void >(self, "calc_locals", p0);
    }

    void default_calc_locals(EMAN::EMData* p0) {
        EMAN::CoordinateFilter::calc_locals(p0);
    }

    bool is_valid() const {
        return call_method< bool >(self, "is_valid");
    }

    bool default_is_valid() const {
        return EMAN::CoordinateFilter::is_valid();
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

    EMAN::TypeDict default_get_param_types() const {
        return EMAN::Filter::get_param_types();
    }

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
    }

    PyObject* self;
};

struct EMAN_FourierFilter_Wrapper: EMAN::FourierFilter
{
    EMAN_FourierFilter_Wrapper(PyObject* self_, const EMAN::FourierFilter& p0):
        EMAN::FourierFilter(p0), self(self_) {}

    EMAN_FourierFilter_Wrapper(PyObject* self_):
        EMAN::FourierFilter(), self(self_) {}

    void process(EMAN::EMData* p0) {
        call_method< void >(self, "process", p0);
    }

    void default_process(EMAN::EMData* p0) {
        EMAN::FourierFilter::process(p0);
    }

    void create_radial_func(std::vector<float,std::allocator<float> >& p0) const {
        call_method< void >(self, "create_radial_func", p0);
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

    EMAN::TypeDict default_get_param_types() const {
        return EMAN::Filter::get_param_types();
    }

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
    }

    PyObject* self;
};

struct EMAN_NormalizeFilter_Wrapper: EMAN::NormalizeFilter
{
    EMAN_NormalizeFilter_Wrapper(PyObject* self_, const EMAN::NormalizeFilter& p0):
        EMAN::NormalizeFilter(p0), self(self_) {}

    EMAN_NormalizeFilter_Wrapper(PyObject* self_):
        EMAN::NormalizeFilter(), self(self_) {}

    void process(EMAN::EMData* p0) {
        call_method< void >(self, "process", p0);
    }

    void default_process(EMAN::EMData* p0) {
        EMAN::NormalizeFilter::process(p0);
    }

    float calc_sigma(EMAN::EMData* p0) const {
        return call_method< float >(self, "calc_sigma", p0);
    }

    float default_calc_sigma(EMAN::EMData* p0) const {
        return EMAN::NormalizeFilter::calc_sigma(p0);
    }

    float calc_mean(EMAN::EMData* p0) const {
        return call_method< float >(self, "calc_mean", p0);
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

    EMAN::TypeDict default_get_param_types() const {
        return EMAN::Filter::get_param_types();
    }

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
    }

    PyObject* self;
};


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyFactory2)
{
    class_< EMAN::Exception, EMAN_Exception_Wrapper >("Exception", init< const EMAN::Exception& >())
        .def(init< optional< const std::string&, int, const std::string&, const std::string& > >())
        .def("set_desc", &EMAN::Exception::set_desc, &EMAN_Exception_Wrapper::default_set_desc)
        .def("get_file", &EMAN::Exception::get_file, &EMAN_Exception_Wrapper::default_get_file)
        .def("get_desc", &EMAN::Exception::get_desc, &EMAN_Exception_Wrapper::default_get_desc)
        .def("get_line_num", &EMAN::Exception::get_line_num, &EMAN_Exception_Wrapper::default_get_line_num)
        .def("set_objname", &EMAN::Exception::set_objname, &EMAN_Exception_Wrapper::default_set_objname)
        .def("get_objname", &EMAN::Exception::get_objname, &EMAN_Exception_Wrapper::default_get_objname)
        .def("what", (const char* (std::exception::*)() const throw())&std::exception::what, (const char* (EMAN_Exception_Wrapper::*)() const)&EMAN_Exception_Wrapper::default_what)
        .def("dump", &EMAN::Exception::dump)
    ;

    scope* EMAN_EMObject_scope = new scope(
    class_< EMAN::EMObject >("EMObject", init<  >())
        .def(init< const EMAN::EMObject& >())
        .def(init< int >())
        .def(init< float >())
        .def(init< double >())
        .def(init< const char* >())
        .def(init< std::string >())
        .def(init< EMAN::EMData* >())
        .def(init< EMAN::XYData* >())
        .def(init< const std::vector<float,std::allocator<float> >& >())
        .def("get_farray", &EMAN::EMObject::get_farray)
        .def("is_null", &EMAN::EMObject::is_null)
        .def("to_str", &EMAN::EMObject::to_str)
        .def("get_object_type_name", &EMAN::EMObject::get_object_type_name)
        .staticmethod("get_object_type_name")
        .def("__int__", &EMAN::EMObject::operator int)
        .def("__float__", &EMAN::EMObject::operator float)
        .def("__float__", &EMAN::EMObject::operator double)
        .def("to_EMAN_EMData", &EMAN::EMObject::operator EMAN::EMData*, return_internal_reference< 1 >())
        .def("to_EMAN_XYData", &EMAN::EMObject::operator EMAN::XYData*, return_internal_reference< 1 >())
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

    class_< EMAN::TypeDict >("TypeDict", init<  >())
        .def(init< const EMAN::TypeDict& >())
        .def("keys", &EMAN::TypeDict::keys)
        .def("size", &EMAN::TypeDict::size)
        .def("put", &EMAN::TypeDict::put, EMAN_TypeDict_put_overloads_2_3())
        .def("get_type", &EMAN::TypeDict::get_type)
        .def("get_desc", &EMAN::TypeDict::get_desc)
        .def("dump", &EMAN::TypeDict::dump)
    ;

    class_< EMAN::Aligner, boost::noncopyable, EMAN_Aligner_Wrapper >("Aligner", init<  >())
        .def("align", pure_virtual(&EMAN::Aligner::align), return_value_policy< manage_new_object >())
        .def("get_params", &EMAN::Aligner::get_params, &EMAN_Aligner_Wrapper::default_get_params)
        .def("set_params", &EMAN::Aligner::set_params, &EMAN_Aligner_Wrapper::default_set_params)
        .def("get_param_types", pure_virtual(&EMAN::Aligner::get_param_types))
        .def("get_name", pure_virtual(&EMAN::Aligner::get_name))
    ;

    class_< EMAN::Factory<EMAN::Aligner>, boost::noncopyable >("Aligners", no_init)
        .def("get", (EMAN::Aligner* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >))&EMAN::Factory<EMAN::Aligner>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Aligner* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&))&EMAN::Factory<EMAN::Aligner>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Aligner>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

    class_< EMAN::Cmp, boost::noncopyable, EMAN_Cmp_Wrapper >("Cmp", init<  >())
        .def("cmp", pure_virtual(&EMAN::Cmp::cmp))
        .def("get_param_types", pure_virtual(&EMAN::Cmp::get_param_types))
        .def("get_name", pure_virtual(&EMAN::Cmp::get_name))
        .def("get_params", &EMAN::Cmp::get_params, &EMAN_Cmp_Wrapper::default_get_params)
        .def("set_params", &EMAN::Cmp::set_params, &EMAN_Cmp_Wrapper::default_set_params)
    ;

    class_< EMAN::Factory<EMAN::Cmp>, boost::noncopyable >("Cmps", no_init)
        .def("get", (EMAN::Cmp* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >))&EMAN::Factory<EMAN::Cmp>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Cmp* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&))&EMAN::Factory<EMAN::Cmp>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Cmp>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

    class_< EMAN::Averager, boost::noncopyable, EMAN_Averager_Wrapper >("Averager", init<  >())
        .def("add_image", pure_virtual(&EMAN::Averager::add_image))
        .def("finish", pure_virtual(&EMAN::Averager::finish), return_value_policy< manage_new_object >())
        .def("add_image_list", &EMAN::Averager::add_image_list, &EMAN_Averager_Wrapper::default_add_image_list)
        .def("get_name", pure_virtual(&EMAN::Averager::get_name))
        .def("set_params", &EMAN::Averager::set_params, &EMAN_Averager_Wrapper::default_set_params)
        .def("get_param_types", &EMAN::Averager::get_param_types, &EMAN_Averager_Wrapper::default_get_param_types)
    ;

    class_< EMAN::Factory<EMAN::Averager>, boost::noncopyable >("Averagers", no_init)
        .def("get", (EMAN::Averager* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >))&EMAN::Factory<EMAN::Averager>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Averager* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&))&EMAN::Factory<EMAN::Averager>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Averager>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

    class_< EMAN::Projector, boost::noncopyable, EMAN_Projector_Wrapper >("Projector", init<  >())
        .def("project3d", pure_virtual(&EMAN::Projector::project3d), return_value_policy< manage_new_object >())
        .def("get_name", pure_virtual(&EMAN::Projector::get_name))
        .def("get_params", &EMAN::Projector::get_params, &EMAN_Projector_Wrapper::default_get_params)
        .def("get_param_types", &EMAN::Projector::get_param_types, &EMAN_Projector_Wrapper::default_get_param_types)
        .def("set_params", &EMAN::Projector::set_params)
    ;

    class_< EMAN::Factory<EMAN::Projector>, boost::noncopyable >("Projectors", no_init)
        .def("get", (EMAN::Projector* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >))&EMAN::Factory<EMAN::Projector>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Projector* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&))&EMAN::Factory<EMAN::Projector>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Projector>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

    class_< EMAN::Reconstructor, boost::noncopyable, EMAN_Reconstructor_Wrapper >("Reconstructor", init<  >())
        .def("setup", pure_virtual(&EMAN::Reconstructor::setup))
        .def("insert_slice", pure_virtual(&EMAN::Reconstructor::insert_slice))
        .def("finish", pure_virtual(&EMAN::Reconstructor::finish), return_value_policy< manage_new_object >())
        .def("get_name", pure_virtual(&EMAN::Reconstructor::get_name))
        .def("get_params", &EMAN::Reconstructor::get_params, &EMAN_Reconstructor_Wrapper::default_get_params)
        .def("set_params", &EMAN::Reconstructor::set_params, &EMAN_Reconstructor_Wrapper::default_set_params)
        .def("get_param_types", pure_virtual(&EMAN::Reconstructor::get_param_types))
    ;

    class_< EMAN::Factory<EMAN::Reconstructor>, boost::noncopyable >("Reconstructors", no_init)
        .def("get", (EMAN::Reconstructor* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >))&EMAN::Factory<EMAN::Reconstructor>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Reconstructor* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&))&EMAN::Factory<EMAN::Reconstructor>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Reconstructor>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

    class_< EMAN::Filter, boost::noncopyable, EMAN_Filter_Wrapper >("Filter", init<  >())
        .def("process", &EMAN::Filter::process, &EMAN_Filter_Wrapper::default_process)
        .def("process_list", &EMAN::Filter::process_list, &EMAN_Filter_Wrapper::default_process_list)
        .def("get_params", &EMAN::Filter::get_params, &EMAN_Filter_Wrapper::default_get_params)
        .def("set_params", &EMAN::Filter::set_params, &EMAN_Filter_Wrapper::default_set_params)
        .def("get_param_types", &EMAN::Filter::get_param_types, &EMAN_Filter_Wrapper::default_get_param_types)
        .def("get_name", pure_virtual(&EMAN::Filter::get_name))
        .def("get_group_desc", &EMAN::Filter::get_group_desc)
        .staticmethod("get_group_desc")
    ;

    class_< EMAN::RealPixelFilter, bases< EMAN::Filter > , boost::noncopyable, EMAN_RealPixelFilter_Wrapper >("RealPixelFilter", init<  >())
        .def("process", (void (EMAN::RealPixelFilter::*)(EMAN::EMData*) )&EMAN::RealPixelFilter::process, (void (EMAN_RealPixelFilter_Wrapper::*)(EMAN::EMData*))&EMAN_RealPixelFilter_Wrapper::default_process)
        .def("set_params", (void (EMAN::RealPixelFilter::*)(const EMAN::Dict&) )&EMAN::RealPixelFilter::set_params, (void (EMAN_RealPixelFilter_Wrapper::*)(const EMAN::Dict&))&EMAN_RealPixelFilter_Wrapper::default_set_params)
        .def("process_list", &EMAN::Filter::process_list, &EMAN_RealPixelFilter_Wrapper::default_process_list)
        .def("get_params", &EMAN::Filter::get_params, &EMAN_RealPixelFilter_Wrapper::default_get_params)
        .def("get_param_types", &EMAN::Filter::get_param_types, &EMAN_RealPixelFilter_Wrapper::default_get_param_types)
        .def("get_name", pure_virtual(&EMAN::Filter::get_name))
        .def("get_group_desc", &EMAN::RealPixelFilter::get_group_desc)
        .staticmethod("get_group_desc")
    ;

    class_< EMAN::BoxStatFilter, bases< EMAN::Filter > , boost::noncopyable, EMAN_BoxStatFilter_Wrapper >("BoxStatFilter", init<  >())
        .def("process", (void (EMAN::BoxStatFilter::*)(EMAN::EMData*) )&EMAN::BoxStatFilter::process, (void (EMAN_BoxStatFilter_Wrapper::*)(EMAN::EMData*))&EMAN_BoxStatFilter_Wrapper::default_process)
        .def("process_list", &EMAN::Filter::process_list, &EMAN_BoxStatFilter_Wrapper::default_process_list)
        .def("get_params", &EMAN::Filter::get_params, &EMAN_BoxStatFilter_Wrapper::default_get_params)
        .def("set_params", &EMAN::Filter::set_params, &EMAN_BoxStatFilter_Wrapper::default_set_params)
        .def("get_param_types", &EMAN::Filter::get_param_types, &EMAN_BoxStatFilter_Wrapper::default_get_param_types)
        .def("get_name", pure_virtual(&EMAN::Filter::get_name))
        .def("get_group_desc", &EMAN::BoxStatFilter::get_group_desc)
        .staticmethod("get_group_desc")
    ;

    class_< EMAN::ComplexPixelFilter, bases< EMAN::Filter > , boost::noncopyable, EMAN_ComplexPixelFilter_Wrapper >("ComplexPixelFilter", init<  >())
        .def("process", (void (EMAN::ComplexPixelFilter::*)(EMAN::EMData*) )&EMAN::ComplexPixelFilter::process, (void (EMAN_ComplexPixelFilter_Wrapper::*)(EMAN::EMData*))&EMAN_ComplexPixelFilter_Wrapper::default_process)
        .def("process_list", &EMAN::Filter::process_list, &EMAN_ComplexPixelFilter_Wrapper::default_process_list)
        .def("get_params", &EMAN::Filter::get_params, &EMAN_ComplexPixelFilter_Wrapper::default_get_params)
        .def("set_params", &EMAN::Filter::set_params, &EMAN_ComplexPixelFilter_Wrapper::default_set_params)
        .def("get_param_types", &EMAN::Filter::get_param_types, &EMAN_ComplexPixelFilter_Wrapper::default_get_param_types)
        .def("get_name", pure_virtual(&EMAN::Filter::get_name))
        .def("get_group_desc", &EMAN::ComplexPixelFilter::get_group_desc)
        .staticmethod("get_group_desc")
    ;

    class_< EMAN::CoordinateFilter, bases< EMAN::Filter > , boost::noncopyable, EMAN_CoordinateFilter_Wrapper >("CoordinateFilter", init<  >())
        .def("process", (void (EMAN::CoordinateFilter::*)(EMAN::EMData*) )&EMAN::CoordinateFilter::process, (void (EMAN_CoordinateFilter_Wrapper::*)(EMAN::EMData*))&EMAN_CoordinateFilter_Wrapper::default_process)
        .def("process_list", &EMAN::Filter::process_list, &EMAN_CoordinateFilter_Wrapper::default_process_list)
        .def("get_params", &EMAN::Filter::get_params, &EMAN_CoordinateFilter_Wrapper::default_get_params)
        .def("set_params", &EMAN::Filter::set_params, &EMAN_CoordinateFilter_Wrapper::default_set_params)
        .def("get_param_types", &EMAN::Filter::get_param_types, &EMAN_CoordinateFilter_Wrapper::default_get_param_types)
        .def("get_name", pure_virtual(&EMAN::Filter::get_name))
        .def("get_group_desc", &EMAN::CoordinateFilter::get_group_desc)
        .staticmethod("get_group_desc")
    ;

    class_< EMAN::FourierFilter, bases< EMAN::Filter > , boost::noncopyable, EMAN_FourierFilter_Wrapper >("FourierFilter", init<  >())
        .def("process", (void (EMAN::FourierFilter::*)(EMAN::EMData*) )&EMAN::FourierFilter::process, (void (EMAN_FourierFilter_Wrapper::*)(EMAN::EMData*))&EMAN_FourierFilter_Wrapper::default_process)
        .def("process_list", &EMAN::Filter::process_list, &EMAN_FourierFilter_Wrapper::default_process_list)
        .def("get_params", &EMAN::Filter::get_params, &EMAN_FourierFilter_Wrapper::default_get_params)
        .def("set_params", &EMAN::Filter::set_params, &EMAN_FourierFilter_Wrapper::default_set_params)
        .def("get_param_types", &EMAN::Filter::get_param_types, &EMAN_FourierFilter_Wrapper::default_get_param_types)
        .def("get_name", pure_virtual(&EMAN::Filter::get_name))
        .def("get_group_desc", &EMAN::FourierFilter::get_group_desc)
        .staticmethod("get_group_desc")
    ;

    class_< EMAN::NormalizeFilter, bases< EMAN::Filter > , boost::noncopyable, EMAN_NormalizeFilter_Wrapper >("NormalizeFilter", init<  >())
        .def("process", (void (EMAN::NormalizeFilter::*)(EMAN::EMData*) )&EMAN::NormalizeFilter::process, (void (EMAN_NormalizeFilter_Wrapper::*)(EMAN::EMData*))&EMAN_NormalizeFilter_Wrapper::default_process)
        .def("process_list", &EMAN::Filter::process_list, &EMAN_NormalizeFilter_Wrapper::default_process_list)
        .def("get_params", &EMAN::Filter::get_params, &EMAN_NormalizeFilter_Wrapper::default_get_params)
        .def("set_params", &EMAN::Filter::set_params, &EMAN_NormalizeFilter_Wrapper::default_set_params)
        .def("get_param_types", &EMAN::Filter::get_param_types, &EMAN_NormalizeFilter_Wrapper::default_get_param_types)
        .def("get_name", pure_virtual(&EMAN::Filter::get_name))
        .def("get_group_desc", &EMAN::NormalizeFilter::get_group_desc)
        .staticmethod("get_group_desc")
    ;

    class_< EMAN::Factory<EMAN::Filter>, boost::noncopyable >("Filters", no_init)
        .def("get", (EMAN::Filter* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >))&EMAN::Factory<EMAN::Filter>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Filter* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&))&EMAN::Factory<EMAN::Filter>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Filter>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

    def("dump_aligners", &EMAN::dump_aligners);
    def("dump_averagers", &EMAN::dump_averagers);
    def("dump_cmps", &EMAN::dump_cmps);
    def("dump_filters", &EMAN::dump_filters);
    def("multi_filters", &EMAN::multi_filters);
    def("group_filters", &EMAN::group_filters);
    def("dump_projectors", &EMAN::dump_projectors);
    def("dump_reconstructors", &EMAN::dump_reconstructors);
}

