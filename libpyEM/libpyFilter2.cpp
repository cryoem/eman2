
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <emdata.h>
#include <emobject.h>
#include <filter.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

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

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
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

    std::string get_desc() const {
        return call_method< std::string >(self, "get_desc");
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

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
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

    std::string get_desc() const {
        return call_method< std::string >(self, "get_desc");
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

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
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

    std::string get_desc() const {
        return call_method< std::string >(self, "get_desc");
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

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
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

    std::string get_desc() const {
        return call_method< std::string >(self, "get_desc");
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

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
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

    std::string get_desc() const {
        return call_method< std::string >(self, "get_desc");
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

    std::string get_name() const {
        return call_method< std::string >(self, "get_name");
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

    std::string get_desc() const {
        return call_method< std::string >(self, "get_desc");
    }

    PyObject* self;
};


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyFilter2)
{
    class_< EMAN::RealPixelFilter, boost::noncopyable, EMAN_RealPixelFilter_Wrapper >("RealPixelFilter", init<  >())
        .def("process", (void (EMAN::RealPixelFilter::*)(EMAN::EMData*) )&EMAN::RealPixelFilter::process, (void (EMAN_RealPixelFilter_Wrapper::*)(EMAN::EMData*))&EMAN_RealPixelFilter_Wrapper::default_process)
        .def("set_params", (void (EMAN::RealPixelFilter::*)(const EMAN::Dict&) )&EMAN::RealPixelFilter::set_params, (void (EMAN_RealPixelFilter_Wrapper::*)(const EMAN::Dict&))&EMAN_RealPixelFilter_Wrapper::default_set_params)
        .def("process_list", (void (EMAN::Filter::*)(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >&) )&EMAN::Filter::process_list, (void (EMAN_RealPixelFilter_Wrapper::*)(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >&))&EMAN_RealPixelFilter_Wrapper::default_process_list)
        .def("get_name", pure_virtual((std::string (EMAN::Filter::*)() const)&EMAN::Filter::get_name))
        .def("get_params", (EMAN::Dict (EMAN::Filter::*)() const)&EMAN::Filter::get_params, (EMAN::Dict (EMAN_RealPixelFilter_Wrapper::*)() const)&EMAN_RealPixelFilter_Wrapper::default_get_params)
        .def("get_param_types", (EMAN::TypeDict (EMAN::Filter::*)() const)&EMAN::Filter::get_param_types, (EMAN::TypeDict (EMAN_RealPixelFilter_Wrapper::*)() const)&EMAN_RealPixelFilter_Wrapper::default_get_param_types)
        .def("get_desc", pure_virtual((std::string (EMAN::Filter::*)() const)&EMAN::Filter::get_desc))
        .def("get_group_desc", &EMAN::RealPixelFilter::get_group_desc)
        .staticmethod("get_group_desc")
    ;

    class_< EMAN::BoxStatFilter, boost::noncopyable, EMAN_BoxStatFilter_Wrapper >("BoxStatFilter", init<  >())
        .def("process", (void (EMAN::BoxStatFilter::*)(EMAN::EMData*) )&EMAN::BoxStatFilter::process, (void (EMAN_BoxStatFilter_Wrapper::*)(EMAN::EMData*))&EMAN_BoxStatFilter_Wrapper::default_process)
        .def("process_list", (void (EMAN::Filter::*)(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >&) )&EMAN::Filter::process_list, (void (EMAN_BoxStatFilter_Wrapper::*)(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >&))&EMAN_BoxStatFilter_Wrapper::default_process_list)
        .def("get_name", pure_virtual((std::string (EMAN::Filter::*)() const)&EMAN::Filter::get_name))
        .def("get_params", (EMAN::Dict (EMAN::Filter::*)() const)&EMAN::Filter::get_params, (EMAN::Dict (EMAN_BoxStatFilter_Wrapper::*)() const)&EMAN_BoxStatFilter_Wrapper::default_get_params)
        .def("set_params", (void (EMAN::Filter::*)(const EMAN::Dict&) )&EMAN::Filter::set_params, (void (EMAN_BoxStatFilter_Wrapper::*)(const EMAN::Dict&))&EMAN_BoxStatFilter_Wrapper::default_set_params)
        .def("get_param_types", (EMAN::TypeDict (EMAN::Filter::*)() const)&EMAN::Filter::get_param_types, (EMAN::TypeDict (EMAN_BoxStatFilter_Wrapper::*)() const)&EMAN_BoxStatFilter_Wrapper::default_get_param_types)
        .def("get_desc", pure_virtual((std::string (EMAN::Filter::*)() const)&EMAN::Filter::get_desc))
        .def("get_group_desc", &EMAN::BoxStatFilter::get_group_desc)
        .staticmethod("get_group_desc")
    ;

    class_< EMAN::ComplexPixelFilter, boost::noncopyable, EMAN_ComplexPixelFilter_Wrapper >("ComplexPixelFilter", init<  >())
        .def("process", (void (EMAN::ComplexPixelFilter::*)(EMAN::EMData*) )&EMAN::ComplexPixelFilter::process, (void (EMAN_ComplexPixelFilter_Wrapper::*)(EMAN::EMData*))&EMAN_ComplexPixelFilter_Wrapper::default_process)
        .def("process_list", (void (EMAN::Filter::*)(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >&) )&EMAN::Filter::process_list, (void (EMAN_ComplexPixelFilter_Wrapper::*)(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >&))&EMAN_ComplexPixelFilter_Wrapper::default_process_list)
        .def("get_name", pure_virtual((std::string (EMAN::Filter::*)() const)&EMAN::Filter::get_name))
        .def("get_params", (EMAN::Dict (EMAN::Filter::*)() const)&EMAN::Filter::get_params, (EMAN::Dict (EMAN_ComplexPixelFilter_Wrapper::*)() const)&EMAN_ComplexPixelFilter_Wrapper::default_get_params)
        .def("set_params", (void (EMAN::Filter::*)(const EMAN::Dict&) )&EMAN::Filter::set_params, (void (EMAN_ComplexPixelFilter_Wrapper::*)(const EMAN::Dict&))&EMAN_ComplexPixelFilter_Wrapper::default_set_params)
        .def("get_param_types", (EMAN::TypeDict (EMAN::Filter::*)() const)&EMAN::Filter::get_param_types, (EMAN::TypeDict (EMAN_ComplexPixelFilter_Wrapper::*)() const)&EMAN_ComplexPixelFilter_Wrapper::default_get_param_types)
        .def("get_desc", pure_virtual((std::string (EMAN::Filter::*)() const)&EMAN::Filter::get_desc))
        .def("get_group_desc", &EMAN::ComplexPixelFilter::get_group_desc)
        .staticmethod("get_group_desc")
    ;

    class_< EMAN::CoordinateFilter, boost::noncopyable, EMAN_CoordinateFilter_Wrapper >("CoordinateFilter", init<  >())
        .def("process", (void (EMAN::CoordinateFilter::*)(EMAN::EMData*) )&EMAN::CoordinateFilter::process, (void (EMAN_CoordinateFilter_Wrapper::*)(EMAN::EMData*))&EMAN_CoordinateFilter_Wrapper::default_process)
        .def("process_list", (void (EMAN::Filter::*)(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >&) )&EMAN::Filter::process_list, (void (EMAN_CoordinateFilter_Wrapper::*)(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >&))&EMAN_CoordinateFilter_Wrapper::default_process_list)
        .def("get_name", pure_virtual((std::string (EMAN::Filter::*)() const)&EMAN::Filter::get_name))
        .def("get_params", (EMAN::Dict (EMAN::Filter::*)() const)&EMAN::Filter::get_params, (EMAN::Dict (EMAN_CoordinateFilter_Wrapper::*)() const)&EMAN_CoordinateFilter_Wrapper::default_get_params)
        .def("set_params", (void (EMAN::Filter::*)(const EMAN::Dict&) )&EMAN::Filter::set_params, (void (EMAN_CoordinateFilter_Wrapper::*)(const EMAN::Dict&))&EMAN_CoordinateFilter_Wrapper::default_set_params)
        .def("get_param_types", (EMAN::TypeDict (EMAN::Filter::*)() const)&EMAN::Filter::get_param_types, (EMAN::TypeDict (EMAN_CoordinateFilter_Wrapper::*)() const)&EMAN_CoordinateFilter_Wrapper::default_get_param_types)
        .def("get_desc", pure_virtual((std::string (EMAN::Filter::*)() const)&EMAN::Filter::get_desc))
        .def("get_group_desc", &EMAN::CoordinateFilter::get_group_desc)
        .staticmethod("get_group_desc")
    ;

    class_< EMAN::FourierFilter, boost::noncopyable, EMAN_FourierFilter_Wrapper >("FourierFilter", init<  >())
        .def("process", (void (EMAN::FourierFilter::*)(EMAN::EMData*) )&EMAN::FourierFilter::process, (void (EMAN_FourierFilter_Wrapper::*)(EMAN::EMData*))&EMAN_FourierFilter_Wrapper::default_process)
        .def("process_list", (void (EMAN::Filter::*)(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >&) )&EMAN::Filter::process_list, (void (EMAN_FourierFilter_Wrapper::*)(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >&))&EMAN_FourierFilter_Wrapper::default_process_list)
        .def("get_name", pure_virtual((std::string (EMAN::Filter::*)() const)&EMAN::Filter::get_name))
        .def("get_params", (EMAN::Dict (EMAN::Filter::*)() const)&EMAN::Filter::get_params, (EMAN::Dict (EMAN_FourierFilter_Wrapper::*)() const)&EMAN_FourierFilter_Wrapper::default_get_params)
        .def("set_params", (void (EMAN::Filter::*)(const EMAN::Dict&) )&EMAN::Filter::set_params, (void (EMAN_FourierFilter_Wrapper::*)(const EMAN::Dict&))&EMAN_FourierFilter_Wrapper::default_set_params)
        .def("get_param_types", (EMAN::TypeDict (EMAN::Filter::*)() const)&EMAN::Filter::get_param_types, (EMAN::TypeDict (EMAN_FourierFilter_Wrapper::*)() const)&EMAN_FourierFilter_Wrapper::default_get_param_types)
        .def("get_desc", pure_virtual((std::string (EMAN::Filter::*)() const)&EMAN::Filter::get_desc))
        .def("get_group_desc", &EMAN::FourierFilter::get_group_desc)
        .staticmethod("get_group_desc")
    ;

    class_< EMAN::NormalizeFilter, boost::noncopyable, EMAN_NormalizeFilter_Wrapper >("NormalizeFilter", init<  >())
        .def("process", (void (EMAN::NormalizeFilter::*)(EMAN::EMData*) )&EMAN::NormalizeFilter::process, (void (EMAN_NormalizeFilter_Wrapper::*)(EMAN::EMData*))&EMAN_NormalizeFilter_Wrapper::default_process)
        .def("process_list", (void (EMAN::Filter::*)(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >&) )&EMAN::Filter::process_list, (void (EMAN_NormalizeFilter_Wrapper::*)(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >&))&EMAN_NormalizeFilter_Wrapper::default_process_list)
        .def("get_name", pure_virtual((std::string (EMAN::Filter::*)() const)&EMAN::Filter::get_name))
        .def("get_params", (EMAN::Dict (EMAN::Filter::*)() const)&EMAN::Filter::get_params, (EMAN::Dict (EMAN_NormalizeFilter_Wrapper::*)() const)&EMAN_NormalizeFilter_Wrapper::default_get_params)
        .def("set_params", (void (EMAN::Filter::*)(const EMAN::Dict&) )&EMAN::Filter::set_params, (void (EMAN_NormalizeFilter_Wrapper::*)(const EMAN::Dict&))&EMAN_NormalizeFilter_Wrapper::default_set_params)
        .def("get_param_types", (EMAN::TypeDict (EMAN::Filter::*)() const)&EMAN::Filter::get_param_types, (EMAN::TypeDict (EMAN_NormalizeFilter_Wrapper::*)() const)&EMAN_NormalizeFilter_Wrapper::default_get_param_types)
        .def("get_desc", pure_virtual((std::string (EMAN::Filter::*)() const)&EMAN::Filter::get_desc))
        .def("get_group_desc", &EMAN::NormalizeFilter::get_group_desc)
        .staticmethod("get_group_desc")
    ;

    def("dump_filters", &EMAN::dump_filters);
    def("multi_filters", &EMAN::multi_filters);
    def("group_filters", &EMAN::group_filters);
    class_< EMAN::Factory<EMAN::Filter>, boost::noncopyable >("Filters", no_init)
        .def("get", (EMAN::Filter* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >))&EMAN::Factory<EMAN::Filter>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Filter* (*)(std::basic_string<char,std::char_traits<char>,std::allocator<char> >, const EMAN::Dict&))&EMAN::Factory<EMAN::Filter>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Filter>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

}

