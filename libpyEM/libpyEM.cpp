
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
#include <emutil.h>
#include <filter.h>
#include <floatstat.h>
#include <geometry.h>
#include <imageio.h>
#include <log.h>
#include <projector.h>
#include <pyem.h>
#include <pylist.h>
#include <reconstructor.h>
#include <transform.h>
#include <util.h>
#include <xydata.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

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

    EMAN::EMData* average(const std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >& p0) const {
        return call_method< EMAN::EMData* >(self, "average", p0);
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

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_to_mass_center_overloads_0_1, to_mass_center, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_fast_translate_overloads_0_1, fast_translate, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_rotate_translate_overloads_0_5, rotate_translate, 0, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_fast_rotate_translate_overloads_0_1, fast_rotate_translate, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_little_big_dot_overloads_1_2, little_big_dot, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_ccf_overloads_1_3, calc_ccf, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_make_rotational_footprint_overloads_0_2, make_rotational_footprint, 0, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_ccfx_overloads_1_4, calc_ccfx, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_mutual_correlation_overloads_1_3, calc_mutual_correlation, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_unwrap_overloads_0_6, unwrap, 0, 6)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_add_random_noise_overloads_4_5, add_random_noise, 4, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_auto_mask_overloads_1_2, auto_mask, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_hist_overloads_1_4, calc_hist, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_dist_overloads_1_2, calc_dist, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_calc_flcf_overloads_1_3, calc_flcf, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_dot_overloads_1_2, dot, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_common_lines_overloads_2_5, common_lines, 2, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_common_lines_real_overloads_2_4, common_lines_real, 2, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_cut_slice_overloads_2_6, cut_slice, 2, 6)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_EMData_uncut_slice_overloads_2_5, uncut_slice, 2, 5)

struct EMAN_Ctf_Wrapper: EMAN::Ctf
{
    EMAN_Ctf_Wrapper(PyObject* self_, const EMAN::Ctf& p0):
        EMAN::Ctf(p0), self(self_) {}

    EMAN_Ctf_Wrapper(PyObject* self_):
        EMAN::Ctf(), self(self_) {}

    bool cmp() const {
        return call_method< bool >(self, "cmp");
    }

    int from_string(std::string p0) {
        return call_method< int >(self, "from_string", p0);
    }

    std::string to_string() const {
        return call_method< std::string >(self, "to_string");
    }

    int from_dict(std::map<std::basic_string<char, std::char_traits<char>, std::allocator<char> >,EMAN::EMObject,std::less<std::basic_string<char, std::char_traits<char>, std::allocator<char> > >,std::allocator<std::pair<const std::basic_string<char, std::char_traits<char>, std::allocator<char> >, EMAN::EMObject> > >& p0) {
        return call_method< int >(self, "from_dict", p0);
    }

    int to_dict(std::map<std::basic_string<char, std::char_traits<char>, std::allocator<char> >,EMAN::EMObject,std::less<std::basic_string<char, std::char_traits<char>, std::allocator<char> > >,std::allocator<std::pair<const std::basic_string<char, std::char_traits<char>, std::allocator<char> >, EMAN::EMObject> > >& p0) const {
        return call_method< int >(self, "to_dict", p0);
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

    bool cmp() const {
        return call_method< bool >(self, "cmp");
    }

    bool default_cmp() const {
        return EMAN::SimpleCtf::cmp();
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

    int from_dict(std::map<std::basic_string<char, std::char_traits<char>, std::allocator<char> >,EMAN::EMObject,std::less<std::basic_string<char, std::char_traits<char>, std::allocator<char> > >,std::allocator<std::pair<const std::basic_string<char, std::char_traits<char>, std::allocator<char> >, EMAN::EMObject> > >& p0) {
        return call_method< int >(self, "from_dict", p0);
    }

    int default_from_dict(std::map<std::basic_string<char, std::char_traits<char>, std::allocator<char> >,EMAN::EMObject,std::less<std::basic_string<char, std::char_traits<char>, std::allocator<char> > >,std::allocator<std::pair<const std::basic_string<char, std::char_traits<char>, std::allocator<char> >, EMAN::EMObject> > >& p0) {
        return EMAN::SimpleCtf::from_dict(p0);
    }

    int to_dict(std::map<std::basic_string<char, std::char_traits<char>, std::allocator<char> >,EMAN::EMObject,std::less<std::basic_string<char, std::char_traits<char>, std::allocator<char> > >,std::allocator<std::pair<const std::basic_string<char, std::char_traits<char>, std::allocator<char> >, EMAN::EMObject> > >& p0) const {
        return call_method< int >(self, "to_dict", p0);
    }

    int default_to_dict(std::map<std::basic_string<char, std::char_traits<char>, std::allocator<char> >,EMAN::EMObject,std::less<std::basic_string<char, std::char_traits<char>, std::allocator<char> > >,std::allocator<std::pair<const std::basic_string<char, std::char_traits<char>, std::allocator<char> >, EMAN::EMObject> > >& p0) const {
        return EMAN::SimpleCtf::to_dict(p0);
    }

    PyObject* self;
};

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_SimpleCtf_compute_map_overloads_1_2, compute_map, 1, 2)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_SimpleCtf_compute_curve_overloads_2_3, EMAN::SimpleCtf::compute_curve, 2, 3)

struct EMAN_ImageIO_Wrapper: EMAN::ImageIO
{
    EMAN_ImageIO_Wrapper(PyObject* self_, const EMAN::ImageIO& p0):
        EMAN::ImageIO(p0), self(self_) {}

    EMAN_ImageIO_Wrapper(PyObject* self_):
        EMAN::ImageIO(), self(self_) {}

    int read_header(std::map<std::basic_string<char, std::char_traits<char>, std::allocator<char> >,EMAN::EMObject,std::less<std::basic_string<char, std::char_traits<char>, std::allocator<char> > >,std::allocator<std::pair<const std::basic_string<char, std::char_traits<char>, std::allocator<char> >, EMAN::EMObject> > >& p0, int p1, const EMAN::Region* p2, bool p3) {
        return call_method< int >(self, "read_header", p0, p1, p2, p3);
    }

    int write_header(std::map<std::basic_string<char, std::char_traits<char>, std::allocator<char> >,EMAN::EMObject,std::less<std::basic_string<char, std::char_traits<char>, std::allocator<char> > >,std::allocator<std::pair<const std::basic_string<char, std::char_traits<char>, std::allocator<char> >, EMAN::EMObject> > >& p0, int p1) {
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


    class_< EMAN::FloatStat >("FloatStat", init< const EMAN::FloatStat& >())
        .def(init< int >())
        .def("clear", &EMAN::FloatStat::clear)
        .def("add", &EMAN::FloatStat::add)
        .def("mean", &EMAN::FloatStat::mean)
        .def("sigma", &EMAN::FloatStat::sigma)
        .def("min", &EMAN::FloatStat::min)
        .def("max", &EMAN::FloatStat::max)
        .def("median", &EMAN::FloatStat::median)
        .def("num_greater", &EMAN::FloatStat::num_greater)
    ;

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

    class_< EMAN::Aligner, boost::noncopyable, EMAN_Aligner_Wrapper >("Aligner", init<  >())
        .def("align", pure_virtual(&EMAN::Aligner::align), return_value_policy< manage_new_object >())
        .def("get_params", &EMAN::Aligner::get_params, &EMAN_Aligner_Wrapper::default_get_params)
        .def("set_params", &EMAN::Aligner::set_params, &EMAN_Aligner_Wrapper::default_set_params)
        .def("get_param_types", pure_virtual(&EMAN::Aligner::get_param_types))
        .def("get_name", pure_virtual(&EMAN::Aligner::get_name))
    ;

    class_< EMAN::Cmp, boost::noncopyable, EMAN_Cmp_Wrapper >("Cmp", init<  >())
        .def("cmp", pure_virtual(&EMAN::Cmp::cmp))
        .def("get_param_types", pure_virtual(&EMAN::Cmp::get_param_types))
        .def("get_name", pure_virtual(&EMAN::Cmp::get_name))
        .def("get_params", &EMAN::Cmp::get_params, &EMAN_Cmp_Wrapper::default_get_params)
        .def("set_params", &EMAN::Cmp::set_params, &EMAN_Cmp_Wrapper::default_set_params)
    ;

    class_< EMAN::Averager, boost::noncopyable, EMAN_Averager_Wrapper >("Averager", init<  >())
        .def("average", pure_virtual(&EMAN::Averager::average), return_value_policy< manage_new_object >())
        .def("get_name", pure_virtual(&EMAN::Averager::get_name))
        .def("set_params", &EMAN::Averager::set_params, &EMAN_Averager_Wrapper::default_set_params)
        .def("get_param_types", &EMAN::Averager::get_param_types, &EMAN_Averager_Wrapper::default_get_param_types)
    ;

    class_< EMAN::Filter, boost::noncopyable, EMAN_Filter_Wrapper >("Filter", init<  >())
        .def("process", &EMAN::Filter::process, &EMAN_Filter_Wrapper::default_process)
        .def("process_list", &EMAN::Filter::process_list, &EMAN_Filter_Wrapper::default_process_list)
        .def("get_params", &EMAN::Filter::get_params, &EMAN_Filter_Wrapper::default_get_params)
        .def("set_params", &EMAN::Filter::set_params, &EMAN_Filter_Wrapper::default_set_params)
        .def("get_param_types", pure_virtual(&EMAN::Filter::get_param_types))
        .def("get_name", pure_virtual(&EMAN::Filter::get_name))
    ;

    class_< EMAN::Projector, boost::noncopyable, EMAN_Projector_Wrapper >("Projector", init<  >())
        .def("project3d", pure_virtual(&EMAN::Projector::project3d), return_value_policy< manage_new_object >())
        .def("get_name", pure_virtual(&EMAN::Projector::get_name))
        .def("get_params", &EMAN::Projector::get_params, &EMAN_Projector_Wrapper::default_get_params)
        .def("get_param_types", &EMAN::Projector::get_param_types, &EMAN_Projector_Wrapper::default_get_param_types)
        .def("set_params", &EMAN::Projector::set_params)
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
        .def("ift_slice", &EMAN::EMData::ift_slice, return_value_policy< manage_new_object >())
        .def("gimme_fft", &EMAN::EMData::gimme_fft)
        .def("normalize", &EMAN::EMData::normalize)
        .def("mask_normalize", &EMAN::EMData::mask_normalize, EMAN_EMData_mask_normalize_overloads_1_2())
        .def("edge_normalize", &EMAN::EMData::edge_normalize, EMAN_EMData_edge_normalize_overloads_0_1())
        .def("row_normalize", &EMAN::EMData::row_normalize)
        .def("normalize_max", &EMAN::EMData::normalize_max)
        .def("normalize_slice", &EMAN::EMData::normalize_slice)
        .def("normalize_to", &EMAN::EMData::normalize_to, EMAN_EMData_normalize_to_overloads_1_5())
        .def("least_square_normalize_to", &EMAN::EMData::least_square_normalize_to, EMAN_EMData_least_square_normalize_to_overloads_1_3())
        .def("ri2ap", &EMAN::EMData::ri2ap)
        .def("ap2ri", &EMAN::EMData::ap2ri)
        .def("to_corner", &EMAN::EMData::to_corner)
        .def("to_mass_center", &EMAN::EMData::to_mass_center, EMAN_EMData_to_mass_center_overloads_0_1())
        .def("rotate_x", &EMAN::EMData::rotate_x)
        .def("rotate_180", &EMAN::EMData::rotate_180)
        .def("fast_translate", &EMAN::EMData::fast_translate, EMAN_EMData_fast_translate_overloads_0_1())
        .def("rotate_translate", &EMAN::EMData::rotate_translate, EMAN_EMData_rotate_translate_overloads_0_5())
        .def("fast_rotate_translate", &EMAN::EMData::fast_rotate_translate, EMAN_EMData_fast_rotate_translate_overloads_0_1())
        .def("dot_rotate_translate", &EMAN::EMData::dot_rotate_translate)
        .def("vertical_flip", &EMAN::EMData::vertical_flip)
        .def("horizontal_flip", &EMAN::EMData::horizontal_flip)
        .def("little_big_dot", &EMAN::EMData::little_big_dot, return_value_policy< manage_new_object >(), EMAN_EMData_little_big_dot_overloads_1_2())
        .def("do_radon", &EMAN::EMData::do_radon, return_value_policy< manage_new_object >())
        .def("calc_ccf", &EMAN::EMData::calc_ccf, return_value_policy< manage_new_object >(), EMAN_EMData_calc_ccf_overloads_1_3())
        .def("make_rotational_footprint", &EMAN::EMData::make_rotational_footprint, return_value_policy< manage_new_object >(), EMAN_EMData_make_rotational_footprint_overloads_0_2())
        .def("calc_ccfx", &EMAN::EMData::calc_ccfx, return_value_policy< manage_new_object >(), EMAN_EMData_calc_ccfx_overloads_1_4())
        .def("calc_mutual_correlation", &EMAN::EMData::calc_mutual_correlation, return_value_policy< manage_new_object >(), EMAN_EMData_calc_mutual_correlation_overloads_1_3())
        .def("unwrap", &EMAN::EMData::unwrap, return_value_policy< manage_new_object >(), EMAN_EMData_unwrap_overloads_0_6())
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
        .def("calc_fourier_shell_correlation", &EMAN::EMData::calc_fourier_shell_correlation)
        .def("calc_hist", &EMAN::EMData::calc_hist, EMAN_EMData_calc_hist_overloads_1_4())
        .def("calc_az_dist", &EMAN::EMData::calc_az_dist)
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
        .def("get_edge_mean", &EMAN::EMData::get_edge_mean)
        .def("get_circle_mean", &EMAN::EMData::get_circle_mean)
        .def("radial_average", &EMAN::EMData::radial_average)
        .def("radial_subtract", &EMAN::EMData::radial_subtract)
        .def("sub_noise", &EMAN::EMData::sub_noise)
        .def("setup_insert_slice", &EMAN::EMData::setup_insert_slice)
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

    scope* EMAN_Ctf_scope = new scope(
    class_< EMAN::Ctf, boost::noncopyable, EMAN_Ctf_Wrapper >("Ctf", init<  >())
        .def("cmp", pure_virtual(&EMAN::Ctf::cmp))
        .def("from_string", pure_virtual(&EMAN::Ctf::from_string))
        .def("to_string", pure_virtual(&EMAN::Ctf::to_string))
        .def("from_dict", pure_virtual(&EMAN::Ctf::from_dict))
        .def("to_dict", pure_virtual(&EMAN::Ctf::to_dict))
    );

    enum_< UniqueInt<0> >("unnamed")
        .value("CTFOS", EMAN::Ctf::CTFOS)
        .export_values()
    ;


    enum_< EMAN::Ctf::CtfMapType >("CtfMapType")
        .value("CTF_MAP_BACKGROUND", EMAN::Ctf::CTF_MAP_BACKGROUND)
        .value("CTF_MAP_B_FACTOR", EMAN::Ctf::CTF_MAP_B_FACTOR)
        .value("CTF_MAP_CTF_NO_B", EMAN::Ctf::CTF_MAP_CTF_NO_B)
        .value("CTF_MAP_AMP_NO_B", EMAN::Ctf::CTF_MAP_AMP_NO_B)
        .value("CTF_MAP_AMP", EMAN::Ctf::CTF_MAP_AMP)
        .value("CTF_MAP_WIENER_CTF_CORRECTION", EMAN::Ctf::CTF_MAP_WIENER_CTF_CORRECTION)
        .value("CTF_MAP_SNR_SIGN", EMAN::Ctf::CTF_MAP_SNR_SIGN)
        .value("CTF_MAP_SIGN", EMAN::Ctf::CTF_MAP_SIGN)
        .value("CTF_MAP_WIENER_FILTER", EMAN::Ctf::CTF_MAP_WIENER_FILTER)
        .value("CTF_MAP_SNR", EMAN::Ctf::CTF_MAP_SNR)
        .value("CTF_MAP_CTF", EMAN::Ctf::CTF_MAP_CTF)
    ;


    enum_< EMAN::Ctf::CtfCurveType >("CtfCurveType")
        .value("CTF_CURVE_WIENER_CTF_CORRECTION1", EMAN::Ctf::CTF_CURVE_WIENER_CTF_CORRECTION1)
        .value("CTF_CURVE_WIENER_CTF_CORRECTION2", EMAN::Ctf::CTF_CURVE_WIENER_CTF_CORRECTION2)
        .value("CTF_CURVE_SNR_WIENER", EMAN::Ctf::CTF_CURVE_SNR_WIENER)
        .value("CTF_CURVE_ABS_AMP_S", EMAN::Ctf::CTF_CURVE_ABS_AMP_S)
        .value("CTF_CURVE_ABS_SNR", EMAN::Ctf::CTF_CURVE_ABS_SNR)
        .value("CTF_CURVE_RELATIVE_SNR", EMAN::Ctf::CTF_CURVE_RELATIVE_SNR)
        .value("CTF_CURVE_NOISE_S", EMAN::Ctf::CTF_CURVE_NOISE_S)
        .value("CTF_CURVE_TOTAL_CURVE", EMAN::Ctf::CTF_CURVE_TOTAL_CURVE)
        .value("CTF_CURVE_AMP_S", EMAN::Ctf::CTF_CURVE_AMP_S)
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
        .def_readwrite("ctfmaptype", &EMAN::SimpleCtf::ctfmaptype)
        .def_readwrite("astig_amp", &EMAN::SimpleCtf::astig_amp)
        .def_readwrite("astig_ang", &EMAN::SimpleCtf::astig_ang)
        .def_readwrite("drift_amp", &EMAN::SimpleCtf::drift_amp)
        .def_readwrite("drift_ang", &EMAN::SimpleCtf::drift_ang)
        .def("cmp", (bool (EMAN::SimpleCtf::*)() const)&EMAN::SimpleCtf::cmp, (bool (EMAN_SimpleCtf_Wrapper::*)() const)&EMAN_SimpleCtf_Wrapper::default_cmp)
        .def("from_string", (int (EMAN::SimpleCtf::*)(std::string) )&EMAN::SimpleCtf::from_string, (int (EMAN_SimpleCtf_Wrapper::*)(std::string))&EMAN_SimpleCtf_Wrapper::default_from_string)
        .def("to_string", (std::string (EMAN::SimpleCtf::*)() const)&EMAN::SimpleCtf::to_string, (std::string (EMAN_SimpleCtf_Wrapper::*)() const)&EMAN_SimpleCtf_Wrapper::default_to_string)
        .def("from_dict", (int (EMAN::SimpleCtf::*)(std::map<std::basic_string<char, std::char_traits<char>, std::allocator<char> >,EMAN::EMObject,std::less<std::basic_string<char, std::char_traits<char>, std::allocator<char> > >,std::allocator<std::pair<const std::basic_string<char, std::char_traits<char>, std::allocator<char> >, EMAN::EMObject> > >&) )&EMAN::SimpleCtf::from_dict, (int (EMAN_SimpleCtf_Wrapper::*)(std::map<std::basic_string<char, std::char_traits<char>, std::allocator<char> >,EMAN::EMObject,std::less<std::basic_string<char, std::char_traits<char>, std::allocator<char> > >,std::allocator<std::pair<const std::basic_string<char, std::char_traits<char>, std::allocator<char> >, EMAN::EMObject> > >&))&EMAN_SimpleCtf_Wrapper::default_from_dict)
        .def("to_dict", (int (EMAN::SimpleCtf::*)(std::map<std::basic_string<char, std::char_traits<char>, std::allocator<char> >,EMAN::EMObject,std::less<std::basic_string<char, std::char_traits<char>, std::allocator<char> > >,std::allocator<std::pair<const std::basic_string<char, std::char_traits<char>, std::allocator<char> >, EMAN::EMObject> > >&) const)&EMAN::SimpleCtf::to_dict, (int (EMAN_SimpleCtf_Wrapper::*)(std::map<std::basic_string<char, std::char_traits<char>, std::allocator<char> >,EMAN::EMObject,std::less<std::basic_string<char, std::char_traits<char>, std::allocator<char> > >,std::allocator<std::pair<const std::basic_string<char, std::char_traits<char>, std::allocator<char> >, EMAN::EMObject> > >&) const)&EMAN_SimpleCtf_Wrapper::default_to_dict)
        .def("get_maptype", &EMAN::SimpleCtf::get_maptype)
        .def("compute_map", &EMAN::SimpleCtf::compute_map, EMAN_SimpleCtf_compute_map_overloads_1_2())
        .def("is_changed", &EMAN::SimpleCtf::is_changed)
        .def("is_set_properly", &EMAN::SimpleCtf::is_set_properly)
        .def("compute_curve", &EMAN::SimpleCtf::compute_curve, EMAN_SimpleCtf_compute_curve_overloads_2_3())
        .staticmethod("compute_curve")
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

    class_< EMAN::Vec3f >("Vec3f", init<  >())
        .def(init< const EMAN::Vec3f& >())
        .def(init< float, float, float >())
        .def(init< const std::vector<float,std::allocator<float> >& >())
        .def("normalize", &EMAN::Vec3f::normalize)
        .def("length", &EMAN::Vec3f::length)
        .def("dot", &EMAN::Vec3f::dot)
        .def("cross", &EMAN::Vec3f::cross)
        .def("negate", &EMAN::Vec3f::negate, return_internal_reference< 1 >())
        .def("get_value", &EMAN::Vec3f::get_value)
        .def("set_value", (void (EMAN::Vec3f::*)(const std::vector<float,std::allocator<float> >&) )&EMAN::Vec3f::set_value)
        .def("set_value", (void (EMAN::Vec3f::*)(float, float, float) )&EMAN::Vec3f::set_value)
        .def( other< float >() * self )
        .def( self - self )
        .def( self * other< float >() )
        .def( other< float >() / self )
        .def( self + self )
        .def( self != self )
        .def( other< EMAN::Matrix3f >() * self )
        .def( self * other< EMAN::Matrix3f >() )
        .def( self / other< float >() )
        .def( self == self )
        .def( self += self )
        .def( self -= self )
        .def( self *= other< float >() )
        .def( self /= other< float >() )
    ;

    class_< EMAN::Matrix3f >("Matrix3f", init<  >())
        .def(init< const EMAN::Matrix3f& >())
        .def(init< float, float, float, float, float, float, float, float, float >())
        .def(init< const std::vector<float,std::allocator<float> >& >())
        .def("make_identity", &EMAN::Matrix3f::make_identity)
        .def("mult_right", &EMAN::Matrix3f::mult_right, return_internal_reference< 1 >())
        .def("mult_left", &EMAN::Matrix3f::mult_left, return_internal_reference< 1 >())
        .def("set_value", &EMAN::Matrix3f::set_value)
        .def("get_value", &EMAN::Matrix3f::get_value)
        .def("inverse", &EMAN::Matrix3f::inverse, return_internal_reference< 1 >())
        .def("transpose", &EMAN::Matrix3f::transpose, return_internal_reference< 1 >())
        .def("create_inverse", &EMAN::Matrix3f::create_inverse)
        .def( self == self )
        .def( self / self )
        .def( self / other< float >() )
        .def( self * other< float >() )
        .def( other< float >() + self )
        .def( other< float >() - self )
        .def( other< float >() * self )
        .def( other< float >() / self )
        .def( self - other< float >() )
        .def( self * other< EMAN::Vec3f >() )
        .def( other< EMAN::Vec3f >() * self )
        .def( self != self )
        .def( self * self )
        .def( self - self )
        .def( self + self )
        .def( self + other< float >() )
        .def( self += other< float >() )
        .def( self -= other< float >() )
        .def( self *= other< float >() )
        .def( self /= other< float >() )
        .def( self += self )
        .def( self -= self )
        .def( self *= self )
        .def( self /= self )
    ;

    class_< EMAN::Matrix4f >("Matrix4f", init<  >())
        .def(init< const EMAN::Matrix4f& >())
        .def(init< const std::vector<float,std::allocator<float> >& >())
        .def(init< const EMAN::Matrix3f& >())
        .def("mult_right", &EMAN::Matrix4f::mult_right, return_internal_reference< 1 >())
        .def("mult_left", &EMAN::Matrix4f::mult_left, return_internal_reference< 1 >())
        .def("make_identity", &EMAN::Matrix4f::make_identity)
        .def("set_value", &EMAN::Matrix4f::set_value)
        .def("get_value", &EMAN::Matrix4f::get_value)
        .def("get_matrix3", &EMAN::Matrix4f::get_matrix3)
        .def("inverse", &EMAN::Matrix4f::inverse, return_internal_reference< 1 >())
        .def("transpose", &EMAN::Matrix4f::transpose, return_internal_reference< 1 >())
        .def("create_inverse", &EMAN::Matrix4f::create_inverse)
        .def( self + other< float >() )
        .def( other< float >() / self )
        .def( other< float >() - self )
        .def( self + self )
        .def( self / other< float >() )
        .def( self * other< float >() )
        .def( self / self )
        .def( self != self )
        .def( self - self )
        .def( self * self )
        .def( other< float >() + self )
        .def( other< float >() * self )
        .def( self - other< float >() )
        .def( self == self )
        .def( self += other< float >() )
        .def( self -= other< float >() )
        .def( self *= other< float >() )
        .def( self /= other< float >() )
        .def( self += self )
        .def( self -= self )
        .def( self *= self )
        .def( self /= self )
    ;

    class_< EMAN::Quaternion >("Quaternion", init<  >())
        .def(init< const EMAN::Quaternion& >())
        .def(init< float, float, float, float >())
        .def(init< float, const EMAN::Vec3f& >())
        .def(init< const EMAN::Vec3f&, float >())
        .def(init< const EMAN::Matrix3f& >())
        .def(init< const EMAN::Matrix4f& >())
        .def("norm", &EMAN::Quaternion::norm)
        .def("conj", &EMAN::Quaternion::conj)
        .def("abs", &EMAN::Quaternion::abs)
        .def("normalize", &EMAN::Quaternion::normalize)
        .def("inverse", &EMAN::Quaternion::inverse, return_internal_reference< 1 >())
        .def("create_inverse", &EMAN::Quaternion::create_inverse)
        .def("rotate", &EMAN::Quaternion::rotate)
        .def("interpolate", (void (EMAN::Quaternion::*)(const EMAN::Quaternion&, float) )&EMAN::Quaternion::interpolate)
        .def("interpolate", (void (EMAN::Quaternion::*)(const EMAN::Quaternion&, const EMAN::Quaternion&, float) )&EMAN::Quaternion::interpolate)
        .def("to_angle", &EMAN::Quaternion::to_angle)
        .def("to_axis", &EMAN::Quaternion::to_axis)
        .def("to_matrix3", &EMAN::Quaternion::to_matrix3)
        .def("to_matrix4", &EMAN::Quaternion::to_matrix4)
        .def("real", &EMAN::Quaternion::real)
        .def("unreal", &EMAN::Quaternion::unreal)
        .def("get_value", &EMAN::Quaternion::get_value)
        .def( self + self )
        .def( self - self )
        .def( self * self )
        .def( self * other< float >() )
        .def( self != self )
        .def( self == self )
        .def( self / self )
        .def( other< float >() * self )
        .def( self += self )
        .def( self -= self )
        .def( self *= self )
        .def( self *= other< float >() )
        .def( self /= self )
        .def( self /= other< float >() )
    ;

    scope* EMAN_Rotation_scope = new scope(
    class_< EMAN::Rotation >("Rotation", init<  >())
        .def(init< const EMAN::Rotation& >())
        .def(init< float, float, float, EMAN::Rotation::Type >())
        .def(init< float, float, float, float, EMAN::Rotation::Type >())
        .def(init< const EMAN::Quaternion& >())
        .def(init< const EMAN::Matrix3f& >())
        .def_readonly("ERR_LIMIT", &EMAN::Rotation::ERR_LIMIT)
        .def("inverse", &EMAN::Rotation::inverse, return_internal_reference< 1 >())
        .def("create_inverse", &EMAN::Rotation::create_inverse)
        .def("diff", &EMAN::Rotation::diff)
        .def("rotate_from_left", &EMAN::Rotation::rotate_from_left, return_internal_reference< 1 >())
        .def("set_sym", &EMAN::Rotation::set_sym)
        .def("get_max_nsym", &EMAN::Rotation::get_max_nsym)
        .def("get_sym", &EMAN::Rotation::get_sym)
        .def("set_angle", (void (EMAN::Rotation::*)(float, float, float, EMAN::Rotation::Type) )&EMAN::Rotation::set_angle)
        .def("set_angle", (void (EMAN::Rotation::*)(float, float, float, float, EMAN::Rotation::Type) )&EMAN::Rotation::set_angle)
        .def("set_angle", (void (EMAN::Rotation::*)(const EMAN::Matrix3f&) )&EMAN::Rotation::set_angle)
        .def("set_angle", (void (EMAN::Rotation::*)(const EMAN::Quaternion&) )&EMAN::Rotation::set_angle)
        .def("is_valid", &EMAN::Rotation::is_valid)
        .def("rectify", &EMAN::Rotation::rectify)
        .def("eman_alt", &EMAN::Rotation::eman_alt)
        .def("eman_az", &EMAN::Rotation::eman_az)
        .def("eman_phi", &EMAN::Rotation::eman_phi)
        .def("mrc_theta", &EMAN::Rotation::mrc_theta)
        .def("mrc_phi", &EMAN::Rotation::mrc_phi)
        .def("mrc_omega", &EMAN::Rotation::mrc_omega)
        .def("imagic_alpha", &EMAN::Rotation::imagic_alpha)
        .def("imagic_beta", &EMAN::Rotation::imagic_beta)
        .def("imagic_gamma", &EMAN::Rotation::imagic_gamma)
        .def("spider_phi", &EMAN::Rotation::spider_phi)
        .def("spider_theta", &EMAN::Rotation::spider_theta)
        .def("spider_psi", &EMAN::Rotation::spider_psi)
        .def("get_spin_axis", (std::vector<float,std::allocator<float> > (EMAN::Rotation::*)() const)&EMAN::Rotation::get_spin_axis)
        .def("get_spin_axis", (void (EMAN::Rotation::*)(float*, float*, float*, float*) const)&EMAN::Rotation::get_spin_axis)
        .def("get_sgi", (std::vector<float,std::allocator<float> > (EMAN::Rotation::*)() const)&EMAN::Rotation::get_sgi)
        .def("get_sgi", (void (EMAN::Rotation::*)(float*, float*, float*, float*) const)&EMAN::Rotation::get_sgi)
        .def("get_quaternion", &EMAN::Rotation::get_quaternion)
        .def("get_matrix3", &EMAN::Rotation::get_matrix3)
        .def("get_matrix4", &EMAN::Rotation::get_matrix4)
        .def( self == self )
        .def( self / self )
        .def( self != self )
        .def( self * self )
        .def( self *= self )
        .def( self /= self )
    );

    enum_< EMAN::Rotation::Type >("Type")
        .value("MATRIX", EMAN::Rotation::MATRIX)
        .value("UNKNOWN", EMAN::Rotation::UNKNOWN)
        .value("IMAGIC", EMAN::Rotation::IMAGIC)
        .value("SPIDER", EMAN::Rotation::SPIDER)
        .value("QUATERNION", EMAN::Rotation::QUATERNION)
        .value("SGIROT", EMAN::Rotation::SGIROT)
        .value("MRC", EMAN::Rotation::MRC)
        .value("SPIN", EMAN::Rotation::SPIN)
        .value("EMAN", EMAN::Rotation::EMAN)
    ;

    delete EMAN_Rotation_scope;

    scope* EMAN_Transform_scope = new scope(
    class_< EMAN::Transform >("Transform", init<  >())
        .def(init< const EMAN::Transform& >())
        .def(init< const EMAN::Matrix4f& >())
        .def(init< const EMAN::Rotation& >())
        .def(init< const EMAN::Rotation&, const EMAN::Vec3f& >())
        .def("set_rotate_instance", (EMAN::Transform& (EMAN::Transform::*)(const EMAN::Rotation&) )&EMAN::Transform::set_rotate_instance, return_internal_reference< 1 >())
        .def("set_rotate_instance", (EMAN::Transform& (EMAN::Transform::*)(const EMAN::Matrix3f&) )&EMAN::Transform::set_rotate_instance, return_internal_reference< 1 >())
        .def("set_translate_instance", &EMAN::Transform::set_translate_instance, return_internal_reference< 1 >())
        .def("set_scale_instance", &EMAN::Transform::set_scale_instance, return_internal_reference< 1 >())
        .def("set_transform_instance", (EMAN::Transform& (EMAN::Transform::*)(const EMAN::Vec3f&, const EMAN::Rotation&, const EMAN::Vec3f&, const EMAN::Rotation&, const EMAN::Vec3f&) )&EMAN::Transform::set_transform_instance, return_internal_reference< 1 >())
        .def("set_transform_instance", (EMAN::Transform& (EMAN::Transform::*)(const EMAN::Vec3f&, const EMAN::Rotation&, const EMAN::Vec3f&) )&EMAN::Transform::set_transform_instance, return_internal_reference< 1 >())
        .def("set_center", &EMAN::Transform::set_center, return_internal_reference< 1 >())
        .def("set_matrix", &EMAN::Transform::set_matrix, return_internal_reference< 1 >())
        .def("set_post_translate", &EMAN::Transform::set_post_translate, return_internal_reference< 1 >())
        .def("inverse", &EMAN::Transform::inverse, return_internal_reference< 1 >())
        .def("create_inverse", &EMAN::Transform::create_inverse)
        .def("transpose", &EMAN::Transform::transpose, return_internal_reference< 1 >())
        .def("post_concatenate", &EMAN::Transform::post_concatenate, return_internal_reference< 1 >())
        .def("pre_concatenate", &EMAN::Transform::pre_concatenate, return_internal_reference< 1 >())
        .def("translate", &EMAN::Transform::translate, return_internal_reference< 1 >())
        .def("rotate", &EMAN::Transform::rotate, return_internal_reference< 1 >())
        .def("rotate_center", &EMAN::Transform::rotate_center, return_internal_reference< 1 >())
        .def("rotate_scale", &EMAN::Transform::rotate_scale, return_internal_reference< 1 >())
        .def("pre_translate_rotate", &EMAN::Transform::pre_translate_rotate, return_internal_reference< 1 >())
        .def("post_translate_rotate", &EMAN::Transform::post_translate_rotate, return_internal_reference< 1 >())
        .def("scale", &EMAN::Transform::scale, return_internal_reference< 1 >())
        .def("transform", &EMAN::Transform::transform)
        .def("inverse_transform", &EMAN::Transform::inverse_transform)
        .def("get_rotation", &EMAN::Transform::get_rotation)
        .def("get_scale", &EMAN::Transform::get_scale)
        .def("get_center", &EMAN::Transform::get_center)
        .def("get_matrix", &EMAN::Transform::get_matrix)
        .def("get_pre_translate", &EMAN::Transform::get_pre_translate)
        .def("get_post_translate", &EMAN::Transform::get_post_translate)
        .def("get_type", &EMAN::Transform::get_type)
        .def("interpolate", &EMAN::Transform::interpolate)
        .staticmethod("interpolate")
        .def( self / other< float >() )
        .def( other< float >() / self )
        .def( other< float >() * self )
        .def( self / self )
        .def( self * self )
        .def( self * other< float >() )
        .def( self + self )
        .def( self - self )
        .def( self += self )
        .def( self -= self )
        .def( self *= self )
        .def( self *= other< float >() )
        .def( self /= self )
        .def( self /= other< float >() )
    );

    enum_< EMAN::Transform::TransformType >("TransformType")
        .value("SCALE", EMAN::Transform::SCALE)
        .value("UNIFORM_SCALE", EMAN::Transform::UNIFORM_SCALE)
        .value("TRANSFORM", EMAN::Transform::TRANSFORM)
        .value("ROTATION", EMAN::Transform::ROTATION)
        .value("TRANSLATION", EMAN::Transform::TRANSLATION)
        .value("IDENTITY", EMAN::Transform::IDENTITY)
    ;

    delete EMAN_Transform_scope;

    scope* EMAN_XYData_scope = new scope(
    class_< EMAN::XYData >("XYData", init<  >())
        .def(init< const EMAN::XYData& >())
        .def("read_file", &EMAN::XYData::read_file)
        .def("write_file", &EMAN::XYData::write_file)
        .def("calc_correlation", &EMAN::XYData::calc_correlation)
        .def("get_yatx", &EMAN::XYData::get_yatx)
        .def("get_x", &EMAN::XYData::get_x)
        .def("get_y", &EMAN::XYData::get_y)
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

    class_< EMAN::Size >("Size", init<  >())
        .def(init< const EMAN::Size& >())
        .def(init< int, int >())
        .def(init< int, int, int >())
        .def_readwrite("xsize", &EMAN::Size::xsize)
        .def_readwrite("ysize", &EMAN::Size::ysize)
        .def_readwrite("zsize", &EMAN::Size::zsize)
        .def("get_ndim", &EMAN::Size::get_ndim)
    ;

    class_< EMAN::Point<int> >("EMAN_Point_int", init<  >())
        .def(init< const EMAN::Point<int>& >())
        .def(init< int, int >())
        .def(init< int, int, int >())
        .def_readwrite("x", &EMAN::Point<int>::x)
        .def_readwrite("y", &EMAN::Point<int>::y)
        .def_readwrite("z", &EMAN::Point<int>::z)
        .def("get_ndim", &EMAN::Point<int>::get_ndim)
    ;

    class_< EMAN::Point<float> >("EMAN_Point_float", init<  >())
        .def(init< const EMAN::Point<float>& >())
        .def(init< float, float >())
        .def(init< float, float, float >())
        .def_readwrite("x", &EMAN::Point<float>::x)
        .def_readwrite("y", &EMAN::Point<float>::y)
        .def_readwrite("z", &EMAN::Point<float>::z)
        .def("get_ndim", &EMAN::Point<float>::get_ndim)
    ;

    class_< EMAN::Region >("Region", init<  >())
        .def(init< const EMAN::Region& >())
        .def(init< float, float, int, int >())
        .def(init< float, float, float, int, int, int >())
        .def(init< const EMAN::Point<float>&, const EMAN::Size& >())
        .def_readwrite("origin", &EMAN::Region::origin)
        .def_readwrite("size", &EMAN::Region::size)
        .def("inside_region", (bool (EMAN::Region::*)() const)&EMAN::Region::inside_region)
        .def("inside_region", (bool (EMAN::Region::*)(const EMAN::Size&) const)&EMAN::Region::inside_region)
        .def("inside_region", (bool (EMAN::Region::*)(float, float) const)&EMAN::Region::inside_region)
        .def("inside_region", (bool (EMAN::Region::*)(float, float, float) const)&EMAN::Region::inside_region)
        .def("get_ndim", &EMAN::Region::get_ndim)
        .def("get_string", &EMAN::Region::get_string)
    ;

    class_< EMAN::Util >("Util", init<  >())
        .def(init< const EMAN::Util& >())
        .def("ap2ri", &EMAN::Util::ap2ri)
        .def("file_lock_wait", &EMAN::Util::file_lock_wait)
        .def("generate_machine_stamp", &EMAN::Util::generate_machine_stamp)
        .def("check_file_by_magic", &EMAN::Util::check_file_by_magic)
        .def("flip_image", &EMAN::Util::flip_image)
        .def("is_sub_string", &EMAN::Util::is_sub_string)
        .def("get_filename_by_ext", &EMAN::Util::get_filename_by_ext)
        .def("calc_least_square_fit", &EMAN::Util::calc_least_square_fit)
        .def("save_data_to_file", (void (*)(const std::vector<float,std::allocator<float> >&, const std::vector<float,std::allocator<float> >&, std::string))&EMAN::Util::save_data_to_file)
        .def("save_data_to_file", (void (*)(float, float, const std::vector<float,std::allocator<float> >&, std::string))&EMAN::Util::save_data_to_file)
        .def("save_data_to_file", (void (*)(float, float, float*, int, std::string))&EMAN::Util::save_data_to_file)
        .def("get_frand", &EMAN::Util::get_frand)
        .def("get_gaussian_rand", &EMAN::Util::get_gaussian_rand)
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
        .staticmethod("square")
        .staticmethod("calc_best_fft_size")
        .staticmethod("square_sum")
        .staticmethod("angle_sub_pi")
        .staticmethod("trilinear_interpolate")
        .staticmethod("file_lock_wait")
        .staticmethod("min")
        .staticmethod("angle_sub_2pi")
        .staticmethod("get_gaussian_rand")
        .staticmethod("find_max")
        .staticmethod("max")
        .staticmethod("fast_floor")
        .staticmethod("ap2ri")
        .staticmethod("save_data_to_file")
        .staticmethod("check_file_by_magic")
        .staticmethod("generate_machine_stamp")
        .staticmethod("agauss")
        .staticmethod("bilinear_interpolate")
        .staticmethod("calc_least_square_fit")
        .staticmethod("find_min_and_max")
        .staticmethod("get_filename_by_ext")
        .staticmethod("get_frand")
        .staticmethod("goodf")
        .staticmethod("hypot3")
        .staticmethod("is_sub_string")
        .staticmethod("flip_image")
        .staticmethod("round")
    ;

}

