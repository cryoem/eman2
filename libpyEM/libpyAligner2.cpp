/*
 * Copyright (c) 2000-2006 Baylor College of Medicine
 *
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * */

#ifdef _WIN32
	#pragma warning(disable:4819)
#endif	//_WIN32

// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <aligner.h>
#include <ctf.h>
#include <emdata.h>
#include <emobject.h>
#include <xydata.h>

#include "emdata_pickle.h"

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct EMAN_Aligner_Wrapper: EMAN::Aligner
{
    EMAN_Aligner_Wrapper(PyObject* py_self_, const EMAN::Aligner& p0):
        EMAN::Aligner(p0), py_self(py_self_) {}

    EMAN_Aligner_Wrapper(PyObject* py_self_):
        EMAN::Aligner(), py_self(py_self_) {}

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
    EMAN_Ctf_Wrapper(PyObject* py_self_, const EMAN::Ctf& p0):
        EMAN::Ctf(p0), py_self(py_self_) {}

    EMAN_Ctf_Wrapper(PyObject* py_self_):
        EMAN::Ctf(), py_self(py_self_) {}

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

    void from_vector(const std::vector<float,std::allocator<float> >& p0) {
        call_method< void >(py_self, "from_vector", p0);
    }

    std::vector<float,std::allocator<float> > to_vector() const {
        return call_method< std::vector<float,std::allocator<float> > >(py_self, "to_vector");
    }

    std::vector<float,std::allocator<float> > compute_1d(int p0,float p1, EMAN::Ctf::CtfType p2, EMAN::XYData* p3) {
        return call_method< std::vector<float,std::allocator<float> > >(py_self, "compute_1d", p0, p1, p2, p3);
    }

    std::vector<float,std::allocator<float> > compute_1d_fromimage(int p0, float p1, EMAN::EMData *p2) {
        return call_method< std::vector<float,std::allocator<float> > >(py_self, "compute_1d_fromimage", p0, p1, p2);
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

    float zero(int n) const {
        return call_method< float >(py_self, "zero", n);
    }

    
    PyObject* py_self;
};

struct EMAN_EMAN1Ctf_Wrapper: EMAN::EMAN1Ctf
{
    EMAN_EMAN1Ctf_Wrapper(PyObject* py_self_, const EMAN::EMAN1Ctf& p0):
        EMAN::EMAN1Ctf(p0), py_self(py_self_) {}

    EMAN_EMAN1Ctf_Wrapper(PyObject* py_self_):
        EMAN::EMAN1Ctf(), py_self(py_self_) {}

    std::vector<float,std::allocator<float> > compute_1d(int p0, float p1,EMAN::Ctf::CtfType p2, EMAN::XYData* p3) {
        return call_method< std::vector<float,std::allocator<float> > >(py_self, "compute_1d", p0, p1, p2, p3);
    }

    std::vector<float,std::allocator<float> > default_compute_1d_2(int p0, float p1, EMAN::Ctf::CtfType p2) {
        return EMAN::EMAN1Ctf::compute_1d(p0, p1, p2);
    }

    std::vector<float,std::allocator<float> > default_compute_1d_3(int p0, float p1, EMAN::Ctf::CtfType p2, EMAN::XYData* p3) {
        return EMAN::EMAN1Ctf::compute_1d(p0, p1, p2, p3);
    }

    std::vector<float,std::allocator<float> > compute_1d_fromimage(int p0, float p1, EMAN::EMData* p2) {
        return call_method< std::vector<float,std::allocator<float> > >(py_self, "compute_1d_fromimage", p0, p1, p2);
    }

    std::vector<float,std::allocator<float> > default_compute_1d_fromimage(int p0, float p1, EMAN::EMData* p2) {
        return EMAN::EMAN1Ctf::compute_1d_fromimage(p0, p1, p2);
    }

    
    void compute_2d_real(EMAN::EMData* p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        call_method< void >(py_self, "compute_2d_real", p0, p1, p2);
    }

    void default_compute_2d_real_2(EMAN::EMData* p0, EMAN::Ctf::CtfType p1) {
        EMAN::EMAN1Ctf::compute_2d_real(p0, p1);
    }

    void default_compute_2d_real_3(EMAN::EMData* p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        EMAN::EMAN1Ctf::compute_2d_real(p0, p1, p2);
    }

    void compute_2d_complex(EMAN::EMData* p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        call_method< void >(py_self, "compute_2d_complex", p0, p1, p2);
    }

    void default_compute_2d_complex_2(EMAN::EMData* p0, EMAN::Ctf::CtfType p1) {
        EMAN::EMAN1Ctf::compute_2d_complex(p0, p1);
    }

    void default_compute_2d_complex_3(EMAN::EMData* p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        EMAN::EMAN1Ctf::compute_2d_complex(p0, p1, p2);
    }

    int from_string(const std::string& p0) {
        return call_method< int >(py_self, "from_string", p0);
    }

    int default_from_string(const std::string& p0) {
        return EMAN::EMAN1Ctf::from_string(p0);
    }

    std::string to_string() const {
        return call_method< std::string >(py_self, "to_string");
    }

    std::string default_to_string() const {
        return EMAN::EMAN1Ctf::to_string();
    }

    void from_dict(const EMAN::Dict& p0) {
        call_method< void >(py_self, "from_dict", p0);
    }

    void default_from_dict(const EMAN::Dict& p0) {
        EMAN::EMAN1Ctf::from_dict(p0);
    }

    EMAN::Dict to_dict() const {
        return call_method< EMAN::Dict >(py_self, "to_dict");
    }

    EMAN::Dict default_to_dict() const {
        return EMAN::EMAN1Ctf::to_dict();
    }

    void from_vector(const std::vector<float,std::allocator<float> >& p0) {
        call_method< void >(py_self, "from_vector", p0);
    }

    void default_from_vector(const std::vector<float,std::allocator<float> >& p0) {
        EMAN::EMAN1Ctf::from_vector(p0);
    }

    std::vector<float,std::allocator<float> > to_vector() const {
        return call_method< std::vector<float,std::allocator<float> > >(py_self, "to_vector");
    }

    std::vector<float,std::allocator<float> > default_to_vector() const {
        return EMAN::EMAN1Ctf::to_vector();
    }

    void copy_from(const EMAN::Ctf* p0) {
        call_method< void >(py_self, "copy_from", p0);
    }

    void default_copy_from(const EMAN::Ctf* p0) {
        EMAN::EMAN1Ctf::copy_from(p0);
    }

    bool equal(const EMAN::Ctf* p0) const {
        return call_method< bool >(py_self, "equal", p0);
    }
    
    bool default_equal(const EMAN::Ctf* p0) const {
        return EMAN::EMAN1Ctf::equal(p0);
    }
    
    float zero(int p0) const {
        return call_method< float >(py_self, "zero", p0);
    }

    float default_zero(int n) const {
        return EMAN::EMAN1Ctf::zero(n);
    }

    PyObject* py_self;
};

struct EMAN_EMAN2Ctf_Wrapper: EMAN::EMAN2Ctf
{
    EMAN_EMAN2Ctf_Wrapper(PyObject* py_self_, const EMAN::EMAN2Ctf& p0):
        EMAN::EMAN2Ctf(p0), py_self(py_self_) {}

    EMAN_EMAN2Ctf_Wrapper(PyObject* py_self_):
        EMAN::EMAN2Ctf(), py_self(py_self_) {}

    std::vector<float,std::allocator<float> > compute_1d(int p0, float p1, EMAN::Ctf::CtfType p2, EMAN::XYData* p3) {
        return call_method< std::vector<float,std::allocator<float> > >(py_self, "compute_1d", p0, p1, p2, p3);
    }

    std::vector<float,std::allocator<float> > default_compute_1d_2(int p0, float p1, EMAN::Ctf::CtfType p2) {
        return EMAN::EMAN2Ctf::compute_1d(p0, p1, p2);
    }

    std::vector<float,std::allocator<float> > default_compute_1d_3(int p0, float p1, EMAN::Ctf::CtfType p2, EMAN::XYData* p3) {
        return EMAN::EMAN2Ctf::compute_1d(p0, p1, p2, p3);
    }

    std::vector<float,std::allocator<float> > compute_1d_fromimage(int p0, float p1, EMAN::Ctf::CtfType p2, EMAN::XYData* p3) {
        return call_method< std::vector<float,std::allocator<float> > >(py_self, "compute_1d_fromimage", p0, p1, p2, p3);
    }

    std::vector<float,std::allocator<float> > default_compute_1d_fromimage(int p0, float p1, EMAN::EMData* p2) {
        return EMAN::EMAN2Ctf::compute_1d_fromimage(p0, p1, p2);
    }

    void compute_2d_real(EMAN::EMData* p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        call_method< void >(py_self, "compute_2d_real", p0, p1, p2);
    }

    void default_compute_2d_real_2(EMAN::EMData* p0, EMAN::Ctf::CtfType p1) {
        EMAN::EMAN2Ctf::compute_2d_real(p0, p1);
    }

    void default_compute_2d_real_3(EMAN::EMData* p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        EMAN::EMAN2Ctf::compute_2d_real(p0, p1, p2);
    }

    void compute_2d_complex(EMAN::EMData* p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        call_method< void >(py_self, "compute_2d_complex", p0, p1, p2);
    }

    void default_compute_2d_complex_2(EMAN::EMData* p0, EMAN::Ctf::CtfType p1) {
        EMAN::EMAN2Ctf::compute_2d_complex(p0, p1);
    }

    void default_compute_2d_complex_3(EMAN::EMData* p0, EMAN::Ctf::CtfType p1, EMAN::XYData* p2) {
        EMAN::EMAN2Ctf::compute_2d_complex(p0, p1, p2);
    }

    int from_string(const std::string& p0) {
        return call_method< int >(py_self, "from_string", p0);
    }

    int default_from_string(const std::string& p0) {
        return EMAN::EMAN2Ctf::from_string(p0);
    }

    std::string to_string() const {
        return call_method< std::string >(py_self, "to_string");
    }

    std::string default_to_string() const {
        return EMAN::EMAN2Ctf::to_string();
    }

    void from_dict(const EMAN::Dict& p0) {
        call_method< void >(py_self, "from_dict", p0);
    }

    void default_from_dict(const EMAN::Dict& p0) {
        EMAN::EMAN2Ctf::from_dict(p0);
    }

    EMAN::Dict to_dict() const {
        return call_method< EMAN::Dict >(py_self, "to_dict");
    }

    EMAN::Dict default_to_dict() const {
        return EMAN::EMAN2Ctf::to_dict();
    }

    void from_vector(const std::vector<float,std::allocator<float> >& p0) {
        call_method< void >(py_self, "from_vector", p0);
    }

    void default_from_vector(const std::vector<float,std::allocator<float> >& p0) {
        EMAN::EMAN2Ctf::from_vector(p0);
    }

    std::vector<float,std::allocator<float> > to_vector() const {
        return call_method< std::vector<float,std::allocator<float> > >(py_self, "to_vector");
    }

    std::vector<float,std::allocator<float> > default_to_vector() const {
        return EMAN::EMAN2Ctf::to_vector();
    }

    void copy_from(const EMAN::Ctf* p0) {
        call_method< void >(py_self, "copy_from", p0);
    }

    void default_copy_from(const EMAN::Ctf* p0) {
        EMAN::EMAN2Ctf::copy_from(p0);
    }

    bool equal(const EMAN::Ctf* p0) const {
        return call_method< bool >(py_self, "equal", p0);
    }

    bool default_equal(const EMAN::Ctf* p0) const {
        return EMAN::EMAN2Ctf::equal(p0);
    }

    float zero(int p0) const {
        return call_method< float >(py_self, "zero", p0);
    }
    
    float default_zero(int n) const {
        return EMAN::EMAN2Ctf::zero(n);
    }

    
    PyObject* py_self;
};

}// namespace


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyAligner2)
{
    def("dump_aligners", &EMAN::dump_aligners);
    def("dump_aligners_list", &EMAN::dump_aligners_list);
    class_< EMAN::Aligner, boost::noncopyable, EMAN_Aligner_Wrapper >("__Aligner", init<  >())
        .def("align", pure_virtual((EMAN::EMData* (EMAN::Aligner::*)(EMAN::EMData*, EMAN::EMData*) const)&EMAN::Aligner::align), return_value_policy< manage_new_object >())
        .def("align", pure_virtual((EMAN::EMData* (EMAN::Aligner::*)(EMAN::EMData*, EMAN::EMData*, const std::string&, const EMAN::Dict&) const)&EMAN::Aligner::align), return_value_policy< manage_new_object >())
		.def("xform_align_nbest", &EMAN::Aligner::xform_align_nbest)
        .def("get_name", pure_virtual(&EMAN::Aligner::get_name))
        .def("get_desc", pure_virtual(&EMAN::Aligner::get_desc))
        .def("get_params", &EMAN::Aligner::get_params, &EMAN_Aligner_Wrapper::default_get_params)
        .def("set_params", &EMAN::Aligner::set_params, &EMAN_Aligner_Wrapper::default_set_params)
        .def("get_param_types", pure_virtual(&EMAN::Aligner::get_param_types))
    ;

#ifdef SPARX_USING_CUDA
    class_< EMAN::CUDA_Aligner, boost::noncopyable>("CUDA_Aligner", init<int>())
    	.def("finish", &EMAN::CUDA_Aligner::finish)
        .def("setup", &EMAN::CUDA_Aligner::setup)
        .def("insert_image", &EMAN::CUDA_Aligner::insert_image)
	.def("filter_stack", &EMAN::CUDA_Aligner::filter_stack)
	.def("sum_oe", &EMAN::CUDA_Aligner::sum_oe)
        .def("alignment_2d", &EMAN::CUDA_Aligner::alignment_2d)
        .def("ali2d_single_iter", &EMAN::CUDA_Aligner::ali2d_single_iter)
    ;

    class_< EMAN::CUDA_multiref_aligner, boost::noncopyable>("CUDA_multiref_aligner", init<int>())
   	.def("finish", &EMAN::CUDA_multiref_aligner::finish)
        .def("setup", &EMAN::CUDA_multiref_aligner::setup)
	.def("setup_params", &EMAN::CUDA_multiref_aligner::setup_params)
        .def("insert_image", &EMAN::CUDA_multiref_aligner::insert_image)
        .def("insert_ref_image", &EMAN::CUDA_multiref_aligner::insert_ref_image)
	.def("multiref_ali2d", &EMAN::CUDA_multiref_aligner::multiref_ali2d)
    ;

#endif


    class_< EMAN::Factory<EMAN::Aligner>, boost::noncopyable >("Aligners", no_init)
        .def("get", (EMAN::Aligner* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&))&EMAN::Factory<EMAN::Aligner>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Aligner* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&, const EMAN::Dict&))&EMAN::Factory<EMAN::Aligner>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Aligner>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

    scope* EMAN_Ctf_scope = new scope(
    class_< EMAN::Ctf, boost::noncopyable, EMAN_Ctf_Wrapper >("Ctf",
    		"Ctf is the base class for all CTF model.\n"
    		"Contrast transfer function (CTF) is the  function that\n"
    		"describes the transfer of information from the object"
    		"to the contrast observed in the image for electron microscopy.",
    		init<  >())
        .def("from_string", pure_virtual(&EMAN::Ctf::from_string))
        .def("to_string", pure_virtual(&EMAN::Ctf::to_string))
        .def("from_dict", pure_virtual(&EMAN::Ctf::from_dict))
        .def("to_dict", pure_virtual(&EMAN::Ctf::to_dict))
        .def("from_vector", pure_virtual(&EMAN::Ctf::from_vector))
        .def("to_vector", pure_virtual(&EMAN::Ctf::to_vector))
        .def("compute_1d", pure_virtual(&EMAN::Ctf::compute_1d))
        .def("compute_1d_fromimage", pure_virtual(&EMAN::Ctf::compute_1d_fromimage))
        .def("compute_2d_real", pure_virtual(&EMAN::Ctf::compute_2d_real))
        .def("compute_2d_complex", pure_virtual(&EMAN::Ctf::compute_2d_complex))
        .def("copy_from", pure_virtual(&EMAN::Ctf::copy_from))
        .def("equal", pure_virtual(&EMAN::Ctf::equal))
        .def("zero", pure_virtual(&EMAN::Ctf::zero))
    );

    enum_< EMAN::Ctf::CtfType >("CtfType")
        .value("CTF_AMP", EMAN::Ctf::CTF_AMP)
        .value("CTF_SIGN", EMAN::Ctf::CTF_SIGN)
        .value("CTF_BACKGROUND", EMAN::Ctf::CTF_BACKGROUND)
        .value("CTF_SNR", EMAN::Ctf::CTF_SNR)
        .value("CTF_SNR_SMOOTH", EMAN::Ctf::CTF_SNR_SMOOTH)
        .value("CTF_WIENER_FILTER", EMAN::Ctf::CTF_WIENER_FILTER)
        .value("CTF_TOTAL", EMAN::Ctf::CTF_TOTAL)
	.value("CTF_FITREF", EMAN::Ctf::CTF_FITREF)
	.value("CTF_POWEVAL", EMAN::Ctf::CTF_POWEVAL)
	.value("CTF_INTEN", EMAN::Ctf::CTF_INTEN)
	.value("CTF_NOISERATIO", EMAN::Ctf::CTF_NOISERATIO)
    ;


//    scope().attr("CTFOS") = (int)EMAN::Ctf::CTFOS;

    delete EMAN_Ctf_scope;

    class_< EMAN::EMAN1Ctf, bases< EMAN::Ctf > , EMAN_EMAN1Ctf_Wrapper >("EMAN1Ctf",
				"EMAN1Ctf is the CTF model used in EMAN1.",
				init<  >())
        .def(init< const EMAN::EMAN1Ctf& >())
        .def(init<const std::vector<float>& >())
        .def_pickle(EMAN1Ctf_pickle_suite())
        .def_readwrite("defocus", &EMAN::EMAN1Ctf::defocus)
        .def_readwrite("bfactor", &EMAN::EMAN1Ctf::bfactor)
        .def_readwrite("amplitude", &EMAN::EMAN1Ctf::amplitude)
        .def_readwrite("ampcont", &EMAN::EMAN1Ctf::ampcont)
        .def_readwrite("noise1", &EMAN::EMAN1Ctf::noise1)
        .def_readwrite("noise2", &EMAN::EMAN1Ctf::noise2)
        .def_readwrite("noise3", &EMAN::EMAN1Ctf::noise3)
        .def_readwrite("noise4", &EMAN::EMAN1Ctf::noise4)
        .def_readwrite("voltage", &EMAN::EMAN1Ctf::voltage)
        .def_readwrite("cs", &EMAN::EMAN1Ctf::cs)
        .def_readwrite("apix", &EMAN::EMAN1Ctf::apix)
        .def("compute_1d", (std::vector<float,std::allocator<float> > (EMAN::EMAN1Ctf::*)(int, float, EMAN::Ctf::CtfType, EMAN::XYData*) )&EMAN::EMAN1Ctf::compute_1d, (std::vector<float,std::allocator<float> > (EMAN_EMAN1Ctf_Wrapper::*)(int, float, EMAN::Ctf::CtfType, EMAN::XYData*))&EMAN_EMAN1Ctf_Wrapper::default_compute_1d_3)
        .def("compute_1d", (std::vector<float,std::allocator<float> > (EMAN_EMAN1Ctf_Wrapper::*)(int, float, EMAN::Ctf::CtfType))&EMAN_EMAN1Ctf_Wrapper::default_compute_1d_2)
        .def("compute_1d_fromimage", (std::vector<float,std::allocator<float> > (EMAN::EMAN1Ctf::*)(int, float, EMAN::EMData*) )&EMAN::EMAN1Ctf::compute_1d_fromimage, (std::vector<float,std::allocator<float> > (EMAN_EMAN1Ctf_Wrapper::*)(int, float, EMAN::EMData*))&EMAN_EMAN1Ctf_Wrapper::default_compute_1d_fromimage)
        .def("compute_2d_real", (void (EMAN::EMAN1Ctf::*)(EMAN::EMData*, EMAN::Ctf::CtfType, EMAN::XYData*) )&EMAN::EMAN1Ctf::compute_2d_real, (void (EMAN_EMAN1Ctf_Wrapper::*)(EMAN::EMData*, EMAN::Ctf::CtfType, EMAN::XYData*))&EMAN_EMAN1Ctf_Wrapper::default_compute_2d_real_3)
        .def("compute_2d_real", (void (EMAN_EMAN1Ctf_Wrapper::*)(EMAN::EMData*, EMAN::Ctf::CtfType))&EMAN_EMAN1Ctf_Wrapper::default_compute_2d_real_2)
        .def("compute_2d_complex", (void (EMAN::EMAN1Ctf::*)(EMAN::EMData*, EMAN::Ctf::CtfType, EMAN::XYData*) )&EMAN::EMAN1Ctf::compute_2d_complex, (void (EMAN_EMAN1Ctf_Wrapper::*)(EMAN::EMData*, EMAN::Ctf::CtfType, EMAN::XYData*))&EMAN_EMAN1Ctf_Wrapper::default_compute_2d_complex_3)
        .def("compute_2d_complex", (void (EMAN_EMAN1Ctf_Wrapper::*)(EMAN::EMData*, EMAN::Ctf::CtfType))&EMAN_EMAN1Ctf_Wrapper::default_compute_2d_complex_2)
        .def("from_string", (int (EMAN::EMAN1Ctf::*)(const std::string&) )&EMAN::EMAN1Ctf::from_string, (int (EMAN_EMAN1Ctf_Wrapper::*)(const std::string&))&EMAN_EMAN1Ctf_Wrapper::default_from_string)
        .def("to_string", (std::string (EMAN::EMAN1Ctf::*)() const)&EMAN::EMAN1Ctf::to_string, (std::string (EMAN_EMAN1Ctf_Wrapper::*)() const)&EMAN_EMAN1Ctf_Wrapper::default_to_string)
        .def("from_dict", (void (EMAN::EMAN1Ctf::*)(const EMAN::Dict&) )&EMAN::EMAN1Ctf::from_dict, (void (EMAN_EMAN1Ctf_Wrapper::*)(const EMAN::Dict&))&EMAN_EMAN1Ctf_Wrapper::default_from_dict)
        .def("to_dict", (EMAN::Dict (EMAN::EMAN1Ctf::*)() const)&EMAN::EMAN1Ctf::to_dict, (EMAN::Dict (EMAN_EMAN1Ctf_Wrapper::*)() const)&EMAN_EMAN1Ctf_Wrapper::default_to_dict)
        .def("from_vector", (void (EMAN::EMAN1Ctf::*)(const std::vector<float,std::allocator<float> >&) )&EMAN::EMAN1Ctf::from_vector, (void (EMAN_EMAN1Ctf_Wrapper::*)(const std::vector<float,std::allocator<float> >&))&EMAN_EMAN1Ctf_Wrapper::default_from_vector)
        .def("to_vector", (std::vector<float,std::allocator<float> > (EMAN::EMAN1Ctf::*)() const)&EMAN::EMAN1Ctf::to_vector, (std::vector<float,std::allocator<float> > (EMAN_EMAN1Ctf_Wrapper::*)() const)&EMAN_EMAN1Ctf_Wrapper::default_to_vector)
        .def("copy_from", (void (EMAN::EMAN1Ctf::*)(const EMAN::Ctf*) )&EMAN::EMAN1Ctf::copy_from, (void (EMAN_EMAN1Ctf_Wrapper::*)(const EMAN::Ctf*))&EMAN_EMAN1Ctf_Wrapper::default_copy_from)
        .def("equal", (bool (EMAN::EMAN1Ctf::*)(const EMAN::Ctf*) const)&EMAN::EMAN1Ctf::equal, (bool (EMAN_EMAN1Ctf_Wrapper::*)(const EMAN::Ctf*) const)&EMAN_EMAN1Ctf_Wrapper::default_equal)
        .def("zero", (float (EMAN::EMAN1Ctf::*)(int) const)&EMAN::EMAN1Ctf::zero, (float (EMAN_EMAN1Ctf_Wrapper::*)(int) const)&EMAN_EMAN1Ctf_Wrapper::default_zero)
    ;

    class_< EMAN::EMAN2Ctf, bases< EMAN::Ctf > , EMAN_EMAN2Ctf_Wrapper >("EMAN2Ctf",
    		"EMAN2Ctf is the default CTF model used in EMAN2",
    		init<  >())
        .def(init< const EMAN::EMAN2Ctf& >())
        .def(init<const std::vector<float>& >())
        .def_pickle(EMAN2Ctf_pickle_suite())
        .def_readwrite("defocus", &EMAN::EMAN2Ctf::defocus)
        .def_readwrite("dfdiff", &EMAN::EMAN2Ctf::dfdiff)
        .def_readwrite("dfang", &EMAN::EMAN2Ctf::dfang)
        .def_readwrite("bfactor", &EMAN::EMAN2Ctf::bfactor)
        .def_readwrite("ampcont", &EMAN::EMAN2Ctf::ampcont)
        .def_readwrite("voltage", &EMAN::EMAN2Ctf::voltage)
        .def_readwrite("cs", &EMAN::EMAN2Ctf::cs)
        .def_readwrite("apix", &EMAN::EMAN2Ctf::apix)
		.def_readwrite("dsbg", &EMAN::EMAN2Ctf::dsbg)
		.add_property("snr", &EMAN::EMAN2Ctf::get_snr, &EMAN::EMAN2Ctf::set_snr)
		.add_property("background", &EMAN::EMAN2Ctf::get_background, &EMAN::EMAN2Ctf::set_background)
        .def("compute_1d", (std::vector<float,std::allocator<float> > (EMAN::EMAN2Ctf::*)(int, float, EMAN::Ctf::CtfType, EMAN::XYData*) )&EMAN::EMAN2Ctf::compute_1d, (std::vector<float,std::allocator<float> > (EMAN_EMAN2Ctf_Wrapper::*)(int, float, EMAN::Ctf::CtfType, EMAN::XYData*))&EMAN_EMAN2Ctf_Wrapper::default_compute_1d_3)
        .def("compute_1d", (std::vector<float,std::allocator<float> > (EMAN_EMAN2Ctf_Wrapper::*)(int, float, EMAN::Ctf::CtfType))&EMAN_EMAN2Ctf_Wrapper::default_compute_1d_2)
        .def("compute_1d_fromimage", (std::vector<float,std::allocator<float> > (EMAN::EMAN2Ctf::*)(int, float, EMAN::EMData*) )&EMAN::EMAN2Ctf::compute_1d, (std::vector<float,std::allocator<float> > (EMAN_EMAN2Ctf_Wrapper::*)(int, float, EMAN::EMData*))&EMAN_EMAN2Ctf_Wrapper::default_compute_1d_fromimage)
        .def("compute_2d_real", (void (EMAN::EMAN2Ctf::*)(EMAN::EMData*, EMAN::Ctf::CtfType, EMAN::XYData*) )&EMAN::EMAN2Ctf::compute_2d_real, (void (EMAN_EMAN2Ctf_Wrapper::*)(EMAN::EMData*, EMAN::Ctf::CtfType, EMAN::XYData*))&EMAN_EMAN2Ctf_Wrapper::default_compute_2d_real_3)
        .def("compute_2d_real", (void (EMAN_EMAN2Ctf_Wrapper::*)(EMAN::EMData*, EMAN::Ctf::CtfType))&EMAN_EMAN2Ctf_Wrapper::default_compute_2d_real_2)
        .def("compute_2d_complex", (void (EMAN::EMAN2Ctf::*)(EMAN::EMData*, EMAN::Ctf::CtfType, EMAN::XYData*) )&EMAN::EMAN2Ctf::compute_2d_complex, (void (EMAN_EMAN2Ctf_Wrapper::*)(EMAN::EMData*, EMAN::Ctf::CtfType, EMAN::XYData*))&EMAN_EMAN2Ctf_Wrapper::default_compute_2d_complex_3)
        .def("compute_2d_complex", (void (EMAN_EMAN2Ctf_Wrapper::*)(EMAN::EMData*, EMAN::Ctf::CtfType))&EMAN_EMAN2Ctf_Wrapper::default_compute_2d_complex_2)
        .def("from_string", (int (EMAN::EMAN2Ctf::*)(const std::string&) )&EMAN::EMAN2Ctf::from_string, (int (EMAN_EMAN2Ctf_Wrapper::*)(const std::string&))&EMAN_EMAN2Ctf_Wrapper::default_from_string)
        .def("to_string", (std::string (EMAN::EMAN2Ctf::*)() const)&EMAN::EMAN2Ctf::to_string, (std::string (EMAN_EMAN2Ctf_Wrapper::*)() const)&EMAN_EMAN2Ctf_Wrapper::default_to_string)
        .def("from_dict", (void (EMAN::EMAN2Ctf::*)(const EMAN::Dict&) )&EMAN::EMAN2Ctf::from_dict, (void (EMAN_EMAN2Ctf_Wrapper::*)(const EMAN::Dict&))&EMAN_EMAN2Ctf_Wrapper::default_from_dict)
        .def("to_dict", (EMAN::Dict (EMAN::EMAN2Ctf::*)() const)&EMAN::EMAN2Ctf::to_dict, (EMAN::Dict (EMAN_EMAN2Ctf_Wrapper::*)() const)&EMAN_EMAN2Ctf_Wrapper::default_to_dict)
        .def("from_vector", (void (EMAN::EMAN2Ctf::*)(const std::vector<float,std::allocator<float> >&) )&EMAN::EMAN2Ctf::from_vector, (void (EMAN_EMAN2Ctf_Wrapper::*)(const std::vector<float,std::allocator<float> >&))&EMAN_EMAN2Ctf_Wrapper::default_from_vector)
        .def("to_vector", (std::vector<float,std::allocator<float> > (EMAN::EMAN2Ctf::*)() const)&EMAN::EMAN2Ctf::to_vector, (std::vector<float,std::allocator<float> > (EMAN_EMAN2Ctf_Wrapper::*)() const)&EMAN_EMAN2Ctf_Wrapper::default_to_vector)
        .def("copy_from", (void (EMAN::EMAN2Ctf::*)(const EMAN::Ctf*) )&EMAN::EMAN2Ctf::copy_from, (void (EMAN_EMAN2Ctf_Wrapper::*)(const EMAN::Ctf*))&EMAN_EMAN2Ctf_Wrapper::default_copy_from)
        .def("equal", (bool (EMAN::EMAN2Ctf::*)(const EMAN::Ctf*) const)&EMAN::EMAN2Ctf::equal, (bool (EMAN_EMAN2Ctf_Wrapper::*)(const EMAN::Ctf*) const)&EMAN_EMAN2Ctf_Wrapper::default_equal)
        .def("zero", (float (EMAN::EMAN2Ctf::*)(int) const)&EMAN::EMAN2Ctf::zero, (float (EMAN_EMAN2Ctf_Wrapper::*)(int) const)&EMAN_EMAN2Ctf_Wrapper::default_zero)
    ;

}

