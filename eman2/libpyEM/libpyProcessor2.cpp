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
#include <emdata.h>
#include <emobject.h>
#include <processor.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct EMAN_Processor_Wrapper: EMAN::Processor
{
    EMAN_Processor_Wrapper(PyObject* py_self_, const EMAN::Processor& p0):
        EMAN::Processor(p0), py_self(py_self_) {}

    EMAN_Processor_Wrapper(PyObject* py_self_):
        EMAN::Processor(), py_self(py_self_) {}

    void process_inplace(EMAN::EMData* p0) {
        call_method< void >(py_self, "process_inplace", p0);
    }

    EMAN::EMData* process(const EMAN::EMData* const p0) {
        return call_method< EMAN::EMData* >(py_self, "process", p0);
    }

    EMAN::EMData* default_process(const EMAN::EMData* const p0) {
        return EMAN::Processor::process(p0);
    }

    void process_list_inplace(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >& p0) {
        call_method< void >(py_self, "process_list_inplace", p0);
    }

    void default_process_list_inplace(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> >& p0) {
        EMAN::Processor::process_list_inplace(p0);
    }

    std::string get_name() const {
        return call_method< std::string >(py_self, "get_name");
    }

    EMAN::Dict get_params() const {
        return call_method< EMAN::Dict >(py_self, "get_params");
    }

    EMAN::Dict default_get_params() const {
        return EMAN::Processor::get_params();
    }

    void set_params(const EMAN::Dict& p0) {
        call_method< void >(py_self, "set_params", p0);
    }

    void default_set_params(const EMAN::Dict& p0) {
        EMAN::Processor::set_params(p0);
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(py_self, "get_param_types");
    }

    EMAN::TypeDict default_get_param_types() const {
        return EMAN::Processor::get_param_types();
    }

    std::string get_desc() const {
        return call_method< std::string >(py_self, "get_desc");
    }

    PyObject* py_self;
};


}// namespace


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyProcessor2)
{
    scope* EMAN_Processor_scope = new scope(
    class_< EMAN::Processor, boost::noncopyable, EMAN_Processor_Wrapper >("Processor", init<  >())
        .def("process_inplace", pure_virtual(&EMAN::Processor::process_inplace))
        .def("process", &EMAN::Processor::process, &EMAN_Processor_Wrapper::default_process, return_value_policy< manage_new_object >())
        .def("process_list_inplace", &EMAN::Processor::process_list_inplace, &EMAN_Processor_Wrapper::default_process_list_inplace)
        .def("get_name", pure_virtual(&EMAN::Processor::get_name))
        .def("get_params", &EMAN::Processor::get_params, &EMAN_Processor_Wrapper::default_get_params)
        .def("set_params", &EMAN::Processor::set_params, &EMAN_Processor_Wrapper::default_set_params)
        .def("get_param_types", &EMAN::Processor::get_param_types, &EMAN_Processor_Wrapper::default_get_param_types)
        .def("get_desc", pure_virtual(&EMAN::Processor::get_desc))
        .def("get_group_desc", &EMAN::Processor::get_group_desc)
        .def("EMFourierFilterInPlace", &EMAN::Processor::EMFourierFilterInPlace)
        .def("EMFourierFilter", &EMAN::Processor::EMFourierFilter, return_value_policy< manage_new_object >())
        .staticmethod("EMFourierFilterInPlace")
        .staticmethod("get_group_desc")
        .staticmethod("EMFourierFilter")
    );

    enum_< EMAN::Processor::fourier_filter_types >("fourier_filter_types")
        .value("GAUSS_HIGH_PASS", EMAN::Processor::GAUSS_HIGH_PASS)
        .value("SHIFT", EMAN::Processor::SHIFT)
        .value("GAUSS_INVERSE", EMAN::Processor::GAUSS_INVERSE)
        .value("BUTTERWORTH_HOMOMORPHIC", EMAN::Processor::BUTTERWORTH_HOMOMORPHIC)
        .value("RADIAL_TABLE", EMAN::Processor::RADIAL_TABLE)
        .value("TANH_HOMOMORPHIC", EMAN::Processor::TANH_HOMOMORPHIC)
        .value("KAISER_SINH_INVERSE", EMAN::Processor::KAISER_SINH_INVERSE)
        .value("TANH_LOW_PASS", EMAN::Processor::TANH_LOW_PASS)
        .value("CTF_", EMAN::Processor::CTF_)
        .value("KAISER_SINH", EMAN::Processor::KAISER_SINH)
        .value("KAISER_I0_INVERSE", EMAN::Processor::KAISER_I0_INVERSE)
        .value("TOP_HAT_BAND_PASS", EMAN::Processor::TOP_HAT_BAND_PASS)
        .value("KAISER_I0", EMAN::Processor::KAISER_I0)
        .value("TOP_HAT_HIGH_PASS", EMAN::Processor::TOP_HAT_HIGH_PASS)
        .value("TANH_HIGH_PASS", EMAN::Processor::TANH_HIGH_PASS)
        .value("BUTTERWORTH_LOW_PASS", EMAN::Processor::BUTTERWORTH_LOW_PASS)
        .value("TOP_HAT_LOW_PASS", EMAN::Processor::TOP_HAT_LOW_PASS)
        .value("GAUSS_HOMOMORPHIC", EMAN::Processor::GAUSS_HOMOMORPHIC)
        .value("GAUSS_LOW_PASS", EMAN::Processor::GAUSS_LOW_PASS)
        .value("GAUSS_BAND_PASS", EMAN::Processor::GAUSS_BAND_PASS)
        .value("TOP_HOMOMORPHIC", EMAN::Processor::TOP_HOMOMORPHIC)
        .value("BUTTERWORTH_HIGH_PASS", EMAN::Processor::BUTTERWORTH_HIGH_PASS)
        .value("TANH_BAND_PASS", EMAN::Processor::TANH_BAND_PASS)
    ;

    delete EMAN_Processor_scope;

    def("dump_processors", &EMAN::dump_processors);
    def("dump_processors_list", &EMAN::dump_processors_list);
    def("multi_processors", &EMAN::multi_processors);
    def("group_processors", &EMAN::group_processors);
    class_< EMAN::Factory<EMAN::Processor>, boost::noncopyable >("Processors", no_init)
        .def("get", (EMAN::Processor* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&))&EMAN::Factory<EMAN::Processor>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Processor* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&, const EMAN::Dict&))&EMAN::Factory<EMAN::Processor>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Processor>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

#ifdef SPARX_USING_CUDA	
// Class to wrap MPI CUDA kmeans code
    class_< EMAN::MPICUDA_kmeans, boost::noncopyable >("MPICUDA_kmeans", init<>())
      .def("setup", &EMAN::MPICUDA_kmeans::setup)
      .def("append_flat_image", &EMAN::MPICUDA_kmeans::append_flat_image)
      .def("init_mem", &EMAN::MPICUDA_kmeans::init_mem)
      .def("compute_im2", &EMAN::MPICUDA_kmeans::compute_im2)
      .def("random_ASG", &EMAN::MPICUDA_kmeans::random_ASG)
      .def("get_ASG", &EMAN::MPICUDA_kmeans::get_ASG)
      .def("get_asg", &EMAN::MPICUDA_kmeans::get_asg)
      .def("compute_NC", &EMAN::MPICUDA_kmeans::compute_NC)
      .def("get_NC", &EMAN::MPICUDA_kmeans::get_NC)
      .def("set_ASG", &EMAN::MPICUDA_kmeans::set_ASG)
      .def("set_NC", &EMAN::MPICUDA_kmeans::set_NC)
      .def("compute_AVE", &EMAN::MPICUDA_kmeans::compute_AVE)
      .def("get_AVE", &EMAN::MPICUDA_kmeans::get_AVE)
      .def("set_AVE", &EMAN::MPICUDA_kmeans::set_AVE)
      .def("get_ct_im_mv", &EMAN::MPICUDA_kmeans::get_ct_im_mv)
      .def("set_T", &EMAN::MPICUDA_kmeans::set_T)
      .def("get_T", &EMAN::MPICUDA_kmeans::get_T)
      .def("one_iter", &EMAN::MPICUDA_kmeans::one_iter)
      //.def("AVE_to_host", &EMAN::MPICUDA_kmeans::AVE_to_host)
      //.def("one_iter_SSE", &EMAN::MPICUDA_kmeans::one_iter_SSE)
      .def("one_iter_SA", &EMAN::MPICUDA_kmeans::one_iter_SA)
      .def("compute_ji", &EMAN::MPICUDA_kmeans::compute_ji)
      .def("compute_criterion", &EMAN::MPICUDA_kmeans::compute_criterion)      
      .def("shutdown", &EMAN::MPICUDA_kmeans::shutdown)
    ;

#endif //SPARX_USING_CUDA

}

