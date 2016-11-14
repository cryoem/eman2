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
#include <analyzer.h>
#include <emdata.h>
#include <emobject.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct EMAN_Analyzer_Wrapper: EMAN::Analyzer
{
    EMAN_Analyzer_Wrapper(PyObject* py_self_, const EMAN::Analyzer& p0):
        EMAN::Analyzer(p0), py_self(py_self_) {}

    EMAN_Analyzer_Wrapper(PyObject* py_self_):
        EMAN::Analyzer(), py_self(py_self_) {}

    int insert_image(EMAN::EMData* p0) {
        return call_method< int >(py_self, "insert_image", p0);
    }

    int insert_images_list(std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> > p0) {
        return call_method< int >(py_self, "insert_images_list", p0);
    }

    std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> > analyze() {
        return call_method< std::vector<EMAN::EMData*,std::allocator<EMAN::EMData*> > >(py_self, "analyze");
    }

    std::string get_name() const {
        return call_method< std::string >(py_self, "get_name");
    }

    std::string get_desc() const {
        return call_method< std::string >(py_self, "get_desc");
    }

    void set_params(const EMAN::Dict& p0) {
        call_method< void >(py_self, "set_params", p0);
    }

    void default_set_params(const EMAN::Dict& p0) {
        EMAN::Analyzer::set_params(p0);
    }

    EMAN::Dict get_params() const {
        return call_method< EMAN::Dict >(py_self, "get_params");
    }

    EMAN::Dict default_get_params() const {
        return EMAN::Analyzer::get_params();
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(py_self, "get_param_types");
    }

    PyObject* py_self;
};


}// namespace


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyAnalyzer2)
{
    def("dump_analyzers", &EMAN::dump_analyzers);
    def("dump_analyzers_list", &EMAN::dump_analyzers_list);
    class_< EMAN::Analyzer, boost::noncopyable, EMAN_Analyzer_Wrapper >("__Analyzer", init<  >())
        .def("insert_image", pure_virtual(&EMAN::Analyzer::insert_image))
        .def("insert_images_list", pure_virtual(&EMAN::Analyzer::insert_images_list))
        .def("analyze", pure_virtual(&EMAN::Analyzer::analyze))
        .def("get_name", pure_virtual(&EMAN::Analyzer::get_name))
        .def("get_desc", pure_virtual(&EMAN::Analyzer::get_desc))
        .def("set_params", &EMAN::Analyzer::set_params, &EMAN_Analyzer_Wrapper::default_set_params)
        .def("get_params", &EMAN::Analyzer::get_params, &EMAN_Analyzer_Wrapper::default_get_params)
        .def("get_param_types", pure_virtual(&EMAN::Analyzer::get_param_types))
    ;

    class_< EMAN::Factory<EMAN::Analyzer>, boost::noncopyable >("Analyzers", no_init)
        .def("get", (EMAN::Analyzer* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&))&EMAN::Factory<EMAN::Analyzer>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Analyzer* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&, const EMAN::Dict&))&EMAN::Factory<EMAN::Analyzer>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Analyzer>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

}

