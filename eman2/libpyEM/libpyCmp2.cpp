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
#include <cmp.h>
#include <emdata.h>
#include <emobject.h>
#include <log.h>
#include <transform.h>
#include <xydata.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct EMAN_Cmp_Wrapper: EMAN::Cmp
{
    EMAN_Cmp_Wrapper(PyObject* py_self_, const EMAN::Cmp& p0):
        EMAN::Cmp(p0), py_self(py_self_) {}

    EMAN_Cmp_Wrapper(PyObject* py_self_):
        EMAN::Cmp(), py_self(py_self_) {}

    float cmp(EMAN::EMData* p0, EMAN::EMData* p1) const {
        return call_method< float >(py_self, "cmp", p0, p1);
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
        return EMAN::Cmp::get_params();
    }

    void set_params(const EMAN::Dict& p0) {
        call_method< void >(py_self, "set_params", p0);
    }

    void default_set_params(const EMAN::Dict& p0) {
        EMAN::Cmp::set_params(p0);
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(py_self, "get_param_types");
    }

    PyObject* py_self;
};

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_Log_end_overloads_1_3, end, 1, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(EMAN_XYData_get_yatx_overloads_1_2, get_yatx, 1, 2)

}// namespace


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyCmp2)
{
    def("dump_cmps", &EMAN::dump_cmps);
    def("dump_cmps_list", &EMAN::dump_cmps_list);
    class_< EMAN::Factory<EMAN::Cmp>, boost::noncopyable >("Cmps", no_init)
        .def("get", (EMAN::Cmp* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&))&EMAN::Factory<EMAN::Cmp>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Cmp* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&, const EMAN::Dict&))&EMAN::Factory<EMAN::Cmp>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Cmp>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;

    class_< EMAN::Cmp, boost::noncopyable, EMAN_Cmp_Wrapper >("__Cmp", init<  >())
        .def("cmp", pure_virtual(&EMAN::Cmp::cmp))
        .def("get_name", pure_virtual(&EMAN::Cmp::get_name))
        .def("get_desc", pure_virtual(&EMAN::Cmp::get_desc))
        .def("get_params", &EMAN::Cmp::get_params, &EMAN_Cmp_Wrapper::default_get_params)
        .def("set_params", &EMAN::Cmp::set_params, &EMAN_Cmp_Wrapper::default_set_params)
        .def("get_param_types", pure_virtual(&EMAN::Cmp::get_param_types))
    ;

    scope* EMAN_Log_scope = new scope(
    class_< EMAN::Log, boost::noncopyable >("Log",
    		"Log defines a way to output logging information.\n"
    		"1) The logs can either go to standard output (default), or go to a user-given file.\n"
    		"2) 4 verbose log levels are defined. by default, ERROR_LOG is used.\n"
    		"3) Typical usage:\n"
    		"Log.logger().set_level(Log.WARNING_LOG)\n",
    		no_init)
        .def("logger", &EMAN::Log::logger, return_value_policy< reference_existing_object >())
        .def("begin", &EMAN::Log::begin)
        .def("end", &EMAN::Log::end, EMAN_Log_end_overloads_1_3())
        .def("set_level", &EMAN::Log::set_level)
        .def("set_logfile", &EMAN::Log::set_logfile)
        .def("loc", &EMAN::Log::loc)
        .staticmethod("logger")
    );

    enum_< EMAN::Log::LogLevel >("LogLevel")
        .value("ERROR_LOG", EMAN::Log::ERROR_LOG)
        .value("VARIABLE_LOG", EMAN::Log::VARIABLE_LOG)
        .value("WARNING_LOG", EMAN::Log::WARNING_LOG)
        .value("DEBUG_LOG", EMAN::Log::DEBUG_LOG)
    ;

    delete EMAN_Log_scope;

    scope* EMAN_XYData_scope = new scope(
    class_< EMAN::XYData >("XYData", "XYData defines a 1D (x,y) data set.", init<  >())
        .enable_pickling()
        .def(init< const EMAN::XYData& >())
        .def("read_file", &EMAN::XYData::read_file)
        .def("write_file", &EMAN::XYData::write_file)
        .def("calc_correlation", &EMAN::XYData::calc_correlation)
        .def("update", &EMAN::XYData::update)
        .def("get_yatx", &EMAN::XYData::get_yatx, EMAN_XYData_get_yatx_overloads_1_2())
        .def("get_yatx_smooth", &EMAN::XYData::get_yatx_smooth)
        .def("insort", &EMAN::XYData::insort)
        .def("dedupx", &EMAN::XYData::dedupx)
        .def("get_x", &EMAN::XYData::get_x)
        .def("set_x", &EMAN::XYData::set_x)
        .def("get_y", &EMAN::XYData::get_y)
        .def("set_y", &EMAN::XYData::set_y)
        .def("get_size", &EMAN::XYData::get_size)
        .def("get_miny", &EMAN::XYData::get_miny)
        .def("get_maxy", &EMAN::XYData::get_maxy)
        .def("is_validx", &EMAN::XYData::is_validx)
        .def("set_xy_list", &EMAN::XYData::set_xy_list)
        .def("set_size", &EMAN::XYData::set_size)
        .def("get_xlist", &EMAN::XYData::get_xlist)
        .def("get_ylist", &EMAN::XYData::get_ylist)
        .def("__getstate__", &EMAN::XYData::get_state)
        .def("__setstate__", &EMAN::XYData::set_state)
    );

    class_< EMAN::XYData::Pair >("Pair", "a pair of float x and y", init< const EMAN::XYData::Pair& >())
        .def(init< float, float >())
        .def_readwrite("x", &EMAN::XYData::Pair::x)
        .def_readwrite("y", &EMAN::XYData::Pair::y)
        .def( self < self )
    ;

    delete EMAN_XYData_scope;

}

