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
#include <reconstructor.h>

// Using =======================================================================
using namespace boost::python;

// debug
#include <iostream>
using std::cout;
using std::endl;

// Declarations ================================================================
using namespace EMAN;

	int reconstructor_insert_slice2(Reconstructor &self, const EMData* slice, const Transform& euler) {
		int ret;
		Py_BEGIN_ALLOW_THREADS
//		printf("wrapper1\n");
		ret=self.insert_slice(slice,euler);
//		ret=call_method< int >(py_self, "insert_slice", slice,euler,1.0f);
		Py_END_ALLOW_THREADS
		return ret;
	}

 	int reconstructor_insert_slice3(Reconstructor &self, const EMData* slice, const Transform& euler,float weight) {
		int ret;
		Py_BEGIN_ALLOW_THREADS
//		printf("wrapper2 %p %p %f\n",&self,slice,weight);
		ret=self.insert_slice(slice,euler,weight);
// 		ret=call_method< int >(py_self, "insert_slice", slice,euler,weight);
		Py_END_ALLOW_THREADS
		return ret;
 	}


struct EMAN_Reconstructor_Wrapper: EMAN::Reconstructor
{
//     EMAN_Reconstructor_Wrapper(PyObject* py_self_, const EMAN::Reconstructor& p0):
//         EMAN::Reconstructor(p0), py_self(py_self_) {}

    EMAN_Reconstructor_Wrapper(PyObject* py_self_):
        EMAN::Reconstructor(), py_self(py_self_) {}

    void setup() {
        call_method< void >(py_self, "setup");
    }

    void setup_seed(const EMAN::EMData* seed,float seed_weight) {
        call_method< void >(py_self, "setup_seed",seed,seed_weight);
    }

    void clear() {
        call_method< void >(py_self, "clear");
    }
    
 	int insert_slice(const EMAN::EMData* const slice, const EMAN::Transform& euler) {
		int ret;
		Py_BEGIN_ALLOW_THREADS
		ret=call_method< int >(py_self, "insert_slice", slice,euler,1.0);
		Py_END_ALLOW_THREADS
		return ret;
	}

 	int insert_slice(const EMAN::EMData* const slice, const EMAN::Transform& euler,float weight) {
		int ret;
		Py_BEGIN_ALLOW_THREADS
 		ret=call_method< int >(py_self, "insert_slice", slice,euler,weight);
		Py_END_ALLOW_THREADS
		return ret;
 	}
 	
// 	int insert_slice2( const EMData* slice, const Transform& euler) {
// 		int ret;
// 		Py_BEGIN_ALLOW_THREADS
// 		printf("wrapper1\n");
// //		ret=py_self->insert_slice(slice,euler);
// 		ret=call_method< int >(py_self, "insert_slice", slice,euler,1.0f);
// 		Py_END_ALLOW_THREADS
// 		return ret;
// 	}
// 
//  	int insert_slice3( const EMData* slice, const Transform& euler,float weight) {
// 		int ret;
// 		Py_BEGIN_ALLOW_THREADS
// 		printf("wrapper2 %p %p %f\n",py_self,slice,weight);
// //		ret=py_self->insert_slice(slice,euler,weight);
//  		ret=call_method< int >(py_self, "insert_slice", slice,euler,weight);
// 		Py_END_ALLOW_THREADS
// 		return ret;
//  	}

    EMAN::EMData* finish(bool doift) {
        return call_method< EMAN::EMData* >(py_self, "finish", doift);
    }

    std::string get_name() const {
        return call_method< std::string >(py_self, "get_name");
    }

    std::string get_desc() const {
        return call_method< std::string >(py_self, "get_desc");
    }

    EMAN::Dict get_params() const {
		printf("call goes here!\n");
        return call_method< EMAN::Dict >(py_self, "get_params");
    }
/*
	void print_params() const {
        call_method< void >(py_self, "print_params");
	}*/

    EMAN::Dict default_get_params() const {
        return EMAN::Reconstructor::get_params();
    }

    void set_params(const EMAN::Dict& p0) {
        call_method< void >(py_self, "set_params", p0);
    }

    void default_set_params(const EMAN::Dict& p0) {
        EMAN::Reconstructor::set_params(p0);
    }

    EMAN::TypeDict get_param_types() const {
        return call_method< EMAN::TypeDict >(py_self, "get_param_types");
    }

    PyObject* py_self;
};


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyReconstructor2)
{
    def("dump_reconstructors", &EMAN::dump_reconstructors);
    def("dump_reconstructors_list", &EMAN::dump_reconstructors_list);
    class_< EMAN::Reconstructor, boost::noncopyable, EMAN_Reconstructor_Wrapper >("__Reconstructor", init<  >())
        .def("setup", pure_virtual(&EMAN::Reconstructor::setup))
        .def("clear", pure_virtual(&EMAN::Reconstructor::clear))
		.def("setup_seed", (int (EMAN::Reconstructor::*)(const EMAN::EMData*, const float))&EMAN::Reconstructor::setup_seed)
//		.def("insert_slice", (int (EMAN::Reconstructor::*)(const EMAN::EMData* const, const EMAN::Transform&, const float))&EMAN::Reconstructor::insert_slice)
//		.def("insert_slice", (int (EMAN::Reconstructor::*)(const EMAN::EMData* const, const EMAN::Transform&))&EMAN::Reconstructor::insert_slice)
// 		.def("insert_slice", (int (EMAN::Reconstructor::*)(const EMAN::EMData* const, const EMAN::Transform&, const float))&EMAN_Reconstructor_Wrapper::insert_slice3)
// 		.def("insert_slice", (int (EMAN::Reconstructor::*)(const EMAN::EMData* const, const EMAN::Transform&))&EMAN_Reconstructor_Wrapper::insert_slice2)
		.def("insert_slice", &reconstructor_insert_slice3)
		.def("insert_slice", &reconstructor_insert_slice2)
		.def("determine_slice_agreement", (int (EMAN::Reconstructor::*)(EMAN::EMData* , const EMAN::Transform&, const float, bool))&EMAN::Reconstructor::determine_slice_agreement)
        .def("preprocess_slice", (EMAN::EMData* (EMAN::Reconstructor::*)(const EMAN::EMData* const, const EMAN::Transform&))&EMAN::Reconstructor::preprocess_slice, return_value_policy< manage_new_object >())
        .def("finish", (EMAN::EMData* (EMAN::Reconstructor::*)(bool))&EMAN::Reconstructor::finish, return_value_policy< manage_new_object >())
        .def("get_name", pure_virtual(&EMAN::Reconstructor::get_name))
        .def("get_desc", pure_virtual(&EMAN::Reconstructor::get_desc))
// 		.def("get_emdata", (&EMAN::Reconstructor::get_emdata),  return_internal_reference< 1 >())
        //.def("get_params", &EMAN::Reconstructor::get_params, &EMAN_Reconstructor_Wrapper::default_get_params)
		.def("get_params", &EMAN::Reconstructor::get_params)
		.def("insert_params", &EMAN::Reconstructor::insert_params)
        .def("set_params", &EMAN::Reconstructor::set_params, &EMAN_Reconstructor_Wrapper::default_set_params)
        .def("set_param", &EMAN::Reconstructor::set_param)
		.def("print_params",  &EMAN::Reconstructor::print_params) // Why is this different to set_params and get_params? Why is the wrapper needed? d.woolford May 2007
        .def("get_param_types", pure_virtual(&EMAN::Reconstructor::get_param_types))
    ;

    class_< EMAN::Factory<EMAN::Reconstructor>, boost::noncopyable >("Reconstructors", no_init)
        .def("get", (EMAN::Reconstructor* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&))&EMAN::Factory<EMAN::Reconstructor>::get, return_value_policy< manage_new_object >())
        .def("get", (EMAN::Reconstructor* (*)(const std::basic_string<char,std::char_traits<char>,std::allocator<char> >&, const EMAN::Dict&))&EMAN::Factory<EMAN::Reconstructor>::get, return_value_policy< manage_new_object >())
        .def("get_list", &EMAN::Factory<EMAN::Reconstructor>::get_list)
        .staticmethod("get_list")
        .staticmethod("get")
    ;


    class_< EMAN::file_store >( "file_store", init< const string&, int, int, bool >() )
        .def( "add_image", &EMAN::file_store::add_image )
        .def( "get_image", &EMAN::file_store::get_image )
        .def( "restart",   &EMAN::file_store::restart )
    ;


    class_< EMAN::newfile_store >( "newfile_store", init< const string&, int, bool >() )
        .def( "add_image", &EMAN::newfile_store::add_image )
        .def( "add_tovol", &EMAN::newfile_store::add_tovol )
	.def( "restart",   &EMAN::newfile_store::restart )
        .def( "read",      &EMAN::newfile_store::read )
    ;


}

